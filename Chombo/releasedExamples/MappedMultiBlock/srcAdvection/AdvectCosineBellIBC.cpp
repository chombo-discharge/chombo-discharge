#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AdvectCosineBellIBC.H"
#include "BoxIterator.H"
#include "CubedSphere2DPanelCS.H"
#include "CylinderEquiangularBlockCS.H"
#include "CylinderSpokesBlockCS.H"
#include "XPointBlockCoordSys.H"
#include "CubedSphereSolidBlockCS.H"
#include "NamespaceHeader.H"

// Null constructor
AdvectCosineBellIBC::AdvectCosineBellIBC()
{
  // m_haveAdvVel = false;
  // CH_assert(false);
  // m_params_are_set = false;
}

// Constructor which defines parameters used by Fortran routines.
// Mapping changes nothing.
AdvectCosineBellIBC::AdvectCosineBellIBC(
                                     const Real&     a_ambientDensity,
                                     const Real&     a_deltaDensity,
                                     const RealVect& a_center,
                                     const Real&     a_size,
                                     const Real&     a_advectVelocity,
                                     const Real&     a_advectAngle,
                                     const Real&     a_evalTime)
{
  m_ambientDensity = a_ambientDensity;
  m_deltaDensity = a_deltaDensity;
  m_center = a_center;
  m_size = a_size;
  m_advectVelocity = a_advectVelocity;
  m_advectAngle = a_advectAngle;
  m_evalTime = a_evalTime;
  m_velType = SOLIDBODY;
  m_useVelUniform = false;
  //  m_haveAdvVel = false;
}

AdvectCosineBellIBC::~AdvectCosineBellIBC()
{
}


void
AdvectCosineBellIBC::setUniformVel()
{
  m_velType = UNIFORM;
}

void
AdvectCosineBellIBC::setSolidBodyRotation()
{
  m_velType = SOLIDBODY;
}

void
AdvectCosineBellIBC::setVelUniform(const RealVect&   a_velUniform)
{
  m_velUniform = a_velUniform;
  m_useVelUniform = true;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysMappedIBC* AdvectCosineBellIBC::new_physIBC()
{
  AdvectCosineBellIBC* retval = new AdvectCosineBellIBC();
/*
  retval->m_isFortranCommonSet = m_isFortranCommonSet;
*/
  if (m_haveTime) retval->setTime(m_time);
  if (m_haveCoordSys) retval->setCoordSys(m_coordSysPtr);
  retval->m_velType = m_velType;
  retval->m_ambientDensity = m_ambientDensity;
  retval->m_deltaDensity = m_deltaDensity;
  retval->m_center = m_center;
  retval->m_size = m_size;
  retval->m_advectVelocity = m_advectVelocity;
  retval->m_advectAngle = m_advectAngle;
  retval->m_evalTime = m_evalTime;
  retval->m_useVelUniform = m_useVelUniform;
  if (m_useVelUniform) retval->m_velUniform = m_velUniform;
  return static_cast<PhysMappedIBC*>(retval);
}

// Set up initial conditions
void AdvectCosineBellIBC::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isDefined == true);
  CH_assert(m_haveCoordSys);
  CH_assert(m_haveTime);

  Real InitialTheta = (m_evalTime + m_time) * m_advectVelocity;

  // const LevelData<FArrayBox>& cellAvgJ = m_coordSysPtr->getJ();
  // int nComp = a_U.nComp();
  const DisjointBoxLayout& layout = a_U.disjointBoxLayout();
  DataIterator dit = layout.dataIterator();

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& baseBox = layout[dit];

      // Storage for current grid
      FArrayBox& UFab = a_U[dit];

      // Box of current grid
      Box uBox = UFab.box();
      // removed by petermc, 9 Feb 2011
      // uBox &= m_domain;

      const NewCoordSys* blockPtr = m_coordSysPtr->getCoordSys(baseBox);
        
      // Xi: mapped space coordinates
      FArrayBox XiFab(uBox, SpaceDim);
      blockPtr->getCenterMappedCoordinates(XiFab, uBox);

      const CubedSphere2DPanelCS* cubedSphereBlockPtr =
        dynamic_cast<const CubedSphere2DPanelCS*>(blockPtr);
      if (cubedSphereBlockPtr != NULL)
        {
          // For each point:
          // set RealVect Xi, which is just linear,
          // and then RealVect X( m_coordSysPtr->realCoord( Xi ) );
          // and then Real J( m_coordSysPtr->pointwiseJ( X ) );

          // lonlat: longitude and latitude
          FArrayBox lonlatFab(uBox, SpaceDim);
          const int LON = 0;
          const int LAT = 1;

          // Evaluate cosine bell
          BoxIterator bit(uBox);
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();

              // Get equiangular coordinates
              RealVect xi;
              xi[0] = XiFab(iv,0);
              xi[1] = XiFab(iv,1);

              Real xyz[3];
              cubedSphereBlockPtr->pointTransformEquiangularToCartesian(xi, xyz);

              // Rotated sample points.
              Real xyzR[3];

              xyzR[0] =
                xyz[0] * (sin(m_advectAngle) * sin(m_advectAngle)
                          + cos(m_advectAngle) * cos(m_advectAngle) * cos(InitialTheta))
                + xyz[1] * cos(m_advectAngle) * sin(InitialTheta)
                - xyz[2] * sin(m_advectAngle) * cos(m_advectAngle)
                * (1.0 - cos(InitialTheta));
              
              xyzR[1] =
                - xyz[0] * cos(m_advectAngle) * sin(InitialTheta)
                + xyz[1] * cos(InitialTheta)
                - xyz[2] * sin(m_advectAngle) * sin(InitialTheta);
              
              xyzR[2] =
                - xyz[0] * sin(m_advectAngle) * cos(m_advectAngle)
                * (1.0 - cos(InitialTheta))
                + xyz[1] * sin(m_advectAngle) * sin(InitialTheta)
                + xyz[2] * (cos(m_advectAngle) * cos(m_advectAngle)
                            + sin(m_advectAngle) * sin(m_advectAngle) * cos(InitialTheta));
              
              // Sample points in RLL coordinates:  longitude, then latitude
              lonlatFab(iv, LON) = atan2(xyzR[1], xyzR[0]); // - M_PI / 2.0;
              lonlatFab(iv, LAT) = asin(xyzR[2]);
            }
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              
              Real dLon = lonlatFab(iv, LON);
              Real dLat = lonlatFab(iv, LAT);
              
              // Great circle distance
              Real dR1 = cos(dLat) * sin(dLon - m_center[LON]);
              Real dR2 = cos(m_center[LAT]) * sin(dLat)
                - sin(m_center[LAT]) * cos(dLat) * cos(dLon - m_center[LON]);
              Real dR3 = sin(m_center[LAT]) * sin(dLat)
                + cos(m_center[LAT]) * cos(dLat) * cos(dLon - m_center[LON]);
              
              Real dR = atan2(sqrt(dR1 * dR1 + dR2 * dR2), dR3);

              UFab(iv, 0) = cosineBell(dR);
            }
        }
      else // NOT cubed-sphere, so we can call blockPtr->realCoord().
        {
          // Find where center of cosine bell is now.
          RealVect centerNow = m_center;
          if (m_velType == UNIFORM)
            {
              if (m_useVelUniform)
                {
                  centerNow += InitialTheta * m_velUniform;
                }
              else
                {
                  // InitialTheta is straight-line distance.
                  Real cosAngle = cos(m_advectAngle);
                  Real sinAngle = sin(m_advectAngle);
                  centerNow[0] += InitialTheta * cosAngle;
                  centerNow[1] += InitialTheta * sinAngle;
#if CH_SPACEDIM == 3
                  centerNow[2] += InitialTheta;
#endif
                }
            }
          else if (m_velType == SOLIDBODY)
            {
              Real thetaOffset = (2.*M_PI) * InitialTheta;
              Real cth = cos(thetaOffset);
              Real sth = sin(thetaOffset);
              centerNow[0] = cth*m_center[0] - sth*m_center[1];
              centerNow[1] = sth*m_center[0] + cth*m_center[1];
            }

          // Not cubed sphere, so call blockPtr->realCoord().
          FArrayBox XYZFab(uBox, SpaceDim);
          blockPtr->realCoord(XYZFab, XiFab, uBox);
          BoxIterator bit(uBox);
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              
              Real dist2 = 0.;
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  Real disp = XYZFab(iv, idir) - centerNow[idir];
                  dist2 += disp*disp;
                }
              Real dist = sqrt(dist2);
              UFab(iv,0) = cosineBell(dist);
            }
        }
    }
}


Real AdvectCosineBellIBC::cosineBell(Real a_dist)
{
  Real U = m_ambientDensity;
  if (a_dist <= m_size)
    {
      // fac(x) is between 0 and 1; fac(0) = 0; fac(m_size) = 1.
      Real fac = (1.0 + cos(M_PI * a_dist / m_size)) / 2.0;
      U += m_deltaDensity * fac*fac*fac;
    }
  return U;
}

// Return the advection velocity at each point
void AdvectCosineBellIBC::advVel(FArrayBox& a_advVelFab,
                                 Real a_time)
{
  CH_TIME("AdvectCosineBellIBC::advVel");
  CH_assert(m_haveCoordSys);

  bool haveAdvVelSet = false;

  if (m_velType == UNIFORM)
    {
      if (m_useVelUniform)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              a_advVelFab.setVal(m_advectVelocity * m_velUniform[idir], idir);
            }
        }
      else
        {
          a_advVelFab.setVal(m_advectVelocity * cos(m_advectAngle), 0);
          a_advVelFab.setVal(m_advectVelocity * sin(m_advectAngle), 1);
#if CH_SPACEDIM == 3
          a_advVelFab.setVal(m_advectVelocity, 2);
#endif
        }
      haveAdvVelSet = true;
    }
  else if (m_velType == SOLIDBODY)
    {
      // Box of current grid
      Box uBox = a_advVelFab.box();
      // removed by petermc, 4 Jan 2011
      // uBox &= m_domain;
      int blockNum = m_coordSysPtr->whichBlockOverlap(uBox);

      const NewCoordSys* blockPtr = m_coordSysPtr->getCoordSys(blockNum);

      // Xi: mapped space coordinates
      FArrayBox XiFab(uBox, SpaceDim);
      blockPtr->getCenterMappedCoordinates(XiFab, uBox);

      const CubedSphere2DPanelCS* cubedSphereBlockPtr =
        dynamic_cast<const CubedSphere2DPanelCS*>(blockPtr);
      if (cubedSphereBlockPtr != NULL)
        {
          // Xi: mapped space coordinates
          FArrayBox XiFab(uBox, SpaceDim);
          blockPtr->getCenterMappedCoordinates(XiFab, uBox);

          const int LON = 0;
          const int LAT = 1;

          // Iterate through all elements
          BoxIterator bit(uBox);
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
      
              // Determine corresponding longitude-latitude point
              RealVect ptXi;
              ptXi[0] = XiFab(iv, 0);
              ptXi[1] = XiFab(iv, 1);
          
              RealVect ptRLL;

              cubedSphereBlockPtr->pointTransformEquiangularToLonLat(ptXi, ptRLL);

              // Calculate longitude-latitude vector field
              // vecRLL[0:1] is vector field in longitude-latitude coordinate basis
              Real vecRLL[2];
              vecRLL[LON] = m_advectVelocity * cos(ptRLL[LAT]) *
                (cos(m_advectAngle) +
                 cos(ptRLL[LON]) * tan(ptRLL[LAT]) * sin(m_advectAngle));
              vecRLL[LAT] = - m_advectVelocity * sin(ptRLL[LON]) * sin(m_advectAngle);
      
              // vecCS[0:1] is vector field in mapped-coordinate basis
              Real vecCS[2];
              cubedSphereBlockPtr->vectorTransformLatLonToEquiangular(ptXi,
                                                                      vecRLL,
                                                                      vecCS);
              /*
                if (cubedSphereBlockPtr->panel() == 1)
                {
                printf("RLL: %1.5e %1.5e\n", vecRLL[0], vecRLL[1]);
                printf("CS: %1.5e %1.5e\n", vecCS[0], vecCS[1]);
                }
              */

              a_advVelFab(iv, 0) = vecCS[0];
              a_advVelFab(iv, 1) = vecCS[1];
            }
          haveAdvVelSet = true;
        }

      const CylinderEquiangularBlockCS* cylinderEquiangularBlockPtr =
        dynamic_cast<const CylinderEquiangularBlockCS*>(blockPtr);
      if (cylinderEquiangularBlockPtr != NULL)
        {
          FArrayBox XYZFab(uBox, SpaceDim);
          blockPtr->realCoord(XYZFab, XiFab, uBox);
          // Set a_advVelFab[0] = -XYZFab[1].
          a_advVelFab.copy(XYZFab, 1, 0);
          a_advVelFab.negate(0);
          // Set a_advVelFab[1] = XYZFab[0].
          a_advVelFab.copy(XYZFab, 0, 1);
          a_advVelFab *= (2.*M_PI) * m_advectVelocity;
          haveAdvVelSet = true;
        }

      const CylinderSpokesBlockCS* cylinderSpokesBlockPtr =
        dynamic_cast<const CylinderSpokesBlockCS*>(blockPtr);
      if (cylinderSpokesBlockPtr != NULL)
        {
          FArrayBox XYZFab(uBox, SpaceDim);
          blockPtr->realCoord(XYZFab, XiFab, uBox);
          // Set a_advVelFab[0] = -XYZFab[1].
          a_advVelFab.copy(XYZFab, 1, 0);
          a_advVelFab.negate(0);
          // Set a_advVelFab[1] = XYZFab[0].
          a_advVelFab.copy(XYZFab, 0, 1);
          a_advVelFab *= (2.*M_PI) * m_advectVelocity;
          haveAdvVelSet = true;
        }

    }
  CH_assert(haveAdvVelSet);
  //  m_advVel[dit].copy(a_advVelFab);

  // m_haveAdvVel = true;

  // convert point values into 4th-order cell averages
  // petermc, 1 Oct 2009:  This is to be done outside, if requested.
  // fourthOrderAverage(a_advVel);

  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}

void AdvectCosineBellIBC::primBC(FArrayBox&            a_WGdnv,
                                 const FArrayBox&      a_Wextrap,
                                 const FArrayBox&      a_W,
                                 const int&            a_dir,
                                 const Side::LoHiSide& a_side,
                                 const Real&           a_time)
{
  const Box& bx = a_Wextrap.box();
  // Set a_WGdnv to zero on boundary faces.
  Real val = 0.;
  int startcomp = 0;
  int ncomp = 1;
  a_WGdnv.setVal(val, bx, startcomp, ncomp);

  // This creates boundary error:
  // a_WGdnv.copy(a_Wextrap, bx);

  // This creates boundary error:
  //  // Set a_WGdnv so that a_W == mean(a_WGdnv(bdry), a_WGdnv(bdry-1)).
  //  // Hence a_WGdnv(bdry) = 2*a_W - a_WGdnv(bdry-1).
  //  // a_WGdnv.copy(a_Wextrap, bx);
  //  int signSide = sign(a_side);
  //  Box bxCells(bx);
  //  bxCells.shiftHalf(a_dir, -signSide);
  //  Box bx1(bx);
  //  bx1.shift(a_dir, -signSide);
  //  a_WGdnv.copy(a_WGdnv, bx1, 0, bx, 0, 1);
  //  a_WGdnv.negate(bx);
  //  // source and dest BaseFabs must be of the same type.
  //  a_WGdnv.shiftHalf(a_dir, -signSide);
  //  a_WGdnv.plus(a_W, bxCells, bxCells, 2., 0, 0, 1);
  //  a_WGdnv.shiftHalf(a_dir, signSide);
}


#include "NamespaceFooter.H"
