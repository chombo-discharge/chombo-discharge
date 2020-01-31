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

#include "NamespaceHeader.H"

// Null constructor
AdvectCosineBellIBC::AdvectCosineBellIBC()
{
  m_gotAdvVel = false;
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
  m_gotAdvVel = false;
}

AdvectCosineBellIBC::~AdvectCosineBellIBC()
{
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
  if (m_gotTime) retval->setTime(m_time);
  if (m_gotCoordSys) retval->setCoordSys(m_coordSysPtr);
*/
  retval->m_ambientDensity = m_ambientDensity;
  retval->m_deltaDensity = m_deltaDensity;
  retval->m_center = m_center;
  retval->m_size = m_size;
  retval->m_advectVelocity = m_advectVelocity;
  retval->m_advectAngle = m_advectAngle;
  retval->m_evalTime = m_evalTime;

  return static_cast<PhysMappedIBC*>(retval);
}

// Set up initial conditions
void AdvectCosineBellIBC::initialize(LevelData<FArrayBox>& a_U)
{
  initializeUnified(a_U, false);
}

void AdvectCosineBellIBC::initializeWithJ(LevelData<FArrayBox>& a_U)
{
  initializeUnified(a_U, true);
}

// Set up initial conditions
void AdvectCosineBellIBC::initializeUnified(LevelData<FArrayBox>& a_U,
                                            bool a_includeJ)
{
  CH_assert(m_isDefined == true);
  CH_assert(m_gotCoordSys);
  CH_assert(m_gotTime);

  const double InitialTheta = (m_evalTime + m_time) * m_advectVelocity;

  // const LevelData<FArrayBox>& cellAvgJ = m_coordSysPtr->getJ();
  int nComp = a_U.nComp();
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

      const CubedSphere2DPanelCS* coordSysBlockPtr =
        dynamic_cast<const CubedSphere2DPanelCS*>(
          m_coordSysPtr->getCoordSys(baseBox));

      CH_assert(coordSysBlockPtr);

      // For each point:
      // set RealVect Xi, which is just linear,
      // and then RealVect X( m_coordSysPtr->realCoord( Xi ) );
      // and then Real J( m_coordSysPtr->pointwiseJ( X ) );

      // Xi: mapped space coordinates
      FArrayBox XiFab(uBox, SpaceDim);

      coordSysBlockPtr->getCellMappedCoordinates(XiFab, uBox);

      // lonlat: longitude and latitude
      FArrayBox lonlatFab(uBox, SpaceDim);

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
          coordSysBlockPtr->pointTransformEquiangularToCartesian(xi, xyz);

          // Rotated sample points
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
          lonlatFab(iv, 0) = atan2(xyzR[1], xyzR[0]); // - M_PI / 2.0;
          lonlatFab(iv, 1) = asin(xyzR[2]);
        }
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();

          Real dLon = lonlatFab(iv, 0);
          Real dLat = lonlatFab(iv, 1);

          // Great circle distance
          Real dR1 = cos(dLat) * sin(dLon - m_center[1]);
          Real dR2 = cos(m_center[0]) * sin(dLat)
                   - sin(m_center[0]) * cos(dLat) * cos(dLon - m_center[1]);
          Real dR3 = sin(m_center[0]) * sin(dLat)
                   + cos(m_center[0]) * cos(dLat) * cos(dLon - m_center[1]);

          Real dR = atan2(sqrt(dR1 * dR1 + dR2 * dR2), dR3);

          UFab(iv,0) = m_ambientDensity;
          if (dR <= m_size)
            {
              // fac(dR) is between 0 and 1; fac(0) = 0; fac(m_size) = 1.
              Real fac = (1.0 + cos(M_PI * dR / m_size)) / 2.0;
              UFab(iv,0) += m_deltaDensity * fac*fac*fac;
              // UFab(iv,0) += m_deltaDensity * fac*fac;
            }
        }

      // Multiply every component of U by J.
      if (a_includeJ)
        {
          // Get pointwise J (not cell-averaged <J>).
          FArrayBox JFab(uBox, 1);

          coordSysBlockPtr->pointwiseJ(JFab, XiFab, uBox);

          for (int comp = 0; comp < nComp; comp++)
            {
              UFab.mult(JFab, 0, comp);
            }
        }
    }
}

// Set boundary fluxes
void AdvectCosineBellIBC::primBC(FArrayBox&            a_WGdnv,
                         const FArrayBox&      a_Wextrap,
                         const FArrayBox&      a_W,
                         const int&            a_dir,
                         const Side::LoHiSide& a_side,
                         const Real&           a_time)
{
}

// Set boundary slopes:
//   The boundary slopes in a_dW are already set to one sided difference
//   approximations.  If this function doesn't change them they will be
//   used for the slopes at the boundaries.
void AdvectCosineBellIBC::setBdrySlopes(FArrayBox&       a_dW,
                                const FArrayBox& a_W,
                                const int&       a_dir,
                                const Real&      a_time)
{
  MayDay::Error("AdvectCosineBellIBC::setBdrySlopes not defined");
  // Cartesian version calls FORT_SLOPEBCSF.
}

// Adjust boundary fluxes to account for artificial viscosity
// Do nothing for this problem
void AdvectCosineBellIBC::artViscBC(FArrayBox&       a_F,
                            const FArrayBox& a_U,
                            const FArrayBox& a_divVel,
                            const int&       a_dir,
                            const Real&      a_time)
{
}

// Return the advection velocity at each point
void AdvectCosineBellIBC::advVel(LevelData<FArrayBox>& a_advVel,
                                 Real a_time)
{
  CH_TIME("AdvectCosineBellIBC::advVel");
  // For this particular IBC, a_advVel will be independent of a_time.
  const DisjointBoxLayout& layout = a_advVel.disjointBoxLayout();
  DataIterator dit = a_advVel.dataIterator();
  if (m_gotAdvVel)
    { // assumes the layout and ghosts of a_advVel are always the same.
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& velFab = a_advVel[dit];
          velFab.copy(m_advVel[dit]);
        }
      return;
    }
  // m_advVel not defined yet.
  int ncomp = a_advVel.nComp();
  const IntVect& ghostVect = a_advVel.ghostVect();
  m_advVel.define(layout, ncomp, ghostVect);
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& baseBox = layout[dit];

      // Storage for current grid
      FArrayBox& velFab = a_advVel[dit];

      // Box of current grid
      Box uBox = velFab.box();
      // removed by petermc, 4 Jan 2011
      // uBox &= m_domain;
      const CubedSphere2DPanelCS* coordSysBlockPtr =
        dynamic_cast<const CubedSphere2DPanelCS*>(
          m_coordSysPtr->getCoordSys(baseBox));

      // Xi: mapped space coordinates
      FArrayBox XiFab(uBox, SpaceDim);
      coordSysBlockPtr->getCellMappedCoordinates(XiFab, uBox);

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

          coordSysBlockPtr->pointTransformEquiangularToLonLat(ptXi, ptRLL);

          // Calculate longitude-latitude vector field
          // vecRLL[0:1] is vector field in longitude-latitude coordinate basis
          Real vecRLL[2];
          vecRLL[0] = m_advectVelocity * cos(ptRLL[1]) *
            (cos(m_advectAngle) +
             cos(ptRLL[0]) * tan(ptRLL[1]) * sin(m_advectAngle));
          vecRLL[1] = - m_advectVelocity * sin(ptRLL[0]) * sin(m_advectAngle);

          // vecCS[0:1] is vector field in mapped-coordinate basis
          Real vecCS[2];
          coordSysBlockPtr->vectorTransformLatLonToEquiangular(ptXi,
                                                               vecRLL,
                                                               vecCS);
/*
          if (coordSysBlockPtr->panel() == 1)
          {
            printf("RLL: %1.5e %1.5e\n", vecRLL[0], vecRLL[1]);
            printf("CS: %1.5e %1.5e\n", vecCS[0], vecCS[1]);
          }
*/

          velFab(iv, 0) = vecCS[0];
          velFab(iv, 1) = vecCS[1];
        }
      // dummy statement in order to get around gdb bug
      int dummy_unused = 0; dummy_unused = 0;
      m_advVel[dit].copy(velFab);
    }
  m_gotAdvVel = true;

  // convert point values into 4th-order cell averages
  // petermc, 1 Oct 2009:  This is to be done outside, if requested.
  // fourthOrderAverage(a_advVel);

  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}

#include "NamespaceFooter.H"
