#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DeformationalFlowIBC.H"
#include "BoxIterator.H"
#include "CubedSphere2DPanelCS.H"

#include "NamespaceHeader.H"

///////////////////////////////////////////////////////////////////////////////

#define COSINEBELLS
//#define GAUSSIANHILLS
//#define SLOTTEDCYLINDERS

///////////////////////////////////////////////////////////////////////////////

// Null constructor
DeformationalFlowIBC::DeformationalFlowIBC()
{
  // CH_assert(false);
  // m_params_are_set = false;
}

// Constructor which defines parameters used by Fortran routines.
// Mapping changes nothing.
DeformationalFlowIBC::DeformationalFlowIBC(
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
}

DeformationalFlowIBC::~DeformationalFlowIBC()
{
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysMappedIBC* DeformationalFlowIBC::new_physIBC()
{
  DeformationalFlowIBC* retval = new DeformationalFlowIBC();
/*
  retval->m_isFortranCommonSet = m_isFortranCommonSet;
*/
  if (m_haveTime) retval->setTime(m_time);
  if (m_haveCoordSys) retval->setCoordSys(m_coordSysPtr);
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
void DeformationalFlowIBC::initialize(LevelData<FArrayBox>& a_U)
{
  // const Real InitialTheta = m_evalTime * m_advectVelocity;

  CH_assert(m_isDefined == true);
  CH_assert(m_haveCoordSys);
  CH_assert(m_haveTime);

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

      // Evaluate cosine bell
      BoxIterator bit(uBox);
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();

          // Get equiangular coordinates
          RealVect xi;
          xi[0] = XiFab(iv,0);
          xi[1] = XiFab(iv,1);

          Real xyzR[3];
          coordSysBlockPtr->pointTransformEquiangularToCartesian(xi, xyzR);

          RealVect ptRLL;
          coordSysBlockPtr->pointTransformEquiangularToLonLat(xi, ptRLL);

          // Sample points in RLL coordinates
          Real dLat = ptRLL[0];
          Real dLon = ptRLL[1];

          // Pointwise height
          Real dH;

#ifdef COSINEBELLS
          Real dLoni[] =
          {
            5.0 * M_PI / 6.0, 7.0 * M_PI / 6.0
          };

          Real dLati[] =
          {
            0.0, 0.0
          };

          const Real dBaseR = 0.5;
          const Real dBaseH = 1.0;
          const Real dBaseB = 0.1;
          const Real dBaseC = 0.9;

          Real dRi[2];

          dRi[0] = acos(sin(dLati[0]) * sin(dLat)
                  + cos(dLati[0]) * cos(dLat) * cos(dLon - dLoni[0]));
          dRi[1] = acos(sin(dLati[1]) * sin(dLat)
                  + cos(dLati[1]) * cos(dLat) * cos(dLon - dLoni[1]));

          if (dRi[0] < dBaseR)
          {
              dH = 0.5 * dBaseH * (1.0 + cos(M_PI * dRi[0] / dBaseR));
              dH = dBaseB * dBaseH + dBaseC * dH;

          }
          else if (dRi[1] < dBaseR)
          {
              dH = 0.5 * dBaseH * (1.0 + cos(M_PI * dRi[1] / dBaseR));
              dH = dBaseB * dBaseH + dBaseC * dH;

          }
          else
          {
              dH = dBaseB;
          }
#endif
#ifdef GAUSSIANHILLS
                        const Real BaseH = 0.95;

                        Real dLoni[] =
                        {
                          5.0 * M_PI / 6.0, 7.0 * M_PI / 6.0
                        };

                        Real dLati[] =
                        {
                          0.0, 0.0
                        };

                        Real dXi[2];
                        Real dYi[2];
                        Real dZi[2];

                        dXi[0] = cos(dLati[0]) * cos(dLoni[0]);
                        dYi[0] = cos(dLati[0]) * sin(dLoni[0]);
                        dZi[0] = sin(dLati[0]);

                        dXi[1] = cos(dLati[1]) * cos(dLoni[1]);
                        dYi[1] = cos(dLati[1]) * sin(dLoni[1]);
                        dZi[1] = sin(dLati[1]);

                        Real dCartX;
                        Real dCartY;
                        Real dCartZ;

                        dCartX = xyzR[0] - dXi[0];
                        dCartY = xyzR[1] - dYi[0];
                        dCartZ = xyzR[2] - dZi[0];

                        dH = BaseH * exp(
                                - 5.0 * (dCartX * dCartX + dCartY * dCartY + dCartZ * dCartZ));

                        dCartX = xyzR[0] - dXi[1];
                        dCartY = xyzR[1] - dYi[1];
                        dCartZ = xyzR[2] - dZi[1];

                        dH += BaseH * exp(
                                - 5.0 * (dCartX * dCartX + dCartY * dCartY + dCartZ * dCartZ));

#endif
#ifdef SLOTTEDCYLINDERS
                        static const Real BaseC = 1.0;
                        static const Real BaseB = 0.1;
                        static const Real BaseR = 0.5;

                        Real dLoni[] =
                        {
                          5.0 * M_PI / 6.0, 7.0 * M_PI / 6.0
                        };

                        Real dLati[] =
                        {
                          0.0, 0.0
                        };

                        Real dRi[2];

                        dRi[0] = acos(sin(dLati[0]) * sin(dLat)
                                        + cos(dLati[0]) * cos(dLat) * cos(dLon - dLoni[0]));
                        dRi[1] = acos(sin(dLati[1]) * sin(dLat)
                                        + cos(dLati[1]) * cos(dLat) * cos(dLon - dLoni[1]));

                        dH = BaseB;

                        if ((dRi[0] <= BaseR) && (fabs(dLon - dLoni[0]) >= BaseR / 6.0))
                        {
                                dH = BaseC;
                        }
                        if ((dRi[1] <= BaseR) && (fabs(dLon - dLoni[1]) >= BaseR / 6.0))
                        {
                                dH = BaseC;
                        }
                        if ((dRi[0] <= BaseR)
                                && (fabs(dLon - dLoni[0]) < BaseR / 6.0)
                                && ((dLat - dLati[0]) < - 5.0 / 12.0 * BaseR)
                        )
                        {
                                dH = BaseC;
                        }
                        if ((dRi[1] <= BaseR)
                                && (fabs(dLon - dLoni[1]) < BaseR / 6.0)
                                && ((dLat - dLati[1]) > 5.0 / 12.0 * BaseR)
                        )
                        {
                                dH = BaseC;
                        }
#endif

          UFab(iv,0) = dH;
        }
    }
}


// Return the advection velocity at each point
void DeformationalFlowIBC::advVel(FArrayBox& a_advVelFab,
                                  Real a_time)
{
  CH_TIME("DeformationalFlowIBC::advVel");
  CH_assert(m_haveCoordSys);

  // For this particular IBC, a_advVel will be independent of a_time.

  // Wind field parameters
  static const Real Period = 5.0; // * 86400.0;
  static const Real Kappa = 10.0 / Period;

  // Box of current grid
  Box uBox = a_advVelFab.box();
  // removed by petermc, 4 Jan 2011
  // uBox &= m_domain;
  int blockNum = m_coordSysPtr->whichBlockOverlap(uBox);
  const CubedSphere2DPanelCS* coordSysBlockPtr =
    dynamic_cast<const CubedSphere2DPanelCS*>(m_coordSysPtr->getCoordSys(blockNum));
                                              
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

      ptRLL[0] = ptRLL[0] - 2.0 * M_PI * a_time / Period;

      // Calculate longitude-latitude vector field
      // vecRLL[0:1] is vector field in longitude-latitude coordinate basis
      Real vecRLL[2];
      vecRLL[0] =
        Kappa * sin(ptRLL[0]) * sin(ptRLL[0]) * sin(2.0 * ptRLL[1])
        * cos(M_PI * a_time / Period)
        + 2.0 * M_PI * cos(ptRLL[1]) / Period;

      vecRLL[1] =
        Kappa * sin(2.0 * ptRLL[0]) * cos(ptRLL[1])
        * cos(M_PI * a_time / Period);

      // vecCS[0:1] is vector field in mapped-coordinate basis
      Real vecCS[2];
      coordSysBlockPtr->vectorTransformLatLonToEquiangular(ptXi,
                                                           vecRLL,
                                                           vecCS);
      a_advVelFab(iv, 0) = vecCS[0];
      a_advVelFab(iv, 1) = vecCS[1];

      // dummy statement in order to get around gdb bug
      int dummy_unused = 0; dummy_unused = 0;
    }

  // convert point values into 4th-order cell averages
  // petermc, 1 Oct 2009:  This is to be done outside, if requested.
  // fourthOrderAverage(a_advVel);

  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}

#include "NamespaceFooter.H"
