#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "GaussianHillsFlowIBC.H"
#include "BoxIterator.H"
#include "CubedSphere2DPanelCS.H"

#include "NamespaceHeader.H"

///////////////////////////////////////////////////////////////////////////////

// Null constructor
GaussianHillsFlowIBC::GaussianHillsFlowIBC()
{
  // CH_assert(false);
  // m_params_are_set = false;
}

// Constructor which defines parameters used by Fortran routines.
// Mapping changes nothing.
GaussianHillsFlowIBC::GaussianHillsFlowIBC(
                     const Real&     a_hmax,
                     const Real&     a_width,
                     const Vector<Real>&  a_longitude,
                     const Vector<Real>&  a_latitude,
                     const Real&     a_period,
                     const Real&     a_kappa,
                     const Real&     a_evalTime) : DefFlowIBC(a_period, a_kappa, a_evalTime)
{
  m_hmax = a_hmax;
  m_width = a_width;
  m_longitude = Vector<Real>(a_longitude);
  m_latitude = Vector<Real>(a_latitude);
  // DefFlowIBC(a_period, a_kappa, a_evalTime);
}

GaussianHillsFlowIBC::~GaussianHillsFlowIBC()
{
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysMappedIBC* GaussianHillsFlowIBC::new_physIBC()
{
  GaussianHillsFlowIBC* retval = new GaussianHillsFlowIBC();
/*
  retval->m_isFortranCommonSet = m_isFortranCommonSet;
*/
  if (m_haveTime) retval->setTime(m_time);
  if (m_haveCoordSys) retval->setCoordSys(m_coordSysPtr);
  retval->m_hmax = m_hmax;
  retval->m_width = m_width;
  retval->m_longitude = Vector<Real>(m_longitude);
  retval->m_latitude = Vector<Real>(m_latitude);
  retval->m_period = m_period;
  retval->m_kappa = m_kappa;
  retval->m_evalTime = m_evalTime;

  return static_cast<PhysMappedIBC*>(retval);
}

// Set up initial conditions
void GaussianHillsFlowIBC::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isDefined == true);
  CH_assert(m_haveCoordSys);
  CH_assert(m_haveTime);

  // const LevelData<FArrayBox>& cellAvgJ = m_coordSysPtr->getJ();
  // int nComp = a_U.nComp();
  const DisjointBoxLayout& layout = a_U.disjointBoxLayout();
  DataIterator dit = layout.dataIterator();

  const int nhills = m_longitude.size();
  Vector<Real> dXi(nhills);
  Vector<Real> dYi(nhills);
  Vector<Real> dZi(nhills);
  for (int i = 0; i < nhills; i++)
    {
      dXi[i] = cos(m_latitude[i]) * cos(m_longitude[i]);
      dYi[i] = cos(m_latitude[i]) * sin(m_longitude[i]);
      dZi[i] = sin(m_latitude[i]);
    }

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

          Real xyz[3];
          coordSysBlockPtr->pointTransformEquiangularToCartesian(xi, xyz);

          // Sample points in RLL coordinates
          // Real dLat = asin(xyz[2]);
          // Real dLon = atan2(xyz[1], xyz[0]) - M_PI / 2.0;

          // Pointwise height
          Real dH = 0.;
          for (int i = 0; i < nhills; i++)
            {
              Real dCartX = xyz[0] - dXi[i];
              Real dCartY = xyz[1] - dYi[i];
              Real dCartZ = xyz[2] - dZi[i];
              Real dist2 = dCartX*dCartX + dCartY*dCartY + dCartZ*dCartZ;
              dH += m_hmax * exp(-m_width * dist2);
            }
          UFab(iv,0) = dH;
        }
    }
}

#include "NamespaceFooter.H"
