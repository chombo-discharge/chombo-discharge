#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BoxIterator.H"

#include "CubedSphere2DPanelCS.H"

#include "BalancedFlowIBC.H"
#include "NamespaceHeader.H"

// Null constructor
BalancedFlowIBC::BalancedFlowIBC()
{
  CH_assert(false);
}

// Constructor which defines parameters used by Fortran routines.
// Mapping changes nothing.
BalancedFlowIBC::BalancedFlowIBC(
                                 const Real&     a_bgVelocity,
                                 const Real&     a_bgHeight)
{
  m_bgVelocity = a_bgVelocity;
  m_bgHeight = a_bgHeight;
}

BalancedFlowIBC::~BalancedFlowIBC()
{
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysMappedIBC* BalancedFlowIBC::new_physIBC()
{
  // BalancedFlowIBC* retval = new BalancedFlowIBC();
  BalancedFlowIBC* retval = new BalancedFlowIBC(m_bgVelocity, m_bgHeight);
  if (m_isFortranCommonSet) retval->setFortranCommon(m_gravity, m_omega, m_alpha);
/*
  if (m_gotTime) retval->setTime(m_time);
*/
  if (m_gotCoordSys) retval->setCoordSys(m_coordSysPtr);
  return static_cast<PhysMappedIBC*>(retval);
}

// Set up initial conditions
void BalancedFlowIBC::initializeUnified(LevelData<FArrayBox>& a_U,
                                        bool a_multiplyJ)
{
  CH_assert(m_isDefined == true);
  CH_assert(m_gotCoordSys);
  CH_assert(m_gotTime);

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

      // Create RLL velocity Fab
      FArrayBox vecRLLFab(uBox, SpaceDim);

      // Create mapped velocity Fab
      FArrayBox velFab(uBox, SpaceDim);

      // Evaluate flow field
      BoxIterator bit(uBox);
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();

          // Get lonlat coordinate
          RealVect xi;
          xi[0] = XiFab(iv,0);
          xi[1] = XiFab(iv,1);
          RealVect lonlat;
          coordSysBlockPtr->pointTransformEquiangularToLonLat(xi, lonlat);
          Real dLon = lonlat[0];
          Real dLat = lonlat[1];

          // Calculate fluid height
          Real dH = (- cos(dLon) * cos(dLat) * sin(m_alpha)
                     + sin(dLat) * cos(m_alpha));

          dH = m_bgHeight - dH * dH * (m_omega * m_bgVelocity
               + 0.5 * m_bgVelocity * m_bgVelocity) / m_gravity;

          UFab(iv,0) = dH;

          // Calculate fluid velocities in longitude-latitude coordinates
          Real dUlon = m_bgVelocity
            * (cos(dLat) * cos(m_alpha) + cos(dLon) * sin(dLat) * sin(m_alpha));

          Real dUlat = - m_bgVelocity * sin(dLon) * sin(m_alpha);

          vecRLLFab(iv,0) = dUlon;
          vecRLLFab(iv,1) = dUlat;
        }

      // Convert velocities from longitude-latitude space to mapped space
      coordSysBlockPtr->fabVectorTransformLatLonToEquiangular(XiFab,
                                                              vecRLLFab,
                                                              velFab);

      // Copy converted velocities
      for (bit.begin(); bit.ok(); ++bit)
        {
          IntVect iv = bit();
          UFab(iv,1) = UFab(iv,0) * velFab(iv,0);
          UFab(iv,2) = UFab(iv,0) * velFab(iv,1);
        }

      // Multiply every component of U by J
      if (a_multiplyJ)
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

// Set up initial conditions
void BalancedFlowIBC::initialize(LevelData<FArrayBox>& a_U)
{
  // Initialize the flow field
  initializeUnified(a_U, false);
}

// Set up initial conditions with J
void BalancedFlowIBC::initializeWithJ(LevelData<FArrayBox>& a_U)
{
  // Initialize the flow field with J
  initializeUnified(a_U, true);
}

// Set boundary fluxes
void BalancedFlowIBC::primBC(FArrayBox&            a_WGdnv,
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
void BalancedFlowIBC::setBdrySlopes(FArrayBox&       a_dW,
                                const FArrayBox& a_W,
                                const int&       a_dir,
                                const Real&      a_time)
{
  MayDay::Error("BalancedFlowIBC::setBdrySlopes not defined");
  // Cartesian version calls FORT_SLOPEBCSF.
}

// Adjust boundary fluxes to account for artificial viscosity
// Do nothing for this problem
void BalancedFlowIBC::artViscBC(FArrayBox&       a_F,
                            const FArrayBox& a_U,
                            const FArrayBox& a_divVel,
                            const int&       a_dir,
                            const Real&      a_time)
{
}

#include "NamespaceFooter.H"
