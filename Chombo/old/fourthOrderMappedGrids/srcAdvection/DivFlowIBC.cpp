#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DivFlowIBC.H"
#include "DivFlowIBCF_F.H"
#include "BoxIterator.H"
#include "CubedSphere2DPanelCS.H"

#include "NamespaceHeader.H"

///////////////////////////////////////////////////////////////////////////////

// Null constructor
DivFlowIBC::DivFlowIBC()
{
  // CH_assert(false);
  // m_params_are_set = false;
}

// Constructor which defines parameters used by Fortran routines.
// Mapping changes nothing.
DivFlowIBC::DivFlowIBC(
                       const Real&  a_period,
                       const Real&  a_k,
                       const Real&  a_evalTime)
{
  m_period = a_period;
  m_k = a_k;
  m_evalTime = a_evalTime;
}

DivFlowIBC::~DivFlowIBC()
{
}

// Set up initial conditions
void DivFlowIBC::initialize(LevelData<FArrayBox>& a_U)
{
  initializeUnified(a_U, false);
}

void DivFlowIBC::initializeWithJ(LevelData<FArrayBox>& a_U)
{
  initializeUnified(a_U, true);
}

// Set boundary fluxes
void DivFlowIBC::primBC(
                        FArrayBox&            a_WGdnv,
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
void DivFlowIBC::setBdrySlopes(
                               FArrayBox&       a_dW,
                               const FArrayBox& a_W,
                               const int&       a_dir,
                               const Real&      a_time)
{
  MayDay::Error("DivFlowIBC::setBdrySlopes not defined");
  // Cartesian version calls FORT_SLOPEBCSF.
}

// Adjust boundary fluxes to account for artificial viscosity
// Do nothing for this problem
void DivFlowIBC::artViscBC(
                           FArrayBox&       a_F,
                           const FArrayBox& a_U,
                           const FArrayBox& a_divVel,
                           const int&       a_dir,
                           const Real&      a_time)
{
}

// Return the advection velocity at each point
void DivFlowIBC::advVel(LevelData<FArrayBox>& a_advVel,
                        Real a_time)
{
  CH_TIME("DivFlowIBC::advVel");
  CH_assert(m_gotCoordSys);

  // For this particular IBC, a_advVel will be independent of a_time.
  const DisjointBoxLayout& layout = a_advVel.disjointBoxLayout();
  DataIterator dit = a_advVel.dataIterator();
  Real time = a_time + m_evalTime;
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

      // rll: longitude-latitude coordinates
      FArrayBox rllFab(uBox, SpaceDim);
      coordSysBlockPtr->fabTransformEquiangularToLonLat(XiFab, rllFab);

      // vecRLL: vector in longitude-latitude basis
      FArrayBox vecRLLFab(uBox, SpaceDim);
      FORT_VECLATLONDIVFLOW(CHF_FRA(vecRLLFab),
                            CHF_CONST_FRA(rllFab),
                            CHF_CONST_REAL(m_k),
                            CHF_CONST_REAL(time),
                            CHF_CONST_REAL(m_period));

      // Conver to velFab: vector in equiangular basis
      coordSysBlockPtr->fabVectorTransformLatLonToEquiangular(XiFab,
                                                              vecRLLFab,
                                                              velFab);

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
