#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DefFlowIBC.H"
#include "DefFlowIBCF_F.H"
#include "BoxIterator.H"
#include "CubedSphere2DPanelCS.H"

#include "NamespaceHeader.H"

///////////////////////////////////////////////////////////////////////////////

// Null constructor
DefFlowIBC::DefFlowIBC()
{
  // CH_assert(false);
  // m_params_are_set = false;
}

// Constructor which defines parameters used by Fortran routines.
// Mapping changes nothing.
DefFlowIBC::DefFlowIBC(
                       const Real&  a_period,
                       const Real&  a_kappa,
                       const Real&  a_evalTime)
{
  m_period = a_period;
  m_kappa = a_kappa;
  m_evalTime = a_evalTime;
}

DefFlowIBC::~DefFlowIBC()
{
}

// Return the advection velocity at each point
void DefFlowIBC::advVel(FArrayBox& a_advVelFab,
                        Real a_time)
{
  CH_TIME("DefFlowIBC::advVel");
  CH_assert(m_haveCoordSys);

  // For this particular IBC, a_advVel will be independent of a_time.
  Real time = a_time + m_evalTime;

  // Box of current grid
  Box uBox = a_advVelFab.box();
  int blockNum = m_coordSysPtr->whichBlockOverlap(uBox);
  // removed by petermc, 4 Jan 2011
  // uBox &= m_domain;
  const CubedSphere2DPanelCS* coordSysBlockPtr =
    dynamic_cast<const CubedSphere2DPanelCS*>(m_coordSysPtr->getCoordSys(blockNum));
                                              
  // Xi: mapped space coordinates
  FArrayBox XiFab(uBox, SpaceDim);
  coordSysBlockPtr->getCellMappedCoordinates(XiFab, uBox);

  // rll: longitude-latitude coordinates
  FArrayBox rllFab(uBox, SpaceDim);
  coordSysBlockPtr->fabTransformEquiangularToLonLat(XiFab, rllFab);

  // vecRLL: vector in longitude-latitude basis
  FArrayBox vecRLLFab(uBox, SpaceDim);
  FORT_VECLATLONDEFFLOW(CHF_FRA(vecRLLFab),
                        CHF_CONST_FRA(rllFab),
                        CHF_CONST_REAL(m_kappa),
                        CHF_CONST_REAL(time),
                        CHF_CONST_REAL(m_period));

  // Conver to a_advVelFab: vector in equiangular basis
  coordSysBlockPtr->fabVectorTransformLatLonToEquiangular(XiFab,
                                                          vecRLLFab,
                                                          a_advVelFab);

  // convert point values into 4th-order cell averages
  // petermc, 1 Oct 2009:  This is to be done outside, if requested.
  // fourthOrderAverage(a_advVel);

  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}

#include "NamespaceFooter.H"
