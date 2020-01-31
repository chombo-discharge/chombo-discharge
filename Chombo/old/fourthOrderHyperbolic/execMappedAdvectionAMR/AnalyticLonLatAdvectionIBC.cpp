#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AnalyticLonLatAdvectionIBC.H"
#include "BoxIterator.H"
#include "SetCentersF_F.H"
#include "CubedSphere2DCS.H"
#include "CubedSphere2DPanelCS.H"
#include "LoHiSide.H"
#include "DebugOut.H"

//-----------------------------------------------------------------------
AnalyticLonLatAdvectionIBC::
AnalyticLonLatAdvectionIBC(RefCountedPtr<ScalarFunction> a_initialState):
  PhysMappedIBC(),
  m_initialState(a_initialState)
{
  CH_assert(!a_initialState.isNull());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
AnalyticLonLatAdvectionIBC::
~AnalyticLonLatAdvectionIBC()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PhysMappedIBC*
AnalyticLonLatAdvectionIBC::
new_physIBC()
{
  AnalyticLonLatAdvectionIBC* retval =
    new AnalyticLonLatAdvectionIBC(m_initialState);
  if (m_haveTime) retval->setTime(m_time);
  if (m_haveCoordSys) retval->setCoordSys(m_coordSysPtr);
  return dynamic_cast<PhysMappedIBC*>(retval);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AnalyticLonLatAdvectionIBC::
initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isDefined == true);
  CH_assert(m_haveCoordSys);
  CH_assert(dynamic_cast<const CubedSphere2DCS*>(m_coordSysPtr) != NULL);

  DataIterator dit = a_U.boxLayout().dataIterator();

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& UFab = a_U[dit()];

    // Box of current grid, including ghost cells.
    Box uBox = UFab.box();

    // Extract the coordinate system for this box.
    Box box = a_U.boxLayout()[dit()];
    const NewCoordSys* coordSys = m_coordSysPtr->getCoordSys(box);
    CH_assert(coordSys != NULL);
    const CubedSphere2DPanelCS* panelCS = dynamic_cast<const CubedSphere2DPanelCS*>(coordSys);
    CH_assert(panelCS != NULL);

    // Xi:  equiangular coordinates.
    FArrayBox XiFab(uBox, 2);
    panelCS->getCellMappedCoordinates(XiFab, uBox);

    // Lon/lat coordinates.
    FArrayBox lonlat(uBox, 2);
    panelCS->fabTransformEquiangularToLonLat(XiFab, lonlat);

    // Set up initial condition in this panel.
    for (BoxIterator bit(uBox); bit.ok(); ++bit)
    {
      IntVect i = bit();
      // Compute lon/lat coordinates from equiangulars.
      RealVect ll(lonlat(i, 0), lonlat(i, 1));
      UFab(i, 0) = (*m_initialState)(ll, m_time);
    }
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AnalyticLonLatAdvectionIBC::
primBC(FArrayBox&             a_WGdnv,
       const FArrayBox&       a_Wextrap,
       const FArrayBox&       a_W,
       const FArrayBox *const a_unitNormalBasisPtr,
       const Interval&        a_velIntv,
       const int&             a_dir,
       const Side::LoHiSide&  a_side,
       const Real&            a_time)
{
//  if (!m_domain.isPeriodic(a_dir))
//    MayDay::Error("AnalyticLonLatAdvectionIBC only supports spherical domains");
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AnalyticLonLatAdvectionIBC::
setBdrySlopes(FArrayBox&       a_dW,
              const FArrayBox& a_W,
              const int&       a_dir,
              const Real&      a_time)
{
//  MayDay::Error("AnalyticLonLatAdvectionIBC only supports spherical domains");
}
//-----------------------------------------------------------------------

