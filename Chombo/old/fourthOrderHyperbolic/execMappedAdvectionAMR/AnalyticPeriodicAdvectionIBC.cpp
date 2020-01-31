#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AnalyticPeriodicAdvectionIBC.H"
#include "BoxIterator.H"
#include "SetCentersF_F.H"
#include "LoHiSide.H"

//-----------------------------------------------------------------------
AnalyticPeriodicAdvectionIBC::
AnalyticPeriodicAdvectionIBC(RefCountedPtr<ScalarFunction> a_initialState):
  PhysMappedIBC(),
  m_initialState(a_initialState)
{
  CH_assert(!a_initialState.isNull());
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
AnalyticPeriodicAdvectionIBC::
~AnalyticPeriodicAdvectionIBC()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PhysMappedIBC*
AnalyticPeriodicAdvectionIBC::
new_physIBC()
{
  AnalyticPeriodicAdvectionIBC* retval =
    new AnalyticPeriodicAdvectionIBC(m_initialState);
  if (m_haveTime) retval->setTime(m_time);
  if (m_haveCoordSys) retval->setCoordSys(m_coordSysPtr);
  return dynamic_cast<PhysMappedIBC*>(retval);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AnalyticPeriodicAdvectionIBC::
initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isDefined == true);
  CH_assert(m_haveCoordSys);

  DataIterator dit = a_U.boxLayout().dataIterator();

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& UFab = a_U[dit()];

    // Box of current grid.
    Box uBox = UFab.box();
//    uBox &= m_domain;

    // Extract the coordinate system for this box.
    Box box = a_U.boxLayout()[dit()];
    const NewCoordSys* coordSys = m_coordSysPtr->getCoordSys(box);
    CH_assert(coordSys != NULL);

    // Xi:  cartesian coordinates
    FArrayBox XiFab(uBox, SpaceDim);
    FORT_SETCELLCENTERS(CHF_FRA(XiFab),
                        CHF_CONST_REAL(m_dx),
                        CHF_BOX(uBox));

    // X:  physical coordinates.
    int dimension = coordSys->dimension();
    FArrayBox XFab(uBox, dimension);
    coordSys->realCoord(XFab, XiFab, uBox);

    // Set up initial condition in this grid
    for (BoxIterator bit(uBox); bit.ok(); ++bit)
    {
      IntVect i = bit();
      RealVect X;
      for (int d = 0; d < SpaceDim; ++d)
        X[d] = XFab(i, d);
      UFab(i, 0) = (*m_initialState)(X, m_time);
    }
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AnalyticPeriodicAdvectionIBC::
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
//    MayDay::Error("AnalyticPeriodicAdvectionIBC only supports periodic BCs");
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AnalyticPeriodicAdvectionIBC::
setBdrySlopes(FArrayBox&       a_dW,
              const FArrayBox& a_W,
              const int&       a_dir,
              const Real&      a_time)
{
//  MayDay::Error("AnalyticPeriodicAdvectionIBC only supports periodic BCs");
}
//-----------------------------------------------------------------------

