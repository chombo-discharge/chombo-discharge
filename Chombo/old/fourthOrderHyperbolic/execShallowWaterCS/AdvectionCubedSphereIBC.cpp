#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AdvectionCubedSphereIBC.H"
#include "BoxIterator.H"
#include "CubedSphere2DPanelCS.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
AdvectionCubedSphereIBC::
AdvectionCubedSphereIBC(RefCountedPtr<ScalarFunction> a_solution,
                        RefCountedPtr<VectorFunction> a_velocity):
  m_solution(a_solution),
  m_velocity(a_velocity),
  m_gotAdvVel(false),
  m_advVel()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
AdvectionCubedSphereIBC::
~AdvectionCubedSphereIBC()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PhysMappedIBC*
AdvectionCubedSphereIBC::
new_physIBC()
{
  AdvectionCubedSphereIBC* retval = new AdvectionCubedSphereIBC(m_solution, m_velocity);
  return retval;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AdvectionCubedSphereIBC::
initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isDefined == true);
  CH_assert(m_gotCoordSys);
  CH_assert(m_gotTime);

  // const LevelData<FArrayBox>& cellAvgJ = m_coordSysPtr->getJ();
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

      // Fill in the solution.
      UFab(iv,0) = (*m_solution)(xi, m_time);
    }
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AdvectionCubedSphereIBC::
primBC(FArrayBox& a_WGdnv,
       const FArrayBox& a_Wextrap,
       const FArrayBox&      a_W,
       const int&            a_dir,
       const Side::LoHiSide& a_side,
       const Real&           a_time)
{
  // Periodic boundaries are handled automagically.
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AdvectionCubedSphereIBC::
setBdrySlopes(FArrayBox&       a_dW,
              const FArrayBox& a_W,
              const int&       a_dir,
              const Real&      a_time)
{
  MayDay::Error("AdvectionCubedSphereIBC::setBdrySlopes not defined");
  // Cartesian version calls FORT_SLOPEBCSF.
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AdvectionCubedSphereIBC::
artViscBC(FArrayBox&       a_F,
          const FArrayBox& a_U,
          const FArrayBox& a_divVel,
          const int&       a_dir,
          const Real&      a_time)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
AdvectionCubedSphereIBC::
advVel(LevelData<FArrayBox>& a_advVel,
       Real a_time)
{
  CH_TIME("AdvectionCubedSphereIBC::advVel");
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

      // Write the values.
      RealVect U = (*m_velocity)(ptXi, a_time);
      velFab(iv, 0) = U[0];
      velFab(iv, 1) = U[1];
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
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
