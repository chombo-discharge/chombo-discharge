#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ShallowWaterCubedSphereIBC.H"
#include "BoxIterator.H"
#include "CubedSphere2DPanelCS.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
ShallowWaterCubedSphereIBC::
ShallowWaterCubedSphereIBC(RefCountedPtr<ScalarFunction> a_height,
                           RefCountedPtr<VectorFunction> a_velocity):
  m_height(a_height),
  m_velocity(a_velocity)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
ShallowWaterCubedSphereIBC::
~ShallowWaterCubedSphereIBC()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
PhysMappedIBC*
ShallowWaterCubedSphereIBC::
new_physIBC()
{
  ShallowWaterCubedSphereIBC* retval = new ShallowWaterCubedSphereIBC(m_height, m_velocity);
  return retval;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
ShallowWaterCubedSphereIBC::
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

      // Fill in the initial height and velocity.
      UFab(iv,0) = (*m_height)(xi, m_time);

      RealVect U = (*m_velocity)(xi, m_time);
      UFab(iv,1) = U[0];
      UFab(iv,2) = U[1];
    }
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
ShallowWaterCubedSphereIBC::
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
ShallowWaterCubedSphereIBC::
setBdrySlopes(FArrayBox&       a_dW,
              const FArrayBox& a_W,
              const int&       a_dir,
              const Real&      a_time)
{
  MayDay::Error("ShallowWaterCubedSphereIBC::setBdrySlopes not defined");
  // Cartesian version calls FORT_SLOPEBCSF.
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
ShallowWaterCubedSphereIBC::
artViscBC(FArrayBox&       a_F,
          const FArrayBox& a_U,
          const FArrayBox& a_divVel,
          const int&       a_dir,
          const Real&      a_time)
{
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
