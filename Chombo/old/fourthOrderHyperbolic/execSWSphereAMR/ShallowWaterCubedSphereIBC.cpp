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
#include "SetCentersF_F.H"
#include "ShallowWaterPhysicsF_F.H"
#include "CubedSphere2DPanelCS.H"
#include "SWintegrator.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
ShallowWaterCubedSphereIBC::
ShallowWaterCubedSphereIBC(RefCountedPtr<ScalarFunction> a_height,
                           RefCountedPtr<VectorFunction> a_velocity,
                           Real a_gravity,
                           Real a_omega,
                           Real a_alpha) :
  m_height(a_height),
  m_velocity(a_velocity),
  m_gravity(a_gravity),
  m_omega(a_omega),
  m_alpha(a_alpha)
{
  FORT_SHALLOWWATERSETF(CHF_CONST_REAL(m_gravity),
                        CHF_CONST_REAL(m_omega),
                        CHF_CONST_REAL(m_alpha));
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
  ShallowWaterCubedSphereIBC* retval =
    new ShallowWaterCubedSphereIBC(m_height, m_velocity,
                                   m_gravity, m_omega, m_alpha);
  // if (m_haveTime) retval->setTime(m_time);
  if (m_haveCoordSys) retval->setCoordSys(m_coordSysPtr);
  return dynamic_cast<PhysMappedIBC*>(retval);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
ShallowWaterCubedSphereIBC::
initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isDefined == true);
  CH_assert(m_haveCoordSys);
  // petermc, 1 Sep 2011: do we need this?
  // CH_assert(m_haveTime);

  // const LevelData<FArrayBox>& cellAvgJ = m_coordSysPtr->getJ();
  const DisjointBoxLayout& layout = a_U.disjointBoxLayout();
  DataIterator dit = layout.dataIterator();

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& UFab = a_U[dit];

    // Box of current grid
    Box uBox = UFab.box();
    // removed by petermc, 9 Feb 2011
    // uBox &= m_domain;

    // Extract the coordinate system for this box.
    const Box& baseBox = layout[dit];
    const NewCoordSys* coordSys = m_coordSysPtr->getCoordSys(baseBox);
    const CubedSphere2DPanelCS* panelCS =
      dynamic_cast<const CubedSphere2DPanelCS*>(coordSys);
    CH_assert(panelCS != NULL);

    // For each point:
    // set RealVect Xi, which is just linear,
    // and then RealVect X( m_coordSysPtr->realCoord( Xi ) );
    // and then Real J( m_coordSysPtr->pointwiseJ( X ) );

    // Xi: mapped space coordinates
    FArrayBox XiFab(uBox, SpaceDim);
    FORT_SETCELLCENTERS(CHF_FRA(XiFab),
                        CHF_CONST_REAL(m_dx),
                        CHF_BOX(uBox));

    // Lon/lat coordinates.
    FArrayBox lonlatFab(uBox, 2);
    panelCS->fabTransformEquiangularToLonLat(XiFab, lonlatFab);

    // longitude-latitude velocity
    FArrayBox vecRLLFab(uBox, SpaceDim);
    // Set up initial condition in this grid
    for (BoxIterator bit(uBox); bit.ok(); ++bit)
    {
      IntVect iv = bit();

      // Compute lat/lon coordinates from equiangulars.
      Real lon = lonlatFab(iv, 0);
      Real lat = lonlatFab(iv, 1);

      // Notice longitude before latitude!
      RealVect lonlatRV(lon, lat);

      // Fill in the initial conserved variables.

      Real h = (*m_height)(lonlatRV);
      RealVect vel = (*m_velocity)(lonlatRV);
      // Seems it might be more elegant to fill in the primitive variables
      // and then convert these to conserved variables,
      // but that would require calling the MOLPhysics class.
      UFab(iv, UHGT) = h;
      vecRLLFab(iv, 0) = vel[0];
      vecRLLFab(iv, 1) = vel[1];
    }
    // mapped velocity
    FArrayBox velFab(uBox, SpaceDim);
    panelCS->fabVectorTransformLatLonToEquiangular(XiFab, vecRLLFab, velFab);
    UFab.copy(UFab, UHGT, UMOMX); // UFab[UMOMX] := UFab[UHGT]
    UFab.copy(UFab, UHGT, UMOMY); // UFab[UMOMY] := UFab[UHGT]
    UFab.mult(velFab, 0, UMOMX); // UFab[UMOMX] *= velFab[0]
    UFab.mult(velFab, 1, UMOMY); // UFab[UMOMY] *= velFab[1]
    // dummy statement in order to get around gdb bug
    int dummy_unused = 0; dummy_unused = 0;
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
ShallowWaterCubedSphereIBC::
primBC(FArrayBox&             a_WGdnv,
       const FArrayBox&       a_Wextrap,
       const FArrayBox&       a_W,
       const FArrayBox *const a_unitNormalBasisPtr,
       const Interval&        a_velIntv,
       const int&             a_dir,
       const Side::LoHiSide&  a_side,
       const Real&            a_time)
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
artViscBC(FArrayBox&                                   a_NtFdir,
          const CHArray<Real, SpaceDim+1, ArRangeCol>& a_Nctg,
          const FArrayBox&                             a_U,
          const FArrayBox&                             a_unitNormalBasis,
          const FArrayBox&                             a_divVel,
          const FArrayBox&                             a_csq,
          const FArrayBox&                             a_dxFace,
          const Interval&                              a_momIntv,
          const Real                                   a_alpha,
          const Real                                   a_beta,
          const Box&                                   a_loFaceBox,
          const int                                    a_hasLo,
          const Box&                                   a_hiFaceBox,
          const int                                    a_hasHi,
          const int                                    a_dir)
{
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
