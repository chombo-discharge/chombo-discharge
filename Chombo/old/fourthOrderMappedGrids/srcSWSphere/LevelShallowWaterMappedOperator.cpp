#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>

#include "LevelShallowWaterMappedOperator.H"
#include "SWintegrator.H"
#include "ShallowWaterPhysicsF_F.H"
#include "CubedSphere2DPanelCS.H"
#include "FourthOrderUtil.H"

#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
LevelShallowWaterMappedOperator::LevelShallowWaterMappedOperator()
{
  m_defined = false;
  m_dx = 0.0;
  m_refineCoarse = 0;
  m_patchConsOperatorPtr = new PatchShallowWaterMappedOperator();
}

//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
LevelShallowWaterMappedOperator::~LevelShallowWaterMappedOperator()
{
}


//////////////////////////////////////////////////////////////////////////////
void
LevelShallowWaterMappedOperator::define(
                                  const DisjointBoxLayout&  a_thisLayout,
                                  const DisjointBoxLayout&  a_coarserLayout,
                                  const ProblemDomain&      a_domain,
                                  const int&                a_refineCoarse,
                                  const Real&               a_dx,
                                  const MOLPhysics* const a_molPhysics,
                                  const int&                a_numFields,
                                  const bool&               a_hasCoarser,
                                  const bool&               a_hasFiner)
{
  LevelConsOperator::define(a_thisLayout, a_coarserLayout,
                            a_domain, a_refineCoarse, a_dx, a_molPhysics,
                            a_numFields, a_hasCoarser, a_hasFiner);
  m_physIBCPtr = (PhysShallowWaterMappedIBC*) a_molPhysics->getPhysIBC();
}

//////////////////////////////////////////////////////////////////////////////
void
LevelShallowWaterMappedOperator::defineMetricStuff()
{
  IntVect ghostVect = m_numGhost * IntVect::Unit;
  // face-centered LevelData<FluxBox>
  m_contravariantMetricFace.define(m_grids, SpaceDim*SpaceDim, ghostVect);
  m_orthoMatrix.define(m_grids, SpaceDim*SpaceDim, ghostVect);
  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& baseBox = m_grids[dit];
      int baseBlockNum = m_coordSysPtr->whichBlock(baseBox);
      const CubedSphere2DPanelCS* panelCS =
        (CubedSphere2DPanelCS*) (m_coordSysPtr->getCoordSys(baseBlockNum));
      FluxBox& orthoMatrixFlub = m_orthoMatrix[dit];
      FluxBox& contravariantMetricFaceFlub =
        m_contravariantMetricFace[dit];
      for (int idirFace = 0; idirFace < SpaceDim; idirFace++)
        {
          FArrayBox& orthoMatrixFab = orthoMatrixFlub[idirFace];
          panelCS->getOrthonormalizingMatrix(orthoMatrixFab,
                                             orthoMatrixFab.box(),
                                             idirFace);

          FArrayBox& contravariantMetricFab =
            contravariantMetricFaceFlub[idirFace];
          int compLo = 0;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              Interval intvl(compLo, compLo+SpaceDim-1);
              FArrayBox metricFab(intvl, contravariantMetricFab);
              panelCS->contravariantMetric(metricFab, idir);
              compLo += SpaceDim;
            }
        }
    }
}


//////////////////////////////////////////////////////////////////////////////
void
LevelShallowWaterMappedOperator::evalRHS(
                                   LevelData<FArrayBox>&       a_LofU,
                                   LevelData<FArrayBox>&       a_U,
                                   LevelFluxRegister&          a_finerFluxRegister,
                                   LevelFluxRegister&          a_coarserFluxRegister,
                                   const LevelData<FArrayBox>& a_UcoarseOld,
                                   const Real&                 a_timeCoarseOld,
                                   const LevelData<FArrayBox>& a_UcoarseNew,
                                   const Real&                 a_timeCoarseNew,
                                   Real                        a_time,
                                   Real                        a_weight)
{
  // Now that we have advection velocities, do what we did before.
  LevelMappedConsOperator::evalRHS(a_LofU, a_U,
                                   a_finerFluxRegister, a_coarserFluxRegister,
                                   a_UcoarseOld, a_timeCoarseOld,
                                   a_UcoarseNew, a_timeCoarseNew,
                                   a_time, a_weight);
}

//////////////////////////////////////////////////////////////////////////////
void
LevelShallowWaterMappedOperator::evalRHSpatches(
                                                LevelData<FArrayBox>&       a_LofU,
                                                const LevelData<FArrayBox>& a_U,
                                                LevelFluxRegister&          a_finerFluxRegister,
                                                LevelFluxRegister&          a_coarserFluxRegister,
                                                Real                        a_weight)
{
  CH_TIME("LevelShallowWaterMappedOperator::evalRHSpatches");
  // Why does FfaceAvgAll have ghost cells?  I'm not sure.
  LevelData<FluxBox> FfaceAvgAll(m_grids, m_numFluxes, IntVect::Unit);
  BlockRegister blockRegister(m_coordSysPtr, m_grids, 1);
  bool setFlattening = (m_useFlattening && (m_evalCount == 1));
  PatchShallowWaterMappedOperator* patchOperatorPtr =
    (PatchShallowWaterMappedOperator*) m_patchConsOperatorPtr;
  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      setPatchIndex(dit());

      FluxBox FfaceAvg;
      FluxBox WfaceOrthoCen; // need this on block interfaces only
      patchOperatorPtr->getNormalFlux(FfaceAvg,
                                      WfaceOrthoCen,
                                      a_U[dit],
                                      a_weight,
                                      setFlattening,
                                      m_flattening[dit]);
      FfaceAvgAll[dit].copy(FfaceAvg, FfaceAvg.box());
      // WfaceOrthoCenAll[dit].copy(WfaceOrthoCen, WfaceOrthoCen.box());
      for (int idir = 0; idir < SpaceDim; idir++)
        for (SideIterator sit; sit.ok(); ++sit)
          {
            Side::LoHiSide side = sit();
            if (blockRegister.hasInterface(dit(), idir, side))
              {
                blockRegister.storeAux(WfaceOrthoCen[idir], dit(), idir, side);
              }
          }
    }
  blockRegister.close();

  /*
    Set common flux on block boundaries.
  */
  setCommonFlux(FfaceAvgAll, blockRegister);

  /*
    Update total fluxes m_fluxes from FfaceAvg:
    m_fluxes[dit] += a_weight * FfaceAvg.

    Also, if there is a finer level,
    update a_finerFluxRegister using FfaceAvg.

    Also, if there is a coarser level,
    update a_coarserFluxRegister using FfaceAvg.
  */
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& FfaceAvg = FfaceAvgAll[dit];
      updateFluxTotalsAndRegisters(FfaceAvg,
                                   a_finerFluxRegister, a_coarserFluxRegister,
                                   dit(), a_weight);

      /*
        Set LofUFab = - div(FfaceAvg).
       */
      FArrayBox& LofUFab = a_LofU[dit];
      // Set LofUFab = div(FfaceAvg).
      setPatchIndex(dit()); // need m_patchConsOperatorPtr->m_currentBox
      m_patchConsOperatorPtr->getFluxDivergence(LofUFab, FfaceAvg);
      // Actually want -div.
      LofUFab.negate();
    }
  // added 22 Oct 2008:  these change nothing
  // a_LofU.exchange();
  // a_U.exchange();
  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}


//////////////////////////////////////////////////////////////////////////////
void
LevelShallowWaterMappedOperator::setCommonFlux(LevelData<FluxBox>&   a_flux,
                                               const BlockRegister&  a_blockRegister) const
{
  CH_TIME("LevelShallowWaterMappedOperator::setCommonFlux");
  // LevelData<FluxBox> FfaceCenAll(m_grids, m_numFluxes, IntVect::Unit);
  // I think I just need to have ghost cells in BlockRegister.
  // a_flux.copyTo(FfaceCenAll);
  const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries =
    m_coordSysPtr->boundaries();
  IntVect csWcomps = IntVect(D_DECL(UMOMX, UMOMY, UMOMZ));
  IntVect orthoWcomps = IntVect(D_DECL(UMOMNORM, UMOMTAN, UMOMTAN2));
  IntVect csVcomps = IntVect(D_DECL(0, 1, 2));
  IntVect orthoVcomps = IntVect(D_DECL(0, 1, 2));

  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& FfaceAvg = a_flux[dit];
      const Box& baseBox = m_grids[dit];
      int baseBlockNum = m_coordSysPtr->whichBlock(baseBox);
      const CubedSphere2DPanelCS* panelCS =
        (CubedSphere2DPanelCS*) (m_coordSysPtr->getCoordSys(baseBlockNum));
      int faceID = 0;
      for (SideIterator sit; sit.ok(); ++sit)
        {
          Side::LoHiSide side = sit();
          Side::LoHiSide sideOther = flip(side);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (a_blockRegister.hasInterface(dit(), idir, side))
                {
                  // maybe better if this is done inside BlockRegister
                  const BlockBoundary& bb = boundaries[baseBlockNum][faceID];
                  int reorientFace = bb.reorientFace(idir);

                  IntVect faceGrowVect = IntVect::Unit;
                  faceGrowVect[idir] = 0;

                  Box faceBaseBox = adjCellBox(baseBox, idir, side, 1);
                  // Shift from cells to faces:
                  // If Lo, then shift +1/2; if Hi, then shift -1/2.
                  faceBaseBox.shiftHalf(idir, -sign(side));

                  Box faceBox = grow(faceBaseBox, faceGrowVect);

                  // Need to define these FABs.
                  FArrayBox WfaceOrthoThisFab(faceBox, WNUM);
                  FArrayBox WfaceOrthoOtherFab(faceBox, WNUM);
                  a_blockRegister.getAux(WfaceOrthoThisFab, dit(),
                                         idir, side, sideOther);
                  a_blockRegister.getAux(WfaceOrthoOtherFab, dit(),
                                         idir, side, side);
                  FArrayBox WfaceOrthoCommonFab(faceBox, WNUM);
                  FORT_SWRIEMANNROEF(CHF_FRA(WfaceOrthoCommonFab),
                                     CHF_CONST_FRA(WfaceOrthoThisFab),
                                     CHF_CONST_FRA(WfaceOrthoOtherFab),
                                     CHF_CONST_INT(idir),
                                     CHF_BOX(faceBox));

                  // WfaceOrthoCommonFab = (h, vnorm, vtan).

                  // Note that in the orthonormal frame, the problem of
                  // computing the flux is 1D (normal component only),
                  // and should be treated as such.

                  // fluxOrthoThisFab =
                  // (h*vnorm, h*vnorm*vnorm+G*h*h/2, h*vnorm*vtan).
                  // These are the normal components of the fluxes.
                  int idirFake = 0; // normal direction (vs. 1 for tangential)
                  FArrayBox fluxOrthoFab(faceBox, m_numFluxes);
                  FORT_SWGETORTHOFLUXF(CHF_FRA(fluxOrthoFab),
                                       CHF_CONST_FRA(WfaceOrthoCommonFab),
                                       CHF_CONST_INT(idirFake),
                                       CHF_BOX(faceBox));

                  // Fill xiFab with mapped coordinates of face centers.
                  FArrayBox xiFab(faceBox, SpaceDim);
                  panelCS->getCenterMappedCoordinates(xiFab, faceBox);
                  // Get J == sqrt(det(covariant metric)) at face centers.
                  FArrayBox JFab(faceBox, 1);
                  panelCS->pointwiseJ(JFab, xiFab, faceBox);

                  // Convert the momentum fluxes to this block.
                  // Note that even after we do this, we'll still need to
                  // deorthonormalize the components separately.
                  FArrayBox fluxMomOrthoFab(faceBox, SpaceDim);
                  // fluxMomOrthoFab[0:SpaceDim-1] := fluxOrthoFab[1:SpaceDim]
                  fluxMomOrthoFab.copy(fluxOrthoFab, 1, 0, SpaceDim);
                  IntVect orthoVmomcomps = IntVect(D_DECL(UMOMX, UMOMY, UMOMZ));
                  panelCS->deorthonormalize(fluxOrthoFab,
                                            fluxMomOrthoFab,
                                            faceBox,
                                            idir, orthoVmomcomps, csVcomps);
                  // fluxOrthoFab[1:SpaceDim] := fluxMomOrthoFab[0:SpaceDim-1]
                  fluxOrthoFab.copy(fluxMomOrthoFab, 0, 1, SpaceDim);
                  // We'll still need to deorthonormalize the components
                  // of fluxOrthoFab separately.

                  FArrayBox fluxCenThisFab(faceBox, m_numFluxes);
                  // For each of the m_numFluxes components of the flux,
                  // deorthonormalize and then find component idir.
                  for (int comp = 0; comp < m_numFluxes; comp++)
                    {
                      FArrayBox fluxCompOrthoFab(faceBox, SpaceDim);
                      // Set fluxCompOrthoFab[0] := fluxOrthoFab[comp]
                      // and other components of fluxCompOrthoFab to zero.
                      fluxCompOrthoFab.copy(fluxOrthoFab, comp, 0);
                      fluxCompOrthoFab.setVal(0., faceBox, 1, SpaceDim-1);
                      FArrayBox fluxCompThisFab(faceBox, SpaceDim);
                      panelCS->deorthonormalize(fluxCompOrthoFab,
                                                fluxCompThisFab,
                                                faceBox,
                                                idir, orthoVcomps, csVcomps);
                      // Set fluxCenThisFab[comp] := fluxCompThisFab[idir].
                      // The other SpaceDim-1 components of fluxCompThisFab
                      // are ignored.
                      fluxCenThisFab.copy(fluxCompThisFab, idir, comp);
                      // fluxCenThisFab[comp] *= JFab[0]
                      fluxCenThisFab.mult(JFab, 0, comp);
                    }
                  // fluxCenThisFab.mult(reorientFace);

                  // Now convolve face-centered fluxCenThisFab
                  // to get face-averaged fluxAvgThisFab.
                  FArrayBox fluxAvgThisFab(faceBaseBox, m_numFluxes);
                  fourthOrderAverageCenterFace(fluxAvgThisFab,
                                               fluxCenThisFab,
                                               idir);

                  FfaceAvg[idir].copy(fluxAvgThisFab);

                  // dummy statement in order to get around gdb bug
                  int dummy_unused = 0; dummy_unused = 0;
                }
              faceID++;
            } // iterate over dimensions
        } // iterate over sides
    } // iterate over patches

  // Now I think we need to multiply a_flux by area of each face.
}



//////////////////////////////////////////////////////////////////////////////
void
LevelShallowWaterMappedOperator::setPatchIndex(const DataIndex&  a_ind) const
{
  // This function must be const because cellUJToCellU must be const.
  // ugh, need to do this for the function to be const
  PatchShallowWaterMappedOperator& patchConsOperator =
    (PatchShallowWaterMappedOperator&) *m_patchConsOperatorPtr;

  const Box& curBox = m_grids[a_ind];
  patchConsOperator.setCurrentBox(curBox);
  patchConsOperator.setContravariantMetricFace(&(m_contravariantMetricFace[a_ind]));
  patchConsOperator.setOrthoMatrix(&(m_orthoMatrix[a_ind]));
}

//////////////////////////////////////////////////////////////////////////////
void
LevelShallowWaterMappedOperator::exchangeGhosts(LevelData<FArrayBox>&   a_U)
{
  LevelConsOperator::exchangeGhosts(a_U);
  Interval intvlScalar(UHGT, UHGT);
  m_mblexPtr->interpGhosts(a_U, intvlScalar);
  Interval intvlVector(UMOMX, UMOMX+SpaceDim-1);
  m_mblexPtr->interpGhostsVector(a_U, intvlVector);
}

#include "NamespaceFooter.H"
