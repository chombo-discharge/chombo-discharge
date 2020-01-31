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

#include "LevelMappedConsOperator.H"
#include "PatchMappedConsOperator.H"
#include "FourthOrderUtil.H" // functions with names beginning "fourthOrder"
#include "CellToEdge.H"

#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
LevelMappedConsOperator::LevelMappedConsOperator()
{
  m_defined = false;
  m_dx = 0.0;
  m_refineCoarse = 0;
  m_patchConsOperatorPtr = new PatchMappedConsOperator();
  m_coordSysPtr = NULL;
  m_mblexPtr = NULL;
}

//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
LevelMappedConsOperator::~LevelMappedConsOperator()
{
  // delete m_patchConsOperatorPtr;
}


//////////////////////////////////////////////////////////////////////////////
void
LevelMappedConsOperator::evalRHS(
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
  CH_TIME("LevelMappedConsOperator::evalRHS");
  const IntVect& ghostVect = a_U.ghostVect();
  int ncomp = a_U.nComp();
  LevelData<FArrayBox> Uonly(m_grids, ncomp, ghostVect);
  a_U.exchange(); // added by petermc, 10 Feb 2011. Exchange <UJ> WITHIN blocks.
  // a_U holds <UJ>; set Uonly to <U>.
  cellUJToCellU(Uonly, a_U);
  LevelConsOperator::evalRHS(a_LofU, Uonly,
                             a_finerFluxRegister, a_coarserFluxRegister,
                             a_UcoarseOld, a_timeCoarseOld,
                             a_UcoarseNew, a_timeCoarseNew,
                             a_time, a_weight);

  if (m_useSourceTerm)
    {
      // a_LofJU += <J * S(Uonly)>
      m_sourceTermPtr->addSourceTerm(a_LofU, Uonly);
    }
}


//////////////////////////////////////////////////////////////////////////////
void
LevelMappedConsOperator::evalRHSpatches(
                                        LevelData<FArrayBox>&       a_LofU,
                                        const LevelData<FArrayBox>& a_U,
                                        LevelFluxRegister&          a_finerFluxRegister,
                                        LevelFluxRegister&          a_coarserFluxRegister,
                                        Real                        a_weight)
{
  CH_TIME("LevelMappedConsOperator::evalRHSpatches");
  // int numFluxesPerField = SpaceDim; // for Cartesian, this is 1
  LevelData<FluxBox> FfaceAvgAll(m_grids,
                                 m_numFluxes, // * numFluxesPerField,
                                 IntVect::Unit);
  BlockRegister blockRegister(m_coordSysPtr, m_grids, 0);
  bool setFlattening = (m_useFlattening && (m_evalCount == 1));
  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      setPatchIndex(dit());

      FluxBox FfaceAvg;
      m_patchConsOperatorPtr->getNormalFlux(FfaceAvg,
                                            a_U[dit],
                                            a_weight,
                                            setFlattening,
                                            m_flattening[dit]);
      FfaceAvgAll[dit].copy(FfaceAvg, FfaceAvg.box());
      for (int idir = 0; idir < SpaceDim; idir++)
        for (SideIterator sit; sit.ok(); ++sit)
          {
            Side::LoHiSide side = sit();
            if (blockRegister.hasInterface(dit(), idir, side))
              {
                blockRegister.storeFlux(FfaceAvg[idir], dit(), idir, side);
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
LevelMappedConsOperator::setCommonFlux(LevelData<FluxBox>&   a_flux,
                                       const BlockRegister&  a_blockRegister) const
{
  CH_TIME("LevelMappedConsOperator::setCommonFlux");
  const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries =
    m_coordSysPtr->boundaries();
  DataIterator dit = m_grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& FfaceAvg = a_flux[dit];
      const Box& baseBox = m_grids[dit];
      int baseBlockNum = m_coordSysPtr->whichBlock(baseBox);
      int faceID = 0;
      for (SideIterator sit; sit.ok(); ++sit)
        {
          Side::LoHiSide side = sit();
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (a_blockRegister.hasInterface(dit(), idir, side))
                {
                  // maybe better if this is done inside BlockRegister
                  const BlockBoundary& bb = boundaries[baseBlockNum][faceID];
                  int reorientFace = bb.reorientFace(idir);
                  Box faceBox = adjCellBox(baseBox, idir, side, 1);
                  // faceBox.grow(faceGrowVect);
                  // if Lo, then shift +1; if Hi, then shift -1
                  faceBox.shiftHalf(idir, -sign(side));
                  Side::LoHiSide sideOther = flip(side);
                  // Need to define these FABs.
                  FArrayBox fluxThisFab(faceBox, m_numFluxes);
                  FArrayBox fluxOtherFab(faceBox, m_numFluxes);
                  a_blockRegister.getFlux(fluxThisFab, dit(),
                                          idir, side, side);
                  fluxThisFab.mult(reorientFace * 0.5);
                  a_blockRegister.getFlux(fluxOtherFab, dit(),
                                          idir, side, sideOther);
                  fluxOtherFab.mult(0.5);
                  fluxThisFab += fluxOtherFab;
                  FfaceAvg[idir].copy(fluxThisFab);
                }
              faceID++;
            } // iterate over dimensions
        } // iterate over sides
    } // iterate over patches
}


//////////////////////////////////////////////////////////////////////////////
void
LevelMappedConsOperator::exchangeGhosts(LevelData<FArrayBox>&   a_U)
{
  LevelConsOperator::exchangeGhosts(a_U);
  m_mblexPtr->interpGhosts(a_U);
}


//////////////////////////////////////////////////////////////////////////////
void
LevelMappedConsOperator::setCoordSys(MultiBlockCoordSys* a_coordSysPtr)
{
  m_coordSysPtr = a_coordSysPtr;
  PatchMappedConsOperator& patchConsOperator =
    (PatchMappedConsOperator&) *m_patchConsOperatorPtr;
  patchConsOperator.setCoordSys(a_coordSysPtr);
}

//////////////////////////////////////////////////////////////////////////////
void
LevelMappedConsOperator::setLevelExchange(MultiBlockLevelExchangeAverage* a_mblexPtr)
{
  m_mblexPtr = a_mblexPtr;
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::useSourceTerm(bool a_useSourceTerm)
{
  m_useSourceTerm = a_useSourceTerm;
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::setSourceTerm(LevelSourceTerm* a_sourceTermPtr)
{
  m_sourceTermPtr = a_sourceTermPtr;
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::computeFaceAverages(
                                             LevelData<FluxBox>& a_face_data,
                                             const LevelData<FArrayBox>& a_cell_data ) const
{
  if (m_spaceOrder == 4)
    {
      fourthOrderCellToFace( a_face_data, a_cell_data );
    }
  else if (m_spaceOrder == 2)
    {
      CellToEdge( a_cell_data, a_face_data );
    }
  else
    {
      MayDay::Error("Bad Space Order in LevelMappedConsOperator");
    }
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::computeFaceCenters(
                                            LevelData<FluxBox>& a_face_data,
                                            const LevelData<FArrayBox>& a_cell_data ) const
{
  if (m_spaceOrder == 4)
    {
      fourthOrderCellToFaceCenters( a_face_data, a_cell_data );
    }
  else if (m_spaceOrder == 2)
    {
      CellToEdge( a_cell_data, a_face_data );
    }
  else
    {
      MayDay::Error("Bad Space Order in LevelMappedConsOperator");
    }
}

//////////////////////////////////////////////////////////////////////////////

inline void
LevelMappedConsOperator::getPhysicalCellVolumes(
                                                LevelData<FArrayBox>& a_volumes ) const
{
  // Do I need this?
  // a_volumes.define( m_coordSysPtr->getCellVolumes() );
}

//////////////////////////////////////////////////////////////////////////////

void
LevelMappedConsOperator::computeCompFaceFluxes( LevelData<FluxBox>& a_uTimesV,
                                          const LevelData<FluxBox>& a_u,
                                          const LevelData<FluxBox>& a_v) const
{
  // Compute the SpaceDim-by-nComp face-averaged fluxes in computational
  // space, where a_v is the SpaceDim-dimensional velocity vector and
  // a_u is the nComp-dim state vector
  int ncomp = a_u.nComp();
  CH_assert(a_v.nComp() == SpaceDim);
  CH_assert(a_uTimesV.nComp() == SpaceDim * ncomp);
  PatchMappedConsOperator& patchConsOperator =
    (PatchMappedConsOperator&) *m_patchConsOperatorPtr;
  // loop over boxes
  DataIterator dit = a_u.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FluxBox& thisUV = a_uTimesV[dit];
      const FluxBox& thisU = a_u[dit];
      const FluxBox& thisV = a_v[dit];
      patchConsOperator.computeCompFaceFluxes(thisUV, thisU, thisV);
    } // end loop over boxes
}

//////////////////////////////////////////////////////////////////////////////
void
LevelMappedConsOperator::cellUJToCellU(LevelData<FArrayBox>& a_Uavg,
                                       const LevelData<FArrayBox>& a_UJavg) const
{
  CH_TIME("LevelMappedConsOperator::cellUJToCellU");
  PatchMappedConsOperator& patchConsOperator =
    (PatchMappedConsOperator&) *m_patchConsOperatorPtr;
  // loop over boxes
  DataIterator dit = a_UJavg.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      setPatchIndex(dit());

      FArrayBox& thisU = a_Uavg[dit];
      const FArrayBox& thisUJ = a_UJavg[dit];
      patchConsOperator.cellUJToCellU(thisU, thisUJ);
    } // end loop over boxes
}

//////////////////////////////////////////////////////////////////////////////
void
LevelMappedConsOperator::cellUToCellUJ(LevelData<FArrayBox>& a_UJavg,
                                       const LevelData<FArrayBox>& a_Uavg) const
{
  CH_TIME("LevelMappedConsOperator::cellUToCellUJ");
  PatchMappedConsOperator& patchConsOperator =
    (PatchMappedConsOperator&) *m_patchConsOperatorPtr;
  // loop over boxes
  DataIterator dit = a_Uavg.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      setPatchIndex(dit());

      FArrayBox& thisUJ = a_UJavg[dit];
      const FArrayBox& thisU = a_Uavg[dit];
      patchConsOperator.cellUToCellUJ(thisUJ, thisU);
    } // end loop over boxes
}

#include "NamespaceFooter.H"
