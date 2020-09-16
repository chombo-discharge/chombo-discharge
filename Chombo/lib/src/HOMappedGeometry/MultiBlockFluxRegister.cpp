#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MultiBlockFluxRegister.H"
#include "LevelFluxRegisterF_F.H"
#include "CubedSphereShellPanelCS.H" // FIXME
#include "MBMiscUtil.H"
#include "MOLUtilities.H"
#include "CH_assert.H"
#include "DebugOut.H"
#include <map>
#include "NamespaceHeader.H"
using namespace std;

//-----------------------------------------------------------------------
MultiBlockFluxRegister::
MultiBlockFluxRegister(MultiBlockCoordSys* a_coordSysCoarse,
                       MultiBlockCoordSys* a_coordSysFine) :
  LevelFluxRegister(),
  m_coordSysCoarse(a_coordSysCoarse),
  m_coordSysFine(a_coordSysFine)
//  m_multiblockFineFlux(),
//  m_neighborFineFlux()
{
  CH_assert(a_coordSysCoarse != NULL);
  CH_assert(a_coordSysFine != NULL);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
MultiBlockFluxRegister::
~MultiBlockFluxRegister()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
MultiBlockFluxRegister::
define(const DisjointBoxLayout& a_dblFine,
       const DisjointBoxLayout& a_dblCoarse,
       const ProblemDomain&     a_dProblem,
       int                      a_nRefine,
       int                      a_nComp,
       bool                     a_scaleFineFluxes,
       int                      a_unrefDir)
{
  CH_TIME("MultiBlockFluxRegister::define");
  // Do all the intrablock stuff.
  CH_assert((a_unrefDir >= -1) && (a_unrefDir < SpaceDim));
  LevelFluxRegister::define(a_dblFine, a_dblCoarse, a_dProblem, a_nRefine,
                            a_nComp, a_scaleFineFluxes, a_unrefDir);

  const Vector<Box>& blocks = m_coordSysCoarse->mappingBlocks();

  // Iterate over all boxes and gather information on those whose boundaries coincide
  // with block boundaries.
  DisjointBoxLayout coarsenedFineGrids;
  coarsen(coarsenedFineGrids, a_dblFine, m_nRefine);
  const Vector<Tuple<BlockBoundary, 2*SpaceDim> >& blockBoundaries =
    m_coordSysCoarse->boundaries();

  m_blockNum.define(coarsenedFineGrids);
  for (DataIterator ditCF(coarsenedFineGrids.dataIterator()); ditCF.ok(); ++ditCF)
    {
      const Box& bxCF = coarsenedFineGrids[ditCF];
      int blockID = m_coordSysCoarse->whichBlock(bxCF);
      // CH_assert(blockID != -1);
      m_blockNum[ditCF] = blockID;
    }

  // We construct data structures for accumulating fine fluxes from remote
  // blocks in each direction.
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      // Don't need them in the unrefined direction
      if (dir == m_unrefinedDirection)
        continue;
      for (DataIterator coarseIter(a_dblCoarse.dataIterator()); coarseIter.ok(); ++coarseIter)
        {
          // Find the block containing this box.
          Box fluxBox = a_dblCoarse[coarseIter];
          //cout << "fluxBox = " << fluxBox << endl;
          int blockID = m_coordSysCoarse->whichBlock(fluxBox);
          CH_assert(blockID != -1);
          const Box& blockBox = blocks[blockID];
          //cout << "  Block box = " << blockBox << endl;
          
          // Now find out whether this flux abuts the block boundary and set up
          // any block-block flux transfers that are needed.
          Box fluxBox1 = fluxBox;
          fluxBox1.grow(dir, 1);
          //**CHECKME Commented out --- why would this ever not be true??
          // if (blockBox.intersectsNotEmpty(fluxBox1))
          // {
          //cout << fluxBox1 << endl;
          if (fluxBox1.smallEnd(dir) < blockBox.smallEnd(dir))
            {
              // Associate this box with the (low side) block boundary.
              const BlockBoundary& boundary =
                blockBoundaries[blockID][dir];
              const int nbrBlockID = boundary.neighbor();
              const IndicesTransformation& transformation =
                boundary.getTransformation();
              
              // Find any fine boxes on the far side of the boundary.
              // FIXME: This is slow and will have to be revisited at some point.
              Vector<Box>& coarseLocations =
                m_coarseLocations[dir+SpaceDim][coarseIter];
              //**FIXME You can probably do some simple sorting optimizations based
              //**      on the corners of the neighbour box
              for (LayoutIterator fineIter = coarsenedFineGrids.layoutIterator();
                   fineIter.ok(); ++fineIter)
                {
                  Box fineBox = coarsenedFineGrids[fineIter];
                  //cout << fineBox << endl;
                  
                  // Only consider boxes from the neighboring block
                  if (m_coordSysCoarse->whichBlock(fineBox) != nbrBlockID) continue;
                  
                  fineBox = transformation.transformBack(fineBox);
                  if (fineBox.intersectsNotEmpty(fluxBox1))
                    {
                      Box intersect = fineBox & fluxBox1;
                      intersect = adjCellHi(intersect, dir);
                      // Add this intersection to the set of coarse flux locations so
                      // that we can increment coarse fluxes here.
                      coarseLocations.push_back(intersect);
                    }
                }
            }
          if (fluxBox1.bigEnd(dir) > blockBox.bigEnd(dir))
            {
              // Associate this box with the (hi side) block boundary.
              const BlockBoundary& boundary =
                blockBoundaries[blockID][dir + SpaceDim];
              const int nbrBlockID = boundary.neighbor();
              const IndicesTransformation& transformation =
                boundary.getTransformation();
              
              // Find any fine boxes on the far side of the boundary.
              // FIXME: This is slow and will have to be revisited at some point.
              Vector<Box>& coarseLocations =
                m_coarseLocations[dir][coarseIter];
              //**FIXME You can probably do some simple sorting optimizations based
              //**      on the corners of the neighbour box
              for (LayoutIterator fineIter = coarsenedFineGrids.layoutIterator();
                   fineIter.ok(); ++fineIter)
                {
                  Box fineBox = coarsenedFineGrids[fineIter];
                  //cout << "finebox = " << fineBox << endl;
                  
                  // Only consider boxes from the neighboring block
                  if (m_coordSysCoarse->whichBlock(fineBox) != nbrBlockID) continue;
                  
                  fineBox = transformation.transformBack(fineBox);
                  //cout << "finebox' = " << fineBox << endl;
                  if (fineBox.intersectsNotEmpty(fluxBox1))
                    {
                      Box intersect = fineBox & fluxBox1;
                      intersect = adjCellLo(intersect, dir);
                      
                      // Add this intersection to the set of coarse flux locations so
                      // that we can increment coarse fluxes here.
                      coarseLocations.push_back(intersect);
                    }
                }
            }
        }
      Vector<Box> sourceBoxes, destBoxes;
      Vector<int> srcProcs, destProcs;
      for (LayoutIterator coarseIter(a_dblCoarse.layoutIterator()); coarseIter.ok(); ++coarseIter)
        {
          // Find the block containing this box.
          Box fluxBox = a_dblCoarse[coarseIter];
          //cout << "fluxBox = " << fluxBox << endl;
          int destProc = a_dblCoarse.procID(coarseIter());
          int blockID = m_coordSysCoarse->whichBlock(fluxBox);
          CH_assert(blockID != -1);
          const Box& blockBox = blocks[blockID];
          //cout << "  Block box = " << blockBox << endl;
          
          // Now find out whether this flux abuts the block boundary and set up
          // any block-block flux transfers that are needed.
          Box fluxBox1 = fluxBox;
          fluxBox1.grow(dir, 1);
          //**CHECKME Commented out --- why would this ever not be true??
          // if (blockBox.intersectsNotEmpty(fluxBox1))
          // {
          //cout << fluxBox1 << endl;
          if (fluxBox1.smallEnd(dir) < blockBox.smallEnd(dir))
            {
              // Associate this box with the (low side) block boundary.
              const BlockBoundary& boundary = blockBoundaries[blockID][dir];
              const int nbrBlockID = boundary.neighbor();
              const IndicesTransformation& transformation = boundary.getTransformation();
              
              // Find any fine boxes on the far side of the boundary.
              // FIXME: This is slow and will have to be revisited at some point.
              //**FIXME You can probably do some simple sorting optimizations based
              //**      on the corners of the neighbour box
              for (LayoutIterator fineIter = coarsenedFineGrids.layoutIterator();
                   fineIter.ok(); ++fineIter)
                {
                  Box fineBox = coarsenedFineGrids[fineIter];
                  int srcProc = coarsenedFineGrids.procID(fineIter());
                  //cout << fineBox << endl;
                  
                  // Only consider boxes from the neighboring block
                  if (m_coordSysCoarse->whichBlock(fineBox) != nbrBlockID) continue;
                  
                  fineBox = transformation.transformBack(fineBox);
                  if (fineBox.intersectsNotEmpty(fluxBox1))
                    {
                      Box intersect = fineBox & fluxBox1;
                      intersect = adjCellHi(intersect, dir);
                      destBoxes.push_back(intersect);
                      //cout << "intersect: " << intersect << endl;
                      Box otherSide = transformation.transformFwd(intersect);
                      //              Box otherSide = transformation.transformFwd(adjCellLo(intersect, dir));
                      //cout << "otherSide: " << otherSide << endl;
                      sourceBoxes.push_back(otherSide);

                      srcProcs.push_back(srcProc);
                      destProcs.push_back(destProc);
                    }
                }
            }
          if (fluxBox1.bigEnd(dir) > blockBox.bigEnd(dir))
            {
              // Associate this box with the (hi side) block boundary.
              const BlockBoundary& boundary = blockBoundaries[blockID][dir + SpaceDim];
              const int nbrBlockID = boundary.neighbor();
              const IndicesTransformation& transformation = boundary.getTransformation();
              
              // Find any fine boxes on the far side of the boundary.
              // FIXME: This is slow and will have to be revisited at some point.
              //**FIXME You can probably do some simple sorting optimizations based
              //**      on the corners of the neighbour box
              for (LayoutIterator fineIter = coarsenedFineGrids.layoutIterator();
                   fineIter.ok(); ++fineIter)
                {
                  Box fineBox = coarsenedFineGrids[fineIter];
                  int srcProc = coarsenedFineGrids.procID(fineIter());
                  //cout << "finebox = " << fineBox << endl;
                  
                  // Only consider boxes from the neighboring block
                  if (m_coordSysCoarse->whichBlock(fineBox) != nbrBlockID) continue;
                  
                  fineBox = transformation.transformBack(fineBox);
                  //cout << "finebox' = " << fineBox << endl;
                  if (fineBox.intersectsNotEmpty(fluxBox1))
                    {
                      Box intersect = fineBox & fluxBox1;
                      intersect = adjCellLo(intersect, dir);
                      destBoxes.push_back(intersect);
                      //cout << "intersect: " << intersect << endl;
                      Box otherSide = transformation.transformFwd(intersect);
                      //              Box otherSide = transformation.transformFwd(adjCellHi(intersect, dir));
                      //cout << "otherSide: " << otherSide << endl;
                      sourceBoxes.push_back(otherSide);

                      srcProcs.push_back(srcProc);
                      destProcs.push_back(destProc);
                    }
                }
            }
        }

      // Now that we have source/destination boxes for coarsened fine fluxes
      // in this direction, we create a copier that will traffick the data
      // from one to the other.
      map<Box, Box> sourceToDestMapping;
      for (int b = 0; b < sourceBoxes.size(); ++b)
        sourceToDestMapping[sourceBoxes[b]] = destBoxes[b];
      DisjointBoxLayout src(sourceBoxes, srcProcs);
      DisjointBoxLayout dest(destBoxes, destProcs);
      m_remoteCopiers[dir].define(m_fineFlux.disjointBoxLayout(), // a_from
                                  m_coarFlux.disjointBoxLayout(), // a_to
                                  src, // src
                                  dest, // dest
                                  sourceToDestMapping,
                                  IntVect::Zero);
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
class MBAddOp: public LDOperator<FArrayBox>
{
public:

//--Constructors/destructors

  // Construct
  MBAddOp(const int                       a_dir,
          const MultiBlockCoordSys* const a_coordSysPtr,
          const Interval&                 a_destVectorIntv,
          const Real                      a_scale = 1.)
    :
    m_dir(a_dir),
    m_coordSysPtr(a_coordSysPtr),
    m_destVectorIntv(a_destVectorIntv),
    m_scale(a_scale)
  { }

  // Use synthesized destructor, copy, and assignment

//--Member functions

  /// Operation for distributed memory (src is the buffer)
  virtual void linearIn(FArrayBox& arg,
                        void* buf,
                        const Box& regionTo,
                        const Interval& comps) const
    {
      Real* buffer = (Real*)buf;

      IndicesTransformation itBack = backTransform(arg, regionTo);
      Box regionFrom = itBack.transformFwd(regionTo);

      FArrayBox tmp(regionFrom, comps.size(), buffer);
      op(arg, regionFrom, comps, regionTo, tmp, comps);
    }

  /// Operation for local memory
  void op(FArrayBox&       a_dest,
          const Box&       a_regionFrom,
          const Interval&  a_Cdest,
          const Box&       a_regionTo,
          const FArrayBox& a_src,
          const Interval&  a_Csrc) const
    {
      CH_assert(a_regionFrom.numPts() == a_regionTo.numPts());
      CH_assert(a_Csrc.size() == a_Cdest.size());

      IndicesTransformation itBack = backTransform(a_dest, a_regionTo);

      // Note: we don't have to worry about the sign of the fluxes because
      // the "fluxes" have already been through the divergence operator.  The
      // data here is ready to be added to cell averages.

      int baseCdest = a_Cdest.begin();
      int baseCsrc  = a_Csrc.begin();
      // Loop over all components
      const int endCdest = a_Cdest.end();
      for (; baseCdest <= endCdest; ++baseCdest, ++baseCsrc)
        {
          for (BoxIterator bit(a_regionTo); bit.ok(); ++bit)
            {
              IntVect ivTo = bit();
              IntVect ivFrom = itBack.transformFwd(ivTo);
              a_dest(ivTo, baseCdest) += m_scale * a_src(ivFrom, baseCsrc);
            }
        }
    }

//--Data members

private:

  IndicesTransformation backTransform(const FArrayBox&  a_destFab,
                                      const Box& a_regionTo) const
  {
    // Call whichBlockOverlap() on the FABs' boxes,
    // instead of whichBlock() on a_regionFrom and a_regionTo,
    // because we're going to have cells OUTSIDE valid blocks.
    // int blockFrom = m_coordSysPtr->whichBlockOverlap(a_src.box());
    int blockTo = m_coordSysPtr->whichBlockOverlap(a_destFab.box());
    // petermc, 14 Nov 2011: Do not obtain the IndicesTransformation
    // this way, because there may be more than one way to get from
    // source block to dest block!
    // const IndicesTransformation& it =
    // m_coordSysPtr->blockTransformation(blockFrom, blockTo);
    // petermc, 14 Nov 2011: Use the fact that a_regionFrom lies just
    // outside the source block.
    Box blockToBox = m_coordSysPtr->mappingBlocks()[blockTo];
    // blockToBoxLo and blockToBoxHi are just inside the dest block.
    Box blockToBoxLo = adjCellLo(blockToBox, m_dir, -1);
    Box blockToBoxHi = adjCellHi(blockToBox, m_dir, -1);

    int faceID = m_dir;
    if (blockToBoxHi.contains(a_regionTo))
      faceID += SpaceDim;
    else
      CH_assert(blockToBoxLo.contains(a_regionTo));

    const BlockBoundary& bb = m_coordSysPtr->boundaries()[blockTo][faceID];
    IndicesTransformation itBack = bb.getTransformation();
    
    return itBack;
  }

  int m_dir;
  const MultiBlockCoordSys* m_coordSysPtr;
  const Interval m_destVectorIntv;

public:
  Real m_scale;                       ///< Extra scaling beyond what is normally
                                      ///< needed to integrated the fine fluxes
};
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
MultiBlockFluxRegister::
reflux(LevelData<FArrayBox>& a_uCoarse,
       const Interval&       a_coarseVectorIntv,
       Real                  a_scale)
{
  Interval interval(0, a_uCoarse.nComp()-1);
  reflux(a_uCoarse, interval, a_coarseVectorIntv, interval, a_scale);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
MultiBlockFluxRegister::
reflux(LevelData<FArrayBox>& a_uCoarse,
       const Interval& a_coarse_interval,
       const Interval& a_coarseVectorIntv,
       const Interval& a_flux_interval,
       Real a_scale)
{
  CH_TIME("MultiBlockFluxRegister::reflux");
  // Do the intrablock reflux and subtract the coarse part of the defect.
  LevelFluxRegister::reflux(a_uCoarse, a_coarse_interval, a_flux_interval, a_scale);

  // Copy all the fine fluxes from their "native" locations to those abutting the coarse fluxes.
  for (int dir = 0; dir < SpaceDim; ++dir)
  {
    if (m_remoteCopiers[dir].isDefined())
    {
      // Copy the coarsened fine fluxes to the corresponding coarse fluxes on their
      // neighboring blocks.
      MBAddOp op(dir, m_coordSysCoarse, a_coarseVectorIntv, -a_scale);
      m_fineFlux.copyTo(a_flux_interval, a_uCoarse, a_coarse_interval,
                        m_remoteCopiers[dir], op);

    }
  }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
MultiBlockFluxRegister::
incrementFineCentered(const FArrayBox& a_fineFlux,
                      Real a_scale,
                      const DataIndex& a_fineDataIndex,
                      const Interval& a_vectorFluxIntv,
                      const Interval& a_srcInterval,
                      const Interval& a_dstInterval,
                      int a_dir)
{
  incrementFineCentered(a_fineFlux, a_scale, a_fineDataIndex,
                        a_vectorFluxIntv, a_srcInterval, a_dstInterval,
                        a_dir, Side::Lo);
  incrementFineCentered(a_fineFlux, a_scale, a_fineDataIndex,
                        a_vectorFluxIntv, a_srcInterval, a_dstInterval,
                        a_dir, Side::Hi);
}

//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
MultiBlockFluxRegister::
incrementFineCentered(const FArrayBox& a_fineFlux,
                      Real a_scale,
                      const DataIndex& a_fineDataIndex,
                      const Interval& a_vectorFluxIntv,
                      const Interval& a_srcInterval,
                      const Interval& a_dstInterval,
                      int a_dir,
                      Side::LoHiSide a_sd)
{
  CH_TIME("MultiBlockFluxRegister::incrementFineCentered");
  FArrayBox fineFluxSaveFab;
  if (a_vectorFluxIntv.size() > 0)
    { // Copy m_fineFlux[a_fineDataIndex] because it is modified by
      // LevelFluxRegister::incrementFine(), and we want
      // to UNDO that modification in certain places.
      const FArrayBox& fineFluxFab = m_fineFlux[a_fineDataIndex];
      fineFluxSaveFab.define(fineFluxFab.box(), fineFluxFab.nComp());
      fineFluxSaveFab.copy(fineFluxFab);
    }
  LevelFluxRegister::incrementFine(a_fineFlux,
                                   a_scale, a_fineDataIndex,
                                   a_srcInterval, a_dstInterval, a_dir, a_sd);

  if (a_vectorFluxIntv.size() > 0)
    {
      int blockID = m_blockNum[a_fineDataIndex];
      int faceID = a_dir;
      if (a_sd == Side::Hi) faceID += SpaceDim;
      const BlockBoundary& bb = m_coordSysFine->boundaries()[blockID][faceID];
      if ( !bb.isDomainBoundary() )
        {
          // Box blockCoarseBox = m_coordSys->mappingBlocks()[blockID];
          // Box blockFineBox = refine(blockCoarseBox, m_nRefine);
          Box blockFineBox = m_coordSysFine->mappingBlocks()[blockID];
          Box adjFineCellBlockBox = adjCellBox(blockFineBox, a_dir, a_sd, 1);

          int shiftOut = sign(a_sd);
          int shiftIn = -sign(a_sd);

          const Box& fineFluxFaceBox = a_fineFlux.box();
          Box fineFluxCellBox(fineFluxFaceBox);
          fineFluxCellBox.shiftHalf(a_dir, shiftOut);
          Box adjFineFluxBox(fineFluxCellBox);
          adjFineFluxBox &= adjFineCellBlockBox;
          if ( !adjFineFluxBox.isEmpty() )
            { // Yes: On block-boundary face.  Need to readjust m_fineFlux.
              
              Box clipBox = m_fineFlux.box(a_fineDataIndex); // NO ghosts
              clipBox.refine(m_nRefine);
              Box fineCellBox = adjCellBox(clipBox, a_dir, a_sd, 1);
              fineCellBox &= fineFluxCellBox;
              // fineCellBox is just outside block, and has no ghosts.

              // cFine lives on the COARSENED cell-centered fine patch.
              FArrayBox& cFine = m_fineFlux[a_fineDataIndex];

              /*
                First replace vector components of
                cFine = m_fineFlux[a_fineDataIndex]
                on block boundary
                with original before calling LevelFluxRegister::incrementFine().
              */
              Real denom = 1.0;
              if (m_scaleFineFluxes)
                {
                  denom *= m_nRefine.product() / m_nRefine[a_dir];
                }
              Real scale = sign(a_sd)*a_scale/denom;
             
              const IntVect& iv = fineCellBox.smallEnd();
              IntVect civ = coarsen(iv, m_nRefine);

              // a_vectorFluxIntv should be within a_srcInterval.
              CH_assert(a_srcInterval.contains(a_vectorFluxIntv.begin()));
              CH_assert(a_srcInterval.contains(a_vectorFluxIntv.end()));
              // first of all comps in source, dest
              int srcComp = a_srcInterval.begin();
              int destComp = a_dstInterval.begin();
              // first of the vector comps in source, dest
              int srcVecComp = a_vectorFluxIntv.begin();
              int destVecComp = (srcVecComp - srcComp) + destComp;
              int ncomp = SpaceDim; // vector components only
              IntVect ivLo = iv;
              // Shift ivLo from cell-centered just outside block
              // to face-centered on block boundary (like a_fineFlux).
              if (a_sd == Side::Lo) ivLo += BASISV(a_dir);
              //              FORT_INCREMENTFINE(CHF_CONST_FRA_SHIFT(a_fineFlux, ivLo),
              //                                 CHF_FRA_SHIFT(cFine, civ),
              //                                 CHF_BOX_SHIFT(fineCellBox, iv),
              //                                 CHF_CONST_INTVECT(m_nRefine),
              //                                 CHF_CONST_REAL(scaleNegated), // -scale
              //                                 CHF_CONST_INT(srcVecComp),
              //                                 CHF_CONST_INT(destVecComp),
              //                                 CHF_CONST_INT(ncomp));
             
              Box coarsenedFineCellBox(fineCellBox);
              coarsenedFineCellBox.shift(-iv);
              coarsenedFineCellBox.coarsen(m_nRefine);
              // Retrieve cFine((fineCellBox+iv)/m_nRefine - civ)
              // on components destVecComp+(0:ncomp-1).
              cFine.shift(-civ);
              fineFluxSaveFab.shift(-civ);
              cFine.copy(fineFluxSaveFab,
                         coarsenedFineCellBox, destVecComp,
                         coarsenedFineCellBox, destVecComp, ncomp);
              fineFluxSaveFab.shift(civ);
              cFine.shift(civ);

              /*
                Now vector components of
                cFine = m_fineFlux[a_fineDataIndex]
                on block boundary are back to what they were before calling
                LevelFluxRegister::incrementFine().
              */

              // Now that we have undone that, let's do what we need to do.

              const NewFourthOrderCoordSys* thisCoordSysFine =
                static_cast<const NewFourthOrderCoordSys*>
                (m_coordSysFine->getCoordSys(blockID));
              // FIXME: Make this work with other coordinate systems.
              const CubedSphereShellPanelCS* thisPanelCS =
                dynamic_cast<const CubedSphereShellPanelCS*>(thisCoordSysFine);
             
              int otherBlockID = bb.neighbor();
              const NewFourthOrderCoordSys* otherCoordSysFine =
                static_cast<const NewFourthOrderCoordSys*>
                (m_coordSysFine->getCoordSys(otherBlockID));
              // FIXME: Make this work with other coordinate systems.
              const CubedSphereShellPanelCS* otherPanelCS =
                dynamic_cast<const CubedSphereShellPanelCS*>(otherCoordSysFine);
             
              Interval baseVectorIntv(0, SpaceDim-1);

              Box fineFaceBox(fineCellBox);
              fineFaceBox.shiftHalf(a_dir, shiftIn);
             
              FArrayBox fineVectorFaceFlux(fineFaceBox, SpaceDim);
              fineVectorFaceFlux.copy(a_fineFlux,
                                      a_vectorFluxIntv.begin(), 0, SpaceDim);

              // Convert fineVectorFaceFlux to other block's frame:

              // 1. Convert face-centered fineVectorFaceFlux
              // from this block's coordinate frame
              // to this block's orthonormal frame on this boundary.
              thisPanelCS->orthonormalizeVectorFluxes(fineVectorFaceFlux,
                                                      fineFaceBox,
                                                      baseVectorIntv);

              // 2. Convert face-centered fineVectorFaceFlux
              // in this block's orthonormal frame on this boundary
              // to face-centered fineVectorFaceOtherFlux
              // in other block's orthonormal frame on this boundary.
              // To do this, we transform indices and possibly reorient.
              IndicesTransformation it = bb.getTransformation();
              Box faceOtherBox = it.transform(fineFaceBox);
              FArrayBox fineVectorFaceOtherFlux(faceOtherBox, SpaceDim);
              it.transformFwd(fineVectorFaceOtherFlux, fineVectorFaceFlux,
                              fineFaceBox, baseVectorIntv, baseVectorIntv);

              // Reorient orthonormal-frame flux components where necessary.
              int reorientNormal = bb.reorientFace(a_dir);
              int orthoDir = 0; // component in orthonormal frame
              if (reorientNormal == -1)
                {
                  fineVectorFaceOtherFlux.negate(orthoDir, 1);
                }
              for (int tanDir = 0; tanDir < SpaceDim; tanDir++)
                {
                  if (tanDir != a_dir)
                    {
                      orthoDir++;
                      int sgnFace = bb.reorientFace(tanDir);
                      if (sgnFace == -1)
                        {
                          fineVectorFaceOtherFlux.negate(orthoDir, 1);
                        }
                    }
                }

              // 3. Convert face-centered fineVectorFaceOtherFlux
              // from other block's orthonormal frame on this boundary
              // to other block's coordinate frame.
              otherPanelCS->deorthonormalizeVectorFluxes(fineVectorFaceOtherFlux,
                                                         faceOtherBox,
                                                         baseVectorIntv);

              // 4. Transform indices
              // from fineVectorFaceOtherFlux face-centered in other block
              // to fineVectorFaceFlux face-centered in this block.
              it.transformBack(fineVectorFaceFlux, fineVectorFaceOtherFlux,
                               faceOtherBox, baseVectorIntv, baseVectorIntv);
             
              // 5. Increment cFine with coarsened fineVectorFaceFlux.
              int baseComp = 0;
              FORT_INCREMENTFINE(CHF_CONST_FRA_SHIFT(fineVectorFaceFlux, ivLo),
                                 CHF_FRA_SHIFT(cFine, civ),
                                 CHF_BOX_SHIFT(fineCellBox, iv),
                                 CHF_CONST_INTVECT(m_nRefine),
                                 CHF_CONST_REAL(scale),
                                 CHF_CONST_INT(baseComp), // comp in fineVectorFaceFlux
                                 CHF_CONST_INT(destVecComp), // comp in cFine
                                 CHF_CONST_INT(ncomp));
              // dummy statement allows us to set a breakpoint here in gdb
              CH_assert(true);
            }
        }
    }
}

//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
MultiBlockFluxRegister::
incrementFineAveraged(const FArrayBox& a_fineFlux,
                      const FArrayBox& a_fineFlux2,
                      Real a_scale,
                      const DataIndex& a_fineDataIndex,
                      const Interval& a_vectorFluxIntv,
                      const Interval& a_srcInterval,
                      const Interval& a_dstInterval,
                      int a_dir)
{
  incrementFineAveraged(a_fineFlux, a_fineFlux2, a_scale, a_fineDataIndex,
                        a_vectorFluxIntv, a_srcInterval, a_dstInterval,
                        a_dir, Side::Lo);
  incrementFineAveraged(a_fineFlux, a_fineFlux2, a_scale, a_fineDataIndex,
                        a_vectorFluxIntv, a_srcInterval, a_dstInterval,
                        a_dir, Side::Hi);
}

//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
MultiBlockFluxRegister::
incrementFineAveraged(const FArrayBox& a_fineFlux,
                      const FArrayBox& a_fineFlux2,
                      Real a_scale,
                      const DataIndex& a_fineDataIndex,
                      const Interval& a_vectorFluxIntv,
                      const Interval& a_srcInterval,
                      const Interval& a_dstInterval,
                      int a_dir,
                      Side::LoHiSide a_sd)
{
  CH_TIME("MultiBlockFluxRegister::incrementFineAveraged");
  FArrayBox fineFluxSaveFab;
  if (a_vectorFluxIntv.size() > 0)
    { // Copy m_fineFlux[a_fineDataIndex] because it is modified by
      // LevelFluxRegister::incrementFine(), and we want
      // to UNDO that modification in certain places.
      const FArrayBox& fineFluxFab = m_fineFlux[a_fineDataIndex];
      fineFluxSaveFab.define(fineFluxFab.box(), fineFluxFab.nComp());
      fineFluxSaveFab.copy(fineFluxFab);
    }
  LevelFluxRegister::incrementFine(a_fineFlux,
                                   a_scale, a_fineDataIndex,
                                   a_srcInterval, a_dstInterval, a_dir, a_sd);

  if (a_vectorFluxIntv.size() > 0)
    {
      int blockID = m_blockNum[a_fineDataIndex];
      int faceID = a_dir;
      if (a_sd == Side::Hi) faceID += SpaceDim;
      const BlockBoundary& bb = m_coordSysFine->boundaries()[blockID][faceID];
      if ( !bb.isDomainBoundary() )
        {
          // Box blockCoarseBox = m_coordSys->mappingBlocks()[blockID];
          // Box blockFineBox = refine(blockCoarseBox, m_nRefine);
          Box blockFineBox = m_coordSysFine->mappingBlocks()[blockID];
          Box adjFineCellBlockBox = adjCellBox(blockFineBox, a_dir, a_sd, 1);

          int shiftOut = sign(a_sd);
          int shiftIn = -sign(a_sd);

          const Box& fineFluxFaceBox = a_fineFlux.box();
          Box fineFluxCellBox(fineFluxFaceBox);
          fineFluxCellBox.shiftHalf(a_dir, shiftOut);
          Box adjFineFluxBox(fineFluxCellBox);
          adjFineFluxBox &= adjFineCellBlockBox;
          if ( !adjFineFluxBox.isEmpty() )
            { // Yes: On block-boundary face.  Need to readjust m_fineFlux.

              Box clipBox = m_fineFlux.box(a_fineDataIndex); // NO ghosts
              clipBox.refine(m_nRefine);
              Box fineCellBox = adjCellBox(clipBox, a_dir, a_sd, 1);
              fineCellBox &= fineFluxCellBox;
              // fineCellBox is just outside block, and has no ghosts.

              // cFine lives on the COARSENED cell-centered fine patch.
              FArrayBox& cFine = m_fineFlux[a_fineDataIndex];

              /*
                First replace vector components of
                cFine = m_fineFlux[a_fineDataIndex]
                on block boundary
                with original before calling LevelFluxRegister::incrementFine().
              */
              Real denom = 1.0;
              if (m_scaleFineFluxes)
                {
                  denom *= m_nRefine.product() / m_nRefine[a_dir];
                }
              Real scale = sign(a_sd)*a_scale/denom;
             
              const IntVect& iv = fineCellBox.smallEnd();
              IntVect civ = coarsen(iv, m_nRefine);

              // a_vectorFluxIntv should be within a_srcInterval.
              CH_assert(a_srcInterval.contains(a_vectorFluxIntv.begin()));
              CH_assert(a_srcInterval.contains(a_vectorFluxIntv.end()));
              // first of all comps in source, dest
              int srcComp = a_srcInterval.begin();
              int destComp = a_dstInterval.begin();
              // first of the vector comps in source, dest
              int srcVecComp = a_vectorFluxIntv.begin();
              int destVecComp = (srcVecComp - srcComp) + destComp;
              int ncomp = SpaceDim; // vector components only
              IntVect ivLo = iv;
              // Shift ivLo from cell-centered just outside block
              // to face-centered on block boundary (like a_fineFlux).
              if (a_sd == Side::Lo) ivLo += BASISV(a_dir);
              //             FORT_INCREMENTFINE(CHF_CONST_FRA_SHIFT(a_fineFlux, ivLo),
              //                                CHF_FRA_SHIFT(cFine, civ),
              //                                CHF_BOX_SHIFT(fineCellBox, iv),
              //                                CHF_CONST_INTVECT(m_nRefine),
              //                                CHF_CONST_REAL(scaleNegated), // -scale
              //                                CHF_CONST_INT(srcVecComp),
              //                                CHF_CONST_INT(destVecComp),
              //                                CHF_CONST_INT(ncomp));
              Box coarsenedFineCellBox(fineCellBox);
              coarsenedFineCellBox.shift(-iv);
              coarsenedFineCellBox.coarsen(m_nRefine);
              // Retrieve cFine((fineCellBox+iv)/m_nRefine - civ)
              // on components destVecComp+(0:ncomp-1).
              cFine.shift(-civ);
              fineFluxSaveFab.shift(-civ);
              cFine.copy(fineFluxSaveFab,
                         coarsenedFineCellBox, destVecComp,
                         coarsenedFineCellBox, destVecComp, ncomp);
              fineFluxSaveFab.shift(civ);
              cFine.shift(civ);
             
              /*
                Now vector components of
                cFine = m_fineFlux[a_fineDataIndex]
                on block boundary are back to what they were before calling
                LevelFluxRegister::incrementFine().
              */

              // Now that we have undone that, let's do what we need to do.
             
              const NewFourthOrderCoordSys* thisCoordSysFine =
                static_cast<const NewFourthOrderCoordSys*>
                (m_coordSysFine->getCoordSys(blockID));
              // FIXME: Make this work with other coordinate systems.
              const CubedSphereShellPanelCS* thisPanelCS =
                dynamic_cast<const CubedSphereShellPanelCS*>(thisCoordSysFine);
             
              int otherBlockID = bb.neighbor();
              const NewFourthOrderCoordSys* otherCoordSysFine =
                static_cast<const NewFourthOrderCoordSys*>
                (m_coordSysFine->getCoordSys(otherBlockID));
              // FIXME: Make this work with other coordinate systems.
              const CubedSphereShellPanelCS* otherPanelCS =
                dynamic_cast<const CubedSphereShellPanelCS*>(otherCoordSysFine);

              Interval baseVectorIntv(0, SpaceDim-1);
             
              // fineCellBox1 is fineCellBox grown in tangential directions
              // and cut off at block boundaries in tangential directions.
              Box fineCellBox1 = grow(fineCellBox,
                                      IntVect::Unit - BASISV(a_dir));

              Box fineFaceDomainBox1(fineCellBox1);
              fineFaceDomainBox1.shiftHalf(a_dir, shiftIn);
              ProblemDomain probDomain1(fineFaceDomainBox1);

              MOLUtilities molUtilThis;
              Real bogusDx = 1.; // not used in deconvolveCenterFace()
              molUtilThis.define(probDomain1, bogusDx);

              Box fineFaceBox(fineCellBox);
              fineFaceBox.shiftHalf(a_dir, shiftIn);
              Box fineFaceBox1(fineCellBox1);
              fineFaceBox1.shiftHalf(a_dir, shiftIn);

              FArrayBox fineVectorFaceFlux(fineFaceBox, SpaceDim);
              FArrayBox fineVectorFaceFlux2(fineFaceBox1, SpaceDim);
              fineVectorFaceFlux.copy(a_fineFlux,
                                      a_vectorFluxIntv.begin(), 0, SpaceDim);
              fineVectorFaceFlux2.copy(a_fineFlux2,
                                       a_vectorFluxIntv.begin(), 0, SpaceDim);

              // Convert fineVectorFaceFlux to face-centered.
              // (This is a copy, not an alias; we do not change a_fineFlux.)
              // int dirFace = faceDimension(faceBox);
              molUtilThis.deconvolveCenterFace(fineVectorFaceFlux,
                                               fineVectorFaceFlux2,
                                               fineFaceBox, a_dir, -1);

              // Convert fineVectorFaceFlux to other block's frame:

              // 1. Convert face-centered fineVectorFaceFlux
              // from this block's coordinate frame
              // to this block's orthonormal frame on this boundary.
              thisPanelCS->orthonormalizeVectorFluxes(fineVectorFaceFlux,
                                                      fineFaceBox,
                                                      baseVectorIntv);
              thisPanelCS->orthonormalizeVectorFluxes(fineVectorFaceFlux2,
                                                      fineFaceBox1,
                                                      baseVectorIntv);

              // 2. Convert face-centered fineVectorFaceFlux
              // in this block's orthonormal frame on this boundary
              // to face-centered fineVectorFaceOtherFlux
              // in other block's orthonormal frame on this boundary.
              // To do this, we transform indices and possibly reorient.
              IndicesTransformation it = bb.getTransformation();
              Box faceOtherBox = it.transform(fineFaceBox);
              Box faceOtherBox1 = it.transform(fineFaceBox1);
              FArrayBox fineVectorFaceOtherFlux(faceOtherBox, SpaceDim);
              FArrayBox fineVectorFaceOtherFlux2(faceOtherBox1, SpaceDim);
              it.transformFwd(fineVectorFaceOtherFlux, fineVectorFaceFlux,
                              fineFaceBox, baseVectorIntv, baseVectorIntv);
              it.transformFwd(fineVectorFaceOtherFlux2, fineVectorFaceFlux2,
                              fineFaceBox1, baseVectorIntv, baseVectorIntv);

              // Reorient orthonormal-frame flux components where necessary.
              int reorientNormal = bb.reorientFace(a_dir);
              int orthoDir = 0; // component in orthonormal frame
              if (reorientNormal == -1)
                {
                  fineVectorFaceOtherFlux.negate(orthoDir, 1);
                  fineVectorFaceOtherFlux2.negate(orthoDir, 1);
                }
              for (int tanDir = 0; tanDir < SpaceDim; tanDir++)
                {
                  if (tanDir != a_dir)
                    {
                      orthoDir++;
                      int sgnFace = bb.reorientFace(tanDir);
                      if (sgnFace == -1)
                        {
                          fineVectorFaceOtherFlux.negate(orthoDir, 1);
                          fineVectorFaceOtherFlux2.negate(orthoDir, 1);
                        }
                    }
                }

              // 3. Convert face-centered fineVectorFaceOtherFlux
              // from other block's orthonormal frame on this boundary
              // to other block's coordinate frame.
              otherPanelCS->deorthonormalizeVectorFluxes(fineVectorFaceOtherFlux,
                                                         faceOtherBox,
                                                         baseVectorIntv);
              otherPanelCS->deorthonormalizeVectorFluxes(fineVectorFaceOtherFlux2,
                                                         faceOtherBox1,
                                                         baseVectorIntv);

              int dirOther = bb.dirOther(a_dir);
              // int shiftOtherIn = (reorientNormal == -1) ? shiftIn : shiftOut;
              int shiftOtherOut = (reorientNormal == -1) ? shiftOut : shiftIn;
              // Like fineCellBox, fineOtherBox is CELL-centered,
              // just outside block, and has no ghosts.
              Box fineOtherBox = it.transform(fineCellBox);
              Box fineOtherBox1 = grow(fineOtherBox,
                                       IntVect::Unit - BASISV(dirOther));

              Box fineFaceDomainOtherBox1(fineOtherBox1);
              fineFaceDomainOtherBox1.shiftHalf(dirOther, shiftOtherOut);
              ProblemDomain probDomainOther1(fineFaceDomainOtherBox1);

              MOLUtilities molUtilOther;
              molUtilOther.define(probDomainOther1, bogusDx);
              molUtilOther.deconvolveCenterFace(fineVectorFaceOtherFlux,
                                                fineVectorFaceOtherFlux2,
                                                faceOtherBox, dirOther, 1);

              // 4. Transform indices
              // from fineVectorFaceOtherFlux face-centered in other block
              // to fineVectorFaceFlux face-centered in this block.
              it.transformBack(fineVectorFaceFlux, fineVectorFaceOtherFlux,
                               faceOtherBox, baseVectorIntv, baseVectorIntv);
             
              // 5. Increment cFine with coarsened fineVectorFaceFlux.
              int baseComp = 0;
              FORT_INCREMENTFINE(CHF_CONST_FRA_SHIFT(fineVectorFaceFlux, ivLo),
                                 CHF_FRA_SHIFT(cFine, civ),
                                 CHF_BOX_SHIFT(fineCellBox, iv),
                                 CHF_CONST_INTVECT(m_nRefine),
                                 CHF_CONST_REAL(scale),
                                 CHF_CONST_INT(baseComp), // comp in fineVectorFaceFlux
                                 CHF_CONST_INT(destVecComp), // comp in cFine
                                 CHF_CONST_INT(ncomp));
              // dummy statement allows us to set a breakpoint here in gdb
              CH_assert(true);
            }
        }
    }
}

//-----------------------------------------------------------------------
#include "NamespaceFooter.H"
