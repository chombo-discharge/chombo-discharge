#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MultiBlockLevelCopier.H"
#include "FourthOrderUtil.H"
#include "MBMiscUtil.H"
#include "MBStencilIterator.H"
#include "MBVectorStencilIterator.H"
#include "CH_Timer.H"

#include "NamespaceHeader.H"

/// constructor
MultiBlockLevelCopier::MultiBlockLevelCopier(const MultiBlockLevelGeom*  a_geomPtr,
                                             const BoxLayout&            a_dstLayout,
                                             int                         a_order)
{
  define(a_geomPtr, a_dstLayout, a_order);
}

MultiBlockLevelCopier::MultiBlockLevelCopier(const MultiBlockLevelGeom*  a_geomPtr,
                                             const DisjointBoxLayout&    a_dstDisjointLayout,
                                             const IntVect&              a_ghostVect,
                                             int                         a_order)
{
  define(a_geomPtr, a_dstDisjointLayout, a_ghostVect, a_order);
}

/// destructor
MultiBlockLevelCopier::~MultiBlockLevelCopier()
{
  undefine();
}

/// destructor
void MultiBlockLevelCopier::undefine()
{
  if (m_isDefined && m_nblocks > 1)
    {
      m_powers.clear();
      DataIterator ditDst = m_dstLayout.dataIterator();
      for (ditDst.begin(); ditDst.ok(); ++ditDst)
        {
          delete m_validIndices[ditDst];
          delete m_validBlock[ditDst];
          delete m_validMappedCenter[ditDst];
        }
      for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
        {
          delete m_validLayout[srcBlock];
          delete m_validFullLayout[srcBlock];
          delete m_copiers[srcBlock];
          delete m_copiersFull[srcBlock];
        }
    }
  m_isVectorDefined = false;
  m_isDefined = false;
}


void
MultiBlockLevelCopier::define(const MultiBlockLevelGeom*  a_geomPtr,
                              const DisjointBoxLayout&    a_dstDisjointLayout,
                              const IntVect&              a_ghostVect,
                              int                         a_order)
{
  undefine();
  BoxLayout dstLayout;
  dstLayout.deepCopy(a_dstDisjointLayout);
  dstLayout.close();
  dstLayout.grow(a_ghostVect);
  define(a_geomPtr, dstLayout, a_order);
}


void
MultiBlockLevelCopier::define(const MultiBlockLevelGeom*  a_geomPtr,
                              const BoxLayout&            a_dstLayout,
                              int                         a_order)
{
  CH_TIME("MultiBlockLevelCopier::define");
  m_geomPtr = (MultiBlockLevelGeom*) a_geomPtr;
  CH_assert(m_geomPtr != NULL);
  CH_assert(m_geomPtr->isDefined());
  m_order = a_order;
  CH_assert((m_order == 1) || (m_order == 2) || (m_order == 4)); // want order 5 eventually

  /*
    Set stuff from m_geomPtr:
    DisjointBoxLayout m_srcLayout;
    MultiBlockCoordSys* m_coordSysPtr;
   */
  m_srcLayout = m_geomPtr->grids();
  m_srcFullLayout = m_geomPtr->gridsFull();

  m_fixedDims = m_geomPtr->fixedDims();
  m_fixedDimsVect = m_geomPtr->fixedDimsVect();
  m_interpDimsVect = m_geomPtr->interpDimsVect();
  m_fixedPt = m_geomPtr->fixedPt();
  m_gridsFixedOff = m_geomPtr->gridsFixedOff();

  m_coordSysPtr = a_geomPtr->coordSysPtr();
  CH_assert(m_coordSysPtr != NULL);
  CH_assert(m_coordSysPtr->isDefined());
  m_nblocks = m_coordSysPtr->numBlocks();
  if (m_nblocks == 1)
    { // only one block; do nothing
      return;
    }

  m_mbUtil = m_geomPtr->mbUtil();

  // FIXME: This could be wrong if the fixed slice is in overlap between boxes.
  m_allGridsHaveFixedPt = m_mbUtil->allGridsHaveFixedPt(a_dstLayout);
  m_mbUtil->getCollapsedLayout(m_dstLayout, a_dstLayout);
  m_mbUtil->getFixedOffLayout(m_dstFixedOff, a_dstLayout);

  /*
    Set basic stuff from m_order:
    int m_radius;
    int m_degree;
    Box m_degreeBox;
    int m_numTaylorCoeffs;
  */
  {
    if (m_order <= 1)
      {
        m_radius = 0;
        m_degree = 0;
      }
    else
      {
        m_radius = 1; // (m_order/2);
        // This radius gives us 2*floor(m_order/2)+1 cells in a row,
        // hence at least m_order+1 cells in a row,
        // which is exactly what we need to find (m_order+1)'th derivatives.
        // NEW, 22 Feb 2011: we have order either 4 or 5 only.
        m_degree = m_order - 1; // either 3 or 4
      }
    // Cleared this in undefine();
    // Vector<IntVect> m_powers;
    m_degreeBox = Box(IntVect::Zero, m_degree * IntVect::Unit);
    // Collapse m_degreeBox in the dimensions in m_fixedDimsVect.
    for (int ind = 0; ind < m_fixedDimsVect.size(); ind++)
      {
        int idir = m_fixedDimsVect[ind];
        m_degreeBox.setRange(idir, 0);
      }
    for (BoxIterator bitPower(m_degreeBox); bitPower.ok(); ++bitPower)
      {
        IntVect pwr = bitPower();
        if (pwr.sum() <= m_degree)
          {
            m_powers.push_back(pwr);
          }
      }
    // m_numTaylorCoeffs = C(m_degree + SpaceDim, SpaceDim)
    m_numTaylorCoeffs = m_powers.size();
  }

  /*
    Everything we have set so far has been independent of the layouts.

    Now define the layout-dependent members of this class.
   */

  /*
    Define LayoutData< IntVectSet > m_extraCellsIVS.
  */
  m_extraCellsIVS.define(m_dstLayout);
  m_blockNumbersDst.define(m_dstLayout);
  DataIterator ditDst = m_dstLayout.dataIterator();
  const Vector<Box>& mappingBlocksAll = m_coordSysPtr->mappingBlocks();
  for (ditDst.begin(); ditDst.ok(); ++ditDst)
    {
      const Box& dstBox = m_dstLayout[ditDst];
      // Take care of the case of dstBox going outside bounds in 
      // one of the fixed dimensions, where it is a slice.
      // This will happen if dstBox contains ghost cells in any
      // fixed dimension.
      Box dstFixedBox = m_mbUtil->boxFixed(dstBox);
      int blockNumberDst = m_coordSysPtr->whichBlockOverlap(dstFixedBox);
      m_blockNumbersDst[ditDst] = blockNumberDst;
      const Box& blockBox = mappingBlocksAll[blockNumberDst];
      if (blockBox.contains(dstBox))
        { // In this case, dstBox has NO extra-block cells.
          m_extraCellsIVS[ditDst] = IntVectSet(); // empty
        }
      else
        { // In this case, dstBox has cells outside blockBox.
          DenseIntVectSet ivs = DenseIntVectSet(dstBox, true);
          ivs -= blockBox;
          // I think we should do this because ivs is DenseIntVectSet.
          ivs.recalcMinBox();
          // convert from DenseIntVectSet
          IntVectSet ivsGen = IntVectSet(ivs);
          m_extraCellsIVS[ditDst] = ivsGen;
        }
    }

  // On thisExtraCell in extraCellsIVSPatch, we'll have
  // MBStencil *m_stencils[ditDst](thisExtraCell, 0)
  m_stencils.define(m_dstLayout);

  for (ditDst.begin(); ditDst.ok(); ++ditDst)
    {
      const IntVectSet& extraCellsIVSPatch = m_extraCellsIVS[ditDst];
        m_stencils[ditDst] =
          RefCountedPtr< IVSFAB<MBStencil> >(new IVSFAB<MBStencil>(extraCellsIVSPatch, 1));
    }

  m_validIndices.define(m_dstLayout);
  m_validBlock.define(m_dstLayout);
  m_validMappedCenter.define(m_dstLayout);
  for (ditDst.begin(); ditDst.ok(); ++ditDst)
    {
      const IntVectSet& extraCellsIVSPatch = m_extraCellsIVS[ditDst];
      m_validIndices[ditDst] = new IVSFAB<IntVect>(extraCellsIVSPatch, 1);
      m_validBlock[ditDst] = new IVSFAB<int>(extraCellsIVSPatch, 1);
      m_validMappedCenter[ditDst] = new IVSFAB<RealVect>(extraCellsIVSPatch, 1);
      m_mbUtil->getValid(*m_validIndices[ditDst],
                         *m_validBlock[ditDst],
                         *m_validMappedCenter[ditDst],
                         extraCellsIVSPatch,
                         m_blockNumbersDst[ditDst]);
    }

  // Vector< RefCountedPtr< LayoutData<Box> > > m_stencilCellsMinBox;
  // Vector< RefCountedPtr< LayoutData<Box> > > m_stencilCellsFullMinBox;
  m_stencilCellsMinBox.resize(m_nblocks);
  m_stencilCellsFullMinBox.resize(m_nblocks);
  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    {
      m_stencilCellsMinBox[srcBlock] = RefCountedPtr< LayoutData<Box> >
        (new LayoutData<Box>(m_dstLayout));
      m_stencilCellsFullMinBox[srcBlock] = RefCountedPtr< LayoutData<Box> >
        (new LayoutData<Box>(m_dstLayout));
    }

  // Find the stencil cells and their source blocks.
  // We will NOT find coefficients now, because those depend on centering;
  // they will be computed in the derived class that specifies centering.
  for (ditDst.begin(); ditDst.ok(); ++ditDst)
    {
      Vector<Box> stencilsMinBox;
      m_mbUtil->getStencilCells(*m_stencils[ditDst],
                                stencilsMinBox,
                                m_extraCellsIVS[ditDst],
                                *m_validIndices[ditDst],
                                *m_validBlock[ditDst],
                                *m_validMappedCenter[ditDst],
                                m_order,
                                m_radius);

      const Box& bxFixedOff = m_dstFixedOff[ditDst];
      for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
        {
          (*m_stencilCellsMinBox[srcBlock])[ditDst] = Box(stencilsMinBox[srcBlock]);

          // If no fixed dimensions, then this will not change.
          Box stencilsFullMinBoxSrc(stencilsMinBox[srcBlock]);
          if ( !stencilsFullMinBoxSrc.isEmpty() )
            { // If stencilsMinBoxSrc empty, keep it empty.
              for (int ind = 0; ind < m_fixedDimsVect.size(); ind++)
                {
                  int idir = m_fixedDimsVect[ind]; 
                  // stencilsFullMinBoxSrc[idir] now has range only 1 cell,
                  // being the minimum, so expand it by resetting big end.
                  int oldHi = stencilsFullMinBoxSrc.bigEnd(idir);
                  int newHi = oldHi + bxFixedOff.bigEnd(idir);
                  stencilsFullMinBoxSrc.setBig(idir, newHi);
                }
            }

          (*m_stencilCellsFullMinBox[srcBlock])[ditDst] = Box(stencilsFullMinBoxSrc);
        }
    } // end loop over dest patches

  m_validLayout.resize(m_nblocks);
  m_validFullLayout.resize(m_nblocks);
  m_copiers.resize(m_nblocks);
  m_copiersFull.resize(m_nblocks);
  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    {
      m_validLayout[srcBlock] =
         new BoxLayout(*m_stencilCellsMinBox[srcBlock]);
      m_validFullLayout[srcBlock] =
         new BoxLayout(*m_stencilCellsFullMinBox[srcBlock]);

      m_copiers[srcBlock] =
         new Copier(m_srcLayout, *m_validLayout[srcBlock]);
      m_copiersFull[srcBlock] =
         new Copier(m_srcFullLayout, *m_validFullLayout[srcBlock]);
    }

  if ( m_allGridsHaveFixedPt )
    {
      for (ditDst.begin(); ditDst.ok(); ++ditDst)
        {
          m_mbUtil->getWeights(*m_stencils[ditDst],
                               m_type,
                               m_extraCellsIVS[ditDst],
                               *m_validIndices[ditDst],
                               *m_validBlock[ditDst],
                               *m_validMappedCenter[ditDst],
                               m_powers,
                               m_blockNumbersDst[ditDst]);
          // dummy statement in order to get around gdb bug
          int dummy_unused = 0; dummy_unused++;
        } // end loop over patches
    }
  else // !m_allGridsHaveFixedPt
    {
      // Set Box blockAllDstBox[dstBlock] to the min Box covering
      // every Box in m_dstLayout in m_interpDimsVect
      // and set to m_fixedPt in m_fixedDimsVect.
      // Note that this is only for boxes on the same processor;
      // we'll need to call minBoxGlobal() on each.
      Vector<Box> blockAllDstBox(m_nblocks);
      for (ditDst.begin(); ditDst.ok(); ++ditDst)
        {
          Box dstFixedBox = m_mbUtil->boxFixed(m_dstLayout[ditDst]);
          int blockNum = m_blockNumbersDst[ditDst];
          blockAllDstBox[blockNum].minBox(dstFixedBox);
        }
      Vector< RefCountedPtr< IVSFAB<MBStencil> > > blockStencils(m_nblocks);
      for (int dstBlock = 0; dstBlock < m_nblocks; dstBlock++)
        {
          Box& blockAllDstBoxDst = blockAllDstBox[dstBlock];
          minBoxGlobal(blockAllDstBoxDst);
          DenseIntVectSet blockAllDstDenseIVS(blockAllDstBoxDst, true);
          // Remove valid cells of block dstBlock,
          // retaining only the extra-block ghost cells.
          blockAllDstDenseIVS -= mappingBlocksAll[dstBlock];
          IntVectSet blockAllDstIVS(blockAllDstDenseIVS);

          IVSFAB<IntVect> blockValidIndices(blockAllDstIVS, 1);
          IVSFAB<int> blockValidBlock(blockAllDstIVS, 1);
          IVSFAB<RealVect> blockValidMappedCenter(blockAllDstIVS, 1);
          m_mbUtil->getValid(blockValidIndices,
                             blockValidBlock,
                             blockValidMappedCenter,
                             blockAllDstIVS,
                             dstBlock);
          
          blockStencils[dstBlock] = RefCountedPtr< IVSFAB<MBStencil> >
            (new IVSFAB<MBStencil>(blockAllDstIVS, 1));
          IVSFAB<MBStencil>& blockStencilsDst = *blockStencils[dstBlock];
          Vector<Box> stencilsMinBox;
          m_mbUtil->getStencilCells(blockStencilsDst,
                                    stencilsMinBox,
                                    blockAllDstIVS,
                                    blockValidIndices,
                                    blockValidBlock,
                                    blockValidMappedCenter,
                                    m_order,
                                    m_radius);

          m_mbUtil->getWeights(blockStencilsDst,
                               m_type,
                               blockAllDstIVS,
                               blockValidIndices,
                               blockValidBlock,
                               blockValidMappedCenter,
                               m_powers,
                               dstBlock);
        }
      for (ditDst.begin(); ditDst.ok(); ++ditDst)
        {
          int dstBlock = m_blockNumbersDst[ditDst];
          const IVSFAB<MBStencil>& blockStencilsDst = *blockStencils[dstBlock];
          m_mbUtil->copyStencilsFromBlock(*m_stencils[ditDst],
                                          m_extraCellsIVS[ditDst],
                                          blockStencilsDst);
        }
    }
}


void
MultiBlockLevelCopier::defineVector()
{ CH_TIME("MultiBlockLevelCopier::defineVector");
  CH_assert(isDefined());
  if (m_coordSysPtr->numBlocks() == 1)
    { // only one block; do nothing
      return;
    }
  DataIterator ditDst = m_dstLayout.dataIterator();
  m_vectorstencils.define(m_dstLayout);
  IntVect growStencilVect = IntVect::Zero;
  if (m_type == IndexType::TheCellType())
    { // need an additional cell in each interpolated dimension.
      int nInterpDims = m_interpDimsVect.size();
      for (int idir = 0; idir < nInterpDims; idir++)
        {
          int interpDir = m_interpDimsVect[idir];
          growStencilVect[interpDir] = 1;
        }
    }
  for (ditDst.begin(); ditDst.ok(); ++ditDst)
    {
      const IntVectSet& extraCellsIVSPatch = m_extraCellsIVS[ditDst];
      int ghostBlockNum = m_blockNumbersDst[ditDst];
      m_vectorstencils[ditDst] = RefCountedPtr< IVSFAB<MBVectorStencil> >
        (new IVSFAB<MBVectorStencil>(extraCellsIVSPatch, 1));
      m_mbUtil->getVectorStencilCells(*m_vectorstencils[ditDst],
                                      *m_stencils[ditDst],
                                      extraCellsIVSPatch);
      Vector< BaseFab<VectorTransformation>* > vectorStencilTransformations(m_nblocks);
      for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
        {
          Box bxStencilSrc = (*m_stencilCellsMinBox[srcBlock])[ditDst];
          bxStencilSrc.grow(growStencilVect);
          vectorStencilTransformations[srcBlock] =
            new BaseFab<VectorTransformation>(bxStencilSrc, 1);
        }
      m_mbUtil->getVectorTransformations(vectorStencilTransformations,
                                         *m_vectorstencils[ditDst],
                                         m_type,
                                         extraCellsIVSPatch,
                                         ghostBlockNum);
      m_mbUtil->getVectorWeights(*m_vectorstencils[ditDst],
                                 vectorStencilTransformations,
                                 m_type,
                                 extraCellsIVSPatch,
                                 *m_validIndices[ditDst],
                                 *m_validBlock[ditDst],
                                 *m_validMappedCenter[ditDst],
                                 m_powers,
                                 ghostBlockNum);
      for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
        {
          delete vectorStencilTransformations[srcBlock];
        }
    }
  m_isVectorDefined = true;
}


void
MultiBlockLevelCopier::copyTo(const LevelData<FArrayBox>&  a_src,
                              BoxLayoutData<FArrayBox>&    a_dst,
                              const Interval&        a_intvl) const
{
  CH_TIME("MultiBlockLevelCopier::copyTo");
  if (m_nblocks == 1)
    { // only one block; do nothing
      return;
    }
  int ncomp = a_intvl.size();
  Interval intvl0(0, ncomp-1);

  // BoxLayoutData *validData[srcBlock] lives on *m_validFullLayout[srcBlock]
  // and will hold
  // all valid data from srcBlock needed to fill m_extraCellsIVS[ditDst].
  Vector< RefCountedPtr< BoxLayoutData<FArrayBox> > > validData(m_nblocks);
  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    {
      validData[srcBlock] = RefCountedPtr< BoxLayoutData<FArrayBox> >
        (new BoxLayoutData<FArrayBox>(*m_validFullLayout[srcBlock], ncomp));
    }

  // COMMUNICATE from a_src to validData.
  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    {
      BoxLayoutData<FArrayBox>& srcData = *validData[srcBlock];
      const Copier& srcCopier = *m_copiersFull[srcBlock];
      { CH_TIME("copyTo communication");
        a_src.copyTo(a_intvl, srcData, intvl0, srcCopier);
      }
    }

  DataIterator ditDst = m_dstLayout.dataIterator();

  // From here on, do LOCAL interpolation.
  for (ditDst.begin(); ditDst.ok(); ++ditDst)
    {
      const IntVectSet& extraCellsIVSPatch = m_extraCellsIVS[ditDst];
      const IVSFAB<MBStencil>& stencilsPatch = *m_stencils[ditDst];

      FArrayBox& dstFab = a_dst[ditDst];
      // Use dst0Fab so that component indices are always 0:ncomp-1.
      FArrayBox dst0Fab(a_intvl, dstFab);
      // Note dstFab and dst0Fab are on the full uncollapsed layout.

      const Box& bxFixedOff = m_dstFixedOff[ditDst];
      for (IVSIterator ivsit(extraCellsIVSPatch); ivsit.ok(); ++ivsit)
        {
          IntVect thisExtraCell = ivsit();

          // Fill in dst0Fab(fullExtraCell, :) where
          // fullExtraCell[m_fixedDimsVect] == thisExtraCell[m_fixedDimsVect].
          dst0Fab.shift(-thisExtraCell);
          dst0Fab.setVal(0., bxFixedOff, 0, ncomp);
          MBStencilIterator stencilit(stencilsPatch(thisExtraCell, 0));
          for (stencilit.begin(); stencilit.ok(); ++stencilit)
            {
              const MBStencilElement& stencilElement = stencilit();
              int srcBlock = stencilElement.block();
              const IntVect& cell = stencilElement.cell();
              Real wt = stencilElement.weight();

              // non-const only because we shift it
              FArrayBox& validDataSrcFab = (*(validData[srcBlock]))[ditDst];
              validDataSrcFab.shift(-cell);

              dst0Fab.plus(validDataSrcFab,
                           bxFixedOff, // source box
                           bxFixedOff, // dest box
                           wt, // scaling
                           0, // start source component index
                           0, // start dest component index
                           ncomp); // number of components
              // shift back to where it started
              validDataSrcFab.shift(cell);
            }
          // shift back
          dst0Fab.shift(thisExtraCell);
        }
    }
}


void
MultiBlockLevelCopier::copyTo(const LevelData<FArrayBox>&  a_src,
                              BoxLayoutData<FArrayBox>&    a_dst) const
{
  if (m_coordSysPtr->numBlocks() == 1)
    { // only one block; do nothing
      return;
    }
  const Interval& intvl = a_src.interval();
  copyTo(a_src, a_dst, intvl);
}


void
MultiBlockLevelCopier::copyToVector(const LevelData<FArrayBox>&  a_src,
                                    BoxLayoutData<FArrayBox>&    a_dst) const
{
  const Interval& intvl = a_dst.interval();
  copyToVector(a_src, a_dst, intvl);
}


void
MultiBlockLevelCopier::copyToVector(const LevelData<FArrayBox>&  a_src,
                                    BoxLayoutData<FArrayBox>&    a_dst,
                                    const Interval&        a_intvl) const
{
  CH_assert(isDefined());
  CH_assert(m_isVectorDefined);
  CH_TIME("MultiBlockLevelCopier::copyToVector");
  if (a_dst.boxLayout().size() == 0)
    {
      return;
    }
  if (m_coordSysPtr->numBlocks() == 1)
    { // only one block; do nothing
      return;
    }
  int ncomp = a_intvl.size();
  // This function is only for vectors of length SpaceDim.
  CH_assert(ncomp == SpaceDim);
  Interval intvl0(0, ncomp-1);

  // BoxLayoutData *validData[srcBlock] lives on *m_validLayout[srcBlock]
  // and will hold
  // all valid data from srcBlock needed to fill m_extraCellsIVS[ditDst].
  Vector< RefCountedPtr< BoxLayoutData<FArrayBox> > > validData(m_nblocks);
  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    {
      validData[srcBlock] = RefCountedPtr< BoxLayoutData<FArrayBox> >
        (new BoxLayoutData<FArrayBox>(*m_validFullLayout[srcBlock], ncomp));
    }

  // COMMUNICATE from a_src to validData.
  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    {
      const Copier& srcCopier = *m_copiersFull[srcBlock];
      BoxLayoutData<FArrayBox>& srcData = *validData[srcBlock];
      { CH_TIME("copyToVector communication");
        a_src.copyTo(a_intvl, srcData, intvl0, srcCopier);
      }
    }
  DataIterator ditDst = m_dstLayout.dataIterator();

  // From here on, do LOCAL interpolation.
  for (ditDst.begin(); ditDst.ok(); ++ditDst)
    {
      const IntVectSet& extraCellsIVSPatch = m_extraCellsIVS[ditDst];
      const IVSFAB<MBVectorStencil>& vectorstencilsPatch = *m_vectorstencils[ditDst];

      FArrayBox& dstFab = a_dst[ditDst];
      // Use dst0Fab so that component indices are always 0:ncomp-1.
      FArrayBox dst0Fab(a_intvl, dstFab);
      // Note dstFab and dst0Fab are on the full uncollapsed layout.

      const Box& bxFixedOff = m_dstFixedOff[ditDst];
      for (IVSIterator ivsit(extraCellsIVSPatch); ivsit.ok(); ++ivsit)
        {
          IntVect thisExtraCell = ivsit();

          // Fill in dst0Fab(fullExtraCell, :) where
          // fullExtraCell[m_fixedDimsVect] == thisExtraCell[m_fixedDimsVect].
          dst0Fab.shift(-thisExtraCell);
          dst0Fab.setVal(0., bxFixedOff, 0, ncomp);
          const MBVectorStencil& thisExtraVectorStencil =
            vectorstencilsPatch(thisExtraCell, 0);
          MBVectorStencilIterator vstencilit(thisExtraVectorStencil);
          for (vstencilit.begin(); vstencilit.ok(); ++vstencilit)
            {
              const MBVectorStencilElement& vstencilElement = vstencilit();
              int srcBlock = vstencilElement.block();
              const IntVect& cell = vstencilElement.cell();
              const Tuple<Real, SpaceDim*SpaceDim>& wt = vstencilElement.weight();
              // not const only because we will shift it
              FArrayBox& validDataSrcFab = (*(validData[srcBlock]))[ditDst];

              validDataSrcFab.shift(-cell);

              int tupleIndex = 0;
              for (int src = 0; src < SpaceDim; src++)
                { // contribution from src'th component of function at stencil cell
                  for (int dst = 0; dst < SpaceDim; dst++)
                    { // add to dst'th component of function at ghost cell
                      // In 2D:
                      // wt[0] is weight of src=0 for dst=0
                      // wt[1] is weight of src=0 for dst=1
                      // wt[2] is weight of src=1 for dst=0
                      // wt[3] is weight of src=1 for dst=1
                      dst0Fab.plus(validDataSrcFab,
                                   bxFixedOff, // source box
                                   bxFixedOff, // dest box
                                   wt[tupleIndex], // scaling
                                   src, // source component index
                                   dst, // dest component index
                                   1); // number of components
                      tupleIndex++;
                    }
                }
              // Shift back to where we were before.
              validDataSrcFab.shift(cell);
            }
          // shift back
          dst0Fab.shift(thisExtraCell);
        }
    }
}

void
MultiBlockLevelCopier::copyToAllWithVector(const LevelData<FArrayBox>&  a_src,
                                           BoxLayoutData<FArrayBox>&    a_dst,
                                           const Interval&         a_vecIntvl) const
{
  if (a_vecIntvl.size() == 0)
    {
      copyTo(a_src, a_dst);
    }
  else
    {
      int compVecLo = a_vecIntvl.begin();
      int compVecHi = a_vecIntvl.end();
      const Interval& srcIntvl = a_src.interval();
      CH_assert(a_dst.interval() == srcIntvl);
      int compDataLo = srcIntvl.begin();
      int compDataHi = srcIntvl.end();
      if (compDataLo < compVecLo)
        { // scalar components before vector components
          Interval scalarIntvl(compDataLo, compVecLo-1);
          copyTo(a_src, a_dst, scalarIntvl);
        }
      // vector components
      copyToVector(a_src, a_dst, a_vecIntvl);
      if (compDataHi > compVecHi)
        { // scalar components after vector components
          Interval scalarIntvl(compVecHi+1, compDataHi);
          copyTo(a_src, a_dst, scalarIntvl);
        }
    }
}

#include "NamespaceFooter.H"
