#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MultiBlockLevelExchange.H"
#include "FourthOrderUtil.H"
#include "MBStencilIterator.H"
#include "MBVectorStencilIterator.H"
#include "CH_Timer.H"

#include "NamespaceHeader.H"

/// constructor
MultiBlockLevelExchange::MultiBlockLevelExchange(const MultiBlockLevelGeom*  a_geomPtr,
                                                 int                         a_ghosts,
                                                 int                         a_order)
{
  define(a_geomPtr, a_ghosts, a_order);
}

/// destructor
MultiBlockLevelExchange::~MultiBlockLevelExchange()
{
  undefine();
}

/// destructor
void MultiBlockLevelExchange::undefine()
{
  CH_TIME("MultiBlockLevelExchange::undefine");
  if (m_isDefined)
    {
      m_powers.clear();
      DataIterator dit = m_grids.dataIterator();
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

void MultiBlockLevelExchange::define(const MultiBlockLevelGeom*  a_geomPtr,
                                     int                         a_ghosts,
                                     int                         a_order)
{
  CH_TIME("MultiBlockLevelExchange::define");
  undefine();

  m_geomPtr = (MultiBlockLevelGeom*) a_geomPtr;
  CH_assert(m_geomPtr != NULL);
  CH_assert(m_geomPtr->isDefined());
  m_ghosts = a_ghosts;
  m_order = a_order;
  CH_assert((m_order == 1) || (m_order == 2) || (m_order == 4)); // want order 5 eventually

  /*
    Set stuff from m_geomPtr:
    DisjointBoxLayout m_grids;
    Vector<int> m_fixedDimsVect, m_interpDimsVect, m_fixedPt;
    MultiBlockCoordSys* m_coordSysPtr;
   */
  m_grids = m_geomPtr->grids();
  m_gridsFull = m_geomPtr->gridsFull();

  m_fixedDims = m_geomPtr->fixedDims();
  m_fixedDimsVect = m_geomPtr->fixedDimsVect();
  m_interpDimsVect = m_geomPtr->interpDimsVect();
  m_fixedPt = m_geomPtr->fixedPt();
  m_gridsFixedOff = m_geomPtr->gridsFixedOff();
  m_allGridsHaveFixedPt = m_geomPtr->allGridsHaveFixedPt();

  m_coordSysPtr = a_geomPtr->coordSysPtr();
  CH_assert(m_coordSysPtr != NULL);
  CH_assert(m_coordSysPtr->isDefined());
  m_nblocks = m_coordSysPtr->numBlocks();
  if (m_nblocks == 1)
    { // only one block; do nothing
      return;
    }

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
    Define m_ghostCells, from m_geomPtr and
    m_ghosts, which is width of destination ghost cell layer.
  */
  m_geomPtr->extraBlockGhosts(m_ghostCells, m_ghosts);

  // On thisGhostCell in ghostCellsIVS, we'll have
  // MBStencil *m_stencils[dit](thisGhostCell, 0)
  DataIterator dit = m_grids.dataIterator();
  { CH_TIME("MultiBlockLevelExchange data allocation");
    m_stencils.define(m_grids);
    for (dit.begin(); dit.ok(); ++dit)
      {
        const IntVectSet& ghostCellsIVS = m_ghostCells[dit];
        m_stencils[dit] = RefCountedPtr< IVSFAB<MBStencil> >
          (new IVSFAB<MBStencil>(ghostCellsIVS, 1));
      }
  }

  m_mbUtil = m_geomPtr->mbUtil();

  const LayoutData< IVSFAB<IntVect>* >& validIndices =
    m_geomPtr->validIndices();
  const LayoutData< IVSFAB<int>* >& validBlock =
    m_geomPtr->validBlock();
  const LayoutData< IVSFAB<RealVect>* >& validMappedCenter =
    m_geomPtr->validMappedCenter();

  // Vector< RefCountedPtr< LayoutData<Box> > > m_stencilCellsMinBox;
  // Vector< RefCountedPtr< LayoutData<Box> > > m_stencilCellsFullMinBox;
  m_stencilCellsMinBox.resize(m_nblocks);
  m_stencilCellsFullMinBox.resize(m_nblocks);
  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    {
      m_stencilCellsMinBox[srcBlock] = RefCountedPtr< LayoutData<Box> >
        (new LayoutData<Box>(m_grids));
      m_stencilCellsFullMinBox[srcBlock] = RefCountedPtr< LayoutData<Box> >
        (new LayoutData<Box>(m_grids));
    }

  // Find the stencil cells and their source blocks.
  // We will NOT find coefficients now, because those depend on centering;
  // they will be computed in the derived class that specifies centering.
  for (dit.begin(); dit.ok(); ++dit)
    {
      Vector<Box> stencilsMinBox;
      m_mbUtil->getStencilCells(*m_stencils[dit],
                                stencilsMinBox,
                                m_ghostCells[dit],
                                *validIndices[dit],
                                *validBlock[dit],
                                *validMappedCenter[dit],
                                m_order,
                                m_radius);

      const Box& bxFixedOff = m_gridsFixedOff[dit];
      { CH_TIME("MultiBlockLevelExchange set stencilCellsMinBoxBlock");
        for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
          {
            (*m_stencilCellsMinBox[srcBlock])[dit] = Box(stencilsMinBox[srcBlock]);

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

            (*m_stencilCellsFullMinBox[srcBlock])[dit] = Box(stencilsFullMinBoxSrc);
          }
      }
    } // end loop over patches

  m_validLayout.resize(m_nblocks);
  m_validFullLayout.resize(m_nblocks);
  m_copiers.resize(m_nblocks);
  m_copiersFull.resize(m_nblocks);
  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    { CH_TIME("MultiBlockLevelExchange copier allocation");
      m_validLayout[srcBlock] =
         new BoxLayout(*m_stencilCellsMinBox[srcBlock]);
      m_validFullLayout[srcBlock] =
         new BoxLayout(*m_stencilCellsFullMinBox[srcBlock]);

      m_copiers[srcBlock] =
         new Copier(m_grids, *m_validLayout[srcBlock]);
      m_copiersFull[srcBlock] =
         new Copier(m_gridsFull, *m_validFullLayout[srcBlock]);
    }

  // Remove from m_ghostCells all cells that do not have complete stencils.
  removeNoValidSource();

  const LayoutData<int>& blockNumbers = m_geomPtr->block();
  if ( m_allGridsHaveFixedPt )
    {
      for (dit.begin(); dit.ok(); ++dit)
        {
          m_mbUtil->getWeights(*m_stencils[dit],
                               m_type,
                               m_ghostCells[dit],
                               *validIndices[dit],
                               *validBlock[dit],
                               *validMappedCenter[dit],
                               m_powers,
                               blockNumbers[dit]);
          // dummy statement in order to get around gdb bug
          int dummy_unused = 0; dummy_unused++;
        } // end loop over patches
    }
  else // !m_allGridsHaveFixedPt
    {
      Vector< RefCountedPtr< IVSFAB<MBStencil> > > blockStencils(m_nblocks);
      for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
        {
          // Sets blockBox[m_fixedDimsVect] = m_fixedPt.
          Box blockBox =
            m_mbUtil->boxFixed(m_coordSysPtr->mappingBlocks()[srcBlock]);

          // This function will grow only in m_interpDims, as we wish.
          IntVectSet blockGhostsSrc =
            m_mbUtil->extraBlockGhosts(blockBox, m_ghosts, srcBlock);

          IVSFAB<IntVect> blockValidIndices(blockGhostsSrc, 1);
          IVSFAB<int> blockValidBlock(blockGhostsSrc, 1);
          IVSFAB<RealVect> blockValidMappedCenter(blockGhostsSrc, 1);
          m_mbUtil->getValid(blockValidIndices,
                             blockValidBlock,
                             blockValidMappedCenter,
                             blockGhostsSrc,
                             srcBlock);

          blockStencils[srcBlock] = RefCountedPtr< IVSFAB<MBStencil> >
            (new IVSFAB<MBStencil>(blockGhostsSrc, 1));
          Vector<Box> stencilsMinBox;
          m_mbUtil->getStencilCells(*blockStencils[srcBlock],
                                    stencilsMinBox,
                                    blockGhostsSrc,
                                    blockValidIndices,
                                    blockValidBlock,
                                    blockValidMappedCenter,
                                    m_order,
                                    m_radius);
          m_mbUtil->getWeights(*blockStencils[srcBlock],
                               m_type,
                               blockGhostsSrc,
                               blockValidIndices,
                               blockValidBlock,
                               blockValidMappedCenter,
                               m_powers,
                               srcBlock);
        }
      for (dit.begin(); dit.ok(); ++dit)
        {
          int srcBlock = blockNumbers[dit];
          m_mbUtil->copyStencilsFromBlock(*m_stencils[dit],
                                          m_ghostCells[dit],
                                          *blockStencils[srcBlock]);
        }
    }
}


void
MultiBlockLevelExchange::defineVector()
{ CH_TIME("MultiBlockLevelExchange::defineVector");
  CH_assert(isDefined());
  if (m_coordSysPtr->numBlocks() == 1)
    { // only one block; do nothing
      return;
    }
  const LayoutData<int>& blockNumbers = m_geomPtr->block();
  DataIterator dit = m_grids.dataIterator();
  m_vectorstencils.define(m_grids);
  const LayoutData< IVSFAB<IntVect>* >& validIndices =
    m_geomPtr->validIndices();
  const LayoutData< IVSFAB<int>* >& validBlock =
    m_geomPtr->validBlock();
  const LayoutData< IVSFAB<RealVect>* >& validMappedCenter =
    m_geomPtr->validMappedCenter();

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
  for (dit.begin(); dit.ok(); ++dit)
    {
      const IntVectSet& ghostCellsIVS = m_ghostCells[dit];
      int ghostBlockNum = blockNumbers[dit];
      m_vectorstencils[dit] = RefCountedPtr< IVSFAB<MBVectorStencil> >
        (new IVSFAB<MBVectorStencil>(ghostCellsIVS, 1));
      m_mbUtil->getVectorStencilCells(*m_vectorstencils[dit],
                                      *m_stencils[dit],
                                      ghostCellsIVS);
      Vector< BaseFab<VectorTransformation>* > vectorStencilTransformations(m_nblocks);
      for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
        {
          Box bxStencilSrc = (*m_stencilCellsMinBox[srcBlock])[dit];
          bxStencilSrc.grow(growStencilVect);
          vectorStencilTransformations[srcBlock] =
            new BaseFab<VectorTransformation>(bxStencilSrc, 1);
        }
      m_mbUtil->getVectorTransformations(vectorStencilTransformations,
                                         *m_vectorstencils[dit],
                                         m_type,
                                         ghostCellsIVS,
                                         ghostBlockNum);
      m_mbUtil->getVectorWeights(*m_vectorstencils[dit],
                                 vectorStencilTransformations,
                                 m_type,
                                 ghostCellsIVS,
                                 *validIndices[dit],
                                 *validBlock[dit],
                                 *validMappedCenter[dit],
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
MultiBlockLevelExchange::removeNoValidSource()
{
  CH_TIME("MultiBlockLevelExchange::removeNoValidSource()");
  DataIterator dit = m_grids.dataIterator();

  LevelData<FArrayBox> zeroData(m_grids, 1);
  for (dit.begin(); dit.ok(); ++dit)
    {
      zeroData[dit].setVal(0.);
    }

  Vector< BoxLayoutData<FArrayBox>* > stencilData(m_nblocks);
  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    {
      stencilData[srcBlock] =
        new BoxLayoutData<FArrayBox>(*m_validLayout[srcBlock], 1);
    }

  // Initialize stencilData with 1.,
  // then overwrite with 0. on valid cells. 
  Interval intvl0(0, 0);
  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    {
      BoxLayoutData<FArrayBox>& stencilDataSrc = *stencilData[srcBlock];
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& stencilDataSrcFab = stencilDataSrc[dit];
          stencilDataSrcFab.setVal(1.);
        }

      const Copier& srcCopier = *m_copiers[srcBlock];
      { CH_TIME("removeNoValidSource communication");
        zeroData.copyTo(intvl0, stencilDataSrc, intvl0, srcCopier);
      }
    }

  for (dit.begin(); dit.ok(); ++dit)
    {
      const IntVectSet& originalGhostsIVS = m_ghostCells[dit];
      // stencilsPatch lives on originalGhostsIVS.
      const IVSFAB<MBStencil>& stencilsPatch = *m_stencils[dit];

      // Loop over cells in originalGhostsIVS,
      // and if any need to be removed,
      // remove them from updatedGhostsIVS.
      IntVectSet updatedGhostsIVS(originalGhostsIVS);
      bool changed = false;
      for (IVSIterator ivsit(originalGhostsIVS); ivsit.ok(); ++ivsit)
        {
          IntVect thisGhostCell = ivsit();
          Real totAbsent = 0.;
          const MBStencil& thisGhostStencil = stencilsPatch(thisGhostCell, 0);
          MBStencilIterator stencilit(thisGhostStencil);
          for (stencilit.begin(); stencilit.ok(); ++stencilit)
            {
              const MBStencilElement& stencilElement = stencilit();
              const IntVect& cell = stencilElement.cell();
              int srcBlock = stencilElement.block();
              const FArrayBox& stencilDataSrc = (*(stencilData[srcBlock]))[dit];
              totAbsent += stencilDataSrc(cell, 0);
            }
          if (totAbsent > 0.)
            {
              updatedGhostsIVS -= thisGhostCell;
              // delete stencilsPatch(thisGhostCell, 0);
              changed = true;
            }
        }
      if (changed)
        {
          m_ghostCells[dit] = IntVectSet(updatedGhostsIVS);
        }
    }

  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    {
      delete stencilData[srcBlock];
    }
}


void
MultiBlockLevelExchange::interpGhosts(LevelData<FArrayBox>&  a_data,
                                      const Interval&        a_intvl) const
{
  CH_TIME("MultiBlockLevelExchange::interpGhosts");
  if (a_data.disjointBoxLayout().size() == 0)
    {
      return;
    }
  if (m_nblocks == 1)
    { // only one block; do nothing
      return;
    }
  int ncomp = a_intvl.size();
  Interval intvl0(0, ncomp-1);

  // BoxLayoutData *validData[srcBlock] lives on *m_validFullLayout[srcBlock]
  // and will hold
  // all valid data from srcBlock needed to fill m_ghostCells[dit].
  Vector< RefCountedPtr< BoxLayoutData<FArrayBox> > > validData(m_nblocks);
  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    {
      validData[srcBlock] = RefCountedPtr< BoxLayoutData<FArrayBox> >
        (new BoxLayoutData<FArrayBox>(*m_validFullLayout[srcBlock], ncomp));
    }

  // COMMUNICATE from data0 to validData.
  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    {
      BoxLayoutData<FArrayBox>& srcData = *validData[srcBlock];
      const Copier& srcCopier = *m_copiersFull[srcBlock];
      { CH_TIME("interpGhosts communication");
        a_data.copyTo(a_intvl, srcData, intvl0, srcCopier);
      }
    }

  DataIterator dit = m_grids.dataIterator();

  // From here on, do LOCAL interpolation.
  for (dit.begin(); dit.ok(); ++dit)
    {
      // ghostCellsIVS and stencilsPatch are on the collapsed layout.
      const IntVectSet& ghostCellsIVS = m_ghostCells[dit];
      const IVSFAB<MBStencil>& stencilsPatch = *m_stencils[dit];

      FArrayBox& dataFab = a_data[dit];
      // Use alias data0Fab so that component indices are always 0:ncomp-1.
      FArrayBox data0Fab(a_intvl, dataFab);
      // Note dataFab and data0Fab are on the full uncollapsed layout.

      const Box& bxFixedOff = m_gridsFixedOff[dit];
      for (IVSIterator ivsit(ghostCellsIVS); ivsit.ok(); ++ivsit)
        {
          // thisGhostCell and stencilHere are on the collapsed layout.
          IntVect thisGhostCell = ivsit();

          // Fill in data0Fab(fullGhostCell, :) where
          // fullGhostCell[m_fixedDimsVect] == thisGhostCell[m_fixedDimsVect].
          data0Fab.shift(-thisGhostCell);
          data0Fab.setVal(0., bxFixedOff, 0, ncomp);
          const MBStencil& thisGhostStencil = stencilsPatch(thisGhostCell, 0);
          MBStencilIterator stencilit(thisGhostStencil);
          for (stencilit.begin(); stencilit.ok(); ++stencilit)
            {
              const MBStencilElement& stencilElement = stencilit();
              int srcBlock = stencilElement.block();
              const IntVect& cell = stencilElement.cell();
              Real wt = stencilElement.weight();

              // non-const only because we shift it
              FArrayBox& validDataSrcFab = (*(validData[srcBlock]))[dit];

              validDataSrcFab.shift(-cell);
              data0Fab.plus(validDataSrcFab,
                            bxFixedOff, // source box
                            bxFixedOff, // dest box
                            wt, // scaling
                            0, // start source component index
                            0, // start dest component index
                            ncomp); // number of components
              // shift back to where it began
              validDataSrcFab.shift(cell);
            }
          // shift back
          data0Fab.shift(thisGhostCell);
        }
    }
}


void
MultiBlockLevelExchange::interpGhosts(LevelData<FArrayBox>&  a_data) const
{
  if (a_data.disjointBoxLayout().size() == 0)
    {
      return;
    }
  if (m_coordSysPtr->numBlocks() == 1)
    { // only one block; do nothing
      return;
    }
  const Interval& intvl = a_data.interval();
  interpGhosts(a_data, intvl);
}

void
MultiBlockLevelExchange::interpGhostsVector(LevelData<FArrayBox>&  a_data) const
{
  const Interval& intvl = a_data.interval();
  interpGhostsVector(a_data, intvl);
}


void
MultiBlockLevelExchange::interpGhostsAllWithVector(LevelData<FArrayBox>&  a_data,
                                                   const Interval&        a_vecIntvl) const
{
  CH_TIME("MultiBlockLevelExchange::interpGhostsAllWithVector");
  if (a_vecIntvl.size() == 0)
    {
      interpGhosts(a_data);
    }
  else
    {
      int compVecLo = a_vecIntvl.begin();
      int compVecHi = a_vecIntvl.end();
      const Interval& dataIntvl = a_data.interval();
      int compDataLo = dataIntvl.begin();
      int compDataHi = dataIntvl.end();
      if (compDataLo < compVecLo)
        { // scalar components before vector components
          Interval scalarIntvl(compDataLo, compVecLo-1);
          interpGhosts(a_data, scalarIntvl);
        }
      // vector components
      interpGhostsVector(a_data, a_vecIntvl);
      if (compDataHi > compVecHi)
        { // scalar components after vector components
          Interval scalarIntvl(compVecHi+1, compDataHi);
          interpGhosts(a_data, scalarIntvl);
        }
    }
}

void
MultiBlockLevelExchange::interpGhostsVector(LevelData<FArrayBox>&  a_data,
                                            const Interval&        a_intvl) const
{
  CH_assert(isDefined());
  CH_assert(m_isVectorDefined);
  CH_TIME("MultiBlockLevelExchange::interpGhostsVector");
  if (a_data.disjointBoxLayout().size() == 0)
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
  // Set data0[intvl0] == a_data[intvl].
  LevelData<FArrayBox> data0;
  aliasLevelData(data0, &a_data, a_intvl);

  // BoxLayoutData *validData[srcBlock] lives on *m_validLayout[srcBlock]
  // and will hold
  // all valid data from srcBlock needed to fill m_ghostCells[dit].
  Vector< RefCountedPtr< BoxLayoutData<FArrayBox> > > validData(m_nblocks);
  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    {
      validData[srcBlock] = RefCountedPtr< BoxLayoutData<FArrayBox> >
        (new BoxLayoutData<FArrayBox>(*m_validFullLayout[srcBlock], ncomp));
    }

  // COMMUNICATE from data0 to validData.
  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    {
      BoxLayoutData<FArrayBox>& srcData = *validData[srcBlock];
      const Copier& srcCopier = *m_copiersFull[srcBlock];
      { CH_TIME("interpGhostsVector communication");
        data0.copyTo(intvl0, srcData, intvl0, srcCopier);
      }
    }

  DataIterator dit = m_grids.dataIterator();
  // From here on, do LOCAL interpolation.
  for (dit.begin(); dit.ok(); ++dit)
    { CH_TIME("interpGhostsVector local interpolation");
      // ghostCellsIVS and stencilsPatch are on the collapsed layout.
      const IntVectSet& ghostCellsIVS = m_ghostCells[dit];
      const IVSFAB<MBVectorStencil>& vectorstencilsPatch = *m_vectorstencils[dit];

      FArrayBox& data0Fab = data0[dit];
      // Note data0Fab is on the full uncollapsed layout.

      const Box& bxFixedOff = m_gridsFixedOff[dit];
      for (IVSIterator ivsit(ghostCellsIVS); ivsit.ok(); ++ivsit)
        {
          // thisGhostCell and stencilsPatch are on the collapsed layout.
          IntVect thisGhostCell = ivsit();

          // Fill in data0Fab(fullGhostCell, :) where
          // fullGhostCell[m_fixedDimsVect] == thisGhostCell[m_fixedDimsVect].
          data0Fab.shift(-thisGhostCell);
          data0Fab.setVal(0., bxFixedOff, 0, ncomp);
          const MBVectorStencil& thisGhostVectorStencil =
            vectorstencilsPatch(thisGhostCell, 0);
          MBVectorStencilIterator vstencilit(thisGhostVectorStencil);
          for (vstencilit.begin(); vstencilit.ok(); ++vstencilit)
            {
              const MBVectorStencilElement& vstencilElement = vstencilit();
              int srcBlock = vstencilElement.block();
              const IntVect& cell = vstencilElement.cell();
              const Tuple<Real, SpaceDim*SpaceDim>& wt = vstencilElement.weight();
              // non-const only because we shift it
              FArrayBox& validDataSrcFab = (*(validData[srcBlock]))[dit];

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
                      data0Fab.plus(validDataSrcFab,
                                    bxFixedOff, // source box
                                    bxFixedOff, // dest box
                                    wt[tupleIndex], // scaling
                                    src, // source component index
                                    dst, // dest component index
                                    1); // number of components
                      tupleIndex++;
                    }
                }
              // shift back to where it began
              validDataSrcFab.shift(cell);
            }
          // shift back
          data0Fab.shift(thisGhostCell);
        }
    }
}

#include "NamespaceFooter.H"
