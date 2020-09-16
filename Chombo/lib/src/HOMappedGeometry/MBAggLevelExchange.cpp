#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MBAggLevelExchange.H"
#include "FourthOrderUtil.H"
#include "MBStencilIterator.H"
#include "CH_Timer.H"
#include "MBAggLevelExchange.H"
#include "NamespaceHeader.H"

/// constructor
MBAggLevelExchange::MBAggLevelExchange(const MultiBlockLevelGeom*  a_geomPtr,
                                       int                         a_ghosts,
                                       int                         a_order)
{
  define(a_geomPtr, a_ghosts, a_order);
}

/// destructor
MBAggLevelExchange::~MBAggLevelExchange()
{
  undefine();
}

/// destructor
void MBAggLevelExchange::undefine()
{
  CH_TIME("MBAggLevelExchange::undefine");
  if (m_isDefined)
    {
      m_powers.clear();

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

void MBAggLevelExchange::define(const MultiBlockLevelGeom*  a_geomPtr,
                                int                         a_ghosts,
                                int                         a_order)
{
  CH_TIME("MBAggLevelExchange::define");
  undefine();

  /// interpolation stencils
  LayoutData< RefCountedPtr< IVSFAB<MBStencil> > > slowStencils;

  m_geomPtr = (MultiBlockLevelGeom*) a_geomPtr;
  CH_assert(m_geomPtr != NULL);
  CH_assert(m_geomPtr->isDefined());
  m_ghosts = a_ghosts;
  m_order = a_order;
  CH_assert((m_order == 2) || (m_order == 4)); // want order 5 eventually

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
    // This radius gives us 2*floor(m_order/2)+1 cells in a row,
    // hence at least m_order+1 cells in a row,
    // which is exactly what we need to find (m_order+1)'th derivatives.
    // NEW, 22 Feb 2011: we have order either 4 or 5 only.
    m_radius = 1; // (m_order/2);

    m_degree = m_order - 1; // either 3 or 4

    // Cleared this in undefine();
    // Vector<int> m_powers;
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

  { 
    CH_TIME("MBAggLevelExchange slowstencil data allocation");
    slowStencils.define(m_grids);

    for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
      {
        const IntVectSet& ghostCellsIVS = m_ghostCells[dit];
        slowStencils[dit] = RefCountedPtr< IVSFAB<MBStencil> >
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

  Vector< RefCountedPtr< LayoutData<Box> > > stencilCellsMinBox(m_nblocks);
  Vector< RefCountedPtr< LayoutData<Box> > > stencilCellsFullMinBox(m_nblocks);
  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    {
      stencilCellsMinBox[srcBlock] = RefCountedPtr< LayoutData<Box> >
        (new LayoutData<Box>(m_grids));
      stencilCellsFullMinBox[srcBlock] = RefCountedPtr< LayoutData<Box> >
        (new LayoutData<Box>(m_grids));
    }

  // Find the stencil cells and their source blocks.
  // We will NOT find coefficients now, because those depend on centering;
  // they will be computed in the derived class that specifies centering.
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      Vector<Box> stencilsMinBox;
      m_mbUtil->getStencilCells(*slowStencils[dit],
                                stencilsMinBox,
                                m_ghostCells[dit],
                                *validIndices[dit],
                                *validBlock[dit],
                                *validMappedCenter[dit],
                                m_order,
                                m_radius);

      const Box& bxFixedOff = m_gridsFixedOff[dit];
      { CH_TIME("MBAggLevelExchange set stencilCellsMinBoxBlock");
        for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
          {
            (*stencilCellsMinBox[srcBlock])[dit] = Box(stencilsMinBox[srcBlock]);

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

            (*stencilCellsFullMinBox[srcBlock])[dit] = Box(stencilsFullMinBoxSrc);
          }
      }
    } // end loop over patches

  m_validLayout.resize(m_nblocks);
  m_validFullLayout.resize(m_nblocks);
  m_copiers.resize(m_nblocks);
  m_copiersFull.resize(m_nblocks);
  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    { CH_TIME("MBAggLevelExchange copier allocation");
      m_validLayout[srcBlock] =
        new BoxLayout(*stencilCellsMinBox[srcBlock]);
      m_validFullLayout[srcBlock] =
        new BoxLayout(*stencilCellsFullMinBox[srcBlock]);

      m_copiers[srcBlock] =
        new Copier(m_grids, *m_validLayout[srcBlock]);
      m_copiersFull[srcBlock] =
        new Copier(m_gridsFull, *m_validFullLayout[srcBlock]);
    }

  // Remove from m_ghostCells all cells that do not have complete stencils.
  removeNoValidSource(slowStencils);

  const LayoutData<int>& blockNumbers = m_geomPtr->block();
  if ( m_allGridsHaveFixedPt )
    {
      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
        {
          m_mbUtil->getWeights(*slowStencils[dit],
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
      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
        {
          int srcBlock = blockNumbers[dit];
          m_mbUtil->copyStencilsFromBlock(*slowStencils[dit],
                                          m_ghostCells[dit],
                                          *blockStencils[srcBlock]);
        }
    }

  //the interpdimsvect interface is more general than the aggstencil
  //thing can support.   we need the to have the dimensionality of the
  //the stencil be the first N components of the stencil 
  //true for all current examples 
  CH_assert(m_interpDimsVect.size() <= SpaceDim);
  CH_assert(m_interpDimsVect.size() > 0);

  if(m_interpDimsVect.size() < SpaceDim)
    {
      for(int idim = 0; idim < m_interpDimsVect.size(); idim++)
        {
          CH_assert(m_interpDimsVect[idim] < m_interpDimsVect.size());
        }
    }
  m_stencils.define(m_grids);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      CH_TIME("MGAggLevelExchange::defining agg stencil");
      const IntVectSet& ghostCellsIVS = m_ghostCells[dit];
      Box ghostBox = m_grids[dit()];
      ghostBox.grow(m_ghosts);
      int stencilDim = m_interpDimsVect.size();
      m_stencils[dit()] = RefCountedPtr<MBAggStencil>
        (new MBAggStencil(slowStencils[dit()], ghostCellsIVS, ghostBox, ghostBox, stencilDim));
    }
                        
}


void
MBAggLevelExchange::defineVector()
{
  MayDay::Error("Missing MBAggLevelExchange::defineVector()");
}


void
MBAggLevelExchange::removeNoValidSource(const LayoutData< RefCountedPtr< IVSFAB<MBStencil> > >& a_slowStencils)
{
  CH_TIME("MBAggLevelExchange::removeNoValidSource()");
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
      const IVSFAB<MBStencil>& stencilsPatch = *a_slowStencils[dit];

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
MBAggLevelExchange::interpGhosts(LevelData<FArrayBox>&  a_data,
                                 const Interval&        a_intvl) const
{
  CH_TIME("MBAggLevelExchange::interpGhosts");
  if (a_data.disjointBoxLayout().size() == 0)
    {
      return;
    }
  if (m_nblocks == 1)
    { // only one block; do nothing
      return;
    }
  int isrc = a_intvl.begin(); 
  int idst = a_intvl.begin();
  int inco = a_intvl.size();
  int ncomp = inco;
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

  // From here on, do LOCAL interpolation.
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      Vector<FArrayBox*> srcVec(m_nblocks, NULL);
      for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
        {
          BoxLayoutData<FArrayBox>& blockData = *(validData[srcBlock]);
          srcVec[srcBlock] = &(blockData[dit()]);
        }
      //the false is for incrementonly
      m_stencils[dit()]->apply(a_data[dit()], srcVec, isrc, idst, inco, false);
    }
}


void
MBAggLevelExchange::interpGhosts(LevelData<FArrayBox>&  a_data) const
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
MBAggLevelExchange::interpGhostsVector(LevelData<FArrayBox>&  a_data) const
{
  const Interval& intvl = a_data.interval();
  interpGhostsVector(a_data, intvl);
}

void
MBAggLevelExchange::interpGhostsVector(LevelData<FArrayBox>&  a_data,
                                       const Interval&        a_intvl) const
{
  if (m_coordSysPtr->numBlocks() == 1)
    { // only one block; do nothing
      return;
    }
  // In derived class.
  MayDay::Error("Missing MBAggLevelExchange::interpGhostsVector()");
}

void
MBAggLevelExchange::interpGhostsAllWithVector(LevelData<FArrayBox>&  a_data,
                                              const Interval&        a_vecIntvl) const
{
  CH_TIME("MBAggLevelExchange::interpGhostsAllWithVector");
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

#include "NamespaceFooter.H"
