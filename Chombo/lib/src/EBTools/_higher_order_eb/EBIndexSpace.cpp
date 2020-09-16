#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "parstream.H"
#include "memtrack.H"
#include "memusage.H"
#include "CH_Attach.H"
#include "IndexTM.H"
#include "BoxIterator.H"
#include "LoadBalance.H"
#include "LayoutIterator.H"
#include "BRMeshRefine.H"
#include "AMRIO.H"

#include "EBCFCopy.H"
#include "EBIndexSpace.H"
#include "EBGraphFactory.H"
#include "EBDataFactory.H"
#include "EBISLayout.H"
#include "VoFIterator.H"
#include "IrregNode.H"
#include "AllRegularService.H"
#include "PolyGeom.H"
#include "EBLevelDataOps.H"

#include "NamespaceHeader.H"

Real EBIndexSpace::s_tolerance = 1.0e-12;
bool EBIndexSpace::s_verbose   = false;
bool EBIndexSpace::s_MFSingleBox=false;
bool EBIndexSpace::s_useMemoryLoadBalance = false;

long long EBIndexSpace::numVoFs(const ProblemDomain& a_domain) const
{
  CH_TIME("EBIndexSpace::numVoFs");
  int ilev = getLevel(a_domain);
  const EBISLevel&  ebisLevel = *m_ebisLevel[ilev];
  long long numPtsLocal = ebisLevel.numVoFsOnProc();
  long long numPtsTot = EBLevelDataOps::parallelSum(numPtsLocal);
  return numPtsTot;
}

Real EBIndexSpace::totalVolFrac(const ProblemDomain& a_domain) const
{
  CH_TIME("EBIndexSpace::totalVolFrac");

  int ilev = getLevel(a_domain);
  const EBISLevel&  ebisLevel = *m_ebisLevel[ilev];

  Real totalVolFracLocal = ebisLevel.totalVolFracOnProc();
  Real totalVolFrac = EBLevelDataOps::parallelSum(totalVolFracLocal);

  return totalVolFrac;
}

bool EBIndexSpace::isDefined() const
{
  return m_isDefined;
}

int EBIndexSpace::getNCellMax() const
{
  return m_nCellMax;
}

#ifdef CH_USE_HDF5

void EBIndexSpace::readInAllLevels(HDF5Handle & a_handle,
                                   ProblemDomain a_finestLevel)
{
  CH_TIME("EBIndexSpace::define_hdf5handle_all_read_no_coarsen");

  BaseIFFAB<FaceData>::setSurroundingNodeSemantic(true);
  clear();
  m_isDefined = true;

  //read in ncellmax
  HDF5HeaderData header;
  header.readFromFile(a_handle);
  m_nCellMax          = header.m_int["EBIS_ncellMax"];
  int numLevelsInFile = header.m_int["EBIS_numLevels"];
  Box finestBoxInFile = header.m_box["EBIS_domain"];
  int startLev = -1;
  Box boxLev = finestBoxInFile;
  for (int ilev = 0; ilev < numLevelsInFile; ilev++)
    {
      if (boxLev == a_finestLevel.domainBox())
        {
          startLev = ilev;
          break;
        }
      boxLev.coarsen(2);
    }
  if (startLev == -1)
    {
      MayDay::Error("EBIS::definefromHDF5::input domain not found");
    }
  m_nlevels = numLevelsInFile - startLev;
  m_domainLevel.resize(m_nlevels);
  m_ebisLevel.resize(m_nlevels);
  ProblemDomain domLev = a_finestLevel;
  for (int ilev = 0; ilev < m_nlevels; ilev++)
    {
      //need this so tags will be right with file 
      int fileint = startLev + ilev; 
      m_ebisLevel[ilev] = new EBISLevel(a_handle, fileint);
      m_domainLevel[ilev] = domLev;
      domLev.coarsen(2);
    }
}

void EBIndexSpace::writeAllLevels(HDF5Handle&   a_handle) const
{
  CH_TIME("EBIndexSpace::write_all_levels");
  //read in ncellmax
  BaseIFFAB<FaceData>::setSurroundingNodeSemantic(true);
  HDF5HeaderData header;
  std::string filedescriptor("EBIndexSpace_allLevels");
  header.m_string ["filetype"]      = filedescriptor;
  header.m_int["EBIS_ncellMax"] = m_nCellMax;
  header.m_int["EBIS_numLevels"] = m_nlevels;
  header.m_box["EBIS_domain"] = m_domainLevel[0].domainBox();
  header.writeToFile(a_handle);
  
  for (int ilev = 0; ilev < m_nlevels; ilev++)
    {
      m_ebisLevel[ilev]->write(a_handle, ilev);
    }
  //write out finest level
  //coarser levels are derived from graph coarsening
}

void EBIndexSpace::define(HDF5Handle & a_handle,
                          int          a_maxCoarsenings)
{
  CH_TIME("EBIndexSpace::define_hdf5handle");

  pout() << "EBIndexSpace::define - From HDF5 file" << endl;

  clear();

  //read in ncellmax
  HDF5HeaderData header;
  header.readFromFile(a_handle);
  int nCellMax = header.m_int["EBIS_ncellMax"];

  pout() << "  Reading first level..." << endl;
  //read in finest level
  //coarser levels are derived from graph coarsening
  EBISLevel* level0 = new EBISLevel(a_handle);

  level0->printGraphSummary("    ");

  pout() << "  Generating coarser levels..." << endl;
  define(level0, nCellMax, a_maxCoarsenings);
}

void EBIndexSpace::write(HDF5Handle&   a_handle,
                         ProblemDomain a_outputLevel) const
{
  CH_TIME("EBIndexSpace::write");
  //read in ncellmax
  HDF5HeaderData header;
  std::string filedescriptor("EBIndexSpace");
  header.m_string ["filetype"]      = filedescriptor;
  header.m_int["EBIS_ncellMax"] = m_nCellMax;
  header.writeToFile(a_handle);

  int ilev = 0;
  if (!a_outputLevel.isEmpty())
    {
      ilev = getLevel(a_outputLevel);
    }
  //write out finest level
  //coarser levels are derived from graph coarsening
  m_ebisLevel[ilev]->write(a_handle);
}
#endif

void EBIndexSpace::define(EBISLevel * a_level0,
                          int         a_nCellMax,
                          int         a_maxCoarsenings)
{
  CH_TIME("EBIndexSpace::define_ebislevel0");

  pout() << "EBIndexSpace::define - Given level 0" << endl;

  m_nCellMax = a_nCellMax;
  m_isDefined = true;

  const ProblemDomain& fineDomain = a_level0->getDomain();

  //figure out how deep we can go
  m_nlevels = 1;
  bool canref = (fineDomain == refine(coarsen(fineDomain,2), 2));
  CH_assert(!fineDomain.isEmpty());
  ProblemDomain refbox = fineDomain;
  while (canref)
    {
      ProblemDomain refcoarbox = coarsen(refbox,2);
      refcoarbox.refine(2);
      if (refcoarbox != refbox)
        {
          canref = false;
        }
      else
        {
          m_nlevels++;
          refbox.coarsen(2);
        }
    }
  if (a_maxCoarsenings != -1)
    {
      CH_assert(a_maxCoarsenings >= 0);
      if (m_nlevels > a_maxCoarsenings+1) m_nlevels = a_maxCoarsenings + 1;
    }

  m_ebisLevel.resize(m_nlevels, NULL);
  m_domainLevel.resize(m_nlevels);

  ProblemDomain  domLevel = fineDomain;
  a_level0->clearMultiBoundaries();
  m_ebisLevel[0] = a_level0;
  m_domainLevel[0] = domLevel;
  AllRegularService dummy;
  for (int ilev = 1; ilev < m_nlevels; ilev++)
    {
      pout() << "  Generating level " << ilev << endl;
      domLevel.coarsen(2);
      m_domainLevel[ilev] = domLevel;
      m_ebisLevel[ilev] = new EBISLevel(*m_ebisLevel[ilev-1], dummy, a_nCellMax);
      m_ebisLevel[ilev]->clearMultiBoundaries();
      m_ebisLevel[ilev]->printGraphSummary("    ");
    }

#ifndef NDEBUG
  for (int ilev = 0; ilev < m_nlevels; ilev++)
    {
      m_ebisLevel[ilev]->sanityCheck(this);
    }
#endif
}

void EBIndexSpace::define(const ProblemDomain    & a_domain,
                          const RealVect         & a_origin,
                          const Real             & a_dx,
                          const GeometryService  & a_geoserver,
                          int                      a_nCellMax,
                          int                      a_maxCoarsenings)
{
  CH_TIME("EBIndexSpace::define_geoserver_domain0");

  pout() << "EBIndexSpace::define - From domain" << endl;

  pout() << "  Building finest level..." << endl;
  if (a_nCellMax > 0)
    {
      m_nCellMax = a_nCellMax;
    }
  else
    {
      if (SpaceDim == 2)
        {
          m_nCellMax = 32;
        }
      else
        {
          m_nCellMax = 32;
        }
    }

  int cellMax = a_nCellMax;
  if((cellMax < 0) && (m_nCellMax > 0))
    {
      // user did not specify a max box size so pull it from member
      cellMax = m_nCellMax;
    }
  buildFirstLevel(a_domain, a_origin, a_dx, a_geoserver, cellMax, a_maxCoarsenings);
  m_ebisLevel[0]->clearMultiBoundaries();

#ifdef CH_MPI
  {
    CH_TIME("EBIndexSpace::done_with_level_0_define_barrier");
    MPI_Barrier(Chombo_MPI::comm);
  }
#endif
  print_memory_line("ebis_after_first_level_build");
  EBISLevel* n =  m_ebisLevel[0];
  n->printGraphSummary("    ");
  pout() << endl;
  int level = 1;
  while (n)
    {
      pout() << "  Building level " << level << "..." << endl;
      n->clearMultiBoundaries();
      n=buildNextLevel(a_geoserver, cellMax);

      if (n)
        {
          n->printGraphSummary("    ");
        }
      else
      {
        pout() << "    Empty" << endl;
      }
      pout() << endl;

      level++;
    }

#ifndef NDEBUG
  for (int ilev = 0; ilev < m_nlevels; ilev++)
    {
      m_ebisLevel[ilev]->sanityCheck(this);
    }
#endif
#ifdef CH_MPI
  {
    CH_TIME("EBIndexSpace::done_with_all_define_barrier");
    MPI_Barrier(Chombo_MPI::comm);
  }
#endif
  print_memory_line("ebis_leaving_define");
}

void EBIndexSpace::define(const ProblemDomain                        & a_entireDomain,
                          const RealVect                             & a_origin,
                          const Real                                 & a_dx,
                          const Vector<RefCountedPtr<EBIndexSpace> > & a_patches,
                          const Vector<IntVect>                      & a_offsets,
                          int                                          a_maxCoarsenings)
{
  CH_TIME("EBIndexSpace::define_stitching");

  pout() << "EBIndexSpace::define - Stitching" << endl;

  CH_assert(a_patches.size() == a_offsets.size());

  int numPatches = a_patches.size();

  if (numPatches > 0)
  {
    Vector<Box> allBoxes;
    Vector<int> allProcIDs;

    pout() << "  Working on finest level..." << endl;

    int numLevels = a_patches[0]->m_nlevels;

    for (int iPatch = 0; iPatch < numPatches; iPatch++)
    {
      const RefCountedPtr<EBIndexSpace> & curPatch  = a_patches[iPatch];
      const IntVect                     & curOffset = a_offsets[iPatch];

      const DisjointBoxLayout & curPatchDBL = curPatch->m_ebisLevel[0]->m_grids;

      Vector<Box> curPatchBoxes   = curPatchDBL.boxArray();
      // Vector<int> curPatchProcIDs = curPatchDBL.procIDs();

      int numBoxes = curPatchBoxes.size();

      for (int iBox = 0; iBox < numBoxes; iBox++)
      {
        Box curBox = curPatchBoxes[iBox] + curOffset;

        if (!a_entireDomain.contains(curBox))
        {
          MayDay::Error("EBIndexSpace::stitchPatches - Entire domain doesn't contain a shifted box from a sub-domain");
        }

        curPatchBoxes[iBox] = curBox;
      }

      allBoxes.append(curPatchBoxes);
    }

    pout() << "  Total boxes: " << allBoxes.size() << endl;

    LoadBalance(allProcIDs,allBoxes);

    DisjointBoxLayout entireDBL(allBoxes,allProcIDs);

    this->m_nCellMax  = a_patches[0]->m_nCellMax;
    this->m_isDefined = a_patches[0]->m_isDefined;

    this->m_distributedData = a_patches[0]->m_distributedData;

    this->m_nlevels = numLevels;

    this->m_domainLevel.resize(numLevels);
    this->m_domainLevel[0] = a_entireDomain;

    this->m_ebisLevel.resize(numLevels);

    for (int iLevel = 0; iLevel < numLevels; iLevel++)
    {
      this->m_ebisLevel[iLevel] = NULL;
    }

    const EBISLevel* patchZeroEBISLevel = a_patches[0]->m_ebisLevel[0];
    EBISLevel* finestEBISLevel = new EBISLevel();

    finestEBISLevel->m_phase = patchZeroEBISLevel->m_phase;

    finestEBISLevel->m_grids = entireDBL;

    finestEBISLevel->m_domain = a_entireDomain;

    finestEBISLevel->m_origin = a_origin;
    finestEBISLevel->m_dx     = a_dx;

    finestEBISLevel->m_tolerance = patchZeroEBISLevel->m_tolerance;

    LevelData<EBGraph> & finestNewGraph = finestEBISLevel->m_graph;

    EBGraphFactory graphFactory(finestEBISLevel->m_domain);
    finestNewGraph.define(finestEBISLevel->m_grids,
                          patchZeroEBISLevel->m_graph.nComp(),
                          patchZeroEBISLevel->m_graph.ghostVect(),
                          graphFactory);

    for (int iPatch = 0; iPatch < numPatches; iPatch++)
    {
      const RefCountedPtr<EBIndexSpace> & curPatch  = a_patches[iPatch];
      const IntVect                     & curOffset = a_offsets[iPatch];

      const EBISLevel* curEBISLevel = curPatch->m_ebisLevel[0];
      const LevelData<EBGraph> & finestCurGraph = curEBISLevel->m_graph;

      EBGraphFactory curGraphFactory(curEBISLevel->m_domain);
      LevelData<EBGraph> offsetCurGraph;

      offsetCurGraph.define(finestCurGraph,curGraphFactory);

      for (DataIterator dit = offsetCurGraph.disjointBoxLayout().dataIterator(); dit.ok(); ++dit)
      {
        RefCountedPtr<EBGraphImplem> curGraphImplem = offsetCurGraph[dit()].m_implem;

        curGraphImplem->m_region += curOffset;
        curGraphImplem->m_domain = a_entireDomain;

        BaseFab<GraphNode> & graph = curGraphImplem->m_graph;

        graph.shift(curOffset);

        for (BoxIterator bit(graph.box()); bit.ok(); ++bit)
        {
          const IntVect & iv = bit();

          GraphNode & graphNode = graph(iv);

          if (graphNode.isIrregular())
          {
            int numVoFs = graphNode.m_cellList->size();

            if (numVoFs > 1)
            {
              MayDay::Error("EBIndexSpace::stitchPatches - EBIndexSpace's being stitched can't have multivalued cells");
            }

            GraphNodeImplem & graphNodeImplem = (*(graphNode.m_cellList))[0];

            for (int idir = 0; idir < SpaceDim; idir++)
            {
              for (SideIterator sit; sit.ok(); ++sit)
              {
                const Side::LoHiSide side = sit();

                IntVect ivOther(iv);
                ivOther.shift(idir,sign(side));

                int arcValue = -1;
                if (a_entireDomain.contains(ivOther))
                {
                  arcValue = 0;
                }

                Vector<int>& arcs = graphNodeImplem.m_arc[graphNodeImplem.index(idir,side)];
                int numArcs = arcs.size();

                if (numArcs > 1)
                {
                  MayDay::Error("EBIndexSpace::stitchPatches - EBIndexSpace's being stitched can't have multiple faces on a side");
                }

                if (numArcs == 1)
                {
                  arcs[0] = arcValue;
                }
              }
            }
          }
        }

        if (curGraphImplem->m_irregIVS != NULL)
        {
          curGraphImplem->m_irregIVS->shift(curOffset);
        }

        if (curGraphImplem->m_multiIVS != NULL)
        {
          curGraphImplem->m_multiIVS->shift(curOffset);
        }

        curGraphImplem->m_mask.shift(curOffset);
      }

      bool exchange = false;
      Copier offsetCopier(offsetCurGraph.disjointBoxLayout(),
                          finestNewGraph.disjointBoxLayout(),
                          exchange,
                          curOffset);
      offsetCurGraph.copyTo(finestNewGraph,offsetCopier);
    }

#if 0
    LevelData<EBGraph> ghostGraph(finestEBISLevel->m_grids,
                                  1,
                                  IntVect::Unit,
                                  graphFactory);

    Interval interval(0,0);
    finestEBISLevel->m_graph.copyTo(interval,ghostGraph,interval);

    EBDataFactory dataFactory;
    finestEBISLevel->m_data.define(finestEBISLevel->m_grids,
                                   1,
                                   IntVect::Zero,
                                   dataFactory);

    for (DataIterator dit = finestEBISLevel->m_grids.dataIterator(); dit.ok(); ++dit)
    {
      finestEBISLevel->m_data[dit()].defineVoFData(ghostGraph[dit()],
                                                   finestEBISLevel->m_grids.get(dit()));
      finestEBISLevel->m_data[dit()].defineFaceData(ghostGraph[dit()],
                                                    finestEBISLevel->m_grids.get(dit()));
    }

    for (int iPatch = 0; iPatch < numPatches; iPatch++)
    {
      const RefCountedPtr<EBIndexSpace> & curPatch  = a_patches[iPatch];
      const IntVect                     & curOffset = a_offsets[iPatch];

      const EBISLevel* curEBISLevel = curPatch->m_ebisLevel[0];
      const LevelData<EBData> & finestCurData = curEBISLevel->m_data;

      EBGraphFactory curGraphFactory(curEBISLevel->m_domain);
      LevelData<EBGraph> offsetCurGraph;

      offsetCurGraph.define(finestCurGraph,curGraphFactory);
    }
#endif

    finestEBISLevel->m_cache.clear();

    finestEBISLevel->m_cacheMisses = 0;
    finestEBISLevel->m_cacheHits   = 0;
    finestEBISLevel->m_cacheStale  = 0;

    this->m_ebisLevel[0] = finestEBISLevel;
  }
}

EBISLevel* EBIndexSpace::buildFirstLevel(const ProblemDomain&   a_domain,
                                         const RealVect&        a_origin,
                                         const Real&            a_dx,
                                         const GeometryService& a_geoserver,
                                         int a_nCellMax,
                                         int a_maxCoarsenings,
                                         bool a_fixRegularNextToMultiValued)
{
  CH_TIME("EBIndexSpace::buildFirstLevel");
  clear();
  m_isDefined = true;
  CH_assert(a_nCellMax > 0);

  m_nlevels = 1;
  bool canref = (a_domain == refine(coarsen(a_domain,2), 2));

  CH_assert(!a_domain.isEmpty());
  ProblemDomain refbox = a_domain;

  while (canref)
    {
      ProblemDomain refcoarbox = coarsen(refbox,2);
      refcoarbox.refine(2);
      if (refcoarbox != refbox)
        {
          canref = false;
        }
      else
        {
          m_nlevels++;
          refbox.coarsen(2);
        }
    }
  if (a_maxCoarsenings != -1)
    {
      CH_assert(a_maxCoarsenings >= 0);
      if (m_nlevels > a_maxCoarsenings+1) m_nlevels = a_maxCoarsenings + 1;
    }

  m_ebisLevel.resize(m_nlevels, NULL);
  m_domainLevel.resize(m_nlevels);

  ProblemDomain  domLevel = a_domain;
  m_ebisLevel[0] = new EBISLevel(domLevel,
                                 a_origin,
                                 a_dx,
                                 a_geoserver,
                                 m_nCellMax,
                                 a_fixRegularNextToMultiValued);
  m_domainLevel[0] = domLevel;
  return m_ebisLevel[0];
}

void EBIndexSpace::resetLevels(int nLevel)
{
  CH_TIME("EBIndexSpace::resetLevels");

  CH_assert(nLevel <= m_nlevels);

  for (int i = nLevel; i < m_nlevels; i++)
    {
      if (m_ebisLevel[i] != NULL) delete m_ebisLevel[i];
      m_ebisLevel[i] = NULL;
    }

  m_nlevels = nLevel;

}

EBISLevel* EBIndexSpace::buildNextLevel(const GeometryService & a_geoserver,
                                        int                     a_nCellMax,
                                        bool                    a_fixRegularNextToMultiValued)
{
  CH_TIME("buildNextLevel");
  int ilev=0;
  for ( ; ilev <m_ebisLevel.size(); ++ilev)
    {
      if (m_ebisLevel[ilev] == NULL) break;
    }
  if (ilev == m_ebisLevel.size()) return NULL;

  {
 
    m_domainLevel[ilev] = m_domainLevel[ilev-1];
    m_domainLevel[ilev].coarsen(2);
    m_ebisLevel[ilev] = new EBISLevel(*m_ebisLevel[ilev-1],
                                      a_geoserver,
                                      a_nCellMax,
                                      a_fixRegularNextToMultiValued);

 
  }

  return m_ebisLevel[ilev];
}

#ifdef CH_USE_HDF5
void EBIndexSpace::writeInfo(HDF5Handle& handle) const
{
  Vector<DisjointBoxLayout> vectGrids(m_nlevels);
  Vector<LevelData<FArrayBox>* >  vectData(m_nlevels);
  Vector<string> vectNames(1,"unknown");
  ProblemDomain domain = m_domainLevel[m_nlevels-1];
  const EBISLevel& level = *(m_ebisLevel[m_nlevels-1]);
  Real dx = level.getDX();
  Vector<int> vectRatio(m_nlevels, 2);

  for (int i = 0; i < m_nlevels; ++i)
    {
      const EBISLevel& level = *(m_ebisLevel[m_nlevels-1-i]);
      vectGrids[i] = level.getGrids();
      vectData[i] = new LevelData<FArrayBox>(level.getGrids(), 1);
      LevelData<FArrayBox>& ld = *(vectData[i]);
      for (DataIterator dit = ld.dataIterator(); dit.ok(); ++dit)
        {
          ld[dit].setVal(0.0);
        }
    }


  WriteAMRHierarchyHDF5(handle, vectGrids, vectData, vectNames, domain.domainBox(),
                        dx, 1.0, 0.0 , vectRatio, m_nlevels);

  for (int i = 0; i < m_nlevels; ++i)
    {

      delete vectData[i];

    }
}
#endif

int EBIndexSpace::numLevels() const
{
  return m_nlevels;
}

void EBIndexSpace::clear()
{
  for (int ilev = 0; ilev < m_ebisLevel.size(); ilev++)
    {
      delete m_ebisLevel[ilev];
      m_ebisLevel[ilev] = NULL;
    }
  m_ebisLevel.resize(0);
  m_domainLevel.resize(0);
  m_nlevels = 0;
  m_isDefined = false;
}

EBIndexSpace::EBIndexSpace()
{
  m_distributedData = false;
  //appalling hack to get stuff working again
  BaseIFFAB<FaceData>::setSurroundingNodeSemantic(false);
}

EBIndexSpace::~EBIndexSpace()
{
  clear();
}

int EBIndexSpace::getLevel(const ProblemDomain& a_domain) const
{
  bool found = false;
  int whichlev = -1;
  for (int ilev = 0; ilev < m_domainLevel.size() && !found; ilev++)
    {
      if (m_domainLevel[ilev].domainBox() == a_domain.domainBox())
        {
          found = true;
          whichlev = ilev;
        }
    }
  return whichlev;
}

void EBIndexSpace::fillEBISLayout(EBISLayout&              a_ebisLayout,
                                  const DisjointBoxLayout& a_grids,
                                  const ProblemDomain&     a_domain,
                                  const int&               a_nghost) const
{
  CH_assert(isDefined());
  CH_TIME("EBIndexSpace::fillEBISLayout");

  //figure out which level we are on
  int whichlev = getLevel(a_domain);
  if (whichlev < 0)
    {
      pout() << "a_domain = " << a_domain
             << " does not correspond to any refinement of EBIS" << endl;
      MayDay::Error("Bad argument to EBIndexSpace::fillEBISLayout");
    }
  m_ebisLevel[whichlev]->fillEBISLayout(a_ebisLayout, a_grids, a_nghost);
  a_ebisLayout.setEBIS(this); //need for mf
}

// Divide the EBIndexSpace into connected components and return a Vector of
// disjoint EBIndexSpace's corresponding to each connected component.
Vector<RefCountedPtr<EBIndexSpace> > EBIndexSpace::findConnectedComponents(int        & a_numComponents,
                                                                           const bool & a_onlyBiggest)
{
  CH_TIME("EBIndexSpace::findConnectedComponents");
  pout() << "    EBIndexSpace::findConnectedComponents" << endl;

  // Begin by setting up some data holders and other infrastructure to number
  // the connected components and remember the renumbering of the VoFs.
  Vector<LevelData<EBCellFAB>* > numberedComponentses(m_nlevels,NULL);
  Vector<LevelData<EBCellFAB>* > newNumberings(m_nlevels,NULL);

  // EBISLayout's for each level
  Vector<EBISLayout> ebisLayouts(m_nlevels);

  // All our temporary data structures have one component
  int nComps = 1;

  // All our temporary data structures have one layer of ghostcells
  int nGhosts = 1;
  IntVect ghostCells = nGhosts * IntVect::Unit;

  // Fill all the EBISLayouts and allocate/define all the local data holders
  pout() << "    EBIndexSpace::connectedComponents-loop1" << endl;
  for (int ilev = 0; ilev < m_nlevels; ilev++)
  {
    CH_TIME("EBIndexSpace::connectedComponents-loop1");
    pout() << "      Level " << ilev << endl;

    EBISLevel* curEBISLevel = m_ebisLevel[ilev];

    DisjointBoxLayout curGrids = curEBISLevel->m_grids;

    EBISLayout& curEBISLayout = ebisLayouts[ilev];
    curEBISLevel->fillEBISLayout(curEBISLayout,curGrids,nGhosts);

    EBCellFactory curFactory(curEBISLayout);

    LevelData<EBCellFAB>* numberedComponents = new LevelData<EBCellFAB>();
    LevelData<EBCellFAB>* newNumbering       = new LevelData<EBCellFAB>();

    numberedComponents->define(curGrids,nComps,ghostCells,curFactory);
    newNumbering      ->define(curGrids,nComps,ghostCells,curFactory);

    numberedComponentses[ilev] = numberedComponents;
    newNumberings       [ilev] = newNumbering;
  }

  // Get the coarsest level and number the connected components at that level.
  // The current implementation assumes there is only one box at the coarsest
  // level - this is checked below.
  EBISLevel* coarEBISLevel = m_ebisLevel[m_nlevels-1];

  DisjointBoxLayout coarGrids = coarEBISLevel->m_grids;

  // This really only needs to be an "int" and not a "Real" but all our
  // infrastructure (e.g., I/O, debugging, library functions) is set up to
  // work with "EBCellFAB".  This can be changed if needed to be an "int" data
  // holder.
  LevelData<EBCellFAB>& coarNumberedComponents = *numberedComponentses[m_nlevels-1];

  // Initialize all the component numbers to -1 (invalid)
  EBLevelDataOps::setVal(coarNumberedComponents,-1.0);

  // Currently this code will only work if there is only one box at the
  // coarsest level of the geometry.  This will need to be fixed by merging
  // different numberings of the same component in different boxes.  This can
  // be done via ghost cells, exchanges, the building of equivalence classes,
  // and, finally, a renumbering operation.
  unsigned int numGrids = coarGrids.size();

  if (numGrids > 1)
  {
    pout() << "Number of coarsest grids: " << numGrids << endl;
    pout() << "Grids: " << coarGrids << endl;
    MayDay::Error("EBIndexSpace::connectedComponents - algorithm currently only works if there is only one box at the coarsest geometry level");
  }

  // For counting the new EBIndexspace's and remember the number assigned to
  // each
  int numEBIS = 0;
  Vector<int> valueEBIS;

  // Used if "a_onlyBiggest" is true to get the biggest connected component
  Real maxVolFrac = 0.0;
  int biggestEBIS = -1;

  // Go through all (one) boxes (see comments and MayDay above).
  pout() << "    EBIndexSpace::connectedComponents-loop2" << endl;
  for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
  {
    CH_TIME("EBIndexSpace::connectedComponents-loop2");

    // Get a bunch of local references
    int boxIndex = coarGrids.index(dit());

    const Box& curBox = coarGrids.get(dit());
    const IntVectSet curIVS(curBox);
    pout() << "      Box: " << curBox << endl;

    EBCellFAB&     curEBCellFAB = coarNumberedComponents[dit()];
    const EBGraph& curEBGraph   = coarEBISLevel->m_graph[dit()];
    const EBISBox& curEBISBox   = curEBCellFAB.getEBISBox();

    unsigned int curNum = boxIndex;

    // Iterate through the box, VoF by VoF and number anything that isn't
    // already numbered.
    for (VoFIterator vofit(curIVS,curEBGraph); vofit.ok(); ++vofit)
    {
      // Current VoF
      const VolIndex& vof = vofit();

      Real totalVolFrac;

      // Examine the connected component starting at "vof".  If it hasn't been
      // numbered yet, return true otherwise false.  Also, return the
      // volume fraction associated with the connect component.
      if (setAllConnectedVoFs(totalVolFrac,curEBCellFAB,curEBGraph,curEBISBox,vof,vof,curNum))
      {
        if (totalVolFrac > 0.0)
        {
          // If a connected component was numbered, check to see for the
          // numbering can be incremented and then increment it.
          if (curNum + numGrids < curNum)
          {
            MayDay::Error("EBIndexSpace::connectedComponents - Component index overflow");
          }

          // Remember the component number and increment the number of
          // components
          valueEBIS.push_back(curNum);
          numEBIS++;

          // Find the biggest connected component
          if (totalVolFrac > maxVolFrac)
          {
            maxVolFrac = totalVolFrac;
            biggestEBIS = curNum;
          }

          // New number of the next component (if there is one)
          curNum += numGrids;
        }
        else
        {
          // If we're not counting the current component (because it has a total
          // volume fraction of 0.0, set the component number back to -1
          resetAllConnectedVoFs(curEBCellFAB,curEBGraph,curEBISBox,vof,vof);
        }
      }
    }
  }

#ifdef CH_MPI
  // Sum up the numEBIS on all processors
  int sumNumEBIS;

  int result = MPI_Allreduce(&numEBIS, &sumNumEBIS, 1, MPI_INT, MPI_SUM, Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
  {
    MayDay::Error("Communication error summing 'numEBIS'");
  }

  numEBIS = sumNumEBIS;

  typedef struct
  {
    double maxVolFrac;
    int    biggestEBIS;
  } CH_MPI_MAXLOC;

  CH_MPI_MAXLOC localMax;
  CH_MPI_MAXLOC globalMax;

  // Get the connected component index of the largest connected component on
  // all processors
  localMax.maxVolFrac  = maxVolFrac;
  localMax.biggestEBIS = biggestEBIS;

  result = MPI_Allreduce(&localMax, &globalMax, 1, MPI_DOUBLE_INT, MPI_MAXLOC, Chombo_MPI::comm);

  if (result != MPI_SUCCESS)
  {
    MayDay::Error("Communication error maximizing 'maxVolFrac'");
  }

  maxVolFrac  = globalMax.maxVolFrac;
  biggestEBIS = globalMax.biggestEBIS;
#endif

  a_numComponents = numEBIS;

  int minEBIS;
  int maxEBIS;

  if (a_onlyBiggest)
  {
    // If we only want the biggest connected componnent, only iterate over
    // that components index
    if (biggestEBIS >= 0)
      {
        minEBIS = biggestEBIS;
        maxEBIS = biggestEBIS;
      }
    else
      {
        // This corresponds to an empty EBIS
        minEBIS =  0;
        maxEBIS = -1;
      }
  }
  else
  {
    // Otherwise, iterate over all the connected component indices
    minEBIS = 0;
    maxEBIS = numEBIS-1;
  }

  // Make sure all ghostcells are correct
  coarNumberedComponents.exchange();

  // At this point, renumber the connected components from 0 to numEBIS-1.
  // Currently, this is already the case so nothing happens here.

  // There is only one component in each temporary data holder
  int comp = 0;

  // The refinement ratio on the geometry size of things is always 2
  int nRef = 2;

  // Pass the numbering for connected components from coarse to fine levels.
  pout() << "    EBIndexSpace::connectedComponents-loop3" << endl;
  for (int ilev = m_nlevels-2; ilev >= 0; ilev--)
  {
    CH_TIME("EBIndexSpace::connectedComponents-loop3");
    pout() << "      Level " << ilev << endl;

    // Get a lot of local references to get things ready to go
    EBISLevel* coarEBISLevel = m_ebisLevel[ilev+1];
    EBISLevel* fineEBISLevel = m_ebisLevel[ilev];

    DisjointBoxLayout coarGrids = coarEBISLevel->m_grids;
    DisjointBoxLayout fineGrids = fineEBISLevel->m_grids;

    const ProblemDomain& coarDomain = coarEBISLevel->m_domain;

    EBISLayout& coarEBISLayout = ebisLayouts[ilev+1];
    EBISLayout& fineEBISLayout = ebisLayouts[ilev];

    LevelData<EBCellFAB>& coarNumberedComponents = *numberedComponentses[ilev+1];
    LevelData<EBCellFAB>& fineNumberedComponents = *numberedComponentses[ilev];


    // Initial all the numberings on this level to -1 (invalid)
    EBLevelDataOps::setVal(fineNumberedComponents,-1.0);

    // Set up a piecewise constant interpolator from the current coarse level
    // to the current fine level
    EBCFCopy constantInterp(fineGrids,
                            coarGrids,
                            fineEBISLayout,
                            coarEBISLayout,
                            coarDomain,
                            nRef,
                            nComps,
                            this,
                            ghostCells);

    // Do the interpolation/copy
    Interval oneComp(comp,comp);
    constantInterp.copy(fineNumberedComponents,
                        coarNumberedComponents,
                        oneComp);

    // Make sure all ghostcells are correct
    fineNumberedComponents.exchange();
  }

  // The set of EBIndexSpace's that are each connected components of current
  // EBIndexSpace.
  Vector<RefCountedPtr<EBIndexSpace> > connectedEBIS(numEBIS);

  // Make copies of the main EBIndexSpace for each connected component.
  pout() << "    EBIndexSpace::connectedComponents-loop4" << endl;
  for (int iEBIS = minEBIS; iEBIS <= maxEBIS; iEBIS++)
  {
    CH_TIME("EBIndexSpace::connectedComponents-loop4");
    pout() << "      EBIS #" << iEBIS << endl;

    // Make a new, empty EBIndexSpace for the current component
    RefCountedPtr<EBIndexSpace> copyEBIS(new EBIndexSpace());

    // Set the member data that is not level dependent
    copyEBIS->m_nCellMax  = m_nCellMax;
    copyEBIS->m_isDefined = m_isDefined;

    copyEBIS->m_distributedData = m_distributedData;

    copyEBIS->m_domainLevel = m_domainLevel;
    copyEBIS->m_nlevels     = m_nlevels;

    // Go through the levels setting the level dependent member data
    for (int ilev = 0; ilev < m_nlevels; ilev++)
    {
      // Original level dependent member data
      EBISLevel* origEBISLevel = m_ebisLevel[ilev];

      // New level dependent member data for this component
      EBISLevel* copyEBISLevel = new EBISLevel();

      // Copy everything except the EBData (this happens at the end)
      copyEBISLevel->m_phase = iEBIS;

      copyEBISLevel->m_grids  = origEBISLevel->m_grids;

      copyEBISLevel->m_domain = origEBISLevel->m_domain;

      copyEBISLevel->m_origin = origEBISLevel->m_origin;
      copyEBISLevel->m_dx     = origEBISLevel->m_dx;

      copyEBISLevel->m_tolerance = origEBISLevel->m_tolerance;

      EBGraphFactory graphFactory(origEBISLevel->m_domain);
      copyEBISLevel->m_graph.define(origEBISLevel->m_graph,graphFactory);

      // Clear this cache and all statistics associated with it
      copyEBISLevel->m_cache.clear();

      copyEBISLevel->m_cacheMisses = 0;
      copyEBISLevel->m_cacheHits   = 0;
      copyEBISLevel->m_cacheStale  = 0;

      // Remember this EBISLevel
      copyEBIS->m_ebisLevel.push_back(copyEBISLevel);
    }

    // Remember this EBIndexSpace
    connectedEBIS[iEBIS] = copyEBIS;
  }

  // Are there any cells which are regular with a multivalued parent
  bool anyRegularWithMultivaluedParent = false;

  // First pass to correct the EBIndexSpace for each conneceted component -
  // Remove AllRegular patches and regular cells that are in other connected
  // components, and mark all multiVoFs that are in other connected components
  pout() << "    EBIndexSpace::connectedComponents-loop5" << endl;
  for (int ilev = 0; ilev < m_nlevels; ilev++)
  {
    CH_TIME("EBIndexSpace::connectedComponents-loop5");
    pout() << "      Level " << ilev << endl;

    // The original EBISLevel
    EBISLevel* origEBISLevel = m_ebisLevel[ilev];

    // All the component numbers for this level
    const LevelData<EBCellFAB>& numberedComponents = *numberedComponentses[ilev];

    // All the EBISLevel's for the connected components
    Vector<EBISLevel*> copyEBISLevel(numEBIS);
    for (int iEBIS = minEBIS; iEBIS <= maxEBIS; iEBIS++)
    {
      copyEBISLevel[iEBIS] = connectedEBIS[iEBIS]->m_ebisLevel[ilev];
    }

    // Iterate through all the boxes
    for (DataIterator dit = origEBISLevel->m_grids.dataIterator(); dit.ok(); ++dit)
    {
      // Get some local references
      const Box& curBox = origEBISLevel->m_grids.get(dit());

      const EBCellFAB& numberedComponentsFAB = numberedComponents[dit()];

      const EBGraphImplem& origEBGraph = *(origEBISLevel->m_graph[dit()].m_implem);

      // If the original graph for this box was all regular
      if (origEBGraph.m_tag == EBGraphImplem::AllRegular)
      {
        const IntVect& oneIVInside = curBox.smallEnd();
        const int vofID = 0;
        const VolIndex vof(oneIVInside,vofID);

        // Get the component number for this (all regular) box
        Real componentNumber = numberedComponentsFAB(vof,comp);

        // Go through all the connected components and mark their graphs for
        // this box as "all covered" if they don't match "componentNumber"
        for (int iEBIS = minEBIS; iEBIS <= maxEBIS; iEBIS++)
        {
          if (iEBIS != componentNumber)
          {
            EBGraphImplem& copyEBGraph = *(copyEBISLevel[iEBIS]->m_graph[dit()].m_implem);
            copyEBGraph.m_tag = EBGraphImplem::AllCovered;
          }
        }
      }
      // If the original graph has irregular cells in this box
      else if (origEBGraph.m_tag == EBGraphImplem::HasIrregular)
      {
        // Iterate through each IntVect in this box
        IntVectSet curIVS(curBox);
        for (IVSIterator ivsit(curIVS); ivsit.ok(); ++ivsit)
        {
          // Current IntVect and it's GraphNode
          const IntVect& iv = ivsit();
          const GraphNode& origGraphNode = origEBGraph.m_graph(iv);

          // If the original GraphNode is all regular with a single valued
          // parent
          if (origGraphNode.isRegularWithSingleValuedParent())
          {
            const int vofID = 0;
            const VolIndex vof(iv,vofID);

            // Get the component number for this "vof"
            Real componentNumber = numberedComponentsFAB(vof,comp);

            // Go through all the connected components and mark their "vof" as
            // covered if they don't match "componentNumber"
            for (int iEBIS = minEBIS; iEBIS <= maxEBIS; iEBIS++)
            {
              if (iEBIS != componentNumber)
              {
                EBGraphImplem& copyEBGraph = *(copyEBISLevel[iEBIS]->m_graph[dit()].m_implem);
                copyEBGraph.m_graph(iv).m_cellList = (Vector<GraphNodeImplem>*)0;
              }
            }
          }
          // If the original GraphNode is anything else except covered
          else if (origGraphNode.hasValidCellList())
          {
            // Iterate through all the cellIndex's in the original GraphNode
            for (int iVoF = 0; iVoF < origGraphNode.size(); iVoF++)
            {
              // The current VoF
              const VolIndex vof(iv,iVoF);

              // Get the component number for this "vof"
              Real componentNumber = numberedComponentsFAB(vof,comp);

              // For each new EBIndexSpace, go to the GraphNodeImplem
              // corresponding to the current VoF and mark it "valid" if it
              // matches the current "componentNumber", otherwise mark it
              // "invalid"
              for (int iEBIS = minEBIS; iEBIS <= maxEBIS; iEBIS++)
              {
                EBGraphImplem&   copyEBGraph         = *(copyEBISLevel[iEBIS]->m_graph[dit()].m_implem);
                GraphNode&       copyGraphNode       = copyEBGraph.m_graph(iv);
                GraphNodeImplem& copyGraphNodeImplem = (*(copyGraphNode.m_cellList))[iVoF];

                if (copyGraphNodeImplem.m_isRegular)
                {
                  anyRegularWithMultivaluedParent = true;
                }

                if (iEBIS != componentNumber)
                {
                  copyGraphNodeImplem.m_isValid = false;
                }
                else
                {
                  copyGraphNodeImplem.m_isValid = true;
                }
              }
            }
          }
        }
      }
    }
  }

  // Second pass to fill the data holder with the new numbering for each
  // multiVoF in each connected component.  These can coexist in one (old)
  // data holder because the connected components are disjoint subsets of
  // the old graph.
  pout() << "    EBIndexSpace::connectedComponents-loop6" << endl;
  for (int ilev = 0; ilev < m_nlevels; ilev++)
  {
    CH_TIME("EBIndexSpace::connectedComponents-loop6");
    pout() << "      Level " << ilev << endl;

    // Data holder for the new numbering on this level
    LevelData<EBCellFAB>& newNumbering = *newNumberings[ilev];

    // Initial all the numberings on this level to 0
    EBLevelDataOps::setVal(newNumbering,0.0);

    // Go through each of the connected components
    for (int iEBIS = minEBIS; iEBIS <= maxEBIS; iEBIS++)
    {
      // Get the current EBIndexSpace and EBISLevel
      EBIndexSpace& curEBIS = *(connectedEBIS[iEBIS]);
      EBISLevel* curEBISLevel = curEBIS.m_ebisLevel[ilev];

      // Iterate through all the boxes on this level
      for (DataIterator dit = curEBISLevel->m_grids.dataIterator(); dit.ok(); ++dit)
      {
        // Get some local references
        const Box& curBox = curEBISLevel->m_grids.get(dit());

        EBCellFAB& newNumberingFAB = newNumbering[dit()];

        const EBGraphImplem& curEBGraph = *(curEBISLevel->m_graph[dit()].m_implem);

        // If the current graph has irregular cells then some renumbering may
        // be needed
        if (curEBGraph.m_tag == EBGraphImplem::HasIrregular)
        {
          // Iterate through all IntVect's in the current box
          IntVectSet curIVS(curBox);
          for (IVSIterator ivsit(curIVS); ivsit.ok(); ++ivsit)
          {
            const IntVect& iv = ivsit();
            const GraphNode& curGraphNode = curEBGraph.m_graph(iv);

            // If the graph has an explicit entry, i.e. a valid cellList, then
            // some changes may be needed
            if (curGraphNode.hasValidCellList())
            {
              // Go through the VoFs, if they are valid then record the
              // current validVoFCount (this VoF's new number in this
              // connected component) in the correct location in the data
              // holder based on the original graph and increment
              // validVoFCount.
              int validVoFCount = 0;
              for (int iVoF = 0; iVoF < curGraphNode.size(); iVoF++)
              {
                GraphNodeImplem& curGraphNodeImplem = (*(curGraphNode.m_cellList))[iVoF];
                if (curGraphNodeImplem.m_isValid == true)
                {
                  const VolIndex vof(iv,iVoF);

                  newNumberingFAB(vof,comp) = validVoFCount;
                  validVoFCount++;
                }
              }
            }
          }
        }
      }
    }

    // Make sure the ghostcells are correct on this level
    newNumbering.exchange();
  }

  // Third pass to correct the EBIndexSpace for each conneceted component on
  // each level (but not between levels) - Go through correcting all
  // references in the connected component using the map created in the second
  // pass.  Then prune away all parts of the graph that are not needed for
  // this connected component.
  pout() << "    EBIndexSpace::connectedComponents-loop7" << endl;
  for (int iEBIS = minEBIS; iEBIS <= maxEBIS; iEBIS++)
  {
    CH_TIME("EBIndexSpace::connectedComponents-loop7");
    pout() << "      EBIS #" << iEBIS << endl;

    EBIndexSpace& curEBIS = *(connectedEBIS[iEBIS]);

    // Go through each level
    for (int ilev = 0; ilev < m_nlevels; ilev++)
    {
      // Get the EBISLevel and new numbering for this lelvel
      EBISLevel* curEBISLevel = curEBIS.m_ebisLevel[ilev];
      LevelData<EBCellFAB>& newNumbering = *newNumberings[ilev];

      // Go through all the boxes in this level
      for (DataIterator dit = curEBISLevel->m_grids.dataIterator(); dit.ok(); ++dit)
      {
        const Box& curBox = curEBISLevel->m_grids.get(dit());
        EBCellFAB& newNumberingFAB = newNumbering[dit()];
        EBGraphImplem& curEBGraph = *(curEBISLevel->m_graph[dit()].m_implem);

        // If the graph in this box contain irregular cells, there may to
        // something to do
        if (curEBGraph.m_tag == EBGraphImplem::HasIrregular)
        {
          // Assume the new graph will be all covered unless we find a valid
          // VoF still in this box
          bool allCovered = true;

          // Iterator through all the IntVect's in this box looking for a
          // valid VoF
          IntVectSet curIVS(curBox);
          for (IVSIterator ivsit(curIVS); ivsit.ok(); ++ivsit)
          {
            const IntVect& iv = ivsit();
            GraphNode& curGraphNode = curEBGraph.m_graph(iv);

            // Nothing to fix but the graph isn't all covered
            if (curGraphNode.isRegularWithSingleValuedParent())
            {
              allCovered = false;
            }
            // If there are VoF's, go through them looking for valid ones and
            // correct the arcs connected that VoF to adjacent VoF's
            else if (curGraphNode.hasValidCellList())
            {
              // Get the old cellList and start building the new cellList
              Vector<GraphNodeImplem>* oldCellList = curGraphNode.m_cellList;
              Vector<GraphNodeImplem>* newCellList = new Vector<GraphNodeImplem>;

              // Go through all the VoF's
              for (int iVoF = 0; iVoF < curGraphNode.size(); iVoF++)
              {
                GraphNodeImplem& oldGraphNodeImplem = (*oldCellList)[iVoF];

                // If this VoF is valid correct the arcs to adjacent VoF's
                if (oldGraphNodeImplem.m_isValid == true)
                {
                  // Look in each direction at the lo and hi sides
                  for (int idir = 0; idir < SpaceDim; idir++)
                  {
                    for (SideIterator sit; sit.ok(); ++sit)
                    {
                      // Correct (in place) all the arcs for all the faces
                      // on all the sides
                      int index = oldGraphNodeImplem.index(idir,sit());

                      IntVect otherIV = iv;
                      otherIV.shift(idir,sign(sit()));

                      Vector<int>& curArcs = oldGraphNodeImplem.m_arc[index];
                      int numArcs = curArcs.size();

                      for (int iFace = 0; iFace < numArcs; iFace++)
                      {
                        if (curArcs[iFace] != -1)
                        {
                          VolIndex otherVoF(otherIV,curArcs[iFace]);
                          curArcs[iFace] = newNumberingFAB(otherVoF,comp);
                        }
                      }
                    }
                  }

                  // Save the new VoF
                  newCellList->push_back(oldGraphNodeImplem);
                }
              }

              // Delete the old cellList and save the new one
              delete oldCellList;
              oldCellList = newCellList;

              // Check the new number of VoF's
              int numCells = oldCellList->size();

              // If there are no VoF's, delete the cellList, mark the IntVect
              // location as covered, and remove the IntVect from the
              // irregular and multivalued IntVectSet's
              if (numCells == 0)
              {
                delete oldCellList;
                oldCellList = ((Vector<GraphNodeImplem>*) 0);

                (*curEBGraph.m_irregIVS) -= iv;
                (*curEBGraph.m_multiIVS) -= iv;
              }
              // If there are some VoF's
              else
              {
                // The graph in this box isn't all covered
                allCovered = false;

                // If there is only one VoF
                if (numCells == 1)
                {
                  // And it's not regular
                  if ((*oldCellList)[0].m_isRegular == false)
                  {
                    // Add it to the irregular IntVectSet and remove it from
                    // the multivalued IntVectSet
                    (*curEBGraph.m_irregIVS) |= iv;
                    (*curEBGraph.m_multiIVS) -= iv;
                  }
                }
              }

              // Save the updated cellList
              curGraphNode.m_cellList = oldCellList;
            }
          }

          // If the graph in this box (for this connected component) has
          // become all covered clean things up
          if (allCovered)
          {
            // Mark it as all covered
            curEBGraph.m_tag = EBGraphImplem::AllCovered;

            // Clear all the GraphNode's
            curEBGraph.m_graph.clear();

            // Delete the irregular and multivalued IntVectSet's
            if (curEBGraph.m_irregIVS != NULL)
            {
              delete curEBGraph.m_irregIVS;
              curEBGraph.m_irregIVS = NULL;
            }

            if (curEBGraph.m_multiIVS != NULL)
            {
              delete curEBGraph.m_multiIVS;
              curEBGraph.m_multiIVS = NULL;
            }
          }
        }
      }
    }
  }

  // Fourth pass to correct all the cell indices for the parent of each VoF.
  // To do this it is necessary to make a coarsened version of the current
  // (fine) level to place the data from the next coarser level (which may
  // have an incompatible DisjointBoxLayout).
  pout() << "    EBIndexSpace::connectedComponents-loop8" << endl;
  for (int ilev = 0; ilev < m_nlevels-1; ilev++)
  {
    CH_TIME("EBIndexSpace::connectedComponents-loop8");
    pout() << "      Level " << ilev << endl;

    // Get the new numberings for the parents of the current (fine) level
    LevelData<EBCellFAB>& newCoarNumbering = *newNumberings[ilev+1];

    // Coarsen the fine DisjointBoxLayout, create a temporary data holder
    // based on that, and copy the new numberings of the parents into this new
    // data holder
    const DisjointBoxLayout& fineGrids = newNumberings[ilev]->disjointBoxLayout();

    DisjointBoxLayout coarGrids;
    coarsen(coarGrids,fineGrids,nRef);

    EBISLayout coarEBISLayout;
    m_ebisLevel[ilev+1]->fillEBISLayout(coarEBISLayout,coarGrids,nGhosts);

    EBCellFactory coarFactory(coarEBISLayout);

    LevelData<EBCellFAB> newCoarsenedFineNumbering(coarGrids,nComps,ghostCells,coarFactory);

    newCoarNumbering.copyTo(newCoarsenedFineNumbering);

    // Go through all the connected components and fix the references to the
    // parents in their graphs
    for (int iEBIS = minEBIS; iEBIS <= maxEBIS; iEBIS++)
    {
      EBIndexSpace& curEBIS = *(connectedEBIS[iEBIS]);
      EBISLevel* fineEBISLevel = curEBIS.m_ebisLevel[ilev];

      // Iterate through all the boxes on this level making the corrections
      for (DataIterator dit = fineEBISLevel->m_grids.dataIterator(); dit.ok(); ++dit)
      {
        const Box& fineBox = fineEBISLevel->m_grids.get(dit());
        EBCellFAB& newNumberingFAB = newCoarsenedFineNumbering[dit()];
        EBGraphImplem& fineEBGraph = *(fineEBISLevel->m_graph[dit()].m_implem);

        // If there are irregular cells, there may be something to correct
        if (fineEBGraph.m_tag == EBGraphImplem::HasIrregular)
        {
          // Go through all the IntVect's in the current box
          IntVectSet fineIVS(fineBox);
          for (IVSIterator ivsit(fineIVS); ivsit.ok(); ++ivsit)
          {
            const IntVect& fineIV = ivsit();
            IntVect coarIV = coarsen(fineIV,nRef);

            GraphNode& fineGraphNode = fineEBGraph.m_graph(fineIV);

            // If there is a valid cellList then there will be references to
            // a parent VoF that need to be corrected
            if (fineGraphNode.hasValidCellList())
            {
              Vector<GraphNodeImplem>* fineCellList = fineGraphNode.m_cellList;

              // Correct the parent references in all the VoF's
              for (int iVoF = 0; iVoF < fineGraphNode.size(); iVoF++)
              {
                GraphNodeImplem& fineGraphNodeImplem = (*fineCellList)[iVoF];

                VolIndex coarVoF(coarIV,fineGraphNodeImplem.m_coarserNode);

                fineGraphNodeImplem.m_coarserNode = newNumberingFAB(coarVoF,comp);
              }
            }
          }
        }
      }
    }
  }

  // Fifth pass to correct all the VolIndex's giving the children of each VoF.
  // To do this it is necessary to make a refined version of the current
  // (coarse) level to place the data from the next finer level (which may
  // have an incompatible DisjointBoxLayout).
  pout() << "    EBIndexSpace::connectedComponents-loop9" << endl;
  for (int ilev = 1; ilev < m_nlevels; ilev++)
  {
    CH_TIME("EBIndexSpace::connectedComponents-loop9");
    pout() << "      Level " << ilev << endl;

    // Get the new numberings for the children of the current (coarse) level
    LevelData<EBCellFAB>& newFineNumbering = *newNumberings[ilev-1];

    // Refine the coarse DisjointBoxLayout, create a temporary data holder
    // based on that, and copy the new numberings of the children into this
    // new data holder
    const DisjointBoxLayout& coarGrids = newNumberings[ilev]->disjointBoxLayout();

    DisjointBoxLayout fineGrids;
    refine(fineGrids,coarGrids,nRef);

    EBISLayout fineEBISLayout;
    m_ebisLevel[ilev-1]->fillEBISLayout(fineEBISLayout,fineGrids,nGhosts);

    EBCellFactory fineFactory(fineEBISLayout);

    LevelData<EBCellFAB> newRefinedCoarNumbering(fineGrids,nComps,ghostCells,fineFactory);

    newFineNumbering.copyTo(newRefinedCoarNumbering);

    // Go through all the connected components and fix the references to the
    // children in their graphs
    for (int iEBIS = minEBIS; iEBIS <= maxEBIS; iEBIS++)
    {
      EBIndexSpace& curEBIS = *(connectedEBIS[iEBIS]);
      EBISLevel* coarEBISLevel = curEBIS.m_ebisLevel[ilev];

      // Iterate through all the boxes on this level making the corrections
      for (DataIterator dit = coarEBISLevel->m_grids.dataIterator(); dit.ok(); ++dit)
      {
        const Box& coarBox = coarEBISLevel->m_grids.get(dit());
        EBCellFAB& newNumberingFAB = newRefinedCoarNumbering[dit()];
        EBGraphImplem& coarEBGraph = *(coarEBISLevel->m_graph[dit()].m_implem);

        // If there are irregular cells, there may be something to correct
        if (coarEBGraph.m_tag == EBGraphImplem::HasIrregular)
        {
          // Go through all the IntVect's in the current box
          IntVectSet coarIVS(coarBox);
          for (IVSIterator ivsit(coarIVS); ivsit.ok(); ++ivsit)
          {
            const IntVect& coarIV = ivsit();

            GraphNode& coarGraphNode = coarEBGraph.m_graph(coarIV);

            // If there is a valid cellList then there will be references to
            // children VoF's that need to be corrected
            if (coarGraphNode.hasValidCellList())
            {
              Vector<GraphNodeImplem>* coarCellList = coarGraphNode.m_cellList;

              // Correct the children references in all the VoF's
              for (int iVoF = 0; iVoF < coarGraphNode.size(); iVoF++)
              {
                GraphNodeImplem& coarGraphNodeImplem = (*coarCellList)[iVoF];


                // Go through all the children making corrections
                int numChildren = coarGraphNodeImplem.m_finerNodes.size();
                for (int iChild = 0; iChild < numChildren; iChild++)
                {
                  VolIndex& curVolIndex = coarGraphNodeImplem.m_finerNodes[iChild];
                  const IntVect& fineIV = curVolIndex.gridIndex();

                  curVolIndex.define(fineIV,newNumberingFAB(curVolIndex,comp));
                }
              }
            }
          }
        }
      }
    }
  }

  // Sixth pass to copy all the relevant EBData into the new EBIndexSpace's.
  // Go from the coarsest level to the finest level.
  pout() << "    EBIndexSpace::connectedComponents-loop10" << endl;
  for (int ilev = m_nlevels-1; ilev >= 0; ilev--)
  {
    CH_TIME("EBIndexSpace::connectedComponents-loop10");
    pout() << "      Level " << ilev << endl;

    // For each connected component, initialize the volume data and face data
    // holders.  This has to be done using a graph with one layer of ghost
    // cells which is created there temporarily.
    for (int iEBIS = minEBIS; iEBIS <= maxEBIS; iEBIS++)
    {
      EBIndexSpace& copyEBIS   = *(connectedEBIS[iEBIS]);
      EBISLevel* copyEBISLevel = copyEBIS.m_ebisLevel[ilev];

      EBGraphFactory graphFactory(copyEBISLevel->m_domain);
      LevelData<EBGraph> ghostGraph(copyEBISLevel->m_grids,
                                    1,
                                    IntVect::Unit,
                                    graphFactory);
      Interval interval(0,0);
      copyEBISLevel->m_graph.copyTo(interval,ghostGraph,interval);

      EBDataFactory dataFactory;
      copyEBISLevel->m_data.define(copyEBISLevel->m_grids,
                                   1,
                                   IntVect::Zero,
                                   dataFactory);

      for (DataIterator dit = copyEBISLevel->m_grids.dataIterator(); dit.ok(); ++dit)
      {
//        copyEBISLevel->m_data[dit()].defineVoFData(ghostGraph[dit()],
//                                                   copyEBISLevel->m_grids.get(dit()));
//        copyEBISLevel->m_data[dit()].defineFaceData(ghostGraph[dit()],
//                                                    copyEBISLevel->m_grids.get(dit()));

        copyEBISLevel->m_data[dit()].define(ghostGraph[dit()],
                                            copyEBISLevel->m_grids.get(dit()),
                                            copyEBISLevel->m_dx, 
                                            copyEBISLevel->m_hasMoments);
      }
    }

    // Get the original EBISLevel to get access the original EBData
    EBISLevel* origEBISLevel = m_ebisLevel[ilev];

    // Get the component numbers and new VoF numbering so that the original
    // VolIndex's can be translated into the new VolIndex's for the correct
    // connected component.
    const LevelData<EBCellFAB>& numberedComponents = *numberedComponentses[ilev];
    const LevelData<EBCellFAB>& newNumbering       = *newNumberings       [ilev];

    // Iterate through all the boxes in this level copying the EBData to the
    // appropriate connected component
    for (DataIterator dit = origEBISLevel->m_grids.dataIterator(); dit.ok(); ++dit)
    {
      RefCountedPtr<EBDataImplem> origEBData = origEBISLevel->m_data[dit()].m_implem;

      const EBCellFAB& numberedComponentsFAB = numberedComponents[dit()];
      const EBCellFAB& newNumberingFAB       = newNumbering      [dit()];

      // Start with the volume data
      const BaseIVFAB<VolData>& origVolData = origEBData->m_volData;

      // Go through all the volume data a VoF at a time
      for (VoFIterator vofit(origVolData.getIVS(),origVolData.getEBGraph()); vofit.ok(); ++vofit)
      {
        // This is the original VolIndex
        VolIndex origVoF = vofit();

        int iEBIS = numberedComponentsFAB(origVoF,comp);

        // Only copy the volume data is the connected component number is
        // valid (if it is invalid then this connected component was discarded
        // because it had a total volume fracion of 0.0)
        if (iEBIS >= minEBIS && iEBIS <= maxEBIS)
        {
          // Construct the new VolIndex for the correct connected component
          int newCellIndex = newNumberingFAB(origVoF,comp);
          VolIndex newVoF(origVoF.gridIndex(),newCellIndex);

          EBIndexSpace&               copyEBIS      = *(connectedEBIS[iEBIS]);
          EBISLevel*                  copyEBISLevel = copyEBIS.m_ebisLevel[ilev];
          RefCountedPtr<EBDataImplem> copyEBData    = copyEBISLevel->m_data[dit()].m_implem;
          BaseIVFAB<VolData>&         copyVolData   = copyEBData->m_volData;

          // Copy the volume data to the correct place
          copyVolData(newVoF,comp) = origVolData(origVoF,comp);
        }
      }

      // Go through all the face data a direction and face at a time
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        const BaseIFFAB<FaceData>& origFaceData = (origEBData->m_faceData)[idir];
        for (FaceIterator faceit(origFaceData.getIVS(),origFaceData.getEBGraph(),idir,FaceStop::SurroundingWithBoundary); faceit.ok(); ++faceit)
        {
          // This is the original FaceIndex and the VolIndex's that define it
          FaceIndex origFace = faceit();
          const VolIndex& origLoVoF = origFace.getVoF(Side::Lo);
          const VolIndex& origHiVoF = origFace.getVoF(Side::Hi);

          // Start out with invalid new values and fill them in from the lo or
          // hi VoF depending on which is valid
          int iEBIS = -1;
          int newLoCellIndex = -1;
          int newHiCellIndex = -1;

          // The lo VoF is valid
          if (origLoVoF.cellIndex() != -1)
          {
            if (iEBIS == -1)
            {
              iEBIS = numberedComponentsFAB(origLoVoF,comp);
            }
            newLoCellIndex = newNumberingFAB(origLoVoF,comp);
          }

          // The hi VoF is valid
          if (origHiVoF.cellIndex() != -1)
          {
            if (iEBIS == -1)
            {
              iEBIS = numberedComponentsFAB(origHiVoF,comp);
            }
            newHiCellIndex = newNumberingFAB(origHiVoF,comp);
          }

          // Only copy the face data is the connected component number is
          // valid (if it is invalid then this connected component was
          // discarded because it had a total volume fracion of 0.0)
          if (iEBIS >= minEBIS && iEBIS <= maxEBIS)
          {
            // Construct the new FaceIndex for the correct connected component
            // (which includes new VolIndex's)
            VolIndex newLoVoF(origLoVoF.gridIndex(),newLoCellIndex);
            VolIndex newHiVoF(origHiVoF.gridIndex(),newHiCellIndex);

            FaceIndex newFace(newLoVoF,newHiVoF,origFace.direction());

            EBIndexSpace&               copyEBIS      = *(connectedEBIS[iEBIS]);
            EBISLevel*                  copyEBISLevel = copyEBIS.m_ebisLevel[ilev];
            RefCountedPtr<EBDataImplem> copyEBData    = copyEBISLevel->m_data[dit()].m_implem;
            BaseIFFAB<FaceData>&        copyFaceData  = (copyEBData->m_faceData)[idir];

            // Copy the face data to the correct place
            copyFaceData(newFace,comp) = origFaceData(origFace,comp);
          }
        }
      }
    }
  }

  // Clean up the temporary data holders defined using the original
  // EBIndexSpace
  pout() << "    EBIndexSpace::connectedComponents-loop11" << endl;
  for (int ilev = 0; ilev < m_nlevels; ilev++)
  {
    CH_TIME("EBIndexSpace::connectedComponents-loop11");
    pout() << "      Level " << ilev << endl;

    delete numberedComponentses[ilev];
    delete newNumberings[ilev];
  }

  if (anyRegularWithMultivaluedParent)
  {
    // TODO: Remove storage for any regular cells which had a multivalued
    // parent in the original graph but don't have a multivalued parent in
    // the connected componenet graphs.  The graphs are all correct without
    // this step but contain extra storage to explicitly store the cellIndex
    // of the parent and result in an extra level of indirection accessing
    // the parent.  On the other hand, these don't occur frequently.
  }

  // Return only the biggest connected component
  if (a_onlyBiggest)
  {
    RefCountedPtr<EBIndexSpace> biggestComponent(NULL);

    // If the EBIS is empty, then biggestEBIS = -1 and connectedEBIS.size() = 0.
    if (biggestEBIS >= 0)
      {
        biggestComponent = connectedEBIS[biggestEBIS];
      }

    connectedEBIS.resize(1);
    connectedEBIS[0] = biggestComponent;
  }

  // Return the Vector of RefCountedPtr<EBIndexSpace> - one for every
  // connected component
  return connectedEBIS;
}

#define CH_LOCAL_CLASS_DATA   4
#define CH_LOCAL_PROC_ID      0
#define CH_LOCAL_BOX_ID       1
#define CH_LOCAL_COMP_NUM     2
#define CH_LOCAL_VOLUME       3

#define CH_GLOBAL_CLASS_DATA  2
#define CH_GLOBAL_COMP_NUM    0
#define CH_GLOBAL_VOLUME      1

#define CH_MAP_KEY_VALUE_DATA (CH_LOCAL_CLASS_DATA - 1 + CH_GLOBAL_CLASS_DATA)
#define CH_MAP_KEY_PROC_ID    (CH_LOCAL_PROC_ID                              )
#define CH_MAP_KEY_BOX_ID     (CH_LOCAL_BOX_ID                               )
#define CH_MAP_KEY_COMP_NUM   (CH_LOCAL_COMP_NUM                             )
#define CH_MAP_VALUE_COMP_NUM (CH_GLOBAL_COMP_NUM + (CH_LOCAL_CLASS_DATA - 1))
#define CH_MAP_VALUE_VOLUME   (CH_GLOBAL_VOLUME   + (CH_LOCAL_CLASS_DATA - 1))

// Divide the EBIndexSpace into connected components and return a Vector of
// disjoint EBIndexSpace's corresponding to each connected component.
Vector<RefCountedPtr<EBIndexSpace> > EBIndexSpace::findConnectedComponentsNew(int        & a_numComponents,
                                                                           const bool & a_onlyBiggest)
{
  CH_TIME("EBIndexSpace::findConnectedComponentsNew");
  pout() << "    EBIndexSpace::findConnectedComponentsNew" << endl;

  // Begin by setting up some data holders and other infrastructure to number
  // the connected components and remember the renumbering of the VoFs.
  Vector<LevelData<EBCellFAB>* > numberedComponentses(m_nlevels,NULL);
  Vector<LevelData<EBCellFAB>* > newNumberings(m_nlevels,NULL);

  // This is only needed at the coarsest level to find connected components
  // over a set of boxes
  LevelData<EBCellFAB> coarComponentInfo;

  // EBISLayout's for each level
  Vector<EBISLayout> ebisLayouts(m_nlevels);

  // All our temporary data structures have one component
  int nComps = 1;

  // All our temporary data structures have one layer of ghostcells
  int nGhosts = 1;
  IntVect ghostCells = nGhosts * IntVect::Unit;

  // Fill all the EBISLayouts and allocate/define all the local data holders
  pout() << "    EBIndexSpace::connectedComponentsNew-loop1" << endl;
  for (int ilev = 0; ilev < m_nlevels; ilev++)
  {
    CH_TIME("EBIndexSpace::connectedComponentsNew-loop1");
    pout() << "      Level " << ilev << endl;

    EBISLevel* curEBISLevel = m_ebisLevel[ilev];

    DisjointBoxLayout curGrids = curEBISLevel->m_grids;

    EBISLayout& curEBISLayout = ebisLayouts[ilev];
    curEBISLevel->fillEBISLayout(curEBISLayout,curGrids,nGhosts);

    EBCellFactory curFactory(curEBISLayout);

    LevelData<EBCellFAB>* numberedComponents = new LevelData<EBCellFAB>();
    LevelData<EBCellFAB>* newNumbering       = new LevelData<EBCellFAB>();

    numberedComponents->define(curGrids,nComps,ghostCells,curFactory);
    newNumbering      ->define(curGrids,nComps,ghostCells,curFactory);

    numberedComponentses[ilev] = numberedComponents;
    newNumberings       [ilev] = newNumbering;

    if (ilev == m_nlevels-1)
    {
      coarComponentInfo.define(curGrids,CH_LOCAL_CLASS_DATA,ghostCells,curFactory);
    }
  }

  // Get the coarsest level and number the connected components at that level.
  // The current implementation assumes there is only one box at the coarsest
  // level - this is checked below.
  EBISLevel* coarEBISLevel = m_ebisLevel[m_nlevels-1];

  DisjointBoxLayout coarGrids = coarEBISLevel->m_grids;

  // This really only needs to be an "int" and not a "Real" but all our
  // infrastructure (e.g., I/O, debugging, library functions) is set up to
  // work with "EBCellFAB".  This can be changed if needed to be an "int" data
  // holder.
  LevelData<EBCellFAB>& coarNumberedComponents = *numberedComponentses[m_nlevels-1];

  // Initialize all the component numbers and box number to -1 (invalid)
  EBLevelDataOps::setVal(coarNumberedComponents,-1.0);
  EBLevelDataOps::setVal(coarComponentInfo,-1.0);

  // Currently this code will only work if there is only one box at the
  // coarsest level of the geometry.  This will need to be fixed by merging
  // different numberings of the same component in different boxes.  This can
  // be done via ghost cells, exchanges, the building of equivalence classes,
  // and, finally, a renumbering operation.
  Real coarDx = coarEBISLevel->m_dx;

#if 0
  unsigned int numGrids = coarGrids.size();

  pout() << "    dx: " << coarDx << endl;
  pout() << endl;

  pout() << "    Number of coarsest grids: " << numGrids << endl;
  pout() << "    Grids: " << endl;
  for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
  {
    pout() << "      " << coarGrids[dit()] << endl;
  }
  pout() << endl;
#endif

  // For counting the new EBIndexspace's and remember the number assigned to
  // each
  int numEBIS = 0;
  Vector<int> valueEBIS;

  int biggestEBIS = -1;

  Vector<IndexTM<Real,CH_LOCAL_CLASS_DATA> > compClasses;
  Vector<int>                                boxClassIndex;

  int myProcID = procID();

  int numCompsLocal = 0;

  // Go through all boxes.
  pout() << "    EBIndexSpace::connectedComponentsNew-loop2" << endl;
  for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
  {
    CH_TIME("EBIndexSpace::connectedComponentsNew-loop2");

    int curBoxClassIndex = compClasses.size();
    boxClassIndex.push_back(curBoxClassIndex);

    // Get a bunch of local references
    int boxIndex = coarGrids.index(dit());

    const Box& curBox = coarGrids.get(dit());
    const IntVectSet curIVS(curBox);
    pout() << "      Box: " << curBox << endl;
#if 0
    pout() << "    Box index: " << boxIndex << ", box: " << curBox << endl;
#endif
    EBCellFAB&     curEBCellFABInfo = coarComponentInfo[dit()];
    const EBGraph& curEBGraph       = coarEBISLevel->m_graph[dit()];
    const EBISBox& curEBISBox       = curEBCellFABInfo.getEBISBox();

    curEBCellFABInfo.setVal(0,myProcID);
    curEBCellFABInfo.setVal(1,boxIndex);

    unsigned int curNum = 0;

    // Iterate through the box, VoF by VoF and number anything that isn't
    // already numbered.
    for (VoFIterator vofit(curIVS,curEBGraph); vofit.ok(); ++vofit)
    {
      // Current VoF
      const VolIndex& vof = vofit();

      Real totalVolFrac;

      // Examine the connected component starting at "vof".  If it hasn't been
      // numbered yet, return true otherwise false.  Also, return the
      // volume fraction associated with the connect component.
      if (setAllConnectedVoFsNew(totalVolFrac,curEBCellFABInfo,curBox,curEBGraph,curEBISBox,vof,vof,curNum))
      {
        Real totalVol = totalVolFrac*pow(coarDx,SpaceDim);
#if 0
        pout() << "        totalVol: " << totalVol << endl;
#endif
        if (totalVolFrac > 0.0)
        {
          // Remember the component number and increment the number of
          // components
          valueEBIS.push_back(curNum);
          numEBIS++;

          IndexTM<Real,CH_LOCAL_CLASS_DATA> cur(myProcID,boxIndex,curNum,totalVol);
          compClasses.push_back(cur);

          curNum++;
          numCompsLocal++;
        }
        else
        {
          // If we're not counting the current component (because it has a total
          // volume fraction of 0.0, set the component number back to -1

          resetAllConnectedVoFsNew(curEBCellFABInfo,curBox,curEBGraph,curEBISBox,vof,vof);
        }
      }
    }

    for (VoFIterator vofit(curIVS,curEBGraph); vofit.ok(); ++vofit)
    {
      const VolIndex& curVoF = vofit();
      const Real comp = curEBCellFABInfo(curVoF,2);

      if (comp >= 0)
      {
        curEBCellFABInfo(curVoF,CH_LOCAL_VOLUME) = compClasses[curBoxClassIndex+comp][CH_LOCAL_VOLUME];
      }
    }

    pout() << "      connected components: " << curNum << endl;
  }

#if 0
  pout() << endl;

  pout() << "    Box class index:" << endl;
  for (int iclass = 0; iclass < boxClassIndex.size(); iclass++)
  {
    pout() << "      " << iclass << ": " << boxClassIndex[iclass] << endl;
  }
  pout() << endl;
#endif

  std::set<std::set<IndexTM<Real,CH_LOCAL_CLASS_DATA> > > subClasses;

#if 0
  pout() << "    Classes:" << endl;
#endif
  Real volume = 0.0;
  for (int iclass = 0; iclass < compClasses.size(); iclass++)
  {
    std::set<IndexTM<Real,CH_LOCAL_CLASS_DATA> > singlet;
    singlet.insert(compClasses[iclass]);

    subClasses.insert(singlet);
#if 0
    pout() << "      " << compClasses[iclass] << endl;
#endif
    volume += compClasses[iclass][CH_LOCAL_VOLUME];
  }
#if 0
  pout() << "    Volume: " << volume << endl;
  pout() << endl;
#endif

  // Make sure all ghostcells are correct
  coarNumberedComponents.exchange();
  coarComponentInfo.exchange();

  pout() << "    EBIndexSpace::connectedComponentsNew-loop2.5" << endl;
  for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
  {
    CH_TIME("EBIndexSpace::connectedComponentsNew-loop2.5");

    // Get a bunch of local references
    const Box& curBox = coarGrids.get(dit());
    const IntVectSet curIVS(curBox);
    pout() << "      Box: " << curBox << endl;
#if 0
    int boxIndex = coarGrids.index(dit());

    pout() << "    Box index: " << boxIndex << ", box: " << curBox << endl;
#endif
    EBCellFAB&     curEBCellFABInfo = coarComponentInfo[dit()];
    const EBGraph& curEBGraph       = coarEBISLevel->m_graph[dit()];
    const EBISBox& curEBISBox       = curEBCellFABInfo.getEBISBox();

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
      {
        const Side::LoHiSide& side = sit();

        Box curSide = adjCellBox(curBox,idir,side,1);
        curSide.shift(idir,-sign(side));

        IntVectSet ivs(curSide);

        for (VoFIterator vofit(ivs,curEBGraph); vofit.ok(); ++vofit)
        {
          const VolIndex& curVoF = vofit();
          const Vector<FaceIndex>& faces = curEBISBox.getFaces(curVoF,idir,side);
          IndexTM<Real,CH_LOCAL_CLASS_DATA> curID(curEBCellFABInfo(curVoF,CH_LOCAL_PROC_ID ),
                                                  curEBCellFABInfo(curVoF,CH_LOCAL_BOX_ID  ),
                                                  curEBCellFABInfo(curVoF,CH_LOCAL_COMP_NUM),
                                                  curEBCellFABInfo(curVoF,CH_LOCAL_VOLUME  ));

          int numFaces = faces.size();

          for (int iface = 0; iface < numFaces; iface++)
          {
            const FaceIndex& face = faces[iface];
            if (!face.isBoundary())
            {
              const VolIndex& neighborVoF = face.getVoF(side);

              IndexTM<Real,CH_LOCAL_CLASS_DATA> neighborID(curEBCellFABInfo(neighborVoF,CH_LOCAL_PROC_ID ),
                                                           curEBCellFABInfo(neighborVoF,CH_LOCAL_BOX_ID  ),
                                                           curEBCellFABInfo(neighborVoF,CH_LOCAL_COMP_NUM),
                                                           curEBCellFABInfo(neighborVoF,CH_LOCAL_VOLUME  ));

              std::set<IndexTM<Real,CH_LOCAL_CLASS_DATA> > curClass;
              curClass.insert(curID);
              curClass.insert(neighborID);

              subClasses.insert(curClass);
            }
          }
        }
      }
    }
  }

  std::list<std::set<IndexTM<Real,CH_LOCAL_CLASS_DATA> > > curClasses;

#if CH_MPI
  {
    int bufferLength = subClasses.size();

    int procs = 0;
    int rank  = 0;
    int root  = 0;

    MPI_Comm_size(Chombo_MPI::comm,&procs);
    MPI_Comm_rank(Chombo_MPI::comm,&rank );

    Vector<int> rBufferSizes1(1);
    Vector<int> rBufferOffsets1(1);
    if (rank == root)
    {
      rBufferSizes1  .resize(procs);
      rBufferOffsets1.resize(procs);
    }

    MPI_Gather(&bufferLength,1,MPI_INT,&rBufferSizes1[0],1,MPI_INT,root,Chombo_MPI::comm);

    Vector<int> rBufferSetSizes(1);
    if (rank == root)
    {
      int offset = 0;
      for (int i = 0; i < procs; i++)
      {
        rBufferOffsets1[i] = offset;
        offset += rBufferSizes1[i];
#if 0
        pout() << "proc "       << i
               << ": offset = " << rBufferOffsets1[i]
               << ", size = "   << rBufferSizes1  [i]
               << endl;
#endif
      }
      rBufferSetSizes.resize(offset);
    }

#if 0
    pout() << "bufferLength: " << bufferLength << endl;
#endif

    Vector<int> sBufferSetSizes(bufferLength);
    if (bufferLength == 0)
    {
      sBufferSetSizes.resize(1);
      sBufferSetSizes[0] = 0;
    }

    int sBufferSize = 0;

    int curCount = 0;
    for (std::set<std::set<IndexTM<Real,CH_LOCAL_CLASS_DATA> > >::iterator it = subClasses.begin();
         it != subClasses.end();
         ++it)
    {
      sBufferSetSizes[curCount]  = it->size();
      sBufferSize               += it->size();
#if 0
      pout() << "  set " << curCount << " size = " << sBufferSetSizes[curCount] << endl;
#endif
      curCount++;
    }

    MPI_Gatherv(&sBufferSetSizes[0],bufferLength,MPI_INT,
                &rBufferSetSizes[0],&rBufferSizes1[0],&rBufferOffsets1[0],
                MPI_INT,root,Chombo_MPI::comm);

    Vector<int> rBufferSizes2(1);
    Vector<int> rBufferOffsets2(1);
    int rSize = 0;

    if (rank == root)
    {
      rBufferSizes2  .resize(procs);
      rBufferOffsets2.resize(procs);

      int totalOffset = 0;
      for (int i = 0; i < procs; i++)
      {
        int offset = rBufferOffsets1[i];
        int totalSize = 0;
        for (int j = 0; j < rBufferSizes1[i]; j++)
        {
          totalSize += rBufferSetSizes[offset+j];
        }

        rBufferSizes2  [i] = totalSize  *CH_LOCAL_CLASS_DATA;
        rBufferOffsets2[i] = totalOffset*CH_LOCAL_CLASS_DATA;
#if 0
        pout() << "proc " << i << ": " << totalSize << ", " << totalOffset << endl;
#endif
        totalOffset += totalSize;
      }

      rSize = totalOffset*CH_LOCAL_CLASS_DATA;
#if 0
      pout() << "rSize = " << rSize << endl;
#endif
    }

    Vector<Real> sBuffer(sBufferSize*CH_LOCAL_CLASS_DATA);
    if (sBufferSize == 0)
    {
      sBuffer.resize(1);
    }
    int sOffset = 0;
#if 0
    pout() << "sBufferSize      = " << sBufferSize << endl;
    pout() << "sBufferSizeReals = " << sBufferSize*CH_LOCAL_CLASS_DATA << endl;
#endif
    for (std::set<std::set<IndexTM<Real,CH_LOCAL_CLASS_DATA> > >::iterator it = subClasses.begin();
         it != subClasses.end();
         ++it)
    {
      for (std::set<IndexTM<Real,CH_LOCAL_CLASS_DATA> >::iterator it2 = it->begin();
           it2 != it->end();
           ++it2)
      {
        sBuffer[sOffset+CH_LOCAL_PROC_ID ] = (*it2)[CH_LOCAL_PROC_ID ];
        sBuffer[sOffset+CH_LOCAL_BOX_ID  ] = (*it2)[CH_LOCAL_BOX_ID  ];
        sBuffer[sOffset+CH_LOCAL_COMP_NUM] = (*it2)[CH_LOCAL_COMP_NUM];
        sBuffer[sOffset+CH_LOCAL_VOLUME  ] = (*it2)[CH_LOCAL_VOLUME  ];
#if 0
        pout() << "sOffset " << sOffset << ":";
        pout() <<  " " << (*it2)[CH_LOCAL_PROC_ID];
        pout() << ", " << (*it2)[CH_LOCAL_BOX_ID];
        pout() << ", " << (*it2)[CH_LOCAL_COMP_NUM];
        pout() << ", " << (*it2)[CH_LOCAL_VOLUME];
        pout() << endl;
#endif
        sOffset += CH_LOCAL_CLASS_DATA;
      }
    }

    Vector<Real> rBuffer(1);
    if (rank == root)
    {
      rBuffer.resize(rSize);
    }

    MPI_Gatherv(&sBuffer[0],sBufferSize*CH_LOCAL_CLASS_DATA,MPI_CH_REAL,
                &rBuffer[0],&rBufferSizes2[0],&rBufferOffsets2[0],
                MPI_CH_REAL,root,Chombo_MPI::comm);

    if (rank == root)
    {
      int rOffset = 0;
      for (int i = 0; i < rBufferSetSizes.size(); i++)
      {
        std::set<IndexTM<Real,CH_LOCAL_CLASS_DATA> > curSet;

        for (int j = 0; j < rBufferSetSizes[i]; j++)
        {
          IndexTM<Real,CH_LOCAL_CLASS_DATA> localComp;

          localComp[CH_LOCAL_PROC_ID ] = rBuffer[rOffset + CH_LOCAL_PROC_ID ];
          localComp[CH_LOCAL_BOX_ID  ] = rBuffer[rOffset + CH_LOCAL_BOX_ID  ];
          localComp[CH_LOCAL_COMP_NUM] = rBuffer[rOffset + CH_LOCAL_COMP_NUM];
          localComp[CH_LOCAL_VOLUME  ] = rBuffer[rOffset + CH_LOCAL_VOLUME  ];

          curSet.insert(localComp);

          rOffset += CH_LOCAL_CLASS_DATA;
        }

        curClasses.push_back(curSet);
      }
    }
  }
#else
  for (std::set<std::set<IndexTM<Real,CH_LOCAL_CLASS_DATA> > >::iterator it = subClasses.begin();
       it != subClasses.end();
       ++it)
  {
    curClasses.push_back(*it);
  }
#endif

  std::list<std::set<IndexTM<Real,CH_LOCAL_CLASS_DATA> > > finalClasses;

  std::list<std::set<IndexTM<Real,CH_LOCAL_CLASS_DATA> > >::iterator it = curClasses.begin();
  while (it != curClasses.end())
  {
#if 0
    pout() << "    Lists size: " << curClasses  .size()
           <<               ", " << finalClasses.size() << endl;
#endif

    std::set<IndexTM<Real,CH_LOCAL_CLASS_DATA> > curSet = *it;

    curClasses.erase(it);

    bool curSetChanged = false;

    std::list<std::set<IndexTM<Real,CH_LOCAL_CLASS_DATA> > >::iterator it2 = curClasses.begin();
    while (it2 != curClasses.end())
    {
      const std::set<IndexTM<Real,CH_LOCAL_CLASS_DATA> >& scanSet = *it2;

      std::set<IndexTM<Real,CH_LOCAL_CLASS_DATA> > intersection;

      std::set_intersection
             (curSet.begin(),curSet.end(),
              scanSet.begin(),scanSet.end(),
              std::inserter(intersection, intersection.end()));

      if (!intersection.empty())
      {
        std::set_union(curSet.begin(),curSet.end(),
                       scanSet.begin(),scanSet.end(),
                       std::inserter(curSet, curSet.end()));

        curSetChanged = true;

        std::list<std::set<IndexTM<Real,CH_LOCAL_CLASS_DATA> > >::iterator it2Next = it2;
        ++it2Next;
        
        curClasses.erase(it2);

        it2 = it2Next;
      }
      else
      {
        ++it2;
      }
    }

    if (curSetChanged)
    {
      curClasses.push_back(curSet);
    }
    else
    {
      finalClasses.push_back(curSet);
    }

    it = curClasses.begin();
  }

#if 0
  pout() << "    Lists size: " << curClasses  .size()
         <<               ", " << finalClasses.size() << endl;

  pout() << endl;
#endif

  Real totalVolume = 0.0;
  int count = 0;

  Real maxVolFrac = -1.0;

  std::map<IndexTM<Real,CH_LOCAL_CLASS_DATA-1>,IndexTM<Real,CH_GLOBAL_CLASS_DATA> > componentMap;

  for (std::list<std::set<IndexTM<Real,CH_LOCAL_CLASS_DATA> > >::const_iterator it = finalClasses.begin();
       it != finalClasses.end();
       ++it, count++)
  {
    pout() << "    Class " << count << ":" << endl;

    Real volume = 0.0;
    for (std::set<IndexTM<Real,CH_LOCAL_CLASS_DATA> >::const_iterator si = it->begin();
         si != it->end();
         ++si)
    {
      const IndexTM<Real,CH_LOCAL_CLASS_DATA>& member = *si;

      pout() << "      " << member << endl;
      volume += member[CH_LOCAL_VOLUME];
    }

    for (std::set<IndexTM<Real,CH_LOCAL_CLASS_DATA> >::const_iterator si = it->begin();
         si != it->end();
         ++si)
    {
      const IndexTM<Real,CH_LOCAL_CLASS_DATA>& member = *si;
      IndexTM<Real,CH_LOCAL_CLASS_DATA-1> key(member[CH_LOCAL_PROC_ID ],
                                              member[CH_LOCAL_BOX_ID  ],
                                              member[CH_LOCAL_COMP_NUM]);

      IndexTM<Real,CH_GLOBAL_CLASS_DATA> value(count,volume);

      componentMap[key] = value;
    }

    if (volume > maxVolFrac)
    {
      maxVolFrac = volume;
      biggestEBIS = count;
    }

    totalVolume += volume;

    pout() << "    Volume: " << volume << endl;
    pout() << endl;
  }

#ifdef CH_MPI
  {
    int procs = 0;
    int rank  = 0;
    int root  = 0;

    MPI_Comm_size(Chombo_MPI::comm,&procs);
    MPI_Comm_rank(Chombo_MPI::comm,&rank );

    MPI_Bcast(&maxVolFrac,1,MPI_CH_REAL,root,Chombo_MPI::comm);
    MPI_Bcast(&biggestEBIS,1,MPI_INT,root,Chombo_MPI::comm);
  }
#endif

  pout() << "    Total volume: " << totalVolume << endl;
  pout() << endl;

  pout() << "    Largest volume: " << maxVolFrac  << endl;
  pout() << "    EBIS #        : " << biggestEBIS << endl;
  pout() << endl;

#ifdef CH_MPI
  int procs = 0;
  int rank  = 0;
  int root  = 0;

  MPI_Comm_size(Chombo_MPI::comm,&procs);
  MPI_Comm_rank(Chombo_MPI::comm,&rank );

  int mapSize = componentMap.size();

  Vector<Real> sMapPacked(1);
  Vector<int>  sNumCompsGlobal(1);
  Vector<int>  sNumCompsOffset(1);

  if (rank == root)
  {
    sMapPacked.resize(mapSize * CH_MAP_KEY_VALUE_DATA);
    sNumCompsGlobal.resize(procs);
    sNumCompsOffset.resize(procs);

    int lastProcID   = -1;
    int curNumComps  =  0;
    int curNumOffset =  0;
    int procIndex    = -1;

    int mapIndex = 0;
    for (std::map<IndexTM<Real,CH_LOCAL_CLASS_DATA-1>,IndexTM<Real,CH_GLOBAL_CLASS_DATA> >::const_iterator mi = componentMap.begin();
         mi != componentMap.end();
         ++mi
        )
    {
      int curProcID = mi->first[CH_LOCAL_PROC_ID ];
      int curBoxID  = mi->first[CH_LOCAL_BOX_ID  ];

      if (curProcID != lastProcID)
      {
        if (procIndex != -1)
        {
          sNumCompsGlobal[procIndex] = CH_MAP_KEY_VALUE_DATA * curNumComps;
          sNumCompsOffset[procIndex] = CH_MAP_KEY_VALUE_DATA * curNumOffset;

          curNumOffset += curNumComps;
        }

        curNumComps =  0;

        for (procIndex++; procIndex < curProcID; procIndex++)
        {
          sNumCompsGlobal[procIndex] = CH_MAP_KEY_VALUE_DATA * curNumComps;
          sNumCompsOffset[procIndex] = CH_MAP_KEY_VALUE_DATA * curNumOffset;
        }
      }

      curNumComps++;

      sMapPacked[mapIndex + CH_MAP_KEY_PROC_ID ] = curProcID;
      sMapPacked[mapIndex + CH_MAP_KEY_BOX_ID  ] = curBoxID;
      sMapPacked[mapIndex + CH_MAP_KEY_COMP_NUM] = mi->first[CH_LOCAL_COMP_NUM];

      sMapPacked[mapIndex + CH_MAP_VALUE_COMP_NUM] = mi->second[CH_GLOBAL_COMP_NUM];
      sMapPacked[mapIndex + CH_MAP_VALUE_VOLUME  ] = mi->second[CH_GLOBAL_VOLUME  ];

      mapIndex += CH_MAP_KEY_VALUE_DATA;

      lastProcID = curProcID;
    }

    sNumCompsGlobal[procIndex] = CH_MAP_KEY_VALUE_DATA * curNumComps;
    sNumCompsOffset[procIndex] = CH_MAP_KEY_VALUE_DATA * curNumOffset;

    componentMap.clear();

#if 0
    pout() << "Computed on root:" << endl;
    for (int i = 0; i < sNumCompsGlobal.size(); i++)
    {
      pout() << "  Proc " << i << ": "
             << sNumCompsGlobal[i] << " comps, "
             << sNumCompsOffset[i] << " offset"
             << endl;
    }
    pout() << endl;
#endif
  }

  int rNumCompsGlobal;

  MPI_Scatter(&sNumCompsGlobal[0],1,MPI_INT,&rNumCompsGlobal,1,MPI_INT,root,Chombo_MPI::comm);

  if (rNumCompsGlobal != CH_MAP_KEY_VALUE_DATA * numCompsLocal)
  {
    MayDay::Error("EBIndexSpace::findConnectedComponentsNew: number of local components don't agree between local and global computation");
  }

#if 0
  pout() << "Total comps locally/globally:" << endl;
  pout() << "  Proc " << rank << ": " 
         << numCompsLocal << " vs. "
         << rNumCompsGlobal << endl;
  pout() << endl;
#endif

  Vector<Real> rMapPacked(CH_MAP_KEY_VALUE_DATA * numCompsLocal);
  if (numCompsLocal == 0)
  {
    rMapPacked.resize(1);
  }

  MPI_Scatterv(&sMapPacked[0],&sNumCompsGlobal[0],&sNumCompsOffset[0],MPI_CH_REAL,
               &rMapPacked[0],CH_MAP_KEY_VALUE_DATA * numCompsLocal,MPI_CH_REAL,
               root,Chombo_MPI::comm);

  for (int i = 0; i < numCompsLocal; i++)
  {
    Real keyProcID  = rMapPacked[CH_MAP_KEY_VALUE_DATA*i + CH_MAP_KEY_PROC_ID ];
    Real keyBoxID   = rMapPacked[CH_MAP_KEY_VALUE_DATA*i + CH_MAP_KEY_BOX_ID  ];
    Real keyCompNum = rMapPacked[CH_MAP_KEY_VALUE_DATA*i + CH_MAP_KEY_COMP_NUM];

    Real valueCompNum = rMapPacked[CH_MAP_KEY_VALUE_DATA*i + CH_MAP_VALUE_COMP_NUM];
    Real valueVolume  = rMapPacked[CH_MAP_KEY_VALUE_DATA*i + CH_MAP_VALUE_VOLUME  ];

#if 0
    pout() << "Local map packed " << i
           << ": (" << keyProcID 
           << " " << keyBoxID
           << " " << keyCompNum
           << ")";

    pout() << " -> (" << valueCompNum
           << " " << valueVolume
           << ")" << endl;
#endif

    IndexTM<Real,CH_LOCAL_CLASS_DATA-1> key(keyProcID,keyBoxID,keyCompNum);
    IndexTM<Real,CH_GLOBAL_CLASS_DATA> value(valueCompNum,valueVolume);

    componentMap[key] = value;
  }

  a_numComponents = finalClasses.size();

  MPI_Bcast(&a_numComponents,1,MPI_INT,root,Chombo_MPI::comm);

#else
  a_numComponents = finalClasses.size();
#endif

  numEBIS = a_numComponents;

  pout() << "Number of components: " << a_numComponents << endl;
  pout() << endl;

  pout() << "    EBIndexSpace::connectedComponentsNew-loop2.75" << endl;
  for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit)
  {
    CH_TIME("EBIndexSpace::connectedComponentsNew-loop2.75");

    int curBoxClassIndex = compClasses.size();
    boxClassIndex.push_back(curBoxClassIndex);

    // Get a bunch of local references
    // int boxIndex = coarGrids.index(dit());

    const Box& curBox = coarGrids.get(dit());
    const IntVectSet curIVS(curBox);
    pout() << "      Box: " << curBox << endl;
#if 0
    pout() << "    Box index: " << boxIndex << ", box: " << curBox << endl;
#endif
    EBCellFAB&     curEBCellFABNum  = coarNumberedComponents[dit()];
    EBCellFAB&     curEBCellFABInfo = coarComponentInfo[dit()];
    const EBGraph& curEBGraph       = coarEBISLevel->m_graph[dit()];
    // const EBISBox& curEBISBox       = curEBCellFABInfo.getEBISBox();

    for (VoFIterator vofit(curIVS,curEBGraph); vofit.ok(); ++vofit)
    {
      const VolIndex& curVoF = vofit();

      int curProcID    = curEBCellFABInfo(curVoF,CH_LOCAL_PROC_ID );
      int curBoxIndex  = curEBCellFABInfo(curVoF,CH_LOCAL_BOX_ID  );
      int curComponent = curEBCellFABInfo(curVoF,CH_LOCAL_COMP_NUM);

      IndexTM<Real,CH_LOCAL_CLASS_DATA-1> key(curProcID,curBoxIndex,curComponent);

      const IndexTM<Real,CH_GLOBAL_CLASS_DATA>& value = componentMap[key];

      curEBCellFABNum(curVoF,0) = value[CH_GLOBAL_COMP_NUM];
    }
  }

  int minEBIS;
  int maxEBIS;

  if (a_onlyBiggest)
  {
    // If we only want the biggest connected componnent, only iterate over
    // that components index
    if (biggestEBIS >= 0)
      {
        minEBIS = biggestEBIS;
        maxEBIS = biggestEBIS;
      }
    else
      {
        // This corresponds to an empty EBIS
        minEBIS =  0;
        maxEBIS = -1;
      }
  }
  else
  {
    // Otherwise, iterate over all the connected component indices
    minEBIS = 0;
    maxEBIS = numEBIS-1;
  }

#if 0
  pout() << "numEBIS:     " << numEBIS     << endl;
  pout() << "biggestEBIS: " << biggestEBIS << endl;
  pout() << "minEBIS:     " << minEBIS     << endl;
  pout() << "maxEBIS:     " << maxEBIS     << endl;

  pout() << "Hi 2..." << endl;
#endif

  // At this point, renumber the connected components from 0 to numEBIS-1.
  // Currently, this is already the case so nothing happens here.

  // There is only one component in each temporary data holder
  int comp = 0;

  // The refinement ratio on the geometry size of things is always 2
  int nRef = 2;

  // Pass the numbering for connected components from coarse to fine levels.
  pout() << "    EBIndexSpace::connectedComponentsNew-loop3" << endl;
  for (int ilev = m_nlevels-2; ilev >= 0; ilev--)
  {
    CH_TIME("EBIndexSpace::connectedComponentsNew-loop3");
    pout() << "      Level " << ilev << endl;

    // Get a lot of local references to get things ready to go
    EBISLevel* coarEBISLevel = m_ebisLevel[ilev+1];
    EBISLevel* fineEBISLevel = m_ebisLevel[ilev];

    DisjointBoxLayout coarGrids = coarEBISLevel->m_grids;
    DisjointBoxLayout fineGrids = fineEBISLevel->m_grids;

    const ProblemDomain& coarDomain = coarEBISLevel->m_domain;

    EBISLayout& coarEBISLayout = ebisLayouts[ilev+1];
    EBISLayout& fineEBISLayout = ebisLayouts[ilev];

    LevelData<EBCellFAB>& coarNumberedComponents = *numberedComponentses[ilev+1];
    LevelData<EBCellFAB>& fineNumberedComponents = *numberedComponentses[ilev];


    // Initial all the numberings on this level to -1 (invalid)
    EBLevelDataOps::setVal(fineNumberedComponents,-1.0);

    // Set up a piecewise constant interpolator from the current coarse level
    // to the current fine level
    EBCFCopy constantInterp(fineGrids,
                            coarGrids,
                            fineEBISLayout,
                            coarEBISLayout,
                            coarDomain,
                            nRef,
                            nComps,
                            this,
                            ghostCells);

    // Do the interpolation/copy
    Interval oneComp(comp,comp);
    constantInterp.copy(fineNumberedComponents,
                        coarNumberedComponents,
                        oneComp);

    // Make sure all ghostcells are correct
    fineNumberedComponents.exchange();
  }

  // The set of EBIndexSpace's that are each connected components of current
  // EBIndexSpace.
  Vector<RefCountedPtr<EBIndexSpace> > connectedEBIS(numEBIS);

  // Make copies of the main EBIndexSpace for each connected component.
  pout() << "    EBIndexSpace::connectedComponentsNew-loop4" << endl;
  for (int iEBIS = minEBIS; iEBIS <= maxEBIS; iEBIS++)
  {
    CH_TIME("EBIndexSpace::connectedComponentsNew-loop4");
    pout() << "      EBIS #" << iEBIS << endl;

    // Make a new, empty EBIndexSpace for the current component
    RefCountedPtr<EBIndexSpace> copyEBIS(new EBIndexSpace());

    // Set the member data that is not level dependent
    copyEBIS->m_nCellMax  = m_nCellMax;
    copyEBIS->m_isDefined = m_isDefined;

    copyEBIS->m_distributedData = m_distributedData;

    copyEBIS->m_domainLevel = m_domainLevel;
    copyEBIS->m_nlevels     = m_nlevels;

    // Go through the levels setting the level dependent member data
    for (int ilev = 0; ilev < m_nlevels; ilev++)
    {
      // Original level dependent member data
      EBISLevel* origEBISLevel = m_ebisLevel[ilev];

      // New level dependent member data for this component
      EBISLevel* copyEBISLevel = new EBISLevel();

      // Copy everything except the EBData (this happens at the end)
      copyEBISLevel->m_phase = iEBIS;

      copyEBISLevel->m_grids  = origEBISLevel->m_grids;

      copyEBISLevel->m_domain = origEBISLevel->m_domain;

      copyEBISLevel->m_origin = origEBISLevel->m_origin;
      copyEBISLevel->m_dx     = origEBISLevel->m_dx;

      copyEBISLevel->m_tolerance = origEBISLevel->m_tolerance;

      EBGraphFactory graphFactory(origEBISLevel->m_domain);
      copyEBISLevel->m_graph.define(origEBISLevel->m_graph,graphFactory);

      // Clear this cache and all statistics associated with it
      copyEBISLevel->m_cache.clear();

      copyEBISLevel->m_cacheMisses = 0;
      copyEBISLevel->m_cacheHits   = 0;
      copyEBISLevel->m_cacheStale  = 0;

      // Remember this EBISLevel
      copyEBIS->m_ebisLevel.push_back(copyEBISLevel);
    }

    // Remember this EBIndexSpace
    connectedEBIS[iEBIS] = copyEBIS;
  }

  // Are there any cells which are regular with a multivalued parent
  bool anyRegularWithMultivaluedParent = false;

  // First pass to correct the EBIndexSpace for each conneceted component -
  // Remove AllRegular patches and regular cells that are in other connected
  // components, and mark all multiVoFs that are in other connected components
  pout() << "    EBIndexSpace::connectedComponentsNew-loop5" << endl;
  for (int ilev = 0; ilev < m_nlevels; ilev++)
  {
    CH_TIME("EBIndexSpace::connectedComponentsNew-loop5");
    pout() << "      Level " << ilev << endl;

    // The original EBISLevel
    EBISLevel* origEBISLevel = m_ebisLevel[ilev];

    // All the component numbers for this level
    const LevelData<EBCellFAB>& numberedComponents = *numberedComponentses[ilev];

    // All the EBISLevel's for the connected components
    Vector<EBISLevel*> copyEBISLevel(numEBIS);
    for (int iEBIS = minEBIS; iEBIS <= maxEBIS; iEBIS++)
    {
      copyEBISLevel[iEBIS] = connectedEBIS[iEBIS]->m_ebisLevel[ilev];
    }

    // Iterate through all the boxes
    for (DataIterator dit = origEBISLevel->m_grids.dataIterator(); dit.ok(); ++dit)
    {
      // Get some local references
      const Box& curBox = origEBISLevel->m_grids.get(dit());

      const EBCellFAB& numberedComponentsFAB = numberedComponents[dit()];

      const EBGraphImplem& origEBGraph = *(origEBISLevel->m_graph[dit()].m_implem);

      // If the original graph for this box was all regular
      if (origEBGraph.m_tag == EBGraphImplem::AllRegular)
      {
        const IntVect& oneIVInside = curBox.smallEnd();
        const int vofID = 0;
        const VolIndex vof(oneIVInside,vofID);

        // Get the component number for this (all regular) box
        Real componentNumber = numberedComponentsFAB(vof,comp);

        // Go through all the connected components and mark their graphs for
        // this box as "all covered" if they don't match "componentNumber"
        for (int iEBIS = minEBIS; iEBIS <= maxEBIS; iEBIS++)
        {
          if (iEBIS != componentNumber)
          {
            EBGraphImplem& copyEBGraph = *(copyEBISLevel[iEBIS]->m_graph[dit()].m_implem);
            copyEBGraph.m_tag = EBGraphImplem::AllCovered;
          }
        }
      }
      // If the original graph has irregular cells in this box
      else if (origEBGraph.m_tag == EBGraphImplem::HasIrregular)
      {
        // Iterate through each IntVect in this box
        IntVectSet curIVS(curBox);
        for (IVSIterator ivsit(curIVS); ivsit.ok(); ++ivsit)
        {
          // Current IntVect and it's GraphNode
          const IntVect& iv = ivsit();
          const GraphNode& origGraphNode = origEBGraph.m_graph(iv);

          // If the original GraphNode is all regular with a single valued
          // parent
          if (origGraphNode.isRegularWithSingleValuedParent())
          {
            const int vofID = 0;
            const VolIndex vof(iv,vofID);

            // Get the component number for this "vof"
            Real componentNumber = numberedComponentsFAB(vof,comp);

            // Go through all the connected components and mark their "vof" as
            // covered if they don't match "componentNumber"
            for (int iEBIS = minEBIS; iEBIS <= maxEBIS; iEBIS++)
            {
              if (iEBIS != componentNumber)
              {
                EBGraphImplem& copyEBGraph = *(copyEBISLevel[iEBIS]->m_graph[dit()].m_implem);
                copyEBGraph.m_graph(iv).m_cellList = (Vector<GraphNodeImplem>*)0;
              }
            }
          }
          // If the original GraphNode is anything else except covered
          else if (origGraphNode.hasValidCellList())
          {
            // Iterate through all the cellIndex's in the original GraphNode
            for (int iVoF = 0; iVoF < origGraphNode.size(); iVoF++)
            {
              // The current VoF
              const VolIndex vof(iv,iVoF);

              // Get the component number for this "vof"
              Real componentNumber = numberedComponentsFAB(vof,comp);

              // For each new EBIndexSpace, go to the GraphNodeImplem
              // corresponding to the current VoF and mark it "valid" if it
              // matches the current "componentNumber", otherwise mark it
              // "invalid"
              for (int iEBIS = minEBIS; iEBIS <= maxEBIS; iEBIS++)
              {
                EBGraphImplem&   copyEBGraph         = *(copyEBISLevel[iEBIS]->m_graph[dit()].m_implem);
                GraphNode&       copyGraphNode       = copyEBGraph.m_graph(iv);
                GraphNodeImplem& copyGraphNodeImplem = (*(copyGraphNode.m_cellList))[iVoF];

                if (copyGraphNodeImplem.m_isRegular)
                {
                  anyRegularWithMultivaluedParent = true;
                }

                if (iEBIS != componentNumber)
                {
                  copyGraphNodeImplem.m_isValid = false;
                }
                else
                {
                  copyGraphNodeImplem.m_isValid = true;
                }
              }
            }
          }
        }
      }
    }
  }

  // Second pass to fill the data holder with the new numbering for each
  // multiVoF in each connected component.  These can coexist in one (old)
  // data holder because the connected components are disjoint subsets of
  // the old graph.
  pout() << "    EBIndexSpace::connectedComponentsNew-loop6" << endl;
  for (int ilev = 0; ilev < m_nlevels; ilev++)
  {
    CH_TIME("EBIndexSpace::connectedComponentsNew-loop6");
    pout() << "      Level " << ilev << endl;

    // Data holder for the new numbering on this level
    LevelData<EBCellFAB>& newNumbering = *newNumberings[ilev];

    // Initial all the numberings on this level to 0
    EBLevelDataOps::setVal(newNumbering,0.0);

    // Go through each of the connected components
    for (int iEBIS = minEBIS; iEBIS <= maxEBIS; iEBIS++)
    {
      // Get the current EBIndexSpace and EBISLevel
      EBIndexSpace& curEBIS = *(connectedEBIS[iEBIS]);
      EBISLevel* curEBISLevel = curEBIS.m_ebisLevel[ilev];

      // Iterate through all the boxes on this level
      for (DataIterator dit = curEBISLevel->m_grids.dataIterator(); dit.ok(); ++dit)
      {
        // Get some local references
        const Box& curBox = curEBISLevel->m_grids.get(dit());

        EBCellFAB& newNumberingFAB = newNumbering[dit()];

        const EBGraphImplem& curEBGraph = *(curEBISLevel->m_graph[dit()].m_implem);

        // If the current graph has irregular cells then some renumbering may
        // be needed
        if (curEBGraph.m_tag == EBGraphImplem::HasIrregular)
        {
          // Iterate through all IntVect's in the current box
          IntVectSet curIVS(curBox);
          for (IVSIterator ivsit(curIVS); ivsit.ok(); ++ivsit)
          {
            const IntVect& iv = ivsit();
            const GraphNode& curGraphNode = curEBGraph.m_graph(iv);

            // If the graph has an explicit entry, i.e. a valid cellList, then
            // some changes may be needed
            if (curGraphNode.hasValidCellList())
            {
              // Go through the VoFs, if they are valid then record the
              // current validVoFCount (this VoF's new number in this
              // connected component) in the correct location in the data
              // holder based on the original graph and increment
              // validVoFCount.
              int validVoFCount = 0;
              for (int iVoF = 0; iVoF < curGraphNode.size(); iVoF++)
              {
                GraphNodeImplem& curGraphNodeImplem = (*(curGraphNode.m_cellList))[iVoF];
                if (curGraphNodeImplem.m_isValid == true)
                {
                  const VolIndex vof(iv,iVoF);

                  newNumberingFAB(vof,comp) = validVoFCount;
                  validVoFCount++;
                }
              }
            }
          }
        }
      }
    }

    // Make sure the ghostcells are correct on this level
    newNumbering.exchange();
  }

  // Third pass to correct the EBIndexSpace for each conneceted component on
  // each level (but not between levels) - Go through correcting all
  // references in the connected component using the map created in the second
  // pass.  Then prune away all parts of the graph that are not needed for
  // this connected component.
  pout() << "    EBIndexSpace::connectedComponentsNew-loop7" << endl;
  for (int iEBIS = minEBIS; iEBIS <= maxEBIS; iEBIS++)
  {
    CH_TIME("EBIndexSpace::connectedComponentsNew-loop7");
    pout() << "      EBIS #" << iEBIS << endl;

    EBIndexSpace& curEBIS = *(connectedEBIS[iEBIS]);

    // Go through each level
    for (int ilev = 0; ilev < m_nlevels; ilev++)
    {
      // Get the EBISLevel and new numbering for this lelvel
      EBISLevel* curEBISLevel = curEBIS.m_ebisLevel[ilev];
      LevelData<EBCellFAB>& newNumbering = *newNumberings[ilev];

      // Go through all the boxes in this level
      for (DataIterator dit = curEBISLevel->m_grids.dataIterator(); dit.ok(); ++dit)
      {
        const Box& curBox = curEBISLevel->m_grids.get(dit());
        EBCellFAB& newNumberingFAB = newNumbering[dit()];
        EBGraphImplem& curEBGraph = *(curEBISLevel->m_graph[dit()].m_implem);

        // If the graph in this box contain irregular cells, there may to
        // something to do
        if (curEBGraph.m_tag == EBGraphImplem::HasIrregular)
        {
          // Assume the new graph will be all covered unless we find a valid
          // VoF still in this box
          bool allCovered = true;

          // Iterator through all the IntVect's in this box looking for a
          // valid VoF
          IntVectSet curIVS(curBox);
          for (IVSIterator ivsit(curIVS); ivsit.ok(); ++ivsit)
          {
            const IntVect& iv = ivsit();
            GraphNode& curGraphNode = curEBGraph.m_graph(iv);

            // Nothing to fix but the graph isn't all covered
            if (curGraphNode.isRegularWithSingleValuedParent())
            {
              allCovered = false;
            }
            // If there are VoF's, go through them looking for valid ones and
            // correct the arcs connected that VoF to adjacent VoF's
            else if (curGraphNode.hasValidCellList())
            {
              // Get the old cellList and start building the new cellList
              Vector<GraphNodeImplem>* oldCellList = curGraphNode.m_cellList;
              Vector<GraphNodeImplem>* newCellList = new Vector<GraphNodeImplem>;

              // Go through all the VoF's
              for (int iVoF = 0; iVoF < curGraphNode.size(); iVoF++)
              {
                GraphNodeImplem& oldGraphNodeImplem = (*oldCellList)[iVoF];

                // If this VoF is valid correct the arcs to adjacent VoF's
                if (oldGraphNodeImplem.m_isValid == true)
                {
                  // Look in each direction at the lo and hi sides
                  for (int idir = 0; idir < SpaceDim; idir++)
                  {
                    for (SideIterator sit; sit.ok(); ++sit)
                    {
                      // Correct (in place) all the arcs for all the faces
                      // on all the sides
                      int index = oldGraphNodeImplem.index(idir,sit());

                      IntVect otherIV = iv;
                      otherIV.shift(idir,sign(sit()));

                      Vector<int>& curArcs = oldGraphNodeImplem.m_arc[index];
                      int numArcs = curArcs.size();

                      for (int iFace = 0; iFace < numArcs; iFace++)
                      {
                        if (curArcs[iFace] != -1)
                        {
                          VolIndex otherVoF(otherIV,curArcs[iFace]);
                          curArcs[iFace] = newNumberingFAB(otherVoF,comp);
                        }
                      }
                    }
                  }

                  // Save the new VoF
                  newCellList->push_back(oldGraphNodeImplem);
                }
              }

              // Delete the old cellList and save the new one
              delete oldCellList;
              oldCellList = newCellList;

              // Check the new number of VoF's
              int numCells = oldCellList->size();

              // If there are no VoF's, delete the cellList, mark the IntVect
              // location as covered, and remove the IntVect from the
              // irregular and multivalued IntVectSet's
              if (numCells == 0)
              {
                delete oldCellList;
                oldCellList = ((Vector<GraphNodeImplem>*) 0);

                (*curEBGraph.m_irregIVS) -= iv;
                (*curEBGraph.m_multiIVS) -= iv;
              }
              // If there are some VoF's
              else
              {
                // The graph in this box isn't all covered
                allCovered = false;

                // If there is only one VoF
                if (numCells == 1)
                {
                  // And it's not regular
                  if ((*oldCellList)[0].m_isRegular == false)
                  {
                    // Add it to the irregular IntVectSet and remove it from
                    // the multivalued IntVectSet
                    (*curEBGraph.m_irregIVS) |= iv;
                    (*curEBGraph.m_multiIVS) -= iv;
                  }
                }
              }

              // Save the updated cellList
              curGraphNode.m_cellList = oldCellList;
            }
          }

          // If the graph in this box (for this connected component) has
          // become all covered clean things up
          if (allCovered)
          {
            // Mark it as all covered
            curEBGraph.m_tag = EBGraphImplem::AllCovered;

            // Clear all the GraphNode's
            curEBGraph.m_graph.clear();

            // Delete the irregular and multivalued IntVectSet's
            if (curEBGraph.m_irregIVS != NULL)
            {
              delete curEBGraph.m_irregIVS;
              curEBGraph.m_irregIVS = NULL;
            }

            if (curEBGraph.m_multiIVS != NULL)
            {
              delete curEBGraph.m_multiIVS;
              curEBGraph.m_multiIVS = NULL;
            }
          }
        }
      }
    }
  }

  // Fourth pass to correct all the cell indices for the parent of each VoF.
  // To do this it is necessary to make a coarsened version of the current
  // (fine) level to place the data from the next coarser level (which may
  // have an incompatible DisjointBoxLayout).
  pout() << "    EBIndexSpace::connectedComponentsNew-loop8" << endl;
  for (int ilev = 0; ilev < m_nlevels-1; ilev++)
  {
    CH_TIME("EBIndexSpace::connectedComponentsNew-loop8");
    pout() << "      Level " << ilev << endl;

    // Get the new numberings for the parents of the current (fine) level
    LevelData<EBCellFAB>& newCoarNumbering = *newNumberings[ilev+1];

    // Coarsen the fine DisjointBoxLayout, create a temporary data holder
    // based on that, and copy the new numberings of the parents into this new
    // data holder
    const DisjointBoxLayout& fineGrids = newNumberings[ilev]->disjointBoxLayout();

    DisjointBoxLayout coarGrids;
    coarsen(coarGrids,fineGrids,nRef);

    EBISLayout coarEBISLayout;
    m_ebisLevel[ilev+1]->fillEBISLayout(coarEBISLayout,coarGrids,nGhosts);

    EBCellFactory coarFactory(coarEBISLayout);

    LevelData<EBCellFAB> newCoarsenedFineNumbering(coarGrids,nComps,ghostCells,coarFactory);

    newCoarNumbering.copyTo(newCoarsenedFineNumbering);

    // Go through all the connected components and fix the references to the
    // parents in their graphs
    for (int iEBIS = minEBIS; iEBIS <= maxEBIS; iEBIS++)
    {
      EBIndexSpace& curEBIS = *(connectedEBIS[iEBIS]);
      EBISLevel* fineEBISLevel = curEBIS.m_ebisLevel[ilev];

      // Iterate through all the boxes on this level making the corrections
      for (DataIterator dit = fineEBISLevel->m_grids.dataIterator(); dit.ok(); ++dit)
      {
        const Box& fineBox = fineEBISLevel->m_grids.get(dit());
        EBCellFAB& newNumberingFAB = newCoarsenedFineNumbering[dit()];
        EBGraphImplem& fineEBGraph = *(fineEBISLevel->m_graph[dit()].m_implem);

        // If there are irregular cells, there may be something to correct
        if (fineEBGraph.m_tag == EBGraphImplem::HasIrregular)
        {
          // Go through all the IntVect's in the current box
          IntVectSet fineIVS(fineBox);
          for (IVSIterator ivsit(fineIVS); ivsit.ok(); ++ivsit)
          {
            const IntVect& fineIV = ivsit();
            IntVect coarIV = coarsen(fineIV,nRef);

            GraphNode& fineGraphNode = fineEBGraph.m_graph(fineIV);

            // If there is a valid cellList then there will be references to
            // a parent VoF that need to be corrected
            if (fineGraphNode.hasValidCellList())
            {
              Vector<GraphNodeImplem>* fineCellList = fineGraphNode.m_cellList;

              // Correct the parent references in all the VoF's
              for (int iVoF = 0; iVoF < fineGraphNode.size(); iVoF++)
              {
                GraphNodeImplem& fineGraphNodeImplem = (*fineCellList)[iVoF];

                VolIndex coarVoF(coarIV,fineGraphNodeImplem.m_coarserNode);

                fineGraphNodeImplem.m_coarserNode = newNumberingFAB(coarVoF,comp);
              }
            }
          }
        }
      }
    }
  }

  // Fifth pass to correct all the VolIndex's giving the children of each VoF.
  // To do this it is necessary to make a refined version of the current
  // (coarse) level to place the data from the next finer level (which may
  // have an incompatible DisjointBoxLayout).
  pout() << "    EBIndexSpace::connectedComponentsNew-loop9" << endl;
  for (int ilev = 1; ilev < m_nlevels; ilev++)
  {
    CH_TIME("EBIndexSpace::connectedComponentsNew-loop9");
    pout() << "      Level " << ilev << endl;

    // Get the new numberings for the children of the current (coarse) level
    LevelData<EBCellFAB>& newFineNumbering = *newNumberings[ilev-1];

    // Refine the coarse DisjointBoxLayout, create a temporary data holder
    // based on that, and copy the new numberings of the children into this
    // new data holder
    const DisjointBoxLayout& coarGrids = newNumberings[ilev]->disjointBoxLayout();

    DisjointBoxLayout fineGrids;
    refine(fineGrids,coarGrids,nRef);

    EBISLayout fineEBISLayout;
    m_ebisLevel[ilev-1]->fillEBISLayout(fineEBISLayout,fineGrids,nGhosts);

    EBCellFactory fineFactory(fineEBISLayout);

    LevelData<EBCellFAB> newRefinedCoarNumbering(fineGrids,nComps,ghostCells,fineFactory);

    newFineNumbering.copyTo(newRefinedCoarNumbering);

    // Go through all the connected components and fix the references to the
    // children in their graphs
    for (int iEBIS = minEBIS; iEBIS <= maxEBIS; iEBIS++)
    {
      EBIndexSpace& curEBIS = *(connectedEBIS[iEBIS]);
      EBISLevel* coarEBISLevel = curEBIS.m_ebisLevel[ilev];

      // Iterate through all the boxes on this level making the corrections
      for (DataIterator dit = coarEBISLevel->m_grids.dataIterator(); dit.ok(); ++dit)
      {
        const Box& coarBox = coarEBISLevel->m_grids.get(dit());
        EBCellFAB& newNumberingFAB = newRefinedCoarNumbering[dit()];
        EBGraphImplem& coarEBGraph = *(coarEBISLevel->m_graph[dit()].m_implem);

        // If there are irregular cells, there may be something to correct
        if (coarEBGraph.m_tag == EBGraphImplem::HasIrregular)
        {
          // Go through all the IntVect's in the current box
          IntVectSet coarIVS(coarBox);
          for (IVSIterator ivsit(coarIVS); ivsit.ok(); ++ivsit)
          {
            const IntVect& coarIV = ivsit();

            GraphNode& coarGraphNode = coarEBGraph.m_graph(coarIV);

            // If there is a valid cellList then there will be references to
            // children VoF's that need to be corrected
            if (coarGraphNode.hasValidCellList())
            {
              Vector<GraphNodeImplem>* coarCellList = coarGraphNode.m_cellList;

              // Correct the children references in all the VoF's
              for (int iVoF = 0; iVoF < coarGraphNode.size(); iVoF++)
              {
                GraphNodeImplem& coarGraphNodeImplem = (*coarCellList)[iVoF];


                // Go through all the children making corrections
                int numChildren = coarGraphNodeImplem.m_finerNodes.size();
                for (int iChild = 0; iChild < numChildren; iChild++)
                {
                  VolIndex& curVolIndex = coarGraphNodeImplem.m_finerNodes[iChild];
                  const IntVect& fineIV = curVolIndex.gridIndex();

                  curVolIndex.define(fineIV,newNumberingFAB(curVolIndex,comp));
                }
              }
            }
          }
        }
      }
    }
  }

  // Sixth pass to copy all the relevant EBData into the new EBIndexSpace's.
  // Go from the coarsest level to the finest level.
  pout() << "    EBIndexSpace::connectedComponentsNew-loop10" << endl;
  for (int ilev = m_nlevels-1; ilev >= 0; ilev--)
  {
    CH_TIME("EBIndexSpace::connectedComponentsNew-loop10");
    pout() << "      Level " << ilev << endl;

    // For each connected component, initialize the volume data and face data
    // holders.  This has to be done using a graph with one layer of ghost
    // cells which is created there temporarily.
    for (int iEBIS = minEBIS; iEBIS <= maxEBIS; iEBIS++)
    {
      EBIndexSpace& copyEBIS   = *(connectedEBIS[iEBIS]);
      EBISLevel* copyEBISLevel = copyEBIS.m_ebisLevel[ilev];

      EBGraphFactory graphFactory(copyEBISLevel->m_domain);
      LevelData<EBGraph> ghostGraph(copyEBISLevel->m_grids,
                                    1,
                                    IntVect::Unit,
                                    graphFactory);
      Interval interval(0,0);
      copyEBISLevel->m_graph.copyTo(interval,ghostGraph,interval);

      EBDataFactory dataFactory;
      copyEBISLevel->m_data.define(copyEBISLevel->m_grids,
                                   1,
                                   IntVect::Zero,
                                   dataFactory);

      for (DataIterator dit = copyEBISLevel->m_grids.dataIterator(); dit.ok(); ++dit)
      {
//        copyEBISLevel->m_data[dit()].defineVoFData(ghostGraph[dit()],
//                                                   copyEBISLevel->m_grids.get(dit()));
//        copyEBISLevel->m_data[dit()].defineFaceData(ghostGraph[dit()],
//                                                    copyEBISLevel->m_grids.get(dit()));
        copyEBISLevel->m_data[dit()].define(ghostGraph[dit()],
                                            copyEBISLevel->m_grids.get(dit()),
                                            copyEBISLevel->m_dx, 
                                            copyEBISLevel->m_hasMoments);
      }
    }

    // Get the original EBISLevel to get access the original EBData
    EBISLevel* origEBISLevel = m_ebisLevel[ilev];

    // Get the component numbers and new VoF numbering so that the original
    // VolIndex's can be translated into the new VolIndex's for the correct
    // connected component.
    const LevelData<EBCellFAB>& numberedComponents = *numberedComponentses[ilev];
    const LevelData<EBCellFAB>& newNumbering       = *newNumberings       [ilev];

    // Iterate through all the boxes in this level copying the EBData to the
    // appropriate connected component
    for (DataIterator dit = origEBISLevel->m_grids.dataIterator(); dit.ok(); ++dit)
    {
      RefCountedPtr<EBDataImplem> origEBData = origEBISLevel->m_data[dit()].m_implem;

      const EBCellFAB& numberedComponentsFAB = numberedComponents[dit()];
      const EBCellFAB& newNumberingFAB       = newNumbering      [dit()];

      // Start with the volume data
      const BaseIVFAB<VolData>& origVolData = origEBData->m_volData;

      // Go through all the volume data a VoF at a time
      for (VoFIterator vofit(origVolData.getIVS(),origVolData.getEBGraph()); vofit.ok(); ++vofit)
      {
        // This is the original VolIndex
        VolIndex origVoF = vofit();

        int iEBIS = numberedComponentsFAB(origVoF,comp);

        // Only copy the volume data is the connected component number is
        // valid (if it is invalid then this connected component was discarded
        // because it had a total volume fracion of 0.0)
        if (iEBIS >= minEBIS && iEBIS <= maxEBIS)
        {
          // Construct the new VolIndex for the correct connected component
          int newCellIndex = newNumberingFAB(origVoF,comp);
          VolIndex newVoF(origVoF.gridIndex(),newCellIndex);

          EBIndexSpace&               copyEBIS      = *(connectedEBIS[iEBIS]);
          EBISLevel*                  copyEBISLevel = copyEBIS.m_ebisLevel[ilev];
          RefCountedPtr<EBDataImplem> copyEBData    = copyEBISLevel->m_data[dit()].m_implem;
          BaseIVFAB<VolData>&         copyVolData   = copyEBData->m_volData;

          // Copy the volume data to the correct place
          copyVolData(newVoF,comp) = origVolData(origVoF,comp);
        }
      }

      // Go through all the face data a direction and face at a time
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        const BaseIFFAB<FaceData>& origFaceData = (origEBData->m_faceData)[idir];
        for (FaceIterator faceit(origFaceData.getIVS(),origFaceData.getEBGraph(),idir,FaceStop::SurroundingWithBoundary); faceit.ok(); ++faceit)
        {
          // This is the original FaceIndex and the VolIndex's that define it
          FaceIndex origFace = faceit();
          const VolIndex& origLoVoF = origFace.getVoF(Side::Lo);
          const VolIndex& origHiVoF = origFace.getVoF(Side::Hi);

          // Start out with invalid new values and fill them in from the lo or
          // hi VoF depending on which is valid
          int iEBIS = -1;
          int newLoCellIndex = -1;
          int newHiCellIndex = -1;

          // The lo VoF is valid
          if (origLoVoF.cellIndex() != -1)
          {
            if (iEBIS == -1)
            {
              iEBIS = numberedComponentsFAB(origLoVoF,comp);
            }
            newLoCellIndex = newNumberingFAB(origLoVoF,comp);
          }

          // The hi VoF is valid
          if (origHiVoF.cellIndex() != -1)
          {
            if (iEBIS == -1)
            {
              iEBIS = numberedComponentsFAB(origHiVoF,comp);
            }
            newHiCellIndex = newNumberingFAB(origHiVoF,comp);
          }

          // Only copy the face data is the connected component number is
          // valid (if it is invalid then this connected component was
          // discarded because it had a total volume fracion of 0.0)
          if (iEBIS >= minEBIS && iEBIS <= maxEBIS)
          {
            // Construct the new FaceIndex for the correct connected component
            // (which includes new VolIndex's)
            VolIndex newLoVoF(origLoVoF.gridIndex(),newLoCellIndex);
            VolIndex newHiVoF(origHiVoF.gridIndex(),newHiCellIndex);

            FaceIndex newFace(newLoVoF,newHiVoF,origFace.direction());

            EBIndexSpace&               copyEBIS      = *(connectedEBIS[iEBIS]);
            EBISLevel*                  copyEBISLevel = copyEBIS.m_ebisLevel[ilev];
            RefCountedPtr<EBDataImplem> copyEBData    = copyEBISLevel->m_data[dit()].m_implem;
            BaseIFFAB<FaceData>&        copyFaceData  = (copyEBData->m_faceData)[idir];

            // Copy the face data to the correct place
            copyFaceData(newFace,comp) = origFaceData(origFace,comp);
          }
        }
      }
    }
  }

  // Clean up the temporary data holders defined using the original
  // EBIndexSpace
  pout() << "    EBIndexSpace::connectedComponentsNew-loop11" << endl;
  for (int ilev = 0; ilev < m_nlevels; ilev++)
  {
    CH_TIME("EBIndexSpace::connectedComponentsNew-loop11");
    pout() << "      Level " << ilev << endl;

    delete numberedComponentses[ilev];
    delete newNumberings[ilev];
  }

  if (anyRegularWithMultivaluedParent)
  {
    // TODO: Remove storage for any regular cells which had a multivalued
    // parent in the original graph but don't have a multivalued parent in
    // the connected componenet graphs.  The graphs are all correct without
    // this step but contain extra storage to explicitly store the cellIndex
    // of the parent and result in an extra level of indirection accessing
    // the parent.  On the other hand, these don't occur frequently.
  }

  // Return only the biggest connected component
  if (a_onlyBiggest)
  {
    RefCountedPtr<EBIndexSpace> biggestComponent(NULL);

    // If the EBIS is empty, then biggestEBIS = -1 and connectedEBIS.size() = 0.
    if (biggestEBIS >= 0)
      {
        biggestComponent = connectedEBIS[biggestEBIS];
      }

    connectedEBIS.resize(1);
    connectedEBIS[0] = biggestComponent;
  }

  // Return the Vector of RefCountedPtr<EBIndexSpace> - one for every
  // connected component
  return connectedEBIS;
}

// Divide the EBIndexSpace into connected components and return a Vector of
// disjoint EBIndexSpace's corresponding to each connected component.
Vector<RefCountedPtr<EBIndexSpace> > EBIndexSpace::connectedComponents()
{
  CH_TIME("EBIndexSpace::connectedComponents");

  int numComponents;
  bool onlyBiggest = false;

  return findConnectedComponents(numComponents,onlyBiggest);
}

// Divide the EBIndexSpace into connected components and return a Vector of
// disjoint EBIndexSpace's corresponding to each connected component.
Vector<RefCountedPtr<EBIndexSpace> > EBIndexSpace::connectedComponentsNew()
{
  CH_TIME("EBIndexSpace::connectedComponents2");

  int numComponents;
  bool onlyBiggest = false;

  return findConnectedComponentsNew(numComponents,onlyBiggest);
}

// Divide the EBIndexSpace into connected components and return the
// EBIndexSpace corresponding to the largest connected component.
RefCountedPtr<EBIndexSpace> EBIndexSpace::biggestConnectedComponent(int & a_numComponents)
{
  CH_TIME("EBIndexSpace::biggestConnectedComponent");

  bool onlyBiggest = true;

  Vector<RefCountedPtr<EBIndexSpace> > biggest = findConnectedComponents(a_numComponents,onlyBiggest);

  if (biggest.size() != 1)
  {
    MayDay::Error("EBIndexSpace::connectedComponents didn't return one EBIndexSpace when asked for only the biggest component");
  }

  return biggest[0];
}

// Divide the EBIndexSpace into connected components and return the
// EBIndexSpace corresponding to the largest connected component.
RefCountedPtr<EBIndexSpace> EBIndexSpace::biggestConnectedComponentNew(int & a_numComponents)
{
  CH_TIME("EBIndexSpace::biggestConnectedComponent2");

  bool onlyBiggest = true;

  Vector<RefCountedPtr<EBIndexSpace> > biggest = findConnectedComponentsNew(a_numComponents,onlyBiggest);

  if (biggest.size() != 1)
  {
    MayDay::Error("EBIndexSpace::connectedComponents didn't return one EBIndexSpace when asked for only the biggest component");
  }

  return biggest[0];
}

// Recursively find all VoFs connected to "a_curVoF" and number them
// "a_curNum".  This works on a single EBCellFAB.  "a_lastVoF" is used to
// avoid recursing back to the VoF that generated the call to this routine.
// For the first call, "a_lastVoF" should equal "a_curVoF" since this allows
// recursion using all faces.
bool EBIndexSpace::setAllConnectedVoFs(Real&               a_totalVolFrac,
                                       EBCellFAB&          a_curEBCellFAB,
                                       const EBGraph&      a_curEBGraph,
                                       const EBISBox&      a_curEBISBox,
                                       const VolIndex&     a_curVoF,
                                       const VolIndex&     a_lastVoF,
                                       const unsigned int& a_curNum)
{
  bool foundNewVoF = false;

  a_totalVolFrac = 0.0;

  int comp = 0;

  if (a_curEBCellFAB(a_curVoF,comp) == -1.0)
  {
    foundNewVoF = true;

    a_curEBCellFAB(a_curVoF,comp) = a_curNum;
    a_totalVolFrac += a_curEBISBox.volFrac(a_curVoF);

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
      {
        const Side::LoHiSide& curSide = sit();

        const Vector<FaceIndex> faces = a_curEBISBox.getFaces(a_curVoF,idir,curSide);

        for (int iface = 0; iface < faces.size(); iface++)
        {
          const FaceIndex& face = faces[iface];

          const VolIndex& nextVoF = face.getVoF(curSide);

          if (nextVoF.cellIndex() >= 0 && nextVoF != a_lastVoF)
          {
            Real totalVolFrac;

            setAllConnectedVoFs(totalVolFrac,
                                a_curEBCellFAB,
                                a_curEBGraph,
                                a_curEBISBox,
                                nextVoF,
                                a_curVoF,
                                a_curNum);

            a_totalVolFrac += totalVolFrac;
          }
        }
      }
    }
  }

  return foundNewVoF;
}

// Recursively find all VoFs connected to "a_curVoF" and number them
// "a_curNum".  This works on a single EBCellFAB.  "a_lastVoF" is used to
// avoid recursing back to the VoF that generated the call to this routine.
// For the first call, "a_lastVoF" should equal "a_curVoF" since this allows
// recursion using all faces.
bool EBIndexSpace::setAllConnectedVoFsNew(Real&               a_totalVolFrac,
                                       EBCellFAB&          a_curEBCellFAB,
                                       const Box&          a_curBox,
                                       const EBGraph&      a_curEBGraph,
                                       const EBISBox&      a_curEBISBox,
                                       const VolIndex&     a_curVoF,
                                       const VolIndex&     a_lastVoF,
                                       const unsigned int& a_curNum)
{
  bool foundNewVoF = false;

  a_totalVolFrac = 0.0;

  int comp = 2;

  if (a_curEBCellFAB(a_curVoF,comp) == -1.0)
  {
    foundNewVoF = true;

    a_curEBCellFAB(a_curVoF,comp) = a_curNum;
    a_totalVolFrac += a_curEBISBox.volFrac(a_curVoF);

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
      {
        const Side::LoHiSide& curSide = sit();

        const Vector<FaceIndex> faces = a_curEBISBox.getFaces(a_curVoF,idir,curSide);

        for (int iface = 0; iface < faces.size(); iface++)
        {
          const FaceIndex& face = faces[iface];

          const VolIndex& nextVoF = face.getVoF(curSide);

          if (a_curBox.contains(nextVoF.gridIndex()) && nextVoF.cellIndex() >= 0 && nextVoF != a_lastVoF)
          {
            Real totalVolFrac;

            setAllConnectedVoFsNew(totalVolFrac,
                                a_curEBCellFAB,
                                a_curBox,
                                a_curEBGraph,
                                a_curEBISBox,
                                nextVoF,
                                a_curVoF,
                                a_curNum);

            a_totalVolFrac += totalVolFrac;
          }
        }
      }
    }
  }

  return foundNewVoF;
}

// Undo everything done by "setAllConnectedVoFs" - this is needed if the total
// volume fraction found by "setAllConnectedVoFs" is 0.0 because then this
// connected component will be ignored.
void EBIndexSpace::resetAllConnectedVoFs(EBCellFAB&          a_curEBCellFAB,
                                         const EBGraph&      a_curEBGraph,
                                         const EBISBox&      a_curEBISBox,
                                         const VolIndex&     a_curVoF,
                                         const VolIndex&     a_lastVoF)
{
  int comp = 0;

  if (a_curEBCellFAB(a_curVoF,comp) != -1.0)
  {
    a_curEBCellFAB(a_curVoF,comp) = -1.0;

    for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
      {
        const Side::LoHiSide& curSide = sit();

        const Vector<FaceIndex> faces = a_curEBISBox.getFaces(a_curVoF,idir,curSide);

        for (int iface = 0; iface < faces.size(); iface++)
        {
          const FaceIndex& face = faces[iface];

          const VolIndex& nextVoF = face.getVoF(curSide);

          if (nextVoF.cellIndex() >= 0 && nextVoF != a_lastVoF)
          {
            resetAllConnectedVoFs(a_curEBCellFAB,
                                  a_curEBGraph,
                                  a_curEBISBox,
                                  nextVoF,
                                  a_curVoF);
          }
        }
      }
    }
  }
}

// Undo everything done by "setAllConnectedVoFsNew" - this is needed if the total
// volume fraction found by "setAllConnectedVoFsNew" is 0.0 because then this
// connected component will be ignored.
void EBIndexSpace::resetAllConnectedVoFsNew(EBCellFAB&          a_curEBCellFAB,
                                         const Box&          a_curBox,
                                         const EBGraph&      a_curEBGraph,
                                         const EBISBox&      a_curEBISBox,
                                         const VolIndex&     a_curVoF,
                                         const VolIndex&     a_lastVoF)
{
  int comp = 2;

  if (a_curEBCellFAB(a_curVoF,comp) != -1.0)
  {
    a_curEBCellFAB(a_curVoF,comp) = -1.0;
  
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      for (SideIterator sit; sit.ok(); ++sit)
      {
        const Side::LoHiSide& curSide = sit();
  
        const Vector<FaceIndex> faces = a_curEBISBox.getFaces(a_curVoF,idir,curSide);
  
        for (int iface = 0; iface < faces.size(); iface++)
        {
          const FaceIndex& face = faces[iface];
  
          const VolIndex& nextVoF = face.getVoF(curSide);
  
          if (nextVoF.cellIndex() >= 0 && nextVoF != a_lastVoF)
          {
            resetAllConnectedVoFsNew(a_curEBCellFAB,
                                  a_curBox,
                                  a_curEBGraph,
                                  a_curEBISBox,
                                  nextVoF,
                                  a_curVoF);
          }
        }
      }
    }
  }
}

#include "NamespaceFooter.H"
