/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Driver.cpp
  @brief  Implementation of CD_Driver.H
  @author Robert Marskar
*/

// Std includes
#include <fstream>
#include <iostream>
#include <sstream>

// Chombo includes
#include <EBArith.H>
#include <PolyGeom.H>
#include <EBAlias.H>
#include <LevelData.H>
#include <EBCellFAB.H>
#include <EBAMRIO.H>
#include <EBAMRDataOps.H>
#include <ParmParse.H>

// Our includes
#include <CD_Driver.H>
#include <CD_Random.H>
#include <CD_VofUtils.H>
#include <CD_DataOps.H>
#include <CD_BoxLoops.H>
#include <CD_MultifluidAlias.H>
#include <CD_Units.H>
#include <CD_MemoryReport.H>
#include <CD_Timer.H>
#include <CD_ParallelOps.H>
#include <CD_DischargeIO.H>
#include <CD_OpenMP.H>
#include <CD_NamespaceHeader.H>

Driver::Driver(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry,
               const RefCountedPtr<TimeStepper>&           a_timeStepper,
               const RefCountedPtr<AmrMesh>&               a_amr,
               const RefCountedPtr<CellTagger>&            a_cellTagger)
{
  CH_TIME("Driver::Driver");

  m_verbosity = -1;

  this->setComputationalGeometry(a_computationalGeometry); // Set computational geometry
  this->setTimeStepper(a_timeStepper);                     // Set time stepper
  this->setAmr(a_amr);                                     // Set amr
  this->setCellTagger(a_cellTagger);                       // Set cell tagger

  // AMR does its thing
  m_amr->sanityCheck();  // Sanity check, make sure everything is set up correctly
  m_amr->buildDomains(); // Build domains and resolutions, nothing else

  // Reset time steps.
  m_timeStep = 0;
  m_time     = 0.0;
  m_dt       = 0.0;
  m_outputDt = -1.0;

  m_profile      = false;
  m_doCoarsening = true;

  // Parse some class options and create the output directories for the simulation.
  this->parseOptions();

  // Always register this Realm and these operators.
  m_realm = Realm::Primal;
  m_amr->registerRealm(m_realm);
  m_amr->registerOperator(s_eb_fill_patch, m_realm, phase::gas);

  // Seed the RNG.
  Random::seed();

  m_time = m_startTime;
}

Driver::~Driver()
{
  CH_TIME("Driver::~Driver()");
}

int
Driver::getNumberOfPlotVariables() const
{
  CH_TIME("Driver::getNumberOfPlotVariables()");
  if (m_verbosity > 5) {
    pout() << "Driver::getNumberOfPlotVariables()" << endl;
  }

  int numPlotVars = 0;

  if (m_plotTags) {
    numPlotVars += 1;
  }
  if (m_plotRanks) {
    numPlotVars += m_amr->getRealms().size();
  }
  if (m_plotLevelset) {
    numPlotVars = numPlotVars + 2;
  }
  if (m_plotLoads) {
    numPlotVars += m_amr->getRealms().size();
  }

  return numPlotVars;
}

Vector<std::string>
Driver::getPlotVariableNames() const
{
  CH_TIME("Driver::getPlotVariableNames()");
  if (m_verbosity > 5) {
    pout() << "Driver::getPlotVariableNames()" << endl;
  }

  // TLDR: Need to write string identifiers for the internal variables that Driver
  //       will write to plot files. These are the cell tags, MPI ranks (for showing patch ownership)
  //       and level-sets. If we use multiple realms, MPI ownshiper of a grid patch will differ among
  //       the ranks.

  Vector<std::string> plotVarNames(0);

  if (m_plotTags) {
    plotVarNames.push_back("cell_tags");
  }
  if (m_plotRanks) {
    const std::string base = "_rank";

    for (const auto& str : m_amr->getRealms()) {
      const std::string id = str + base;
      plotVarNames.push_back(id);
    }
  }
  if (m_plotLevelset) {
    plotVarNames.push_back("levelset_gas");
    plotVarNames.push_back("levelset_solid");
  }
  if (m_plotLoads) {
    const std::string base = "_load";

    for (const auto& str : m_amr->getRealms()) {
      const std::string id = str + base;

      plotVarNames.push_back(id);
    }
  }

  return plotVarNames;
}

void
Driver::allocateInternals()
{
  CH_TIME("Driver::allocateInternals()");
  if (m_verbosity > 5) {
    pout() << "Driver::allocateInternals()" << endl;
  }

  // We define a DenseIntVectSet (because it's fast) to hold the tags that CellTagger creates. This is defined
  // over the primal realm.

  const int finestLevel = m_amr->getFinestLevel();

  m_tags.resize(1 + finestLevel);

  for (int lvl = 0; lvl <= finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    m_tags[lvl] = RefCountedPtr<LayoutData<DenseIntVectSet>>(new LayoutData<DenseIntVectSet>(dbl));

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      DenseIntVectSet& ivs = (*m_tags[lvl])[din];

      ivs = DenseIntVectSet(dbl[din], false);
    }
  }
}

void
Driver::cacheTags(const EBAMRTags& a_tags)
{
  CH_TIME("Driver::cacheTags(EBAMRTags)");
  if (m_verbosity > 5) {
    pout() << "Driver::cacheTags(EBAMRTags)" << endl;
  }

  constexpr int comp  = 0;
  constexpr int ncomp = 1;

  const int finestLevel = m_amr->getFinestLevel();
  const int ghost       = 0;

  m_amr->allocate(m_cachedTags, m_realm, ncomp, ghost);

  m_cachedTags.resize(1 + finestLevel);

  for (int lvl = 0; lvl <= finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    // Copy tags onto boolean mask. Unlike DenseIntVectSet, BaseFab<bool> can be put in LevelData<BaseFab<bool> > so
    // we can copy it to the new grid after regridding.
    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      BaseFab<bool>& cachedTags = (*m_cachedTags[lvl])[din];

      // Default is false, but we set cells to true if they were tagged earlier.
      cachedTags.setVal(false);

      const IntVectSet ivs = IntVectSet((*m_tags[lvl])[din]);
      for (IVSIterator ivsIt(ivs); ivsIt.ok(); ++ivsIt) {
        cachedTags(ivsIt(), comp) = true;
      }
    }
  }
}

void
Driver::getGeometryTags()
{
  CH_TIME("Driver::getGeometryTags()");
  if (m_verbosity > 5) {
    pout() << "Driver::getGeometryTags()" << endl;
  }

  // TLDR: This routine fetches cut-cell indexes using various criteria (supplied through the input script).

  const int maxAmrDepth = m_amr->getMaxAmrDepth();

  m_geomTags.resize(maxAmrDepth);

  const RefCountedPtr<EBIndexSpace>& ebisGas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebisSol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  // Note that we only need tags up to maxAmrDepth-1 since the grid on maxAmrDepth are generated from tags on maxAmrDepth-1
  for (int lvl = 0; lvl < maxAmrDepth; lvl++) {

    // Our level indexing disagrees with EBIS level indexing (where the finest level is on index 0). This is
    // a way to get the EBIndexSpace level from the current AMR Level.
    const ProblemDomain& curDomain = m_amr->getDomains()[lvl];     // Domain we're looking at right now.
    const int            curLevel  = ebisGas->getLevel(curDomain); // Corresponding level index in EBIndexSpace.

    IntVectSet condTags; // Irregular cells on electrodes
    IntVectSet dielTags; // Irregular cells on dielectrics

    // Conductor cells
    if (m_conductorTagsDepth > lvl) {
      condTags = ebisGas->irregCells(curLevel);
      if (!ebisSol.isNull()) {
        condTags |= ebisSol->irregCells(curLevel);
        condTags -= m_multifluidIndexSpace->interfaceRegion(curDomain);
      }
    }

    // Dielectric cells
    if (m_dielectricTagsDepth > lvl) {
      if (!ebisSol.isNull()) {
        dielTags = ebisSol->irregCells(curLevel);
      }
    }

    IntVectSet& geomTags = m_geomTags[lvl];
    geomTags.makeEmpty();

    // Evaluate angles between cut-cells and refine based on that.
    // Need one ghost cell because we fetch normal vectors from neighboring cut-cells.
    RefCountedPtr<EBIndexSpace> indexSpaces[2];
    indexSpaces[0] = ebisGas;
    indexSpaces[1] = ebisSol;

    for (int is = 0; is < 2; is++) {
      const RefCountedPtr<EBIndexSpace>& ebis = indexSpaces[is];

      if (!(ebis.isNull())) {
        DisjointBoxLayout irregGrids = ebis->getIrregGrids(curDomain);
        EBISLayout        ebisl;
#warning "Beware of this change -- I've filled two ghost cells to narrow gaps and find stencils that way"
        ebis->fillEBISLayout(ebisl, irregGrids, curDomain, 2);

        const RealVect     probLo = m_amr->getProbLo();
        const DataIterator dit    = irregGrids.dataIterator();

        const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(+ : geomTags)
        for (int mybox = 0; mybox < nbox; mybox++) {
          const DataIndex& din = dit[mybox];

          const Box        box     = irregGrids[din];
          const EBISBox&   ebisbox = ebisl[din];
          const EBGraph&   ebgraph = ebisbox.getEBGraph();
          const IntVectSet irreg   = ebisbox.getIrregIVS(box);

          VoFIterator vofit(irreg, ebgraph);

          auto kernel = [&](const VolIndex& vof) -> void {
            const IntVect  iv     = vof.gridIndex();
            const RealVect normal = ebisbox.normal(vof);

            // Check the angle between the normal vector in this irregular cell and neighboring irregular cells. If the
            // angle exceeds a specified threshold the cell is refined.
            const Vector<VolIndex> otherVofs = VofUtils::getVofsInRadius(vof,
                                                                         ebisbox,
                                                                         2,
                                                                         VofUtils::Connectivity::MonotonePath,
                                                                         false);

            for (int i = 0; i < otherVofs.size(); i++) {
              const VolIndex& curVof = otherVofs[i];

              if (ebisbox.isIrregular(curVof.gridIndex())) {
                const RealVect curNormal        = ebisbox.normal(curVof);
                const Real     degreesPerRadian = 180.0 / Units::pi;
                const Real     normLength       = normal.vectorLength();
                const Real     curNormLength    = curNormal.vectorLength();
                const Real     cosAngle         = PolyGeom::dot(normal, curNormal) / (normLength * curNormLength);
                const Real     theta            = acos(cosAngle) * degreesPerRadian;

                // Refine if angle exceeds threshold
                if (std::abs(theta) > m_refineAngle) {
                  geomTags |= iv;
                }
              }
            }
          };

          BoxLoops::loop(vofit, kernel);

          // Always refine multi-valued cells.
          geomTags |= ebisbox.getMultiCells(box);
        }
      }
    }

    // Things from depth specifications
    geomTags |= dielTags;
    geomTags |= condTags;
  }

  // Grow tags with specified factor.
  for (int lvl = 0; lvl < maxAmrDepth; lvl++) {
    m_geomTags[lvl].grow(m_irregTagGrowth);
  }

  // Processes may not agree what is the maximum tag depth. Make sure they're all on the same page.
  int deepestTagLevel = 0;
  for (int lvl = 0; lvl < m_geomTags.size(); lvl++) {
    if (!m_geomTags[lvl].isEmpty())
      deepestTagLevel = lvl;
  }

  m_geometricTagsDepth = ParallelOps::max(deepestTagLevel);
}

void
Driver::getCellsAndBoxes(long long&                       a_numLocalCells,
                         long long&                       a_numLocalCellsGhosts,
                         long long&                       a_numLocalBoxes,
                         long long&                       a_numTotalCells,
                         long long&                       a_numTotalCellsGhosts,
                         long long&                       a_numTotalBoxes,
                         Vector<long long>&               a_numLocalLevelBoxes,
                         Vector<long long>&               a_numTotalLevelBoxes,
                         Vector<long long>&               a_numLocalLevelCells,
                         Vector<long long>&               a_numTotalLevelCells,
                         const int&                       a_finestLevel,
                         const Vector<DisjointBoxLayout>& a_grids)
{
  CH_TIME("Driver::getCellsAndBoxes(long long x6, Vector<long long> x4, int, Vector<DisjointBoxLayout>)");
  if (m_verbosity > 5) {
    pout() << "Driver::getCellsAndBoxes(long long x6, Vector<long long> x4, int, Vector<DisjointBoxLayout>)" << endl;
  }

  a_numLocalCells       = 0;
  a_numLocalCellsGhosts = 0;
  a_numLocalBoxes       = 0;

  a_numTotalCells       = 0;
  a_numTotalCellsGhosts = 0;
  a_numTotalBoxes       = 0;

  a_numLocalLevelBoxes.resize(1 + a_finestLevel);
  a_numTotalLevelBoxes.resize(1 + a_finestLevel);
  a_numLocalLevelCells.resize(1 + a_finestLevel);
  a_numTotalLevelCells.resize(1 + a_finestLevel);

  const int ghost = m_amr->getNumberOfGhostCells();

  for (int lvl = 0; lvl <= a_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl  = a_grids[lvl];
    const DataIterator&      dit  = dbl.dataIterator();
    const int                nbox = dit.size();

    // Number of cells/boxes -- local to each MPI rank.
    long long numCellsNoGhosts   = 0;
    long long numCellsWithGhosts = 0;
    long long numBoxes           = 0;

#pragma omp parallel for schedule(runtime) reduction(+ : numCellsNoGhosts, numCellsWithGhosts, numBoxes)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const Box box      = dbl[din];
      const Box grownBox = grow(dbl[din], ghost);

      numCellsNoGhosts += box.numPts();
      numCellsWithGhosts += grownBox.numPts();
      numBoxes += 1;
    }

    const long long cellsThisLevel       = ParallelOps::sum(numCellsNoGhosts);
    const long long cellsThisLevelGhosts = ParallelOps::sum(numCellsWithGhosts);
    const long long boxesThisLevel       = ParallelOps::sum(numBoxes);

    // Total for this level
    a_numLocalCells += numCellsNoGhosts;
    a_numLocalCellsGhosts += numCellsWithGhosts;
    a_numLocalBoxes += numBoxes;

    a_numTotalCells += cellsThisLevel;
    a_numTotalCellsGhosts += cellsThisLevelGhosts;
    a_numTotalBoxes += boxesThisLevel;

    a_numLocalLevelBoxes[lvl] = numBoxes;
    a_numTotalLevelBoxes[lvl] = boxesThisLevel;
    a_numLocalLevelCells[lvl] = numCellsNoGhosts;
    a_numTotalLevelCells[lvl] = cellsThisLevel;
  }
}

void
Driver::gridReport()
{
  CH_TIME("Driver::gridReport()");
  if (m_verbosity > 5) {
    pout() << "Driver::gridReport()" << endl;
  }

  pout() << endl;

  const int                    finestLevel = m_amr->getFinestLevel();
  const Vector<ProblemDomain>& domains     = m_amr->getDomains();
  const Vector<Real>           dx          = m_amr->getDx();

  // Grid stuff goes into here
  long long totalCells;
  long long totalCellsGhosts;
  long long totalBoxes;
  long long localCells;
  long long localCellsGhosts;
  long long localBoxes;

  Vector<long long> localLevelBoxes;
  Vector<long long> totalLevelBoxes;
  Vector<long long> localLevelCells;
  Vector<long long> totalLevelCells;
  Vector<long long> validLevelCells;

  // Total number of grid points for a Cartesian grid covering entire finest domain. Used for "grid sparsity".
  const long long uniformPoints = (domains[finestLevel].domainBox()).numPts();

  // Some stuff
  const ProblemDomain coarsest_domain = m_amr->getDomains()[0];
  const ProblemDomain finest_domain   = m_amr->getDomains()[finestLevel];
  const Box           finestBox       = finest_domain.domainBox();
  const Box           coarsestBox     = coarsest_domain.domainBox();
  Vector<int>         refRat          = m_amr->getRefinementRatios();
  Vector<int>         ref_rat(finestLevel);
  for (int lvl = 0; lvl < finestLevel; lvl++) {
    ref_rat[lvl] = refRat[lvl];
  }

  this->getCellsAndBoxes(localCells,
                         localCellsGhosts,
                         localBoxes,
                         totalCells,
                         totalCellsGhosts,
                         totalBoxes,
                         localLevelBoxes,
                         totalLevelBoxes,
                         localLevelCells,
                         totalLevelCells,
                         finestLevel,
                         m_amr->getGrids(m_realm));

  long long totalValidCells = totalLevelCells[finestLevel];

  validLevelCells = totalLevelCells;

  for (int lvl = 0; lvl < finestLevel; lvl++) {
    validLevelCells[lvl] -= totalLevelCells[lvl + 1] / std::pow(refRat[lvl], SpaceDim);

    totalValidCells += validLevelCells[lvl];
  }

  int numThreads = 0;
#pragma omp parallel reduction(+ : numThreads)
  {
    numThreads += 1;
  }

  // Begin writing a report.
  pout() << "=======================================================================" << endl
         << "Driver::Grid report - timestep = " << m_timeStep << endl
         << "---------------------------------------------------------" << endl
#ifdef CH_MPI
         << "\tMPI ranks................ = " << numProc() << endl
#endif
#ifdef _OPENMP
         << "\tOpenMP threads........... = " << numThreads << endl
#endif
         << "\tFinest level............. = " << finestLevel << endl
#if CH_SPACEDIM == 2
         << "\tFinest AMR domain........ = " << finestBox.size()[0] << " x " << finestBox.size()[1] << endl
         << "\tCoarsest AMR domain...... = " << coarsestBox.size()[0] << " x " << coarsestBox.size()[1] << endl
#elif CH_SPACEDIM == 3
         << "\tFinest AMR domain........ = " << finestBox.size()[0] << " x " << finestBox.size()[1] << " x "
         << finestBox.size()[2] << endl
         << "\tCoarsest AMR domain...... = " << coarsestBox.size()[0] << " x " << coarsestBox.size()[1] << " x "
         << coarsestBox.size()[2] << endl
#endif
         << "\tRefinement ratios........ = " << ref_rat << endl
         << "\tGrid sparsity............ = " << 1.0 * totalCells / uniformPoints << endl
         << "\tFinest dx................ = " << dx[finestLevel] << endl
         << "\tTotal number boxes....... = " << DischargeIO::numberFmt(totalBoxes) << endl
         << "\tNumber of cells.......... = " << DischargeIO::numberFmt(totalCells) << endl
         << "\tIncluding ghost cells.... = " << DischargeIO::numberFmt(totalCellsGhosts) << endl
         << "\tValid # of cells......... = " << DischargeIO::numberFmt(totalValidCells) << endl
         << "\tTotal # of boxes (lvl)... = " << DischargeIO::numberFmt(totalLevelBoxes) << endl
         << "\tTotal # of cells (lvl)... = " << DischargeIO::numberFmt(totalLevelCells) << endl
         << "\tValid # of cells (lvl)... = " << DischargeIO::numberFmt(validLevelCells) << endl;

  // Do a local report for each Realm
  const std::vector<std::string> realms = m_amr->getRealms();

  for (const auto& str : realms) {
    this->getCellsAndBoxes(localCells,
                           localCellsGhosts,
                           localBoxes,
                           totalCells,
                           totalCellsGhosts,
                           totalBoxes,
                           localLevelBoxes,
                           totalLevelBoxes,
                           localLevelCells,
                           totalLevelCells,
                           finestLevel,
                           m_amr->getGrids(str));

    pout() << "\t**************" << endl
           << "\tRealm = " << str << endl
           << "\t...Proc. # of valid cells... = " << DischargeIO::numberFmt(localCells) << endl
           << "\t...Including ghost cells.... = " << DischargeIO::numberFmt(localCellsGhosts) << endl
           << "\t...Proc. # of boxes......... = " << DischargeIO::numberFmt(localBoxes) << endl
           << "\t...Proc. # of boxes (lvl)... = " << DischargeIO::numberFmt(localLevelBoxes) << endl
           << "\t...Proc. # of cells (lvl)... = " << DischargeIO::numberFmt(localLevelCells) << endl;
  }

  // Write a memory report if Chombo was to compiled to use memory tracking.
#ifdef CH_USE_MEMORY_TRACKING
  constexpr Real BytesPerMB = 1024.0 * 1024.0;

  long long localUnfreedMemory = 0LL;
  long long localPeakMemory    = 0LL;

  overallMemoryUsage(localUnfreedMemory, localPeakMemory);

  pout() << "\tUnfreed memory        = " << std::ceil(localUnfreedMemory / BytesPerMB) << " (MB)" << endl
         << "\tPeak memory usage     = " << std::ceil(localPeakMemory / BytesPerMB) << " (MB)" << endl;
#ifdef CH_MPI

  // If this is an MPI run we want to include the maximum consum memory in the report as well. We compute the
  // smallest/largest memory consumptions.
  const long long minUnfreedMemory = ParallelOps::min(localUnfreedMemory);
  const long long minPeakMemory    = ParallelOps::min(localPeakMemory);
  const long long maxUnfreedMemory = ParallelOps::max(localUnfreedMemory);
  const long long maxPeakMemory    = ParallelOps::max(localPeakMemory);

  pout() << "\tMin unfreed memory    = " << std::ceil(minUnfreedMemory / BytesPerMB) << " (MB)" << endl
         << "\tMin peak memory       = " << std::ceil(minPeakMemory / BytesPerMB) << " (MB)" << endl
         << "\tMax unfreed memory    = " << std::ceil(maxUnfreedMemory / BytesPerMB) << " (MB)" << endl
         << "\tMax peak memory       = " << std::ceil(maxPeakMemory / BytesPerMB) << " (MB)" << endl;
#endif
#endif
  pout() << "=======================================================================" << endl;

  pout() << endl;
}

void
Driver::regrid(const int a_lmin, const int a_lmax, const bool a_useInitialData)
{
  CH_TIMERS("Driver::regrid");
  CH_TIMER("Driver::regrid::compact_tags", t1);
  if (m_verbosity > 2) {
    pout() << "Driver::regrid" << endl;
  }

  // Use a timer here because I want to be able to put some diagnostics into this function.
  Timer timer("Driver::regrid(int, int, bool)");

  // We are allowing geometric tags to change under the hood, but we need a method for detecting if they changed. If they did,
  // we certainly have to regrid.
  timer.startEvent("Get geometry tags");
  Vector<IntVectSet> tags;

  if (m_needsNewGeometricTags) {
    this->getGeometryTags();
  }
  timer.stopEvent("Get geometry tags");

  timer.startEvent("Tag cells");
  const bool newCellTags = this->tagCells(tags, m_tags); // Tag cells using the cell tagger
  timer.stopEvent("Tag cells");

  if (!newCellTags && !m_needsNewGeometricTags) {
    if (a_useInitialData) {
      m_timeStepper->initialData();
    }

    if (m_verbosity > 1) {
      pout() << "\nDriver::regrid - Didn't find any new cell tags. Skipping the regrid step\n" << endl;
    }
    return;
  }
  else { // Compact tags
    CH_START(t1);
    for (int i = 0; i < tags.size(); i++) {
      tags[i].compact();
    }
    CH_STOP(t1);
  }

  // Store things that need to be regridded
  timer.startEvent("Pre-regrid");
  this->cacheTags(m_tags); // Cache m_tags because after regrid, ownership will change
  m_timeStepper->preRegrid(a_lmin, m_amr->getFinestLevel());
  if (!(m_cellTagger.isNull())) {
    m_cellTagger->preRegrid();
  }
  m_amr->preRegrid();
  timer.stopEvent("Pre-regrid");

  // Regrid AMR. Only levels [lmin, lmax] are allowed to change.
  timer.startEvent("Regrid AmrMesh");
  const int oldFinestLevel = m_amr->getFinestLevel();
  m_amr->regridAmr(tags, a_lmin, a_lmax);
  const int newFinestLevel = m_amr->getFinestLevel();
  timer.stopEvent("Regrid AmrMesh");

  // Load balance and regrid the various Realms
  timer.startEvent("Load balancing");
  const std::vector<std::string>& Realms = m_amr->getRealms();
  for (const auto& str : Realms) {
    if (m_timeStepper->loadBalanceThisRealm(str)) {

      Vector<Vector<int>> procs;
      Vector<Vector<Box>> boxes;

      m_timeStepper->loadBalanceBoxes(procs, boxes, str, m_amr->getProxyGrids(), a_lmin, newFinestLevel);

      m_amr->regridRealm(str, procs, boxes, a_lmin);
    }
  }
  timer.stopEvent("Load balancing");

  // Regrid the operators
  timer.startEvent("Regrid operators");
  m_amr->regridOperators(a_lmin);
  m_amr->postRegrid();
  timer.stopEvent("Regrid operators");

  // Regrid Driver, timestepper, and celltagger
  timer.startEvent("Regrid timestepper");
  this->regridInternals(oldFinestLevel, newFinestLevel);         // Regrid internals for Driver
  m_timeStepper->regrid(a_lmin, oldFinestLevel, newFinestLevel); // Regrid solvers
  if (a_useInitialData) {
    m_timeStepper->initialData();
  }
  timer.stopEvent("Regrid timestepper");

  // Regrid cell tagger if we have one.
  if (!m_cellTagger.isNull()) {
    timer.startEvent("Regrid celltagger");
    m_cellTagger->regrid();
    timer.stopEvent("Regrid celltagger");
  }

  // If it wants to, TimeStepper can do a postRegrid operation.
  timer.startEvent("Post-regrid");
  m_timeStepper->postRegrid();
  timer.stopEvent("Post-regrid");

  m_needsNewGeometricTags = false;

  if (m_profile) {
    timer.eventReport(pout(), true);
  }
}

void
Driver::regridInternals(const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("Driver::regridInternals(int, int)");
  if (m_verbosity > 2) {
    pout() << "Driver::regridInternals(int, int)" << endl;
  }

  constexpr int comp  = 0;
  constexpr int nComp = 1;

  // Allocate internals -- this allocates m_tags
  this->allocateInternals();

  // Copy cached tags back over to m_tags
  for (int lvl = 0; lvl <= std::min(a_oldFinestLevel, a_newFinestLevel); lvl++) {
    const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit  = dbl.dataIterator();
    const int                nbox = dit.size();

    LevelData<BaseFab<bool>> tmp;
    tmp.define(dbl, nComp, IntVect::Zero);

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      tmp[din].setVal(false);
    }

    // Ideally we'd use a pre-defined Copier here. But since this is only called once per regrid it should be fine.
    m_cachedTags[lvl]->copyTo(tmp);

    // m_cachedTags was allocated on the grids while tmp is allocated on the new. Look through
    // the bools in tmp and reconstruct the tags.
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const BaseFab<bool>& tmpFab = tmp[din];
      const Box&           box    = dbl[din];

      DenseIntVectSet& tags = (*m_tags[lvl])[din];

      auto kernel = [&](const IntVect& iv) -> void {
        if (tmpFab(iv, comp)) {
          tags |= iv;
        }
      };

      BoxLoops::loop(box, kernel);
    }
  }
}

void
Driver::run(const Real a_startTime, const Real a_endTime, const int a_maxSteps)
{
  CH_TIME("Driver::run(Real, Real, int)");
  if (m_verbosity > 1) {
    pout() << "Driver::run(Real, Real, int)" << endl;
  }

  constexpr Real abortThreshold = 1.E-5;

  if (m_verbosity > 0) {
    pout() << "=================================" << endl;
    if (!m_restart) {
      pout() << "Driver::run -- starting run" << endl;
    }
    else {
      pout() << "Driver::run -- restarting run" << endl;
    }
  }

  if (a_maxSteps > 0) {
    if (!m_restart) {
      m_time     = a_startTime;
      m_timeStep = 0;
    }

    m_dt = m_timeStepper->computeDt();
    m_timeStepper->synchronizeSolverTimes(m_timeStep, m_time, m_dt);

    bool isLastStep  = false;
    bool isFirstStep = true;

    // Store the initial dt in case we have to abort (in case the time step became too small for some reason).
    const Real initDt = m_dt;

    if (m_verbosity > 0) {
      this->gridReport();
    }

    // Store the initial time before we starting stepping. This is used for giving an estimate
    // for how long the simulation will run.
    m_wallClockStart = Timer::wallClock();

    while (m_time < a_endTime && m_timeStep < a_maxSteps && !isLastStep) {
      const int maxSimDepth = m_amr->getMaxSimulationDepth();
      const int maxAmrDepth = m_amr->getMaxAmrDepth();

      // Check if we should regrid. This can be due to regular intervals, or the TimeStepper can call for it.
      const bool canRegrid         = maxSimDepth > 0 && maxAmrDepth > 0 && m_regridInterval > 0;
      const bool regridStep        = m_timeStep % m_regridInterval == 0;
      const bool regridTimeStepper = m_timeStepper->needToRegrid();
      if (canRegrid && (regridStep || regridTimeStepper)) {
        if (!isFirstStep) {

          // Regrid all levels, but restrict the addition to one at a time. As always, new grids on level l are generated through tags
          // on levels (l-1);
          // This means that if we refine, we can only add one level at a time.
          const int lmin = 0;
          const int lmax = m_amr->getFinestLevel() + 1;

          if (m_writeRegridFiles) {
            this->writePreRegridFile();
          }

          this->regrid(lmin, lmax, false);

          // Write a grid report, displaying information about the new grids. Can also write
          // a regrid file if necessary.
          if (m_verbosity > 0) {
            this->gridReport();
          }
          if (m_writeRegridFiles) {
            this->writePostRegridFile();
          }
        }
      }

      // Compute a time step for the TimeStepper::advance(...) method.
      if (!isFirstStep) {
        m_dt = m_timeStepper->computeDt();
      }
      else { // In this case we already had one.
        isFirstStep = false;
      }

      // Abort if the time step became too small.
      if (m_dt < abortThreshold * initDt) {
        m_timeStep++;

        if (m_writeMemory) {
          this->writeMemoryUsage();
        }
        if (m_writeLoads) {
          this->writeComputationalLoads();
        }
        this->writeCrashFile();

        MayDay::Error("Driver::run(Real, Real, int) - the time step became too small.");
      }

      // Adjust time step to conform with the specified output intervals
      bool writePltFile = false;

      if (m_outputDt > 0.0) {
        const int k = std::floor(m_time / m_outputDt);

        Real lastOutputTime = k * m_outputDt;
        Real nextOutputTime = (k + 1) * m_outputDt;

        const Real thresh = 1.E-10 * m_outputDt;

        // Weird, but can happen due to flooring when m_time is an integer multiple of m_outputDt
        if (std::abs(m_time - nextOutputTime) < thresh) {
          lastOutputTime += m_outputDt;
          nextOutputTime += m_outputDt;
        }

        // Adjust dt so that we land on the next output interval
        if ((m_time + m_dt) >= nextOutputTime) {
          m_dt = nextOutputTime - m_time;
        }

        // Write plot file if we're landing on the next output time.
        if (std::abs((m_time + m_dt) - nextOutputTime) <= thresh) {
          writePltFile = true;
        }
      }

      // Adjust last time step -- it can be smaller than the one we computed because we want to
      // end the simulation at a_endTime.
      if (m_time + m_dt > a_endTime) {
        m_dt = a_endTime - m_time;

        isLastStep = true;
      }

      // Time stepper advances solutions. Note that the time stepper can choose to use a time step different
      // from the one we computed (because some time-steppers use adaptive time-stepping).
      m_wallClockOne      = Timer::wallClock();
      const Real actualDt = m_timeStepper->advance(m_dt);
      m_wallClockTwo      = Timer::wallClock();

      // Synchronize times
      m_dt = actualDt;
      m_time += actualDt;
      m_timeStep += 1;
      m_timeStepper->synchronizeSolverTimes(m_timeStep, m_time, m_dt);

      // Check if this was the last step.
      if (std::abs(m_time - a_endTime) < m_dt * 1.E-5) {
        isLastStep = true;
      }
      if (m_timeStep == m_maxSteps) {
        isLastStep = true;
      }

      // Print step report
      if (m_verbosity > 0) {
        this->stepReport(a_startTime, a_endTime, a_maxSteps);
      }

#ifdef CH_USE_HDF5
      // Write plot files, memory files, loads, checkpoint etc.
      if (m_plotInterval > 0) {

        // Aux data
        if (m_writeMemory) {
          this->writeMemoryUsage();
        }
        if (m_writeLoads) {
          this->writeComputationalLoads();
        }

        // Check if we should output this time step anyways.
        if (m_outputDt <= 0.0) {
          writePltFile = (m_timeStep % m_plotInterval == 0) || isLastStep == true;
        }

        if (writePltFile) {
          if (m_verbosity > 2) {
            pout() << "Driver::run -- Writing plot file" << endl;
          }

          // Write plot file.
          this->writePlotFile();
        }
      }

      // Write checkpoint file
      if (m_checkpointInterval > 0) {
        if (m_timeStep % m_checkpointInterval == 0 || isLastStep == true) {
          if (m_verbosity > 2) {
            pout() << "Driver::run -- Writing checkpoint file" << endl;
          }
          this->writeCheckpointFile();
        }
      }
#endif

      // Rebuild the ParmParse table and read input parameters again. Some parameters are allowed to change during runtime.
      this->rebuildParmParse();

      this->parseRuntimeOptions();
      m_amr->parseRuntimeOptions();
      m_timeStepper->parseRuntimeOptions();
      if (!m_cellTagger.isNull()) {
        m_cellTagger->parseRuntimeOptions();
      }
    }
  }

  if (m_verbosity > 0) {
    pout() << "==================================" << endl;
    pout() << "Driver::run -- ending run  " << endl;
    pout() << "==================================" << endl;

    MemoryReport::getMaxMinMemoryUsage();
  }
}

void
Driver::setupAndRun(const std::string a_inputFile)
{
  CH_TIME("Driver::setupAndRun(std::string)");
  if (m_verbosity > 0) {
    pout() << "Driver::setupAndRun(std::string)" << endl;
  }

  // TLDR: Call setup(...), which will select among the various setup functions.

  char iter_str[100];
  sprintf(iter_str, ".check%07d.%dd.hdf5", m_restartStep, SpaceDim);
  const std::string restartFile = m_outputDirectory + "/chk/" + m_outputFileNames + std::string(iter_str);

  this->setup(a_inputFile, m_initialRegrids, m_restart, restartFile);

  // Run the simulation.
  if (!m_geometryOnly) {
    this->run(m_startTime, m_stopTime, m_maxSteps);
  }
}

void
Driver::rebuildParmParse() const
{
  CH_TIME("Driver::rebuildParmParse()");
  if (m_verbosity > 5) {
    pout() << "Driver::rebuildParmParse()" << endl;
  }

  ParmParse pp("Driver");

  pp.redefine(m_inputFile.c_str());
}

void
Driver::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry)
{
  CH_TIME("Driver::setComputationalGeometry(RefCountedPtr<ComputationalGeometry>");
  if (m_verbosity > 5) {
    pout() << "Driver::setComputationalGeometry(RefCountedPtr<ComputationalGeometry>)" << endl;
  }

  CH_assert(!a_computationalGeometry.isNull());

  m_computationalGeometry = a_computationalGeometry;
  m_multifluidIndexSpace  = a_computationalGeometry->getMfIndexSpace();
}

void
Driver::setTimeStepper(const RefCountedPtr<TimeStepper>& a_timeStepper)
{
  CH_TIME("Driver::setTimeStepper(RefCountedPtr<TimeStepper>)");
  if (m_verbosity > 5) {
    pout() << "Driver::setTimeStepper(RefCountedPtr<TimeStepper>)" << endl;
  }

  CH_assert(!a_timeStepper.isNull());

  m_timeStepper = a_timeStepper;
}

void
Driver::setCellTagger(const RefCountedPtr<CellTagger>& a_cellTagger)
{
  CH_TIME("Driver::setCellTagger(RefCountedPtr<CellTagger>)");
  if (m_verbosity > 5) {
    pout() << "Driver::setCellTagger(RefCountedPtr<CellTagger>)" << endl;
  }

  m_cellTagger = a_cellTagger;

  if (!a_cellTagger.isNull()) {
    m_cellTagger->parseOptions();
  }
}

void
Driver::parseOptions()
{
  CH_TIME("Driver::parseOptions()");

  ParmParse pp("Driver");

  pp.get("verbosity", m_verbosity);
  pp.query("profile", m_profile);
  if (m_verbosity > 5) {
    pout() << "Driver::parseOptions()" << endl;
  }

  pp.get("regrid_interval", m_regridInterval);
  pp.get("initial_regrids", m_initialRegrids);
  pp.get("restart", m_restartStep);
  pp.get("write_memory", m_writeMemory);
  pp.get("write_loads", m_writeLoads);
  pp.get("output_directory", m_outputDirectory);
  pp.get("output_names", m_outputFileNames);
  pp.get("plot_interval", m_plotInterval);
  pp.get("checkpoint_interval", m_checkpointInterval);
  pp.get("write_regrid_files", m_writeRegridFiles);
  pp.get("write_restart_files", m_writeRestartFiles);
  pp.get("num_plot_ghost", m_numPlotGhost);
  pp.get("allow_coarsening", m_allowCoarsening);
  pp.get("geometry_only", m_geometryOnly);
  pp.get("ebis_memory_load_balance", m_ebisMemoryLoadBalance);
  pp.get("max_steps", m_maxSteps);
  pp.get("start_time", m_startTime);
  pp.get("stop_time", m_stopTime);
  pp.get("max_plot_depth", m_maxPlotLevel);
  pp.get("max_chk_depth", m_maxCheckpointDepth);
  pp.get("do_init_load_balance", m_doInitLoadBalancing);
  pp.get("output_dt", m_outputDt);

  m_restart = (m_restartStep > 0) ? true : false;

  // Stuff that's a little too verbose to include here directly.
  this->parsePlotVariables();
  this->parseGeometryGeneration();
  this->parseGeometryRefinement();
  this->parseIrregTagGrowth();

  // Not a required thing.
  pp.query("coarsening", m_doCoarsening);
}

void
Driver::parseRuntimeOptions()
{
  CH_TIME("Driver::parseRuntimeOptions()");

  ParmParse pp("Driver");

  pp.get("verbosity", m_verbosity);
  if (m_verbosity > 5) {
    pout() << "Driver::parseRuntimeOptions()" << endl;
  }
  pp.get("write_memory", m_writeMemory);
  pp.get("write_loads", m_writeLoads);
  pp.get("plot_interval", m_plotInterval);
  pp.get("regrid_interval", m_regridInterval);
  pp.get("checkpoint_interval", m_checkpointInterval);
  pp.get("write_regrid_files", m_writeRegridFiles);
  pp.get("write_restart_files", m_writeRestartFiles);
  pp.get("num_plot_ghost", m_numPlotGhost);
  pp.get("allow_coarsening", m_allowCoarsening);
  pp.get("max_plot_depth", m_maxPlotLevel);
  pp.get("max_steps", m_maxSteps);
  pp.get("stop_time", m_stopTime);
  pp.get("output_dt", m_outputDt);

  this->parseGeometryRefinement();
  this->parsePlotVariables();
  this->parseIrregTagGrowth();
}

void
Driver::parsePlotVariables()
{
  CH_TIME("Driver::parsePlotVariables()");
  if (m_verbosity > 5) {
    pout() << "Driver::parsePlotVariables()" << endl;
  }

  // False by default.
  m_plotTags     = false;
  m_plotRanks    = false;
  m_plotLevelset = false;
  m_plotLoads    = false;

  // Check input script to see which variables to include in plot files.
  ParmParse pp("Driver");
  const int num = pp.countval("plt_vars");

  Vector<std::string> str(num);
  pp.getarr("plt_vars", str, 0, num);

  for (int i = 0; i < num; i++) {
    if (str[i] == "tags") {
      m_plotTags = true;
    }
    else if (str[i] == "mpi_rank") {
      m_plotRanks = true;
    }
    else if (str[i] == "levelset") {
      m_plotLevelset = true;
    }
    else if (str[i] == "loads") {
      m_plotLoads = true;
    }
  }
}

void
Driver::parseIrregTagGrowth()
{
  CH_TIME("Driver::parseIrregTagGrowth()");
  if (m_verbosity > 5) {
    pout() << "Driver::parseIrregTagGrowth()" << endl;
  }

  ParmParse pp("Driver");

  pp.get("grow_geo_tags", m_irregTagGrowth);

  m_irregTagGrowth = std::max(0, m_irregTagGrowth);
}

void
Driver::parseGeometryGeneration()
{
  CH_TIME("Driver::parseGeometryGeneration()");
  if (m_verbosity > 5) {
    pout() << "Driver::parseGeometryGeneration()" << endl;
  }

  ParmParse pp("Driver");

  pp.get("geometry_generation", m_geometryGeneration);
  pp.get("geometry_scan_level", m_geoScanLevel);

  if (!(m_geometryGeneration == "chombo-discharge" || m_geometryGeneration == "chombo")) {
    MayDay::Abort("Driver:parseGeometryGeneration - unsupported argument requested");
  }
}

void
Driver::parseGeometryRefinement()
{
  CH_TIME("Driver::parseGeometryRefinement()");
  if (m_verbosity > 5) {
    pout() << "Driver::parseGeometryRefinement()" << endl;
  }

  // This routine parses the default depth at which we refine EBs, and also the refinement angle. This is used for setting up grids where only flags on the
  // EB are involved.
  //
  ParmParse pp("Driver");

  const auto c1 = m_refineAngle;
  const auto c2 = m_conductorTagsDepth;
  const auto c3 = m_dielectricTagsDepth;

  pp.get("refine_angles", m_refineAngle);
  pp.get("refine_electrodes", m_conductorTagsDepth);
  pp.get("refine_dielectrics", m_dielectricTagsDepth);

  if (m_conductorTagsDepth < 0) {
    m_conductorTagsDepth = m_amr->getMaxAmrDepth();
  }
  if (m_dielectricTagsDepth < 0) {
    m_dielectricTagsDepth = m_amr->getMaxAmrDepth();
  }

  // We are also allowing geometry refinement criteria to change as simulations progress (i.e. m_timeStep > 0) where we call this routine again. But we need
  // something to tell us that we got new refinement criteria for geometric tags so we can avoid regrid if they didn't change. This is my clunky way of doing that.
  if (m_timeStep >
      0) { // Simulation is already running, and we need to check if we need new geometric tags for regridding.
    if (c1 != m_refineAngle || c2 != m_conductorTagsDepth || c3 != m_dielectricTagsDepth) {
      m_needsNewGeometricTags = true;
    }
  }
  else { // Fresh simulation -- we should have gotten geometric tags before entering regrid so we can set this to false.
    m_needsNewGeometricTags = false;
  }
}

void
Driver::createOutputDirectories()
{
  CH_TIME("Driver::createOutputDirectories()");
  if (m_verbosity > 5) {
    pout() << "Driver::createOutputDirectories()" << endl;
  }

  // TLDR: Master rank creates the output directories for everyone.
  if (procID() == 0) {
    int success;

    std::string cmd;

    cmd     = "mkdir -p " + m_outputDirectory;
    success = system(cmd.c_str());
    if (success != 0) {
      std::cout << "Driver::createOutputDirectories - master could not create directory" << std::endl;
    }

    cmd     = "mkdir -p " + m_outputDirectory + "/plt";
    success = system(cmd.c_str());
    if (success != 0) {
      std::cout << "Driver::createOutputDirectories - master could not create plot directory" << std::endl;
    }

    cmd     = "mkdir -p " + m_outputDirectory + "/geo";
    success = system(cmd.c_str());
    if (success != 0) {
      std::cout << "Driver::createOutputDirectories - master could not create geo directory" << std::endl;
    }

    cmd     = "mkdir -p " + m_outputDirectory + "/chk";
    success = system(cmd.c_str());
    if (success != 0) {
      std::cout << "Driver::createOutputDirectories - master could not create checkpoint directory" << std::endl;
    }

    cmd     = "mkdir -p " + m_outputDirectory + "/mpi";
    success = system(cmd.c_str());
    if (success != 0) {
      std::cout << "Driver::createOutputDirectories - master could not create mpi directory" << std::endl;
    }

    cmd     = "mkdir -p " + m_outputDirectory + "/mpi/memory";
    success = system(cmd.c_str());
    if (success != 0) {
      std::cout << "Driver::createOutputDirectories - master could not create mpi/memory directory" << std::endl;
    }

    cmd     = "mkdir -p " + m_outputDirectory + "/mpi/loads";
    success = system(cmd.c_str());
    if (success != 0) {
      std::cout << "Driver::createOutputDirectories - master could not create mpi/loads directory" << std::endl;
    }

    cmd     = "mkdir -p " + m_outputDirectory + "/regrid";
    success = system(cmd.c_str());
    if (success != 0) {
      std::cout << "Driver::createOutputDirectories - master could not create regrid directory" << std::endl;
    }

    cmd     = "mkdir -p " + m_outputDirectory + "/restart";
    success = system(cmd.c_str());
    if (success != 0) {
      std::cout << "Driver::createOutputDirectories - master could not create restart directory" << std::endl;
    }

    cmd     = "mkdir -p " + m_outputDirectory + "/crash";
    success = system(cmd.c_str());
    if (success != 0) {
      std::cout << "Driver::createOutputDirectories - master could not create crash directory" << std::endl;
    }
  }

  ParallelOps::barrier();
}

void
Driver::setAmr(const RefCountedPtr<AmrMesh>& a_amrMesh)
{
  CH_TIME("Driver::setAmr(RefCountedPtr<AmrMesh>)");
  if (m_verbosity > 5) {
    pout() << "Driver::setAmr(RefCountedPtr<AmrMesh>)" << endl;
  }

  CH_assert(!a_amrMesh.isNull());

  m_amr = a_amrMesh;
  m_amr->setMultifluidIndexSpace(m_computationalGeometry->getMfIndexSpace());
}

void
Driver::setup(const std::string a_inputFile,
              const int         a_initialRegrids,
              const bool        a_restart,
              const std::string a_restartFile)
{
  CH_TIME("Driver::setup(string, int, bool, string)");
  if (m_verbosity > 5) {
    pout() << "Driver::setup(string, int, bool, string)" << endl;
  }

  // We store the input file because we want ParmParse to read parameters
  // at every time step (for run-time configuration purposes).
  m_inputFile = a_inputFile;

  this->createOutputDirectories();

  if (m_geometryOnly) {
    this->setupGeometryOnly();
  }
  else {
    if (!a_restart) {
      this->setupFresh(a_initialRegrids);
#ifdef CH_USE_HDF5
      if (m_plotInterval > 0) {
        if (m_writeMemory) {
          this->writeMemoryUsage();
        }
        if (m_writeLoads) {
          this->writeComputationalLoads();
        }

        this->writePlotFile();
      }
#endif
    }
    else {
#ifdef CH_USE_HDF5
      this->setupForRestart(a_initialRegrids, a_restartFile);
#else
      this->setupFresh(a_initialRegrids);
#endif
    }
  }
}

void
Driver::setupGeometryOnly()
{
  CH_TIME("Driver::setupGeometryOnly");
  if (m_verbosity > 5) {
    pout() << "Driver::setupGeometryOnly" << endl;
  }

  this->sanityCheck();

  const Real t0 = Timer::wallClock();

  // Need to activate some flags that trigger Chombo or chombo-discharge geo-generation method.
  if (m_geometryGeneration == "chombo-discharge") {

    // We will run with ScanShop, which builds the EBIndexSpace map using recursive pruning of
    // regions that don't contain cut-cells. The user asks for a coarsening/refinement of the
    // base AMR level, which we make here.
    ProblemDomain scanDomain;
    if (m_geoScanLevel < 0) {
      // Make a coarsening of the base AMR level, but don't coarsen beyond 2x2 domains.
      scanDomain  = m_amr->getDomains()[0];
      int numCoar = 0;
      while (scanDomain.domainBox().shortside() >= 4 && numCoar < std::abs(m_geoScanLevel)) {
        scanDomain.coarsen(2);
        numCoar++;
      }
    }
    else {
      const int amrLevel = std::min(m_geoScanLevel, m_amr->getMaxAmrDepth());
      scanDomain         = m_amr->getDomains()[amrLevel];
    }

    EBISLevel::s_distributedData = true;

    m_computationalGeometry->useScanShop(scanDomain);
  }
  else if (m_geometryGeneration == "chombo") {
    m_computationalGeometry->useChomboShop();

    if (m_ebisMemoryLoadBalance) {
      EBIndexSpace::s_useMemoryLoadBalance = true;
    }
    else {
      EBIndexSpace::s_useMemoryLoadBalance = false;
    }
  }

  const int numCoarsenings = m_doCoarsening ? -1 : m_amr->getMaxAmrDepth();

  m_computationalGeometry->buildGeometries(m_amr->getFinestDomain(),
                                           m_amr->getProbLo(),
                                           m_amr->getFinestDx(),
                                           m_amr->getMaxEbisBoxSize(),
                                           m_amr->getNumberOfEbGhostCells(),
                                           numCoarsenings);
  const Real t1 = Timer::wallClock();
  if (procID() == 0)
    std::cout << "geotime = " << t1 - t0 << std::endl;

  // Set implicit functions now.
  m_amr->setBaseImplicitFunction(phase::gas, m_computationalGeometry->getGasImplicitFunction());
  m_amr->setBaseImplicitFunction(phase::solid, m_computationalGeometry->getSolidImplicitFunction());

  if (m_writeMemory) {
    this->writeMemoryUsage();
  }

  this->getGeometryTags(); // Get geometric tags.

  if (m_writeMemory) {
    this->writeMemoryUsage();
  }

  // Regrid using geometric tags only.
  m_amr->regridAmr(m_geomTags, 0);

  if (m_verbosity > 0) {
    this->gridReport();
  }

  //  this->writeMemoryUsage();
  if (m_plotInterval > 0) {
    this->writeGeometry(); // Write geometry only
  }
}

void
Driver::setupFresh(const int a_initialRegrids)
{
  CH_TIME("Driver::setupFresh");
  if (m_verbosity > 5) {
    pout() << "Driver::setupFresh" << endl;
  }

  this->sanityCheck(); // Sanity check before doing anything expensive

  // Need to specify geometry generation method.
  if (m_geometryGeneration == "chombo-discharge") {
    // We will run with ScanShop, which builds the EBIndexSpace map using recursive pruning of
    // regions that don't contain cut-cells. The user asks for a coarsening/refinement of the
    // base AMR level, which we make here.
    ProblemDomain scanDomain;
    if (m_geoScanLevel < 0) {
      // Make a coarsening of the base AMR level, but don't coarsen beyond 2x2 domains.
      scanDomain  = m_amr->getDomains()[0];
      int numCoar = 0;
      while (scanDomain.domainBox().shortside() >= 4 && numCoar < std::abs(m_geoScanLevel)) {
        scanDomain.coarsen(2);
        numCoar++;
      }
    }
    else {
      const int amrLevel = std::min(m_geoScanLevel, m_amr->getMaxAmrDepth());
      scanDomain         = m_amr->getDomains()[amrLevel];
    }

    EBISLevel::s_distributedData = true;

    m_computationalGeometry->useScanShop(scanDomain);
  }
  else if (m_geometryGeneration == "chombo") {
    if (m_ebisMemoryLoadBalance) {
      EBIndexSpace::s_useMemoryLoadBalance = true;
    }
    else {
      EBIndexSpace::s_useMemoryLoadBalance = false;
    }
    m_computationalGeometry->useChomboShop();
  }

  const int numCoarsenings = m_doCoarsening ? -1 : m_amr->getMaxAmrDepth();
  m_computationalGeometry->buildGeometries(m_amr->getFinestDomain(),
                                           m_amr->getProbLo(),
                                           m_amr->getFinestDx(),
                                           m_amr->getMaxEbisBoxSize(),
                                           m_amr->getNumberOfEbGhostCells(),
                                           numCoarsenings);

  // Register Realms
  m_timeStepper->setAmr(m_amr);
  m_timeStepper->registerRealms();

  // Set implicit functions now.
  m_amr->setBaseImplicitFunction(phase::gas, m_computationalGeometry->getGasImplicitFunction());
  m_amr->setBaseImplicitFunction(phase::solid, m_computationalGeometry->getSolidImplicitFunction());

  // Get geometry tags
  this->getGeometryTags();

  // TLDR:
  // -----
  // The stuff below here might seem a bit convoluted, but we are first letting AmrMesh compute a set of grids which
  // is load balanced by the patch volume. All realms and operators are set up with this initial set of grids. We then
  // use those grids to let the time stepper predict computational loads, which we use to regrid both the realms and the time
  // stepper.

  // When we're setting up fresh, we need to regrid everything from the
  // base level and upwards, so no hardcap on the permitted grids.
  const int lmin    = 0;
  const int hardcap = -1;
  m_amr->regridAmr(m_geomTags, lmin, hardcap);
  const int lmax = m_amr->getFinestLevel();

  // Allocate internal storage
  this->allocateInternals();

  // Provide TimeStepper with geometry in case it needs it.
  m_timeStepper->setComputationalGeometry(m_computationalGeometry);

  // TimeStepper setup. This instantiatse solvers (but does not necessarily fill them with data).
  m_timeStepper->setupSolvers();
  m_timeStepper->synchronizeSolverTimes(m_timeStep, m_time, m_dt);

  // Set up the AMR operators
  m_timeStepper->registerOperators();
  m_amr->regridOperators(lmin);

  // Fill solves with initial data
  m_timeStepper->allocate();
  m_timeStepper->initialData();

  // If called for -- we can perform a
  if (m_doInitLoadBalancing) {
    this->cacheTags(m_tags);
    m_timeStepper->preRegrid(lmin, lmax);
    if (!(m_cellTagger.isNull())) {
      m_cellTagger->preRegrid();
    }
    for (const auto& str : m_amr->getRealms()) {
      if (m_timeStepper->loadBalanceThisRealm(str)) {

        Vector<Vector<int>> procs;
        Vector<Vector<Box>> boxes;

        const int lmin = 0;
        const int lmax = m_amr->getFinestLevel();

        m_timeStepper->loadBalanceBoxes(procs, boxes, str, m_amr->getProxyGrids(), lmin, lmax);

        m_amr->regridRealm(str, procs, boxes, lmin);
      }
    }
    m_amr->regridOperators(lmin);            // Regrid operators again.
    this->regridInternals(lmax, lmax);       // Regrid internals for Driver.
    m_timeStepper->regrid(lmin, lmax, lmax); // Regrid solvers.
    m_timeStepper->initialData();            // Need to fill with initial data again.
  }

  // CellTagger
  if (!m_cellTagger.isNull()) {
    m_cellTagger->regrid();
  }

  // Do a grid report of the initial grid
  if (m_verbosity > 0) {
    this->gridReport();
  }

  // Initial regrids
  for (int i = 0; i < a_initialRegrids; i++) {
    if (m_verbosity > 5) {
      pout() << "Driver -- initial regrid # " << i + 1 << endl;
    }

    const int lmin = 0;
    const int lmax = m_amr->getFinestLevel() + 1;
    this->regrid(lmin, lmax, true);

    if (m_verbosity > 0) {
      this->gridReport();
    }
  }

  // Do post initialize stuff
  m_timeStepper->postInitialize();
}

#ifdef CH_USE_HDF5
void
Driver::setupForRestart(const int a_initialRegrids, const std::string a_restartFile)
{
  CH_TIME("Driver::setupForRestart");
  if (m_verbosity > 5) {
    pout() << "Driver::setupForRestart" << endl;
  }

  this->checkRestartFile(a_restartFile);

  this->sanityCheck(); // Sanity check before doing anything expensive

  // Need to activate some flags that trigger Chombo or chombo-discharge geo-generation method.
  if (m_geometryGeneration == "chombo-discharge") {
    EBISLevel::s_distributedData = true;
    m_computationalGeometry->useScanShop(m_amr->getDomains()[m_geoScanLevel]);
  }
  else if (m_geometryGeneration == "chombo") {
    m_computationalGeometry->useChomboShop();

    if (m_ebisMemoryLoadBalance) {
      EBIndexSpace::s_useMemoryLoadBalance = true;
    }
    else {
      EBIndexSpace::s_useMemoryLoadBalance = false;
    }
  }

  const int numCoarsenings = m_doCoarsening ? -1 : m_amr->getMaxAmrDepth();

  m_computationalGeometry->buildGeometries(m_amr->getFinestDomain(),
                                           m_amr->getProbLo(),
                                           m_amr->getFinestDx(),
                                           m_amr->getMaxEbisBoxSize(),
                                           m_amr->getNumberOfEbGhostCells(),
                                           numCoarsenings);

  this->getGeometryTags(); // Get geometric tags.

  m_timeStepper->setAmr(m_amr);                                     // Set amr
  m_timeStepper->registerRealms();                                  // Register Realms
  m_timeStepper->setComputationalGeometry(m_computationalGeometry); // Set computational geometry

  // Set implicit functions now.
  m_amr->setBaseImplicitFunction(phase::gas, m_computationalGeometry->getGasImplicitFunction());
  m_amr->setBaseImplicitFunction(phase::solid, m_computationalGeometry->getSolidImplicitFunction());

  // Read checkpoint file
  this->readCheckpointFile(a_restartFile);

  // Time stepper does post checkpoint setup
  m_timeStepper->postCheckpointSetup();

  // Prepare storage for CellTagger
  if (!m_cellTagger.isNull()) {
    m_cellTagger->regrid();
  }

  if (m_writeRestartFiles) {
    this->writeRestartFile();
  }

  // Initial regrids
  for (int i = 0; i < a_initialRegrids; i++) {
    if (m_verbosity > 0) {
      pout() << "Driver -- initial regrid # " << i + 1 << endl;
    }

    if (m_writeRegridFiles) {
      this->writePreRegridFile();
    }

    const int lmin = 0;
    const int lmax = m_amr->getFinestLevel() + 1;

    this->regrid(lmin, lmax, false);

    if (m_writeRegridFiles) {
      this->writePostRegridFile();
    }

    if (m_verbosity > 0) {
      this->gridReport();
    }
  }
}
#endif

void
Driver::checkRestartFile(const std::string a_restartFile) const
{
  CH_TIME("Driver::checkRestartFile");
  if (m_verbosity > 4) {
    pout() << "Driver::checkRestartFile" << endl;
  }

  ifstream f(a_restartFile.c_str());
  if (!f.good()) {
    pout() << "Driver::checkRestartFile - could not find file = " << a_restartFile << endl;
    MayDay::Abort("Driver::checkRestartFile - abort, could not find file");
  }
}

void
Driver::sanityCheck()
{
  CH_TIME("Driver::sanityCheck");
  if (m_verbosity > 4) {
    pout() << "Driver::sanityCheck" << endl;
  }

  CH_assert(!m_timeStepper.isNull());
}

void
Driver::stepReport(const Real a_startTime, const Real a_endTime, const int a_maxSteps)
{
  CH_TIME("Driver::stepReport");
  if (m_verbosity > 5) {
    pout() << "Driver::stepReport" << endl;
  }

  pout() << endl;
  pout() << "Driver::Time step report -- Time step #" << m_timeStep << endl
         << "                                   Time  = " << m_time << endl
         << "                                   dt    = " << m_dt << endl;

  m_timeStepper->printStepReport();

  char metrics[300];

  // Percentage completed of time steps
  const Real percentStep = (1.0 * m_timeStep / a_maxSteps) * 100.;
  sprintf(metrics, "%31c -- %5.2f percentage of time steps completed", ' ', percentStep);
  pout() << metrics << endl;

  const Real percentTime = ((m_time - a_startTime) / (a_endTime - a_startTime)) * 100.;
  sprintf(metrics, "%31c -- %5.2f percentage of simulation time completed", ' ', percentTime);
  pout() << metrics << endl;

  // Hours, minutes, seconds and millisecond of the previous iteration
  const Real elapsed    = m_wallClockTwo - m_wallClockStart;
  const int  elapsedHrs = floor(elapsed / 3600);
  const int  elapsedMin = floor((elapsed - 3600 * elapsedHrs) / 60);
  const int  elapsedSec = floor(elapsed - 3600 * elapsedHrs - 60 * elapsedMin);
  const int  elapsedMs  = floor((elapsed - 3600 * elapsedHrs - 60 * elapsedMin - elapsedSec) * 1000);

  // Write a string with total elapsed time
  sprintf(metrics,
          "%31c -- Elapsed time          : %3.3ih %2.2im %2.2is %3.3ims",
          ' ',
          elapsedHrs,
          elapsedMin,
          elapsedSec,
          elapsedMs);
  pout() << metrics << endl;

  // Hours, minutes, seconds and millisecond of the previous iteration
  const Real lastadv = m_wallClockTwo - m_wallClockOne;
  const int  advHrs  = floor(lastadv / 3600);
  const int  advMin  = floor((lastadv - 3600 * advHrs) / 60);
  const int  advSec  = floor(lastadv - 3600 * advHrs - 60 * advMin);
  const int  advMs   = floor((lastadv - 3600 * advHrs - 60 * advMin - advSec) * 1000);

  // Write a string with the previous iteration metrics
  sprintf(metrics, "%31c -- Last time step        : %3.3ih %2.2im %2.2is %3.3ims", ' ', advHrs, advMin, advSec, advMs);
  pout() << metrics << endl;

  // Hours, minutes, seconds and millisecond of the previous iteration
  const Real wt_factor = std::pow(10.0, 3.0 * (int(std::floor(log10(m_dt))) / 3));
  const Real wt_ns     = (m_wallClockTwo - m_wallClockOne) * wt_factor / m_dt;
  const int  wt_Hrs    = floor(wt_ns / 3600);
  const int  wt_Min    = floor((wt_ns - 3600 * wt_Hrs) / 60);
  const int  wt_Sec    = floor(wt_ns - 3600 * wt_Hrs - 60 * wt_Min);
  const int  wt_Ms     = floor((wt_ns - 3600 * wt_Hrs - 60 * wt_Min - wt_Sec) * 1000);
  sprintf(metrics, "%31c -- Wall time per ns      : %3.3ih %2.2im %2.2is %3.3ims", ' ', wt_Hrs, wt_Min, wt_Sec, wt_Ms);
  sprintf(metrics,
          "%31c -- Wall time per/%1.0E   : %3.3ih %2.2im %2.2is %3.3ims",
          ' ',
          wt_factor,
          wt_Hrs,
          wt_Min,
          wt_Sec,
          wt_Ms);
  pout() << metrics << endl;

  // This is the time remaining
  const Real maxPercent = Max(percentTime, percentStep);
  const Real remaining  = 100. * elapsed / maxPercent - elapsed;
  const int  remHrs     = floor(remaining / 3600);
  const int  remMin     = floor((remaining - 3600 * remHrs) / 60);
  const int  remSec     = floor(remaining - 3600 * remHrs - 60 * remMin);
  const int  remMs      = floor((remaining - 3600 * remHrs - 60 * remMin - remSec) * 1000);

  // Write a string with the previous iteration metrics
  sprintf(metrics, "%31c -- Estimated remaining   : %3.3ih %2.2im %2.2is %3.3ims", ' ', remHrs, remMin, remSec, remMs);
  pout() << metrics << endl;

  // Write memory usage
#ifdef CH_USE_MEMORY_TRACKING
  const Real bytesPerMB = 1024. * 1024.;

  long long unfreedMem = 0LL;
  long long peakMem    = 0LL;

  //  overallMemoryUsage(unfreedMem, peakMem);

  pout() << "                                -- Unfreed memory        : " << std::ceil(unfreedMem / bytesPerMB)
         << "(MB)" << endl;
  pout() << "                                -- Peak memory usage     : " << std::ceil(peakMem / bytesPerMB) << "(MB)"
         << endl;

#ifdef CH_MPI
  const long long maxUnfreedMem = ParallelOps::max(unfreedMem);
  const long long maxPeakMem    = ParallelOps::max(peakMem);

  pout() << "                                -- Max unfreed memory    : " << std::ceil(maxUnfreedMem / bytesPerMB)
         << "(MB)" << endl;
  pout() << "                                -- Max peak memory usage : " << std::ceil(maxPeakMem / bytesPerMB)
         << "(MB)" << endl;
#endif
#endif
}

int
Driver::getFinestTagLevel(const EBAMRTags& a_cellTags) const
{
  CH_TIME("Driver::getFinestTagLevel");
  if (m_verbosity > 5) {
    pout() << "Driver::getFinestTagLevel" << endl;
  }

  int finestTagLevel = -1;
  for (int lvl = 0; lvl < a_cellTags.size(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime) reduction(max : finestTagLevel)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const DenseIntVectSet& tags = (*a_cellTags[lvl])[din];

      if (!tags.isEmpty()) {
        finestTagLevel = std::max(finestTagLevel, lvl);
      }
    }
  }

  return ParallelOps::max(finestTagLevel);
}

bool
Driver::tagCells(Vector<IntVectSet>& a_allTags, EBAMRTags& a_cellTags)
{
  CH_TIME("Driver::tagCells");
  if (m_verbosity > 5) {
    pout() << "Driver::tagCells" << endl;
  }

  // TLDR: This routine collects tags from both the cell tagger (which is user-provided) and from geometric tags in
  //       the Driver class. These are aggregated on an IntVectSet on each level.

  bool gotNewTags = false;

  // Note that when we regrid we add at most one level at a time. This means that if we have a
  // simulation with AMR depth l and we want to add a level l+1, we need tags on levels 0 through l.
  const int finestLevel = m_amr->getFinestLevel();
  a_allTags.resize(1 + finestLevel, IntVectSet());

  if (!m_cellTagger.isNull()) {
    gotNewTags = m_cellTagger->tagCells(a_cellTags);
  }

  // Gather tags from the cell tagger.
  for (int lvl = 0; lvl <= finestLevel; lvl++) {
    const DataIterator& dit  = a_cellTags[lvl]->dataIterator();
    const int           nbox = dit.size();

    IntVectSet& tags = a_allTags[lvl];

#pragma omp parallel for schedule(runtime) reduction(+ : tags)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      // Implicity converts to TreeIntVecSet
      tags |= IntVectSet((*a_cellTags[lvl])[din]);
    }

    // Grow tags with cell taggers buffer
    if (!m_cellTagger.isNull()) {
      const int buf = m_cellTagger->getBuffer();

      tags.grow(buf);
    }
  }

  // Add geometric tags.
  if (m_allowCoarsening) {
    const int finestTagLevel = this->getFinestTagLevel(a_cellTags);

    for (int lvl = 0; lvl <= finestLevel; lvl++) {
      if (lvl <= finestTagLevel) {
        a_allTags[lvl] |= m_geomTags[lvl];
      }
    }
  }
  else {
    // Loop only goes to the current finest level because we only add one level at a time
    for (int lvl = 0; lvl <= finestLevel; lvl++) {
      if (lvl < m_amr->getMaxAmrDepth()) { // Geometric tags don't exist on AmrMesh.m_maxAmrDepth
        a_allTags[lvl] |= m_geomTags[lvl];
      }
    }
  }

  return gotNewTags;
}

void
Driver::writeMemoryUsage()
{
  CH_TIME("Driver::writeMemoryUsage()");
  if (m_verbosity > 3) {
    pout() << "Driver::writeMemoryUsage()" << endl;
  }

#ifdef CH_USE_MEMORY_TRACKING
  char              file_char[1000];
  const std::string prefix = m_outputDirectory + "/mpi/memory/" + m_outputFileNames;
  sprintf(file_char, "%s.memory.step%07d.%dd.dat", prefix.c_str(), m_timeStep, SpaceDim);
  std::string fname(file_char);

  // Get memory stuff
  Vector<Real> peakMemory;
  Vector<Real> unfreedMemory;
  MemoryReport::getMemoryUsage(peakMemory, unfreedMemory);

  // Begin writing output
  if (procID() == 0) {
    std::ofstream f;
    f.open(fname, std::ios_base::trunc);
    const int width = 12;

    // Write header
    f << std::left << std::setw(width) << "# MPI rank" << "\t" << std::left << std::setw(width) << "Peak memory" << "\t"
      << std::left << std::setw(width) << "Unfreed memory" << "\t" << endl;

    // Write memory
    for (int i = 0; i < numProc(); i++) {
      f << std::left << std::setw(width) << i << "\t" << std::left << std::setw(width) << peakMemory[i] << "\t"
        << std::left << std::setw(width) << unfreedMemory[i] << "\t" << endl;
    }
  }
#endif
}

void
Driver::writeComputationalLoads()
{
  CH_TIME("Driver::writeComputationalLoads()");
  if (m_verbosity > 3) {
    pout() << "Driver::writeComputationalLoads()" << endl;
  }

  // TLDR: This routine writes loads to files. The file format is in the form
  //
  // rank realm1_loads realm2_loads
  // 0       X             Y
  // 1       XX            YY

  const int nProc = numProc();

  // Filename for output.
  char              file_char[1000];
  const std::string prefix = m_outputDirectory + "/mpi/loads/" + m_outputFileNames;
  sprintf(file_char, "%s.loads.step%07d.%dd.dat", prefix.c_str(), m_timeStep, SpaceDim);
  std::string fname(file_char);

  // Get sum of all loads on all realms
  std::map<std::string, Vector<long int>> realmLoads;
  for (const auto& r : m_amr->getRealms()) {

    // Compute total loads on each rank. This does a call to m_timeStepper to fetch the loads
    // on each realm, which we put in sumLoads, the loads for realm 'r'. Note that each rank
    // only fills the sum of the loads for the patches the it own, so we need to accumulate the result.
    Vector<long int> sumLoads(nProc, 0L);
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const Vector<long int> boxLoads = m_timeStepper->getCheckpointLoads(r, lvl);

      const DisjointBoxLayout& dbl  = m_amr->getGrids(r)[lvl];
      const DataIterator&      dit  = dbl.dataIterator();
      const int                nbox = dit.size();

      long int& load = sumLoads[procID()];

#pragma omp parallel for schedule(runtime) reduction(+ : load)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din     = dit[mybox];
        const int        intCode = din.intCode();

        load += boxLoads[intCode];
      }
    }

    ParallelOps::vectorSum(sumLoads);

    realmLoads.emplace(r, sumLoads);
  }

  // When we're here we have a global view of all the loads for all the realms.
#ifdef CH_MPI
  if (procID() == 0) {
#endif
    const int width = 16;

    std::ofstream f;
    f.open(fname, std::ios_base::trunc);

    // Write header
    std::stringstream ss;
    ss << std::left << std::setw(width) << "# Rank";
    for (auto r : realmLoads) {
      ss << std::left << std::setw(width) << r.first;
    }
    f << ss.str() << endl;

    // Write data.
    for (int irank = 0; irank < nProc; irank++) {
      std::stringstream ds;

      ds << std::left << std::setw(width) << irank;
      for (auto r : realmLoads) {
        ds << std::left << std::setw(width) << r.second[irank];
      }
      f << ds.str() << std::endl;
    }

    f.close();
#ifdef CH_MPI
  }
#endif
}

void
Driver::writeGeometry()
{
  CH_TIME("Driver::writeGeometry()");
  if (m_verbosity > 3) {
    pout() << "Driver::writeGeometry()" << endl;
  }

  // This is a special routine that writes a plot file containing only the level-set function. Very useful when adjusting the
  // geometry.

  const int ncomp = 2;

  EBAMRCellData output;
  m_amr->allocate(output, m_realm, phase::gas, ncomp);
  DataOps::setValue(output, 0.0);

  // Names
  Vector<std::string> names(2);
  names[0] = "levelset_1";
  names[1] = "levelset_2";

  // Write levelsets
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    int icomp = 0;
    this->writeLevelset(*output[lvl], icomp, lvl);
  }

  // Use raw pointers because Chombo IO is not too smart.
  Vector<LevelData<EBCellFAB>*> outputPtr(1 + m_amr->getFinestLevel());
  m_amr->alias(outputPtr, output);

  // Dummy file name
  char              file_char[1000];
  const std::string prefix = m_outputDirectory + "/geo/" + m_outputFileNames;
  sprintf(file_char, "%s.geometry.%dd.hdf5", prefix.c_str(), SpaceDim);
  string fname(file_char);

#ifdef CH_USE_HDF5
  DischargeIO::writeEBHDF5(fname,
                           names,
                           m_amr->getGrids(m_realm),
                           outputPtr,
                           m_amr->getDomains(),
                           m_amr->getDx(),
                           m_amr->getRefinementRatios(),
                           m_dt,
                           m_time,
                           m_amr->getProbLo(),
                           1 + m_amr->getFinestLevel(),
                           m_numPlotGhost);
#endif
}

void
Driver::writePlotFile()
{
  CH_TIME("Driver::writePlotFile()");
  if (m_verbosity > 3) {
    pout() << "Driver::writePlotFile()" << endl;
  }

  // TLDR: This writes a plot file to the /plt/ folder

  // Filename
  char              file_char[1000];
  const std::string prefix = m_outputDirectory + "/plt/" + m_outputFileNames;
  sprintf(file_char, "%s.step%07d.%dd.hdf5", prefix.c_str(), m_timeStep, SpaceDim);
  string fname(file_char);

  // Write.
  this->writePlotFile(fname);
}

void
Driver::writePreRegridFile()
{
  CH_TIME("Driver::writePreRegridFile()");
  if (m_verbosity > 3) {
    pout() << "Driver::writePreRegridFile()" << endl;
  }

  // Filename
  char              file_char[1000];
  const std::string prefix = m_outputDirectory + "/regrid/" + m_outputFileNames;
  sprintf(file_char, "%s.preRegrid%07d.%dd.hdf5", prefix.c_str(), m_timeStep, SpaceDim);
  string fname(file_char);

  this->writePlotFile(fname);
}

void
Driver::writePostRegridFile()
{
  CH_TIME("Driver::writePostRegridFile()");
  if (m_verbosity > 3) {
    pout() << "Driver::writePostRegridFile()" << endl;
  }

  // Filename
  char              file_char[1000];
  const std::string prefix = m_outputDirectory + "/regrid/" + m_outputFileNames;
  sprintf(file_char, "%s.postRegrid%07d.%dd.hdf5", prefix.c_str(), m_timeStep, SpaceDim);
  string fname(file_char);

  this->writePlotFile(fname);
}

void
Driver::writeRestartFile()
{
  CH_TIME("Driver::writeRestartFile()");
  if (m_verbosity > 3) {
    pout() << "Driver::writeRestartFile()" << endl;
  }

  // Filename
  char              file_char[1000];
  const std::string prefix = m_outputDirectory + "/restart/" + m_outputFileNames;
  sprintf(file_char, "%s.restart%07d.%dd.hdf5", prefix.c_str(), m_timeStep, SpaceDim);
  string fname(file_char);

  this->writePlotFile(fname);
}

void
Driver::writeCrashFile()
{
  CH_TIME("Driver::writeCrashFile()");
  if (m_verbosity > 3) {
    pout() << "Driver::writeCrashFile()" << endl;
  }

  // Filename
  char              file_char[1000];
  const std::string prefix = m_outputDirectory + "/crash/" + m_outputFileNames;
  sprintf(file_char, "%s.crash%07d.%dd.hdf5", prefix.c_str(), m_timeStep, SpaceDim);
  string fname(file_char);

  this->writePlotFile(fname);
}

void
Driver::writePlotFile(const std::string a_filename)
{
  CH_TIMERS("Driver::writePlotFile(string)");
  CH_TIMER("Driver::writePlotFile::allocate", t1);
  CH_TIMER("Driver::writePlotFile::assemble", t2);
  CH_TIMER("Driver::writePlotFile::interp_exchange", t3);
  CH_TIMER("Driver::writePlotFile::copy_internal", t4);
  CH_TIMER("Driver::writePlotFile::hdf5_write", t5);

  if (m_verbosity >= 1) {
    pout() << "Driver::writePlotFile(string)" << endl;
  }

  // Time stepper does pre-plot operations
  m_timeStepper->prePlot();
  if (!(m_cellTagger.isNull())) {
    m_cellTagger->prePlot();
  }

  // Get total number of components for output
  int numOutputComp = m_timeStepper->getNumberOfPlotVariables();
  if (!m_cellTagger.isNull()) {
    numOutputComp += m_cellTagger->getNumberOfPlotVariables();
  }
  numOutputComp += this->getNumberOfPlotVariables();

  if (numOutputComp > 0) {
    // Users get to restrict the maximum plot depth.
    const int finestLevel   = m_amr->getFinestLevel();
    const int maxPlotLevel  = (m_maxPlotLevel < 0) ? finestLevel : std::min(m_maxPlotLevel, finestLevel);
    const int numPlotLevels = maxPlotLevel + 1;

    // Get plot variable names.
    Vector<std::string> plotVariableNames;
    plotVariableNames.append(m_timeStepper->getPlotVariableNames());
    if (!(m_cellTagger.isNull())) {
      plotVariableNames.append(m_cellTagger->getPlotVariableNames());
    }
    plotVariableNames.append(this->getPlotVariableNames());

    // Write HDF5 header.
#ifdef CH_USE_HDF5
    HDF5Handle handle(a_filename.c_str(), HDF5Handle::CREATE);
    DischargeIO::writeEBHDF5Header(handle, numPlotLevels, m_amr->getProbLo(), plotVariableNames);
    handle.close();
#endif

    Timer timer("Driver::writePlotFile");
    for (int lvl = 0; lvl <= maxPlotLevel; lvl++) {
      timer.startEvent("Allocate");
      LevelData<EBCellFAB> outputData;
      m_amr->allocate(outputData, m_realm, phase::gas, lvl, numOutputComp);
      DataOps::setValue(outputData, 0.0);
      timer.stopEvent("Allocate");

      // Relevant components collect data for IO.
      int comp = 0;
      timer.startEvent("Assemble data");
      if (m_verbosity > 2) {
        pout() << "Driver::writePlotFile -- assembling data on level = " << lvl << endl;
      }
      m_timeStepper->writePlotData(outputData, comp, m_realm, lvl);
      if (!(m_cellTagger.isNull())) {
        m_cellTagger->writePlotData(outputData, comp, m_realm, lvl);
      }
      this->writePlotData(outputData, comp, lvl);
      timer.stopEvent("Assemble data");

      // Do the HDF5 write.
#ifdef CH_USE_HDF5
      if (m_verbosity > 2) {
        pout() << "Driver::writePlotFile -- HDF5 write on level = " << lvl << endl;
        MemoryReport::getMaxMinMemoryUsage();
      }

      timer.startEvent("HDF5 write");
      HDF5Handle handle(a_filename.c_str(), HDF5Handle::OPEN_RDWR);
      const int  refRat = (lvl < m_amr->getFinestLevel()) ? m_amr->getRefinementRatios()[lvl] : 1;
      DischargeIO::writeEBHDF5Level(handle,
                                    outputData,
                                    m_amr->getDomains()[lvl],
                                    m_amr->getDx()[lvl],
                                    m_dt,
                                    m_time,
                                    lvl,
                                    refRat,
                                    m_numPlotGhost);
      handle.close();
      timer.stopEvent("HDF5 write");

      if (m_verbosity > 2) {
        MemoryReport::getMaxMinMemoryUsage();
      }
#endif
    }

    if (m_profile) {
      timer.eventReport(pout(), true);
    }
  }
  else {
    const std::string msg1 = "Driver::writePlotFile - skipping file '";
    const std::string msg2 = "' since there is nothing to write";

    pout() << msg1 << a_filename.c_str() << msg2 << endl;
  }

  // TimeStepper does post-plot operations.
  m_timeStepper->postPlot();
}

void
Driver::writePlotData(LevelData<EBCellFAB>& a_output, int& a_comp, const int a_level) const noexcept
{
  CH_TIME("Driver::writePlotData(LevelData<EBCellFAB>, int, int)");
  if (m_verbosity > 3) {
    pout() << "Driver::writePlotData(LevelData<EBCellFAB>, int, int)" << endl;
  }

  CH_assert(a_level >= 0);
  CH_assert(a_level <= m_amr->getFinestLevel());

  if (m_plotTags) {
    this->writeTags(a_output, a_comp, a_level);
  }
  if (m_plotRanks) {
    this->writeRanks(a_output, a_comp, a_level);
  }
  if (m_plotLevelset) {
    this->writeLevelset(a_output, a_comp, a_level);
  }
  if (m_plotLoads) {
    this->writeLoads(a_output, a_comp, a_level);
  }
}

void
Driver::writeTags(LevelData<EBCellFAB>& a_output, int& a_comp, const int a_level) const noexcept
{
  CH_TIME("Driver::writeTags");
  if (m_verbosity > 3) {
    pout() << "Driver::writeTags" << endl;
  }

  CH_assert(a_level >= 0);
  CH_assert(a_level <= m_amr->getFinestLevel());
  CH_assert(a_output.nComp() > a_comp);

  // Need some scratch storage for the tags because DenseIntVectSet can't be written to file. We just run through
  // the cells and set data == 1 if a cell had been flagged for refinement.
  LevelData<EBCellFAB> scratch;
  m_amr->allocate(scratch, m_realm, phase::gas, a_level, 1);
  DataOps::setValue(scratch, 0.0);

  // Set tagged cells = 1
  const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[a_level];
  const DataIterator&      dit  = dbl.dataIterator();
  const int                nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    FArrayBox&             regTags = scratch[din].getFArrayBox();
    const DenseIntVectSet& ivs     = (*m_tags[a_level])[din];

    auto kernel = [&](const IntVect& iv) -> void {
      if (ivs[iv]) {
        regTags(iv, 0) = 1.0;
      }
    };

    BoxLoops::loop(dbl[din], kernel);
  }

  // Copy 'tags' over to 'a_output', starting on component a_comp.
  const Interval srcInterv(0, 0);
  const Interval dstInterv(a_comp, a_comp);

  m_amr->copyData(scratch, a_output, a_level, m_realm, m_realm, dstInterv, srcInterv);

  a_comp++;
}

void
Driver::writeRanks(LevelData<EBCellFAB>& a_output, int& a_comp, const int a_level) const noexcept
{
  CH_TIME("Driver::writeRanks");
  if (m_verbosity > 3) {
    pout() << "Driver::writeRanks" << endl;
  }

  CH_assert(a_level >= 0);
  CH_assert(a_level <= m_amr->getFinestLevel());

  for (const auto& r : m_amr->getRealms()) {
    CH_assert(a_output.nComp() > a_comp);

    LevelData<EBCellFAB> scratch;
    m_amr->allocate(scratch, r, phase::gas, a_level, 1);

    DataOps::setValue(scratch, 1.0 * procID());

    const Interval srcInterv(0, 0);
    const Interval dstInterv(a_comp, a_comp);

    m_amr->copyData(a_output, scratch, a_level, m_realm, r, dstInterv, srcInterv);

    a_comp++;
  }
}

void
Driver::writeLevelset(LevelData<EBCellFAB>& a_output, int& a_comp, const int a_level) const noexcept
{
  CH_TIME("Driver::writeLevelset");
  if (m_verbosity > 3) {
    pout() << "Driver::writeLevelset" << endl;
  }

  CH_assert(a_level >= 0);
  CH_assert(a_level <= m_amr->getFinestLevel());
  CH_assert(a_output.nComp() > a_comp + 1);

  const RefCountedPtr<BaseIF>& lsf1 = m_computationalGeometry->getGasImplicitFunction();
  const RefCountedPtr<BaseIF>& lsf2 = m_computationalGeometry->getSolidImplicitFunction();

  const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[a_level];
  const Real               dx     = m_amr->getDx()[a_level];
  const RealVect           probLo = m_amr->getProbLo();
  const DataIterator&      dit    = dbl.dataIterator();

  const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    FArrayBox& fab = a_output[din].getFArrayBox();

    fab.setVal(std::numeric_limits<Real>::max(), a_comp);
    fab.setVal(std::numeric_limits<Real>::max(), a_comp + 1);

    auto kernel = [&](const IntVect& iv) -> void {
      const RealVect pos = probLo + (RealVect(iv) + 0.5 * RealVect::Unit) * dx;

      if (!lsf1.isNull()) {
        fab(iv, a_comp) = lsf1->value(pos);
      }
      if (!lsf2.isNull()) {
        fab(iv, a_comp + 1) = lsf2->value(pos);
      }
    };

    BoxLoops::loop(fab.box(), kernel);
  }

  a_comp = a_comp + 2;
}

void
Driver::writeLoads(LevelData<EBCellFAB>& a_output, int& a_comp, const int a_level) const noexcept
{
  CH_TIME("Driver::writeLoads");
  if (m_verbosity > 3) {
    pout() << "Driver::writeLoads" << endl;
  }

  CH_assert(a_level >= 0);
  CH_assert(a_level <= m_amr->getFinestLevel());

  for (const auto& r : m_amr->getRealms()) {
    CH_assert(a_output.nComp() > a_comp);

    LevelData<EBCellFAB> scratch;
    m_amr->allocate(scratch, r, phase::gas, a_level, 1);

    const Vector<long int> loads = m_timeStepper->getCheckpointLoads(r, a_level);

    const DisjointBoxLayout& dbl = m_amr->getGrids(r)[a_level];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      scratch[din].setVal(loads[din.intCode()]);
    }

    const Interval srcInterv(0, 0);
    const Interval dstInterv(a_comp, a_comp);

    m_amr->copyData(a_output, scratch, a_level, m_realm, r, dstInterv, srcInterv);

    a_comp++;
  }
}

void
Driver::writeCheckpointFile()
{
  CH_TIME("Driver::writeCheckpointFile()");
  if (m_verbosity >= 1) {
    pout() << "Driver::writeCheckpointFile()" << endl;
  }

#ifdef CH_USE_HDF5
  const int finestLevel      = m_amr->getFinestLevel();
  int       finestCheckLevel = Min(m_maxCheckpointDepth, finestLevel);
  if (m_maxCheckpointDepth < 0) {
    finestCheckLevel = finestLevel;
  }

  // Write header containing time step information, grid resolution, etc.
  HDF5HeaderData header;
  header.m_int["finest_eb_lvl"] = m_amr->getMaxAmrDepth();
  header.m_real["coarsest_dx"]  = m_amr->getDx()[0];
  header.m_real["time"]         = m_time;
  header.m_real["dt"]           = m_dt;
  header.m_int["step"]          = m_timeStep;
  header.m_int["finestLevel"]   = finestLevel;

  // Write realm names -- these are needed because we also write computational loads to checkpoint files
  // so we can load balance on immediately restart, using the checkpointed loads.
  for (auto r : m_amr->getRealms()) {
    header.m_string[r] = r;
  }

  // Time stepper writes necessary meta-data to header.
  m_timeStepper->writeCheckpointHeader(header);

  // Create the output file name.
  char              str[100];
  const std::string prefix = m_outputDirectory + "/chk/" + m_outputFileNames;
  sprintf(str, "%s.check%07d.%dd.hdf5", prefix.c_str(), m_timeStep, SpaceDim);

  // Output file
  HDF5Handle handleOut(str, HDF5Handle::CREATE);
  header.writeToFile(handleOut);

  Timer timer("Driver::writeCheckpointFile");
  if (m_verbosity >= 3) {
    pout() << "Driver::writeCheckpointFile - writing checkpoint file..." << endl;
  }
  if (m_profile) {
    timer.startEvent("Write data");
  }

  for (int lvl = 0; lvl <= finestCheckLevel; lvl++) {
    if (m_verbosity >= 4) {
      pout() << "Driver::writeCheckpointFile -- writing level = " << lvl << endl;
    }

    handleOut.setGroupToLevel(lvl);

    // Write amr grids
    write(handleOut, m_amr->getGrids(m_realm)[lvl]); // write AMR grids

    // Time stepper checkpoints solver data.
    m_timeStepper->writeCheckpointData(handleOut, lvl);

    // Driver checkpoints internal data. This involves the cell tags and the computational loads on the various realms.
    this->writeCheckpointLevel(handleOut, lvl);
  }

  if (m_profile) {
    timer.stopEvent("Write data");
    timer.eventReport(pout(), true);
  }

  handleOut.close();
#endif
}

#ifdef CH_USE_HDF5
void
Driver::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level)
{
  CH_TIME("Driver::writeCheckpointLevel(HDF5Handle, int)");
  if (m_verbosity > 5) {
    pout() << "Driver::writeCheckpointLevel(HDF5Handle, int)" << endl;
  }

  this->writeCheckpointTags(a_handle, a_level);
  this->writeCheckpointRealmLoads(a_handle, a_level);
}
#endif

#ifdef CH_USE_HDF5
void
Driver::writeCheckpointTags(HDF5Handle& a_handle, const int a_level)
{
  CH_TIME("Driver::writeCheckpointTags(HDF5Handle, int)");
  if (m_verbosity > 5) {
    pout() << "Driver::writeCheckpointTags(HDF5Handle, int)" << endl;
  }

  // We have a set of tags in DenseIntVect that we want to checkpoint. Chombo doesn't do I/O of IntVect
  // so we create some mesh data on level a_lvl which we fill with grid tags. Then just checkpoint that file instead.

  EBCellFactory        fact(m_amr->getEBISLayout(m_realm, phase::gas)[a_level]);
  LevelData<EBCellFAB> scratch(m_amr->getGrids(m_realm)[a_level], 1, 3 * IntVect::Unit, fact);
  DataOps::setValue(scratch, 0.0);

  // Set tags = 1
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_level];
  const DataIterator&      dit = dbl.dataIterator();

  const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    FArrayBox&             fab  = scratch[din].getFArrayBox();
    const DenseIntVectSet& tags = (*m_tags[a_level])[din];

    auto kernel = [&](const IntVect& iv) -> void {
      if (tags[iv]) {
        fab(iv, 0) = 1.0;
      }
    };

    BoxLoops::loop(dbl[din], kernel);

    DataOps::setCoveredValue(scratch, 0, 0.0);
  }

  // Write tags
  write(a_handle, scratch, "tagged_cells");
}
#endif

#ifdef CH_USE_HDF5
void
Driver::writeCheckpointRealmLoads(HDF5Handle& a_handle, const int a_level)
{
  CH_TIME("Driver::writeCheckpointRealmLoads(HDF5Handle, int)");
  if (m_verbosity > 5) {
    pout() << "Driver::writeCheckpointRealmLoads(HDF5Handle, int)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_level];
  const DataIterator&      dit = dbl.dataIterator();

  const int nbox = dit.size();

  // Make some storage and set the computational load to be the same in every patch. Later when
  // we read the file we can just fetch the computational load from one of the cells.
  LevelData<FArrayBox> scratch(dbl, 1, IntVect::Zero);

  // Get loads.
  for (auto r : m_amr->getRealms()) {
    const Vector<long int> loads = m_timeStepper->getCheckpointLoads(r, a_level);

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      scratch[din].setVal(loads[din.intCode()]);
    }

    // String identifier in HDF file.
    const std::string str = r + "_loads";

    // Write
    write(a_handle, scratch, str);
  }
}
#endif

#ifdef CH_USE_HDF5
void
Driver::readCheckpointFile(const std::string& a_restartFile)
{
  CH_TIME("Driver::readCheckpointFile(string)");
  if (m_verbosity > 0) {
    pout() << "Driver::readCheckpointFile(string)" << endl;
  }

  // Time stepper can register realms immediately.
  m_timeStepper->registerRealms();

  // Read the header that was written by new_readCheckpointFile
  HDF5Handle     handle_in(a_restartFile, HDF5Handle::OPEN_RDONLY);
  HDF5HeaderData header;
  header.readFromFile(handle_in);

  m_time     = header.m_real["time"];
  m_dt       = header.m_real["dt"];
  m_timeStep = header.m_int["step"];

  const Real coarsestDx  = header.m_real["coarsest_dx"];
  const int  baseLevel   = 0;
  const int  finestLevel = header.m_int["finestLevel"];
  const int  prevMaxEB   = header.m_int["finest_eb_lvl"];

  if (prevMaxEB != m_amr->getMaxAmrDepth()) {
    const std::string err = "Driver::readCheckpointFile -- max EB level changed. Proceed at your own peril";

    MayDay::Warning(err.c_str());
  }

  // Get the names of the realms that were checkpointed. This is a part of the HDF header.
  std::map<std::string, Vector<Vector<long int>>> checkpointedLoads;
  for (auto s : header.m_string) {
    checkpointedLoads.emplace(s.second, Vector<Vector<long int>>());
  }

  // Get then names of the realms that will be used for simulations. These are not necessarily the same because users
  // can restart with a different number of realms. If the user uses a new realm, the computational load is taken as
  // the primal realm load, i.e. volume-based.
  std::map<std::string, Vector<Vector<long int>>> simulationLoads;
  for (const auto& iRealm : m_amr->getRealms()) {
    simulationLoads.emplace(iRealm, Vector<Vector<long int>>());
  }

  // Print checkpointed Realm names.
  if (m_verbosity > 2) {
    pout() << "Driver::readCheckpointFile - checked Realms are: ";
    for (auto r : checkpointedLoads) {
      pout() << '"' << r.first << '"' << "\t";
    }
    pout() << endl;
  }

  m_timeStepper->readCheckpointHeader(header);

  // Abort if base resolution has changed.
  if (!(coarsestDx == m_amr->getDx()[0])) {
    MayDay::Error(
      "Driver::readCheckpointFile - coarsestDx != dx[0], did you change the base level resolution when restarting?!?");
  }

  // Read in grids. If the file has no grids then something has gone wrong.
  Vector<Vector<Box>> boxes(1 + finestLevel);
  for (int lvl = 0; lvl <= finestLevel; lvl++) {
    handle_in.setGroupToLevel(lvl);

    const int status = read(handle_in, boxes[lvl]);

    if (status != 0) {
      MayDay::Error("Driver::readCheckpointFile - file has no grids");
    }
  }

  // Read in the computational loads from the HDF5 file.
  for (auto& r : checkpointedLoads) {
    const std::string&        realmName  = r.first;
    Vector<Vector<long int>>& realmLoads = r.second;

    realmLoads.resize(1 + finestLevel);
    for (int lvl = 0; lvl <= finestLevel; lvl++) {
      realmLoads[lvl].resize(boxes[lvl].size(), 0L);

      this->readCheckpointRealmLoads(realmLoads[lvl], handle_in, realmName, lvl);
    }
  }

  // In case we restart with more or fewer realms, we need to decide how to assign the computational loads. If the realm was a new realm we may
  // not have the computational loads for that. in that case we take the computational loads from the primal realm which we know is always there.
  for (auto& s : simulationLoads) {

    const std::string&        curRealm = s.first;
    Vector<Vector<long int>>& curLoads = s.second;

    // For every realm in the restarted simulation, check if there was an equivalent checkpointed realm. If there was,
    // fetch the loads from that realm.
    bool foundCheckedLoads = false;
    for (const auto& c : checkpointedLoads) {
      if (curRealm == c.first) {
        foundCheckedLoads = true;
      }
    }

    curLoads = (foundCheckedLoads) ? checkpointedLoads.at(curRealm) : checkpointedLoads.at(Realm::Primal);
  }

  // Define AmrMesh and Realms.
  m_amr->setFinestLevel(finestLevel);
  m_amr->setGrids(boxes, simulationLoads);

  // Instantiate solvers and register operators
  m_timeStepper->setupSolvers();
  m_timeStepper->registerOperators();
  m_timeStepper->synchronizeSolverTimes(m_timeStep, m_time, m_dt);
  m_amr->regridOperators(baseLevel);
  m_timeStepper->allocate();

  // Allocate internal stuff (e.g. space for tags)
  this->allocateInternals();

  // Go through level by level and have solvers extract their data
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    handle_in.setGroupToLevel(lvl);

    // Time stepper reads in data
    m_timeStepper->readCheckpointData(handle_in, lvl);

    // Read in internal data
    readCheckpointLevel(handle_in, lvl);
  }

  // Close input file
  handle_in.close();

  if (m_verbosity > 0) {
    pout() << "Driver::readCheckpointFile(string) -- DONE!" << endl;
  }
}
#endif

#ifdef CH_USE_HDF5
void
Driver::readCheckpointLevel(HDF5Handle& a_handle, const int a_level)
{
  CH_TIME("Driver::readCheckpointLevel(HDF5Handle, int)");
  if (m_verbosity > 5) {
    pout() << "Driver::readCheckpointLevel(HDF5Handle, int)" << endl;
  }

  const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[a_level];
  const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, phase::gas)[a_level];
  const DataIterator&      dit   = dbl.dataIterator();
  const int                nbox  = dit.size();

  // Some scratch data we can use
  LevelData<EBCellFAB> scratch(dbl, 1, IntVect::Zero, EBCellFactory(ebisl));
  DataOps::setValue(scratch, 0.0);

  // Read in tags
  read<EBCellFAB>(a_handle, scratch, "tagged_cells", dbl, Interval(0, 0), false);

  // Instantiate m_tags. When we wrote the tags we put a floating point value of 0 where we didn't have tags
  // and a floating point value of 1 where we had.
#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    DenseIntVectSet& taggedCells = (*m_tags[a_level])[din];
    FArrayBox&       scratchFAB  = scratch[din].getFArrayBox();

    auto kernel = [&](const IntVect& iv) -> void {
      if (scratchFAB(iv, 0) > 0.9999) {
        taggedCells |= iv;
      }
    };

    BoxLoops::loop(dbl[din], kernel);
  }
}
#endif

#ifdef CH_USE_HDF5
void
Driver::readCheckpointRealmLoads(Vector<long int>& a_loads,
                                 HDF5Handle&       a_handle,
                                 const std::string a_realm,
                                 const int         a_level)
{
  CH_TIME("Driver::readCheckpointRealmLoads(Vector<long int>, HDF5Handle, string, int)");
  if (m_verbosity > 5) {
    pout() << "Driver::readCheckpointRealmLoads((Vector<long int>, HDF5Handle, string, int))" << endl;
  }

  // This reads computational loads from a file. When we wrote them we set the load in each grid patch, so we can just fetch using
  // FArrayBox::max() function for getting it.

  // Identifier in HDF5 file.
  const std::string str = a_realm + "_loads";

#ifdef CH_MPI
  const int                 nBoxes   = a_loads.size();
  const std::pair<int, int> beginEnd = ParallelOps::partition(nBoxes);
  for (int i = 0; i < nBoxes; i++) {
    a_loads[i] = 0L;
  }

  for (int ibox = beginEnd.first; ibox <= beginEnd.second; ibox++) {
    FArrayBox fab;

    readFArrayBox(a_handle, fab, a_level, ibox, Interval(0, 0), str);

    a_loads[ibox] = lround(fab.max());
  }

  ParallelOps::vectorSum(a_loads);
#else
  // Read into an FArrayBox.
  FArrayBox fab;
  for (int ibox = 0; ibox < a_loads.size(); ibox++) {
    readFArrayBox(a_handle, fab, a_level, ibox, Interval(0, 0), str);

    a_loads[ibox] = lround(fab.max());
  }
#endif
}
#endif

#include <CD_NamespaceFooter.H>
