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
#include <CD_VofUtils.H>
#include <CD_DataOps.H>
#include <CD_MultifluidAlias.H>
#include <CD_Units.H>
#include <CD_MemoryReport.H>
#include <CD_NamespaceHeader.H>

Driver::Driver(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry,
	       const RefCountedPtr<TimeStepper>&           a_timeStepper,
	       const RefCountedPtr<AmrMesh>&               a_amr,
	       const RefCountedPtr<CellTagger>&            a_cellTagger,
	       const RefCountedPtr<GeoCoarsener>&          a_geoCoarsen){
  CH_TIME("Driver::Driver(full)");

  setComputationalGeometry(a_computationalGeometry); // Set computational geometry
  setTimeStepper(a_timeStepper);                     // Set time stepper
  setAmr(a_amr);                                      // Set amr
  setCellTagger(a_cellTagger);                       // Set cell tagger
  setGeoCoarsen(a_geoCoarsen);                       // Set geo coarsener

  // AMR does its thing
  m_amr->sanityCheck();                 // Sanity check, make sure everything is set up correctly
  m_amr->buildDomains();                // Build domains and resolutions, nothing else

  // Ok we're ready to go. 
  m_timeStep      = 0;
  m_time          = 0.0;
  m_dt            = 0.0;  

  // Parse some class options
  parseOptions();

  // Create output directories. 
  createOutputDirectories();

  // Always register this Realm and these operators. 
  m_realm = Realm::Primal;
  m_amr->registerRealm(m_realm);
  m_amr->registerOperator(s_eb_fill_patch, m_realm, phase::gas); 
}

Driver::~Driver(){
  CH_TIME("Driver::~Driver");
}

int Driver::getNumberOfPlotVariables() const {
  CH_TIME("Driver::getNumberOfPlotVariables");
  if(m_verbosity > 5){
    pout() << "Driver::getNumberOfPlotVariables" << endl;
  }

  int num_output = 0;

  if(m_plotTags)     num_output = num_output + 1;
  if(m_plotRanks)    {
    const int num_realms = m_amr->getRealms().size();
    num_output = num_output + num_realms;
  }
  if(m_plotLevelset) num_output = num_output + 2;

  return num_output;
}

Vector<std::string> Driver::getPlotVariableNames() const {
  CH_TIME("Driver::getPlotVariableNames");
  if(m_verbosity > 5){
    pout() << "Driver::getPlotVariableNames" << endl;
  }
  
  Vector<std::string> names(0);
  
  if(m_plotTags) names.push_back("cell_tags");
  if(m_plotRanks) {
    const std::string base = "_rank";
    for (const auto& str : m_amr->getRealms()){
      const std::string id = str + base;
      names.push_back(id);
    }
  }
  if(m_plotLevelset){
    names.push_back("levelset_1");
    names.push_back("levelset_2");
  }
  
  return names;
}

void Driver::allocateInternals(){
  CH_TIME("Driver::allocateInternals");
  if(m_verbosity > 5){
    pout() << "Driver::allocateInternals" << endl;
  }
  
  const int ncomp        = 1;
  const int finest_level = m_amr->getFinestLevel();
  const IntVect ghost    = IntVect::Zero;

  m_tags.resize(1 + finest_level);
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    m_tags[lvl] = RefCountedPtr<LayoutData<DenseIntVectSet> > (new LayoutData<DenseIntVectSet>(dbl));

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      DenseIntVectSet& divs = (*m_tags[lvl])[dit()];
      divs = DenseIntVectSet(dbl.get(dit()), false);
    }
  }
}

void Driver::cacheTags(const EBAMRTags& a_tags){
  CH_TIME("Driver::cacheTags");
  if(m_verbosity > 5){
    pout() << "Driver::cacheTags" << endl;
  }

  const int ncomp         = 1;
  const int finest_level  = m_amr->getFinestLevel();
  const int ghost         = 0;

  m_amr->allocate(m_cachedTags, m_realm, ncomp, ghost);
  m_cachedTags.resize(1+finest_level);
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    // Copy tags onto boolean mask
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      BaseFab<bool>& cached_tags  = (*m_cachedTags[lvl])[dit()];
      cached_tags.setVal(false);

      const IntVectSet divs = IntVectSet((*m_tags[lvl])[dit()]);
      for (IVSIterator ivsIt(divs); ivsIt.ok(); ++ivsIt){
	const IntVect iv = ivsIt();
	cached_tags(iv,0) = true;
      }
    }
  }
}

void Driver::deallocateInternals(){
  CH_TIME("Driver::deallocateInternals");
  if(m_verbosity > 5){
    pout() << "Driver::deallocateInternals" << endl;
  }

  //  m_amr->deallocate(m_tags);
}

void Driver::getGeometryTags(){
  CH_TIME("Driver::getGeometryTags");
  if(m_verbosity > 5){
    pout() << "Driver::getGeometryTags" << endl;
  }

  const int maxdepth = m_amr->getMaxAmrDepth();

  m_geomTags.resize(maxdepth);

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  CH_assert(ebis_gas != NULL);

  for (int lvl = 0; lvl < maxdepth; lvl++){ // Don't need tags on maxdepth, we will never generate grids below that.
    const ProblemDomain& cur_dom = m_amr->getDomains()[lvl];
    const int which_level = ebis_gas->getLevel(cur_dom);

    IntVectSet cond_tags;
    IntVectSet diel_tags;
    IntVectSet gas_tags;
    IntVectSet solid_tags;
    IntVectSet gas_diel_tags;
    IntVectSet gas_solid_tags;

    // Conductor cells
    if(m_conductorTagsDepth > lvl){ 
      cond_tags = ebis_gas->irregCells(which_level);
      if(!ebis_sol.isNull()){
	cond_tags |= ebis_sol->irregCells(which_level);
	cond_tags -= m_multifluidIndexSpace->interfaceRegion(cur_dom);
      }
    }

    // dielectric cells
    if(m_dielectricTagsDepth > lvl){ 
      if(!ebis_sol.isNull()){
	diel_tags = ebis_sol->irregCells(which_level);
      }
    }

    // Gas-solid interface cells
    if(m_gasSolidInterfaceTagDepth > lvl){ 
      if(!ebis_sol.isNull()){
    	gas_tags = ebis_gas->irregCells(which_level);
      }
    }

    // Gas-Dielectric interface cells
    if(m_gasDielectricInterfaceTagDepth > lvl){
      if(!ebis_sol.isNull()){
    	gas_diel_tags = m_multifluidIndexSpace->interfaceRegion(cur_dom);
      }
    }

    // Gas-conductor interface cells
    if(m_gasConductorInterfaceTagDepth > lvl){ 
      gas_solid_tags = ebis_gas->irregCells(which_level);
      if(!ebis_sol.isNull()){
    	gas_solid_tags -= m_multifluidIndexSpace->interfaceRegion(cur_dom);
      }
    }

    // Solid-solid interfaces
    if(m_solidSolidInterfaceTagDepth > lvl){ 
      if(!ebis_sol.isNull()){
    	solid_tags = ebis_sol->irregCells(which_level);

    	// Do the intersection with the conductor cells
    	IntVectSet tmp = ebis_gas->irregCells(which_level);
    	tmp |= ebis_sol->irregCells(which_level);
    	tmp -= m_multifluidIndexSpace->interfaceRegion(cur_dom);

    	solid_tags &= tmp;
      }
    }

    m_geomTags[lvl].makeEmpty();

    // Things from depth specifications
    m_geomTags[lvl] |= diel_tags;
    m_geomTags[lvl] |= cond_tags;
    m_geomTags[lvl] |= gas_diel_tags;
    m_geomTags[lvl] |= gas_solid_tags;
    m_geomTags[lvl] |= gas_tags;
    m_geomTags[lvl] |= solid_tags;


    // Evaluate angles between cut-cells and refine based on that. 
    DisjointBoxLayout irregGrids = ebis_gas->getIrregGrids(cur_dom);
    EBISLayout ebisl;
    ebis_gas->fillEBISLayout(ebisl, irregGrids, cur_dom, 1); // Need one ghost cell because we fetch neighboring vofs. 

    const Real dx         = m_amr->getDx()[lvl];
    const RealVect probLo = m_amr->getProbLo();
    
    for (DataIterator dit(irregGrids); dit.ok(); ++dit){
      const Box box = irregGrids[dit()];

      const EBISBox& ebisbox = ebisl[dit()];
      const EBGraph& ebgraph = ebisbox.getEBGraph();
      const IntVectSet irreg = ebisbox.getIrregIVS(box);

      for (VoFIterator vofit(irreg, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof    = vofit();
	const IntVect   iv     = vof.gridIndex();
	const RealVect  normal = ebisbox.normal(vof);
	
	// Check the angle between the normal vector in this irregular cell and neighboring irregular cells. If the
	// angle exceeds a specified threshold the cell is refined. 
	const Vector<VolIndex> otherVofs = VofUtils::getVofsInRadius(vof, ebisbox, 1, VofUtils::Connectivity::MonotonePath, false);

	for (int i = 0; i < otherVofs.size(); i++){
	  const VolIndex& curVof = otherVofs[i];

	  if(ebisbox.isIrregular(curVof.gridIndex())){ // Only check irregular cells
	    constexpr Real degreesPerRadian = 180.0/Units::pi;
	    const RealVect curNormal = ebisbox.normal(curVof);
	    const Real cosAngle      = PolyGeom::dot(normal, curNormal)/(normal.vectorLength()*curNormal.vectorLength());
	    const Real theta         = acos(cosAngle)*degreesPerRadian;

	    // Refine if angle exceeds threshold
	    if(std::abs(theta) > m_refineAngle) m_geomTags[lvl] |= iv;
	  }
	}
      }
    }
  }

  // Grow tags. This is an ad-hoc fix that prevents ugly grid near EBs (i.e. cases where only ghost cells are used
  // for elliptic equations)
  for (int lvl = 0; lvl < maxdepth; lvl++){
    m_geomTags[lvl].grow(m_irregTagGrowth);
  }

  // Remove tags using the geocoarsener if we have it
  if(!m_geoCoarsen.isNull()){
    m_geoCoarsen->coarsenTags(m_geomTags, m_amr->getDx(), m_amr->getProbLo());
  }

#ifdef CH_MPI
  // Processes may not agree what is the maximum tag depth. Make sure they're all on the same page. 
  int deepestTagLevel = 0;
  for (int lvl = 0; lvl < m_geomTags.size(); lvl++){
    if(!m_geomTags[lvl].isEmpty()) deepestTagLevel = lvl;
  }

  int tmp = -1;
  MPI_Allreduce(&deepestTagLevel, &m_geometricTagsDepth, 1, MPI_INT, MPI_MAX, Chombo_MPI::comm);
#endif
}

void Driver::getLoadsAndBoxes(long long& a_myPoints,
			      long long& a_myPointsGhosts,
			      long long& a_myBoxes,
			      long long& a_totalPoints,
			      long long& a_totalPointsGhosts,
			      long long& a_totalBoxes,
			      Vector<long long>& a_localLevelBoxes,
			      Vector<long long>& a_allLevelBoxes,
			      Vector<long long>& a_localLevelPoints,
			      Vector<long long>& a_totalLevelPoints,
			      const int& a_finestLevel,
			      const Vector<DisjointBoxLayout>& a_grids){
  CH_TIME("Driver::getLoadsAndBoxes");
  if(m_verbosity > 5){
    pout() << "Driver::getLoadsAndBoxes" << endl;
  }

  a_myPoints          = 0;
  a_myPointsGhosts    = 0;
  a_myBoxes           = 0;
  a_totalPoints       = 0;
  a_totalPointsGhosts = 0;
  a_totalBoxes        = 0;

  a_localLevelBoxes.resize(1 + a_finestLevel);
  a_allLevelBoxes.resize(1 + a_finestLevel);
  a_localLevelPoints.resize(1 + a_finestLevel);
  a_totalLevelPoints.resize(1 + a_finestLevel);

  const int ghost = m_amr->getNumberOfGhostCells();

  for (int lvl = 0; lvl <= a_finestLevel; lvl++){
    const DisjointBoxLayout& dbl = a_grids[lvl];
    const Vector<Box> boxes      = dbl.boxArray();
    const Vector<int> procs      = dbl.procIDs();
    
    // Find the total number of points and boxes for this level
    long long pointsThisLevel       = 0;
    long long pointsThisLevelGhosts = 0;
    long long boxesThisLevel        = 0;
    for (LayoutIterator lit = dbl.layoutIterator(); lit.ok(); ++lit){
      Box box      = dbl[lit()];
      Box grownBox = dbl[lit()];
      grownBox.grow(ghost);
      
      //
      pointsThisLevel       += box.numPts();
      pointsThisLevelGhosts += grownBox.numPts();
      boxesThisLevel        += 1;
    }
    a_totalLevelPoints[lvl] = pointsThisLevel;
    a_allLevelBoxes[lvl]  = boxesThisLevel;


    // Find the total number of points and boxes that this processor owns
    long long myPointsLevel       = 0;
    long long myPointsLevelGhosts = 0;
    long long myBoxesLevel        = 0;
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      Box box      = dbl[dit()];
      Box grownBox = dbl[dit()];
      grownBox.grow(3);
      
      myPointsLevel       += box.numPts();
      myPointsLevelGhosts += grownBox.numPts();
      myBoxesLevel        += 1;
    }


    // Total for this level
    a_totalPoints           += pointsThisLevel;
    a_totalPointsGhosts     += pointsThisLevelGhosts;
    a_totalBoxes            += boxesThisLevel;
    a_myPoints              += myPointsLevel;
    a_myPointsGhosts        += myPointsLevelGhosts;
    a_myBoxes               += myBoxesLevel;
    a_localLevelBoxes[lvl]    = myBoxesLevel;
    a_allLevelBoxes[lvl] = boxesThisLevel;
    a_localLevelPoints[lvl]   = myPointsLevel;
    a_localLevelBoxes[lvl]    = myBoxesLevel;
  }
}

const std::string Driver::numberFmt(const long long n, char sep) const {
  stringstream fmt;
  fmt << n;
  string s = fmt.str();
  s.reserve(s.length() + s.length() / 3);

  for (int i = 0, j = 3 - s.length() % 3; i < s.length(); ++i, ++j)
    if (i != 0 && j % 3 == 0)
      s.insert(i++, 1, sep);

  return s;
}

const Vector<std::string> Driver::numberFmt(const Vector<long long> a_number, char a_sep) const{
  Vector<std::string> ret(a_number.size());
  for (int i = 0; i < a_number.size(); i++){
    ret[i] = numberFmt(a_number[i], a_sep) + " ";
  }

  return ret;
}

void Driver::gridReport(){
  CH_TIME("Driver::gridReport");
  if(m_verbosity > 5){
    pout() << "Driver::gridReport" << endl;
  }

  pout() << endl;

  const int finest_level                 = m_amr->getFinestLevel();
  const Vector<DisjointBoxLayout>& grids = m_amr->getGrids(m_realm);
  const Vector<ProblemDomain>& domains   = m_amr->getDomains();
  const Vector<Real> dx                  = m_amr->getDx();

  // Grid stuff goes into here
  long long totPoints;
  long long totPointsGhosts;
  long long totBoxes;
  long long myPoints;
  long long myPointsGhosts;
  long long myBoxes;
  Vector<long long> my_level_boxes;
  Vector<long long> total_level_boxes;
  Vector<long long> my_level_points;
  Vector<long long> total_level_points;

  //
  const long long uniformPoints = (domains[finest_level].domainBox()).numPts();

  // Track memory
#ifdef CH_USE_MEMORY_TRACKING
  int BytesPerMB = 1024*1024;
  long long curMem;
  long long peakMem;
  overallMemoryUsage(curMem, peakMem);

#ifdef CH_MPI
  int unfreed_mem = curMem;
  int peak_mem    = peakMem;

  int max_unfreed_mem;
  int max_peak_mem;

  int result1 = MPI_Allreduce(&unfreed_mem, &max_unfreed_mem, 1, MPI_INT, MPI_MAX, Chombo_MPI::comm);
  int result2 = MPI_Allreduce(&peak_mem,    &max_peak_mem,    1, MPI_INT, MPI_MAX, Chombo_MPI::comm);
#endif

  //  ReportUnfreedMemory(pout());
#endif

  // Some stuff
  const ProblemDomain coarsest_domain = m_amr->getDomains()[0];
  const ProblemDomain finest_domain   = m_amr->getDomains()[finest_level];
  const Box finestBox   = finest_domain.domainBox();
  const Box coarsestBox = coarsest_domain.domainBox();
  Vector<int> refRat = m_amr->getRefinementRatios();
  Vector<int> ref_rat(finest_level);
  for (int lvl = 0; lvl < finest_level; lvl++){
    ref_rat[lvl] = refRat[lvl];
  }

  // Get boxes for each Realm
  const std::vector<std::string> Realms = m_amr->getRealms();
  for (auto str : Realms){
    this->getLoadsAndBoxes(myPoints,
			   myPointsGhosts,
			   myBoxes,
			   totPoints,
			   totPointsGhosts,
			   totBoxes,
			   my_level_boxes,
			   total_level_boxes,
			   my_level_points,
			   total_level_points,
			   finest_level,
			   m_amr->getGrids(str));
  }

  // Begin writing a report. 
  pout() << "-----------------------------------------------------------------------" << endl
	 << "Driver::Grid report - timestep = " << m_timeStep << endl
	 << "\t\t\t        Finest level           = " << finest_level << endl
	 << "\t\t\t        Finest domain          = " << finestBox.size()[0] << " x " << finestBox.size()[1] <<
#if CH_SPACEDIM==2
    endl
#elif CH_SPACEDIM==3
    " x " << finestBox.size()[2] << endl
#endif
	 << "\t\t\t        Coarsest domain        = " << coarsestBox.size()[0] << " x " << coarsestBox.size()[1] <<
#if CH_SPACEDIM==2
    endl
#elif CH_SPACEDIM==3
    " x " << coarsestBox.size()[2] << endl
#endif
	 << "\t\t\t        Refinement ratios      = " << ref_rat << endl
	 << "\t\t\t        Grid sparsity          = " << 1.0*totPoints/uniformPoints << endl
	 << "\t\t\t        Finest dx              = " << dx[finest_level] << endl
	 << "\t\t\t        Total number boxes     = " << numberFmt(totBoxes) << endl
	 << "\t\t\t        Number of valid cells  = " << numberFmt(totPoints) << endl
	 << "\t\t\t        Including ghost cells  = " << numberFmt(totPointsGhosts)  << endl
	 << "\t\t\t        Total # of boxes (lvl) = " << numberFmt(total_level_boxes) << endl
	 << "\t\t\t        Total # of cells (lvl) = " << numberFmt(total_level_points) << endl;

  // Do a local report for each Realm
  for (auto str : Realms){
    this->getLoadsAndBoxes(myPoints,
			   myPointsGhosts,
			   myBoxes,
			   totPoints,
			   totPointsGhosts,
			   totBoxes,
			   my_level_boxes,
			   total_level_boxes,
			   my_level_points,
			   total_level_points,
			   finest_level,
			   m_amr->getGrids(str));
    pout() << "\t\t\t        Realm = " << str << endl
	   << "\t\t\t\t        Proc. # of valid cells = " << numberFmt(myPoints) << endl
	   << "\t\t\t\t        Including ghost cells  = " << numberFmt(myPointsGhosts) << endl
	   << "\t\t\t\t        Proc. # of boxes       = " << numberFmt(myBoxes) << endl
	   << "\t\t\t\t        Proc. # of boxes (lvl) = " << numberFmt(my_level_boxes) << endl
	   << "\t\t\t\t        Proc. # of cells (lvl) = " << numberFmt(my_level_points) << endl;
  }
  
  pout()
#ifdef CH_USE_MEMORY_TRACKING
    << "\t\t\t        Unfreed memory        = " << curMem/BytesPerMB << " (MB)" << endl
    << "\t\t\t        Peak memory usage     = " << peakMem/BytesPerMB << " (MB)" << endl
#ifdef CH_MPI
    << "\t\t\t        Max unfreed memory    = " << max_unfreed_mem/BytesPerMB << " (MB)" << endl
    << "\t\t\t        Max peak memory       = " << max_peak_mem/BytesPerMB << " (MB)" << endl
#endif
    << "-----------------------------------------------------------------------" << endl
#endif
    << endl;

  pout() << endl;
}

void Driver::memoryReport(const MemoryReportMode a_mode){
#ifdef CH_USE_MEMORY_TRACKING
  CH_TIME("Driver::gridReport");
  if(m_verbosity > 5){
    pout() << "Driver::gridReport" << endl;
  }

  if(a_mode == MemoryReportMode::overall){
    overallMemoryUsage();
  }
  else if(a_mode == MemoryReportMode::unfreed){
    ReportUnfreedMemory(pout());
  }
  else if(a_mode == MemoryReportMode::allocated){
    ReportAllocatedMemory(pout());
  }
  pout() << endl;
#endif
}

void Driver::regrid(const int a_lmin, const int a_lmax, const bool a_useInitialData){
  CH_TIME("Driver::regrid");
  if(m_verbosity > 2){
    pout() << "Driver::regrid" << endl;
  }

  // We need to be careful with memory allocations here. Therefore we do:
  // --------------------------------------------------------------------
  // 1.  Tag cells, this calls CellTagger which allocates and deallocate its own storage so 
  //     there's a peak in memory consumption here. We have to eat this one because we
  //     potentially need all the solver data for tagging, so that data can't be touched.
  //     If we don't get new tags, we exit this routine already here. 
  // 2.  Deallocate internal storage for the TimeStepper - this frees up a bunch of memory that
  //     we don't need since we won't advance until after regridding anyways. 
  // 3.  Cache tags, this doubles up on the memory for m_tags but that shouldn't matter.
  // 4.  Free up m_tags for safety because it will be regridded anyways. 
  // 5.  Cache solver states
  // 6.  Deallocate internal storage for solver. This should free up a bunch of memory. 
  // 7.  Regrid AmrMesh - this shouldn't cause any extra memory issues
  // 8.  Regrid Driver - this
  // 9.  Regrid the cell tagger. I'm not explicitly releasing storage from here since it's so small.
  // 10. Solve elliptic equations and fill solvers


  Vector<IntVectSet> tags;

  const Real start_time = MPI_Wtime();   // Timer

  // We are allowing geometric tags to change under the hood, but we need a method for detecting if they changed. If they did,
  // we certainly have to regrid.
  if(m_needsNewGeometricTags){
    this->getGeometryTags();
  }

  const bool newCellTags = this->tagCells(tags, m_tags); // Tag cells using the cell tagger

  if(!newCellTags && !m_needsNewGeometricTags){
    if(a_useInitialData){
      m_timeStepper->initialData();
    }

    if(m_verbosity > 1){
      pout() << "\nDriver::regrid - Didn't find any new cell tags. Skipping the regrid step\n" << endl;
    }
    return;
  }
  else{ // Compact tags
    for (int i = 0; i < tags.size(); i++){
      tags[i].compact();
    }
  }

  // Store things that need to be regridded
  this->cacheTags(m_tags);              // Cache m_tags because after regrid, ownership will change
  m_timeStepper->preRegrid(a_lmin, m_amr->getFinestLevel());

  // Deallocate unnecessary storage
  this->deallocateInternals();          // Deallocate internal storage for Driver
  
  const Real cell_tags = MPI_Wtime();    // Timer

  // Regrid AMR. Only levels [lmin, lmax] are allowed to change. 
  const int old_finestLevel = m_amr->getFinestLevel();
  m_amr->regridAmr(tags, a_lmin, a_lmax);
  const int new_finestLevel = m_amr->getFinestLevel();

  // Load balance and regrid the various Realms
  const std::vector<std::string>& Realms = m_amr->getRealms();
  for (const auto& str : Realms){
    if(m_timeStepper->loadBalanceThisRealm(str)){
      
      Vector<Vector<int> > procs;
      Vector<Vector<Box> > boxes;
      
      m_timeStepper->loadBalanceBoxes(procs, boxes, str, m_amr->getProxyGrids(), a_lmin, new_finestLevel);

      m_amr->regridRealm(str, procs, boxes, a_lmin);
    }
  }


  // Regrid the operators
  m_amr->regridOperators(a_lmin);
  const Real base_regrid = MPI_Wtime(); // Base regrid time

  // Regrid Driver, timestepper, and celltagger
  this->regridInternals(old_finestLevel, new_finestLevel);          // Regrid internals for Driver
  m_timeStepper->regrid(a_lmin, old_finestLevel, new_finestLevel);   // Regrid solvers
  if(a_useInitialData){
    m_timeStepper->initialData();
  }

  // Regrid cell tagger if we have one. 
  if(!m_cellTagger.isNull()){
    m_cellTagger->regrid();             
  }

  // If it wants to, TimeStepper can do a postRegrid operation. 
  m_timeStepper->postRegrid();

  const Real solver_regrid = MPI_Wtime(); // Timer

  if(m_verbosity > 1){
    this->regridReport(solver_regrid - start_time,
		       cell_tags - start_time,
		       base_regrid - cell_tags,
		       solver_regrid - base_regrid);
  }

  m_needsNewGeometricTags = false;
}

void Driver::regridInternals(const int a_oldFinestLevel, const int a_newFinestLevel){
  CH_TIME("Driver::regridInternals");
  if(m_verbosity > 2){
    pout() << "Driver::regridInternals" << endl;
  }

  this->allocateInternals();

  // Copy cached tags back over to m_tags
  for (int lvl = 0; lvl <= Min(a_oldFinestLevel, a_newFinestLevel); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    
    // Copy mask
    LevelData<BaseFab<bool> > tmp;
    tmp.define(dbl, 1, IntVect::Zero);
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      tmp[dit()].setVal(false);
    }

    m_cachedTags[lvl]->copyTo(tmp);
    
    for(DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const BaseFab<bool>& tmpFab = tmp[dit()];
      const Box& box = dbl.get(dit());

      DenseIntVectSet& tags = (*m_tags[lvl])[dit()];
      
      for (BoxIterator bit(box); bit.ok(); ++bit){
	const IntVect iv = bit();
	if(tmpFab(iv,0)){
	  tags |= iv;
	}
      }
    }
  }
}

void Driver::regridReport(const Real a_totalTime,
			  const Real a_tagTime,
			  const Real a_baseRegridTime,
			  const Real a_solverRegridTime){
  CH_TIME("Driver::regridReport");
  if(m_verbosity > 5){
    pout() << "Driver::regridReport" << endl;
  }

  const Real elapsed    = a_totalTime;
  const int elapsed_hrs = floor(elapsed/3600);
  const int elapsed_min = floor((elapsed - 3600*elapsed_hrs)/60);
  const int elapsed_sec = floor( elapsed - 3600*elapsed_hrs - 60*elapsed_min);
  const int elapsed_ms  = floor((elapsed - 3600*elapsed_hrs - 60*elapsed_min - elapsed_sec)*1000.);

  char metrics[30];
  sprintf(metrics, "%3.3ih %2.2im %2.2is %3.3ims",
	  elapsed_hrs, 
	  elapsed_min, 
	  elapsed_sec, 
	  elapsed_ms);

  pout() << "-----------------------------------------------------------------------" << endl
	 << "Driver::regridReport breakdown - Time step #" << m_timeStep << endl
	 << "\t\t\t" << "Total regrid time : " << metrics << endl
	 << "\t\t\t" << "Cell tagging      : " << 100.*(a_tagTime/a_totalTime) << "%" << endl
	 << "\t\t\t" << "Base regrid       : " << 100.*(a_baseRegridTime/a_totalTime) << "%" << endl
	 << "\t\t\t" << "Solver regrid     : " << 100.*(a_solverRegridTime/a_totalTime) << "%" << endl
	 << "-----------------------------------------------------------------------" << endl;
}

void Driver::run(const Real a_startTime, const Real a_endTime, const int a_maxSteps){
  CH_TIME("Driver::run");
  if(m_verbosity > 1){
    pout() << "Driver::run" << endl;
  }

  if(m_verbosity > 0){
    pout() << "=================================" << endl;
    if(!m_restart){
      pout() << "Driver::run -- starting run" << endl;
    }
    else{
      pout() << "Driver::run -- restarting run" << endl;
    }
  }

  if(a_maxSteps > 0){
    if(!m_restart){
      m_time = a_startTime;
      m_timeStep = 0;
    }

    m_timeStepper->computeDt(m_dt, m_timeCode);
    m_timeStepper->synchronizeSolverTimes(m_timeStep, m_time, m_dt);

    bool last_step     = false;
    bool first_step    = true;
    const Real init_dt = m_dt;

    if(m_verbosity > 0){
      this->gridReport();
    }

    m_wallClockStart = MPI_Wtime();

    while(m_time < a_endTime && m_timeStep < a_maxSteps && !last_step){
      const int max_sim_depth = m_amr->getMaxSimulationDepth();
      const int max_amr_depth = m_amr->getMaxAmrDepth();

      // This is the regrid test. We do some dummy tests first and then do the recursive/non-recursive stuff
      // inside the loop. 
      const bool can_regrid        = max_sim_depth > 0 && max_amr_depth > 0;
      const bool check_step        = m_timeStep%m_regridInterval == 0 && m_regridInterval > 0;
      const bool check_timeStepper = m_timeStepper->needToRegrid() && m_regridInterval > 0;
      if(can_regrid && (check_step || check_timeStepper)){
	if(!first_step){

	  // We'll regrid levels lmin through lmax. As always, new grids on level l are generated through tags
	  // on levels (l-1);
	  const int lmin = 0;
	  const int lmax = m_amr->getFinestLevel() + 1; // This means that if we refine, we can only add one level at a time. 

	  this->regrid(lmin, lmax, false);
	  if(m_verbosity > 0){
	    this->gridReport();
	  }
	  if(m_writeRegridFiles){
	    this->writeRegridFile();
	  }
	}
      }

      if(!first_step){
	m_timeStepper->computeDt(m_dt, m_timeCode);
      }

      if(first_step){
	first_step = false;
      }

      // Did the time step become too small?
      if(m_dt < 1.0E-5*init_dt){
	m_timeStep++;

	if(m_writeMemory){
	  this->writeMemoryUsage();
	}
	if(m_writeLoads){
	  this->writeComputationalLoads();
	}
#ifdef CH_USE_HDF5
	this->writeCrashFile();
	//	this->writeCheckpointFile();
#endif

	MayDay::Abort("Driver::run - the time step became too small");
      }

      // Last time step can be smaller than m_dt so that we end on a_endTime
      if(m_time + m_dt > a_endTime){
	m_dt = a_endTime - m_time;
	last_step = true;
      }


      // Time stepper advances solutions
      m_wallClockOne = MPI_Wtime();
      const Real actual_dt = m_timeStepper->advance(m_dt);
      m_wallClockTwo = MPI_Wtime();

      // Synchronize times
      m_dt    = actual_dt;
      m_time += actual_dt;
      m_timeStep += 1;
      m_timeStepper->synchronizeSolverTimes(m_timeStep, m_time, m_dt);

      if(Abs(m_time - a_endTime) < m_dt*1.E-5){
	last_step = true;
      }
      if(m_timeStep == m_maxSteps){
	last_step = true;
      }

      // Print step report
      if(m_verbosity > 0){
	this->stepReport(a_startTime, a_endTime, a_maxSteps);
      }

#ifdef CH_USE_HDF5
      if(m_plotInterval > 0){

	// Aux data
	if(m_writeMemory){
	  this->writeMemoryUsage();
	}
	if(m_writeLoads){
	  this->writeComputationalLoads();
	}
	
	// Plot file
	if(m_timeStep%m_plotInterval == 0 || last_step == true){
	  if(m_verbosity > 2){
	    pout() << "Driver::run -- Writing plot file" << endl;
	  }

	  this->writePlotFile();
	}
      }

      // Write checkpoint file
      if(m_timeStep % m_checkpointInterval == 0 && m_checkpointInterval > 0 || last_step == true && m_checkpointInterval > 0){
	if(m_verbosity > 2){
	  pout() << "Driver::run -- Writing checkpoint file" << endl;
	}
	this->writeCheckpointFile();
      }
#endif

      // Rebuild input parameters
      this->rebuildParmParse();
      this->parseRuntimeOptions();
      m_amr->parseRuntimeOptions();
      m_timeStepper->parseRuntimeOptions();
      if(!m_cellTagger.isNull()){
	m_cellTagger->parseRuntimeOptions();
      }
    }
  }

  if(m_verbosity > 0){
    this->gridReport();
  }

  if(m_verbosity > 0){
    pout() << "==================================" << endl;
    pout() << "Driver::run -- ending run  " << endl;
    pout() << "==================================" << endl;
  }
}

void Driver::setupAndRun(const std::string a_inputFile){
  CH_TIME("Driver::setupAndRun");
  if(m_verbosity > 0){
    pout() << "Driver::setupAndRun" << endl;
  }

  char iter_str[100];
  sprintf(iter_str, ".check%07d.%dd.hdf5", m_restartStep, SpaceDim);
  const std::string restart_file = m_outputDirectory + "/chk/" + m_outputFileNames + std::string(iter_str);

  this->setup(a_inputFile, m_initialRegrids, m_restart, restart_file);

  if(!m_geometryOnly){
    this->run(m_startTime, m_stopTime, m_maxSteps);
  }
}

void Driver::rebuildParmParse() const {
  ParmParse pp("Driver");

  pp.redefine(m_inputFile.c_str());
}

void Driver::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry){
  CH_TIME("Driver::setComputationalGeometry");
  if(m_verbosity > 5){
    pout() << "Driver::setComputationalGeometry" << endl;
  }
  m_computationalGeometry = a_computationalGeometry;
  m_multifluidIndexSpace     = a_computationalGeometry->getMfIndexSpace();
}

void Driver::setTimeStepper(const RefCountedPtr<TimeStepper>& a_timeStepper){
  CH_TIME("Driver::setTimeStepper");
  if(m_verbosity > 5){
    pout() << "Driver::setTimeStepper" << endl;
  }
  m_timeStepper = a_timeStepper;
}

void Driver::setCellTagger(const RefCountedPtr<CellTagger>& a_cellTagger){
  CH_TIME("Driver::setCellTagger");
  if(m_verbosity > 5){
    pout() << "Driver::setCellTagger" << endl;
  }

  m_cellTagger = a_cellTagger;
  if(!a_cellTagger.isNull()){
    m_cellTagger->parseOptions();
  }
}

void Driver::setGeoCoarsen(const RefCountedPtr<GeoCoarsener>& a_geoCoarsen){
  CH_TIME("Driver::setGeoCoarsen");
  if(m_verbosity > 5){
    pout() << "Driver::setGeoCoarsen" << endl;
  }
  m_geoCoarsen = a_geoCoarsen;
}

void Driver::parseOptions(){
  CH_TIME("Driver::parseOptions");

  ParmParse pp("Driver");

  pp.get("verbosity",                m_verbosity);
  if(m_verbosity > 5){
    pout() << "Driver::parseOptions" << endl;
  }
  
  pp.get("regrid_interval",          m_regridInterval);
  pp.get("initial_regrids",          m_initialRegrids);
  pp.get("restart",                  m_restartStep); 
  pp.get("write_memory",             m_writeMemory);
  pp.get("write_loads",              m_writeLoads);
  pp.get("output_directory",         m_outputDirectory);
  pp.get("output_names",             m_outputFileNames);
  pp.get("plot_interval",            m_plotInterval);
  pp.get("checkpoint_interval",      m_checkpointInterval);
  pp.get("write_regrid_files",       m_writeRegridFiles);
  pp.get("write_restart_files",      m_writeRestartFiles);
  pp.get("num_plot_ghost",           m_numPlotGhost);
  pp.get("allow_coarsening",         m_allowCoarsening);
  pp.get("geometry_only",            m_geometryOnly);
  pp.get("ebis_memory_load_balance", m_ebisMemoryLoadBalance);
  pp.get("max_steps",                m_maxSteps);
  pp.get("start_time",               m_startTime);
  pp.get("stop_time",                m_stopTime);
  pp.get("max_plot_depth",           m_maxPlotDepth);
  pp.get("max_chk_depth",            m_maxCheckpointDepth);

  m_restart   = (m_restartStep > 0) ? true : false;

  // Stuff that's a little too verbose to include here directly. 
  this->parsePlotVariables();
  this->parseGeometryGeneration();
  this->parseGeometryRefinement();
  this->parseIrregTagGrowth();
}

void Driver::parseRuntimeOptions(){
  CH_TIME("Driver::parseRuntimeOptions");

  ParmParse pp("Driver");

  pp.get("verbosity",                m_verbosity);
  if(m_verbosity > 5){
    pout() << "Driver::parseRuntimeOptions" << endl;
  }
  pp.get("write_memory",             m_writeMemory);
  pp.get("write_loads",              m_writeLoads);
  pp.get("plot_interval",            m_plotInterval);
  pp.get("regrid_interval",          m_regridInterval);
  pp.get("checkpoint_interval",      m_checkpointInterval);
  pp.get("write_regrid_files",       m_writeRegridFiles);
  pp.get("write_restart_files",      m_writeRestartFiles);
  pp.get("num_plot_ghost",           m_numPlotGhost);
  pp.get("allow_coarsening",         m_allowCoarsening);
  pp.get("max_steps",                m_maxSteps);
  pp.get("stop_time",                m_stopTime);

  this->parseGeometryRefinement();  
  this->parsePlotVariables();
  this->parseIrregTagGrowth();  
}

void Driver::parsePlotVariables(){
  ParmParse pp("Driver");
  const int num = pp.countval("plt_vars");
  Vector<std::string> str(num);
  pp.getarr("plt_vars", str, 0, num);

  m_plotTags     = false;
  m_plotRanks    = false;
  m_plotLevelset = false;
  
  for (int i = 0; i < num; i++){
    if(     str[i] == "tags")     m_plotTags     = true;
    else if(str[i] == "mpi_rank") m_plotRanks    = true;
    else if(str[i] == "levelset") m_plotLevelset = true;
  }
}

void Driver::parseIrregTagGrowth(){
  CH_TIME("Driver::parseIrregTagGrowth()");
  if(m_verbosity > 5){
    pout() << "Driver::parseIrregTagGrowth()" << endl;
  }

  ParmParse pp("Driver");

  pp.get("grow_geo_tags", m_irregTagGrowth);

  m_irregTagGrowth = std::max(0, m_irregTagGrowth);
}

void Driver::parseGeometryGeneration(){
  CH_TIME("Driver::parseGeometryGeneration");
  if(m_verbosity > 5){
    pout() << "Driver::parseGeometryGeneration" << endl;
  }

  ParmParse pp("Driver");
  pp.get("geometry_generation", m_geometryGeneration);
  pp.get("geometry_scan_level", m_geoScanLevel);
  

  if(m_geometryGeneration == "chombo-discharge"){ // Need to activate some flags that trigger the correct code. 
    ComputationalGeometry::s_use_new_gshop = true;
    EBISLevel::s_distributedData            = true;
    ComputationalGeometry::s_ScanDomain = m_amr->getDomains()[m_geoScanLevel];
  }
  else if(m_geometryGeneration == "chombo"){
  }
  else{
    MayDay::Abort("Driver:parseGeometryGeneration - unsupported argument requested");
  }
}

void Driver::parseGeometryRefinement(){
  CH_TIME("Driver::set_geom_refinement_depth");
  if(m_verbosity > 5){
    pout() << "Driver::set_geom_refinement_depth" << endl;
  }

  ParmParse pp("Driver");

  const int max_depth = m_amr->getMaxAmrDepth();

  // These are special options. We are allowing geometry refinement criteria to change as simulations progress (i.e. m_timeStep > 0) where we call
  // this routine again. But we need something to tell us that we got new refinement criteria for geometric tags so we can avoid regrid if they didn't change.
  // This is my clunky way of doing that. 
  
  const auto c1 = m_refineAngle;
  const auto c2 = m_conductorTagsDepth;
  const auto c3 = m_dielectricTagsDepth;
  const auto c4 = m_gasConductorInterfaceTagDepth;
  const auto c5 = m_gasDielectricInterfaceTagDepth;
  const auto c6 = m_gasSolidInterfaceTagDepth;
  const auto c7 = m_solidSolidInterfaceTagDepth;


  int depth0;

  pp.get("refine_angles",                   m_refineAngle);
  pp.get("refine_geometry",                 depth0);
  pp.get("refine_electrodes",               m_conductorTagsDepth);
  pp.get("refine_dielectrics",              m_dielectricTagsDepth);
  pp.get("refine_electrode_gas_interface",  m_gasConductorInterfaceTagDepth);
  pp.get("refine_dielectric_gas_interface", m_gasDielectricInterfaceTagDepth);
  pp.get("refine_solid_gas_interface",      m_gasSolidInterfaceTagDepth);
  pp.get("refine_solid_solid_interface",    m_solidSolidInterfaceTagDepth);

  depth0                            = (depth0                           < 0) ? max_depth : depth0;
  m_conductorTagsDepth              = (m_conductorTagsDepth             < 0) ? depth0 : m_conductorTagsDepth;
  m_dielectricTagsDepth             = (m_dielectricTagsDepth            < 0) ? depth0 : m_dielectricTagsDepth;
  m_gasConductorInterfaceTagDepth   = (m_gasConductorInterfaceTagDepth  < 0) ? depth0 : m_gasConductorInterfaceTagDepth;
  m_gasDielectricInterfaceTagDepth  = (m_gasDielectricInterfaceTagDepth < 0) ? depth0 : m_gasDielectricInterfaceTagDepth;
  m_gasSolidInterfaceTagDepth       = (m_gasSolidInterfaceTagDepth      < 0) ? depth0 : m_gasSolidInterfaceTagDepth;
  m_solidSolidInterfaceTagDepth     = (m_solidSolidInterfaceTagDepth    < 0) ? depth0 : m_solidSolidInterfaceTagDepth;

  if(m_timeStep > 0){ // Simulation is already running, and we need to check if we need new geometric tags for regridding. 
    if(c1 != m_refineAngle                    ||
       c2 != m_conductorTagsDepth             ||
       c3 != m_dielectricTagsDepth            ||
       c4 != m_gasConductorInterfaceTagDepth  ||
       c5 != m_gasDielectricInterfaceTagDepth ||
       c6 != m_gasSolidInterfaceTagDepth      ||
       c7 != m_solidSolidInterfaceTagDepth){
      
      m_needsNewGeometricTags = true;
    }
  }
  else{ // Fresh simulation -- we should have gotten geometric tags before entering regrid so we can set this to false. 
    m_needsNewGeometricTags = false;
  }
}

void Driver::createOutputDirectories(){
  CH_TIME("Driver::createOutputDirectories");
  if(m_verbosity > 5){
    pout() << "Driver::createOutputDirectories" << endl;
  }

  // If directory does not exist, create it
  int success = 0;
  if(procID() == 0){
    std::string cmd;

    cmd = "mkdir -p " + m_outputDirectory;
    success = system(cmd.c_str());
    if(success != 0){
      std::cout << "Driver::set_outputDirectoryectory - master could not create directory" << std::endl;
    }

    cmd = "mkdir -p " + m_outputDirectory + "/plt";
    success = system(cmd.c_str());
    if(success != 0){
      std::cout << "Driver::set_outputDirectoryectory - master could not create plot directory" << std::endl;
    }

    cmd = "mkdir -p " + m_outputDirectory + "/geo";
    success = system(cmd.c_str());
    if(success != 0){
      std::cout << "Driver::set_outputDirectoryectory - master could not create geo directory" << std::endl;
    }

    cmd = "mkdir -p " + m_outputDirectory + "/chk";
    success = system(cmd.c_str());
    if(success != 0){
      std::cout << "Driver::set_outputDirectoryectory - master could not create checkpoint directory" << std::endl;
    }

    cmd = "mkdir -p " + m_outputDirectory + "/mpi";
    success = system(cmd.c_str());
    if(success != 0){
      std::cout << "Driver::set_outputDirectoryectory - master could not create mpi directory" << std::endl;
    }

    cmd = "mkdir -p " + m_outputDirectory + "/mpi/memory";
    success = system(cmd.c_str());
    if(success != 0){
      std::cout << "Driver::set_outputDirectoryectory - master could not create mpi/memory directory" << std::endl;
    }

    cmd = "mkdir -p " + m_outputDirectory + "/mpi/loads";
    success = system(cmd.c_str());
    if(success != 0){
      std::cout << "Driver::set_outputDirectoryectory - master could not create mpi/loads directory" << std::endl;
    }

    cmd = "mkdir -p " + m_outputDirectory + "/regrid";
    success = system(cmd.c_str());
    if(success != 0){
      std::cout << "Driver::set_outputDirectoryectory - master could not create regrid directory" << std::endl;
    }

    cmd = "mkdir -p " + m_outputDirectory + "/restart";
    success = system(cmd.c_str());
    if(success != 0){
      std::cout << "Driver::set_outputDirectoryectory - master could not create restart directory" << std::endl;
    }

    cmd = "mkdir -p " + m_outputDirectory + "/crash";
    success = system(cmd.c_str());
    if(success != 0){
      std::cout << "Driver::set_outputDirectoryectory - master could not create crash directory" << std::endl;
    }    
  }
  
  MPI_Barrier(Chombo_MPI::comm);
  if(success != 0){
    MayDay::Abort("Driver::set_outputDirectoryectory - could not create directories for output");
  }
}



void Driver::setAmr(const RefCountedPtr<AmrMesh>& a_amr){
  CH_TIME("Driver::setAmr");
  if(m_verbosity > 5){
    pout() << "Driver::setAmr" << endl;
  }

  m_amr = a_amr;
  m_amr->setMultifluidIndexSpace(m_computationalGeometry->getMfIndexSpace());

}

void Driver::setup(const std::string a_inputFile, const int a_initialRegrids, const bool a_restart, const std::string a_restartFile){
  CH_TIME("Driver::setup");
  if(m_verbosity > 5){
    pout() << "Driver::setup" << endl;
  }

  m_inputFile = a_inputFile;

  if(m_geometryOnly){
    this->setupGeometryOnly();
  }
  else{
    if(!a_restart){
      this->setupFresh(a_initialRegrids);
#ifdef CH_USE_HDF5
      if(m_plotInterval > 0){
	if(m_writeMemory){
	  this->writeMemoryUsage();
	}
	if(m_writeLoads){
	  this->writeComputationalLoads();
	}
	this->writePlotFile();
      }
#endif
    }
    else{
      this->setupForRestart(a_initialRegrids, a_restartFile);
    }
  }
}

void Driver::setupGeometryOnly(){
  CH_TIME("Driver::setupGeometryOnly");
  if(m_verbosity > 5){
    pout() << "Driver::setupGeometryOnly" << endl;
  }

  this->sanityCheck();

  if(m_ebisMemoryLoadBalance){
    EBIndexSpace::s_useMemoryLoadBalance = true;
  }
  else {
    EBIndexSpace::s_useMemoryLoadBalance = false;
  }

  const Real t0 = MPI_Wtime();
  m_computationalGeometry->buildGeometries(m_amr->getFinestDomain(),
					    m_amr->getProbLo(),
					    m_amr->getFinestDx(),
					    m_amr->getMaxEbisBoxSize());
  const Real t1 = MPI_Wtime();
  if(procID() == 0) std::cout << "geotime = " << t1 - t0 << std::endl;

  // Set implicit functions now. 
  m_amr->setBaseImplicitFunction(phase::gas,   m_computationalGeometry->getGasImplicitFunction());
  m_amr->setBaseImplicitFunction(phase::solid, m_computationalGeometry->getSolidImplicitFunction());

  if(m_writeMemory){
    this->writeMemoryUsage();
  }

  this->getGeometryTags();       // Get geometric tags.

  if(m_writeMemory){
    this->writeMemoryUsage();
  }

  // Regrid using geometric tags only.
  m_amr->regridAmr(m_geomTags, 0);   

  if(m_verbosity > 0){
    this->gridReport();
  }

  //  this->writeMemoryUsage();
  if(m_plotInterval > 0){
    this->writeGeometry();                             // Write geometry only
  }
}

void Driver::setupFresh(const int a_initialRegrids){
  CH_TIME("Driver::setupFresh");
  if(m_verbosity > 5){
    pout() << "Driver::setupFresh" << endl;
  }

  this->sanityCheck();                                    // Sanity check before doing anything expensive

  if(m_ebisMemoryLoadBalance){
    EBIndexSpace::s_useMemoryLoadBalance = true;
  }
  else {
    EBIndexSpace::s_useMemoryLoadBalance = false;
  }

  m_computationalGeometry->buildGeometries(m_amr->getFinestDomain(),
					    m_amr->getProbLo(),
					    m_amr->getFinestDx(),
					    m_amr->getMaxEbisBoxSize());


  // Register Realms
  m_timeStepper->setAmr(m_amr);
  m_timeStepper->registerRealms();

  // Set implicit functions now. 
  m_amr->setBaseImplicitFunction(phase::gas,   m_computationalGeometry->getGasImplicitFunction());
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
  m_timeStepper->setComputationalGeometry(m_computationalGeometry);       // Set computational geometry

  // TimeStepper setup
  m_timeStepper->setupSolvers();                                 // Instantiate solvers
  m_timeStepper->synchronizeSolverTimes(m_timeStep, m_time, m_dt);  // Sync solver times
  m_timeStepper->registerOperators();
  m_amr->regridOperators(lmin);
  m_timeStepper->allocate();

  // Fill solves with initial data
  m_timeStepper->initialData();                                  // Fill solvers with initial data

  // We now load balance and define operators and stuff like that. 
  this->cacheTags(m_tags);
  m_timeStepper->preRegrid(lmin, lmax);
  for (const auto& str : m_amr->getRealms()){
    if(m_timeStepper->loadBalanceThisRealm(str)){
      
      Vector<Vector<int> > procs;
      Vector<Vector<Box> > boxes;

      const int lmin   = 0;
      const int lmax = m_amr->getFinestLevel(); 
      
      m_timeStepper->loadBalanceBoxes(procs, boxes, str, m_amr->getProxyGrids(), lmin, lmax);

      m_amr->regridRealm(str, procs, boxes, lmin);
    }
  }
  m_amr->regridOperators(lmin);             // Regrid operators again.
  this->regridInternals(lmax, lmax);        // Regrid internals for Driver.
  m_timeStepper->regrid(lmin, lmax, lmax);  // Regrid solvers.
  m_timeStepper->initialData();             // Need to fill with initial data again. 
  
  // Do post initialize stuff
  m_timeStepper->postInitialize();

  // CellTagger
  if(!m_cellTagger.isNull()){
    m_cellTagger->regrid();
  }

  // Do a grid report of the initial grid
  if(m_verbosity > 0){
    this->gridReport();
  }

  // Initial regrids
  for (int i = 0; i < a_initialRegrids; i++){
    if(m_verbosity > 5){
      pout() << "Driver -- initial regrid # " << i + 1 << endl;
    }

    const int lmin = 0;
    const int lmax = m_amr->getFinestLevel() + 1;
    this->regrid(lmin, lmax, true);

    if(m_verbosity > 0){
      this->gridReport();
    }
  }
}

void Driver::setupForRestart(const int a_initialRegrids, const std::string a_restartFile){
  CH_TIME("Driver::setupForRestart");
  if(m_verbosity > 5){
    pout() << "Driver::setupForRestart" << endl;
  }

  this->checkRestartFile(a_restartFile);

  this->sanityCheck();                                    // Sanity check before doing anything expensive

  m_computationalGeometry->buildGeometries(m_amr->getFinestDomain(),
					    m_amr->getProbLo(),
					    m_amr->getFinestDx(),
					    m_amr->getMaxEbisBoxSize());

  this->getGeometryTags();       // Get geometric tags.

  m_timeStepper->setAmr(m_amr);                         // Set amr
  m_timeStepper->registerRealms();                      // Register Realms
  m_timeStepper->setComputationalGeometry(m_computationalGeometry); // Set computational geometry

  // Set implicit functions now. 
  m_amr->setBaseImplicitFunction(phase::gas,   m_computationalGeometry->getGasImplicitFunction());
  m_amr->setBaseImplicitFunction(phase::solid, m_computationalGeometry->getSolidImplicitFunction());

  // Read checkpoint file
  this->readCheckpointFile(a_restartFile); // Read checkpoint file - this sets up amr, instantiates solvers and fills them

  // Time stepper does post checkpoint setup
  m_timeStepper->postCheckpointSetup();
  
  // Prepare storage for CellTagger
  if(!m_cellTagger.isNull()){
    m_cellTagger->regrid();         
  }

  if(m_writeRestartFiles){
    this->writeRestartFile();
  }

  // Initial regrids
  for (int i = 0; i < a_initialRegrids; i++){
    if(m_verbosity > 0){
      pout() << "Driver -- initial regrid # " << i + 1 << endl;
    }

    const int lmin = 0;
    const int lmax = m_amr->getFinestLevel() + 1;
    this->regrid(lmin, lmax, false);

    if(m_verbosity > 0){
      this->gridReport();
    }
  }


}

void Driver::checkRestartFile(const std::string a_restartFile) const {
  CH_TIME("Driver::checkRestartFile");
  if(m_verbosity > 4){
    pout() << "Driver::checkRestartFile" << endl;
  }

  ifstream f(a_restartFile.c_str());
  if(!f.good()){
    pout() << "Driver::checkRestartFile - could not find file = " << a_restartFile << endl;
    MayDay::Abort("Driver::checkRestartFile - abort, could not find file");
  }
}

void Driver::sanityCheck(){
  CH_TIME("Driver::sanityCheck");
  if(m_verbosity > 4){
    pout() << "Driver::sanityCheck" << endl;
  }

  CH_assert(!m_timeStepper.isNull());
}

void Driver::stepReport(const Real a_startTime, const Real a_endTime, const int a_maxSteps){
  CH_TIME("Driver::stepReport");
  if(m_verbosity > 5){
    pout() << "Driver::stepReport" << endl;
  }

  pout() << endl;
  pout() << "Driver::Time step report -- Time step #" << m_timeStep << endl
	 << "                                   Time  = " << m_time << endl
	 << "                                   dt    = " << m_dt << endl;

  m_timeStepper->printStepReport();

  // Get the total number of poitns across all levels
  const int finest_level                 = m_amr->getFinestLevel();
  const Vector<DisjointBoxLayout>& grids = m_amr->getGrids(m_realm);
  const Vector<ProblemDomain>& domains   = m_amr->getDomains();
  const Vector<Real>& dx                 = m_amr->getDx();
  long long totalPoints = 0;
  long long uniformPoints = (domains[finest_level].domainBox()).numPts();
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    long long pointsThisLevel = 0;
    for (LayoutIterator lit = grids[lvl].layoutIterator(); lit.ok(); ++lit){
      pointsThisLevel += grids[lvl][lit()].numPts();
    }
    totalPoints += pointsThisLevel;
  }

  char metrics[300];

  // Percentage completed of time steps
  const Real percentStep = (1.0*m_timeStep/a_maxSteps)*100.;
  sprintf(metrics,"%31c -- %5.2f percentage of time steps completed",' ', percentStep);
  pout() << metrics << endl;

  const Real percentTime = ((m_time - a_startTime)/(a_endTime - a_startTime))*100.;
  sprintf(metrics,"%31c -- %5.2f percentage of simulation time completed",' ', percentTime);
  pout() << metrics << endl;


  // Hours, minutes, seconds and millisecond of the previous iteration
  const Real elapsed   = m_wallClockTwo - m_wallClockStart;
  const int elapsedHrs = floor(elapsed/3600);
  const int elapsedMin = floor((elapsed - 3600*elapsedHrs)/60);
  const int elapsedSec = floor( elapsed - 3600*elapsedHrs - 60*elapsedMin);
  const int elapsedMs  = floor((elapsed - 3600*elapsedHrs - 60*elapsedMin - elapsedSec)*1000);

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
  const int advHrs = floor(lastadv/3600);
  const int advMin = floor((lastadv - 3600*advHrs)/60);
  const int advSec = floor( lastadv - 3600*advHrs - 60*advMin);
  const int advMs  = floor((lastadv - 3600*advHrs - 60*advMin - advSec)*1000);

  // Write a string with the previous iteration metrics
  sprintf(metrics, 
	  "%31c -- Last time step        : %3.3ih %2.2im %2.2is %3.3ims",
	  ' ',
	  advHrs, 
	  advMin, 
	  advSec, 
	  advMs);
  pout() << metrics << endl;

  // Hours, minutes, seconds and millisecond of the previous iteration
  const Real wt_ns = (m_wallClockTwo - m_wallClockOne)*1.E-9/m_dt;
  const int wt_Hrs = floor(wt_ns/3600);
  const int wt_Min = floor((wt_ns - 3600*wt_Hrs)/60);
  const int wt_Sec = floor( wt_ns - 3600*wt_Hrs - 60*wt_Min);
  const int wt_Ms  = floor((wt_ns - 3600*wt_Hrs - 60*wt_Min - wt_Sec)*1000);
  sprintf(metrics, 
	  "%31c -- Wall time per ns      : %3.3ih %2.2im %2.2is %3.3ims",
	  ' ',
	  wt_Hrs, 
	  wt_Min, 
	  wt_Sec, 
	  wt_Ms);
  pout() << metrics << endl;


  // This is the time remaining
  const Real maxPercent = Max(percentTime, percentStep);
  const Real remaining  = 100.*elapsed/maxPercent - elapsed;
  const int remHrs = floor(remaining/3600);
  const int remMin = floor((remaining - 3600*remHrs)/60);
  const int remSec = floor( remaining - 3600*remHrs - 60*remMin);
  const int remMs  = floor((remaining - 3600*remHrs - 60*remMin - remSec)*1000);

  // Write a string with the previous iteration metrics
  sprintf(metrics, 
	  "%31c -- Estimated remaining   : %3.3ih %2.2im %2.2is %3.3ims",
	  ' ',
	  remHrs, 
	  remMin, 
	  remSec, 
	  remMs);
  pout() << metrics << endl;

  // Write memory usage
#ifdef CH_USE_MEMORY_TRACKING
  const int BytesPerMB = 1024*1024;
  long long curMem;
  long long peakMem;
  overallMemoryUsage(curMem, peakMem);

  pout() << "                                -- Unfreed memory        : " << curMem/BytesPerMB << "(MB)" << endl;
  pout() << "                                -- Peak memory usage     : " << peakMem/BytesPerMB << "(MB)" << endl;

#ifdef CH_MPI
  int unfreed_mem = curMem;
  int peak_mem    = peakMem;

  int max_unfreed_mem;
  int max_peak_mem;

  int result1 = MPI_Allreduce(&unfreed_mem, &max_unfreed_mem, 1, MPI_INT, MPI_MAX, Chombo_MPI::comm);
  int result2 = MPI_Allreduce(&peak_mem,    &max_peak_mem,    1, MPI_INT, MPI_MAX, Chombo_MPI::comm);
  pout() << "                                -- Max unfreed memory    : " << max_unfreed_mem/BytesPerMB << "(MB)" << endl;
  pout() << "                                -- Max peak memory usage : " << max_peak_mem/BytesPerMB << "(MB)" << endl;
#endif
#endif



}

int Driver::getFinestTagLevel(const EBAMRTags& a_cellTags) const{
  CH_TIME("Driver::getFinestTagLevel");
  if(m_verbosity > 5){
    pout() << "Driver::getFinestTagLevel" << endl;
  }

  int finest_tag_level = -1;
  for (int lvl = 0; lvl < a_cellTags.size(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const DenseIntVectSet& tags = (*a_cellTags[lvl])[dit()];

      if(!tags.isEmpty()){
	finest_tag_level = Max(finest_tag_level, lvl);
      }
    }
  }

#ifdef CH_MPI
  int finest;
  MPI_Allreduce(&finest_tag_level, &finest, 1, MPI_INT, MPI_MAX, Chombo_MPI::comm);

  finest_tag_level = finest;
#endif

  return finest_tag_level;
}

bool Driver::tagCells(Vector<IntVectSet>& a_allTags, EBAMRTags& a_cellTags){
  CH_TIME("Driver::tagCells");
  if(m_verbosity > 5){
    pout() << "Driver::tagCells" << endl;
  }

  bool got_new_tags = false;

  // Note that when we regrid we add at most one level at a time. This means that if we have a
  // simulation with AMR depth l and we want to add a level l+1, we need tags on levels 0 through l.
  const int finest_level  = m_amr->getFinestLevel();
  a_allTags.resize(1 + finest_level, IntVectSet());

  if(!m_cellTagger.isNull()){
    got_new_tags = m_cellTagger->tagCells(a_cellTags);
  }


  // Gather tags from a_tags
  for (int lvl = 0; lvl <= finest_level; lvl++){
    for (DataIterator dit = a_cellTags[lvl]->dataIterator(); dit.ok(); ++dit){
      a_allTags[lvl] |= IntVectSet((*a_cellTags[lvl])[dit()]);// This should become a TreeIntVectSet
    }

    // Grow tags with cell taggers buffer
    if(!m_cellTagger.isNull()){
      const int buf = m_cellTagger->getBuffer();
      a_allTags[lvl].grow(buf);
    }
  }

  // Add geometric tags.
  int tag_level = this->getFinestTagLevel(a_cellTags);
  if(m_allowCoarsening){
    for (int lvl = 0; lvl <= finest_level; lvl++){
      if(lvl <= tag_level){
	a_allTags[lvl] |= m_geomTags[lvl];
      }
    }
  }
  else{
    // Loop only goes to the current finest level because we only add one level at a time
    for (int lvl = 0; lvl <= finest_level; lvl++){
      if(lvl < m_amr->getMaxAmrDepth()){ // Geometric tags don't exist on AmrMesh.m_maxAmrDepth
	a_allTags[lvl] |= m_geomTags[lvl];
      }
    }
  }

#if 0 // Debug - if this fails, you have tags on m_amr->m_maxAmrDepth and something has gone wrong. 
  if(finest_level == m_amr->getMaxAmrDepth()){
    for (int lvl = 0; lvl <= finest_level; lvl++){
      pout() << "level = " << lvl << "\t num_pts = " << a_allTags[lvl].numPts() << endl;
    }
    CH_assert(a_allTags[finest_level].isEmpty());
  }
#endif

  // Get the total number of tags
  Vector<int> num_local_tags(1+finest_level);
  for (int lvl = 0; lvl <= finest_level; lvl++){
    num_local_tags[lvl] = a_allTags[lvl].numPts();
  }

  return got_new_tags;
}

void Driver::writeMemoryUsage(){
  CH_TIME("Driver::writeMemoryUsage");
  if(m_verbosity > 3){
    pout() << "Driver::writeMemoryUsage" << endl;
  }

  char file_char[1000];
  const std::string prefix = m_outputDirectory + "/mpi/memory/" + m_outputFileNames;
  sprintf(file_char, "%s.memory.step%07d.%dd.dat", prefix.c_str(), m_timeStep, SpaceDim);
  std::string fname(file_char);

  // Get memory stuff
  Vector<Real> peak, unfreed;
  MemoryReport::getMemoryUsage(peak, unfreed);
  
  // Begin writing output
  if(procID() == 0){
    std::ofstream f;
    f.open(fname, std::ios_base::trunc);
    const int width = 12;

    // Write header
    f << std::left << std::setw(width) << "# MPI rank" << "\t"
      << std::left << std::setw(width) << "Peak memory" << "\t"
      << std::left << std::setw(width) << "Unfreed memory" << "\t"
      << endl;

    // Write memory 
    for (int i = 0; i < numProc(); i++){
      f << std::left << std::setw(width) << i << "\t"
	<< std::left << std::setw(width) << peak[i] << "\t"
	<< std::left << std::setw(width) << unfreed[i] << "\t"
	<< endl;
    }
  }
}

void Driver::writeComputationalLoads(){
  CH_TIME("Driver::writeComputationalLoads");
  if(m_verbosity > 3){
    pout() << "Driver::writeComputationalLoads" << endl;
  }

  const int nProc = numProc();

  // Filename for output. 
  char file_char[1000];
  const std::string prefix = m_outputDirectory + "/mpi/loads/" + m_outputFileNames;
  sprintf(file_char, "%s.loads.step%07d.%dd.dat", prefix.c_str(), m_timeStep, SpaceDim);
  std::string fname(file_char);

  // Get sum of all loads on all Realms
  std::map<std::string, Vector<long int> > RealmLoads;
  for (const auto& r : m_amr->getRealms()){

    // Compute total loads on each rank.
    Vector<long int> sumLoads(nProc, 0L);
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const Vector<long int> boxLoads = m_timeStepper->getCheckpointLoads(r, lvl);

      const DisjointBoxLayout& dbl = m_amr->getGrids(r)[lvl];
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
	sumLoads[procID()] += boxLoads[dit().intCode()];
      }
    }

    // Reduce onto output rank.
#ifdef CH_MPI
    Vector<long int> tmp(nProc, 0L);
    MPI_Allreduce(&(sumLoads[0]), &(tmp[0]), nProc, MPI_LONG, MPI_SUM, Chombo_MPI::comm);
    sumLoads = tmp;
#endif

    RealmLoads.emplace(r, sumLoads);
  }

  // Write header
  if(procID() == 0){
    const int width = 12;
    
    std::ofstream f;
    f.open(fname, std::ios_base::trunc);

    // Write header
    std::stringstream ss;
    ss << std::left << std::setw(width) << "# Rank";
    for (auto r : RealmLoads){
      ss << std::left << std::setw(width) << r.first;
    }
    f << ss.str() << endl;

    // Write data.
    for (int irank = 0; irank < nProc; irank++){
      std::stringstream ds;

      ds << std::left << std::setw(width) << irank;
      for (auto r : RealmLoads){
	ds << std::left << std::setw(width) << r.second[irank];
      }
      f << ds.str() << std::endl;
    }

    
    f.close();
  }
}

void Driver::writeGeometry(){
  CH_TIME("Driver::writeGeometry");
  if(m_verbosity > 3){
    pout() << "Driver::writeGeometry" << endl;
  }

  const int ncomp = 2;

  EBAMRCellData output;
  m_amr->allocate(output, m_realm, phase::gas, ncomp);
  DataOps::setValue(output, 0.0);

  // Names
  Vector<std::string> names(2);
  names[0] = "levelset_1";
  names[1] = "levelset_2";

  // Write levelsets
  int icomp = 0;
  this->writeLevelset(output, icomp);

  //
  const int finest_level                 = m_amr->getFinestLevel();
  const Vector<DisjointBoxLayout>& grids = m_amr->getGrids(m_realm);
  const Vector<ProblemDomain>& domains   = m_amr->getDomains();
  const Vector<Real>& dx                 = m_amr->getDx();
  const Vector<int>& ref_rat             = m_amr->getRefinementRatios();

  Vector<LevelData<EBCellFAB>*> output_ptr(1 + finest_level);
  m_amr->alias(output_ptr, output);

  bool replace_covered = false;
  Vector<Real> covered_values;

  // Dummy file name
  char file_char[1000];
  const std::string prefix = m_outputDirectory + "/geo/" + m_outputFileNames;
  sprintf(file_char, "%s.geometry.%dd.hdf5", prefix.c_str(), SpaceDim);
  string fname(file_char);

  writeEBHDF5(fname, 
	      grids,
	      output_ptr,
	      names, 
	      domains[0],
	      dx[0], 
	      m_dt,
	      m_time,
	      ref_rat,
	      finest_level + 1,
	      replace_covered,
	      covered_values,
	      m_numPlotGhost*IntVect::Unit);
}

void Driver::writePlotFile(){
  CH_TIME("Driver::writePlotFile");
  if(m_verbosity > 3){
    pout() << "Driver::writePlotFile" << endl;
  }

  // Filename
  char file_char[1000];
  const std::string prefix = m_outputDirectory + "/plt/" + m_outputFileNames;
  sprintf(file_char, "%s.step%07d.%dd.hdf5", prefix.c_str(), m_timeStep, SpaceDim);
  string fname(file_char);

  this->writePlotFile(fname);

}

void Driver::writeRegridFile(){
  CH_TIME("Driver::writeRegridFile");
  if(m_verbosity > 3){
    pout() << "Driver::writeRegridFile" << endl;
  }

  // Filename
  char file_char[1000];
  const std::string prefix = m_outputDirectory + "/regrid/" + m_outputFileNames;
  sprintf(file_char, "%s.regrid%07d.%dd.hdf5", prefix.c_str(), m_timeStep, SpaceDim);
  string fname(file_char);

  this->writePlotFile(fname);
}

void Driver::writeRestartFile(){
  CH_TIME("Driver::writeRestartFile");
  if(m_verbosity > 3){
    pout() << "Driver::writeRestartFile" << endl;
  }

  // Filename
  char file_char[1000];
  const std::string prefix = m_outputDirectory + "/restart/" + m_outputFileNames;
  sprintf(file_char, "%s.restart%07d.%dd.hdf5", prefix.c_str(), m_timeStep, SpaceDim);
  string fname(file_char);

  this->writePlotFile(fname);
}

void Driver::writeCrashFile(){
  CH_TIME("Driver::writeCrashFile");
  if(m_verbosity > 3){
    pout() << "Driver::writeCrashFile" << endl;
  }

  // Filename
  char file_char[1000];
  const std::string prefix = m_outputDirectory + "/crash/" + m_outputFileNames;
  sprintf(file_char, "%s.crash%07d.%dd.hdf5", prefix.c_str(), m_timeStep, SpaceDim);
  string fname(file_char);

  this->writePlotFile(fname);
}

void Driver::writePlotFile(const std::string a_filename){
  CH_TIME("Driver::writePlotFile");
  if(m_verbosity > 3){
    pout() << "Driver::writePlotFile" << endl;
  }

  // Output file
  EBAMRCellData output;
  EBAMRCellData scratch;

  // Names for output variables  
  Vector<std::string> names(0);

  // Get total number of components for output
  int ncomp = m_timeStepper->getNumberOfPlotVariables();
  if(!m_cellTagger.isNull()) {
    ncomp += m_cellTagger->getNumberOfPlotVariables();
  }
  ncomp += this->getNumberOfPlotVariables();

  // Allocate storage
  m_amr->allocate(output,  m_realm, phase::gas, ncomp);
  m_amr->allocate(scratch, m_realm, phase::gas, 1);
  DataOps::setValue(output, 0.0);
  DataOps::setValue(scratch, 0.0);

  // Assemble data
  int icomp = 0;             // Used as reference for output components
  Real t_assemble = -MPI_Wtime();
  if(m_verbosity >= 3){
    pout() << "Driver::writePlotFile - assembling data..." << endl;
  }
  
  // Time stepper writes its data
  m_timeStepper->writePlotData(output, names, icomp);

  // Cell tagger writes data
  if(!m_cellTagger.isNull()){
    m_cellTagger->writePlotData(output, names, icomp);
  }
									       
  // Data file aliasing, because Chombo IO wants dumb pointers. 
  Vector<LevelData<EBCellFAB>* > output_ptr(1 + m_amr->getFinestLevel());
  m_amr->alias(output_ptr, output);

  // Restrict plot depth if need be
  int plot_depth;
  if(m_maxPlotDepth < 0){
    plot_depth = m_amr->getFinestLevel();
  }
  else{
    plot_depth = Min(m_maxPlotDepth, m_amr->getFinestLevel());
  }

  // Interpolate ghost cells. This might be important if we use multiple Realms. 
  for (int icomp = 0; icomp < ncomp; icomp++){
    const Interval interv(icomp, icomp);
    
    for (int lvl = 1; lvl <= m_amr->getFinestLevel(); lvl++){

      LevelData<EBCellFAB> fineAlias;
      LevelData<EBCellFAB> coarAlias;

      aliasLevelData(fineAlias, output_ptr[lvl],   interv);
      aliasLevelData(coarAlias, output_ptr[lvl-1], interv);

      m_amr->interpGhost(fineAlias, coarAlias, lvl, m_realm, phase::gas);
    }
  }


  // Write internal data
  names.append(this->getPlotVariableNames());
  this->writePlotData(output, icomp);
  t_assemble += MPI_Wtime();


  // Write HDF5 file
  if(m_verbosity >= 3){
    pout() << "Driver::writePlotFile - writing plot file..." << endl;
  }
  Real t_write = -MPI_Wtime();

  // Write. 
  writeEBHDF5(a_filename, 
	      m_amr->getGrids(m_realm),
	      output_ptr,
	      names, 
	      m_amr->getDomains()[0],
	      m_amr->getDx()[0], 
	      m_dt,
	      m_time,
	      m_amr->getRefinementRatios(),
	      plot_depth + 1,
	      false,
	      Vector<Real>(),
	      m_numPlotGhost*IntVect::Unit);
  t_write += MPI_Wtime();

  const Real t_tot = t_write + t_assemble;
  if(m_verbosity >= 3){
    pout() << "Driver::writePlotFile - writing plot file... DONE!. " << endl
	   << "\t Total time    = " << t_tot << " seconds" << endl
	   << "\t Assemble data = " << 100.*t_assemble/t_tot << "%" << endl
	   << "\t Write time    = " << 100.*t_write/t_tot << "%" << endl;
  }
}

void Driver::writePlotData(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("Driver::writePlotData");
  if(m_verbosity > 3){
    pout() << "Driver::writePlotData" << endl;
  }

  if(m_plotTags)     writeTags(a_output, a_comp);
  if(m_plotRanks)    writeRanks(a_output, a_comp);
  if(m_plotLevelset) writeLevelset(a_output, a_comp);
}

void Driver::writeTags(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("Driver::writeTags");
  if(m_verbosity > 3){
    pout() << "Driver::writeTags" << endl;
  }

  
  // Alloc some temporary storage
  EBAMRCellData tags;
  m_amr->allocate(tags, m_realm, phase::gas, 1);
  DataOps::setValue(tags, 0.0);
    
  // Set tagged cells = 1
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, phase::gas)[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const DenseIntVectSet& ivs = (*m_tags[lvl])[dit()];
      const Box box              = dbl.get(dit());

      // Do regular cells only.
      BaseFab<Real>& tags_fab = (*tags[lvl])[dit()].getSingleValuedFAB();
      for (BoxIterator bit(box); bit.ok(); ++bit){
	const IntVect iv = bit();
	if(ivs[iv]){
	  tags_fab(iv, 0) = 1.0;
	}
      }
    }
  }

  DataOps::setCoveredValue(tags, 0, 0.0);

  const Interval src_interv(0, 0);
  const Interval dst_interv(a_comp, a_comp);
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    tags[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
  }

  a_comp++;
}

void Driver::writeRanks(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("Driver::writeRanks");
  if(m_verbosity > 3){
    pout() << "Driver::writeRanks" << endl;
  }

  for (const auto& r : m_amr->getRealms()){
    EBAMRCellData scratch;
    m_amr->allocate(scratch, r, phase::gas, 1);

    DataOps::setValue(scratch, 1.0*procID());

    const Interval src(0,0);
    const Interval dst(a_comp, a_comp);
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      scratch[lvl]->copyTo(src, *a_output[lvl], dst);
    }

    a_comp++;
  }
}

void Driver::writeLevelset(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("Driver::writeLevelset");
  if(m_verbosity > 3){
    pout() << "Driver::writeLevelset" << endl;
  }

  const RefCountedPtr<BaseIF>& lsf1 = m_computationalGeometry->getGasImplicitFunction();
  const RefCountedPtr<BaseIF>& lsf2 = m_computationalGeometry->getSolidImplicitFunction();

  const RealVect prob_lo = m_amr->getProbLo();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(a_output.getRealm())[lvl];
    const Real dx = m_amr->getDx()[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      FArrayBox& fab = (*a_output[lvl])[dit()].getFArrayBox();

      fab.setVal(0.0, a_comp);
      fab.setVal(0.0, a_comp+1);

      const Box box          = fab.box();

      for (BoxIterator bit(box); bit.ok(); ++bit){
	const IntVect iv = bit();
	
	const RealVect pos = prob_lo + (RealVect(iv)+ 0.5*RealVect::Unit)*dx;

	if(!lsf1.isNull()) fab(iv, a_comp  ) = lsf1->value(pos);
	if(!lsf2.isNull()) fab(iv, a_comp+1) = lsf2->value(pos);
      }
    }
  }

  a_comp = a_comp + 2;
}

void Driver::writeCheckpointFile(){
  CH_TIME("Driver::writeCheckpointFile");
  if(m_verbosity > 3){
    pout() << "Driver::writeCheckpointFile" << endl;
  }
  
  const int finest_level = m_amr->getFinestLevel();
  int finest_chk_level  = Min(m_maxCheckpointDepth, finest_level);
  if(m_maxCheckpointDepth < 0){
    finest_chk_level = finest_level;
  }

  // Write header. 
  HDF5HeaderData header;
  header.m_real["coarsest_dx"] = m_amr->getDx()[0];
  header.m_real["time"]        = m_time;
  header.m_real["dt"]          = m_dt;
  header.m_int["step"]         = m_timeStep;
  header.m_int["finest_level"] = finest_level;

  // Write Realm names.
  for (auto r : m_amr->getRealms()){
    header.m_string[r] = r;
  }

  // Output file name
  char str[100];
  const std::string prefix = m_outputDirectory + "/chk/" + m_outputFileNames;
  sprintf(str, "%s.check%07d.%dd.hdf5", prefix.c_str(), m_timeStep, SpaceDim);

  // Output file
  HDF5Handle handle_out(str, HDF5Handle::CREATE);
  header.writeToFile(handle_out);

  // Write stuff level by level
  const Real t0 = MPI_Wtime();
  if(m_verbosity >= 3){
    pout() << "Driver::writeCheckpointFile - writing checkpoint file..." << endl;
  }
  
  for (int lvl = 0; lvl <= finest_chk_level; lvl++){
    handle_out.setGroupToLevel(lvl);

    // write amr grids
    write(handle_out, m_amr->getGrids(m_realm)[lvl]); // write AMR grids

    // time stepper checkpoints data
    m_timeStepper->writeCheckpointData(handle_out, lvl); 

    // Driver checkpoints internal data
    this->writeCheckpointLevel(handle_out, lvl); 
  }
  const Real t1 = MPI_Wtime();

  if(m_verbosity >= 3){
    pout() << "Driver::writeCheckpointFile - writing checkpoint file... DONE! " << endl
	   << "\t Total time    = " << t1 - t0 << " seconds" << endl;
  }
  
  handle_out.close();
}

void Driver::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level){
  CH_TIME("Driver::writeCheckpointLevel");
  if(m_verbosity > 5){
    pout() << "Driver::writeCheckpointLevel" << endl;
  }

  this->writeCheckpointTags(a_handle, a_level);
  this->writeCheckpointRealmLoads(a_handle, a_level);

}

void Driver::writeCheckpointTags(HDF5Handle& a_handle, const int a_level){
  CH_TIME("Driver::writeCheckpointTags");
  if(m_verbosity > 5){
    pout() << "Driver::writeCheckpointTags" << endl;
  }

  // Create some scratch data = 0 which can grok
  EBCellFactory fact(m_amr->getEBISLayout(m_realm, phase::gas)[a_level]);
  LevelData<EBCellFAB> scratch(m_amr->getGrids(m_realm)[a_level], 1, 3*IntVect::Unit, fact);
  DataOps::setValue(scratch, 0.0);

  // Set tags = 1
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_level];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, phase::gas)[a_level];
    
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box = dbl.get(dit());
    const DenseIntVectSet& tags = (*m_tags[a_level])[dit()];

    BaseFab<Real>& scratch_fab = scratch[dit()].getSingleValuedFAB();
    for (BoxIterator bit(box); bit.ok(); ++bit){
      const IntVect iv = bit();
      if(tags[iv]){
	scratch_fab(iv, 0) = 1.0;
      }
    }

    DataOps::setCoveredValue(scratch, 0, 0.0);
  }

  // Write tags
  write(a_handle, scratch, "tagged_cells");
}

void Driver::writeCheckpointRealmLoads(HDF5Handle& a_handle, const int a_level){
  CH_TIME("Driver::writeCheckpointRealmLoads");
  if(m_verbosity > 5){
    pout() << "Driver::writeCheckpointRealmLoads" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_level];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, phase::gas)[a_level];

  // Make some storage. 
  LevelData<FArrayBox> scratch(dbl, 1, 3*IntVect::Unit);

  // Get loads. 
  for (auto r : m_amr->getRealms()){
    const Vector<long int> loads = m_timeStepper->getCheckpointLoads(r, a_level);

    // Set loads on an FArrayBox
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      scratch[dit()].setVal(loads[dit().intCode()]);
    }

    // String identifier in HDF file.
    const std::string str = r + "_loads";

    // Write
    write(a_handle, scratch, str);
  }
}

void Driver::readCheckpointFile(const std::string& a_restartFile){
  CH_TIME("Driver::readCheckpointFile");
  if(m_verbosity > 3){
    pout() << "Driver::readCheckpointFile" << endl;
  }

  // Time stepper can register Realms immediately. 
  m_timeStepper->registerRealms();

  // Read the header that was written by new_readCheckpointFile
  HDF5Handle handle_in(a_restartFile, HDF5Handle::OPEN_RDONLY);
  HDF5HeaderData header;
  header.readFromFile(handle_in);

  m_time        = header.m_real["time"];
  m_dt          = header.m_real["dt"];
  m_timeStep        = header.m_int["step"];

  const Real coarsest_dx = header.m_real["coarsest_dx"];
  const int base_level   = 0;
  const int finest_level = header.m_int["finest_level"];

  // Get the names of the Realms that were checkpointed. This is a part of the HDF header. 
  std::map<std::string, Vector<Vector<long int > > > chk_loads;
  for (auto s : header.m_string){
    chk_loads.emplace(s.second, Vector<Vector<long int> >());
  }

  // Get then names of the Realms that will be used for simulations.
  std::map<std::string, Vector<Vector<long int > > > sim_loads;
  for (const auto& iRealm : m_amr->getRealms()){
    sim_loads.emplace(iRealm, Vector<Vector<long int> >());
  }

  // Print checkpointed Realm names. 
  if(m_verbosity > 2){
    pout() << "Driver::readCheckpointFile - checked Realms are: ";
    for (auto r : chk_loads){
      pout() << '"' << r.first << '"' << "\t";
    }
    pout() << endl;
  }

  // Abort if base resolution has changed. 
  if(!coarsest_dx == m_amr->getDx()[0]){
    MayDay::Abort("Driver::readCheckpointFile - coarsest_dx != dx[0], did you change the base level resolution?!?");
  }

  // Read in grids. If the file has no grids we must abort. 
  Vector<Vector<Box> > boxes(1 + finest_level);
  for (int lvl = 0; lvl <= finest_level; lvl++){
    handle_in.setGroupToLevel(lvl);
    
    const int status = read(handle_in, boxes[lvl]);
    
    if(status != 0) {
      MayDay::Error("Driver::readCheckpointFile - file has no grids");
    }
  }

  // Read in the computational loads from the HDF5 file. 
  for (auto& r : chk_loads){
    const std::string& Realm_name = r.first;
    Vector<Vector<long int> >& Realm_loads = r.second;

    Realm_loads.resize(1 + finest_level);
    for (int lvl = 0; lvl <= finest_level; lvl++){
      Realm_loads[lvl].resize(boxes[lvl].size(), 0L);

      this->readCheckpointRealmLoads(Realm_loads[lvl], handle_in, Realm_name, lvl);
    }
  }


  // In case we restart with more or fewer Realms, we need to decide how to assign the computational loads. If the Realm was a new Realm we may
  // not have the computational loads for that. In that case we take the computational loads from the primal Realm. 
  for (auto& s : sim_loads){

    const std::string&         cur_Realm = s.first;
    Vector<Vector<long int> >& cur_loads = s.second;

#if 0 // Original code
    for (const auto& c : chk_loads){
      if(cur_Realm == c.first){
	cur_loads = c.second;
      }
      else{
	cur_loads = chk_loads.at(Realm::Primal);
      }
    }
#else

    bool found_checked_loads = false;
    for (const auto& c : chk_loads){
      if(cur_Realm == c.first){
	found_checked_loads = true;
      }
    }

    cur_loads = (found_checked_loads) ? chk_loads.at(cur_Realm) : chk_loads.at(Realm::Primal);
#endif
  }

  // Define AmrMesh and Realms. 
  m_amr->setFinestLevel(finest_level); 
  m_amr->setGrids(boxes, sim_loads);
  
  // Instantiate solvers and register operators
  m_timeStepper->setupSolvers();
  m_timeStepper->registerOperators();
  m_timeStepper->synchronizeSolverTimes(m_timeStep, m_time, m_dt);  
  m_amr->regridOperators(base_level);
  m_timeStepper->allocate();

  // Allocate internal stuff (e.g. space for tags)
  this->allocateInternals();

  // Go through level by level and have solvers extract their data
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    handle_in.setGroupToLevel(lvl);

    // time stepper reads in data
    m_timeStepper->readCheckpointData(handle_in, lvl);

    // Read in internal data
    readCheckpointLevel(handle_in, lvl);
  }

  // Close input file
  handle_in.close();
}

void Driver::readCheckpointLevel(HDF5Handle& a_handle, const int a_level){
  CH_TIME("Driver::readCheckpointLevel");
  if(m_verbosity > 5){
    pout() << "Driver::readCheckpointLevel" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_level];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, phase::gas)[a_level];

  // Some scratch data we can use
  EBCellFactory fact(ebisl);
  LevelData<EBCellFAB> scratch(dbl, 1, 3*IntVect::Unit, fact);
  DataOps::setValue(scratch, 0.0);

  // Read in tags
  read<EBCellFAB>(a_handle, scratch, "tagged_cells", dbl, Interval(0,0), false);

  // Instantiate m_tags
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box          = dbl.get(dit());
    const EBCellFAB& tmp   = scratch[dit()];
    const EBISBox& ebisbox = tmp.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs(box);

    DenseIntVectSet& tagged_cells = (*m_tags[a_level])[dit()];

    BaseFab<Real>& scratch_fab = scratch[dit()].getSingleValuedFAB();
    for (BoxIterator bit(box); bit.ok(); ++bit){
      const IntVect iv = bit();
      if(scratch_fab(iv, 0) > 0.9999){
	tagged_cells |= iv;
      }
    }
  }
}

void Driver::readCheckpointRealmLoads(Vector<long int>& a_loads, HDF5Handle& a_handle, const std::string a_realm, const int a_level){
  CH_TIME("Driver::readCheckpointRealmLoads(...)");
  if(m_verbosity > 5){
    pout() << "Driver::readCheckpointRealmLoads(...)" << endl;
  }

  // HDF identifier.
  const std::string str = a_realm + "_loads";

  // Read into an FArrayBox.
  FArrayBox fab;
  for (int ibox = 0; ibox < a_loads.size(); ibox++){
    readFArrayBox(a_handle, fab, a_level, ibox, Interval(0,0), str);

    a_loads[ibox] = lround(fab.max());
  }
}

void Driver::writeVectorData(HDF5HeaderData&     a_header,
			     const Vector<Real>& a_data,
			     const std::string   a_name,
			     const int           a_elements){
  CH_TIME("Driver::writeVectorData");
  if(m_verbosity > 3){
    pout() << "Driver::writeVectorData" << endl;
  }

  char step[100];

  const int last = a_data.size() < a_elements ? a_data.size() : a_elements;
  for (int i = 0; i < last; i++){
    sprintf(step, "%07d", i);

    const std::string identifier = a_name + std::string(step);
    a_header.m_real[identifier] = a_data[i];
  }
}

void Driver::readVectorData(HDF5HeaderData& a_header,
			    Vector<Real>&         a_data,
			    const std::string     a_name,
			    const int             a_elements){
  CH_TIME("Driver::readVectorData");
  if(m_verbosity > 3){
    pout() << "Driver::readVectorData" << endl;
  }

  char step[100];

  const int last = a_data.size() < a_elements ? a_data.size() : a_elements;
  for (int i = 0; i < last; i++){
    sprintf(step, "%07d", i);

    const std::string identifier = a_name + std::string(step);
    a_data[i] = a_header.m_real[identifier];
  }
}

#include <CD_NamespaceFooter.H>
