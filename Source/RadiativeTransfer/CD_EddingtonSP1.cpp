/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_EddingtonSP1.cpp
  @brief  Implementation of CD_EddingtonSP1.H
  @author Robert Marskar
*/

// Std includes
#include <chrono>
#include <time.h>

// Chombo includes
#include <ParmParse.H>
#include <EBAMRIO.H>
#include <BRMeshRefine.H>

// Our includes
#include <CD_EddingtonSP1.H>
#include <CD_DataOps.H>
#include <CD_Units.H>
#include <CD_ConductivityDomainBcWrapper.H>
#include <CD_NamespaceHeader.H>

#define EddingtonSP1_feature 1 // Comment Feb. 14 2018: I think we can keep this - it appears to produce the correct physics.

EddingtonSP1::EddingtonSP1() : RtSolver() {
  m_name = "EddingtonSP1";
  m_className = "EddingtonSP1";

  m_verbosity  = -1;
  m_needsMultigridSetup = true;
  m_hasDeeperMultigridLevels = false;
}

EddingtonSP1::~EddingtonSP1(){
}

int EddingtonSP1::queryGhost() const{
  return 3;
}

void EddingtonSP1::parseOptions(){
  CH_TIME("EddingtonSP1::parseOptions");
  if(m_verbosity > 5){
    pout() << m_name + "::parseOptions" << endl;
  }
  
  parseDomainBc();          // Parses domain BC options
  parseStationary();        // Parse stationary solver
  parsePlotVariables();     // Parses plot variables
  parseMultigridSettings(); // Parses solver parameters for geometric multigrid
}

void EddingtonSP1::parseRuntimeOptions(){
  CH_TIME("EddingtonSP1::parseRuntimeOptions");
  if(m_verbosity > 5){
    pout() << m_name + "::parseRuntimeOptions" << endl;
  }
  
  parsePlotVariables();     // Parses plot variables
  parseMultigridSettings(); // Parses solver parameters for geometric multigrid
}

void EddingtonSP1::parseDomainBc(){
  CH_TIME("EddingtonSP1::parseDomainBc");
  if(m_verbosity > 5){
    pout() << m_name + "::parseDomainBc" << endl;
  }

  allocateWallBc();

  // Get BC from input script
  ParmParse pp(m_className.c_str());
  
  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const Side::LoHiSide side = sit();
	
      std::string str_dir;
      if(dir == 0){
	str_dir = "x";
      }
      else if(dir == 1){
	str_dir = "y";
      }
      else if(dir == 2){
	str_dir = "z";
      }

      if(side == Side::Lo){
	std::string type;
	std::string bc_string = "bc_" + str_dir + "_low";
	if(pp.contains(bc_string.c_str())){
	  pp.get(bc_string.c_str(), type);
	  if(type == "neumann"){
	    this->setNeumannWallBc(dir, Side::Lo, 0.0);
	  }
	  else if(type == "robin"){
	    this->setRobinWallBc(dir, Side::Lo, 0.0);
	  }
	  else {
	    std::string error = "EddingtonSP1::EddingtonSP1 - unknown bc requested for " + bc_string;
	    MayDay::Abort(error.c_str());
	  }
	}
      }
      else if(side == Side::Hi){
	std::string type;
	std::string bc_string = "bc_" + str_dir + "_high";
	if(pp.contains(bc_string.c_str())){
	  pp.get(bc_string.c_str(), type);
	  if(type == "neumann"){
	    this->setNeumannWallBc(dir, Side::Hi, 0.0);
	  }
	  else if(type == "robin"){
	    this->setRobinWallBc(dir, Side::Hi, 0.0);
	  }
	  else {
	    std::string error = "EddingtonSP1::EddingtonSP1 - unknown bc requested for " + bc_string;
	    MayDay::Abort(error.c_str());
	  }
	}
      }
    }
  }
}

void EddingtonSP1::parseStationary(){
  CH_TIME("EddingtonSP1::parseStationary");
  if(m_verbosity > 5){
    pout() << m_name + "::parseStationary" << endl;
  }

  ParmParse pp(m_className.c_str());
  std::string str;
  
  pp.get("stationary", m_stationary);
  pp.get("use_tga",    m_useTGA);
}

void EddingtonSP1::parsePlotVariables(){
  CH_TIME("EddingtonSP1::parsePlotVariables");
  if(m_verbosity > 5){
    pout() << m_name + "::parsePlotVariables" << endl;
  }

  m_plotPhi    = false;
  m_plotSource = false;

  ParmParse pp(m_className.c_str());
  const int num = pp.countval("plt_vars");
  Vector<std::string> str(num);
  pp.getarr("plt_vars", str, 0, num);

  for (int i = 0; i < num; i++){
    if(     str[i] == "phi") m_plotPhi    = true;
    else if(str[i] == "src") m_plotSource = true;
  }
}

void EddingtonSP1::parseMultigridSettings(){
  ParmParse pp(m_className.c_str());

  std::string str;

  pp.get("gmg_verbosity",   m_multigridVerbosity);
  pp.get("gmg_coarsen",     m_numCoarseningsBeforeAggregation);
  pp.get("gmg_pre_smooth",  m_multigridPreSmooth);
  pp.get("gmg_post_smooth", m_multigridPostSmooth);
  pp.get("gmg_bott_smooth", m_multigridBottomSmooth);
  pp.get("gmg_max_iter",    m_multigridMaxIterations);
  pp.get("gmg_min_iter",    m_multigridMinIterations);
  pp.get("gmg_tolerance",   m_multigridTolerance);
  pp.get("gmg_hang",        m_multigridHang);
  pp.get("gmg_bottom_drop", m_numCellsBottomDrop);

  // Bottom solver
  pp.get("gmg_bottom_solver", str);
  if(str == "simple"){
    m_bottomSolver = BottomSolver::Simple;
  }
  else if(str == "bicgstab"){
    m_bottomSolver = BottomSolver::BiCGStab;
  }
  else{
    MayDay::Abort("EddingtonSP1::parseMultigridSettings - unknown bottom solver requested");
  }

  // Relaxation type
  pp.get("gmg_relax_type", str);
  if(str == "gsrb"){
    m_multigridRelaxMethod = RelaxationMethod::GSRBFast;
  }
  else if(str == "jacobi"){
    m_multigridRelaxMethod = RelaxationMethod::Jacobi;
  }
  else if(str == "gauss_seidel"){
    m_multigridRelaxMethod = RelaxationMethod::GaussSeidel;
  }
  else{
    MayDay::Abort("EddingtonSP1::parseMultigridSettings - unknown relaxation method requested");
  }

  // Cycle type
  pp.get("gmg_cycle", str);
  if(str == "vcycle"){
    m_multigridType = MultigridType::VCycle;
  }
  else{
    MayDay::Abort("EddingtonSP1::parseMultigridSettings - unknown cycle type requested");
  }

  // No lower than 2. 
  if(m_numCellsBottomDrop < 2){
    m_numCellsBottomDrop = 2;
  }
}

void EddingtonSP1::parseReflection(){
  CH_TIME("EddingtonSP1::parse_reflectivity");
  if(m_verbosity > 5){
    pout() << m_name + "::parse_reflectivity" << endl;
  }
  
  ParmParse pp(m_className.c_str());

  Real r;
  pp.get("reflectivity", r);

  m_reflectionCoefficientOne = r/(2.0);
  m_reflectionCoefficientTwo = r/(3.0);
}

void EddingtonSP1::allocateWallBc(){
  CH_TIME("EddingtonSP1::allocateWallBc");
  if(m_verbosity > 5){
    pout() << "EddingtonSP1::allocateWallBc" << endl;
  }
  m_wallBc.resize(2*SpaceDim);
  for (int i = 0; i < 2*SpaceDim; i++){
    m_wallBc[i] = RefCountedPtr<WallBc> (NULL);
  }
}

void EddingtonSP1::preRegrid(const int a_base, const int a_oldFinestLevel){
  CH_TIME("EddingtonSP1::preRegrid");
  if(m_verbosity > 5){
    pout() << m_name + "::preRegrid" << endl;
  }

  const int ncomp = 1;
  const int finest_level = m_amr->getFinestLevel();

  m_amr->allocate(m_cache, m_realm, m_phase, ncomp);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    m_phi[lvl]->localCopyTo(*m_cache[lvl]);
  }
}

void EddingtonSP1::setReflectionCoefficients(const Real a_r1, const Real a_r2){
  CH_TIME("EddingtonSP1::setReflectionCoefficients");
  if(m_verbosity > 5){
    pout() << m_name + "::setReflectionCoefficients" << endl;
  }
  
  m_reflectionCoefficientOne = a_r1;
  m_reflectionCoefficientTwo = a_r2;
}

void EddingtonSP1::allocateInternals(){
  CH_TIME("EddingtonSP1::allocateInternals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocateInternals" << endl;
  }
  
  const int ncomp = 1;

  m_amr->allocate(m_aCoef,       m_realm, m_phase, ncomp);
  m_amr->allocate(m_bco,       m_realm, m_phase, ncomp);
  m_amr->allocate(m_bco_irreg, m_realm, m_phase, ncomp);
  m_amr->allocate(m_phi,     m_realm, m_phase, ncomp);
  m_amr->allocate(m_source,    m_realm, m_phase, ncomp);
  m_amr->allocate(m_resid,     m_realm, m_phase, ncomp);

  DataOps::setValue(m_resid,  0.0);
  DataOps::setValue(m_phi,  0.0);
  DataOps::setValue(m_source, 0.0);

  this->setACoefAndBCoef();
}

void EddingtonSP1::deallocateInternals(){
  m_amr->deallocate(m_aCoef);
  m_amr->deallocate(m_bco);
  m_amr->deallocate(m_bco_irreg);
  m_amr->deallocate(m_phi);
  m_amr->deallocate(m_source);
  m_amr->deallocate(m_resid);
}

void EddingtonSP1::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) {
  CH_TIME("EddingtonSP1::regrid");
  if(m_verbosity > 5){
    pout() << m_name + "::regrid" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const Interval interv(comp, comp);

  this->allocateInternals();

  Vector<RefCountedPtr<EBPWLFineInterp> >& interpolator = m_amr->getPwlInterpolator(m_realm, m_phase);

  // These levels have not changed
  for (int lvl = 0; lvl <= Max(0, a_lmin-1); lvl++){
    m_cache[lvl]->copyTo(*m_phi[lvl]); // Base level should never change, but ownership might
  }

  // These levels have changed
  for (int lvl = Max(1,a_lmin); lvl <= a_newFinestLevel; lvl++){
    interpolator[lvl]->interpolate(*m_phi[lvl], *m_phi[lvl-1], interv);

    if(lvl <= a_oldFinestLevel){
      m_cache[lvl]->copyTo(*m_phi[lvl]);
    }
  }

  m_needsMultigridSetup = true;
}

void EddingtonSP1::registerOperators(){
  CH_TIME("EddingtonSP1::registerOperators");
  if(m_verbosity > 5){
    pout() << m_name + "::registerOperators" << endl;
  }

  if(m_amr.isNull()){
    MayDay::Abort("EddingtonSP1::registerOperators - need to set AmrMesh!");
  }
  else{
    m_amr->registerOperator(s_eb_coar_ave,     m_realm, m_phase);
    m_amr->registerOperator(s_eb_fill_patch,   m_realm, m_phase);
    m_amr->registerOperator(s_eb_flux_reg,     m_realm, m_phase);
    m_amr->registerOperator(s_eb_quad_cfi,     m_realm, m_phase);
    m_amr->registerOperator(s_eb_gradient,     m_realm, m_phase);
    m_amr->registerOperator(s_eb_irreg_interp, m_realm, m_phase);
  }
}

bool EddingtonSP1::advance(const Real a_dt, EBAMRCellData& a_phi, const EBAMRCellData& a_source, const bool a_zerophi){
  CH_TIME("EddingtonSP1::advance(ebamrcell, ebamrcell)");
  if(m_verbosity > 5){
    pout() << m_name + "::advance(ebamrcell, ebamrcell)" << endl;
  }

  const int ncomp        = 1;
  const int finest_level = m_amr->getFinestLevel();

  if(m_needsMultigridSetup){
    this->setupMultigrid();
  }

  bool converged;

  // Must have a dummy for checking initial residual
  EBAMRCellData dummy;
  EBAMRCellData source;
  m_amr->allocate(dummy,  m_realm, m_phase, ncomp);
  m_amr->allocate(source, m_realm, m_phase, ncomp);
  DataOps::setValue(dummy, 0.0);

  // Various source term manipulations. 
  DataOps::setValue(source, 0.0);
  DataOps::incr(source, a_source, 1.0);
#if EddingtonSP1_feature
  DataOps::scale(source, 1./Units::c); // Source should be scaled by 1./c0
#endif
  if(m_stationary){ // Should kappa-scale for transient solvres
    DataOps::kappaScale(source);
  }

  Vector<LevelData<EBCellFAB>* > phi, rhs, res, zero;
  m_amr->alias(phi,  a_phi);
  m_amr->alias(rhs,  source);
  m_amr->alias(res,  m_resid);
  m_amr->alias(zero, dummy);

  if(m_stationary){
    const Real phi_resid  = m_multigridSolver->computeAMRResidual(phi,  rhs, finest_level, 0); // Incoming residual
    const Real zero_resid = m_multigridSolver->computeAMRResidual(zero, rhs, finest_level, 0); // Zero residual

    if(phi_resid > zero_resid*m_multigridTolerance){ // Residual is too large
      m_multigridSolver->m_convergenceMetric = zero_resid;
      m_multigridSolver->solveNoInitResid(phi, res, rhs, finest_level, 0, a_zerophi);
      
      const int status = m_multigridSolver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
      if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
	converged = true;
      }
    }
    else{ // Solution is already good enough
      converged = true;
    }
    m_multigridSolver->revert(phi, rhs, finest_level, 0);

    DataOps::setCoveredValue(a_phi, 0, 0.0);
  }
  else{
    if(m_useTGA){
#if EddingtonSP1_feature
      m_tgaSolver->oneStep(res, phi, rhs, Units::c*a_dt, 0, finest_level, m_time);
#else
      m_tgaSolver->oneStep(res, phi, rhs, a_dt, 0, finest_level, m_time);
#endif
    }
    else{
#if EddingtonSP1_feature
      m_eulerSolver->oneStep(res, phi, rhs, Units::c*a_dt, 0, finest_level, false);
#else
      m_eulerSolver->oneStep(res, phi, rhs, a_dt, 0, finest_level, false);
#endif
    }
    
    const int status = m_multigridSolver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }

    // We solve onto res, copy back to state
    DataOps::copy(a_phi, m_resid);
  }

  m_amr->averageDown(a_phi, m_realm, m_phase);
  m_amr->interpGhost(a_phi, m_realm, m_phase);

  DataOps::floor(a_phi, 0.0);

  return converged;
}

void EddingtonSP1::setupMultigrid(){
  CH_TIME("EddingtonSP1::setupMultigrid");
  if(m_verbosity > 5){
    pout() << m_name + "::setupMultigrid" << endl;
  }

  if(!m_hasDeeperMultigridLevels){
    this->defineDeeperMultigridLevels();
    m_hasDeeperMultigridLevels = true;
  }
  
  this->setMultigridCoefficients();       // Set coefficients, kappa, aco, bco
  this->setupOperatorFactory(); // Set the operator factory
  this->setupMultigridSolver();        // Set up the AMR multigrid solver

  if(!m_stationary){
    if(m_useTGA){
      this->setupTGA();
    }
    else{
      this->setupEuler();
    }
  }

  m_needsMultigridSetup = false;
}

void EddingtonSP1::setMultigridCoefficients(){
  CH_TIME("EddingtonSP1::setMultigridCoefficients");
  if(m_verbosity > 5){
    pout() << m_name + "::setMultigridCoefficients" << endl;
  }

  const int ncomp = 1;
  const int ghost = 3;

  m_amr->allocate(m_aCoef,        m_realm, m_phase, ncomp, ghost);
  m_amr->allocate(m_bco,        m_realm, m_phase, ncomp, ghost);
  m_amr->allocate(m_bco_irreg,  m_realm, m_phase, ncomp, ghost);

  this->setACoefAndBCoef();
}

void EddingtonSP1::setACoefAndBCoef(){
  CH_TIME("EddingtonSP1::setACoefAndBCoef");
  if(m_verbosity > 5){
    pout() << m_name + "::setACoefAndBCoef" << endl;
  }

  // This loop fills aco with kappa and bco_irreg with 1./kappa
  if(m_RtSpecies->isKappaConstant()){
    const Real kap = m_RtSpecies->getKappa(RealVect::Zero);
    DataOps::setValue(m_aCoef, kap);
    DataOps::setValue(m_bco, 1./kap);
    DataOps::setValue(m_bco_irreg, 1./kap);
  }
  else{ // If kappa is not constant, we need to go through each cell to determine it
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const RealVect origin = m_amr->getProbLo();
      const Real dx         = m_amr->getDx()[lvl];
    
      LevelData<EBCellFAB>& aco            = *m_aCoef[lvl];
      LevelData<EBFluxFAB>& bco            = *m_bco[lvl];
      LevelData<BaseIVFAB<Real> >& bco_irr = *m_bco_irreg[lvl];

      for (DataIterator dit = aco.dataIterator(); dit.ok(); ++dit){
	const Box box = (m_amr->getGrids(m_realm)[lvl]).get(dit());
	this->setACoefAndBCoef_box(aco[dit()], bco_irr[dit()], box, origin, dx, lvl, dit());
      }
    }

    m_amr->averageDown(m_aCoef, m_realm, m_phase);
    m_amr->interpGhost(m_aCoef, m_realm, m_phase);
    DataOps::averageCellToFaceAllComps(m_bco, m_aCoef, m_amr->getDomains()); // Average aco onto face
    DataOps::invert(m_bco); // Make m_bco = 1./kappa
  }

#if EddingtonSP1_feature // Different scaling for the RTE
  DataOps::scale(m_aCoef,       1.0);       // aco = c*kappa
  DataOps::scale(m_bco,       1.0/(3.0)); // bco = c/(3*kappa)
  DataOps::scale(m_bco_irreg, 1.0/(3.0)); // bco = c/(3*kappa)
#else // Original code before different scaling
  DataOps::scale(m_aCoef,       Units::c);       // aco = c*kappa
  DataOps::scale(m_bco,       Units::c/(3.0)); // bco = c/(3*kappa)
  DataOps::scale(m_bco_irreg, Units::c/(3.0)); // bco = c/(3*kappa)
#endif
}

void EddingtonSP1::setACoefAndBCoef_box(EBCellFAB&       a_aco,
					BaseIVFAB<Real>& a_bco,
					const Box        a_box,
					const RealVect   a_origin,
					const Real       a_dx,
					const int        a_lvl,
					const DataIndex& a_dit){
  CH_TIME("EddingtonSP1::setACoefAndBCoef_box");
  if(m_verbosity > 10){
    pout() << m_name + "::setACoefAndBCoef_box" << endl;
  }
  
  const int comp  = 0;
  const int ncomp = 1;
  
  const EBISBox& ebisbox = a_aco.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();

  // Regular aco
  BaseFab<Real>& aco_fab = a_aco.getSingleValuedFAB();
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv = bit();

    const RealVect pos = a_origin + iv*a_dx*RealVect::Unit;
    aco_fab(iv, comp) = m_RtSpecies->getKappa(pos);
  }


  // Irregular stuff
  VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();

    const RealVect pos  = EBArith::getVofLocation(vof, a_dx*RealVect::Unit, a_origin);
    const Real tmp = m_RtSpecies->getKappa(pos);
    a_aco(vof, comp) = tmp;
    a_bco(vof, comp) = 1./tmp;
  }
}

void EddingtonSP1::defineDeeperMultigridLevels(){
  CH_TIME("EddingtonSP1::defineDeeperMultigridLevels");
  if(m_verbosity > 5){
    pout() << m_name + "::defineDeeperMultigridLevels" << endl;
  }

  const int coar_ref = 2;

  Vector<ProblemDomain> m_mg_domains(0);
  Vector<DisjointBoxLayout> m_mg_grids(0);
  
  m_mg_domains.resize(0);
  m_mg_grids.resize(0);
  m_mg_levelgrids.resize(0);

  // Get some stuff from AmrMesh on how to decompose the levels
  const Vector<ProblemDomain>& domains = m_amr->getDomains();
  const int max_box_size               = m_amr->getMaxBoxSize();
  const int blocking_factor            = m_amr->getBlockingFactor();
  const int num_numEbGhostsCells       = m_amr->getNumberOfEbGhostCells();

  int num_coar       = 0;
  bool has_coar      = true;
  ProblemDomain fine = domains[0];

  // Coarsen problem domains and create grids
  while(num_coar < m_numCoarseningsBeforeAggregation || !has_coar){

    // Check if we can coarsen
    const ProblemDomain coar = fine.coarsen(coar_ref);
    const Box coar_box       = coar.domainBox();
    for (int dir = 0; dir < SpaceDim; dir++){
      if(coar_box.size()[dir] < max_box_size || coar_box.size()[dir]%max_box_size != 0){
	has_coar = false;
      }
    }

    if(has_coar){
      // Split the domain into pieces, then order and load balance them
      Vector<Box> boxes;
      Vector<int> proc_assign;
      domainSplit(coar, boxes, max_box_size, blocking_factor);
      mortonOrdering(boxes);
      LoadBalancing::makeBalance(proc_assign, boxes);

      // Add problem domain and grid
      m_mg_domains.push_back(coar);
      m_mg_grids.push_back(DisjointBoxLayout(boxes, proc_assign, coar));

      // Define the EBLevelGrids
      const int idx = m_mg_grids.size() - 1; // Last element added
      m_mg_levelgrids.push_back(EBLevelGrid(m_mg_grids[idx],
					    m_mg_domains[idx],
					    num_numEbGhostsCells,
					    m_ebis));

      // Next iterate
      fine = coar;
      num_coar++;
    }
    else{
      break;
    }
  }
}

void EddingtonSP1::setupOperatorFactory(){
  CH_TIME("EddingtonSP1::setupOperatorFactory");
  if(m_verbosity > 5){
    pout() << m_name + "::setupOperatorFactory" << endl;
  }

  const int finest_level                 = m_amr->getFinestLevel();
  const int ghost                        = m_amr->getNumberOfGhostCells();
  const Vector<DisjointBoxLayout>& grids = m_amr->getGrids(m_realm);
  const Vector<int>& refinement_ratios   = m_amr->getRefinementRatios();
  const Vector<ProblemDomain>& domains   = m_amr->getDomains();
  const Vector<Real>& dx                 = m_amr->getDx();
  const RealVect& origin                 = m_amr->getProbLo();
  const Vector<EBISLayout>& ebisl        = m_amr->getEBISLayout(m_realm, m_phase);
  
  const Vector<RefCountedPtr<EBQuadCFInterp> >& quadcfi  = m_amr->getEBQuadCFInterp(m_realm, m_phase);
  const Vector<RefCountedPtr<EBFluxRegister> >& fastFR   = m_amr->getFluxRegister(m_realm, m_phase);

  Vector<EBLevelGrid> levelgrids;

  for (int lvl = 0; lvl <= finest_level; lvl++){ 
    levelgrids.push_back(*(m_amr->getEBLevelGrid(m_realm, m_phase)[lvl])); // AmrMesh uses RefCounted levelgrids. EBConductivityOp does not. 
  }

#if 0
  Vector<EBLevelGrid> mg_levelgrids;
  Vector<RefCountedPtr<EBLevelGrid> >& mg_eblg = m_amr->get_mg_eblg(m_realm, m_phase);
  for (int lvl = 0; lvl < mg_eblg.size(); lvl++){
    mg_levelgrids.push_back(*mg_eblg[lvl]);
  }
#endif

  // Appropriate coefficients. 
  const Real alpha =  1.0;
  const Real beta  = -1.0;

  // Appropriate coefficients for this type of Robin BC
  m_robinCoefficients = RefCountedPtr<LarsenCoefficients> (new LarsenCoefficients(m_RtSpecies, m_reflectionCoefficientOne, m_reflectionCoefficientTwo));

  // Domain BC
#if 0
  m_robinDomainBcFactory = RefCountedPtr<RobinConductivityDomainBcFactory> (new RobinConductivityDomainBcFactory());
  m_robinDomainBcFactory->setCoefficients(m_robinCoefficients);
#else
  RefCountedPtr<BaseDomainBCFactory> domfact = RefCountedPtr<BaseDomainBCFactory>(NULL);
  ConductivityDomainBcWrapperFactory* bcfact = new ConductivityDomainBcWrapperFactory();
  Vector<RefCountedPtr<RobinCoefficients> > coefs(2*SpaceDim, m_robinCoefficients);
  bcfact->setWallBc(m_wallBc);
  bcfact->setRobinCoefficients(coefs);
  domfact = RefCountedPtr<BaseDomainBCFactory> (bcfact);
#endif

  // EBBC
  m_robinEbBcFactory  = RefCountedPtr<RobinConductivityEbBcFactory> (new RobinConductivityEbBcFactory(origin));
  m_robinEbBcFactory->setCoefficients(m_robinCoefficients);
  m_robinEbBcFactory->setStencilType(IrregStencil::StencilType::LeastSquares);

  //  Make relaxation type into int code
  int relax_type = 0;
  if(m_multigridRelaxMethod == RelaxationMethod::Jacobi){
    relax_type = 0;
  }
  else if(m_multigridRelaxMethod == RelaxationMethod::GaussSeidel){
    relax_type = 1;
  }
  else if(m_multigridRelaxMethod == RelaxationMethod::GSRBFast){
    relax_type = 2;
  }

  // Create operator factory.
  m_operatorFactory = RefCountedPtr<EbHelmholtzOpFactory> (new EbHelmholtzOpFactory(levelgrids,
										    quadcfi,
										    fastFR,
										    alpha,
										    beta,
										    m_aCoef.getData(),
										    m_bco.getData(),
										    m_bco_irreg.getData(),
										    dx[0],
										    refinement_ratios,
										    domfact,
										    m_robinEbBcFactory,
										    ghost*IntVect::Unit,
										    ghost*IntVect::Unit,
										    relax_type,
										    m_numCellsBottomDrop,
										    -1,
										    m_mg_levelgrids));
}

void EddingtonSP1::setupMultigridSolver(){
  CH_TIME("EddingtonSP1::setupMultigridSolver");
  if(m_verbosity > 5){
    pout() << m_name + "::setupMultigridSolver" << endl;
  }

  const int finest_level       = m_amr->getFinestLevel();
  const ProblemDomain coar_dom = m_amr->getDomains()[0];

  // Select bottom solver
  LinearSolver<LevelData<EBCellFAB> >* botsolver = NULL;
  if(m_bottomSolver == BottomSolver::Simple){
    m_simpleSolver.setNumSmooths(m_numSmoothingsForSimpleSolver);
    botsolver = &m_simpleSolver;
  }
  else if(m_bottomSolver == BottomSolver::BiCGStab){
    botsolver = &m_bicgstab;
  }

  // Make m_multigridType into an int for multigrid
  int gmg_type;
  if(m_multigridType == MultigridType::FAS){
    gmg_type = 0;
  }
  else if(m_multigridType == MultigridType::VCycle){
    gmg_type = 1;
  }
  else if(m_multigridType == MultigridType::FCycle){
    gmg_type = 2;
  }

  m_multigridSolver = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > (new AMRMultiGrid<LevelData<EBCellFAB> >());
  m_multigridSolver->define(coar_dom, *m_operatorFactory, botsolver, 1 + finest_level);
  m_multigridSolver->setSolverParameters(m_multigridPreSmooth,
					 m_multigridPostSmooth,
					 m_multigridBottomSmooth,
					 gmg_type,
					 m_multigridMaxIterations,
					 m_multigridTolerance,
					 m_multigridHang,
					 1.E-99); // Residue set through other means
  m_multigridSolver->m_imin    = m_multigridMinIterations;
  m_multigridSolver->m_verbosity = m_multigridVerbosity;

  // Dummies for init
  const int ncomp = 1;
  EBAMRCellData dummy1, dummy2;
  m_amr->allocate(dummy1, m_realm, m_phase, ncomp);
  m_amr->allocate(dummy2, m_realm, m_phase, ncomp);
  DataOps::setValue(dummy1, 0.0);
  DataOps::setValue(dummy2, 0.0);

  // Aliasing
  Vector<LevelData<EBCellFAB>* > phi, rhs;
  m_amr->alias(phi, dummy1);
  m_amr->alias(rhs, dummy2);

  // Init solver
  m_multigridSolver->init(phi, rhs, finest_level, 0);
}

void EddingtonSP1::setupTGA(){
  CH_TIME("EddingtonSP1::setupTGA");
  if(m_verbosity > 5){
    pout() << m_name + "::setupTGA" << endl;
  }
  
  const int finest_level       = m_amr->getFinestLevel();
  const ProblemDomain coar_dom = m_amr->getDomains()[0];
  const Vector<int> ref_rat    = m_amr->getRefinementRatios();

  m_tgaSolver = RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > >
    (new AMRTGA<LevelData<EBCellFAB> > (m_multigridSolver, *m_operatorFactory, coar_dom, ref_rat, 1 + finest_level, m_multigridSolver->m_verbosity));

  // Must init gmg for TGA
  Vector<LevelData<EBCellFAB>* > phi, rhs;
  m_amr->alias(phi, m_phi);
  m_amr->alias(rhs, m_source);
  m_multigridSolver->init(phi, rhs, finest_level, 0);
}

void EddingtonSP1::setupEuler(){
  CH_TIME("EddingtonSP1::setupEuler");
  if(m_verbosity > 5){
    pout() << m_name + "::setupEuler" << endl;
  }

  const int finest_level       = m_amr->getFinestLevel();
  const ProblemDomain coar_dom = m_amr->getDomains()[0];
  const Vector<int> ref_rat    = m_amr->getRefinementRatios();

  m_eulerSolver = RefCountedPtr<EBBackwardEuler> 
    (new EBBackwardEuler (m_multigridSolver, *m_operatorFactory, coar_dom, ref_rat, 1 + finest_level, m_multigridSolver->m_verbosity));

  // Note: If this crashes, try to init gmg first
}

void EddingtonSP1::computeBoundaryFlux(EBAMRIVData& a_ebFlux, const EBAMRCellData& a_phi){
  CH_TIME("EddingtonSP1::computeBoundaryFlux");
  if(m_verbosity > 5){
    pout() << m_name + "::computeBoundaryFlux" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();
  
  const IrregAmrStencil<EbCentroidInterpolationStencil>& sten = m_amr->getEbCentroidInterpolationStencilStencils(m_realm, m_phase);
  for(int lvl = 0; lvl <= finest_level; lvl++){
    sten.apply(*a_ebFlux[lvl], *a_phi[lvl], lvl, true);
  }

  m_amr->averageDown(a_ebFlux, m_realm, m_phase);

  DataOps::scale(a_ebFlux, 0.5*Units::c);
}

void EddingtonSP1::computeDomainFlux(EBAMRIFData& a_domainflux, const EBAMRCellData& a_data){
  CH_TIME("EddingtonSP1::computeDomainFlux");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDomainFlux" << endl;
  }


  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const int ncomp = a_data[lvl]->nComp();
      
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const EBCellFAB& data         = (*a_data[lvl])[dit()];
      const EBISBox& ebisbox        = ebisl[dit()];
      const BaseFab<Real>& data_fab = data.getSingleValuedFAB();
      
      for (int dir = 0; dir < SpaceDim; dir++){
	for (SideIterator sit; sit.ok(); ++sit){
	  BaseIFFAB<Real>& extrap = (*a_domainflux[lvl])[dit()](dir, sit());

	  const IntVectSet& ivs  = extrap.getIVS();
	  const EBGraph& ebgraph = extrap.getEBGraph();

	  // Extrapolate to the boundary. Use face-centered stuff for all faces (also multivalued ones)
	  const FaceStop::WhichFaces crit = FaceStop::AllBoundaryOnly;
	  for (FaceIterator faceit(ivs, ebgraph, dir, crit); faceit.ok(); ++faceit){
	    const FaceIndex& face = faceit();

	    const int sgn = sign(sit()); // Lo = -1, Hi = 1
	    
	    const VolIndex& vof = face.getVoF(flip(sit()));
	    const IntVect iv0   = vof.gridIndex();
	    const IntVect iv1   = iv0 - sgn*BASISV(dir);

	    if(ebisbox.isCovered(iv0)){ // Just provide some bogus data because the face
	      for (int comp = 0; comp < ncomp; comp++){
		extrap(face, comp) = 0.0;
	      }
	    }
	    else{
	      if(!ebisbox.isCovered(iv1)){ // linear extrapolation
		for (int comp = 0; comp < ncomp; comp++){
		  extrap(face, comp) = 1.5*data_fab(iv0, comp) - 0.5*data_fab(iv1, comp); // Should be ok
		}
	      }
	      else{ // Not enough cells available, use cell-centered only
		for (int comp = 0; comp < ncomp; comp++){
		  extrap(face, comp) = data_fab(iv0, comp);
		}
	      }
	    }

	    // Necessary scaling
	    for (int comp = 0; comp < ncomp; comp++){
	      extrap(face, comp) = 0.5*Units::c*extrap(face, comp);
	    }
	  }
	}
      }
    }
  }
}

void EddingtonSP1::computeFlux(EBAMRCellData& a_flux, const EBAMRCellData& a_phi){
  CH_TIME("EddingtonSP1::computeFlux");
  if(m_verbosity > 5){
    pout() << m_name + "::computeFlux" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();

  m_amr->computeGradient(a_flux, a_phi, m_realm, m_phase); // flux = grad(phi)
  for (int lvl = 0; lvl <= finest_level; lvl++){
    DataOps::divideByScalar(*a_flux[lvl], *m_aCoef[lvl]);   // flux = grad(phi)/(c*kappa)
    DataOps::scale(*a_flux[lvl], -Units::c*Units::c/3.0);  // flux = -c*grad(phi)/3.
  }

  m_amr->averageDown(a_flux, m_realm, m_phase);
  m_amr->interpGhost(a_flux, m_realm, m_phase);
}


void EddingtonSP1::computeDensity(EBAMRCellData& a_isotropic, const EBAMRCellData& a_phi){
  CH_TIME("EddingtonSP1::computeDensity");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDensity" << endl;
  }

  const int finest_level = m_amr->getFinestLevel();
  
  for (int lvl = 0; lvl <= finest_level; lvl++){
    const Interval interv(0,0);
    a_phi[lvl]->localCopyTo(interv, *a_isotropic[lvl], interv);
  }
}

void EddingtonSP1::writePlotFile(){
  CH_TIME("EddingtonSP1::writePlotFile");
  if(m_verbosity > 5){
    pout() << m_name + "::writePlotFile" << endl;
  }

  char file_char[1000];
  sprintf(file_char, "%s.step%07d.%dd.hdf5", m_name.c_str(), m_timeStep, SpaceDim);

  const int ncomps = 3 + SpaceDim;
  Vector<string> names(ncomps);
  names[0] = "density";
  names[1] = "x-flux";
  names[2] = "y-flux";
  if(SpaceDim == 3){
    names[3] = "z-flux";
    names[4] = "isotropic source";
    names[5] = "residue";
  }
  else{
    names[3] = "isotropic source";
    names[4] = "residue";
  }

  // Compute the flux
  EBAMRCellData flux;
  m_amr->allocate(flux, m_realm, m_phase, SpaceDim);
  this->computeFlux(flux, m_phi);

  // Allocate output storage
  EBAMRCellData output;
  m_amr->allocate(output, m_realm, m_phase, ncomps, 1);


  for (int lvl = 0; lvl < output.size(); lvl++){
    LevelData<EBCellFAB>& state  = *m_phi[lvl];
    LevelData<EBCellFAB>& source = *m_source[lvl];
    LevelData<EBCellFAB>& flx    = *flux[lvl];
    LevelData<EBCellFAB>& res    = *m_resid[lvl];

    state.localCopyTo(Interval(0,0),          *output[lvl], Interval(0,0));
    flx.localCopyTo(Interval(0,SpaceDim - 1), *output[lvl], Interval(1, SpaceDim));
    source.localCopyTo(Interval(0,0),         *output[lvl], Interval(1 + SpaceDim, 1 + SpaceDim));
    res.localCopyTo(Interval(0,0),            *output[lvl], Interval(2+SpaceDim, 2+SpaceDim));
  }

  // Transform to centroid-centered
  const IrregAmrStencil<CentroidInterpolationStencil>& sten = m_amr->getCentroidInterpolationStencils(m_realm, phase::gas);
  sten.apply(output);

  // Alias this stuff
  Vector<LevelData<EBCellFAB>* > output_ptr;
  m_amr->alias(output_ptr, output);

  Vector<Real> covered_values(ncomps, 0.0);
  string fname(file_char);
  writeEBHDF5(fname,
	      m_amr->getGrids(m_realm),
	      output_ptr,
	      names,
	      m_amr->getDomains()[0].domainBox(),
	      m_amr->getDx()[0],
	      m_dt,
	      m_time,
	      m_amr->getRefinementRatios(),
	      m_amr->getFinestLevel() + 1,
	      true,
	      covered_values);
}

void EddingtonSP1::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const {
  CH_TIME("EddingtonSP1::writeCheckpointLevel");
  if(m_verbosity > 5){
    pout() << m_name + "::writeCheckpointLevel" << endl;
  }

  // Write state vector
  write(a_handle, *m_phi[a_level], m_name);
}

void EddingtonSP1::readCheckpointLevel(HDF5Handle& a_handle, const int a_level){
  CH_TIME("EddingtonSP1::readCheckpointLevel");
  if(m_verbosity > 5){
    pout() << m_name + "::readCheckpointLevel" << endl;
  }

  read<EBCellFAB>(a_handle, *m_phi[a_level], m_name, m_amr->getGrids(m_realm)[a_level], Interval(0,0), false);
}

void EddingtonSP1::setNeumannWallBc(const int a_dir, Side::LoHiSide a_side, const Real a_value){
  CH_TIME("EddingtonSP1::setNeumannWallBc");
  if(m_verbosity > 5){
    pout() << "EddingtonSP1::setNeumannWallBc" << endl;
  }

  const int idx = WallBc::map_bc(a_dir, a_side);
  m_wallBc[idx] = RefCountedPtr<WallBc> (new WallBc(a_dir, a_side, wallbc::neumann));
  m_wallBc[idx]->setValue(a_value);
}

void EddingtonSP1::setRobinWallBc(const int a_dir, Side::LoHiSide a_side, const Real a_value){
  CH_TIME("EddingtonSP1::setRobinWallBc");
  if(m_verbosity > 5){
    pout() << "EddingtonSP1::setRobinWallBc" << endl;
  }

  const int idx = WallBc::map_bc(a_dir, a_side);
  m_wallBc[idx] = RefCountedPtr<WallBc> (new WallBc(a_dir, a_side, wallbc::robin));
  m_wallBc[idx]->setValue(a_value);
}

#include <CD_NamespaceFooter.H>
