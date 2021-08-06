/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
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
#if 1 // All of these should be removed
#include <EBAMRIO.H>
#include <BRMeshRefine.H>
#include <NeumannConductivityEBBC.H>
#include <DirichletConductivityEBBC.H>
#endif

// Our includes
#include <CD_EddingtonSP1.H>
#include <CD_DataOps.H>
#include <CD_Units.H>
#if 1 // Should be removed
#include <CD_ConductivityEddingtonSP1DomainBcFactory.H>
#include <CD_RobinConductivityEbBcFactory.H>
#endif
#include <CD_EBHelmholtzDirichletEBBCFactory.H>      
#include <CD_EBHelmholtzNeumannEBBCFactory.H>
#include <CD_EBHelmholtzLarsenEBBCFactory.H>
#include <CD_EBHelmholtzLarsenDomainBCFactory.H>
#include <CD_EBHelmholtzDirichletDomainBCFactory.H> // Domain stuff should be collapsed since we need a wrapper anyways. 
#include <CD_EBHelmholtzNeumannDomainBCFactory.H>
#include <CD_EBHelmholtzRobinDomainBCFactory.H>
#include <CD_NamespaceHeader.H>

constexpr Real EddingtonSP1::m_alpha;
constexpr Real EddingtonSP1::m_beta;

Real EddingtonSP1::s_defaultDomainBcFunction(const RealVect a_position, const Real a_time){
  return 1.0;
}

EddingtonSP1::EddingtonSP1() : RtSolver() {

  // Default settings
  m_name         = "EddingtonSP1";
  m_className    = "EddingtonSP1";
  
  m_verbosity                = -1;
  m_isSolverSetup            = false;
  m_hasDeeperMultigridLevels = false;
  m_dataLocation             = Location::Cell::Center;  
  
  this->setDefaultDomainBcFunctions(); // This fills m_domainBcFunctions with s_defaultDomainBcFunction on every domain side. 
}

EddingtonSP1::~EddingtonSP1(){
}

void EddingtonSP1::parseOptions(){
  CH_TIME("EddingtonSP1::parseOptions");
  if(m_verbosity > 5){
    pout() << m_name + "::parseOptions" << endl;
  }

  this->parseDomainBc();          // Parses domain bc
  this->parseEbBc();              // Parse ebbc
  this->parseReflection();        // Parses "reflection coefficients"
  this->parseStationary();        // Parse stationary solver
  this->parsePlotVariables();     // Parses plot variables
  this->parseMultigridSettings(); // Parses solver parameters for geometric multigrid
}

void EddingtonSP1::parseRuntimeOptions(){
  CH_TIME("EddingtonSP1::parseRuntimeOptions");
  if(m_verbosity > 5){
    pout() << m_name + "::parseRuntimeOptions" << endl;
  }
  
  this->parsePlotVariables();     // Parses plot variables
  this->parseMultigridSettings(); // Parses solver parameters for geometric multigrid
}

void EddingtonSP1::setDefaultDomainBcFunctions(){
  CH_TIME("EddingtonSP1::setDefaultDomainBcFunctions");
  if(m_verbosity > 5){
    pout() << m_name + "::setDefaultDomainBcFunctions" << endl;
  }

  m_domainBcFunctions.clear();
  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      m_domainBcFunctions.emplace(std::make_pair(dir, sit()), EddingtonSP1::s_defaultDomainBcFunction);
    }
  }
}

void EddingtonSP1::setDomainBcWallFunction(const int a_dir, const Side::LoHiSide a_side, const EddingtonSP1DomainBc::BcFunction& a_function){
  CH_TIME("EddingtonSP1::setDomainBcWallFunction");
  if(m_verbosity > 4){
    pout() << m_name + "::setDomainBcWallFunction" << endl;
  }

  const EddingtonSP1DomainBc::Wall curWall = std::make_pair(a_dir, a_side);

  m_domainBcFunctions.at(curWall) = a_function;
}

std::string EddingtonSP1::makeBcString(const int a_dir, const Side::LoHiSide a_side) const {
  CH_TIME("EddingtonSP1::makeBcString");
  if(m_verbosity > 5){
    pout() << m_name + "::makeBcString" << endl;
  }

  // TLDR: Returns string in format "bc.x.lo"
  std::string strDir;
  std::string strSide;
  
  if(a_dir == 0){
    strDir = "x";
  }
  else if(a_dir == 1){
    strDir = "y";
  }
  else if(a_dir == 2){
    strDir = "z";
  }

  if(a_side == Side::Lo){
    strSide = "lo";
  }
  else if(a_side == Side::Hi){
    strSide = "hi";
  }

  const std::string ret = std::string("bc.") + strDir + std::string(".") + strSide;

  return ret;
}

EddingtonSP1DomainBc::BcType EddingtonSP1::parseBcString(const std::string a_str) const {
  CH_TIME("EddingtonSP1::parseBcString");
  if(m_verbosity > 5){
    pout() << m_name + "::parseBcString" << endl;
  }

  EddingtonSP1DomainBc::BcType ret;

  if(a_str == "dirichlet"){
    ret = EddingtonSP1DomainBc::BcType::Dirichlet;
  }
  else if(a_str == "neumann"){
    ret = EddingtonSP1DomainBc::BcType::Neumann;
  }
  else if(a_str == "larsen"){
    ret = EddingtonSP1DomainBc::BcType::Larsen;
  }
  else{
    MayDay::Abort("EddingtonSP1DomainBc::BcType - unknown BC type!");
  }

  return ret;
}

void EddingtonSP1::parseDomainBc(){
  CH_TIME("EddingtonSP1::parseDomainBc");
  if(m_verbosity > 5){
    pout() << m_name + "::parseDomainBc" << endl;
  }

  // TLDR: This routine might seem big and complicated. What we are doing is that we are creating one function object which returns some value
  //       anywhere in space and time on a domain edge (face). The code below simply creates those functions and associates them with an edge (face in 3D)
  //
  //       Because of the stationary interface of EBHelmholtzOp and the transient interface of the radiative transfer solvers, we capture the time by solver
  //       time (RtSolver::m_time) by reference in these functions. That allows us to pass the functions in as stationary-looking functions in EBHelmholtzOp,
  //       yet retain transient BCs.
  //
  //       For flexibility we want to be able to specify the BC functions in two forms, using a multiplier or not. The specification for this is e.g.
  //       "dirichlet <number>" in which case the bc function is multiplied by <number>. Use "dirichlet_custom" to use the bc functions directly.

  ParmParse pp(m_className.c_str());

  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const EddingtonSP1DomainBc::Wall curWall = std::make_pair(dir, sit());
      const std::string bcString = this->makeBcString(dir, sit());
      const int num = pp.countval(bcString.c_str());

      std::string str;

      EddingtonSP1DomainBc::BcType      bcType;
      EddingtonSP1DomainBc::BcFunction& bcFunc = m_domainBcFunctions.at(curWall);

      std::function<Real(const RealVect, const Real)> curFunc;

      switch(num){
      case 1:{ // If only providing dirichlet_custom, neumann_custom, larsen_custom, the second value is overridden and only the function is used as BC. 
	pp.get(bcString.c_str(), str, 0);

	// Set the function. Capture solver time by reference. 
	curFunc = [&bcFunc, &time=this->m_time] (const RealVect a_pos, const Real a_time) {
	  return bcFunc(a_pos, time);
	};

	// Ensure that we actually have a valid BC. 
	if(str == "dirichlet_custom"){
	  bcType  = EddingtonSP1DomainBc::BcType::Dirichlet;
	}
	else if(str == "neumann_custom"){
	  bcType  = EddingtonSP1DomainBc::BcType::Neumann;
	}
	else if(str == "larsen_custom"){
	  bcType  = EddingtonSP1DomainBc::BcType::Larsen;
	}
	else{
	  MayDay::Abort("EddingtonSP1::parseDomainBc -- got only one argument but this argument was not in the form of e.g. 'neumann_custom'!");
	}

	break;
      }
      case 2:{ // Had two arguments in input, e.g. neumann 0.0. In this case we set the BC to be the specified function times the value.
	Real val;
	
	pp.get(bcString.c_str(), str, 0);
	pp.get(bcString.c_str(), val, 1);

	bcType = this->parseBcString(str);

	curFunc = [&bcFunc, &time=this->m_time, val] (const RealVect a_pos, const Real a_time) {
	  return bcFunc(a_pos, time)*val;
	};

	break;
      }
      default: {
	MayDay::Abort("EddingtonSP1::parseDomainBc -- unknown argument to domain bc");
	break;
      }
      }

      m_domainBc.setBc(curWall, std::make_pair(bcType, curFunc));
    }
  }
}

void EddingtonSP1::parseEbBc(){
  CH_TIME("EddingtonSP1::parseEbBc");
  if(m_verbosity > 5){
    pout() << m_name + "::parseEbBc" << endl;
  }

  ParmParse pp(m_className.c_str());

  const std::string bcString = "ebbc";
  const int num = pp.countval(bcString.c_str());

  if(num == 2){
    Real val;
    std::string str;

    pp.get(bcString.c_str(), str, 0);
    pp.get(bcString.c_str(), val, 1);

    EbBcType bcType;

    if(str == "dirichlet"){
      bcType = EbBcType::Dirichlet;
    }
    else if(str == "neumann"){
      bcType = EbBcType::Neumann;
    }
    else if(str == "larsen"){
      bcType = EbBcType::Larsen;
    }
    else{
      MayDay::Abort("EddingtonSP1::parseEbBc -- not dirichlet/neumann/larsen!");
    }

    m_ebbc = std::make_pair(bcType, val);
  }
  else{
    MayDay::Abort("EddingtonSP1::parseEbBc -- bad input argument. Should be in the form 'larsen 0.0', 'dirichlet 0.0', etc");
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

  pp.get("gmg_verbosity",    m_multigridVerbosity);
  pp.get("gmg_coarsen",      m_numCoarseningsBeforeAggregation);
  pp.get("gmg_pre_smooth",   m_multigridPreSmooth);
  pp.get("gmg_post_smooth",  m_multigridPostSmooth);
  pp.get("gmg_bott_smooth",  m_multigridBottomSmooth);
  pp.get("gmg_max_iter",     m_multigridMaxIterations);
  pp.get("gmg_min_iter",     m_multigridMinIterations);
  pp.get("gmg_exit_tol",     m_multigridExitTolerance);
  pp.get("gmg_exit_hang",    m_multigridExitHang);
  pp.get("gmg_min_cells",    m_minCellsBottom);
  pp.get("gmg_bc_order",     m_multigridBcOrder);
  pp.get("gmg_bc_weight",    m_multigridBcWeight);

  // Fetch the desired bottom solver from the input script. We look for things like EddingtonSP1.gmg_bottom_solver = bicgstab or '= simple <number>'
  // where <number> is the number of relaxation for the smoothing solver. 
  const int num = pp.countval("gmg_bottom_solver");
  if(num == 1){
    pp.get("gmg_bottom_solver", str);
    if(str == "bicgstab"){
      m_bottomSolverType = BottomSolverType::BiCGStab;
    }
    else{
      MayDay::Error("FieldSolverMultigrid::parseMultigridSettings - logic bust, you've specified one parameter and I expected either 'bicgstab' or 'gmres'");
    }
  }
  else if(num == 2){
    int numSmooth;
    pp.get("gmg_bottom_solver", str,       0);
    pp.get("gmg_bottom_solver", numSmooth, 1);
    if(str == "simple"){
      m_bottomSolverType = BottomSolverType::Simple;
      m_simpleSolver.setNumSmooths(numSmooth);
    }
    else{
      MayDay::Error("FieldSolverMultigrid::parseMultigridSettings - logic bust, you've specified two parameters and I expected 'simple <number>'");
    }
  }
  else{
    MayDay::Error("FieldSolverMultigrid::parseMultigridSettings - logic bust in bottom solver. You must specify ' = bicgstab', ' = gmres', or ' = simple <number>'");
  }

  // Relaxation type
  pp.get("gmg_smoother", str);
  if(str == "jacobi"){
    m_multigridRelaxMethod = EBHelmholtzOp::Smoother::PointJacobi;    
  }
  else if(str == "red_black"){
    m_multigridRelaxMethod = EBHelmholtzOp::Smoother::GauSaiRedBlack;    
  }
  else if(str == "multi_color"){
    m_multigridRelaxMethod = EBHelmholtzOp::Smoother::GauSaiMultiColor;
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
  if(m_minCellsBottom < 2){
    m_minCellsBottom = 2;
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

  m_reflectCoefOne = r/(2.0);
  m_reflectCoefTwo = r/(3.0);
}

void EddingtonSP1::preRegrid(const int a_base, const int a_oldFinestLevel){
  CH_TIME("EddingtonSP1::preRegrid");
  if(m_verbosity > 5){
    pout() << m_name + "::preRegrid" << endl;
  }

  const int finestLevel = m_amr->getFinestLevel();

  m_amr->allocate(m_cache, m_realm, m_phase, m_nComp);

  for (int lvl = 0; lvl <= finestLevel; lvl++){
    m_phi[lvl]->localCopyTo(*m_cache[lvl]);
  }
}

void EddingtonSP1::allocateInternals(){
  CH_TIME("EddingtonSP1::allocateInternals");
  if(m_verbosity > 5){
    pout() << m_name + "::allocateInternals" << endl;
  }
  
  m_amr->allocate(m_aco,       m_realm, m_phase, m_nComp);
  m_amr->allocate(m_bco,       m_realm, m_phase, m_nComp);
  m_amr->allocate(m_bco_irreg, m_realm, m_phase, m_nComp);
  m_amr->allocate(m_phi,       m_realm, m_phase, m_nComp);
  m_amr->allocate(m_source,    m_realm, m_phase, m_nComp);
  m_amr->allocate(m_resid,     m_realm, m_phase, m_nComp);

  DataOps::setValue(m_resid,  0.0);
  DataOps::setValue(m_phi,    0.0);
  DataOps::setValue(m_source, 0.0);

  this->setACoefAndBCoef();
}

void EddingtonSP1::deallocateInternals(){
  m_amr->deallocate(m_aco);
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

  const Interval interv(m_comp, m_comp);

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

  m_isSolverSetup = false;
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
    m_amr->registerOperator(s_eb_pwl_interp,   m_realm, m_phase);
    m_amr->registerOperator(s_eb_multigrid,    m_realm, m_phase);
  }
}

bool EddingtonSP1::advance(const Real a_dt, EBAMRCellData& a_phi, const EBAMRCellData& a_source, const bool a_zerophi){
  CH_TIME("EddingtonSP1::advance(ebamrcell, ebamrcell)");
  if(m_verbosity > 5){
    pout() << m_name + "::advance(ebamrcell, ebamrcell)" << endl;
  }

  const int finestLevel = m_amr->getFinestLevel();

  if(!m_isSolverSetup){
    this->setupSolver();
  }

  bool converged;

  // Must have a dummy for checking initial residual
  EBAMRCellData dummy;
  EBAMRCellData source;
  m_amr->allocate(dummy,  m_realm, m_phase, m_nComp);
  m_amr->allocate(source, m_realm, m_phase, m_nComp);
  DataOps::setValue(dummy, 0.0);

  // Various source term manipulations. 
  DataOps::setValue(source, 0.0);
  DataOps::incr(source, a_source, 1.0);
  DataOps::scale(source, 1./Units::c); // Source should be scaled by 1./c0

  if(m_stationary){ // Should kappa-scale for transient solvres
    DataOps::kappaScale(source);
  }

  Vector<LevelData<EBCellFAB>* > phi, rhs, res, zero;
  m_amr->alias(phi,  a_phi);
  m_amr->alias(rhs,  source);
  m_amr->alias(res,  m_resid);
  m_amr->alias(zero, dummy);

  if(m_stationary){
    const Real phi_resid  = m_multigridSolver->computeAMRResidual(phi,  rhs, finestLevel, 0); // Incoming residual
    const Real zero_resid = m_multigridSolver->computeAMRResidual(zero, rhs, finestLevel, 0); // Zero residual

    if(phi_resid > zero_resid*m_multigridExitTolerance){ // Residual is too large
      m_multigridSolver->m_convergenceMetric = zero_resid;
      m_multigridSolver->solveNoInitResid(phi, res, rhs, finestLevel, 0, a_zerophi);
      
      const int status = m_multigridSolver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
      if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
	converged = true;
      }
    }
    else{ // Solution is already good enough
      converged = true;
    }
    m_multigridSolver->revert(phi, rhs, finestLevel, 0);

    DataOps::setCoveredValue(a_phi, 0, 0.0);
  }
  else{
    if(m_useTGA){

      m_tgaSolver->oneStep(res, phi, rhs, Units::c*a_dt, 0, finestLevel, m_time);

    }
    else{

      m_eulerSolver->oneStep(res, phi, rhs, Units::c*a_dt, 0, finestLevel, false);

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

void EddingtonSP1::setupSolver(){
  CH_TIME("EddingtonSP1::setupSolver");
  if(m_verbosity > 5){
    pout() << m_name + "::setupSolver" << endl;
  }

  if(!m_hasDeeperMultigridLevels){
    this->defineDeeperMultigridLevels();
    m_hasDeeperMultigridLevels = true;
  }
  
  this->setMultigridCoefficients(); // Set coefficients, kappa, aco, bco
  this->setupOperatorFactory();     // Set the operator factory
  this->setupHelmholtzFactory();    // Set up the Helmholtz operator factory
  this->setupMultigrid();     // Set up the AMR multigrid solver

  if(!m_stationary){
    if(m_useTGA){
      this->setupTGA();
    }
    else{
      this->setupEuler();
    }
  }

  m_isSolverSetup = true;
}

void EddingtonSP1::setMultigridCoefficients(){
  CH_TIME("EddingtonSP1::setMultigridCoefficients");
  if(m_verbosity > 5){
    pout() << m_name + "::setMultigridCoefficients" << endl;
  }

  m_amr->allocate(m_aco,        m_realm, m_phase, m_nComp);
  m_amr->allocate(m_bco,        m_realm, m_phase, m_nComp);
  m_amr->allocate(m_bco_irreg,  m_realm, m_phase, m_nComp);

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
    DataOps::setValue(m_aco, kap);
    DataOps::setValue(m_bco, 1./kap);
    DataOps::setValue(m_bco_irreg, 1./kap);
  }
  else{ // If kappa is not constant, we need to go through each cell to determine it
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
      const RealVect origin = m_amr->getProbLo();
      const Real dx         = m_amr->getDx()[lvl];
    
      LevelData<EBCellFAB>& aco            = *m_aco[lvl];
      LevelData<EBFluxFAB>& bco            = *m_bco[lvl];
      LevelData<BaseIVFAB<Real> >& bco_irr = *m_bco_irreg[lvl];

      for (DataIterator dit = aco.dataIterator(); dit.ok(); ++dit){
	const Box box = (m_amr->getGrids(m_realm)[lvl]).get(dit());
	this->setACoefAndBCoefBox(aco[dit()], bco_irr[dit()], box, origin, dx, lvl, dit());
      }
    }

    m_amr->averageDown(m_aco, m_realm, m_phase);
    m_amr->interpGhost(m_aco, m_realm, m_phase);
    DataOps::averageCellToFaceAllComps(m_bco, m_aco, m_amr->getDomains()); // Average aco onto face
    DataOps::invert(m_bco); // Make m_bco = 1./kappa
  }


  DataOps::scale(m_aco,       1.0);       // aco = c*kappa
  DataOps::scale(m_bco,       1.0/(3.0)); // bco = c/(3*kappa)
  DataOps::scale(m_bco_irreg, 1.0/(3.0)); // bco = c/(3*kappa)

}

void EddingtonSP1::setACoefAndBCoefBox(EBCellFAB&       a_aco,
				       BaseIVFAB<Real>& a_bco,
				       const Box        a_box,
				       const RealVect   a_origin,
				       const Real       a_dx,
				       const int        a_lvl,
				       const DataIndex& a_dit){
  CH_TIME("EddingtonSP1::setACoefAndBCoefBox");
  if(m_verbosity > 10){
    pout() << m_name + "::setACoefAndBCoefBox" << endl;
  }
  
  const EBISBox& ebisbox = a_aco.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();

  // Regular aco
  BaseFab<Real>& aco_fab = a_aco.getSingleValuedFAB();
  for (BoxIterator bit(a_box); bit.ok(); ++bit){
    const IntVect iv = bit();

    const RealVect pos = a_origin + iv*a_dx*RealVect::Unit;
    aco_fab(iv, m_comp) = m_RtSpecies->getKappa(pos);
  }


  // Irregular stuff
  VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_lvl])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof = vofit();

    const RealVect pos  = EBArith::getVofLocation(vof, a_dx*RealVect::Unit, a_origin);
    const Real tmp = m_RtSpecies->getKappa(pos);
    a_aco(vof, m_comp) = tmp;
    a_bco(vof, m_comp) = 1./tmp;
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

void EddingtonSP1::setupHelmholtzFactory(){
  CH_TIME("EddingtonSP1::setupHelmholtzFactory");
  if(m_verbosity > 5){
    pout() << m_name + "::setupHelmholtzFactory" << endl;
  }

  const int finestLevel = m_amr->getFinestLevel();

  const Vector<RefCountedPtr<EBLevelGrid> >&             levelGrids   = m_amr->getEBLevelGrid          (m_realm, m_phase);
  const Vector<RefCountedPtr<EbCoarAve> >&               coarAve      = m_amr->getCoarseAverage        (m_realm, m_phase);
  const Vector<RefCountedPtr<EBFluxRegister> >&          fluxReg      = m_amr->getFluxRegister         (m_realm, m_phase);
  const Vector<RefCountedPtr<EBMultigridInterpolator> >& interpolator = m_amr->getMultigridInterpolator(m_realm, m_phase);

  // Coarsest domain used for multigrid. The user specifies the minimum number of cells in any
  // coordinate direction, and we coarsen until we have a domain which satisfies that constraint. 
  ProblemDomain bottomDomain = m_amr->getDomains()[0];
  while(bottomDomain.domainBox().shortside() >= 2*m_minCellsBottom){
    bottomDomain.coarsen(2);
  }

  // Number of ghost cells in data holders
  const IntVect ghostPhi   = m_amr->getNumberOfGhostCells()*IntVect::Unit;
  const IntVect ghostRhs   = m_amr->getNumberOfGhostCells()*IntVect::Unit;
  
  // Just run with Larsen for now. These need to be updated with proper source terms. 
  auto domainBcFactory = RefCountedPtr<EBHelmholtzDomainBCFactory> (new EBHelmholtzLarsenDomainBCFactory(m_RtSpecies, m_reflectCoefOne, m_reflectCoefTwo));
  auto ebBcFactory     = RefCountedPtr<EBHelmholtzEBBCFactory>     (new EBHelmholtzLarsenEBBCFactory    (m_RtSpecies, m_reflectCoefOne, m_reflectCoefTwo));

  // Set up the operator
  m_helmholtzOpFactory = RefCountedPtr<EBHelmholtzOpFactory> (new EBHelmholtzOpFactory(m_dataLocation,
										       m_alpha,
										       m_beta,
										       m_amr->getProbLo(),
										       levelGrids,
										       interpolator,
										       fluxReg,
										       coarAve,
										       m_amr->getRefinementRatios(),
										       m_amr->getDx(),
										       m_aco.getData(),
										       m_bco.getData(),
										       m_bco_irreg.getData(),
										       domainBcFactory,
										       ebBcFactory,
										       ghostPhi,
										       ghostRhs,
										       m_multigridRelaxMethod,
										       bottomDomain,
										       m_amr->getMaxBoxSize()));
										       

  
}

void EddingtonSP1::setupOperatorFactory(){
  CH_TIME("EddingtonSP1::setupOperatorFactory");
  if(m_verbosity > 5){
    pout() << m_name + "::setupOperatorFactory" << endl;
  }

  // Set the domain bc factory. This can use arbitrary neumann/dirichelt/larsen with functions. This the same larsen coefficients on all domain sides.
  // If you want to use different coefficients on different domain sides, here is where you would do it. 
  std::map<EddingtonSP1DomainBc::Wall, RefCountedPtr<RobinCoefficients> > larsenCoeffs;
  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      EddingtonSP1DomainBc::Wall             curWall = std::make_pair(dir, sit());
      const EddingtonSP1DomainBc::BcFunction curFunc = m_domainBc.getBc(curWall).second;

      auto larsen = RefCountedPtr<LarsenCoefficients>(new LarsenCoefficients(m_RtSpecies, m_reflectCoefOne, m_reflectCoefTwo, curFunc));

      larsenCoeffs.emplace(curWall, larsen);
    }
  }
  auto fact = RefCountedPtr<ConductivityEddingtonSP1DomainBcFactory> (new ConductivityEddingtonSP1DomainBcFactory(m_domainBc, larsenCoeffs, m_amr->getProbLo()));

  // Set the embedded boundary bc factory. This uses a constant value for the right-hand side on the EB. 
  RefCountedPtr<BaseEBBCFactory> ebbcFactory;
  switch(m_ebbc.first){
  case EbBcType::Dirichlet: {
    auto dirichletEbBcFactory = RefCountedPtr<DirichletConductivityEBBCFactory> (new DirichletConductivityEBBCFactory());
    dirichletEbBcFactory->setValue(m_ebbc.second);
    dirichletEbBcFactory->setOrder(1);
    ebbcFactory = RefCountedPtr<BaseEBBCFactory> (dirichletEbBcFactory);
    break;
  }
  case EbBcType::Neumann: {
    auto neumannEbBcFactory = RefCountedPtr<NeumannConductivityEBBCFactory> (new NeumannConductivityEBBCFactory());
    neumannEbBcFactory->setValue(m_ebbc.second);
    ebbcFactory = RefCountedPtr<BaseEBBCFactory> (neumannEbBcFactory);
    break;
  }
  case EbBcType::Larsen: {
    // Make the appropriate coefficients (Larsen) with right-hand side. 
    std::function<Real(const RealVect, const Real) > ebFunc = [&ebbc=this->m_ebbc](const RealVect a_pos, const Real a_time) {
      return ebbc.second;
    };
    
    auto robinCoefficients = RefCountedPtr<LarsenCoefficients> (new LarsenCoefficients(m_RtSpecies, m_reflectCoefOne, m_reflectCoefTwo, ebFunc));

    auto robinEbBcFactory = RefCountedPtr<RobinConductivityEbBcFactory> (new RobinConductivityEbBcFactory(m_amr->getProbLo()));
    robinEbBcFactory->setCoefficients(robinCoefficients);
    robinEbBcFactory->setStencilType(IrregStencil::StencilType::LeastSquares);
    ebbcFactory = RefCountedPtr<BaseEBBCFactory> (robinEbBcFactory);
    break;
  }
  default:
    MayDay::Abort("EddingtonSP1::setupOperatorFactory -- bad EBBC!");
    break;
  }


  // alpha and beta-coefficients for the Helmholtz equation. 
  const Real alpha =  1.0;
  const Real beta  = -1.0;

  // Number of ghost cells. 
  const int ghost = m_amr->getNumberOfGhostCells();


  // AmrMesh uses RefCounted levelgrids. EbHelmholtzOpFactory does not. Convert them there. 
  Vector<EBLevelGrid> levelgrids;
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){ 
    levelgrids.push_back(*(m_amr->getEBLevelGrid(m_realm, m_phase)[lvl])); 
  }

  // Create operator factory.
  m_operatorFactory = RefCountedPtr<EbHelmholtzOpFactory> (new EbHelmholtzOpFactory(levelgrids,
										    m_amr->getEBQuadCFInterp(m_realm, m_phase),
										    m_amr->getFluxRegister(m_realm, m_phase),
										    m_alpha,
										    m_beta,
										    m_aco.getData(),
										    m_bco.getData(),
										    m_bco_irreg.getData(),
										    m_amr->getDx()[0],
										    m_amr->getRefinementRatios(),
										    fact,
										    ebbcFactory,
										    ghost*IntVect::Unit,
										    ghost*IntVect::Unit,
										    2,
										    m_minCellsBottom,
										    -1,
										    m_mg_levelgrids));
}

void EddingtonSP1::setupMultigrid(){
  CH_TIME("EddingtonSP1::setupMultigrid");
  if(m_verbosity > 5){
    pout() << m_name + "::setupMultigrid" << endl;
  }

  // Select bottom solver
  LinearSolver<LevelData<EBCellFAB> >* botsolver = NULL;
  if(m_bottomSolverType == BottomSolverType::Simple){
    botsolver = &m_simpleSolver;
  }
  else if(m_bottomSolverType == BottomSolverType::BiCGStab){
    botsolver = &m_bicgstab;
  }

  // Make m_multigridType into an int for multigrid
  int gmgType;
  switch(m_multigridType){
  case MultigridType::VCycle:
    gmgType = 1;
    break;
  case MultigridType::WCycle:
    gmgType = 2;
    break;
  default:
    MayDay::Error("EddingtonSP1::setupMultigrid -- logic bust in multigrid type selection");
  }

  const int finestLevel              = m_amr->getFinestLevel();
  const ProblemDomain coarsestDomain = m_amr->getDomains()[0];

  // Define AMRMultiGrid
  m_multigridSolver = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > (new AMRMultiGrid<LevelData<EBCellFAB> >());
    m_multigridSolver->define(coarsestDomain, *m_operatorFactory, botsolver, 1 + finestLevel);
    //m_multigridSolver->define(coarsestDomain, *m_helmholtzOpFactory, botsolver, 1 + finestLevel);
  m_multigridSolver->setSolverParameters(m_multigridPreSmooth,
					 m_multigridPostSmooth,
					 m_multigridBottomSmooth,
					 gmgType,
					 m_multigridMaxIterations,
					 m_multigridExitTolerance,
					 m_multigridExitHang,
					 1.E-99); // Residue set through other means
  
  m_multigridSolver->m_imin      = m_multigridMinIterations;
  m_multigridSolver->m_verbosity = m_multigridVerbosity;

  // Dummies for init
  EBAMRCellData dummy1;
  EBAMRCellData dummy2;
  
  m_amr->allocate(dummy1, m_realm, m_phase, m_nComp);
  m_amr->allocate(dummy2, m_realm, m_phase, m_nComp);
  
  DataOps::setValue(dummy1, 0.0);
  DataOps::setValue(dummy2, 0.0);

  // Aliasing
  Vector<LevelData<EBCellFAB>* > phi;
  Vector<LevelData<EBCellFAB>* > rhs;

  m_amr->alias(phi, dummy1);
  m_amr->alias(rhs, dummy2);

  // Init solver. This instantiates all the operators in AMRMultiGrid so we can just call "solve"
  m_multigridSolver->init(phi, rhs, finestLevel, 0);
}

void EddingtonSP1::setupTGA(){
  CH_TIME("EddingtonSP1::setupTGA");
  if(m_verbosity > 5){
    pout() << m_name + "::setupTGA" << endl;
  }
  
  const int finestLevel              = m_amr->getFinestLevel();
  const ProblemDomain coarsestDomain = m_amr->getDomains()[0];
  const Vector<int> refRat           = m_amr->getRefinementRatios();

  m_tgaSolver = RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > >
    (new AMRTGA<LevelData<EBCellFAB> > (m_multigridSolver, *m_operatorFactory, coarsestDomain, refRat, 1 + finestLevel, m_multigridSolver->m_verbosity));
}

void EddingtonSP1::setupEuler(){
  CH_TIME("EddingtonSP1::setupEuler");
  if(m_verbosity > 5){
    pout() << m_name + "::setupEuler" << endl;
  }

  const int finestLevel              = m_amr->getFinestLevel();
  const ProblemDomain coarsestDomain = m_amr->getDomains()[0];
  const Vector<int> refRat           = m_amr->getRefinementRatios();

  m_eulerSolver = RefCountedPtr<EBBackwardEuler> 
    (new EBBackwardEuler (m_multigridSolver, *m_operatorFactory, coarsestDomain, refRat, 1 + finestLevel, m_multigridSolver->m_verbosity));
}

void EddingtonSP1::computeBoundaryFlux(EBAMRIVData& a_ebFlux, const EBAMRCellData& a_phi){
  CH_TIME("EddingtonSP1::computeBoundaryFlux");
  if(m_verbosity > 5){
    pout() << m_name + "::computeBoundaryFlux" << endl;
  }

  const int finestLevel = m_amr->getFinestLevel();
  
  const IrregAmrStencil<EbCentroidInterpolationStencil>& sten = m_amr->getEbCentroidInterpolationStencilStencils(m_realm, m_phase);
  for(int lvl = 0; lvl <= finestLevel; lvl++){
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
	      extrap(face, m_comp) = 0.0;
	    }
	    else{
	      if(!ebisbox.isCovered(iv1)){ // linear extrapolation
		extrap(face, m_comp) = 1.5*data_fab(iv0, m_comp) - 0.5*data_fab(iv1, m_comp); // Should be ok
	      }
	      else{ // Not enough cells available, use cell-centered only
		  extrap(face, m_comp) = data_fab(iv0, m_comp);
	      }
	    }

	    // Necessary scaling
	    extrap(face, m_comp) = 0.5*Units::c*extrap(face, m_comp);
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

  m_amr->computeGradient(a_flux, a_phi, m_realm, m_phase); // flux = grad(phi)
  
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    DataOps::divideByScalar(*a_flux[lvl], *m_aco[lvl]);    // flux = grad(phi)/(c*kappa)
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
  
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const Interval interv(m_comp, m_comp);
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

#include <CD_NamespaceFooter.H>
