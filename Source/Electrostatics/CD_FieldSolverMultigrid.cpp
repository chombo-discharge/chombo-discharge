/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   FieldSolverMultigrid.cpp
  @brief  Implementation of FieldSolverMultigrid.H
  @author Robert Marskar
  @todo   See if we can find a way to monkey the setup functions so that we don't have to explicitly cast in e.g. ItoPlasmaGodunovStepper. Maybe we can have
          FieldSolver take extra arguments representing the conductivity or somesuch, and add a setupSolver() virtual routine...?
  @todo   Replace vof iterators with iterations from AmrMesh
  @todo   When specifying the bottom solver, use syntax "simple 400" to indicate the number of smoothings. 
  @todo   gmg parameters should be specified as "FieldSolverMultigrid.gmg.<setting>"
  @todo   B-coefficient => face permittivities
  @todo   A-coefficient -> should be invisible for users
  @todo   Purge the pre-coarsening stuff. 
  @todo   Register the new multigrid interpolators
*/

// Chombo includes
#include <ParmParse.H>
#include <BRMeshRefine.H>

// Our includes
#include <CD_FieldSolverMultigrid.H>
#include <CD_DataOps.H>
#include <CD_MFQuadCFInterp.H>
#include <CD_MFInterfaceFAB.H>
#include <CD_JumpBc.H>
#include <CD_AmrMesh.H>
#include <CD_ConductivityElectrostaticDomainBcFactory.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

constexpr Real FieldSolverMultigrid::m_alpha;
constexpr Real FieldSolverMultigrid::m_beta;

#define POISSON_MF_GMG_TIMER 0

FieldSolverMultigrid::FieldSolverMultigrid() : FieldSolver() {
  m_needsMultigridSetup      = true;
  m_hasDeeperMultigridLevels = false;
  m_className                = "FieldSolverMultigrid";
}

FieldSolverMultigrid::~FieldSolverMultigrid(){

}

void FieldSolverMultigrid::parseOptions(){
  this->parseDomainBc();
  this->parsePlotVariables();
  this->parseMultigridSettings();
  this->parseKappaSource();
}

void FieldSolverMultigrid::parseRuntimeOptions(){
  this->parseMultigridSettings();
  this->parseKappaSource();
  this->parsePlotVariables();
}

void FieldSolverMultigrid::parseKappaSource(){
  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("kappa_source", m_kappaSource);
}

void FieldSolverMultigrid::parseMultigridSettings(){
  ParmParse pp(m_className.c_str());

  std::string str;

  pp.get("gmg_coarsen",     m_numCoarseningsBeforeAggregation);
  pp.get("gmg_verbosity",   m_multigridVerbosity);
  pp.get("gmg_pre_smooth",  m_multigridPreSmoothg);
  pp.get("gmg_post_smooth", m_multigridPostSmooth);
  pp.get("gmg_bott_smooth", m_multigridBottomSmooth);
  pp.get("gmg_max_iter",    m_multigridMaxIterations);
  pp.get("gmg_min_iter",    m_multigridMinIterations);
  pp.get("gmg_tolerance",   m_multigridTolerance);
  pp.get("gmg_hang",        m_multigridHang);
  pp.get("gmg_bottom_drop", m_numCellsBottomDrop);
  pp.get("gmg_bc_order",    m_multigridBcOrder);

  if(!(m_multigridBcOrder == 1 || m_multigridBcOrder == 2)){
    MayDay::Abort("FieldSolverMultigrid::parseMultigridSettings - boundary condition order must be 1 or 2");
  }

  // Bottom solver
  pp.get("gmg_bottom_solver", str);
  if(str == "simple"){
    m_bottomSolver = BottomSolver::Simple;
  }
  else if(str == "bicgstab"){
    m_bottomSolver = BottomSolver::BiCGStab;
  }
  else if(str == "gmres"){
    m_bottomSolver = BottomSolver::GMRES;
  }
  else{
    MayDay::Abort("FieldSolverMultigrid::parseMultigridSettings - unknown bottom solver requested");
  }

  // Relaxation type
  pp.get("gmg_relax_type", str);
  if(str == "gsrb"){
    m_multigridRelaxMethod = RelaxationMethod::GSRBFast;
  }
  else if( str == "gauss_seidel"){
    m_multigridRelaxMethod = RelaxationMethod::GaussSeidel;
  }
  else{
    MayDay::Abort("FieldSolverMultigrid::parseMultigridSettings - unsupported relaxation method requested");
  }

  // Cycle type
  pp.get("gmg_cycle", str);
  if(str == "vcycle"){
    m_multigridType = MultigridType::VCycle;
  }
  else{
    MayDay::Abort("FieldSolverMultigrid::parseMultigridSettings - unsupported multigrid cycle type requested. Only vcycle supported for now. ");
  }

  // No lower than 2. 
  if(m_numCellsBottomDrop < 2){
    m_numCellsBottomDrop = 2;
  }
}

void FieldSolverMultigrid::allocateInternals(){
  CH_TIME("FieldSolverMultigrid::allocateInternals");
  if(m_verbosity > 5){
    pout() << "FieldSolverMultigrid::allocateInternals" << endl;
  }

  FieldSolver::allocateInternals();

  const int ncomp = 1;

  m_amr->allocate(m_zero,         m_realm, ncomp);
  m_amr->allocate(m_scaledSource, m_realm, ncomp);
  m_amr->allocate(m_scaledSigma,  m_realm, phase::gas, ncomp);

  DataOps::setValue(m_zero, 0.0);
}



bool FieldSolverMultigrid::solve(MFAMRCellData&       a_phi,
				 const MFAMRCellData& a_source,
				 const EBAMRIVData&   a_sigma,
				 const bool           a_zerophi){
  CH_TIME("FieldSolverMultigrid::solve(mfamrcell, mfamrcell");
  if(m_verbosity > 5){
    pout() << "FieldSolverMultigrid::solve(mfamrcell, mfamrcell)" << endl;
  }

  bool converged = false;

  if(m_needsMultigridSetup){
    this->setupSolver(); // This does everything, allocates coefficients, gets bc stuff and so on
  }

  const Real t0 = MPI_Wtime();

  const int ncomp        = 1;
  const int finest_level = m_amr->getFinestLevel();


  const Real t1 = MPI_Wtime();

  // Do the scaled space charge density
  DataOps::copy(m_scaledSource, a_source);
  DataOps::scale(m_scaledSource, 1./(Units::eps0));
  DataOps::scale(m_scaledSource, 1./(m_lengthScale*m_lengthScale));

  if(m_kappaSource){ // Scale source by kappa
    DataOps::kappaScale(m_scaledSource);
  }

  // Do the scaled surface charge
  DataOps::copy(m_scaledSigma, a_sigma);
  DataOps::scale(m_scaledSigma, 1./(m_lengthScale*m_lengthScale));
  m_operatorFactory->setJump(m_scaledSigma, 1.0/Units::eps0);

  const Real t2 = MPI_Wtime();
  
#if 0 // Debug
  MayDay::Warning("FieldSolverMultigrid::solve - debug mode");
  m_operatorFactory->setJump(0.0, 1.0);
  DataOps::setValue(source, 0.0);
#endif

#if 0 // Check NaN/Inf input
  EBAMRCellData phiGas = m_amr->alias(phase::gas,   a_phi);
  EBAMRCellData phiSol = m_amr->alias(phase::solid, a_phi);

  EBAMRCellData srcGas = m_amr->alias(phase::gas,   m_scaledSource);
  EBAMRCellData srcSol = m_amr->alias(phase::solid, m_scaledSource);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    EBLevelDataOps::checkNANINF(*phiGas[lvl]);
    EBLevelDataOps::checkNANINF(*phiSol[lvl]);

    EBLevelDataOps::checkNANINF(*srcGas[lvl]);
    EBLevelDataOps::checkNANINF(*srcSol[lvl]);
  }
#endif

  // Aliasing
  Vector<LevelData<MFCellFAB>* > phi, cpy, rhs, res, zero;
  m_amr->alias(phi,     a_phi);
  m_amr->alias(rhs,     m_scaledSource);
  m_amr->alias(res,     m_residue);
  m_amr->alias(zero,    m_zero);

  // GMG solve. Use phi = zero as initial metric. Want to reduce this by m_multigridTolerance
  //  m_multigridSolver.init(phi, rhs, finest_level, 0);
  const Real phi_resid  = m_multigridSolver.computeAMRResidual(phi,  rhs, finest_level, 0);
  const Real zero_resid = m_multigridSolver.computeAMRResidual(zero, rhs, finest_level, 0);

  m_convergedResidue = zero_resid*m_multigridTolerance;

  const Real t3 = MPI_Wtime();

  if(phi_resid > m_convergedResidue){ // Residual is too large, recompute solution
    m_multigridSolver.m_convergenceMetric = zero_resid;
    m_multigridSolver.solveNoInitResid(phi, res, rhs, finest_level, 0, a_zerophi);

    const int status = m_multigridSolver.m_exitStatus;   // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8){  // 8 => Norm sufficiently small
      converged = true;
    }
  }
  else{ // Solution is already converged
    converged = true;
  }

  // Solver hang. Try again. 
  if(!converged){
    this->setupSolver();
    DataOps::setValue(a_phi, 0.0);
    m_multigridSolver.solveNoInitResid(phi, res, rhs, finest_level, 0, a_zerophi);

    const int status = m_multigridSolver.m_exitStatus;   // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }
  }

  const Real t4 = MPI_Wtime();

#if 0 // Why is this required??? Is it because of op->zeroCovered()????
  Real new_resid = m_multigridSolver.computeAMRResidual(phi, rhs, finest_level, 0);
  new_resid = m_multigridSolver.computeAMRResidual(phi, rhs, finest_level, 0);
#endif

  
#if 0 // Debug. Solve again
  m_multigridSolver.solveNoInitResid(phi, res, rhs, finest_level, 0, a_zerophi);
#endif


  m_multigridSolver.revert(phi, rhs, finest_level, 0);

  m_amr->averageDown(a_phi, m_realm);
  m_amr->interpGhost(a_phi, m_realm);

  const Real t5 = MPI_Wtime();

#if POISSON_MF_GMG_TIMER
  const Real T = t5 - t0;
  pout() << endl;
  pout() << "FieldSolverMultigrid::solve breakdown" << endl;
  pout() << "alloc :     " << 100.*(t1-t0)/T << "%" << endl;
  pout() << "set jump:   " << 100.*(t2-t1)/T << "%" << endl;
  pout() << "resid:      " << 100.*(t3-t2)/T << "%" << endl;
  pout() << "solve:      " << 100.*(t4-t3)/T << "%" << " = " << t4 - t3 << endl;
  pout() << "revert/avg: " << 100.*(t5-t4)/T << "%" << endl;
  pout() << "Total time: " << T << endl;
#endif

  this->computeElectricField();

  return converged;
}

void FieldSolverMultigrid::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){
  CH_TIME("FieldSolverMultigrid::regrid");
  if(m_verbosity > 5){
    pout() << "FieldSolverMultigrid::regrid" << endl;
  }
  
  FieldSolver::regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  
  m_needsMultigridSetup = true;
}

void FieldSolverMultigrid::registerOperators(){
  CH_TIME("FieldSolverMultigrid::registerOperators");
  if(m_verbosity > 5){
    pout() << "FieldSolverMultigrid::registerOperators" << endl;
  }

  if(m_amr.isNull()){
    MayDay::Abort("FieldSolverMultigrid::registerOperators - need to set AmrMesh!");
  }
  else{
    m_amr->registerOperator(s_eb_coar_ave,     m_realm, phase::gas);
    m_amr->registerOperator(s_eb_coar_ave,     m_realm, phase::solid);
    m_amr->registerOperator(s_eb_fill_patch,   m_realm, phase::gas);
    m_amr->registerOperator(s_eb_fill_patch,   m_realm, phase::solid);
    m_amr->registerOperator(s_eb_pwl_interp,   m_realm, phase::gas);
    m_amr->registerOperator(s_eb_pwl_interp,   m_realm, phase::solid);
    m_amr->registerOperator(s_eb_quad_cfi,     m_realm, phase::gas);
    m_amr->registerOperator(s_eb_quad_cfi,     m_realm, phase::solid);
    m_amr->registerOperator(s_eb_irreg_interp, m_realm, phase::gas);
    m_amr->registerOperator(s_eb_irreg_interp, m_realm, phase::solid);
    m_amr->registerOperator(s_eb_flux_reg,     m_realm, phase::gas);
    m_amr->registerOperator(s_eb_flux_reg,     m_realm, phase::solid);
  }
}

void FieldSolverMultigrid::defineDeeperMultigridLevels(){
  CH_TIME("FieldSolverMultigrid::define_mg_level");
  if(m_verbosity > 5){
    pout() << "FieldSolverMultigrid::defineDeeperMultigridLevels" << endl;
  }

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  const int coar_ref = 2;
  
  m_deeperMultigridDomains.resize(0);
  m_deeperMultigridGrids.resize(0);
  m_mg_mflg.resize(0);
  
  m_mg_eblg.resize(phase::numPhases);

  m_mg_eblg[phase::gas].resize(0);
  m_mg_eblg[phase::solid].resize(0);

  // Get some stuff from AmrMesh on how to decompose the levels
  const Vector<ProblemDomain>& domains = m_amr->getDomains();
  const int max_box_size               = m_amr->getMaxBoxSize();
  const int blocking_factor            = m_amr->getBlockingFactor();
  const int num_numEbGhostsCells                = m_amr->getNumberOfEbGhostCells();

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
      m_deeperMultigridDomains.push_back(coar);
      m_deeperMultigridGrids.push_back(DisjointBoxLayout(boxes, proc_assign, coar));

      // Define the EBLevelGrids
      const int idx = m_deeperMultigridGrids.size() - 1; // Last element added
      if(!ebis_gas.isNull()){
	m_mg_eblg[phase::gas].push_back(RefCountedPtr<EBLevelGrid> (new EBLevelGrid(m_deeperMultigridGrids[idx],
										    m_deeperMultigridDomains[idx],
										    num_numEbGhostsCells,
										    ebis_gas)));
      }
      if(!ebis_sol.isNull()){
	m_mg_eblg[phase::solid].push_back(RefCountedPtr<EBLevelGrid> (new EBLevelGrid(m_deeperMultigridGrids[idx],
										      m_deeperMultigridDomains[idx],
										      num_numEbGhostsCells,
										      ebis_sol)));
      }

      // Define the MFLevelGrid object
      Vector<EBLevelGrid> eblgs;
      if(!ebis_gas.isNull()) eblgs.push_back(*m_mg_eblg[phase::gas][idx]);
      if(!ebis_sol.isNull()) eblgs.push_back(*m_mg_eblg[phase::solid][idx]);
      
      m_mg_mflg.push_back(RefCountedPtr<MFLevelGrid> (new MFLevelGrid(m_multifluidIndexSpace, eblgs)));

      // Next iterate
      fine = coar;
      num_coar++;
    }
    else{
      break;
    }
  }
  
  m_hasDeeperMultigridLevels = true;
}

void FieldSolverMultigrid::setupSolver(){
  CH_TIME("FieldSolverMultigrid::setupSolver()");
  if(m_verbosity > 5){
    pout() << "FieldSolverMultigrid::setupSolver()" << endl;
  }

  if(!m_hasDeeperMultigridLevels){
    this->defineDeeperMultigridLevels();     // Define MG levels. These don't change during regrids so we only need to set them once. 
  }
  this->setupOperatorFactory(); // Set the operator factory
  this->setupMultigrid();           // Set up the AMR multigrid solver

  m_needsMultigridSetup = false;

}

void FieldSolverMultigrid::setupOperatorFactory(){
  CH_TIME("FieldSolverMultigrid::setupOperatorFactory");
  if(m_verbosity > 5){
    pout() << "FieldSolverMultigrid::setupOperatorFactory" << endl;
  }

  const int nphases                      = m_multifluidIndexSpace->numPhases();
  const int finest_level                 = m_amr->getFinestLevel();
  const Vector<DisjointBoxLayout>& grids = m_amr->getGrids(m_realm);
  const Vector<int>& refinement_ratios   = m_amr->getRefinementRatios();
  const Vector<ProblemDomain>& domains   = m_amr->getDomains();
  const Vector<Real>& dx                 = m_amr->getDx();
  const RealVect& origin                 = m_amr->getProbLo();

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  // This stuff is needed for the operator factory
  Vector<MFLevelGrid>    mflg(1 + finest_level);
  Vector<MFQuadCFInterp> mfquadcfi(1 + finest_level);
  Vector<MFFluxReg>  mffluxreg(1 + finest_level);
  for (int lvl = 0; lvl <= finest_level; lvl++){
    Vector<EBLevelGrid>                    eblg_phases   (nphases);
    Vector<RefCountedPtr<EBQuadCFInterp> > quadcfi_phases(nphases);
    Vector<RefCountedPtr<EBFluxRegister> > fluxreg_phases(nphases);

    if(!ebis_gas.isNull()) eblg_phases[phase::gas]      = *(m_amr->getEBLevelGrid(m_realm, phase::gas)[lvl]);
    if(!ebis_sol.isNull()) eblg_phases[phase::solid]    = *(m_amr->getEBLevelGrid(m_realm, phase::solid)[lvl]);

    if(!ebis_gas.isNull()) quadcfi_phases[phase::gas]   = (m_amr->getEBQuadCFInterp(m_realm, phase::gas)[lvl]);
    if(!ebis_sol.isNull()) quadcfi_phases[phase::solid] = (m_amr->getEBQuadCFInterp(m_realm, phase::solid)[lvl]);

    if(!ebis_gas.isNull()) fluxreg_phases[phase::gas]   = (m_amr->getFluxRegister(m_realm, phase::gas)[lvl]);
    if(!ebis_sol.isNull()) fluxreg_phases[phase::solid] = (m_amr->getFluxRegister(m_realm, phase::solid)[lvl]);
    
    mflg[lvl].define(m_multifluidIndexSpace, eblg_phases);
    mfquadcfi[lvl].define(quadcfi_phases);
    mffluxreg[lvl].define(fluxreg_phases);
  }

  // Appropriate coefficients for poisson equation

  const IntVect ghost_phi = m_amr->getNumberOfGhostCells()*IntVect::Unit;
  const IntVect ghost_rhs = m_amr->getNumberOfGhostCells()*IntVect::Unit;

  Vector<MFLevelGrid> mg_levelgrids;
  for (int lvl = 0; lvl < m_mg_mflg.size(); lvl++){
    mg_levelgrids.push_back(*m_mg_mflg[lvl]);
  }

  // Set the length scale for the Poisson equation. This is equivalent to solving the Poisson equation
  // on the domain [-1,1] in the x-direction
  m_lengthScale = 2./(domains[0].size(0)*dx[0]);

  int relax_type;
  if(m_multigridRelaxMethod == RelaxationMethod::Jacobi){
    relax_type = 0;
  }
  else if(m_multigridRelaxMethod == RelaxationMethod::GaussSeidel){
    relax_type = 1;
  }
  else if(m_multigridRelaxMethod == RelaxationMethod::GSRBFast){
    relax_type = 2;
  }

  auto domfact = RefCountedPtr<ConductivityElectrostaticDomainBcFactory> (new ConductivityElectrostaticDomainBcFactory(m_domainBc, m_amr->getProbLo()));

  // Create factory and set potential
  m_operatorFactory = RefCountedPtr<MfHelmholtzOpFactory> (new MfHelmholtzOpFactory(m_multifluidIndexSpace,
										 mflg,
										 mfquadcfi,
										 mffluxreg,
										 refinement_ratios,
										 grids,
										 m_permittivityCell,
										 m_permittivityFace,
										 m_permittivityEB,
										 m_alpha,
										 m_beta,
										 m_lengthScale,
										 dx[0]*m_lengthScale,
										 domains[0],
										 domfact,
										 origin,
										 ghost_phi,
										 ghost_rhs,
										 m_multigridBcOrder,
										 relax_type,
										 m_numCellsBottomDrop,
										 1 + finest_level,
										 mg_levelgrids));

  // Parse Dirichlet boundary conditions on electrodes to the operator factory. 
  m_operatorFactory->setDirichletEbBc(m_ebBc);
}

void FieldSolverMultigrid::setupMultigrid(){
  CH_TIME("FieldSolverMultigrid::setupMultigrid");
  if(m_verbosity > 5){
    pout() << "FieldSolverMultigrid::setupMultigrid" << endl;
  }

  const int finest_level       = m_amr->getFinestLevel();
  const ProblemDomain coar_dom = m_amr->getDomains()[0];

  // Select bottom solver
  LinearSolver<LevelData<MFCellFAB> >* botsolver = NULL;
  
  if(m_bottomSolver == BottomSolver::Simple){
    m_mfsolver.setNumSmooths(m_numSmoothingsForSimpleSolver);
    botsolver = &m_mfsolver;
  }
  else if(m_bottomSolver == BottomSolver::BiCGStab){
    botsolver = &m_bicgstab;
    if(m_multifluidIndexSpace->numPhases() == 2){ // BiCGStab doesn't work with multifluid (yet)
      botsolver = &m_mfsolver;
      
      if(m_verbosity > 0){
	pout() << "FieldSolverMultigrid::FieldSolverMultigrid - BiCGStab not supported for multifluid" << endl;
      }
    }
  }
  else if(m_bottomSolver == BottomSolver::GMRES){
    botsolver = &m_gmres;
  }

  int gmg_type;
  if(m_multigridType == MultigridType::VCycle){
    gmg_type = 1;
  }
  else if(m_multigridType == MultigridType::WCycle){
    gmg_type = 2;
  }
  else{
    MayDay::Abort("FieldSolverMultigrid::setupMultigrid - unknown cycle type requested");
  }

  m_multigridSolver.define(coar_dom, *m_operatorFactory, botsolver, 1 + finest_level);
  m_multigridSolver.setSolverParameters(m_multigridPreSmoothg,
					m_multigridPostSmooth,
					m_multigridBottomSmooth,
					gmg_type,
					m_multigridMaxIterations,
					m_multigridTolerance,
					m_multigridHang,
					1.E-99); // Norm thresh will be set via eps
  
  m_multigridSolver.m_imin      = m_multigridMinIterations;
  m_multigridSolver.m_verbosity = m_multigridVerbosity;


  // Dummies for init
  const int ncomp = 1;
  MFAMRCellData dummy1, dummy2;
  m_amr->allocate(dummy1, m_realm, ncomp);
  m_amr->allocate(dummy2, m_realm, ncomp);
  DataOps::setValue(dummy1, 0.0);
  DataOps::setValue(dummy2, 0.0);

  // Aliasing
  Vector<LevelData<MFCellFAB>* > phi, rhs;
  m_amr->alias(phi, dummy1);
  m_amr->alias(rhs, dummy2);

  // Init solver
  m_multigridSolver.init(phi, rhs, finest_level, 0);
}

Vector<long long> FieldSolverMultigrid::computeLoads(const DisjointBoxLayout& a_dbl, const int a_level) {
  CH_TIME("FieldSolverMultigrid::computeLoads");
  if(m_verbosity > 5){
    pout() << "FieldSolverMultigrid::computeLoads" << endl;
  }

  constexpr int numApply = 50;

  if(m_needsMultigridSetup){
    this->setupSolver();
  }

  // Dummy storage where the multigrid operator will do relaxations. 
  MFAMRCellData dummy1, dummy2;

  m_amr->allocate(dummy1, m_realm, 1);
  m_amr->allocate(dummy2, m_realm, 1);

  // Make an operator
  auto oper = (MfHelmholtzOp*) m_operatorFactory->MGnewOp(m_amr->getDomains()[a_level], 0, false);

  // Reset time
  TimedDataIterator dit = a_dbl.timedDataIterator();
  dit.clearTime();
  dit.enableTime();

  // Do some relaxations on each level. This includes BCs.
  for (int k = 0; k < numApply; k++){
    oper->applyOp(*dummy1[a_level], *dummy2[a_level], dit, true);
  }

  // Merge times
  dit.disableTime();
  dit.mergeTime();

  // Now do the load balancing. When we do this, the boxes are standard-sorted! I.e. the new mesh ignores the desired sorting from AmrMesh
  Vector<unsigned long long> loads = dit.getTime();

  Vector<long long> ret(loads.size());
  for (int i = 0; i < loads.size(); i++){
    ret[i] = (long long) loads[i];
  }

  delete oper;

  return ret;
}

#include <CD_NamespaceFooter.H>
