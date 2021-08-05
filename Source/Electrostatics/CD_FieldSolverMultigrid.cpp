/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   FieldSolverMultigrid.cpp
  @brief  Implementation of FieldSolverMultigrid.H
  @author Robert Marskar
  @todo   When specifying the bottom solver, use syntax "simple 400" to indicate the number of smoothings. 
  @todo   gmg parameters should be specified as "FieldSolverMultigrid.gmg.<setting>"
  @todo   Purge the pre-coarsening stuff, e.g. defineDeeperMultigridLevels
  @todo   Register the new multigrid interpolators
  @todo   Once the new operator is in, check the computeLoads routine. 
*/

// Chombo includes
#include <ParmParse.H>
#include <BRMeshRefine.H>

// Our includes
#include <CD_FieldSolverMultigrid.H>
#include <CD_DataOps.H>
#include <CD_MFQuadCFInterp.H>
#include <CD_MFInterfaceFAB.H>
#include <CD_MFMultigridInterpolator.H>
#include <CD_MFHelmholtzDirichletEBBCFactory.H>
#include <CD_MFHelmholtzDirichletDomainBCFactory.H>
#include <CD_JumpBc.H>
#include <CD_AmrMesh.H>
#include <CD_ConductivityElectrostaticDomainBcFactory.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

constexpr Real FieldSolverMultigrid::m_alpha;
constexpr Real FieldSolverMultigrid::m_beta;

FieldSolverMultigrid::FieldSolverMultigrid() : FieldSolver() {

  // Default settings
  m_isSolverSetup            = false;
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

  pp.get("gmg_coarsen",      m_numCoarseningsBeforeAggregation);
  pp.get("gmg_verbosity",    m_multigridVerbosity);
  pp.get("gmg_pre_smooth",   m_multigridPreSmooth);
  pp.get("gmg_post_smooth",  m_multigridPostSmooth);
  pp.get("gmg_bott_smooth",  m_multigridBottomSmooth);
  pp.get("gmg_max_iter",     m_multigridMaxIterations);
  pp.get("gmg_min_iter",     m_multigridMinIterations);
  pp.get("gmg_exit_tol",     m_multigridExitTolerance);
  pp.get("gmg_exit_hang",    m_multigridExitHang);
  pp.get("gmg_min_cells",    m_minCellsBottom);
  pp.get("gmg_bc_order",     m_multigridBcOrder);

  if(!(m_multigridBcOrder == 1 || m_multigridBcOrder == 2)){
    MayDay::Abort("FieldSolverMultigrid::parseMultigridSettings - boundary condition order must be 1 or 2");
  }

  // Fetch the desired bottom solver from the input script. We look for things like FieldSolverMultigrid.gmg_bottom_solver = bicgstab or '= simple <number>'
  // where <number> is the number of relaxation for the smoothing solver. 
  const int num = pp.countval("gmg_bottom_solver");
  if(num == 1){
    pp.get("gmg_bottom_solver", str);
    if(str == "bicgstab"){
      m_bottomSolverType = BottomSolverType::BiCGStab;
    }
    else if(str == "gmres"){
      m_bottomSolverType = BottomSolverType::GMRES;
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
      m_mfsolver.setNumSmooths(numSmooth);
    }
    else{
      MayDay::Error("FieldSolverMultigrid::parseMultigridSettings - logic bust, you've specified two parameters and I expected 'simple <number>'");
    }
  }
  else{
    MayDay::Error("FieldSolverMultigrid::parseMultigridSettings - logic bust in bottom solver. You must specify ' = bicgstab', ' = gmres', or ' = simple <number>'");
  }

  // Set the multigrid relaxation type. 
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
  if(m_minCellsBottom < 2){
    m_minCellsBottom = 2;
  }
}

void FieldSolverMultigrid::allocateInternals(){
  CH_TIME("FieldSolverMultigrid::allocateInternals");
  if(m_verbosity > 5){
    pout() << "FieldSolverMultigrid::allocateInternals" << endl;
  }

  FieldSolver::allocateInternals();

  m_amr->allocate(m_zero,           m_realm,             m_nComp);
  m_amr->allocate(m_kappaRhoByEps0, m_realm,             m_nComp);
  m_amr->allocate(m_sigmaByEps0,    m_realm, phase::gas, m_nComp);

  DataOps::setValue(m_zero, 0.0);
}

bool FieldSolverMultigrid::solve(MFAMRCellData&       a_phi,
				 const MFAMRCellData& a_rho,
				 const EBAMRIVData&   a_sigma,
				 const bool           a_zerophi){
  CH_TIME("FieldSolverMultigrid::solve(mfamrcell, mfamrcell");
  if(m_verbosity > 5){
    pout() << "FieldSolverMultigrid::solve(mfamrcell, mfamrcell)" << endl;
  }

  // TLDR: The operator factory was set up as kappa*L(phi) = -kappa*div(eps*grad(phi)) which we use to solve the Poisson equation
  //       div(eps*grad(phi)) = -rho/eps0. 
  //
  //       So, we must scale rho (which is defined on centroids) by kappa and divide by eps0. The minus-sign is in the operator. 

  bool converged = false;

  if(!m_isSolverSetup){
    this->setupSolver(); // This does everything, allocates coefficients, gets bc stuff and so on
  }

  // Do the scaled space charge density
  DataOps::copy (m_kappaRhoByEps0, a_rho);
  DataOps::scale(m_kappaRhoByEps0, 1./(Units::eps0));
  DataOps::scale(m_kappaRhoByEps0, 1./(m_lengthScale*m_lengthScale));

  // Special flag for when a_rho is on the centroid but was not scaled by kappa on input. The multigrid operator solves
  // kappa*L(phi) = kappa*rho so the right-hand side must be kappa-weighted. 
  if(m_kappaSource){ 
    DataOps::kappaScale(m_kappaRhoByEps0);
  }
  
  // Do the scaled surface charge
  DataOps::copy (m_sigmaByEps0, a_sigma);
  DataOps::scale(m_sigmaByEps0, 1./(Units::eps0));
  DataOps::scale(m_sigmaByEps0, 1./(m_lengthScale*m_lengthScale));
  m_operatorFactory->setJump(m_sigmaByEps0, 1.0);

  // Aliasing
  Vector<LevelData<MFCellFAB>* > phi, cpy, rhs, res, zero;
  m_amr->alias(phi,     a_phi);
  m_amr->alias(rhs,     m_kappaRhoByEps0);
  m_amr->alias(res,     m_residue);
  m_amr->alias(zero,    m_zero);

  // GMG solve. Use phi = zero as initial metric. Want to reduce this by m_multigridExitTolerance
  const int finestLevel = m_amr->getFinestLevel();
  
  const Real phi_resid  = m_multigridSolver.computeAMRResidual(phi,  rhs, finestLevel, 0);
  const Real zero_resid = m_multigridSolver.computeAMRResidual(zero, rhs, finestLevel, 0);

  m_convergedResidue = zero_resid*m_multigridExitTolerance;

  // Do a multigrid solve if the residual is too large
  if(phi_resid > m_convergedResidue){ 
    m_multigridSolver.m_convergenceMetric = zero_resid;
    m_multigridSolver.solveNoInitResid(phi, res, rhs, finestLevel, 0, a_zerophi);

    const int status = m_multigridSolver.m_exitStatus;   // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8){                      // 8 => Norm sufficiently small
      converged = true;
    }
  }
  else{ // Solution is already converged
    converged = true;
  }

#if 1 // Solver hang. I don't know why this would work but OK.
  if(!converged){
    this->setupSolver();
    DataOps::setValue(a_phi, 0.0);
    m_multigridSolver.solveNoInitResid(phi, res, rhs, finestLevel, 0, a_zerophi);

    const int status = m_multigridSolver.m_exitStatus;   // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }
  }
#endif

  m_multigridSolver.revert(phi, rhs, finestLevel, 0);

  m_amr->averageDown(a_phi, m_realm);
  m_amr->interpGhost(a_phi, m_realm);

  this->computeElectricField();

  return converged;
}

void FieldSolverMultigrid::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){
  CH_TIME("FieldSolverMultigrid::regrid");
  if(m_verbosity > 5){
    pout() << "FieldSolverMultigrid::regrid" << endl;
  }
  
  FieldSolver::regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  
  m_isSolverSetup = false;
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
    m_amr->registerOperator(s_eb_coar_ave,     m_realm, phase::gas  );
    m_amr->registerOperator(s_eb_coar_ave,     m_realm, phase::solid);
    m_amr->registerOperator(s_eb_fill_patch,   m_realm, phase::gas  );
    m_amr->registerOperator(s_eb_fill_patch,   m_realm, phase::solid);
    m_amr->registerOperator(s_eb_pwl_interp,   m_realm, phase::gas  );
    m_amr->registerOperator(s_eb_pwl_interp,   m_realm, phase::solid);
    m_amr->registerOperator(s_eb_quad_cfi,     m_realm, phase::gas  );
    m_amr->registerOperator(s_eb_quad_cfi,     m_realm, phase::solid);
    m_amr->registerOperator(s_eb_irreg_interp, m_realm, phase::gas  );
    m_amr->registerOperator(s_eb_irreg_interp, m_realm, phase::solid);
    m_amr->registerOperator(s_eb_flux_reg,     m_realm, phase::gas  );
    m_amr->registerOperator(s_eb_flux_reg,     m_realm, phase::solid);
    m_amr->registerOperator(s_eb_multigrid,    m_realm, phase::gas  );
    m_amr->registerOperator(s_eb_multigrid,    m_realm, phase::solid);    
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
  this->setupHelmholtzFactory(); // Set up the operator factory
  this->setupOperatorFactory();  // Set the operator factory
  this->setupMultigrid();        // Set up the AMR multigrid solver

  m_isSolverSetup = true;
}

void FieldSolverMultigrid::setupHelmholtzFactory(){
  CH_TIME("FieldSolverMultigrid::setupHelmholtzFactory");
  if(m_verbosity > 5){
    pout() << "FieldSolverMultigrid::setupHelmholtzFactory" << endl;
  }

  pout() << "setting up helmholtz" << endl;

  // 
  const RefCountedPtr<EBIndexSpace>& ebisGas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebisSol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  // Need to set up multifluid wrappers
  const int finestLevel = m_amr->getFinestLevel();
  const int numPhases   = m_multifluidIndexSpace->numPhases();
  
  Vector<MFLevelGrid>             mflg     (1 + finestLevel);
  Vector<MFMultigridInterpolator> mfInterp (1 + finestLevel);
  Vector<MFFluxReg>               mfFluxReg(1 + finestLevel);
  Vector<MFCoarAve>               mfCoarAve(1 + finestLevel);

  for (int lvl = 0; lvl <= finestLevel; lvl++){
    Vector<EBLevelGrid>                             eblgPhases(numPhases);
    Vector<RefCountedPtr<EBMultigridInterpolator> > interpPhases(numPhases);
    Vector<RefCountedPtr<EBFluxRegister> >          fluxRegPhases(numPhases);
    Vector<RefCountedPtr<EbCoarAve> >               avePhases(numPhases);

    if(!ebisGas.isNull()) eblgPhases[phase::gas  ]    = *(m_amr->getEBLevelGrid(m_realm, phase::gas)  [lvl]);
    if(!ebisSol.isNull()) eblgPhases[phase::solid]    = *(m_amr->getEBLevelGrid(m_realm, phase::solid)[lvl]);

    if(!ebisGas.isNull()) interpPhases[phase::gas  ]  = (m_amr->getMultigridInterpolator(m_realm, phase::gas  )[lvl]);
    if(!ebisGas.isNull()) interpPhases[phase::solid]  = (m_amr->getMultigridInterpolator(m_realm, phase::solid)[lvl]);

    if(!ebisGas.isNull()) fluxRegPhases[phase::gas  ] = (m_amr->getFluxRegister(m_realm, phase::gas)  [lvl]);
    if(!ebisSol.isNull()) fluxRegPhases[phase::solid] = (m_amr->getFluxRegister(m_realm, phase::solid)[lvl]);

    if(!ebisGas.isNull()) avePhases[phase::gas  ]     = (m_amr->getCoarseAverage(m_realm, phase::gas)  [lvl]);
    if(!ebisSol.isNull()) avePhases[phase::solid]     = (m_amr->getCoarseAverage(m_realm, phase::solid)[lvl]);

    mflg[lvl].define(m_multifluidIndexSpace, eblgPhases);
    mfInterp[lvl].define(interpPhases);
    mfFluxReg[lvl].define(fluxRegPhases);
    mfCoarAve[lvl].define(avePhases);
  }

  // Coarsest domain used for multigrid
  ProblemDomain bottomDomain = m_amr->getDomains()[0];
  while(bottomDomain.domainBox().shortside() >= 2*m_minCellsBottom){
    bottomDomain.coarsen(2);
  }

  // Number of ghost cells in data holders
  const IntVect ghostPhi = m_amr->getNumberOfGhostCells()*IntVect::Unit;
  const IntVect ghostRhs = m_amr->getNumberOfGhostCells()*IntVect::Unit;

  // Set up Dirichlet boundary conditions for now
#if 1
  RefCountedPtr<MFHelmholtzEBBCFactory> ebbcFactory;
  RefCountedPtr<MFHelmholtzDomainBCFactory> domainBcFactory;

  int  jumpOrder = 2;  // Needs to be input parameter
  int  jumpWeight = 2; // Needs to be input parameter
  
  int  eb_order  = 2; // Needs to be input parameter
  int  eb_weight = 2; // Needs to be input parameter

  auto relaxType = MFHelmholtzOp::RelaxationMethod::GauSaiRedBlack; // Needs to be input parameter

  // dummies
  Real eb_value  = 1.0;
  auto bcFunction = [](const RealVect& a_pos) -> Real {
    return -1.0;
  };

  ebbcFactory     = RefCountedPtr<MFHelmholtzEBBCFactory>     (new MFHelmholtzDirichletEBBCFactory(eb_order, eb_weight, eb_value));
  domainBcFactory = RefCountedPtr<MFHelmholtzDomainBCFactory> (new MFHelmholtzDirichletDomainBCFactory(bcFunction));
#endif
  

#if 1
  m_helmholtzOpFactory = RefCountedPtr<MFHelmholtzOpFactory>(new MFHelmholtzOpFactory(m_multifluidIndexSpace,
										      m_dataLocation,
										      m_alpha,
										      m_beta,
										      m_amr->getProbLo(),
										      mflg,
										      mfInterp,
										      mfFluxReg,
										      mfCoarAve,
										      m_amr->getRefinementRatios(),
										      m_amr->getDx(),
										      m_permittivityCell.getData(), // Dummy argument (recall m_alpha = 0.0)
										      m_permittivityFace.getData(),
										      m_permittivityEB.getData(),
										      domainBcFactory,
										      ebbcFactory,
										      ghostPhi,
										      ghostRhs,
										      relaxType,
										      bottomDomain,
										      jumpOrder,
										      jumpWeight,
										      m_amr->getMaxBoxSize()));
										    
										    
#endif
}

void FieldSolverMultigrid::setupOperatorFactory(){
  CH_TIME("FieldSolverMultigrid::setupOperatorFactory");
  if(m_verbosity > 5){
    pout() << "FieldSolverMultigrid::setupOperatorFactory" << endl;
  }

  const int nphases                      = m_multifluidIndexSpace->numPhases();
  const int finestLevel                 = m_amr->getFinestLevel();
  const Vector<DisjointBoxLayout>& grids = m_amr->getGrids(m_realm);
  const Vector<int>& refinement_ratios   = m_amr->getRefinementRatios();
  const Vector<ProblemDomain>& domains   = m_amr->getDomains();
  const Vector<Real>& dx                 = m_amr->getDx();
  const RealVect& origin                 = m_amr->getProbLo();

  const RefCountedPtr<EBIndexSpace>& ebis_gas = m_multifluidIndexSpace->getEBIndexSpace(phase::gas);
  const RefCountedPtr<EBIndexSpace>& ebis_sol = m_multifluidIndexSpace->getEBIndexSpace(phase::solid);

  // This stuff is needed for the operator factory
  Vector<MFLevelGrid>    mflg(1 + finestLevel);
  Vector<MFQuadCFInterp> mfquadcfi(1 + finestLevel);
  Vector<MFFluxReg>  mffluxreg(1 + finestLevel);
  for (int lvl = 0; lvl <= finestLevel; lvl++){
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
										 m_minCellsBottom,
										 1 + finestLevel,
										 mg_levelgrids));

  // Parse Dirichlet boundary conditions on electrodes to the operator factory. 
  m_operatorFactory->setDirichletEbBc(m_ebBc);
}

void FieldSolverMultigrid::setupMultigrid(){
  CH_TIME("FieldSolverMultigrid::setupMultigrid");
  if(m_verbosity > 5){
    pout() << "FieldSolverMultigrid::setupMultigrid" << endl;
  }

  // Select bottom solver
  LinearSolver<LevelData<MFCellFAB> >* bottomSolver = NULL;
  switch(m_bottomSolverType){
  case BottomSolverType::Simple:
    bottomSolver = &m_mfsolver;
    break;
  case BottomSolverType::BiCGStab:
    bottomSolver = &m_bicgstab;
    break;
  case BottomSolverType::GMRES:
    bottomSolver = &m_gmres;
    break;
  default:
    MayDay::Error("FieldSolverMultigrid::setupMultigrid - logic bust in bottom solver");
    break;
  }

  // Select the multigrid type
  int gmgType;
  switch(m_multigridType){
  case MultigridType::VCycle:
    gmgType = 1;
    break;
  case MultigridType::WCycle:
    gmgType = 2;
    break;
  default:
    MayDay::Error("FieldSolverMultigrid::setupMultigrid - logic bust in multigrid type selection");
  }

#if 1   // This code should be removed later when the new operator is functional
  if(m_bottomSolverType == BottomSolverType::Simple){
    m_mfsolver.setNumSmooths(m_numSmoothingsForSimpleSolver);
    bottomSolver = &m_mfsolver;
  }
  else if(m_bottomSolverType == BottomSolverType::BiCGStab){
    bottomSolver = &m_bicgstab;
    if(m_multifluidIndexSpace->numPhases() == 2){ // BiCGStab doesn't work with multifluid (yet)
      bottomSolver = &m_mfsolver;
      
      if(m_verbosity > 0){
	pout() << "FieldSolverMultigrid::FieldSolverMultigrid - BiCGStab not supported for multifluid" << endl;
      }
    }
  }
  else if(m_bottomSolverType == BottomSolverType::GMRES){
    bottomSolver = &m_gmres;
  }
#endif

  const int finestLevel              = m_amr->getFinestLevel();
  const ProblemDomain coarsestDomain = m_amr->getDomains()[0];

  // Define the Chombo multigrid solver
  m_multigridSolver.define(coarsestDomain, *m_operatorFactory, bottomSolver, 1 + finestLevel);
  m_multigridSolver.setSolverParameters(m_multigridPreSmooth,
					m_multigridPostSmooth,
					m_multigridBottomSmooth,
					gmgType,
					m_multigridMaxIterations,
					m_multigridExitTolerance,
					m_multigridExitHang,
					1.E-99); // Norm thresh will be set via eps
  
  m_multigridSolver.m_imin      = m_multigridMinIterations;
  m_multigridSolver.m_verbosity = m_multigridVerbosity;


  // Create some dummy storage for multigrid initialization. 
  MFAMRCellData dummy1;
  MFAMRCellData dummy2;
  
  m_amr->allocate(dummy1, m_realm, m_nComp);
  m_amr->allocate(dummy2, m_realm, m_nComp);
  
  DataOps::setValue(dummy1, 0.0);
  DataOps::setValue(dummy2, 0.0);

  // Aliasing because Chombo is not too clever about smart pointers. 
  Vector<LevelData<MFCellFAB>* > phi;
  Vector<LevelData<MFCellFAB>* > rhs;

  m_amr->alias(phi, dummy1);
  m_amr->alias(rhs, dummy2);

  // Init the solver. This instantiates the all the operators. 
  m_multigridSolver.init(phi, rhs, finestLevel, 0);
}

Vector<long long> FieldSolverMultigrid::computeLoads(const DisjointBoxLayout& a_dbl, const int a_level) {
  CH_TIME("FieldSolverMultigrid::computeLoads");
  if(m_verbosity > 5){
    pout() << "FieldSolverMultigrid::computeLoads" << endl;
  }

  constexpr int numApply = 50;

  if(!m_isSolverSetup){
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
