/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrTGA.cpp
  @brief  Implementation of CD_CdrTGA.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>
#include <EBAMRIO.H>

// Our includes
#include <CD_CdrTGA.H>
#include <CD_DataOps.H>
#include <CD_EBHelmholtzNeumannDomainBCFactory.H>
#include <CD_EBHelmholtzNeumannEBBCFactory.H>
#include <CD_NamespaceHeader.H>

CdrTGA::CdrTGA() : CdrSolver() {

  // Default settings
  m_name         = "CdrTGA";
  m_className    = "CdrTGA";
}

CdrTGA::~CdrTGA(){

}

void CdrTGA::registerOperators() {
  if(m_amr.isNull()){
    MayDay::Abort("CdrGA::registerOperators - need to set AmrMesh!");
  }
  else{

    // Need to register everything that the base class registered, plus the amr interpolator
    CdrSolver::registerOperators();
    
    m_amr->registerOperator(s_eb_multigrid, m_realm, m_phase);
  }
}



void CdrTGA::advanceEuler(EBAMRCellData&       a_newPhi,
			  const EBAMRCellData& a_oldPhi,
			  const EBAMRCellData& a_source,
			  const Real           a_dt){
  CH_TIME("CdrTGA::advanceEuler");
  if(m_verbosity > 5){
    pout() << m_name + "::advanceEuler" << endl;
  }
  
  if(m_isDiffusive){
    this->setupDiffusionSolver(); // Set up gmg again since diffusion coefficients might change 
    
    bool converged = false;

    const int comp         = 0;
    const int ncomp        = 1;
    const int finestLevel = m_amr->getFinestLevel();

    // Set convergence metric. Make zero = 0 and m_scratch = phi^k + dt*S^k. Then
    // compute the residue L(zero) - m_scratch which sets the convergence metric. 
    EBAMRCellData zero;
    m_amr->allocate(zero, m_realm, m_phase, m_nComp);
    DataOps::setValue(zero, 0.0);
    DataOps::copy(m_scratch, a_oldPhi);
    DataOps::incr(m_scratch, a_source, a_dt);

    Vector<LevelData<EBCellFAB>* > orez;
    Vector<LevelData<EBCellFAB>* > shr;
    m_amr->alias(orez, zero);
    m_amr->alias(shr,  m_scratch);

    m_eulerSolver->resetAlphaAndBeta(1.0, -a_dt);
    
    const Real zero_resid = m_multigridSolver->computeAMRResidual(orez, shr, finestLevel, 0);
    const Real stopcrit   = zero_resid;
    m_multigridSolver->m_convergenceMetric = zero_resid;


    // Now do the solve. 
    Vector<LevelData<EBCellFAB>* > new_state, old_state, source;
    m_amr->alias(new_state, a_newPhi);
    m_amr->alias(old_state, a_oldPhi);
    m_amr->alias(source,    a_source);


    // Euler solve
    DataOps::setValue(m_ebCenteredDiffusionCoefficient, 0.0);
    m_eulerSolver->oneStep(new_state, old_state, source, a_dt, 0, finestLevel, false);
    const int status = m_multigridSolver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }
  }
  else{
    DataOps::copy(a_newPhi, a_oldPhi);
  }
}



void CdrTGA::advanceTGA(EBAMRCellData&       a_newPhi,
			const EBAMRCellData& a_oldPhi,
			const EBAMRCellData& a_source,
			const Real           a_dt){
  CH_TIME("CdrTGA::advanceTGA(full)");
  if(m_verbosity > 5){
    pout() << m_name + "::advanceTGA(full)" << endl;
  }
  
  if(m_isDiffusive){
    this->setupDiffusionSolver(); // Set up gmg again since diffusion coefficients might change
    
    bool converged = false;

    const int comp         = 0;
    const int ncomp        = 1;
    const int finestLevel = m_amr->getFinestLevel();

    // Do the aliasing stuff
    Vector<LevelData<EBCellFAB>* > new_state;
    Vector<LevelData<EBCellFAB>* > old_state;
    Vector<LevelData<EBCellFAB>* > source;

    m_amr->alias(new_state, a_newPhi);
    m_amr->alias(old_state, a_oldPhi);
    m_amr->alias(source,    a_source);

    const Real alpha = 0.0;
    const Real beta  = 1.0;

    DataOps::setValue(m_ebCenteredDiffusionCoefficient, 0.0);

    // TGA solve
    m_tgaSolver->resetAlphaAndBeta(alpha, beta); // Needed...?
    m_tgaSolver->oneStep(new_state, old_state, source, a_dt, 0, finestLevel, false);

    const int status = m_multigridSolver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }
  }
  else{
    DataOps::copy(a_newPhi, a_oldPhi);
  }
}

void CdrTGA::setupDiffusionSolver(){
  CH_TIME("CdrTGA::setupDiffusionSolver");
  if(m_verbosity > 5){
    pout() << m_name + "::setupDiffusionSolver" << endl;
  }

  this->setupHelmholtzFactory();
  this->setupMultigrid();
  this->setupTGA();
  this->setupEuler();
}

void CdrTGA::setupHelmholtzFactory(){
  CH_TIME("CdrTGA::setupHelmholtzFactory");
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

  // Handle to boundary condition factories for domain and EB in an EBHelmholtzOp context.
  auto domainBcFactory = RefCountedPtr<EBHelmholtzNeumannDomainBCFactory> (new EBHelmholtzNeumannDomainBCFactory(0.0));
  auto ebbcFactory     = RefCountedPtr<EBHelmholtzNeumannEBBCFactory>     (new EBHelmholtzNeumannEBBCFactory    (0.0));

  // Temp alpha and beta. Diffusion solvers will reset these later. 
  const Real alpha = 1.0;
  const Real beta  = 1.0;

  // Transient storage which will be removed when we transition to Poisson operator
  EBAMRCellData m_helmAcoef;
  m_amr->allocate(m_helmAcoef, m_realm, m_phase, 1);
  DataOps::setValue(m_helmAcoef, 1.0);

  // Set up the operator
  m_poissonOpFactory = RefCountedPtr<EBHelmholtzOpFactory> (new EBHelmholtzOpFactory(Location::Cell::Center,
										       alpha,
										       beta,
										       m_amr->getProbLo(),
										       levelGrids,
										       interpolator,
										       fluxReg,
										       coarAve,
										       m_amr->getRefinementRatios(),
										       m_amr->getDx(),
										       m_helmAcoef.getData(),
										       m_faceCenteredDiffusionCoefficient.getData(),
										       m_ebCenteredDiffusionCoefficient.getData(),
										       domainBcFactory,
										       ebbcFactory,
										       ghostPhi,
										       ghostRhs,
										       m_smoother,
										       bottomDomain,
										       m_amr->getMaxBoxSize()));
}

void CdrTGA::setupMultigrid(){
  CH_TIME("CdrTGA::setupMultigrid");
  if(m_verbosity > 5){
    pout() << m_name + "::setupMultigrid" << endl;
  }

  // Select the bottom solver
  LinearSolver<LevelData<EBCellFAB> >* botsolver = NULL;
  switch(m_bottomSolverType){
  case BottomSolverType::Simple:
    botsolver = &m_simpleSolver;
    break;
  case BottomSolverType::BiCGStab:
    botsolver = &m_bicgstab;
    break;
  case BottomSolverType::GMRES:
    botsolver = &m_gmres;
    m_gmres.m_verbosity = 0; // Shut up. 
  default:
    MayDay::Error("CdrTGA::setupMultigrid - logic bust in bottom solver setup");
    break;
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
    MayDay::Error("CdrTGA::setupMultigrid -- logic bust in multigrid type selection");
  }

  const int finestLevel              = m_amr->getFinestLevel();
  const ProblemDomain coarsestDomain = m_amr->getDomains()[0];

  // Define AMRMultiGrid
  m_multigridSolver = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > (new AMRMultiGrid<LevelData<EBCellFAB> >());
  m_multigridSolver->define(coarsestDomain, *m_poissonOpFactory, botsolver, 1 + finestLevel);
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

void CdrTGA::setupTGA(){
  CH_TIME("CdrTGA::setupTGA");
  if(m_verbosity > 5){
    pout() << m_name + "::setupTGA" << endl;
  }
  
  const int finestLevel       = m_amr->getFinestLevel();
  const ProblemDomain coar_dom = m_amr->getDomains()[0];
  const Vector<int> ref_rat    = m_amr->getRefinementRatios();

  m_tgaSolver = RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > >
    (new AMRTGA<LevelData<EBCellFAB> > (m_multigridSolver, *m_poissonOpFactory, coar_dom, ref_rat, 1 + finestLevel, m_multigridSolver->m_verbosity));
}

void CdrTGA::setupEuler(){
  CH_TIME("CdrTGA::setupEuler");
  if(m_verbosity > 5){
    pout() << m_name + "::setupEuler" << endl;
  }

  const int finestLevel       = m_amr->getFinestLevel();
  const ProblemDomain coar_dom = m_amr->getDomains()[0];
  const Vector<int> ref_rat    = m_amr->getRefinementRatios();

  m_eulerSolver = RefCountedPtr<EBBackwardEuler> 
    (new EBBackwardEuler (m_multigridSolver, *m_poissonOpFactory, coar_dom, ref_rat, 1 + finestLevel, m_multigridSolver->m_verbosity));
}

void CdrTGA::computeDivJ(EBAMRCellData& a_divJ, EBAMRCellData& a_phi, const Real a_extrapDt, const bool a_ebFlux){
  CH_TIME("CdrTGA::computeDivJ(divF, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDivJ(divF, state)" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;

  // Fill ghost cells
  m_amr->interpGhostPwl(a_phi, m_realm, m_phase);

  if(m_isMobile || m_isDiffusive){

    if(m_useMassWeightedRedistribution){
      this->setRedistWeights(a_phi);
    }

    // We will let m_scratchFluxOne hold the total flux = advection + diffusion fluxes
    DataOps::setValue(m_scratchFluxOne, 0.0);

    // Compute advection flux. This is mostly the same as computeDivF. If we can, add domain fluxes here. 
    if(m_isMobile){
      m_amr->interpGhostPwl(m_cellVelocity, m_realm, m_phase);
      
      this->averageVelocityToFaces(); // Update m_faceVelocity from m_cellVelocity
      this->advectToFaces(m_faceStates, a_phi, a_extrapDt); // Advect to faces
      this->computeAdvectionFlux(m_scratchFluxTwo, m_faceVelocity, m_faceStates, m_domainBc);

      DataOps::incr(m_scratchFluxOne, m_scratchFluxTwo, 1.0);
    }

    // Compute diffusion flux. If we don't have advection, add the domain flux here. 
    if(m_isDiffusive){
      if(m_isMobile){
	this->computeDiffusionFlux(m_scratchFluxTwo, a_phi, CdrBc::Wall); // Domain flux already in advective flux
      }
      else{
	this->computeDiffusionFlux(m_scratchFluxTwo, a_phi, m_domainBc); // No advective deriv. Put flux here. 
      }
      DataOps::incr(m_scratchFluxOne, m_scratchFluxTwo, -1.0);
    }

    // General divergence computation. Also inject charge. Domain fluxes came in above but eb fluxes come in here. 
    EBAMRIVData* ebflux;
    if(a_ebFlux){
      ebflux = &m_ebFlux;
    }
    else{
      ebflux = &m_ebZero;
    }
    this->computeDivG(a_divJ, m_scratchFluxOne, *ebflux);
  }
  else{ 
    DataOps::setValue(a_divJ, 0.0);
  }

  return;
}

void CdrTGA::computeDivF(EBAMRCellData& a_divF, EBAMRCellData& a_phi, const CdrBc a_whichBc, const Real a_extrapDt, const bool a_ebFlux){
  CH_TIME("CdrTGA::computeDivF(divF, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDivF(divF, state)" << endl;
  }

  if(m_isMobile){

    // Fill ghost cells
    m_amr->interpGhostPwl(a_phi,     m_realm, m_phase);
    m_amr->interpGhostPwl(m_cellVelocity, m_realm, m_phase);

    if(m_useMassWeightedRedistribution){
      this->setRedistWeights(a_phi);
    }
    this->averageVelocityToFaces();                                                         // Cell-centered velocities become face-centered velocities. 
    this->advectToFaces(m_faceStates, a_phi, a_extrapDt);                                   // Face extrapolation to cell-centered faces
    this->computeAdvectionFlux(m_scratchFluxOne, m_faceVelocity, m_faceStates, a_whichBc);  // Compute face-centered fluxes

    EBAMRIVData* ebflux;
    if(a_ebFlux){
      ebflux = &m_ebFlux;
    }
    else{
      ebflux = &m_ebZero;
    }
    this->computeDivG(a_divF, m_scratchFluxOne, *ebflux); 
  }
  else{
    DataOps::setValue(a_divF, 0.0);
  }

#if 0 // Debug
  Real volume;
  const Real sum = EBLevelDataOps::noKappaSumLevel(volume, *a_divF[0], 0, m_amr->getDomains()[0]);

  std::cout << sum << std::endl;
#endif

}

void CdrTGA::computeDivD(EBAMRCellData& a_divD, EBAMRCellData& a_phi, const CdrBc a_domainBc, const bool a_ebFlux){
  CH_TIME("CdrTGA::computeDivD");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDivD" << endl;
  }

  if(m_isDiffusive){
    const int ncomp = 1;

    // Fill ghost cells
    m_amr->interpGhostPwl(a_phi, m_realm, m_phase);

    if(m_useMassWeightedRedistribution){
      this->setRedistWeights(a_phi);
    }

    this->computeDiffusionFlux(m_scratchFluxOne, a_phi, a_domainBc);  // Compute the face-centered diffusion flux

    EBAMRIVData* ebflux;
    if(a_ebFlux){
      ebflux = &m_ebFlux;
    }
    else{
      ebflux = &m_ebZero;
    }
    this->computeDivG(a_divD, m_scratchFluxOne, *ebflux); // General face-centered flux to divergence magic.

    m_amr->averageDown(a_divD, m_realm, m_phase);
    m_amr->interpGhost(a_divD, m_realm, m_phase);
  }
  else{
    DataOps::setValue(a_divD, 0.0);
  }
}

void CdrTGA::parseMultigridSettings(){
  ParmParse pp(m_className.c_str());

  std::string str;

  pp.get("gmg_verbosity",    m_multigridVerbosity);
  pp.get("gmg_pre_smooth",   m_multigridPreSmooth);
  pp.get("gmg_post_smooth",  m_multigridPostSmooth);
  pp.get("gmg_bott_smooth",  m_multigridBottomSmooth);
  pp.get("gmg_max_iter",     m_multigridMaxIterations);
  pp.get("gmg_min_iter",     m_multigridMinIterations);
  pp.get("gmg_exit_tol",     m_multigridExitTolerance);
  pp.get("gmg_exit_hang",    m_multigridExitHang);
  pp.get("gmg_min_cells",    m_minCellsBottom);

  // Fetch the desired bottom solver from the input script. We look for things like CdrTGA.gmg_bottom_solver = bicgstab or '= simple <number>'
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
      MayDay::Error("CdrTGA::parseMultigridSettings - logic bust, you've specified one parameter and I expected either 'bicgstab' or 'gmres'");
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
      MayDay::Error("CdrTGA::parseMultigridSettings - logic bust, you've specified two parameters and I expected 'simple <number>'");
    }
  }
  else{
    MayDay::Error("CdrTGA::parseMultigridSettings - logic bust in bottom solver. You must specify ' = bicgstab', ' = gmres', or ' = simple <number>'");
  }

  // Relaxation type
  pp.get("gmg_smoother", str);
  if(str == "jacobi"){
    m_smoother = EBHelmholtzOp::Smoother::PointJacobi;    
  }
  else if(str == "red_black"){
    m_smoother = EBHelmholtzOp::Smoother::GauSaiRedBlack;    
  }
  else if(str == "multi_color"){
    m_smoother = EBHelmholtzOp::Smoother::GauSaiMultiColor;
  }
  else{
    MayDay::Error("CdrTGA::parseMultigridSettings - unknown relaxation method requested");
  }

  // Cycle type
  pp.get("gmg_cycle", str);
  if(str == "vcycle"){
    m_multigridType = MultigridType::VCycle;
  }
  else{
    MayDay::Error("CdrTGA::parseMultigridSettings - unknown cycle type requested");
  }

  // No lower than 2. 
  if(m_minCellsBottom < 2){
    m_minCellsBottom = 2;
  }
}

void CdrTGA::writePlotData(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("CdrTGA::writePlotData");
  if(m_verbosity > 5){
    pout() << m_name + "::writePlotData" << endl;
  }

  if(m_plotPhi) {
    if(!m_plotNumbers){ // Regular write
      this->writeData(a_output, a_comp, m_phi, true);
    }
    else{ // Scale, write, and scale back
     
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
	DataOps::scale(*m_phi[lvl], (pow(m_amr->getDx()[lvl], 3)));
      }
      this->writeData(a_output, a_comp, m_phi,    false);
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
	DataOps::scale(*m_phi[lvl], 1./(pow(m_amr->getDx()[lvl], 3)));
      }
    }
  }

  if(m_plotDiffusionCoefficient && m_isDiffusive) { // Need to compute the cell-centerd stuff first
    DataOps::setValue(m_scratch, 0.0);
    DataOps::averageFaceToCell(m_scratch, m_faceCenteredDiffusionCoefficient, m_amr->getDomains());
    
    this->writeData(a_output, a_comp, m_scratch,false);
  }

  if(m_plotSource) {
    if(!m_plotNumbers){
      this->writeData(a_output, a_comp, m_source, false);
    }
    else { // Scale, write, and scale back
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
	DataOps::scale(*m_source[lvl], (pow(m_amr->getDx()[lvl], 3)));
      }
      writeData(a_output, a_comp, m_source,    false);
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
	DataOps::scale(*m_source[lvl], 1./(pow(m_amr->getDx()[lvl], 3)));
      }
    }
  }

  if(m_plotVelocity && m_isMobile) {
    this->writeData(a_output, a_comp, m_cellVelocity, false);
  }

  // Plot EB fluxes
  if(m_plotEbFlux && m_isMobile){
    DataOps::setValue(m_scratch, 0.0);
    DataOps::incr(m_scratch, m_ebFlux, 1.0);
    
    this->writeData(a_output, a_comp, m_scratch, false);
  }
}

#include <CD_NamespaceFooter.H>
