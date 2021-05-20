/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_CdrTGA.cpp
  @brief  Implementation of CD_CdrTGA.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>
#include <EBArith.H>
#include <EBAMRIO.H>
#include <NeumannConductivityEBBC.H>
#include <NeumannConductivityDomainBC.H>
#include <BRMeshRefine.H>

// Our includes
#include <CD_CdrTGA.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

CdrTGA::CdrTGA() : CdrSolver() {
  m_name         = "CdrTGA";
  m_className    = "CdrTGA";
  m_hasDeeperMultigridLevels = false;
}

CdrTGA::~CdrTGA(){

}

void CdrTGA::advanceEuler(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const Real a_dt){
  CH_TIME("CdrTGA::advanceEuler(no source)");
  if(m_verbosity > 5){
    pout() << m_name + "::advanceEuler(no source)" << endl;
  }
  
  if(m_isDiffusive){
    // Create a source term = S = 0.0;
    const int ncomp = 1;
    EBAMRCellData src;
    m_amr->allocate(src, m_realm, m_phase, ncomp);
    DataOps::setValue(src, 0.0);

    // Call version with source term
    advanceEuler(a_newPhi, a_oldPhi, src, a_dt);
  }
  else{
    DataOps::copy(a_newPhi, a_oldPhi);
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
    this->setupMultigrid(); // Set up gmg again since diffusion coefficients might change 
    
    bool converged = false;

    const int comp         = 0;
    const int ncomp        = 1;
    const int finest_level = m_amr->getFinestLevel();

    // Set convergence metric. Make zero = 0 and m_scratch = phi^k + dt*S^k. Then
    // compute the residue L(zero) - m_scratch which sets the convergence metric. 
    EBAMRCellData zero;
    m_amr->allocate(zero, m_realm, m_phase, ncomp);
    DataOps::setValue(zero, 0.0);
    DataOps::copy(m_scratch, a_oldPhi);
    DataOps::incr(m_scratch, a_source, a_dt);

    Vector<LevelData<EBCellFAB>* > orez;
    Vector<LevelData<EBCellFAB>* > shr;
    m_amr->alias(orez, zero);
    m_amr->alias(shr,  m_scratch);

    m_eulerSolver->resetAlphaAndBeta(1.0, -a_dt);
    
    const Real zero_resid = m_multigridSolver->computeAMRResidual(orez, shr, finest_level, 0);
    const Real stopcrit   = zero_resid;
    m_multigridSolver->m_convergenceMetric = zero_resid;


    // Now do the solve. 
    Vector<LevelData<EBCellFAB>* > new_state, old_state, source;
    m_amr->alias(new_state, a_newPhi);
    m_amr->alias(old_state, a_oldPhi);
    m_amr->alias(source,    a_source);


    // Euler solve
    DataOps::setValue(m_ebCenteredDiffusionCoefficient, 0.0);
    m_eulerSolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, false);
    const int status = m_multigridSolver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }
  }
  else{
    DataOps::copy(a_newPhi, a_oldPhi);
  }
}

void CdrTGA::advanceTGA(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const Real a_dt){
  CH_TIME("CdrTGA::advanceTGA(no source)");
  if(m_verbosity > 5){
    pout() << m_name + "::advanceTGA(no source)" << endl;
  }

  if(m_isDiffusive){

    // Dummy source term
    const int ncomp = 1;
    EBAMRCellData src;
    m_amr->allocate(src, m_realm, m_phase, ncomp);
    DataOps::setValue(src, 0.0);

    // Call other version
    advanceTGA(a_newPhi, a_oldPhi, src, a_dt);
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
    this->setupMultigrid(); // Set up gmg again since diffusion coefficients might change
    
    bool converged = false;

    const int comp         = 0;
    const int ncomp        = 1;
    const int finest_level = m_amr->getFinestLevel();

    // Do the aliasing stuff
    Vector<LevelData<EBCellFAB>* > new_state, old_state, source;
    m_amr->alias(new_state, a_newPhi);
    m_amr->alias(old_state, a_oldPhi);
    m_amr->alias(source,    a_source);

    const Real alpha = 0.0;
    const Real beta  = 1.0;

    DataOps::setValue(m_ebCenteredDiffusionCoefficient, 0.0);

    // TGA solve
    m_tgaSolver->resetAlphaAndBeta(alpha, beta);
    m_tgaSolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, false);

    const int status = m_multigridSolver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }
  }
  else{
    DataOps::copy(a_newPhi, a_oldPhi);
  }
}

void CdrTGA::setupMultigrid(){
  CH_TIME("CdrTGA::setupMultigrid");
  if(m_verbosity > 5){
    pout() << m_name + "::setupMultigrid" << endl;
  }

  if(!m_hasDeeperMultigridLevels){
    this->defineDeeperMultigridLevels();
    m_hasDeeperMultigridLevels = true;
  }
  this->setupOperatorFactory();
  this->setupMultigridSolver();

  this->setupTGA();
  this->setupEuler();
}

void CdrTGA::defineDeeperMultigridLevels(){
  CH_TIME("CdrTGA::defineDeeperMultigridLevels");
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
  const int numEbGhostsCells                = m_amr->getNumberOfEbGhostCells();

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
					    numEbGhostsCells,
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

void CdrTGA::setupOperatorFactory(){
  CH_TIME("CdrTGA::setupOperatorFactory");
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
    levelgrids.push_back(*(m_amr->getEBLevelGrid(m_realm, m_phase)[lvl])); // EBConductivityOp does not not refcounted operators
  }

  // Appropriate coefficients. These don't matter right now. 
  const Real alpha =  1.0;
  const Real beta  =  1.0;

  // Default is Neumann. This might change in the future. 
  RefCountedPtr<NeumannConductivityDomainBCFactory> domfact = RefCountedPtr<NeumannConductivityDomainBCFactory>
    (new NeumannConductivityDomainBCFactory());
  RefCountedPtr<NeumannConductivityEBBCFactory> ebfact  = RefCountedPtr<NeumannConductivityEBBCFactory>
    (new NeumannConductivityEBBCFactory());

  domfact->setValue(0.0);
  ebfact->setValue(0.0);

  // Set the relaxation type
  int relax_type = 0;
  if(m_multigridRelaxMethod == RelaxationMethod::Jacobi){
    relax_type = 0;
  }
  else if(m_multigridRelaxMethod == RelaxationMethod::GaussSeidel){
    relax_type = 1;
  }
  else if(m_multigridRelaxMethod == RelaxationMethod::GSRB){
    relax_type = 2;
  }

      
  
  // Create operator factory.
  DataOps::setValue(m_aCoef, 1.0); // We're usually solving (1 - dt*nabla^2)*phi^(k+1) = phi^k + dt*S^k so aco=1
  DataOps::setValue(m_ebCenteredDiffusionCoefficient, 0.0);
  m_operatorFactory = RefCountedPtr<EbHelmholtzOpFactory> (new EbHelmholtzOpFactory(levelgrids,
										 quadcfi,
										 fastFR,
										 alpha,
										 beta,
										 m_aCoef.get_data(),
										 m_faceCenteredDiffusionCoefficient.get_data(),
										 m_ebCenteredDiffusionCoefficient.get_data(),
										 dx[0],
										 refinement_ratios,
										 domfact,
										 ebfact,
										 ghost*IntVect::Unit,
										 ghost*IntVect::Unit,
										 relax_type,
										 m_numCellsBottomDrop,
										 -1,
										 m_mg_levelgrids));
}

void CdrTGA::setupMultigridSolver(){
  CH_TIME("CdrTGA::setupMultigrid");
  if(m_verbosity > 5){
    pout() << m_name + "::setupMultigrid" << endl;
  }

  const int finest_level       = m_amr->getFinestLevel();
  const ProblemDomain coar_dom = m_amr->getDomains()[0];

  // Select bottom solver
  LinearSolver<LevelData<EBCellFAB> >* botsolver = NULL;
  if(m_bottomSolver == BottomSolver::Simple){
    m_simpleSolver.setNumSmooths(m_numSmoothingsForSimpleSolver);
    botsolver = &m_simpleSolver;
  }
  else if (m_bottomSolver == BottomSolver::BiCGStab){
    botsolver = &m_bicgstab;
  }
  m_multigridSolver = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > (new AMRMultiGrid<LevelData<EBCellFAB> >());
  m_multigridSolver->m_imin = m_multigridMinIterations;
  m_multigridSolver->m_verbosity = m_multigridVerbosity;
  m_multigridSolver->define(coar_dom, *m_operatorFactory, botsolver, 1 + finest_level);

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
  
  m_multigridSolver->setSolverParameters(m_multigridPreSmooth,
				    m_multigridPostSmooth,
				    m_multigridBottomSmooth,
				    gmg_type,
				    m_multigridMaxIterations,
				    m_multigridTolerance,
				    m_multigridHang,
				    1.E-90); // Residue set through other means

}

void CdrTGA::setupTGA(){
  CH_TIME("CdrTGA::setupTGA");
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

void CdrTGA::setupEuler(){
  CH_TIME("CdrTGA::setupEuler");
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

    // Compute advection flux. This is mostly the same as computeDivF
    if(m_isMobile){
      m_amr->interpGhostPwl(m_cellVelocity, m_realm, m_phase);
      
      this->averageVelocityToFaces(); // Update m_faceVelocity from m_cellVelocity
      this->advectToFaces(m_faceStates, a_phi, a_extrapDt); // Advect to faces
      this->computeFlux(m_scratchFluxTwo, m_faceVelocity, m_faceStates, m_domainFlux);

      DataOps::incr(m_scratchFluxOne, m_scratchFluxTwo, 1.0);
    }

    // Compute diffusion flux. 
    if(m_isDiffusive){
      this->computeDiffusionFlux(m_scratchFluxTwo, a_phi);
      DataOps::incr(m_scratchFluxOne, m_scratchFluxTwo, -1.0);
    }

    // General divergence computation. Also inject charge. Domain fluxes came in through the compute
    // advective flux function and eb fluxes come in through the divergence computation.
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

void CdrTGA::computeDivF(EBAMRCellData& a_divF, EBAMRCellData& a_phi, const Real a_extrapDt, const bool a_ebFlux){
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
    this->averageVelocityToFaces();
    this->advectToFaces(m_faceStates, a_phi, a_extrapDt);          // Face extrapolation to cell-centered faces
    this->computeFlux(m_scratchFluxOne, m_faceVelocity, m_faceStates, m_domainFlux);  // Compute face-centered fluxes

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

void CdrTGA::computeDivD(EBAMRCellData& a_divD, EBAMRCellData& a_phi, const bool a_ebFlux){
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

    this->computeDiffusionFlux(m_scratchFluxOne, a_phi);  // Compute the face-centered diffusion flux

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
  
  pp.get("gmg_coarsen",     m_numCoarseningsBeforeAggregation);
  pp.get("gmg_verbosity",   m_multigridVerbosity);
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
    MayDay::Abort("CdrTGA::parseMultigridSettings - unknown bottom solver requested");
  }

  // Relaxation type
  pp.get("gmg_relax_type", str);
  if(str == "gsrb"){
    m_multigridRelaxMethod = RelaxationMethod::GSRB;
  }
  else if(str == "jacobi"){
    m_multigridRelaxMethod = RelaxationMethod::Jacobi;
  }
  else if(str == "GaussSeidel"){
    m_multigridRelaxMethod = RelaxationMethod::GaussSeidel;
  }
  else{
    MayDay::Abort("CdrTGA::parseMultigridSettings - unknown relaxation method requested");
  }

  // Cycle type
  pp.get("gmg_cycle", str);
  if(str == "vcycle"){
    m_multigridType = MultigridType::VCycle;
  }
  else{
    MayDay::Abort("CdrTGA::parseMultigridSettings - unknown cycle type requested");
  }

  // No lower than 2. 
  if(m_numCellsBottomDrop < 2){
    m_numCellsBottomDrop = 2;
  }
}

void CdrTGA::writePlotData(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("CdrTGA::writePlotData");
  if(m_verbosity > 5){
    pout() << m_name + "::writePlotData" << endl;
  }

  if(m_plotPhi) {
    if(!m_plotNumbers){ // Regular write
      writeData(a_output, a_comp, m_phi, true);
    }
    else{ // Scale, write, and scale back
     
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
	DataOps::scale(*m_phi[lvl], (pow(m_amr->getDx()[lvl], 3)));
      }
      writeData(a_output, a_comp, m_phi,    false);
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
	DataOps::scale(*m_phi[lvl], 1./(pow(m_amr->getDx()[lvl], 3)));
      }
    }
  }

  if(m_plotDiffusionCoefficient && m_isDiffusive) { // Need to compute the cell-centerd stuff first
    DataOps::setValue(m_scratch, 0.0);
    DataOps::averageFaceToCell(m_scratch, m_faceCenteredDiffusionCoefficient, m_amr->getDomains());
    writeData(a_output, a_comp, m_scratch,   false);
  }

  if(m_plotSource) {
    if(!m_plotNumbers){
      writeData(a_output, a_comp, m_source,    false);
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
    writeData(a_output, a_comp, m_cellVelocity, false);
  }

  // Plot EB fluxes
  if(m_plotEbFlux && m_isMobile){
    DataOps::setValue(m_scratch, 0.0);
    DataOps::incr(m_scratch, m_ebFlux, 1.0);
    writeData(a_output, a_comp, m_scratch, false);
  }
}

#include <CD_NamespaceFooter.H>
