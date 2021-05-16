/*!
  @file   cdr_tga.cpp
  @brief  Implementation of cdr_tga.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "cdr_tga.H"
#include "data_ops.H"

#include <ParmParse.H>
#include <EBArith.H>
#include <EBAMRIO.H>
#include <NeumannConductivityEBBC.H>
#include <NeumannConductivityDomainBC.H>
#include <BRMeshRefine.H>

#include "CD_NamespaceHeader.H"

cdr_tga::cdr_tga() : CdrSolver() {
  m_name         = "cdr_tga";
  m_className   = "cdr_tga";
  m_has_mg_stuff = false;
}

cdr_tga::~cdr_tga(){

}

void cdr_tga::advanceEuler(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const Real a_dt){
  CH_TIME("cdr_tga::advanceEuler(no source)");
  if(m_verbosity > 5){
    pout() << m_name + "::advanceEuler(no source)" << endl;
  }
  
  if(m_isDiffusive){
    // Create a source term = S = 0.0;
    const int ncomp = 1;
    EBAMRCellData src;
    m_amr->allocate(src, m_Realm, m_phase, ncomp);
    data_ops::set_value(src, 0.0);

    // Call version with source term
    advanceEuler(a_newPhi, a_oldPhi, src, a_dt);
  }
  else{
    data_ops::copy(a_newPhi, a_oldPhi);
  }
}

void cdr_tga::advanceEuler(EBAMRCellData&       a_newPhi,
			    const EBAMRCellData& a_oldPhi,
			    const EBAMRCellData& a_source,
			    const Real           a_dt){
  CH_TIME("cdr_tga::advanceEuler");
  if(m_verbosity > 5){
    pout() << m_name + "::advanceEuler" << endl;
  }
  
  if(m_isDiffusive){
    this->setup_gmg(); // Set up gmg again since diffusion coefficients might change 
    
    bool converged = false;

    const int comp         = 0;
    const int ncomp        = 1;
    const int finest_level = m_amr->getFinestLevel();

    // Set convergence metric. Make zero = 0 and m_scratch = phi^k + dt*S^k. Then
    // compute the residue L(zero) - m_scratch which sets the convergence metric. 
    EBAMRCellData zero;
    m_amr->allocate(zero, m_Realm, m_phase, ncomp);
    data_ops::set_value(zero, 0.0);
    data_ops::copy(m_scratch, a_oldPhi);
    data_ops::incr(m_scratch, a_source, a_dt);

    Vector<LevelData<EBCellFAB>* > orez;
    Vector<LevelData<EBCellFAB>* > shr;
    m_amr->alias(orez, zero);
    m_amr->alias(shr,  m_scratch);

    m_eulersolver->resetAlphaAndBeta(1.0, -a_dt);
    
    const Real zero_resid = m_gmg_solver->computeAMRResidual(orez, shr, finest_level, 0);
    const Real stopcrit   = zero_resid;
    m_gmg_solver->m_convergenceMetric = zero_resid;


    // Now do the solve. 
    Vector<LevelData<EBCellFAB>* > new_state, old_state, source;
    m_amr->alias(new_state, a_newPhi);
    m_amr->alias(old_state, a_oldPhi);
    m_amr->alias(source,    a_source);


    // Euler solve
    data_ops::set_value(m_faceCenteredDiffusionCoefficient_eb, 0.0);
    m_eulersolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, false);
    const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }
  }
  else{
    data_ops::copy(a_newPhi, a_oldPhi);
  }
}

void cdr_tga::advanceTGA(EBAMRCellData& a_newPhi, const EBAMRCellData& a_oldPhi, const Real a_dt){
  CH_TIME("cdr_tga::advanceTGA(no source)");
  if(m_verbosity > 5){
    pout() << m_name + "::advanceTGA(no source)" << endl;
  }

  if(m_isDiffusive){

    // Dummy source term
    const int ncomp = 1;
    EBAMRCellData src;
    m_amr->allocate(src, m_Realm, m_phase, ncomp);
    data_ops::set_value(src, 0.0);

    // Call other version
    advanceTGA(a_newPhi, a_oldPhi, src, a_dt);
  }
  else{
    data_ops::copy(a_newPhi, a_oldPhi);
  }
}

void cdr_tga::advanceTGA(EBAMRCellData&       a_newPhi,
			  const EBAMRCellData& a_oldPhi,
			  const EBAMRCellData& a_source,
			  const Real           a_dt){
  CH_TIME("cdr_tga::advanceTGA(full)");
  if(m_verbosity > 5){
    pout() << m_name + "::advanceTGA(full)" << endl;
  }
  
  if(m_isDiffusive){
    this->setup_gmg(); // Set up gmg again since diffusion coefficients might change
    
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

    data_ops::set_value(m_faceCenteredDiffusionCoefficient_eb, 0.0);

    // TGA solve
    m_tgasolver->resetAlphaAndBeta(alpha, beta);
    m_tgasolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, false);

    const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }
  }
  else{
    data_ops::copy(a_newPhi, a_oldPhi);
  }
}

void cdr_tga::set_bottom_solver(const int a_whichsolver){
  CH_TIME("cdr_tga::set_bottom_solver");
  if(m_verbosity > 5){
    pout() << m_name + "::set_bottom_solver" << endl;
  }

  if(a_whichsolver == 0 || a_whichsolver == 1){
    m_bottomsolver = a_whichsolver;

    std::string str;
    ParmParse pp("cdr_tga");
    pp.get("gmg_bottom_solver", str);
    if(str == "simple"){
      m_bottomsolver = 0;
    }
    else if(str == "bicgstab"){
      m_bottomsolver = 1;
    }
  }
  else{
    MayDay::Abort("cdr_tga::set_bottom_solver - Unsupported solver type requested");
  }
}

void cdr_tga::set_botsolver_smooth(const int a_numsmooth){
  CH_TIME("cdr_tga::set_botsolver_smooth");
  if(m_verbosity > 5){
    pout() << m_name + "::set_botsolver_smooth" << endl;
  }
  CH_assert(a_numsmooth > 0);


  ParmParse pp("cdr_tga");
  pp.query("gmg_bottom_relax", m_numsmooth);

  if(m_numsmooth < 0){
    m_numsmooth = 0;
  }
}

void cdr_tga::set_bottom_drop(const int a_bottom_drop){
  CH_TIME("cdr_tga::set_bottom_drop");
  if(m_verbosity > 5){
    pout() << m_name + "::set_bottom_drop" << endl;
  }

  ParmParse pp("cdr_tga");
  pp.get("gmg_bottom_drop", m_bottom_drop);

  if(m_bottom_drop < 2){
    m_bottom_drop = 2;
  }
}

void cdr_tga::set_tga(const bool a_use_tga){
  CH_TIME("cdr_tga::set_tga");
  if(m_verbosity > 5){
    pout() << m_name + "::set_tga" << endl;
  }
  
  m_use_tga = a_use_tga;

  ParmParse pp("cdr_tga");
  std::string str;
  if(pp.contains("use_tga")){
    pp.get("use_tga", str);
    if(str == "true"){
      m_use_tga = true;
    }
    else if(str == "false"){
      m_use_tga = false;
    }
  }
}

void cdr_tga::set_gmg_solver_parameters(relax      a_relax_type,
					amrmg      a_gmg_type,      
					const int  a_verbosity,          
					const int  a_pre_smooth,         
					const int  a_post_smooth,       
					const int  a_bot_smooth,         
					const int  a_max_iter,
					const int  a_min_iter,
					const Real a_eps,               
					const Real a_hang){
  CH_TIME("cdr_tga::set_gmg_solver_parameters");
  if(m_verbosity > 5){
    pout() << m_name + "::set_gmg_solver_parameters" << endl;
  }

  m_gmg_relax_type  = a_relax_type;
  m_gmg_type        = a_gmg_type;
  m_gmg_verbosity   = a_verbosity;
  m_gmg_pre_smooth  = a_pre_smooth;
  m_gmg_post_smooth = a_post_smooth;
  m_gmg_bot_smooth  = a_bot_smooth;
  m_gmg_max_iter    = a_max_iter;
  m_gmg_min_iter    = a_min_iter;
  m_gmg_eps         = a_eps;
  m_gmg_hang        = a_hang;

  ParmParse pp("cdr_tga");
  std::string str;


  pp.get("gmg_verbosity",   m_gmg_verbosity);
  pp.get("gmg_pre_smooth",  m_gmg_pre_smooth);
  pp.get("gmg_post_smooth", m_gmg_post_smooth);
  pp.get("gmg_bott_smooth", m_gmg_bot_smooth);
  pp.get("gmg_max_iter",    m_gmg_max_iter);
  pp.get("gmg_min_iter",    m_gmg_min_iter);
  pp.get("gmg_tolerance",   m_gmg_eps);
  pp.get("gmg_hang",        m_gmg_hang);

  // Get the relaxation type
  pp.get("gmg_relax",        str);
  if(str == "jacobi"){
    m_gmg_relax_type = relax::jacobi;
  }
  else if(str == "gauss_seidel"){
    m_gmg_relax_type = relax::gauss_seidel;
  }
  else if(str == "gsrb"){
    m_gmg_relax_type = relax::gsrb_fast;
  }
  else{
    MayDay::Abort("cdr_tga::set_gmg_solver_parameters - unknown relaxation method requested");
  }

  // Get the MG cycle
  pp.get("gmg_cycle",        str);
  if(str == "full"){
    m_gmg_type = amrmg::full;
  }
  else if(str == "vvcycle"){
    m_gmg_type = amrmg::vcycle;
  }
  else if(str == "fcycle"){
    m_gmg_type = amrmg::fcycle;
  }
  else{
    MayDay::Abort("cdr_tga::set_gmg_solver_parameters - unknown MG cycle method requested");
  }
}

void cdr_tga::setup_gmg(){
  CH_TIME("cdr_tga::setup_gmg");
  if(m_verbosity > 5){
    pout() << m_name + "::setup_gmg" << endl;
  }

  if(!m_has_mg_stuff){
    this->define_mg_levels();
    m_has_mg_stuff = true;
  }
  this->setup_operator_factory();
  this->setup_multigrid();

  this->setup_tga();
  this->setup_euler();
}

void cdr_tga::define_mg_levels(){
  CH_TIME("cdr_tga::define_mg_levels");
  if(m_verbosity > 5){
    pout() << m_name + "::define_mg_levels" << endl;
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
  const int num_numEbGhostsCells                = m_amr->getNumberOfEbGhostCells();

  int num_coar       = 0;
  bool has_coar      = true;
  ProblemDomain fine = domains[0];

  // Coarsen problem domains and create grids
  while(num_coar < m_gmg_coarsen || !has_coar){

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

void cdr_tga::setup_operator_factory(){
  CH_TIME("cdr_tga::setup_operator_factory");
  if(m_verbosity > 5){
    pout() << m_name + "::setup_operator_factory" << endl;
  }

  const int finest_level                 = m_amr->getFinestLevel();
  const int ghost                        = m_amr->getNumberOfGhostCells();
  const Vector<DisjointBoxLayout>& grids = m_amr->getGrids(m_Realm);
  const Vector<int>& refinement_ratios   = m_amr->getRefinementRatios();
  const Vector<ProblemDomain>& domains   = m_amr->getDomains();
  const Vector<Real>& dx                 = m_amr->getDx();
  const RealVect& origin                 = m_amr->getProbLo();
  const Vector<EBISLayout>& ebisl        = m_amr->getEBISLayout(m_Realm, m_phase);
  
  const Vector<RefCountedPtr<EBQuadCFInterp> >& quadcfi  = m_amr->getEBQuadCFInterp(m_Realm, m_phase);
  const Vector<RefCountedPtr<EBFluxRegister> >& fastFR   = m_amr->getFluxRegister(m_Realm, m_phase);

  Vector<EBLevelGrid> levelgrids;
  for (int lvl = 0; lvl <= finest_level; lvl++){ 
    levelgrids.push_back(*(m_amr->getEBLevelGrid(m_Realm, m_phase)[lvl])); // EBConductivityOp does not not refcounted operators
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
  if(m_gmg_relax_type == relax::jacobi){
    relax_type = 0;
  }
  else if(m_gmg_relax_type == relax::gauss_seidel){
    relax_type = 1;
  }
  else if(m_gmg_relax_type == relax::gsrb_fast){
    relax_type = 2;
  }

      
  
  // Create operator factory.
  data_ops::set_value(m_aCoefficient, 1.0); // We're usually solving (1 - dt*nabla^2)*phi^(k+1) = phi^k + dt*S^k so aco=1
  data_ops::set_value(m_faceCenteredDiffusionCoefficient_eb, 0.0);
  m_opfact = RefCountedPtr<ebconductivityopfactory> (new ebconductivityopfactory(levelgrids,
										 quadcfi,
										 fastFR,
										 alpha,
										 beta,
										 m_aCoefficient.get_data(),
										 m_faceCenteredDiffusionCoefficient.get_data(),
										 m_faceCenteredDiffusionCoefficient_eb.get_data(),
										 dx[0],
										 refinement_ratios,
										 domfact,
										 ebfact,
										 ghost*IntVect::Unit,
										 ghost*IntVect::Unit,
										 relax_type,
										 m_bottom_drop,
										 -1,
										 m_mg_levelgrids));
}

void cdr_tga::setup_multigrid(){
  CH_TIME("cdr_tga::setup_multigrid");
  if(m_verbosity > 5){
    pout() << m_name + "::setup_multigrid" << endl;
  }

  const int finest_level       = m_amr->getFinestLevel();
  const ProblemDomain coar_dom = m_amr->getDomains()[0];

  // Select bottom solver
  LinearSolver<LevelData<EBCellFAB> >* botsolver = NULL;
  if(m_bottomsolver == 0){
    m_simple_solver.setNumSmooths(m_numsmooth);
    botsolver = &m_simple_solver;
  }
  else{
    botsolver = &m_bicgstab;
  }
  m_gmg_solver = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > (new AMRMultiGrid<LevelData<EBCellFAB> >());
  m_gmg_solver->m_imin = m_gmg_min_iter;
  m_gmg_solver->m_verbosity = m_gmg_verbosity;
  m_gmg_solver->define(coar_dom, *m_opfact, botsolver, 1 + finest_level);

  // Make m_gmg_type into an int for multigrid
  int gmg_type;
  if(m_gmg_type == amrmg::full){
    gmg_type = 0;
  }
  else if(m_gmg_type == amrmg::vcycle){
    gmg_type = 1;
  }
  else if(m_gmg_type == amrmg::fcycle){
    gmg_type = 2;
  }
  
  m_gmg_solver->setSolverParameters(m_gmg_pre_smooth,
				    m_gmg_post_smooth,
				    m_gmg_bot_smooth,
				    gmg_type,
				    m_gmg_max_iter,
				    m_gmg_eps,
				    m_gmg_hang,
				    1.E-90); // Residue set through other means

}

void cdr_tga::setup_tga(){
  CH_TIME("cdr_tga::setup_tga");
  if(m_verbosity > 5){
    pout() << m_name + "::setup_tga" << endl;
  }
  
  const int finest_level       = m_amr->getFinestLevel();
  const ProblemDomain coar_dom = m_amr->getDomains()[0];
  const Vector<int> ref_rat    = m_amr->getRefinementRatios();

  m_tgasolver = RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > >
    (new AMRTGA<LevelData<EBCellFAB> > (m_gmg_solver, *m_opfact, coar_dom, ref_rat, 1 + finest_level, m_gmg_solver->m_verbosity));

  // Must init gmg for TGA
  Vector<LevelData<EBCellFAB>* > phi, rhs;
  m_amr->alias(phi, m_phi);
  m_amr->alias(rhs, m_source);
  m_gmg_solver->init(phi, rhs, finest_level, 0);
}

void cdr_tga::setup_euler(){
  CH_TIME("cdr_tga::setup_euler");
  if(m_verbosity > 5){
    pout() << m_name + "::setup_euler" << endl;
  }

  const int finest_level       = m_amr->getFinestLevel();
  const ProblemDomain coar_dom = m_amr->getDomains()[0];
  const Vector<int> ref_rat    = m_amr->getRefinementRatios();

  m_eulersolver = RefCountedPtr<EBBackwardEuler> 
    (new EBBackwardEuler (m_gmg_solver, *m_opfact, coar_dom, ref_rat, 1 + finest_level, m_gmg_solver->m_verbosity));

  // Note: If this crashes, try to init gmg first
}

void cdr_tga::computeDivJ(EBAMRCellData& a_divJ, EBAMRCellData& a_phi, const Real a_extrapDt, const bool a_ebFlux){
  CH_TIME("cdr_tga::computeDivJ(divF, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDivJ(divF, state)" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;

  // Fill ghost cells
  m_amr->interpGhostPwl(a_phi, m_Realm, m_phase);

  if(m_isMobile || m_isDiffusive){

    if(m_useMassWeightedRedistribution){
      this->setRedistWeights(a_phi);
    }

    // We will let m_scratchFluxOne hold the total flux = advection + diffusion fluxes
    data_ops::set_value(m_scratchFluxOne, 0.0);

    // Compute advection flux. This is mostly the same as computeDivF
    if(m_isMobile){
      m_amr->interpGhostPwl(m_cellVelocity, m_Realm, m_phase);
      
      this->averageVelocityToFaces(); // Update m_faceVelocity from m_cellVelocity
      this->advect_to_faces(m_faceStates, a_phi, a_extrapDt); // Advect to faces
      this->computeFlux(m_scratchFluxTwo, m_faceVelocity, m_faceStates, m_domainFlux);

      data_ops::incr(m_scratchFluxOne, m_scratchFluxTwo, 1.0);
    }

    // Compute diffusion flux. 
    if(m_isDiffusive){
      this->computeDiffusionFlux(m_scratchFluxTwo, a_phi);
      data_ops::incr(m_scratchFluxOne, m_scratchFluxTwo, -1.0);
    }

    // General divergence computation. Also inject charge. Domain fluxes came in through the compute
    // advective flux function and eb fluxes come in through the divergence computation.
    EBAMRIVData* ebflux;
    if(a_ebFlux){
      ebflux = &m_ebFlux;
    }
    else{
      ebflux = &m_EbZero;
    }
    this->computeDivG(a_divJ, m_scratchFluxOne, *ebflux);
  }
  else{ 
    data_ops::set_value(a_divJ, 0.0);
  }

  return;
}

void cdr_tga::computeDivF(EBAMRCellData& a_divF, EBAMRCellData& a_phi, const Real a_extrapDt, const bool a_ebFlux){
  CH_TIME("cdr_tga::computeDivF(divF, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDivF(divF, state)" << endl;
  }

  if(m_isMobile){

    // Fill ghost cells
    m_amr->interpGhostPwl(a_phi,     m_Realm, m_phase);
    m_amr->interpGhostPwl(m_cellVelocity, m_Realm, m_phase);

    if(m_useMassWeightedRedistribution){
      this->setRedistWeights(a_phi);
    }
    this->averageVelocityToFaces();
    this->advect_to_faces(m_faceStates, a_phi, a_extrapDt);          // Face extrapolation to cell-centered faces
    this->computeFlux(m_scratchFluxOne, m_faceVelocity, m_faceStates, m_domainFlux);  // Compute face-centered fluxes

    EBAMRIVData* ebflux;
    if(a_ebFlux){
      ebflux = &m_ebFlux;
    }
    else{
      ebflux = &m_EbZero;
    }
    this->computeDivG(a_divF, m_scratchFluxOne, *ebflux); 
  }
  else{
    data_ops::set_value(a_divF, 0.0);
  }

#if 0 // Debug
  Real volume;
  const Real sum = EBLevelDataOps::noKappaSumLevel(volume, *a_divF[0], 0, m_amr->getDomains()[0]);

  std::cout << sum << std::endl;
#endif

}

void cdr_tga::computeDivD(EBAMRCellData& a_divD, EBAMRCellData& a_phi, const bool a_ebFlux){
  CH_TIME("cdr_tga::computeDivD");
  if(m_verbosity > 5){
    pout() << m_name + "::computeDivD" << endl;
  }

  if(m_isDiffusive){
    const int ncomp = 1;

    // Fill ghost cells
    m_amr->interpGhostPwl(a_phi, m_Realm, m_phase);

    if(m_useMassWeightedRedistribution){
      this->setRedistWeights(a_phi);
    }

    this->computeDiffusionFlux(m_scratchFluxOne, a_phi);  // Compute the face-centered diffusion flux

    EBAMRIVData* ebflux;
    if(a_ebFlux){
      ebflux = &m_ebFlux;
    }
    else{
      ebflux = &m_EbZero;
    }
    this->computeDivG(a_divD, m_scratchFluxOne, *ebflux); // General face-centered flux to divergence magic.

    m_amr->averageDown(a_divD, m_Realm, m_phase);
    m_amr->interpGhost(a_divD, m_Realm, m_phase);
  }
  else{
    data_ops::set_value(a_divD, 0.0);
  }
}

void cdr_tga::parse_gmg_settings(){
  ParmParse pp(m_className.c_str());

  std::string str;
  
  pp.get("gmg_coarsen",     m_gmg_coarsen);
  pp.get("gmg_verbosity",   m_gmg_verbosity);
  pp.get("gmg_pre_smooth",  m_gmg_pre_smooth);
  pp.get("gmg_post_smooth", m_gmg_post_smooth);
  pp.get("gmg_bott_smooth", m_gmg_bot_smooth);
  pp.get("gmg_max_iter",    m_gmg_max_iter);
  pp.get("gmg_min_iter",    m_gmg_min_iter);
  pp.get("gmg_tolerance",   m_gmg_eps);
  pp.get("gmg_hang",        m_gmg_hang);
  pp.get("gmg_bottom_drop", m_bottom_drop);

  // Bottom solver
  pp.get("gmg_bottom_solver", str);
  if(str == "simple"){
    m_bottomsolver = 0;
  }
  else if(str == "bicgstab"){
    m_bottomsolver = 1;
  }
  else{
    MayDay::Abort("cdr_tga::parse_gmg_settings - unknown bottom solver requested");
  }

  // Relaxation type
  pp.get("gmg_relax_type", str);
  if(str == "gsrb"){
    m_gmg_relax_type = relax::gsrb_fast;
  }
  else if(str == "jacobi"){
    m_gmg_relax_type = relax::jacobi;
  }
  else if(str == "gauss_seidel"){
    m_gmg_relax_type = relax::gauss_seidel;
  }
  else{
    MayDay::Abort("cdr_tga::parse_gmg_settings - unknown relaxation method requested");
  }

  // Cycle type
  pp.get("gmg_cycle", str);
  if(str == "vcycle"){
    m_gmg_type = amrmg::vcycle;
  }
  else{
    MayDay::Abort("cdr_tga::parse_gmg_settings - unknown cycle type requested");
  }

  // No lower than 2. 
  if(m_bottom_drop < 2){
    m_bottom_drop = 2;
  }
}

void cdr_tga::writePlotData(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("cdr_tga::writePlotData");
  if(m_verbosity > 5){
    pout() << m_name + "::writePlotData" << endl;
  }

  if(m_plotPhi) {
    if(!m_plotNumbers){ // Regular write
      writeData(a_output, a_comp, m_phi, true);
    }
    else{ // Scale, write, and scale back
     
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
	data_ops::scale(*m_phi[lvl], (pow(m_amr->getDx()[lvl], 3)));
      }
      writeData(a_output, a_comp, m_phi,    false);
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
	data_ops::scale(*m_phi[lvl], 1./(pow(m_amr->getDx()[lvl], 3)));
      }
    }
  }

  if(m_plotDiffusionCoefficient && m_isDiffusive) { // Need to compute the cell-centerd stuff first
    data_ops::set_value(m_scratch, 0.0);
    data_ops::average_face_to_cell(m_scratch, m_faceCenteredDiffusionCoefficient, m_amr->getDomains());
    writeData(a_output, a_comp, m_scratch,   false);
  }

  if(m_plotSource) {
    if(!m_plotNumbers){
      writeData(a_output, a_comp, m_source,    false);
    }
    else { // Scale, write, and scale back
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
	data_ops::scale(*m_source[lvl], (pow(m_amr->getDx()[lvl], 3)));
      }
      writeData(a_output, a_comp, m_source,    false);
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
	data_ops::scale(*m_source[lvl], 1./(pow(m_amr->getDx()[lvl], 3)));
      }
    }
  }

  if(m_plotVelocity && m_isMobile) {
    writeData(a_output, a_comp, m_cellVelocity, false);
  }

  // Plot EB fluxes
  if(m_plotEbFlux && m_isMobile){
    data_ops::set_value(m_scratch, 0.0);
    data_ops::incr(m_scratch, m_ebFlux, 1.0);
    writeData(a_output, a_comp, m_scratch, false);
  }
}
#include "CD_NamespaceFooter.H"
