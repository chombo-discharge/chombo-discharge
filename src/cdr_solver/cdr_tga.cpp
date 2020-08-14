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

cdr_tga::cdr_tga() : cdr_solver() {
  m_name         = "cdr_tga";
  m_class_name   = "cdr_tga";
  m_has_mg_stuff = false;
}

cdr_tga::~cdr_tga(){

}

void cdr_tga::advance_euler(EBAMRCellData& a_new_state, const EBAMRCellData& a_old_state, const Real a_dt){
  CH_TIME("cdr_tga::advance_euler(no source)");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_euler(no source)" << endl;
  }
  
  if(m_diffusive){
    // Create a source term = S = 0.0;
    const int ncomp = 1;
    EBAMRCellData src;
    m_amr->allocate(src, m_realm, m_phase, ncomp);
    data_ops::set_value(src, 0.0);

    // Call version with source term
    advance_euler(a_new_state, a_old_state, src, a_dt);
  }
  else{
    data_ops::copy(a_new_state, a_old_state);
  }
}

void cdr_tga::advance_euler(EBAMRCellData&       a_new_state,
			    const EBAMRCellData& a_old_state,
			    const EBAMRCellData& a_source,
			    const Real           a_dt){
  CH_TIME("cdr_tga::advance_euler");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_euler" << endl;
  }
  
  if(m_diffusive){
    this->setup_gmg(); // Set up gmg again since diffusion coefficients might change 
    
    bool converged = false;

    const int comp         = 0;
    const int ncomp        = 1;
    const int finest_level = m_amr->get_finest_level();

    // Do the aliasing stuff
    Vector<LevelData<EBCellFAB>* > new_state, old_state, source;
    m_amr->alias(new_state, a_new_state);
    m_amr->alias(old_state, a_old_state);
    m_amr->alias(source,    a_source);

    const Real alpha = 0.0;
    const Real beta  = 1.0;

    m_eulersolver->resetAlphaAndBeta(alpha, beta);
    data_ops::set_value(m_diffco_eb, 0.0);
    
    // Euler solve
    m_eulersolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, false);
    const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }
  }
  else{
    data_ops::copy(a_new_state, a_old_state);
  }
}

void cdr_tga::advance_tga(EBAMRCellData& a_new_state, const EBAMRCellData& a_old_state, const Real a_dt){
  CH_TIME("cdr_tga::advance_tga(no source)");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_tga(no source)" << endl;
  }

  if(m_diffusive){

    // Dummy source term
    const int ncomp = 1;
    EBAMRCellData src;
    m_amr->allocate(src, m_realm, m_phase, ncomp);
    data_ops::set_value(src, 0.0);

    // Call other version
    advance_tga(a_new_state, a_old_state, src, a_dt);
  }
  else{
    data_ops::copy(a_new_state, a_old_state);
  }
}

void cdr_tga::advance_tga(EBAMRCellData&       a_new_state,
			  const EBAMRCellData& a_old_state,
			  const EBAMRCellData& a_source,
			  const Real           a_dt){
  CH_TIME("cdr_tga::advance_tga(full)");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_tga(full)" << endl;
  }
  
  if(m_diffusive){
    this->setup_gmg(); // Set up gmg again since diffusion coefficients might change
    
    bool converged = false;

    const int comp         = 0;
    const int ncomp        = 1;
    const int finest_level = m_amr->get_finest_level();

    // Do the aliasing stuff
    Vector<LevelData<EBCellFAB>* > new_state, old_state, source;
    m_amr->alias(new_state, a_new_state);
    m_amr->alias(old_state, a_old_state);
    m_amr->alias(source,    a_source);

    const Real alpha = 0.0;
    const Real beta  = 1.0;

    data_ops::set_value(m_diffco_eb, 0.0);

    // TGA solve
    m_tgasolver->resetAlphaAndBeta(alpha, beta);
    m_tgasolver->oneStep(new_state, old_state, source, a_dt, 0, finest_level, false);

    const int status = m_gmg_solver->m_exitStatus;  // 1 => Initial norm sufficiently reduced
    if(status == 1 || status == 8 || status == 9){  // 8 => Norm sufficiently small
      converged = true;
    }
  }
  else{
    data_ops::copy(a_new_state, a_old_state);
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

  // Get some stuff from amr_mesh on how to decompose the levels
  const Vector<ProblemDomain>& domains = m_amr->get_domains();
  const int max_box_size               = m_amr->get_max_box_size();
  const int blocking_factor            = m_amr->get_blocking_factor();
  const int num_ebghost                = m_amr->get_eb_ghost();

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
      load_balance::balance_volume(proc_assign, boxes);

      // Add problem domain and grid
      m_mg_domains.push_back(coar);
      m_mg_grids.push_back(DisjointBoxLayout(boxes, proc_assign, coar));

      // Define the EBLevelGrids
      const int idx = m_mg_grids.size() - 1; // Last element added
      m_mg_levelgrids.push_back(EBLevelGrid(m_mg_grids[idx],
					    m_mg_domains[idx],
					    num_ebghost,
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

  const int finest_level                 = m_amr->get_finest_level();
  const int ghost                        = m_amr->get_num_ghost();
  const Vector<DisjointBoxLayout>& grids = m_amr->get_grids(m_realm);
  const Vector<int>& refinement_ratios   = m_amr->get_ref_rat();
  const Vector<ProblemDomain>& domains   = m_amr->get_domains();
  const Vector<Real>& dx                 = m_amr->get_dx();
  const RealVect& origin                 = m_amr->get_prob_lo();
  const Vector<EBISLayout>& ebisl        = m_amr->get_ebisl(m_realm, m_phase);
  
  const Vector<RefCountedPtr<EBQuadCFInterp> >& quadcfi  = m_amr->get_old_quadcfi(m_realm, m_phase);
  const Vector<RefCountedPtr<EBFasterFR> >& fastFR       = m_amr->get_eb_fast_fr(m_realm, m_phase);

  Vector<EBLevelGrid> levelgrids;
  for (int lvl = 0; lvl <= finest_level; lvl++){ 
    levelgrids.push_back(*(m_amr->get_eblg(m_realm, m_phase)[lvl])); // EBConductivityOp does not not refcounted operators
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
  data_ops::set_value(m_aco, 1.0); // We're usually solving (1 - dt*nabla^2)*phi^(k+1) = phi^k + dt*S^k so aco=1
  data_ops::set_value(m_diffco_eb, 0.0);
  m_opfact = RefCountedPtr<ebconductivityopfactory> (new ebconductivityopfactory(levelgrids,
										 quadcfi,
										 fastFR,
										 alpha,
										 beta,
										 m_aco.get_data(),
										 m_diffco.get_data(),
										 m_diffco_eb.get_data(),
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

  const int finest_level       = m_amr->get_finest_level();
  const ProblemDomain coar_dom = m_amr->get_domains()[0];

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
				    1.E-99); // Residue set through other means

}

void cdr_tga::setup_tga(){
  CH_TIME("cdr_tga::setup_tga");
  if(m_verbosity > 5){
    pout() << m_name + "::setup_tga" << endl;
  }
  
  const int finest_level       = m_amr->get_finest_level();
  const ProblemDomain coar_dom = m_amr->get_domains()[0];
  const Vector<int> ref_rat    = m_amr->get_ref_rat();

  m_tgasolver = RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > >
    (new AMRTGA<LevelData<EBCellFAB> > (m_gmg_solver, *m_opfact, coar_dom, ref_rat, 1 + finest_level, m_gmg_solver->m_verbosity));

  // Must init gmg for TGA
  Vector<LevelData<EBCellFAB>* > phi, rhs;
  m_amr->alias(phi, m_state);
  m_amr->alias(rhs, m_source);
  m_gmg_solver->init(phi, rhs, finest_level, 0);
}

void cdr_tga::setup_euler(){
  CH_TIME("cdr_tga::setup_euler");
  if(m_verbosity > 5){
    pout() << m_name + "::setup_euler" << endl;
  }

  const int finest_level       = m_amr->get_finest_level();
  const ProblemDomain coar_dom = m_amr->get_domains()[0];
  const Vector<int> ref_rat    = m_amr->get_ref_rat();

  m_eulersolver = RefCountedPtr<EBBackwardEuler> 
    (new EBBackwardEuler (m_gmg_solver, *m_opfact, coar_dom, ref_rat, 1 + finest_level, m_gmg_solver->m_verbosity));

  // Note: If this crashes, try to init gmg first
}

void cdr_tga::compute_divJ(EBAMRCellData& a_divJ, EBAMRCellData& a_state, const Real a_extrap_dt, const bool a_ebflux){
  CH_TIME("cdr_tga::compute_divJ(divF, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_divJ(divF, state)" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;

  // Fill ghost cells
  m_amr->interp_ghost_pwl(a_state, m_realm, m_phase);

  if(m_mobile || m_diffusive){

    if(m_redist_mass_weighted){
      this->reset_redist_weights(a_state);
    }

    // We will let m_scratchFluxOne hold the total flux = advection + diffusion fluxes
    data_ops::set_value(m_scratchFluxOne, 0.0);

    // Compute advection flux. This is mostly the same as compute_divF
    if(m_mobile){
      m_amr->interp_ghost_pwl(m_velo_cell, m_realm, m_phase);
      
      this->average_velo_to_faces(); // Update m_velo_face from m_velo_cell
      this->advect_to_faces(m_face_states, a_state, a_extrap_dt); // Advect to faces
      this->compute_flux(m_scratchFluxTwo, m_velo_face, m_face_states, m_domainflux);

      data_ops::incr(m_scratchFluxOne, m_scratchFluxTwo, 1.0);
    }

    // Compute diffusion flux. 
    if(m_diffusive){
      this->compute_diffusion_flux(m_scratchFluxTwo, a_state);
      data_ops::incr(m_scratchFluxOne, m_scratchFluxTwo, -1.0);
    }

    // General divergence computation. Also inject charge. Domain fluxes came in through the compute
    // advective flux function and eb fluxes come in through the divergence computation.
    EBAMRIVData* ebflux;
    if(a_ebflux){
      ebflux = &m_ebflux;
    }
    else{
      ebflux = &m_eb_zero;
    }
    this->compute_divG(a_divJ, m_scratchFluxOne, *ebflux);
  }
  else{ 
    data_ops::set_value(a_divJ, 0.0);
  }

  return;
}

void cdr_tga::compute_divF(EBAMRCellData& a_divF, EBAMRCellData& a_state, const Real a_extrap_dt, const bool a_ebflux){
  CH_TIME("cdr_tga::compute_divF(divF, state)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_divF(divF, state)" << endl;
  }

  if(m_mobile){

    // Fill ghost cells
    m_amr->interp_ghost_pwl(a_state,     m_realm, m_phase);
    m_amr->interp_ghost_pwl(m_velo_cell, m_realm, m_phase);

    if(m_redist_mass_weighted){
      this->reset_redist_weights(a_state);
    }
    this->average_velo_to_faces();
    this->advect_to_faces(m_face_states, a_state, a_extrap_dt);          // Face extrapolation to cell-centered faces
    this->compute_flux(m_scratchFluxOne, m_velo_face, m_face_states, m_domainflux);  // Compute face-centered fluxes

    EBAMRIVData* ebflux;
    if(a_ebflux){
      ebflux = &m_ebflux;
    }
    else{
      ebflux = &m_eb_zero;
    }
    this->compute_divG(a_divF, m_scratchFluxOne, *ebflux); 
  }
  else{
    data_ops::set_value(a_divF, 0.0);
  }

#if 0 // Debug
  Real volume;
  const Real sum = EBLevelDataOps::noKappaSumLevel(volume, *a_divF[0], 0, m_amr->get_domains()[0]);

  std::cout << sum << std::endl;
#endif

}

void cdr_tga::compute_divD(EBAMRCellData& a_divD, EBAMRCellData& a_state, const bool a_ebflux){
  CH_TIME("cdr_tga::compute_divD");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_divD" << endl;
  }

  if(m_diffusive){
    const int ncomp = 1;

    // Fill ghost cells
    m_amr->interp_ghost_pwl(a_state, m_realm, m_phase);

    if(m_redist_mass_weighted){
      this->reset_redist_weights(a_state);
    }

    this->compute_diffusion_flux(m_scratchFluxOne, a_state);  // Compute the face-centered diffusion flux

    EBAMRIVData* ebflux;
    if(a_ebflux){
      ebflux = &m_ebflux;
    }
    else{
      ebflux = &m_eb_zero;
    }
    this->compute_divG(a_divD, m_scratchFluxOne, *ebflux); // General face-centered flux to divergence magic.

    m_amr->average_down(a_divD, m_realm, m_phase);
    m_amr->interp_ghost(a_divD, m_realm, m_phase);
  }
  else{
    data_ops::set_value(a_divD, 0.0);
  }
}

void cdr_tga::parse_gmg_settings(){
  ParmParse pp(m_class_name.c_str());

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

void cdr_tga::write_plot_data(EBAMRCellData& a_output, int& a_comp){
  CH_TIME("cdr_tga::write_plot_data");
  if(m_verbosity > 5){
    pout() << m_name + "::write_plot_data" << endl;
  }

  if(m_plot_phi) {
    if(!m_plot_numbers){ // Regular write
      write_data(a_output, a_comp, m_state, true);
    }
    else{ // Scale, write, and scale back
     
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	data_ops::scale(*m_state[lvl], (pow(m_amr->get_dx()[lvl], 3)));
      }
      write_data(a_output, a_comp, m_state,    false);
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	data_ops::scale(*m_state[lvl], 1./(pow(m_amr->get_dx()[lvl], 3)));
      }
    }
  }

  if(m_plot_dco && m_diffusive) { // Need to compute the cell-centerd stuff first
    data_ops::set_value(m_scratch, 0.0);
    data_ops::average_face_to_cell(m_scratch, m_diffco, m_amr->get_domains());
    write_data(a_output, a_comp, m_scratch,   false);
  }

  if(m_plot_src) {
    if(!m_plot_numbers){
      write_data(a_output, a_comp, m_source,    false);
    }
    else { // Scale, write, and scale back
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	data_ops::scale(*m_source[lvl], (pow(m_amr->get_dx()[lvl], 3)));
      }
      write_data(a_output, a_comp, m_source,    false);
      for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
	data_ops::scale(*m_source[lvl], 1./(pow(m_amr->get_dx()[lvl], 3)));
      }
    }
  }

  if(m_plot_vel && m_mobile) {
    write_data(a_output, a_comp, m_velo_cell, false);
  }

  // Plot EB fluxes
  if(m_plot_ebf && m_mobile){
    data_ops::set_value(m_scratch, 0.0);
    data_ops::incr(m_scratch, m_ebflux, 1.0);
    write_data(a_output, a_comp, m_scratch, false);
  }
}
