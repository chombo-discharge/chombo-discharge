/*!
  @file   sisdc.cpp
  @brief  Implementation of sisdc.H
  @author Robert Marskar
  @date   Sept. 2018
*/

#include "sisdc.H"
#include "sisdc_storage.H"
#include "data_ops.H"
#include "units.H"
#include "cdr_tga.H"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <ParmParse.H>

typedef sisdc::cdr_storage     cdr_storage;
typedef sisdc::poisson_storage poisson_storage;
typedef sisdc::rte_storage     rte_storage;
typedef sisdc::sigma_storage   sigma_storage;

sisdc::sisdc(){
  m_order        = 2;
  m_min_order    = 2;
  m_max_order    = 2;
  m_error_norm   = 2;
  m_minCFL       = 0.1;
  m_maxCFL       = 0.9;   
  m_err_thresh   = 1.E-4;
  m_safety       = 0.9;

  m_adaptive_dt    = false;
  m_adaptive_order = false;
  m_have_dtf       = false;

  // Basically only for debugging
  m_compute_v      = true;
  m_compute_S      = true;
  m_compute_D      = true;
  m_consistent_E   = true;
  m_consistent_rte = true;
  m_do_advec_src   = true;
  m_do_diffusion   = true;
  m_do_rte         = true;
  m_do_poisson     = true;
  m_print_diagno   = false;
  m_write_diagno   = false;

  // Get parameters from input script
  {
    ParmParse pp("sisdc");

    std::string str;

    pp.query("order",           m_order);
    pp.query("min_order",       m_min_order);
    pp.query("max_order",       m_max_order);
    pp.query("error_norm",      m_error_norm);
    pp.query("min_cfl",         m_minCFL);
    pp.query("max_cfl",         m_maxCFL);
    pp.query("max_error",       m_err_thresh);
    pp.query("safety",          m_safety);

    if(pp.contains("adaptive_dt")){
      pp.get("adaptive_dt", str);
      if(str == "true"){
	m_adaptive_dt = true;
      }
    }
    if(pp.contains("adaptive_order")){
      pp.get("adaptive_order", str);
      if(str == "true"){
	m_adaptive_order = true;
      }
    }
    if(pp.contains("compute_D")){
      pp.get("compute_D", str);
      if(str == "false"){
	m_compute_D = false;
      }
    }
    if(pp.contains("compute_v")){
      pp.get("compute_v", str);
      if(str == "false"){
	m_compute_v = false;
      }
    }
    if(pp.contains("compute_S")){
      pp.get("compute_S", str);
      if(str == "false"){
	m_compute_S = false;
      }
    }
    if(pp.contains("consistent_E")){
      pp.get("consistent_E", str);
      if(str == "false"){
	m_consistent_E = false;
      }
    }
    if(pp.contains("consistent_rte")){
      pp.get("consistent_rte", str);
      if(str == "false"){
	m_consistent_rte = false;
      }
    }
    if(pp.contains("print_diagnostics")){
      pp.get("print_diagnostics", str);
      if(str == "true"){
	m_print_diagno = true;
      }
    }
    if(pp.contains("write_diagnostics")){
      pp.get("write_diagnostics", str);
      if(str == "true"){
	m_write_diagno = true;
      }
    }
    if(pp.contains("do_advec_src")){
      pp.get("do_advec_src", str);
      if(str == "false"){
	m_do_advec_src = false;
	if(m_verbosity > 2){
	  pout() << "sisdc::sisdc - Turning off advection & source" << endl;
	}
      }
    }
    if(pp.contains("do_diffusion")){
      pp.get("do_diffusion", str);
      if(str == "false"){
	m_do_diffusion = false;
	if(m_verbosity > 2){
	  pout() << "sisdc::sisdc - Turning off diffusion" << endl;
	}
      }
    }
    if(pp.contains("do_rte")){
      pp.get("do_rte", str);
      if(str == "false"){
	m_do_rte = false;

	if(m_verbosity > 2){
	  pout() << "sisdc::sisdc - Turning off rte" << endl;
	}
      }
    }
    if(pp.contains("do_poisson")){
      pp.get("do_poisson", str);
      if(str == "false"){
	m_do_poisson = false;

	if(m_verbosity > 2){
	  pout() << "sisdc::sisdc - Turning off poisson" << endl;
	}
      }
    }
  }
}

sisdc::~sisdc(){

}

RefCountedPtr<cdr_storage>& sisdc::get_cdr_storage(const cdr_iterator& a_solverit){
  return m_cdr_scratch[a_solverit.get_solver()];
}

RefCountedPtr<rte_storage>& sisdc::get_rte_storage(const rte_iterator& a_solverit){
  return m_rte_scratch[a_solverit.get_solver()];
}

Real sisdc::restrict_dt(){
  return 1.E99;
}

Real sisdc::get_max_error(){
  CH_TIME("sisdc::get_max_error");
  if(m_verbosity > 2){
    pout() << "sisdc::get_max_error" << endl;
  }

  Real cur_err = m_sigma_error;
  for (int i = 0; i < m_cdr_error.size(); i++){
    cur_err = Max(cur_err, m_cdr_error[i]);
  }

  return cur_err;
}

Real sisdc::advance(const Real a_dt){
  CH_TIME("sisdc::advance");
  if(m_verbosity > 2){
    pout() << "sisdc::advance" << endl;
  }

  // TLDR: This timestepper can do multiple steps per a_dt. If you're using adaptive time integration,
  //       the number of steps are adjusted such that the invididual step errors remain below a specified
  //       threshold. A side effect of this is that the final substep may cause integration beyond a_dt.
  //
  //       PS: When we enter this routine, solvers SHOULD have been filled with valid ready and be ready 
  //           advancement. If you think that this may not be the case, activate the debugging below
#if 0 // Debugging
  this->compute_E_into_scratch();
  this->compute_cdr_gradients();
  this->compute_cdr_velo(m_time);
  this->compute_cdr_sources(m_time);
  time_stepper::compute_cdr_diffusion(m_poisson_scratch->get_E_cell(), m_poisson_scratch->get_E_eb());
#endif

  this->store_previous_solutions(); // Store old solution. This stores the old solutions in storage->m_backup

  MayDay::Abort("stop");
  
  return 0.0;
}

void sisdc::setup_gauss_lobatto(const int a_order){
  CH_TIME("sisdc::regrid_internals");
  if(m_verbosity > 5){
    pout() << "sisdc::regrid_internals" << endl;
  }

  // TLDR: The nodes and weights are hardcoded. A better programmer would compute these
  //       recursively with Legendre polynomials. 

  if(a_order == 1){
    m_gl_weights.resize(2);
    m_gl_nodes.resize(2);
  }
  else{
    m_gl_weights.resize(a_order);
    m_gl_nodes.resize(a_order);
  }

  if(a_order == 1 || a_order == 2){
    m_gl_weights[0] = 1.0;
    m_gl_weights[1] = 1.0;
    m_gl_nodes[0]   = -1.0;
    m_gl_nodes[1]   =  1.0;
  }
  else if(a_order == 3){
    m_gl_nodes[0] = -1.0;
    m_gl_nodes[1] =  0.0;
    m_gl_nodes[2] =  1.0;

    m_gl_weights[0] = 1./3.;
    m_gl_weights[1] = 4./3.;
    m_gl_weights[2] = 1./3.;
  }
  else if(a_order == 4){
    m_gl_nodes[0] = -1.0;
    m_gl_nodes[1] = -1./sqrt(5.);
    m_gl_nodes[2] =  1./sqrt(5.);
    m_gl_nodes[3] =  1.0;

    m_gl_weights[0] = 1./6.;
    m_gl_weights[1] = 5./6.;
    m_gl_weights[2] = 5./6.;
    m_gl_weights[3] = 1./6.;
  }
  else if(a_order == 5){
    m_gl_nodes[0] = -1.0;
    m_gl_nodes[1] = -sqrt(3./7);
    m_gl_nodes[2] =  0.0;
    m_gl_nodes[3] =  sqrt(3./7);
    m_gl_nodes[4] =  1.0;

    m_gl_weights[0] = 1./10.;
    m_gl_weights[1] = 49./90.;
    m_gl_weights[2] = 32./45;
    m_gl_weights[3] = 49./90.;
    m_gl_weights[4] = 1./10.;
  }
  else if(a_order == 6){
    m_gl_nodes[0] = -1.0;
    m_gl_nodes[1] = -0.76595532;
    m_gl_nodes[2] = -0.28532152;
    m_gl_nodes[3] =  0.28532152;
    m_gl_nodes[4] =  0.76595532;
    m_gl_nodes[5] =  1.0;

    m_gl_weights[0] = 1./15;
    m_gl_weights[1] = 0.37847496;
    m_gl_weights[2] = 0.55485838;
    m_gl_weights[3] = 0.55485838;
    m_gl_weights[4] = 0.37847496;
    m_gl_weights[5] = 1./15;
  }
  else if(a_order == 7){
    m_gl_nodes[0] = -1.0;
    m_gl_nodes[1] = -0.83022390;
    m_gl_nodes[2] = -0.46884879;
    m_gl_nodes[3] =  0.0;
    m_gl_nodes[4] =  0.46884879;
    m_gl_nodes[5] =  0.83022390;
    m_gl_nodes[6] =  1.0;

    m_gl_weights[0] = 0.04761904;
    m_gl_weights[1] = 0.27682604;
    m_gl_weights[2] = 0.43174538;
    m_gl_weights[3] = 0.48761904;
    m_gl_weights[4] = 0.43174538;
    m_gl_weights[5] = 0.27682604;
    m_gl_weights[6] = 0.04761904;
  }
  else{
    MayDay::Abort("sisdc::compute_gauss_lobatto - requested order exceeds. If you want order > 7, compute your own damn weights!");
  }
}

void sisdc::setup_subintervals(const int a_order, const Real a_time, const Real a_dt){
  CH_TIME("sisdc::setup_subintervals");
  if(m_verbosity > 5){
    pout() << "sisdc::setup_subintervals" << endl;
  }

  // m_gl_nodes are Gauss-Lobatto nodes on [-1,1]. These must
  // be shifted to [t_n,t_n + a_dt]
  m_tm.resize(a_order);
  Vector<Real> shifted_gl_nodes = m_gl_nodes;
  for (int m = 0; m < shifted_gl_nodes.size(); m++){

    shifted_gl_nodes[m] += 1.0;    // [0,2]
    shifted_gl_nodes[m] *= 0.5;    // [0,1]
    shifted_gl_nodes[m] *= a_dt;   // [0, a_dt]
    shifted_gl_nodes[m] += a_time; // [a_time, a_time + a_dt]
  }

  // dtm = t_{m+1} - t_m. Order 1 is special since we only use the SISDC predictor from a second order formulation
  if(a_order == 1){
    m_dtm.resize(1);
    m_dtm = m_tm[1] - m_tm[0];
  }
  else{
    m_dtm.resize(a_order - 1);
    for (int m = 0; m < m_tm.size()-1; m++){
      m_dtm[m] = m_tm[m+1] - m_tm[m];
    }
  }
}

void sisdc::compute_dt(Real& a_dt, time_code::which_code& a_timecode){
  CH_TIME("sisdc::compute_dt");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_dt" << endl;
  }

  Real dt = 1.E99;

  m_dt_cfl = m_cdr->compute_cfl_dt();
  const Real dt_cfl = m_cfl*m_dt_cfl;
  if(dt_cfl < dt){
    dt = dt_cfl;
    a_timecode = time_code::cfl;
  }

  const Real dt_src = m_src_growth*m_cdr->compute_source_dt(m_src_tolerance, m_src_elec_only);
  if(dt_src < dt){
    dt = dt_src;
    a_timecode = time_code::source;
  }

  const Real dt_relax = m_relax_time*this->compute_relaxation_time();
  if(dt_relax < dt){
    dt = dt_relax;
    a_timecode = time_code::relaxation_time;
  }

  const Real dt_restrict = this->restrict_dt();
  if(dt_restrict < dt){
    dt = dt_restrict;
    a_timecode = time_code::restricted;
  }

  if(dt < m_min_dt){
    dt = m_min_dt;
    a_timecode = time_code::hardcap;
  }

  if(dt > m_max_dt){
    dt = m_max_dt;
    a_timecode = time_code::hardcap;
  }

  a_dt = dt;

  // Copy the time code, it is needed for diagnostics
  m_timecode = a_timecode;
}

void sisdc::regrid_internals(){
  CH_TIME("sisdc::regrid_internals");
  if(m_verbosity > 5){
    pout() << "sisdc::regrid_internals" << endl;
  }

  m_cdr_error.resize(m_plaskin->get_num_species());
  
  this->allocate_cdr_storage();
  this->allocate_poisson_storage();
  this->allocate_rte_storage();
  this->allocate_sigma_storage();

  //  this->build_gauss_lobatto(m_order);
}

void sisdc::allocate_cdr_storage(){
  const int ncomp       = 1;
  const int num_species = m_plaskin->get_num_species();

  m_cdr_scratch.resize(num_species);
  
  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx] = RefCountedPtr<cdr_storage> (new cdr_storage(m_amr, m_cdr->get_phase(), ncomp));
    m_cdr_scratch[idx]->allocate_storage(m_order);
  }
}

void sisdc::allocate_poisson_storage(){
  const int ncomp = 1;
  m_poisson_scratch = RefCountedPtr<poisson_storage> (new poisson_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_poisson_scratch->allocate_storage(m_order);
}

void sisdc::allocate_rte_storage(){
  const int ncomp       = 1;
  const int num_photons = m_plaskin->get_num_photons();
  m_rte_scratch.resize(num_photons);
  
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx] = RefCountedPtr<rte_storage> (new rte_storage(m_amr, m_rte->get_phase(), ncomp));
    m_rte_scratch[idx]->allocate_storage(m_order);
  }
}

void sisdc::allocate_sigma_storage(){
  const int ncomp = 1;
  m_sigma_scratch = RefCountedPtr<sigma_storage> (new sigma_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_sigma_scratch->allocate_storage(m_order);
}

void sisdc::deallocate_internals(){
  CH_TIME("sisdc::deallocate_internals");
  if(m_verbosity > 5){
    pout() << "sisdc::deallocate_internals" << endl;
  }
  
  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx]->deallocate_storage();
    m_cdr_scratch[idx] = RefCountedPtr<cdr_storage>(0);
  }

  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx]->deallocate_storage();
    m_rte_scratch[idx] = RefCountedPtr<rte_storage>(0);
  }

  m_cdr_scratch.resize(0);
  m_rte_scratch.resize(0);

  m_poisson_scratch->deallocate_storage();
  m_poisson_scratch = RefCountedPtr<poisson_storage>(0);
  
  m_sigma_scratch->deallocate_storage();
  m_sigma_scratch = RefCountedPtr<sigma_storage>(0);
}

void sisdc::store_previous_solutions(){
  CH_TIME("sisdc::store_previous_solutions");
  if(m_verbosity > 5){
    pout() << "sisdc::store_previous_solutions" << endl;
  }
  
  // Backup cdr solutions
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();

    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    EBAMRCellData& previous = storage->get_previous();

    data_ops::copy(previous, solver->get_state());
  }

  {// Backup Poisson solution
    MFAMRCellData& previous = m_poisson_scratch->get_previous();
    data_ops::copy(previous, m_poisson->get_state());
  }

  // Backup RTE solutions
  for (rte_iterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<rte_solver>& solver = solver_it();

    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);
    EBAMRCellData& previous = storage->get_previous();

    data_ops::copy(previous, solver->get_state());
  }

  { // Backup sigma
    EBAMRIVData& previous = m_sigma_scratch->get_previous();
    data_ops::copy(previous, m_sigma->get_state());
  }
}

void sisdc::restore_previous_solutions(){
  CH_TIME("sisdc::restore_previous_solutions");
  if(m_verbosity > 5){
    pout() << "sisdc::restore_previous_solutions" << endl;
  }
  
  // Revert cdr solutions
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();

    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    const EBAMRCellData& previous = storage->get_previous();

    data_ops::copy(solver->get_state(), previous);
  }

  {// Revert Poisson solution
    const MFAMRCellData& previous = m_poisson_scratch->get_previous();
    data_ops::copy(m_poisson->get_state(), previous);
  }

  // Revert RTE solutions
  for (rte_iterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<rte_solver>& solver = solver_it();

    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);
    const EBAMRCellData& previous = storage->get_previous();

    data_ops::copy(solver->get_state(), previous);
  }

  { // Revert sigma solution
    const EBAMRIVData& previous = m_sigma_scratch->get_previous();
    data_ops::copy(m_sigma->get_state(), previous);
  }
}

void sisdc::compute_E_into_scratch(){
  CH_TIME("sisdc::compute_E_into_scratch");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_E_into_scratch" << endl;
  }
  
  EBAMRCellData& E_cell = m_poisson_scratch->get_E_cell();
  EBAMRFluxData& E_face = m_poisson_scratch->get_E_face();
  EBAMRIVData&   E_eb   = m_poisson_scratch->get_E_eb();
  EBAMRIFData&   E_dom  = m_poisson_scratch->get_E_domain();

  const MFAMRCellData& phi = m_poisson->get_state();
  
  this->compute_E(E_cell, m_cdr->get_phase(), phi);     // Compute cell-centered field
  this->compute_E(E_face, m_cdr->get_phase(), E_cell);  // Compute face-centered field
  this->compute_E(E_eb,   m_cdr->get_phase(), E_cell);  // EB-centered field

  time_stepper::extrapolate_to_domain_faces(E_dom, m_cdr->get_phase(), E_cell);
}

void sisdc::compute_cdr_gradients(){
  CH_TIME("sisdc::compute_cdr_gradients");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_cdr_gradients" << endl;
  }

  sisdc::compute_cdr_gradients(m_cdr->get_states());
}

void sisdc::compute_cdr_gradients(const Vector<EBAMRCellData*>& a_states){
  CH_TIME("sisdc::compute_cdr_gradients");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_cdr_gradients" << endl;
  }

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    EBAMRCellData& grad = storage->get_gradient();

    m_amr->compute_gradient(grad, *a_states[idx]);
    m_amr->average_down(grad, m_cdr->get_phase());
    m_amr->interp_ghost(grad, m_cdr->get_phase());
  }
}

void sisdc::compute_cdr_velo(const Real a_time){
  CH_TIME("sisdc::compute_cdr_velo");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_cdr_velo" << endl;
  }

  this->compute_cdr_velo(m_cdr->get_states(), a_time);
}

void sisdc::compute_cdr_velo(const Vector<EBAMRCellData*>& a_states, const Real a_time){
  CH_TIME("sisdc::compute_cdr_velo(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_cdr_velo(Vector<EBAMRCellData*>, Real)" << endl;
  }

  Vector<EBAMRCellData*> velocities = m_cdr->get_velocities();
  this->compute_cdr_velocities(velocities, a_states, m_poisson_scratch->get_E_cell(), a_time);
}

void sisdc::compute_cdr_eb_states(){
  CH_TIME("sisdc::compute_cdr_eb_states");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_cdr_eb_states" << endl;
  }

  Vector<EBAMRIVData*>   eb_gradients;
  Vector<EBAMRIVData*>   eb_states;
  Vector<EBAMRCellData*> cdr_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    cdr_states.push_back(&(solver->get_state()));
    eb_states.push_back(&(storage->get_eb_state()));
    eb_gradients.push_back(&(storage->get_eb_grad()));
    cdr_gradients.push_back(&(storage->get_gradient())); // Should already have been computed
  }

  // Extrapolate states to the EB and floor them so we cannot get negative values on the boundary. This
  // won't hurt mass conservation because the mass hasn't been injected yet
  this->extrapolate_to_eb(eb_states, m_cdr->get_phase(), cdr_states);
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    data_ops::floor(*eb_states[idx], 0.0);
  }

  // We should already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIVData eb_gradient;
  m_amr->allocate(eb_gradient, m_cdr->get_phase(), SpaceDim);
  for (int i = 0; i < cdr_states.size(); i++){
    this->extrapolate_to_eb(eb_gradient, m_cdr->get_phase(), *cdr_gradients[i]);
    this->project_flux(*eb_gradients[i], eb_gradient);
  }
}

void sisdc::compute_cdr_eb_states(const Vector<EBAMRCellData*>& a_states){
  CH_TIME("sisdc::compute_cdr_eb_states(vec)");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_cdr_eb_states(vec)" << endl;
  }

  Vector<EBAMRIVData*>   eb_gradients;
  Vector<EBAMRIVData*>   eb_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    eb_states.push_back(&(storage->get_eb_state()));
    eb_gradients.push_back(&(storage->get_eb_grad()));
    cdr_gradients.push_back(&(storage->get_gradient())); // Should already have been computed
  }

  // Extrapolate states to the EB and floor them so we cannot get negative values on the boundary. This
  // won't hurt mass conservation because the mass hasn't been injected yet
  this->extrapolate_to_eb(eb_states, m_cdr->get_phase(), a_states);
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    data_ops::floor(*eb_states[idx], 0.0);
  }

  // We should already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIVData eb_gradient;
  m_amr->allocate(eb_gradient, m_cdr->get_phase(), SpaceDim);
  for (int i = 0; i < a_states.size(); i++){
    this->extrapolate_to_eb(eb_gradient, m_cdr->get_phase(), *cdr_gradients[i]);
    this->project_flux(*eb_gradients[i], eb_gradient);
  }
}

void sisdc::compute_cdr_domain_states(){
  CH_TIME("sisdc::compute_cdr_domain_states");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_cdr_domain_states" << endl;
  }

  Vector<EBAMRIFData*>   domain_gradients;
  Vector<EBAMRIFData*>   domain_states;
  Vector<EBAMRCellData*> cdr_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    cdr_states.push_back(&(solver->get_state()));
    domain_states.push_back(&(storage->get_domain_state()));
    domain_gradients.push_back(&(storage->get_domain_grad()));
    cdr_gradients.push_back(&(storage->get_gradient())); // Should already be computed
  }

  // Extrapolate states to the domain faces
  this->extrapolate_to_domain_faces(domain_states, m_cdr->get_phase(), cdr_states);

  // We already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIFData grad;
  m_amr->allocate(grad, m_cdr->get_phase(), SpaceDim);
  for (int i = 0; i < cdr_states.size(); i++){
    this->extrapolate_to_domain_faces(grad, m_cdr->get_phase(), *cdr_gradients[i]);
    this->project_domain(*domain_gradients[i], grad);
  }
}

void sisdc::compute_cdr_domain_states(const Vector<EBAMRCellData*>& a_states){
  CH_TIME("sisdc::compute_cdr_domain_states");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_cdr_domain_states" << endl;
  }

  Vector<EBAMRIFData*>   domain_gradients;
  Vector<EBAMRIFData*>   domain_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    domain_states.push_back(&(storage->get_domain_state()));
    domain_gradients.push_back(&(storage->get_domain_grad()));
    cdr_gradients.push_back(&(storage->get_gradient()));
  }

  // Extrapolate states to the domain faces
  this->extrapolate_to_domain_faces(domain_states, m_cdr->get_phase(), a_states);

  // We already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIFData grad;
  m_amr->allocate(grad, m_cdr->get_phase(), SpaceDim);
  for (int i = 0; i < a_states.size(); i++){
    this->extrapolate_to_domain_faces(grad, m_cdr->get_phase(), *cdr_gradients[i]);
    this->project_domain(*domain_gradients[i], grad);
  }
}

void sisdc::compute_cdr_fluxes(const Real a_time){
  CH_TIME("sisdc::compute_cdr_fluxes");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_cdr_fluxes" << endl;
  }

  this->compute_cdr_fluxes(m_cdr->get_states(), a_time);
}

void sisdc::compute_cdr_fluxes(const Vector<EBAMRCellData*>& a_states, const Real a_time){
  CH_TIME("sisdc::compute_cdr_fluxes(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_cdr_fluxes(Vector<EBAMRCellData*>, Real)" << endl;
  }

  Vector<EBAMRIVData*> cdr_fluxes;
  Vector<EBAMRIVData*> extrap_cdr_fluxes;
  Vector<EBAMRIVData*> extrap_cdr_densities;
  Vector<EBAMRIVData*> extrap_cdr_velocities;
  Vector<EBAMRIVData*> extrap_cdr_gradients;
  Vector<EBAMRIVData*> extrap_rte_fluxes;

  cdr_fluxes = m_cdr->get_ebflux();

  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRIVData& dens_eb = storage->get_eb_state();
    EBAMRIVData& velo_eb = storage->get_eb_velo();
    EBAMRIVData& flux_eb = storage->get_eb_flux();
    EBAMRIVData& grad_eb = storage->get_eb_grad();

    extrap_cdr_densities.push_back(&dens_eb);  // Computed in compute_cdr_eb_states
    extrap_cdr_velocities.push_back(&velo_eb); // Not yet computed
    extrap_cdr_fluxes.push_back(&flux_eb);     // Not yet computed
    extrap_cdr_gradients.push_back(&grad_eb);  // Computed in compute_cdr_eb_states
  }


  // Extrapolate densities, velocities, and fluxes
  Vector<EBAMRCellData*> cdr_velocities = m_cdr->get_velocities();
  this->compute_extrapolated_fluxes(extrap_cdr_fluxes, a_states, cdr_velocities, m_cdr->get_phase());
  this->extrapolate_to_eb(extrap_cdr_velocities, m_cdr->get_phase(), cdr_velocities);

  // Compute RTE flux on the boundary
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);

    EBAMRIVData& flux_eb = storage->get_eb_flux();
    solver->compute_boundary_flux(flux_eb, solver->get_state());
    extrap_rte_fluxes.push_back(&flux_eb);
  }

  const EBAMRIVData& E = m_poisson_scratch->get_E_eb();
  time_stepper::compute_cdr_fluxes(cdr_fluxes,
				   extrap_cdr_fluxes,
				   extrap_cdr_densities,
				   extrap_cdr_velocities,
				   extrap_cdr_gradients,
				   extrap_rte_fluxes,
				   E,
				   a_time);
}

void sisdc::compute_cdr_domain_fluxes(const Real a_time){
  CH_TIME("sisdc::compute_cdr_domain_fluxes");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_cdr_domain_fluxes" << endl;
  }

  this->compute_cdr_domain_fluxes(m_cdr->get_states(), a_time);
}

void sisdc::compute_cdr_domain_fluxes(const Vector<EBAMRCellData*>& a_states, const Real a_time){
  CH_TIME("sisdc::compute_cdr_domain_fluxes(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_cdr_domain_fluxes(Vector<EBAMRCellData*>, Real)" << endl;
  }

  Vector<EBAMRIFData*>   cdr_fluxes;
  Vector<EBAMRIFData*>   extrap_cdr_fluxes;
  Vector<EBAMRIFData*>   extrap_cdr_densities;
  Vector<EBAMRIFData*>   extrap_cdr_velocities;
  Vector<EBAMRIFData*>   extrap_cdr_gradients;
  Vector<EBAMRIFData*>   extrap_rte_fluxes;

  Vector<EBAMRCellData*> cdr_velocities;
  Vector<EBAMRCellData*> cdr_gradients;

  cdr_fluxes = m_cdr->get_domainflux();
  cdr_velocities = m_cdr->get_velocities();
  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRIFData& dens_domain = storage->get_domain_state();
    EBAMRIFData& velo_domain = storage->get_domain_velo();
    EBAMRIFData& flux_domain = storage->get_domain_flux();
    EBAMRIFData& grad_domain = storage->get_domain_grad();
    EBAMRCellData& gradient  = storage->get_gradient();

    extrap_cdr_densities.push_back(&dens_domain);  // Has not been computed
    extrap_cdr_velocities.push_back(&velo_domain); // Has not been computed
    extrap_cdr_fluxes.push_back(&flux_domain);     // Has not been computed
    extrap_cdr_gradients.push_back(&grad_domain);  // Has not been computed
    cdr_gradients.push_back(&gradient);
  }

  // Compute extrapolated velocities and fluxes at the domain faces
  this->extrapolate_to_domain_faces(extrap_cdr_densities,         m_cdr->get_phase(), a_states);
  this->extrapolate_vector_to_domain_faces(extrap_cdr_velocities, m_cdr->get_phase(), cdr_velocities);
  this->compute_extrapolated_domain_fluxes(extrap_cdr_fluxes,     a_states,           cdr_velocities, m_cdr->get_phase());
  this->extrapolate_vector_to_domain_faces(extrap_cdr_gradients,  m_cdr->get_phase(), cdr_gradients);

  // Compute RTE flux on domain faces
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);

    EBAMRIFData& domain_flux = storage->get_domain_flux();
    solver->compute_domain_flux(domain_flux, solver->get_state());
    extrap_rte_fluxes.push_back(&domain_flux);
  }

  const EBAMRIFData& E = m_poisson_scratch->get_E_domain();

  // This fills the solvers' domain fluxes
  time_stepper::compute_cdr_domain_fluxes(cdr_fluxes,
					  extrap_cdr_fluxes,
					  extrap_cdr_densities,
					  extrap_cdr_velocities,
					  extrap_cdr_gradients,
					  extrap_rte_fluxes,
					  E,
					  a_time);
}

void sisdc::compute_sigma_flux(){
  CH_TIME("sisdc::compute_sigma_flux");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_sigma_flux" << endl;
  }

  EBAMRIVData& flux = m_sigma->get_flux();
  data_ops::set_value(flux, 0.0);

  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    const RefCountedPtr<species>& spec      = solver_it.get_species();
    const EBAMRIVData& solver_flux          = solver->get_ebflux();

    data_ops::incr(flux, solver_flux, spec->get_charge()*units::s_Qe);
  }

  m_sigma->reset_cells(flux);
}

void sisdc::compute_cdr_sources(const Real a_time){
  CH_TIME("sisdc::compute_cdr_sources_into_scratch");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_cdr_sources_into_scratch" << endl;
  }

  this->compute_cdr_sources(m_cdr->get_states(), a_time);
}

void sisdc::compute_cdr_sources(const Vector<EBAMRCellData*>& a_states, const Real a_time){
  CH_TIME("sisdc::compute_cdr_sources(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_cdr_sources(Vector<EBAMRCellData*>, Real)" << endl;
  }
  
  Vector<EBAMRCellData*> cdr_sources = m_cdr->get_sources();
  Vector<EBAMRCellData*> rte_states  = m_rte->get_states();
  EBAMRCellData& E                   = m_poisson_scratch->get_E_cell();

  Vector<EBAMRCellData*> cdr_gradients;
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    cdr_gradients.push_back(&(storage->get_gradient())); // These should already have been computed
  }

  time_stepper::compute_cdr_sources(cdr_sources, a_states, cdr_gradients, rte_states, E, a_time, centering::cell_center);
}

void sisdc::update_poisson(){
  if(m_verbosity > 5){
    pout() << "sisdc::update_poisson" << endl;
  }
  
  if(m_do_poisson){ // Solve Poisson equation
    if((m_step +1) % m_fast_poisson == 0){
      time_stepper::solve_poisson();
      this->compute_E_into_scratch();
    }
  }
}

void sisdc::update_rte(const Real a_time){
  if(m_verbosity > 5){
    pout() << "sisdc::update_rte" << endl;
  }
  
  if(m_do_rte){
    if((m_step + 1) % m_fast_rte == 0){
      Vector<EBAMRCellData*> rte_states  = m_rte->get_states();
      Vector<EBAMRCellData*> rte_sources = m_rte->get_sources();
      Vector<EBAMRCellData*> cdr_states  = m_cdr->get_states();

      EBAMRCellData& E = m_poisson_scratch->get_E_cell();

      const Real dummy_dt = 0.0;
      this->solve_rte(rte_states, rte_sources, cdr_states, E, a_time, dummy_dt, centering::cell_center);
    }
  }
}
