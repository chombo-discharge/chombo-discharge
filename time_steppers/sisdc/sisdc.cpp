/*!
  @file   sisdc.cpp
  @brief  Implementation of sisdc.H
  @author Robert Marskar
  @date   Feb. 2019
*/

#include "sisdc.H"
#include "sisdcF_F.H"
#include "sisdc_storage.H"
#include "data_ops.H"
#include "units.H"
#include "cdr_tga.H"
#include "cdr_fhd.H"
#include "cdr_gdnv.H"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <ParmParse.H>

typedef sisdc::cdr_storage     cdr_storage;
typedef sisdc::poisson_storage poisson_storage;
typedef sisdc::rte_storage     rte_storage;
typedef sisdc::sigma_storage   sigma_storage;

sisdc::sisdc(){
  m_class_name = "sisdc";
}

sisdc::~sisdc(){

}

void sisdc::parse_options(){
  CH_TIME("sisdc::parse_options");
  if(m_verbosity > 5){
    pout() << "sisdc::parse_options" << endl;
  }

  // Regular stuff from time_stepper that we almost always need
  parse_verbosity();
  parse_solver_verbosity();
  parse_cfl();
  parse_relax_time();
  parse_fast_rte();
  parse_fast_poisson();
  parse_min_dt();
  parse_max_dt();
  parse_source_comp();

  // Specific to this class
  parse_nodes();
  parse_diffusion_coupling();
  parse_adaptive_options();
  parse_subcycle_options();
  parse_debug_options();
  parse_advection_options();

  // Setup nodes
  sisdc::setup_quadrature_nodes(m_p);
  sisdc::setup_qmj(m_p);
}

void sisdc::parse_nodes(){
  CH_TIME("sisdc::parse_nodes");
  if(m_verbosity > 5){
    pout() << "sisdc::parse_nodes" << endl;
  }

  ParmParse pp(m_class_name.c_str());
  std::string str;
  
  pp.get("quad_nodes", str);
  if(str == "lobatto"){
    m_which_nodes = "lobatto";
  }
  else if(str == "uniform"){
    m_which_nodes = "uniform";
  }
  else if(str == "chebyshev"){
    m_which_nodes = "chebyshev";
  }
  else {
    MayDay::Abort("sisdc::parse_nodes - unknown node type requested");
  }

  pp.get("subintervals",     m_p);
  pp.get("corr_iter",        m_k);

  if(m_p < 1){
    MayDay::Abort("sisdc::parse_nodes - sisdc.subintervals cannot be < 1");
  }
  if(m_k < 0){
    MayDay::Abort("sisdc::parse_nodes - sisdc.corr_iter cannot be < 0");
  }
}

void sisdc::parse_diffusion_coupling(){
  CH_TIME("sisdc::parse_diffusion_coupling");
  if(m_verbosity > 5){
    pout() << "sisdc::parse_diffusion_coupling" << endl;
  }

  ParmParse pp(m_class_name.c_str());
  
  std::string str;

  pp.get("use_tga", str);
  m_use_tga = (str == "true") ? true : false;

  pp.get("diffusive_coupling", str);
  m_strong_diffu = (str == "true") ? true : false;
  
  pp.get("num_diff_corr",  m_num_diff_corr);
  if(m_num_diff_corr < 0){
    MayDay::Abort("sisdc::parse_diffusion_coupling - option 'sisdc.num_diff_corr' cannot be negative");
  }
}

void sisdc::parse_adaptive_options(){
  CH_TIME("sisdc::parse_adaptive_options");
  if(m_verbosity > 5){
    pout() << "sisdc::parse_adaptive_options" << endl;
  }

  ParmParse pp(m_class_name.c_str());
  std::string str;

  pp.get("error_norm",       m_error_norm);
  pp.get("min_corr",         m_min_corr);
  pp.get("max_retries",      m_max_retries);
  pp.get("max_growth",       m_max_growth);
  pp.get("decrease_safety",  m_decrease_safe);
  pp.get("min_cfl",          m_minCFL);
  pp.get("max_cfl",          m_maxCFL);
  pp.get("max_error",        m_err_thresh);
  pp.get("error_index",      m_error_idx);
  pp.get("safety",           m_safety);

  pp.get("print_report", str);
  m_print_report = (str == "true") ? true : false;
  
  pp.get("adaptive_dt", str);
  m_adaptive_dt = (str == "true") ? true : false;
}

void sisdc::parse_subcycle_options(){
  CH_TIME("sisdc::parse_subcycle_options");
  if(m_verbosity > 5){
    pout() << "sisdc::parse_subcycle_options" << endl;
  }

  ParmParse pp(m_class_name.c_str());
  std::string str;

  pp.get("regrid_cfl",       m_regrid_cfl);
  pp.get("subcycle_cfl",     m_cycleCFL);
  pp.get("min_cycle_cfl",    m_min_cycle_cfl);
  pp.get("max_cycle_cfl",    m_max_cycle_cfl);

  pp.get("subcycle", str);
  if(str == "false"){
    m_subcycle = false;
  }
  else if(str == "standard"){
    m_subcycle = true;
    m_optimal_subcycling = false;
  }
  else if(str == "optimal"){
    m_subcycle = true;
    m_optimal_subcycling = true;
  }
  else if(str == "multistep"){
    m_subcycle  = true;
    m_multistep = true;
  }
  else{
    MayDay::Abort("sisdc::sisdc - unknown sisdc.subcycle = ???");
  }
  
  pp.get("cycle_sources", str);
  m_cycle_sources = (str == "true") ? true : false;
}

void sisdc::parse_debug_options(){
  CH_TIME("sisdc::parse_debug_options");
  if(m_verbosity > 5){
    pout() << "sisdc::parse_debug_options" << endl;
  }

  ParmParse pp(m_class_name.c_str());
  std::string str;

  pp.get("compute_D", str); m_compute_D = (str == "true") ? true : false;
  pp.get("compute_v", str); m_compute_v = (str == "true") ? true : false;
  pp.get("compute_S", str); m_compute_S = (str == "true") ? true : false;
  
  pp.get("consistent_E",   str); m_consistent_E   = (str == "true") ? true : false;
  pp.get("consistent_rte", str); m_consistent_rte = (str == "true") ? true : false;
  
  pp.get("do_advec_src", str); 	m_do_advec_src = (str == "true") ? true : false;
  pp.get("do_diffusion", str); 	m_do_diffusion = (str == "true") ? true : false;
  pp.get("do_rte", str); 	m_do_rte       = (str == "true") ? true : false;
  pp.get("do_poisson", str); 	m_do_poisson   = (str == "true") ? true : false;
  
  pp.get("profile_steps", str); m_profile_steps = (str == "true") ? true : false;
}

void sisdc::parse_advection_options(){
  CH_TIME("sisdc::parse_advection_options");
  if(m_verbosity > 5){
    pout() << "sisdc::parse_advection_options" << endl;
  }

  ParmParse pp(m_class_name.c_str());
  std::string str;

  m_extrap_dt = 0.5;  // Relic of an ancient past. I don't see any reason why extrapolating to anything but the half interval
                      // would make sense. 

  pp.get("extrap_advect", str); m_extrap_advect = (str == "true") ? true : false;
}

RefCountedPtr<cdr_storage>& sisdc::get_cdr_storage(const cdr_iterator& a_solverit){
  return m_cdr_scratch[a_solverit.get_solver()];
}

RefCountedPtr<rte_storage>& sisdc::get_rte_storage(const rte_iterator& a_solverit){
  return m_rte_scratch[a_solverit.get_solver()];
}

bool sisdc::need_to_regrid(){
  CH_TIME("sisdc::need_to_regrid");
  if(m_verbosity > 5){
    pout() << "sisdc::need_to_regrid" << endl;
  }
  const bool regrid = m_accum_cfl > m_regrid_cfl;
  if(regrid) m_accum_cfl = 0.0;
  
  return regrid;
}

Real sisdc::restrict_dt(){
  return 1.E99;
}

Real sisdc::get_max_node_distance(){
  CH_TIME("sisdc::get_max_node_distance");
  if(m_verbosity > 5){
    pout() << "sisdc::get_max_node_distance" << endl;
  }

  Real max_dist = 0.0;
  for (int m = 0; m < m_p; m++){
    max_dist = Max(max_dist, m_nodes[m+1] - m_nodes[m]);
  }

  return max_dist;
}

void sisdc::init_source_terms(){
  CH_TIME("sisdc::init_source_terms");
  if(m_verbosity > 5){
    pout() << "sisdc::init_source_terms" << endl;
  }
  
  time_stepper::init_source_terms();
}

void sisdc::setup_quadrature_nodes(const int a_p){
  CH_TIME("sisdc::setup_quadrature_nodes");
  if(m_verbosity > 5){
    pout() << "sisdc::setup_quadrature_nodes" << endl;
  }

  if(m_which_nodes == "uniform"){
    sisdc::setup_uniform_nodes(a_p);
  }
  else if(m_which_nodes == "lobatto"){
    sisdc::setup_lobatto_nodes(a_p);
  }
  else if(m_which_nodes == "chebyshev"){
    sisdc::setup_chebyshev_nodes(a_p);
  }
  else {
    MayDay::Abort("sisdc::setup_quadrature_nodes - unknown nodes requested");
  }
}

void sisdc::setup_uniform_nodes(const int a_p){
  CH_TIME("sisdc::setup_uniform_nodes");
  if(m_verbosity > 5){
    pout() << "sisdc::setup_uniform_nodes" << endl;
  }

  // TLDR: The nodes and weights are hardcoded. A better programmer would compute these
  //       recursively with Legendre polynomials. 
  m_nodes.resize(1+a_p);

  const Real delta = 2./a_p;
  for (int m = 0; m <= a_p; m++){
    m_nodes[m] = m*delta;
  }
}

void sisdc::setup_lobatto_nodes(const int a_p){
  CH_TIME("sisdc::setup_lobatto_nodes");
  if(m_verbosity > 5){
    pout() << "sisdc::setup_lobatto_nodes" << endl;
  }

  // TLDR: The nodes and weights are hardcoded. A better programmer would compute these
  //       recursively with Legendre polynomials. 
  m_nodes.resize(1+a_p);

  if(a_p == 1){
    m_nodes[0]   = -1.0;
    m_nodes[1]   =  1.0;
  }
  else if(a_p == 2){
    m_nodes[0] = -1.0;
    m_nodes[1] =  0.0;
    m_nodes[2] =  1.0;
  }
  else if(a_p == 3){
    m_nodes[0] = -1.0;
    m_nodes[1] = -1./sqrt(5.);
    m_nodes[2] =  1./sqrt(5.);
    m_nodes[3] =  1.0;
  }
  else if(a_p == 4){
    m_nodes[0] = -1.0;
    m_nodes[1] = -sqrt(3./7);
    m_nodes[2] =  0.0;
    m_nodes[3] =  sqrt(3./7);
    m_nodes[4] =  1.0;
  }
  else if(a_p == 5){
    m_nodes[0] = -1.0;
    m_nodes[1] = -0.76595532;
    m_nodes[2] = -0.28532152;
    m_nodes[3] =  0.28532152;
    m_nodes[4] =  0.76595532;
    m_nodes[5] =  1.0;
  }
  else if(a_p == 6){
    m_nodes[0] = -1.0;
    m_nodes[1] = -0.83022390;
    m_nodes[2] = -0.46884879;
    m_nodes[3] =  0.0;
    m_nodes[4] =  0.46884879;
    m_nodes[5] =  0.83022390;
    m_nodes[6] =  1.0;
  }
  else{
    MayDay::Abort("sisdc::setup_lobatto_nodes - requested order exceeds 7. Compute your own damn nodes!");
  }
}

void sisdc::setup_chebyshev_nodes(const int a_p){
  CH_TIME("sisdc::setup_chebyshev_nodes");
  if(m_verbosity > 5){
    pout() << "sisdc::setup_chebyshev_nodes" << endl;
  }

  // TLDR: The nodes and weights are hardcoded. A better programmer would compute these
  //       recursively with Legendre polynomials. 
  m_nodes.resize(1+a_p);
  m_nodes[0] = -1.0;
  for (int m = 1; m < a_p; m++){
    m_nodes[m] = -cos((2*m-1)*units::s_pi/(2*(a_p-1)));
  }
  m_nodes[a_p] = 1.0;
}

void sisdc::setup_qmj(const int a_p){
  CH_TIME("sisdc::setup_qmj");
  if(m_verbosity > 5){
    pout() << "sisdc::setup_qmj" << endl;
  }

  const int nnodes = 1 + a_p;

  // Resize the integration matrix, it should be p x (p+1)
  m_qmj.resize(a_p);
  for (int m = 0; m < a_p; m++){
    m_qmj[m].resize(nnodes, 0.0);
  }

  // Generate the qmj matrix
  for (int j=0; j < nnodes; j++){

    // Set up the Vandermonde matrix (in Fortran order since we will call LaPack)
    double V[nnodes*nnodes];
    for (int j=0; j<nnodes; j++){
      for (int i=0; i<nnodes; i++){
	const int k = j*nnodes + i;
	V[k] = pow(m_nodes[i],j);
      }
    }

    // Setup f = delta_kj. When we solve, this becomes the solution vector
    double cj[nnodes];
    for (int k=0; k<nnodes; k++){
      cj[k] = (k==j) ? 1.0 : 0.0;
    }

    // Solve V*c = f. This calls LAPACK
    int N    = nnodes;
    int NRHS = 1;
    int LDA  = nnodes;
    int IPIV[nnodes];
    int LDB  = nnodes;
    int INFO = 10;
    dgesv_(&N, &NRHS, V, &LDA, IPIV, cj, &LDB, &INFO);
    if(INFO != 0) MayDay::Abort("sisdc::setup_qmj - could not compute weights");
    
    // Now construct qmj
    for (int m = 0; m < a_p; m++){
      m_qmj[m][j] = 0.0;
      for (int k = 0; k < nnodes; k++){
	m_qmj[m][j] += cj[k]*(pow(m_nodes[m+1], k+1) - pow(m_nodes[m], k+1))/(k+1);
      }
    }
  }
}

void sisdc::setup_subintervals(const Real a_time, const Real a_dt){
  CH_TIME("sisdc::setup_subintervals");
  if(m_verbosity > 5){
    pout() << "sisdc::setup_subintervals" << endl;
  }

  // m_nodes are Gauss-Lobatto nodes on [-1,1]. These must
  // be shifted to [t_n,t_n + a_dt]
  m_tm.resize(m_nodes.size());
  Vector<Real> shifted_nodes = m_nodes;
  for (int m = 0; m < shifted_nodes.size(); m++){
    shifted_nodes[m] += 1.0;    // [0,2]
    shifted_nodes[m] *= 0.5;    // [0,1]
    shifted_nodes[m] *= a_dt;   // [0, a_dt]
    shifted_nodes[m] += a_time; // [a_time, a_time + a_dt]

    m_tm[m] = shifted_nodes[m];
  }

  // dtm = t_{m+1} - t_m. Order 1 is special since we only use the SISDC predictor from a second order formulation
  m_dtm.resize(m_tm.size() - 1);
  for (int m = 0; m < m_tm.size()-1; m++){
    m_dtm[m] = m_tm[m+1] - m_tm[m];
  }
}

void sisdc::quad(EBAMRCellData& a_quad, const Vector<EBAMRCellData>& a_integrand, const int a_m){
  CH_TIME("sisdc::quad");
  if(m_verbosity > 5){
    pout() << "sisdc::quad" << endl;
  }

  if(a_m < 0)     MayDay::Abort("sisdc::quad - bad index a_m < 0");
  if(a_m >= m_p)  MayDay::Abort("sisdc::quad - bad index a_m >= m_p");

  data_ops::set_value(a_quad, 0.0);
  for (int j = 0; j <= m_p; j++){
    data_ops::incr(a_quad, a_integrand[j], m_qmj[a_m][j]);
  }
}

void sisdc::quad(EBAMRIVData& a_quad, const Vector<EBAMRIVData>& a_integrand, const int a_m){
  CH_TIME("sisdc::quad");
  if(m_verbosity > 5){
    pout() << "sisdc::quad" << endl;
  }

  if(a_m < 0)     MayDay::Abort("sisdc::quad - bad index a_m < 0");
  if(a_m >= m_p)  MayDay::Abort("sisdc::quad - bad index a_m >= m_p");

  data_ops::set_value(a_quad, 0.0);
  for (int j = 0; j <= m_p; j++){
    data_ops::incr(a_quad, a_integrand[j], m_qmj[a_m][j]);
  }
}
  
void sisdc::copy_phi_p_to_cdr(){
  CH_TIME("sisdc::copy_phi_p_to_cdr");
  if(m_verbosity > 5){
    pout() << "sisdc::copy_phi_p_to_cdr" << endl;
  }

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>&  solver  = solver_it();
    RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);

    EBAMRCellData& phi = solver->get_state();
    const EBAMRCellData& phip = storage->get_phi()[m_p];
    data_ops::copy(phi, phip);
  }
}

void sisdc::copy_sigma_p_to_sigma(){
  CH_TIME("sisdc::copy_sigma_p_to_sigma");
  if(m_verbosity > 5){
    pout() << "sisdc::copy_sigma_p_to_sigma" << endl;
  }

  EBAMRIVData& sigma        = m_sigma->get_state();
  const EBAMRIVData& sigmap = m_sigma_scratch->get_sigma()[m_p];
  data_ops::copy(sigma, sigmap);
}

Real sisdc::advance(const Real a_dt){
  CH_TIME("sisdc::advance");
  if(m_verbosity > 2){
    pout() << "sisdc::advance" << endl;
  }

  // ---------------------------------------------------------------------------------------------------
  // TLDR:  When we enter this routine, solvers SHOULD have been filled with valid ready and be ready 
  //        advancement. If you think that this may not be the case, activate the debugging below
  // ---------------------------------------------------------------------------------------------------

  // Initialize integrations. If we do corrections, we need FD(phi_0) since this is implicit
  sisdc::copy_cdr_to_phi_m0();
  sisdc::copy_sigma_to_sigma_m0();
  if(m_k > 0) sisdc::compute_FD_0();
  if(m_k > 0 && m_adaptive_dt) {
    sisdc::store_solvers();
  }

  // SISDC advance
  Real first_dt       = a_dt;
  Real actual_dt      = a_dt;
  int num_reject      = 0;
  int num_corrections = 0;
  bool accept_step    = false;
  bool retry_step     = true;

  m_max_error = 0.1234E5;
  Real t = 0.0;
  while(!accept_step && retry_step){
    num_corrections = 0;
    sisdc::setup_subintervals(m_time, actual_dt);
    sisdc::integrate(a_dt, m_time, false);
    for(int icorr = 0; icorr < Max(m_k, m_min_corr); icorr++){
      num_corrections++;

      sisdc::initialize_errors();
      sisdc::reconcile_integrands();
      sisdc::integrate(a_dt, m_time, true);
      sisdc::finalize_errors();
      if(m_max_error < m_err_thresh && m_adaptive_dt && icorr >= m_min_corr) break; // No need in going beyond
    }

    // Compute a new time step. If it is smaller than the minimum allowed CFL step, accept the step anyways
    if(m_adaptive_dt){
      sisdc::compute_new_dt(accept_step, actual_dt, num_corrections);
      
      if(!accept_step){  // Step rejection, use the new dt for next step.
#if 0 // Debug
	if(procID() == 0){
	  std::cout << "Rejecting step = " << m_step
		    << "\t dt = " << actual_dt
		    << "\t Err = " << m_max_error
		    << "\t New dt = " << m_new_dt
		    << std::endl;
	}
#endif
	actual_dt = m_new_dt;
	num_reject++;

	retry_step  = num_reject <= m_max_retries;
	
	if(retry_step){
	  sisdc::restore_solvers();
	  sisdc::compute_E_into_scratch();
	  sisdc::compute_cdr_gradients();
	  sisdc::compute_cdr_velo(m_time);
	  time_stepper::compute_cdr_diffusion(m_poisson_scratch->get_E_cell(), m_poisson_scratch->get_E_eb());
	  sisdc::compute_cdr_sources(m_time);
	}
      }
    }
    else{
      m_new_dt = 1.234567E89;
      accept_step = true;
    }
  }


  // Copy results back to solvers, and update the Poisson and radiative transfer equations
  sisdc::copy_phi_p_to_cdr();
  sisdc::copy_sigma_p_to_sigma();

  const Real t0 = MPI_Wtime();
  sisdc::update_poisson();
  sisdc::update_stationary_rte(m_time + actual_dt); // Only triggers if m_rte->is_stationar() == true

  // Always recompute source terms and velocities for the next time step. These were computed ea
  const Real t1 = MPI_Wtime();
  sisdc::compute_cdr_gradients();
  const Real t2 = MPI_Wtime();
  sisdc::compute_cdr_velo(m_time + actual_dt);
  const Real t3 = MPI_Wtime();
  time_stepper::compute_cdr_diffusion(m_poisson_scratch->get_E_cell(), m_poisson_scratch->get_E_eb());
  const Real t4 = MPI_Wtime();

  // In case we're using FHD, we need to tell the kinetics module about the time step before computign sources
  Real next_dt;
  time_code::which_code dummy;
  compute_dt(next_dt, dummy);
  m_plaskin->set_dt(next_dt);
  sisdc::compute_cdr_sources(m_time + actual_dt);
  if(!m_rte->is_stationary()){
    
  }
  const Real t5 = MPI_Wtime();

  // Profile step
  if(m_print_report)  sisdc::adaptive_report(first_dt, actual_dt, m_new_dt, num_corrections, num_reject, m_max_error);
  if(m_profile_steps) sisdc::write_step_profile(actual_dt, m_max_error, m_p, num_corrections, num_reject);

#if 0 // Debug
  pout() << "poisson time = " << t1 - t0 << endl;
  pout() << "gradient = " << t2 - t1 << endl;
  pout() << "velo = " << t3 - t2 << endl;
  pout() << "diffco = " << t4 - t3 << endl;
  pout() << "source = " << t5 - t4 << endl;
#endif

  // Store current error. 
  m_have_err  = true;
  m_pre_error = m_max_error;

  //
  m_accum_cfl += actual_dt/m_dt_cfl;
  
  return actual_dt;
}

void sisdc::copy_cdr_to_phi_m0(){
  CH_TIME("sisdc::copy_cdr_to_phi_m0");
  if(m_verbosity > 5){
    pout() << "sisdc::copy_cdr_to_phi_m0" << endl;
  }

  // CDR solvers
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>&  solver  = solver_it();
    RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);
    
    EBAMRCellData& phi0 = storage->get_phi()[0];
    const EBAMRCellData& phi = solver->get_state();
    data_ops::copy(phi0, phi);
  }
}

void sisdc::copy_sigma_to_sigma_m0(){
  CH_TIME("sisdc::copy_sigma_sigma_m0");
  if(m_verbosity > 5){
    pout() << "sisdc::copy_sigma_to_sigma_m0" << endl;
  }

  // Copy sigma to starting state
  EBAMRIVData& sigma0      = m_sigma_scratch->get_sigma()[0];
  const EBAMRIVData& sigma = m_sigma->get_state();
  data_ops::copy(sigma0, sigma);
}

void sisdc::compute_FD_0(){
  CH_TIME("sisdc::compute_FD_0");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_FD_0" << endl;
  }

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);
    
    const EBAMRCellData& phi_0 = storage->get_phi()[0]; // phi_0
    EBAMRCellData& FD_0        = storage->get_FD()[0];  // FD(phi_0)
    
    if(solver->is_diffusive()){
      cdr_tga* tgasolver = (cdr_tga*) (&(*solver));
      tgasolver->compute_divD(FD_0, phi_0);

      // Shouldn't be necesary
      // m_amr->average_down(FD_0, m_cdr->get_phase());
      // m_amr->interp_ghost(FD_0, m_cdr->get_phase());
    }
    else{
      data_ops::set_value(FD_0, 0.0);
    }
  }
}

void sisdc::integrate(const Real a_dt, const Real a_time, const bool a_corrector){
  CH_TIME("sisdc::integrate");
  if(m_verbosity > 5){
    pout() << "sisdc::integrate" << endl;
  }

  // 1. The first time we enter this routine, source terms and velocities were updated.
  // 2. For further calls, source terms and velocities have been overwritten, but since the explicit
  //    operator slopes do not change, this is perfectly fine. We just increment with the lagged terms. 

  // Always update boundary conditions on the way in. All of these calls use the stuff that reside in the solvers,
  // which is what we need to do at the start of the time step.
  Real t0, t1;
  Real total_time = 0.0;
  Real setup_time  = 0.0;
  Real advect_time = 0.0;
  Real diffusive_time = 0.0;

  t0 = MPI_Wtime();
  sisdc::compute_E_into_scratch();
  sisdc::compute_cdr_eb_states();
  sisdc::compute_cdr_fluxes(a_time);
  sisdc::compute_cdr_domain_states();
  sisdc::compute_cdr_domain_fluxes(a_time);
  sisdc::compute_sigma_flux();
  t1 = MPI_Wtime();

  total_time = -t0;
  setup_time = t1-t0;

  // We begin with phi[0] = phi(t_n). Then update phi[m+1].
  for(int m = 0; m < m_p; m++){
    
    // This does the transient rte advance. Stationary solves are done after computing the E-field.
    // The transient solve needs to happen BEFORE the reaction solve in case there is a tight coupling
    // between the RTE and CDR equations (relaxations of excited states). The integrate_rte routine compues
    // source terms just before advancing, and since source terms for the CDR equations are not updated
    // between this routine and the integrate_advection_reaction call, we are ensured that the source terms
    // are consistent. 
    t0 = MPI_Wtime();
    if(m_consistent_rte) {
      sisdc::integrate_rte(a_dt, m, a_corrector);
    }

    // This computes phi_(m+1) = phi_m + dtm*FAR_m(phi_m) + lagged quadrature and lagged advection-reaction
    t0 = MPI_Wtime();
    sisdc::integrate_advection_reaction(a_dt, m, a_corrector);
    t1 = MPI_Wtime();
    advect_time += t1-t0;

    // This does the diffusion advance. It also adds in the remaining lagged diffusion terms before the implicit diffusion solve
    t0 = MPI_Wtime();
    sisdc::integrate_diffusion(a_dt, m, a_corrector);
    t1 = MPI_Wtime();
    diffusive_time += t1-t0;



    // After the diffusion step we should update source terms and boundary conditions for the next step. We don't
    // do this on the last step. This is done either in the reconcile_integrands routine, or after SISDC is done
    // with its substebs. 
    const bool last = (m == m_p-1);
    if(!last){
      Vector<EBAMRCellData*> cdr_densities_mp1 = sisdc::get_cdr_phik(m+1);
      EBAMRIVData& sigma_mp1 = sisdc::get_sigmak(m+1);
      const Real t_mp1 = m_tm[m+1];

      // Update electric field, RTE equations, source terms, and velocities. 
      if(m_consistent_E)   sisdc::update_poisson(cdr_densities_mp1, sigma_mp1);
      if(m_consistent_rte) sisdc::update_stationary_rte(cdr_densities_mp1, t_mp1);
      if(m_compute_S)      sisdc::compute_cdr_gradients(cdr_densities_mp1);
      if(m_compute_v)      sisdc::compute_cdr_velo(cdr_densities_mp1, t_mp1);
      if(m_compute_S)      sisdc::compute_cdr_sources(cdr_densities_mp1, t_mp1);
      if(m_compute_D)      sisdc::update_diffusion_coefficients();

      // Update boundary conditions for cdr and sigma equations. 
      sisdc::compute_cdr_eb_states(cdr_densities_mp1);
      sisdc::compute_cdr_fluxes(cdr_densities_mp1, t_mp1);
      sisdc::compute_cdr_domain_states(cdr_densities_mp1);
      sisdc::compute_cdr_domain_fluxes(cdr_densities_mp1, t_mp1);
      sisdc::compute_sigma_flux();
    }
  }
  t1 = MPI_Wtime();

  total_time += t1;

#if 0
  pout() << endl
	 << "setup time = " << setup_time << endl
	 << "advect_time = " << advect_time << endl
	 << "diffusive_time = " << diffusive_time << endl
	 << "total time = " << total_time << endl
	 << endl;
#endif
}

void sisdc::integrate_advection_reaction(const Real a_dt, const int a_m, const bool a_corrector){
  CH_TIME("sisdc::integrate_advection_reaction");
  if(m_verbosity > 5){
    pout() << "sisdc::integrate_advection_reaction" << endl;
  }

  // Advance phi_(m+1) = phi_m + dtm*F_A using either subcyling or not. These routines do nothing
  // with the operator slopes for phi, but they do adjust the slopes m_Fsig (but not m_Fsum) for sigma. Incidentally,
  // if m=0 and a_corrector=true, we can increment directly with the precomputed advection-reaction. This means that
  // we can skip the advective advance. The sigma advance is accordingly also skipped.
  const bool skip = (a_m == 0 && a_corrector);
  const Real t0 = MPI_Wtime();
  if(!skip){
    if(m_subcycle){
      if(m_multistep){
	sisdc::integrate_advection_multistep(a_dt, a_m, a_corrector);
      }
      else{
	sisdc::integrate_advection_subcycle(a_dt, a_m, a_corrector);
      }
    }
    else{
      sisdc::integrate_advection_nosubcycle(a_dt, a_m, a_corrector);
    }
  }
  const Real t1 = MPI_Wtime();

  // Add in the reaction term and then compute the operator slopes.
  // If this is the corrector and m=0, we skipped the advection advance because we can use the precomputed
  // advection-reaction operator slope. In this case phi_(m+1) is bogus and we need to recompute it. Otherwise,
  // phi_(m+1) = phi_m + dtm*FA_m, and we just increment with the reaction operator. 
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);

    // phi_(m+1) = phi_M
    EBAMRCellData& phi_m1      = storage->get_phi()[a_m+1];
    EBAMRCellData& scratch     = storage->get_scratch();
    const EBAMRCellData& phi_m = storage->get_phi()[a_m];
    
    // Increment with operator slopes. m=0 and corrector is a special case where we skipped the advective advance,
    // choosing instead to use the old slopes (which did not change)
    if(skip){ // Can use the old slopes
      const EBAMRCellData& FAR_m = storage->get_FAR()[a_m]; // Slope, doesn't require recomputation. 
      data_ops::copy(phi_m1, phi_m);
      data_ops::incr(phi_m1, FAR_m, m_dtm[a_m]);
      if(a_corrector) {
	data_ops::copy(scratch, FAR_m);
      }
    }
    else{ // If we made it here, phi_(m+1) = phi_m + dtm*FA(phi_m) through the integrate_advection_subcycle routine
      EBAMRCellData& FAR_m     = storage->get_FAR()[a_m]; // Currently the old slope
      EBAMRCellData& src = solver->get_source();    // Updated source

      // Increment swith source and then compute slope. This has already been done 
      if(!(m_cycle_sources && m_subcycle)){
	data_ops::incr(phi_m1, src, m_dtm[a_m]);  // phi_(m+1) = phi_m + dtm*(FA_m + FR_m)
      }

      // This shouldn't be necessary
      // m_amr->average_down(phi_m1, m_cdr->get_phase());
      // m_amr->interp_ghost(phi_m1, m_cdr->get_phase());

      if(a_corrector){ // Back up the old slope first, we will need it for the lagged term
	data_ops::copy(scratch, FAR_m);
      }
      data_ops::copy(FAR_m, phi_m1);            // FAR_m = (phi_(m+1) - phi_m)/dtm
      data_ops::incr(FAR_m, phi_m, -1.0);       // :
      data_ops::scale(FAR_m, 1./m_dtm[a_m]);    // :

      // Shouldn't be necessary
      // m_amr->average_down(FAR_m, m_cdr->get_phase());
      // m_amr->interp_ghost(FAR_m, m_cdr->get_phase());
    }

    // Now add in the lagged advection-reaction and quadrature terms. This is a bit weird, but we did overwrite
    // FAR_m above after the advection-reaction advance, but we also backed up the old term into scratch. 
    if(a_corrector){
      data_ops::incr(phi_m1, scratch, -m_dtm[a_m]); // phi_(m+1)^(k+1) = phi_m^(k+1) + dtm*(FAR_m^(k+1) - FAR_m^k)
      sisdc::quad(scratch, storage->get_F(), a_m);  // Does the quadrature of the lagged operator slopes. 
      data_ops::incr(phi_m1, scratch, 0.5*a_dt);    // phi_(m+1)^(k+1) = phi_m^(k+1) + dtm*(FAR_m^(k+1) - FAR_m^k) + I_m^(m+1)
    }
  }
  const Real t2 = MPI_Wtime();

  // Add in the lagged terms for sigma. As above, m=0 and corrector is a special case where we just use the old slopes.
  EBAMRIVData& sigma_m1      = m_sigma_scratch->get_sigma()[a_m+1];
  const EBAMRIVData& sigma_m = m_sigma_scratch->get_sigma()[a_m];
  if(skip){
    const EBAMRIVData& Fsig_m = m_sigma_scratch->get_Fold()[a_m]; // Here, we should be able to use either Fold or Fnew
    data_ops::copy(sigma_m1, sigma_m);                            // since Fsig_0 is only computed once. 
    data_ops::incr(sigma_m1, Fsig_m, m_dtm[a_m]);
  }

  const Real t3 = MPI_Wtime();
  if(a_corrector){ // Add in the lagged terms. When we make it here, sigma_(m+1) = sigma_m + dtm*Fsig_m. 
    EBAMRIVData& Fsig_lag = m_sigma_scratch->get_Fold()[a_m];
    data_ops::incr(sigma_m1, Fsig_lag, -m_dtm[a_m]);

    // Add in the quadrature term
    EBAMRIVData& scratch = m_sigma_scratch->get_scratch();
    sisdc::quad(scratch, m_sigma_scratch->get_Fold(), a_m);
    data_ops::incr(sigma_m1, scratch, 0.5*a_dt); // Mult by 0.5*a_dt due to scaling on [-1,1] for quadrature
  }
  const Real t4 = MPI_Wtime();

#if 0
  pout() << "integrate_advection_reaction::" << endl;
  pout() << "t1-t0 = " << t1-t0 << endl;
  pout() << "t2-t1 = " << t2-t2 << endl;
  pout() << "t3-t2 = " << t3-t2 << endl;
  pout() << "t4-t3 = " << t4-t3 << endl;
#endif
}

void sisdc::integrate_advection_nosubcycle(const Real a_dt, const int a_m, const bool a_corrector){
  CH_TIME("sisdc::integrate_advection_nosubcycle");
  if(m_verbosity > 5){
    pout() << "sisdc::integrate_advection_nosubcycle" << endl;
  }

  // TLDR; This routine should do phi_(m+1) = phi_m + dtm*FA_m, and sigma_(m+1) = sigma_m + dt*Fsig_m.
  //       It also computes the sigma slope.
  //
  //       The lagged terms are not a part of this routine. 

  if(a_m == 0 && a_corrector){
    MayDay::Abort("sisdc::integrate_advection_nosubcycle - (m==0 && corrector==true) should never happen");
  }

  // Advance cdr equations
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);

    EBAMRCellData& phi_m1      = storage->get_phi()[a_m+1];
    EBAMRCellData& scratch     = storage->get_scratch();
    const EBAMRCellData& phi_m = storage->get_phi()[a_m];

    if(solver->is_mobile()){
      const Real extrap_dt = m_extrap_advect ? 2.0*m_extrap_dt*m_dtm[a_m] : 0.0; // Factor of 2 due to EBPatchAdvect
      solver->compute_divF(scratch, phi_m, extrap_dt, true);                     // scratch =  Div(v_m*phi_m^(k+1))
    
      data_ops::copy(phi_m1, phi_m);
      data_ops::incr(phi_m1, scratch, -m_dtm[a_m]);
      // m_amr->average_down(phi_m1, m_cdr->get_phase());
      // m_amr->interp_ghost(phi_m1, m_cdr->get_phase());
      data_ops::floor(phi_m1, 0.0);
    }
    else{
      data_ops::copy(phi_m1, phi_m);
    }
  }

  // Update sigma. Also compute the new slope.
  EBAMRIVData& sigma_m1      = m_sigma_scratch->get_sigma()[a_m+1];
  EBAMRIVData& Fsig_new      = m_sigma_scratch->get_Fnew()[a_m];
  const EBAMRIVData& sigma_m = m_sigma_scratch->get_sigma()[a_m];
  m_sigma->compute_rhs(Fsig_new); // Fills Fsig_new with BCs from CDR solvers
  data_ops::copy(sigma_m1, sigma_m);
  data_ops::incr(sigma_m1, Fsig_new, m_dtm[a_m]);
}

void sisdc::integrate_advection_multistep(const Real a_dt, const int a_m, const bool a_corrector){
  CH_TIME("sisdc::integrate_advection_nosubcycle");
  if(m_verbosity > 5){
    pout() << "sisdc::integrate_advection_nosubcycle" << endl;
  }

  // TLDR; This routine should do phi_(m+1) = phi_m + dtm*FA_m, and sigma_(m+1) = sigma_m + dt*Fsig_m.
  //       It also computes the sigma slope.
  //
  //       The lagged terms are not a part of this routine. 

  if(a_m == 0 && a_corrector){
    MayDay::Abort("sisdc::integrate_advection_multistep - (m==0 && corrector==true) should never happen");
  }

  // Copy into phi[a_m+1]
  sisdc::subcycle_copy_states(a_m);

  const int nsteps = ceil(m_dtm[a_m]/(m_cycleCFL*m_dt_cfl));
  const Real dt    = m_dtm[a_m]/nsteps;

#if 0
  MayDay::Abort("sisdc::integrate_advection_multistep - multistep is not done");
#endif
  
  // Advance cdr equations. Use a Heun's method
  for (int istep = 0; istep < nsteps; istep++){

    // Always need to update boundary conditions and source terms (if we substep with the sources)
    if (istep > 0){
      Vector<EBAMRCellData*> cdr_densities_mp1 = sisdc::get_cdr_phik(a_m+1);
      EBAMRIVData& sigma_mp1 = sisdc::get_sigmak(a_m+1);
      const Real t_mp1 = m_tm[a_m+1];

      if(m_cycle_sources){
	if(m_compute_S)      sisdc::compute_cdr_gradients(cdr_densities_mp1);
	if(m_compute_S)      sisdc::compute_cdr_sources(cdr_densities_mp1, t_mp1);
      }
      
      sisdc::compute_cdr_eb_states(cdr_densities_mp1);
      sisdc::compute_cdr_fluxes(cdr_densities_mp1, t_mp1);
      sisdc::compute_cdr_domain_states(cdr_densities_mp1);
      sisdc::compute_cdr_domain_fluxes(cdr_densities_mp1, t_mp1);
      sisdc::compute_sigma_flux();
    }

    // First Heun method stage
    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);

      EBAMRCellData& phi_m1      = storage->get_phi()[a_m+1];
      EBAMRCellData& scratch     = storage->get_scratch();
      EBAMRCellData& scratch2    = storage->get_scratch2();
      const EBAMRCellData& src   = solver->get_source();

      const Real extrap_dt = m_extrap_advect ? 2.0*m_extrap_dt*dt : 0.0; // Factor of 2 due to EBPatchAdvect
      solver->compute_divF(scratch, phi_m1, extrap_dt, true);            // scratch =  Div(v_m*phi_m^(k+1))

      data_ops::copy(scratch2, phi_m1);
      data_ops::incr(phi_m1, scratch, -dt);
      if(m_cycle_sources){
	data_ops::incr(phi_m1, src, dt);
      }
      m_amr->average_down(phi_m1, m_cdr->get_phase());
      m_amr->interp_ghost(phi_m1, m_cdr->get_phase());
      data_ops::floor(phi_m1, 0.0);
    }

    // Update sigma. Also compute the new slope.
#if 0 // Bogus code, please, please fix this!
  EBAMRIVData& sigma_m1      = m_sigma_scratch->get_sigma()[a_m+1];
  EBAMRIVData& Fsig_new      = m_sigma_scratch->get_Fnew()[a_m];
  const EBAMRIVData& sigma_m = m_sigma_scratch->get_sigma()[a_m];
  m_sigma->compute_rhs(Fsig_new); // Fills Fsig_new with BCs from CDR solvers
  data_ops::copy(sigma_m1, sigma_m);
  data_ops::incr(sigma_m1, Fsig_new, m_dtm[a_m]);
#endif

    // Always update BC between stages
    Vector<EBAMRCellData*> cdr_densities_mp1 = sisdc::get_cdr_phik(a_m+1);
    EBAMRIVData& sigma_mp1 = sisdc::get_sigmak(a_m+1);
    const Real t_mp1 = m_tm[a_m+1];

    if(m_cycle_sources){
      if(m_compute_S)      sisdc::compute_cdr_gradients(cdr_densities_mp1);
      if(m_compute_S)      sisdc::compute_cdr_sources(cdr_densities_mp1, t_mp1);
    }
      
    sisdc::compute_cdr_eb_states(cdr_densities_mp1);
    sisdc::compute_cdr_fluxes(cdr_densities_mp1, t_mp1);
    sisdc::compute_cdr_domain_states(cdr_densities_mp1);
    sisdc::compute_cdr_domain_fluxes(cdr_densities_mp1, t_mp1);
    sisdc::compute_sigma_flux();

    // Second Heun method stage
    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);

      EBAMRCellData& phi_m1      = storage->get_phi()[a_m+1];
      EBAMRCellData& scratch     = storage->get_scratch();
      EBAMRCellData& scratch2    = storage->get_scratch2();
      const EBAMRCellData& src   = solver->get_source();

      const Real extrap_dt = m_extrap_advect ? 2.0*m_extrap_dt*dt : 0.0; // Factor of 2 due to EBPatchAdvect
      solver->compute_divF(scratch, phi_m1, extrap_dt, true);            // scratch =  Div(v_m*phi_m^(k+1))
    
      data_ops::incr(phi_m1, scratch, -dt);
      if(m_cycle_sources){
	data_ops::incr(phi_m1, src, dt);
      }
      data_ops::incr(phi_m1, scratch2, 1.0);
      data_ops::scale(phi_m1, 0.5);
      m_amr->average_down(phi_m1, m_cdr->get_phase());
      m_amr->interp_ghost(phi_m1, m_cdr->get_phase());
      data_ops::floor(phi_m1, 0.0);
    }
  }


}

void sisdc::integrate_diffusion(const Real a_dt, const int a_m, const bool a_corrector){
  CH_TIME("sisdc::integrate_diffusion");
  if(m_verbosity > 5){
    pout() << "sisdc::integrate_diffusion" << endl;
  }

  // TLDR: We're solving
  //
  // phi_(m+1)^(k+1) = phi_(m)^(k+1,\ast) + dtm*FD_(m+1)^(k+1) + sources. 
  //
  // This routine does not modify FD_(m+1)^k. This is replaced by FD_(m+1)^(k+1) later on. 
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = sisdc::get_cdr_storage(solver_it);
    
    if(solver->is_diffusive()){
      EBAMRCellData& phi_m1      = storage->get_phi()[a_m+1]; // Advected solution. Possibly with lagged terms. 
      const EBAMRCellData& phi_m = storage->get_phi()[a_m];

      // Build the diffusion source term
      EBAMRCellData& source   = storage->get_scratch();
      EBAMRCellData& init_soln = storage->get_scratch2();
      data_ops::set_value(source, 0.0); // No source term
      
      data_ops::copy(init_soln, phi_m1);      // Copy initial solutions
      if(a_corrector){
	const EBAMRCellData& FD_m1k = storage->get_FD()[a_m+1];      // FD_(m+1)^k. Lagged term.
	data_ops::incr(init_soln, FD_m1k, -m_dtm[a_m]);
      }
      m_amr->average_down(init_soln, m_cdr->get_phase());
      m_amr->interp_ghost(init_soln, m_cdr->get_phase());
#if 1 // Original code
      data_ops::copy(phi_m1, phi_m);
#else // Debug code
      data_ops::copy(phi_m1, init_soln);
#endif

      // Solve
      if(m_use_tga){
	solver->advance_tga(phi_m1, init_soln, source, m_dtm[a_m]); // No source. 
      }
      else{
	solver->advance_euler(phi_m1, init_soln, source, m_dtm[a_m]); // No source. 
      }
      m_amr->average_down(phi_m1, m_cdr->get_phase());
      m_amr->interp_ghost(phi_m1, m_cdr->get_phase());
      data_ops::floor(phi_m1, 0.0);

      // Update the operator slope
      EBAMRCellData& FD_m1k = storage->get_FD()[a_m+1];
      data_ops::set_value(FD_m1k, 0.0);
      data_ops::incr(FD_m1k, phi_m1, 1.0);
      data_ops::incr(FD_m1k, init_soln, -1.0);
      data_ops::scale(FD_m1k, 1./m_dtm[a_m]);

      m_amr->average_down(FD_m1k, m_cdr->get_phase());
      m_amr->interp_ghost(FD_m1k, m_cdr->get_phase());
    }
    else{
      EBAMRCellData& FD_m1k = storage->get_FD()[a_m+1];
      data_ops::set_value(FD_m1k, 0.0);
    }
  }
}

void sisdc::reconcile_integrands(){
  CH_TIME("sisdc::reconcile_integrands");
  if(m_verbosity > 5){
    pout() << "sisdc::reconcile_integrands" << endl;
  }

  Vector<EBAMRCellData*> cdr_densities_p = sisdc::get_cdr_phik(m_p);
  EBAMRIVData& sigma_p = sisdc::get_sigmak(m_p);
  const Real t_p = m_tm[m_p];

//  Update electric field, RTE equations, source terms, and velocities
  if(m_consistent_E)   sisdc::update_poisson(cdr_densities_p, sigma_p);
  if(m_consistent_rte) sisdc::update_stationary_rte(cdr_densities_p, t_p);
  if(m_compute_S)      sisdc::compute_cdr_gradients(cdr_densities_p);
  if(m_compute_v)      sisdc::compute_cdr_velo(cdr_densities_p, t_p);
  if(m_compute_S)      sisdc::compute_cdr_sources(cdr_densities_p, t_p);

  // Update boundary conditions for cdr and sigma equations
  sisdc::compute_cdr_eb_states(cdr_densities_p);
  sisdc::compute_cdr_fluxes(cdr_densities_p, t_p);
  sisdc::compute_cdr_domain_states(cdr_densities_p);
  sisdc::compute_cdr_domain_fluxes(cdr_densities_p, t_p);
  sisdc::compute_sigma_flux();

  // Now compute FAR_p - that wasn't done when we integrated
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = sisdc::get_cdr_storage(solver_it);
    const int idx = solver_it.get_solver();

    // This has not been computed yet. Do it.
    EBAMRCellData& FAR_p       = storage->get_FAR()[m_p];
    const EBAMRCellData& phi_p = *cdr_densities_p[idx] ;
    const EBAMRCellData& src   = solver->get_source();

    if(solver->is_mobile()){
      const Real extrap_dt = m_extrap_advect ? 2.0*m_extrap_dt*m_dtm[m_p-1] : 0.0; // Factor of 2 because of EBPatchAdvect
      solver->compute_divF(FAR_p, phi_p, extrap_dt, true); // FAR_p =  Div(v_p*phi_p)
      data_ops::scale(FAR_p, -1.0);                        // FAR_p = -Div(v_p*phi_p)
    }
    else{
      data_ops::set_value(FAR_p, 0.0);
    }
    data_ops::incr(FAR_p, src, 1.0);                     // RHS = -Div(v_m*phi_m) + S_m = FAR(phi_m)

    // Build the integrand
    for (int m = 0; m <= m_p; m++){
      EBAMRCellData& F_m   = storage->get_F()[m];
      EBAMRCellData& FD_m  = storage->get_FD()[m];
      EBAMRCellData& FAR_m = storage->get_FAR()[m];

      data_ops::copy(F_m, FAR_m);
      if(solver->is_diffusive()){
	data_ops::incr(F_m, FD_m, 1.0);
      }

      // Shouldn't be necessary
      // m_amr->average_down(F_m, m_cdr->get_phase());
      // m_amr->interp_ghost(F_m, m_cdr->get_phase());
    }
  }

  // Compute Fsig_p - that wasn't done either
  EBAMRIVData& Fnew_p = m_sigma_scratch->get_Fnew()[m_p];
  m_sigma->compute_rhs(Fnew_p);
  for (int m = 0; m <= m_p; m++){
    EBAMRIVData& Fold_m = m_sigma_scratch->get_Fold()[m];
    EBAMRIVData& Fnew_m = m_sigma_scratch->get_Fnew()[m];
    data_ops::copy(Fold_m, Fnew_m);
  }
}

void sisdc::initialize_errors(){
  CH_TIME("sisdc::corrector_initialize_errors");
  if(m_verbosity > 5){
    pout() << "sisdc::corrector_initialize_errors" << endl;
  }

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = sisdc::get_cdr_storage(solver_it);
    const int idx = solver_it.get_solver();
    
    // These should be zero
    if(idx == m_error_idx || m_error_idx < 0){
      EBAMRCellData& error = storage->get_error();
      const EBAMRCellData& phi_final = storage->get_phi()[m_p];

      data_ops::set_value(error, 0.0);
      data_ops::incr(error, phi_final, -1.0);
    }
  }

  EBAMRIVData& error = m_sigma_scratch->get_error();
  const EBAMRIVData& sigma_final = m_sigma_scratch->get_sigma()[m_p];
  data_ops::set_value(error, 0.0);
  data_ops::incr(error, sigma_final, -1.0);
}

void sisdc::finalize_errors(){
  CH_TIME("sisdc::corrector_finalize_errors");
  if(m_verbosity > 5){
    pout() << "sisdc::corrector_finalize_errors" << endl;
  }

  const Real safety = 1.E-20;

  m_max_error = 0.0;
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = sisdc::get_cdr_storage(solver_it);
    const int idx = solver_it.get_solver();

    // Compute error
    if(idx == m_error_idx || m_error_idx < 0){
      EBAMRCellData& error       = storage->get_error();
      const EBAMRCellData& phi_p = storage->get_phi()[m_p];
      data_ops::incr(error, phi_p, 1.0);

      // Compute norms. Only coarsest level
      Real Lerr, Lphi;
      const int lvl = 0;
      data_ops::norm(Lerr, *error[lvl], m_amr->get_domains()[lvl], m_error_norm);
      data_ops::norm(Lphi, *phi_p[lvl], m_amr->get_domains()[lvl], m_error_norm);

      if(Lphi > 0.0){
	m_cdr_error[idx] = Lerr/Lphi;

	m_max_error = Max(m_cdr_error[idx], m_max_error);
      }
#if 0 // Debug
      if(procID() == 0){
	std::cout << "Lerr = " << Lerr << "\t Lphi = " << Lphi << "\t Lerr/Lphi = " << Lerr/Lphi << std::endl;
      }
#endif
    }
  }

  // Override if
  if(m_error_idx >= 0){
    m_max_error = m_cdr_error[m_error_idx];
  }

  // Compute the surface charge conservation error
  EBAMRIVData& error = m_sigma_scratch->get_error();
  const EBAMRIVData& sigma_final = m_sigma_scratch->get_sigma()[m_p];
  data_ops::incr(error, sigma_final, 1.0);
  m_sigma_error = 0.0; // I don't think this is ever used...


}

void sisdc::compute_new_dt(bool& a_accept_step, const Real a_dt, const int a_num_corrections){
  CH_TIME("sisdc::compute_new_dt");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_new_dt" << endl;
  }

  // If a_dt was the smallest possible CFL or hardcap time step, we just have to accept it
  const Real max_gl_dist = sisdc::get_max_node_distance();
  Real dt_cfl = 2.0*m_dt_cfl/max_gl_dist; // This is the smallest time step ON THE FINEST LEVEL

  int Nref = 1;
  if(m_subcycle){
    for (int lvl = 0; lvl < m_amr->get_finest_level(); lvl++){
      Nref = Nref*m_amr->get_ref_rat()[lvl];
    }
    dt_cfl = dt_cfl*Nref;
  }

  // Try time step
  const Real rel_err    = (m_safety*m_err_thresh)/m_max_error;
  const Real dt_adapt   = (m_max_error > 0.0) ? a_dt*pow(rel_err, 1.0/(a_num_corrections+1)) : m_max_dt;
  const Real min_dt_cfl = dt_cfl*m_minCFL;
  const Real max_dt_cfl = dt_cfl*m_maxCFL;

  if(m_max_error <= m_err_thresh){ // Always accept, and compute new step
    a_accept_step = true;

    // Do not grow step too fast
    if(rel_err > 1.0){ // rel_err > 1 => dt_adapt > a_dt
      m_new_dt = Min(m_max_growth*a_dt, dt_adapt);
    }
    else{ // rel_err > 1 => dt_adapt < a_dt. This shrinks the error down to the safety factor. 
      m_new_dt = dt_adapt;
    }
    m_new_dt = Max(m_new_dt, m_min_dt);                   // Don't go below hardcap
    m_new_dt = Min(m_new_dt, m_max_dt);                   // Don't go above other hardcap
    if(!m_subcycle){
      m_new_dt = Max(m_new_dt, dt_cfl*m_minCFL);            // Don't drop below minimum CFL
      m_new_dt = Min(m_new_dt, dt_cfl*m_maxCFL);            // Don't go above maximum CFL
    }
  }
  else{
    a_accept_step = false;

    m_new_dt = m_decrease_safe*dt_adapt; // Decrease time step a little bit extra to avoid another rejection
    if(a_dt <= min_dt_cfl/Nref || a_dt < m_min_dt){ // Step already at minimum. Accept it anyways.
      a_accept_step = true;
    }
    
    if(!m_subcycle){
      m_new_dt = Max(m_new_dt, dt_cfl*m_minCFL);            // Don't drop below minimum CFL
      m_new_dt = Min(m_new_dt, dt_cfl*m_maxCFL);            // Don't go above maximum CFL
    }
    m_new_dt = Max(m_new_dt, m_min_dt);                   // Don't go below hardcap
    m_new_dt = Min(m_new_dt, m_max_dt);                   // Don't go above other hardcap
  }

#if 0 // Debug
  if(procID() == 0) std::cout << "accept = " << a_accept_step
			      << " dt = " << a_dt
			      << " new_dt = " << m_new_dt
			      << " fraction = " << m_new_dt/a_dt
			      << std::endl;
#endif

  m_have_dt_err = true;
}

void sisdc::adaptive_report(const Real a_first_dt, const Real a_dt, const Real a_new_dt, const int a_corr, const int a_rej, const Real a_max_err){
  CH_TIME("sisdc::adaptive_report");
  if(m_verbosity > 5){
    pout() << "sisdc::adaptive_report" << endl;
  }

  pout() << "\n";
  pout() << "sisdc::adaptive_report breakdown" << endl;
  pout() << "--------------------------------\n";
  pout() << "\t Try dt       = " << a_first_dt << endl;
  pout() << "\t Advanced dt  = " << a_dt << endl;
  pout() << "\t New dt       = " << a_new_dt << endl;
  pout() << "\t Subintervals = " << m_p << endl;
  pout() << "\t Corrections  = " << a_corr << endl;
  pout() << "\t Rejections   = " << a_rej << endl;
  pout() << "\t Max error    = " << a_max_err << endl;
  pout() << "\n";
}

void sisdc::compute_dt(Real& a_dt, time_code::which_code& a_timecode){
  CH_TIME("sisdc::compute_dt");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_dt" << endl;
  }

  Real dt = 1.E99;

  int Nref = 1;
  for (int lvl = 0; lvl < m_amr->get_finest_level(); lvl++){
    Nref = Nref*m_amr->get_ref_rat()[lvl];
  }
  const Real max_gl_dist = sisdc::get_max_node_distance();
  m_dt_cfl = m_cdr->compute_cfl_dt();

  //  Real dt_cfldt_cfl = (m_subcycle) ? m_dt_cfl*Nref : m_dt_cfl;
  Real dt_cfl = 2.0*m_dt_cfl/max_gl_dist;
  dt_cfl = m_subcycle ? Nref*dt_cfl : dt_cfl;
  
  // Time step selection for non-adaptive stepping
  if(!m_adaptive_dt){
    if(dt_cfl < dt){
      dt = m_cfl*dt_cfl;
      a_timecode = time_code::cfl;
    }
  }
  else{
    Real new_dt;

    // Step should not exceed m_new_dt. Also, it shoul
    if(m_have_dt_err){
      new_dt = m_new_dt;
      if(!m_subcycle){
	new_dt = Max(new_dt, dt_cfl*m_minCFL);
	new_dt = Min(new_dt, dt_cfl*m_maxCFL);
      }
      else{
	//	new_dt = Min(new_dt, dt_cfl*m_maxCFL);
      }
    }
    else{
      if(!m_subcycle){
	new_dt = m_minCFL*dt_cfl;
      }
      else{
	new_dt = 0.5*max_gl_dist*m_maxCFL*dt_cfl/Nref; // Default start step for subcycling advances
      }
    }

    if(new_dt < dt){
      dt = new_dt;
      a_timecode = time_code::error;
    }
  }

  // Truncate
  if(m_subcycle){
    dt = Max(dt, m_min_cycle_cfl*m_dt_cfl);
    dt = Min(dt, m_max_cycle_cfl*m_dt_cfl);
  }

  // EVERYTHING BELOW HERE IS "STANDARD"
  // -----------------------------------
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

#if 0 // Debug
  if(procID() == 0){
    std::cout << "compute_dt = " << a_dt << "\t m_new_dt = " << m_new_dt << std::endl; 
  }
#endif
}

void sisdc::regrid_internals(){
  CH_TIME("sisdc::regrid_internals");
  if(m_verbosity > 5){
    pout() << "sisdc::regrid_internals" << endl;
  }

  m_accum_cfl = 0.0;
  m_cdr_error.resize(m_plaskin->get_num_species());
  
  sisdc::allocate_cdr_storage();
  sisdc::allocate_poisson_storage();
  sisdc::allocate_rte_storage();
  sisdc::allocate_sigma_storage();

  sisdc::setup_quadrature_nodes(m_p);
  sisdc::setup_qmj(m_p);
}

void sisdc::allocate_cdr_storage(){
  const int ncomp       = 1;
  const int num_species = m_plaskin->get_num_species();

  m_cdr_scratch.resize(num_species);
  
  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx] = RefCountedPtr<cdr_storage> (new cdr_storage(m_amr, m_cdr->get_phase(), ncomp));
    m_cdr_scratch[idx]->allocate_storage(m_p);
  }
}

void sisdc::allocate_poisson_storage(){
  const int ncomp = 1;
  m_poisson_scratch = RefCountedPtr<poisson_storage> (new poisson_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_poisson_scratch->allocate_storage(m_p);
}

void sisdc::allocate_rte_storage(){
  const int ncomp       = 1;
  const int num_photons = m_plaskin->get_num_photons();
  m_rte_scratch.resize(num_photons);
  
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx] = RefCountedPtr<rte_storage> (new rte_storage(m_amr, m_rte->get_phase(), ncomp));
    m_rte_scratch[idx]->allocate_storage(m_p);
  }
}

void sisdc::allocate_sigma_storage(){
  const int ncomp = 1;
  m_sigma_scratch = RefCountedPtr<sigma_storage> (new sigma_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_sigma_scratch->allocate_storage(m_p);
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
  
  sisdc::compute_E(E_cell, m_cdr->get_phase(), phi);     // Compute cell-centered field
  sisdc::compute_E(E_face, m_cdr->get_phase(), E_cell);  // Compute face-centered field
  sisdc::compute_E(E_eb,   m_cdr->get_phase(), E_cell);  // EB-centered field

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
    RefCountedPtr<cdr_storage>& storage = sisdc::get_cdr_storage(solver_it);
    EBAMRCellData& grad = storage->get_gradient();
    m_amr->compute_gradient(grad, *a_states[idx]);
    //    m_amr->average_down(grad, m_cdr->get_phase());
    m_amr->interp_ghost(grad, m_cdr->get_phase());
  }
}

void sisdc::compute_cdr_velo(const Real a_time){
  CH_TIME("sisdc::compute_cdr_velo");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_cdr_velo" << endl;
  }

  sisdc::compute_cdr_velo(m_cdr->get_states(), a_time);
}

void sisdc::compute_cdr_velo(const Vector<EBAMRCellData*>& a_states, const Real a_time){
  CH_TIME("sisdc::compute_cdr_velo(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_cdr_velo(Vector<EBAMRCellData*>, Real)" << endl;
  }

  Vector<EBAMRCellData*> velocities = m_cdr->get_velocities();
  sisdc::compute_cdr_velocities(velocities, a_states, m_poisson_scratch->get_E_cell(), a_time);
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
    RefCountedPtr<cdr_storage>& storage = sisdc::get_cdr_storage(solver_it);

    cdr_states.push_back(&(solver->get_state()));
    eb_states.push_back(&(storage->get_eb_state()));
    eb_gradients.push_back(&(storage->get_eb_grad()));
    cdr_gradients.push_back(&(storage->get_gradient())); // Should already have been computed
  }

  // Extrapolate states to the EB and floor them so we cannot get negative values on the boundary. This
  // won't hurt mass conservation because the mass hasn't been injected yet
  sisdc::extrapolate_to_eb(eb_states, m_cdr->get_phase(), cdr_states);
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    data_ops::floor(*eb_states[idx], 0.0);
  }

  // We should already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIVData eb_gradient;
  m_amr->allocate(eb_gradient, m_cdr->get_phase(), SpaceDim);
  for (int i = 0; i < cdr_states.size(); i++){
    sisdc::extrapolate_to_eb(eb_gradient, m_cdr->get_phase(), *cdr_gradients[i]);
    sisdc::project_flux(*eb_gradients[i], eb_gradient);
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
    RefCountedPtr<cdr_storage>& storage = sisdc::get_cdr_storage(solver_it);

    eb_states.push_back(&(storage->get_eb_state()));
    eb_gradients.push_back(&(storage->get_eb_grad()));
    cdr_gradients.push_back(&(storage->get_gradient())); // Should already have been computed
  }

  // Extrapolate states to the EB and floor them so we cannot get negative values on the boundary. This
  // won't hurt mass conservation because the mass hasn't been injected yet
  sisdc::extrapolate_to_eb(eb_states, m_cdr->get_phase(), a_states);
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    data_ops::floor(*eb_states[idx], 0.0);
  }

  // We should already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIVData eb_gradient;
  m_amr->allocate(eb_gradient, m_cdr->get_phase(), SpaceDim);
  for (int i = 0; i < a_states.size(); i++){
    sisdc::extrapolate_to_eb(eb_gradient, m_cdr->get_phase(), *cdr_gradients[i]);
    sisdc::project_flux(*eb_gradients[i], eb_gradient);
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
    RefCountedPtr<cdr_storage>& storage = sisdc::get_cdr_storage(solver_it);

    cdr_states.push_back(&(solver->get_state()));
    domain_states.push_back(&(storage->get_domain_state()));
    domain_gradients.push_back(&(storage->get_domain_grad()));
    cdr_gradients.push_back(&(storage->get_gradient())); // Should already be computed
  }

  // Extrapolate states to the domain faces
  sisdc::extrapolate_to_domain_faces(domain_states, m_cdr->get_phase(), cdr_states);

  // We already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIFData grad;
  m_amr->allocate(grad, m_cdr->get_phase(), SpaceDim);
  for (int i = 0; i < cdr_states.size(); i++){
    sisdc::extrapolate_to_domain_faces(grad, m_cdr->get_phase(), *cdr_gradients[i]);
    sisdc::project_domain(*domain_gradients[i], grad);
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
  time_stepper::compute_extrapolated_fluxes(extrap_cdr_fluxes, a_states, cdr_velocities, m_cdr->get_phase());
  time_stepper::compute_extrapolated_velocities(extrap_cdr_velocities, cdr_velocities, m_cdr->get_phase());

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
  CH_TIME("sisdc::update_poisson(solver)");
  if(m_verbosity > 5){
    pout() << "sisdc::update_poisson(solver)" << endl;
  }
  
  if(m_do_poisson){ // Solve Poisson equation
    if((m_step +1) % m_fast_poisson == 0){
      time_stepper::solve_poisson();
      this->compute_E_into_scratch();
    }
  }
}

void sisdc::update_poisson(const Vector<EBAMRCellData*>& a_densities, const EBAMRIVData& a_sigma){
  CH_TIME("sisdc::update_poisson(full)");
  if(m_verbosity > 5){
    pout() << "sisdc::update_poisson(full)" << endl;
  }
  
  if(m_do_poisson){ // Solve Poisson equation
    if((m_step +1) % m_fast_poisson == 0){
      time_stepper::solve_poisson(m_poisson->get_state(),
				  m_poisson->get_source(),
				  a_densities,
				  a_sigma,
				  centering::cell_center);
      this->compute_E_into_scratch();
    }
  }
}

void sisdc::update_stationary_rte(const Real a_time){
  CH_TIME("sisdc::update_stationary_rte(solver)");
  if(m_verbosity > 5){
    pout() << "sisdc::update_stationary_rte(solver)" << endl;
  }
  
  if(m_do_rte){
    if((m_step + 1) % m_fast_rte == 0){
      if(m_rte->is_stationary()){
	Vector<EBAMRCellData*> rte_states  = m_rte->get_states();
	Vector<EBAMRCellData*> rte_sources = m_rte->get_sources();
	Vector<EBAMRCellData*> cdr_states  = m_cdr->get_states();

	EBAMRCellData& E = m_poisson_scratch->get_E_cell();

	const Real dummy_dt = 0.0;
	this->solve_rte(rte_states, rte_sources, cdr_states, E, a_time, dummy_dt, centering::cell_center);
      }
    }
  }
}

void sisdc::update_stationary_rte(const Vector<EBAMRCellData*>& a_cdr_states, const Real a_time){
  CH_TIME("sisdc::update_stationary_rte(full)");
  if(m_verbosity > 5){
    pout() << "sisdc::update_stationary_rte(full)" << endl;
  }
  
  if(m_do_rte){
    if((m_step + 1) % m_fast_rte == 0){
      if(m_rte->is_stationary()){
	Vector<EBAMRCellData*> rte_states  = m_rte->get_states();
	Vector<EBAMRCellData*> rte_sources = m_rte->get_sources();

	EBAMRCellData& E = m_poisson_scratch->get_E_cell();

	const Real dummy_dt = 0.0;
	this->solve_rte(rte_states, rte_sources, a_cdr_states, E, a_time, dummy_dt, centering::cell_center);
      }
    }
  }
}

void sisdc::integrate_rte(const Real a_dt, const int a_m, const bool a_corrector){
  CH_TIME("sisdc::integrate_rte(full)");
  if(m_verbosity > 5){
    pout() << "sisdc::integrate_rte(full)" << endl;
  }

  if(m_do_rte){
    if((m_step + 1) % m_fast_rte == 0){
      const Real time = m_time + m_dtm[a_m]; // This is the current time
      
      Vector<EBAMRCellData*>  rte_states  = m_rte->get_states();
      Vector<EBAMRCellData*>  rte_sources = m_rte->get_sources();
      Vector<EBAMRCellData*>  cdr_states  = sisdc::get_cdr_phik(a_m);
      EBAMRCellData& E = m_poisson_scratch->get_E_cell();
      this->solve_rte(rte_states, rte_sources, cdr_states, E, time, a_dt, centering::cell_center);
    }
  }
}

void sisdc::update_diffusion_coefficients(){
  CH_TIME("sisdc::update_diffusion_coefficients");
  if(m_verbosity > 5){
    pout() << "sisdc::update_diffusion_coefficients" << endl;
  }
  time_stepper::compute_cdr_diffusion(m_poisson_scratch->get_E_cell(), m_poisson_scratch->get_E_eb());
}

Vector<EBAMRCellData*> sisdc::get_cdr_errors(){
  CH_TIME("sisdc::get_cdr_errors");
  if(m_verbosity > 5){
    pout() << "sisdc::get_cdr_errors" << endl;
  }
  
  Vector<EBAMRCellData*> ret;
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = sisdc::get_cdr_storage(solver_it);
    ret.push_back(&(storage->get_error()));
  }

  return ret;
}

Vector<EBAMRCellData*> sisdc::get_cdr_phik(const int a_m){
  CH_TIME("sisdc::get_cdr_phik");
  if(m_verbosity > 5){
    pout() << "sisdc::get_cdr_phik" << endl;
  }
  
  Vector<EBAMRCellData*> ret;
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = sisdc::get_cdr_storage(solver_it);
    ret.push_back(&(storage->get_phi()[a_m]));
  }

  return ret;
}

EBAMRIVData& sisdc::get_sigmak(const int a_m){
  CH_TIME("sisdc::get_sigmak");
  if(m_verbosity > 5){
    pout() << "sisdc::get_sigmak)" << endl;
  }
  return m_sigma_scratch->get_sigma()[a_m];
}

void sisdc::write_step_profile(const Real a_dt,
			       const Real a_error,
			       const int  a_substeps,
			       const int  a_corrections,
			       const int  a_rejections){
  CH_TIME("sissdc::write_step_profile");
  if(m_verbosity > 5){
    pout() << "sisdc::write_step_profile" << endl;
  }

  if(procID() == 0 ){

    const std::string fname("sisdc_step_profile.txt");
    
    bool write_header;
    { // Write header if we must
      std::ifstream infile(fname);
      write_header = infile.peek() == std::ifstream::traits_type::eof() ? true : false;
    }

    // Write output
    std::ofstream f;
    f.open(fname, std::ios_base::app);
    const int width = 12;

    if(write_header){
      f << std::left << std::setw(width) << "# Step" << "\t"
	<< std::left << std::setw(width) << "Time" << "\t"
	<< std::left << std::setw(width) << "dt" << "\t"
	<< std::left << std::setw(width) << "Substeps" << "\t"
	<< std::left << std::setw(width) << "Corrections" << "\t"
	<< std::left << std::setw(width) << "Rejections" << "\t"
	<< std::left << std::setw(width) << "Error" << "\t"
	<< std::left << std::setw(width) << "CFL" << "\t"
	<< endl;
    }

    f << std::left << std::setw(width) << m_step << "\t"
      << std::left << std::setw(width) << m_time << "\t"            
      << std::left << std::setw(width) << a_dt << "\t"        
      << std::left << std::setw(width) << a_substeps << "\t"
      << std::left << std::setw(width) << a_corrections << "\t"
      << std::left << std::setw(width) << a_rejections << "\t"
      << std::left << std::setw(width) << m_max_error << "\t"
      << std::left << std::setw(width) << m_dt_cfl << "\t"
      << endl;
  }
}

void sisdc::reset_finer_flux_registers_level(const int a_lvl,
					     const int a_coarsest_level,
					     const int a_finest_level){
  CH_TIME("sisdc::reset_flux_registers_level");
  if(m_verbosity > 5){
    pout() << "sisdc::reset_flux_registers_level" << endl;
  }

  const phase::which_phase phase = m_cdr->get_phase();
  const bool has_fine = a_lvl < a_finest_level;
  if(has_fine){
    EBFluxRegister* fluxreg_fine = m_amr->get_flux_reg(phase)[a_lvl];
    fluxreg_fine->setToZero();
  }
}

void sisdc::reset_redist_registers_level(const int a_lvl, const int a_coarsest_level, const int a_finest_level){
  CH_TIME("sisdc::reset_redist_registers_level");
  if(m_verbosity > 5){
    pout() << "sisdc::reset_redist_registers_level" << endl;
  }

  const phase::which_phase phase = m_cdr->get_phase();
  EBLevelRedist& level_redist = *(m_amr->get_level_redist(phase)[a_lvl]);
  level_redist.setToZero();

  if(m_amr->get_ebcf()){
    const bool has_fine = a_lvl < a_finest_level;
    const bool has_coar = a_lvl > a_coarsest_level;

    RefCountedPtr<EBFineToCoarRedist>& fine2coar_redist = m_amr->get_fine_to_coar_redist(m_cdr->get_phase())[a_lvl];
    RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->get_coar_to_fine_redist(m_cdr->get_phase())[a_lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->get_coar_to_coar_redist(m_cdr->get_phase())[a_lvl];

    if(has_coar){
      fine2coar_redist->setToZero();
    }

    if(has_fine){
      coar2fine_redist->setToZero();
      coar2coar_redist->setToZero();
    }
  }
}

void sisdc::update_flux_registers(LevelData<EBFluxFAB>& a_flux,
				  const int             a_solver,
				  const int             a_lvl,
				  const int             a_coarsest_level,
				  const int             a_finest_level,
				  const Real            a_dt){

  // Increment the coarser flux register and initialize the finer flux register. a_flux holds phi*vel which we can use
  const bool has_fine = a_lvl < a_finest_level;
  const bool has_coar = a_lvl > a_coarsest_level;

  const phase::which_phase phase = m_cdr->get_phase();
  const Interval interv(a_solver, a_solver);
  EBFluxRegister* fluxreg_fine = NULL;
  EBFluxRegister* fluxreg_coar = NULL;

  // Remember, register on a_lvl holds flux between level a_lvl and a_lvl+1
  if(has_fine) fluxreg_fine = m_amr->get_flux_reg(phase)[a_lvl];   
  if(has_coar) fluxreg_coar = m_amr->get_flux_reg(phase)[a_lvl-1]; 

  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];
  const EBISLayout& ebisl      = m_amr->get_ebisl(phase)[a_lvl];

  // This is a bit stupid but in order to get the correct flux in the correct place in the register, we need some storage
  // with enough variables
  EBFluxFactory fact(ebisl);
  LevelData<EBFluxFAB> scratch(dbl, m_plaskin->get_num_species(), a_flux.ghostVect(), fact);
  EBLevelDataOps::setVal(scratch, 0.0);
  a_flux.localCopyTo(Interval(0,0), scratch, interv);
  scratch.exchange(interv);

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const EBISBox& ebisbox = ebisl[dit()];
    const Box box          = dbl.get(dit());
      
    for (int dir = 0; dir < SpaceDim; dir++){

      // Initialize finer flux register with flux leaving this level and into the finer level (or invalid region of that level)
      if(has_fine) {
	for (SideIterator sit; sit.ok(); ++sit){
	  //fluxreg_fine->incrementCoarseBoth(a_flux[dit()][dir], a_dt, dit(), interv, dir, sit());
	  fluxreg_fine->incrementCoarseBoth(scratch[dit()][dir], a_dt, dit(), interv, dir, sit());
	}
      }

      // Increment coarser flux register with flux entering that level from this level 
      if(has_coar){ // The coarser level has already been initialized with the coarse side flux
	for (SideIterator sit; sit.ok(); ++sit){
	  //fluxreg_coar->incrementFineBoth(a_flux[dit()][dir], a_dt, dit(), interv, dir, sit());
	  fluxreg_coar->incrementFineBoth(scratch[dit()][dir], a_dt, dit(), interv, dir, sit());
	}
      }
    }
  }
}

void sisdc::update_redist_register(const LevelData<BaseIVFAB<Real> >& a_mass_diff, const int a_solver, const int a_lvl){
  CH_TIME("cdr_solver::update_redist_register");
  if(m_verbosity > 5){
    pout() << "sisdc::::update_redist_register" << endl;
  }

  const Interval interv(a_solver, a_solver);
  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];

  const phase::which_phase phase = m_cdr->get_phase();
  EBLevelRedist& level_redist = *(m_amr->get_level_redist(phase)[a_lvl]);

  // Again, this is a bit stupid but to get the correct data from the correct interval, we have to do this
  EBAMRIVData diff;
  m_amr->allocate(diff, m_cdr->get_phase(), m_plaskin->get_num_species());
  a_mass_diff.localCopyTo(Interval(0,0), *diff[a_lvl], interv);
  diff[a_lvl]->exchange(interv);

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    //    level_redist.increment(a_mass_diff[dit()], dit(), interv);
    level_redist.increment((*diff[a_lvl])[dit()], dit(), interv);
  }
}

void sisdc::update_coarse_fine_register(const LevelData<BaseIVFAB<Real> >& a_mass_diff,
					const int a_solver,
					const int a_lvl,
					const int a_coarsest_level,
					const int a_finest_level){
  CH_TIME("cdr_solver::update_redist_register");
  if(m_verbosity > 5){
    pout() << "sisdc::::update_redist_register" << endl;
  }

  const Interval interv(a_solver, a_solver);
  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];
  const bool has_coar = a_lvl > a_coarsest_level;
  const bool has_fine = a_lvl < a_finest_level;

  if(has_coar || has_fine){

    const Real dx = m_amr->get_dx()[a_lvl];
    
    // Again, this is a bit stupid but to get the correct data from the correct interval, we have to do this
    EBAMRIVData diff;
    m_amr->allocate(diff, m_cdr->get_phase(), m_plaskin->get_num_species());
    data_ops::set_value(*diff[a_lvl], 0.0);
    a_mass_diff.localCopyTo(Interval(0,0), *diff[a_lvl], interv);
    diff[a_lvl]->exchange();

    RefCountedPtr<EBFineToCoarRedist>& fine2coar_redist = m_amr->get_fine_to_coar_redist(m_cdr->get_phase())[a_lvl];
    RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->get_coar_to_fine_redist(m_cdr->get_phase())[a_lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->get_coar_to_coar_redist(m_cdr->get_phase())[a_lvl];


    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      if(has_coar){
	fine2coar_redist->increment((*diff[a_lvl])[dit()], dit(), interv);
      }

      if(has_fine){
	coar2fine_redist->increment((*diff[a_lvl])[dit()], dit(), interv); 
	coar2coar_redist->increment((*diff[a_lvl])[dit()], dit(), interv); 
      }
    }

    // Tell the flux register about what is going on with EBCF. 
    if(has_fine){
      RefCountedPtr<EBFluxRegister>& fluxreg = m_amr->get_flux_reg(m_cdr->get_phase())[a_lvl];
      fluxreg->incrementRedistRegister(*coar2fine_redist, interv, -dx);
      fluxreg->incrementRedistRegister(*coar2coar_redist, interv, -dx);
    }
  }
  
}

void sisdc::reflux_level(EBAMRCellData& a_state,
			 const int      a_solver,
			 const int      a_lvl,
			 const int      a_coarsest_level,
			 const int      a_finest_level,
			 const Real     a_scale){
  CH_TIME("sisdc::reflux_level");
  if(m_verbosity > 5){
    pout() << "sisdc::::reflux_level" << endl;
  }
  
  const phase::which_phase phase = m_cdr->get_phase();
  const Interval soln_interv(0, 0);
  const Interval flux_interv(a_solver, a_solver);
  const Real dx = m_amr->get_dx()[a_lvl];
  RefCountedPtr<EBFluxRegister >& fluxreg = m_amr->get_flux_reg(phase)[a_lvl];
  fluxreg->reflux(*a_state[a_lvl], soln_interv, flux_interv, 1./dx);
}

void sisdc::redist_level(LevelData<EBCellFAB>&       a_state,
			 const int                   a_solver,   
			 const LevelData<EBCellFAB>& a_weights,
			 const int                   a_lvl){
  CH_TIME("sisdc::redist_level");
  if(m_verbosity > 5){
    pout() << "sisdc::redist_level" << endl;
  }
  const phase::which_phase phase = m_cdr->get_phase();
  const Interval solver_interv(0, 0);
  const Interval redist_interv(a_solver, a_solver);
  EBLevelRedist& level_redist = *(m_amr->get_level_redist(phase)[a_lvl]);
  if(m_cdr->get_mass_redist()){
    level_redist.resetWeights(a_weights, 0);
  }
  level_redist.redistribute(a_state, redist_interv, solver_interv);
}

void sisdc::integrate_advection_subcycle(const Real a_dt, const int a_m, const bool a_corrector){
  CH_TIME("sisdc::integrate_advection_subcycle");
  if(m_verbosity > 5){
    pout() << "sisdc::integrate_advection_subcycle" << endl;
  }

  // Required storages for advection outside of cdr_solver. Don't need all of these for every species, fortunately. 
  const int comp       = 0;
  const int ncomp      = 1;
  const int redist_rad = m_amr->get_redist_rad();

  const phase::which_phase phase = m_cdr->get_phase();

  // These are all temporaries that are required for evaluating the advective derivative. 
  EBAMRFluxData flux;
  EBAMRFluxData face_state;
  EBAMRCellData divF_c;
  EBAMRCellData weights;
  EBAMRIVData   divF_nc;
  EBAMRIVData   mass_diff;
    
  m_amr->allocate(face_state, phase, ncomp);
  m_amr->allocate(flux,       phase, ncomp);
  m_amr->allocate(divF_nc,    phase, ncomp);
  m_amr->allocate(mass_diff,  phase, ncomp);
  m_amr->allocate(divF_c,     phase, ncomp);
  m_amr->allocate(weights,    phase, ncomp, 2*redist_rad);

  // Compute advection velocities. These don't change. Then copy phi_m to storage for phi_(m+1) and start solving from there. 
  sisdc::subcycle_compute_advection_velocities();
  sisdc::subcycle_copy_states(a_m);

  const int coar_lvl = 0;
  const int fine_lvl = m_amr->get_finest_level();

  Vector<Real> tnew(1 + fine_lvl, m_tm[a_m]);
  Vector<Real> told(1 + fine_lvl, m_tm[a_m]);

  // Advance phi[a_m] to phi[a_m+1] using subcycling. First, compute a dt that is below the CFL limit
  int tref = 1;
  for (int lvl = coar_lvl; lvl < fine_lvl; lvl++){
    tref = tref*m_amr->get_ref_rat()[lvl];
  }
  const Real coar_dt_cfl = m_dt_cfl*tref;

  const int nsteps = ceil(m_dtm[a_m]/(coar_dt_cfl*m_cycleCFL));
  const Real dt    = m_dtm[a_m]/nsteps;
  for (int i = 0; i < nsteps; i++){
    sisdc::subcycle_advect_amr(flux, face_state, divF_c, weights, divF_nc, mass_diff, tnew, told,
			       a_m, coar_lvl, coar_lvl, fine_lvl, dt);
  }
}

void sisdc::subcycle_compute_advection_velocities(){
  CH_TIME("sisdc::subcycle_compute_advection_velocities");
  if(m_verbosity > 5){
    pout() << "sisdc::subcycle_compute_advection_velocities" << endl;
  }

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    cdr_gdnv* gdnv = dynamic_cast<cdr_gdnv*>(&(*solver));
    if(gdnv == NULL){
      MayDay::Abort("sisdc::subcycle_compute_advection_velocities - Only cdr_gdnv can subcycle these days...");
    }
    else{
      if(solver->is_mobile()){
	gdnv->average_velo_to_faces();
      }
    }
  }
}

void sisdc::subcycle_copy_states(const int a_m){
  CH_TIME("sisdc::subcycle_copy_states");
  if(m_verbosity > 5){
    pout() << "sisdc::subcycle_copy_states" << endl;
  }

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = sisdc::get_cdr_storage(solver_it);
    EBAMRCellData& phi_m1      = storage->get_phi()[a_m+1];
    const EBAMRCellData& phi_m = storage->get_phi()[a_m];
    
    data_ops::copy(phi_m1, phi_m);
  }

  EBAMRIVData& sigma_m1      = m_sigma_scratch->get_sigma()[a_m+1];
  const EBAMRIVData& sigma_m = m_sigma_scratch->get_sigma()[a_m];
  data_ops::copy(sigma_m1, sigma_m);
}
    
void sisdc::subcycle_advect_amr(EBAMRFluxData& a_flux,
				EBAMRFluxData& a_face_states,
				EBAMRCellData& a_divF_c,
				EBAMRCellData& a_weights,
				EBAMRIVData&   a_divF_nc,
				EBAMRIVData&   a_mass_diff,
				Vector<Real>&  a_tnew,
				Vector<Real>&  a_told,
				const int      a_m,
				const int      a_lvl,
				const int      a_coarsest_level,
				const int      a_finest_level,
				const Real     a_dt){
  CH_TIME("sisdc::subcycle_advect_amr");
  if(m_verbosity > 5){
    pout() << "sisdc::subcycle_advect_amr" << endl;
  }

  Real coar_time_old = 0.0;
  Real coar_time_new = 0.0;

  const bool has_fine = a_lvl < a_finest_level;
  const bool has_coar = a_lvl > a_coarsest_level;

  if(has_coar){
    coar_time_old = a_told[a_lvl-1];
    coar_time_new = a_tnew[a_lvl-1];
  }

  // Prepare level solve
  sisdc::subcycle_copy_current_to_old_states(a_m, a_lvl);
  sisdc::reset_finer_flux_registers_level(a_lvl, a_coarsest_level, a_finest_level);
  sisdc::reset_redist_registers_level(a_lvl, a_coarsest_level, a_finest_level);

  // Level solve. Update boundary conditions and source terms on this level
  sisdc::subcycle_update_transport_bc(a_m, a_lvl, a_tnew[a_lvl]);
  if(m_cycle_sources){
    sisdc::subcycle_update_sources(a_m, a_lvl, a_tnew[a_lvl]);
  }
  sisdc::subcycle_integrate_level(*a_flux[a_lvl],
				  *a_face_states[a_lvl],
				  *a_divF_c[a_lvl],
				  *a_weights[a_lvl],
				  *a_mass_diff[a_lvl],
				  *a_divF_nc[a_lvl],
				  a_m,
				  a_lvl,
				  a_coarsest_level,
				  a_finest_level,
				  coar_time_old,
				  coar_time_new,
				  a_tnew[a_lvl],
				  a_dt);

  // We have advance this level. Updates new times
  a_told[a_lvl] = a_tnew[a_lvl];
  a_tnew[a_lvl] = a_tnew[a_lvl] + a_dt;

  // If there is a finer level, advance it nref times so that we synchronize.
  if(has_fine){

    int nref;
    Real dt_ref;
    if(!m_optimal_subcycling){ // Standard, use refinement ratio
      nref   = m_amr->get_ref_rat()[a_lvl];
      dt_ref = a_dt/nref;
    }
    else{ // Now advance a_lvl+1
      int tref = 1;
      for (int lvl = a_coarsest_level; lvl < a_finest_level; lvl++){
      	tref *= m_amr->get_ref_rat()[lvl];
      }
      
      Real dt_cfl_level = m_dt_cfl*tref;
      for (int lvl = a_coarsest_level; lvl <= a_lvl; lvl++){
      	dt_cfl_level /= m_amr->get_ref_rat()[lvl];
      }
    
      // Try to take as few steps as possible. We will ensure that local cfl = [minCFL, maxCFL]
      nref   = 1;
      dt_ref = a_dt;
      while(dt_ref > m_cycleCFL*dt_cfl_level){
	nref++;
	dt_ref = a_dt/nref;
      }

#if 0 // Debug
    if(procID() == 0 && a_lvl == 0) std::cout << std::endl;
    if(procID() == 0) std::cout << "lvl = " << a_lvl
				<< " dt_cfl = " << m_dt_cfl
				<< " dt = " << a_dt
				<< " .....On next level I'm taking nref = " << nref << " steps"
				<< " with time = " << dt_ref 
				<< " for total time = " << dt_ref*nref
				<< " with level cfl = " << dt_ref/dt_cfl_level
				<< endl;
#endif
    }

    for (int i = 0; i < nref; i++){
      sisdc::subcycle_advect_amr(a_flux, a_face_states, a_divF_c, a_weights, a_divF_nc, a_mass_diff,
				 a_tnew, a_told, a_m, a_lvl+1, a_coarsest_level, a_finest_level, dt_ref);
    }


    // Finer level has caught up. Sync levels
    sisdc::subcycle_sync_levels(a_m, a_lvl, a_coarsest_level, a_finest_level);
  }
}

void sisdc::subcycle_update_transport_bc(const int a_m, const int a_lvl, const Real a_time){
  CH_TIME("sisdc::subcycle_update_transport_bc");
  if(m_verbosity > 5){
    pout() << "sisdc::subcycle_update_transport_bc" << endl;
  }

  Vector<LevelData<EBCellFAB>* >        cell_states;      
  Vector<LevelData<EBCellFAB>* >        cell_gradients;   
  Vector<LevelData<EBCellFAB>* >        cell_velocities;
  Vector<LevelData<BaseIVFAB<Real> >* > solver_eb_fluxes;
  Vector<LevelData<BaseIVFAB<Real> >* > eb_states;
  Vector<LevelData<BaseIVFAB<Real> >* > eb_velocities;
  Vector<LevelData<BaseIVFAB<Real> >* > eb_fluxes;
  Vector<LevelData<BaseIVFAB<Real> >* > eb_gradients;
  Vector<LevelData<BaseIVFAB<Real> >* > rte_fluxes;
  Vector<LevelData<DomainFluxIFFAB>* >  domain_states;
  Vector<LevelData<DomainFluxIFFAB>* >  domain_gradients;
  
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = sisdc::get_cdr_storage(solver_it);

    cell_states.push_back(storage->get_phi()[a_m+1][a_lvl]); // This is what we are updating
    cell_gradients.push_back(storage->get_gradient()[a_lvl]);
    cell_velocities.push_back(solver->get_velo_cell()[a_lvl]);
    
    eb_states.push_back(storage->get_eb_state()[a_lvl]);
    eb_velocities.push_back(storage->get_eb_velo()[a_lvl]);
    eb_gradients.push_back(storage->get_eb_grad()[a_lvl]);
    eb_fluxes.push_back(storage->get_eb_flux()[a_lvl]);
    
    domain_states.push_back(storage->get_domain_state()[a_lvl]);
    domain_gradients.push_back(storage->get_domain_grad()[a_lvl]);

    solver_eb_fluxes.push_back(solver->get_ebflux()[a_lvl]);
  }

  // Compute all the things that are necessary for BCs
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    RefCountedPtr<cdr_storage>& storage = sisdc::get_cdr_storage(solver_it);

    // Storage we can use for computations
    LevelData<EBCellFAB>&        scratchD    = *storage->get_scratchD()[a_lvl];
    LevelData<BaseIVFAB<Real> >& scratchIV_1 = *storage->get_eb_scratch1()[a_lvl]; // Scalar
    LevelData<BaseIVFAB<Real> >& scratchIV_D = *storage->get_eb_scratchD()[a_lvl]; // SpaceDim comps

    // 1. Compute cell-centered gradients
    cell_states[idx]->exchange();
    m_amr->compute_gradient(*cell_gradients[idx], *cell_states[idx], a_lvl);
    
    // 2. Extrapolate cell-centered gradient to the EB
    time_stepper::extrapolate_to_eb(scratchIV_D, m_cdr->get_phase(), *cell_gradients[idx], a_lvl);

    // 3. Dot EB-centered gradient with normal vector
    time_stepper::project_flux(*eb_gradients[idx], scratchIV_D, a_lvl);
    
    // 4. Extrapolate the cell-centered states to the EB
    time_stepper::extrapolate_to_eb(*eb_states[idx], m_cdr->get_phase(), *cell_states[idx], a_lvl);

    // 5. Extrapolate cell-centered velocities to the EB
    time_stepper::extrapolate_to_eb(scratchIV_D, m_cdr->get_phase(), *cell_velocities[idx], a_lvl);

    // 6. Project normal velocity
    time_stepper::project_flux(*eb_velocities[idx], scratchIV_D, a_lvl);

    // 7. Compute the extrapolated flux at the boundary
    cell_velocities[idx]->localCopyTo(scratchD);
    data_ops::multiply_scalar(scratchD, *cell_states[idx]);
    time_stepper::extrapolate_to_eb(scratchIV_D, m_cdr->get_phase(), scratchD, a_lvl);
    time_stepper::project_flux(*eb_fluxes[idx], scratchIV_D, a_lvl);
  }

  // Radiative transfer fluxes at boundary
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);

    EBAMRIVData& flux_eb = storage->get_eb_flux();
    solver->compute_boundary_flux(flux_eb, solver->get_state());
    rte_fluxes.push_back(flux_eb[a_lvl]);
  }

  // Electric field at boundary
  const EBAMRIVData& E = m_poisson_scratch->get_E_eb();
  
  // Update the stinking EB fluxes
  time_stepper::compute_cdr_fluxes(solver_eb_fluxes,
				   eb_fluxes,
				   eb_states,
				   eb_velocities,
				   eb_gradients,
				   rte_fluxes,
				   *E[a_lvl],
				   a_time,
				   a_lvl);

#if 0 
  MayDay::Warning("sisdc::update_transport_bc - we have not yet done the domain boundary conditiosn!!!");
#endif
}

void sisdc::subcycle_update_sources(const int a_m, const int a_lvl, const Real a_time){
  CH_TIME("sisdc::subcycle_update_sources");
  if(m_verbosity > 5){
    pout() << "sisdc::subcycle_update_sources" << endl;
  }

#if 0 //
  // Actually, we DO need to interpolate the coarse grid data to a_time in order to get update ghost cells.
  MayDay::Warning("sisdc::subcycle_update_sources - Coarse grid data should be interpolated to current time");
#endif

  // This should have been called AFTER update_transport_bc, which computed the cell centered gradients. So
  // at least we don't have to do that again...
  const int num_species        = m_plaskin->get_num_species();
  const int num_photons        = m_plaskin->get_num_photons();
  
  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_cdr->get_phase())[a_lvl];
  const RealVect origin        = m_physdom->get_prob_lo();
  const Real dx                = m_amr->get_dx()[a_lvl];

  // Stencils for extrapolating things to cell centroids
  const irreg_amr_stencil<centroid_interp>& interp_stencils = m_amr->get_centroid_interp_stencils(m_cdr->get_phase());

  // We must have the gradient of E. This block of code does that. 
  EBAMRCellData grad_E, E_norm;
  m_amr->allocate(grad_E, m_cdr->get_phase(), SpaceDim);  // Allocate storage for grad(|E|)
  m_amr->allocate(E_norm, m_cdr->get_phase(), 1);         // Allocate storage for |E|
  const EBAMRCellData& E = m_poisson_scratch->get_E_cell();
  data_ops::vector_length(*E_norm[a_lvl], *E[a_lvl]);            // Compute |E| on this level
  m_amr->compute_gradient(*grad_E[a_lvl], *E_norm[a_lvl], a_lvl);// Compute grad(|E|) on this level

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      Vector<EBCellFAB*> sources(num_species);
      Vector<EBCellFAB*> cdr_densities(num_species);
      Vector<EBCellFAB*> cdr_gradients(num_species);
      Vector<EBCellFAB*> rte_densities(num_photons);
      for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	RefCountedPtr<cdr_solver>& solver   = solver_it();
	RefCountedPtr<cdr_storage>& storage = sisdc::get_cdr_storage(solver_it);
	const int idx = solver_it.get_solver();
	EBAMRCellData& src   = solver->get_source();
	EBAMRCellData& phim  = storage->get_phi()[a_m+1];
	EBAMRCellData& gradm = storage->get_gradient();
	
	sources[idx]       = &((*src[a_lvl])[dit()]);
	cdr_densities[idx] = &((*phim[a_lvl])[dit()]);
	cdr_gradients[idx] = &((*gradm[a_lvl])[dit()]);
      }
      for (rte_iterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      	const int idx = solver_it.get_solver();
	RefCountedPtr<rte_solver>& solver = solver_it();
	EBAMRCellData& state = solver->get_state();
      	rte_densities[idx] = &((*state[a_lvl])[dit()]);
      }
      const EBCellFAB& e  = (*E[a_lvl])[dit()];
      const EBCellFAB& gE = (*grad_E[a_lvl])[dit()];

      // This does all cells
      time_stepper::compute_cdr_sources_reg(sources,
      					    cdr_densities,
      					    cdr_gradients,
      					    rte_densities,
      					    e,
      					    gE,
      					    dbl.get(dit()),
      					    a_time,
      					    dx);

      // Have to redo irregular cells
      time_stepper::compute_cdr_sources_irreg(sources,
					      cdr_densities,
					      cdr_gradients,
					      rte_densities,
					      e,
					      gE,
					      interp_stencils[a_lvl][dit()],
					      dbl.get(dit()),
					      a_time,
					      dx);
  }
}

void sisdc::subcycle_sync_levels(const int a_m, const int a_lvl, const int a_coarsest_level, const int a_finest_level){
  CH_TIME("sisdc::subcycle_sync_levels");
  if(m_verbosity > 5){
    pout() << "sisdc::subcycle_sync_levels" << endl;
  }

  // Finer level has reached this level. Average down solution on this level and reflux mass.
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int solver_idx                = solver_it.get_solver();
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = sisdc::get_cdr_storage(solver_it);

    EBAMRCellData& state = storage->get_phi()[a_m+1]; // This is the one that we update

    // Since the redistribution registers don't allow components, do some transient storage stuff
    
    
    // Reflux state
    if(solver->is_mobile()){
      m_amr->average_down(state, m_cdr->get_phase(), a_lvl);
      sisdc::reflux_level(state, solver_idx, a_lvl, a_coarsest_level, a_finest_level, 1.0);
      // EBCF related code. 
      if(m_amr->get_ebcf()){
	const Interval inter0(0,0);
	const Interval interv(solver_idx, solver_idx);
	
	const bool has_fine = a_lvl < a_finest_level;
	const bool has_coar = a_lvl > a_coarsest_level;

	// Bah, extra storage becase redistribution registers don't let me use different for mass diffs and target FAB
	EBAMRCellData dummy;
	m_amr->allocate(dummy, m_cdr->get_phase(), m_plaskin->get_num_species());
	data_ops::set_value(dummy, 0.0);

	RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->get_coar_to_fine_redist(m_cdr->get_phase())[a_lvl];
	RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->get_coar_to_coar_redist(m_cdr->get_phase())[a_lvl];
	RefCountedPtr<EBFineToCoarRedist>& fine2coar_redist = m_amr->get_fine_to_coar_redist(m_cdr->get_phase())[a_lvl];


	LevelData<EBCellFAB>* dummy_lev  = dummy[a_lvl];
	LevelData<EBCellFAB>* dummy_coar = has_coar ? &(*dummy[a_lvl-1]) : NULL;
	LevelData<EBCellFAB>* dummy_fine = has_fine ? &(*dummy[a_lvl+1]) : NULL;

	// Level copies
	if(has_coar || has_fine){
	  state[a_lvl]->localCopyTo(inter0, *dummy_lev, interv);

	  if(has_coar) state[a_lvl-1]->localCopyTo(inter0, *dummy_coar, interv);
	  if(has_fine) state[a_lvl+1]->localCopyTo(inter0, *dummy_fine, interv);
	
	  if(has_coar){
	    state[a_lvl-1]->localCopyTo(Interval(0,0), *dummy_coar, interv);
	    fine2coar_redist->redistribute(*dummy_coar, interv);
	  }

	  if(has_fine){
	    coar2fine_redist->redistribute(*dummy_fine, interv);
	    coar2coar_redist->redistribute(*dummy_lev,  interv);
	  }

	  // Copy back
	  dummy_lev->localCopyTo(interv, *state[a_lvl], inter0);
	  if(has_coar) dummy_coar->localCopyTo(interv, *state[a_lvl-1], inter0);
	  if(has_fine) dummy_fine->localCopyTo(interv, *state[a_lvl+1], inter0);
	}
      }
    }
    m_amr->average_down(state, m_cdr->get_phase(), a_lvl);
    state[a_lvl]->exchange();
  }
}

void sisdc::subcycle_copy_current_to_old_states(const int a_m, const int a_lvl){
  CH_TIME("sisdc::subcycle_copy_current_to_old_states");
  if(m_verbosity > 5){
    pout() << "sisdc::subcycle_copy_current_to_old_states" << endl;
  }

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = sisdc::get_cdr_storage(solver_it);

    EBAMRCellData& old = storage->get_old();
    const EBAMRCellData& current = storage->get_phi()[a_m+1];

    current[a_lvl]->localCopyTo(*old[a_lvl]);
    old[a_lvl]->exchange();
  }
}
  
void sisdc::subcycle_integrate_level(LevelData<EBFluxFAB>&        a_flux,
				     LevelData<EBFluxFAB>&        a_face_states,
				     LevelData<EBCellFAB>&        a_divF_c,
				     LevelData<EBCellFAB>&        a_weights,
				     LevelData<BaseIVFAB<Real> >& a_mass_diff,
				     LevelData<BaseIVFAB<Real> >& a_divF_nc,
				     const int                    a_m,
				     const int                    a_lvl,
				     const int                    a_coarsest_level,
				     const int                    a_finest_level,
				     const Real                   a_coar_time_old,
				     const Real                   a_coar_time_new,
				     const Real                   a_time,
				     const Real                   a_dt){
  CH_TIME("sisdc::subcycle_integrate_level");
  if(m_verbosity > 5){
    pout() << "sisdc::subcycle_integrate_level" << endl;
  }

#if 0 //
  if(procID() == 0) std::cout << "Integrating level = " << a_lvl << std::endl;
  if(procID() == 0) std::cout << "time debug: " << a_coar_time_old
			      << "\t" << a_time << "\t" << a_coar_time_new << std::endl;
#endif

  const bool ebcf     = m_amr->get_ebcf();
  const bool has_coar = a_lvl > a_coarsest_level;
  const bool has_fine = a_lvl < a_finest_level;
  
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int solver_idx = solver_it.get_solver();
    RefCountedPtr<cdr_storage>& storage = sisdc::get_cdr_storage(solver_it);
    RefCountedPtr<cdr_solver>& solver   = solver_it();

    LevelData<EBCellFAB>& state_m1 = (*storage->get_phi()[a_m+1][a_lvl]); // We will update this one.
    LevelData<EBCellFAB>& source   = *solver->get_source()[a_lvl];

    if(solver->is_mobile()){
      
      cdr_gdnv* gdnv = (cdr_gdnv*) (&(*solver));
      LevelData<EBCellFAB>* coar_old = NULL;
      LevelData<EBCellFAB>* coar_new = NULL;

      if(has_coar){
	coar_old = storage->get_old()[a_lvl-1];
	coar_new = storage->get_phi()[a_m+1][a_lvl-1];
      }
      
      // Advect to faces and compute fluxes on face centers, and compute the conservative divergence on regular cells
      const Real extr_dt = m_extrap_advect ? 2.0*m_extrap_dt*a_dt : 0.0;
      gdnv->advect_to_faces(a_face_states, state_m1, coar_old, coar_new, a_time, a_coar_time_old,a_coar_time_new, a_lvl, extr_dt);
      gdnv->new_compute_flux(a_flux, a_face_states, a_lvl);
      gdnv->consdiv_regular(a_divF_c, a_flux, a_lvl);

      // Set up flux interpolant and compute conservative flux on irregular cells
      LevelData<BaseIFFAB<Real> > flux_interp[SpaceDim];
      const LevelData<BaseIVFAB<Real> >& ebflux = *solver->get_ebflux()[a_lvl];
      gdnv->setup_flux_interpolant(flux_interp, a_flux, a_lvl); // Compute interpolant
      gdnv->interpolate_flux_to_centroids(flux_interp, a_lvl);  // Interpolant now holds face centroid-centered fluxes
      gdnv->compute_divF_irreg(a_divF_c, flux_interp, ebflux, a_lvl); // Update the conservative divergence in the cut cells

      // Recess: So far the conservative divergence is scaled by 1/dx but not yet divided by the volume fraction. This
      // means that the actual advance without the hybrid stuff would be
      //
      // new_state -= dt*a_divF_c/kappa
      //
      // The stuff below was originally written for d(phi)/dt = Div(F) rather than d(phi)/dt = -Div(F) so that's why
      // there's a (-a_dt) in all the stuff below. This design choice was made because I am, in fact, an ass.

      // Compute the nonconservative and hybrid divergences (hybrid put on storage for divF_c, which is lost)
      gdnv->nonconservative_divergence(a_divF_nc, a_face_states, a_lvl);
      gdnv->hybrid_divergence(a_divF_c, a_mass_diff, a_divF_nc, a_lvl); // Puts hybrid in a_divF_c. mass_diff as usual without dt,
      data_ops::scale(a_mass_diff, -a_dt);                              // Sign convention

      // Update flux and redistribution registers
      gdnv->new_compute_flux(a_flux, a_face_states, a_lvl);
      sisdc::update_flux_registers(a_flux, solver_idx, a_lvl, a_coarsest_level, a_finest_level, a_dt);
      sisdc::update_redist_register(a_mass_diff, solver_idx, a_lvl);

      // If we have EBCF, we must update those as well
      if(ebcf){
	sisdc::update_coarse_fine_register(a_mass_diff, solver_idx, a_lvl, a_coarsest_level, a_finest_level);
      }
      
      if(m_cdr->get_mass_redist()){
	data_ops::set_value(a_weights, 0.0);
	data_ops::incr(a_weights, state_m1, 1.0);
      }

      // Euler advance with redistribution
      data_ops::incr(state_m1, a_divF_c, -a_dt);
      sisdc::redist_level(state_m1, solver_idx, a_weights, a_lvl);
    }

    // Add in the source term
    if(m_cycle_sources){
      data_ops::incr(state_m1, source, a_dt);
    }
    state_m1.exchange();
  }
}

void sisdc::store_solvers(){
  CH_TIME("sisdc::store_solvers");
  if(m_verbosity > 5){
    pout() << "sisdc::store_solvers" << endl;
  }

  // SISDC does not manipulate cdr and sigma solvers until the end of the time step. Only need to do
  // Poisson and RTE here. 

  // Poisson
  MFAMRCellData& previous    = m_poisson_scratch->get_previous();
  const MFAMRCellData& state = m_poisson->get_state();
  data_ops::copy(previous, state);

  // RTE
  for (rte_iterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_storage>& storage     = sisdc::get_rte_storage(solver_it);
    const RefCountedPtr<rte_solver>& solver = solver_it();

    EBAMRCellData& previous = storage->get_previous();
    const EBAMRCellData& state = solver->get_state();

    data_ops::copy(previous, state);
  }
}

void sisdc::restore_solvers(){
  CH_TIME("sisdc::restore_solvers");
  if(m_verbosity > 5){
    pout() << "sisdc::restore_solvers" << endl;
  }

  // SISDC does not manipulate cdr and sigma solvers until the end of the time step. Only need to do
  // Poisson and RTE here. 

  // Poisson
  MFAMRCellData& state = m_poisson->get_state();
  const MFAMRCellData& previous    = m_poisson_scratch->get_previous();

  data_ops::copy(state, previous);

  // RTE
  for (rte_iterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_storage>& storage     = sisdc::get_rte_storage(solver_it);
    RefCountedPtr<rte_solver>& solver = solver_it();

    EBAMRCellData& previous = storage->get_previous();
    EBAMRCellData& state = solver->get_state();

    data_ops::copy(state, previous);
  }
}
