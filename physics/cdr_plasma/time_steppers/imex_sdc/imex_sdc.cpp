/*!
  @file   imex_sdc.cpp
  @brief  Implementation of imex_sdc.H
  @author Robert Marskar
  @date   Feb. 2019
*/

#include "imex_sdc.H"
#include "imex_sdc_storage.H"
#include "data_ops.H"
#include "units.H"
#include "cdr_gdnv.H"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <ParmParse.H>

using namespace physics::cdr_plasma;

typedef imex_sdc::cdr_storage     cdr_storage;
typedef imex_sdc::poisson_storage poisson_storage;
typedef imex_sdc::rte_storage     rte_storage;
typedef imex_sdc::sigma_storage   sigma_storage;

imex_sdc::imex_sdc(){
  m_class_name = "imex_sdc";
  m_subcycle   = false;
}

imex_sdc::imex_sdc(RefCountedPtr<cdr_plasma_physics>& a_physics) : imex_sdc() {
  m_physics = a_physics;
}

imex_sdc::~imex_sdc(){

}

void imex_sdc::parse_options(){
  CH_TIME("imex_sdc::parse_options");
  if(m_verbosity > 5){
    pout() << "imex_sdc::parse_options" << endl;
  }

  // Regular stuff from cdr_plasma_stepper that we almost always need
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
  parse_debug_options();
  parse_advection_options();

  // Setup nodes
  imex_sdc::setup_quadrature_nodes(m_p);
  imex_sdc::setup_qmj(m_p);
}

void imex_sdc::parse_runtime_options(){
  CH_TIME("imex_sdc::parse_runtime_options");
  if(m_verbosity > 5){
    pout() << "imex_sdc::parse_runtime_options" << endl;
  }

  // Regular stuff from cdr_plasma_stepper that we almost always need
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
  parse_debug_options();
  parse_advection_options();

  // Setup nodes
  imex_sdc::setup_quadrature_nodes(m_p);
  imex_sdc::setup_qmj(m_p);

  m_cdr->parse_runtime_options();
  m_rte->parse_runtime_options();
  m_poisson->parse_runtime_options();
}

void imex_sdc::parse_nodes(){
  CH_TIME("imex_sdc::parse_nodes");
  if(m_verbosity > 5){
    pout() << "imex_sdc::parse_nodes" << endl;
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
    MayDay::Abort("imex_sdc::parse_nodes - unknown node type requested");
  }

  pp.get("subintervals",     m_p);
  pp.get("corr_iter",        m_k);

  if(m_p < 1){
    MayDay::Abort("imex_sdc::parse_nodes - imex_sdc.subintervals cannot be < 1");
  }
  if(m_k < 0){
    MayDay::Abort("imex_sdc::parse_nodes - imex_sdc.corr_iter cannot be < 0");
  }
}

void imex_sdc::parse_diffusion_coupling(){
  CH_TIME("imex_sdc::parse_diffusion_coupling");
  if(m_verbosity > 5){
    pout() << "imex_sdc::parse_diffusion_coupling" << endl;
  }

  ParmParse pp(m_class_name.c_str());
  
  std::string str;

  pp.get("use_tga", str);
  m_use_tga = (str == "true") ? true : false;
}

void imex_sdc::parse_adaptive_options(){
  CH_TIME("imex_sdc::parse_adaptive_options");
  if(m_verbosity > 5){
    pout() << "imex_sdc::parse_adaptive_options" << endl;
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

void imex_sdc::parse_debug_options(){
  CH_TIME("imex_sdc::parse_debug_options");
  if(m_verbosity > 5){
    pout() << "imex_sdc::parse_debug_options" << endl;
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

void imex_sdc::parse_advection_options(){
  CH_TIME("imex_sdc::parse_advection_options");
  if(m_verbosity > 5){
    pout() << "imex_sdc::parse_advection_options" << endl;
  }

  ParmParse pp(m_class_name.c_str());
  std::string str;

  m_extrap_dt = 0.5;  // Relic of an ancient past. I don't see any reason why extrapolating to anything but the half interval
                      // would make sense. 

  pp.get("extrap_advect", str); m_extrap_advect = (str == "true") ? true : false;
}

RefCountedPtr<cdr_storage>& imex_sdc::get_cdr_storage(const cdr_iterator<cdr_solver>& a_solverit){
  return m_cdr_scratch[a_solverit.index()];
}

RefCountedPtr<rte_storage>& imex_sdc::get_rte_storage(const rte_iterator<rte_solver>& a_solverit){
  return m_rte_scratch[a_solverit.index()];
}

bool imex_sdc::need_to_regrid(){
  CH_TIME("imex_sdc::need_to_regrid");
  if(m_verbosity > 5){
    pout() << "imex_sdc::need_to_regrid" << endl;
  }

  return false;
}

Real imex_sdc::restrict_dt(){
  return 1.E99;
}

Real imex_sdc::get_max_node_distance(){
  CH_TIME("imex_sdc::get_max_node_distance");
  if(m_verbosity > 5){
    pout() << "imex_sdc::get_max_node_distance" << endl;
  }

  Real max_dist = 0.0;
  for (int m = 0; m < m_p; m++){
    max_dist = Max(max_dist, m_nodes[m+1] - m_nodes[m]);
  }

  return max_dist;
}

void imex_sdc::init(){
  CH_TIME("imex_sdc::init");
  if(m_verbosity > 5){
    pout() << "imex_sdc::init" << endl;
  }

  advance_reaction_network(m_time, m_dt);
}

void imex_sdc::setup_quadrature_nodes(const int a_p){
  CH_TIME("imex_sdc::setup_quadrature_nodes");
  if(m_verbosity > 5){
    pout() << "imex_sdc::setup_quadrature_nodes" << endl;
  }

  if(m_which_nodes == "uniform"){
    imex_sdc::setup_uniform_nodes(a_p);
  }
  else if(m_which_nodes == "lobatto"){
    imex_sdc::setup_lobatto_nodes(a_p);
  }
  else if(m_which_nodes == "chebyshev"){
    imex_sdc::setup_chebyshev_nodes(a_p);
  }
  else {
    MayDay::Abort("imex_sdc::setup_quadrature_nodes - unknown nodes requested");
  }
}

void imex_sdc::setup_uniform_nodes(const int a_p){
  CH_TIME("imex_sdc::setup_uniform_nodes");
  if(m_verbosity > 5){
    pout() << "imex_sdc::setup_uniform_nodes" << endl;
  }

  // TLDR: The nodes and weights are hardcoded. A better programmer would compute these
  //       recursively with Legendre polynomials. 
  m_nodes.resize(1+a_p);

  const Real delta = 2./a_p;
  for (int m = 0; m <= a_p; m++){
    m_nodes[m] = m*delta;
  }
}

void imex_sdc::setup_lobatto_nodes(const int a_p){
  CH_TIME("imex_sdc::setup_lobatto_nodes");
  if(m_verbosity > 5){
    pout() << "imex_sdc::setup_lobatto_nodes" << endl;
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
    MayDay::Abort("imex_sdc::setup_lobatto_nodes - requested order exceeds 7. Compute your own damn nodes!");
  }
}

void imex_sdc::setup_chebyshev_nodes(const int a_p){
  CH_TIME("imex_sdc::setup_chebyshev_nodes");
  if(m_verbosity > 5){
    pout() << "imex_sdc::setup_chebyshev_nodes" << endl;
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

void imex_sdc::setup_qmj(const int a_p){
  CH_TIME("imex_sdc::setup_qmj");
  if(m_verbosity > 5){
    pout() << "imex_sdc::setup_qmj" << endl;
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
    if(INFO != 0) MayDay::Abort("imex_sdc::setup_qmj - could not compute weights");
    
    // Now construct qmj
    for (int m = 0; m < a_p; m++){
      m_qmj[m][j] = 0.0;
      for (int k = 0; k < nnodes; k++){
	m_qmj[m][j] += cj[k]*(pow(m_nodes[m+1], k+1) - pow(m_nodes[m], k+1))/(k+1);
      }
    }
  }
}

void imex_sdc::setup_subintervals(const Real a_time, const Real a_dt){
  CH_TIME("imex_sdc::setup_subintervals");
  if(m_verbosity > 5){
    pout() << "imex_sdc::setup_subintervals" << endl;
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

  // dtm = t_{m+1} - t_m. Order 1 is special since we only use the IMEX_SDC predictor from a second order formulation
  m_dtm.resize(m_tm.size() - 1);
  for (int m = 0; m < m_tm.size()-1; m++){
    m_dtm[m] = m_tm[m+1] - m_tm[m];
  }
}

void imex_sdc::quad(EBAMRCellData& a_quad, const Vector<EBAMRCellData>& a_integrand, const int a_m){
  CH_TIME("imex_sdc::quad");
  if(m_verbosity > 5){
    pout() << "imex_sdc::quad" << endl;
  }

  if(a_m < 0)     MayDay::Abort("imex_sdc::quad - bad index a_m < 0");
  if(a_m >= m_p)  MayDay::Abort("imex_sdc::quad - bad index a_m >= m_p");

  data_ops::set_value(a_quad, 0.0);
  for (int j = 0; j <= m_p; j++){
    data_ops::incr(a_quad, a_integrand[j], m_qmj[a_m][j]);
  }
}

void imex_sdc::quad(EBAMRIVData& a_quad, const Vector<EBAMRIVData>& a_integrand, const int a_m){
  CH_TIME("imex_sdc::quad");
  if(m_verbosity > 5){
    pout() << "imex_sdc::quad" << endl;
  }

  if(a_m < 0)     MayDay::Abort("imex_sdc::quad - bad index a_m < 0");
  if(a_m >= m_p)  MayDay::Abort("imex_sdc::quad - bad index a_m >= m_p");

  data_ops::set_value(a_quad, 0.0);
  for (int j = 0; j <= m_p; j++){
    data_ops::incr(a_quad, a_integrand[j], m_qmj[a_m][j]);
  }
}
  
void imex_sdc::copy_phi_p_to_cdr(){
  CH_TIME("imex_sdc::copy_phi_p_to_cdr");
  if(m_verbosity > 5){
    pout() << "imex_sdc::copy_phi_p_to_cdr" << endl;
  }

  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>&  solver  = solver_it();
    RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);

    EBAMRCellData& phi = solver->get_state();
    const EBAMRCellData& phip = storage->get_phi()[m_p];
    data_ops::copy(phi, phip);
  }
}

void imex_sdc::copy_sigma_p_to_sigma(){
  CH_TIME("imex_sdc::copy_sigma_p_to_sigma");
  if(m_verbosity > 5){
    pout() << "imex_sdc::copy_sigma_p_to_sigma" << endl;
  }

  EBAMRIVData& sigma        = m_sigma->get_state();
  const EBAMRIVData& sigmap = m_sigma_scratch->get_sigma()[m_p];
  data_ops::copy(sigma, sigmap);
}

Real imex_sdc::advance(const Real a_dt){
  CH_TIME("imex_sdc::advance");
  if(m_verbosity > 2){
    pout() << "imex_sdc::advance" << endl;
  }

  // ---------------------------------------------------------------------------------------------------
  // TLDR:  When we enter this routine, solvers SHOULD have been filled with valid ready and be ready 
  //        advancement. If you think that this may not be the case, activate the debugging below
  // ---------------------------------------------------------------------------------------------------

  // Initialize integrations. If we do corrections, we need FD(phi_0) since this is implicit. If we do adaptive_dt, we should
  // also a
  imex_sdc::copy_cdr_to_phi_m0();
  imex_sdc::copy_sigma_to_sigma_m0();
  imex_sdc::compute_FD_0();
  imex_sdc::store_solvers();

  // IMEX_SDC advance
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
    imex_sdc::setup_subintervals(m_time, actual_dt);

    // First SDC sweep. No lagged slopes here. 
    imex_sdc::integrate(actual_dt, m_time, false);

    // SDC correction sweeps. Need to take care of lagged terms. 
    for(int icorr = 0; icorr < Max(m_k, m_min_corr); icorr++){
      num_corrections++;

      // Initialize error and reconcile integrands (i.e. make them quadrature-ready)
      imex_sdc::initialize_errors();
      imex_sdc::reconcile_integrands();

      // SDC correction along whole interval
      imex_sdc::integrate(actual_dt, m_time, true);

      // Compute error and check if we need to keep iterating
      imex_sdc::finalize_errors();
      if(m_max_error < m_err_thresh && m_adaptive_dt && icorr >= m_min_corr) break; // No need in going beyond
    }

    // Compute a new time step. If it is smaller than the minimum allowed CFL step, accept the step anyways
    if(m_adaptive_dt){
      imex_sdc::compute_new_dt(accept_step, actual_dt, num_corrections);
      
      if(!accept_step){  // Step rejection, use the new dt for next step.
	actual_dt = m_new_dt;
	num_reject++;

	retry_step  = num_reject <= m_max_retries;
	
	if(retry_step){
	  imex_sdc::restore_solvers();
	  imex_sdc::compute_E_into_scratch();
	  imex_sdc::compute_cdr_gradients();
	  imex_sdc::compute_cdr_velo(m_time);
	  cdr_plasma_stepper::compute_cdr_diffusion(m_poisson_scratch->get_E_cell(), m_poisson_scratch->get_E_eb());
	}
      }
    }
    else{
      m_new_dt = 1.234567E89;
      accept_step = true;
    }
  }

  // Copy results back to solvers
  imex_sdc::copy_phi_p_to_cdr();
  imex_sdc::copy_sigma_p_to_sigma();

  // Always recompute velocities and diffusion coefficients before the next time step. The Poisson and RTE equations
  // have been updated when we come in here. 
  imex_sdc::compute_cdr_velo(m_time + actual_dt);
  cdr_plasma_stepper::compute_cdr_diffusion(m_poisson_scratch->get_E_cell(), m_poisson_scratch->get_E_eb());


  // Profile step
  if(m_print_report)  imex_sdc::adaptive_report(first_dt, actual_dt, m_new_dt, num_corrections, num_reject, m_max_error);
  if(m_profile_steps) imex_sdc::write_step_profile(actual_dt, m_max_error, m_p, num_corrections, num_reject);

  // Store current error. 
  m_have_err  = true;
  m_pre_error = m_max_error;
  
  return actual_dt;
}

void imex_sdc::copy_cdr_to_phi_m0(){
  CH_TIME("imex_sdc::copy_cdr_to_phi_m0");
  if(m_verbosity > 5){
    pout() << "imex_sdc::copy_cdr_to_phi_m0" << endl;
  }

  // CDR solvers
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>&  solver  = solver_it();
    RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);
    
    EBAMRCellData& phi0 = storage->get_phi()[0];
    const EBAMRCellData& phi = solver->get_state();
    data_ops::copy(phi0, phi);
  }
}

void imex_sdc::copy_sigma_to_sigma_m0(){
  CH_TIME("imex_sdc::copy_sigma_sigma_m0");
  if(m_verbosity > 5){
    pout() << "imex_sdc::copy_sigma_to_sigma_m0" << endl;
  }

  // Copy sigma to starting state
  EBAMRIVData& sigma0      = m_sigma_scratch->get_sigma()[0];
  const EBAMRIVData& sigma = m_sigma->get_state();
  data_ops::copy(sigma0, sigma);
}

void imex_sdc::compute_FD_0(){
  CH_TIME("imex_sdc::compute_FD_0");
  if(m_verbosity > 5){
    pout() << "imex_sdc::compute_FD_0" << endl;
  }

  if(m_k > 0){ // We only need this if we're actually doing any corrections....
    for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);
    
      EBAMRCellData& phi_0 = storage->get_phi()[0]; // phi_0
      EBAMRCellData& FD_0  = storage->get_FD()[0];  // FD(phi_0)
    
      if(solver->is_diffusive()){
	solver->compute_divD(FD_0, phi_0);

	// Shouldn't be necesary
	m_amr->average_down(FD_0, m_realm, m_cdr->get_phase());
	m_amr->interp_ghost(FD_0, m_realm, m_cdr->get_phase());
      }
      else{
	data_ops::set_value(FD_0, 0.0);
      }
    }
  }
}

void imex_sdc::integrate(const Real a_dt, const Real a_time, const bool a_lagged_terms){
  CH_TIME("imex_sdc::integrate");
  if(m_verbosity > 5){
    pout() << "imex_sdc::integrate" << endl;
  }


  // 1. The first time we enter this routine, velocities were updated.
  // 2. For further calls, source terms and velocities have been overwritten, but since the explicit
  //    operator slopes do not change, this is perfectly fine. We just increment with the lagged terms. 
  Real t0, t1;
  Real total_time = 0.0;
  Real setup_time  = 0.0;
  Real advect_time = 0.0;
  Real diffusive_time = 0.0;

  // We begin with phi[0] = phi(t_n). Then update phi[m+1].
  Real time = a_time;
  for(int m = 0; m < m_p; m++){

    // Update source terms every time we go through this
    imex_sdc::compute_E_into_scratch();
    imex_sdc::compute_reaction_network(m, a_time, m_dtm[m]); // Ppdate the CDR and RTE source terms using the correct step size

    // Always update boundary conditions on the way in. All of these calls use the stuff that reside in the solvers,
    // which is what we need to do at the start of the time step. In principle, these things do not change
    // and so we could probably store them somewhere for increased performance. 
    if(m == 0 && !a_lagged_terms){ // This updates the CDR boundary conditions; since these are used to compute the slopes,
      t0 = MPI_Wtime();            // and the m=0 hyperbolic slopes do not change, we only need to do this for the predictor. 
      imex_sdc::compute_cdr_eb_states();
      imex_sdc::compute_cdr_fluxes(a_time);
      imex_sdc::compute_cdr_domain_states();
      imex_sdc::compute_cdr_domain_fluxes(a_time);
      imex_sdc::compute_sigma_flux();
      t1 = MPI_Wtime();

      total_time = -t0;
      setup_time = t1-t0;
    }
    
    // This does the transient rte advance. Source terms were uåpdated in the compute_reaction_network routine above. 
    t0 = MPI_Wtime();
    if(!(m_rte->is_stationary())) imex_sdc::integrate_rte_transient(a_dt);

    // This computes phi_(m+1) = phi_m + dtm*FAR_m(phi_m) + lagged quadrature and lagged advection-reaction
    t0 = MPI_Wtime();
    imex_sdc::integrate_advection_reaction(a_dt, m, a_lagged_terms);
    t1 = MPI_Wtime();
    advect_time += t1-t0;

    // This does the diffusion advance. It also adds in the remaining lagged diffusion terms before the implicit diffusion solve
    t0 = MPI_Wtime();
    imex_sdc::integrate_diffusion(a_dt, m, a_lagged_terms);
    t1 = MPI_Wtime();
    diffusive_time += t1-t0;

    // After the diffusion step we update the Poisson and *stationary* RTE equations
    Vector<EBAMRCellData*> cdr_densities_mp1 = imex_sdc::get_cdr_phik(m+1);
    EBAMRIVData& sigma_mp1 = imex_sdc::get_sigmak(m+1);
    const Real t_mp1 = m_tm[m+1];

    // Update electric field and stationary RTE equations
    if(m_consistent_E)   imex_sdc::update_poisson(cdr_densities_mp1, sigma_mp1);
    if(m_consistent_rte) {
      if(m_rte->is_stationary()){
	imex_sdc::compute_reaction_network(m+1, time + m_dtm[m], m_dtm[m]);
	imex_sdc::integrate_rte_stationary();
      }
    }

    // If we need another step, we should update boundary conditions agains. We DONT do this on the last step
    // because this was also done on the way INTO this routine. If we've updated m=m_p, we either recompute
    // boundary conditions in the next SDC sweep, or we allow the next time step to take care of this. 
    const int last = m == m_p-1;
    if(!last){
      if(m_compute_S)      imex_sdc::compute_cdr_gradients(cdr_densities_mp1);
      if(m_compute_v)      imex_sdc::compute_cdr_velo(cdr_densities_mp1, t_mp1);
      if(m_compute_D)      imex_sdc::update_diffusion_coefficients();

      // Update boundary conditions for cdr and sigma equations. We need them for the next step
      imex_sdc::compute_cdr_eb_states(cdr_densities_mp1);
      imex_sdc::compute_cdr_fluxes(cdr_densities_mp1, t_mp1);
      imex_sdc::compute_cdr_domain_states(cdr_densities_mp1);
      imex_sdc::compute_cdr_domain_fluxes(cdr_densities_mp1, t_mp1);
      imex_sdc::compute_sigma_flux();
    }

    time += m_dtm[m];
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

void imex_sdc::integrate_advection_reaction(const Real a_dt, const int a_m, const bool a_lagged_terms){
  CH_TIME("imex_sdc::integrate_advection_reaction");
  if(m_verbosity > 5){
    pout() << "imex_sdc::integrate_advection_reaction" << endl;
  }

  // Advance phi_(m+1) = phi_m + dtm*F_A. These routines do nothing
  // with the operator slopes for phi, but they do adjust the slopes m_Fsig (but not m_Fsum) for sigma. Incidentally,
  // if m=0 and a_lagged_terms=true, we can increment directly with the precomputed advection-reaction. This means that
  // we can skip the advective advance. The sigma advance is accordingly also skipped.
  const bool skip = (a_m == 0 && a_lagged_terms);
  const Real t0 = MPI_Wtime();
  if(!skip){
    imex_sdc::integrate_advection(a_dt, a_m, a_lagged_terms);
  }
  const Real t1 = MPI_Wtime();

  // Add in the reaction term and then compute the new operator slopes.
  // If this is the corrector and m=0, we skipped the advection advance because we can use the precomputed
  // advection-reaction operator slope. In this case phi_(m+1) is bogus and we need to recompute it. Otherwise,
  // phi_(m+1) = phi_m + dtm*FA_m, and we just increment with the reaction operator. 
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
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
      if(a_lagged_terms) {
	data_ops::copy(scratch, FAR_m);
      }
    }
    else{ // If we made it here, phi_(m+1) = phi_m + dtm*FA(phi_m) through the integrate_advection routine
      EBAMRCellData& FAR_m     = storage->get_FAR()[a_m]; // Currently the old slope
      EBAMRCellData& src = solver->get_source();    // Updated source

      // Increment swith source and then compute slope. This has already been done 
      data_ops::incr(phi_m1, src, m_dtm[a_m]);  // phi_(m+1) = phi_m + dtm*(FA_m + FR_m)

      // This shouldn't be necessary
      m_amr->average_down(phi_m1, m_realm, m_cdr->get_phase());
      m_amr->interp_ghost(phi_m1, m_realm, m_cdr->get_phase());

      if(a_lagged_terms){ // Back up the old slope first, we will need it for the lagged term
	data_ops::copy(scratch, FAR_m);
      }

      // Re-compute the advection-reaction slope for node t_m
      data_ops::copy(FAR_m, phi_m1);            // FAR_m = (phi_(m+1) - phi_m)/dtm
      data_ops::incr(FAR_m, phi_m, -1.0);       // :
      data_ops::scale(FAR_m, 1./m_dtm[a_m]);    // :

      // Shouldn't be necessary
      m_amr->average_down(FAR_m, m_realm, m_cdr->get_phase());
      m_amr->interp_ghost(FAR_m, m_realm, m_cdr->get_phase());
    }

    // Now add in the lagged advection-reaction and quadrature terms. This is a bit weird, but we did overwrite
    // FAR_m above after the advection-reaction advance, but we also backed up the old term into scratch. 
    if(a_lagged_terms){
      data_ops::incr(phi_m1, scratch, -m_dtm[a_m]); // phi_(m+1)^(k+1) = phi_m^(k+1) + dtm*(FAR_m^(k+1) - FAR_m^k)
      imex_sdc::quad(scratch, storage->get_F(), a_m);  // Does the quadrature of the lagged operator slopes. 
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
  if(a_lagged_terms){ // Add in the lagged terms. When we make it here, sigma_(m+1) = sigma_m + dtm*Fsig_m. 
    EBAMRIVData& Fsig_lag = m_sigma_scratch->get_Fold()[a_m];
    data_ops::incr(sigma_m1, Fsig_lag, -m_dtm[a_m]);

    // Add in the quadrature term
    EBAMRIVData& scratch = m_sigma_scratch->get_scratch();
    imex_sdc::quad(scratch, m_sigma_scratch->get_Fold(), a_m);
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

void imex_sdc::integrate_advection(const Real a_dt, const int a_m, const bool a_lagged_terms){
  CH_TIME("imex_sdc::integrate_advection");
  if(m_verbosity > 5){
    pout() << "imex_sdc::integrate_advection" << endl;
  }

  // TLDR; This routine should do phi_(m+1) = phi_m + dtm*FA_m, and sigma_(m+1) = sigma_m + dt*Fsig_m.
  //       It also computes the sigma slope.
  //
  //       The lagged terms are not a part of this routine. 

  if(a_m == 0 && a_lagged_terms){
    MayDay::Abort("imex_sdc::integrate_advection - (m==0 && corrector==true) logic bust which should never happen");
  }

  // Advance cdr equations
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);

    EBAMRCellData& phi_m1  = storage->get_phi()[a_m+1];
    EBAMRCellData& scratch = storage->get_scratch();
    EBAMRCellData& phi_m   = storage->get_phi()[a_m];

    if(solver->is_mobile()){
      const Real extrap_dt = m_extrap_advect ? 2.0*m_extrap_dt*m_dtm[a_m] : 0.0; // Factor of 2 due to EBPatchAdvect
      solver->compute_divF(scratch, phi_m, extrap_dt);                           // scratch =  Div(v_m*phi_m^(k+1))
    
      data_ops::copy(phi_m1, phi_m);
      data_ops::incr(phi_m1, scratch, -m_dtm[a_m]);
      data_ops::floor(phi_m1, 0.0);
      m_amr->average_down(phi_m1, m_realm, m_cdr->get_phase());
      m_amr->interp_ghost(phi_m1, m_realm, m_cdr->get_phase());

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

void imex_sdc::integrate_diffusion(const Real a_dt, const int a_m, const bool a_lagged_terms){
  CH_TIME("imex_sdc::integrate_diffusion");
  if(m_verbosity > 5){
    pout() << "imex_sdc::integrate_diffusion" << endl;
  }

  // TLDR: We're solving
  //
  // phi_(m+1)^(k+1) = phi_(m)^(k+1,\ast) + dtm*FD_(m+1)^(k+1) + sources. 
  //
  // This routine does not modify FD_(m+1)^k. This is replaced by FD_(m+1)^(k+1) later on. 
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = imex_sdc::get_cdr_storage(solver_it);
    
    if(solver->is_diffusive()){
      EBAMRCellData& phi_m1      = storage->get_phi()[a_m+1]; // Advected solution. Possibly with lagged terms. 
      const EBAMRCellData& phi_m = storage->get_phi()[a_m];

      // Build the diffusion source term
      EBAMRCellData& source   = storage->get_scratch();
      EBAMRCellData& init_soln = storage->get_scratch2();
      data_ops::set_value(source, 0.0); // No source term
      
      data_ops::copy(init_soln, phi_m1);      // Copy initial solutions
      if(a_lagged_terms){
	const EBAMRCellData& FD_m1k = storage->get_FD()[a_m+1];      // FD_(m+1)^k. Lagged term.
	data_ops::incr(init_soln, FD_m1k, -m_dtm[a_m]);
      }
      m_amr->average_down(init_soln, m_realm, m_cdr->get_phase());
      m_amr->interp_ghost(init_soln, m_realm, m_cdr->get_phase());
      data_ops::copy(phi_m1, phi_m);

      // Solve
      if(m_use_tga){
	solver->advance_tga(phi_m1, init_soln, source, m_dtm[a_m]); // No source. 
      }
      else{
	solver->advance_euler(phi_m1, init_soln, source, m_dtm[a_m]); // No source. 
      }
      m_amr->average_down(phi_m1, m_realm, m_cdr->get_phase());
      m_amr->interp_ghost(phi_m1, m_realm, m_cdr->get_phase());
      data_ops::floor(phi_m1, 0.0);

      // Update the operator slope
      EBAMRCellData& FD_m1k = storage->get_FD()[a_m+1];
      data_ops::set_value(FD_m1k, 0.0);
      data_ops::incr(FD_m1k, phi_m1, 1.0);
      data_ops::incr(FD_m1k, init_soln, -1.0);
      data_ops::scale(FD_m1k, 1./m_dtm[a_m]);

      m_amr->average_down(FD_m1k, m_realm, m_cdr->get_phase());
      m_amr->interp_ghost(FD_m1k, m_realm, m_cdr->get_phase());
    }
    else{
      EBAMRCellData& FD_m1k = storage->get_FD()[a_m+1];
      data_ops::set_value(FD_m1k, 0.0);
    }
  }
}

void imex_sdc::reconcile_integrands(){
  CH_TIME("imex_sdc::reconcile_integrands");
  if(m_verbosity > 5){
    pout() << "imex_sdc::reconcile_integrands" << endl;
  }

  // TLDR: When we come in here, all solutions (Poisson, CDR, RTE, Sigma) are known at node m_p. But we
  //       do need the extra slopes for the explicit operators

  Vector<EBAMRCellData*> cdr_densities_p = imex_sdc::get_cdr_phik(m_p);
  EBAMRIVData& sigma_p = imex_sdc::get_sigmak(m_p);
  const Real t_p = m_tm[m_p];

  // Update boundary conditions for cdr and sigma equations before getting the slope at the final node
  imex_sdc::compute_cdr_eb_states(cdr_densities_p);
  imex_sdc::compute_cdr_fluxes(cdr_densities_p, t_p);
  imex_sdc::compute_cdr_domain_states(cdr_densities_p);
  imex_sdc::compute_cdr_domain_fluxes(cdr_densities_p, t_p);
  imex_sdc::compute_sigma_flux();

  // Now compute FAR_p - that wasn't done when we integrated
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = imex_sdc::get_cdr_storage(solver_it);
    const int idx = solver_it.index();

    // This has not been computed yet. Do it.
    EBAMRCellData& FAR_p     = storage->get_FAR()[m_p];
    EBAMRCellData& phi_p     = *cdr_densities_p[idx] ;
    const EBAMRCellData& src = solver->get_source();

    if(solver->is_mobile()){
      const Real extrap_dt = m_extrap_advect ? 2.0*m_extrap_dt*m_dtm[m_p-1] : 0.0; // Factor of 2 because of EBPatchAdvect
      solver->compute_divF(FAR_p, phi_p, extrap_dt);       // FAR_p =  Div(v_p*phi_p)
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
      m_amr->average_down(F_m, m_realm, m_cdr->get_phase());
      m_amr->interp_ghost(F_m, m_realm, m_cdr->get_phase());
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

void imex_sdc::initialize_errors(){
  CH_TIME("imex_sdc::corrector_initialize_errors");
  if(m_verbosity > 5){
    pout() << "imex_sdc::corrector_initialize_errors" << endl;
  }

  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = imex_sdc::get_cdr_storage(solver_it);
    const int idx = solver_it.index();
    
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

void imex_sdc::finalize_errors(){
  CH_TIME("imex_sdc::corrector_finalize_errors");
  if(m_verbosity > 5){
    pout() << "imex_sdc::corrector_finalize_errors" << endl;
  }

  const Real safety = 1.E-20;

  m_max_error = 0.0;
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = imex_sdc::get_cdr_storage(solver_it);
    const int idx = solver_it.index();

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

void imex_sdc::compute_new_dt(bool& a_accept_step, const Real a_dt, const int a_num_corrections){
  CH_TIME("imex_sdc::compute_new_dt");
  if(m_verbosity > 5){
    pout() << "imex_sdc::compute_new_dt" << endl;
  }

  // If a_dt was the smallest possible CFL or hardcap time step, we just have to accept it
  const Real max_gl_dist = imex_sdc::get_max_node_distance();
  Real dt_cfl = 2.0*m_dt_cfl/max_gl_dist; // This is the smallest time step ON THE FINEST LEVEL

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

    m_new_dt = Max(m_new_dt, dt_cfl*m_minCFL);            // Don't drop below minimum CFL
    m_new_dt = Min(m_new_dt, dt_cfl*m_maxCFL);            // Don't go above maximum CFL

  }
  else{
    a_accept_step = false;

    m_new_dt = m_decrease_safe*dt_adapt; // Decrease time step a little bit extra to avoid another rejection
    if(a_dt <= min_dt_cfl || a_dt < m_min_dt){ // Step already at minimum. Accept it anyways.
      a_accept_step = true;
    }
    
    m_new_dt = Max(m_new_dt, dt_cfl*m_minCFL);            // Don't drop below minimum CFL
    m_new_dt = Min(m_new_dt, dt_cfl*m_maxCFL);            // Don't go above maximum CFL

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

void imex_sdc::adaptive_report(const Real a_first_dt, const Real a_dt, const Real a_new_dt, const int a_corr, const int a_rej, const Real a_max_err){
  CH_TIME("imex_sdc::adaptive_report");
  if(m_verbosity > 5){
    pout() << "imex_sdc::adaptive_report" << endl;
  }

  pout() << "\n";
  pout() << "imex_sdc::adaptive_report breakdown" << endl;
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

void imex_sdc::compute_dt(Real& a_dt, time_code& a_timecode){
  CH_TIME("imex_sdc::compute_dt");
  if(m_verbosity > 5){
    pout() << "imex_sdc::compute_dt" << endl;
  }

  Real dt = 1.E99;

  int Nref = 1;
  for (int lvl = 0; lvl < m_amr->get_finest_level(); lvl++){
    Nref = Nref*m_amr->get_ref_rat()[lvl];
  }
  const Real max_gl_dist = imex_sdc::get_max_node_distance();
  m_dt_cfl = m_cdr->compute_advection_dt();

  Real dt_cfl = 2.0*m_dt_cfl/max_gl_dist;
  
  // Time step selection for non-adaptive stepping
  if(!m_adaptive_dt){
    if(dt_cfl < dt){
      dt = m_cfl*dt_cfl;
      a_timecode = time_code::advection;
    }
  }
  else{
    Real new_dt;

    // Step should not exceed m_new_dt. Also, it shoul
    if(m_have_dt_err){
      new_dt = m_new_dt;

      new_dt = Max(new_dt, dt_cfl*m_minCFL);
      new_dt = Min(new_dt, dt_cfl*m_maxCFL);
    }
    else{
      new_dt = m_minCFL*dt_cfl;
    }

    if(new_dt < dt){
      dt = new_dt;
      a_timecode = time_code::error;
    }
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

void imex_sdc::regrid_internals(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
  CH_TIME("imex_sdc::regrid_internals");
  if(m_verbosity > 5){
    pout() << "imex_sdc::regrid_internals" << endl;
  }

  // Nothing to see here
}

void imex_sdc::allocate_internals(){
  CH_TIME("imex_sdc::allocate_internals");
  if(m_verbosity > 5){
    pout() << "imex_sdc::allocate_internals" << endl;
  }

  m_cdr_error.resize(m_physics->get_num_cdr_species());
  
  imex_sdc::allocate_cdr_storage();
  imex_sdc::allocate_poisson_storage();
  imex_sdc::allocate_rte_storage();
  imex_sdc::allocate_sigma_storage();

  imex_sdc::setup_quadrature_nodes(m_p);
  imex_sdc::setup_qmj(m_p);
}

void imex_sdc::allocate_cdr_storage(){
  const int ncomp       = 1;
  const int num_species = m_physics->get_num_cdr_species();

  m_cdr_scratch.resize(num_species);
  
  for (cdr_iterator<cdr_solver> solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    m_cdr_scratch[idx] = RefCountedPtr<cdr_storage> (new cdr_storage(m_amr, m_realm, m_cdr->get_phase(), ncomp));
    m_cdr_scratch[idx]->allocate_storage(m_p);
  }
}

void imex_sdc::allocate_poisson_storage(){
  const int ncomp = 1;
  m_poisson_scratch = RefCountedPtr<poisson_storage> (new poisson_storage(m_amr, m_realm, m_cdr->get_phase(), ncomp));
  m_poisson_scratch->allocate_storage(m_p);
}

void imex_sdc::allocate_rte_storage(){
  const int ncomp       = 1;
  const int num_photons = m_physics->get_num_rte_species();
  m_rte_scratch.resize(num_photons);
  
  for (rte_iterator<rte_solver> solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    m_rte_scratch[idx] = RefCountedPtr<rte_storage> (new rte_storage(m_amr, m_realm, m_rte->get_phase(), ncomp));
    m_rte_scratch[idx]->allocate_storage(m_p);
  }
}

void imex_sdc::allocate_sigma_storage(){
  const int ncomp = 1;
  m_sigma_scratch = RefCountedPtr<sigma_storage> (new sigma_storage(m_amr, m_realm, m_cdr->get_phase(), ncomp));
  m_sigma_scratch->allocate_storage(m_p);
}

void imex_sdc::deallocate_internals(){
  CH_TIME("imex_sdc::deallocate_internals");
  if(m_verbosity > 5){
    pout() << "imex_sdc::deallocate_internals" << endl;
  }

  for (cdr_iterator<cdr_solver> solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    m_cdr_scratch[idx]->deallocate_storage();
    m_cdr_scratch[idx] = RefCountedPtr<cdr_storage>(0);
  }

  for (rte_iterator<rte_solver> solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
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

void imex_sdc::compute_E_into_scratch(){
  CH_TIME("imex_sdc::compute_E_into_scratch");
  if(m_verbosity > 5){
    pout() << "imex_sdc::compute_E_into_scratch" << endl;
  }
  
  EBAMRCellData& E_cell = m_poisson_scratch->get_E_cell();
  EBAMRFluxData& E_face = m_poisson_scratch->get_E_face();
  EBAMRIVData&   E_eb   = m_poisson_scratch->get_E_eb();
  EBAMRIFData&   E_dom  = m_poisson_scratch->get_E_domain();

  const MFAMRCellData& phi = m_poisson->get_state();
  
  imex_sdc::compute_E(E_cell, m_cdr->get_phase(), phi);     // Compute cell-centered field
  imex_sdc::compute_E(E_face, m_cdr->get_phase(), E_cell);  // Compute face-centered field
  imex_sdc::compute_E(E_eb,   m_cdr->get_phase(), E_cell);  // EB-centered field

  cdr_plasma_stepper::extrapolate_to_domain_faces(E_dom, m_cdr->get_phase(), E_cell);
}

void imex_sdc::compute_cdr_gradients(){
  CH_TIME("imex_sdc::compute_cdr_gradients");
  if(m_verbosity > 5){
    pout() << "imex_sdc::compute_cdr_gradients" << endl;
  }

  imex_sdc::compute_cdr_gradients(m_cdr->get_states());
}

void imex_sdc::compute_cdr_gradients(const Vector<EBAMRCellData*>& a_states){
  CH_TIME("imex_sdc::compute_cdr_gradients");
  if(m_verbosity > 5){
    pout() << "imex_sdc::compute_cdr_gradients" << endl;
  }

  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    RefCountedPtr<cdr_storage>& storage = imex_sdc::get_cdr_storage(solver_it);
    EBAMRCellData& grad = storage->get_gradient();
    m_amr->compute_gradient(grad, *a_states[idx], m_realm, m_cdr->get_phase());
    //    m_amr->average_down(grad, m_realm, m_cdr->get_phase());
    m_amr->interp_ghost(grad, m_realm, m_cdr->get_phase());
  }
}

void imex_sdc::compute_cdr_velo(const Real a_time){
  CH_TIME("imex_sdc::compute_cdr_velo");
  if(m_verbosity > 5){
    pout() << "imex_sdc::compute_cdr_velo" << endl;
  }

  imex_sdc::compute_cdr_velo(m_cdr->get_states(), a_time);
}

void imex_sdc::compute_cdr_velo(const Vector<EBAMRCellData*>& a_states, const Real a_time){
  CH_TIME("imex_sdc::compute_cdr_velo(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "imex_sdc::compute_cdr_velo(Vector<EBAMRCellData*>, Real)" << endl;
  }

  Vector<EBAMRCellData*> velocities = m_cdr->get_velocities();
  imex_sdc::compute_cdr_velocities(velocities, a_states, m_poisson_scratch->get_E_cell(), a_time);
}

void imex_sdc::compute_cdr_eb_states(){
  CH_TIME("imex_sdc::compute_cdr_eb_states");
  if(m_verbosity > 5){
    pout() << "imex_sdc::compute_cdr_eb_states" << endl;
  }

  Vector<EBAMRIVData*>   eb_gradients;
  Vector<EBAMRIVData*>   eb_states;
  Vector<EBAMRCellData*> cdr_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    RefCountedPtr<cdr_storage>& storage = imex_sdc::get_cdr_storage(solver_it);

    cdr_states.push_back(&(solver->get_state()));
    eb_states.push_back(&(storage->get_eb_state()));
    eb_gradients.push_back(&(storage->get_eb_grad()));
    cdr_gradients.push_back(&(storage->get_gradient())); // Should already have been computed
  }

  // Extrapolate states to the EB and floor them so we cannot get negative values on the boundary. This
  // won't hurt mass conservation because the mass hasn't been injected yet
  imex_sdc::extrapolate_to_eb(eb_states, m_cdr->get_phase(), cdr_states);
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    data_ops::floor(*eb_states[idx], 0.0);
  }

  // We should already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIVData eb_gradient;
  m_amr->allocate(eb_gradient, m_realm, m_cdr->get_phase(), SpaceDim);
  for (int i = 0; i < cdr_states.size(); i++){
    imex_sdc::extrapolate_to_eb(eb_gradient, m_cdr->get_phase(), *cdr_gradients[i]);
    imex_sdc::project_flux(*eb_gradients[i], eb_gradient);
  }
}

void imex_sdc::compute_cdr_eb_states(const Vector<EBAMRCellData*>& a_states){
  CH_TIME("imex_sdc::compute_cdr_eb_states(vec)");
  if(m_verbosity > 5){
    pout() << "imex_sdc::compute_cdr_eb_states(vec)" << endl;
  }

  Vector<EBAMRIVData*>   eb_gradients;
  Vector<EBAMRIVData*>   eb_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = imex_sdc::get_cdr_storage(solver_it);

    eb_states.push_back(&(storage->get_eb_state()));
    eb_gradients.push_back(&(storage->get_eb_grad()));
    cdr_gradients.push_back(&(storage->get_gradient())); // Should already have been computed
  }

  // Extrapolate states to the EB and floor them so we cannot get negative values on the boundary. This
  // won't hurt mass conservation because the mass hasn't been injected yet
  imex_sdc::extrapolate_to_eb(eb_states, m_cdr->get_phase(), a_states);
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    data_ops::floor(*eb_states[idx], 0.0);
  }

  // We should already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIVData eb_gradient;
  m_amr->allocate(eb_gradient, m_realm, m_cdr->get_phase(), SpaceDim);
  for (int i = 0; i < a_states.size(); i++){
    imex_sdc::extrapolate_to_eb(eb_gradient, m_cdr->get_phase(), *cdr_gradients[i]);
    imex_sdc::project_flux(*eb_gradients[i], eb_gradient);
  }
}

void imex_sdc::compute_cdr_domain_states(){
  CH_TIME("imex_sdc::compute_cdr_domain_states");
  if(m_verbosity > 5){
    pout() << "imex_sdc::compute_cdr_domain_states" << endl;
  }

  Vector<EBAMRIFData*>   domain_gradients;
  Vector<EBAMRIFData*>   domain_states;
  Vector<EBAMRCellData*> cdr_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    RefCountedPtr<cdr_storage>& storage = imex_sdc::get_cdr_storage(solver_it);

    cdr_states.push_back(&(solver->get_state()));
    domain_states.push_back(&(storage->get_domain_state()));
    domain_gradients.push_back(&(storage->get_domain_grad()));
    cdr_gradients.push_back(&(storage->get_gradient())); // Should already be computed
  }

  // Extrapolate states to the domain faces
  imex_sdc::extrapolate_to_domain_faces(domain_states, m_cdr->get_phase(), cdr_states);

  // We already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIFData grad;
  m_amr->allocate(grad, m_realm, m_cdr->get_phase(), SpaceDim);
  for (int i = 0; i < cdr_states.size(); i++){
    imex_sdc::extrapolate_to_domain_faces(grad, m_cdr->get_phase(), *cdr_gradients[i]);
    imex_sdc::project_domain(*domain_gradients[i], grad);
  }
}

void imex_sdc::compute_cdr_domain_states(const Vector<EBAMRCellData*>& a_states){
  CH_TIME("imex_sdc::compute_cdr_domain_states");
  if(m_verbosity > 5){
    pout() << "imex_sdc::compute_cdr_domain_states" << endl;
  }

  Vector<EBAMRIFData*>   domain_gradients;
  Vector<EBAMRIFData*>   domain_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
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
  m_amr->allocate(grad, m_realm, m_cdr->get_phase(), SpaceDim);
  for (int i = 0; i < a_states.size(); i++){
    this->extrapolate_to_domain_faces(grad, m_cdr->get_phase(), *cdr_gradients[i]);
    this->project_domain(*domain_gradients[i], grad);
  }
}

void imex_sdc::compute_cdr_fluxes(const Real a_time){
  CH_TIME("imex_sdc::compute_cdr_fluxes");
  if(m_verbosity > 5){
    pout() << "imex_sdc::compute_cdr_fluxes" << endl;
  }

  this->compute_cdr_fluxes(m_cdr->get_states(), a_time);
}

void imex_sdc::compute_cdr_fluxes(const Vector<EBAMRCellData*>& a_states, const Real a_time){
  CH_TIME("imex_sdc::compute_cdr_fluxes(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "imex_sdc::compute_cdr_fluxes(Vector<EBAMRCellData*>, Real)" << endl;
  }

  Vector<EBAMRIVData*> cdr_fluxes;
  Vector<EBAMRIVData*> extrap_cdr_fluxes;
  Vector<EBAMRIVData*> extrap_cdr_densities;
  Vector<EBAMRIVData*> extrap_cdr_velocities;
  Vector<EBAMRIVData*> extrap_cdr_gradients;
  Vector<EBAMRIVData*> extrap_rte_fluxes;

  cdr_fluxes = m_cdr->get_ebflux();

  for (cdr_iterator<cdr_solver> solver_it(*m_cdr); solver_it.ok(); ++solver_it){
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
  cdr_plasma_stepper::compute_extrapolated_fluxes(extrap_cdr_fluxes, a_states, cdr_velocities, m_cdr->get_phase());
  cdr_plasma_stepper::compute_extrapolated_velocities(extrap_cdr_velocities, cdr_velocities, m_cdr->get_phase());

  // Compute RTE flux on the boundary
  for (rte_iterator<rte_solver> solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);

    EBAMRIVData& flux_eb = storage->get_eb_flux();
    solver->compute_boundary_flux(flux_eb, solver->get_state());
    extrap_rte_fluxes.push_back(&flux_eb);
  }

  const EBAMRIVData& E = m_poisson_scratch->get_E_eb();
  cdr_plasma_stepper::compute_cdr_fluxes(cdr_fluxes,
				   extrap_cdr_fluxes,
				   extrap_cdr_densities,
				   extrap_cdr_velocities,
				   extrap_cdr_gradients,
				   extrap_rte_fluxes,
				   E,
				   a_time);
}

void imex_sdc::compute_cdr_domain_fluxes(const Real a_time){
  CH_TIME("imex_sdc::compute_cdr_domain_fluxes");
  if(m_verbosity > 5){
    pout() << "imex_sdc::compute_cdr_domain_fluxes" << endl;
  }

  this->compute_cdr_domain_fluxes(m_cdr->get_states(), a_time);
}

void imex_sdc::compute_cdr_domain_fluxes(const Vector<EBAMRCellData*>& a_states, const Real a_time){
  CH_TIME("imex_sdc::compute_cdr_domain_fluxes(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "imex_sdc::compute_cdr_domain_fluxes(Vector<EBAMRCellData*>, Real)" << endl;
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
  for (cdr_iterator<cdr_solver> solver_it(*m_cdr); solver_it.ok(); ++solver_it){
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
  this->extrapolate_velo_to_domain_faces(extrap_cdr_velocities, m_cdr->get_phase(), cdr_velocities);
  this->compute_extrapolated_domain_fluxes(extrap_cdr_fluxes,     a_states,           cdr_velocities, m_cdr->get_phase());
  this->extrapolate_vector_to_domain_faces(extrap_cdr_gradients,  m_cdr->get_phase(), cdr_gradients);

  // Compute RTE flux on domain faces
  for (rte_iterator<rte_solver> solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);

    EBAMRIFData& domain_flux = storage->get_domain_flux();
    solver->compute_domain_flux(domain_flux, solver->get_state());
    extrap_rte_fluxes.push_back(&domain_flux);
  }

  const EBAMRIFData& E = m_poisson_scratch->get_E_domain();

  // This fills the solvers' domain fluxes
  cdr_plasma_stepper::compute_cdr_domain_fluxes(cdr_fluxes,
					  extrap_cdr_fluxes,
					  extrap_cdr_densities,
					  extrap_cdr_velocities,
					  extrap_cdr_gradients,
					  extrap_rte_fluxes,
					  E,
					  a_time);
}

void imex_sdc::compute_sigma_flux(){
  CH_TIME("imex_sdc::compute_sigma_flux");
  if(m_verbosity > 5){
    pout() << "imex_sdc::compute_sigma_flux" << endl;
  }

  EBAMRIVData& flux = m_sigma->get_flux();
  data_ops::set_value(flux, 0.0);

  for (cdr_iterator<cdr_solver> solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    const RefCountedPtr<cdr_species>& spec  = solver_it.get_species();
    const EBAMRIVData& solver_flux          = solver->get_ebflux();

    data_ops::incr(flux, solver_flux, spec->get_charge()*units::s_Qe);
  }

  m_sigma->reset_cells(flux);
}

void imex_sdc::compute_reaction_network(const int a_m, const Real a_time, const Real a_dt){
  CH_TIME("imex_sdc::compute_reaction_network");
  if(m_verbosity > 5){
    pout() << "imex_sdc::compute_reaction_network";
  }

  Vector<EBAMRCellData*> cdr_sources = m_cdr->get_sources();
  Vector<EBAMRCellData*> rte_sources = m_rte->get_sources();

  const Vector<EBAMRCellData*> cdr_densities = get_cdr_phik(a_m);
  const Vector<EBAMRCellData*> rte_densities = m_rte->get_states();
  const EBAMRCellData& E = m_poisson_scratch->get_E_cell();

  cdr_plasma_stepper::advance_reaction_network(cdr_sources, rte_sources, cdr_densities, rte_densities, E, a_time, a_dt);
}

void imex_sdc::update_poisson(){
  CH_TIME("imex_sdc::update_poisson(solver)");
  if(m_verbosity > 5){
    pout() << "imex_sdc::update_poisson(solver)" << endl;
  }
  
  if(m_do_poisson){ // Solve Poisson equation
    if((m_step +1) % m_fast_poisson == 0){
      cdr_plasma_stepper::solve_poisson();
      this->compute_E_into_scratch();
    }
  }
}

void imex_sdc::update_poisson(const Vector<EBAMRCellData*>& a_densities, const EBAMRIVData& a_sigma){
  CH_TIME("imex_sdc::update_poisson(full)");
  if(m_verbosity > 5){
    pout() << "imex_sdc::update_poisson(full)" << endl;
  }
  
  if(m_do_poisson){ // Solve Poisson equation
    if((m_step +1) % m_fast_poisson == 0){
      cdr_plasma_stepper::solve_poisson(m_poisson->get_state(),
				  m_poisson->get_source(),
				  a_densities,
				  a_sigma,
				  centering::cell_center);
      this->compute_E_into_scratch();
    }
  }
}

void imex_sdc::integrate_rte_transient(const Real a_dt){
  CH_TIME("imex_sdc::integrate_rte_transient");
  if(m_verbosity > 5){
    pout() << "imex_sdc::integrate_rte_transient" << endl;
  }

  if(m_do_rte){
    if((m_step + 1) % m_fast_rte == 0){
      if(!(m_rte->is_stationary())){
	for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
	  RefCountedPtr<rte_solver>& solver = solver_it();
	  solver->advance(a_dt);
	}
      }
    }
  }
}

void imex_sdc::integrate_rte_stationary(){
  CH_TIME("imex_sdc::integrate_rte_transient");
  if(m_verbosity > 5){
    pout() << "imex_sdc::integrate_rte_transient" << endl;
  }

  if(m_do_rte){
    if((m_step + 1) % m_fast_rte == 0){
      if((m_rte->is_stationary())){
	for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
	  RefCountedPtr<rte_solver>& solver = solver_it();
	  solver->advance(0.0);
	}
      }
    }
  }
}

void imex_sdc::update_diffusion_coefficients(){
  CH_TIME("imex_sdc::update_diffusion_coefficients");
  if(m_verbosity > 5){
    pout() << "imex_sdc::update_diffusion_coefficients" << endl;
  }
  cdr_plasma_stepper::compute_cdr_diffusion(m_poisson_scratch->get_E_cell(), m_poisson_scratch->get_E_eb());
}

Vector<EBAMRCellData*> imex_sdc::get_cdr_errors(){
  CH_TIME("imex_sdc::get_cdr_errors");
  if(m_verbosity > 5){
    pout() << "imex_sdc::get_cdr_errors" << endl;
  }
  
  Vector<EBAMRCellData*> ret;
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = imex_sdc::get_cdr_storage(solver_it);
    ret.push_back(&(storage->get_error()));
  }

  return ret;
}

Vector<EBAMRCellData*> imex_sdc::get_cdr_phik(const int a_m){
  CH_TIME("imex_sdc::get_cdr_phik");
  if(m_verbosity > 5){
    pout() << "imex_sdc::get_cdr_phik" << endl;
  }
  
  Vector<EBAMRCellData*> ret;
  for (cdr_iterator<cdr_solver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = imex_sdc::get_cdr_storage(solver_it);
    ret.push_back(&(storage->get_phi()[a_m]));
  }

  return ret;
}

EBAMRIVData& imex_sdc::get_sigmak(const int a_m){
  CH_TIME("imex_sdc::get_sigmak");
  if(m_verbosity > 5){
    pout() << "imex_sdc::get_sigmak)" << endl;
  }
  return m_sigma_scratch->get_sigma()[a_m];
}

void imex_sdc::write_step_profile(const Real a_dt,
			       const Real a_error,
			       const int  a_substeps,
			       const int  a_corrections,
			       const int  a_rejections){
  CH_TIME("sissdc::write_step_profile");
  if(m_verbosity > 5){
    pout() << "imex_sdc::write_step_profile" << endl;
  }

  if(procID() == 0 ){

    const std::string fname("imex_sdc_step_profile.txt");
    
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

void imex_sdc::store_solvers(){
  CH_TIME("imex_sdc::store_solvers");
  if(m_verbosity > 5){
    pout() << "imex_sdc::store_solvers" << endl;
  }

  if(m_k > 0 && m_adaptive_dt){
    // IMEX_SDC does not manipulate cdr and sigma solvers until the end of the time step. Only need to do
    // Poisson and RTE here.

    // Poisson
    MFAMRCellData& previous    = m_poisson_scratch->get_previous();
    const MFAMRCellData& state = m_poisson->get_state();
    data_ops::copy(previous, state);

    // RTE
    for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<rte_storage>& storage     = imex_sdc::get_rte_storage(solver_it);
      const RefCountedPtr<rte_solver>& solver = solver_it();

      EBAMRCellData& previous = storage->get_previous();
      const EBAMRCellData& state = solver->get_state();

      data_ops::copy(previous, state);
    }
  }
}

void imex_sdc::restore_solvers(){
  CH_TIME("imex_sdc::restore_solvers");
  if(m_verbosity > 5){
    pout() << "imex_sdc::restore_solvers" << endl;
  }

  // IMEX_SDC does not manipulate cdr and sigma solvers until the end of the time step. Only need to do
  // Poisson and RTE here. 

  // Poisson
  MFAMRCellData& state = m_poisson->get_state();
  const MFAMRCellData& previous    = m_poisson_scratch->get_previous();

  data_ops::copy(state, previous);

  // RTE
  for (rte_iterator<rte_solver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_storage>& storage     = imex_sdc::get_rte_storage(solver_it);
    RefCountedPtr<rte_solver>& solver = solver_it();

    EBAMRCellData& previous = storage->get_previous();
    EBAMRCellData& state = solver->get_state();

    data_ops::copy(state, previous);
  }
}
