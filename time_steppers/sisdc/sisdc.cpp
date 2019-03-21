/*!
  @file   sisdc.cpp
  @brief  Implementation of sisdc.H
  @author Robert Marskar
n  @date   Feb. 2019
*/

#include "sisdc.H"
#include "sisdcF_F.H"
#include "sisdc_storage.H"
#include "data_ops.H"
#include "units.H"
#include "cdr_tga.H"
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
  m_p             = 2;
  m_k             = 1;
  m_error_norm    = 2;
  m_minCFL        = 0.1;
  m_maxCFL        = 0.9;   
  m_err_thresh    = 1.E-4;
  m_safety        = 0.9;
  m_num_diff_corr = 0;
  m_new_dt        = 1.234567E89;
  m_max_tries     = 1;
  m_min_corr      = 1;
  m_upwd_idx      = 0;
  m_extrap_dt     = 0.5;

  m_which_nodes   = "lobatto";

  m_extrap_advect  = false;
  m_subcycle       = false;
  m_print_report   = false;
  m_adaptive_dt    = false;
  m_strong_diffu   = false;
  m_have_dt_err    = false;

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

  // Profiling
  m_profile_steps  = false;

  // Get parameters from input script
  {
    ParmParse pp("sisdc");

    std::string str;

    pp.query("subintervals",    m_p);
    pp.query("corr_iter",       m_k);
    pp.query("error_norm",      m_error_norm);
    pp.query("min_corr",        m_min_corr);
    pp.query("min_cfl",         m_minCFL);
    pp.query("max_cfl",         m_maxCFL);
    pp.query("max_error",       m_err_thresh);
    pp.query("safety",          m_safety);
    pp.query("num_corrections", m_num_diff_corr);
    pp.query("max_tries",       m_max_tries);
    pp.query("extrap_dt",       m_extrap_dt);

    if(pp.contains("subcycle")){
      pp.get("subcycle", str);
      if(str == "true"){
	m_subcycle = true;
      }
      else{
	m_subcycle = false;
      }
    }
    if(pp.contains("extrap_advect")){
      pp.get("extrap_advect", str);
      if(str == "true"){
	m_extrap_advect = true;
      }
      else{
	m_extrap_advect = false;
      }
    }
    if(pp.contains("quad_nodes")){
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
    }
    if(pp.contains("print_report")){
      pp.get("print_report", str);
      if(str == "true"){
	m_print_report = true;
      }
    }
    if(pp.contains("adaptive_dt")){
      pp.get("adaptive_dt", str);
      if(str == "true"){
	m_adaptive_dt = true;
      }
    }
    if(pp.contains("diffusive_coupling")){
      pp.get("diffusive_coupling", str);
      if(str == "strong"){
	m_strong_diffu = true;
      }
      if(str == "weak"){
	m_strong_diffu = false;
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
    if(pp.contains("profile_steps")){
      pp.get("profile_steps", str);
      if(str == "true"){
	m_profile_steps = true;
      }
      else {
	m_profile_steps = false;
      }
    }
  }

  // Setup nodes
  sisdc::setup_quadrature_nodes(m_p);
  sisdc::setup_qmj(m_p);

  if(m_subcycle){
    pout() << "sisdc::sisdc - subcycling in time not (yet) supported. Turning it off" << endl;
    m_subcycle = false;
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

  // SISDC advance
  Real first_dt       = a_dt;
  Real actual_dt      = a_dt;
  int num_reject      = 0;
  int num_corrections = 0;
  bool accept_step    = false;
  bool retry_step     = true;
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
      // This restricts new_dt to > min_cfl and > min_hardcap. If actual_dt is equal these bounds, we accept the step
      sisdc::compute_new_dt(accept_step, actual_dt, num_corrections);
      
      if(!accept_step){  // Step rejection, use the new dt for next step. 
	actual_dt = m_new_dt;
	num_reject++;

	retry_step  = num_reject <= m_max_tries;

	if(retry_step){
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
  sisdc::update_poisson();
  sisdc::update_rte(m_time + actual_dt);

  // Always recompute source terms and velocities for the next time step. These were computed ea
  sisdc::compute_cdr_gradients();
  sisdc::compute_cdr_velo(m_time + actual_dt);
  time_stepper::compute_cdr_diffusion(m_poisson_scratch->get_E_cell(), m_poisson_scratch->get_E_eb());
  sisdc::compute_cdr_sources(m_time + actual_dt);

  // Profile step
  if(m_print_report)  sisdc::adaptive_report(first_dt, actual_dt, m_new_dt, num_corrections, num_reject, m_max_error);
  if(m_profile_steps) sisdc::write_step_profile(actual_dt, m_max_error, m_p, num_corrections, num_reject);

  
  return actual_dt;
}

void sisdc::copy_cdr_to_phi_m0(){
  CH_TIME("sisdc::copy_cdr_to_phi_m0");
  if(m_verbosity > 5){
    pout() << "sisdc::copy_cdr_to_phi_m0" << endl;
  }

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
  CH_TIME("sisdc::predictor_compute_FD_0");
  if(m_verbosity > 5){
    pout() << "sisdc::predictor_compute_FD_0" << endl;
  }

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);
  
    const EBAMRCellData& phi_0 = storage->get_phi()[0]; // phi_0
    EBAMRCellData& FD_0        = storage->get_FD()[0];  // FD(phi_0)

    cdr_tga* tgasolver = (cdr_tga*) (&(*solver));
    tgasolver->compute_divD(FD_0, phi_0);
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
  sisdc::compute_E_into_scratch();
  sisdc::compute_cdr_eb_states();
  sisdc::compute_cdr_fluxes(a_time);
  sisdc::compute_cdr_domain_states();
  sisdc::compute_cdr_domain_fluxes(a_time);
  sisdc::compute_sigma_flux();

  // We begin with phi[0] = phi(t_n). Then update phi[m+1].
  for(int m = 0; m < m_p; m++){

    // This computes phi_(m+1) = phi_m + dtm*FAR_m(phi_m) + lagged quadrature and lagged advection-reaction
    sisdc::integrate_advection_reaction(a_dt, m, a_corrector);

    // This does the diffusion advance. It also adds in the remaining lagged diffusion terms before the implicit diffusion solve
    sisdc::integrate_diffusion(a_dt, m, a_corrector);

    // After the diffusion step we should update source terms and boundary conditions for the next step. We don't
    // do this on the last step. This is done either in the reconcile_integrands routine, or after SISDC is done
    // with its substebs. 
    const bool last = (m == m_p-1);
    if(!last){
      Vector<EBAMRCellData*> cdr_densities_mp1 = sisdc::get_cdr_phik(m+1);
      EBAMRIVData& sigma_mp1 = sisdc::get_sigmak(m+1);
      const Real t_mp1 = m_tm[m+1];

      // Update electric field, RTE equations, source terms, and velocities. No need to update diffusion since
      // that is done in the predictor_diffusion routine
      if(m_consistent_E)   sisdc::update_poisson(cdr_densities_mp1, sigma_mp1);
      if(m_consistent_rte) sisdc::update_rte(cdr_densities_mp1, t_mp1);
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
  if(!skip){
    if(m_subcycle){
      sisdc::integrate_advection_subcycle(a_dt, a_m, a_corrector);    
    }
    else{
      sisdc::integrate_advection_nosubcycle(a_dt, a_m, a_corrector);
    }
  }

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
      const EBAMRCellData& src = solver->get_source();    // Updated source

      // Increment swith source and then compute slope. 
      data_ops::incr(phi_m1, src, m_dtm[a_m]);  // phi_(m+1) = phi_m + dtm*(FA_m + FR_m)
      if(a_corrector){ // Back up the old slope first, we will need it for the lagged term
	data_ops::copy(scratch, FAR_m);
      }
      data_ops::copy(FAR_m, phi_m1);            // FAR_m = (phi_(m+1) - phi_m)/dtm
      data_ops::incr(FAR_m, phi_m, -1.0);       // :
      data_ops::scale(FAR_m, 1./m_dtm[a_m]);    // :
    }

    // Now add in the lagged advection-reaction and quadrature terms. This is a bit weird, but we did overwrite
    // FAR_m above after the advection-reaction advance, but we also backed up the old term into scratch. 
    if(a_corrector){
      data_ops::incr(phi_m1, scratch, -m_dtm[a_m]); // phi_(m+1)^(k+1) = phi_m^(k+1) + dtm*(FAR_m^(k+1) - FAR_m^k)
      sisdc::quad(scratch, storage->get_F(), a_m);  // Does the quadrature of the lagged operator slopes. 
      data_ops::incr(phi_m1, scratch, 0.5*a_dt);    // phi_(m+1)^(k+1) = phi_m^(k+1) + dtm*(FAR_m^(k+1) - FAR_m^k) + I_m^(m+1)
    }
  }

  // Add in the lagged terms for sigma. As above, m=0 and corrector is a special case where we just use the old slopes.
  EBAMRIVData& sigma_m1      = m_sigma_scratch->get_sigma()[a_m+1];
  const EBAMRIVData& sigma_m = m_sigma_scratch->get_sigma()[a_m];
  if(a_m == 0 && a_corrector){ 
    const EBAMRIVData& Fsig_m = m_sigma_scratch->get_Fold()[a_m]; // Here, we should be able to use either Fold or Fnew
    data_ops::copy(sigma_m1, sigma_m);                            // since Fsig_0 is only computed once. 
    data_ops::incr(sigma_m1, Fsig_m, m_dtm[a_m]);
  }

  if(a_corrector){ // Add in the lagged terms. When we make it here, sigma_(m+1) = sigma_m + dtm*Fsig_m. 
    EBAMRIVData& Fsig_lag = m_sigma_scratch->get_Fold()[a_m];
    data_ops::incr(sigma_m1, Fsig_lag, -m_dtm[a_m]);

    // Add in the quadrature term
    EBAMRIVData& scratch = m_sigma_scratch->get_scratch();
    sisdc::quad(scratch, m_sigma_scratch->get_Fold(), a_m);
    data_ops::incr(sigma_m1, scratch, 0.5*a_dt); // Mult by 0.5*a_dt due to scaling on [-1,1] for quadrature
  }      
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
    MayDay::Abort("sisdc::integrate_advection_nosubcycle - m=0 and corrector=true should never happen");
  }

  // Advance cdr equations
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);

    EBAMRCellData& phi_m1      = storage->get_phi()[a_m+1];
    EBAMRCellData& scratch     = storage->get_scratch();
    const EBAMRCellData& phi_m = storage->get_phi()[a_m];

    const Real extrap_dt = m_extrap_advect ? 2.0*m_extrap_dt*m_dtm[a_m] : 0.0; // Factor of 2 due to EBPatchAdvect
    solver->compute_divF(scratch, phi_m, extrap_dt, true);                     // scratch =  Div(v_m*phi_m^(k+1))

    data_ops::copy(phi_m1, phi_m);
    data_ops::incr(phi_m1, scratch, -m_dtm[a_m]);

    m_amr->average_down(phi_m1, m_cdr->get_phase());
    m_amr->interp_ghost(phi_m1, m_cdr->get_phase());
    data_ops::floor(phi_m1, 0.0);
  }

  // Update sigma. Also compute the new slope.
  EBAMRIVData& sigma_m1      = m_sigma_scratch->get_sigma()[a_m+1];
  EBAMRIVData& Fsig_new      = m_sigma_scratch->get_Fnew()[a_m];
  const EBAMRIVData& sigma_m = m_sigma_scratch->get_sigma()[a_m];
  m_sigma->compute_rhs(Fsig_new); // Fills Fsig_new with BCs from CDR solvers
  data_ops::copy(sigma_m1, sigma_m);
  data_ops::incr(sigma_m1, Fsig_new, m_dtm[a_m]);
}

void sisdc::integrate_advection_subcycle(const Real a_dt, const int a_m, const bool a_corrector){
  CH_TIME("sisdc::integrate_advection_nosubcycle");
  if(m_verbosity > 5){
    pout() << "sisdc::integrate_advection_nosubcycle" << endl;
  }

  MayDay::Abort("sisdc::integrate_advection_subcycle - not implemented");
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
      EBAMRCellData& phi_m1      = storage->get_phi()[a_m+1]; // 
      const EBAMRCellData& phi_m = storage->get_phi()[a_m];

      // Build the diffusion source term
      EBAMRCellData& source  = storage->get_scratch();
      data_ops::copy(source, phi_m1);      // Copy initial solutions
      data_ops::incr(source, phi_m, -1.0); // Copy initial solutions
      data_ops::scale(source, 1./m_dtm[a_m]);
      if(a_corrector){
	const EBAMRCellData& FD_m1k = storage->get_FD()[a_m+1];      // FD_(m+1)^k. Lagged term.
	data_ops::incr(source, FD_m1k, -1.0);
      }

      // Solve
      cdr_tga* tgasolver = (cdr_tga*) (&(*solver));
      tgasolver->advance_euler(phi_m1, phi_m, source, m_dtm[a_m]); // Source is -Fd_(m+1)^k
      m_amr->average_down(phi_m1, m_cdr->get_phase());
      m_amr->interp_ghost(phi_m1, m_cdr->get_phase());
      data_ops::floor(phi_m1, 0.0);


      // Update the operator slope
      EBAMRCellData& FD_m1k = storage->get_FD()[a_m+1];
      tgasolver->compute_divD(FD_m1k, phi_m1);
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

  // Update electric field, RTE equations, source terms, and velocities
  if(m_consistent_E)   sisdc::update_poisson(cdr_densities_p, sigma_p);
  if(m_consistent_rte) sisdc::update_rte(cdr_densities_p, t_p);
  if(m_compute_S)      sisdc::compute_cdr_gradients(cdr_densities_p);
  if(m_compute_v)      sisdc::compute_cdr_velo(cdr_densities_p, t_p);
  if(m_compute_S)      sisdc::compute_cdr_sources(cdr_densities_p, t_p);

  // Update boundary conditions for cdr and sigma equations
  sisdc::compute_cdr_eb_states(cdr_densities_p);
  sisdc::compute_cdr_fluxes(cdr_densities_p, t_p);
  sisdc::compute_cdr_domain_states(cdr_densities_p);
  sisdc::compute_cdr_domain_fluxes(cdr_densities_p, t_p);
  sisdc::compute_sigma_flux();

  // Now compute FAR_p - that wasn't done in the predictor or the corrector. 
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = sisdc::get_cdr_storage(solver_it);
    const int idx = solver_it.get_solver();

    // This has not been computed yet. Do it.
    EBAMRCellData& FAR_p       = storage->get_FAR()[m_p];
    const EBAMRCellData& phi_p = *cdr_densities_p[idx] ;
    const EBAMRCellData& src   = solver->get_source();

    const Real extrap_dt = m_extrap_advect ? 2.0*m_extrap_dt*m_dtm[m_p-1] : 0.0; // Factor of 2 because of EBPatchAdvect
    solver->compute_divF(FAR_p, phi_p, extrap_dt, true); // FAR_p =  Div(v_p*phi_p)
    data_ops::scale(FAR_p, -1.0);                        // FAR_p = -Div(v_p*phi_p)
    data_ops::incr(FAR_p, src, 1.0);                     // RHS = -Div(v_m*phi_m) + S_m = FAR(phi_m)

    // Build the integrand
    for (int m = 0; m <= m_p; m++){
      EBAMRCellData& F_m         = storage->get_F()[m];
      const EBAMRCellData& FD_m  = storage->get_FD()[m];
      const EBAMRCellData& FAR_m = storage->get_FAR()[m];

      data_ops::copy(F_m, FAR_m);
      data_ops::incr(F_m, FD_m, 1.0);

      m_amr->average_down(F_m, m_cdr->get_phase());
      m_amr->interp_ghost(F_m, m_cdr->get_phase());
    }
  }

  // Compute Fsig_p - that wasn't done in the predictor either
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

    // These should be zero
    EBAMRCellData& error = storage->get_error();
    const EBAMRCellData& phi_final = storage->get_phi()[m_p];

    data_ops::set_value(error, 0.0);
    data_ops::incr(error, phi_final, -1.0);
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
    EBAMRCellData& error       = storage->get_error();
    const EBAMRCellData& phi_p = storage->get_phi()[m_p];
    data_ops::incr(error, phi_p, 1.0);

    // Compute norms
    Real Lerr, Lphi;
    data_ops::norm(Lerr, *error[0], m_amr->get_domains()[0], m_error_norm);
    data_ops::norm(Lphi, *phi_p[0], m_amr->get_domains()[0], m_error_norm);
    m_cdr_error[idx] = (Lphi > 0.0) ? (Lerr/Lphi) : 0.0;

    m_max_error = Max(m_cdr_error[idx], m_max_error);
  }

  // Compute the surface charge conservation error
  EBAMRIVData& error = m_sigma_scratch->get_error();
  const EBAMRIVData& sigma_final = m_sigma_scratch->get_sigma()[m_p];
  data_ops::incr(error, sigma_final, 1.0);
  m_sigma_error = 0.0;

#if 1 // Debug
  if(procID() == 0) std::cout << m_cdr_error << std::endl;
#endif
}

void sisdc::compute_new_dt(bool& a_accept_step, const Real a_dt, const int a_num_corrections){
  CH_TIME("sisdc::compute_new_dt");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_new_dt" << endl;
  }

  // If a_dt was the smallest possible CFL or hardcap time step, we just have to accept it
  const Real max_gl_dist = sisdc::get_max_node_distance();
  const Real dt_cfl = 2.0*m_dt_cfl/max_gl_dist;
  if(a_dt <= dt_cfl*m_minCFL && m_max_error > m_err_thresh){ // No choice but to accept
    a_accept_step = true;
    m_new_dt = dt_cfl*m_minCFL;
    return;
  }
  if(a_dt <= m_min_dt && m_max_error > m_err_thresh){ // No choice but to accept
    a_accept_step = true;
    m_new_dt = m_min_dt;
    return;
  }

  // If we made it here, we should be able to decrease or increase the time step as we see fit.
  const Real rel_err  = (m_safety*m_err_thresh)/m_max_error;
  const Real dt_adapt = (m_max_error > 0.0) ? a_dt*pow(rel_err, 1.0/(a_num_corrections+1)) : m_max_dt;

  if(m_max_error <= m_err_thresh){ // Increase time step, but only if we're sufficiently far away form the error threshold
    if(m_max_error < m_safety*m_err_thresh){ 
      m_new_dt = dt_adapt;
    }
    else{
      //      m_new_dt = m_safety*a_dt;  // If we're too close, reduce the new time step with safety margin
      m_new_dt = m_safety*dt_adapt;
    }
    a_accept_step = true;
  }
  else{ // Decrease time step
    m_new_dt = m_safety*dt_adapt;
    a_accept_step = false;
  }

  // New time step can't go below minimum permitted CFL
  if(m_new_dt < dt_cfl*m_minCFL){ // Can't decrease above minimum CFL step, accept the step
    m_new_dt = Max(m_new_dt, dt_cfl*m_minCFL);
  }
  if(m_new_dt < m_min_dt){ // Can't decrease below hardcap, accept the step
    m_new_dt = Max(m_new_dt, m_min_dt);
  }
  if(m_new_dt > dt_cfl*m_maxCFL){
    m_new_dt = Min(m_new_dt, dt_cfl*m_maxCFL);
  }

  // Also, can't go beyond what is limited by dt_m...
  if(m_new_dt > dt_cfl){
    m_new_dt = Min(dt_cfl, m_new_dt);
  }

  { // Debug
#if 0 // Debug
    pout() << m_max_error << "\t"
	   << a_num_corrections << "\t"
	   << a_dt << "\t"
	   << m_new_dt << "\t"
	   << dt_adapt << "\t"
	   << endl;
#endif
  }

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

  const Real max_gl_dist = sisdc::get_max_node_distance();

  m_dt_cfl = m_cdr->compute_cfl_dt();
#if 0 // Debug
  if(procID() == 0) std::cout << m_dt_cfl << std::endl;
#endif
  if(!m_adaptive_dt){
    const Real dt_cfl = 2.0*m_cfl*m_dt_cfl/max_gl_dist;
    //const Real dt_cfl = m_cfl*m_dt_cfl;
    if(dt_cfl < dt){
      dt = dt_cfl;
      a_timecode = time_code::cfl;
    }
  }
  else{
    Real new_dt;
    const Real dt_cfl = 2.0*m_dt_cfl/max_gl_dist;
    if(m_have_dt_err){
      new_dt = m_new_dt;
      new_dt = Max(new_dt, dt_cfl*m_minCFL);
      new_dt = Min(new_dt, dt_cfl*m_maxCFL);
    }
    else{
      //      new_dt = 0.5*(m_maxCFL+m_minCFL)*dt_cfl;
      new_dt = m_maxCFL*dt_cfl;
    }

    if(new_dt < dt){
      dt = new_dt;
      a_timecode = time_code::error;
    }
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
    m_amr->average_down(grad, m_cdr->get_phase());
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

void sisdc::update_rte(const Real a_time){
  CH_TIME("sisdc::update_rte(solver)");
  if(m_verbosity > 5){
    pout() << "sisdc::update_rte(solver)" << endl;
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

void sisdc::update_rte(const Vector<EBAMRCellData*>& a_cdr_states, const Real a_time){
  CH_TIME("sisdc::update_rte(full)");
  if(m_verbosity > 5){
    pout() << "sisdc::update_rte(full)" << endl;
  }
  
  if(m_do_rte){
    if((m_step + 1) % m_fast_rte == 0){
      Vector<EBAMRCellData*> rte_states  = m_rte->get_states();
      Vector<EBAMRCellData*> rte_sources = m_rte->get_sources();

      EBAMRCellData& E = m_poisson_scratch->get_E_cell();

      const Real dummy_dt = 0.0;
      this->solve_rte(rte_states, rte_sources, a_cdr_states, E, a_time, dummy_dt, centering::cell_center);
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
