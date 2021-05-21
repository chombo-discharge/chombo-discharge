/*!
  @file   sdc.cpp
  @brief  Implementation of sdc.H
  @author Robert Marskar
  @date   Feb. 2019
*/

#include "sdc.H"
#include "sdc_storage.H"
#include <CD_DataOps.H>
#include <CD_Units.H>
#include <CD_CdrGodunov.H>
#include "CD_FieldSolverMultigrid.H"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <ParmParse.H>

typedef sdc::cdr_storage     cdr_storage;
typedef sdc::poisson_storage poisson_storage;
typedef sdc::rte_storage     rte_storage;
typedef sdc::sigma_storage   sigma_storage;

sdc::sdc(){
  m_className = "sdc";
  m_subcycle = false;
}

sdc::~sdc(){

}

void sdc::parseOptions(){
  CH_TIME("sdc::parseOptions");
  if(m_verbosity > 5){
    pout() << "sdc::parseOptions" << endl;
  }

  // Regular stuff from TimeStepper that we almost always need
  parseVerbosity();
  parse_solver_verbosity();
  parse_cfl();
  parse_relax_time();
  parse_fast_rte();
  parse_fast_poisson();
  parse_min_dt();
  parse_max_dt();
  parse_semi_implicit();

  // Specific to this class
  parse_nodes();
  parse_adaptive_options();
  parse_debug_options();
  parse_advection_options();
  parse_photoi_options();

  // Setup nodes
  sdc::setup_quadrature_nodes(m_p);
  sdc::setup_qmj(m_p);
}

void sdc::parse_nodes(){
  CH_TIME("sdc::parse_nodes");
  if(m_verbosity > 5){
    pout() << "sdc::parse_nodes" << endl;
  }

  ParmParse pp(m_className.c_str());
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
    MayDay::Abort("sdc::parse_nodes - unknown node type requested");
  }

  pp.get("subintervals",     m_p);
  pp.get("corr_iter",        m_k);

  if(m_p < 1){
    MayDay::Abort("sdc::parse_nodes - sdc.subintervals cannot be < 1");
  }
  if(m_k < 0){
    MayDay::Abort("sdc::parse_nodes - sdc.corr_iter cannot be < 0");
  }
}

void sdc::parse_adaptive_options(){
  CH_TIME("sdc::parse_adaptive_options");
  if(m_verbosity > 5){
    pout() << "sdc::parse_adaptive_options" << endl;
  }

  ParmParse pp(m_className.c_str());
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

void sdc::parse_debug_options(){
  CH_TIME("sdc::parse_debug_options");
  if(m_verbosity > 5){
    pout() << "sdc::parse_debug_options" << endl;
  }

  ParmParse pp(m_className.c_str());
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

void sdc::parse_advection_options(){
  CH_TIME("sdc::parse_advection_options");
  if(m_verbosity > 5){
    pout() << "sdc::parse_advection_options" << endl;
  }

  ParmParse pp(m_className.c_str());
  std::string str;

  m_extrap_dt = 0.5;  // Relic of an ancient past. I don't see any reason why extrapolating to anything but the half interval
                      // would make sense. 

  pp.get("extrap_advect", str); m_extrap_advect = (str == "true") ? true : false;
}

void sdc::parse_semi_implicit(){
  CH_TIME("sdc::parse_semi_implicit");
  if(m_verbosity > 5){
    pout() << "sdc::parse_semi_implicit" << endl;
  }

  ParmParse pp(m_className.c_str());
  std::string str;

  pp.get("semi_implicit", str); m_semi_implicit = (str == "true") ? true : false;
}

void sdc::parse_photoi_options(){
  CH_TIME("sdc::parse_photoi_options");
  if(m_verbosity > 5){
    pout() << "sdc::parse_photoi_options" << endl;
  }

  ParmParse pp(m_className.c_str());
  std::string str;

  pp.get("freeze_photoi", str); m_freeze_photoi = (str == "true") ? true : false;
}

RefCountedPtr<cdr_storage>& sdc::get_cdr_storage(const CdrIterator& a_solverit){
  return m_cdr_scratch[a_solverit.get_solver()];
}

RefCountedPtr<rte_storage>& sdc::get_rte_storage(const RtIterator& a_solverit){
  return m_rte_scratch[a_solverit.get_solver()];
}

bool sdc::needToRegrid(){
  CH_TIME("sdc::needToRegrid");
  if(m_verbosity > 5){
    pout() << "sdc::needToRegrid" << endl;
  }

  return false;
}

Real sdc::restrict_dt(){
  return 1.E99;
}

Real sdc::get_max_node_distance(){
  CH_TIME("sdc::get_max_node_distance");
  if(m_verbosity > 5){
    pout() << "sdc::get_max_node_distance" << endl;
  }

  Real max_dist = 0.0;
  for (int m = 0; m < m_p; m++){
    max_dist = Max(max_dist, m_nodes[m+1] - m_nodes[m]);
  }

  return max_dist;
}

void sdc::init(){
  CH_TIME("sdc::init");
  if(m_verbosity > 5){
    pout() << "sdc::init" << endl;
  }

  advance_reaction_network(m_time, m_dt);
}

void sdc::setup_quadrature_nodes(const int a_p){
  CH_TIME("sdc::setup_quadrature_nodes");
  if(m_verbosity > 5){
    pout() << "sdc::setup_quadrature_nodes" << endl;
  }

  if(m_which_nodes == "uniform"){
    sdc::setup_uniform_nodes(a_p);
  }
  else if(m_which_nodes == "lobatto"){
    sdc::setup_lobatto_nodes(a_p);
  }
  else if(m_which_nodes == "chebyshev"){
    sdc::setup_chebyshev_nodes(a_p);
  }
  else {
    MayDay::Abort("sdc::setup_quadrature_nodes - unknown nodes requested");
  }
}

void sdc::setup_uniform_nodes(const int a_p){
  CH_TIME("sdc::setup_uniform_nodes");
  if(m_verbosity > 5){
    pout() << "sdc::setup_uniform_nodes" << endl;
  }

  // TLDR: The nodes and weights are hardcoded. A better programmer would compute these
  //       recursively with Legendre polynomials. 
  m_nodes.resize(1+a_p);

  const Real delta = 2./a_p;
  for (int m = 0; m <= a_p; m++){
    m_nodes[m] = m*delta;
  }
}

void sdc::setup_lobatto_nodes(const int a_p){
  CH_TIME("sdc::setup_lobatto_nodes");
  if(m_verbosity > 5){
    pout() << "sdc::setup_lobatto_nodes" << endl;
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
    MayDay::Abort("sdc::setup_lobatto_nodes - requested order exceeds 7. Compute your own damn nodes!");
  }
}

void sdc::setup_chebyshev_nodes(const int a_p){
  CH_TIME("sdc::setup_chebyshev_nodes");
  if(m_verbosity > 5){
    pout() << "sdc::setup_chebyshev_nodes" << endl;
  }

  // TLDR: The nodes and weights are hardcoded. A better programmer would compute these
  //       recursively with Legendre polynomials. 
  m_nodes.resize(1+a_p);
  m_nodes[0] = -1.0;
  for (int m = 1; m < a_p; m++){
    m_nodes[m] = -cos((2*m-1)*Units::pi/(2*(a_p-1)));
  }
  m_nodes[a_p] = 1.0;
}

void sdc::setup_qmj(const int a_p){
  CH_TIME("sdc::setup_qmj");
  if(m_verbosity > 5){
    pout() << "sdc::setup_qmj" << endl;
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
    if(INFO != 0) MayDay::Abort("sdc::setup_qmj - could not compute weights");
    
    // Now construct qmj
    for (int m = 0; m < a_p; m++){
      m_qmj[m][j] = 0.0;
      for (int k = 0; k < nnodes; k++){
	m_qmj[m][j] += cj[k]*(pow(m_nodes[m+1], k+1) - pow(m_nodes[m], k+1))/(k+1);
      }
    }
  }
}

void sdc::setup_subintervals(const Real a_time, const Real a_dt){
  CH_TIME("sdc::setup_subintervals");
  if(m_verbosity > 5){
    pout() << "sdc::setup_subintervals" << endl;
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

  // dtm = t_{m+1} - t_m. Order 1 is special since we only use the SDC predictor from a second order formulation
  m_dtm.resize(m_tm.size() - 1);
  for (int m = 0; m < m_tm.size()-1; m++){
    m_dtm[m] = m_tm[m+1] - m_tm[m];
  }
}

void sdc::quad(EBAMRCellData& a_quad, const Vector<EBAMRCellData>& a_integrand, const int a_m){
  CH_TIME("sdc::quad");
  if(m_verbosity > 5){
    pout() << "sdc::quad" << endl;
  }

  if(a_m < 0)     MayDay::Abort("sdc::quad - bad index a_m < 0");
  if(a_m >= m_p)  MayDay::Abort("sdc::quad - bad index a_m >= m_p");

  DataOps::setValue(a_quad, 0.0);
  for (int j = 0; j <= m_p; j++){
    DataOps::incr(a_quad, a_integrand[j], m_qmj[a_m][j]);
  }
}

void sdc::quad(EBAMRIVData& a_quad, const Vector<EBAMRIVData>& a_integrand, const int a_m){
  CH_TIME("sdc::quad");
  if(m_verbosity > 5){
    pout() << "sdc::quad" << endl;
  }

  if(a_m < 0)     MayDay::Abort("sdc::quad - bad index a_m < 0");
  if(a_m >= m_p)  MayDay::Abort("sdc::quad - bad index a_m >= m_p");

  DataOps::setValue(a_quad, 0.0);
  for (int j = 0; j <= m_p; j++){
    DataOps::incr(a_quad, a_integrand[j], m_qmj[a_m][j]);
  }
}
  
void sdc::copy_phi_p_to_cdr(){
  CH_TIME("sdc::copy_phi_p_to_cdr");
  if(m_verbosity > 5){
    pout() << "sdc::copy_phi_p_to_cdr" << endl;
  }

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>&  solver  = solver_it();
    RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);

    EBAMRCellData& phi = solver->getPhi();
    const EBAMRCellData& phip = storage->get_phi()[m_p];
    DataOps::copy(phi, phip);
  }
}

void sdc::copy_phi_p_to_rte(){
  CH_TIME("sdc::copy_phi_p_to_rte");
  if(m_verbosity > 5){
    pout() << "sdc::copy_phi_p_to_rte" << endl;
  }

  for (RtIterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>&  solver  = solver_it();
    RefCountedPtr<rte_storage>& storage = get_rte_storage(solver_it);

    EBAMRCellData& phi = solver->getPhi();
    const EBAMRCellData& phip = storage->get_phi()[m_p];
    DataOps::copy(phi, phip);
  }
}

void sdc::copy_sigma_p_to_sigma(){
  CH_TIME("sdc::copy_sigma_p_to_sigma");
  if(m_verbosity > 5){
    pout() << "sdc::copy_sigma_p_to_sigma" << endl;
  }

  EBAMRIVData& sigma        = m_sigma->getPhi();
  const EBAMRIVData& sigmap = m_sigma_scratch->get_sigma()[m_p];
  DataOps::copy(sigma, sigmap);
}

void sdc::copy_FRp_to_FR0(){
  CH_TIME("sdc::copy_FRp_to_FR0");
  if(m_verbosity > 5){
    pout() << "sdc::copy_FRp_to_FR0" << endl;
  }

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);

    EBAMRCellData& FR0       = storage->get_FR()[0];
    const EBAMRCellData& FRp = storage->get_FR()[m_p];
    DataOps::copy(FR0, FRp);
  }
}

Real sdc::advance(const Real a_dt){
  CH_TIME("sdc::advance");
  if(m_verbosity > 2){
    pout() << "sdc::advance" << endl;
  }

  // ---------------------------------------------------------------------------------------------------
  // TLDR:  When we enter this routine, solvers SHOULD have been filled with valid ready and be ready 
  //        advancement. If you think that this may not be the case, activate the debugging below
  // ---------------------------------------------------------------------------------------------------

  // Backup solvers
  sdc::store_solvers();
  
  // Initialize integrations. If we do corrections
  sdc::copy_cdr_to_phi_m0();
  sdc::copy_sigma_to_sigma_m0();
  sdc::copy_rte_to_phi_m0();

  // SDC advance
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
    sdc::setup_subintervals(m_time, actual_dt);

    // First SDC sweep. No lagged slopes here.
    sdc::sweep(a_dt, m_time, false);

    // SDC correction sweeps. Need to take care of lagged terms. 
    for(int icorr = 0; icorr < Max(m_k, m_min_corr); icorr++){
      num_corrections++;

      // Initialize error and reconcile integrands (i.e. make them quadrature-ready)
      sdc::initialize_errors();
      sdc::reconcile_integrands();

      // SDC correction along whole interval
      sdc::sweep(a_dt, m_time, true);

      // Compute error and check if we need to keep iterating
      sdc::finalize_errors();
      if(m_max_error < m_err_thresh && m_adaptive_dt && icorr >= m_min_corr) break; // No need in going beyond
    }

    // Compute a new time step. If it is smaller than the minimum allowed CFL step, accept the step anyways
    if(m_adaptive_dt){
      sdc::compute_new_dt(accept_step, actual_dt, num_corrections);
      
      if(!accept_step){  // Step rejection, use the new dt for next step.
	actual_dt = m_new_dt;
	num_reject++;

	retry_step  = num_reject <= m_max_retries;
	
	if(retry_step){
	  sdc::restore_solvers();
	  sdc::compute_E_into_scratch();
	  sdc::compute_cdr_gradients();
	  sdc::compute_cdr_velo(m_time);
	  TimeStepper::compute_cdr_diffusion(m_fieldSolver_scratch->get_E_cell(), m_fieldSolver_scratch->get_E_eb());
	}
      }
    }
    else{
      m_new_dt = 1.234567E89;
      accept_step = true;
    }
  }

  // Copy results back to solvers. Also add the reactive slope for the next iteration. This is needed because I_m needs slopes
  // at all nodes for the high-order quadrature so we take the reactive slope at m=p at the end of this step and put it
  // on m=0 for the next step. If this step is the last step before a regrid, this class will regrid the FR_0 slope and
  // put it back where it belongs. 
  sdc::copy_phi_p_to_rte();
  sdc::copy_phi_p_to_cdr();
  sdc::copy_sigma_p_to_sigma();
  sdc::copy_FRp_to_FR0();       

  // Always recompute velocities and diffusion coefficients before the next time step. The Poisson and RTE equations
  // have been updated when we come in here. 
  sdc::compute_cdr_velo(m_time + actual_dt);
  TimeStepper::compute_cdr_diffusion(m_fieldSolver_scratch->get_E_cell(), m_fieldSolver_scratch->get_E_eb());


  // Profile step
  if(m_print_report)  sdc::adaptive_report(first_dt, actual_dt, m_new_dt, num_corrections, num_reject, m_max_error);
  if(m_profile_steps) sdc::write_step_profile(actual_dt, m_max_error, m_p, num_corrections, num_reject);

  // Store current error. 
  m_have_err  = true;
  m_pre_error = m_max_error;
  
  return actual_dt;
}

void sdc::copy_cdr_to_phi_m0(){
  CH_TIME("sdc::copy_cdr_to_phi_m0");
  if(m_verbosity > 5){
    pout() << "sdc::copy_cdr_to_phi_m0" << endl;
  }

  // CDR solvers
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>&  solver  = solver_it();
    RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);
    
    EBAMRCellData& phi0 = storage->get_phi()[0];
    const EBAMRCellData& phi = solver->getPhi();
    DataOps::copy(phi0, phi);
  }
}

void sdc::copy_rte_to_phi_m0(){
  CH_TIME("sdc::copy_rte_to_phi_m0");
  if(m_verbosity > 5){
    pout() << "sdc::copy_rte_to_phi_m0" << endl;
  }

  // RTE solvers
  for (RtIterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>&  solver  = solver_it();
    RefCountedPtr<rte_storage>& storage = get_rte_storage(solver_it);
    
    EBAMRCellData& phi0 = storage->get_phi()[0];
    const EBAMRCellData& phi = solver->getPhi();
    DataOps::copy(phi0, phi);
  }
}

void sdc::copy_sigma_to_sigma_m0(){
  CH_TIME("sdc::copy_sigma_sigma_m0");
  if(m_verbosity > 5){
    pout() << "sdc::copy_sigma_to_sigma_m0" << endl;
  }

  // Copy sigma to starting state
  EBAMRIVData& sigma0      = m_sigma_scratch->get_sigma()[0];
  const EBAMRIVData& sigma = m_sigma->getPhi();
  DataOps::copy(sigma0, sigma);
}

void sdc::sweep(const Real a_dt, const Real a_time, const bool a_lagged_terms){
  CH_TIME("sdc::sweep");
  if(m_verbosity > 5){
    pout() << "sdc::sweep" << endl;
  }

  if(m_semi_implicit){
    sdc::sweep_semi_implicit(a_dt, a_time, a_lagged_terms);
  }
  else{
    sdc::sweep_explicit(a_dt, a_time, a_lagged_terms);
  }
}

void sdc::sweep_semi_implicit(const Real a_dt, const Real a_time, const bool a_corrector){
  CH_TIME("sdc::sweep_semi_implicit");
  if(m_verbosity > 5){
    pout() << "sdc::sweep_semi_implicit" << endl;
  }

  // For the semi-implicit coupling we are solving
  //
  //    phi_(m+1)^(k+1) = phi_m^(k+1) + dt*[ div(D*grad(phi_m^(k+1))) - div(mu*E^(k+1)_(m+1)*phi_m^(k+1))]
  //                                  - dt*[ div(D*grad(phi_m^k)) - div(mu*E^k_(m+1)*phi_m^k)]
  //                                  + I_m^(m+1)(phi^k)
  //
  // We then obtain E^(k+1)_(m+1) through the Poisson equation by inserting phi_(m+1)^(k+1) as rho^(k+1)_(m+1). Charge
  // injection onto dielectrics is handled explicitly so that the CDR BCs are updated BEFORE computing E. Here are all the
  // stages for one SDC substep:
  //
  // 1.  Update boundary conditions
  //
  // 2.  Compute the explicit diffusion operator div(D*grad(phi_m^(k+1))) = DivD
  //
  // 3.  Make the right-hand side for the Poisson equation. This is
  //
  //        rho^(k+1)_(m+1) = Sum[q*(phi_m^(k+1) + dtm*divD - dtm*(FA(phi^k) + FD(phi^k) + I(phi^k))
  //
  // 4.  Compute mobilities as mu = |v|/|E|. Beware division by zero. 
  //
  // 5.  Adjust the permittivity, i.e. re-initilizae multigrid, for the Poisson equation as
  //
  //        bco = eps + dtm/eps0*Sum[q*mu_phi*phi_m^(k+1)]
  //
  // 6.  Advance sigma equation to sigma_(m+1)^(k+1)
  // 
  // 7.  Solve for E_(m+1)^(k+1)
  //
  // 8.  Advance the reaction network to get the source
  //
  // 9.  Compute new cdr velocities
  //
  // 10. Compute DivF
  //
  // 11. Solve the equation above for phi_(m+1)^(k+1)
  //
  // 12. Solve the radiative transfer equations

  Real time = a_time;

  // This loop solves for phi_(m+1)^(k+1). 
  for (int m = 0; m < m_p; m++){

    // We must have this for doing mu = |v|/|E|
    sdc::compute_E_into_scratch();
      
    // 1. Update boundary conditions and charge fluxes
    sdc::compute_cdr_gradients();
    sdc::compute_cdr_eb_states(); 
    sdc::compute_cdr_fluxes(a_time);
    sdc::compute_cdr_domain_states();
    sdc::compute_cdr_domain_fluxes(a_time);
    sdc::compute_sigma_flux();

    // 2. Compute the explicit divergence. The operator div(D*grad(phi_m^(k+1))) goes onto m_divD storage
    sdc::computeDivD(m, a_corrector);

    // 3. Make semi-implicit rho
    sdc::compute_semi_implicit_rho(m, a_corrector);

    // 4. Compute mobilities as mu = |v|/|E|
    sdc::compute_semi_implicit_mobilities(m, a_corrector);

    // 5. Adjust permittivities and setup poisson solver again
    sdc::set_semi_implicit_permittivities();

    // 6. Advance sigma equation
    sdc::advance_sigma(m, a_corrector);

    // 7. Solve the Poisson equation and compute the electric field E_(m+1)^(k+1)
    sdc::solve_semi_implicit_poisson(m);

    // 8. Advance the reaction network to get source terms
    sdc::compute_reaction_network(m, time, m_dtm[m], a_corrector);

    // 9. Compute the new cdr velocities
    sdc::compute_semi_implicit_cdr_velocities(m, time);

    // 9b. Compute new boundary conditions using the new field and velocities. Really need to rework this one. 
#if 1
    sdc::compute_cdr_gradients();
    sdc::compute_cdr_eb_states();
    sdc::compute_cdr_fluxes(a_time);
    sdc::compute_cdr_domain_states();
    sdc::compute_cdr_domain_fluxes(a_time);
    sdc::compute_sigma_flux();
#endif

    // 10. Compute div(v*phi)
    sdc::computeDivF(m, a_corrector);

    // 11. Solve the stinking equation
    sdc::substep_cdr(m, a_corrector);

    // 12. Solve the RTE equation
    sdc::substep_rte(m, a_corrector);

    // Increment time
    time = time + m_dtm[m];
  }
}

void sdc::computeDivD(const int a_m, const bool a_corrector){
  CH_TIME("sdc::computeDivD");
  if(m_verbosity > 5){
    pout() << "sdc::computeDivD" << endl;
  }

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);

    EBAMRCellData& divD = storage->get_divD();
    const EBAMRCellData& phi = storage->get_phi()[a_m];
    const EBAMRCellData& FD0 = storage->get_FD()[0];
    
    // If a_m==0 we should only have to compute divD once, and this happens during the initial SDC sweep. For the corrector
    // we should always be able to use the m_FD[0] slope so we simply copy that over
    if(solver->isDiffusive()){
      if(a_corrector && a_m > 0){
	solver->computeDivD(divD, phi);
      }
      else if(a_m == 0 && !a_corrector){
	solver->computeDivD(divD, phi);
      }
      else{
	DataOps::copy(divD, FD0);
      }
    }
  }
}

void sdc::computeDivF(const int a_m, const bool a_corrector){
  CH_TIME("sdc::computeDivF");
  if(m_verbosity > 5){
    pout() << "sdc::computeDivF" << endl;
  }

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);

    EBAMRCellData& divF      = storage->get_divF();
    const EBAMRCellData& phi = storage->get_phi()[a_m];
    const EBAMRCellData& FA0 = storage->get_FA()[0];

    const Real extrap_dt = (m_extrap_advect) ? m_dtm[a_m] : 0.0;
    
    // If a_m==0 we should only have to compute divD once, and this happens during the initial SDC sweep. For the corrector
    // we should always be able to use the m_FD[0] slope so we simply copy that over
    if(solver->isMobile()){
      if(a_corrector && a_m > 0){
	solver->computeDivF(divF, phi, extrap_dt);
      }
      else if(a_m == 0 && !a_corrector){
	solver->computeDivF(divF, phi, extrap_dt);
      }
      else{
	DataOps::copy(divF, FA0);
      }
    }
    else{
      DataOps::setValue(divF, 0.0);
    }
  }
}

void sdc::compute_semi_implicit_mobilities(const int a_m, const bool a_corrector){
  CH_TIME("sdc::compute_semi_implicit_mobilities");
  if(m_verbosity > 5){
    pout() << "sdc::compute_semi_implicit_mobilities" << endl;
  }

  // TLDR: This routine computes the mobility = |v|/|E| and then modifies it as
  //
  //    mobility = |v|/|E|*q*dt/eps0*phi

  // Compute |E| first
  const EBAMRCellData& E = m_fieldSolver_scratch->get_E_cell();
  DataOps::vectorLength(m_scratch1, E);

  // This is the dt for m->m+1
  const Real dtm = m_dtm[a_m];
  
  // This iterates over solvers. If they are mobile and have charge != 0 we extact a mobility
  // by computing mobility = |v|/|E|
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver   = solver_it();
    const RefCountedPtr<cdr_storage>& storage = sdc::get_cdr_storage(solver_it);

    const int q      = (solver_it.getSpecies())->getChargeNumber();
    const int q_sign = q > 0 ? 1 : -1;

    if(q != 0 && solver->isMobile()){
      EBAMRCellData& cell_mob   = storage->get_cell_mob();
      EBAMRFluxData& face_mob   = storage->get_face_mob();
      EBAMRIVData&   eb_mob     = storage->get_eb_mob();
      const EBAMRCellData& velo = solver->getCellCenteredVelocity();
      const EBAMRCellData& phi  = storage->get_phi()[a_m];

      DataOps::vectorLength(cell_mob, velo);
      DataOps::divideByScalar(cell_mob, m_scratch1); // This gives the magnitude of the cell-centered mobility
      DataOps::scale(cell_mob, q_sign);             // Set the correct sign, i.e. v = sgn(q)*mobility*E

      // Now make the cell-centered mobility equal to
      DataOps::multiply(cell_mob, phi);
      DataOps::scale(cell_mob, q*dtm*Units::Qe/Units::eps0);

      // Get valid ghost cells before averaging
      m_amr->averageDown(cell_mob, phase::gas);
      m_amr->interpGhost(cell_mob, phase::gas);

      // Compute face-centered "mobility"
      DataOps::averageCellToFaceAllComps(face_mob, cell_mob, m_amr->getDomains());

      // Compute EB-centered mobility
      TimeStepper::extrapolate_to_eb(eb_mob, phase::gas, cell_mob);
    }
  }
  
}

void sdc::compute_semi_implicit_rho(const int a_m,  const bool a_corrector){
  CH_TIME("sdc::compute_semi_implicit_rho");
  if(m_verbosity > 5){
    pout() << "sdc::compute_semi_implicit_rho" << endl;
  }

  // TLDR: This routine computes
  //
  //    rho = sum q*[phi_m^k + dtm*div(D*grad(phi_m^k)) - dtm*FD(phi_m^k) - dtm*FA(phi_m^k) + I_m(phi^k)]
  // 
  
  // Get gas-side rho
  MFAMRCellData& src = m_fieldSolver->getRho();
  DataOps::setValue(src, 0.0);

  // Get handle to gas-side rho
  EBAMRCellData rho_gas;
  m_amr->allocatePointer(rho_gas); 
  m_amr->alias(rho_gas, phase::gas, src);

  // This is the time step from m->m+1
  const Real dtm = m_dtm[a_m];

  // This iterates over solvers and add the terms if q > 0
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver   = solver_it();
    const RefCountedPtr<cdr_storage>& storage = sdc::get_cdr_storage(solver_it);

    const int q = (solver_it.getSpecies())->getChargeNumber();
    
    if(q != 0){

      // This adds the term q*phi_m^k
      const EBAMRCellData& phi = storage->get_phi()[a_m];
      DataOps::incr(rho_gas, phi, 1.0*q);

      // Add diffusive term. This is not a lagged term
      if(solver->isDiffusive()){
	const EBAMRCellData& divD = storage->get_divD();
	DataOps::incr(rho_gas, divD, dtm*q);
      }

      // This adds the other terms. The quadrature required for I_m(phi^k) is performed in place. 
      if(a_corrector){
	const EBAMRCellData& FA   = storage->get_FA()[a_m];
	const EBAMRCellData& FD   = storage->get_FD()[a_m];

	// Compute the quadrature
	EBAMRCellData& scratch = storage->get_scratch();
	sdc::quad(scratch, storage->get_F(), a_m);

	// Increment
	if(solver->isMobile()){
	  DataOps::incr(rho_gas, FA,     -dtm*q);
	}
	if(solver->isDiffusive()){
	  DataOps::incr(rho_gas, FD,     -dtm*q);
	}
	DataOps::incr(rho_gas, scratch,  0.5*q*dtm); // Scaled by 0.5*dt since quadrature takes place on [-1,1]
      }
    }
  }

  // Now do the scaling
  DataOps::scale(rho_gas, Units::Qe);
  m_amr->interpToCentroids(rho_gas, phase::gas);
}

void sdc::set_semi_implicit_permittivities(){
  CH_TIME("sdc::set_semi_implicit_permittivities");
  if(m_verbosity > 5){
    pout() << "sdc::set_semi_implicit_permittivities" << endl;
  }

  FieldSolverMultigrid* poisson = (FieldSolverMultigrid*) (&(*m_fieldSolver));

  // Set coefficients as usual
  poisson->setMultigridCoefficients();

  // Get bco and increment with mobilities
  MFAMRFluxData& bco   = poisson->getAvgBco();
  MFAMRIVData& bco_irr = poisson->getAvgBco_irreg();
  
  EBAMRFluxData bco_gas;
  EBAMRIVData   bco_irr_gas;
  
  m_amr->allocatePointer(bco_gas);
  m_amr->allocatePointer(bco_irr_gas);
  
  m_amr->alias(bco_gas,     phase::gas, bco);
  m_amr->alias(bco_irr_gas, phase::gas, bco_irr);
  
  // Increment that shit. 
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver   = solver_it();
    const RefCountedPtr<cdr_storage>& storage = sdc::get_cdr_storage(solver_it);

    const int q = (solver_it.getSpecies())->getChargeNumber();
    
    if(q != 0 && solver->isMobile()){
      const EBAMRFluxData& face_mob = storage->get_face_mob();
      const EBAMRIVData& eb_mob     = storage->get_eb_mob();

      DataOps::incr(bco_gas,     face_mob, 1.0);
      DataOps::incr(bco_irr_gas, eb_mob,   1.0);
    }
  }

  // Set up the multigrid solver
  poisson->setupOperatorFactory();
  poisson->setup_solver();
  poisson->set_needs_setup(false);
}

void sdc::advance_sigma(const int a_m, const bool a_corrector){
  CH_TIME("sdc::advance_sigma");
  if(m_verbosity > 5){
    pout() << "sdc::advance_sigma" << endl;
  }

  // Add in the lagged terms for sigma. As above, m=0 and corrector is a special case where we just use the old slopes.
  EBAMRIVData& sigma_m1      = m_sigma_scratch->get_sigma()[a_m+1];
  const EBAMRIVData& sigma_m = m_sigma_scratch->get_sigma()[a_m];
  const EBAMRIVData& Fsig_m = m_sigma_scratch->get_Fold()[a_m]; // Here, we should be able to use either Fold or Fnew

  const Real dtm = m_dtm[a_m];
  
  // 
  DataOps::copy(sigma_m1, sigma_m); 
  DataOps::incr(sigma_m1, Fsig_m, dtm);

  // Corrector needs lagged terms
  if(a_corrector){
    EBAMRIVData& Fsig_lag = m_sigma_scratch->get_Fold()[a_m];
    DataOps::incr(sigma_m1, Fsig_lag, -dtm);

    // Add in the quadrature term
    EBAMRIVData& scratch = m_sigma_scratch->get_scratch();
    sdc::quad(scratch, m_sigma_scratch->get_Fold(), a_m);
    DataOps::incr(sigma_m1, scratch, 0.5*dtm); // Mult by 0.5*a_dt due to scaling on [-1,1] for quadrature
  }
}

void sdc::solve_semi_implicit_poisson(const int a_m){
  CH_TIME("sdc::solve_semi_implicit_poisson");
  if(m_verbosity > 5){
    pout() << "sdc::solve_semi_implicit_poisson" << endl;
  }

  MFAMRCellData& phi       = m_fieldSolver->getPotential();
  const MFAMRCellData& rho = m_fieldSolver->getRho();
  const EBAMRIVData& sigma = m_sigma_scratch->get_sigma()[a_m+1];

  m_fieldSolver->solve(phi, rho, sigma);

  sdc::compute_E_into_scratch();
}

void sdc::compute_semi_implicit_cdr_velocities(const int a_m, const Real a_time){
  CH_TIME("sdc::compute_semi_implicit_cdr_velocities");
  if(m_verbosity > 5){
    pout() << "sdc::compute_semi_implicit_cdr_velocities" << endl;
  }

  Vector<EBAMRCellData*> states = get_cdr_phik(a_m);
  sdc::compute_cdr_velo(states, a_time);
}

void sdc::substep_cdr(const int a_m, const bool a_corrector){
  CH_TIME("sdc::substep_cdr");
  if(m_verbosity > 5){
    pout() << "sdc::substep_cdr" << endl;
  }

  const Real dtm = m_dtm[a_m];

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);

    EBAMRCellData& phi_m1      = storage->get_phi()[a_m+1];
    EBAMRCellData& I_m         = storage->get_scratch();  
    const EBAMRCellData& phi_m = storage->get_phi()[a_m];
    const EBAMRCellData& divF  = storage->get_divF();     // Beware the sign, FA = -DivF
    const EBAMRCellData& divD  = storage->get_divD();  
    const EBAMRCellData& R     = solver->getSource();

    EBAMRCellData& FA_lag    = storage->get_FA()[a_m];
    EBAMRCellData& FD_lag    = storage->get_FD()[a_m];
    EBAMRCellData& FR_lag    = storage->get_FR()[a_m+1]; // Beware the centering monster for implicit reactions
    const Vector<EBAMRCellData>& F = storage->get_F();

    // Do the lagged term high-order quadrature
    if(a_corrector){
      sdc::quad(I_m, F, a_m);
      DataOps::scale(I_m, 0.5*dtm);
    }

    // phi_(m+1) = phi_m - dt*divF - dt*FA^k
    DataOps::copy(phi_m1, phi_m);
    if(solver->isMobile()){
      DataOps::incr(phi_m1, divF, -dtm);
      
      if(a_corrector){ // Do lagged term if this is the corrector
	DataOps::incr(phi_m1, FA_lag, -dtm);
      }
    }

    // phi_(m+1) = phi_m - dt*divF - dt*FA^k + dt*divD - dt*FD^k
    if(solver->isDiffusive()){
      DataOps::incr(phi_m1, divD, dtm);
      
      if(a_corrector){ // Do lagged term if this is the corrector
	DataOps::incr(phi_m1, FD_lag, -dtm);
      }
    }

    // phi_(m+1) = phi_m - dt*divF - dt*FA^k + dt*divD -dt*FD^k + dt*FR - dt*FR^k
    DataOps::incr(phi_m1, R, dtm);
    if(a_corrector){
      DataOps::incr(phi_m1, FR_lag, -dtm);
    }

    // phi_(m+1) = phi_m - dt*divF - dt*FA^k + dt*divD -dt*FD^k + dt*FR - dt*FR^k + I_m(phi^k)
    if(a_corrector){
      DataOps::incr(phi_m1, I_m, 1.0);
    }

    // Non-negative magic stencil
    //    solver->make_non_negative(phi_m1);
    DataOps::floor(phi_m1, 0.0);

    // Overwrite the slopes. Could probably save some m=0 computations here since divF and divD doesn't change there
    DataOps::copy(FA_lag,   divF); 
    DataOps::scale(FA_lag, -1.0); 
    DataOps::copy(FD_lag,   divD);
    DataOps::copy(FR_lag,   R);
  }
}

void sdc::substep_rte(const int a_m, const bool a_corrector){
  CH_TIME("sdc::substep_rte");
  if(m_verbosity > 5){
    pout() << "sdc::substep_rte" << endl;
  }

  const Real dtm = m_dtm[a_m];

  for (RtIterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = get_rte_storage(solver_it);


    EBAMRCellData& phi_m1      = storage->get_phi()[a_m+1];
    const EBAMRCellData& phi_m = storage->get_phi()[a_m];
    
    if(!m_freeze_photoi){
      DataOps::copy(phi_m1, phi_m); // Faster MG solve (if MG is used)
      solver->advance(dtm, phi_m1);
    }
    else if(m_freeze_photoi && !a_corrector){
      DataOps::copy(phi_m1, phi_m); // Faster MG solve (if MG is used)
      solver->advance(dtm, phi_m1);
    }
  }
}

void sdc::reconcile_integrands(){
  CH_TIME("sdc::reconcile_integrands");
  if(m_verbosity > 5){
    pout() << "sdc::reconcile_integrands" << endl;
  }

  // TLDR: When we come in here, all solutions (Poisson, CDR, RTE, Sigma) are known at node m_p. But we
  //       do need the extra slopes for the explicit operators so we have to compute divF and divD again. 

  Vector<EBAMRCellData*> cdr_densities_p = sdc::get_cdr_phik(m_p);
  EBAMRIVData& sigma_p = sdc::get_sigmak(m_p);
  const Real t_p = m_tm[m_p];

  // Update boundary conditions for cdr and sigma equations before getting the slope at the final node
  sdc::compute_cdr_eb_states(cdr_densities_p);
  sdc::compute_cdr_fluxes(cdr_densities_p, t_p);
  sdc::compute_cdr_domain_states(cdr_densities_p);
  sdc::compute_cdr_domain_fluxes(cdr_densities_p, t_p);
  sdc::compute_sigma_flux();

  // Now compute FAR_p - that wasn't done when we integrated
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = sdc::get_cdr_storage(solver_it);
    const int idx = solver_it.get_solver();

    // This has not been computed yet. Do it.
    EBAMRCellData& FA_p       = storage->get_FA()[m_p];
    EBAMRCellData& FD_p       = storage->get_FD()[m_p];
    const EBAMRCellData& phi_p = *cdr_densities_p[idx] ;

    if(solver->isMobile()){
      const Real extrap_dt = m_extrap_advect ? m_dtm[m_p-1] : 0.0; // Factor of 2 because of EBPatchAdvect
      solver->computeDivF(FA_p, phi_p, extrap_dt);                // FAR_p =  Div(v_p*phi_p)
      DataOps::scale(FA_p, -1.0);                                 // FAR_p = -Div(v_p*phi_p)
    }
    else{
      DataOps::setValue(FA_p, 0.0);
    }

    if(solver->isDiffusive()){
      solver->computeDivD(FD_p, phi_p);
    }
    else{
      DataOps::setValue(FD_p, 0.0);
    }

    // Build the integrand required for the lagged high-order quadrature
    for (int m = 0; m <= m_p; m++){
      EBAMRCellData& F_m        = storage->get_F()[m];
      const EBAMRCellData& FA_m = storage->get_FA()[m];
      const EBAMRCellData& FD_m = storage->get_FD()[m];
      const EBAMRCellData& FR_m = storage->get_FR()[m];

      DataOps::copy(F_m, FR_m);
      if(solver->isDiffusive()){
	DataOps::incr(F_m, FD_m, 1.0);
      }
      if(solver->isMobile()){
	DataOps::incr(F_m, FA_m, 1.0);
      }

      // Shouldn't be necessary
      m_amr->averageDown(F_m, m_cdr->getPhase());
      m_amr->interpGhost(F_m, m_cdr->getPhase());
    }
  }

  // Compute Fsig_p - that wasn't done either
  EBAMRIVData& Fnew_p = m_sigma_scratch->get_Fnew()[m_p];
  m_sigma->computeRHS(Fnew_p);
  for (int m = 0; m <= m_p; m++){
    EBAMRIVData& Fold_m = m_sigma_scratch->get_Fold()[m];
    EBAMRIVData& Fnew_m = m_sigma_scratch->get_Fnew()[m];
    DataOps::copy(Fold_m, Fnew_m);
  }
}

void sdc::sweep_explicit(const Real a_dt, const Real a_time, const bool a_lagged_terms){
  CH_TIME("sdc::sweep_explicit");
  if(m_verbosity > 5){
    pout() << "sdc::sweep_explicit" << endl;
  }

  MayDay::Abort("sdc::sweep_explicit - not implemented yet");
}

void sdc::initialize_errors(){
  CH_TIME("sdc::corrector_initialize_errors");
  if(m_verbosity > 5){
    pout() << "sdc::corrector_initialize_errors" << endl;
  }

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = sdc::get_cdr_storage(solver_it);
    const int idx = solver_it.get_solver();
    
    // These should be zero
    if(idx == m_error_idx || m_error_idx < 0){
      EBAMRCellData& error = storage->get_error();
      const EBAMRCellData& phi_final = storage->get_phi()[m_p];

      DataOps::setValue(error, 0.0);
      DataOps::incr(error, phi_final, -1.0);
    }
  }

  EBAMRIVData& error = m_sigma_scratch->get_error();
  const EBAMRIVData& sigma_final = m_sigma_scratch->get_sigma()[m_p];
  DataOps::setValue(error, 0.0);
  DataOps::incr(error, sigma_final, -1.0);
}

void sdc::finalize_errors(){
  CH_TIME("sdc::corrector_finalize_errors");
  if(m_verbosity > 5){
    pout() << "sdc::corrector_finalize_errors" << endl;
  }

  const Real safety = 1.E-20;

  m_max_error = 0.0;
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = sdc::get_cdr_storage(solver_it);
    const int idx = solver_it.get_solver();

    // Compute error
    if(idx == m_error_idx || m_error_idx < 0){
      EBAMRCellData& error       = storage->get_error();
      const EBAMRCellData& phi_p = storage->get_phi()[m_p];
      DataOps::incr(error, phi_p, 1.0);

      // Compute norms. Only coarsest level
      Real Lerr, Lphi;
      const int lvl = 0;
      DataOps::norm(Lerr, *error[lvl], m_amr->getDomains()[lvl], m_error_norm);
      DataOps::norm(Lphi, *phi_p[lvl], m_amr->getDomains()[lvl], m_error_norm);

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
  DataOps::incr(error, sigma_final, 1.0);
  m_sigma_error = 0.0; // I don't think this is ever used...


}

void sdc::compute_new_dt(bool& a_accept_step, const Real a_dt, const int a_num_corrections){
  CH_TIME("sdc::compute_new_dt");
  if(m_verbosity > 5){
    pout() << "sdc::compute_new_dt" << endl;
  }

  // If a_dt was the smallest possible CFL or hardcap time step, we just have to accept it
  const Real max_gl_dist = sdc::get_max_node_distance();
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

void sdc::adaptive_report(const Real a_first_dt, const Real a_dt, const Real a_new_dt, const int a_corr, const int a_rej, const Real a_max_err){
  CH_TIME("sdc::adaptive_report");
  if(m_verbosity > 5){
    pout() << "sdc::adaptive_report" << endl;
  }

  pout() << "\n";
  pout() << "sdc::adaptive_report breakdown" << endl;
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

void sdc::computeDt(Real& a_dt, TimeCode::which_code& a_timeCode){
  CH_TIME("sdc::computeDt");
  if(m_verbosity > 5){
    pout() << "sdc::computeDt" << endl;
  }

  Real dt = 1.E99;

  int Nref = 1;
  for (int lvl = 0; lvl < m_amr->getFinestLevel(); lvl++){
    Nref = Nref*m_amr->getRefinementRatios()[lvl];
  }
  const Real max_gl_dist = sdc::get_max_node_distance();
  m_dt_cfl = m_cdr->compute_cfl_dt();

  Real dt_cfl = 2.0*m_dt_cfl/max_gl_dist;
  
  // Time step selection for non-adaptive stepping
  if(!m_adaptive_dt){
    if(dt_cfl < dt){
      dt = m_cfl*dt_cfl;
      a_timeCode = TimeCode::cfl;
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
      a_timeCode = TimeCode::Error;
    }
  }

  // EVERYTHING BELOW HERE IS "STANDARD"
  // -----------------------------------
  const Real dt_relax = m_relax_time*this->compute_relaxation_time();
  if(dt_relax < dt){
    dt = dt_relax;
    a_timeCode = TimeCode::RelaxationTime;
  }

  const Real dt_restrict = this->restrict_dt();
  if(dt_restrict < dt){
    dt = dt_restrict;
    a_timeCode = TimeCode::Restricted;
  }

  if(dt < m_min_dt){
    dt = m_min_dt;
    a_timeCode = TimeCode::Hardcap;
  }

  if(dt > m_max_dt){
    dt = m_max_dt;
    a_timeCode = TimeCode::Hardcap;
  }

  a_dt = dt;

  // Copy the time code, it is needed for diagnostics
  m_timeCode = a_timeCode;

#if 0 // Debug
  if(procID() == 0){
    std::cout << "computeDt = " << a_dt << "\t m_new_dt = " << m_new_dt << std::endl; 
  }
#endif
}

void sdc::cache_internals(){
  CH_TIME("sdc::cache_internals");
  if(m_verbosity > 5){
    pout() << "sdc::cache_internals" << endl;
  }

  m_cache_FR0.resize(m_plaskin->get_num_species());

  //
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_storage>& storage = sdc::get_cdr_storage(solver_it);
    const int idx = solver_it.get_solver();
    m_amr->allocate(m_cache_FR0[idx], phase::gas, 1);

    // Storage the reactive slope
    const EBAMRCellData& FR0 = storage->get_FR()[0];
    DataOps::copy(m_cache_FR0[idx], FR0);
  }
}

void sdc::allocateInternals(){
  CH_TIME("sdc::allocateInternals");
  if(m_verbosity > 5){
    pout() << "sdc::allocateInternals" << endl;
  }
  
  m_cdr_error.resize(m_plaskin->get_num_species());
  m_dummy_rte.resize(m_plaskin->get_num_Photons());

  for (RtIterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_dummy_rte[idx] = new EBAMRCellData();
    m_amr->allocate(*m_dummy_rte[idx], phase::gas, 1);
  }

  m_amr->allocate(m_scratch1,  phase::gas, 1);
  m_amr->allocate(m_scratchD,  phase::gas, SpaceDim);
  
  sdc::allocate_cdr_storage();
  sdc::allocate_poisson_storage();
  sdc::allocate_rte_storage();
  sdc::allocate_sigma_storage();

  sdc::setup_quadrature_nodes(m_p);
  sdc::setup_qmj(m_p);
}

void sdc::regridInternals(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){
  CH_TIME("sdc::regridInternals");
  if(m_verbosity > 5){
    pout() << "sdc::regridInternals" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;
  const Interval interv(comp, comp);

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_storage>& storage = sdc::get_cdr_storage(solver_it);
    const int idx = solver_it.get_solver();

    // This is data on the new grids which needs to be filled
    EBAMRCellData& FR0 = storage->get_FR()[0];


    // Get the interpolator from amr
    Vector<RefCountedPtr<EBPWLFineInterp> >& interpolator = m_amr->getPwlInterpolator(phase::gas);
    

    // These levels have not changed and can be copied
    for (int lvl = 0; lvl <= Max(0,a_lmin-1); lvl++){
      m_cache_FR0[idx][lvl]->copyTo(*FR0[lvl]); // Base level should never change, but ownership can.
    }

    // These levels have changed and need to be interpolated. End by copying regions that didn't change
    for (int lvl = a_lmin; lvl <= a_newFinestLevel; lvl++){
      interpolator[lvl]->interpolate(*FR0[lvl], *FR0[lvl-1], interv);
      if(lvl <= Min(a_oldFinestLevel, a_newFinestLevel)){
	m_cache_FR0[idx][lvl]->copyTo(*FR0[lvl]);
      }
    }

    m_amr->averageDown(FR0, phase::gas);
    m_amr->interpGhost(FR0, phase::gas);


    m_amr->deallocate(m_cache_FR0[idx]);
  }


}

void sdc::allocate_cdr_storage(){
  const int ncomp       = 1;
  const int num_species = m_plaskin->get_num_species();

  m_cdr_scratch.resize(num_species);
  
  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx] = RefCountedPtr<cdr_storage> (new cdr_storage(m_amr, m_cdr->getPhase(), ncomp));
    m_cdr_scratch[idx]->allocate_storage(m_p);
  }
}

void sdc::allocate_poisson_storage(){
  const int ncomp = 1;
  m_fieldSolver_scratch = RefCountedPtr<poisson_storage> (new poisson_storage(m_amr, m_cdr->getPhase(), ncomp));
  m_fieldSolver_scratch->allocate_storage(m_p);
}

void sdc::allocate_rte_storage(){
  const int ncomp       = 1;
  const int num_Photons = m_plaskin->get_num_Photons();
  m_rte_scratch.resize(num_Photons);
  
  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx] = RefCountedPtr<rte_storage> (new rte_storage(m_amr, m_rte->getPhase(), ncomp));
    m_rte_scratch[idx]->allocate_storage(m_p);
  }
}

void sdc::allocate_sigma_storage(){
  const int ncomp = 1;
  m_sigma_scratch = RefCountedPtr<sigma_storage> (new sigma_storage(m_amr, m_cdr->getPhase(), ncomp));
  m_sigma_scratch->allocate_storage(m_p);
}

void sdc::deallocateInternals(){
  CH_TIME("sdc::deallocateInternals");
  if(m_verbosity > 5){
    pout() << "sdc::deallocateInternals" << endl;
  }

  m_amr->deallocate(m_scratch1);
  m_amr->deallocate(m_scratchD);

  for (RtIterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_amr->deallocate(*m_dummy_rte[idx]);
  }

  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx]->deallocate_storage();
    m_cdr_scratch[idx] = RefCountedPtr<cdr_storage>(0);
  }

  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx]->deallocate_storage();
    m_rte_scratch[idx] = RefCountedPtr<rte_storage>(0);
  }

  m_dummy_rte.resize(0);
  m_cdr_scratch.resize(0);
  m_rte_scratch.resize(0);

  m_fieldSolver_scratch->deallocate_storage();
  m_fieldSolver_scratch = RefCountedPtr<poisson_storage>(0);
  
  m_sigma_scratch->deallocate_storage();
  m_sigma_scratch = RefCountedPtr<sigma_storage>(0);
}

void sdc::compute_E_into_scratch(){
  CH_TIME("sdc::compute_E_into_scratch");
  if(m_verbosity > 5){
    pout() << "sdc::compute_E_into_scratch" << endl;
  }
  
  EBAMRCellData& E_cell = m_fieldSolver_scratch->get_E_cell();
  EBAMRFluxData& E_face = m_fieldSolver_scratch->get_E_face();
  EBAMRIVData&   E_eb   = m_fieldSolver_scratch->get_E_eb();
  EBAMRIFData&   E_dom  = m_fieldSolver_scratch->get_E_domain();

  const MFAMRCellData& phi = m_fieldSolver->getPotential();
  
  sdc::compute_E(E_cell, m_cdr->getPhase(), phi);     // Compute cell-centered field
  sdc::compute_E(E_face, m_cdr->getPhase(), E_cell);  // Compute face-centered field
  sdc::compute_E(E_eb,   m_cdr->getPhase(), E_cell);  // EB-centered field

  TimeStepper::extrapolate_to_domain_faces(E_dom, m_cdr->getPhase(), E_cell);
}

void sdc::compute_cdr_gradients(){
  CH_TIME("sdc::compute_cdr_gradients");
  if(m_verbosity > 5){
    pout() << "sdc::compute_cdr_gradients" << endl;
  }

  sdc::compute_cdr_gradients(m_cdr->getPhis());
}

void sdc::compute_cdr_gradients(const Vector<EBAMRCellData*>& a_phis){
  CH_TIME("sdc::compute_cdr_gradients");
  if(m_verbosity > 5){
    pout() << "sdc::compute_cdr_gradients" << endl;
  }

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    RefCountedPtr<cdr_storage>& storage = sdc::get_cdr_storage(solver_it);
    EBAMRCellData& grad = storage->get_gradient();
    m_amr->computeGradient(grad, *a_phis[idx], m_cdr->getPhase());
    //    m_amr->averageDown(grad, m_cdr->getPhase());
    m_amr->interpGhost(grad, m_cdr->getPhase());
  }
}

void sdc::compute_cdr_velo(const Real a_time){
  CH_TIME("sdc::compute_cdr_velo");
  if(m_verbosity > 5){
    pout() << "sdc::compute_cdr_velo" << endl;
  }

  sdc::compute_cdr_velo(m_cdr->getPhis(), a_time);
}

void sdc::compute_cdr_velo(const Vector<EBAMRCellData*>& a_phis, const Real a_time){
  CH_TIME("sdc::compute_cdr_velo(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "sdc::compute_cdr_velo(Vector<EBAMRCellData*>, Real)" << endl;
  }

  Vector<EBAMRCellData*> velocities = m_cdr->getVelocities();
  sdc::compute_cdr_velocities(velocities, a_phis, m_fieldSolver_scratch->get_E_cell(), a_time);
}

void sdc::compute_cdr_eb_states(){
  CH_TIME("sdc::compute_cdr_eb_states");
  if(m_verbosity > 5){
    pout() << "sdc::compute_cdr_eb_states" << endl;
  }

  Vector<EBAMRIVData*>   eb_gradients;
  Vector<EBAMRIVData*>   eb_states;
  Vector<EBAMRCellData*> cdr_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    RefCountedPtr<cdr_storage>& storage = sdc::get_cdr_storage(solver_it);

    cdr_states.push_back(&(solver->getPhi()));
    eb_states.push_back(&(storage->get_eb_state()));
    eb_gradients.push_back(&(storage->get_eb_grad()));
    cdr_gradients.push_back(&(storage->get_gradient())); // Should already have been computed
  }

  // Extrapolate states to the EB and floor them so we cannot get negative values on the boundary. This
  // won't hurt mass conservation because the mass hasn't been injected yet
  sdc::extrapolate_to_eb(eb_states, m_cdr->getPhase(), cdr_states);
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    DataOps::floor(*eb_states[idx], 0.0);
  }

  // We should already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIVData eb_gradient;
  m_amr->allocate(eb_gradient, m_cdr->getPhase(), SpaceDim);
  for (int i = 0; i < cdr_states.size(); i++){
    sdc::extrapolate_to_eb(eb_gradient, m_cdr->getPhase(), *cdr_gradients[i]);
    sdc::project_flux(*eb_gradients[i], eb_gradient);
  }
}

void sdc::compute_cdr_eb_states(const Vector<EBAMRCellData*>& a_phis){
  CH_TIME("sdc::compute_cdr_eb_states(vec)");
  if(m_verbosity > 5){
    pout() << "sdc::compute_cdr_eb_states(vec)" << endl;
  }

  Vector<EBAMRIVData*>   eb_gradients;
  Vector<EBAMRIVData*>   eb_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = sdc::get_cdr_storage(solver_it);

    eb_states.push_back(&(storage->get_eb_state()));
    eb_gradients.push_back(&(storage->get_eb_grad()));
    cdr_gradients.push_back(&(storage->get_gradient())); // Should already have been computed
  }

  // Extrapolate states to the EB and floor them so we cannot get negative values on the boundary. This
  // won't hurt mass conservation because the mass hasn't been injected yet
  sdc::extrapolate_to_eb(eb_states, m_cdr->getPhase(), a_phis);
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    DataOps::floor(*eb_states[idx], 0.0);
  }

  // We should already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIVData eb_gradient;
  m_amr->allocate(eb_gradient, m_cdr->getPhase(), SpaceDim);
  for (int i = 0; i < a_phis.size(); i++){
    sdc::extrapolate_to_eb(eb_gradient, m_cdr->getPhase(), *cdr_gradients[i]);
    sdc::project_flux(*eb_gradients[i], eb_gradient);
  }
}

void sdc::compute_cdr_domain_states(){
  CH_TIME("sdc::compute_cdr_domain_states");
  if(m_verbosity > 5){
    pout() << "sdc::compute_cdr_domain_states" << endl;
  }

  Vector<EBAMRIFData*>   domain_gradients;
  Vector<EBAMRIFData*>   domain_states;
  Vector<EBAMRCellData*> cdr_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    RefCountedPtr<cdr_storage>& storage = sdc::get_cdr_storage(solver_it);

    cdr_states.push_back(&(solver->getPhi()));
    domain_states.push_back(&(storage->getDomain_state()));
    domain_gradients.push_back(&(storage->getDomain_grad()));
    cdr_gradients.push_back(&(storage->get_gradient())); // Should already be computed
  }

  // Extrapolate states to the domain faces
  sdc::extrapolate_to_domain_faces(domain_states, m_cdr->getPhase(), cdr_states);

  // We already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIFData grad;
  m_amr->allocate(grad, m_cdr->getPhase(), SpaceDim);
  for (int i = 0; i < cdr_states.size(); i++){
    sdc::extrapolate_to_domain_faces(grad, m_cdr->getPhase(), *cdr_gradients[i]);
    sdc::project_domain(*domain_gradients[i], grad);
  }
}

void sdc::compute_cdr_domain_states(const Vector<EBAMRCellData*>& a_phis){
  CH_TIME("sdc::compute_cdr_domain_states");
  if(m_verbosity > 5){
    pout() << "sdc::compute_cdr_domain_states" << endl;
  }

  Vector<EBAMRIFData*>   domain_gradients;
  Vector<EBAMRIFData*>   domain_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    domain_states.push_back(&(storage->getDomain_state()));
    domain_gradients.push_back(&(storage->getDomain_grad()));
    cdr_gradients.push_back(&(storage->get_gradient()));
  }

  // Extrapolate states to the domain faces
  this->extrapolate_to_domain_faces(domain_states, m_cdr->getPhase(), a_phis);

  // We already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIFData grad;
  m_amr->allocate(grad, m_cdr->getPhase(), SpaceDim);
  for (int i = 0; i < a_phis.size(); i++){
    this->extrapolate_to_domain_faces(grad, m_cdr->getPhase(), *cdr_gradients[i]);
    this->project_domain(*domain_gradients[i], grad);
  }
}

void sdc::compute_cdr_fluxes(const Real a_time){
  CH_TIME("sdc::compute_cdr_fluxes");
  if(m_verbosity > 5){
    pout() << "sdc::compute_cdr_fluxes" << endl;
  }

  this->compute_cdr_fluxes(m_cdr->getPhis(), a_time);
}

void sdc::compute_cdr_fluxes(const Vector<EBAMRCellData*>& a_phis, const Real a_time){
  CH_TIME("sdc::compute_cdr_fluxes(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "sdc::compute_cdr_fluxes(Vector<EBAMRCellData*>, Real)" << endl;
  }

  Vector<EBAMRIVData*> cdr_fluxes;
  Vector<EBAMRIVData*> extrap_cdr_fluxes;
  Vector<EBAMRIVData*> extrap_cdr_densities;
  Vector<EBAMRIVData*> extrap_cdr_velocities;
  Vector<EBAMRIVData*> extrap_cdr_gradients;
  Vector<EBAMRIVData*> extrap_rte_fluxes;

  cdr_fluxes = m_cdr->getEbFlux();

  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
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
  Vector<EBAMRCellData*> cdr_velocities = m_cdr->getVelocities();
  TimeStepper::compute_extrapolated_fluxes(extrap_cdr_fluxes, a_phis, cdr_velocities, m_cdr->getPhase());
  TimeStepper::compute_extrapolated_velocities(extrap_cdr_velocities, cdr_velocities, m_cdr->getPhase());

  // Compute RTE flux on the boundary
  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);

    EBAMRIVData& flux_eb = storage->get_eb_flux();
    solver->computeBoundaryFlux(flux_eb, solver->getPhi());
    extrap_rte_fluxes.push_back(&flux_eb);
  }

  const EBAMRIVData& E = m_fieldSolver_scratch->get_E_eb();
  TimeStepper::compute_cdr_fluxes(cdr_fluxes,
				   extrap_cdr_fluxes,
				   extrap_cdr_densities,
				   extrap_cdr_velocities,
				   extrap_cdr_gradients,
				   extrap_rte_fluxes,
				   E,
				   a_time);
}

void sdc::compute_cdr_domain_fluxes(const Real a_time){
  CH_TIME("sdc::compute_cdr_domain_fluxes");
  if(m_verbosity > 5){
    pout() << "sdc::compute_cdr_domain_fluxes" << endl;
  }

  this->compute_cdr_domain_fluxes(m_cdr->getPhis(), a_time);
}

void sdc::compute_cdr_domain_fluxes(const Vector<EBAMRCellData*>& a_phis, const Real a_time){
  CH_TIME("sdc::compute_cdr_domain_fluxes(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "sdc::compute_cdr_domain_fluxes(Vector<EBAMRCellData*>, Real)" << endl;
  }

  Vector<EBAMRIFData*>   cdr_fluxes;
  Vector<EBAMRIFData*>   extrap_cdr_fluxes;
  Vector<EBAMRIFData*>   extrap_cdr_densities;
  Vector<EBAMRIFData*>   extrap_cdr_velocities;
  Vector<EBAMRIFData*>   extrap_cdr_gradients;
  Vector<EBAMRIFData*>   extrap_rte_fluxes;

  Vector<EBAMRCellData*> cdr_velocities;
  Vector<EBAMRCellData*> cdr_gradients;

  cdr_fluxes = m_cdr->getDomainFlux();
  cdr_velocities = m_cdr->getVelocities();
  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRIFData& dens_domain = storage->getDomain_state();
    EBAMRIFData& velo_domain = storage->getDomain_velo();
    EBAMRIFData& flux_domain = storage->getDomain_flux();
    EBAMRIFData& grad_domain = storage->getDomain_grad();
    EBAMRCellData& gradient  = storage->get_gradient();

    extrap_cdr_densities.push_back(&dens_domain);  // Has not been computed
    extrap_cdr_velocities.push_back(&velo_domain); // Has not been computed
    extrap_cdr_fluxes.push_back(&flux_domain);     // Has not been computed
    extrap_cdr_gradients.push_back(&grad_domain);  // Has not been computed
    cdr_gradients.push_back(&gradient);
  }

  // Compute extrapolated velocities and fluxes at the domain faces
  this->extrapolate_to_domain_faces(extrap_cdr_densities,         m_cdr->getPhase(), a_phis);
  this->extrapolate_velo_to_domain_faces(extrap_cdr_velocities, m_cdr->getPhase(), cdr_velocities);
  this->compute_extrapolated_domain_fluxes(extrap_cdr_fluxes,     a_phis,           cdr_velocities, m_cdr->getPhase());
  this->extrapolate_vector_to_domain_faces(extrap_cdr_gradients,  m_cdr->getPhase(), cdr_gradients);

  // Compute RTE flux on domain faces
  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);

    EBAMRIFData& domain_flux = storage->getDomain_flux();
    solver->computeDomainFlux(domain_flux, solver->getPhi());
    extrap_rte_fluxes.push_back(&domain_flux);
  }

  const EBAMRIFData& E = m_fieldSolver_scratch->get_E_domain();

  // This fills the solvers' domain fluxes
  TimeStepper::compute_cdr_domain_fluxes(cdr_fluxes,
					  extrap_cdr_fluxes,
					  extrap_cdr_densities,
					  extrap_cdr_velocities,
					  extrap_cdr_gradients,
					  extrap_rte_fluxes,
					  E,
					  a_time);
}

void sdc::compute_sigma_flux(){
  CH_TIME("sdc::compute_sigma_flux");
  if(m_verbosity > 5){
    pout() << "sdc::compute_sigma_flux" << endl;
  }

  EBAMRIVData& flux = m_sigma->getFlux();
  DataOps::setValue(flux, 0.0);

  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const RefCountedPtr<species>& spec      = solver_it.getSpecies();
    const EBAMRIVData& solver_flux          = solver->getEbFlux();

    DataOps::incr(flux, solver_flux, spec->getChargeNumber()*Units::Qe);
  }

  m_sigma->resetCells(flux);
}

void sdc::compute_reaction_network(const int a_m, const Real a_time, const Real a_dt, const bool a_corrector){
  CH_TIME("sdc::compute_reaction_network");
  if(m_verbosity > 5){
    pout() << "sdc::compute_reaction_network";
  }

  Vector<EBAMRCellData*> cdr_sources = m_cdr->getSources();
  Vector<EBAMRCellData*> rte_sources;

  // If we freeze the photoionization profile in the first sweep, avoid filling the solver source terms
  if(m_freeze_photoi && a_corrector){
    if(procID() == 0) std::cout << "filling dummy" << std::endl;
    rte_sources = m_dummy_rte;
  }
  else{
    if(procID() == 0) std::cout << "filling sources" << std::endl;
    rte_sources = m_rte->getSources();
  }

  const Vector<EBAMRCellData*> cdr_densities = get_cdr_phik(a_m);
  const Vector<EBAMRCellData*> rte_densities = get_rte_phik(a_m);
  const EBAMRCellData& E = m_fieldSolver_scratch->get_E_cell();


  TimeStepper::advance_reaction_network(cdr_sources, rte_sources, cdr_densities, rte_densities, E, a_time, a_dt);
}

void sdc::update_poisson(){
  CH_TIME("sdc::update_poisson(solver)");
  if(m_verbosity > 5){
    pout() << "sdc::update_poisson(solver)" << endl;
  }
  
  if(m_do_poisson){ // Solve Poisson equation
    if((m_timeStep +1) % m_fast_poisson == 0){
      TimeStepper::solve_poisson();
      this->compute_E_into_scratch();
    }
  }
}

void sdc::update_poisson(const Vector<EBAMRCellData*>& a_densities, const EBAMRIVData& a_sigma){
  CH_TIME("sdc::update_poisson(full)");
  if(m_verbosity > 5){
    pout() << "sdc::update_poisson(full)" << endl;
  }
  
  if(m_do_poisson){ // Solve Poisson equation
    if((m_timeStep +1) % m_fast_poisson == 0){
      TimeStepper::solve_poisson(m_fieldSolver->getPotential(),
				  m_fieldSolver->getRho(),
				  a_densities,
				  a_sigma,
				  centering::cell_center);
      this->compute_E_into_scratch();
    }
  }
}

void sdc::integrate_rte_transient(const Real a_dt){
  CH_TIME("sdc::integrate_rte_transient");
  if(m_verbosity > 5){
    pout() << "sdc::integrate_rte_transient" << endl;
  }

  if(m_do_rte){
    if((m_timeStep + 1) % m_fast_rte == 0){
      if(!(m_rte->isStationary())){
	for (RtIterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
	  RefCountedPtr<RtSolver>& solver = solver_it();
	  solver->advance(a_dt);
	}
      }
    }
  }
}

void sdc::integrate_rte_stationary(){
  CH_TIME("sdc::integrate_rte_transient");
  if(m_verbosity > 5){
    pout() << "sdc::integrate_rte_transient" << endl;
  }

  if(m_do_rte){
    if((m_timeStep + 1) % m_fast_rte == 0){
      if((m_rte->isStationary())){
	for (RtIterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
	  RefCountedPtr<RtSolver>& solver = solver_it();
	  solver->advance(0.0);
	}
      }
    }
  }
}

void sdc::update_diffusion_coefficients(){
  CH_TIME("sdc::update_diffusion_coefficients");
  if(m_verbosity > 5){
    pout() << "sdc::update_diffusion_coefficients" << endl;
  }
  TimeStepper::compute_cdr_diffusion(m_fieldSolver_scratch->get_E_cell(), m_fieldSolver_scratch->get_E_eb());
}

Vector<EBAMRCellData*> sdc::get_cdr_errors(){
  CH_TIME("sdc::get_cdr_errors");
  if(m_verbosity > 5){
    pout() << "sdc::get_cdr_errors" << endl;
  }
  
  Vector<EBAMRCellData*> ret;
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = sdc::get_cdr_storage(solver_it);
    ret.push_back(&(storage->get_error()));
  }

  return ret;
}

Vector<EBAMRCellData*> sdc::get_cdr_phik(const int a_m){
  CH_TIME("sdc::get_cdr_phik");
  if(m_verbosity > 5){
    pout() << "sdc::get_cdr_phik" << endl;
  }
  
  Vector<EBAMRCellData*> ret;
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = sdc::get_cdr_storage(solver_it);
    ret.push_back(&(storage->get_phi()[a_m]));
  }

  return ret;
}

Vector<EBAMRCellData*> sdc::get_rte_phik(const int a_m){
  CH_TIME("sdc::get_rte_phik");
  if(m_verbosity > 5){
    pout() << "sdc::get_rte_phik" << endl;
  }
  
  Vector<EBAMRCellData*> ret;
  for (RtIterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_storage>& storage = sdc::get_rte_storage(solver_it);
    ret.push_back(&(storage->get_phi()[a_m]));
  }

  return ret;
}

EBAMRIVData& sdc::get_sigmak(const int a_m){
  CH_TIME("sdc::get_sigmak");
  if(m_verbosity > 5){
    pout() << "sdc::get_sigmak)" << endl;
  }
  return m_sigma_scratch->get_sigma()[a_m];
}

void sdc::write_step_profile(const Real a_dt,
			     const Real a_error,
			     const int  a_substeps,
			     const int  a_corrections,
			     const int  a_rejections){
  CH_TIME("sissdc::write_step_profile");
  if(m_verbosity > 5){
    pout() << "sdc::write_step_profile" << endl;
  }

  if(procID() == 0 ){

    const std::string fname("sdc_step_profile.txt");
    
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

    f << std::left << std::setw(width) << m_timeStep << "\t"
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

void sdc::store_solvers(){
  CH_TIME("sdc::store_solvers");
  if(m_verbosity > 5){
    pout() << "sdc::store_solvers" << endl;
  }

  if(m_k > 0 && m_adaptive_dt){
    // SDC does not manipulate cdr and sigma solvers until the end of the time step. Only need to do
    // Poisson and RTE here.

    // Poisson
    MFAMRCellData& previous    = m_fieldSolver_scratch->get_previous();
    const MFAMRCellData& state = m_fieldSolver->getPotential();
    DataOps::copy(previous, state);

    // RTE
    for (RtIterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<rte_storage>& storage     = sdc::get_rte_storage(solver_it);
      const RefCountedPtr<RtSolver>& solver = solver_it();

      EBAMRCellData& previous = storage->get_previous();
      const EBAMRCellData& state = solver->getPhi();

      DataOps::copy(previous, state);
    }
  }
}

void sdc::restore_solvers(){
  CH_TIME("sdc::restore_solvers");
  if(m_verbosity > 5){
    pout() << "sdc::restore_solvers" << endl;
  }

  // SDC does not manipulate cdr and sigma solvers until the end of the time step. Only need to do
  // Poisson and RTE here. 

  // Poisson
  MFAMRCellData& state = m_fieldSolver->getPotential();
  const MFAMRCellData& previous    = m_fieldSolver_scratch->get_previous();

  DataOps::copy(state, previous);

  // RTE
  for (RtIterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_storage>& storage     = sdc::get_rte_storage(solver_it);
    RefCountedPtr<RtSolver>& solver = solver_it();

    EBAMRCellData& previous = storage->get_previous();
    EBAMRCellData& state = solver->getPhi();

    DataOps::copy(state, previous);
  }
#include "CD_NamespaceFooter.H"
