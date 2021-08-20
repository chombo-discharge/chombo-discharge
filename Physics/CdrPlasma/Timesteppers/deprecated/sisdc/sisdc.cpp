/*!
  @file   sisdc.cpp
  @brief  Implementation of sisdc.H
  @author Robert Marskar
  @date   Feb. 2019
*/

#include "sisdc.H"
#include "sisdcF_F.H"
#include "sisdc_storage.H"
#include <CD_DataOps.H>
#include <CD_Units.H>
#include <CD_CdrGodunov.H>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <ParmParse.H>

typedef sisdc::CdrStorage     CdrStorage;
typedef sisdc::FieldStorage FieldStorage;
typedef sisdc::RtStorage     RtStorage;
typedef sisdc::SigmaStorage   SigmaStorage;

sisdc::sisdc(){
  m_className = "sisdc";
}

sisdc::~sisdc(){

}

void sisdc::parseOptions(){
  CH_TIME("sisdc::parseOptions");
  if(m_verbosity > 5){
    pout() << "sisdc::parseOptions" << endl;
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
  parse_source_comp();

  // Specific to this class
  parseNodes();
  parseAP();
  parseDiffusionCoupling();
  parseAdaptiveOptions();
  parseSubcycleOptions();
  parseDebugOptions();
  parseAdvectionOptions();

  // Setup nodes
  sisdc::setupQuadratureNodes(m_p);
  sisdc::setupQmj(m_p);
}

void sisdc::parseNodes(){
  CH_TIME("sisdc::parseNodes");
  if(m_verbosity > 5){
    pout() << "sisdc::parseNodes" << endl;
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
    MayDay::Abort("sisdc::parseNodes - unknown node type requested");
  }

  pp.get("subintervals",     m_p);
  pp.get("corr_iter",        m_k);

  if(m_p < 1){
    MayDay::Abort("sisdc::parseNodes - sisdc.subintervals cannot be < 1");
  }
  if(m_k < 0){
    MayDay::Abort("sisdc::parseNodes - sisdc.corr_iter cannot be < 0");
  }
}

void sisdc::parseAP(){
  CH_TIME("sisdc::parseAP");
  if(m_verbosity > 5){
    pout() << "sisdc::parseAP" << endl;
  }

  ParmParse pp(m_className.c_str());
  
  std::string str;

  pp.get("use_AP", str);
  m_use_AP = (str == "true") ? true : false;
}

void sisdc::parseDiffusionCoupling(){
  CH_TIME("sisdc::parseDiffusionCoupling");
  if(m_verbosity > 5){
    pout() << "sisdc::parseDiffusionCoupling" << endl;
  }

  ParmParse pp(m_className.c_str());
  
  std::string str;

  pp.get("use_tga", str);
  m_useTGA = (str == "true") ? true : false;

  pp.get("diffusive_coupling", str);
  m_strong_diffu = (str == "true") ? true : false;
  
  pp.get("num_diff_corr",  m_num_diff_corr);
  if(m_num_diff_corr < 0){
    MayDay::Abort("sisdc::parseDiffusionCoupling - option 'sisdc.num_diff_corr' cannot be negative");
  }
}

void sisdc::parseAdaptiveOptions(){
  CH_TIME("sisdc::parseAdaptiveOptions");
  if(m_verbosity > 5){
    pout() << "sisdc::parseAdaptiveOptions" << endl;
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

void sisdc::parseSubcycleOptions(){
  CH_TIME("sisdc::parseSubcycleOptions");
  if(m_verbosity > 5){
    pout() << "sisdc::parseSubcycleOptions" << endl;
  }

  ParmParse pp(m_className.c_str());
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

void sisdc::parseDebugOptions(){
  CH_TIME("sisdc::parseDebugOptions");
  if(m_verbosity > 5){
    pout() << "sisdc::parseDebugOptions" << endl;
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

void sisdc::parseAdvectionOptions(){
  CH_TIME("sisdc::parseAdvectionOptions");
  if(m_verbosity > 5){
    pout() << "sisdc::parseAdvectionOptions" << endl;
  }

  ParmParse pp(m_className.c_str());
  std::string str;

  m_extrap_dt = 0.5;  // Relic of an ancient past. I don't see any reason why extrapolating to anything but the half interval
                      // would make sense. 

  pp.get("extrap_advect", str); m_extrap_advect = (str == "true") ? true : false;
}

RefCountedPtr<CdrStorage>& sisdc::get_CdrStorage(const CdrIterator& a_solverit){
  return m_cdr_scratch[a_solverit.get_solver()];
}

RefCountedPtr<RtStorage>& sisdc::get_RtStorage(const RtIterator& a_solverit){
  return m_rte_scratch[a_solverit.get_solver()];
}

bool sisdc::needToRegrid(){
  CH_TIME("sisdc::needToRegrid");
  if(m_verbosity > 5){
    pout() << "sisdc::needToRegrid" << endl;
  }
  const bool regrid = m_accum_cfl > m_regrid_cfl;
  if(regrid) m_accum_cfl = 0.0;
  
  return regrid;
}

Real sisdc::restrict_dt(){
  return 1.E99;
}

Real sisdc::getMaxNodeDistance(){
  CH_TIME("sisdc::getMaxNodeDistance");
  if(m_verbosity > 5){
    pout() << "sisdc::getMaxNodeDistance" << endl;
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
  
  TimeStepper::init_source_terms();
}

void sisdc::setupQuadratureNodes(const int a_p){
  CH_TIME("sisdc::setupQuadratureNodes");
  if(m_verbosity > 5){
    pout() << "sisdc::setupQuadratureNodes" << endl;
  }

  if(m_which_nodes == "uniform"){
    sisdc::setupUniformNodes(a_p);
  }
  else if(m_which_nodes == "lobatto"){
    sisdc::setupLobattoNodes(a_p);
  }
  else if(m_which_nodes == "chebyshev"){
    sisdc::setupChebyshevNodes(a_p);
  }
  else {
    MayDay::Abort("sisdc::setupQuadratureNodes - unknown nodes requested");
  }
}

void sisdc::setupUniformNodes(const int a_p){
  CH_TIME("sisdc::setupUniformNodes");
  if(m_verbosity > 5){
    pout() << "sisdc::setupUniformNodes" << endl;
  }

  // TLDR: The nodes and weights are hardcoded. A better programmer would compute these
  //       recursively with Legendre polynomials. 
  m_nodes.resize(1+a_p);

  const Real delta = 2./a_p;
  for (int m = 0; m <= a_p; m++){
    m_nodes[m] = m*delta;
  }
}

void sisdc::setupLobattoNodes(const int a_p){
  CH_TIME("sisdc::setupLobattoNodes");
  if(m_verbosity > 5){
    pout() << "sisdc::setupLobattoNodes" << endl;
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
    MayDay::Abort("sisdc::setupLobattoNodes - requested order exceeds 7. Compute your own damn nodes!");
  }
}

void sisdc::setupChebyshevNodes(const int a_p){
  CH_TIME("sisdc::setupChebyshevNodes");
  if(m_verbosity > 5){
    pout() << "sisdc::setupChebyshevNodes" << endl;
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

void sisdc::setupQmj(const int a_p){
  CH_TIME("sisdc::setupQmj");
  if(m_verbosity > 5){
    pout() << "sisdc::setupQmj" << endl;
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
    if(INFO != 0) MayDay::Abort("sisdc::setupQmj - could not compute weights");
    
    // Now construct qmj
    for (int m = 0; m < a_p; m++){
      m_qmj[m][j] = 0.0;
      for (int k = 0; k < nnodes; k++){
	m_qmj[m][j] += cj[k]*(pow(m_nodes[m+1], k+1) - pow(m_nodes[m], k+1))/(k+1);
      }
    }
  }
}

void sisdc::setupSubintervals(const Real a_time, const Real a_dt){
  CH_TIME("sisdc::setupSubintervals");
  if(m_verbosity > 5){
    pout() << "sisdc::setupSubintervals" << endl;
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

  DataOps::setValue(a_quad, 0.0);
  for (int j = 0; j <= m_p; j++){
    DataOps::incr(a_quad, a_integrand[j], m_qmj[a_m][j]);
  }
}

void sisdc::quad(EBAMRIVData& a_quad, const Vector<EBAMRIVData>& a_integrand, const int a_m){
  CH_TIME("sisdc::quad");
  if(m_verbosity > 5){
    pout() << "sisdc::quad" << endl;
  }

  if(a_m < 0)     MayDay::Abort("sisdc::quad - bad index a_m < 0");
  if(a_m >= m_p)  MayDay::Abort("sisdc::quad - bad index a_m >= m_p");

  DataOps::setValue(a_quad, 0.0);
  for (int j = 0; j <= m_p; j++){
    DataOps::incr(a_quad, a_integrand[j], m_qmj[a_m][j]);
  }
}
  
void sisdc::copyPhiPToCdr(){
  CH_TIME("sisdc::copyPhiPToCdr");
  if(m_verbosity > 5){
    pout() << "sisdc::copyPhiPToCdr" << endl;
  }

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>&  solver  = solver_it();
    RefCountedPtr<CdrStorage>& storage = get_CdrStorage(solver_it);

    EBAMRCellData& phi = solver->getPhi();
    const EBAMRCellData& phip = storage->getPhi()[m_p];
    DataOps::copy(phi, phip);
  }
}

void sisdc::copySigmaPToSigma(){
  CH_TIME("sisdc::copySigmaPToSigma");
  if(m_verbosity > 5){
    pout() << "sisdc::copySigmaPToSigma" << endl;
  }

  EBAMRIVData& sigma        = m_sigma->getPhi();
  const EBAMRIVData& sigmap = m_sigma_scratch->get_sigma()[m_p];
  DataOps::copy(sigma, sigmap);
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
  sisdc::copyCdrToPhiM0();
  sisdc::copySigmaToM0();
  if(m_k > 0) sisdc::computeFD0();
  if(m_k > 0 && m_adaptive_dt) {
    sisdc::storeSolvers();
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
    sisdc::setupSubintervals(m_time, actual_dt);

    // First SDC sweep. No lagged slopes here. 
    if(m_use_AP){
      sisdc::integrate_AP(a_dt, m_time, false);
    }
    else{
      sisdc::integrate(a_dt, m_time, false);
    }

    // SDC correction sweeps. Need to take care of lagged terms. 
    for(int icorr = 0; icorr < Max(m_k, m_min_corr); icorr++){
      num_corrections++;

      // Initialize error and reconcile integrands (i.e. make them quadrature-ready)
      sisdc::initializeErrors();
      sisdc::reconcileIntegrands();

      // SDC correction along whole interval
      if(m_use_AP){
	sisdc::integrate_AP(a_dt, m_time, true);
      }
      else{
	sisdc::integrate(a_dt, m_time, true);
      }

      // Compute error and check if we need to keep iterating
      sisdc::finalizeErrors();
      if(m_max_error < m_err_thresh && m_adaptive_dt && icorr >= m_min_corr) break; // No need in going beyond
    }

    // Compute a new time step. If it is smaller than the minimum allowed CFL step, accept the step anyways
    if(m_adaptive_dt){
      sisdc::computeNewDt(accept_step, actual_dt, num_corrections);
      
      if(!accept_step){  // Step rejection, use the new dt for next step.
#if 0 // Debug
	if(procID() == 0){
	  std::cout << "Rejecting step = " << m_timeStep
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
	  sisdc::restoreSolvers();
	  sisdc::compute_E_into_scratch();
	  sisdc::computeCdrGradients();
	  sisdc::computeCdrVelo(m_time);
	  TimeStepper::compute_cdr_diffusion(m_fieldSolver_scratch->getElectricFieldCell(), m_fieldSolver_scratch->getElectricFieldEb());
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
  sisdc::copyPhiPToCdr();
  sisdc::copySigmaPToSigma();

  const Real t0 = Timer::wallClock();
  sisdc::updateField();
  sisdc::update_stationary_rte(m_time + actual_dt); // Only triggers if m_rte->is_stationar() == true

  // Always recompute source terms and velocities for the next time step. These were computed ea
  const Real t1 = Timer::wallClock();
  sisdc::computeCdrGradients();
  const Real t2 = Timer::wallClock();
  sisdc::computeCdrVelo(m_time + actual_dt);
  const Real t3 = Timer::wallClock();
  TimeStepper::compute_cdr_diffusion(m_fieldSolver_scratch->getElectricFieldCell(), m_fieldSolver_scratch->getElectricFieldEb());
  const Real t4 = Timer::wallClock();

  // In case we're using FHD, we need to tell the kinetics module about the time step before computign sources
  Real next_dt;
  TimeCode::which_code dummy;
  computeDt(next_dt, dummy);
  m_plaskin->setDt(next_dt);
  sisdc::compute_cdr_sources(m_time + actual_dt);
  if(!m_rte->isStationary()){
    
  }
  const Real t5 = Timer::wallClock();

  // Profile step
  if(m_print_report)  sisdc::adaptiveReport(first_dt, actual_dt, m_new_dt, num_corrections, num_reject, m_max_error);
  if(m_profile_steps) sisdc::writeStepProfile(actual_dt, m_max_error, m_p, num_corrections, num_reject);

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

void sisdc::copyCdrToPhiM0(){
  CH_TIME("sisdc::copyCdrToPhiM0");
  if(m_verbosity > 5){
    pout() << "sisdc::copyCdrToPhiM0" << endl;
  }

  // CDR solvers
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>&  solver  = solver_it();
    RefCountedPtr<CdrStorage>& storage = get_CdrStorage(solver_it);
    
    EBAMRCellData& phi0 = storage->getPhi()[0];
    const EBAMRCellData& phi = solver->getPhi();
    DataOps::copy(phi0, phi);
  }
}

void sisdc::copySigmaToM0(){
  CH_TIME("sisdc::copy_sigma_sigma_m0");
  if(m_verbosity > 5){
    pout() << "sisdc::copySigmaToM0" << endl;
  }

  // Copy sigma to starting state
  EBAMRIVData& sigma0      = m_sigma_scratch->get_sigma()[0];
  const EBAMRIVData& sigma = m_sigma->getPhi();
  DataOps::copy(sigma0, sigma);
}

void sisdc::computeFD0(){
  CH_TIME("sisdc::computeFD0");
  if(m_verbosity > 5){
    pout() << "sisdc::computeFD0" << endl;
  }

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = get_CdrStorage(solver_it);
    
    const EBAMRCellData& phi_0 = storage->getPhi()[0]; // phi_0
    EBAMRCellData& FD_0        = storage->getFD()[0];  // FD(phi_0)
    
    if(solver->isDiffusive()){
      //      CdrTGA* tgasolver = (CdrTGA*) (&(*solver));
      //tgasolver->computeDivD(FD_0, phi_0);
      solver->computeDivD(FD_0, phi_0);

      // Shouldn't be necesary
      // m_amr->averageDown(FD_0, m_cdr->getPhase());
      // m_amr->interpGhost(FD_0, m_cdr->getPhase());
    }
    else{
      DataOps::setValue(FD_0, 0.0);
    }
  }
}

void sisdc::integrate(const Real a_dt, const Real a_time, const bool a_lagged_terms){
  CH_TIME("sisdc::integrate");
  if(m_verbosity > 5){
    pout() << "sisdc::integrate" << endl;
  }


  // 1. The first time we enter this routine, source terms and velocities were updated.
  // 2. For further calls, source terms and velocities have been overwritten, but since the explicit
  //    operator slopes do not change, this is perfectly fine. We just increment with the lagged terms. 


  Real t0, t1;
  Real total_time = 0.0;
  Real setup_time  = 0.0;
  Real advect_time = 0.0;
  Real diffusive_time = 0.0;

  // We begin with phi[0] = phi(t_n). Then update phi[m+1].
  for(int m = 0; m < m_p; m++){

    // Always update boundary conditions on the way in. All of these calls use the stuff that reside in the solvers,
    // which is what we need to do at the start of the time step. In principle, these things do not change
    // and so we could probably store them somewhere for increased performance. 
    if(m == 0 && !a_lagged_terms){ // This updates the boundary conditions; since these are used to compute the slopes,
      t0 = Timer::wallClock();         // and the m=0 slopes do not change, we only need to do this for the predictor. 
      sisdc::compute_E_into_scratch();
      sisdc::computeCdrEbStates();
      sisdc::compute_cdr_fluxes(a_time);
      sisdc::computeCdrDomainStates();
      sisdc::computeCdrDomainFluxes(a_time);
      sisdc::computeSigmaFlux();
      t1 = Timer::wallClock();

      total_time = -t0;
      setup_time = t1-t0;
    }
    
    // This does the transient rte advance. Stationary solves are done after computing the E-field.
    // The transient solve needs to happen BEFORE the reaction solve in case there is a tight coupling
    // between the RTE and CDR equations (relaxations of excited states). The integrate_rte routine compues
    // source terms just before advancing, and since source terms for the CDR equations are not updated
    // between this routine and the integrateAdvectionReaction call, we are ensured that the source terms
    // are consistent.
    t0 = Timer::wallClock();
    if(m_consistent_rte) {
      sisdc::integrate_rte(a_dt, m, a_lagged_terms);
    }

    // This computes phi_(m+1) = phi_m + dtm*FAR_m(phi_m) + lagged quadrature and lagged advection-reaction
    t0 = Timer::wallClock();
    sisdc::integrateAdvectionReaction(a_dt, m, a_lagged_terms);
    t1 = Timer::wallClock();
    advect_time += t1-t0;

    // This does the diffusion advance. It also adds in the remaining lagged diffusion terms before the implicit diffusion solve
    t0 = Timer::wallClock();
    sisdc::integrateDiffusion(a_dt, m, a_lagged_terms);
    t1 = Timer::wallClock();
    diffusive_time += t1-t0;

    // After the diffusion step we should update source terms and boundary conditions for the next step. We don't
    // do this on the last step. This is done either in the reconcileIntegrands routine, or after SISDC is done
    // with its substebs. 
    const bool last = (m == m_p-1);
    if(!last){
      Vector<EBAMRCellData*> cdr_densities_mp1 = sisdc::get_cdr_phik(m+1);
      EBAMRIVData& sigma_mp1 = sisdc::get_sigmak(m+1);
      const Real t_mp1 = m_tm[m+1];

      // Update electric field, RTE equations, source terms, and velocities. 
      if(m_consistent_E)   sisdc::updateField(cdr_densities_mp1, sigma_mp1);
      if(m_consistent_rte) sisdc::update_stationary_rte(cdr_densities_mp1, t_mp1);
      if(m_compute_S)      sisdc::computeCdrGradients(cdr_densities_mp1);
      if(m_compute_v)      sisdc::computeCdrVelo(cdr_densities_mp1, t_mp1);
      if(m_compute_S)      sisdc::compute_cdr_sources(cdr_densities_mp1, t_mp1);
      if(m_compute_D)      sisdc::updateDiffusionCoefficients();

      // Update boundary conditions for cdr and sigma equations. 
      sisdc::computeCdrEbStates(cdr_densities_mp1);
      sisdc::compute_cdr_fluxes(cdr_densities_mp1, t_mp1);
      sisdc::computeCdrDomainStates(cdr_densities_mp1);
      sisdc::computeCdrDomainFluxes(cdr_densities_mp1, t_mp1);
      sisdc::computeSigmaFlux();
    }
  }
  t1 = Timer::wallClock();

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

void sisdc::integrate_AP(const Real a_dt, const Real a_time, const bool a_lagged_terms){
  CH_TIME("sisdc::integrate_AP");
  if(m_verbosity > 5){
    pout() << "sisdc::integrate_AP" << endl;
  }


  // TLDR: This is the asymptotic preserving method for the advance routine. It can be optimized
  //  MayDay::Warning("sisdc::integrate_AP - routine is not yet done");


  // 1. The first time we enter this routine, source terms and velocities were updated.
  // 2. For further calls, source terms and velocities have been overwritten, but since the explicit
  //    operator slopes do not change, this is perfectly fine. We just increment with the lagged terms.

  // Not supported with subcycling
  if(m_subcycle == true) MayDay::Abort("sisdc::integrate_AP - AP method not supported with subcycling in time (yet)");


  Real t0, t1;
  Real total_time = 0.0;
  Real setup_time  = 0.0;
  Real advect_time = 0.0;
  Real diffusive_time = 0.0;

  // We begin with phi[0] = phi(t_n). Then update phi[m+1].
  for(int m = 0; m < m_p; m++){

    // Always update boundary conditions on the way in. All of these calls use the stuff that reside in the solvers,
    // which is what we need to do at the start of the time step. In principle, these things do not change
    // and so we could probably store them somewhere for increased performance. 
    if(m == 0 && !a_lagged_terms){ // This updates the boundary conditions; since these are used to compute the slopes,
      t0 = Timer::wallClock();         // and the m=0 slopes do not change, we only need to do this for the first SDC sweep
      sisdc::compute_E_into_scratch();
      sisdc::computeCdrEbStates();
      sisdc::compute_cdr_fluxes(a_time);
      sisdc::computeCdrDomainStates();
      sisdc::computeCdrDomainFluxes(a_time);
      sisdc::computeSigmaFlux();
      t1 = Timer::wallClock();

      total_time = -t0;
      setup_time = t1-t0;
    }
    
    // This does the transient rte advance. Stationary solves are done after computing the E-field.
    // The transient solve needs to happen BEFORE the reaction solve in case there is a tight coupling
    // between the RTE and CDR equations (relaxations of excited states). The integrate_rte routine compues
    // source terms just before advancing, and since source terms for the CDR equations are not updated
    // between this routine and the integrateAdvectionReaction call, we are ensured that the source terms
    // are consistent.
    t0 = Timer::wallClock();
    if(m_consistent_rte) {
      sisdc::integrate_rte(a_dt, m, a_lagged_terms);
    }

    // Predictor + corrector. Can definitely optimize a bunch of stuff here. 
    sisdc::integrate_AP_advection_reaction(a_dt, m, a_lagged_terms, true);
    sisdc::integrateDiffusion(a_dt, m, a_lagged_terms);

    updateField(get_cdr_phik(m+1), get_sigmak(m+1));
    compute_cdr_sources(get_cdr_phik(m), m_tm[m+1]);

    sisdc::integrate_AP_advection_reaction(a_dt, m, a_lagged_terms, false);
    sisdc::integrateDiffusion(a_dt, m, a_lagged_terms);

    // After the diffusion step we should update source terms and boundary conditions for the next step. We don't
    // do this on the last step. This is done either in the reconcileIntegrands routine, or after SISDC is done
    // with its substebs. 
    const bool last = (m == m_p-1);
    if(!last){
      Vector<EBAMRCellData*> cdr_densities_mp1 = sisdc::get_cdr_phik(m+1);
      EBAMRIVData& sigma_mp1 = sisdc::get_sigmak(m+1);
      const Real t_mp1 = m_tm[m+1];

      // Update electric field, RTE equations, source terms, and velocities. 
      if(m_consistent_E)   sisdc::updateField(cdr_densities_mp1, sigma_mp1);
      if(m_consistent_rte) sisdc::update_stationary_rte(cdr_densities_mp1, t_mp1);
      if(m_compute_S)      sisdc::computeCdrGradients(cdr_densities_mp1);
      if(m_compute_v)      sisdc::computeCdrVelo(cdr_densities_mp1, t_mp1);
      if(m_compute_S)      sisdc::compute_cdr_sources(cdr_densities_mp1, t_mp1);
      if(m_compute_D)      sisdc::updateDiffusionCoefficients();

      // Update boundary conditions for cdr and sigma equations. 
      sisdc::computeCdrEbStates(cdr_densities_mp1);
      sisdc::compute_cdr_fluxes(cdr_densities_mp1, t_mp1);
      sisdc::computeCdrDomainStates(cdr_densities_mp1);
      sisdc::computeCdrDomainFluxes(cdr_densities_mp1, t_mp1);
      sisdc::computeSigmaFlux();
    }
  }
  t1 = Timer::wallClock();

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

void sisdc::integrateAdvectionReaction(const Real a_dt, const int a_m, const bool a_lagged_terms){
  CH_TIME("sisdc::integrateAdvectionReaction");
  if(m_verbosity > 5){
    pout() << "sisdc::integrateAdvectionReaction" << endl;
  }

  // Advance phi_(m+1) = phi_m + dtm*F_A using either subcyling or not. These routines do nothing
  // with the operator slopes for phi, but they do adjust the slopes m_Fsig (but not m_Fsum) for sigma. Incidentally,
  // if m=0 and a_lagged_terms=true, we can increment directly with the precomputed advection-reaction. This means that
  // we can skip the advective advance. The sigma advance is accordingly also skipped.
  const bool skip = (a_m == 0 && a_lagged_terms);
  const Real t0 = Timer::wallClock();
  if(!skip){
    if(m_subcycle){
      if(m_multistep){
	sisdc::integrateAdvection_multistep(a_dt, a_m, a_lagged_terms);
      }
      else{
	sisdc::integrateAdvection_subcycle(a_dt, a_m, a_lagged_terms);
      }
    }
    else{
      sisdc::integrateAdvection_nosubcycle(a_dt, a_m, a_lagged_terms);
    }
  }
  const Real t1 = Timer::wallClock();

  // Add in the reaction term and then compute the new operator slopes.
  // If this is the corrector and m=0, we skipped the advection advance because we can use the precomputed
  // advection-reaction operator slope. In this case phi_(m+1) is bogus and we need to recompute it. Otherwise,
  // phi_(m+1) = phi_m + dtm*FA_m, and we just increment with the reaction operator. 
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = get_CdrStorage(solver_it);

    // phi_(m+1) = phi_M
    EBAMRCellData& phi_m1      = storage->getPhi()[a_m+1];
    EBAMRCellData& scratch     = storage->getScratch();
    const EBAMRCellData& phi_m = storage->getPhi()[a_m];
    
    // Increment with operator slopes. m=0 and corrector is a special case where we skipped the advective advance,
    // choosing instead to use the old slopes (which did not change)
    if(skip){ // Can use the old slopes
      const EBAMRCellData& FAR_m = storage->getFAR()[a_m]; // Slope, doesn't require recomputation. 
      DataOps::copy(phi_m1, phi_m);
      DataOps::incr(phi_m1, FAR_m, m_dtm[a_m]);
      if(a_lagged_terms) {
	DataOps::copy(scratch, FAR_m);
      }
    }
    else{ // If we made it here, phi_(m+1) = phi_m + dtm*FA(phi_m) through the integrateAdvection_subcycle routine
      EBAMRCellData& FAR_m     = storage->getFAR()[a_m]; // Currently the old slope
      EBAMRCellData& src = solver->getSource();    // Updated source

      // Increment swith source and then compute slope. This has already been done 
      if(!(m_cycle_sources && m_subcycle)){
	DataOps::incr(phi_m1, src, m_dtm[a_m]);  // phi_(m+1) = phi_m + dtm*(FA_m + FR_m)
      }

      // This shouldn't be necessary
      m_amr->averageDown(phi_m1, m_cdr->getPhase());
      m_amr->interpGhost(phi_m1, m_cdr->getPhase());

      if(a_lagged_terms){ // Back up the old slope first, we will need it for the lagged term
	DataOps::copy(scratch, FAR_m);
      }

      // Re-compute the advection-reaction slope for node t_m
      DataOps::copy(FAR_m, phi_m1);            // FAR_m = (phi_(m+1) - phi_m)/dtm
      DataOps::incr(FAR_m, phi_m, -1.0);       // :
      DataOps::scale(FAR_m, 1./m_dtm[a_m]);    // :

      // Shouldn't be necessary
      m_amr->averageDown(FAR_m, m_cdr->getPhase());
      m_amr->interpGhost(FAR_m, m_cdr->getPhase());
    }

    // Now add in the lagged advection-reaction and quadrature terms. This is a bit weird, but we did overwrite
    // FAR_m above after the advection-reaction advance, but we also backed up the old term into scratch. 
    if(a_lagged_terms){
      DataOps::incr(phi_m1, scratch, -m_dtm[a_m]); // phi_(m+1)^(k+1) = phi_m^(k+1) + dtm*(FAR_m^(k+1) - FAR_m^k)
      sisdc::quad(scratch, storage->getF(), a_m);  // Does the quadrature of the lagged operator slopes. 
      DataOps::incr(phi_m1, scratch, 0.5*a_dt);    // phi_(m+1)^(k+1) = phi_m^(k+1) + dtm*(FAR_m^(k+1) - FAR_m^k) + I_m^(m+1)
    }
  }
  const Real t2 = Timer::wallClock();

  // Add in the lagged terms for sigma. As above, m=0 and corrector is a special case where we just use the old slopes.
  EBAMRIVData& sigma_m1      = m_sigma_scratch->get_sigma()[a_m+1];
  const EBAMRIVData& sigma_m = m_sigma_scratch->get_sigma()[a_m];
  if(skip){
    const EBAMRIVData& Fsig_m = m_sigma_scratch->getFold()[a_m]; // Here, we should be able to use either Fold or Fnew
    DataOps::copy(sigma_m1, sigma_m);                            // since Fsig_0 is only computed once. 
    DataOps::incr(sigma_m1, Fsig_m, m_dtm[a_m]);
  }

  const Real t3 = Timer::wallClock();
  if(a_lagged_terms){ // Add in the lagged terms. When we make it here, sigma_(m+1) = sigma_m + dtm*Fsig_m. 
    EBAMRIVData& Fsig_lag = m_sigma_scratch->getFold()[a_m];
    DataOps::incr(sigma_m1, Fsig_lag, -m_dtm[a_m]);

    // Add in the quadrature term
    EBAMRIVData& scratch = m_sigma_scratch->getScratch();
    sisdc::quad(scratch, m_sigma_scratch->getFold(), a_m);
    DataOps::incr(sigma_m1, scratch, 0.5*a_dt); // Mult by 0.5*a_dt due to scaling on [-1,1] for quadrature
  }
  const Real t4 = Timer::wallClock();

#if 0
  pout() << "integrateAdvectionReaction::" << endl;
  pout() << "t1-t0 = " << t1-t0 << endl;
  pout() << "t2-t1 = " << t2-t2 << endl;
  pout() << "t3-t2 = " << t3-t2 << endl;
  pout() << "t4-t3 = " << t4-t3 << endl;
#endif
}

void sisdc::integrateAdvection_nosubcycle(const Real a_dt, const int a_m, const bool a_lagged_terms){
  CH_TIME("sisdc::integrateAdvection_nosubcycle");
  if(m_verbosity > 5){
    pout() << "sisdc::integrateAdvection_nosubcycle" << endl;
  }

  // TLDR; This routine should do phi_(m+1) = phi_m + dtm*FA_m, and sigma_(m+1) = sigma_m + dt*Fsig_m.
  //       It also computes the sigma slope.
  //
  //       The lagged terms are not a part of this routine. 

  if(a_m == 0 && a_lagged_terms){
    MayDay::Abort("sisdc::integrateAdvection_nosubcycle - (m==0 && corrector==true) should never happen");
  }

  // Advance cdr equations
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = get_CdrStorage(solver_it);

    EBAMRCellData& phi_m1      = storage->getPhi()[a_m+1];
    EBAMRCellData& scratch     = storage->getScratch();
    const EBAMRCellData& phi_m = storage->getPhi()[a_m];

    if(solver->isMobile()){
      const Real extrap_dt = m_extrap_advect ? 2.0*m_extrap_dt*m_dtm[a_m] : 0.0; // Factor of 2 due to EBPatchAdvect
      solver->computeDivF(scratch, phi_m, extrap_dt, true);                     // scratch =  Div(v_m*phi_m^(k+1))
    
      DataOps::copy(phi_m1, phi_m);
      DataOps::incr(phi_m1, scratch, -m_dtm[a_m]);
      //      DataOps::floor(phi_m1, 0.0);
      m_amr->averageDown(phi_m1, m_cdr->getPhase());
      m_amr->interpGhost(phi_m1, m_cdr->getPhase());

    }
    else{
      DataOps::copy(phi_m1, phi_m);
    }
  }

  // Update sigma. Also compute the new slope.
  EBAMRIVData& sigma_m1      = m_sigma_scratch->get_sigma()[a_m+1];
  EBAMRIVData& Fsig_new      = m_sigma_scratch->getFnew()[a_m];
  const EBAMRIVData& sigma_m = m_sigma_scratch->get_sigma()[a_m];
  m_sigma->computeRHS(Fsig_new); // Fills Fsig_new with BCs from CDR solvers
  DataOps::copy(sigma_m1, sigma_m);
  DataOps::incr(sigma_m1, Fsig_new, m_dtm[a_m]);
}

void sisdc::integrateAdvection_multistep(const Real a_dt, const int a_m, const bool a_lagged_terms){
  CH_TIME("sisdc::integrateAdvection_nosubcycle");
  if(m_verbosity > 5){
    pout() << "sisdc::integrateAdvection_nosubcycle" << endl;
  }

  // TLDR; This routine should do phi_(m+1) = phi_m + dtm*FA_m, and sigma_(m+1) = sigma_m + dt*Fsig_m.
  //       It also computes the sigma slope.
  //
  //       The lagged terms are not a part of this routine. 

  if(a_m == 0 && a_lagged_terms){
    MayDay::Abort("sisdc::integrateAdvection_multistep - (m==0 && corrector==true) should never happen");
  }

  // Copy into phi[a_m+1]
  sisdc::subcycle_copy_states(a_m);

  const int nsteps = ceil(m_dtm[a_m]/(m_cycleCFL*m_dt_cfl));
  const Real dt    = m_dtm[a_m]/nsteps;

#if 0
  MayDay::Abort("sisdc::integrateAdvection_multistep - multistep is not done");
#endif
  
  // Advance cdr equations. Use a Heun's method
  for (int istep = 0; istep < nsteps; istep++){

    // Always need to update boundary conditions and source terms (if we substep with the sources)
    if (istep > 0){
      Vector<EBAMRCellData*> cdr_densities_mp1 = sisdc::get_cdr_phik(a_m+1);
      EBAMRIVData& sigma_mp1 = sisdc::get_sigmak(a_m+1);
      const Real t_mp1 = m_tm[a_m+1];

      if(m_cycle_sources){
	if(m_compute_S)      sisdc::computeCdrGradients(cdr_densities_mp1);
	if(m_compute_S)      sisdc::compute_cdr_sources(cdr_densities_mp1, t_mp1);
      }
      
      sisdc::computeCdrEbStates(cdr_densities_mp1);
      sisdc::compute_cdr_fluxes(cdr_densities_mp1, t_mp1);
      sisdc::computeCdrDomainStates(cdr_densities_mp1);
      sisdc::computeCdrDomainFluxes(cdr_densities_mp1, t_mp1);
      sisdc::computeSigmaFlux();
    }

    // First Heun method stage
    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<CdrStorage>& storage = get_CdrStorage(solver_it);

      EBAMRCellData& phi_m1      = storage->getPhi()[a_m+1];
      EBAMRCellData& scratch     = storage->getScratch();
      EBAMRCellData& scratch2    = storage->getScratch2();
      const EBAMRCellData& src   = solver->getSource();

      const Real extrap_dt = m_extrap_advect ? 2.0*m_extrap_dt*dt : 0.0; // Factor of 2 due to EBPatchAdvect
      solver->computeDivF(scratch, phi_m1, extrap_dt, true);            // scratch =  Div(v_m*phi_m^(k+1))

      DataOps::copy(scratch2, phi_m1);
      DataOps::incr(phi_m1, scratch, -dt);
      if(m_cycle_sources){
	DataOps::incr(phi_m1, src, dt);
      }
      m_amr->averageDown(phi_m1, m_cdr->getPhase());
      m_amr->interpGhost(phi_m1, m_cdr->getPhase());
      DataOps::floor(phi_m1, 0.0);
    }

    // Update sigma. Also compute the new slope.
#if 0 // Bogus code, please, please fix this!
    EBAMRIVData& sigma_m1      = m_sigma_scratch->get_sigma()[a_m+1];
    EBAMRIVData& Fsig_new      = m_sigma_scratch->getFnew()[a_m];
    const EBAMRIVData& sigma_m = m_sigma_scratch->get_sigma()[a_m];
    m_sigma->computeRHS(Fsig_new); // Fills Fsig_new with BCs from CDR solvers
    DataOps::copy(sigma_m1, sigma_m);
    DataOps::incr(sigma_m1, Fsig_new, m_dtm[a_m]);
#endif

    // Always update BC between stages
    Vector<EBAMRCellData*> cdr_densities_mp1 = sisdc::get_cdr_phik(a_m+1);
    EBAMRIVData& sigma_mp1 = sisdc::get_sigmak(a_m+1);
    const Real t_mp1 = m_tm[a_m+1];

    if(m_cycle_sources){
      if(m_compute_S)      sisdc::computeCdrGradients(cdr_densities_mp1);
      if(m_compute_S)      sisdc::compute_cdr_sources(cdr_densities_mp1, t_mp1);
    }
      
    sisdc::computeCdrEbStates(cdr_densities_mp1);
    sisdc::compute_cdr_fluxes(cdr_densities_mp1, t_mp1);
    sisdc::computeCdrDomainStates(cdr_densities_mp1);
    sisdc::computeCdrDomainFluxes(cdr_densities_mp1, t_mp1);
    sisdc::computeSigmaFlux();

    // Second Heun method stage
    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<CdrStorage>& storage = get_CdrStorage(solver_it);

      EBAMRCellData& phi_m1      = storage->getPhi()[a_m+1];
      EBAMRCellData& scratch     = storage->getScratch();
      EBAMRCellData& scratch2    = storage->getScratch2();
      const EBAMRCellData& src   = solver->getSource();

      const Real extrap_dt = m_extrap_advect ? 2.0*m_extrap_dt*dt : 0.0; // Factor of 2 due to EBPatchAdvect
      solver->computeDivF(scratch, phi_m1, extrap_dt, true);            // scratch =  Div(v_m*phi_m^(k+1))
    
      DataOps::incr(phi_m1, scratch, -dt);
      if(m_cycle_sources){
	DataOps::incr(phi_m1, src, dt);
      }
      DataOps::incr(phi_m1, scratch2, 1.0);
      DataOps::scale(phi_m1, 0.5);
      m_amr->averageDown(phi_m1, m_cdr->getPhase());
      m_amr->interpGhost(phi_m1, m_cdr->getPhase());
      DataOps::floor(phi_m1, 0.0);
    }
  }


}

void sisdc::integrateDiffusion(const Real a_dt, const int a_m, const bool a_lagged_terms){
  CH_TIME("sisdc::integrateDiffusion");
  if(m_verbosity > 5){
    pout() << "sisdc::integrateDiffusion" << endl;
  }

  // TLDR: We're solving
  //
  // phi_(m+1)^(k+1) = phi_(m)^(k+1,\ast) + dtm*FD_(m+1)^(k+1) + sources. 
  //
  // This routine does not modify FD_(m+1)^k. This is replaced by FD_(m+1)^(k+1) later on. 
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = sisdc::get_CdrStorage(solver_it);
    
    if(solver->isDiffusive()){
      EBAMRCellData& phi_m1      = storage->getPhi()[a_m+1]; // Advected solution. Possibly with lagged terms. 
      const EBAMRCellData& phi_m = storage->getPhi()[a_m];

      // Build the diffusion source term
      EBAMRCellData& source   = storage->getScratch();
      EBAMRCellData& init_soln = storage->getScratch2();
      DataOps::setValue(source, 0.0); // No source term
      
      DataOps::copy(init_soln, phi_m1);      // Copy initial solutions
      if(a_lagged_terms){
	const EBAMRCellData& FD_m1k = storage->getFD()[a_m+1];      // FD_(m+1)^k. Lagged term.
	DataOps::incr(init_soln, FD_m1k, -m_dtm[a_m]);
      }
      m_amr->averageDown(init_soln, m_cdr->getPhase());
      m_amr->interpGhost(init_soln, m_cdr->getPhase());
#if 1 // Original code
      DataOps::copy(phi_m1, phi_m);
#else // Debug code
      DataOps::copy(phi_m1, init_soln);
#endif

      // Solve
      if(m_useTGA){
	solver->advanceTGA(phi_m1, init_soln, source, m_dtm[a_m]); // No source. 
      }
      else{
	solver->advanceEuler(phi_m1, init_soln, source, m_dtm[a_m]); // No source. 
      }
      m_amr->averageDown(phi_m1, m_cdr->getPhase());
      m_amr->interpGhost(phi_m1, m_cdr->getPhase());
      DataOps::floor(phi_m1, 0.0);

      // Update the operator slope
      EBAMRCellData& FD_m1k = storage->getFD()[a_m+1];
      DataOps::setValue(FD_m1k, 0.0);
      DataOps::incr(FD_m1k, phi_m1, 1.0);
      DataOps::incr(FD_m1k, init_soln, -1.0);
      DataOps::scale(FD_m1k, 1./m_dtm[a_m]);

      m_amr->averageDown(FD_m1k, m_cdr->getPhase());
      m_amr->interpGhost(FD_m1k, m_cdr->getPhase());
    }
    else{
      EBAMRCellData& FD_m1k = storage->getFD()[a_m+1];
      DataOps::setValue(FD_m1k, 0.0);
    }
  }
}

void sisdc::integrate_AP_advection_reaction(const Real a_dt, const int a_m, const bool a_lagged_terms, const bool a_predictor){
  CH_TIME("sisdc::integrate_AP_advection_reaction");
  if(m_verbosity > 5){
    pout() << "sisdc::integrate_AP_advection_reaction" << endl;
  }

  // Advance phi_(m+1) = phi_m + dtm*F_A using either subcyling or not. These routines do nothing
  // with the operator slopes for phi, but they do adjust the slopes m_Fsig (but not m_Fsum) for sigma. Incidentally,
  // if m=0 and a_lagged_terms=true, we can increment directly with the precomputed advection-reaction. This means that
  // we can skip the advective advance. The sigma advance is accordingly also skipped.
  const bool skip = (a_m == 0 && a_lagged_terms); // First step in corrector => true. False otherwise. 

  // Compute phi_(m+1) = phi_m + dt_m*DivF. If this is the AP predictor, we can store divF. If this is the AP corrector,
  // we can fetch divF from the predictor. 
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = get_CdrStorage(solver_it);

    EBAMRCellData& phi_m1      = storage->getPhi()[a_m+1]; 
    const EBAMRCellData& phi_m = storage->getPhi()[a_m];


    if(solver->isMobile()){
      EBAMRCellData& divF = storage->getDivF();
      if(a_predictor && !skip){ // Need to compute divF for predictor. For the corrector, it has already been computed. Yay!
	const Real extrap_dt = m_extrap_advect ? 2.0*m_extrap_dt*m_dtm[a_m] : 0.0; // Factor of 2 due to EBPatchAdvect
	solver->computeDivF(divF, phi_m, extrap_dt, true);                        // divF =  Div(v_m*phi_m^(k+1))
      }

      // phi_(m+1) = phi_m - dt*div(F)
      DataOps::copy(phi_m1, phi_m);
      DataOps::incr(phi_m1, divF, -m_dtm[a_m]);
      DataOps::floor(phi_m1, 0.0);
    }
    else{
      DataOps::copy(phi_m1, phi_m);
    }
  }

  // Update sigma. Also compute the new slope. This codes computes the new slope and does
  // sigma_(m+1) = sigma_m + dt*J_m(phi_m^(k+1)). If this is the first step in the SDC correcting sweeps, we
  // can skip this update. Since the boundary flux comes in through div(F), we only need to do this in the AP predictor
  if(a_predictor){
    EBAMRIVData& sigma_m1      = m_sigma_scratch->get_sigma()[a_m+1];
    EBAMRIVData& Fsig_new      = m_sigma_scratch->getFnew()[a_m];
    const EBAMRIVData& sigma_m = m_sigma_scratch->get_sigma()[a_m];
    m_sigma->computeRHS(Fsig_new); // Fills Fsig_new with BC data in solver
    DataOps::copy(sigma_m1, sigma_m);
    DataOps::incr(sigma_m1, Fsig_new, m_dtm[a_m]);
  }

  // Add in the reaction operator
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = get_CdrStorage(solver_it);

    EBAMRCellData& phi_m1      = storage->getPhi()[a_m+1]; 
    const EBAMRCellData& src   = solver->getSource();

    DataOps::incr(phi_m1, src, m_dtm[a_m]);
  }

  // Add in the lagged advection-reaction terms for the CDR equations. 
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = get_CdrStorage(solver_it);

    EBAMRCellData& phi_m1  = storage->getPhi()[a_m+1]; // This contains phi_m - dt*div(F) + dt*R
    EBAMRCellData& phi_m   = storage->getPhi()[a_m];   // This contains phi_m - dt*div(F) + dt*R
    EBAMRCellData& FAR_m   = storage->getFAR()[a_m];   // Old operator slope
    EBAMRCellData& scratch = storage->getScratch();    // Scratch storage


    // Back up the old slope first
    if(a_predictor){ // For the predictor, we shouldn't monkey with the slopes
      if(a_lagged_terms){
	DataOps::incr(phi_m1, FAR_m, -m_dtm[a_m]);   // phi_(m+1)^(k+1) = phi_m^(k+1) + dtm*(-div(F) + R - FAR_m^k)
	sisdc::quad(scratch, storage->getF(), a_m);  // Does the quadrature of the lagged operator slopes. 
	DataOps::incr(phi_m1, scratch, 0.5*a_dt);    // phi_(m+1)^(k+1) = phi_m^(k+1) + dtm*(FAR_m^(k+1) - FAR_m^k) + I_m^(m+1)
      }
    }
    else{ // For the corrector, we do need to update the slopes as we move along
      if(a_lagged_terms) { 
	DataOps::copy(scratch, FAR_m); // Put the old operator slope in scratch; we need a backup
      }

      // Update the slopes
      DataOps::copy(FAR_m, phi_m1);            // FAR_m = (phi_(m+1) - phi_m)/dtm
      DataOps::incr(FAR_m, phi_m, -1.0);       // :
      DataOps::scale(FAR_m, 1./m_dtm[a_m]);    // :

      //
      if(a_lagged_terms){
	DataOps::incr(phi_m1, scratch, -m_dtm[a_m]); // phi_(m+1)^(k+1) = phi_m^(k+1) + dtm*(FAR_m^(k+1) - FAR_m^k)
	sisdc::quad(scratch, storage->getF(), a_m);  // Does the quadrature of the lagged operator slopes. 
	DataOps::incr(phi_m1, scratch, 0.5*a_dt);    // phi_(m+1)^(k+1) = phi_m^(k+1) + dtm*(FAR_m^(k+1) - FAR_m^k) + I_m^(m+1)
      }
    }
  }

  // Add in the lagged terms for sigma. We have already done sigma_(m+1) = sigma_m + dt*J_m. Now add the lagged terms
  if(a_predictor){
    EBAMRIVData& sigma_m1      = m_sigma_scratch->get_sigma()[a_m+1];
    const EBAMRIVData& sigma_m = m_sigma_scratch->get_sigma()[a_m];

    if(a_lagged_terms){
      EBAMRIVData& Fsig_lag = m_sigma_scratch->getFold()[a_m]; // Add in the lagged term
      DataOps::incr(sigma_m1, Fsig_lag, -m_dtm[a_m]);

      // Add in the quadrature term
      EBAMRIVData& scratch = m_sigma_scratch->getScratch();
      sisdc::quad(scratch, m_sigma_scratch->getFold(), a_m);
      DataOps::incr(sigma_m1, scratch, 0.5*a_dt); // Mult by 0.5*a_dt due to scaling on [-1,1] for quadrature
    }
  }
}

void sisdc::reconcileIntegrands(){
  CH_TIME("sisdc::reconcileIntegrands");
  if(m_verbosity > 5){
    pout() << "sisdc::reconcileIntegrands" << endl;
  }

  Vector<EBAMRCellData*> cdr_densities_p = sisdc::get_cdr_phik(m_p);
  EBAMRIVData& sigma_p = sisdc::get_sigmak(m_p);
  const Real t_p = m_tm[m_p];

  //  Update electric field, RTE equations, source terms, and velocities
  if(m_consistent_E)   sisdc::updateField(cdr_densities_p, sigma_p);
  if(m_consistent_rte) sisdc::update_stationary_rte(cdr_densities_p, t_p);
  if(m_compute_S)      sisdc::computeCdrGradients(cdr_densities_p);
  if(m_compute_v)      sisdc::computeCdrVelo(cdr_densities_p, t_p);
  if(m_compute_S)      sisdc::compute_cdr_sources(cdr_densities_p, t_p);

  // Update boundary conditions for cdr and sigma equations
  sisdc::computeCdrEbStates(cdr_densities_p);
  sisdc::compute_cdr_fluxes(cdr_densities_p, t_p);
  sisdc::computeCdrDomainStates(cdr_densities_p);
  sisdc::computeCdrDomainFluxes(cdr_densities_p, t_p);
  sisdc::computeSigmaFlux();

  // Now compute FAR_p - that wasn't done when we integrated
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = sisdc::get_CdrStorage(solver_it);
    const int idx = solver_it.get_solver();

    // This has not been computed yet. Do it.
    EBAMRCellData& FAR_p       = storage->getFAR()[m_p];
    const EBAMRCellData& phi_p = *cdr_densities_p[idx] ;
    const EBAMRCellData& src   = solver->getSource();

    if(solver->isMobile()){
      const Real extrap_dt = m_extrap_advect ? 2.0*m_extrap_dt*m_dtm[m_p-1] : 0.0; // Factor of 2 because of EBPatchAdvect
      solver->computeDivF(FAR_p, phi_p, extrap_dt, true); // FAR_p =  Div(v_p*phi_p)
      DataOps::scale(FAR_p, -1.0);                        // FAR_p = -Div(v_p*phi_p)
    }
    else{
      DataOps::setValue(FAR_p, 0.0);
    }
    DataOps::incr(FAR_p, src, 1.0);                     // RHS = -Div(v_m*phi_m) + S_m = FAR(phi_m)

    // Build the integrand
    for (int m = 0; m <= m_p; m++){
      EBAMRCellData& F_m   = storage->getF()[m];
      EBAMRCellData& FD_m  = storage->getFD()[m];
      EBAMRCellData& FAR_m = storage->getFAR()[m];

      DataOps::copy(F_m, FAR_m);
      if(solver->isDiffusive()){
	DataOps::incr(F_m, FD_m, 1.0);
      }

      // Shouldn't be necessary
      m_amr->averageDown(F_m, m_cdr->getPhase());
      m_amr->interpGhost(F_m, m_cdr->getPhase());
    }
  }

  // Compute Fsig_p - that wasn't done either
  EBAMRIVData& Fnew_p = m_sigma_scratch->getFnew()[m_p];
  m_sigma->computeRHS(Fnew_p);
  for (int m = 0; m <= m_p; m++){
    EBAMRIVData& Fold_m = m_sigma_scratch->getFold()[m];
    EBAMRIVData& Fnew_m = m_sigma_scratch->getFnew()[m];
    DataOps::copy(Fold_m, Fnew_m);
  }
}

void sisdc::initializeErrors(){
  CH_TIME("sisdc::corrector_initializeErrors");
  if(m_verbosity > 5){
    pout() << "sisdc::corrector_initializeErrors" << endl;
  }

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = sisdc::get_CdrStorage(solver_it);
    const int idx = solver_it.get_solver();
    
    // These should be zero
    if(idx == m_error_idx || m_error_idx < 0){
      EBAMRCellData& error = storage->getError();
      const EBAMRCellData& phi_final = storage->getPhi()[m_p];

      DataOps::setValue(error, 0.0);
      DataOps::incr(error, phi_final, -1.0);
    }
  }

  EBAMRIVData& error = m_sigma_scratch->getError();
  const EBAMRIVData& sigma_final = m_sigma_scratch->get_sigma()[m_p];
  DataOps::setValue(error, 0.0);
  DataOps::incr(error, sigma_final, -1.0);
}

void sisdc::finalizeErrors(){
  CH_TIME("sisdc::corrector_finalizeErrors");
  if(m_verbosity > 5){
    pout() << "sisdc::corrector_finalizeErrors" << endl;
  }

  const Real safety = 1.E-20;

  m_max_error = 0.0;
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = sisdc::get_CdrStorage(solver_it);
    const int idx = solver_it.get_solver();

    // Compute error
    if(idx == m_error_idx || m_error_idx < 0){
      EBAMRCellData& error       = storage->getError();
      const EBAMRCellData& phi_p = storage->getPhi()[m_p];
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
  EBAMRIVData& error = m_sigma_scratch->getError();
  const EBAMRIVData& sigma_final = m_sigma_scratch->get_sigma()[m_p];
  DataOps::incr(error, sigma_final, 1.0);
  m_sigma_error = 0.0; // I don't think this is ever used...


}

void sisdc::computeNewDt(bool& a_accept_step, const Real a_dt, const int a_num_corrections){
  CH_TIME("sisdc::computeNewDt");
  if(m_verbosity > 5){
    pout() << "sisdc::computeNewDt" << endl;
  }

  // If a_dt was the smallest possible CFL or hardcap time step, we just have to accept it
  const Real max_gl_dist = sisdc::getMaxNodeDistance();
  Real dt_cfl = 2.0*m_dt_cfl/max_gl_dist; // This is the smallest time step ON THE FINEST LEVEL

  int Nref = 1;
  if(m_subcycle){
    for (int lvl = 0; lvl < m_amr->getFinestLevel(); lvl++){
      Nref = Nref*m_amr->getRefinementRatios()[lvl];
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

void sisdc::adaptiveReport(const Real a_first_dt, const Real a_dt, const Real a_new_dt, const int a_corr, const int a_rej, const Real a_max_err){
  CH_TIME("sisdc::adaptiveReport");
  if(m_verbosity > 5){
    pout() << "sisdc::adaptiveReport" << endl;
  }

  pout() << "\n";
  pout() << "sisdc::adaptiveReport breakdown" << endl;
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

void sisdc::computeDt(Real& a_dt, TimeCode::which_code& a_timeCode){
  CH_TIME("sisdc::computeDt");
  if(m_verbosity > 5){
    pout() << "sisdc::computeDt" << endl;
  }

  Real dt = 1.E99;

  int Nref = 1;
  for (int lvl = 0; lvl < m_amr->getFinestLevel(); lvl++){
    Nref = Nref*m_amr->getRefinementRatios()[lvl];
  }
  const Real max_gl_dist = sisdc::getMaxNodeDistance();
  m_dt_cfl = m_cdr->compute_cfl_dt();

  //  Real dt_cfldt_cfl = (m_subcycle) ? m_dt_cfl*Nref : m_dt_cfl;
  Real dt_cfl = 2.0*m_dt_cfl/max_gl_dist;
  dt_cfl = m_subcycle ? Nref*dt_cfl : dt_cfl;
  
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
      a_timeCode = TimeCode::Error;
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

void sisdc::regridInternals(){
  CH_TIME("sisdc::regridInternals");
  if(m_verbosity > 5){
    pout() << "sisdc::regridInternals" << endl;
  }

  m_accum_cfl = 0.0;
  m_cdr_error.resize(m_plaskin->get_num_species());
  
  sisdc::allocateCdrStorage();
  sisdc::allocateFieldStorage();
  sisdc::allocateRtStorage();
  sisdc::allocateSigmaStorage();

  sisdc::setupQuadratureNodes(m_p);
  sisdc::setupQmj(m_p);
}

void sisdc::allocateCdrStorage(){
  const int ncomp       = 1;
  const int num_species = m_plaskin->get_num_species();

  m_cdr_scratch.resize(num_species);
  
  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx] = RefCountedPtr<CdrStorage> (new CdrStorage(m_amr, m_cdr->getPhase(), ncomp));
    m_cdr_scratch[idx]->allocate_storage(m_p);
  }
}

void sisdc::allocateFieldStorage(){
  const int ncomp = 1;
  m_fieldSolver_scratch = RefCountedPtr<FieldStorage> (new FieldStorage(m_amr, m_cdr->getPhase(), ncomp));
  m_fieldSolver_scratch->allocate_storage(m_p);
}

void sisdc::allocateRtStorage(){
  const int ncomp       = 1;
  const int num_Photons = m_plaskin->get_num_Photons();
  m_rte_scratch.resize(num_Photons);
  
  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx] = RefCountedPtr<RtStorage> (new RtStorage(m_amr, m_rte->getPhase(), ncomp));
    m_rte_scratch[idx]->allocate_storage(m_p);
  }
}

void sisdc::allocateSigmaStorage(){
  const int ncomp = 1;
  m_sigma_scratch = RefCountedPtr<SigmaStorage> (new SigmaStorage(m_amr, m_cdr->getPhase(), ncomp));
  m_sigma_scratch->allocate_storage(m_p);
}

void sisdc::deallocateInternals(){
  CH_TIME("sisdc::deallocateInternals");
  if(m_verbosity > 5){
    pout() << "sisdc::deallocateInternals" << endl;
  }

  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx]->deallocate_storage();
    m_cdr_scratch[idx] = RefCountedPtr<CdrStorage>(0);
  }

  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx]->deallocate_storage();
    m_rte_scratch[idx] = RefCountedPtr<RtStorage>(0);
  }

  m_cdr_scratch.resize(0);
  m_rte_scratch.resize(0);

  m_fieldSolver_scratch->deallocate_storage();
  m_fieldSolver_scratch = RefCountedPtr<FieldStorage>(0);
  
  m_sigma_scratch->deallocate_storage();
  m_sigma_scratch = RefCountedPtr<SigmaStorage>(0);
}

void sisdc::compute_E_into_scratch(){
  CH_TIME("sisdc::compute_E_into_scratch");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_E_into_scratch" << endl;
  }
  
  EBAMRCellData& E_cell = m_fieldSolver_scratch->getElectricFieldCell();
  EBAMRFluxData& E_face = m_fieldSolver_scratch->getElectricFieldFace();
  EBAMRIVData&   E_eb   = m_fieldSolver_scratch->getElectricFieldEb();
  EBAMRIFData&   E_dom  = m_fieldSolver_scratch->getElectricFieldDomain();

  const MFAMRCellData& phi = m_fieldSolver->getPotential();
  
  sisdc::compute_E(E_cell, m_cdr->getPhase(), phi);     // Compute cell-centered field
  sisdc::compute_E(E_face, m_cdr->getPhase(), E_cell);  // Compute face-centered field
  sisdc::compute_E(E_eb,   m_cdr->getPhase(), E_cell);  // EB-centered field

  TimeStepper::extrapolate_to_domain_faces(E_dom, m_cdr->getPhase(), E_cell);
}

void sisdc::computeCdrGradients(){
  CH_TIME("sisdc::computeCdrGradients");
  if(m_verbosity > 5){
    pout() << "sisdc::computeCdrGradients" << endl;
  }

  sisdc::computeCdrGradients(m_cdr->getPhis());
}

void sisdc::computeCdrGradients(const Vector<EBAMRCellData*>& a_phis){
  CH_TIME("sisdc::computeCdrGradients");
  if(m_verbosity > 5){
    pout() << "sisdc::computeCdrGradients" << endl;
  }

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    RefCountedPtr<CdrStorage>& storage = sisdc::get_CdrStorage(solver_it);
    EBAMRCellData& grad = storage->getGradient();
    m_amr->computeGradient(grad, *a_phis[idx]);
    //    m_amr->averageDown(grad, m_cdr->getPhase());
    m_amr->interpGhost(grad, m_cdr->getPhase());
  }
}

void sisdc::computeCdrVelo(const Real a_time){
  CH_TIME("sisdc::computeCdrVelo");
  if(m_verbosity > 5){
    pout() << "sisdc::computeCdrVelo" << endl;
  }

  sisdc::computeCdrVelo(m_cdr->getPhis(), a_time);
}

void sisdc::computeCdrVelo(const Vector<EBAMRCellData*>& a_phis, const Real a_time){
  CH_TIME("sisdc::computeCdrVelo(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "sisdc::computeCdrVelo(Vector<EBAMRCellData*>, Real)" << endl;
  }

  Vector<EBAMRCellData*> velocities = m_cdr->getVelocities();
  sisdc::computeCdrDriftVelocities(velocities, a_phis, m_fieldSolver_scratch->getElectricFieldCell(), a_time);
}

void sisdc::computeCdrEbStates(){
  CH_TIME("sisdc::computeCdrEbStates");
  if(m_verbosity > 5){
    pout() << "sisdc::computeCdrEbStates" << endl;
  }

  Vector<EBAMRIVData*>   eb_gradients;
  Vector<EBAMRIVData*>   eb_states;
  Vector<EBAMRCellData*> cdr_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    RefCountedPtr<CdrStorage>& storage = sisdc::get_CdrStorage(solver_it);

    cdr_states.push_back(&(solver->getPhi()));
    eb_states.push_back(&(storage->getEbState()));
    eb_gradients.push_back(&(storage->getEbGrad()));
    cdr_gradients.push_back(&(storage->getGradient())); // Should already have been computed
  }

  // Extrapolate states to the EB and floor them so we cannot get negative values on the boundary. This
  // won't hurt mass conservation because the mass hasn't been injected yet
  sisdc::extrapolate_to_eb(eb_states, m_cdr->getPhase(), cdr_states);
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    DataOps::floor(*eb_states[idx], 0.0);
  }

  // We should already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIVData eb_gradient;
  m_amr->allocate(eb_gradient, m_cdr->getPhase(), SpaceDim);
  for (int i = 0; i < cdr_states.size(); i++){
    sisdc::extrapolate_to_eb(eb_gradient, m_cdr->getPhase(), *cdr_gradients[i]);
    sisdc::project_flux(*eb_gradients[i], eb_gradient);
  }
}

void sisdc::computeCdrEbStates(const Vector<EBAMRCellData*>& a_phis){
  CH_TIME("sisdc::computeCdrEbStates(vec)");
  if(m_verbosity > 5){
    pout() << "sisdc::computeCdrEbStates(vec)" << endl;
  }

  Vector<EBAMRIVData*>   eb_gradients;
  Vector<EBAMRIVData*>   eb_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrStorage>& storage = sisdc::get_CdrStorage(solver_it);

    eb_states.push_back(&(storage->getEbState()));
    eb_gradients.push_back(&(storage->getEbGrad()));
    cdr_gradients.push_back(&(storage->getGradient())); // Should already have been computed
  }

  // Extrapolate states to the EB and floor them so we cannot get negative values on the boundary. This
  // won't hurt mass conservation because the mass hasn't been injected yet
  sisdc::extrapolate_to_eb(eb_states, m_cdr->getPhase(), a_phis);
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    DataOps::floor(*eb_states[idx], 0.0);
  }

  // We should already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIVData eb_gradient;
  m_amr->allocate(eb_gradient, m_cdr->getPhase(), SpaceDim);
  for (int i = 0; i < a_phis.size(); i++){
    sisdc::extrapolate_to_eb(eb_gradient, m_cdr->getPhase(), *cdr_gradients[i]);
    sisdc::project_flux(*eb_gradients[i], eb_gradient);
  }
}

void sisdc::computeCdrDomainStates(){
  CH_TIME("sisdc::computeCdrDomainStates");
  if(m_verbosity > 5){
    pout() << "sisdc::computeCdrDomainStates" << endl;
  }

  Vector<EBAMRIFData*>   domain_gradients;
  Vector<EBAMRIFData*>   domain_states;
  Vector<EBAMRCellData*> cdr_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    RefCountedPtr<CdrStorage>& storage = sisdc::get_CdrStorage(solver_it);

    cdr_states.push_back(&(solver->getPhi()));
    domain_states.push_back(&(storage->getDomainState()));
    domain_gradients.push_back(&(storage->getDomainGrad()));
    cdr_gradients.push_back(&(storage->getGradient())); // Should already be computed
  }

  // Extrapolate states to the domain faces
  sisdc::extrapolate_to_domain_faces(domain_states, m_cdr->getPhase(), cdr_states);

  // We already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIFData grad;
  m_amr->allocate(grad, m_cdr->getPhase(), SpaceDim);
  for (int i = 0; i < cdr_states.size(); i++){
    sisdc::extrapolate_to_domain_faces(grad, m_cdr->getPhase(), *cdr_gradients[i]);
    sisdc::project_domain(*domain_gradients[i], grad);
  }
}

void sisdc::computeCdrDomainStates(const Vector<EBAMRCellData*>& a_phis){
  CH_TIME("sisdc::computeCdrDomainStates");
  if(m_verbosity > 5){
    pout() << "sisdc::computeCdrDomainStates" << endl;
  }

  Vector<EBAMRIFData*>   domain_gradients;
  Vector<EBAMRIFData*>   domain_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);

    domain_states.push_back(&(storage->getDomainState()));
    domain_gradients.push_back(&(storage->getDomainGrad()));
    cdr_gradients.push_back(&(storage->getGradient()));
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

void sisdc::compute_cdr_fluxes(const Real a_time){
  CH_TIME("sisdc::compute_cdr_fluxes");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_cdr_fluxes" << endl;
  }

  this->compute_cdr_fluxes(m_cdr->getPhis(), a_time);
}

void sisdc::compute_cdr_fluxes(const Vector<EBAMRCellData*>& a_phis, const Real a_time){
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

  cdr_fluxes = m_cdr->getEbFlux();

  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);

    EBAMRIVData& dens_eb = storage->getEbState();
    EBAMRIVData& velo_eb = storage->getEbVelo();
    EBAMRIVData& flux_eb = storage->getEbFlux();
    EBAMRIVData& grad_eb = storage->getEbGrad();

    extrap_cdr_densities.push_back(&dens_eb);  // Computed in computeCdrEbStates
    extrap_cdr_velocities.push_back(&velo_eb); // Not yet computed
    extrap_cdr_fluxes.push_back(&flux_eb);     // Not yet computed
    extrap_cdr_gradients.push_back(&grad_eb);  // Computed in computeCdrEbStates
  }


  // Extrapolate densities, velocities, and fluxes
  Vector<EBAMRCellData*> cdr_velocities = m_cdr->getVelocities();
  TimeStepper::compute_extrapolated_fluxes(extrap_cdr_fluxes, a_phis, cdr_velocities, m_cdr->getPhase());
  TimeStepper::compute_extrapolated_velocities(extrap_cdr_velocities, cdr_velocities, m_cdr->getPhase());

  // Compute RTE flux on the boundary
  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver   = solver_it();
    RefCountedPtr<RtStorage>& storage = this->get_RtStorage(solver_it);

    EBAMRIVData& flux_eb = storage->getEbFlux();
    solver->computeBoundaryFlux(flux_eb, solver->getPhi());
    extrap_rte_fluxes.push_back(&flux_eb);
  }

  const EBAMRIVData& E = m_fieldSolver_scratch->getElectricFieldEb();
  TimeStepper::compute_cdr_fluxes(cdr_fluxes,
				   extrap_cdr_fluxes,
				   extrap_cdr_densities,
				   extrap_cdr_velocities,
				   extrap_cdr_gradients,
				   extrap_rte_fluxes,
				   E,
				   a_time);
}

void sisdc::computeCdrDomainFluxes(const Real a_time){
  CH_TIME("sisdc::computeCdrDomainFluxes");
  if(m_verbosity > 5){
    pout() << "sisdc::computeCdrDomainFluxes" << endl;
  }

  this->computeCdrDomainFluxes(m_cdr->getPhis(), a_time);
}

void sisdc::computeCdrDomainFluxes(const Vector<EBAMRCellData*>& a_phis, const Real a_time){
  CH_TIME("sisdc::computeCdrDomainFluxes(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "sisdc::computeCdrDomainFluxes(Vector<EBAMRCellData*>, Real)" << endl;
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
    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);

    EBAMRIFData& dens_domain = storage->getDomainState();
    EBAMRIFData& velo_domain = storage->getDomainVelo();
    EBAMRIFData& flux_domain = storage->getDomainFlux();
    EBAMRIFData& grad_domain = storage->getDomainGrad();
    EBAMRCellData& gradient  = storage->getGradient();

    extrap_cdr_densities.push_back(&dens_domain);  // Has not been computed
    extrap_cdr_velocities.push_back(&velo_domain); // Has not been computed
    extrap_cdr_fluxes.push_back(&flux_domain);     // Has not been computed
    extrap_cdr_gradients.push_back(&grad_domain);  // Has not been computed
    cdr_gradients.push_back(&gradient);
  }

  // Compute extrapolated velocities and fluxes at the domain faces
  this->extrapolate_to_domain_faces(extrap_cdr_densities,         m_cdr->getPhase(), a_phis);
  this->extrapolate_vector_to_domain_faces(extrap_cdr_velocities, m_cdr->getPhase(), cdr_velocities);
  this->compute_extrapolated_domain_fluxes(extrap_cdr_fluxes,     a_phis,           cdr_velocities, m_cdr->getPhase());
  this->extrapolate_vector_to_domain_faces(extrap_cdr_gradients,  m_cdr->getPhase(), cdr_gradients);

  // Compute RTE flux on domain faces
  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver   = solver_it();
    RefCountedPtr<RtStorage>& storage = this->get_RtStorage(solver_it);

    EBAMRIFData& domain_flux = storage->getDomainFlux();
    solver->computeDomainFlux(domain_flux, solver->getPhi());
    extrap_rte_fluxes.push_back(&domain_flux);
  }

  const EBAMRIFData& E = m_fieldSolver_scratch->getElectricFieldDomain();

  // This fills the solvers' domain fluxes
  TimeStepper::computeCdrDomainFluxes(cdr_fluxes,
					  extrap_cdr_fluxes,
					  extrap_cdr_densities,
					  extrap_cdr_velocities,
					  extrap_cdr_gradients,
					  extrap_rte_fluxes,
					  E,
					  a_time);
}

void sisdc::computeSigmaFlux(){
  CH_TIME("sisdc::computeSigmaFlux");
  if(m_verbosity > 5){
    pout() << "sisdc::computeSigmaFlux" << endl;
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

void sisdc::compute_cdr_sources(const Real a_time){
  CH_TIME("sisdc::compute_cdr_sources_into_scratch");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_cdr_sources_into_scratch" << endl;
  }

  this->compute_cdr_sources(m_cdr->getPhis(), a_time);
}

void sisdc::compute_cdr_sources(const Vector<EBAMRCellData*>& a_phis, const Real a_time){
  CH_TIME("sisdc::compute_cdr_sources(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "sisdc::compute_cdr_sources(Vector<EBAMRCellData*>, Real)" << endl;
  }
  
  Vector<EBAMRCellData*> cdr_sources = m_cdr->getSources();
  Vector<EBAMRCellData*> rte_states  = m_rte->getPhis();
  EBAMRCellData& E                   = m_fieldSolver_scratch->getElectricFieldCell();

  Vector<EBAMRCellData*> cdr_gradients;
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);
    cdr_gradients.push_back(&(storage->getGradient())); // These should already have been computed
  }

  TimeStepper::compute_cdr_sources(cdr_sources, a_phis, cdr_gradients, rte_states, E, a_time, centering::cell_center);
}

void sisdc::updateField(){
  CH_TIME("sisdc::updateField(solver)");
  if(m_verbosity > 5){
    pout() << "sisdc::updateField(solver)" << endl;
  }
  
  if(m_do_poisson){ // Solve Poisson equation
    if((m_timeStep +1) % m_fast_poisson == 0){
      TimeStepper::solve_poisson();
      this->compute_E_into_scratch();
    }
  }
}

void sisdc::updateField(const Vector<EBAMRCellData*>& a_densities, const EBAMRIVData& a_sigma){
  CH_TIME("sisdc::updateField(full)");
  if(m_verbosity > 5){
    pout() << "sisdc::updateField(full)" << endl;
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

void sisdc::update_stationary_rte(const Real a_time){
  CH_TIME("sisdc::update_stationary_rte(solver)");
  if(m_verbosity > 5){
    pout() << "sisdc::update_stationary_rte(solver)" << endl;
  }
  
  if(m_do_rte){
    if((m_timeStep + 1) % m_fast_rte == 0){
      if(m_rte->isStationary()){
	Vector<EBAMRCellData*> rte_states  = m_rte->getPhis();
	Vector<EBAMRCellData*> rte_sources = m_rte->getSources();
	Vector<EBAMRCellData*> cdr_states  = m_cdr->getPhis();

	EBAMRCellData& E = m_fieldSolver_scratch->getElectricFieldCell();

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
    if((m_timeStep + 1) % m_fast_rte == 0){
      if(m_rte->isStationary()){
	Vector<EBAMRCellData*> rte_states  = m_rte->getPhis();
	Vector<EBAMRCellData*> rte_sources = m_rte->getSources();

	EBAMRCellData& E = m_fieldSolver_scratch->getElectricFieldCell();

	const Real dummy_dt = 0.0;
	this->solve_rte(rte_states, rte_sources, a_cdr_states, E, a_time, dummy_dt, centering::cell_center);
      }
    }
  }
}

void sisdc::integrate_rte(const Real a_dt, const int a_m, const bool a_lagged_terms){
  CH_TIME("sisdc::integrate_rte(full)");
  if(m_verbosity > 5){
    pout() << "sisdc::integrate_rte(full)" << endl;
  }

  if(m_do_rte){
    if((m_timeStep + 1) % m_fast_rte == 0){
      if(!(m_rte->isStationary())){
	const Real time = m_time + m_dtm[a_m]; // This is the current time
      
	Vector<EBAMRCellData*>  rte_states  = m_rte->getPhis();
	Vector<EBAMRCellData*>  rte_sources = m_rte->getSources();
	Vector<EBAMRCellData*>  cdr_states  = sisdc::get_cdr_phik(a_m);
	EBAMRCellData& E = m_fieldSolver_scratch->getElectricFieldCell();
	this->solve_rte(rte_states, rte_sources, cdr_states, E, time, a_dt, centering::cell_center);
#if 1 // Debug
	MayDay::Abort("sisdc::integrate_rte - shouldn't happen");
#endif
      }
    }
  }
}

void sisdc::updateDiffusionCoefficients(){
  CH_TIME("sisdc::updateDiffusionCoefficients");
  if(m_verbosity > 5){
    pout() << "sisdc::updateDiffusionCoefficients" << endl;
  }
  TimeStepper::compute_cdr_diffusion(m_fieldSolver_scratch->getElectricFieldCell(), m_fieldSolver_scratch->getElectricFieldEb());
}

Vector<EBAMRCellData*> sisdc::get_cdr_errors(){
  CH_TIME("sisdc::get_cdr_errors");
  if(m_verbosity > 5){
    pout() << "sisdc::get_cdr_errors" << endl;
  }
  
  Vector<EBAMRCellData*> ret;
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrStorage>& storage = sisdc::get_CdrStorage(solver_it);
    ret.push_back(&(storage->getError()));
  }

  return ret;
}

Vector<EBAMRCellData*> sisdc::get_cdr_phik(const int a_m){
  CH_TIME("sisdc::get_cdr_phik");
  if(m_verbosity > 5){
    pout() << "sisdc::get_cdr_phik" << endl;
  }
  
  Vector<EBAMRCellData*> ret;
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrStorage>& storage = sisdc::get_CdrStorage(solver_it);
    ret.push_back(&(storage->getPhi()[a_m]));
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

void sisdc::writeStepProfile(const Real a_dt,
			       const Real a_error,
			       const int  a_substeps,
			       const int  a_corrections,
			       const int  a_rejections){
  CH_TIME("sissdc::writeStepProfile");
  if(m_verbosity > 5){
    pout() << "sisdc::writeStepProfile" << endl;
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

void sisdc::reset_finer_flux_registers_level(const int a_lvl,
					     const int a_coarsest_level,
					     const int a_finestLevel){
  CH_TIME("sisdc::resetFluxRegisters_level");
  if(m_verbosity > 5){
    pout() << "sisdc::resetFluxRegisters_level" << endl;
  }

  const phase::which_phase phase = m_cdr->getPhase();
  const bool has_fine = a_lvl < a_finestLevel;
  if(has_fine){
    EBFluxRegister* fluxreg_fine = m_amr->getFluxRegister(phase)[a_lvl];
    fluxreg_fine->setToZero();
  }
}

void sisdc::reset_redist_registers_level(const int a_lvl, const int a_coarsest_level, const int a_finestLevel){
  CH_TIME("sisdc::reset_redist_registers_level");
  if(m_verbosity > 5){
    pout() << "sisdc::reset_redist_registers_level" << endl;
  }

  const phase::which_phase phase = m_cdr->getPhase();
  EBLevelRedist& level_redist = *(m_amr->getLevelRedist(phase)[a_lvl]);
  level_redist.setToZero();

  if(m_amr->getEbCf()){
    const bool has_fine = a_lvl < a_finestLevel;
    const bool has_coar = a_lvl > a_coarsest_level;

    RefCountedPtr<EBFineToCoarRedist>& fine2coar_redist = m_amr->getFineToCoarRedist(m_cdr->getPhase())[a_lvl];
    RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->getCoarToFineRedist(m_cdr->getPhase())[a_lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->getCoarToCoarRedist(m_cdr->getPhase())[a_lvl];

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
				  const int             a_finestLevel,
				  const Real            a_dt){

  // Increment the coarser flux register and initialize the finer flux register. a_flux holds phi*vel which we can use
  const bool has_fine = a_lvl < a_finestLevel;
  const bool has_coar = a_lvl > a_coarsest_level;

  const phase::which_phase phase = m_cdr->getPhase();
  const Interval interv(a_solver, a_solver);
  EBFluxRegister* fluxreg_fine = NULL;
  EBFluxRegister* fluxreg_coar = NULL;

  // Remember, register on a_lvl holds flux between level a_lvl and a_lvl+1
  if(has_fine) fluxreg_fine = m_amr->getFluxRegister(phase)[a_lvl];   
  if(has_coar) fluxreg_coar = m_amr->getFluxRegister(phase)[a_lvl-1]; 

  const DisjointBoxLayout& dbl = m_amr->getGrids()[a_lvl];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(phase)[a_lvl];

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

void sisdc::update_redist_register(const LevelData<BaseIVFAB<Real> >& a_massDifference, const int a_solver, const int a_lvl){
  CH_TIME("CdrSolver::update_redist_register");
  if(m_verbosity > 5){
    pout() << "sisdc::::update_redist_register" << endl;
  }

  const Interval interv(a_solver, a_solver);
  const DisjointBoxLayout& dbl = m_amr->getGrids()[a_lvl];

  const phase::which_phase phase = m_cdr->getPhase();
  EBLevelRedist& level_redist = *(m_amr->getLevelRedist(phase)[a_lvl]);

  // Again, this is a bit stupid but to get the correct data from the correct interval, we have to do this
  EBAMRIVData diff;
  m_amr->allocate(diff, m_cdr->getPhase(), m_plaskin->get_num_species());
  a_massDifference.localCopyTo(Interval(0,0), *diff[a_lvl], interv);
  diff[a_lvl]->exchange(interv);

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    //    level_redist.increment(a_massDifference[dit()], dit(), interv);
    level_redist.increment((*diff[a_lvl])[dit()], dit(), interv);
  }
}

void sisdc::update_coarse_fine_register(const LevelData<BaseIVFAB<Real> >& a_massDifference,
					const int a_solver,
					const int a_lvl,
					const int a_coarsest_level,
					const int a_finestLevel){
  CH_TIME("CdrSolver::update_redist_register");
  if(m_verbosity > 5){
    pout() << "sisdc::::update_redist_register" << endl;
  }

  const Interval interv(a_solver, a_solver);
  const DisjointBoxLayout& dbl = m_amr->getGrids()[a_lvl];
  const bool has_coar = a_lvl > a_coarsest_level;
  const bool has_fine = a_lvl < a_finestLevel;

  if(has_coar || has_fine){

    const Real dx = m_amr->getDx()[a_lvl];
    
    // Again, this is a bit stupid but to get the correct data from the correct interval, we have to do this
    EBAMRIVData diff;
    m_amr->allocate(diff, m_cdr->getPhase(), m_plaskin->get_num_species());
    DataOps::setValue(*diff[a_lvl], 0.0);
    a_massDifference.localCopyTo(Interval(0,0), *diff[a_lvl], interv);
    diff[a_lvl]->exchange();

    RefCountedPtr<EBFineToCoarRedist>& fine2coar_redist = m_amr->getFineToCoarRedist(m_cdr->getPhase())[a_lvl];
    RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->getCoarToFineRedist(m_cdr->getPhase())[a_lvl];
    RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->getCoarToCoarRedist(m_cdr->getPhase())[a_lvl];


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
      RefCountedPtr<EBFluxRegister>& fluxreg = m_amr->getFluxRegister(m_cdr->getPhase())[a_lvl];
      fluxreg->incrementRedistRegister(*coar2fine_redist, interv, -dx);
      fluxreg->incrementRedistRegister(*coar2coar_redist, interv, -dx);
    }
  }
  
}

void sisdc::reflux_level(EBAMRCellData& a_phi,
			 const int      a_solver,
			 const int      a_lvl,
			 const int      a_coarsest_level,
			 const int      a_finestLevel,
			 const Real     a_scale){
  CH_TIME("sisdc::reflux_level");
  if(m_verbosity > 5){
    pout() << "sisdc::::reflux_level" << endl;
  }
  
  const phase::which_phase phase = m_cdr->getPhase();
  const Interval soln_interv(0, 0);
  const Interval flux_interv(a_solver, a_solver);
  const Real dx = m_amr->getDx()[a_lvl];
  RefCountedPtr<EBFluxRegister >& fluxreg = m_amr->getFluxRegister(phase)[a_lvl];
  fluxreg->reflux(*a_phi[a_lvl], soln_interv, flux_interv, 1./dx);
}

void sisdc::redist_level(LevelData<EBCellFAB>&       a_phi,
			 const int                   a_solver,   
			 const LevelData<EBCellFAB>& a_weights,
			 const int                   a_lvl){
  CH_TIME("sisdc::redist_level");
  if(m_verbosity > 5){
    pout() << "sisdc::redist_level" << endl;
  }
  const phase::which_phase phase = m_cdr->getPhase();
  const Interval solver_interv(0, 0);
  const Interval redist_interv(a_solver, a_solver);
  EBLevelRedist& level_redist = *(m_amr->getLevelRedist(phase)[a_lvl]);
  if(m_cdr->get_mass_redist()){
    level_redist.resetWeights(a_weights, 0);
  }
  level_redist.redistribute(a_phi, redist_interv, solver_interv);
}

void sisdc::integrateAdvection_subcycle(const Real a_dt, const int a_m, const bool a_lagged_terms){
  CH_TIME("sisdc::integrateAdvection_subcycle");
  if(m_verbosity > 5){
    pout() << "sisdc::integrateAdvection_subcycle" << endl;
  }

  // Required storages for advection outside of CdrSolver. Don't need all of these for every species, fortunately. 
  const int comp       = 0;
  const int ncomp      = 1;
  const int redist_rad = m_amr->getRedistributionRadius();

  const phase::which_phase phase = m_cdr->getPhase();

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
  const int fine_lvl = m_amr->getFinestLevel();

  Vector<Real> tnew(1 + fine_lvl, m_tm[a_m]);
  Vector<Real> told(1 + fine_lvl, m_tm[a_m]);

  // Advance phi[a_m] to phi[a_m+1] using subcycling. First, compute a dt that is below the CFL limit
  int tref = 1;
  for (int lvl = coar_lvl; lvl < fine_lvl; lvl++){
    tref = tref*m_amr->getRefinementRatios()[lvl];
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

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    CdrGodunov* gdnv = dynamic_cast<CdrGodunov*>(&(*solver));
    if(gdnv == NULL){
      MayDay::Abort("sisdc::subcycle_compute_advection_velocities - Only CdrGodunov can subcycle these days...");
    }
    else{
      if(solver->isMobile()){
	solver->averageVelocityToFaces();
      }
    }
  }
}

void sisdc::subcycle_copy_states(const int a_m){
  CH_TIME("sisdc::subcycle_copy_states");
  if(m_verbosity > 5){
    pout() << "sisdc::subcycle_copy_states" << endl;
  }

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrStorage>& storage = sisdc::get_CdrStorage(solver_it);
    EBAMRCellData& phi_m1      = storage->getPhi()[a_m+1];
    const EBAMRCellData& phi_m = storage->getPhi()[a_m];
    
    DataOps::copy(phi_m1, phi_m);
  }

  EBAMRIVData& sigma_m1      = m_sigma_scratch->get_sigma()[a_m+1];
  const EBAMRIVData& sigma_m = m_sigma_scratch->get_sigma()[a_m];
  DataOps::copy(sigma_m1, sigma_m);
}
    
void sisdc::subcycle_advect_amr(EBAMRFluxData& a_flux,
				EBAMRFluxData& a_facePhi,
				EBAMRCellData& a_divF_c,
				EBAMRCellData& a_weights,
				EBAMRIVData&   a_nonConservativeDivergence,
				EBAMRIVData&   a_massDifference,
				Vector<Real>&  a_tnew,
				Vector<Real>&  a_told,
				const int      a_m,
				const int      a_lvl,
				const int      a_coarsest_level,
				const int      a_finestLevel,
				const Real     a_dt){
  CH_TIME("sisdc::subcycle_advect_amr");
  if(m_verbosity > 5){
    pout() << "sisdc::subcycle_advect_amr" << endl;
  }

  Real coar_time_old = 0.0;
  Real coar_time_new = 0.0;

  const bool has_fine = a_lvl < a_finestLevel;
  const bool has_coar = a_lvl > a_coarsest_level;

  if(has_coar){
    coar_time_old = a_told[a_lvl-1];
    coar_time_new = a_tnew[a_lvl-1];
  }

  // Prepare level solve
  sisdc::subcycle_copy_current_to_old_states(a_m, a_lvl);
  sisdc::reset_finer_flux_registers_level(a_lvl, a_coarsest_level, a_finestLevel);
  sisdc::reset_redist_registers_level(a_lvl, a_coarsest_level, a_finestLevel);

  // Level solve. Update boundary conditions and source terms on this level
  sisdc::subcycle_update_transport_bc(a_m, a_lvl, a_tnew[a_lvl]);
  if(m_cycle_sources){
    sisdc::subcycle_update_sources(a_m, a_lvl, a_tnew[a_lvl]);
  }
  sisdc::subcycle_integrate_level(*a_flux[a_lvl],
				  *a_facePhi[a_lvl],
				  *a_divF_c[a_lvl],
				  *a_weights[a_lvl],
				  *a_massDifference[a_lvl],
				  *a_nonConservativeDivergence[a_lvl],
				  a_m,
				  a_lvl,
				  a_coarsest_level,
				  a_finestLevel,
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
      nref   = m_amr->getRefinementRatios()[a_lvl];
      dt_ref = a_dt/nref;
    }
    else{ // Now advance a_lvl+1
      int tref = 1;
      for (int lvl = a_coarsest_level; lvl < a_finestLevel; lvl++){
      	tref *= m_amr->getRefinementRatios()[lvl];
      }
      
      Real dt_cfl_level = m_dt_cfl*tref;
      for (int lvl = a_coarsest_level; lvl <= a_lvl; lvl++){
      	dt_cfl_level /= m_amr->getRefinementRatios()[lvl];
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
      sisdc::subcycle_advect_amr(a_flux, a_facePhi, a_divF_c, a_weights, a_nonConservativeDivergence, a_massDifference,
				 a_tnew, a_told, a_m, a_lvl+1, a_coarsest_level, a_finestLevel, dt_ref);
    }


    // Finer level has caught up. Sync levels
    sisdc::subcycle_sync_levels(a_m, a_lvl, a_coarsest_level, a_finestLevel);
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
  
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = sisdc::get_CdrStorage(solver_it);

    cell_states.push_back(storage->getPhi()[a_m+1][a_lvl]); // This is what we are updating
    cell_gradients.push_back(storage->getGradient()[a_lvl]);
    cell_velocities.push_back(solver->getCellCenteredVelocity()[a_lvl]);
    
    eb_states.push_back(storage->getEbState()[a_lvl]);
    eb_velocities.push_back(storage->getEbVelo()[a_lvl]);
    eb_gradients.push_back(storage->getEbGrad()[a_lvl]);
    eb_fluxes.push_back(storage->getEbFlux()[a_lvl]);
    
    domain_states.push_back(storage->getDomainState()[a_lvl]);
    domain_gradients.push_back(storage->getDomainGrad()[a_lvl]);

    solver_eb_fluxes.push_back(solver->getEbFlux()[a_lvl]);
  }

  // Compute all the things that are necessary for BCs
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    RefCountedPtr<CdrStorage>& storage = sisdc::get_CdrStorage(solver_it);

    // Storage we can use for computations
    LevelData<EBCellFAB>&        scratchD    = *storage->getScratchD()[a_lvl];
    LevelData<BaseIVFAB<Real> >& scratchIV_1 = *storage->getEbScratch1()[a_lvl]; // Scalar
    LevelData<BaseIVFAB<Real> >& scratchIV_D = *storage->getEbScratchD()[a_lvl]; // SpaceDim comps

    // 1. Compute cell-centered gradients
    cell_states[idx]->exchange();
    m_amr->computeGradient(*cell_gradients[idx], *cell_states[idx], a_lvl);
    
    // 2. Extrapolate cell-centered gradient to the EB
    TimeStepper::extrapolate_to_eb(scratchIV_D, m_cdr->getPhase(), *cell_gradients[idx], a_lvl);

    // 3. Dot EB-centered gradient with normal vector
    TimeStepper::project_flux(*eb_gradients[idx], scratchIV_D, a_lvl);
    
    // 4. Extrapolate the cell-centered states to the EB
    TimeStepper::extrapolate_to_eb(*eb_states[idx], m_cdr->getPhase(), *cell_states[idx], a_lvl);

    // 5. Extrapolate cell-centered velocities to the EB
    TimeStepper::extrapolate_to_eb(scratchIV_D, m_cdr->getPhase(), *cell_velocities[idx], a_lvl);

    // 6. Project normal velocity
    TimeStepper::project_flux(*eb_velocities[idx], scratchIV_D, a_lvl);

    // 7. Compute the extrapolated flux at the boundary
    cell_velocities[idx]->localCopyTo(scratchD);
    DataOps::multiplyScalar(scratchD, *cell_states[idx]);
    TimeStepper::extrapolate_to_eb(scratchIV_D, m_cdr->getPhase(), scratchD, a_lvl);
    TimeStepper::project_flux(*eb_fluxes[idx], scratchIV_D, a_lvl);
  }

  // Radiative transfer fluxes at boundary
  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver   = solver_it();
    RefCountedPtr<RtStorage>& storage = this->get_RtStorage(solver_it);

    EBAMRIVData& flux_eb = storage->getEbFlux();
    solver->computeBoundaryFlux(flux_eb, solver->getPhi());
    rte_fluxes.push_back(flux_eb[a_lvl]);
  }

  // Electric field at boundary
  const EBAMRIVData& E = m_fieldSolver_scratch->getElectricFieldEb();
  
  // Update the stinking EB fluxes
  TimeStepper::compute_cdr_fluxes(solver_eb_fluxes,
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
  const int num_Photons        = m_plaskin->get_num_Photons();
  
  const DisjointBoxLayout& dbl = m_amr->getGrids()[a_lvl];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_cdr->getPhase())[a_lvl];
  const RealVect origin        = m_physdom->getProbLo();
  const Real dx                = m_amr->getDx()[a_lvl];

  // Stencils for extrapolating things to cell centroids
  const IrregAmrStencil<CentroidInterpolationStencil>& interp_stencils = m_amr->getCentroidInterpolationStencils(m_cdr->getPhase());

  // We must have the gradient of E. This block of code does that. 
  EBAMRCellData grad_E, E_norm;
  m_amr->allocate(grad_E, m_cdr->getPhase(), SpaceDim);  // Allocate storage for grad(|E|)
  m_amr->allocate(E_norm, m_cdr->getPhase(), 1);         // Allocate storage for |E|
  const EBAMRCellData& E = m_fieldSolver_scratch->getElectricFieldCell();
  DataOps::vectorLength(*E_norm[a_lvl], *E[a_lvl]);            // Compute |E| on this level
  m_amr->computeGradient(*grad_E[a_lvl], *E_norm[a_lvl], a_lvl);// Compute grad(|E|) on this level

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    Vector<EBCellFAB*> sources(num_species);
    Vector<EBCellFAB*> cdr_densities(num_species);
    Vector<EBCellFAB*> cdr_gradients(num_species);
    Vector<EBCellFAB*> rte_densities(num_Photons);
    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<CdrStorage>& storage = sisdc::get_CdrStorage(solver_it);
      const int idx = solver_it.get_solver();
      EBAMRCellData& src   = solver->getSource();
      EBAMRCellData& phim  = storage->getPhi()[a_m+1];
      EBAMRCellData& gradm = storage->getGradient();
	
      sources[idx]       = &((*src[a_lvl])[dit()]);
      cdr_densities[idx] = &((*phim[a_lvl])[dit()]);
      cdr_gradients[idx] = &((*gradm[a_lvl])[dit()]);
    }
    for (RtIterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
      const int idx = solver_it.get_solver();
      RefCountedPtr<RtSolver>& solver = solver_it();
      EBAMRCellData& state = solver->getPhi();
      rte_densities[idx] = &((*state[a_lvl])[dit()]);
    }
    const EBCellFAB& e  = (*E[a_lvl])[dit()];
    const EBCellFAB& gE = (*grad_E[a_lvl])[dit()];

    // This does all cells
    TimeStepper::compute_cdr_sources_reg(sources,
					  cdr_densities,
					  cdr_gradients,
					  rte_densities,
					  e,
					  gE,
					  dbl.get(dit()),
					  a_time,
					  dx);

    // Have to redo irregular cells
    TimeStepper::compute_cdr_sources_irreg(sources,
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

void sisdc::subcycle_sync_levels(const int a_m, const int a_lvl, const int a_coarsest_level, const int a_finestLevel){
  CH_TIME("sisdc::subcycle_sync_levels");
  if(m_verbosity > 5){
    pout() << "sisdc::subcycle_sync_levels" << endl;
  }

  // Finer level has reached this level. Average down solution on this level and reflux mass.
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int solver_idx                = solver_it.get_solver();
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = sisdc::get_CdrStorage(solver_it);

    EBAMRCellData& state = storage->getPhi()[a_m+1]; // This is the one that we update

    // Since the redistribution registers don't allow components, do some transient storage stuff
    
    
    // Reflux state
    if(solver->isMobile()){
      m_amr->averageDown(state, m_cdr->getPhase(), a_lvl);
      sisdc::reflux_level(state, solver_idx, a_lvl, a_coarsest_level, a_finestLevel, 1.0);
      // EBCF related code. 
      if(m_amr->getEbCf()){
	const Interval inter0(0,0);
	const Interval interv(solver_idx, solver_idx);
	
	const bool has_fine = a_lvl < a_finestLevel;
	const bool has_coar = a_lvl > a_coarsest_level;

	// Bah, extra storage becase redistribution registers don't let me use different for mass diffs and target FAB
	EBAMRCellData dummy;
	m_amr->allocate(dummy, m_cdr->getPhase(), m_plaskin->get_num_species());
	DataOps::setValue(dummy, 0.0);

	RefCountedPtr<EBCoarToFineRedist>& coar2fine_redist = m_amr->getCoarToFineRedist(m_cdr->getPhase())[a_lvl];
	RefCountedPtr<EBCoarToCoarRedist>& coar2coar_redist = m_amr->getCoarToCoarRedist(m_cdr->getPhase())[a_lvl];
	RefCountedPtr<EBFineToCoarRedist>& fine2coar_redist = m_amr->getFineToCoarRedist(m_cdr->getPhase())[a_lvl];


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
    m_amr->averageDown(state, m_cdr->getPhase(), a_lvl);
    state[a_lvl]->exchange();
  }
}

void sisdc::subcycle_copy_current_to_old_states(const int a_m, const int a_lvl){
  CH_TIME("sisdc::subcycle_copy_current_to_old_states");
  if(m_verbosity > 5){
    pout() << "sisdc::subcycle_copy_current_to_old_states" << endl;
  }

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrStorage>& storage = sisdc::get_CdrStorage(solver_it);

    EBAMRCellData& old = storage->getOld();
    const EBAMRCellData& current = storage->getPhi()[a_m+1];

    current[a_lvl]->localCopyTo(*old[a_lvl]);
    old[a_lvl]->exchange();
  }
}
  
void sisdc::subcycle_integrate_level(LevelData<EBFluxFAB>&        a_flux,
				     LevelData<EBFluxFAB>&        a_facePhi,
				     LevelData<EBCellFAB>&        a_divF_c,
				     LevelData<EBCellFAB>&        a_weights,
				     LevelData<BaseIVFAB<Real> >& a_massDifference,
				     LevelData<BaseIVFAB<Real> >& a_nonConservativeDivergence,
				     const int                    a_m,
				     const int                    a_lvl,
				     const int                    a_coarsest_level,
				     const int                    a_finestLevel,
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

  const bool ebcf     = m_amr->getEbCf();
  const bool has_coar = a_lvl > a_coarsest_level;
  const bool has_fine = a_lvl < a_finestLevel;
  
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int solver_idx = solver_it.get_solver();
    RefCountedPtr<CdrStorage>& storage = sisdc::get_CdrStorage(solver_it);
    RefCountedPtr<CdrSolver>& solver   = solver_it();

    LevelData<EBCellFAB>& state_m1 = (*storage->getPhi()[a_m+1][a_lvl]); // We will update this one.
    LevelData<EBCellFAB>& source   = *solver->getSource()[a_lvl];

    if(solver->isMobile()){
      
      CdrGodunov* gdnv = (CdrGodunov*) (&(*solver));
      LevelData<EBCellFAB>* coar_old = NULL;
      LevelData<EBCellFAB>* coar_new = NULL;

      if(has_coar){
	coar_old = storage->getOld()[a_lvl-1];
	coar_new = storage->getPhi()[a_m+1][a_lvl-1];
      }
      
      // Advect to faces and compute fluxes on face centers, and compute the conservative divergence on regular cells
      const Real extr_dt = m_extrap_advect ? 2.0*m_extrap_dt*a_dt : 0.0;
      gdnv->advectToFaces(a_facePhi, state_m1, coar_old, coar_new, a_time, a_coar_time_old,a_coar_time_new, a_lvl, extr_dt);
      gdnv->new_computeFlux(a_flux, a_facePhi, a_lvl);
      gdnv->conservativeDivergenceRegular(a_divF_c, a_flux, a_lvl);

      // Set up flux interpolant and compute conservative flux on irregular cells
      LevelData<BaseIFFAB<Real> > flux_interp[SpaceDim];
      const LevelData<BaseIVFAB<Real> >& ebflux = *solver->getEbFlux()[a_lvl];
      gdnv->setupFluxInterpolant(flux_interp, a_flux, a_lvl); // Compute interpolant
      gdnv->interpolateFluxToFaceCentroids(flux_interp, a_lvl);  // Interpolant now holds face centroid-centered fluxes
      gdnv->computeDivF_irreg(a_divF_c, flux_interp, ebflux, a_lvl); // Update the conservative divergence in the cut cells

      // Recess: So far the conservative divergence is scaled by 1/dx but not yet divided by the volume fraction. This
      // means that the actual advance without the hybrid stuff would be
      //
      // new_state -= dt*a_divF_c/kappa
      //
      // The stuff below was originally written for d(phi)/dt = Div(F) rather than d(phi)/dt = -Div(F) so that's why
      // there's a (-a_dt) in all the stuff below. This design choice was made because I am, in fact, an ass.

      // Compute the nonconservative and hybrid divergences (hybrid put on storage for divF_c, which is lost)
      gdnv->nonConservativeDivergence(a_nonConservativeDivergence, a_facePhi, a_lvl);
      gdnv->hybridDivergence(a_divF_c, a_massDifference, a_nonConservativeDivergence, a_lvl); // Puts hybrid in a_divF_c. mass_diff as usual without dt,
      DataOps::scale(a_massDifference, -a_dt);                              // Sign convention

      // Update flux and redistribution registers
      gdnv->new_computeFlux(a_flux, a_facePhi, a_lvl);
      sisdc::update_flux_registers(a_flux, solver_idx, a_lvl, a_coarsest_level, a_finestLevel, a_dt);
      sisdc::update_redist_register(a_massDifference, solver_idx, a_lvl);

      // If we have EBCF, we must update those as well
      if(ebcf){
	sisdc::update_coarse_fine_register(a_massDifference, solver_idx, a_lvl, a_coarsest_level, a_finestLevel);
      }
      
      if(m_cdr->get_mass_redist()){
	DataOps::setValue(a_weights, 0.0);
	DataOps::incr(a_weights, state_m1, 1.0);
      }

      // Euler advance with redistribution
      DataOps::incr(state_m1, a_divF_c, -a_dt);
      sisdc::redist_level(state_m1, solver_idx, a_weights, a_lvl);
    }

    // Add in the source term
    if(m_cycle_sources){
      DataOps::incr(state_m1, source, a_dt);
    }
    state_m1.exchange();
  }
}

void sisdc::storeSolvers(){
  CH_TIME("sisdc::storeSolvers");
  if(m_verbosity > 5){
    pout() << "sisdc::storeSolvers" << endl;
  }

  // SISDC does not manipulate cdr and sigma solvers until the end of the time step. Only need to do
  // Poisson and RTE here. 

  // Poisson
  MFAMRCellData& previous    = m_fieldSolver_scratch->getPrevious();
  const MFAMRCellData& state = m_fieldSolver->getPotential();
  DataOps::copy(previous, state);

  // RTE
  for (RtIterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtStorage>& storage     = sisdc::get_RtStorage(solver_it);
    const RefCountedPtr<RtSolver>& solver = solver_it();

    EBAMRCellData& previous = storage->getPrevious();
    const EBAMRCellData& state = solver->getPhi();

    DataOps::copy(previous, state);
  }
}

void sisdc::restoreSolvers(){
  CH_TIME("sisdc::restoreSolvers");
  if(m_verbosity > 5){
    pout() << "sisdc::restoreSolvers" << endl;
  }

  // SISDC does not manipulate cdr and sigma solvers until the end of the time step. Only need to do
  // Poisson and RTE here. 

  // Poisson
  MFAMRCellData& state = m_fieldSolver->getPotential();
  const MFAMRCellData& previous    = m_fieldSolver_scratch->getPrevious();

  DataOps::copy(state, previous);

  // RTE
  for (RtIterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtStorage>& storage     = sisdc::get_RtStorage(solver_it);
    RefCountedPtr<RtSolver>& solver = solver_it();

    EBAMRCellData& previous = storage->getPrevious();
    EBAMRCellData& state = solver->getPhi();

    DataOps::copy(state, previous);
  }
#include "CD_NamespaceFooter.H"
