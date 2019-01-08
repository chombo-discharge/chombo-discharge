/*!
  @file   strang2.cpp
  @brief  Implementation of strang2.H
  @author Robert Marskar
  @date   Sept. 2018
*/

#include "strang2.H"
#include "strang2_storage.H"
#include "data_ops.H"
#include "units.H"
#include "cdr_tga.H"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <ParmParse.H>

typedef strang2::cdr_storage     cdr_storage;
typedef strang2::poisson_storage poisson_storage;
typedef strang2::rte_storage     rte_storage;
typedef strang2::sigma_storage   sigma_storage;

strang2::strang2(){
  m_advective_order     = 2;
  m_minCFL              = 0.1;
  m_maxCFL              = 0.9;   // Overriden below since Strang and Godunov splittings have different CFLs
  m_err_thresh          = 1.E-4;
  m_safety              = 0.9;
  m_error_norm          = 2;

  // This is maximum order supported right now
  m_max_advection_order = 3;

  m_adaptive_dt    = true;
  m_compute_error  = true;
  m_have_dtf       = false;
  m_fixed_order    = false;

  m_splitting      = "strang";

  // Basically only for debugging
  m_use_embedded   = false;
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

  // Get parameters from input script
  {
    ParmParse pp("strang2");

    std::string str;

    if(pp.contains("splitting")){
      pp.get("splitting", str);
      if(str == "godunov"){
	m_splitting = "godunov";
      }
      else if(str == "strang"){
	m_splitting = "strang";
      }
      else {
	MayDay::Abort("strang2::strang2 - unknown argument 'splitting'. This must be 'godunov' or 'strang'");
      }
    }

    if(m_splitting == "godunov"){
      m_maxCFL = 0.9;
    }
    if(m_splitting == "strang"){
      m_maxCFL = 1.8;
    }

    pp.query("advective_order", m_advective_order);
    pp.query("min_cfl",         m_minCFL);
    pp.query("max_cfl",         m_maxCFL);
    pp.query("max_error",       m_err_thresh);
    pp.query("safety",          m_safety);
    pp.query("error_norm",      m_error_norm);

    if(pp.contains("fixed_order")){
      pp.get("fixed_order", str);
      if(str == "true"){
	m_fixed_order = true;
      }
    }
    if(pp.contains("accept_error")){
      pp.get("accept_error", str);
      if(str == "true"){
	m_use_embedded = true;
      }
    }
    if(pp.contains("compute_error")){
      pp.get("compute_error", str);
      if(str == "true"){
	m_compute_error = true;
      }
      else{
	m_compute_error = false;
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
    if(pp.contains("adaptive_dt")){
      pp.get("adaptive_dt", str);
      if(str == "false"){
	m_adaptive_dt = false;
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
	  pout() << "strang2::strang2 - Turning off advection & source" << endl;
	}
      }
    }
    if(pp.contains("do_diffusion")){
      pp.get("do_diffusion", str);
      if(str == "false"){
	m_do_diffusion = false;
	if(m_verbosity > 2){
	  pout() << "strang2::strang2 - Turning off diffusion" << endl;
	}
      }
    }
    if(pp.contains("do_rte")){
      pp.get("do_rte", str);
      if(str == "false"){
	m_do_rte = false;

	if(m_verbosity > 2){
	  pout() << "strang2::strang2 - Turning off rte" << endl;
	}
      }
    }
    if(pp.contains("do_poisson")){
      pp.get("do_poisson", str);
      if(str == "false"){
	m_do_poisson = false;

	if(m_verbosity > 2){
	  pout() << "strang2::strang2 - Turning off poisson" << endl;
	}
      }
    }
  }

  // Print error if user tries to exceed supported orders
  if(m_advective_order < 2 || m_advective_order > m_max_advection_order){
    pout() << "strang2 ::strang2 - order < 1 or order > 3 requested!." << endl;
    MayDay::Abort("strang2 ::strang2 - order < 1 or order > 3 requested!.");
  }

  // If we don't use adaptive, we can't go beyond CFL
  if(!m_adaptive_dt){
    m_cfl = m_maxCFL;
  }

  // If we're doing adaptive stepping, errors must always be comptued
  if(m_adaptive_dt){
    m_compute_error = true;
  }

  // Strang splitting limited by CFL of 2, Godunov of 1
  if(m_splitting == "strang"){
    m_maxCFL = Min(m_maxCFL, 2.0);
  }
  else if(m_splitting == "godunov"){
    m_maxCFL = Min(m_maxCFL, 1.0);
  }
}

strang2::~strang2(){

}

RefCountedPtr<cdr_storage>& strang2::get_cdr_storage(const cdr_iterator& a_solverit){
  return m_cdr_scratch[a_solverit.get_solver()];
}

RefCountedPtr<rte_storage>& strang2::get_rte_storage(const rte_iterator& a_solverit){
  return m_rte_scratch[a_solverit.get_solver()];
}

Real strang2::restrict_dt(){
  return 1.E99;
}

Real strang2::get_max_error(){
  CH_TIME("strang2::get_max_error");
  if(m_verbosity > 2){
    pout() << "strang2::get_max_error" << endl;
  }

  Real cur_err = m_sigma_error;
  for (int i = 0; i < m_cdr_error.size(); i++){
    cur_err = Max(cur_err, m_cdr_error[i]);
  }

  return cur_err;
}

Real strang2::advance(const Real a_dt){
  CH_TIME("strang2::advance");
  if(m_verbosity > 2){
    pout() << "strang2::advance" << endl;
  }

  // TLDR: This timestepper can do multiple steps per a_dt. If you're using adaptive time integration,
  //       the number of steps are adjusted such that the invididual step errors remain below a specified
  //       threshold. A side effect of this is that the final substep may cause integration beyond a_dt.
  //
  //       PS: When we enter this routine, solvers SHOULD have been filled with valid ready and be ready 
  //           advancement. If you think that this may not be the case, activate the debugging below

  {
#if 0
  this->compute_E_into_scratch();
  this->compute_cdr_gradients();
  this->compute_cdr_velo(m_time);
  this->compute_cdr_sources(m_time);
  time_stepper::compute_cdr_diffusion(m_poisson_scratch->get_E_cell(), m_poisson_scratch->get_E_eb());
#endif
  }

  this->cache_solutions(); // Store old solution. This stores the old solutions in storage->m_cache
  
  Real cfl;              // Adaptive or fixed CFL for advective+source term advancements
  Real actual_dt = a_dt; // This is the actual time step taken. It is equal to or larger than a_dt
  int substeps = 0;      // Number of individual steps taken
  
  // Advection and source term advancements
  if(m_do_advec_src){
    if(m_adaptive_dt){ // Loop for adaptive time stepping

      // If we don't have a good estimate for the substep size, estimate it by equal division. This is
      // probably a bad time step for the first step but it should stabilize fairly quickly.
      if(!m_have_dtf){
	substeps   = ceil(a_dt/(m_maxCFL*m_dt_cfl));
	cfl        = a_dt/(substeps*m_dt_cfl);
	m_dt_adapt = cfl*m_dt_cfl;
	m_have_dtf = true;
      }

      // Start adaptive substepping. The actual_dt is the total time integrated (may be longer than a_dt)
      actual_dt = this->advance_adaptive(substeps, m_dt_adapt, m_time, a_dt);
    }
    else{ // Non-adaptive time-stepping scheme loop
      // Do equal division of the time steps. This integration always lands on a_dt
      substeps   = ceil(a_dt/(m_maxCFL*m_dt_cfl));
      cfl        = a_dt/(substeps*m_dt_cfl);
      m_dt_adapt = cfl*m_dt_cfl;
      actual_dt = this->advance_fixed(substeps, m_dt_adapt);
    }
  }

  // End of the large time step, resolve the Poisson and RTE equations
  this->update_poisson();
  this->update_rte(m_time + actual_dt);
  
  // Recompute source terms and velocities. This puts the cdr solvers back in useable state so that
  // we can reliably use them in the next time step (or anywhere else).
  this->compute_cdr_gradients();
  this->compute_cdr_sources(m_time + a_dt);
  this->compute_cdr_velo(m_time + a_dt);
  time_stepper::compute_cdr_diffusion(m_poisson_scratch->get_E_cell(), m_poisson_scratch->get_E_eb());

  // Print diagnostics
  if(m_print_diagno){
    pout() << "\n";
    pout() << "\t strang2::advance(Real a_dt) breakdown" << endl;
    pout() << "\t ==============================================\n" << endl;
    
    pout() << "\t Convection-source advance: \n ";
    pout() << "\t --------------------------\n";
    pout() << "\t\t Steps     = " << substeps << endl;
    if(m_adaptive_dt){
      pout() << "\t\t Avg. cfl  = " << actual_dt/(substeps*m_dt_cfl) << endl;
      pout() << "\t\t Avg. dt   = " << actual_dt/substeps << endl;
    }
    else{
      pout() << "\t\t Local cfl = " << cfl        << endl;
      pout() << "\t\t Local dt  = " << m_dt_adapt << endl;
    }
  }
  
  return actual_dt;
}

void strang2::advance_diffusion(const Real a_time, const Real a_dt){
  CH_TIME("strang2::advance_diffusion");
  if(m_verbosity > 2){
    pout() << "strang2::advance_diffusion" << endl;
  }

  if(m_do_diffusion){
    this->advance_tga_diffusion(a_time, a_dt);   // This is the 2nd order update
    this->advance_euler_diffusion(a_time, a_dt); // This is the embedded formula 1st order update, it uses the solution from
                                                 // advance_tga_diffusion as initial guess in order to optimize.
  }
}

void strang2::advance_tga_diffusion(const Real a_time, const Real a_dt){
  CH_TIME("strang2::advance_tga_diffusion");
  if(m_verbosity > 2){
    pout() << "strang2::advance_tga_diffusion" << endl;
  }

  const int ncomp = 1;

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
      
    EBAMRCellData& phi = solver->get_state();
    if(solver->is_diffusive()){
      EBAMRCellData old_phi;
      m_amr->allocate(old_phi, phase::gas, ncomp);
      data_ops::copy(old_phi, phi);
      
      cdr_tga* tgasolver = (cdr_tga*) (&(*solver));
      tgasolver->advance_tga(phi, old_phi, a_dt);
    }
  }
}

void strang2::advance_euler_diffusion(const Real a_time, const Real a_dt){
  CH_TIME("strang2::advance_euler_diffusion");
  if(m_verbosity > 2){
    pout() << "strang2::advance_euler_diffusion" << endl;
  }

  const int ncomp = 1;

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
      
    EBAMRCellData& exact_phi = solver->get_state();
    EBAMRCellData& wrong_phi = storage->get_error();
    if(solver->is_diffusive()){

      // We advance the error but take the 2nd order solution as the initial guess
      // since we anticipate that it will be a "close" solution.
      EBAMRCellData old_phi;
      m_amr->allocate(old_phi, phase::gas, ncomp);
      data_ops::copy(old_phi, wrong_phi);
      data_ops::copy(wrong_phi, exact_phi);

      cdr_tga* tgasolver = (cdr_tga*) (&(*solver));
      tgasolver->advance_euler(wrong_phi, old_phi, a_dt);
    }
  }
}

Real strang2::advance_fixed(const int a_substeps, const Real a_dt){
  CH_TIME("strang2::advance_fixed");
  if(m_verbosity > 2){
    pout() << "strang2::advance_fixed" << endl;
  }

  this->compute_E_into_scratch(); // Compute the electric field
  this->compute_cdr_gradients();  // Precompute gradients

  for (int step = 0; step < a_substeps; step++){
    const Real time      = m_time + step*a_dt;
    const bool last_step = (step == a_substeps - 1);

    // If we're doing Strang splitting, we're actually taking two
    // advective steps with 0.5*dt. Since a_dt is the fine time step,
    // we adjust here. 
    Real diffusive_dt = a_dt;
    Real advective_dt;
    if(m_splitting == "godunov"){
      advective_dt = a_dt;
    }
    else if(m_splitting == "strang"){
      advective_dt = 0.5*a_dt;
    }
    else {
      MayDay::Abort("strang2::advance_fixed - unknown splitting");
    }

    // Advection-reaction integration. We must have a place to put u^n (the current solution) before
    // we advance, so it is put in storage->m_previous. We do this because it's easier to fill
    // the solvers with the intermediate Runge-Kutta stage (it lets us use a simpler interface).
    // The 'true' flag indicates that we should compute the embedded formula error, which we only do
    // on the first time step (the embedded splitting uses a Lie splitting)
    this->store_solvers();
    if(m_advective_order == 2){
      this->advance_rk2(time, advective_dt, m_compute_error);
    }
    else if(m_advective_order == 3){
      this->advance_rk3(time, advective_dt, m_compute_error);
    }
    else{
      MayDay::Abort("strang2::advance_fixed - unsupported order requested");
    }

    // In principle, we should update E before computing new diffusion coefficients. Do that, then update
    // the diffusion coefficients and then advance both solutions
    if(m_do_diffusion){
      if(m_consistent_E){
	this->update_poisson();
      }
      if(m_compute_D){
	time_stepper::compute_cdr_diffusion(m_poisson_scratch->get_E_cell(), m_poisson_scratch->get_E_eb());
      }
      this->advance_diffusion(time, diffusive_dt);
    }

    // If we're doing a Strang splitting, we need to advance the final advection-reaction stuff. But for
    // Godunov splitting we don't need to do this. 
    if(m_splitting == "strang"){
      // Strictly speaking, we should update the electric field and RTE equations again
      if(m_consistent_E) this->update_poisson();
      if(m_consistent_rte) this->update_rte(m_time + 0.5*a_dt);

      // We're now doing a full time step, so solvers need to be filled again
      this->store_solvers();
      this->compute_cdr_gradients();
      if(m_compute_v) this->compute_cdr_velo(time);
      if(m_compute_S) this->compute_cdr_sources(time);

      // Now do the second advance. No need to monkey with more error computations since the embedded formula
      // is a Lie split
      if(m_advective_order == 2){
	this->advance_rk2(time, advective_dt, false);
      }
      else if(m_advective_order == 3){
	this->advance_rk3(time, advective_dt, false);
      }
      else{
	MayDay::Abort("strang2::advance_fixed - unsupported order requested");
      }
    }

    // If we're doing one more substep, we need to make sure that the solvers are up to date. If not, let advance()
    // take care of the rest of the synchronization
    if(!last_step){
      // Update the electric field and RTE equations
      if(m_consistent_E)   this->update_poisson();
      if(m_consistent_rte) this->update_rte(time);

      // Compute new velocities and source terms for the next advective step. No need to compute the
      // diffusion coefficients because they are automateically filled before advance_diffusion() call
      if(m_compute_S) this->compute_cdr_gradients();
      if(m_compute_v) this->compute_cdr_velo(time);
      if(m_compute_S) this->compute_cdr_sources(time);
    }
  }

  return a_dt;
}

Real strang2::advance_adaptive(int& a_substeps, Real& a_dt, const Real a_time, const Real a_dtc){
  CH_TIME("strang2::advance_adaptive");
  if(m_verbosity > 2){
    pout() << "strang2::advance_adaptive" << endl;
  }

  MayDay::Abort("strang2::advance_adaptive - not implemented (yet)");
  return 0.0;
}

void strang2::advance_rk2(const Real a_time, const Real a_dt, const bool a_compute_err){
  CH_TIME("strang2::advance_rk2");
  if(m_verbosity > 2){
    pout() << "strang2::advance_rk2" << endl;
  }

  // NOTE: If we use strang splitting, the time step that comes in here has been halved. Furthermore,
  //       the solvers have been filled with source terms and velocities, but not much else.

  { // u^1 = u^n + dt*L(u^n)
    this->compute_E_into_scratch();
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();

    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& rhs       = storage->get_scratch();
      EBAMRCellData& phi       = solver->get_state();
      const EBAMRCellData& src = solver->get_source();

      solver->compute_divF(rhs, phi, 0.0, true); // RHS =  Div(u^n*v^n)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^n*v^n)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^n - Div(u^n*v^n)
      data_ops::incr(phi, rhs, a_dt);            // u^1 = u^n + dt*L(u^n)

      m_amr->average_down(phi, m_cdr->get_phase());
      m_amr->interp_ghost(phi, m_cdr->get_phase());

      data_ops::floor(phi, 0.0);

       // For Heun's method, the embedded error is just u^1 (the forward Euler)
      if(a_compute_err){
	EBAMRCellData& err = storage->get_error();
	data_ops::copy(err, phi);         // err = u^n + dt*L(u^n). Correct for Godunov splitting. 
	if(m_splitting == "strang"){
	  data_ops::incr(err, rhs, a_dt); // err = u^n + 2*dt*L(u^n). Correct for Strang splitting. 
	}
      }
    }

    EBAMRIVData& sigma = m_sigma->get_state();
    EBAMRIVData& rhs   = m_sigma_scratch->get_scratch();
    m_sigma->compute_rhs(rhs);
    data_ops::incr(sigma, rhs, a_dt); // sigma^1 = sigma^n + dt*F^n

    // For Heun's method, the embedded error is just u^1 (the forward Euler)
    if(a_compute_err){
      EBAMRIVData& err = m_sigma_scratch->get_error();
      data_ops::set_value(err, 0.0);
      data_ops::incr(err, sigma, 1.0);  // err = u^n + dt*L(u^n). Correct for Godunov splitting.
      if(m_splitting == "strang"){
	data_ops::incr(err, rhs, a_dt); // err = u^n + 2*dt*L(u^n). Correct for Strang splitting.
      }
    }
  }

  if(m_consistent_E)   this->update_poisson();
  if(m_consistent_rte) this->update_rte(a_time + a_dt);


  { // u^(n+1) = u^n + 0.5*dt*[L(u^n) + L(u^1)]
    this->compute_cdr_gradients();
    if(m_compute_v) this->compute_cdr_velo(a_time + a_dt);
    if(m_compute_S) this->compute_cdr_sources(a_time + a_dt);
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();

    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& rhs       = storage->get_scratch();
      EBAMRCellData& phi       = solver->get_state();
      const EBAMRCellData& src = solver->get_source();
      const EBAMRCellData& pre = storage->get_previous();

      solver->compute_divF(rhs, phi, 0.0, true); // RHS =  Div(u^n*v^n)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^n*v^n)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^n - Div(u^n*v^n)
      
      data_ops::incr(phi, rhs, a_dt);            // u^(n+1) = u^1 + dt*L(u^1)
      data_ops::incr(phi, pre, 1.0);             // u^(n+1) = u^n + u^1 + dt*L(u^1);
      data_ops::scale(phi, 0.5);                 // u^(n+1) = 0.5*[u^n + u^1 + dt*L(u^1)];

      m_amr->average_down(phi, m_cdr->get_phase());
      m_amr->interp_ghost(phi, m_cdr->get_phase());

      data_ops::floor(phi, 0.0);
    }

    EBAMRIVData& sigma     = m_sigma->get_state();
    EBAMRIVData& rhs       = m_sigma_scratch->get_scratch();
    const EBAMRIVData& pre = m_sigma_scratch->get_previous();
    m_sigma->compute_rhs(rhs);
    data_ops::incr(sigma, rhs, a_dt); // sigma^(n+1) = sigma^1 + dt*F^1
    data_ops::incr(sigma, pre, 1.0);  // sigma^(n+1) = sigma^n + sigma^1 + dt*F^1
    data_ops::scale(sigma, 0.5);      // sigma^(n+1) = 0.5*[sigma^n + sigma^1 + dt*F^1]
  }
}

void strang2::advance_rk3(const Real a_time, const Real a_dt, const bool a_compute_err){
  CH_TIME("strang2::advance_rk3");
  if(m_verbosity > 2){
    pout() << "strang2::advance_rk3" << endl;
  }

  // u^1 = u^n + dt*L(u^n)
  { 
    this->compute_E_into_scratch();
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();
  
    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    
      EBAMRCellData& phi       = solver->get_state();    // u^n
      EBAMRCellData& rhs       = storage->get_scratch();
      const EBAMRCellData& src = solver->get_source();
    
      solver->compute_divF(rhs, phi, 0.0, true); // RHS =  Div(u^n*v^n)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^n*v^n)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^n - Div(u^n*v^n)
      data_ops::incr(phi, rhs, a_dt);            // u^1 = u^n + dt*L(u^n)
    
      m_amr->average_down(phi, m_cdr->get_phase());
      m_amr->interp_ghost(phi, m_cdr->get_phase());
    
      data_ops::floor(phi, 0.0);

      // For RK3, the embedded method is u^(n+1) = 0.5*[u^n + u^1 + dt*L(u^1)]. We must
      // make err = u^1 because we will later need to compute 
      if(a_compute_err){
	EBAMRCellData& err = storage->get_error();
	data_ops::copy(err, phi);         // err = u^n + dt*L(u^n). Correct so far for the Godunov splitting. 
	if(m_splitting == "strang"){
	  data_ops::incr(err, rhs, a_dt); // err = u^n + 2*dt*L(u^n). Correct so far for Strang splitting. 
	}
      }
    }

    EBAMRIVData& sigma = m_sigma->get_state();
    EBAMRIVData& rhs   = m_sigma_scratch->get_scratch();
    m_sigma->compute_rhs(rhs);
    data_ops::incr(sigma, rhs, a_dt); // sigma^1 = sigma^n + dt*F^n

    // For Heun's method, the embedded error is just u^1 (the forward Euler)
    if(a_compute_err){
      EBAMRIVData& err = m_sigma_scratch->get_error();
      data_ops::set_value(err, 0.0);
      data_ops::incr(err, sigma, 1.0);  // err = u^n + dt*F^n. Correct so far for the Godunov splitting. 
      if(m_splitting == "strang"){      
	data_ops::incr(err, rhs, a_dt); // err = u^n + 2*dt*F^n. Correct so far for the Strang splitting. 
      }
    }
  }

  if(m_consistent_E)   this->update_poisson();
  if(m_consistent_rte) this->update_rte(a_time + a_dt);

  // Final Heun update for the error in case of Strang splitting. In principle, it should use the fields
  // obtained from the embedded formula update, but we skip that here to save some computations. It's just the error
  // and we don't care if it is exact.
  //
  // This loop is only for Strang splitting because the Godunov splitting shares function evaluations with
  // the full RK3 scheme and therefore has no overhead. 
  if(a_compute_err && m_splitting == "strang"){
    Vector<EBAMRCellData*> err_states;
    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);
      err_states.push_back(&(storage->get_error()));
    }

    this->compute_cdr_gradients(err_states);
    if(m_compute_v) this->compute_cdr_velo(err_states, a_time + a_dt);
    if(m_compute_S) this->compute_cdr_sources(err_states, a_time + a_dt);
    this->compute_cdr_eb_states(err_states);
    this->compute_cdr_fluxes(err_states,a_time);
    this->compute_cdr_domain_fluxes(err_states,a_time);
    this->compute_sigma_flux();

    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);
      
      EBAMRCellData& err       = storage->get_error();    // err = u^n + 2*dt*L(u^n) for Strang splitting. 
      EBAMRCellData& rhs       = storage->get_scratch();
      const EBAMRCellData& pre = storage->get_previous();
      const EBAMRCellData& src = solver->get_source();

      solver->compute_divF(rhs, err, 0.0, true); // RHS = Div(u^1*v^1)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^n*v^n)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^n - Div(u^n*v^n)
      
      data_ops::incr(err, rhs, 2*a_dt);            // err = u^1 + 2*dt*L(u^1). Correct for Strang splitting.
      data_ops::incr(err, pre, 1.0);               // err = u^n + u^1 + 2*dt*L(u^1); 
      data_ops::scale(err, 0.5);                   // err = 0.5*[u^n + u^1 + 2*dt*L(u^1)];

      m_amr->average_down(err, m_cdr->get_phase());
      m_amr->interp_ghost(err, m_cdr->get_phase());
    }

    EBAMRIVData& err = m_sigma_scratch->get_error();     // err = sigma^n + 2*dt*F^n = sigma^1
    EBAMRIVData& rhs = m_sigma_scratch->get_scratch();
    EBAMRIVData& pre = m_sigma_scratch->get_previous();
    m_sigma->compute_rhs(rhs);
    data_ops::incr(err, rhs, 2*a_dt); // err = sigma^1 + 2*dt*F^n. Correct for Strang splitting
    data_ops::incr(err, pre, 2*a_dt); // err = sigma^n + sigma^1 + 2*dt*F^1
    data_ops::scale(err, 0.5);        // err = 0.5*(sigma^n + sigma^1 + 2*dt*F^1)
  }

  // u^2 = 0.25*(3*u^n + u^1 + dt*L(u^1)). For Godunov splitting, the embedded formula is
  // err = 0.5*(u^n + u^1 + dt*L(u^1)). The Strang splitting error is handled in the loop above. 
  { 
    this->compute_cdr_gradients();
    if(m_compute_v) this->compute_cdr_velo(a_time + a_dt);
    if(m_compute_S) this->compute_cdr_sources(a_time + a_dt);
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();

    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    
      EBAMRCellData& rhs       = storage->get_scratch();
      EBAMRCellData& phi       = solver->get_state();     // phi = u^1 = u^n + dt*L(u^n)
      const EBAMRCellData& src = solver->get_source();
      const EBAMRCellData& pre = storage->get_previous(); // pre = u^n
    
      solver->compute_divF(rhs, phi, 0.0, true); // RHS =  Div(u^1*v^1)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^1*v^1)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^1 - Div(u^1*v^1)
      data_ops::incr(phi,  rhs, a_dt);           // u^2 = u^1 + dt*L(u^1)
      data_ops::incr(phi,  pre, 3.0);            // u^2 = 3*u^n + u^1 + dt*L(u^1)
      data_ops::scale(phi, 0.25);                // u^2 = 0.25*[3*u^n + u^1 + dt*L(u^1)]
    
      m_amr->average_down(phi, m_cdr->get_phase());
      m_amr->interp_ghost(phi, m_cdr->get_phase());
    
      data_ops::floor(phi, 0.0);

      // For RK3 Strang splitting, the error has already been computed. For Godunov splitting,
      // the embedded method is u^(n+1) = 0.5*[u^n + u^1 + dt*L(u^1)] and we computed L(u^1) above into rhs
      // so just increment the error with what we need. 
      if(a_compute_err && m_splitting == "godunov"){
	EBAMRCellData& err = storage->get_error(); // err = u^1
	data_ops::incr(err,  rhs, a_dt);           // err = u^1 + dt*L(u^1)
	data_ops::incr(err,  pre, 1.0);            // err = u^n + u^1 + dt*L(u^1)
	data_ops::scale(err, 0.5);                 // err = 0.5*[u^n + u^1 + dt*L(u^1)]
      }
    }

    EBAMRIVData& sigma = m_sigma->get_state();           // u^1
    EBAMRIVData& rhs = m_sigma_scratch->get_scratch();   // Storage for right hand side
    EBAMRIVData& pre = m_sigma_scratch->get_previous();  // u^n
    m_sigma->compute_rhs(rhs);
    data_ops::incr(sigma, rhs, a_dt);  // sigma = u^1 + dt*L(u^1)
    data_ops::incr(sigma, pre, 3.0);   // sigma = 3*u^n + u^1 + dt*L(u^1)
    data_ops::scale(sigma, 0.25);      // sigma = 0.25*[3*u^n + u^1 + dt*L(u^1)]

    if(a_compute_err && m_splitting == "strang"){
      EBAMRIVData& err = m_sigma_scratch->get_error(); // u^1
      data_ops::incr(err, rhs, a_dt);                  // err = u^1 + dt*L(u^1)
      data_ops::incr(err, pre, 1.0);                   // err = u^n + u^1 + dt*L(u^1)
      data_ops::scale(err, 0.5);                       // err = 0.5*[u^n + u^1 + dt*L(u^1)]
    }
  }

  if(m_consistent_E)   this->update_poisson();
  if(m_consistent_rte) this->update_rte(a_time + a_dt);

  // u^3 = u^n + 2*u^2 + 2*dt*L(u^2). Embedded errors have already been computed so this only advances
  // the final stage of the solution. 
  {
    this->compute_cdr_gradients();
    if(m_compute_v) this->compute_cdr_velo(a_time + a_dt);
    if(m_compute_S) this->compute_cdr_sources(a_time + a_dt);
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();

    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    
      EBAMRCellData& rhs       = storage->get_scratch();
      EBAMRCellData& phi       = solver->get_state();     // phi = u^2
      const EBAMRCellData& src = solver->get_source();
      const EBAMRCellData& pre = storage->get_previous(); // pre = u^n
    
      solver->compute_divF(rhs, phi, 0.0, true); // RHS =  Div(u^2*v^2)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^2*v^2)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^2 - Div(u^2*v^2)
      data_ops::incr(phi, rhs, a_dt);            // u^(n+1) = u^2 + dt*L(u^2)
      data_ops::scale(phi, 2.0);                 // u^(n+1) = 2*u^2 + 2*dt*L(u^2)
      data_ops::incr(phi, pre, 1.0);             // u^(n+1) = u^n + 2*u^2 + dt*L(u^2)
      data_ops::scale(phi, 1./3.);               // u^(n+1) = (1/3)*[u^n + 2*u^2 + 2*dt*L(u^2)]
    
      m_amr->average_down(phi, m_cdr->get_phase());
      m_amr->interp_ghost(phi, m_cdr->get_phase());
    
      data_ops::floor(phi, 0.0);

    }

    EBAMRIVData& sigma = m_sigma->get_state();           // u^2
    EBAMRIVData& rhs = m_sigma_scratch->get_scratch();   // Storage for right hand side
    EBAMRIVData& pre = m_sigma_scratch->get_previous();  // u^n
    m_sigma->compute_rhs(rhs);
    data_ops::incr(sigma, rhs, a_dt);  // sigma = u^2 + dt*L(u^2)
    data_ops::scale(sigma, 2.0);       // sigma = 2*u^2 + 2*dt*L(u^2)
    data_ops::incr(sigma, pre, 1.0);   // sigma = u^n + 2*u^2 + 2*dt*L(u^2)
    data_ops::scale(sigma, 1./3.);     // sigma = (1/3)*[u^n + 2*u^2 + 2*dt*L(u^2)]
  }
}

void strang2::compute_cdr_gradients(){
  CH_TIME("time_stepper::compute_cdr_gradients");
  if(m_verbosity > 5){
    pout() << "time_stepper::compute_cdr_gradients" << endl;
  }

  strang2::compute_cdr_gradients(m_cdr->get_states());
}

void strang2::compute_cdr_gradients(const Vector<EBAMRCellData*>& a_states){
  CH_TIME("time_stepper::compute_cdr_gradients");
  if(m_verbosity > 5){
    pout() << "time_stepper::compute_cdr_gradients" << endl;
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

void strang2::compute_errors(){
  CH_TIME("time_stepper::compute_errors");
  if(m_verbosity > 5){
    pout() << "time_stepper::compute_errors" << endl;
  }

  const Real safety = 1.E-20;
  const int comp = 0;
  Real max, min, emax, emin;

  // CDR errors
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    const int which = solver_it.get_solver();
    
    const EBAMRCellData& phi = solver->get_state();
    EBAMRCellData& err       = storage->get_error();

    // So far 'err' contains the embedded formula and 'phi' is the numerical solution
    data_ops::scale(err, -1.0);    // err -> -err
    data_ops::incr(err, phi, 1.0); // err -> (phi-err)

    m_amr->average_down(err, m_cdr->get_phase());
    
    Real Lerr, Lphi;
    data_ops::norm(Lerr, *err[0], m_amr->get_domains()[0], m_error_norm);
    data_ops::norm(Lphi, *phi[0], m_amr->get_domains()[0], m_error_norm);

    // Don't want to divide by zero...
    Lphi = Max(Lphi, safety);
    
    m_cdr_error[which] = Lerr/Lphi;
  }

  // Sigma error. So far, 'err' contains the embedded formula and 'phi' is the numerical solution
  EBAMRIVData& phi = m_sigma->get_state();
  EBAMRIVData& err = m_sigma_scratch->get_error();
  data_ops::scale(err, -1.0);
  data_ops::incr(err, phi, 1.0);
  data_ops::get_max_min_norm(max, min, phi);
  data_ops::get_max_min_norm(emax, emin, err);

  // DOn't want to divide by zero
  max = Max(max, safety);
  m_sigma_error = emax/max;

  // Maximum error
  m_max_error = this->get_max_error();
}

void strang2::compute_dt(Real& a_dt, time_code::which_code& a_timecode){
  CH_TIME("time_stepper::compute_dt");
  if(m_verbosity > 5){
    pout() << "time_stepper::compute_dt" << endl;
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

void strang2::regrid_internals(){
  CH_TIME("strang2::regrid_internals");
  if(m_verbosity > 5){
    pout() << "strang2::regrid_internals" << endl;
  }

  m_cdr_error.resize(m_plaskin->get_num_species());
  
  this->allocate_cdr_storage();
  this->allocate_poisson_storage();
  this->allocate_rte_storage();
  this->allocate_sigma_storage();
}

void strang2::allocate_cdr_storage(){
  const int ncomp       = 1;
  const int num_species = m_plaskin->get_num_species();

  m_cdr_scratch.resize(num_species);
  
  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx] = RefCountedPtr<cdr_storage> (new cdr_storage(m_advective_order, m_amr, m_cdr->get_phase(), ncomp));
    m_cdr_scratch[idx]->allocate_storage();
  }
}

void strang2::allocate_poisson_storage(){
  const int ncomp = 1;
  m_poisson_scratch = RefCountedPtr<poisson_storage> (new poisson_storage(m_advective_order, m_amr, m_cdr->get_phase(), ncomp));
  m_poisson_scratch->allocate_storage();
}

void strang2::allocate_rte_storage(){
  const int ncomp       = 1;
  const int num_photons = m_plaskin->get_num_photons();
  m_rte_scratch.resize(num_photons);
  
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx] = RefCountedPtr<rte_storage> (new rte_storage(m_advective_order, m_amr, m_rte->get_phase(), ncomp));
    m_rte_scratch[idx]->allocate_storage();
  }
}

void strang2::allocate_sigma_storage(){
  const int ncomp = 1;
  m_sigma_scratch = RefCountedPtr<sigma_storage> (new sigma_storage(m_advective_order, m_amr, m_cdr->get_phase(), ncomp));
  m_sigma_scratch->allocate_storage();
}

void strang2::deallocate_internals(){
  CH_TIME("time_stepper::deallocate_internals");
  if(m_verbosity > 5){
    pout() << "time_stepper::deallocate_internals" << endl;
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

void strang2::cache_solutions(){
  CH_TIME("strang2::cache_solutions");
  if(m_verbosity > 5){
    pout() << "strang2::cache_solutions" << endl;
  }
  
  // Cache cdr solutions
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();

    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    EBAMRCellData& cache = storage->get_cache();

    data_ops::copy(cache, solver->get_state());
  }

  {// Cache Poisson solution
    MFAMRCellData& cache = m_poisson_scratch->get_cache();
    data_ops::copy(cache, m_poisson->get_state());
  }

  // Cache RTE solutions
  for (rte_iterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<rte_solver>& solver = solver_it();

    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);
    EBAMRCellData& cache = storage->get_cache();

    data_ops::copy(cache, solver->get_state());
  }

  { // Cache sigma
    EBAMRIVData& cache = m_sigma_scratch->get_cache();
    data_ops::copy(cache, m_sigma->get_state());
  }
}

void strang2::uncache_solutions(){
  CH_TIME("strang2::uncache_solutions");
  if(m_verbosity > 5){
    pout() << "strang2::uncache_solutions" << endl;
  }
  
  // Uncache cdr solutions
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();

    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    const EBAMRCellData& cache = storage->get_cache();

    data_ops::copy(solver->get_state(), cache);
  }

  {// Uncache Poisson solution
    const MFAMRCellData& cache = m_poisson_scratch->get_cache();
    data_ops::copy(m_poisson->get_state(), cache);
  }

  // Uncache RTE solutions
  for (rte_iterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<rte_solver>& solver = solver_it();

    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);
    const EBAMRCellData& cache = storage->get_cache();

    data_ops::copy(solver->get_state(), cache);
  }

  { // Uncache sigma
    const EBAMRIVData& cache = m_sigma_scratch->get_cache();
    data_ops::copy(m_sigma->get_state(), cache);
  }
}

void strang2::compute_E_into_scratch(){
  CH_TIME("strang2::compute_E_into_scratch");
  if(m_verbosity > 5){
    pout() << "strang2::compute_E_into_scratch" << endl;
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

void strang2::compute_cdr_velo(const Real a_time){
  CH_TIME("strang2::compute_cdr_velo");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_velo" << endl;
  }

  this->compute_cdr_velo(m_cdr->get_states(), a_time);
}

void strang2::compute_cdr_velo(const Vector<EBAMRCellData*>& a_states, const Real a_time){
  CH_TIME("strang2::compute_cdr_velo(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_velo(Vector<EBAMRCellData*>, Real)" << endl;
  }

  Vector<EBAMRCellData*> velocities = m_cdr->get_velocities();
  this->compute_cdr_velocities(velocities, a_states, m_poisson_scratch->get_E_cell(), a_time);
}

void strang2::compute_cdr_eb_states(){
  CH_TIME("strang2::compute_cdr_eb_states");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_eb_states" << endl;
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

void strang2::compute_cdr_eb_states(const Vector<EBAMRCellData*>& a_states){
  CH_TIME("strang2::compute_cdr_eb_states(vec)");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_eb_states(vec)" << endl;
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

void strang2::compute_cdr_domain_states(){
  CH_TIME("strang2::compute_cdr_domain_states");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_domain_states" << endl;
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

void strang2::compute_cdr_domain_states(const Vector<EBAMRCellData*>& a_states){
  CH_TIME("strang2::compute_cdr_domain_states");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_domain_states" << endl;
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

void strang2::compute_cdr_fluxes(const Real a_time){
  CH_TIME("strang2::compute_cdr_fluxes");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_fluxes" << endl;
  }

  this->compute_cdr_fluxes(m_cdr->get_states(), a_time);
}

void strang2::compute_cdr_fluxes(const Vector<EBAMRCellData*>& a_states, const Real a_time){
  CH_TIME("strang2::compute_cdr_fluxes(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_fluxes(Vector<EBAMRCellData*>, Real)" << endl;
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

void strang2::compute_cdr_domain_fluxes(const Real a_time){
  CH_TIME("strang2::compute_cdr_domain_fluxes");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_domain_fluxes" << endl;
  }

  this->compute_cdr_domain_fluxes(m_cdr->get_states(), a_time);
}

void strang2::compute_cdr_domain_fluxes(const Vector<EBAMRCellData*>& a_states, const Real a_time){
  CH_TIME("strang2::compute_cdr_domain_fluxes(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_domain_fluxes(Vector<EBAMRCellData*>, Real)" << endl;
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

void strang2::compute_sigma_flux(){
  CH_TIME("strang2::compute_sigma_flux");
  if(m_verbosity > 5){
    pout() << "strang2::compute_sigma_flux" << endl;
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

void strang2::compute_cdr_sources(const Real a_time){
  CH_TIME("strang2::compute_cdr_sources_into_scratch");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_sources_into_scratch" << endl;
  }

  this->compute_cdr_sources(m_cdr->get_states(), a_time);
}

void strang2::compute_cdr_sources(const Vector<EBAMRCellData*>& a_states, const Real a_time){
  CH_TIME("strang2::compute_cdr_sources(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_sources(Vector<EBAMRCellData*>, Real)" << endl;
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

void strang2::advance_rte_stationary(const Real a_time){
  CH_TIME("strang2::compute_rte_k1_stationary");
  if(m_verbosity > 5){
    pout() << "strang2::compute_k1_stationary" << endl;
  }

  if((m_step + 1) % m_fast_rte == 0){
    Vector<EBAMRCellData*> rte_states  = m_rte->get_states();
    Vector<EBAMRCellData*> rte_sources = m_rte->get_sources();
    Vector<EBAMRCellData*> cdr_states  = m_cdr->get_states();

    EBAMRCellData& E = m_poisson_scratch->get_E_cell();

    const Real dummy_dt = 0.0;
    this->solve_rte(rte_states, rte_sources, cdr_states, E, a_time, dummy_dt, centering::cell_center);
  }
}

void strang2::store_solvers(){
  CH_TIME("strang2::store_solvers");
  if(m_verbosity > 5){
    pout() << "strang2::store_solvers" << endl;
  }

  // Copy solver states into storage->m_previous
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& state = solver->get_state();
    EBAMRCellData& prev  = storage->get_previous();

    data_ops::copy(prev, state);
  }

  // Copy solver state into storage->m_previous
  EBAMRIVData& phi = m_sigma->get_state();
  EBAMRIVData& pre = m_sigma_scratch->get_previous();
  data_ops::set_value(pre, 0.0);
  data_ops::incr(pre, phi, 1.0);
}

void strang2::restore_solvers(){
  CH_TIME("strang2::restore_solvers");
  if(m_verbosity > 5){
    pout() << "strang2::restore_solvers" << endl;
  }

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& state = solver->get_state();
    EBAMRCellData& prev  = storage->get_previous();

    data_ops::copy(state, prev);
  }

  EBAMRIVData& phi = m_sigma->get_state();
  EBAMRIVData& pre = m_sigma_scratch->get_previous();
  data_ops::set_value(phi, 0.0);
  data_ops::incr(phi, pre, 1.0);
}

void strang2::update_poisson(){
  if(m_verbosity > 5){
    pout() << "strang2::update_poisson" << endl;
  }
  
  if(m_do_poisson){ // Solve Poisson equation
    if((m_step +1) % m_fast_poisson == 0){
      time_stepper::solve_poisson();
      this->compute_E_into_scratch();
    }
  }
}

void strang2::update_rte(const Real a_time){
  if(m_verbosity > 5){
    pout() << "strang2::update_rte" << endl;
  }
  
  if(m_do_rte){
    this->advance_rte_stationary(a_time);
  }
}
