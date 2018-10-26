/*!
  @file   multirate_rkSSP.cpp
  @brief  Implementation of multirate_rkSSP.H
  @author Robert Marskar
  @date   Sept. 2018
*/

#include "multirate_rkSSP.H"
#include "multirate_rkSSP_storage.H"
#include "data_ops.H"
#include "units.H"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <ParmParse.H>

typedef multirate_rkSSP::cdr_storage     cdr_storage;
typedef multirate_rkSSP::poisson_storage poisson_storage;
typedef multirate_rkSSP::rte_storage     rte_storage;
typedef multirate_rkSSP::sigma_storage   sigma_storage;

multirate_rkSSP::multirate_rkSSP(){
  m_order      = 2;
  m_maxCFL     = 0.5;
  m_err_thresh = 1.E-2;;
  m_safety     = 0.95;
  m_error_norm = 2;

  m_compute_v    = true;
  m_compute_S    = true;
  m_adaptive_dt  = true;
  m_multirate    = true;
  m_consistent   = true;
  m_have_dtf     = false;

  // Basically only for debugging
  m_do_advec_src = true;
  m_do_diffusion = true;
  m_do_rte       = true;
  m_do_poisson   = true;
  m_print_diagno = false;
  m_write_diagno = false;
  
  {
    ParmParse pp("multirate_rkSSP");

    std::string str;

    pp.query("order",      m_order);
    pp.query("max_cfl",    m_maxCFL);
    pp.query("max_error",  m_err_thresh);
    pp.query("safety",     m_safety);
    pp.query("error_norm", m_error_norm);

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
    if(pp.contains("consistent")){
      pp.get("consistent", str);
      if(str == "false"){
	m_consistent = false;
      }
    }
    if(pp.contains("use_multirate")){
      pp.get("use_multirate", str);
      if(str == "false"){
	m_multirate = false;
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
    if(pp.contains("turn_off_advection")){
      pp.get("turn_off_advection_source", str);
      if(str == "true"){
	m_do_advec_src = false;
	if(m_verbosity > 2){
	  pout() << "multirate_rkSSP::multirate_rkSSP - Turning off advection & source" << endl;
	}
      }
    }
    if(pp.contains("turn_off_diffusion")){
      pp.get("turn_off_diffusion", str);
      if(str == "true"){
	m_do_diffusion = false;
	if(m_verbosity > 2){
	  pout() << "multirate_rkSSP::multirate_rkSSP - Turning off diffusion" << endl;
	}
      }
    }
    if(pp.contains("turn_off_rte")){
      pp.get("turn_off_rte", str);
      if(str == "true"){
	m_do_rte = false;

	if(m_verbosity > 2){
	  pout() << "multirate_rkSSP::multirate_rkSSP - Turning off rte" << endl;
	}
      }
    }
    if(pp.contains("turn_off_poisson")){
      pp.get("turn_off_poisson", str);
      if(str == "true"){
	m_do_poisson = false;

	if(m_verbosity > 2){
	  pout() << "multirate_rkSSP::multirate_rkSSP - Turning off poisson" << endl;
	}
      }
    }
  }

  if(m_order < 1 || m_order > 3){
    pout() << "multirate_rkSSP ::multirate_rkSSP - order < 1 or order > 3 requested!." << endl;
    MayDay::Abort("multirate_rkSSP ::multirate_rkSSP - order < 1 or order > 3 requested!.");
  }

  // No embedded schemes lower than first order...
  if(m_order == 1){
    m_adaptive_dt = false;
  }

  // If we don't use multirate, we can't go beyond CFL
  if(!m_multirate){
    m_cfl = m_maxCFL;
  }
}

multirate_rkSSP::~multirate_rkSSP(){

}

RefCountedPtr<cdr_storage>& multirate_rkSSP::get_cdr_storage(const cdr_iterator& a_solverit){
  return m_cdr_scratch[a_solverit.get_solver()];
}

RefCountedPtr<rte_storage>& multirate_rkSSP::get_rte_storage(const rte_iterator& a_solverit){
  return m_rte_scratch[a_solverit.get_solver()];
}

Real multirate_rkSSP::restrict_dt(){
  return 1.E99;
}

Real multirate_rkSSP::get_max_error(){
  CH_TIME("multirate_rkSSP::get_max_error");
  if(m_verbosity > 2){
    pout() << "multirate_rkSSP::get_max_error" << endl;
  }

  Real cur_err = 0.0;//m_sigma_error;
  for (int i = 0; i < m_cdr_error.size(); i++){
    cur_err = Max(cur_err, m_cdr_error[i]);
  }

  return cur_err;
}

Real multirate_rkSSP::advance(const Real a_dt){
  CH_TIME("multirate_rkSSP::advance");
  if(m_verbosity > 2){
    pout() << "multirate_rkSSP::advance" << endl;
  }
  
  // When we enter this routine, we assume that velocities, source terms, and diffusion coefficients have already
  // been computed. 
  const Real t0 = MPI_Wtime();
  this->cache_solutions(); // Cache old solutions. We can safely manipulate directly into solver states.

  int advective_substeps = 0; // Number of advective_substeps
  int diffusive_substeps = 0; // Number of advective_substeps
  
  Real cfl;     // Average or fixed CFL for advective+source term advancements
  Real diff_dt;
  
  // Advection and source term advancements
  if(m_do_advec_src){
    if(m_adaptive_dt){
      if(!m_have_dtf){
	// Estimate dtf by substepping. This is probably a bad time step but it will be adjusted
	advective_substeps   = ceil(a_dt/(m_maxCFL*m_dt_cfl));
	cfl        = a_dt/(advective_substeps*m_dt_cfl);
	m_dt_advec = cfl*m_dt_cfl;
	m_have_dtf = true;
      }

      // Start substepping
      advective_substeps = 0;
      this->advance_multirate_adaptive(advective_substeps, m_dt_advec, m_time, a_dt);
    }
    else{
      advective_substeps   = ceil(a_dt/(m_maxCFL*m_dt_cfl));
      cfl        = a_dt/(advective_substeps*m_dt_cfl);
      m_dt_advec = cfl*m_dt_cfl;
    
      this->advance_multirate_fixed(advective_substeps, m_dt_advec);
    }
  }
  const Real t1 = MPI_Wtime();

#if 0 // Debug
  if(m_adaptive_dt){
    if(m_order >= 2){
      if(procID() == 0){
	std::cout << "Advection:" << "\t cdr_errors = " << m_cdr_error << "\t sigma error = " << m_sigma_error << std::endl;
      }
    }
  }
#endif


  // Diffusion advance
  if(m_do_diffusion){
    diff_dt = a_dt;
    diffusive_substeps = 0;
    this->advance_diffusion(diffusive_substeps, diff_dt, a_dt);
  }
#if 0 // Debug
  if(m_adaptive_dt){
    if(m_order >= 2){
      if(procID() == 0){
	std::cout << "Diffusion:" << "\t cdr_errors = " << m_cdr_error << "\t sigma error = " << m_sigma_error << std::endl;
      }
    }
  }
#endif
  const Real t2 = MPI_Wtime();


  // Solve Poisson equation
  if(m_do_poisson){ // Solve Poisson equation
    if((m_step +1) % m_fast_poisson == 0){
      time_stepper::solve_poisson();
      this->compute_E_into_scratch();
    }
  }
  const Real t3 = MPI_Wtime();

  // Solve RTE equations
  if(m_do_rte){ // Solve for final RTE stage. Poisson equation should already have been solved at the end of the advection
    this->advance_rte_stationary(m_time + a_dt); // and diffusion stages
  }
  const Real t4 = MPI_Wtime();

  
  // Put cdr solvers back in useable state so that we can reliably compute the next time step.
  this->compute_cdr_gradients();
  multirate_rkSSP::compute_cdr_velo(m_time + a_dt);
  time_stepper::compute_cdr_diffusion(m_poisson_scratch->get_E_cell(), m_poisson_scratch->get_E_eb());
  multirate_rkSSP::compute_cdr_sources(m_time + a_dt);
  const Real t5 = MPI_Wtime();

  // Print diagnostics
  if(m_print_diagno){
    pout() << "\n";
    pout() << "\t multirate_rkSSP::advance(Real a_dt) breakdown" << endl;
    pout() << "\t ==============================================\n" << endl;
    
    pout() << "\t Convection-source advance: \n ";
    pout() << "\t --------------------------\n";
    pout() << "\t\t Steps     = " << advective_substeps << endl;
    if(m_adaptive_dt){
      pout() << "\t\t Avg. cfl  = " << a_dt/(advective_substeps*m_dt_cfl) << endl;
      pout() << "\t\t Avg. dt   = " << a_dt/advective_substeps << endl;
    }
    else{
      pout() << "\t\t Local cfl = " << cfl        << endl;
      pout() << "\t\t Local dt  = " << m_dt_advec << endl;
    }
    pout() << endl;
    pout() << "\t\t Time      = " << 100.*(t1-t0)/(t5-t0) << "%\n" << endl;

    pout() << "\n";
    pout() << "\t Diffusion advance:\n ";
    pout() << "\t -----------------\n";
    pout() << "\t\t Steps     = " << diffusive_substeps << endl;
    if(m_adaptive_dt){
      pout() << "\t\t Avg. dt   = " << a_dt/diffusive_substeps << endl;
    }
    else{
      pout() << "\t\t Local dt  = " << a_dt << endl;
    }
    pout() << endl;
    pout() << "\t\t Time      = " << 100.*(t2-t1)/(t5-t0) << "%\n" << endl;

    pout() << "\n";
    pout() << "\t Poisson solve:\n ";
    pout() << "\t --------------\n";
    pout() << "\t\t Time      = " << 100.*(t3-t2)/(t5-t0) << "%\n" << endl;

    pout() << "\n";
    pout() << "\t RTE solve:\n ";
    pout() << "\t ----------\n";
    pout() << "\t\t Time      = " << 100.*(t4-t3)/(t5-t0) << "%\n" << endl;

    pout() << "\n";
    pout() << "\t Fill solvers at end:\n ";
    pout() << "\t --------------------\n";
    pout() << "\t\t Time      = " << 100.*(t5-t4)/(t5-t0) << "%\n" << endl;

    pout() << "\n";
    pout() << "\t Total time = " << t5-t0 << endl;
    pout() << "\t ==============================================\n" << endl;


  }
  if(m_write_diagno){
    this->write_diagnostics(advective_substeps, m_dt_advec, a_dt, cfl, t1-t0, t2-t1, t3-t2, t4-t3, t5-t4, t5-t0);
  }

  return a_dt;
}

void multirate_rkSSP::advance_diffusion(int& a_substeps, Real& a_dt, const Real a_dtc){
  CH_TIME("multirate_rkSSP::advance_diffusion");
  if(m_verbosity > 5){
    pout() << "multirate_rkSSP::advance_diffusion" << endl;
  }

  // Diffusive advance for all cdr equations
  bool diffusive_states = false;
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();

    diffusive_states = solver->is_diffusive() ? true : diffusive_states;
  }

  // Do the diffusion advance
  if(diffusive_states){
    if(m_consistent){ // Consistent => Update E before computing diffusion coefficients
      if(m_do_poisson){ // Solve Poisson equation
	if((m_step +1) % m_fast_poisson == 0){
	  time_stepper::solve_poisson();
	  this->compute_E_into_scratch();
	}
      }
    }

    // Necessary stuff
    m_cdr->set_source(0.0);
    this->compute_cdr_diffusion(m_poisson_scratch->get_E_cell(), m_poisson_scratch->get_E_eb());

    Real sum_dt = 0.0;
    while(sum_dt < a_dtc){

      // Last time step can be smaller
      Real actual_dt = a_dt;
      bool last_step = false;
      if(sum_dt + a_dt >= a_dtc){
	actual_dt = a_dtc - sum_dt;
	last_step = true;
      }

      // Storage solvers
      this->store_solvers();
      
      for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	RefCountedPtr<cdr_solver>& solver   = solver_it();
	RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
      
	EBAMRCellData& phi   = solver->get_state();
	EBAMRCellData& error = storage->get_error();
	if(solver->is_diffusive()){
	  solver->advance_diffusion(phi, error, actual_dt);
	}
	else{
	  data_ops::set_value(error, 0.0);
	}
      }

      // Compute maximum error
      m_sigma_error = 0.0; // No sigma error now
      this->compute_errors();

      const Real new_dt = a_dt*sqrt(m_err_thresh/m_max_error);
      const bool accept = m_max_error <= m_err_thresh;

      if(accept){
	sum_dt += actual_dt;
	a_substeps += 1;

	// Set new time step, but only if we are sufficiently far from the error threshold
	// Note: No special actual is taken if the error was estimated from the last time step
	// This will give a smaller time step for the first step in the next multirate loop,
	// but this should be ok, it will be increased anyways.
	const Real rel_err = m_max_error/m_err_thresh;
	if(rel_err < m_safety){
	  a_dt  = m_safety*new_dt; // Errors are far from the threshold, use a larger time step
	}
	else{ 
	  a_dt  = m_safety*a_dt;   // If errors get too close, reduce the time step a little bit so we don̈́t get rejected steps
	}
      }
      else{
	a_dt = m_safety*new_dt;
	this->restore_solvers();
      }
    }
  }
}

void multirate_rkSSP::advance_multirate_fixed(const int a_substeps, const Real a_dt){
  CH_TIME("multirate_rkSSP::advance_multirate_fixed");
  if(m_verbosity > 2){
    pout() << "multirate_rkSSP::advance_multirate_fixed" << endl;
  }

  this->compute_E_into_scratch();
  this->compute_cdr_gradients(); // Precompute gradients

  for (int step = 0; step < a_substeps; step++){
    const Real time = m_time + step*a_dt;
    const bool last_step = (step == a_substeps - 1);

    // Integrate
    if(m_order == 1){
      this->advance_rk1(time, a_dt);
    }
    else if(m_order == 2){
      this->store_solvers();
      this->advance_rk2(time, a_dt);
    }
    else if(m_order == 3){
      this->store_solvers();
      this->advance_rk3(time, a_dt);
    }
    else{
      MayDay::Abort("multirate_rkSSP::advance_multirate_fixed - unsupported order requested");
    }



    // Update for next iterate. This should generally be done because the solution might move many grid cells.
    if(!last_step){
      if(m_consistent){ // Consistent => Update E and RTE before the next fine step
	if(m_do_poisson){ // Solve Poisson equation
	  if((m_step +1) % m_fast_poisson == 0){
	    time_stepper::solve_poisson();
	    this->compute_E_into_scratch();
	  }
	}
	if(m_do_rte){this->advance_rte_stationary(m_time + a_dt);}
      }

      // Compute new velocities and source terms
      this->compute_cdr_gradients();
      if(m_compute_v) this->compute_cdr_velo(time);
      if(m_compute_S) this->compute_cdr_sources(time);
    }
  }
}

void multirate_rkSSP::advance_multirate_adaptive(int& a_substeps, Real& a_dt, const Real a_time, const Real a_dtc){
  CH_TIME("multirate_rkSSP::advance_multirate_adaptive");
  if(m_verbosity > 2){
    pout() << "multirate_rkSSP::advance_multirate_adaptive" << endl;
  }

  this->compute_E_into_scratch();

  Real sum_dt = 0.0;
  while(sum_dt < a_dtc){

    // Adjust last time step so that we land on a_dtc
    Real actual_dt = a_dt;
    bool last_step = false;
    if(sum_dt + a_dt >= a_dtc){
      actual_dt = a_dtc - sum_dt;
      last_step = true;
    }
    const Real cur_time  = a_time + sum_dt;

    // Store solver states
    this->store_solvers();

    // Advance with errors estimated with embedded formulas
    if(m_order == 2){
      this->advance_rk2(cur_time, actual_dt);
    }
    else if(m_order == 3){
      this->advance_rk3(cur_time, actual_dt);
    }

    // Compute errors and new time step
    this->compute_errors();
    const Real new_dt = a_dt*pow(m_err_thresh/m_max_error, 1.0/m_order);
    
    // Accept or discard time step
    const bool accept = m_max_error < m_err_thresh;
    if(accept){
      sum_dt += a_dt;
      a_substeps += 1;

      // Set new time step, but only if we are sufficiently far from the error threshold
      // Note: No special actual is taken if the error was estimated from the last time step
      // This will give a smaller time step for the first step in the next multirate loop,
      // but this should be ok, it will be increased anyways.
      const Real rel_err = m_max_error/m_err_thresh;
      if(rel_err < m_safety){
	a_dt  = m_safety*new_dt; // Errors are far from the threshold, use a larger time step
      }
      else{ 
	a_dt  = m_safety*a_dt;   // If errors get too close, reduce the time step a little bit so we don̈́t get rejected steps
      }
    }
    else{
      a_dt = m_safety*new_dt;
      this->restore_solvers();
    }

    // New time step should never exceed CFL constraints
    a_dt = Min(a_dt, m_maxCFL*m_dt_cfl);

    // Update for next iterate. This should generally be done because the solution might move many grid cells.
    if(!last_step){

      if(m_consistent){ // Consistent => Update E and RTE before the next fine step
	if(m_do_poisson){ // Solve Poisson equation
	  if((m_step +1) % m_fast_poisson == 0){
	    time_stepper::solve_poisson();
	    this->compute_E_into_scratch();
	  }
	}
	if(m_do_rte){this->advance_rte_stationary(m_time + a_dt);}
      }

      this->compute_cdr_gradients();
      if(m_compute_v) this->compute_cdr_velo(cur_time);
      if(m_compute_S) this->compute_cdr_sources(cur_time);
    } 
  }
}

void multirate_rkSSP::advance_rk1(const Real a_time, const Real a_dt){
  CH_TIME("multirate_rkSSP::advance_rk1");
  if(m_verbosity > 2){
    pout() << "multirate_rkSSP::advance_rk1" << endl;
  }

  this->compute_cdr_domain_states();
  this->compute_cdr_eb_states();
  this->compute_cdr_domain_fluxes(a_time);
  this->compute_cdr_fluxes(a_time);
  this->compute_sigma_flux();

  // Advance advection-reaction
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& phi = solver->get_state();
    EBAMRCellData& src = solver->get_source();
      
    EBAMRCellData& rhs = storage->get_scratch();
    data_ops::set_value(rhs, 0.0);              
    solver->compute_divF(rhs, phi, 0.0, true);  // RHS =  div(u*v)
    data_ops::scale(rhs, -1.0);                 // RHS = -div(u*v)
    data_ops::incr(rhs, src, 1.0);              // RHS = S - div(u*v)
    data_ops::incr(phi, rhs, a_dt);             // u^(n+1) = u^n + (S - div(u*v))*dt

    m_amr->average_down(phi, m_cdr->get_phase());
    m_amr->interp_ghost(phi, m_cdr->get_phase());
  }

  // Advance sigma
  EBAMRIVData& phi = m_sigma->get_state();
  EBAMRIVData& rhs = m_sigma_scratch->get_scratch();

  m_sigma->compute_rhs(rhs);
  m_sigma->reset_cells(rhs);
  data_ops::incr(phi, rhs, a_dt);

  m_amr->average_down(phi, m_cdr->get_phase());
  m_sigma->reset_cells(phi);
}

void multirate_rkSSP::advance_rk2(const Real a_time, const Real a_dt){
  CH_TIME("multirate_rkSSP::advance_rk2");
  if(m_verbosity > 2){
    pout() << "multirate_rkSSP::advance_rk2" << endl;
  }

  // u^1 = u^n + dt*L(u^n)
  // u^n resides in cache, put u^1 in the solver. Use scratch for computing L(u^n)
  {
    // CDR gradients already computed so proceed as usual
    this->compute_cdr_domain_states();
    this->compute_cdr_eb_states();
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_cdr_fluxes(a_time);
    this->compute_sigma_flux();
      
    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& phi  = solver->get_state();
      EBAMRCellData& rhs  = storage->get_scratch();
      EBAMRCellData& src  = solver->get_source();
      EBAMRCellData& pre  = storage->get_previous();
	
      solver->compute_divF(rhs, phi, 0.0, true); // RHS =  div(u*v)
      data_ops::scale(rhs, -1.0);                // RHS = -div(u*v)
      data_ops::incr(rhs, src, 1.0);             // RHS = S - div(u*v)
      data_ops::incr(phi, rhs, a_dt);            // u^(n+1) = u^n + dt*L(u^n)

      m_amr->average_down(phi, m_cdr->get_phase());
      m_amr->interp_ghost(phi, m_cdr->get_phase());

      // Embedded Euler method for error computation
      EBAMRCellData& err = storage->get_error();
      data_ops::copy(err, phi);   // err =  (u^n + dt*L(u^n))
      data_ops::scale(err, -1.0); // err = -(u^n + dt*L(u^n))
    }

    // Advance sigma
    EBAMRIVData& phi = m_sigma->get_state();
    EBAMRIVData& rhs = m_sigma_scratch->get_scratch();
    EBAMRIVData& pre = m_sigma_scratch->get_previous();
      
    m_sigma->compute_rhs(rhs);
    m_sigma->reset_cells(rhs);
    data_ops::incr(phi, rhs, a_dt); // phi = u^1 = u^n + dt*L(u^n)

    m_amr->average_down(phi, m_cdr->get_phase());
    m_sigma->reset_cells(phi);

    // Embedded Euler method
    EBAMRIVData& err = m_sigma_scratch->get_error();
    data_ops::set_value(err, 0.0);
    data_ops::incr(err, phi, -1.0);  // error = -(u^n + dt*L(u^n))
  }

  // For consistent computations, we must do this
  if(m_consistent){ // Must update E and RTE
    if(m_do_poisson){ // Solve Poisson equation
      if((m_step +1) % m_fast_poisson == 0){
	time_stepper::solve_poisson();
	this->compute_E_into_scratch();
      }
    }
    if(m_do_rte){ // Update RTE equations
      this->advance_rte_stationary(m_time + a_dt);
    }
  }

  // u^(n+1) = 0.5*(u^n + u^1 + dt*L(u^(1)))
  // u^n resides in temp storage, u^1 resides in solver. Put u^2 in solver at end and use scratch fo computing L(u^1)
  {
    this->compute_cdr_gradients();
    this->compute_cdr_domain_states();
    this->compute_cdr_eb_states();
    this->compute_cdr_domain_fluxes(a_time + a_dt);
    this->compute_cdr_fluxes(a_time + a_dt);
    this->compute_sigma_flux();

    if(m_compute_v) this->compute_cdr_velo(a_time + a_dt);
    if(m_compute_S) this->compute_cdr_sources(a_time + a_dt);

    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& phi = solver->get_state();     // u^1
      EBAMRCellData& rhs = storage->get_scratch();  // Storage for RHS
      EBAMRCellData& src = solver->get_source();    // Source
      EBAMRCellData& pre = storage->get_previous(); // u^n

      solver->compute_divF(rhs, phi, 0.0, true);   // RHS =  div(u^1*v)
      data_ops::scale(rhs, -1.0);                  // RHS = -div(u^1*v)
      data_ops::incr(rhs, src, 1.0);               // RHS = S - div(u^1*v)
      data_ops::incr(phi, rhs, a_dt);              // phi = u^1 + dt*L(u^1)
      data_ops::incr(phi, pre, 1.0);               // phi = u^n + u^1 + dt*L(u^1)
      data_ops::scale(phi, 0.5);                   // phi = 0.5*(u^n + u^1 + dt*L(u^1))

      m_amr->average_down(phi, m_cdr->get_phase());
      m_amr->interp_ghost(phi, m_cdr->get_phase());

      // Compute cdr error
      EBAMRCellData& err = storage->get_error(); // err  = -(u^n + dt*L(u^n)) // 1st. order
      data_ops::incr(err, phi, 1.0);             // err += (second order)
    }

    // Advance sigma
    EBAMRIVData& phi = m_sigma->get_state();
    EBAMRIVData& rhs = m_sigma_scratch->get_scratch();
    EBAMRIVData& pre = m_sigma_scratch->get_previous();

    m_sigma->compute_rhs(rhs);
    m_sigma->reset_cells(rhs);
    data_ops::incr(phi, rhs, a_dt); // phi = u^1 + dt*L(u^1)
    data_ops::incr(phi, pre, 1.0);  // phi = u^n + u^1 + dt*L(u^1)
    data_ops::scale(phi, 0.5);      // phi = 0.5*(u^n + u^1 + dt*L(u^1))

    m_amr->average_down(phi, m_cdr->get_phase());
    m_sigma->reset_cells(phi);

    // Compute error sigma error
    EBAMRIVData& err = m_sigma_scratch->get_error(); // err = -(u^n + dt*L(u^n))
    data_ops::incr(err, phi, 1.0);                   // err = phi - (u^n + dt*L(u^n))
    m_sigma->reset_cells(err);
  }
}

void multirate_rkSSP::advance_rk3(const Real a_time, const Real a_dt){
  CH_TIME("multirate_rkSSP::advance_rk3");
  if(m_verbosity > 2){
    pout() << "multirate_rkSSP::advance_rk3" << endl;
  }

  // u^1 = u^n + dt*L(u^n)
  // u^n resides in solver. Use this and put u^1 in the solver instead. Use scratch for computing L(u^n)
  {
    this->compute_cdr_domain_states();
    this->compute_cdr_eb_states();
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_cdr_fluxes(a_time);
    this->compute_sigma_flux();

    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& phi        = solver->get_state();
      EBAMRCellData& rhs        = storage->get_scratch();
      EBAMRCellData& pre        = storage->get_previous();
      const EBAMRCellData& src  = solver->get_source();
	
      solver->compute_divF(rhs, phi, 0.0, true); // RHS =  div(u*v)
      data_ops::scale(rhs, -1.0);                // RHS = -div(u*v)
      data_ops::incr(rhs, src, 1.0);             // RHS = S - div(u*v)
      data_ops::incr(phi, rhs, a_dt);            // u^1 = u^n + dt*L(u^n)

      m_amr->average_down(phi, m_cdr->get_phase());
      m_amr->interp_ghost(phi, m_cdr->get_phase());

      // Partial embedded formula
      EBAMRCellData& err = storage->get_error();
      data_ops::copy(err, pre);      // err = u^n
      data_ops::incr(err, phi, 1.0); // err = u^n + u^1
    }

    // Advance sigma
    EBAMRIVData& phi = m_sigma->get_state();
    EBAMRIVData& rhs = m_sigma_scratch->get_scratch();
    EBAMRIVData& pre = m_sigma_scratch->get_previous();
      
    m_sigma->compute_rhs(rhs);
    m_sigma->reset_cells(rhs);
    data_ops::incr(phi, rhs, a_dt); // phi = u^n + dt*L(u^n)

    m_amr->average_down(phi, m_cdr->get_phase());
    m_sigma->reset_cells(phi);

    // Partial embedded formula for 
    EBAMRIVData& err = m_sigma_scratch->get_error();
    data_ops::set_value(err, 0.0);
    data_ops::incr(err, pre, 1.0); // err = u^n
    data_ops::incr(err, phi, 1.0); // err = u^n + u^1
  }

  // u^2 = (3*u^n + u^1 + dt*L(u^(1)))/4
  // u^n is stored, u^1 resides in solver. Put u^2 in solver at end and use scratch fo computing L(u^1)
  {
    this->compute_cdr_gradients();
    this->compute_cdr_domain_states();
    this->compute_cdr_eb_states();
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_cdr_fluxes(a_time);
    this->compute_sigma_flux();

    if(m_compute_v) this->compute_cdr_velo(a_time);
    if(m_compute_S) this->compute_cdr_sources(a_time);

    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& phi       = solver->get_state();     // u^1
      EBAMRCellData& rhs       = storage->get_scratch();  // Storage for RHS
      const EBAMRCellData& src = solver->get_source();    // Source
      const EBAMRCellData& pre = storage->get_previous(); // u^n

      solver->compute_divF(rhs, phi, 0.0, true);   // RHS =  div(u^1*v^1)
      data_ops::scale(rhs, -1.0);                  // RHS = -div(u^1*v^1)
      data_ops::incr(rhs, src, 1.0);               // RHS = S^1 - div(u^1*v^1)

      data_ops::incr(phi, pre, 3);                 // u^2 = 3*u^n + u^1 
      data_ops::incr(phi, rhs, a_dt);              // u^2 = 3^u^n + u^1 + dt*L(u^1)
      data_ops::scale(phi, 0.25);                  // phi = (3*u^n + u^1 + dt*L(u^1))/4

      m_amr->average_down(phi, m_cdr->get_phase());
      m_amr->interp_ghost(phi, m_cdr->get_phase());

      // Embedded formula for error
      EBAMRCellData& err = storage->get_error(); // err =  (u^n + u^1)
      data_ops::incr(err, rhs, a_dt);            // err =  (u^n + u^1 + dt*L(u^1))
      data_ops::scale(err, -0.5);                // err = -(u^n + u^1 + dt*L(u^1))/2
    }

    // Advance sigma
    EBAMRIVData& phi = m_sigma->get_state();
    EBAMRIVData& rhs = m_sigma_scratch->get_scratch();
    EBAMRIVData& pre = m_sigma_scratch->get_previous();

    m_sigma->compute_rhs(rhs);
    m_sigma->reset_cells(rhs);
    data_ops::incr(phi, rhs, a_dt); // phi = u^1 + dt*L(u^1)
    data_ops::incr(phi, pre, 3.0);  // phi = 3*u^n + u^1 + dt*L(u^1)
    data_ops::scale(phi, 0.25);     // phi = (3*u^n + u^1 + dt*L(u^1))/4

    m_amr->average_down(phi, m_cdr->get_phase());
    m_sigma->reset_cells(phi);

    // Embedded formula for error
    EBAMRIVData& err = m_sigma_scratch->get_error(); // err =  (u^n + u^1)
    data_ops::incr(err, rhs, a_dt);                  // err =  (u^n + u^1 + dt*L(u^1))
    data_ops::scale(err, -0.5);                      // err = -(u^n + u^1 + dt*L(u^1))/2
  }

  // u^(n+1) = (u^n + 2*u^2 + 2*dt*L(u^(1)))/3
  // u^n is stored, u^2 resides in solver. Put u^(n+1) in solver at end and use scratch for computing L(u^1)
  {
    this->compute_cdr_gradients();
    this->compute_cdr_domain_states();
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_sigma_flux();

    if(m_compute_v) this->compute_cdr_velo(a_time);
    if(m_compute_S) this->compute_cdr_sources(a_time);

    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& phi       = solver->get_state();     // u^(2)
      EBAMRCellData& rhs       = storage->get_scratch();  // Storage for RHS
      const EBAMRCellData& src = solver->get_source();    // Source
      const EBAMRCellData& pre = storage->get_previous(); // u^n

      solver->compute_divF(rhs, phi, 0.0, true);   // RHS =  div(u*v)
      data_ops::scale(rhs, -1.0);                  // RHS = -div(u*v)
      data_ops::incr(rhs, src, 1.0);               // RHS = S - div(u*v)

      data_ops::scale(phi, 2.0);                   // u^(n+1) = 2*u^2
      data_ops::incr(phi, pre, 1);                 // u^(n+1) = u^n + 2*u^2
      data_ops::incr(phi, rhs, 2.0*a_dt);          // u^(n+1) = u^n + 2*u^2 + 2*dt*L(u^2)
      data_ops::scale(phi, 1./3.);                 // u^(n+1) = (u^n + 2*u^2 + 2*dt*L(u^2))/3

      m_amr->average_down(phi, m_cdr->get_phase());
      m_amr->interp_ghost(phi, m_cdr->get_phase());

      // Error computation
      EBAMRCellData& err = storage->get_error(); // err  = -0.5*(u^n + u^1 + dt*L(u^1)) // second order
      data_ops::incr(err, phi, 1.0);             // err += (third order approximation)
    }

    // Advance sigma
    EBAMRIVData& phi = m_sigma->get_state();
    EBAMRIVData& rhs = m_sigma_scratch->get_scratch();
    EBAMRIVData& pre = m_sigma_scratch->get_previous();

    m_sigma->compute_rhs(rhs);
    m_sigma->reset_cells(rhs);
    data_ops::incr(phi, rhs, a_dt); // phi = u^2 + dt*L(u^2)
    data_ops::scale(phi, 2.0);      // phi = 2*u^2 + 2*dt*L(u^2)
    data_ops::incr(phi, pre, 1.0);  // phi = u^n + 2*u^2 + 2*dt*L(u^2)
    data_ops::scale(phi, 1./3.);    // phi = (u^n + 2*u^2 + 2*dt*L(u^2))/3

    m_amr->average_down(phi, m_cdr->get_phase());
    m_sigma->reset_cells(phi);

    // Error computation
    EBAMRIVData& err = m_sigma_scratch->get_error(); // err = -(u^n + dt*L(u^n))
    data_ops::incr(err, phi, 1.0);                   // err = phi - (u^n + dt*L(u^n))
  }
}

void multirate_rkSSP::regrid_internals(){
  CH_TIME("time_stepper::regrid_internals");
  if(m_verbosity > 5){
    pout() << "time_stepper::regrid_internals" << endl;
  }

  m_cdr_error.resize(m_plaskin->get_num_species());
  
  this->allocate_cdr_storage();
  this->allocate_poisson_storage();
  this->allocate_rte_storage();
  this->allocate_sigma_storage();
}

void multirate_rkSSP::compute_cdr_gradients(){
  CH_TIME("time_stepper::compute_cdr_gradients");
  if(m_verbosity > 5){
    pout() << "time_stepper::compute_cdr_gradients" << endl;
  }

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    
    const EBAMRCellData& phi = solver->get_state();
    EBAMRCellData& grad = storage->get_gradient();

    m_amr->compute_gradient(grad, phi);
    m_amr->average_down(grad, m_cdr->get_phase());
    m_amr->interp_ghost(grad, m_cdr->get_phase());
  }
}

void multirate_rkSSP::compute_errors(){
  CH_TIME("time_stepper::compute_errors");
  if(m_verbosity > 5){
    pout() << "time_stepper::compute_errors" << endl;
  }

  const int comp = 0;
  Real max, min, emax, emin;

  // CDR errors
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    const int which = solver_it.get_solver();
    
    const EBAMRCellData& phi = solver->get_state();
    EBAMRCellData& err = storage->get_error();

    m_amr->average_down(err, m_cdr->get_phase());


    // const int comp = 0;

    // data_ops::get_max_min(max, min, phi, 0);
    // data_ops::get_max_min_norm(emax, emin, err);
    // m_cdr_error[which] = emax/max;;
    Real Lerr, Lphi;
    data_ops::norm(Lerr, *err[0], m_amr->get_domains()[0], m_error_norm);
    data_ops::norm(Lphi, *phi[0], m_amr->get_domains()[0], m_error_norm);

    m_cdr_error[which] = Lerr/Lphi;
  }

  // Sigma error
  EBAMRIVData& phi = m_sigma->get_state();
  EBAMRIVData& err = m_sigma_scratch->get_error();
  data_ops::get_max_min_norm(max, min, phi);
  data_ops::get_max_min_norm(emax, emin, err);
  m_sigma_error = emax/max;

  // Maximum error
  m_max_error = this->get_max_error();
}

void multirate_rkSSP::compute_dt(Real& a_dt, time_code::which_code& a_timecode){
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

void multirate_rkSSP::allocate_cdr_storage(){
  const int ncomp       = 1;
  const int num_species = m_plaskin->get_num_species();
  m_cdr_scratch.resize(num_species);
  
  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx] = RefCountedPtr<cdr_storage> (new cdr_storage(m_order, m_amr, m_cdr->get_phase(), ncomp));
    m_cdr_scratch[idx]->allocate_storage();
  }
}

void multirate_rkSSP::allocate_poisson_storage(){
  const int ncomp = 1;
  m_poisson_scratch = RefCountedPtr<poisson_storage> (new poisson_storage(m_order, m_amr, m_cdr->get_phase(), ncomp));
  m_poisson_scratch->allocate_storage();
}

void multirate_rkSSP::allocate_rte_storage(){
  const int ncomp       = 1;
  const int num_photons = m_plaskin->get_num_photons();
  m_rte_scratch.resize(num_photons);
  
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx] = RefCountedPtr<rte_storage> (new rte_storage(m_order, m_amr, m_rte->get_phase(), ncomp));
    m_rte_scratch[idx]->allocate_storage();
  }
}

void multirate_rkSSP::allocate_sigma_storage(){
  const int ncomp = 1;
  m_sigma_scratch = RefCountedPtr<sigma_storage> (new sigma_storage(m_order, m_amr, m_cdr->get_phase(), ncomp));
  m_sigma_scratch->allocate_storage();
}

void multirate_rkSSP::deallocate_internals(){
  CH_TIME("time_stepper::deallocate_internals");
  if(m_verbosity > 5){
    pout() << "time_stepper::deallocate_internals" << endl;
  }
  
  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx]->deallocate_storage();
  }

  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx]->deallocate_storage();
  }

  m_poisson_scratch->deallocate_storage();
  m_sigma_scratch->deallocate_storage();
}

void multirate_rkSSP::cache_solutions(){
  CH_TIME("multirate_rkSSP::cache_solutions");
  if(m_verbosity > 5){
    pout() << "multirate_rkSSP::cache_solutions" << endl;
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

void multirate_rkSSP::compute_E_into_scratch(){
  CH_TIME("multirate_rkSSP::compute_E_into_scratch");
  if(m_verbosity > 5){
    pout() << "multirate_rkSSP::compute_E_into_scratch" << endl;
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

void multirate_rkSSP::compute_cdr_velo(const Real a_time){
  CH_TIME("multirate_rkSSP::compute_cdr_velo");
  if(m_verbosity > 5){
    pout() << "multirate_rkSSP::compute_cdr_velo" << endl;
  }

  this->compute_cdr_velo(m_cdr->get_states(), a_time);
}

void multirate_rkSSP::compute_cdr_velo(const Vector<EBAMRCellData*>& a_states, const Real a_time){
  CH_TIME("multirate_rkSSP::compute_cdr_velo(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "multirate_rkSSP::compute_cdr_velo(Vector<EBAMRCellData*>, Real)" << endl;
  }

  Vector<EBAMRCellData*> velocities = m_cdr->get_velocities();
  this->compute_cdr_velocities(velocities, a_states, m_poisson_scratch->get_E_cell(), a_time);
}

void multirate_rkSSP::compute_cdr_eb_states(){
  CH_TIME("multirate_rkSSP::compute_cdr_eb_states");
  if(m_verbosity > 5){
    pout() << "multirate_rkSSP::compute_cdr_eb_states" << endl;
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
    cdr_gradients.push_back(&(storage->get_gradient()));
  }

  // Extrapolate states to the EB
  this->extrapolate_to_eb(eb_states, m_cdr->get_phase(), cdr_states);

  // We already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIVData eb_gradient;
  m_amr->allocate(eb_gradient, m_cdr->get_phase(), SpaceDim);
  for (int i = 0; i < cdr_states.size(); i++){
    this->extrapolate_to_eb(eb_gradient, m_cdr->get_phase(), *cdr_gradients[i]);
    this->project_flux(*eb_gradients[i], eb_gradient);
  }
}

void multirate_rkSSP::compute_cdr_domain_states(){
  CH_TIME("multirate_rkSSP::compute_cdr_domain_states");
  if(m_verbosity > 5){
    pout() << "multirate_rkSSP::compute_cdr_domain_states" << endl;
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
    cdr_gradients.push_back(&(storage->get_gradient()));
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

void multirate_rkSSP::compute_cdr_fluxes(const Real a_time){
  CH_TIME("multirate_rkSSP::compute_cdr_fluxes");
  if(m_verbosity > 5){
    pout() << "multirate_rkSSP::compute_cdr_fluxes" << endl;
  }

  this->compute_cdr_fluxes(m_cdr->get_states(), a_time);
}

void multirate_rkSSP::compute_cdr_fluxes(const Vector<EBAMRCellData*>& a_states, const Real a_time){
  CH_TIME("multirate_rkSSP::compute_cdr_fluxes(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "multirate_rkSSP::compute_cdr_fluxes(Vector<EBAMRCellData*>, Real)" << endl;
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
    extrap_cdr_velocities.push_back(&velo_eb);
    extrap_cdr_fluxes.push_back(&flux_eb);
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

void multirate_rkSSP::compute_cdr_domain_fluxes(const Real a_time){
  CH_TIME("multirate_rkSSP::compute_cdr_domain_fluxes");
  if(m_verbosity > 5){
    pout() << "multirate_rkSSP::compute_cdr_domain_fluxes" << endl;
  }

  this->compute_cdr_domain_fluxes(m_cdr->get_states(), a_time);
}

void multirate_rkSSP::compute_cdr_domain_fluxes(const Vector<EBAMRCellData*>& a_states, const Real a_time){
  CH_TIME("multirate_rkSSP::compute_cdr_domain_fluxes(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "multirate_rkSSP::compute_cdr_domain_fluxes(Vector<EBAMRCellData*>, Real)" << endl;
  }

  Vector<EBAMRIFData*> cdr_fluxes;
  Vector<EBAMRIFData*> extrap_cdr_fluxes;
  Vector<EBAMRIFData*> extrap_cdr_densities;
  Vector<EBAMRIFData*> extrap_cdr_velocities;
  Vector<EBAMRIFData*> extrap_cdr_gradients;
  Vector<EBAMRIFData*> extrap_rte_fluxes;


  cdr_fluxes = m_cdr->get_domainflux();


  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRIFData& dens_domain = storage->get_domain_state();
    EBAMRIFData& velo_domain = storage->get_domain_velo();
    EBAMRIFData& flux_domain = storage->get_domain_flux();
    EBAMRIFData& grad_domain = storage->get_domain_grad();

    extrap_cdr_densities.push_back(&dens_domain);  // Computed in compute_cdr_eb_states
    extrap_cdr_velocities.push_back(&velo_domain); // Has not been computed
    extrap_cdr_fluxes.push_back(&flux_domain);     // Has not been computed
    extrap_cdr_gradients.push_back(&grad_domain);  // Computed in compute_cdr_eb_states
  }

  // Compute extrapolated velocities and fluxes at the domain faces
  Vector<EBAMRCellData*> cdr_velocities = m_cdr->get_velocities();
  this->compute_extrapolated_domain_fluxes(extrap_cdr_fluxes, a_states, cdr_velocities, m_cdr->get_phase());
  this->extrapolate_to_domain_faces(extrap_cdr_velocities, m_cdr->get_phase(), a_states);

  // Compute RTE flux on domain faces
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);

    EBAMRIFData& domain_flux = storage->get_domain_flux();
    solver->compute_domain_flux(domain_flux, solver->get_state());
    extrap_rte_fluxes.push_back(&domain_flux);
  }

  const EBAMRIFData& E = m_poisson_scratch->get_E_domain();
  
  time_stepper::compute_cdr_domain_fluxes(cdr_fluxes,
				   extrap_cdr_fluxes,
				   extrap_cdr_densities,
				   extrap_cdr_velocities,
				   extrap_cdr_gradients,
				   extrap_rte_fluxes,
				   E,
				   a_time);
  
#if 0
  // Extrapolate densities, velocities, and fluxes

  this->compute_extrapolated_fluxes(extrap_cdr_fluxes, a_states, cdr_velocities, m_cdr->get_phase());
  this->extrapolate_to_eb(extrap_cdr_velocities, m_cdr->get_phase(), cdr_velocities);

  // Compute RTE flux on the boundary
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);

    EBAMRIVData& flux_eb = storage->get_eb_flux();
    solver->compute_domain_flux(flux_eb, solver->get_state());
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
#endif
}

void multirate_rkSSP::compute_sigma_flux(){
  CH_TIME("multirate_rkSSP::compute_sigma_flux");
  if(m_verbosity > 5){
    pout() << "multirate_rkSSP::compute_sigma_flux" << endl;
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

void multirate_rkSSP::compute_cdr_sources(const Real a_time){
  CH_TIME("multirate_rkSSP::compute_cdr_sources_into_scratch");
  if(m_verbosity > 5){
    pout() << "multirate_rkSSP::compute_cdr_sources_into_scratch" << endl;
  }

  this->compute_cdr_sources(m_cdr->get_states(), a_time);
}

void multirate_rkSSP::compute_cdr_sources(const Vector<EBAMRCellData*>& a_states, const Real a_time){
  CH_TIME("multirate_rkSSP::compute_cdr_sources(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "multirate_rkSSP::compute_cdr_sources(Vector<EBAMRCellData*>, Real)" << endl;
  }
  
  Vector<EBAMRCellData*> cdr_sources = m_cdr->get_sources();
  Vector<EBAMRCellData*> rte_states  = m_rte->get_states();
  EBAMRCellData& E                   = m_poisson_scratch->get_E_cell();

  Vector<EBAMRCellData*> cdr_gradients;
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    cdr_gradients.push_back(&(storage->get_gradient()));
  }

  time_stepper::compute_cdr_sources(cdr_sources, a_states, cdr_gradients, rte_states, E, a_time, centering::cell_center);
}

void multirate_rkSSP::advance_rte_stationary(const Real a_time){
  CH_TIME("multirate_rkSSP::compute_rte_k1_stationary");
  if(m_verbosity > 5){
    pout() << "multirate_rkSSP::compute_k1_stationary" << endl;
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

void multirate_rkSSP::write_diagnostics(const int  a_substeps,
					const Real a_sub_dt,
					const Real a_glob_dt,
					const Real a_sub_cfl,
					const Real a_convection_time,
					const Real a_diffusion_time,
					const Real a_poisson_time,
					const Real a_rte_time,
					const Real a_misc_time,
					const Real a_total_time){

  // Compute the number of cells in the amr hierarchy
  long long num_cells = 0;
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    num_cells += (m_amr->get_grids()[lvl]).numCells();
  }

  long long uniform_points = (m_amr->get_domains()[m_amr->get_finest_level()]).domainBox().numPts();
  const Real compression = 1.0*num_cells/uniform_points;
  
  if(procID() == 0 ){

    const std::string fname("multirate_rkSSP_diagnostics.txt");
    
    bool write_header;
    { // Write header if we must
      std::ifstream infile(fname);
      write_header = infile.peek() == std::ifstream::traits_type::eof() ? true : false;
    }

    // Write output
    std::ofstream f;
    f.open("multirate_rkSSP_diagnostics.txt", std::ios_base::app);
    const int width = 12;


    if(write_header){
      f << std::left << std::setw(width) << "# Step" << "\t"
	<< std::left << std::setw(width) << "Time" << "\t"
	<< std::left << std::setw(width) << "Time code" << "\t"
	<< std::left << std::setw(width) << "Substeps" << "\t"
	<< std::left << std::setw(width) << "Global dt" << "\t"
	<< std::left << std::setw(width) << "Local dt" << "\t"
	<< std::left << std::setw(width) << "Local cfl" << "\t"
	<< std::left << std::setw(width) << "Convection time" << "\t"
	<< std::left << std::setw(width) << "Diffusion time" << "\t"
	<< std::left << std::setw(width) << "Poisson time" << "\t"
	<< std::left << std::setw(width) << "RTE time" << "\t"
	<< std::left << std::setw(width) << "Misc time" << "\t"
	<< std::left << std::setw(width) << "Total time" << "\t"
	<< std::left << std::setw(width) << "Mesh cells" << "\t"
	<< std::left << std::setw(width) << "Mesh compression" << "\t"
	<< endl;
    }

    f << std::left << std::setw(width) << m_step << "\t"
      << std::left << std::setw(width) << m_time << "\t"            
      << std::left << std::setw(width) << m_timecode << "\t"        
      << std::left << std::setw(width) << a_substeps << "\t"        
      << std::left << std::setw(width) << a_glob_dt << "\t"         
      << std::left << std::setw(width) << a_sub_dt << "\t"
      << std::left << std::setw(width) << a_sub_cfl << "\t"          
      << std::left << std::setw(width) << a_convection_time << "\t" 
      << std::left << std::setw(width) << a_diffusion_time  << "\t"
      << std::left << std::setw(width) << a_poisson_time  << "\t"    
      << std::left << std::setw(width) << a_rte_time  << "\t"        
      << std::left << std::setw(width) << a_misc_time  << "\t"       
      << std::left << std::setw(width) << a_total_time  << "\t"
      << std::left << std::setw(width) << num_cells  << "\t"      
      << std::left << std::setw(width) << compression  << "\t"      
      << std::left << std::setw(width) << endl;
    
  }
}

void multirate_rkSSP::store_solvers(){
  CH_TIME("multirate_rkSSP::store_solvers");
  if(m_verbosity > 5){
    pout() << "multirate_rkSSP::store_solvers" << endl;
  }

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& state = solver->get_state();
    EBAMRCellData& prev  = storage->get_previous();

    data_ops::copy(prev, state);
  }

  EBAMRIVData& phi = m_sigma->get_state();
  EBAMRIVData& pre = m_sigma_scratch->get_previous();
  data_ops::set_value(pre, 0.0);
  data_ops::incr(pre, phi, 1.0);
}

void multirate_rkSSP::restore_solvers(){
  CH_TIME("multirate_rkSSP::restore_solvers");
  if(m_verbosity > 5){
    pout() << "multirate_rkSSP::restore_solvers" << endl;
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
