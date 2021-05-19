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
#include <CD_CdrTGA.H>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <ParmParse.H>

typedef strang2::cdr_storage     cdr_storage;
typedef strang2::poisson_storage poisson_storage;
typedef strang2::rte_storage     rte_storage;
typedef strang2::sigma_storage   sigma_storage;

strang2::strang2(){
  m_rk_order     = 3;
  m_rk_stages    = 3;
  m_minCFL       = 0.1;
  m_maxCFL       = 1.8;   
  m_err_thresh   = 1.E-4;
  m_safety       = 0.9;
  m_error_norm   = 2;
  m_alpha        = 0.25;
  m_min_alpha    = 0.2;
  m_max_alpha    = 1.5;
  
  m_auto_safety_cfl = 0.8;
  m_auto_min_stages = 2;
  m_auto_max_stages = 5;

  m_auto_stages    = false;
  m_adaptive_dt    = true;
  m_compute_error  = true;
  m_have_dtf       = false;
  m_fixed_order    = false;


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
    Vector<Real> rk_method(2);

    if(pp.contains("rk_method")){
      pp.getarr("rk_method", rk_method, 0, 2);
      m_rk_stages = rk_method[0];
      m_rk_order  = rk_method[1];
    }
    pp.query("min_cfl",         m_minCFL);
    pp.query("max_cfl",         m_maxCFL);
    pp.query("max_error",       m_err_thresh);
    pp.query("safety",          m_safety);
    pp.query("error_norm",      m_error_norm);
    pp.query("alpha",           m_alpha);
    pp.query("min_alpha",       m_min_alpha);
    pp.query("max_alpha",       m_max_alpha);
    pp.query("auto_max_stages", m_auto_max_stages);
    pp.query("auto_min_stages", m_auto_min_stages);
    pp.query("auto_safety_cfl", m_auto_safety_cfl);

    if(pp.contains("fixed_order")){
      pp.get("fixed_order", str);
      if(str == "true"){
	m_fixed_order = true;
      }
    }
    if(pp.contains("auto_stages")){
      pp.get("auto_stages", str);
      if(str == "true"){
	m_auto_stages = true;
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

  // If we don't use adaptive, we can't go beyond CFL
  if(!m_adaptive_dt){
    m_cfl = m_maxCFL;
  }

  // If we're doing adaptive stepping, errors must always be comptued
  if(m_adaptive_dt){
    m_compute_error = true;
  }

  // If we do automatic selection of # of stages, order must be do
  if(m_auto_stages){
    m_rk_order  = 2;
    m_rk_stages = 2;
  }
}

strang2::~strang2(){

}

RefCountedPtr<cdr_storage>& strang2::get_cdr_storage(const CdrIterator& a_solverit){
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
#if 0 // Debugging
  this->compute_E_into_scratch();
  this->compute_cdr_gradients();
  this->compute_cdr_velo(m_time);
  this->compute_cdr_sources(m_time);
  TimeStepper::compute_cdr_diffusion(m_fieldSolver_scratch->get_E_cell(), m_fieldSolver_scratch->get_E_eb());
#endif

  this->backup_solutions(); // Store old solution. This stores the old solutions in storage->m_backup
  
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
      // Do equal division of the time steps. This integration always lands on a_dt.

      // Automatic stage selection for second order
      if(m_auto_stages){
	const Real safe_cfl = m_dt_cfl;
	const int s = 1 + ceil(a_dt/(2*m_dt_cfl)); // Because advection is done with a_dt/2. 
	m_rk_stages = s;
	m_rk_stages = Min(m_rk_stages, m_auto_max_stages);
	m_rk_stages = Max(m_rk_stages, m_auto_min_stages);

#if 1 // Debug
	if(procID() == 0){
	  std::cout << m_rk_stages << std::endl;
	}
#endif

	// Given order s. Use more time steps if we must. 
	substeps = ceil(a_dt/(2*(m_rk_stages-1)*m_auto_safety_cfl*m_dt_cfl));
      }
      else{
	substeps = ceil(a_dt/(m_maxCFL*m_dt_cfl)); // Just do equal division
      }

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
  TimeStepper::compute_cdr_diffusion(m_fieldSolver_scratch->get_E_cell(), m_fieldSolver_scratch->get_E_eb());

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
    this->advanceTGA_diffusion(a_time, a_dt);   // This is the 2nd order update
    //    this->advanceEuler_diffusion(a_time, a_dt); // This is the embedded formula 1st order update, it uses the solution from
    // advanceTGA_diffusion as initial guess in order to optimize.
  }
}

void strang2::advanceTGA_diffusion(const Real a_time, const Real a_dt){
  CH_TIME("strang2::advanceTGA_diffusion");
  if(m_verbosity > 2){
    pout() << "strang2::advanceTGA_diffusion" << endl;
  }

  const int ncomp = 1;

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
      
    EBAMRCellData& phi = solver->getPhi();
    if(solver->isDiffusive()){
      EBAMRCellData old_phi;
      m_amr->allocate(old_phi, phase::gas, ncomp);
      data_ops::copy(old_phi, phi);
      
      CdrTGA* tgasolver = (CdrTGA*) (&(*solver));
      tgasolver->advanceTGA(phi, old_phi, a_dt);
    }
  }
}

void strang2::advanceEuler_diffusion(const Real a_time, const Real a_dt){
  CH_TIME("strang2::advanceEuler_diffusion");
  if(m_verbosity > 2){
    pout() << "strang2::advanceEuler_diffusion" << endl;
  }

  const int ncomp = 1;

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
      
    EBAMRCellData& exact_phi = solver->getPhi();
    EBAMRCellData& wrong_phi = storage->get_error();
    if(solver->isDiffusive()){

      // We advance the error but take the 2nd order solution as the initial guess
      // since we anticipate that it will be a "close" solution.
      EBAMRCellData old_phi;
      m_amr->allocate(old_phi, phase::gas, ncomp);
      data_ops::copy(old_phi, wrong_phi);
      data_ops::copy(wrong_phi, exact_phi);

      CdrTGA* tgasolver = (CdrTGA*) (&(*solver));
      tgasolver->advanceTGA(wrong_phi, old_phi, a_dt);
    }
  }
}

Real strang2::advance_one_step(const Real a_time, const Real a_dt){
  CH_TIME("strang2::advance_one_step");
  if(m_verbosity > 2){
    pout() << "strang2::advance_one_step" << endl;
  }

  // Advection-reaction integration. We must have a place to put u^n (the current solution) before
  // we advance, so it is put in storage->m_previous. We do this because it's easier to fill
  // the solvers with the intermediate Runge-Kutta stage (it lets us use a simpler interface).
  // The 'true' flag indicates that we should compute the embedded formula error, which we only do
  // on the first time step (the embedded splitting uses a Lie splitting)
  this->store_solvers();
  this->advance_rk(a_time, 0.5*a_dt);

  // Update elliptic equations. 
  if(m_consistent_E) this->update_poisson();
  if(m_consistent_rte) this->update_rte(a_time + a_dt);

  // Embedded formula is the same for the first half step, but it should use another advective half step
  // followed by diffusion. The embedded formulas will mess with the elliptic solutions, so we should back
  // up those solutions first, and then revert later. 
  if(m_compute_error){
    this->copy_solvers_to_cache();                       // scratch -> e^(0.5*A*dt)u. Also does elliptic equations
    this->store_solvers();                               // Store solvers
    this->compute_cdr_gradients();                       // Compute gradients, then do velo and source
    if(m_compute_v) this->compute_cdr_velo(a_time + a_dt);
    if(m_compute_S) this->compute_cdr_sources(a_time + a_dt);
    this->advance_rk(a_time, 0.5*a_dt);                  // Solvers now contain e^(0.5*A*dt)*e^(0.5*A*dt)*u

    // Last embedded diffusive step tep
    if(m_do_diffusion){
      if(m_consistent_E) this->update_poisson();
      TimeStepper::compute_cdr_diffusion(m_fieldSolver_scratch->get_E_cell(), m_fieldSolver_scratch->get_E_eb());
      this->advance_diffusion(a_time, a_dt);
    }

    // Solvers now contain e^(D*dt)*e^(0.5*A*dt)*e^(0.5*A*dt)*u, but this is the error so we need
    // to copy that to the correct place, and revert the solvers to their original states
    this->copy_solvers_to_error();   // error -> e^(D*dt)*e^(0.5*A*dt)*e^(0.5*A*dt)*u
    this->copy_cache_to_solvers();   // solvers -> e^(0.5*A*dt)*u, this also reverts the potential and RTE solutions
    this->compute_E_into_scratch();  // Field were not reverted (for storage reasons), recompute them here. 
  }

  // Diffusion advance
  if(m_do_diffusion){
    if(m_compute_D){
      TimeStepper::compute_cdr_diffusion(m_fieldSolver_scratch->get_E_cell(), m_fieldSolver_scratch->get_E_eb());
    }
    this->advance_diffusion(a_time, a_dt);
  }

  // Update the electric field and RTE equations again
  if(m_consistent_E) this->update_poisson();
  if(m_consistent_rte) this->update_rte(a_time + a_dt);

  // Solvers need to be refilled with stuff
  this->store_solvers();
  this->compute_cdr_gradients();
  if(m_compute_v) this->compute_cdr_velo(a_time + a_dt);
  if(m_compute_S) this->compute_cdr_sources(a_time + a_dt);
  this->advance_rk(a_time, 0.5*a_dt);
}

Real strang2::advance_fixed(const int a_substeps, const Real a_dt){
  CH_TIME("strang2::advance_fixed");
  if(m_verbosity > 2){
    pout() << "strang2::advance_fixed" << endl;
  }

  this->compute_E_into_scratch(); // Compute the electric field
  this->compute_cdr_gradients();  // Precompute gradients

  Real sum_dt = 0.0;
  for (int step = 0; step < a_substeps; step++){
    const Real time      = m_time + step*a_dt;
    const bool last_step = (step == a_substeps - 1);

    this->advance_one_step(time, a_dt);

    if(m_compute_error){
      this->compute_errors();
    }
    if(m_use_embedded){
      this->copy_error_to_solvers();
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

    sum_dt += a_dt;
  }

  return sum_dt;
}

Real strang2::advance_adaptive(int& a_substeps, Real& a_dt, const Real a_time, const Real a_dtc){
  CH_TIME("strang2::advance_adaptive");
  if(m_verbosity > 2){
    pout() << "strang2::advance_adaptive" << endl;
  }

  this->compute_E_into_scratch();
  this->compute_cdr_gradients();  // Precompute gradients

  Real sum_dt = 0.0;
  const Real fudge = 1E-6;
  while (sum_dt - a_dtc < -fudge*a_dt){
    const bool last_step = (sum_dt + a_dt >= a_dtc - a_dt*1.E-3) ? true : false;
    const Real cur_time  = a_time + sum_dt;

    // Split step advance and compute errors
    this->advance_one_step(cur_time, a_dt);
    if(m_compute_error) {
      this->compute_errors();
    }
    else{
      MayDay::Abort("strang2::advance_adaptive - m_compute_errors = false but this shouldn't happen!");
    }
    if(m_use_embedded){
      this->copy_error_to_solvers();
    }

    const Real old_dt  = a_dt;
    const Real rel_err = m_err_thresh/m_max_error;
#if 1 // Simple way
    const Real new_dt  = (m_max_error > 0.0) ? a_dt*sqrt(rel_err) : m_max_dt;
#else // Better way?
    const Real new_dt  = (m_max_error > 0.0) ? a_dt*Min(m_max_alpha, Max(m_min_alpha, sqrt(m_alpha*rel_err))) : m_max_dt;
#endif
    
    const bool accept_err   = m_max_error <= m_err_thresh;
    const bool decr_order   = accept_err  && (a_dt >= m_maxCFL*m_dt_cfl);
    const bool incr_order   = !accept_err && (a_dt <= m_minCFL*m_dt_cfl);
    const bool accept_order = !decr_order && !accept_order;
    const bool accept_step  = accept_err;

    if(accept_step){
      sum_dt     += a_dt;
      a_substeps += 1;

      // Set the new time step, but only if we are sufficiently far from the error threshold.
      const Real err_factor = 1./rel_err;
      if(err_factor < m_safety){
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

#if 1 // Debug
    if(procID() == 0){
      std::cout << "accept = " << accept_step
		<< "\t accept_order = " << accept_order
		<< "\t RK order = " << m_rk_order
		<< "\t last = " << last_step
		<< "\t rel error = " << m_err_thresh/m_max_error
		<< "\tmax error = " << m_max_error
		<< "\t a_dt = " << old_dt
		<< "\t new_dt = " << a_dt
		<< "\t dt_cfl = " << m_dt_cfl
		<< "\t cfl = " << old_dt/m_dt_cfl 
		<< std::endl;
    }
#endif

    // New time step should never exceed CFL or hardcap constraints
    a_dt = Min(a_dt, m_maxCFL*m_dt_cfl);
    a_dt = Max(a_dt, m_minCFL*m_dt_cfl);
    a_dt = Min(a_dt, m_max_dt);
    a_dt = Max(a_dt, m_min_dt);

    // If we're doing one more substep, we need to make sure that the solvers are up to date. If not, let advance()
    // take care of the rest of the synchronization
    if(!last_step){
      // Update the electric field and RTE equations
      if(m_consistent_E)   this->update_poisson();
      if(m_consistent_rte) this->update_rte(cur_time);

      // Compute new velocities and source terms for the next advective step. No need to compute the
      // diffusion coefficients because they are automateically filled before advance_diffusion() call
      if(m_compute_S) this->compute_cdr_gradients();
      if(m_compute_v) this->compute_cdr_velo(cur_time);
      if(m_compute_S) this->compute_cdr_sources(cur_time);
    }

  }

  return sum_dt;
}

void strang2::advance_rk(const Real a_time, const Real a_dt){
  CH_TIME("strang2::advance_rk");
  if(m_verbosity > 2){
    pout() << "strang2::advance_rk" << endl;
  }

  if(m_rk_order == 2){
    if(m_rk_stages >= 2){
      this->advance_rkN2(a_time, a_dt, m_rk_stages);
    }
    else{
      MayDay::Abort("strang2::advance_rk - unknown second RK method requested");
    }
  }
  else if(m_rk_order == 3){
    if(m_rk_stages == 3){
      this->advance_rk33(a_time, a_dt);
    }
    else if(m_rk_stages == 4){
      this->advance_rk43(a_time, a_dt);
    }
    else if(m_rk_stages == 5){
      this->advance_rk53(a_time, a_dt);
    }
    else{
      MayDay::Abort("strang2::advance_rk - unknown third order RK method requested");
    }
  }
  else if(m_rk_order == 4){
    if(m_rk_stages == 5){
      this->advance_rk54(a_time, a_dt);
    }
    else{
      MayDay::Abort("strang2::advance_rk - unknown fourth order RK method requested");
    }
  }
  else{
    MayDay::Abort("strang2::advance_rk - unknown RK order requested");
  }
}

void strang2::advance_rkN2(const Real a_time, const Real a_dt, const int a_stages){
  CH_TIME("strang2::advance_rkN2");
  if(m_verbosity > 5){
    pout() << "strang2::advance_rkN2" << endl;
  }

  // The tables for this scheme are
  // alpha            
  // -----------------
  //  1               
  //  0    1
  //  0    0   1
  //  :    :   :   1
  //  :    :   :   :   1
  // 1/s   0   0   0   0  (s-1)/s
  //
  // beta 
  // -----------------
  //  b
  //  0   b 
  //  :   0   b
  //  :   :   0   b
  //  0   0   0   0  1/s
  //
  // with b = 1/(s-1)
  // I.e. the first three stages are equal (we put the result into the solver). Put them in a loop. 

  const Real s     = a_stages;
  const Real alpha = (s-1)/s;
  const Real beta  = 1/(s-1);
  const Real sinv  = 1./s;

  // u^(i) = u^(i-1) + third*dt*L(u^(i-1))
  for (int i = 0; i < a_stages-1; i++){
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();

    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& rhs       = storage->get_scratch();
      EBAMRCellData& phi       = solver->getPhi();  // u^i
      const EBAMRCellData& src = solver->getSource(); // S^i

      solver->computeDivF(rhs, phi, 0.0, true); // RHS =  Div(u^i*v^i)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^i*v^i)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^i - Div(u^i*v^i)
      data_ops::incr(phi, rhs, beta*a_dt);       // u^(i+1) = u^i + [1/(s-1)]*dt*L(u^i)

      m_amr->averageDown(phi, m_cdr->get_phase());
      m_amr->interpGhost(phi, m_cdr->get_phase());

      data_ops::floor(phi, 0.0);
    }

    EBAMRIVData& sigma = m_sigma->getPhi(); // sigma^i
    EBAMRIVData& rhs   = m_sigma_scratch->get_scratch();
    m_sigma->computeRHS(rhs);
    data_ops::incr(sigma, rhs, beta*a_dt); // sigma^(i+1) = sigma^i + [1/(s-1)]*dt*F^i

    if(m_consistent_E)   this->update_poisson();
    if(m_consistent_rte) this->update_rte(a_time + a_dt);
    this->compute_cdr_gradients();
    if(m_compute_v) this->compute_cdr_velo(a_time + a_dt);
    if(m_compute_S) this->compute_cdr_sources(a_time + a_dt);
  }


  { // u^(n+1) = (1/s)*u^n + {(s-1)/s]*u^(s-1) + (1/s)*dt*L[u^(s-1)]
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();

    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& rhs       = storage->get_scratch();
      EBAMRCellData& phi       = solver->getPhi();     // u^i
      const EBAMRCellData& src = solver->getSource();
      const EBAMRCellData& pre = storage->get_previous(); // u^n

      solver->computeDivF(rhs, phi, 0.0, true); // RHS =  Div(u^i*v^i)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^i*v^i)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^3 - Div(u^i*v^i)
      data_ops::scale(phi, alpha);               // u^(n+1) = (s-1/s)*u^i
      data_ops::incr(phi, pre, sinv);            // u^(n+1) = (1/s)*u^n + (s-1/s)*u^i
      data_ops::incr(phi, rhs, sinv*a_dt);       // u^(n+1) = (1/s)*u^n + (1/s)*dt*L(u^i)

      m_amr->averageDown(phi, m_cdr->get_phase());
      m_amr->interpGhost(phi, m_cdr->get_phase());

      data_ops::floor(phi, 0.0);
    }

    EBAMRIVData& sigma     = m_sigma->getPhi();            // sigma^3
    EBAMRIVData& rhs       = m_sigma_scratch->get_scratch();
    const EBAMRIVData& pre = m_sigma_scratch->get_previous(); // sigma^n
    m_sigma->computeRHS(rhs);
    data_ops::scale(sigma, alpha);           // sigma^(n+1) = (s-1/s)*sigma^i
    data_ops::incr(sigma, pre, sinv);        // sigma^(n+1) = (1/s)*sigma^n + (s-1/s)*sigma^i
    data_ops::incr(sigma, rhs, sinv*a_dt);   // sigma^(n+1) = (1/s)*sigma^n + (s-1/s)*sigma^i + (1/s)*dt*L(sigma^i)
  }
}

void strang2::advance_rk33(const Real a_time, const Real a_dt){
  CH_TIME("strang2::advance_rk33");
  if(m_verbosity > 2){
    pout() << "strang2::advance_rk33" << endl;
  }

  // The tables for this scheme are
  // alpha            
  // -----------------
  //  1               
  // 3/4  1/4
  // 1/3   0   2/3
  //
  // beta            
  // -----------------
  //  1               
  //  0   1/4
  //  0    0    2/3

  // u^1 = u^n + dt*L(u^n)
  { 
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();
  
    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    
      EBAMRCellData& phi       = solver->getPhi();    // u^n
      EBAMRCellData& rhs       = storage->get_scratch();
      const EBAMRCellData& src = solver->getSource();
    
      solver->computeDivF(rhs, phi, 0.0, true); // RHS =  Div(u^n*v^n)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^n*v^n)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^n - Div(u^n*v^n)
      data_ops::incr(phi, rhs, a_dt);            // u^1 = u^n + dt*L(u^n)
    
      m_amr->averageDown(phi, m_cdr->get_phase());
      m_amr->interpGhost(phi, m_cdr->get_phase());
    
      data_ops::floor(phi, 0.0);

    }

    EBAMRIVData& sigma = m_sigma->getPhi();
    EBAMRIVData& rhs   = m_sigma_scratch->get_scratch();
    m_sigma->computeRHS(rhs);
    data_ops::incr(sigma, rhs, a_dt); // sigma^1 = sigma^n + dt*F^n
  }

  if(m_consistent_E)   this->update_poisson();
  if(m_consistent_rte) this->update_rte(a_time + a_dt);

  // u^2 = 0.25*(3*u^n + u^1 + dt*L(u^1)). 
  { 
    this->compute_cdr_gradients();
    if(m_compute_v) this->compute_cdr_velo(a_time + a_dt);
    if(m_compute_S) this->compute_cdr_sources(a_time + a_dt);
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();

    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    
      EBAMRCellData& rhs       = storage->get_scratch();
      EBAMRCellData& phi       = solver->getPhi();     // phi = u^1 = u^n + dt*L(u^n)
      const EBAMRCellData& src = solver->getSource();
      const EBAMRCellData& pre = storage->get_previous(); // pre = u^n
    
      solver->computeDivF(rhs, phi, 0.0, true); // RHS =  Div(u^1*v^1)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^1*v^1)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^1 - Div(u^1*v^1)
      data_ops::incr(phi,  rhs, a_dt);           // u^2 = u^1 + dt*L(u^1)
      data_ops::incr(phi,  pre, 3.0);            // u^2 = 3*u^n + u^1 + dt*L(u^1)
      data_ops::scale(phi, 0.25);                // u^2 = 0.25*[3*u^n + u^1 + dt*L(u^1)]
    
      m_amr->averageDown(phi, m_cdr->get_phase());
      m_amr->interpGhost(phi, m_cdr->get_phase());
    
      data_ops::floor(phi, 0.0);
    }

    EBAMRIVData& sigma = m_sigma->getPhi();           // u^1
    EBAMRIVData& rhs = m_sigma_scratch->get_scratch();   // Storage for right hand side
    EBAMRIVData& pre = m_sigma_scratch->get_previous();  // u^n
    m_sigma->computeRHS(rhs);
    data_ops::incr(sigma, rhs, a_dt);  // sigma = u^1 + dt*L(u^1)
    data_ops::incr(sigma, pre, 3.0);   // sigma = 3*u^n + u^1 + dt*L(u^1)
    data_ops::scale(sigma, 0.25);      // sigma = 0.25*[3*u^n + u^1 + dt*L(u^1)]
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

    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    
      EBAMRCellData& rhs       = storage->get_scratch();
      EBAMRCellData& phi       = solver->getPhi();     // phi = u^2
      const EBAMRCellData& src = solver->getSource();
      const EBAMRCellData& pre = storage->get_previous(); // pre = u^n
    
      solver->computeDivF(rhs, phi, 0.0, true); // RHS =  Div(u^2*v^2)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^2*v^2)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^2 - Div(u^2*v^2)
      data_ops::incr(phi, rhs, a_dt);            // u^(n+1) = u^2 + dt*L(u^2)
      data_ops::scale(phi, 2.0);                 // u^(n+1) = 2*u^2 + 2*dt*L(u^2)
      data_ops::incr(phi, pre, 1.0);             // u^(n+1) = u^n + 2*u^2 + dt*L(u^2)
      data_ops::scale(phi, 1./3.);               // u^(n+1) = (1/3)*[u^n + 2*u^2 + 2*dt*L(u^2)]
    
      m_amr->averageDown(phi, m_cdr->get_phase());
      m_amr->interpGhost(phi, m_cdr->get_phase());
    
      data_ops::floor(phi, 0.0);

    }

    EBAMRIVData& sigma = m_sigma->getPhi();           // u^2
    EBAMRIVData& rhs = m_sigma_scratch->get_scratch();   // Storage for right hand side
    EBAMRIVData& pre = m_sigma_scratch->get_previous();  // u^n
    m_sigma->computeRHS(rhs);
    data_ops::incr(sigma, rhs, a_dt);  // sigma = u^2 + dt*L(u^2)
    data_ops::scale(sigma, 2.0);       // sigma = 2*u^2 + 2*dt*L(u^2)
    data_ops::incr(sigma, pre, 1.0);   // sigma = u^n + 2*u^2 + 2*dt*L(u^2)
    data_ops::scale(sigma, 1./3.);     // sigma = (1/3)*[u^n + 2*u^2 + 2*dt*L(u^2)]
  }
}

void strang2::advance_rk43(const Real a_time, const Real a_dt){
  CH_TIME("strang2::advance_rk43");
  if(m_verbosity > 5){
    pout() << "strang2::advance_rk43" << endl;
  }

  // The tables for this scheme are
  // alpha            
  // -----------------
  //  1               
  //  0    1
  // 2/3   0   1/3
  //  0    0    0    1
  //
  // beta            
  // -----------------
  // 1/2               
  //  0   1/2
  //  0    0   1/6
  //  0    0    0   1/2
  //
  // The first two stages are equal, so they can be loop'ed.

  const Real third = 1./3.;
  const Real sixth = 1./6;

  // u^(i) = u^(i-1) + third*dt*L(u^(i-1))
  for (int i = 0; i <= 1; i++){
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();

    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& rhs       = storage->get_scratch();
      EBAMRCellData& phi       = solver->getPhi();  // u^i
      const EBAMRCellData& src = solver->getSource(); // S^i

      solver->computeDivF(rhs, phi, 0.0, true); // RHS =  Div(u^i*v^i)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^i*v^i)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^i - Div(u^i*v^i)
      data_ops::incr(phi, rhs, 0.5*a_dt);      // u^(i+1) = u^i + 0.5*dt*L(u^i)

      m_amr->averageDown(phi, m_cdr->get_phase());
      m_amr->interpGhost(phi, m_cdr->get_phase());

      data_ops::floor(phi, 0.0);
    }

    EBAMRIVData& sigma = m_sigma->getPhi(); // sigma^i
    EBAMRIVData& rhs   = m_sigma_scratch->get_scratch();
    m_sigma->computeRHS(rhs);
    data_ops::incr(sigma, rhs, 0.5*a_dt); // sigma^(i+1) = sigma^i + 0.5*dt*F^i

    if(m_consistent_E)   this->update_poisson();
    if(m_consistent_rte) this->update_rte(a_time + a_dt);
    this->compute_cdr_gradients();
    if(m_compute_v) this->compute_cdr_velo(a_time + a_dt);
    if(m_compute_S) this->compute_cdr_sources(a_time + a_dt);
  }


  // u^3 = (2/3)*u^n + (1/3)*u^2 + (1/6)*dt*L(u^2)]
  { 
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();

    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& rhs       = storage->get_scratch();
      EBAMRCellData& phi       = solver->getPhi();     // u^2
      const EBAMRCellData& src = solver->getSource();
      const EBAMRCellData& pre = storage->get_previous(); // u^n

      solver->computeDivF(rhs, phi, 0.0, true); // RHS =  Div(u^2*v^2)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^2*v^2)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^2 - Div(u^2*v^2)
      data_ops::scale(phi, 2.0);                 // u^(n+1) = 2*u^2
      data_ops::incr(phi, pre, 4.0);             // u^(n+1) = 4*u^n + 2*u^2
      data_ops::incr(phi, rhs, a_dt);            // u^(n+1) = 4*u^n + 2*u^2 + dt*L(u^2)
      data_ops::scale(phi, sixth);               // u^(n+1) = (2/3)u^n + (1/3)*u^3 + (1/6)*dt*L(u^3)]

      m_amr->averageDown(phi, m_cdr->get_phase());
      m_amr->interpGhost(phi, m_cdr->get_phase());

      data_ops::floor(phi, 0.0);
    }

    EBAMRIVData& sigma     = m_sigma->getPhi();            // sigma^2
    EBAMRIVData& rhs       = m_sigma_scratch->get_scratch();
    const EBAMRIVData& pre = m_sigma_scratch->get_previous(); // sigma^n
    m_sigma->computeRHS(rhs);
    data_ops::scale(sigma, 2.0);            // sigma^(n+1) = 2*sigma^2
    data_ops::incr(sigma, pre, 4.0);        // sigma^(n+1) = 4*sigma^n + 2*sigma^2
    data_ops::incr(sigma, rhs, a_dt);       // sigma^(n+1) = 4*sigma^n + 2*sigma^2 + dt*F^2
    data_ops::scale(sigma, sixth);          // sigma^(n+1) = (2/3)/*sigma^n + (1/3)*sigma^2 + (1/6)*dt*F^3

    if(m_consistent_E)   this->update_poisson();
    if(m_consistent_rte) this->update_rte(a_time + a_dt);
    this->compute_cdr_gradients();
    if(m_compute_v) this->compute_cdr_velo(a_time + a_dt);
    if(m_compute_S) this->compute_cdr_sources(a_time + a_dt);
  }

  // u^(n+1) = u^3 + 0.5*dt*L(u^3)
  {
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();

    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& rhs       = storage->get_scratch();
      EBAMRCellData& phi       = solver->getPhi();  // u^i
      const EBAMRCellData& src = solver->getSource(); // S^i

      solver->computeDivF(rhs, phi, 0.0, true); // RHS =  Div(u^i*v^i)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^i*v^i)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^i - Div(u^i*v^i)
      data_ops::incr(phi, rhs, 0.5*a_dt);      // u^(i+1) = u^i + 0.5*dt*L(u^i)

      m_amr->averageDown(phi, m_cdr->get_phase());
      m_amr->interpGhost(phi, m_cdr->get_phase());

      data_ops::floor(phi, 0.0);
    }

    EBAMRIVData& sigma = m_sigma->getPhi(); // sigma^i
    EBAMRIVData& rhs   = m_sigma_scratch->get_scratch();
    m_sigma->computeRHS(rhs);
    data_ops::incr(sigma, rhs, 0.5*a_dt); // sigma^(i+1) = sigma^i + 0.5*dt*F^i
  }
}

void strang2::advance_rk53(const Real a_time, const Real a_dt){
  CH_TIME("strang2::advance_rk53");
  if(m_verbosity > 5){
    pout() << "strang2::advance_rk53" << endl;
  }

  // The tables for this scheme are
  // alpha            
  // -----------------
  //  1
  //  0    1
  // a20   0   a22
  // a30  a31   0    a33
  // a40  a41  a42    0   a44
  //
  // beta            
  // -----------------
  //  1
  //  0   b11
  //  0    0   b22
  // b30   0    0    b33
  // b40  b41   0     0   b44
  //
  // This scheme requires extra storage for u^1, u^2, L(u^1), L(u^2) from the final update. We will put 
  //
  // u^1    -> extra[0]
  // u^2    -> extra[1]
  // L(u^1) -> extra[2]
  // L(u^2  -> extra[3]

  const Real a00 = 1.0;
  const Real a20 = 0.56656131914033;
  const Real a30 = 0.09299483444413;
  const Real a40 = 0.00736132260920;
  const Real a11 = 1.0;
  const Real a31 = 0.00002090369620;
  const Real a41 = 0.20127980325145;
  const Real a22 = 0.43343868085967;
  const Real a42 = 0.00182955389682;
  const Real a33 = 0.90698426185967;
  const Real a44 = 0.78952932024253;

  const Real b00 = 0.37726891511710;
  const Real b30 = 0.00071997378654;
  const Real b40 = 0.00277719819460;
  const Real b11 = 0.37726891511710;
  const Real b41 = 0.00001567934613;
  const Real b22 = 0.16352294089771;
  const Real b33 = 0.34217696850008;
  const Real b44 = 0.29786487010104;

  // Allocate extra storage on way in
  {
    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);
      storage->allocate_extra_storage(4);
    }
    m_sigma_scratch->allocate_extra_storage(4);
  }

  // u^1 = u^n + dt*b00*L(u^n)
  // 
  { 
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();

    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& rhs       = *(storage->get_extra_storage()[2]); // Will become L(u^n)
      EBAMRCellData& phi       = solver->getPhi();                // u^n
      EBAMRCellData& u1        = *(storage->get_extra_storage()[0]); // Becomes u^1
      const EBAMRCellData& src = solver->getSource();               // S^n

      solver->computeDivF(rhs, phi, 0.0, true); // RHS =  Div(u^n*v^n)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^n*v^n)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^n - Div(u^n*v^n) = L(u^n)
      data_ops::incr(phi, rhs, a_dt*b00);        // u^1 = u^n + dt*b00*L(u^n)

      m_amr->averageDown(phi, m_cdr->get_phase());
      m_amr->interpGhost(phi, m_cdr->get_phase());

      data_ops::floor(phi, 0.0);
      data_ops::copy(u1, phi); // Backup of u^1
    }

    EBAMRIVData& sigma  = m_sigma->getPhi(); // sigma^n
    EBAMRIVData& rhs    = *(m_sigma_scratch->get_extra_storage()[2]); // L(sigma^n)
    EBAMRIVData& sigma1 = *(m_sigma_scratch->get_extra_storage()[0]); // Will become sigma^1
    m_sigma->computeRHS(rhs);
    data_ops::incr(sigma, rhs, a_dt*b00);   // sigma^1 = sigma^n + dt*b00*L(sigma^n)
    data_ops::set_value(sigma1, 0.0);
    data_ops::incr(sigma1, sigma, 1.0);     // Backup of sigma^1

    if(m_consistent_E)   this->update_poisson();
    if(m_consistent_rte) this->update_rte(a_time + a_dt);
    this->compute_cdr_gradients();
    if(m_compute_v) this->compute_cdr_velo(a_time + a_dt);
    if(m_compute_S) this->compute_cdr_sources(a_time + a_dt);
  }

  // u^2 = u^1 + b11*L(u^1)
  { 
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();

    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& rhs       = *(storage->get_extra_storage()[3]); // Will become L(u^1)
      EBAMRCellData& phi       = solver->getPhi();                // u^1
      EBAMRCellData& u2        = *(storage->get_extra_storage()[1]); // Becomes u^2
      const EBAMRCellData& src = solver->getSource();               // S^1

      solver->computeDivF(rhs, phi, 0.0, true); // RHS =  Div(u^1*v^1)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^1*v^1)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^1 - Div(u^1*v^1) = L(u^1)
      data_ops::incr(phi, rhs, a_dt*b11);        // u^2 = u^1 + dt*b11*L(u^1)

      m_amr->averageDown(phi, m_cdr->get_phase());
      m_amr->interpGhost(phi, m_cdr->get_phase());

      data_ops::floor(phi, 0.0);
      data_ops::copy(u2, phi);   // Backup of u^2
    }

    EBAMRIVData& sigma  = m_sigma->getPhi(); // sigma^1
    EBAMRIVData& rhs    = *(m_sigma_scratch->get_extra_storage()[3]); // Will become L(sigma^1)
    EBAMRIVData& sigma2 = *(m_sigma_scratch->get_extra_storage()[1]); // Will become sigma^1
    m_sigma->computeRHS(rhs);
    data_ops::incr(sigma, rhs, a_dt*b11);   // sigma^2 = sigma^1 + dt*b11*L(sigma^n)
    data_ops::set_value(sigma2, 0.0);
    data_ops::incr(sigma2, sigma, 1.0);     // Backup of sigma^1

    if(m_consistent_E)   this->update_poisson();
    if(m_consistent_rte) this->update_rte(a_time + a_dt);
    this->compute_cdr_gradients();
    if(m_compute_v) this->compute_cdr_velo(a_time + a_dt);
    if(m_compute_S) this->compute_cdr_sources(a_time + a_dt);
  }

  // u^3 = a20*u^n + a22*u^2 + b22*L(u^2)
  { 
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();

    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& phi       = solver->getPhi();     // u^2
      EBAMRCellData& rhs       = storage->get_scratch();  // RHS     
      EBAMRCellData& pre       = storage->get_previous(); // u^n
      const EBAMRCellData& src = solver->getSource();    // S^3

      solver->computeDivF(rhs, phi, 0.0, true); // RHS =  Div(u^2*v^2)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^2*v^2)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^2 - Div(u^2*v^2) = L(u^2)
      data_ops::scale(phi, a22);                 // u^3 = a22*u^2
      data_ops::incr(phi, pre, a20);             // u^3 = a20*u^n + a22*u^2
      data_ops::incr(phi, rhs, a_dt*b22);        // u^3 = a20*u^n + a22*u^2 + dt*b22*L(u^2)

      m_amr->averageDown(phi, m_cdr->get_phase());
      m_amr->interpGhost(phi, m_cdr->get_phase());

      data_ops::floor(phi, 0.0);
    }

    EBAMRIVData& sigma  = m_sigma->getPhi();            // sigma^2
    EBAMRIVData& rhs    = m_sigma_scratch->get_scratch();  // RHS
    EBAMRIVData& pre    = m_sigma_scratch->get_previous(); // sigma^n
    m_sigma->computeRHS(rhs);
    data_ops::scale(sigma, a22);                           // sigma^3 = a2*sigma^2
    data_ops::incr(sigma, pre, a20);                       // sigma^3 = a20*sigma^n + a22*sigma^2
    data_ops::incr(sigma, rhs, a_dt*b22);                  // sigma^3 = a20*sigma^n + a22*sigma^2 + b22*dt*L(sigma^2)

    if(m_consistent_E)   this->update_poisson();
    if(m_consistent_rte) this->update_rte(a_time + a_dt);
    this->compute_cdr_gradients();
    if(m_compute_v) this->compute_cdr_velo(a_time + a_dt);
    if(m_compute_S) this->compute_cdr_sources(a_time + a_dt);
  }

  // u^4 = a30*u^n + a31*u^1 + a33*u^3 + b30*dt*L(u^n) + b33*dt*L(u^3)
  { 
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();

    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& phi       = solver->getPhi();                // u^3
      EBAMRCellData& rhs       = storage->get_scratch();             // RHS     
      const EBAMRCellData& pre = storage->get_previous();            // u^n
      const EBAMRCellData& u1  = *(storage->get_extra_storage()[0]); // u^1
      const EBAMRCellData& Lun = *(storage->get_extra_storage()[2]); // L(u^n)
      const EBAMRCellData& src = solver->getSource();               // S^3

      solver->computeDivF(rhs, phi, 0.0, true); // RHS =  Div(u^3*v^3)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^3*v^3)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^3 - Div(u^3*v^3)
      data_ops::scale(phi, a33);                 // u^4 = a33*u^3
      data_ops::incr(phi, pre, a30);             // u^4 = a30*u^n + a33*u^3
      data_ops::incr(phi, u1,  a31);             // u^4 = a30*u^n + a31*u^1 + a33*u^3
      data_ops::incr(phi, Lun, a_dt*b30);        // u^4 = a30*u^n + a31*u^1 + a33*u^3 + dt*b30*L(u^n)
      data_ops::incr(phi, rhs, a_dt*b33);        // u^4 = a30*u^n + a31*u^1 + a33*u^3 + dt*b30*L(u^n) + dt*b33*L(u^3)

      m_amr->averageDown(phi, m_cdr->get_phase());
      m_amr->interpGhost(phi, m_cdr->get_phase());

      data_ops::floor(phi, 0.0);
    }

    EBAMRIVData& sigma         = m_sigma->getPhi();               // sigma^3
    EBAMRIVData& rhs           = m_sigma_scratch->get_scratch();     // RHS
    const EBAMRIVData& pre     = m_sigma_scratch->get_previous();    // sigma^n
    const EBAMRIVData& sigma1  = *(m_sigma_scratch->get_extra_storage()[0]); // sigma^1
    const EBAMRIVData& Lsigman = *(m_sigma_scratch->get_extra_storage()[2]); // L(sigma^n)
    m_sigma->computeRHS(rhs);
    data_ops::scale(sigma, a33);                // sigma^4 = a33*sigma^3
    data_ops::incr(sigma,  pre,     a30);       // sigma^4 = a30*sigma^n + a33*sigma^3
    data_ops::incr(sigma,  sigma1,  a31);       // sigma^4 = a30*sigma^n + a31*sigma^1 + a33*sigma^3
    data_ops::incr(sigma,  Lsigman, a_dt*b30);  // sigma^3 = a30*sigma^n + a31*sigma^1 + a33*sigma^3 + b30*dt*L(sigma^n)
    data_ops::incr(sigma,  rhs,     a_dt*b33);  // sigma^3 = a30*sigma^n + a31*sigma^1 + a33*sigma^3 + b30*dt*L(sigma^n) + a33*dt

    if(m_consistent_E)   this->update_poisson();
    if(m_consistent_rte) this->update_rte(a_time + a_dt);
    this->compute_cdr_gradients();
    if(m_compute_v) this->compute_cdr_velo(a_time + a_dt);
    if(m_compute_S) this->compute_cdr_sources(a_time + a_dt);
  }

  // u^(n+1) = a40*u^n + a41*u^1 + a42*u^2 + a44*u^4 + b40*dt*L(u^n) + b41*dt*L(u^1) + b44*dt*L(u^4)
  { 
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();

    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& phi       = solver->getPhi();                // u^4
      EBAMRCellData& rhs       = storage->get_scratch();             // RHS     
      const EBAMRCellData& pre = storage->get_previous();            // u^n
      const EBAMRCellData& u1  = *(storage->get_extra_storage()[0]); // u^1
      const EBAMRCellData& u2  = *(storage->get_extra_storage()[1]); // u^2
      const EBAMRCellData& Lun = *(storage->get_extra_storage()[2]); // L(u^n)
      const EBAMRCellData& Lu1 = *(storage->get_extra_storage()[3]); // L(u^1)
      const EBAMRCellData& src = solver->getSource();               // S^4

      solver->computeDivF(rhs, phi, 0.0, true); // RHS =  Div(u^4*v^4)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^4*v^4)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^4 - Div(u^4*v^4)
      data_ops::scale(phi, a44);                 // u^(n+1)  = a44*u^4
      data_ops::incr(phi, pre, a40);             // u^(n+1) += a40*u^n
      data_ops::incr(phi, u1,  a41);             // u^(n+1) += a41*u^1 
      data_ops::incr(phi, u2,  a42);             // u^(n+1) += a42*u^2
      data_ops::incr(phi, Lun, a_dt*b40);        // u^(n+1) += b40*dt*L(u^n)
      data_ops::incr(phi, Lu1, a_dt*b41);        // u^(n+1) += b41*dt*L(u^1)
      data_ops::incr(phi, rhs, a_dt*b44);        // u^(n+1) += b44*dt*L(u^4)

      m_amr->averageDown(phi, m_cdr->get_phase());
      m_amr->interpGhost(phi, m_cdr->get_phase());

      data_ops::floor(phi, 0.0);
    }

    EBAMRIVData& sigma         = m_sigma->getPhi();               // sigma^4
    EBAMRIVData& rhs           = m_sigma_scratch->get_scratch();     // RHS
    const EBAMRIVData& pre     = m_sigma_scratch->get_previous();    // sigma^n
    const EBAMRIVData& sigma1  = *(m_sigma_scratch->get_extra_storage()[0]); // sigma^1
    const EBAMRIVData& sigma2  = *(m_sigma_scratch->get_extra_storage()[1]); // sigma^2
    const EBAMRIVData& Lsigman = *(m_sigma_scratch->get_extra_storage()[2]); // L(sigma^n)
    const EBAMRIVData& Lsigma1 = *(m_sigma_scratch->get_extra_storage()[3]); // L(sigma^n)
    m_sigma->computeRHS(rhs);
    data_ops::scale(sigma, a44);                // sigma^(n+1) = a44*sigma^4
    data_ops::incr(sigma,  pre,     a40);       // sigma^(n+1) += a40*sigma^n
    data_ops::incr(sigma,  sigma1,  a41);       // sigma^(n+1) += a41*sigma^1
    data_ops::incr(sigma,  sigma2,  a42);       // sigma^(n+1) += a42*sigma^2
    data_ops::incr(sigma,  Lsigman, a_dt*b40);  // sigma^(n+1) += b40*dt*L(sigma^n)
    data_ops::incr(sigma,  Lsigma1, a_dt*b41);  // sigma^(n+1) += b41*dt*L(sigma^1)
    data_ops::incr(sigma,  rhs,     a_dt*b44);  // sigma^(n+1) += b44*dt*L(sigma^4)

    if(m_consistent_E)   this->update_poisson();
    if(m_consistent_rte) this->update_rte(a_time + a_dt);
    this->compute_cdr_gradients();
    if(m_compute_v) this->compute_cdr_velo(a_time + a_dt);
    if(m_compute_S) this->compute_cdr_sources(a_time + a_dt);
  }

  // Deallocate extra storage on way out
  {
    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);
      storage->deallocate_extra_storage();
    }
    m_sigma_scratch->deallocate_extra_storage();
  }
}

void strang2::advance_rk54(const Real a_time, const Real a_dt){
  CH_TIME("strang2::advance_rk54");
  if(m_verbosity > 5){
    pout() << "strang2::advance_rk54" << endl;
  }

  // The tables for this scheme are
  // alpha            
  // -----------------
  // a00
  // a10  a11
  // a20   0   a22
  // a30   0   0    a33
  // a40   0   a42  a43  a44
  //
  // beta            
  // -----------------
  //  1
  //  0   b11
  //  0    0   b22
  //  0    0    0    b33
  //  0    0    0    b43   b44
  //
  // This scheme requires extra storage for u^2, u^3, L(u^3)
  //
  // u^2    -> extra[0]
  // u^3    -> extra[1]
  // L(u^3) -> extra[2]

  const Real a00 = 1.0;
  const Real a10 = 0.44437049406734;
  const Real a20 = 0.62010185138540;
  const Real a30 = 0.17807995410773;
  const Real a40 = 0.00683325884039;
  const Real a11 = 0.55562950593266;
  const Real a22 = 0.37989814861469;
  const Real a42 = 0.51723167208978;
  const Real a33 = 0.82192004589227;
  const Real a43 = 0.12759831133288;
  const Real a44 = 0.34833675773694;

  const Real b00 = 0.39175222700392;
  const Real b11 = 0.36841029262959;
  const Real b22 = 0.25189177424738;
  const Real b33 = 0.54497475021237;
  const Real b43 = 0.08460416338212;
  const Real b44 = 0.22600748319395;

  // Allocate extra storage on way in
  {
    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);
      storage->allocate_extra_storage(3);
    }
    m_sigma_scratch->allocate_extra_storage(3);
  }

  // u^1 = u^n + dt*b00*L(u^n)
  { 
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();

    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& rhs       = storage->get_scratch(); // rhs
      EBAMRCellData& phi       = solver->getPhi();    // u^n
      const EBAMRCellData& src = solver->getSource();   // S^n

      solver->computeDivF(rhs, phi, 0.0, true); // RHS =  Div(u^n*v^n)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^n*v^n)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^n - Div(u^n*v^n) = L(u^n)
      data_ops::incr(phi, rhs, a_dt*b00);        // u^1 = u^n + dt*b00*L(u^n)

      m_amr->averageDown(phi, m_cdr->get_phase());
      m_amr->interpGhost(phi, m_cdr->get_phase());

      data_ops::floor(phi, 0.0);
    }

    EBAMRIVData& sigma  = m_sigma->getPhi();           // sigma^n
    EBAMRIVData& rhs    = m_sigma_scratch->get_scratch(); // rhs
    m_sigma->computeRHS(rhs);
    data_ops::incr(sigma, rhs, a_dt*b00);   // sigma^1 = sigma^n + dt*b00*L(sigma^n)

    if(m_consistent_E)   this->update_poisson();
    if(m_consistent_rte) this->update_rte(a_time + a_dt);
    this->compute_cdr_gradients();
    if(m_compute_v) this->compute_cdr_velo(a_time + a_dt);
    if(m_compute_S) this->compute_cdr_sources(a_time + a_dt);
  }

  // u^2 = a10*u^n + a11*u^1 + b11*L(u^1)
  { 
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();

    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& rhs       = storage->get_scratch();  // rhs
      EBAMRCellData& phi       = solver->getPhi();     // u^1
      EBAMRCellData& u2        = *(storage->get_extra_storage()[0]); // Becomes u^2
      const EBAMRCellData& pre = storage->get_previous(); // u^n
      const EBAMRCellData& src = solver->getSource();    // S^1

      solver->computeDivF(rhs, phi, 0.0, true); // RHS =  Div(u^1*v^1)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^1*v^1)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^1 - Div(u^1*v^1) = L(u^1)
      data_ops::scale(phi, a11);                 // u^2 = a11*u^1
      data_ops::incr(phi, pre, a10);             // u^2 = a10*u^n + a11*u^1
      data_ops::incr(phi, rhs, a_dt*b11);        // u^2 = a10*u^n + a11*u^1 + dt*b11*L(u^1)

      m_amr->averageDown(phi, m_cdr->get_phase());
      m_amr->interpGhost(phi, m_cdr->get_phase());

      data_ops::floor(phi, 0.0);
      data_ops::copy(u2, phi);   // Backup of u^2
    }

    EBAMRIVData& sigma     = m_sigma->getPhi();           // sigma^1
    EBAMRIVData& rhs       = m_sigma_scratch->get_scratch(); // RHS
    EBAMRIVData& sigma2    = *(m_sigma_scratch->get_extra_storage()[0]); // Will become sigma^2
    const EBAMRIVData& pre = m_sigma_scratch->get_previous();
    m_sigma->computeRHS(rhs);
    data_ops::scale(sigma, a11);            // sigma^2 = a11*sigma^1
    data_ops::incr(sigma, pre, a10);        // sigma^2 = a10*sigma^n + a11*sigma^1
    data_ops::incr(sigma, rhs, a_dt*b11);   // sigma^2 = a10*sigma^n + a11*sigma^1 + dt*b11*L(sigma^n)
    data_ops::set_value(sigma2, 0.0);
    data_ops::incr(sigma2, sigma, 1.0);     // Backup of sigma^2

    if(m_consistent_E)   this->update_poisson();
    if(m_consistent_rte) this->update_rte(a_time + a_dt);
    this->compute_cdr_gradients();
    if(m_compute_v) this->compute_cdr_velo(a_time + a_dt);
    if(m_compute_S) this->compute_cdr_sources(a_time + a_dt);
  }

  // u^3 = a20*u^n + a22*u^2 + b22*L(u^2)
  { 
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();

    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& phi       = solver->getPhi();     // u^2
      EBAMRCellData& rhs       = storage->get_scratch();  // RHS
      EBAMRCellData& u3        = *(storage->get_extra_storage()[1]); // Becomes u^3
      const EBAMRCellData& pre = storage->get_previous(); // u^n
      const EBAMRCellData& src = solver->getSource();    // S^3

      solver->computeDivF(rhs, phi, 0.0, true); // RHS =  Div(u^2*v^2)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^2*v^2)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^2 - Div(u^2*v^2) = L(u^2)
      data_ops::scale(phi, a22);                 // u^3 = a22*u^2
      data_ops::incr(phi, pre, a20);             // u^3 = a20*u^n + a22*u^2
      data_ops::incr(phi, rhs, a_dt*b22);        // u^3 = a20*u^n + a22*u^2 + dt*b22*L(u^2)

      m_amr->averageDown(phi, m_cdr->get_phase());
      m_amr->interpGhost(phi, m_cdr->get_phase());

      data_ops::floor(phi, 0.0);
      data_ops::copy(u3, phi);   // Backup of u^3
    }

    EBAMRIVData& sigma     = m_sigma->getPhi();            // sigma^2
    EBAMRIVData& rhs       = m_sigma_scratch->get_scratch();  // RHS
    EBAMRIVData& sigma3    = *(m_sigma_scratch->get_extra_storage()[1]); // Becomes sigma^3
    const EBAMRIVData& pre = m_sigma_scratch->get_previous(); // sigma^n
    m_sigma->computeRHS(rhs);
    data_ops::scale(sigma, a22);                           // sigma^3 = a2*sigma^2
    data_ops::incr(sigma, pre, a20);                       // sigma^3 = a20*sigma^n + a22*sigma^2
    data_ops::incr(sigma, rhs, a_dt*b22);                  // sigma^3 = a20*sigma^n + a22*sigma^2 + b22*dt*L(sigma^2)
    data_ops::set_value(sigma3, 0.0);
    data_ops::incr(sigma3, sigma, 1.0);

    if(m_consistent_E)   this->update_poisson();
    if(m_consistent_rte) this->update_rte(a_time + a_dt);
    this->compute_cdr_gradients();
    if(m_compute_v) this->compute_cdr_velo(a_time + a_dt);
    if(m_compute_S) this->compute_cdr_sources(a_time + a_dt);
  }

  // u^4 = a30*u^n + a33*u^3 + b33*dt*L(u^3)
  { 
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();

    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& phi       = solver->getPhi();                // u^3
      EBAMRCellData& rhs       = *(storage->get_extra_storage()[2]); // RHS, becomes L(u^3) after RHS computation
      const EBAMRCellData& pre = storage->get_previous();            // u^n
      const EBAMRCellData& src = solver->getSource();               // S^3

      solver->computeDivF(rhs, phi, 0.0, true); // RHS =  Div(u^3*v^3)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^3*v^3)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^3 - Div(u^3*v^3)
      data_ops::scale(phi, a33);                 // u^4 = a33*u^3
      data_ops::incr(phi, pre, a30);             // u^4 = a30*u^n + a33*u^3
      data_ops::incr(phi, rhs, a_dt*b33);        // u^4 = a30*u^n + a33*u^3 + dt*b33*L(u^3)

      m_amr->averageDown(phi, m_cdr->get_phase());
      m_amr->interpGhost(phi, m_cdr->get_phase());

      data_ops::floor(phi, 0.0);
    }

    EBAMRIVData& sigma         = m_sigma->getPhi();                       // sigma^3
    EBAMRIVData& rhs           = *(m_sigma_scratch->get_extra_storage()[2]); // Becomes L(sigma^3) after RHS computation
    const EBAMRIVData& pre     = m_sigma_scratch->get_previous();            // sigma^n
    m_sigma->computeRHS(rhs);
    data_ops::scale(sigma, a33);            // sigma^4 = a33*sigma^3
    data_ops::incr(sigma,  pre, a30);       // sigma^4 = a30*sigma^n + a33*sigma^3
    data_ops::incr(sigma,  rhs, a_dt*b33);  // sigma^3 = a30*sigma^n + a33*sigma^3 + b33*dt*L(u^3)

    if(m_consistent_E)   this->update_poisson();
    if(m_consistent_rte) this->update_rte(a_time + a_dt);
    this->compute_cdr_gradients();
    if(m_compute_v) this->compute_cdr_velo(a_time + a_dt);
    if(m_compute_S) this->compute_cdr_sources(a_time + a_dt);
  }

  // u^(n+1) = a40*u^n + a42*u^2 + a43*u^3 + a44*u^4 + b43*L(u^3) + b44*L(u^4)
  { 
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(a_time);
    this->compute_cdr_domain_fluxes(a_time);
    this->compute_sigma_flux();

    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver   = solver_it();
      RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

      EBAMRCellData& phi       = solver->getPhi();                // u^4
      EBAMRCellData& rhs       = storage->get_scratch();             // RHS     
      const EBAMRCellData& pre = storage->get_previous();            // u^n
      const EBAMRCellData& u2  = *(storage->get_extra_storage()[0]); // u^2
      const EBAMRCellData& u3  = *(storage->get_extra_storage()[1]); // u^3
      const EBAMRCellData& Lu3 = *(storage->get_extra_storage()[2]); // L(u^3)
      const EBAMRCellData& src = solver->getSource();               // S^4

      solver->computeDivF(rhs, phi, 0.0, true); // RHS =  Div(u^4*v^4)
      data_ops::scale(rhs, -1.0);                // RHS = -Div(u^4*v^4)
      data_ops::incr(rhs,  src, 1.0);            // RHS = S^4 - Div(u^4*v^4)
      
      data_ops::scale(phi, a44);                 // u^(n+1) = a44*u^4
      data_ops::incr(phi, pre, a40);             // u^(n+1) = a40*u^n + a44*u^4
      data_ops::incr(phi, u2,  a42);             // u^(n+1) = a40*u^n + a42*u^2 + a44*u^4
      data_ops::incr(phi, u3,  a43);             // u^(n+1) = a40*u^n + a42*u^2 + a43*u^3 a44*u^4
      data_ops::incr(phi, Lu3, a_dt*b43);        // u^(n+1) = a40*u^n + a42*u^2 + a43*u^3 a44*u^4 + b43*dt*L(u^3)
      data_ops::incr(phi, rhs, a_dt*b44);        // u^(n+1) = a40*u^n + a42*u^2 + a43*u^3 a44*u^4 + b43*dt*L(u^3) + b44*dt*L(u^4)

      m_amr->averageDown(phi, m_cdr->get_phase());
      m_amr->interpGhost(phi, m_cdr->get_phase());

      data_ops::floor(phi, 0.0);
    }

    EBAMRIVData& sigma         = m_sigma->getPhi();               // sigma^4
    EBAMRIVData& rhs           = m_sigma_scratch->get_scratch();     // RHS
    const EBAMRIVData& pre     = m_sigma_scratch->get_previous();    // sigma^n
    const EBAMRIVData& sigma2  = *(m_sigma_scratch->get_extra_storage()[0]); // sigma^1
    const EBAMRIVData& sigma3  = *(m_sigma_scratch->get_extra_storage()[1]); // sigma^2
    const EBAMRIVData& Lsigma3 = *(m_sigma_scratch->get_extra_storage()[2]); // L(sigma^3)
    m_sigma->computeRHS(rhs);
    data_ops::scale(sigma, a44);                // sigma^(n+1)  = a44*sigma^4
    data_ops::incr(sigma,  pre,     a40);       // sigma^(n+1)  = a40*sigma^n + a44*sigma^4
    data_ops::incr(sigma,  sigma2,  a42);       // sigma^(n+1)  = a40*sigma^n + a42*sigma^2 + a44*sigma^4
    data_ops::incr(sigma,  sigma3,  a43);       // sigma^(n+1)  = a40*sigma^n + a42*sigma^2 + a43*sigma^3 + a44*sigma^4
    data_ops::incr(sigma,  Lsigma3, a_dt*b43);  // sigma^(n+1) += b43*dt*L(sigma^3)
    data_ops::incr(sigma,  rhs,     a_dt*b44);  // sigma^(n+1) += b44*dt*L(sigma^4)

    if(m_consistent_E)   this->update_poisson();
    if(m_consistent_rte) this->update_rte(a_time + a_dt);
    this->compute_cdr_gradients();
    if(m_compute_v) this->compute_cdr_velo(a_time + a_dt);
    if(m_compute_S) this->compute_cdr_sources(a_time + a_dt);
  }

  // Deallocate extra storage on way out
  {
    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);
      storage->deallocate_extra_storage();
    }
    m_sigma_scratch->deallocate_extra_storage();
  }
}

void strang2::compute_cdr_gradients(){
  CH_TIME("strang2::compute_cdr_gradients");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_gradients" << endl;
  }

  strang2::compute_cdr_gradients(m_cdr->getPhis());
}

void strang2::compute_cdr_gradients(const Vector<EBAMRCellData*>& a_phis){
  CH_TIME("strang2::compute_cdr_gradients");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_gradients" << endl;
  }

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    EBAMRCellData& grad = storage->get_gradient();

    m_amr->computeGradient(grad, *a_phis[idx]);
    m_amr->averageDown(grad, m_cdr->get_phase());
    m_amr->interpGhost(grad, m_cdr->get_phase());
  }
}

void strang2::compute_errors(){
  CH_TIME("strang2::compute_errors");
  if(m_verbosity > 5){
    pout() << "strang2::compute_errors" << endl;
  }

  const Real safety = 1.E-20;
  const int comp = 0;
  Real max, min, emax, emin;

  // CDR errors
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    const int which = solver_it.get_solver();
    
    const EBAMRCellData& phi = solver->getPhi();
    EBAMRCellData& err       = storage->get_error();

    // So far 'err' contains the embedded formula and 'phi' is the numerical solution
    data_ops::incr(err, phi, -1.0); // err -> (err-phi), this is opposite, but the norm takes the magnitude anyways

    m_amr->averageDown(err, m_cdr->get_phase());
    
    Real Lerr, Lphi;
    data_ops::norm(Lerr, *err[0], m_amr->getDomains()[0], m_error_norm);
    data_ops::norm(Lphi, *phi[0], m_amr->getDomains()[0], m_error_norm);

    // Don't want to divide by zero...
    Lphi = Max(Lphi, safety);
    
    m_cdr_error[which] = Lerr/Lphi;

    // Done computing errors, make err into the embedded formula again
    data_ops::incr(err, phi, 1.0);
  }

  // Sigma error. So far, 'err' contains the embedded formula and 'phi' is the numerical solution
  EBAMRIVData& phi = m_sigma->getPhi();
  EBAMRIVData& err = m_sigma_scratch->get_error();
  data_ops::incr(err, phi, -1.0);
  data_ops::get_max_min_norm(max, min, phi);
  data_ops::get_max_min_norm(emax, emin, err);

  // DOn't want to divide by zero
  max = Max(max, safety);
  m_sigma_error = emax/max;

  // Done with sigma, make it the embedded formula again
  data_ops::incr(err, phi, 1.0);

  // Maximum error
  m_max_error = this->get_max_error();
}

void strang2::computeDt(Real& a_dt, TimeCode::which_code& a_timeCode){
  CH_TIME("strang2::computeDt");
  if(m_verbosity > 5){
    pout() << "strang2::computeDt" << endl;
  }

  Real dt = 1.E99;

  m_dt_cfl = m_cdr->compute_cfl_dt();
  const Real dt_cfl = m_cfl*m_dt_cfl;
  if(dt_cfl < dt){
    dt = dt_cfl;
    a_timeCode = TimeCode::cfl;
  }

  const Real dt_src = m_src_growth*m_cdr->computeSourceDt(m_src_tolerance, m_src_elec_only);
  if(dt_src < dt){
    dt = dt_src;
    a_timeCode = TimeCode::Source;
  }

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
}

void strang2::regridInternals(){
  CH_TIME("strang2::regridInternals");
  if(m_verbosity > 5){
    pout() << "strang2::regridInternals" << endl;
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
  
  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx] = RefCountedPtr<cdr_storage> (new cdr_storage(m_rk_order, m_amr, m_cdr->get_phase(), ncomp));
    m_cdr_scratch[idx]->allocate_storage();
  }
}

void strang2::allocate_poisson_storage(){
  const int ncomp = 1;
  m_fieldSolver_scratch = RefCountedPtr<poisson_storage> (new poisson_storage(m_rk_order, m_amr, m_cdr->get_phase(), ncomp));
  m_fieldSolver_scratch->allocate_storage();
}

void strang2::allocate_rte_storage(){
  const int ncomp       = 1;
  const int num_photons = m_plaskin->get_num_photons();
  m_rte_scratch.resize(num_photons);
  
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx] = RefCountedPtr<rte_storage> (new rte_storage(m_rk_order, m_amr, m_rte->get_phase(), ncomp));
    m_rte_scratch[idx]->allocate_storage();
  }
}

void strang2::allocate_sigma_storage(){
  const int ncomp = 1;
  m_sigma_scratch = RefCountedPtr<sigma_storage> (new sigma_storage(m_rk_order, m_amr, m_cdr->get_phase(), ncomp));
  m_sigma_scratch->allocate_storage();
}

void strang2::deallocateInternals(){
  CH_TIME("strang2::deallocateInternals");
  if(m_verbosity > 5){
    pout() << "strang2::deallocateInternals" << endl;
  }
  
  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
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

  m_fieldSolver_scratch->deallocate_storage();
  m_fieldSolver_scratch = RefCountedPtr<poisson_storage>(0);
  
  m_sigma_scratch->deallocate_storage();
  m_sigma_scratch = RefCountedPtr<sigma_storage>(0);
}

void strang2::backup_solutions(){
  CH_TIME("strang2::backup_solutions");
  if(m_verbosity > 5){
    pout() << "strang2::backup_solutions" << endl;
  }
  
  // Backup cdr solutions
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();

    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    EBAMRCellData& backup = storage->get_backup();

    data_ops::copy(backup, solver->getPhi());
  }

  {// Backup Poisson solution
    MFAMRCellData& backup = m_fieldSolver_scratch->get_backup();
    data_ops::copy(backup, m_fieldSolver->getPotential());
  }

  // Backup RTE solutions
  for (rte_iterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<rte_solver>& solver = solver_it();

    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);
    EBAMRCellData& backup = storage->get_backup();

    data_ops::copy(backup, solver->getPhi());
  }

  { // Backup sigma
    EBAMRIVData& backup = m_sigma_scratch->get_backup();
    data_ops::copy(backup, m_sigma->getPhi());
  }
}

void strang2::revert_backup(){
  CH_TIME("strang2::revert_backup");
  if(m_verbosity > 5){
    pout() << "strang2::revert_backup" << endl;
  }
  
  // Revert cdr solutions
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();

    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    const EBAMRCellData& backup = storage->get_backup();

    data_ops::copy(solver->getPhi(), backup);
  }

  {// Revert Poisson solution
    const MFAMRCellData& backup = m_fieldSolver_scratch->get_backup();
    data_ops::copy(m_fieldSolver->getPotential(), backup);
  }

  // Revert RTE solutions
  for (rte_iterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<rte_solver>& solver = solver_it();

    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);
    const EBAMRCellData& backup = storage->get_backup();

    data_ops::copy(solver->getPhi(), backup);
  }

  { // Revert sigma solution
    const EBAMRIVData& backup = m_sigma_scratch->get_backup();
    data_ops::copy(m_sigma->getPhi(), backup);
  }
}

void strang2::copy_error_to_solvers(){
  CH_TIME("strang2::copy_error_to_solvers");
  if(m_verbosity > 5){
    pout() << "strang2::copy_error_to_solvers" << endl;
  }

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& state      = solver->getPhi();
    const EBAMRCellData& err  = storage->get_error();

    data_ops::copy(state, err);
  }

  EBAMRIVData& phi       = m_sigma->getPhi();
  const EBAMRIVData& err = m_sigma_scratch->get_error();
  data_ops::set_value(phi, 0.0);
  data_ops::incr(phi, err, 1.0);
}

void strang2::copy_solvers_to_error(){
  CH_TIME("strang2::copy_error_to_solvers");
  if(m_verbosity > 5){
    pout() << "strang2::copy_error_to_solvers" << endl;
  }

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    const EBAMRCellData& state = solver->getPhi();
    EBAMRCellData& err   = storage->get_error();

    data_ops::copy(err, state);
  }

  const EBAMRIVData& phi = m_sigma->getPhi();
  EBAMRIVData& err = m_sigma_scratch->get_error();
  data_ops::set_value(err, 0.0);
  data_ops::incr(err, phi, 1.0);
}

void strang2::copy_solvers_to_cache(){
  CH_TIME("strang2::copy_solvers_to_cache");
  if(m_verbosity > 5){
    pout() << "strang2::copy_solvers_to_cache" << endl;
  }
  
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    const EBAMRCellData& state = solver->getPhi();
    EBAMRCellData& cache       = storage->get_cache();

    data_ops::copy(cache, state);
  }

  for (rte_iterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);

    const EBAMRCellData& state = solver->getPhi();
    EBAMRCellData& cache       = storage->get_cache();

    data_ops::copy(cache, state);
  }

  { // Cache the Poisson solution
    const MFAMRCellData& state = m_fieldSolver->getPotential();
    MFAMRCellData& cache       = m_fieldSolver_scratch->get_cache();
    data_ops::copy(cache, state);
  }

  { // Cache the surface charge density
    const EBAMRIVData& phi = m_sigma->getPhi();
    EBAMRIVData& cache = m_sigma_scratch->get_cache();
    data_ops::set_value(cache, 0.0);
    data_ops::incr(cache, phi, 1.0);
  }
}

void strang2::copy_cache_to_solvers(){
  CH_TIME("strang2::copy_cache_to_solvers");
  if(m_verbosity > 5){
    pout() << "strang2::copy_cache_to_solvers" << endl;
  }
  
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& state       = solver->getPhi();
    const EBAMRCellData& cache = storage->get_cache();

    data_ops::copy(state, cache);
  }

  for (rte_iterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);

    EBAMRCellData& state       = solver->getPhi();
    const EBAMRCellData& cache = storage->get_cache();

    data_ops::copy(state, cache);
  }

  { // Uncache the Poisson solution
    MFAMRCellData& state       = m_fieldSolver->getPotential();
    const MFAMRCellData& cache = m_fieldSolver_scratch->get_cache();
    data_ops::copy(state, cache);
  }

  { // Uncache the surface charge density
    EBAMRIVData& sigma       = m_sigma->getPhi();
    const EBAMRIVData& cache = m_sigma_scratch->get_cache();
    data_ops::set_value(sigma, 0.0);
    data_ops::incr(sigma, cache, 1.0);
  }
}

void strang2::compute_E_into_scratch(){
  CH_TIME("strang2::compute_E_into_scratch");
  if(m_verbosity > 5){
    pout() << "strang2::compute_E_into_scratch" << endl;
  }
  
  EBAMRCellData& E_cell = m_fieldSolver_scratch->get_E_cell();
  EBAMRFluxData& E_face = m_fieldSolver_scratch->get_E_face();
  EBAMRIVData&   E_eb   = m_fieldSolver_scratch->get_E_eb();
  EBAMRIFData&   E_dom  = m_fieldSolver_scratch->get_E_domain();

  const MFAMRCellData& phi = m_fieldSolver->getPotential();
  
  this->compute_E(E_cell, m_cdr->get_phase(), phi);     // Compute cell-centered field
  this->compute_E(E_face, m_cdr->get_phase(), E_cell);  // Compute face-centered field
  this->compute_E(E_eb,   m_cdr->get_phase(), E_cell);  // EB-centered field

  TimeStepper::extrapolate_to_domain_faces(E_dom, m_cdr->get_phase(), E_cell);
}

void strang2::compute_cdr_velo(const Real a_time){
  CH_TIME("strang2::compute_cdr_velo");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_velo" << endl;
  }

  this->compute_cdr_velo(m_cdr->getPhis(), a_time);
}

void strang2::compute_cdr_velo(const Vector<EBAMRCellData*>& a_phis, const Real a_time){
  CH_TIME("strang2::compute_cdr_velo(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_velo(Vector<EBAMRCellData*>, Real)" << endl;
  }

  Vector<EBAMRCellData*> velocities = m_cdr->get_velocities();
  this->compute_cdr_velocities(velocities, a_phis, m_fieldSolver_scratch->get_E_cell(), a_time);
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
  
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    cdr_states.push_back(&(solver->getPhi()));
    eb_states.push_back(&(storage->get_eb_state()));
    eb_gradients.push_back(&(storage->get_eb_grad()));
    cdr_gradients.push_back(&(storage->get_gradient())); // Should already have been computed
  }

  // Extrapolate states to the EB and floor them so we cannot get negative values on the boundary. This
  // won't hurt mass conservation because the mass hasn't been injected yet
  this->extrapolate_to_eb(eb_states, m_cdr->get_phase(), cdr_states);
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
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

void strang2::compute_cdr_eb_states(const Vector<EBAMRCellData*>& a_phis){
  CH_TIME("strang2::compute_cdr_eb_states(vec)");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_eb_states(vec)" << endl;
  }

  Vector<EBAMRIVData*>   eb_gradients;
  Vector<EBAMRIVData*>   eb_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    eb_states.push_back(&(storage->get_eb_state()));
    eb_gradients.push_back(&(storage->get_eb_grad()));
    cdr_gradients.push_back(&(storage->get_gradient())); // Should already have been computed
  }

  // Extrapolate states to the EB and floor them so we cannot get negative values on the boundary. This
  // won't hurt mass conservation because the mass hasn't been injected yet
  this->extrapolate_to_eb(eb_states, m_cdr->get_phase(), a_phis);
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    data_ops::floor(*eb_states[idx], 0.0);
  }

  // We should already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIVData eb_gradient;
  m_amr->allocate(eb_gradient, m_cdr->get_phase(), SpaceDim);
  for (int i = 0; i < a_phis.size(); i++){
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
  
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    cdr_states.push_back(&(solver->getPhi()));
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

void strang2::compute_cdr_domain_states(const Vector<EBAMRCellData*>& a_phis){
  CH_TIME("strang2::compute_cdr_domain_states");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_domain_states" << endl;
  }

  Vector<EBAMRIFData*>   domain_gradients;
  Vector<EBAMRIFData*>   domain_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    domain_states.push_back(&(storage->get_domain_state()));
    domain_gradients.push_back(&(storage->get_domain_grad()));
    cdr_gradients.push_back(&(storage->get_gradient()));
  }

  // Extrapolate states to the domain faces
  this->extrapolate_to_domain_faces(domain_states, m_cdr->get_phase(), a_phis);

  // We already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIFData grad;
  m_amr->allocate(grad, m_cdr->get_phase(), SpaceDim);
  for (int i = 0; i < a_phis.size(); i++){
    this->extrapolate_to_domain_faces(grad, m_cdr->get_phase(), *cdr_gradients[i]);
    this->project_domain(*domain_gradients[i], grad);
  }
}

void strang2::compute_cdr_fluxes(const Real a_time){
  CH_TIME("strang2::compute_cdr_fluxes");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_fluxes" << endl;
  }

  this->compute_cdr_fluxes(m_cdr->getPhis(), a_time);
}

void strang2::compute_cdr_fluxes(const Vector<EBAMRCellData*>& a_phis, const Real a_time){
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
  Vector<EBAMRCellData*> cdr_velocities = m_cdr->get_velocities();
  this->compute_extrapolated_fluxes(extrap_cdr_fluxes, a_phis, cdr_velocities, m_cdr->get_phase());
  this->extrapolate_to_eb(extrap_cdr_velocities, m_cdr->get_phase(), cdr_velocities);

  // Compute RTE flux on the boundary
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
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

void strang2::compute_cdr_domain_fluxes(const Real a_time){
  CH_TIME("strang2::compute_cdr_domain_fluxes");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_domain_fluxes" << endl;
  }

  this->compute_cdr_domain_fluxes(m_cdr->getPhis(), a_time);
}

void strang2::compute_cdr_domain_fluxes(const Vector<EBAMRCellData*>& a_phis, const Real a_time){
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

  cdr_fluxes = m_cdr->getDomainFlux();
  cdr_velocities = m_cdr->get_velocities();
  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
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
  this->extrapolate_to_domain_faces(extrap_cdr_densities,         m_cdr->get_phase(), a_phis);
  this->extrapolate_vector_to_domain_faces(extrap_cdr_velocities, m_cdr->get_phase(), cdr_velocities);
  this->compute_extrapolated_domain_fluxes(extrap_cdr_fluxes,     a_phis,           cdr_velocities, m_cdr->get_phase());
  this->extrapolate_vector_to_domain_faces(extrap_cdr_gradients,  m_cdr->get_phase(), cdr_gradients);

  // Compute RTE flux on domain faces
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);

    EBAMRIFData& domain_flux = storage->get_domain_flux();
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

void strang2::compute_sigma_flux(){
  CH_TIME("strang2::compute_sigma_flux");
  if(m_verbosity > 5){
    pout() << "strang2::compute_sigma_flux" << endl;
  }

  EBAMRIVData& flux = m_sigma->getFlux();
  data_ops::set_value(flux, 0.0);

  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const RefCountedPtr<species>& spec      = solver_it.getSpecies();
    const EBAMRIVData& solver_flux          = solver->getEbFlux();

    data_ops::incr(flux, solver_flux, spec->getChargeNumber()*units::s_Qe);
  }

  m_sigma->resetCells(flux);
}

void strang2::compute_cdr_sources(const Real a_time){
  CH_TIME("strang2::compute_cdr_sources_into_scratch");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_sources_into_scratch" << endl;
  }

  this->compute_cdr_sources(m_cdr->getPhis(), a_time);
}

void strang2::compute_cdr_sources(const Vector<EBAMRCellData*>& a_phis, const Real a_time){
  CH_TIME("strang2::compute_cdr_sources(Vector<EBAMRCellData*>, Real)");
  if(m_verbosity > 5){
    pout() << "strang2::compute_cdr_sources(Vector<EBAMRCellData*>, Real)" << endl;
  }
  
  Vector<EBAMRCellData*> cdr_sources = m_cdr->getSources();
  Vector<EBAMRCellData*> rte_states  = m_rte->getPhis();
  EBAMRCellData& E                   = m_fieldSolver_scratch->get_E_cell();

  Vector<EBAMRCellData*> cdr_gradients;
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    cdr_gradients.push_back(&(storage->get_gradient())); // These should already have been computed
  }

  TimeStepper::compute_cdr_sources(cdr_sources, a_phis, cdr_gradients, rte_states, E, a_time, centering::cell_center);
}

void strang2::advance_rte_stationary(const Real a_time){
  CH_TIME("strang2::compute_rte_k1_stationary");
  if(m_verbosity > 5){
    pout() << "strang2::compute_k1_stationary" << endl;
  }

  if((m_timeStep + 1) % m_fast_rte == 0){
    Vector<EBAMRCellData*> rte_states  = m_rte->getPhis();
    Vector<EBAMRCellData*> rte_sources = m_rte->getSources();
    Vector<EBAMRCellData*> cdr_states  = m_cdr->getPhis();

    EBAMRCellData& E = m_fieldSolver_scratch->get_E_cell();

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
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& state = solver->getPhi();
    EBAMRCellData& prev  = storage->get_previous();

    data_ops::copy(prev, state);
  }

  // Copy solver state into storage->m_previous
  EBAMRIVData& phi = m_sigma->getPhi();
  EBAMRIVData& pre = m_sigma_scratch->get_previous();
  data_ops::set_value(pre, 0.0);
  data_ops::incr(pre, phi, 1.0);
}

void strang2::restore_solvers(){
  CH_TIME("strang2::restore_solvers");
  if(m_verbosity > 5){
    pout() << "strang2::restore_solvers" << endl;
  }

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& state = solver->getPhi();
    EBAMRCellData& prev  = storage->get_previous();

    data_ops::copy(state, prev);
  }

  EBAMRIVData& phi = m_sigma->getPhi();
  EBAMRIVData& pre = m_sigma_scratch->get_previous();
  data_ops::set_value(phi, 0.0);
  data_ops::incr(phi, pre, 1.0);
}

void strang2::update_poisson(){
  if(m_verbosity > 5){
    pout() << "strang2::update_poisson" << endl;
  }
  
  if(m_do_poisson){ // Solve Poisson equation
    if((m_timeStep +1) % m_fast_poisson == 0){
      TimeStepper::solve_poisson();
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
#include "CD_NamespaceFooter.H"
