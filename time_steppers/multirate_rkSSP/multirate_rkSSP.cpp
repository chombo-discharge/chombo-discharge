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
  m_order  = 2;
  m_maxCFL = 0.5;

  if(procID() == 0){
    std::cout << "multirate_rkSSP::multirate_rkSSP - this class might be in error. Talk to Hans Kristian about this" << std::endl;
  }

  // Basically only for debugging
  m_print_diagno = false;
  m_write_diagno = false;
  m_do_advec_src = true;
  m_do_diffusion = true;
  m_do_rte       = true;
  m_do_poisson   = true;
  m_compute_v    = true;
  m_compute_S    = true;
  m_compute_err  = true;
  
  {
    ParmParse pp("multirate_rkSSP");

    std::string str;

    pp.get("order",   m_order);
    pp.get("max_cfl", m_maxCFL);

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
    if(pp.contains("compute_error")){
      pp.get("compute_error", str);
      if(str == "false"){
	m_compute_err = false;
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
    m_compute_err = false;
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

Real multirate_rkSSP::advance(const Real a_dt){
  CH_TIME("multirate_rkSSP::advance");
  if(m_verbosity > 2){
    pout() << "multirate_rkSSP::advance" << endl;
  }

  const Real t0 = MPI_Wtime();

  // When we enter this routine, we assume that velocities, source terms, and diffusion coefficients have already
  // been computed. If you change the order of the integration, you must ensure that appropriate velocities, diffusion
  // coefficients, and source terms, are updated where they should. So please don't do that...

  this->cache_solutions(); // Cache old solutions. We can safely manipulate directly into solver states.

  int substeps; // Number of substeps
  Real cfl;     // Substep CFL
  Real sub_dt;  // Substep time


  if(m_do_advec_src){
    substeps = ceil(a_dt/(m_maxCFL*m_dt_cfl));
    cfl      = a_dt/(substeps*m_dt_cfl);
    sub_dt   = cfl*m_dt_cfl;

    this->advance_multirate_advec_src(substeps, sub_dt);
  }
  const Real t1 = MPI_Wtime();
  
  if(m_do_diffusion){
    this->advance_diffusion(a_dt);
  }
  const Real t2 = MPI_Wtime();

  if(m_do_poisson){ // Solve Poisson equation
    if((m_step +1) % m_fast_poisson == 0){
      time_stepper::solve_poisson();
      this->compute_E_into_scratch();
    }
  }
  const Real t3 = MPI_Wtime();
  
  if(m_do_rte){ // Solve for final RTE stage. Poisson equation should already have been solved at the end of the advection
    this->advance_rte_stationary(m_time + a_dt); // and diffusion stages
  }
  const Real t4 = MPI_Wtime();
  // Put cdr solvers back in useable state so that we can reliably compute the next time step.

  multirate_rkSSP::compute_cdr_velo(m_time + a_dt);
  time_stepper::compute_cdr_diffusion(m_poisson_scratch->get_E_cell(), m_poisson_scratch->get_E_eb());
  multirate_rkSSP::compute_cdr_sources(m_time + a_dt);
  const Real t5 = MPI_Wtime();

  if(m_print_diagno){
    pout() << "\n";
    pout() << "\t multirate_rkSSP::advance(Real a_dt) breakdown" << endl;
    pout() << "\t ==============================================\n" << endl;
    
    pout() << "\t Convection-source advance: \n ";
    pout() << "\t --------------------------\n";
    pout() << "\t\t Steps     = " << substeps << endl;
    pout() << "\t\t Local cfl = " << cfl      << endl;
    pout() << "\t\t Local dt  = " << sub_dt   << endl;
    pout() << endl;
    pout() << "\t\t Time      = " << 100.*(t1-t0)/(t5-t0) << "%\n" << endl;

    pout() << "\n";
    pout() << "\t Diffusion advance:\n ";
    pout() << "\t -----------------\n";
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
    this->write_diagnostics(substeps, sub_dt, a_dt, cfl, t1-t0, t2-t1, t3-t2, t4-t3, t5-t4, t5-t0);
  }

  
  return a_dt;
}

void multirate_rkSSP::advance_diffusion(const Real a_dt){
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
    m_cdr->set_source(0.0); // This is necessary because advance_diffusion also works with source terms
    this->compute_cdr_diffusion(m_poisson_scratch->get_E_cell(), m_poisson_scratch->get_E_eb());
    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver>& solver = solver_it();

      solver->advance_diffusion(a_dt);
    }
  }
}

void multirate_rkSSP::advance_multirate_advec_src(const int a_substeps, const Real a_dt){
  CH_TIME("multirate_rkSSP::advance_multirate_advec_src");
  if(m_verbosity > 2){
    pout() << "multirate_rkSSP::advance_multirate_advec_src" << endl;
  }

  this->compute_E_into_scratch(); 

  if(m_order == 1){
    this->advance_advec_rk1(a_substeps, a_dt);
  }
  else if(m_order == 2){
    this->advance_advec_rk2(a_substeps, a_dt);
  }
  else if(m_order == 3){
    this->advance_advec_rk3(a_substeps, a_dt);
  }
}

void multirate_rkSSP::advance_advec_rk1(const int a_substeps, const Real a_dt){
  CH_TIME("multirate_rkSSP::advance_advec_rk1");
  if(m_verbosity > 2){
    pout() << "multirate_rkSSP::advance_advec_rk1" << endl;
  }

  for (int step = 0; step < a_substeps; step++){
    const Real time = m_time + step*a_dt;
    this->compute_cdr_eb_states();
    this->compute_cdr_fluxes(time);
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
    
    // Update for next iterate. This should generally be done because the solution might move many grid cells.
    if(step < a_substeps -1){ // No need on last step
      if(m_compute_v) this->compute_cdr_velo(time);
      if(m_compute_S) this->compute_cdr_sources(time);
    }
  }
}

void multirate_rkSSP::advance_advec_rk2(const int a_substeps, const Real a_dt){
  CH_TIME("multirate_rkSSP::advance_advec_rk2");
  if(m_verbosity > 2){
    pout() << "multirate_rkSSP::advance_advec_rk2" << endl;
  }

  for (int step = 0; step < a_substeps; step++){
    const Real time = m_time + step*a_dt;

    // u^1 = u^n + dt*L(u^n)
    // u^n resides in cache, put u^1 in the solver. Use scratch for computing L(u^n)
    {
      // Compute fluxes using the solver state
      this->compute_cdr_eb_states();
      this->compute_cdr_fluxes(time);
      this->compute_sigma_flux();
      
      for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	RefCountedPtr<cdr_solver>& solver   = solver_it();
	RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

	EBAMRCellData& phi  = solver->get_state();
	EBAMRCellData& rhs  = storage->get_scratch();
	EBAMRCellData& src  = solver->get_source();
	EBAMRCellData& pre  = storage->get_previous();

	data_ops::copy(pre, phi); // Store u^n
	
	solver->compute_divF(rhs, phi, 0.0, true); // RHS =  div(u*v)
	data_ops::scale(rhs, -1.0);                // RHS = -div(u*v)
	data_ops::incr(rhs, src, 1.0);             // RHS = S - div(u*v)
	data_ops::incr(phi, rhs, a_dt);            // u^(n+1) = u^n + dt*L(u^n)

	m_amr->average_down(phi, m_cdr->get_phase());
	m_amr->interp_ghost(phi, m_cdr->get_phase());

	if(m_compute_err){
	  EBAMRCellData& err = storage->get_error();

	  data_ops::copy(err, phi);   // err =  (u^n + dt*L(u^n))
	  data_ops::scale(err, -1.0); // err = -(u^n + dt*L(u^n))
	}
      }

      // Advance sigma
      EBAMRIVData& phi = m_sigma->get_state();
      EBAMRIVData& rhs = m_sigma_scratch->get_scratch();
      EBAMRIVData& pre = m_sigma_scratch->get_previous();

      data_ops::copy(pre, phi); // Store u^n before doing anything
      
      m_sigma->compute_rhs(rhs);
      m_sigma->reset_cells(rhs);
      data_ops::incr(phi, rhs, a_dt); // phi = u^1 = u^n + dt*L(u^n)

      m_amr->average_down(phi, m_cdr->get_phase());
      m_sigma->reset_cells(phi);

      if(m_compute_err){
	EBAMRIVData& err = m_sigma_scratch->get_error();
	data_ops::set_value(err, 0.0);
	data_ops::incr(err, phi, -1.0);  // error = -(u^n + dt*L(u^n))
      }
    }

    // u^(n+1) = 0.5*(u^n + u^1 + dt*L(u^(1)))
    // u^n resides in temp storage, u^1 resides in solver. Put u^2 in solver at end and use scratch fo computing L(u^1)
    {

      // Compute fluxes from u^1, which resides in solvers. 
      this->compute_cdr_eb_states();
      this->compute_cdr_fluxes(time);
      this->compute_sigma_flux();

      if(m_compute_v) this->compute_cdr_velo(time);
      if(m_compute_S) this->compute_cdr_sources(time);

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

	if(m_compute_err){
	  EBAMRCellData& err = storage->get_error(); // err  = -(u^n + dt*L(u^n)) // 1st. order
	  data_ops::incr(err, phi, 1.0);             // err += (second order)

	  const int comp = 0;
	  Real max, min;
	  data_ops::get_max_min(max, min, phi, 0);
	  
	  Real emax, emin;
	  data_ops::get_max_min_norm(emax, emin, err);
	  emax = emax/max;

	  const int which = solver_it.get_solver();
	  m_cdr_error[which] = emax;
	}
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

      if(m_compute_err){
	EBAMRIVData& err = m_sigma_scratch->get_error(); // err = -(u^n + dt*L(u^n))
	data_ops::incr(err, phi, 1.0);                   // err = phi - (u^n + dt*L(u^n))

	m_sigma->reset_cells(err);

	const int comp = 0;
	Real max, min;
	data_ops::get_max_min_norm(max, min, phi);

	Real emax, emin;
	data_ops::get_max_min_norm(emax, emin, err);
	emax = emax/max;

	m_sigma_error = emax;
      }
    }
    
    // Update for next iterate. This should generally be done because the solution might move many grid cells.
    if(step < a_substeps -1){ // No need on last step
      if(m_compute_v) this->compute_cdr_velo(time);
      if(m_compute_S) this->compute_cdr_sources(time);
    }

#if 1 // Debug
    if(procID() == 0){
      std::cout << "cdr_err = " << m_cdr_error << "\t sigma_err = " << m_sigma_error << std::endl;
    }
#endif

  }
}

void multirate_rkSSP::advance_advec_rk3(const int a_substeps, const Real a_dt){
  CH_TIME("multirate_rkSSP::advance_advec_rk3");
  if(m_verbosity > 2){
    pout() << "multirate_rkSSP::advance_advec_rk3" << endl;
  }

#if 1 // Debug warning 
  if(procID() == 0){
    std::cout << "multirate_rkSSP::advance_advec_rk3 - There are indications of implementation errors in this routine; the computed error is larger than for the second order embedded scheme" << std::endl;
  }
#endif

  for (int step = 0; step < a_substeps; step++){
    const Real time = m_time + step*a_dt;

    // u^1 = u^n + dt*L(u^n)
    // u^n resides in solver. Store this and put u^1 in the solver instead. Use scratch for computing L(u^n)
    {

      // Compute fluxes from u^n
      this->compute_cdr_eb_states();
      this->compute_cdr_fluxes(time);
      this->compute_sigma_flux();

      for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
	RefCountedPtr<cdr_solver>& solver   = solver_it();
	RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

	EBAMRCellData& phi        = solver->get_state();
	EBAMRCellData& rhs        = storage->get_scratch();
	EBAMRCellData& pre        = storage->get_previous();
	const EBAMRCellData& src  = solver->get_source();

	data_ops::copy(pre, phi); // Store u^n
	
	solver->compute_divF(rhs, phi, 0.0, true); // RHS =  div(u*v)
	data_ops::scale(rhs, -1.0);                // RHS = -div(u*v)
	data_ops::incr(rhs, src, 1.0);             // RHS = S - div(u*v)
	data_ops::incr(phi, rhs, a_dt);            // u^(n+1) = u^n + dt*L(u^n)

	m_amr->average_down(phi, m_cdr->get_phase());
	m_amr->interp_ghost(phi, m_cdr->get_phase());

	if(m_compute_err){
	  EBAMRCellData& err = storage->get_error();

	  data_ops::copy(err, pre);      // err = u^n
	  data_ops::incr(err, phi, 1.0); // err = u^n + u^1
	}
      }

      // Advance sigma
      EBAMRIVData& phi = m_sigma->get_state();
      EBAMRIVData& rhs = m_sigma_scratch->get_scratch();
      EBAMRIVData& pre = m_sigma_scratch->get_previous();

      data_ops::copy(pre, phi); // Store u^n
      
      m_sigma->compute_rhs(rhs);
      m_sigma->reset_cells(rhs);
      data_ops::incr(phi, rhs, a_dt); // phi = u^n + dt*L(u^n)

      m_amr->average_down(phi, m_cdr->get_phase());
      m_sigma->reset_cells(phi);

      if(m_compute_err){
	EBAMRIVData& err = m_sigma_scratch->get_error();

	data_ops::set_value(err, 0.0);
	data_ops::incr(err, pre, 1.0); // err = u^n
	data_ops::incr(err, phi, 1.0); // err = u^n + u^1
      }
    }

    // u^2 = (3*u^n + u^1 + dt*L(u^(1)))/4
    // u^n is stored, u^1 resides in solver. Put u^2 in solver at end and use scratch fo computing L(u^1)
    {
      // Compute fluxes from u^1, which reside in solver
      this->compute_cdr_eb_states();
      this->compute_cdr_fluxes(time);
      this->compute_sigma_flux();

      if(m_compute_v) this->compute_cdr_velo(time);
      if(m_compute_S) this->compute_cdr_sources(time);

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

	if(m_compute_err){
	  EBAMRCellData& err = storage->get_error(); // err =  (u^n + u^1)
	  data_ops::incr(err, rhs, a_dt);            // err =  (u^n + u^1 + dt*L(u^1))
	  data_ops::scale(err, -0.5);                // err = -(u^n + u^1 + dt*L(u^1))/2
	}
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

      if(m_compute_err){
	EBAMRIVData& err = m_sigma_scratch->get_error(); // err =  (u^n + u^1)
	data_ops::incr(err, rhs, a_dt);                  // err =  (u^n + u^1 + dt*L(u^1))
	data_ops::scale(err, -0.5);                      // err =  (u^n + u^1 + dt*L(u^1))/2
      }
    }

    // u^(n+1) = (u^n + 2*u^2 + 2*dt*L(u^(1)))/3
    // u^n is stored, u^2 resides in solver. Put u^(n+1) in solver at end and use scratch for computing L(u^1)
    {
      this->compute_cdr_eb_states();
      this->compute_cdr_fluxes(time);
      this->compute_sigma_flux();

      if(m_compute_v) this->compute_cdr_velo(time);
      if(m_compute_S) this->compute_cdr_sources(time);

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

	if(m_compute_err){
	  EBAMRCellData& err = storage->get_error(); // err  = -0.5*(u^n + u^1 + dt*L(u^1)) // second order
	  data_ops::incr(err, phi, 1.0);             // err += (third order approximation)

	  const int comp = 0;
	  Real max, min;
	  data_ops::get_max_min(max, min, phi, 0);
	  
	  Real emax, emin;
	  data_ops::get_max_min_norm(emax, emin, err);
	  emax = emax/max;

	  const int which = solver_it.get_solver();
	  m_cdr_error[which] = emax;
	}
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

      if(m_compute_err){
	EBAMRIVData& err = m_sigma_scratch->get_error(); // err = -(u^n + dt*L(u^n))
	data_ops::incr(err, phi, 1.0);                   // err = phi - (u^n + dt*L(u^n))

	m_sigma->reset_cells(err);

	const int comp = 0;
	Real max, min;
	data_ops::get_max_min_norm(max, min, phi);

	Real emax, emin;
	data_ops::get_max_min_norm(emax, emin, err);
	emax = emax/max;

	m_sigma_error = emax;
      }
    }
    
    // Update for next iterate. This should generally be done because the solution might move many grid cells.
    if(step < a_substeps -1){ // No need on last step
      if(m_compute_v) this->compute_cdr_velo(time);
      if(m_compute_S) this->compute_cdr_sources(time);
    }

#if 1 // Debug
    if(procID() == 0){
      std::cout << "cdr_err = " << m_cdr_error << "\t sigma_err = " << m_sigma_error << std::endl;
    }
#endif
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

  const MFAMRCellData& phi = m_poisson->get_state();
  
  this->compute_E(E_cell, m_cdr->get_phase(), phi);     // Compute cell-centered field
  this->compute_E(E_face, m_cdr->get_phase(), E_cell);  // Compute face-centered field
  this->compute_E(E_eb,   m_cdr->get_phase(), E_cell);  // EB-centered field
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
  
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    cdr_states.push_back(&(solver->get_state()));
    eb_states.push_back(&(storage->get_eb_state()));
    eb_gradients.push_back(&(storage->get_eb_grad()));
  }

  this->extrapolate_to_eb(eb_states,          m_cdr->get_phase(), cdr_states);
  this->compute_gradients_at_eb(eb_gradients, m_cdr->get_phase(), cdr_states);
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

  time_stepper::compute_cdr_sources(cdr_sources, a_states, rte_states, E, a_time, centering::cell_center);
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
