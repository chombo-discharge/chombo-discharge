/*!
  @file   multirate_eulerf.cpp
  @brief  Implementation of multirate_eulerf.H
  @author Robert Marskar
  @date   Sept. 2018
*/


#include "multirate_eulerf.H"
#include "multirate_eulerf_storage.H"
#include "data_ops.H"
#include "units.H"

#include <ParmParse.H>

typedef multirate_eulerf::cdr_storage     cdr_storage;
typedef multirate_eulerf::poisson_storage poisson_storage;
typedef multirate_eulerf::rte_storage     rte_storage;
typedef multirate_eulerf::sigma_storage   sigma_storage;

multirate_eulerf::multirate_eulerf(){
  m_maxCFL = 0.5;

  // Basically only for debugging
  m_diagnostics  = false;
  m_do_advec_src = true;
  m_do_diffusion = true;
  m_do_rte       = true;
  m_do_poisson   = true;
  
  {
    ParmParse pp("multirate_eulerf");

    std::string str;
    
    pp.get("max_cfl", m_maxCFL);

    if(pp.contains("diagnostics")){
      pp.get("diagnostics", str);
      if(str == "true"){
	m_diagnostics = true;
      }
    }

    if(pp.contains("turn_off_advection")){
      pp.get("turn_off_advection_source", str);
      if(str == "true"){
	m_do_advec_src = false;
	if(m_verbosity > 2){
	  pout() << "multirate_eulerf::multirate_eulerf - Turning off advection & source" << endl;
	}
      }
    }
    if(pp.contains("turn_off_diffusion")){
      pp.get("turn_off_diffusion", str);
      if(str == "true"){
	m_do_diffusion = false;
	if(m_verbosity > 2){
	  pout() << "multirate_eulerf::multirate_eulerf - Turning off diffusion" << endl;
	}
      }
    }
    if(pp.contains("turn_off_rte")){
      pp.get("turn_off_rte", str);
      if(str == "true"){
	m_do_rte = false;

	if(m_verbosity > 2){
	  pout() << "multirate_eulerf::multirate_eulerf - Turning off rte" << endl;
	}
      }
    }
    if(pp.contains("turn_off_poisson")){
      pp.get("turn_off_poisson", str);
      if(str == "true"){
	m_do_poisson = false;

	if(m_verbosity > 2){
	  pout() << "multirate_eulerf::multirate_eulerf - Turning off poisson" << endl;
	}
      }
    }
  }
}

multirate_eulerf::~multirate_eulerf(){

}

RefCountedPtr<cdr_storage>& multirate_eulerf::get_cdr_storage(const cdr_iterator& a_solverit){
  return m_cdr_scratch[a_solverit.get_solver()];
}

RefCountedPtr<rte_storage>& multirate_eulerf::get_rte_storage(const rte_iterator& a_solverit){
  return m_rte_scratch[a_solverit.get_solver()];
}

Real multirate_eulerf::restrict_dt(){
  return 1.E99;
}

Real multirate_eulerf::advance(const Real a_dt){
  CH_TIME("multirate_eulerf::advance");
  if(m_verbosity > 2){
    pout() << "multirate_eulerf::advance" << endl;
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

  multirate_eulerf::compute_cdr_velo(m_time + a_dt);
  time_stepper::compute_cdr_diffusion(m_poisson_scratch->get_E_cell(), m_poisson_scratch->get_E_eb());
  multirate_eulerf::compute_cdr_sources(m_time + a_dt);
  const Real t5 = MPI_Wtime();

  if(m_diagnostics){
    pout() << "\t multirate_eulerf::advance(Real a_dt) breakdown" << endl;
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
  
  return a_dt;
}

void multirate_eulerf::advance_diffusion(const Real a_dt){
  CH_TIME("multirate_eulerf::advance_diffusion");
  if(m_verbosity > 5){
    pout() << "multirate_eulerf::advance_diffusion" << endl;
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

void multirate_eulerf::advance_multirate_advec_src(const int a_substeps, const Real a_dt){
  CH_TIME("multirate_eulerf::advance_multirate_advec_src");
  if(m_verbosity > 2){
    pout() << "multirate_eulerf::advance_multirate_advec_src" << endl;
  }

  this->compute_E_into_scratch(); 
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
      solver->compute_divF(rhs, phi, 0.0, true);
      data_ops::scale(rhs, -1.0);
      data_ops::incr(rhs, src, 1.0);
      data_ops::incr(phi, rhs, a_dt);

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
    this->compute_cdr_velo(time); 
    this->compute_cdr_sources(time);
  }
}

void multirate_eulerf::regrid_internals(){
  CH_TIME("time_stepper::regrid_internals");
  if(m_verbosity > 5){
    pout() << "time_stepper::regrid_internals" << endl;
  }
  
  this->allocate_cdr_storage();
  this->allocate_poisson_storage();
  this->allocate_rte_storage();
  this->allocate_sigma_storage();
}

void multirate_eulerf::compute_dt(Real& a_dt, time_code::which_code& a_timecode){
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
}

void multirate_eulerf::allocate_cdr_storage(){
  const int ncomp       = 1;
  const int num_species = m_plaskin->get_num_species();
  m_cdr_scratch.resize(num_species);
  
  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx] = RefCountedPtr<cdr_storage> (new cdr_storage(m_amr, m_cdr->get_phase(), ncomp));
    m_cdr_scratch[idx]->allocate_storage();
  }
}

void multirate_eulerf::allocate_poisson_storage(){
  const int ncomp = 1;
  m_poisson_scratch = RefCountedPtr<poisson_storage> (new poisson_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_poisson_scratch->allocate_storage();
}

void multirate_eulerf::allocate_rte_storage(){
  const int ncomp       = 1;
  const int num_photons = m_plaskin->get_num_photons();
  m_rte_scratch.resize(num_photons);
  
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx] = RefCountedPtr<rte_storage> (new rte_storage(m_amr, m_rte->get_phase(), ncomp));
    m_rte_scratch[idx]->allocate_storage();
  }
}

void multirate_eulerf::allocate_sigma_storage(){
  const int ncomp = 1;
  m_sigma_scratch = RefCountedPtr<sigma_storage> (new sigma_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_sigma_scratch->allocate_storage();
}

void multirate_eulerf::deallocate_internals(){
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

void multirate_eulerf::cache_solutions(){
  CH_TIME("multirate_eulerf::cache_solutions");
  if(m_verbosity > 5){
    pout() << "multirate_eulerf::cache_solutions" << endl;
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

void multirate_eulerf::compute_E_into_scratch(){
  CH_TIME("multirate_eulerf::compute_E_into_scratch");
  if(m_verbosity > 5){
    pout() << "multirate_eulerf::compute_E_into_scratch" << endl;
  }
  
  EBAMRCellData& E_cell = m_poisson_scratch->get_E_cell();
  EBAMRFluxData& E_face = m_poisson_scratch->get_E_face();
  EBAMRIVData&   E_eb   = m_poisson_scratch->get_E_eb();

  const MFAMRCellData& phi = m_poisson->get_state();
  
  this->compute_E(E_cell, m_cdr->get_phase(), phi);     // Compute cell-centered field
  this->compute_E(E_face, m_cdr->get_phase(), E_cell);  // Compute face-centered field
  this->compute_E(E_eb,   m_cdr->get_phase(), E_cell);  // EB-centered field
}

void multirate_eulerf::compute_cdr_velo(const Real a_time){
  CH_TIME("multirate_eulerf::compute_cdr_velo");
  if(m_verbosity > 5){
    pout() << "multirate_eulerf::compute_cdr_velo" << endl;
  }

  Vector<EBAMRCellData*> states     = m_cdr->get_states();
  Vector<EBAMRCellData*> velocities = m_cdr->get_velocities();
  this->compute_cdr_velocities(velocities, states, m_poisson_scratch->get_E_cell(), a_time);
}

void multirate_eulerf::compute_cdr_eb_states(){
  CH_TIME("multirate_eulerf::compute_cdr_eb_states");
  if(m_verbosity > 5){
    pout() << "multirate_eulerf::compute_cdr_eb_states" << endl;
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

void multirate_eulerf::compute_cdr_fluxes(const Real a_time){
  CH_TIME("multirate_eulerf::compute_cdr_fluxes");
  if(m_verbosity > 5){
    pout() << "multirate_eulerf::compute_cdr_fluxes" << endl;
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
  Vector<EBAMRCellData*> cdr_densities  = m_cdr->get_states();
  Vector<EBAMRCellData*> cdr_velocities = m_cdr->get_velocities();
  this->compute_extrapolated_fluxes(extrap_cdr_fluxes, cdr_densities, cdr_velocities, m_cdr->get_phase());
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

void multirate_eulerf::compute_sigma_flux(){
  CH_TIME("multirate_eulerf::compute_sigma_flux");
  if(m_verbosity > 5){
    pout() << "multirate_eulerf::compute_sigma_flux" << endl;
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


void multirate_eulerf::compute_cdr_sources(const Real a_time){
  CH_TIME("multirate_eulerf::compute_cdr_sources_into_scratch");
  if(m_verbosity > 5){
    pout() << "multirate_eulerf::compute_cdr_sources_into_scratch" << endl;
  }
  
  Vector<EBAMRCellData*> cdr_sources = m_cdr->get_sources();
  Vector<EBAMRCellData*> cdr_states  = m_cdr->get_states();
  Vector<EBAMRCellData*> rte_states  = m_rte->get_states();
  EBAMRCellData& E                   = m_poisson_scratch->get_E_cell();

  time_stepper::compute_cdr_sources(cdr_sources, cdr_states, rte_states, E, a_time, centering::cell_center);
}

void multirate_eulerf::advance_rte_stationary(const Real a_time){
  CH_TIME("multirate_eulerf::compute_rte_k1_stationary");
  if(m_verbosity > 5){
    pout() << "multirate_eulerf::compute_k1_stationary" << endl;
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
