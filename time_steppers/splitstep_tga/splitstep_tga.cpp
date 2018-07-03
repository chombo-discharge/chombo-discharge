/*!
  @file   splitstep_tga.cpp
  @brief  Implementation of splitstep_tga.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "splitstep_tga.H"
#include "splitstep_tga_storage.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"
#include "data_ops.H"
#include "units.H"

#include <ParmParse.H>

typedef splitstep_tga::cdr_storage     cdr_storage;
typedef splitstep_tga::poisson_storage poisson_storage;
typedef splitstep_tga::rte_storage     rte_storage;
typedef splitstep_tga::sigma_storage   sigma_storage;

splitstep_tga::splitstep_tga(){
  m_alpha = 1.0;
  m_simpi = false;

  // Basically only for debugging
  m_do_advection = true;
  m_do_diffusion = true;
  m_do_source    = true;
  m_do_rte       = true;
  m_do_poisson   = true;
  
  {
    ParmParse pp("splitstep_tga");

    std::string str;
    
    pp.query("rk2_alpha", m_alpha);
    
    if(pp.contains("coupling")){

      pp.query("coupling", str);
      if(str == "semi_implicit"){
	//	m_simpi = true; Not yet supportd
      }
      else if(str == "explicit"){
	m_simpi = false;
      }
    }

    if(pp.contains("turn_off_advection")){
      pp.get("turn_off_advection", str);
      if(str == "true"){
	m_do_advection = false;
	if(m_verbosity > 2){
	  pout() << "splitstep_tga::splitstep_tga - Turning off advection" << endl;
	}
      }
    }
    
    if(pp.contains("turn_off_diffusion")){
      pp.get("turn_off_diffusion", str);
      if(str == "true"){
	m_do_diffusion = false;
	if(m_verbosity > 2){
	  pout() << "splitstep_tga::splitstep_tga - Turning off diffusion" << endl;
	}
      }
    }
    
    if(pp.contains("turn_off_source")){
      pp.get("turn_off_source", str);
      if(str == "true"){
	m_do_source = false;

	if(m_verbosity > 2){
	  pout() << "splitstep_tga::splitstep_tga - Turning off source" << endl;
	}
      }
    }

    if(pp.contains("turn_off_rte")){
      pp.get("turn_off_rte", str);
      if(str == "true"){
	m_do_rte = false;

	if(m_verbosity > 2){
	  pout() << "splitstep_tga::splitstep_tga - Turning off rte" << endl;
	}
      }
    }

    if(pp.contains("turn_off_poisson")){
      pp.get("turn_off_poisson", str);
      if(str == "true"){
	m_do_poisson = false;

	if(m_verbosity > 2){
	  pout() << "splitstep_tga::splitstep_tga - Turning off poisson" << endl;
	}
      }
    }
  }
}

splitstep_tga::~splitstep_tga(){

}

RefCountedPtr<cdr_storage>& splitstep_tga::get_cdr_storage(const cdr_iterator& a_solverit){
  return m_cdr_scratch[a_solverit.get_solver()];
}

RefCountedPtr<rte_storage>& splitstep_tga::get_rte_storage(const rte_iterator& a_solverit){
  return m_rte_scratch[a_solverit.get_solver()];
}

Real splitstep_tga::restrict_dt(){
  return 1.E99;
}

Real splitstep_tga::advance(const Real a_dt){
  CH_TIME("splitstep_tga::advance");
  if(m_verbosity > 2){
    pout() << "splitstep_tga::advance" << endl;
  }

  this->cache_solutions(); // Cache old solutions. Not used for anything (yet), but derived classes might. 

  if(m_do_advection){
    this->advance_advection(a_dt); // Advective advance. After this, solvers contain the advected states. The poisson 
  }                                // solver contains the potential after advection
  if(m_do_diffusion){
    this->advance_diffusion(a_dt); // Diffusion advance. After this, solvers contain the diffused advected states.
  }
  if(m_do_source){
    this->advance_sources(a_dt);   // Source term advance. 
  }

  // Put solver back in useable state so that we can reliably compute the next time step. 
  this->compute_cdr_velocities();
  this->compute_cdr_diffusion();
  this->compute_cdr_sources();
  
  return a_dt;
}

void splitstep_tga::regrid_internals(){
  CH_TIME("time_stepper::regrid_internals");
  if(m_verbosity > 5){
    pout() << "time_stepper::regrid_internals" << endl;
  }
  
  this->allocate_cdr_storage();
  this->allocate_poisson_storage();
  this->allocate_rte_storage();
  this->allocate_sigma_storage();
}

void splitstep_tga::compute_dt(Real& a_dt, time_code::which_code& a_timecode){
  CH_TIME("time_stepper::compute_dt");
  if(m_verbosity > 5){
    pout() << "time_stepper::compute_dt" << endl;
  }

  Real dt = 1.E99;

  const Real dt_cfl = m_cfl*m_cdr->compute_cfl_dt();
  if(dt_cfl < dt){
    dt = dt_cfl;
    a_timecode = time_code::cfl;
  }

  const Real dt_src = m_src_growth*m_cdr->compute_source_dt();
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

void splitstep_tga::allocate_cdr_storage(){
  const int ncomp       = 1;
  const int num_species = m_plaskin->get_num_species();
  m_cdr_scratch.resize(num_species);
  
  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx] = RefCountedPtr<cdr_storage> (new cdr_storage(m_amr, m_cdr->get_phase(), ncomp));
    m_cdr_scratch[idx]->allocate_storage();
  }
}

void splitstep_tga::allocate_poisson_storage(){
  const int ncomp = 1;
  m_poisson_scratch = RefCountedPtr<poisson_storage> (new poisson_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_poisson_scratch->allocate_storage();
}

void splitstep_tga::allocate_rte_storage(){
  const int ncomp       = 1;
  const int num_photons = m_plaskin->get_num_photons();
  m_rte_scratch.resize(num_photons);
  
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx] = RefCountedPtr<rte_storage> (new rte_storage(m_amr, m_rte->get_phase(), ncomp));
    m_rte_scratch[idx]->allocate_storage();
  }
}

void splitstep_tga::allocate_sigma_storage(){
  const int ncomp = 1;
  m_sigma_scratch = RefCountedPtr<sigma_storage> (new sigma_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_sigma_scratch->allocate_storage();
}

void splitstep_tga::deallocate_internals(){
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

void splitstep_tga::cache_solutions(){
  CH_TIME("splitstep_tga::cache_solutions");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::cache_solutions" << endl;
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

void splitstep_tga::advance_advection(const Real a_dt){
  CH_TIME("splitstep_tga::advance_advection");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::advance_advection" << endl;
  }

  // Set sources and diffusion coefficients to zero
  m_cdr->set_source(0.0);
  m_cdr->set_diffco(0.0);

  // Compute necessary things for k1 advance
  this->compute_E_at_start_of_time_step();
  this->compute_cdr_velo_at_start_of_time_step();
  this->compute_cdr_eb_states_at_start_of_time_step();
  this->compute_cdr_fluxes_at_start_of_time_step();
  this->compute_sigma_flux_at_start_of_time_step();

  // Do k1 advance
  this->advance_advection_cdr_k1(a_dt);
  this->advance_advection_sigma_k1(a_dt);
  if(m_do_poisson){
    this->solve_poisson_k1();
  }
  this->compute_E_after_k1();

  // Recompute things in order to do k2 advance
  this->compute_cdr_eb_states_after_k1();
  this->compute_cdr_velo_after_k1(a_dt);
  this->compute_cdr_fluxes_after_k1(a_dt);
  this->compute_sigma_flux_after_k1();

  // Do k2 advance
  this->advance_advection_cdr_k2(a_dt);
  this->advance_advection_sigma_k2(a_dt);
  if(m_do_poisson){
    this->solve_poisson_k2();
  }
}

void splitstep_tga::compute_E_at_start_of_time_step(){
  CH_TIME("splitstep_tga::compute_E_at_start_of_time_step");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::compute_E_at_start_of_time_step" << endl;
  }
  
  EBAMRCellData& E_cell = m_poisson_scratch->get_E_cell();
  EBAMRFluxData& E_face = m_poisson_scratch->get_E_face();
  EBAMRIVData&   E_eb   = m_poisson_scratch->get_E_eb();

  const MFAMRCellData& phi = m_poisson->get_state();
  
  this->compute_E(E_cell, m_cdr->get_phase(), phi);     // Compute cell-centered field
  this->compute_E(E_face, m_cdr->get_phase(), E_cell);  // Compute face-centered field
  this->compute_E(E_eb,   m_cdr->get_phase(), E_cell);  // EB-centered field
}

void splitstep_tga::compute_cdr_velo_at_start_of_time_step(){
  CH_TIME("splitstep_tga::compute_cdr_velo_at_start_of_time_step");
  if(m_verbosity > 5){
    pout() << "splitstep_::compute_cdr_velo_at_start_of_time_step" << endl;
  }

  Vector<EBAMRCellData*> states     = m_cdr->get_states();
  Vector<EBAMRCellData*> velocities = m_cdr->get_velocities();
  this->compute_cdr_velocities(velocities, states, m_poisson_scratch->get_E_cell(), m_time);
}

void splitstep_tga::compute_cdr_eb_states_at_start_of_time_step(){
  CH_TIME("splitstep_tga::compute_cdr_eb_states_at_start_of_time_step");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::compute_cdr_eb_states_at_start_of_time_step" << endl;
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

void splitstep_tga::compute_cdr_fluxes_at_start_of_time_step(){
  CH_TIME("splitstep_tga::compute_cdr_fluxes_at_start_of_time_step");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::compute_cdr_fluxes_at_start_of_time_step" << endl;
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

    extrap_cdr_densities.push_back(&dens_eb);  // Already been computed
    extrap_cdr_velocities.push_back(&velo_eb);
    extrap_cdr_fluxes.push_back(&flux_eb);
    extrap_cdr_gradients.push_back(&grad_eb);  // Already been computed
  }

  // Extrapolate densities, velocities, and fluxes
  Vector<EBAMRCellData*> cdr_densities = m_cdr->get_states();
  Vector<EBAMRCellData*> cdr_velocities = m_cdr->get_velocities();
  this->compute_extrapolated_fluxes(extrap_cdr_fluxes, cdr_densities, cdr_velocities, m_cdr->get_phase());
  this->extrapolate_to_eb(extrap_cdr_velocities, m_cdr->get_phase(), cdr_velocities);
  //  this->extrapolate_to_eb(extrap_cdr_densities,  m_cdr->get_phase(), cdr_densities); // Already been done, no?

  // Compute RTE flux on the boundary
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);

    EBAMRIVData& flux_eb = storage->get_eb_flux();
    solver->compute_boundary_flux(flux_eb, solver->get_state());
    extrap_rte_fluxes.push_back(&flux_eb);
  }

  const EBAMRIVData& E = m_poisson_scratch->get_E_eb();

  this->compute_cdr_fluxes(cdr_fluxes,
			   extrap_cdr_fluxes,
			   extrap_cdr_densities,
			   extrap_cdr_velocities,
			   extrap_cdr_gradients,
			   extrap_rte_fluxes,
			   E,
			   m_time);
}

void splitstep_tga::compute_sigma_flux_at_start_of_time_step(){
  CH_TIME("splitstep_tga::compute_sigma_flux_at_start_of_time_step");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::compute_sigma_flux_at_start_of_time_step" << endl;
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

void splitstep_tga::advance_advection_cdr_k1(const Real a_dt){
  CH_TIME("splitstep_tga::advance_advection_cdr_k1");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::advance_advection_cdr_k1" << endl;
  }

  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& k1  = storage->get_k1();
    EBAMRCellData& phi = storage->get_phi();

    const EBAMRCellData& state = solver->get_state();

    // Compute rhs
    data_ops::set_value(k1, 0.0);
    solver->compute_divF(k1, state, 0.0, true);
    data_ops::scale(k1, -1.0);

    // Compute phi
    data_ops::set_value(phi, 0.0);
    data_ops::incr(phi, state, 1.0);
    data_ops::incr(phi, k1,    m_alpha*a_dt);

    m_amr->average_down(phi, m_cdr->get_phase());
    m_amr->interp_ghost(phi, m_cdr->get_phase());

    data_ops::floor(phi, 0.0);
  }
}

void splitstep_tga::advance_advection_sigma_k1(const Real a_dt){
  CH_TIME("splitstep_tga::advance_advection_sigma_k1");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::advance_advection_sigma_k1" << endl;
  }

  const EBAMRIVData& state = m_sigma->get_state();
  
  EBAMRIVData& k1  = m_sigma_scratch->get_k1();
  EBAMRIVData& phi = m_sigma_scratch->get_phi();
  m_sigma->compute_rhs(k1);
  data_ops::set_value(phi,   0.0);
  data_ops::incr(phi, state, 1.0);
  data_ops::incr(phi, k1,    m_alpha*a_dt);

  m_amr->average_down(phi, m_cdr->get_phase());
  
  m_sigma->reset_cells(k1);
  m_sigma->reset_cells(phi);
}

void splitstep_tga::solve_poisson_k1(){
  CH_TIME("splitstep_tga::solve_poisson_k1");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::solve_poisson_k1" << endl;
  }

  MFAMRCellData& scratch_pot = m_poisson_scratch->get_phi();
  EBAMRIVData& sigma         = m_sigma_scratch->get_phi();
  Vector<EBAMRCellData*> cdr_densities;
  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    cdr_densities.push_back(&(storage->get_phi()));
  }

  data_ops::set_value(scratch_pot, 0.0);
  data_ops::incr(scratch_pot, m_poisson->get_state(), 1.0);

  if((m_step + 1) % m_fast_poisson == 0){
    bool converged = this->solve_poisson(scratch_pot, m_poisson->get_source(), cdr_densities, sigma, centering::cell_center);
    if(!converged){
      pout() << "splitstep_tga::solve_poisson_k1 - solver did not converge at step = " << m_step << endl;
    }
  }
}

void splitstep_tga::compute_E_after_k1(){
  CH_TIME("splitstep_tga::compute_E_after_k1");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::compute_E_after_k1" << endl;
  }
  
  EBAMRCellData& E_cell = m_poisson_scratch->get_E_cell();
  EBAMRFluxData& E_face = m_poisson_scratch->get_E_face();
  EBAMRIVData&   E_eb   = m_poisson_scratch->get_E_eb();

  const MFAMRCellData& phi = m_poisson_scratch->get_phi();
  
  this->compute_E(E_cell, m_cdr->get_phase(), phi);     // Compute cell-centered field
  this->compute_E(E_face, m_cdr->get_phase(), E_cell);  // Compute face-centered field
  this->compute_E(E_eb,   m_cdr->get_phase(), E_cell);  // EB-centered field
}

void splitstep_tga::compute_cdr_eb_states_after_k1(){
  CH_TIME("splitstep_tga::compute_cdr_eb_states_after_k1");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::compute_cdr_eb_states_after_k1" << endl;
  }

  Vector<EBAMRIVData*>   eb_gradients;
  Vector<EBAMRIVData*>   eb_states;
  Vector<EBAMRCellData*> cdr_states;
  
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    cdr_states.push_back(&(storage->get_phi()));
    eb_states.push_back(&(storage->get_eb_state()));
    eb_gradients.push_back(&(storage->get_eb_grad()));
  }

  this->extrapolate_to_eb(eb_states,          m_cdr->get_phase(), cdr_states);
  this->compute_gradients_at_eb(eb_gradients, m_cdr->get_phase(), cdr_states);
}

void splitstep_tga::compute_cdr_velo_after_k1(const Real a_dt){
  CH_TIME("splitstep_tga::compute_cdr_velo_after_k1");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::compute_cdr_velo_after_k1";
  }
  
  const int num_species = m_plaskin->get_num_species();
  
  const Real time = m_time + m_alpha*a_dt;
  
  Vector<EBAMRCellData*> states(num_species);
  Vector<EBAMRCellData*> velocities = m_cdr->get_velocities();

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    const int idx = solver_it.get_solver();
    states[idx] = &(storage->get_phi());
  }
  
  this->compute_cdr_velocities(velocities, states, m_poisson_scratch->get_E_cell(), time);
}

void splitstep_tga::compute_cdr_fluxes_after_k1(const Real a_dt){
  CH_TIME("splitstep_tga::compute_cdr_fluxes_after_k1");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::compute_cdr_fluxes_after_k1" << endl;
  }

  const Real time = m_time + m_alpha*a_dt;
  
  Vector<EBAMRIVData*> cdr_fluxes;
  Vector<EBAMRIVData*> extrap_cdr_fluxes;
  Vector<EBAMRIVData*> extrap_cdr_densities;
  Vector<EBAMRIVData*> extrap_cdr_velocities;
  Vector<EBAMRIVData*> extrap_cdr_gradients;
  Vector<EBAMRIVData*> extrap_rte_fluxes;
  
  Vector<EBAMRCellData*> cdr_densities;
  Vector<EBAMRCellData*> rte_densities;

  cdr_fluxes = m_cdr->get_ebflux();

  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& dens  = storage->get_phi();
    EBAMRIVData& dens_eb = storage->get_eb_state();
    EBAMRIVData& velo_eb = storage->get_eb_velo();
    EBAMRIVData& flux_eb = storage->get_eb_flux();
    EBAMRIVData& grad_eb = storage->get_eb_grad();

    cdr_densities.push_back(&dens);
    extrap_cdr_densities.push_back(&dens_eb);  // This has already been extrapolated to the EB
    extrap_cdr_velocities.push_back(&velo_eb); // This has not.
    extrap_cdr_fluxes.push_back(&flux_eb);     // This hasn't either. 
    extrap_cdr_gradients.push_back(&grad_eb);  // This has already been extrapolated to the EB
  }

  // Extrapolate the flux and the velocity
  Vector<EBAMRCellData*> cdr_velocities = m_cdr->get_velocities();
  this->compute_extrapolated_fluxes(extrap_cdr_fluxes, cdr_densities, cdr_velocities, m_cdr->get_phase());
  this->extrapolate_to_eb(extrap_cdr_velocities, m_cdr->get_phase(), cdr_velocities);
  //  this->extrapolate_to_eb(extrap_cdr_densities,  m_cdr->get_phase(), cdr_densities); // This has already been done, no?

  // Compute RTE flux on the boundary
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);

    EBAMRIVData& flux_eb = storage->get_eb_flux();
    solver->compute_boundary_flux(flux_eb, storage->get_phi());
    extrap_rte_fluxes.push_back(&flux_eb);
  }

  const EBAMRIVData& E = m_poisson_scratch->get_E_eb();

  this->compute_cdr_fluxes(cdr_fluxes,
			   extrap_cdr_fluxes,
			   extrap_cdr_densities,
			   extrap_cdr_velocities,
			   extrap_cdr_gradients, 
			   extrap_rte_fluxes,
			   E,
			   time);
}

void splitstep_tga::compute_sigma_flux_after_k1(){
  CH_TIME("splitstep_tga::compute_sigma_flux_after_k1");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::compute_sigma_flux_after_k1" << endl;
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

void splitstep_tga::advance_advection_cdr_k2(const Real a_dt){
  CH_TIME("splitstep_tga::advance_advection_cdr_k2");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::advance_advection_cdr_k2" << endl;
  }

  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& state   = solver->get_state();
    EBAMRCellData& k1      = storage->get_k1();
    EBAMRCellData& k2      = storage->get_k2();
    EBAMRCellData& phi     = storage->get_phi();

    // Compute RHS
    data_ops::set_value(k2, 0.0);
    solver->compute_divF(k2, phi, 0.0, true);
    data_ops::scale(k2, -1.0);

    // RK2 advance
    data_ops::incr(state, k1, a_dt*(1 - 1./(2.*m_alpha)));
    data_ops::incr(state, k2, a_dt*1./(2.*m_alpha));

    m_amr->average_down(state, m_cdr->get_phase());
    m_amr->interp_ghost(state, m_cdr->get_phase());

    data_ops::floor(state, 0.0);
  }
}

void splitstep_tga::advance_advection_sigma_k2(const Real a_dt){
  CH_TIME("splitstep_tga::advance_advection_sigma_k2");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::advance_advection_sigma_k2" << endl;
  }
  
  EBAMRIVData& k1  = m_sigma_scratch->get_k1();
  EBAMRIVData& k2  = m_sigma_scratch->get_k2();
  m_sigma->compute_rhs(k2);

  EBAMRIVData& state = m_sigma->get_state();
  data_ops::incr(state, k1, a_dt*(1 - 1./(2.*m_alpha)));
  data_ops::incr(state, k2, a_dt*1./(2.*m_alpha));

  m_sigma->reset_cells(state);
}

void splitstep_tga::solve_poisson_k2(){
  CH_TIME("splitstep_tga::solve_poisson_k2");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::solve_poisson_k2" << endl;
  }

  // We computed the intermediate potential at time t_k + alpha*dt. Linearly extrapolate that result to the end of
  // the time step. The result for this is y_extrap = y_alpha/alpha - y_0*(1-alpha)/alpha. This (usually) brings the
  // initial guess closer to the true solution.
  MFAMRCellData& pot     = m_poisson->get_state();
  MFAMRCellData& phi     = m_poisson_scratch->get_phi();

  if((m_step + 1) % m_fast_poisson == 0){
    data_ops::scale(pot, -(1.0 - m_alpha)/m_alpha);
    data_ops::incr(pot, phi, 1./m_alpha);

    const bool converged = this->solve_poisson(pot,
					       m_poisson->get_source(),
					       m_cdr->get_states(),
					       m_sigma->get_state(),
					       centering::cell_center);
    if(!converged){
      pout() << "rk2::solve_poisson_k2 - solver did not converge at step = " << m_step << endl;
    }
  }
  else{
    data_ops::set_value(pot, 0.0);
    data_ops::incr(pot, phi, 1.0);
  }
}

void splitstep_tga::advance_diffusion(const Real a_dt){
  CH_TIME("splitstep_tga::advance_diffusion");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::advance_diffusion" << endl;
  }

  // Set source terms and velocity to zero. Then compute diffusion coefficients

  //  m_cdr->set_velocity(RealVect::Zero); // This breaks sometimes, but doing this is not necessary for diffusion advance
  m_cdr->set_source(0.0); // This is necessary because advance_diffusion also takes a source term. 

  // Diffusive advance for all cdr equations
  this->compute_cdr_diffusion();
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();

    solver->advance_diffusion(a_dt);
  }

  // Update poisson equation afterwards
  this->solve_poisson();
}

void splitstep_tga::advance_sources(const Real a_dt){
  CH_TIME("splitstep_tga::advance_sources");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::advance_sources" << endl;
  }

  // Set sources and diffusion coefficients to zero
  //  m_cdr->set_velocity(RealVect::Zero);
  m_cdr->set_diffco(0.0);

  // Compute necessary things for k1 advance
  this->compute_E_at_start_of_time_step();
  this->compute_cdr_sources_at_start_of_time_step();

  // Do k1 advance
  this->advance_source_cdr_k1(a_dt);
  if(m_do_poisson){
    this->solve_poisson_k1();
  }
  if(m_do_rte){
    if(m_rte->is_stationary()){
      this->advance_rte_k1_stationary(a_dt);
    }
    else{
      this->advance_rte_k1_transient(a_dt);
    }
  }


  // Recompute necessary things for k2 advance
  this->compute_E_after_k1();
  this->compute_cdr_sources_after_k1(a_dt);

  // Do k2 advance
  this->advance_source_cdr_k2(a_dt);
  if(m_do_poisson){
    this->solve_poisson_k2();
  }
  if(m_do_rte){
    if(m_rte->is_stationary()){
      this->advance_rte_k2_stationary(a_dt);
    }
    else{
      this->advance_rte_k2_transient(a_dt);
    }
  }
}

void splitstep_tga::compute_cdr_sources_at_start_of_time_step(){
  CH_TIME("splitstep_tga::compute_cdr_sources_at_start_of_time_step");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::compute_cdr_sources_at_start_of_time_step" << endl;
  }
  
  Vector<EBAMRCellData*> cdr_sources = m_cdr->get_sources();
  Vector<EBAMRCellData*> cdr_states  = m_cdr->get_states();
  Vector<EBAMRCellData*> rte_states  = m_rte->get_states();
  EBAMRCellData& E                   = m_poisson_scratch->get_E_cell();

  this->compute_cdr_sources(cdr_sources, cdr_states, rte_states, E, m_time, centering::cell_center);
}

void splitstep_tga::advance_source_cdr_k1(const Real a_dt){
  CH_TIME("splitstep_tga::advance_source_cdr_k1");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::advance_source_cdr_k1" << endl;
  }

  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& k1  = storage->get_k1();
    EBAMRCellData& phi = storage->get_phi();

    const EBAMRCellData& state = solver->get_state();

    data_ops::copy(k1, solver->get_source());

    data_ops::set_value(phi, 0.0);
    data_ops::incr(phi, state, 1.0);
    data_ops::incr(phi, k1,    m_alpha*a_dt);

    m_amr->average_down(phi, m_cdr->get_phase());
    m_amr->interp_ghost(phi, m_cdr->get_phase());

    data_ops::floor(phi, 0.0);
  }
}

void splitstep_tga::advance_rte_k1_stationary(const Real a_dt){
  CH_TIME("splitstep_tga::compute_rte_k1_stationary");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::compute_k1_stationary" << endl;
  }

  const Real time = m_time + m_alpha*a_dt;

  Vector<EBAMRCellData*> rte_states;
  Vector<EBAMRCellData*> rte_sources;
  Vector<EBAMRCellData*> cdr_states;

  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);
    
    EBAMRCellData& phi    = storage->get_phi();
    EBAMRCellData& state  = solver->get_state();
    EBAMRCellData& source = solver->get_source();

    data_ops::set_value(phi, 0.0);
    data_ops::incr(phi, state, 1.0);

    rte_states.push_back(&(phi));
    rte_sources.push_back(&(source));
  }

  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    cdr_states.push_back(&(storage->get_phi()));
  }

  EBAMRCellData& E = m_poisson_scratch->get_E_cell();

  if((m_step + 1) % m_fast_rte == 0){
    const Real dummy_dt = 0.0;
    this->solve_rte(rte_states, rte_sources, cdr_states, E, time, dummy_dt, centering::cell_center);
  }
}

void splitstep_tga::advance_rte_k1_transient(const Real a_dt){
  CH_TIME("splitstep_tga::compute_rte_k1_transient");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::compute_k1_transient" << endl;
  }

  MayDay::Abort("splitstep_tga::advance_rte_k1_transient - transient RTE not (yet) supported for this timestepper");
}

void splitstep_tga::compute_cdr_sources_after_k1(const Real a_dt){
  CH_TIME("splitstep_tga::compute_cdr_sources_after_k1");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::compute_cdr_sources_after_k1" << endl;
  }

  const Real time = m_time + m_alpha*a_dt;

  Vector<EBAMRCellData*> cdr_sources;
  Vector<EBAMRCellData*> cdr_states;
  Vector<EBAMRCellData*> rte_states;

  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);
    rte_states.push_back(&(storage->get_phi()));
  }

  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    
    cdr_states.push_back(&(storage->get_phi()));
    cdr_sources.push_back(&(solver->get_source()));
  }

  EBAMRCellData& E = m_poisson_scratch->get_E_cell();

  this->compute_cdr_sources(cdr_sources, cdr_states, rte_states, E, time, centering::cell_center);
}

void splitstep_tga::advance_source_cdr_k2(const Real a_dt){
  CH_TIME("splitstep_tga::advance_source_cdr_k2");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::advance_source_cdr_k2" << endl;
  }

  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& state   = solver->get_state();
    EBAMRCellData& k1      = storage->get_k1();
    EBAMRCellData& k2      = storage->get_k2();
    EBAMRCellData& phi     = storage->get_phi();

    // Source term
    data_ops::copy(k2, solver->get_source());

    // 
    data_ops::incr(state, k1, a_dt*(1 - 1./(2.*m_alpha)));
    data_ops::incr(state, k2, a_dt*1./(2.*m_alpha));

    m_amr->average_down(state, m_cdr->get_phase());
    m_amr->interp_ghost(state, m_cdr->get_phase());

    data_ops::floor(state, 0.0);
  }
}

void splitstep_tga::advance_rte_k2_stationary(const Real a_dt){
  CH_TIME("splitstep_tga::compute_rte_k2_stationary");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::compute_rte_k2_stationary" << endl;
  }

  const Real time = m_time + a_dt; // Source terms centered on the end of the time step

  Vector<EBAMRCellData*> rte_states;
  Vector<EBAMRCellData*> rte_sources;
  Vector<EBAMRCellData*> cdr_states;

  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);
    
    EBAMRCellData& phi    = storage->get_phi();
    EBAMRCellData& state  = solver->get_state();
    EBAMRCellData& source = solver->get_source();

    
    // We computed the intermediate potential at time t_k + alpha*dt. Linearly extrapolate that result to the end of
    // the time step. The result for this is y_extrap = y_alpha/alpha - y_0*(1-alpha)/alpha. This (usually) brings the
    // initial guess closer to the true solution.
    if((m_step + 1) % m_fast_rte == 0){
      data_ops::scale(state, -(1.0 - m_alpha)/(m_alpha));
      data_ops::incr(state, phi, 1.0/m_alpha);
    }
    else {
      data_ops::set_value(state, 0.0);
      data_ops::incr(state, phi, 1.0);
    }
    
    rte_states.push_back(&(state));
    rte_sources.push_back(&(source));
  }

  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    cdr_states.push_back(&(solver->get_state()));
  }

  EBAMRCellData& E = m_poisson_scratch->get_E_cell();

  if((m_step + 1) % m_fast_rte == 0){
    const Real dummy_dt = 0.0;
    this->solve_rte(rte_states, rte_sources, cdr_states, E, time, dummy_dt, centering::cell_center);
  }
}

void splitstep_tga::advance_rte_k2_transient(const Real a_dt){
  CH_TIME("splitstep_tga::compute_rte_k1_transient");
  if(m_verbosity > 5){
    pout() << "splitstep_tga::compute_k1_transient" << endl;
  }

  MayDay::Abort("splitstep_tga::advance_rte_k2_transient - transient RTE not (yet) supported for this timestepper");

}
