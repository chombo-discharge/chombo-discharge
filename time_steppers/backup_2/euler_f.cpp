/*!
  @file   euler_f.cpp
  @brief  Implementation of euler_f.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "euler_f.H"
#include "euler_f_storage.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"
#include "data_ops.H"
#include "units.H"

typedef euler_f::cdr_storage     cdr_storage;
typedef euler_f::poisson_storage poisson_storage;
typedef euler_f::rte_storage     rte_storage;
typedef euler_f::sigma_storage   sigma_storage;

euler_f::euler_f() : time_stepper(){

}

euler_f::~euler_f(){

}

RefCountedPtr<cdr_storage>& euler_f::get_cdr_storage(const cdr_iterator& a_solverit){
  return m_cdr_scratch[a_solverit.get_solver()];
}

RefCountedPtr<rte_storage>& euler_f::get_rte_storage(const rte_iterator& a_solverit){
  return m_rte_scratch[a_solverit.get_solver()];
}

void euler_f::allocate_cdr_storage(){

  const int ncomp       = 1;
  const int num_species = m_plaskin->get_num_species();
  m_cdr_scratch.resize(num_species);
  
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx] = RefCountedPtr<cdr_storage> (new cdr_storage(m_amr, m_cdr->get_phase(), ncomp));
    m_cdr_scratch[idx]->allocate_storage();
  }
}

void euler_f::allocate_poisson_storage(){
  const int ncomp = 1;
  m_poisson_scratch = RefCountedPtr<poisson_storage> (new poisson_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_poisson_scratch->allocate_storage();
}

void euler_f::allocate_rte_storage(){
  const int ncomp       = 1;
  const int num_photons = m_plaskin->get_num_photons();
  m_rte_scratch.resize(num_photons);
  
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx] = RefCountedPtr<rte_storage> (new rte_storage(m_amr, m_rte->get_phase(), ncomp));
    m_rte_scratch[idx]->allocate_storage();
  }
}

void euler_f::allocate_sigma_storage(){
  const int ncomp = 1;
  m_sigma_scratch = RefCountedPtr<sigma_storage> (new sigma_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_sigma_scratch->allocate_storage();
}

void euler_f::regrid_internals(){
  CH_TIME("euler_f::regrid_internals");
  if(m_verbosity > 5){
    pout() << "euler_f::regrid_internals" << endl;
  }
  
  this->allocate_cdr_storage();
  this->allocate_poisson_storage();
  this->allocate_rte_storage();
  this->allocate_sigma_storage();
}

Real euler_f::advance(const Real a_dt){
  CH_TIME("euler_f::advance");
  if(m_verbosity > 2){
    pout() << "euler_f::advance" << endl;
  }

  this->compute_E_at_start_of_time_step();
  this->compute_cdr_eb_states_at_start_of_time_step();
  this->compute_cdr_velo_at_start_of_time_step();
  this->compute_cdr_diffco_at_start_of_time_step();
  this->compute_cdr_sources_at_start_of_time_step();
  this->compute_cdr_fluxes_at_start_of_time_step();
  this->compute_sigma_flux_at_start_of_time_step();

  // Advance
  this->advance_cdr(a_dt);
  this->advance_sigma(a_dt);
  this->advance_poisson();
  this->compute_E_after_poisson();
  if(m_rte->is_stationary()){
    this->advance_rte_stationary(a_dt);
  }
  else{
    this->advance_rte_transient(a_dt);
  }

  return a_dt;
}

void euler_f::compute_E_at_start_of_time_step(){
  CH_TIME("euler_f::compute_E_at_start_of_time_step");
  if(m_verbosity > 2){
    pout() << "euler_f::compute_E_at_start_of_time_step" << endl;
  }

  EBAMRCellData& E_cell = m_poisson_scratch->get_E_cell();
  EBAMRFluxData& E_face = m_poisson_scratch->get_E_face();
  EBAMRIVData&   E_eb   = m_poisson_scratch->get_E_eb();

  const MFAMRCellData& phi = m_poisson->get_state();
  
  time_stepper::compute_E(E_cell, m_cdr->get_phase(), phi);     // Compute cell-centered field
  time_stepper::compute_E(E_face, m_cdr->get_phase(), E_cell);  // Compute face-centered field
  time_stepper::compute_E(E_eb,   m_cdr->get_phase(), E_cell);  // EB-centered field
}

void euler_f::compute_cdr_eb_states_at_start_of_time_step(){
  CH_TIME("euler_f::compute_cdr_eb_states_at_start_of_time_step");
  if(m_verbosity > 2){
    pout() << "euler_f::compute_cdr_eb_states_at_start_of_time_step" << endl;
  }

  const irreg_amr_stencil<eb_centroid_interp>& stencil = m_amr->get_eb_centroid_interp_stencils(m_cdr->get_phase());
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    const EBAMRCellData& state              = solver->get_state();

    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    EBAMRIVData& cdr_eb = storage->get_eb_state();

    stencil.apply(cdr_eb, state);
  }
}

void euler_f::compute_cdr_velo_at_start_of_time_step(){
  CH_TIME("euler_f::compute_cdr_velo_at_start_of_time_step");
  if(m_verbosity > 2){
    pout() << "euler_f::compute_cdr_velo_at_start_of_time_step" << endl;
  }

  Vector<EBAMRCellData*> states     = m_cdr->get_states();
  Vector<EBAMRCellData*> velocities = m_cdr->get_velocities();
  this->compute_cdr_velocities(velocities, states, m_poisson_scratch->get_E_cell(), m_time);
}

void euler_f::compute_cdr_diffco_at_start_of_time_step(){
  CH_TIME("euler_f::compute_cdr_diffco_at_start_of_time_step");
  if(m_verbosity > 2){
    pout() << "euler_f::compute_cdr_diffco_at_start_of_time_step" << endl;
  }

  const int num_species = m_plaskin->get_num_species();

  Vector<EBAMRCellData*> cdr_states  = m_cdr->get_states();
  Vector<EBAMRFluxData*> diffco_face = m_cdr->get_diffco_face();
  Vector<EBAMRIVData*> diffco_eb     = m_cdr->get_diffco_eb();

  const EBAMRCellData& E_cell = m_poisson_scratch->get_E_cell();
  const EBAMRIVData& E_eb     = m_poisson_scratch->get_E_eb();

  // Get extrapolated states
  Vector<EBAMRIVData*> eb_states(num_species);
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    eb_states[idx] = &(storage->get_eb_state());
  }
  
  this->compute_cdr_diffco_face(diffco_face, cdr_states, E_cell, m_time);
  this->compute_cdr_diffco_eb(diffco_eb,     eb_states,  E_eb,   m_time);
}

void euler_f::compute_cdr_sources_at_start_of_time_step(){
  CH_TIME("euler_f::compute_cdr_sources_at_start_of_time_step");
  if(m_verbosity > 2){
    pout() << "euler_f::compute_cdr_sources_at_start_of_time_step" << endl;
  }

  Vector<EBAMRCellData*> cdr_sources = m_cdr->get_sources();
  Vector<EBAMRCellData*> cdr_states  = m_cdr->get_states();
  Vector<EBAMRCellData*> rte_states  = m_rte->get_states();
  EBAMRCellData& E                   = m_poisson_scratch->get_E_cell();

  this->compute_cdr_sources(cdr_sources, cdr_states, rte_states, E, m_time, centering::cell_center);
}

void euler_f::compute_cdr_fluxes_at_start_of_time_step(){
  CH_TIME("euler_f::compute_cdr_fluxes_at_start_of_time_step");
  if(m_verbosity > 2){
    pout() << "euler_f::compute_cdr_fluxes_at_start_of_time_step" << endl;
  }

  Vector<EBAMRIVData*> cdr_fluxes;
  Vector<EBAMRIVData*> extrap_cdr_fluxes;
  Vector<EBAMRIVData*> extrap_cdr_densities;
  Vector<EBAMRIVData*> extrap_cdr_velocities;
  Vector<EBAMRIVData*> extrap_rte_fluxes;

  cdr_fluxes = m_cdr->get_ebflux();

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRIVData& dens_eb = storage->get_eb_state();
    EBAMRIVData& velo_eb = storage->get_eb_velo();
    EBAMRIVData& flux_eb = storage->get_eb_flux();

    extrap_cdr_densities.push_back(&dens_eb);
    extrap_cdr_velocities.push_back(&velo_eb);
    extrap_cdr_fluxes.push_back(&flux_eb);
  }

  // Extrapolate densities, velocities, and fluxes
  Vector<EBAMRCellData*> cdr_densities = m_cdr->get_states();
  Vector<EBAMRCellData*> cdr_velocities = m_cdr->get_velocities();
  this->compute_extrapolated_fluxes(extrap_cdr_fluxes, cdr_densities, cdr_velocities, m_cdr->get_phase());
  this->extrapolate_to_eb(extrap_cdr_densities,  m_cdr->get_phase(), cdr_densities); // Already been done, no?
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

  this->compute_cdr_fluxes(cdr_fluxes,
			   extrap_cdr_fluxes,
			   extrap_cdr_densities,
			   extrap_cdr_velocities,
			   extrap_rte_fluxes,
			   E,
			   m_time);
}

void euler_f::compute_sigma_flux_at_start_of_time_step(){
  CH_TIME("euler_f::compute_sigma_flux_at_start_of_time_step");
  if(m_verbosity > 2){
    pout() << "euler_f::compute_sigma_flux_at_start_of_time_step" << endl;
  }

  EBAMRIVData& flux = m_sigma->get_flux();
  data_ops::set_value(flux, 0.0);

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    const RefCountedPtr<species>& spec      = solver_it.get_species();
    const EBAMRIVData& solver_flux          = solver->get_ebflux();

    data_ops::incr(flux, solver_flux, spec->get_charge()*units::s_Qe);
  }

  m_sigma->reset_cells(flux);
}

void euler_f::advance_cdr(const Real a_dt){
  CH_TIME("euler_f::advance_cdr");
  if(m_verbosity > 2){
    pout() << "euler_f::advance_cdr" << endl;
  }

  const int ncomp = 1;

  EBAMRCellData rhs;
  m_amr->allocate(rhs, m_cdr->get_phase(), ncomp);

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();

    EBAMRCellData& state = solver->get_state();
    solver->compute_rhs(rhs, state, a_dt);

    data_ops::incr(state, rhs, a_dt);

    m_amr->average_down(state, m_cdr->get_phase());
    m_amr->interp_ghost(state, m_cdr->get_phase());

    data_ops::floor(state, 0.0);
  }
}

void euler_f::advance_sigma(const Real a_dt){
  CH_TIME("euler_f::advance_sigma");
  if(m_verbosity > 2){
    pout() << "euler_f::advance_sigma" << endl;
  }

  EBAMRIVData& state = m_sigma->get_state();
  EBAMRIVData& rhs   = m_sigma_scratch->get_rhs();

  m_sigma->compute_rhs(rhs);
  data_ops::incr(state, rhs, a_dt);

  m_amr->average_down(state, m_cdr->get_phase());
  
  m_sigma->reset_cells(state);
}

void euler_f::advance_poisson(){
  CH_TIME("euler_f::advance_poisson");
  if(m_verbosity > 2){
    pout() << "euler_f::advance_poisson" << endl;
  }

  if((m_step + 1) % m_fast_poisson == 0){
    const bool converged = this->solve_poisson();
    if(!converged){
      pout() << "euler_f::advance_poisson - solver did not converge at step = " << m_step << endl;
    }
  }
}

void euler_f::compute_E_after_poisson(){
  CH_TIME("euler_f::compute_E_after_poisson");
  if(m_verbosity > 2){
    pout() << "euler_f::compute_E_after_poisson" << endl;
  }

  this->compute_E_at_start_of_time_step();
}

void euler_f::advance_rte_stationary(const Real a_dt){
  CH_TIME("euler_f::advance_rte_stationary");
  if(m_verbosity > 2){
    pout() << "euler_f::advance_rte_stationary" << endl;
  }

  Vector<EBAMRCellData*> rte_states = m_rte->get_states();
  Vector<EBAMRCellData*> rte_source = m_rte->get_sources();
  Vector<EBAMRCellData*> cdr_states = m_cdr->get_states();

  const EBAMRCellData& E = m_poisson_scratch->get_E_cell();
  const Real time        = m_time + m_alpha*a_dt;

  if((m_step + 1) % m_fast_rte == 0){
    const Real dummy_dt = 0.0;
    this->solve_rte(rte_states, rte_source, cdr_states, E, time, dummy_dt, centering::cell_center);
  }
}

void euler_f::advance_rte_transient(const Real a_dt){
  CH_TIME("euler_f::advance_rte_transient");
  if(m_verbosity > 2){
    pout() << "euler_f::advance_rte_transient" << endl;
  }

  MayDay::Abort("euler_f::advance_rte_transient - not implemented (yet). Please use stationary RTE solvers");
}

Real euler_f::restrict_dt(){
  return 1.E99;
}
