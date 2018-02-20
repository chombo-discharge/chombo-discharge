/*!
  @file   splitstep_euler_f.cpp
  @brief  Implementation of spltistep_euler_f.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "splitstep_euler_f.H"
#include "splitstep_euler_f_storage.H"
#include "data_ops.H"
#include "units.H"

typedef splitstep_euler_f::cdr_storage     cdr_storage;
typedef splitstep_euler_f::poisson_storage poisson_storage;
typedef splitstep_euler_f::rte_storage     rte_storage;
typedef splitstep_euler_f::sigma_storage   sigma_storage;

splitstep_euler_f::splitstep_euler_f(){

}

splitstep_euler_f::~splitstep_euler_f(){

}

RefCountedPtr<cdr_storage>& splitstep_euler_f::get_cdr_storage(const cdr_iterator& a_solverit){
  return m_cdr_scratch[a_solverit.get_solver()];
}

RefCountedPtr<rte_storage>& splitstep_euler_f::get_rte_storage(const rte_iterator& a_solverit){
  return m_rte_scratch[a_solverit.get_solver()];
}

Real splitstep_euler_f::restrict_dt(){
  return 1.E99;
}

Real splitstep_euler_f::advance(const Real a_dt){
  CH_TIME("splitstep_euler_f::advance");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::advance" << endl;
  }

  // We don't support this (yet)
  if(!m_rte->is_stationary()){
    MayDay::Abort("splitstep_euler_f::advance - transient RTE currently not supported. Please use a different time stepper");
  }

  // Transport step
  this->compute_E_at_start_of_time_step();               // Compute E from previous solution
  this->compute_cdr_velo_at_start_of_time_step();        // Compute cdr velocities
  this->compute_cdr_diffco_at_start_of_time_step();      // Compute cdr diffusion coefficients
  this->compute_cdr_fluxes_at_start_of_time_step();      // Compute cdr fluxes
  this->compute_sigma_flux_at_start_of_time_step();      // Compute sigma fluxes
  this->set_cdr_source_to_zero_at_start_of_time_step();  // Set source term to zero
  
  this->advance_cdr_transport(a_dt);             // Forward Euler method for transport step
  this->advance_sigma_transport(a_dt);           // Forward Euler method for the transport step
  this->solve_poisson_transport();               // Poisson solve after transport step
  this->compute_E_after_transport();             // Compute a new E for the source term advance
  if(m_rte->is_stationary()){
    this->advance_rte_transport_stationary();    // Stationary RTE solve
  }
  else{
    this->advance_rte_transport_transient(a_dt); // Advance transiently using a zero source
  }

  // Source step
  this->compute_cdr_source_after_transport();     // Compute source term
  this->advance_cdr_source(a_dt);                 // Advance with source terms
  this->advance_sigma_source(a_dt);               // Advance with source terms
  this->solve_poisson_source();                   // Recompute Poisson equation. 
  this->compute_E_after_source();                 // Recompute electric fields
  if(m_rte->is_stationary()){
    this->advance_rte_transport_stationary();     // Stationary RTE solve
  }
  else{
    this->advance_rte_transport_transient(a_dt);  // Advance transiently using a zero source
  }

  return a_dt;
}

void splitstep_euler_f::regrid_internals(){
  CH_TIME("splitstep_euler_f::regrid_internals");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::regrid_internals" << endl;
  }
  
  this->allocate_cdr_storage();
  this->allocate_poisson_storage();
  this->allocate_rte_storage();
  this->allocate_sigma_storage();
}

void splitstep_euler_f::allocate_cdr_storage(){
  const int ncomp       = 1;
  const int num_species = m_plaskin->get_num_species();
  m_cdr_scratch.resize(num_species);
  
  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx] = RefCountedPtr<cdr_storage> (new cdr_storage(m_amr, m_cdr->get_phase(), ncomp));
    m_cdr_scratch[idx]->allocate_storage();
  }
}

void splitstep_euler_f::allocate_poisson_storage(){
  const int ncomp = 1;
  m_poisson_scratch = RefCountedPtr<poisson_storage> (new poisson_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_poisson_scratch->allocate_storage();
}

void splitstep_euler_f::allocate_rte_storage(){
  const int ncomp       = 1;
  const int num_photons = m_plaskin->get_num_photons();
  m_rte_scratch.resize(num_photons);
  
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx] = RefCountedPtr<rte_storage> (new rte_storage(m_amr, m_rte->get_phase(), ncomp));
    m_rte_scratch[idx]->allocate_storage();
  }
}

void splitstep_euler_f::allocate_sigma_storage(){
  const int ncomp = 1;
  m_sigma_scratch = RefCountedPtr<sigma_storage> (new sigma_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_sigma_scratch->allocate_storage();
}

void splitstep_euler_f::compute_E_at_start_of_time_step(){
  CH_TIME("splitstep_euler_f::compute_E_at_start_of_time_step");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::compute_E_at_start_of_time_step" << endl;
  }

  EBAMRCellData& E_cell = m_poisson_scratch->get_E_cell();
  EBAMRFluxData& E_face = m_poisson_scratch->get_E_face();
  EBAMRIVData&   E_eb   = m_poisson_scratch->get_E_eb();

  const MFAMRCellData& phi = m_poisson->get_state();
  
  this->compute_E(E_cell, m_cdr->get_phase(), phi);     // Compute cell-centered field
  this->compute_E(E_face, m_cdr->get_phase(), E_cell);  // Compute face-centered field
  this->compute_E(E_eb,   m_cdr->get_phase(), E_cell);  // EB-centered field
}

void splitstep_euler_f::compute_cdr_velo_at_start_of_time_step(){
  CH_TIME("splitstep_euler_f::compute_cdr_velo_at_start_of_time_step");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::compute_cdr_velo_at_start_of_time_step" << endl;
  }

  Vector<EBAMRCellData*> velocities = m_cdr->get_velocities();
  this->compute_cdr_velocities(velocities, m_poisson_scratch->get_E_cell());
}

void splitstep_euler_f::compute_cdr_diffco_at_start_of_time_step(){
  CH_TIME("splitstep_euler_f::compute_cdr_diffco_at_start_of_time_step");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::compute_cdr_diffco_at_start_of_time_step" << endl;
  }

  Vector<EBAMRFluxData*> diffco_face = m_cdr->get_diffco_face();
  Vector<EBAMRIVData*> diffco_eb     = m_cdr->get_diffco_eb();

  const EBAMRFluxData& E_face = m_poisson_scratch->get_E_face();
  const EBAMRIVData& E_eb     = m_poisson_scratch->get_E_eb();
  
  this->compute_cdr_diffco_face(diffco_face, E_face);
  this->compute_cdr_diffco_eb(diffco_eb,     E_eb);
}

void splitstep_euler_f::compute_cdr_fluxes_at_start_of_time_step(){
  CH_TIME("splitstep_euler_f::compute_cdr_fluxes_at_start_of_time_step");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::compute_cdr_fluxes_at_start_of_time_step" << endl;
  }

  Vector<EBAMRIVData*> cdr_fluxes;
  Vector<EBAMRIVData*> extrap_cdr_fluxes;
  Vector<EBAMRIVData*> extrap_cdr_densities;
  Vector<EBAMRIVData*> extrap_cdr_velocities;
  Vector<EBAMRIVData*> extrap_rte_fluxes;

  cdr_fluxes = m_cdr->get_ebflux();

  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
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
  this->extrapolate_to_eb(extrap_cdr_densities,  m_cdr->get_phase(), cdr_densities);
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
			   E);
}

void splitstep_euler_f::compute_sigma_flux_at_start_of_time_step(){
  CH_TIME("splitstep_euler_f::compute_sigma_flux_at_start_of_time_step");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::compute_sigma_flux_at_start_of_time_step" << endl;
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

void splitstep_euler_f::set_cdr_source_to_zero_at_start_of_time_step(){
  CH_TIME("splitstep_euler_f::set_cdr_source_to_zero_at_start_of_time_step");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::set_cdr_source_to_zero_at_start_of_time_step" << endl;
  }

  m_cdr->set_source(0.0);
}

void splitstep_euler_f::advance_cdr_transport(const Real a_dt){
  CH_TIME("splitstep_euler_f::advance_cdr_transport");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::advance_cdr_transport" << endl;
  }

  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& phi         = storage->get_phi();
    EBAMRCellData& rhs         = storage->get_scratch();
    const EBAMRCellData& state = solver->get_state();


    solver->compute_rhs(rhs, state, a_dt);

    data_ops::set_value(phi, 0.0);
    data_ops::incr(phi, state,   1.0);
    data_ops::incr(phi, rhs, a_dt);

    m_amr->average_down(phi, m_cdr->get_phase());
    m_amr->interp_ghost(phi, m_cdr->get_phase());

    data_ops::floor(phi, 0.0);
  }
}

void splitstep_euler_f::advance_sigma_transport(const Real a_dt){
  CH_TIME("splitstep_euler_f::advance_sigma_transport");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::advance_sigma_transport" << endl;
  }

  const EBAMRIVData& state = m_sigma->get_state();
  
  EBAMRIVData& rhs  = m_sigma_scratch->get_scratch();
  EBAMRIVData& phi  = m_sigma_scratch->get_phi();
  m_sigma->compute_rhs(rhs);
  data_ops::set_value(phi,   0.0);
  data_ops::incr(phi, state,   1.0);
  data_ops::incr(phi, rhs, a_dt);

  m_amr->average_down(phi, m_cdr->get_phase());
  
  m_sigma->reset_cells(phi);
}

void splitstep_euler_f::solve_poisson_transport(){
  CH_TIME("splitstep_euler_f::solve_poisson_transport");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::solve_poisson_transport" << endl;
  }

  MFAMRCellData& scratch_pot = m_poisson_scratch->get_phi();

  // Sigma and rho's to solve for
  const EBAMRIVData& sigma         = m_sigma_scratch->get_phi();
  Vector<EBAMRCellData*> cdr_densities;
  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    cdr_densities.push_back(&(storage->get_phi()));
  }

  // Copy so we don't iterate from scratch
  data_ops::set_value(scratch_pot, 0.0);
  data_ops::incr(scratch_pot, m_poisson->get_state(), 1.0);


  // Solve
  if((m_step + 1) % m_fast_poisson == 0){
    const bool converged = this->solve_poisson(scratch_pot, m_poisson->get_source(), cdr_densities, sigma, centering::cell_center);
    
    if(!converged){
      pout() << "splitstep_euler_f::solve_poisson_transport - solver did not converge at step = " << m_step << endl;
    }
  }
}

void splitstep_euler_f::advance_rte_transport_stationary(){
  CH_TIME("splitstep_euler_f::advance_rte_transport_stationary");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::advance_rte_transport_stationary" << endl;
  }

  // TLDR: The transport-advanced states lies in the storage containers. Use those in the computation of the source
  //       terms and solve the RTE

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
    this->solve_rte(rte_states, rte_sources, cdr_states, E, dummy_dt, centering::cell_center);
  }
}

void splitstep_euler_f::advance_rte_transport_transient(const Real a_dt){
  CH_TIME("splitstep_euler_f::advance_rte_transport_transient");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::advance_rte_transport_transient" << endl;
  }

  MayDay::Abort("splitstep_euler_f::advance_rte_transport_transient - not implemented (yet)");
}

void splitstep_euler_f::compute_E_after_transport(){
  CH_TIME("splitstep_euler_f::compute_E_after_transport");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::compute_E_after_transport";
  }

  EBAMRCellData& E_cell = m_poisson_scratch->get_E_cell();
  EBAMRFluxData& E_face = m_poisson_scratch->get_E_face();
  EBAMRIVData&   E_eb   = m_poisson_scratch->get_E_eb();

  const MFAMRCellData& phi = m_poisson_scratch->get_phi();
  
  this->compute_E(E_cell, m_cdr->get_phase(), phi);     // Compute cell-centered field
  this->compute_E(E_face, m_cdr->get_phase(), E_cell);  // Compute face-centered field
  this->compute_E(E_eb,   m_cdr->get_phase(), E_cell);  // EB-centered field
}

void splitstep_euler_f::compute_cdr_source_after_transport(){
  CH_TIME("splitstep_euler_f::compute_cdr_velo_after_transport");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::compute_cdr_velo_after_transport";
  }

  // TLDR; The transport-advanced states lie in the scratch data holders for cdr, rte, and poisson. Use those
  //       for computing the source term
  
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

  this->compute_cdr_sources(cdr_sources, cdr_states, rte_states, E, centering::cell_center);
}

void splitstep_euler_f::advance_cdr_source(const Real a_dt){
  CH_TIME("splitstep_euler_f::advance_cdr_source");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::advance_cdr_source";
  }

  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& state   = solver->get_state();
    EBAMRCellData& source  = solver->get_source();
    EBAMRCellData& phi     = storage->get_phi();

    data_ops::set_value(state,    0.0);
    data_ops::incr(state, phi,    1.0);
    data_ops::incr(state, source, a_dt);

    m_amr->average_down(state, m_cdr->get_phase());
    m_amr->interp_ghost(state, m_cdr->get_phase());

    data_ops::floor(state, 0.0);
  }
}

void splitstep_euler_f::advance_sigma_source(const Real a_dt){
  CH_TIME("splitstep_euler_f::advance_sigma_source");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::advance_sigma_source";
  }

  EBAMRIVData& state = m_sigma->get_state();
  const EBAMRIVData& phi = m_sigma_scratch->get_phi();

  data_ops::set_value(state, 0.0);
  data_ops::incr(state, phi, 1.0);

  m_amr->average_down(state, m_cdr->get_phase());
  
  m_sigma->reset_cells(state);
}

void splitstep_euler_f::solve_poisson_source(){
  CH_TIME("splitstep_euler_f::solve_poisson_source");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::solve_poisson_source";
  }

  // We've computed an intermediate potential using the advective stuff. Copy that into the solver so that we can shave
  // off a few MG iterations
  MFAMRCellData& pot     = m_poisson->get_state();
  MFAMRCellData& phi     = m_poisson_scratch->get_phi();
  MFAMRCellData& scratch = m_poisson_scratch->get_scratch_phi();
  
  // For transient RTE solvers I need the source term at half time steps. Since the internal state
  // inside the solver will be overwritten, I take a backup into poisson_storage.scratch_phi
  if(!m_rte->is_stationary()){
    data_ops::set_value(scratch, 0.0);
    data_ops::incr(scratch, pot, 1.0); 
  }

  if((m_step + 1) % m_fast_poisson == 0){
    data_ops::set_value(pot, 0.0);
    data_ops::incr(pot, phi, 1.0);
    
    const bool converged = this->solve_poisson(pot,
					       m_poisson->get_source(),
					       m_cdr->get_states(),
					       m_sigma->get_state(),
					       centering::cell_center);
    if(!converged){
      pout() << "splitstep_euler_f::solve_poisson_source - solver did not converge at step = " << m_step << endl;
    }
  }
  else{
    data_ops::set_value(pot, 0.0);
    data_ops::incr(pot, phi, 1.0);
  }
}

void splitstep_euler_f::compute_E_after_source(){
  CH_TIME("splitstep_euler_f::compute_E_after_source");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::compute_E_after_source" << endl;
  }

  EBAMRCellData& E_cell = m_poisson_scratch->get_E_cell();
  EBAMRFluxData& E_face = m_poisson_scratch->get_E_face();
  EBAMRIVData&   E_eb   = m_poisson_scratch->get_E_eb();

  const MFAMRCellData& phi = m_poisson->get_state();
  
  this->compute_E(E_cell, m_cdr->get_phase(), phi);     // Compute cell-centered field
  this->compute_E(E_face, m_cdr->get_phase(), E_cell);  // Compute face-centered field
  this->compute_E(E_eb,   m_cdr->get_phase(), E_cell);  // EB-centered field
}

void splitstep_euler_f::advance_rte_source_stationary(){
  CH_TIME("splitstep_euler_f::advance_rte_source_stationary");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::advance_rte_source_stationary";
  }

  // When we arrive here, all the solvers have been filled.

  Vector<EBAMRCellData*> rte_states;
  Vector<EBAMRCellData*> rte_sources;
  Vector<EBAMRCellData*> cdr_states;

  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);
    
    EBAMRCellData& phi    = storage->get_phi();
    EBAMRCellData& state  = solver->get_state();
    EBAMRCellData& source = solver->get_source();

    
    // We computed the intermediate RTE solution in the transport step. Copy that solution over here so our iterative methods
    // can bring it closer to 
    // the time step. The result for this is y_extrap = y_alpha/alpha - y_0*(1-alpha)/alpha. This (usually) brings the
    // initial guess closer to the true solution.
    data_ops::set_value(state, 0.0);
    data_ops::incr(state, phi, 1.0);
    
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
    this->solve_rte(rte_states, rte_sources, cdr_states, E, dummy_dt, centering::cell_center);
  }
}

void splitstep_euler_f::advance_rte_source_transient(const Real a_dt){
  CH_TIME("splitstep_euler_f::advance_rte_source_transient");
  if(m_verbosity > 5){
    pout() << "splitstep_euler_f::avance_rte_source_transient";
  }

  MayDay::Abort("splitstep_euler_f::advance_rte_transport_transient - not implemented (yet)");
}
