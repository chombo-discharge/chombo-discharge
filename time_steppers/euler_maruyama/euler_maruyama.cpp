/*!
  @file   euler_maruyama.cpp
  @brief  Implementation of euler_maruyama.H
  @author Robert Marskar
  @date   Aug. 2019
*/

#include "euler_maruyama.H"
#include "euler_maruyama_storage.H"
#include "data_ops.H"
#include "units.H"

typedef euler_maruyama::cdr_storage     cdr_storage;
typedef euler_maruyama::poisson_storage poisson_storage;
typedef euler_maruyama::rte_storage     rte_storage;
typedef euler_maruyama::sigma_storage   sigma_storage;

euler_maruyama::euler_maruyama(){

}

euler_maruyama::~euler_maruyama(){

}

bool euler_maruyama::need_to_regrid(){
  CH_TIME("euler_maruyama::deallocate_internals");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::need_to_regrid" << endl;
  }

  return false;
}

RefCountedPtr<cdr_storage>& euler_maruyama::get_cdr_storage(const cdr_iterator& a_solverit){
  return m_cdr_scratch[a_solverit.get_solver()];
}

RefCountedPtr<rte_storage>& euler_maruyama::get_rte_storage(const rte_iterator& a_solverit){
  return m_rte_scratch[a_solverit.get_solver()];
}

Real euler_maruyama::restrict_dt(){
  CH_TIME("euler_maruyama::euler_maruyama");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::euler_maruyama" << endl;
  }

  return 1.E99;
}

Real euler_maruyama::advance(const Real a_dt){
  if(m_verbosity > 5){
    pout() << "euler_maruyama::advance" << endl;
  }

  // INFO: Solvers should have been filled with velocities and diffusion coefficients. We must still do:
  // 1. Compute E
  // 2. Extrapolate everything to the EB
  // 3. Compute fluxes at the EB and domain
  // 4. Advance the reaction network. This provides source terms for CDR and RTE equations
  // 5. Compute the hyperbolic terms
  // 6. Solve the semi-implicit discretization
  // 7. Advance the RTE equations
  // 8. Update the Poisson equation
  // 9. Recompute solver velocities and diffusion coefficients

  euler_maruyama::compute_E_into_scratch();
  euler_maruyama::compute_cdr_eb_states(); // Extrapolate cell-centered stuff to EB centroids
  euler_maruyama::compute_cdr_eb_fluxes(); // Extrapolate cell-centered fluxes to EB centroids
  
  return a_dt;
}

void euler_maruyama::init_source_terms(){
  CH_TIME("euler_maruyama::init_source_terms");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::init_source_terms" << endl;
  }

  // No need to do anything in this routine yet
}

void euler_maruyama::regrid_internals(){
  CH_TIME("euler_maruyama::regrid_internals");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::regrid_internals" << endl;
  }

  const int ncomp       = 1;
  const int num_species = m_plaskin->get_num_species();
  const int num_photons = m_plaskin->get_num_photons();

  // Allocate cdr storage
  m_cdr_scratch.resize(num_species);
  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx] = RefCountedPtr<cdr_storage> (new cdr_storage(m_amr, m_cdr->get_phase(), ncomp));
    m_cdr_scratch[idx]->allocate_storage();
  }

  // Allocate RTE storage
  m_rte_scratch.resize(num_photons);
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx] = RefCountedPtr<rte_storage> (new rte_storage(m_amr, m_rte->get_phase(), ncomp));
    m_rte_scratch[idx]->allocate_storage();
  }

  // Allocate Poisson storage
  m_poisson_scratch = RefCountedPtr<poisson_storage> (new poisson_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_poisson_scratch->allocate_storage();
  
  // Allocate sigma storage
  m_sigma_scratch = RefCountedPtr<sigma_storage> (new sigma_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_sigma_scratch->allocate_storage();
}

void euler_maruyama::deallocate_internals(){
  CH_TIME("euler_maruyama::deallocate_internals");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::deallocate_internals" << endl;
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

void euler_maruyama::compute_E_into_scratch(){
  CH_TIME("euler_maruyama::compute_E_into_scratch");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::compute_E_into_scratch" << endl;
  }
  
  EBAMRCellData& E_cell = m_poisson_scratch->get_E_cell();
  EBAMRFluxData& E_face = m_poisson_scratch->get_E_face();
  EBAMRIVData&   E_eb   = m_poisson_scratch->get_E_eb();
  EBAMRIFData&   E_dom  = m_poisson_scratch->get_E_domain();

  const MFAMRCellData& phi = m_poisson->get_state();
  
  time_stepper::compute_E(E_cell, m_cdr->get_phase(), phi);     // Compute cell-centered field
  time_stepper::compute_E(E_face, m_cdr->get_phase(), E_cell);  // Compute face-centered field
  time_stepper::compute_E(E_eb,   m_cdr->get_phase(), E_cell);  // EB-centered field
  time_stepper::extrapolate_to_domain_faces(E_dom, m_cdr->get_phase(), E_cell); // Domain centered field
}

void euler_maruyama::compute_cdr_gradients(){
  CH_TIME("euler_maruyama::compute_cdr_gradients");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::compute_cdr_gradients" << endl;
  }

  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    RefCountedPtr<cdr_solver>& solver = solver_it();
    RefCountedPtr<cdr_storage>& storage = euler_maruyama::get_cdr_storage(solver_it);

    EBAMRCellData& grad = storage->get_gradient();
    m_amr->compute_gradient(grad, solver->get_state());
    m_amr->average_down(grad, m_cdr->get_phase());
    m_amr->interp_ghost(grad, m_cdr->get_phase());
  }
}

void euler_maruyama::compute_cdr_eb_states(){
  CH_TIME("euler_maruyama::compute_cdr_eb_states");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::compute_cdr_eb_states" << endl;
  }

  Vector<EBAMRIVData*>   eb_gradients;
  Vector<EBAMRIVData*>   eb_states;
  Vector<EBAMRCellData*> cdr_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    RefCountedPtr<cdr_storage>& storage = euler_maruyama::get_cdr_storage(solver_it);

    cdr_states.push_back(&(solver->get_state()));
    eb_states.push_back(&(storage->get_eb_state()));
    eb_gradients.push_back(&(storage->get_eb_grad()));
    cdr_gradients.push_back(&(storage->get_gradient())); // Should already have been computed
  }

  // Extrapolate states to the EB and floor them so we cannot get negative values on the boundary. This
  // won't hurt mass conservation because the mass hasn't been injected yet
  time_stepper::extrapolate_to_eb(eb_states, m_cdr->get_phase(), cdr_states);
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    data_ops::floor(*eb_states[idx], 0.0);
  }

  // We should already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIVData eb_gradient;
  m_amr->allocate(eb_gradient, m_cdr->get_phase(), SpaceDim);
  for (int i = 0; i < cdr_states.size(); i++){
    time_stepper::extrapolate_to_eb(eb_gradient, m_cdr->get_phase(), *cdr_gradients[i]);
    time_stepper::project_flux(*eb_gradients[i], eb_gradient);
  }
}

void euler_maruyama::compute_cdr_eb_fluxes(){
  CH_TIME("euler_maruyama::compute_cdr_eb_fluxes()");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::compute_cdr_eb_fluxes()";
  }

  Vector<EBAMRCellData*> states = m_cdr->get_states();

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


  // Compute extrapolated fluxes and velocities at the EB
  Vector<EBAMRCellData*> cdr_velocities = m_cdr->get_velocities();
  time_stepper::compute_extrapolated_fluxes(extrap_cdr_fluxes, states, cdr_velocities, m_cdr->get_phase());
  time_stepper::compute_extrapolated_velocities(extrap_cdr_velocities, cdr_velocities, m_cdr->get_phase());

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
				   m_time);
}


