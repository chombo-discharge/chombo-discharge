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

#define EULER_MARUYAMA_TIMER 1

euler_maruyama::euler_maruyama(){
  m_class_name = "euler_maruyama";
  m_extrap_advect = true;
}

euler_maruyama::~euler_maruyama(){

}

void euler_maruyama::parse_options(){
  CH_TIME("euler_maruyama::parse_options");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::parse_options" << endl;
  }

  parse_verbosity();
  parse_solver_verbosity();
  parse_fast_poisson();
  parse_cfl();
  parse_relax_time();
  parse_min_dt();
  parse_max_dt();
  parse_source_comp();
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

  Real t_fil1 = 0.0;
  Real t_reac = 0.0;
  Real t_cdr  = 0.0;
  Real t_rte  = 0.0;
  Real t_sig  = 0.0;
  Real t_pois = 0.0;
  Real t_fil2 = 0.0;
  Real t_tot  = 0.0;

  Real t0, t1;
  t0 = MPI_Wtime();
  t_tot  = -t0;

  euler_maruyama::compute_E_into_scratch();       // Compute the electric field
  euler_maruyama::compute_cdr_eb_states();        // Extrapolate cell-centered stuff to EB centroids
  euler_maruyama::compute_cdr_eb_fluxes();        // Extrapolate cell-centered fluxes to EB centroids
  euler_maruyama::compute_cdr_domain_states();    // Extrapolate cell-centered states to domain edges
  euler_maruyama::compute_cdr_domain_fluxes();    // Extrapolate cell-centered fluxes to domain edges
  euler_maruyama::compute_sigma_flux();           // Update charge flux for sigma solver
  t1 = MPI_Wtime();
  t_fil1 = t1 - t0;

  t0 = MPI_Wtime();
  euler_maruyama::compute_reaction_network(a_dt); // Advance the reaction network
  t1 = MPI_Wtime();
  t_reac = t1-t0;

  t0 = MPI_Wtime();
  euler_maruyama::advance_cdr(a_dt);              // Update cdr equations
  t1 = MPI_Wtime();
  t_cdr = t1 - t0;

  t0 = MPI_Wtime();
  euler_maruyama::advance_rte(a_dt);              // Update RTE equations
  t1 = MPI_Wtime();
  t_rte = t1-t0;

  t0 = MPI_Wtime();
  euler_maruyama::advance_sigma(a_dt);            // Update sigma equation
  t1 = MPI_Wtime();
  t_sig = t1 - t0;
  
  t0 = MPI_Wtime();
  if((m_step +1) % m_fast_poisson == 0){
    time_stepper::solve_poisson();                  // Update the Poisson equation
  }
  t1 = MPI_Wtime();
  t_pois = t1 - t0;

  t0 = MPI_Wtime();
  euler_maruyama::compute_E_into_scratch();       // Update electric fields too

  // Update velocities and diffusion coefficients. We don't do sources here. 
  euler_maruyama::compute_cdr_velo(m_time + a_dt);
  euler_maruyama::compute_cdr_diffco(m_time + a_dt);
  t1 = MPI_Wtime();
  t_fil2 = t1 - t0;
  t_tot += t1;

#if EULER_MARUYAMA_TIMER
  pout() << endl;
  pout() << "euler_maruyama::advance breakdown:" << endl
	 << "BC fill   = " << 100.0*t_fil1/t_tot << "%" << endl
	 << "Reactions = " << 100.*t_reac/t_tot << "%" << endl
	 << "CDR adv.  = " << 100.*t_cdr/t_tot << "%" << endl
	 << "RTE adv.  = " << 100.*t_rte/t_tot << "%" << endl
	 << "Poisson   = " << 100.*t_pois/t_tot << "%" << endl
	 << "Vel/Dco   = " << 100.*t_fil2/t_tot << "%" << endl
	 << "TOTAL = " << t_tot << "seconds" << endl;
  pout() << endl;
#endif
  
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

  Vector<EBAMRCellData*> cdr_states = m_cdr->get_states();
  
  Vector<EBAMRIVData*>   eb_gradients;
  Vector<EBAMRIVData*>   eb_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    RefCountedPtr<cdr_storage>& storage = euler_maruyama::get_cdr_storage(solver_it);

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

void euler_maruyama::compute_cdr_domain_states(){
  CH_TIME("euler_maruyama::compute_cdr_domain_states");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::compute_cdr_domain_states" << endl;
  }

  Vector<EBAMRIFData*>   domain_gradients;
  Vector<EBAMRIFData*>   domain_states;
  Vector<EBAMRCellData*> cdr_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (auto solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    RefCountedPtr<cdr_storage>& storage = euler_maruyama::get_cdr_storage(solver_it);

    cdr_states.push_back(&(solver->get_state()));
    domain_states.push_back(&(storage->get_domain_state()));
    domain_gradients.push_back(&(storage->get_domain_grad()));
    cdr_gradients.push_back(&(storage->get_gradient())); // Should already be computed
  }

  // Extrapolate states to the domain faces
  time_stepper::extrapolate_to_domain_faces(domain_states, m_cdr->get_phase(), cdr_states);

  // We already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIFData grad;
  m_amr->allocate(grad, m_cdr->get_phase(), SpaceDim);
  for (int i = 0; i < cdr_states.size(); i++){
    time_stepper::extrapolate_to_domain_faces(grad, m_cdr->get_phase(), *cdr_gradients[i]);
    time_stepper::project_domain(*domain_gradients[i], grad);
  }
}

void euler_maruyama::compute_cdr_domain_fluxes(){
  CH_TIME("euler_maruyama::compute_cdr_domain_fluxes()");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::compute_cdr_domain_fluxes()" << endl;
  }

  Vector<EBAMRCellData*> states = m_cdr->get_states();

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
  this->extrapolate_to_domain_faces(extrap_cdr_densities,         m_cdr->get_phase(), states);
  this->extrapolate_vector_to_domain_faces(extrap_cdr_velocities, m_cdr->get_phase(), cdr_velocities);
  this->compute_extrapolated_domain_fluxes(extrap_cdr_fluxes,     states,             cdr_velocities, m_cdr->get_phase());
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
					  m_time);
}

void euler_maruyama::compute_sigma_flux(){
  CH_TIME("euler_maruyama::compute_sigma_flux");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::compute_sigma_flux" << endl;
  }

  EBAMRIVData& flux = m_sigma->get_flux();
  data_ops::set_value(flux, 0.0);

  for (auto solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();
    const RefCountedPtr<species>& spec      = solver_it.get_species();
    const EBAMRIVData& solver_flux          = solver->get_ebflux();

    data_ops::incr(flux, solver_flux, spec->get_charge()*units::s_Qe);
  }

  m_sigma->reset_cells(flux);
}

void euler_maruyama::compute_reaction_network(const Real a_dt){
  CH_TIME("euler_maruyama::compute_reaction_network");
  if(m_verbosity > 5){
    pout() << "euler_maruaya::compute_reaction_network" << endl;
  }

  time_stepper::advance_reaction_network(m_time, a_dt);
}

void euler_maruyama::advance_cdr(const Real a_dt){
  CH_TIME("euler_maruyama::advance_cdr");
  if(m_verbosity > 5){
    pout() << "euler_maruaya::advance_cdr" << endl;
  }

  Real t0, t1;

  Real t_divF = 0.0;
  Real t_sour = 0.0;
  Real t_diff = 0.0;
  Real t_sync = 0.0;
  Real t_tot  = -MPI_Wtime();
  for (auto solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = euler_maruyama::get_cdr_storage(solver_it);

    EBAMRCellData& phi = solver->get_state();
    EBAMRCellData& src = solver->get_source();
    
    EBAMRCellData& scratch  = storage->get_scratch();
    EBAMRCellData& scratch2 = storage->get_scratch2();

    // Compute hyperbolic term into scratch
    t0 = MPI_Wtime();
    if(solver->is_mobile()){
      const Real extrap_dt = m_extrap_advect ? a_dt : 0.0;
      solver->compute_divF(scratch, phi, extrap_dt, true);
      data_ops::scale(scratch, -1.0);
    }
    else{
      data_ops::set_value(scratch, 0.0);
    }
    t1 = MPI_Wtime();
    t_divF += t1-t0;

    // Also include the source term
    t0 = MPI_Wtime();
    data_ops::incr(scratch, src, 1.0);  // scratch = [-div(F) + R]
    data_ops::scale(scratch, a_dt);     // scratch = [-div(F) + R]*dt
    data_ops::incr(phi, scratch, 1.0);  // Make phi = phi^k - dt*div(F) + dt*R
    t1 = MPI_Wtime();
    t_sour += t1 - t0;
    data_ops::floor(phi, 0.0);


    // Solve diffusion equation. This looks weird but we're solving
    //
    // phi^(k+1) = phi^k - dt*div(F) + dt*R + dt*div(D*div(phi^k+1))
    //
    // This discretization is equivalent to a diffusion-only discretization with phi^k -dt*div(F) + dt*R as initial solution
    // so we just use that for simplicity
    t0 = MPI_Wtime();
    if(solver->is_diffusive()){
      data_ops::copy(scratch, phi); // Weird-ass initial solution, as explained above
      data_ops::set_value(scratch2, 0.0); // No source, those are a part of the initial solution
      solver->advance_euler(phi, scratch, scratch2, a_dt); 
    }
    t1 = MPI_Wtime();
    t_diff += t1 - t0;

    t0 = MPI_Wtime();
    data_ops::floor(phi, 0.0);
    m_amr->average_down(phi, m_cdr->get_phase());
    m_amr->interp_ghost(phi, m_cdr->get_phase());
    t1 = MPI_Wtime();
    t_sync += t1 - t0;

    
  }
  t_tot += MPI_Wtime();

#if EULER_MARUYAMA_TIMER
  pout() << endl;
  pout() << "euler_maruayama::advance_cdr breakdown:" << endl
	 << "divF       = " << 100.*t_divF/t_tot << "%" << endl
	 << "source     = " << 100.*t_sour/t_tot << "%" << endl
	 << "diffusion  = " << 100.*t_diff/t_tot << "%" << endl
	 << "avg/interp = " << 100.*t_sync/t_tot << "%" << endl
	 << "TOTAL = " << t_tot << "seconds" << endl;
  pout() << endl;
#endif
}

void euler_maruyama::advance_rte(const Real a_dt){
  CH_TIME("euler_maruyama::advance_rte");
  if(m_verbosity > 5){
    pout() << "euler_maruaya::advance_rte" << endl;
  }

  // Source terms should already be in place so we can solve directly.
  for (rte_iterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    solver->advance(a_dt);
  }
}

void euler_maruyama::advance_sigma(const Real a_dt){
  CH_TIME("euler_maruyama::advance_sigma");
  if(m_verbosity > 5){
    pout() << "euler_maruaya::advance_sigma" << endl;
  }

  // Advance the sigma equation
  EBAMRIVData& sigma = m_sigma->get_state();
  const EBAMRIVData& rhs = m_sigma->get_flux();
  data_ops::incr(sigma, rhs, a_dt);
}

void euler_maruyama::compute_cdr_velo(const Real a_time){
  CH_TIME("euler_maruyama::compute_cdr_velo");
  if(m_verbosity > 5){
    pout() << "euler_maruaya::compute_cdr_velo" << endl;
  }

  Vector<EBAMRCellData*> velocities = m_cdr->get_velocities();
  time_stepper::compute_cdr_velocities(velocities, m_cdr->get_states(), m_poisson_scratch->get_E_cell(), a_time);
}

void euler_maruyama::compute_cdr_diffco(const Real a_time){
  CH_TIME("euler_maruyama::compute_cdr_diffco");
  if(m_verbosity > 5){
    pout() << "euler_maruaya::compute_cdr_diffco" << endl;
  }

  time_stepper::compute_cdr_diffusion(m_poisson_scratch->get_E_cell(), m_poisson_scratch->get_E_eb());
}

void euler_maruyama::compute_dt(Real& a_dt, time_code::which_code& a_timecode){
  CH_TIME("euler_maruyama::compute_dt");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::compute_dt" << endl;
  }

  Real dt = 1.E99;

  m_dt_cfl          = m_cdr->compute_cfl_dt();
  const Real dt_cfl = m_cfl*m_dt_cfl;
  if(dt_cfl < dt){
    dt = dt_cfl;
    a_timecode = time_code::cfl;
  }

  const Real dt_relax = m_relax_time*this->compute_relaxation_time();
  if(dt_relax < dt){
    dt = dt_relax;
    a_timecode = time_code::relaxation_time;
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
