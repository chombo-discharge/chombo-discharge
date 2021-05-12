/*!
  @file   rk2_tga.cpp
  @brief  Implementation of rk2_tga.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "rk2_tga.H"
#include "rk2_tga_storage.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"
#include "data_ops.H"
#include "units.H"

#include <ParmParse.H>

typedef rk2_tga::cdr_storage     cdr_storage;
typedef rk2_tga::poisson_storage poisson_storage;
typedef rk2_tga::rte_storage     rte_storage;
typedef rk2_tga::sigma_storage   sigma_storage;

rk2_tga::rk2_tga(){
  m_alpha = 1.0;

  // Basically only for debugging
  m_do_advec_src = true;
  m_do_diffusion = true;
  m_do_rte       = true;
  m_do_poisson   = true;
  
  {
    ParmParse pp("rk2_tga");

    std::string str;
    
    pp.query("rk2_tga_alpha", m_alpha);

    if(pp.contains("turn_off_advection")){
      pp.get("turn_off_advection_source", str);
      if(str == "true"){
	m_do_advec_src = false;
	if(m_verbosity > 2){
	  pout() << "rk2_tga::rk2_tga - Turning off advection & source" << endl;
	}
      }
    }
    if(pp.contains("turn_off_diffusion")){
      pp.get("turn_off_diffusion", str);
      if(str == "true"){
	m_do_diffusion = false;
	if(m_verbosity > 2){
	  pout() << "rk2_tga::rk2_tga - Turning off diffusion" << endl;
	}
      }
    }
    if(pp.contains("turn_off_rte")){
      pp.get("turn_off_rte", str);
      if(str == "true"){
	m_do_rte = false;

	if(m_verbosity > 2){
	  pout() << "rk2_tga::rk2_tga - Turning off rte" << endl;
	}
      }
    }
    if(pp.contains("turn_off_poisson")){
      pp.get("turn_off_poisson", str);
      if(str == "true"){
	m_do_poisson = false;

	if(m_verbosity > 2){
	  pout() << "rk2_tga::rk2_tga - Turning off poisson" << endl;
	}
      }
    }
  }
}

rk2_tga::~rk2_tga(){

}

RefCountedPtr<cdr_storage>& rk2_tga::get_cdr_storage(const cdr_iterator& a_solverit){
  return m_cdr_scratch[a_solverit.get_solver()];
}

RefCountedPtr<rte_storage>& rk2_tga::get_rte_storage(const rte_iterator& a_solverit){
  return m_rte_scratch[a_solverit.get_solver()];
}

Real rk2_tga::restrict_dt(){
  return 1.E99;
}

Real rk2_tga::advance(const Real a_dt){
  CH_TIME("rk2_tga::advance");
  if(m_verbosity > 2){
    pout() << "rk2_tga::advance" << endl;
  }

  // When we enter this routine, we assume that velocities, source terms, and diffusion coefficients have already
  // been computed. If you change the order of the integration, you must ensure that appropriate velocities, diffusion
  // coefficients, and source terms, are updated where they should. So please don't do that...

  this->cache_solutions(); // Cache old solutions. Not used for anything (yet), but derived classes might.
                           // This means that everything we do here is done directly in the solver states. 

  if(m_do_advec_src){
    this->advance_advection_source(a_dt);
  }
  if(m_do_diffusion){
    this->advance_diffusion(a_dt);
  }

  if(m_do_rte){ // Solve for final RTE stage. Poisson equation should already have been solved at the end of the advection
    this->advance_rte_stationary(m_time + a_dt); // and diffusion stages
  }

  // Put cdr solvers back in useable state so that we can reliably compute the next time step. 
  this->compute_cdr_velocities();
  this->compute_cdr_diffusion();
  this->compute_cdr_sources();
  
  return a_dt;
}

void rk2_tga::regrid_internals(){
  CH_TIME("time_stepper::regrid_internals");
  if(m_verbosity > 5){
    pout() << "time_stepper::regrid_internals" << endl;
  }
  
  this->allocate_cdr_storage();
  this->allocate_poisson_storage();
  this->allocate_rte_storage();
  this->allocate_sigma_storage();
}

void rk2_tga::compute_dt(Real& a_dt, time_code::which_code& a_timecode){
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

void rk2_tga::allocate_cdr_storage(){
  const int ncomp       = 1;
  const int num_species = m_plaskin->get_num_species();
  m_cdr_scratch.resize(num_species);
  
  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx] = RefCountedPtr<cdr_storage> (new cdr_storage(m_amr, m_cdr->get_phase(), ncomp));
    m_cdr_scratch[idx]->allocate_storage();
  }
}

void rk2_tga::allocate_poisson_storage(){
  const int ncomp = 1;
  m_fieldSolver_scratch = RefCountedPtr<poisson_storage> (new poisson_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_fieldSolver_scratch->allocate_storage();
}

void rk2_tga::allocate_rte_storage(){
  const int ncomp       = 1;
  const int num_photons = m_plaskin->get_num_photons();
  m_rte_scratch.resize(num_photons);
  
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx] = RefCountedPtr<rte_storage> (new rte_storage(m_amr, m_rte->get_phase(), ncomp));
    m_rte_scratch[idx]->allocate_storage();
  }
}

void rk2_tga::allocate_sigma_storage(){
  const int ncomp = 1;
  m_sigma_scratch = RefCountedPtr<sigma_storage> (new sigma_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_sigma_scratch->allocate_storage();
}

void rk2_tga::deallocate_internals(){
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

  m_fieldSolver_scratch->deallocate_storage();
  m_sigma_scratch->deallocate_storage();
}

void rk2_tga::cache_solutions(){
  CH_TIME("rk2_tga::cache_solutions");
  if(m_verbosity > 5){
    pout() << "rk2_tga::cache_solutions" << endl;
  }
  
  // Cache cdr solutions
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();

    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);
    EBAMRCellData& cache = storage->get_cache();

    data_ops::copy(cache, solver->get_state());
  }

  {// Cache Poisson solution
    MFAMRCellData& cache = m_fieldSolver_scratch->get_cache();
    data_ops::copy(cache, m_fieldSolver->getPotential());
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

void rk2_tga::advance_advection_source(const Real a_dt){
  CH_TIME("rk2_tga::advance_advection_source");
  if(m_verbosity > 5){
    pout() << "rk2_tga::advance_advection_source" << endl;
  }

  const Real t0 = m_time;
  const Real t1 = m_time + m_alpha*a_dt;

  // Compute necessary things for k1 advance
  this->compute_E_into_scratch();             // Electric field
  this->compute_cdr_eb_states();              // Compute extrapolation n and grad(n) on the EB
  this->compute_cdr_fluxes(t0);               // Compute EB fluxes
  this->compute_sigma_flux_into_scratch();    // Compute sum of EB fluxes

  // Do k1 advance
  this->advance_advection_source_cdr_k1(a_dt);// First RK stage advance. Make phi = phi + k1*alpha*dt, phi being the solver state
  this->advance_advection_sigma_k1(a_dt);     // First RK stage advance. Make phi = phi + k1*alpha*dt, phi being the solver state
  if(m_do_poisson){
    this->solve_poisson();                    // Solvers contain the intermediate states, resolve Poisson
  }
  this->compute_E_into_scratch();             // Recompute E
  if(m_do_rte){                               // Do the RTE solve in order to get new source terms
    this->advance_rte_stationary(t1);
  }

  // Recompute things in order to do k2 advance
  this->compute_cdr_sources_into_scratch(t1);
  this->compute_cdr_velo(t1);
  this->compute_cdr_eb_states();
  this->compute_cdr_fluxes(t1);
  this->compute_sigma_flux_into_scratch();

  // Do k2 advance
  this->advance_advection_source_cdr_k2(a_dt);
  this->advance_advection_sigma_k2(a_dt);
  if(m_do_poisson){ // Solve Poisson, but don't solve RTE because we should do that after implicit diffusion
    this->solve_poisson();
  }
}

void rk2_tga::compute_E_into_scratch(){
  CH_TIME("rk2_tga::compute_E_into_scratch");
  if(m_verbosity > 5){
    pout() << "rk2_tga::compute_E_into_scratch" << endl;
  }
  
  EBAMRCellData& E_cell = m_fieldSolver_scratch->get_E_cell();
  EBAMRFluxData& E_face = m_fieldSolver_scratch->get_E_face();
  EBAMRIVData&   E_eb   = m_fieldSolver_scratch->get_E_eb();

  const MFAMRCellData& phi = m_fieldSolver->getPotential();
  
  this->compute_E(E_cell, m_cdr->get_phase(), phi);     // Compute cell-centered field
  this->compute_E(E_face, m_cdr->get_phase(), E_cell);  // Compute face-centered field
  this->compute_E(E_eb,   m_cdr->get_phase(), E_cell);  // EB-centered field
}

void rk2_tga::compute_cdr_velo(const Real a_time){
  CH_TIME("rk2_tga::compute_cdr_velo");
  if(m_verbosity > 5){
    pout() << "splitstep_::compute_cdr_velo" << endl;
  }

  Vector<EBAMRCellData*> states     = m_cdr->get_states();
  Vector<EBAMRCellData*> velocities = m_cdr->get_velocities();
  this->compute_cdr_velocities(velocities, states, m_fieldSolver_scratch->get_E_cell(), a_time);
}

void rk2_tga::compute_cdr_eb_states(){
  CH_TIME("rk2_tga::compute_cdr_eb_states");
  if(m_verbosity > 5){
    pout() << "rk2_tga::compute_cdr_eb_states" << endl;
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
  this->computeGradients_at_eb(eb_gradients, m_cdr->get_phase(), cdr_states);
}

void rk2_tga::compute_cdr_fluxes(const Real a_time){
  CH_TIME("rk2_tga::compute_cdr_fluxes");
  if(m_verbosity > 5){
    pout() << "rk2_tga::compute_cdr_fluxes" << endl;
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

  const EBAMRIVData& E = m_fieldSolver_scratch->get_E_eb();

  time_stepper::compute_cdr_fluxes(cdr_fluxes,
				   extrap_cdr_fluxes,
				   extrap_cdr_densities,
				   extrap_cdr_velocities,
				   extrap_cdr_gradients,
				   extrap_rte_fluxes,
				   E,
				   a_time);
}

void rk2_tga::compute_sigma_flux_into_scratch(){
  CH_TIME("rk2_tga::compute_sigma_flux_into_scratch");
  if(m_verbosity > 5){
    pout() << "rk2_tga::compute_sigma_flux_into_scratch" << endl;
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

void rk2_tga::advance_advection_source_cdr_k1(const Real a_dt){
  CH_TIME("rk2_tga::advance_advection_source_cdr_k1");
  if(m_verbosity > 5){
    pout() << "rk2_tga::advance_advection_source_cdr_k1" << endl;
  }

  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& k1  = storage->get_k1();
    EBAMRCellData& phi = solver->get_state();
    EBAMRCellData& src = solver->get_source();

    // Compute rhs
    data_ops::set_value(k1, 0.0);
    solver->compute_divF(k1, phi, 0.0, true);
    data_ops::scale(k1, -1.0); 
    data_ops::incr(k1, src, 1.0);

    // Make phi = phi + k1*alpha*dt
    data_ops::incr(phi, k1, m_alpha*a_dt);

    m_amr->averageDown(phi, m_cdr->get_phase());
    m_amr->interpGhost(phi, m_cdr->get_phase());

    data_ops::floor(phi, 0.0);
  }
}

void rk2_tga::advance_advection_sigma_k1(const Real a_dt){
  CH_TIME("rk2_tga::advance_advection_sigma_k1");
  if(m_verbosity > 5){
    pout() << "rk2_tga::advance_advection_sigma_k1" << endl;
  }

  EBAMRIVData& phi = m_sigma->get_state();
  EBAMRIVData& k1  = m_sigma_scratch->get_k1();
  
  m_sigma->compute_rhs(k1);

  // Make phi = phi + k1*alpha*dt
  data_ops::incr(phi, k1,    m_alpha*a_dt);

  m_amr->averageDown(phi, m_cdr->get_phase());
  
  m_sigma->reset_cells(k1);
  m_sigma->reset_cells(phi);
}

void rk2_tga::advance_advection_source_cdr_k2(const Real a_dt){
  CH_TIME("rk2_tga::advance_advection_source_cdr_k2");
  if(m_verbosity > 5){
    pout() << "rk2_tga::advance_advection_source_cdr_k2" << endl;
  }

  for (cdr_iterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = this->get_cdr_storage(solver_it);

    EBAMRCellData& state       = solver->get_state();
    EBAMRCellData& k2          = storage->get_k2();
    const EBAMRCellData& k1    = storage->get_k1();
    EBAMRCellData& src         = solver->get_source();

    // Compute RHS
    data_ops::set_value(k2, 0.0);
    solver->compute_divF(k2, state, 0.0, true);
    data_ops::scale(k2, -1.0);
    data_ops::incr(k2, src, 1.0);

    // RK2 advance. The extract m_alpha subtraction is because when we came here, the solver state (which we update in place)
    // contained the intermediate state phi + k1*alpha_dt. But we want phi = phi + k1*a_dt*(1-1/(2*alpha)) + k2*dt/(2*alpha),
    // so we just adjust the factor directly.
    const Real k1_factor = a_dt*(1.0 - 1.0/(2.0*m_alpha) - m_alpha);
    const Real k2_factor = a_dt/(2.0*m_alpha);
    
    data_ops::incr(state, k1, k1_factor);
    data_ops::incr(state, k2, k2_factor);

    m_amr->averageDown(state, m_cdr->get_phase());
    m_amr->interpGhost(state, m_cdr->get_phase());

    data_ops::floor(state, 0.0);
  }
}

void rk2_tga::advance_advection_sigma_k2(const Real a_dt){
  CH_TIME("rk2_tga::advance_advection_sigma_k2");
  if(m_verbosity > 5){
    pout() << "rk2_tga::advance_advection_sigma_k2" << endl;
  }
  
  const EBAMRIVData& k1  = m_sigma_scratch->get_k1();
  EBAMRIVData& k2        = m_sigma_scratch->get_k2();
  m_sigma->compute_rhs(k2);

  EBAMRIVData& state       = m_sigma->get_state();

  // RK2 advance. The extract m_alpha subtraction is because when we came here, the solver state (which we update in place)
  // contained the intermediate state phi + k1*alpha_dt. But we want phi = phi + k1*a_dt*(1-1/(2*alpha)) + k2*dt/(2*alpha),
  // so we just adjust the factor directly.
  const Real k1_factor = a_dt*(1.0 - 1.0/(2.0*m_alpha) - m_alpha);
  const Real k2_factor = a_dt/(2.0*m_alpha);
    
  data_ops::incr(state, k1, k1_factor);
  data_ops::incr(state, k2, k2_factor);

  m_amr->averageDown(state, m_cdr->get_phase());
  m_sigma->reset_cells(state);
}

void rk2_tga::advance_diffusion(const Real a_dt){
  CH_TIME("rk2_tga::advance_diffusion");
  if(m_verbosity > 5){
    pout() << "rk2_tga::advance_diffusion" << endl;
  }

  // Diffusive advance for all cdr equations
  bool diffusive_states = false;
  for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<cdr_solver>& solver = solver_it();

    if(solver->is_diffusive()){
      diffusive_states = true;
    }
  }

  // Do the diffusion advance
  if(diffusive_states){
    m_cdr->set_source(0.0); // This is necessary because advance_diffusion also works with source terms
    this->compute_cdr_diffusion();
    for (cdr_iterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<cdr_solver>& solver = solver_it();

      solver->advance_diffusion(a_dt);
    }

    // Update poisson equation afterwards
    if(m_do_poisson){
      this->solve_poisson();
    }
  }
}

void rk2_tga::compute_cdr_sources_into_scratch(const Real a_time){
  CH_TIME("rk2_tga::compute_cdr_sources_into_scratch");
  if(m_verbosity > 5){
    pout() << "rk2_tga::compute_cdr_sources_into_scratch" << endl;
  }
  
  Vector<EBAMRCellData*> cdr_sources = m_cdr->get_sources();
  Vector<EBAMRCellData*> cdr_states  = m_cdr->get_states();
  Vector<EBAMRCellData*> rte_states  = m_rte->get_states();
  EBAMRCellData& E                   = m_fieldSolver_scratch->get_E_cell();

  this->compute_cdr_sources(cdr_sources, cdr_states, rte_states, E, a_time, centering::cell_center);
}

void rk2_tga::advance_rte_stationary(const Real a_time){
  CH_TIME("rk2_tga::compute_rte_k1_stationary");
  if(m_verbosity > 5){
    pout() << "rk2_tga::compute_k1_stationary" << endl;
  }

  if((m_step + 1) % m_fast_rte == 0){
    Vector<EBAMRCellData*> rte_states  = m_rte->get_states();
    Vector<EBAMRCellData*> rte_sources = m_rte->get_sources();
    Vector<EBAMRCellData*> cdr_states  = m_cdr->get_states();

    EBAMRCellData& E = m_fieldSolver_scratch->get_E_cell();

    const Real dummy_dt = 0.0;
    this->solve_rte(rte_states, rte_sources, cdr_states, E, a_time, dummy_dt, centering::cell_center);
  }
#include "CD_NamespaceFooter.H"
