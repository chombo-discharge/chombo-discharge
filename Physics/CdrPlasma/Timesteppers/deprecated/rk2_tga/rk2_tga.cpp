/*!
  @file   rk2_tga.cpp
  @brief  Implementation of rk2_tga.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "rk2_tga.H"
#include "rk2_tga_storage.H"
#include <CD_CdrIterator.H>
#include <CD_RtIterator.H>
#include <CD_DataOps.H>
#include <CD_Units.H>

#include <ParmParse.H>

typedef rk2_tga::CdrStorage     CdrStorage;
typedef rk2_tga::FieldStorage FieldStorage;
typedef rk2_tga::RtStorage     RtStorage;
typedef rk2_tga::SigmaStorage   SigmaStorage;

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

RefCountedPtr<CdrStorage>& rk2_tga::get_CdrStorage(const CdrIterator& a_solverit){
  return m_cdr_scratch[a_solverit.get_solver()];
}

RefCountedPtr<RtStorage>& rk2_tga::get_RtStorage(const RtIterator& a_solverit){
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
  this->computeCdrDriftVelocities();
  this->compute_cdr_diffusion();
  this->compute_cdr_sources();
  
  return a_dt;
}

void rk2_tga::regridInternals(){
  CH_TIME("TimeStepper::regridInternals");
  if(m_verbosity > 5){
    pout() << "TimeStepper::regridInternals" << endl;
  }
  
  this->allocateCdrStorage();
  this->allocateFieldStorage();
  this->allocateRtStorage();
  this->allocateSigmaStorage();
}

void rk2_tga::computeDt(Real& a_dt, TimeCode::which_code& a_timeCode){
  CH_TIME("TimeStepper::computeDt");
  if(m_verbosity > 5){
    pout() << "TimeStepper::computeDt" << endl;
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
}

void rk2_tga::allocateCdrStorage(){
  const int ncomp       = 1;
  const int num_species = m_plaskin->get_num_species();
  m_cdr_scratch.resize(num_species);
  
  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx] = RefCountedPtr<CdrStorage> (new CdrStorage(m_amr, m_cdr->getPhase(), ncomp));
    m_cdr_scratch[idx]->allocate_storage();
  }
}

void rk2_tga::allocateFieldStorage(){
  const int ncomp = 1;
  m_fieldSolver_scratch = RefCountedPtr<FieldStorage> (new FieldStorage(m_amr, m_cdr->getPhase(), ncomp));
  m_fieldSolver_scratch->allocate_storage();
}

void rk2_tga::allocateRtStorage(){
  const int ncomp       = 1;
  const int num_Photons = m_plaskin->get_num_Photons();
  m_rte_scratch.resize(num_Photons);
  
  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx] = RefCountedPtr<RtStorage> (new RtStorage(m_amr, m_rte->getPhase(), ncomp));
    m_rte_scratch[idx]->allocate_storage();
  }
}

void rk2_tga::allocateSigmaStorage(){
  const int ncomp = 1;
  m_sigma_scratch = RefCountedPtr<SigmaStorage> (new SigmaStorage(m_amr, m_cdr->getPhase(), ncomp));
  m_sigma_scratch->allocate_storage();
}

void rk2_tga::deallocateInternals(){
  CH_TIME("TimeStepper::deallocateInternals");
  if(m_verbosity > 5){
    pout() << "TimeStepper::deallocateInternals" << endl;
  }
  
  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx]->deallocate_storage();
  }

  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
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
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();

    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);
    EBAMRCellData& cache = storage->get_cache();

    DataOps::copy(cache, solver->getPhi());
  }

  {// Cache Poisson solution
    MFAMRCellData& cache = m_fieldSolver_scratch->get_cache();
    DataOps::copy(cache, m_fieldSolver->getPotential());
  }

  // Cache RTE solutions
  for (RtIterator solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<RtSolver>& solver = solver_it();

    RefCountedPtr<RtStorage>& storage = this->get_RtStorage(solver_it);
    EBAMRCellData& cache = storage->get_cache();

    DataOps::copy(cache, solver->getPhi());
  }

  { // Cache sigma
    EBAMRIVData& cache = m_sigma_scratch->get_cache();
    DataOps::copy(cache, m_sigma->getPhi());
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
  this->computeCdrEbStates();              // Compute extrapolation n and grad(n) on the EB
  this->compute_cdr_fluxes(t0);               // Compute EB fluxes
  this->computeSigmaFlux_into_scratch();    // Compute sum of EB fluxes

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
  this->computeCdrVelo(t1);
  this->computeCdrEbStates();
  this->compute_cdr_fluxes(t1);
  this->computeSigmaFlux_into_scratch();

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
  
  EBAMRCellData& E_cell = m_fieldSolver_scratch->getElectricFieldCell();
  EBAMRFluxData& E_face = m_fieldSolver_scratch->getElectricFieldFace();
  EBAMRIVData&   E_eb   = m_fieldSolver_scratch->getElectricFieldEb();

  const MFAMRCellData& phi = m_fieldSolver->getPotential();
  
  this->compute_E(E_cell, m_cdr->getPhase(), phi);     // Compute cell-centered field
  this->compute_E(E_face, m_cdr->getPhase(), E_cell);  // Compute face-centered field
  this->compute_E(E_eb,   m_cdr->getPhase(), E_cell);  // EB-centered field
}

void rk2_tga::computeCdrVelo(const Real a_time){
  CH_TIME("rk2_tga::computeCdrVelo");
  if(m_verbosity > 5){
    pout() << "splitstep_::computeCdrVelo" << endl;
  }

  Vector<EBAMRCellData*> states     = m_cdr->getPhis();
  Vector<EBAMRCellData*> velocities = m_cdr->getVelocities();
  this->computeCdrDriftVelocities(velocities, states, m_fieldSolver_scratch->getElectricFieldCell(), a_time);
}

void rk2_tga::computeCdrEbStates(){
  CH_TIME("rk2_tga::computeCdrEbStates");
  if(m_verbosity > 5){
    pout() << "rk2_tga::computeCdrEbStates" << endl;
  }

  Vector<EBAMRIVData*>   eb_gradients;
  Vector<EBAMRIVData*>   eb_states;
  Vector<EBAMRCellData*> cdr_states;
  
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);

    cdr_states.push_back(&(solver->getPhi()));
    eb_states.push_back(&(storage->getEbState()));
    eb_gradients.push_back(&(storage->getEbGrad()));
  }

  this->extrapolate_to_eb(eb_states,          m_cdr->getPhase(), cdr_states);
  this->computeGradients_at_eb(eb_gradients, m_cdr->getPhase(), cdr_states);
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

  cdr_fluxes = m_cdr->getEbFlux();

  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);

    EBAMRIVData& dens_eb = storage->getEbState();
    EBAMRIVData& velo_eb = storage->getEbVelo();
    EBAMRIVData& flux_eb = storage->getEbFlux();
    EBAMRIVData& grad_eb = storage->getEbGrad();

    extrap_cdr_densities.push_back(&dens_eb);  // Computed in computeCdrEbStates
    extrap_cdr_velocities.push_back(&velo_eb);
    extrap_cdr_fluxes.push_back(&flux_eb);
    extrap_cdr_gradients.push_back(&grad_eb);  // Computed in computeCdrEbStates
  }

  // Extrapolate densities, velocities, and fluxes
  Vector<EBAMRCellData*> cdr_densities  = m_cdr->getPhis();
  Vector<EBAMRCellData*> cdr_velocities = m_cdr->getVelocities();
  this->compute_extrapolated_fluxes(extrap_cdr_fluxes, cdr_densities, cdr_velocities, m_cdr->getPhase());
  this->extrapolate_to_eb(extrap_cdr_velocities, m_cdr->getPhase(), cdr_velocities);

  // Compute RTE flux on the boundary
  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver   = solver_it();
    RefCountedPtr<RtStorage>& storage = this->get_RtStorage(solver_it);

    EBAMRIVData& flux_eb = storage->getEbFlux();
    solver->computeBoundaryFlux(flux_eb, solver->getPhi());
    extrap_rte_fluxes.push_back(&flux_eb);
  }

  const EBAMRIVData& E = m_fieldSolver_scratch->getElectricFieldEb();

  TimeStepper::compute_cdr_fluxes(cdr_fluxes,
				   extrap_cdr_fluxes,
				   extrap_cdr_densities,
				   extrap_cdr_velocities,
				   extrap_cdr_gradients,
				   extrap_rte_fluxes,
				   E,
				   a_time);
}

void rk2_tga::computeSigmaFlux_into_scratch(){
  CH_TIME("rk2_tga::computeSigmaFlux_into_scratch");
  if(m_verbosity > 5){
    pout() << "rk2_tga::computeSigmaFlux_into_scratch" << endl;
  }

  EBAMRIVData& flux = m_sigma->getFlux();
  DataOps::setValue(flux, 0.0);

  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const RefCountedPtr<species>& spec      = solver_it.getSpecies();
    const EBAMRIVData& solver_flux          = solver->getEbFlux();

    DataOps::incr(flux, solver_flux, spec->getChargeNumber()*Units::Qe);
  }

  m_sigma->resetCells(flux);
}

void rk2_tga::advance_advection_source_cdr_k1(const Real a_dt){
  CH_TIME("rk2_tga::advance_advection_source_cdr_k1");
  if(m_verbosity > 5){
    pout() << "rk2_tga::advance_advection_source_cdr_k1" << endl;
  }

  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);

    EBAMRCellData& k1  = storage->get_k1();
    EBAMRCellData& phi = solver->getPhi();
    EBAMRCellData& src = solver->getSource();

    // Compute rhs
    DataOps::setValue(k1, 0.0);
    solver->computeDivF(k1, phi, 0.0, true);
    DataOps::scale(k1, -1.0); 
    DataOps::incr(k1, src, 1.0);

    // Make phi = phi + k1*alpha*dt
    DataOps::incr(phi, k1, m_alpha*a_dt);

    m_amr->averageDown(phi, m_cdr->getPhase());
    m_amr->interpGhost(phi, m_cdr->getPhase());

    DataOps::floor(phi, 0.0);
  }
}

void rk2_tga::advance_advection_sigma_k1(const Real a_dt){
  CH_TIME("rk2_tga::advance_advection_sigma_k1");
  if(m_verbosity > 5){
    pout() << "rk2_tga::advance_advection_sigma_k1" << endl;
  }

  EBAMRIVData& phi = m_sigma->getPhi();
  EBAMRIVData& k1  = m_sigma_scratch->get_k1();
  
  m_sigma->computeRHS(k1);

  // Make phi = phi + k1*alpha*dt
  DataOps::incr(phi, k1,    m_alpha*a_dt);

  m_amr->averageDown(phi, m_cdr->getPhase());
  
  m_sigma->resetCells(k1);
  m_sigma->resetCells(phi);
}

void rk2_tga::advance_advection_source_cdr_k2(const Real a_dt){
  CH_TIME("rk2_tga::advance_advection_source_cdr_k2");
  if(m_verbosity > 5){
    pout() << "rk2_tga::advance_advection_source_cdr_k2" << endl;
  }

  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);

    EBAMRCellData& state       = solver->getPhi();
    EBAMRCellData& k2          = storage->get_k2();
    const EBAMRCellData& k1    = storage->get_k1();
    EBAMRCellData& src         = solver->getSource();

    // Compute RHS
    DataOps::setValue(k2, 0.0);
    solver->computeDivF(k2, state, 0.0, true);
    DataOps::scale(k2, -1.0);
    DataOps::incr(k2, src, 1.0);

    // RK2 advance. The extract m_alpha subtraction is because when we came here, the solver state (which we update in place)
    // contained the intermediate state phi + k1*alpha_dt. But we want phi = phi + k1*a_dt*(1-1/(2*alpha)) + k2*dt/(2*alpha),
    // so we just adjust the factor directly.
    const Real k1_factor = a_dt*(1.0 - 1.0/(2.0*m_alpha) - m_alpha);
    const Real k2_factor = a_dt/(2.0*m_alpha);
    
    DataOps::incr(state, k1, k1_factor);
    DataOps::incr(state, k2, k2_factor);

    m_amr->averageDown(state, m_cdr->getPhase());
    m_amr->interpGhost(state, m_cdr->getPhase());

    DataOps::floor(state, 0.0);
  }
}

void rk2_tga::advance_advection_sigma_k2(const Real a_dt){
  CH_TIME("rk2_tga::advance_advection_sigma_k2");
  if(m_verbosity > 5){
    pout() << "rk2_tga::advance_advection_sigma_k2" << endl;
  }
  
  const EBAMRIVData& k1  = m_sigma_scratch->get_k1();
  EBAMRIVData& k2        = m_sigma_scratch->get_k2();
  m_sigma->computeRHS(k2);

  EBAMRIVData& state       = m_sigma->getPhi();

  // RK2 advance. The extract m_alpha subtraction is because when we came here, the solver state (which we update in place)
  // contained the intermediate state phi + k1*alpha_dt. But we want phi = phi + k1*a_dt*(1-1/(2*alpha)) + k2*dt/(2*alpha),
  // so we just adjust the factor directly.
  const Real k1_factor = a_dt*(1.0 - 1.0/(2.0*m_alpha) - m_alpha);
  const Real k2_factor = a_dt/(2.0*m_alpha);
    
  DataOps::incr(state, k1, k1_factor);
  DataOps::incr(state, k2, k2_factor);

  m_amr->averageDown(state, m_cdr->getPhase());
  m_sigma->resetCells(state);
}

void rk2_tga::advance_diffusion(const Real a_dt){
  CH_TIME("rk2_tga::advance_diffusion");
  if(m_verbosity > 5){
    pout() << "rk2_tga::advance_diffusion" << endl;
  }

  // Diffusive advance for all cdr equations
  bool diffusive_states = false;
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();

    if(solver->isDiffusive()){
      diffusive_states = true;
    }
  }

  // Do the diffusion advance
  if(diffusive_states){
    m_cdr->setSource(0.0); // This is necessary because advance_diffusion also works with source terms
    this->compute_cdr_diffusion();
    for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<CdrSolver>& solver = solver_it();

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
  
  Vector<EBAMRCellData*> cdr_sources = m_cdr->getSources();
  Vector<EBAMRCellData*> cdr_states  = m_cdr->getPhis();
  Vector<EBAMRCellData*> rte_states  = m_rte->getPhis();
  EBAMRCellData& E                   = m_fieldSolver_scratch->getElectricFieldCell();

  this->compute_cdr_sources(cdr_sources, cdr_states, rte_states, E, a_time, centering::cell_center);
}

void rk2_tga::advance_rte_stationary(const Real a_time){
  CH_TIME("rk2_tga::compute_rte_k1_stationary");
  if(m_verbosity > 5){
    pout() << "rk2_tga::compute_k1_stationary" << endl;
  }

  if((m_timeStep + 1) % m_fast_rte == 0){
    Vector<EBAMRCellData*> rte_states  = m_rte->getPhis();
    Vector<EBAMRCellData*> rte_sources = m_rte->getSources();
    Vector<EBAMRCellData*> cdr_states  = m_cdr->getPhis();

    EBAMRCellData& E = m_fieldSolver_scratch->getElectricFieldCell();

    const Real dummy_dt = 0.0;
    this->solve_rte(rte_states, rte_sources, cdr_states, E, a_time, dummy_dt, centering::cell_center);
  }
#include "CD_NamespaceFooter.H"
