/*!
  @file rk2.cpp
  @brief Implementation of rk2.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "rk2.H"
#include "rk2_storage.H"
#include <CD_CdrIterator.H>
#include <CD_RtIterator.H>
#include <CD_DataOps.H>
#include <CD_Units.H>

#include <ParmParse.H>

#define RK2_DEBUG_TIMER 0
#define RK2_DEBUG_TIMER_STOP 0

typedef rk2::CdrStorage     CdrStorage;
typedef rk2::FieldStorage FieldStorage;
typedef rk2::RtStorage     RtStorage;
typedef rk2::SigmaStorage   SigmaStorage;

rk2::rk2() : TimeStepper() {
  m_alpha = 1.0;

  {
    ParmParse pp("rk2");
    pp.query("alpha", m_alpha);
  }
}

rk2::~rk2(){
  
}

RefCountedPtr<CdrStorage>& rk2::get_CdrStorage(const CdrIterator& a_solverit){
  return m_cdr_scratch[a_solverit.get_solver()];
}

RefCountedPtr<RtStorage>& rk2::get_RtStorage(const RtIterator& a_solverit){
  return m_rte_scratch[a_solverit.get_solver()];
}

void rk2::allocateCdrStorage(){

  const int ncomp       = 1;
  const int num_species = m_plaskin->get_num_species();
  m_cdr_scratch.resize(num_species);
  
  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx] = RefCountedPtr<CdrStorage> (new CdrStorage(m_amr, m_cdr->getPhase(), ncomp));
    m_cdr_scratch[idx]->allocate_storage();
  }
}

void rk2::allocateFieldStorage(){
  const int ncomp = 1;
  m_fieldSolver_scratch = RefCountedPtr<FieldStorage> (new FieldStorage(m_amr, m_cdr->getPhase(), ncomp));
  m_fieldSolver_scratch->allocate_storage();
}

void rk2::allocateRtStorage(){
  const int ncomp       = 1;
  const int num_Photons = m_plaskin->get_num_Photons();
  m_rte_scratch.resize(num_Photons);
  
  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx] = RefCountedPtr<RtStorage> (new RtStorage(m_amr, m_rte->getPhase(), ncomp));
    m_rte_scratch[idx]->allocate_storage();
  }
}

void rk2::allocateSigmaStorage(){
  const int ncomp = 1;
  m_sigma_scratch = RefCountedPtr<SigmaStorage> (new SigmaStorage(m_amr, m_cdr->getPhase(), ncomp));
  m_sigma_scratch->allocate_storage();
}

void rk2::deallocateInternals(){
  CH_TIME("rk2::deallocateInternals");
  if(m_verbosity > 5){
    pout() << "rk2::deallocateInternals" << endl;
  }

  m_fieldSolver_scratch->deallocate_storage();
  m_sigma_scratch->deallocate_storage();

  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_cdr_scratch[idx]->deallocate_storage();
  }

  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    m_rte_scratch[idx]->deallocate_storage();
  }
}

void rk2::regridInternals(){
  CH_TIME("rk2::regridInternals");
  if(m_verbosity > 5){
    pout() << "rk2::regridInternals" << endl;
  }
  
  this->allocateCdrStorage();
  this->allocateFieldStorage();
  this->allocateRtStorage();
  this->allocateSigmaStorage();
}

Real rk2::advance(const Real a_dt){
  CH_TIME("rk2::advance");
  if(m_verbosity > 2){
    pout() << "rk2::advance" << endl;
  }

  // Prepare for k1 advance
  const Real t0 = MPI_Wtime();
  this->compute_E_at_start_of_time_step();
  const Real t00 = MPI_Wtime();
  this->computeCdrVelo_at_start_of_time_step();
  const Real t01 = MPI_Wtime();
  this->computeCdrEbStates_at_start_of_time_step();
  const Real t02 = MPI_Wtime();
  this->compute_cdr_diffco_at_start_of_time_step();
  const Real t03 = MPI_Wtime();
  this->compute_cdr_sources_at_start_of_time_step();
  const Real t04 = MPI_Wtime();
  this->compute_cdr_fluxes_at_start_of_time_step();
  const Real t05 = MPI_Wtime();
  this->computeSigmaFlux_at_start_of_time_step();

  // Do k1 advance
  const Real t1 = MPI_Wtime();
  this->advance_cdr_k1(a_dt);
  const Real t10 = MPI_Wtime();
  this->advance_sigma_k1(a_dt);
  const Real t11 = MPI_Wtime();
  this->solve_poisson_k1();
  const Real t12 = MPI_Wtime();
  this->compute_E_after_k1();
  const Real t13 = MPI_Wtime();
  if(m_rte->isStationary()){
    this->advance_rte_k1_stationary(a_dt);
  }
  else{
    this->advance_rte_k1_transient(a_dt);
  }
  const Real t14 = MPI_Wtime();

  // Recompute things in order to do k2 advance
  const Real t2 = MPI_Wtime();
  this->computeCdrEbStates_after_k1();
  const Real t20 = MPI_Wtime();
  this->computeCdrVelo_after_k1(a_dt);
  const Real t21 = MPI_Wtime();
  this->compute_cdr_diffco_after_k1(a_dt);
  const Real t22 = MPI_Wtime();
  this->compute_cdr_sources_after_k1(a_dt);
  const Real t23 = MPI_Wtime();
  this->compute_cdr_fluxes_after_k1(a_dt);
  const Real t24 = MPI_Wtime();
  this->computeSigmaFlux_after_k1();
  const Real t25 = MPI_Wtime();
  
  // Do k2 advance
  const Real t3 = MPI_Wtime();
  const Real t30 = MPI_Wtime();
  this->advance_cdr_k2(a_dt);
  const Real t31 = MPI_Wtime();
  this->advance_sigma_k2(a_dt);
  const Real t32 = MPI_Wtime();
  this->solve_poisson_k2();
  const Real t33 = MPI_Wtime();
  this->compute_E_after_k2();
  const Real t34 = MPI_Wtime();
  if(m_rte->isStationary()){
    this->advance_rte_k2_stationary(a_dt);
  }
  else{
    this->advance_rte_k2_transient(a_dt);
  }
  const Real t4 = MPI_Wtime();

#if RK2_DEBUG_TIMER
  pout() << endl;
  pout() << "rk2::advance breakdown" << endl;

  pout() << "t1 - t0 % = " << 100.*(t1 - t0)/(t4-t0) << "%" << endl;
  pout() << "t2 - t1 % = " << 100.*(t2 - t1)/(t4-t0) << "%" << endl;
  pout() << "t3 - t2 % = " << 100.*(t3 - t2)/(t4-t0) << "%" << endl;
  pout() << "t4 - t3 % = " << 100.*(t4 - t3)/(t4-t0) << "%" << endl;
  pout() << "Total time = " << t4 - t0 << endl;
  pout() << endl;
  pout() << "t00 - t0  = " << 100.*(t00 - t0)/(t4-t0) << "%"  << endl;
  pout() << "t01 - t00 = " << 100.*(t01 - t00)/(t4-t0) << "%" << endl;
  pout() << "t02 - t01 = " << 100.*(t02 - t01)/(t4-t0) << "%" << endl;
  pout() << "t03 - t02 = " << 100.*(t03 - t02)/(t4-t0) << "%" << endl;
  pout() << "t04 - t03 = " << 100.*(t04 - t03)/(t4-t0) << "%" << endl;
  pout() << "t05 - t04 = " << 100.*(t05 - t04)/(t4-t0) << "%" << endl;
  pout() << "t1  - t05 = " << 100.*(t1  - t05)/(t4-t0) << "%" << endl;
  pout() << "Total = " << 100.*(t1-t0)/(t4-t0) << "%" << endl;
  pout() << endl;
  pout() << "t10 - t1  = " << 100.*(t10 - t1)/(t4-t0) << "%"  << endl;
  pout() << "t11 - t10 = " << 100.*(t11 - t10)/(t4-t0) << "%" << endl;
  pout() << "t12 - t11 = " << 100.*(t12 - t11)/(t4-t0) << "%" << endl;
  pout() << "t13 - t12 = " << 100.*(t13 - t12)/(t4-t0) << "%" << endl;
  pout() << "t14 - t13 = " << 100.*(t14 - t13)/(t4-t0) << "%" << endl;
  pout() << "Total = " << 100.*(t2-t1)/(t4-t0) << "%" << endl;
  pout() << endl;
  pout() << "t20 - t2  = " << 100.*(t20 - t2)/(t4-t0) << "%"  << endl;
  pout() << "t21 - t20 = " << 100.*(t21 - t20)/(t4-t0) << "%" << endl;
  pout() << "t22 - t21 = " << 100.*(t22 - t21)/(t4-t0) << "%" << endl;
  pout() << "t23 - t22 = " << 100.*(t23 - t22)/(t4-t0) << "%" << endl;
  pout() << "t24 - t23 = " << 100.*(t24 - t23)/(t4-t0) << "%" << endl;
  pout() << "t25 - t24 = " << 100.*(t25 - t24)/(t4-t0) << "%" << endl;
  pout() << "Total = " << 100.*(t3-t2)/(t4-t0) << "%" << endl;
  pout() << endl;
  pout() << "t31 - t30 = " << 100.*(t31 - t30)/(t4-t0) << "%" << endl;
  pout() << "t32 - t31 = " << 100.*(t32 - t31)/(t4-t0) << "%" << endl;
  pout() << "t33 - t32 = " << 100.*(t33 - t32)/(t4-t0) << "%" << endl;
  pout() << "t34 - t33 = " << 100.*(t34 - t33)/(t4-t0) << "%" << endl;
  pout() << "t4  - t34 = " << 100.*(t4  - t34)/(t4-t0) << "%" << endl;
  pout() << "Total = " << 100.*(t4-t3)/(t4-t0) << "%" << endl;
  pout() << endl;
#if RK2_DEBUG_TIMER_STOP
  MayDay::Abort("rk2::advance - debug timer stop");
#endif
#endif

  return a_dt;
}

void rk2::compute_E_at_start_of_time_step(){
  CH_TIME("rk2::compute_E_at_start_of_time_step");
  if(m_verbosity > 5){
    pout() << "rk2::compute_E_at_start_of_time_step" << endl;
  }

  EBAMRCellData& E_cell = m_fieldSolver_scratch->getElectricFieldCell();
  EBAMRFluxData& E_face = m_fieldSolver_scratch->getElectricFieldFace();
  EBAMRIVData&   E_eb   = m_fieldSolver_scratch->getElectricFieldEb();

  const MFAMRCellData& phi = m_fieldSolver->getPotential();
  
  this->compute_E(E_cell, m_cdr->getPhase(), phi);     // Compute cell-centered field
  this->compute_E(E_face, m_cdr->getPhase(), E_cell);  // Compute face-centered field
  this->compute_E(E_eb,   m_cdr->getPhase(), E_cell);  // EB-centered field
}

void rk2::computeCdrVelo_at_start_of_time_step(){
  CH_TIME("rk2::computeCdrVelo_at_start_of_time_step");
  if(m_verbosity > 5){
    pout() << "rk2::computeCdrVelo_at_start_of_time_step" << endl;
  }

  Vector<EBAMRCellData*> states     = m_cdr->getPhis();
  Vector<EBAMRCellData*> velocities = m_cdr->getVelocities();
  this->computeCdrDriftVelocities(velocities, states, m_fieldSolver_scratch->getElectricFieldCell(), m_time);
}

void rk2::computeCdrEbStates_at_start_of_time_step(){
  CH_TIME("rk2::computeCdrEbStates_at_start_of_time_step");
  if(m_verbosity > 5){
    pout() << "rk2::computeCdrEbStates_at_start_of_time_step" << endl;
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

void rk2::compute_cdr_diffco_at_start_of_time_step(){
  CH_TIME("rk2::compute_cdr_diffco_at_start_of_time_step");
  if(m_verbosity > 5){
    pout() << "rk2::compute_cdr_diffco_at_start_of_time_step" << endl;
  }

  const int num_species = m_plaskin->get_num_species();

  Vector<EBAMRCellData*> cdr_states  = m_cdr->getPhis();
  Vector<EBAMRFluxData*> diffco_face = m_cdr->getFaceCenteredDiffusionCoefficient();
  Vector<EBAMRIVData*> diffco_eb     = m_cdr->getEbCenteredDiffusionCoefficient();

  const EBAMRCellData& E_cell = m_fieldSolver_scratch->getElectricFieldCell();
  const EBAMRIVData& E_eb     = m_fieldSolver_scratch->getElectricFieldEb();

  // Get extrapolated states
  Vector<EBAMRIVData*> eb_states(num_species);
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.get_solver();
    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);
    eb_states[idx] = &(storage->getEbState());
  }
  
  this->compute_cdr_diffco_face(diffco_face, cdr_states, E_cell, m_time);
  this->compute_cdr_diffco_eb(diffco_eb,     eb_states,  E_eb,   m_time);
}

void rk2::compute_cdr_sources_at_start_of_time_step(){
  CH_TIME("rk2::compute_cdr_sources_at_start_of_time_step");
  if(m_verbosity > 5){
    pout() << "rk2::compute_cdr_sources_at_start_of_time_step" << endl;
  }
  
  Vector<EBAMRCellData*> cdr_sources = m_cdr->getSources();
  Vector<EBAMRCellData*> cdr_states  = m_cdr->getPhis();
  Vector<EBAMRCellData*> rte_states  = m_rte->getPhis();
  EBAMRCellData& E                   = m_fieldSolver_scratch->getElectricFieldCell();

  this->compute_cdr_sources(cdr_sources, cdr_states, rte_states, E, m_time, centering::cell_center);
}

void rk2::compute_cdr_fluxes_at_start_of_time_step(){
  CH_TIME("rk2::compute_cdr_fluxes_at_start_of_time_step");
  if(m_verbosity > 5){
    pout() << "rk2::compute_cdr_fluxes_at_start_of_time_step" << endl;
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

    extrap_cdr_densities.push_back(&dens_eb);  // Already been computed
    extrap_cdr_velocities.push_back(&velo_eb);
    extrap_cdr_fluxes.push_back(&flux_eb);
    extrap_cdr_gradients.push_back(&grad_eb);  // Already been computed
  }
  
  // Extrapolate densities, velocities, and fluxes
  Vector<EBAMRCellData*> cdr_densities = m_cdr->getPhis();
  Vector<EBAMRCellData*> cdr_velocities = m_cdr->getVelocities();
  this->compute_extrapolated_fluxes(extrap_cdr_fluxes, cdr_densities, cdr_velocities, m_cdr->getPhase());
  this->extrapolate_to_eb(extrap_cdr_velocities, m_cdr->getPhase(), cdr_velocities);
  //  this->extrapolate_to_eb(extrap_cdr_densities,  m_cdr->getPhase(), cdr_densities); // Already been done, no?

  // Compute RTE flux on the boundary
  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver   = solver_it();
    RefCountedPtr<RtStorage>& storage = this->get_RtStorage(solver_it);

    EBAMRIVData& flux_eb = storage->getEbFlux();
    solver->computeBoundaryFlux(flux_eb, solver->getPhi());
    extrap_rte_fluxes.push_back(&flux_eb);
  }

  const EBAMRIVData& E = m_fieldSolver_scratch->getElectricFieldEb();

  this->compute_cdr_fluxes(cdr_fluxes,
			   extrap_cdr_fluxes,
			   extrap_cdr_densities,
			   extrap_cdr_velocities,
			   extrap_cdr_gradients,
			   extrap_rte_fluxes,
			   E,
			   m_time);
}

void rk2::computeSigmaFlux_at_start_of_time_step(){
  CH_TIME("rk2::computeSigmaFlux_at_start_of_time_step");
  if(m_verbosity > 5){
    pout() << "rk2::computeSigmaFlux_at_start_of_time_step" << endl;
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

void rk2::advance_cdr_k1(const Real a_dt){
  CH_TIME("rk2::advance_cdr_k1");
  if(m_verbosity > 5){
    pout() << "rk2::advance_cdr_k1" << endl;
  }

  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);

    EBAMRCellData& k1  = storage->get_k1();
    EBAMRCellData& phi = storage->getPhi();

    const EBAMRCellData& state = solver->getPhi();

    DataOps::setValue(k1, 0.0);
    solver->computeRHS(k1, state, a_dt);

    DataOps::setValue(phi, 0.0);
    DataOps::incr(phi, state, 1.0);
    DataOps::incr(phi, k1,    m_alpha*a_dt);

    m_amr->averageDown(phi, m_cdr->getPhase());
    m_amr->interpGhost(phi, m_cdr->getPhase());

    DataOps::floor(phi, 0.0);
  }
}

void rk2::advance_sigma_k1(const Real a_dt){
  CH_TIME("rk2::advance_sigma_k1");
  if(m_verbosity > 5){
    pout() << "rk2::advance_sigma_k1" << endl;
  }

  const EBAMRIVData& state = m_sigma->getPhi();
  
  EBAMRIVData& k1  = m_sigma_scratch->get_k1();
  EBAMRIVData& phi = m_sigma_scratch->getPhi();
  m_sigma->computeRHS(k1);
  DataOps::setValue(phi,   0.0);
  DataOps::incr(phi, state, 1.0);
  DataOps::incr(phi, k1,    m_alpha*a_dt);

  m_amr->averageDown(phi, m_cdr->getPhase());
  
  m_sigma->resetCells(k1);
  m_sigma->resetCells(phi);
}

void rk2::solve_poisson_k1(){
  CH_TIME("rk2::solve_poisson_k1");
  if(m_verbosity > 5){
    pout() << "rk2::solve_poisson_k1" << endl;
  }

  MFAMRCellData& scratch_pot = m_fieldSolver_scratch->getPhi();
  EBAMRIVData& sigma         = m_sigma_scratch->getPhi();
  Vector<EBAMRCellData*> cdr_densities;
  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);
    cdr_densities.push_back(&(storage->getPhi()));
  }

  DataOps::setValue(scratch_pot, 0.0);
  DataOps::incr(scratch_pot, m_fieldSolver->getPotential(), 1.0);

  if((m_timeStep + 1) % m_fast_poisson == 0){
    const bool converged = this->solve_poisson(scratch_pot, m_fieldSolver->getRho(), cdr_densities, sigma, centering::cell_center);
    if(!converged){
      pout() << "rk2::solve_poisson_k1 - solver did not converge at step = " << m_timeStep << endl;
    }
  }
}

void rk2::compute_E_after_k1(){
  CH_TIME("rk2::compute_E_after_k1");
  if(m_verbosity > 5){
    pout() << "rk2::compute_E_after_k1" << endl;
  }
  
  EBAMRCellData& E_cell = m_fieldSolver_scratch->getElectricFieldCell();
  EBAMRFluxData& E_face = m_fieldSolver_scratch->getElectricFieldFace();
  EBAMRIVData&   E_eb   = m_fieldSolver_scratch->getElectricFieldEb();

  const MFAMRCellData& phi = m_fieldSolver_scratch->getPhi();
  
  this->compute_E(E_cell, m_cdr->getPhase(), phi);     // Compute cell-centered field
  this->compute_E(E_face, m_cdr->getPhase(), E_cell);  // Compute face-centered field
  this->compute_E(E_eb,   m_cdr->getPhase(), E_cell);  // EB-centered field
}

void rk2::advance_rte_k1_stationary(const Real a_dt){
  CH_TIME("rk2::compute_rte_k1_stationary");
  if(m_verbosity > 5){
    pout() << "rk2::compute_k1_stationary" << endl;
  }

  const Real time = m_time + m_alpha*a_dt;

  Vector<EBAMRCellData*> rte_states;
  Vector<EBAMRCellData*> rte_sources;
  Vector<EBAMRCellData*> cdr_states;

  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver   = solver_it();
    RefCountedPtr<RtStorage>& storage = this->get_RtStorage(solver_it);
    
    EBAMRCellData& phi    = storage->getPhi();
    EBAMRCellData& state  = solver->getPhi();
    EBAMRCellData& source = solver->getSource();

    DataOps::setValue(phi, 0.0);
    DataOps::incr(phi, state, 1.0);

    rte_states.push_back(&(phi));
    rte_sources.push_back(&(source));
  }

  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);
    cdr_states.push_back(&(storage->getPhi()));
  }

  EBAMRCellData& E = m_fieldSolver_scratch->getElectricFieldCell();

  if((m_timeStep + 1) % m_fast_rte == 0){
    const Real dummy_dt = 0.0;
    this->solve_rte(rte_states, rte_sources, cdr_states, E, time, dummy_dt, centering::cell_center);
  }
}

void rk2::advance_rte_k1_transient(const Real a_dt){
  CH_TIME("rk2::compute_rte_k1_transient");
  if(m_verbosity > 5){
    pout() << "rk2::compute_k1_transient" << endl;
  }

  const Real time = m_time + 0.5*m_alpha*a_dt; // Source terms are now centered on the half time step

  Vector<EBAMRCellData*> rte_states;
  Vector<EBAMRCellData*> rte_sources;
  Vector<EBAMRCellData*> cdr_states;

  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver   = solver_it();
    RefCountedPtr<RtStorage>& storage = this->get_RtStorage(solver_it);
    
    EBAMRCellData& phi    = storage->getPhi();
    EBAMRCellData& state  = solver->getPhi();
    EBAMRCellData& source = solver->getSource();

    DataOps::setValue(phi, 0.0);
    DataOps::incr(phi, state, 1.0);
    
    rte_states.push_back(&(phi));
    rte_sources.push_back(&(source));
  }

  // Source term must be centered between time tn and the intermediate time
  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);
    
    const EBAMRCellData& state = solver->getPhi();
    const EBAMRCellData& phi   = storage->getPhi();

    EBAMRCellData& scratch = storage->getScratch();

    DataOps::setValue(scratch, 0.0);
    DataOps::incr(scratch, state, 0.5);
    DataOps::incr(scratch, phi,   0.5);

    cdr_states.push_back(&(scratch));
  }


  if((m_timeStep + 1) % m_fast_rte == 0){ // Actual solve
    const MFAMRCellData& state = m_fieldSolver->getPotential();
    const MFAMRCellData& phi   = m_fieldSolver_scratch->getPhi();
  
    MFAMRCellData& scratch_phi = m_fieldSolver_scratch->getScratch_phi();
    EBAMRCellData& scratch_E   = m_fieldSolver_scratch->getScratch_E();

    DataOps::setValue(scratch_phi, 0.0);
    DataOps::incr(scratch_phi, state, 0.5);
    DataOps::incr(scratch_phi, phi, 0.5);

    m_amr->averageDown(scratch_phi);
    m_amr->interpGhost(scratch_phi);

    this->compute_E(scratch_E, m_cdr->getPhase(), scratch_phi);
  

    this->solve_rte(rte_states, rte_sources, cdr_states, scratch_E, time, m_alpha*a_dt, centering::cell_center);
  }
}

void rk2::computeCdrEbStates_after_k1(){
  CH_TIME("rk2::computeCdrEbStates_after_k1");
  if(m_verbosity > 5){
    pout() << "rk2::computeCdrEbStates_after_k1" << endl;
  }

  Vector<EBAMRIVData*>   eb_gradients;
  Vector<EBAMRIVData*>   eb_states;
  Vector<EBAMRCellData*> cdr_states;
  
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);
    cdr_states.push_back(&(storage->getPhi()));
    eb_states.push_back(&(storage->getEbState()));
    eb_gradients.push_back(&(storage->getEbGrad()));
  }

  this->extrapolate_to_eb(eb_states,          m_cdr->getPhase(), cdr_states);
  this->computeGradients_at_eb(eb_gradients, m_cdr->getPhase(), cdr_states);
}

void rk2::computeCdrVelo_after_k1(const Real a_dt){
  CH_TIME("rk2::computeCdrVelo_after_k1");
  if(m_verbosity > 5){
    pout() << "rk2::computeCdrVelo_after_k1" << endl;
  }
  
  const int num_species = m_plaskin->get_num_species();
  
  const Real time = m_time + m_alpha*a_dt;
  
  Vector<EBAMRCellData*> states(num_species);
  Vector<EBAMRCellData*> velocities = m_cdr->getVelocities();

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);
    const int idx = solver_it.get_solver();
    states[idx] = &(storage->getPhi());
  }
  
  this->computeCdrDriftVelocities(velocities, states, m_fieldSolver_scratch->getElectricFieldCell(), time);
}

void rk2::compute_cdr_diffco_after_k1(const Real a_dt){
  CH_TIME("rk2::compute_cdr_diffco_after_k1");
  if(m_verbosity > 5){
    pout() << "rk2::compute_cdr_diffco_after_k1" << endl;
  }

  const int num_species = m_plaskin->get_num_species();

  const Real time = m_time + m_alpha*a_dt;
  
  Vector<EBAMRFluxData*> diffco_face = m_cdr->getFaceCenteredDiffusionCoefficient();
  Vector<EBAMRIVData*> diffco_eb     = m_cdr->getEbCenteredDiffusionCoefficient();

  const EBAMRCellData& E_cell = m_fieldSolver_scratch->getElectricFieldCell();
  const EBAMRIVData& E_eb     = m_fieldSolver_scratch->getElectricFieldEb();

  Vector<EBAMRCellData*> states(num_species);
  Vector<EBAMRIVData*> eb_states(num_species);
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);
    const int idx  = solver_it.get_solver();
    states[idx]    = &(storage->getPhi());
    eb_states[idx] = &(storage->getEbState());
  }
  
  this->compute_cdr_diffco_face(diffco_face, states,    E_cell, time);
  this->compute_cdr_diffco_eb(diffco_eb,     eb_states, E_eb,   time);
}

void rk2::compute_cdr_sources_after_k1(const Real a_dt){
  CH_TIME("rk2::compute_cdr_sources_after_k1");
  if(m_verbosity > 5){
    pout() << "rk2::compute_cdr_sources_after_k1" << endl;
  }

  const Real time = m_time + m_alpha*a_dt;

  Vector<EBAMRCellData*> cdr_sources;
  Vector<EBAMRCellData*> cdr_states;
  Vector<EBAMRCellData*> rte_states;

  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtStorage>& storage = this->get_RtStorage(solver_it);
    rte_states.push_back(&(storage->getPhi()));
  }

  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);
    
    cdr_states.push_back(&(storage->getPhi()));
    cdr_sources.push_back(&(solver->getSource()));
  }

  EBAMRCellData& E = m_fieldSolver_scratch->getElectricFieldCell();

  this->compute_cdr_sources(cdr_sources, cdr_states, rte_states, E, time, centering::cell_center);
}

void rk2::compute_cdr_fluxes_after_k1(const Real a_dt){
  CH_TIME("rk2::compute_cdr_fluxes_after_k1");
  if(m_verbosity > 5){
    pout() << "rk2::compute_cdr_fluxes_after_k1" << endl;
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

  cdr_fluxes = m_cdr->getEbFlux();

  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);

    EBAMRCellData& dens  = storage->getPhi();
    EBAMRIVData& dens_eb = storage->getEbState();
    EBAMRIVData& velo_eb = storage->getEbVelo();
    EBAMRIVData& flux_eb = storage->getEbFlux();
    EBAMRIVData& grad_eb = storage->getEbGrad();

    cdr_densities.push_back(&dens);
    extrap_cdr_densities.push_back(&dens_eb);  // This has already been extrapolated to the EB
    extrap_cdr_velocities.push_back(&velo_eb); // This has not.
    extrap_cdr_fluxes.push_back(&flux_eb);     // This hasn't either. 
    extrap_cdr_gradients.push_back(&grad_eb);  // This has already been extrapolated to the EB
  }

  // Extrapolate the flux and the velocity
  Vector<EBAMRCellData*> cdr_velocities = m_cdr->getVelocities();
  this->compute_extrapolated_fluxes(extrap_cdr_fluxes, cdr_densities, cdr_velocities, m_cdr->getPhase());
  this->extrapolate_to_eb(extrap_cdr_velocities, m_cdr->getPhase(), cdr_velocities);
  //  this->extrapolate_to_eb(extrap_cdr_densities,  m_cdr->getPhase(), cdr_densities); // This has already been done, no?

  // Compute RTE flux on the boundary
  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver   = solver_it();
    RefCountedPtr<RtStorage>& storage = this->get_RtStorage(solver_it);

    EBAMRIVData& flux_eb = storage->getEbFlux();
    solver->computeBoundaryFlux(flux_eb, storage->getPhi());
    extrap_rte_fluxes.push_back(&flux_eb);
  }


  const EBAMRIVData& E = m_fieldSolver_scratch->getElectricFieldEb();

  this->compute_cdr_fluxes(cdr_fluxes,
			   extrap_cdr_fluxes,
			   extrap_cdr_densities,
			   extrap_cdr_velocities,
			   extrap_cdr_gradients, 
			   extrap_rte_fluxes,
			   E,
			   time);
}

void rk2::computeSigmaFlux_after_k1(){
  CH_TIME("rk2::computeSigmaFlux_after_k1");
  if(m_verbosity > 5){
    pout() << "rk2::computeSigmaFlux_after_k1" << endl;
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

void rk2::advance_cdr_k2(const Real a_dt){
  CH_TIME("rk2::advance_cdr_k2");
  if(m_verbosity > 5){
    pout() << "rk2::advance_cdr_k2" << endl;
  }

  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);

    EBAMRCellData& state   = solver->getPhi();
    EBAMRCellData& k1      = storage->get_k1();
    EBAMRCellData& k2      = storage->get_k2();
    EBAMRCellData& phi     = storage->getPhi();
    EBAMRCellData& scratch = storage->getScratch();

    solver->computeRHS(k2, phi, a_dt);

    // For transient RTE solvers I need the source term at half time steps. Since the internal state
    // inside the solver will be overwritten, I take a backup into rk2_storage.scratch
    if(!m_rte->isStationary()){
      DataOps::setValue(scratch,   0.0);
      DataOps::incr(scratch, state, 1.0);
    }

    
    DataOps::incr(state, k1, a_dt*(1 - 1./(2.*m_alpha)));
    DataOps::incr(state, k2, a_dt*1./(2.*m_alpha));

    m_amr->averageDown(state, m_cdr->getPhase());
    m_amr->interpGhost(state, m_cdr->getPhase());

    DataOps::floor(state, 0.0);
  }
}

void rk2::advance_sigma_k2(const Real a_dt){
  CH_TIME("rk2::advance_sigma_k2");
  if(m_verbosity > 5){
    pout() << "rk2::advance_sigma_k2" << endl;
  }
  
  EBAMRIVData& k1  = m_sigma_scratch->get_k1();
  EBAMRIVData& k2  = m_sigma_scratch->get_k2();
  m_sigma->computeRHS(k2);

  EBAMRIVData& state = m_sigma->getPhi();
  DataOps::incr(state, k1, a_dt*(1 - 1./(2.*m_alpha)));
  DataOps::incr(state, k2, a_dt*1./(2.*m_alpha));

  m_amr->averageDown(state, m_cdr->getPhase());
  m_sigma->resetCells(state);
}

void rk2::solve_poisson_k2(){
  CH_TIME("rk2::solve_poisson_k2");
  if(m_verbosity > 5){
    pout() << "rk2::solve_poisson_k2" << endl;
  }

  // We computed the intermediate potential at time t_k + alpha*dt. Linearly extrapolate that result to the end of
  // the time step. The result for this is y_extrap = y_alpha/alpha - y_0*(1-alpha)/alpha. This (usually) brings the
  // initial guess closer to the true solution.
  MFAMRCellData& pot     = m_fieldSolver->getPotential();
  MFAMRCellData& phi     = m_fieldSolver_scratch->getPhi();
  MFAMRCellData& scratch = m_fieldSolver_scratch->getScratch_phi();

  // For transient RTE solvers I need the source term at half time steps. Since the internal state
  // inside the solver will be overwritten, I take a backup into FieldStorage.scratch_phi
  if(!m_rte->isStationary()){
    DataOps::setValue(scratch, 0.0);
    DataOps::incr(scratch, pot, 1.0); 
  }

  if((m_timeStep + 1) % m_fast_poisson == 0){
    DataOps::scale(pot, -(1.0 - m_alpha)/m_alpha);
    DataOps::incr(pot, phi, 1./m_alpha);

    const bool converged = this->solve_poisson(pot,
					       m_fieldSolver->getRho(),
					       m_cdr->getPhis(),
					       m_sigma->getPhi(),
					       centering::cell_center);

    if(!converged){
      pout() << "rk2::solve_poisson_k2 - solver did not converge at step = " << m_timeStep << endl;
    }
  }
  else{
    DataOps::setValue(pot, 0.0);
    DataOps::incr(pot, phi, 1.0);
  }
}

void rk2::compute_E_after_k2(){
  CH_TIME("rk2::compute_E_at_after_k2");
  if(m_verbosity > 5){
    pout() << "rk2::compute_E_after_k2" << endl;
  }
  
  EBAMRCellData& E_cell = m_fieldSolver_scratch->getElectricFieldCell();
  EBAMRFluxData& E_face = m_fieldSolver_scratch->getElectricFieldFace();
  EBAMRIVData&   E_eb   = m_fieldSolver_scratch->getElectricFieldEb();

  const MFAMRCellData& phi = m_fieldSolver->getPotential();
  
  this->compute_E(E_cell, m_cdr->getPhase(), phi);     // Compute cell-centered field
  this->compute_E(E_face, m_cdr->getPhase(), E_cell);  // Compute face-centered field
  this->compute_E(E_eb,   m_cdr->getPhase(), E_cell);  // EB-centered field
}

void rk2::advance_rte_k2_stationary(const Real a_dt){
  CH_TIME("rk2::compute_rte_k2_stationary");
  if(m_verbosity > 5){
    pout() << "rk2::compute_rte_k2_stationary" << endl;
  }

  const Real time = m_time + a_dt; // Source terms centered on the end of the time step

  Vector<EBAMRCellData*> rte_states;
  Vector<EBAMRCellData*> rte_sources;
  Vector<EBAMRCellData*> cdr_states;

  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver   = solver_it();
    RefCountedPtr<RtStorage>& storage = this->get_RtStorage(solver_it);
    
    EBAMRCellData& phi    = storage->getPhi();
    EBAMRCellData& state  = solver->getPhi();
    EBAMRCellData& source = solver->getSource();

    
    // We computed the intermediate potential at time t_k + alpha*dt. Linearly extrapolate that result to the end of
    // the time step. The result for this is y_extrap = y_alpha/alpha - y_0*(1-alpha)/alpha. This (usually) brings the
    // initial guess closer to the true solution.
    if((m_timeStep + 1) % m_fast_rte == 0){
      DataOps::scale(state, -(1.0 - m_alpha)/(m_alpha));
      DataOps::incr(state, phi, 1.0/m_alpha);
    }
    else {
      DataOps::setValue(state, 0.0);
      DataOps::incr(state, phi, 1.0);
    }
    
    rte_states.push_back(&(state));
    rte_sources.push_back(&(source));
  }

  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    cdr_states.push_back(&(solver->getPhi()));
  }

  EBAMRCellData& E = m_fieldSolver_scratch->getElectricFieldCell();

  if((m_timeStep + 1) % m_fast_rte == 0){
    const Real dummy_dt = 0.0;
    this->solve_rte(rte_states, rte_sources, cdr_states, E, time, dummy_dt, centering::cell_center);
  }
}

void rk2::advance_rte_k2_transient(const Real a_dt){
  CH_TIME("rk2::compute_rte_k1_transient");
  if(m_verbosity > 5){
    pout() << "rk2::compute_k1_transient" << endl;
  }

  const Real time = m_time + 0.5*a_dt; // Source terms centered on the half time step

  // If we made it here, the old potential lies in m_fieldSolver_scratch->scratch and the old cdr solutions
  // lie in CdrStorage->m_scratch. The internal state inside the solver is unaffected by anything done previously

  Vector<EBAMRCellData*> rte_states;
  Vector<EBAMRCellData*> rte_sources;
  Vector<EBAMRCellData*> cdr_states;

  for (RtIterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver   = solver_it();
    RefCountedPtr<RtStorage>& storage = this->get_RtStorage(solver_it);
    
    EBAMRCellData& state  = solver->getPhi();  // This has been unaffected so far because we solved onto RtStorage.phi
    EBAMRCellData& source = solver->getSource(); // in the k1-stage. 

    rte_states.push_back(&(state));
    rte_sources.push_back(&(source));
  }

    
  // Source term must be centered between time tn and the intermediate time
  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = this->get_CdrStorage(solver_it);
    
    const EBAMRCellData& state   = solver->getPhi();
    const EBAMRCellData& scratch = storage->getScratch();

    EBAMRCellData& phi = storage->getPhi();

    DataOps::setValue(phi, 0.0);
    DataOps::incr(phi, scratch, 0.5);
    DataOps::incr(phi, phi,     0.5);

    cdr_states.push_back(&(phi));
  }
  
  if((m_timeStep + 1) % m_fast_rte == 0){ // Actual solve
    const MFAMRCellData& state       = m_fieldSolver->getPotential();
    const MFAMRCellData& scratch_phi = m_fieldSolver_scratch->getScratch_phi();
    
    MFAMRCellData& phi   = m_fieldSolver_scratch->getPhi();
    EBAMRCellData& scratch_E   = m_fieldSolver_scratch->getScratch_E();
    
    DataOps::setValue(phi, 0.0);
    DataOps::incr(phi, state, 0.5);
    DataOps::incr(phi, scratch_phi, 0.5);
    
    m_amr->averageDown(phi);
    m_amr->interpGhost(phi);
    
    this->compute_E(scratch_E, m_cdr->getPhase(), phi);
    
    this->solve_rte(rte_states, rte_sources, cdr_states, scratch_E, time, a_dt, centering::cell_center);
  }
}

Real rk2::restrict_dt(){
  return 1.E99;
#include "CD_NamespaceFooter.H"
