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
#include <CD_CdrGodunov.H>

#include <ParmParse.H>

typedef euler_maruyama::cdr_storage     cdr_storage;
typedef euler_maruyama::poisson_storage poisson_storage;
typedef euler_maruyama::rte_storage     rte_storage;
typedef euler_maruyama::sigma_storage   sigma_storage;

euler_maruyama::euler_maruyama(){
  m_className = "euler_maruyama";
  m_extrap_advect = true;
}

euler_maruyama::~euler_maruyama(){

}

void euler_maruyama::parseOptions(){
  CH_TIME("euler_maruyama::parseOptions");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::parseOptions" << endl;
  }

  parseVerbosity();
  parse_solver_verbosity();
  parse_fast_poisson();
  parse_cfl();
  parse_relax_time();
  parse_min_dt();
  parse_max_dt();
  parse_source_comp();
  parse_diffusion();
  parse_advection();
  parse_floor();
  parse_debug();
}

void euler_maruyama::parseRuntimeOptions(){
  CH_TIME("euler_maruyama::parseRuntimeOptions");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::parseRuntimeOptions" << endl;
  }

  parseVerbosity();
  parse_solver_verbosity();
  parse_fast_poisson();
  parse_cfl();
  parse_relax_time();
  parse_min_dt();
  parse_max_dt();
  parse_source_comp();
  parse_diffusion();
  parse_advection();
  parse_floor();
  parse_debug();

  m_cdr->parseRuntimeOptions();
  m_rte->parseRuntimeOptions();
  m_fieldSolver->parseRuntimeOptions();
}

void euler_maruyama::parse_diffusion(){
  CH_TIME("euler_maruyama::parse_diffusion");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::parse_diffusion" << endl;
  }

  ParmParse pp("euler_maruyama");

  std::string str;
  pp.get("diffusion", str);
  if(str == "explicit"){
    m_implicit_diffusion = false;
    m_whichDiffusion = whichDiffusion::Explicit;
  }
  else if(str == "implicit"){
    m_implicit_diffusion = true;
    m_whichDiffusion = whichDiffusion::Implicit;
  }
  else if(str == "auto"){
    m_implicit_diffusion = true;
    m_whichDiffusion = whichDiffusion::Automatic;
  }
  else{
    MayDay::Abort("euler_maruayama::parse_diffusion - unknown diffusion type requested");
  }
}

void euler_maruyama::parse_advection(){
  CH_TIME("euler_maruyama::parse_advection");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::parse_advection" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("extrap_advect", str);
  if(str == "true"){
    m_extrap_advect = true;
  }
  else if(str == "false"){
    m_extrap_advect = false;
  }
  else{
    MayDay::Abort("euler_maruyama::parse_advection - unknown argument");
  }
}

void euler_maruyama::parse_floor(){
  CH_TIME("euler_maruyama::parse_floor");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::parse_floor" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("floor_cdr", str);
  if(str == "true"){
    m_floor = true;
  }
  else if(str == "false"){
    m_floor = false;
  }
  else{
    MayDay::Abort("euler_maruayma::parse_floor - unknown argument requested.");
  }
}

void euler_maruyama::parse_debug(){
  CH_TIME("euler_maruyama::parse_debug");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::parse_debug" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("debug", str);
  if(str == "true"){
    m_debug = true;
  }
  else if(str == "false"){
    m_debug = false;
  }
  else{
    MayDay::Abort("euler_maruayma::parse_debug - unknown argument requested.");
  }
}

bool euler_maruyama::needToRegrid(){
  CH_TIME("euler_maruyama::deallocateInternals");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::needToRegrid" << endl;
  }

  return false;
}

RefCountedPtr<cdr_storage>& euler_maruyama::get_cdr_storage(const CdrIterator& a_solverit){
  return m_cdr_scratch[a_solverit.index()];
}

RefCountedPtr<rte_storage>& euler_maruyama::get_rte_storage(const rte_iterator& a_solverit){
  return m_rte_scratch[a_solverit.index()];
}

Real euler_maruyama::restrict_dt(){
  CH_TIME("euler_maruyama::euler_maruyama");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::euler_maruyama" << endl;
  }

  return 1.E99;
}

Real euler_maruyama::advance(const Real a_dt){
  CH_TIME("euler_maruyama::advance");
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

  Real t_grad = 0.0;
  Real t_filBC= 0.0;
  Real t_reac = 0.0;
  Real t_cdr  = 0.0;
  Real t_rte  = 0.0;
  Real t_sig  = 0.0;
  Real t_pois = 0.0;
  Real t_filE = 0.0;
  Real t_filV = 0.0;
  Real t_filD = 0.0;
  Real t_tot  = 0.0;

  Real t0, t1;
  t0 = MPI_Wtime();
  t_tot  = -t0;

  // These calls are responsible for filling CDR and sigma solver boundary conditions
  // on the EB and on the domain walls
  t0 = MPI_Wtime();
  euler_maruyama::compute_E_into_scratch();       // Compute the electric field
  euler_maruyama::compute_cdr_gradients();        // Extrapolate cell-centered stuff to EB centroids
  t1 = MPI_Wtime();
  t_grad = t1 - t0;

  t0 = MPI_Wtime();
  euler_maruyama::compute_cdr_eb_states();        // Extrapolate cell-centered stuff to EB centroids
  euler_maruyama::compute_cdr_eb_fluxes();        // Extrapolate cell-centered fluxes to EB centroids
  euler_maruyama::compute_cdr_domain_states();    // Extrapolate cell-centered states to domain edges
  euler_maruyama::compute_cdr_domain_fluxes();    // Extrapolate cell-centered fluxes to domain edges
  euler_maruyama::compute_sigma_flux();           // Update charge flux for sigma solver
  t1 = MPI_Wtime();
  t_filBC = t1 - t0;

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
  if((m_timeStep +1) % m_fast_poisson == 0){
    TimeStepper::solve_poisson();                  // Update the Poisson equation
  }
  t1 = MPI_Wtime();
  t_pois = t1 - t0;

  t0 = MPI_Wtime();
  euler_maruyama::compute_E_into_scratch();       // Update electric fields too
  t1 = MPI_Wtime();
  t_filE = t1-t0;

  // Update velocities and diffusion coefficients. We don't do sources here.
  t0 = MPI_Wtime();
  euler_maruyama::compute_cdr_velo(m_time + a_dt);
  t1 = MPI_Wtime();
  t_filV = t1 - t0;
  t0 = MPI_Wtime();
  euler_maruyama::compute_cdr_diffco(m_time + a_dt);
  t1 = MPI_Wtime();
  t_filD = t1 - t0;
  t_tot += t1;

  if(m_debug){
    pout() << endl;
    pout() << "euler_maruyama::advance breakdown:" << endl
	   << "E & grad  = " << 100.0*t_grad/t_tot << "%" << endl
	   << "BC fill   = " << 100.0*t_filBC/t_tot << "%" << endl
	   << "Reactions = " << 100.*t_reac/t_tot << "%" << endl
	   << "CDR adv.  = " << 100.*t_cdr/t_tot << "%" << endl
	   << "RTE adv.  = " << 100.*t_rte/t_tot << "%" << endl
	   << "Poisson   = " << 100.*t_pois/t_tot << "%" << endl
	   << "Ecomp     = " << 100.*t_filE/t_tot << "%" << endl
	   << "Vel       = " << 100.*t_filV/t_tot << "%" << endl
	   << "Dco       = " << 100.*t_filD/t_tot << "%" << endl
	   << "TOTAL = " << t_tot << "seconds" << endl;
    pout() << endl;
  }
  
  return a_dt;
}

void euler_maruyama::init(){
  CH_TIME("euler_maruyama::init");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::init" << endl;
  }

  // No need to do anything in this routine yet
}

void euler_maruyama::regridInternals(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){
  CH_TIME("euler_maruyama::regridInternals");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::regridInternals" << endl;
  }

  // These aren't the droids you're looking for. 
}

void euler_maruyama::allocateInternals(){
  CH_TIME("euler_maruyama::allocateInternals");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::allocateInternals" << endl;
  }

  const int ncomp       = 1;
  const int num_species = m_plaskin->get_num_species();
  const int num_photons = m_plaskin->get_num_photons();

  // Allocate cdr storage
  m_cdr_scratch.resize(num_species);
  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    m_cdr_scratch[idx] = RefCountedPtr<cdr_storage> (new cdr_storage(m_amr, m_cdr->get_phase(), ncomp));
    m_cdr_scratch[idx]->allocate_storage();
  }

  // Allocate RTE storage
  m_rte_scratch.resize(num_photons);
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    m_rte_scratch[idx] = RefCountedPtr<rte_storage> (new rte_storage(m_amr, m_rte->get_phase(), ncomp));
    m_rte_scratch[idx]->allocate_storage();
  }

  // Allocate Poisson storage
  m_fieldSolver_scratch = RefCountedPtr<poisson_storage> (new poisson_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_fieldSolver_scratch->allocate_storage();
  
  // Allocate sigma storage
  m_sigma_scratch = RefCountedPtr<sigma_storage> (new sigma_storage(m_amr, m_cdr->get_phase(), ncomp));
  m_sigma_scratch->allocate_storage();
}

void euler_maruyama::deallocateInternals(){
  CH_TIME("euler_maruyama::deallocateInternals");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::deallocateInternals" << endl;
  }

  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    m_cdr_scratch[idx]->deallocate_storage();
    m_cdr_scratch[idx] = RefCountedPtr<cdr_storage>(0);
  }

  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    m_rte_scratch[idx]->deallocate_storage();
    m_rte_scratch[idx] = RefCountedPtr<rte_storage>(0);
  }

  m_cdr_scratch.resize(0);
  m_rte_scratch.resize(0);

  m_fieldSolver_scratch->deallocate_storage();
  m_fieldSolver_scratch = RefCountedPtr<poisson_storage>(0);
  
  m_sigma_scratch->deallocate_storage();
  m_sigma_scratch = RefCountedPtr<sigma_storage>(0);
}

void euler_maruyama::compute_E_into_scratch(){
  CH_TIME("euler_maruyama::compute_E_into_scratch");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::compute_E_into_scratch" << endl;
  }
  
  EBAMRCellData& E_cell = m_fieldSolver_scratch->get_E_cell();
  EBAMRIVData&   E_eb   = m_fieldSolver_scratch->get_E_eb();
  EBAMRIFData&   E_dom  = m_fieldSolver_scratch->get_E_domain();

  const MFAMRCellData& phi = m_fieldSolver->getPotential();

  TimeStepper::compute_E(E_cell, m_cdr->get_phase(), phi);     // Compute cell-centered field
  TimeStepper::compute_E(E_eb,   m_cdr->get_phase(), E_cell);  // EB-centered field
  TimeStepper::extrapolate_to_domain_faces(E_dom, m_cdr->get_phase(), E_cell); // Domain centered field
}

void euler_maruyama::compute_cdr_gradients(){
  CH_TIME("euler_maruyama::compute_cdr_gradients");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::compute_cdr_gradients" << endl;
  }

  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    RefCountedPtr<CdrSolver>& solver = solver_it();
    RefCountedPtr<cdr_storage>& storage = euler_maruyama::get_cdr_storage(solver_it);

    EBAMRCellData& grad = storage->get_gradient();
    m_amr->computeGradient(grad, solver->getPhi(), phase::gas);
    m_amr->averageDown(grad, m_cdr->get_phase());
    m_amr->interpGhost(grad, m_cdr->get_phase());
  }
}

void euler_maruyama::compute_cdr_eb_states(){
  CH_TIME("euler_maruyama::compute_cdr_eb_states");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::compute_cdr_eb_states" << endl;
  }

  Vector<EBAMRCellData*> cdr_states = m_cdr->getPhis();
  
  Vector<EBAMRIVData*>   eb_gradients;
  Vector<EBAMRIVData*>   eb_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    RefCountedPtr<cdr_storage>& storage = euler_maruyama::get_cdr_storage(solver_it);

    eb_states.push_back(&(storage->get_eb_state()));
    eb_gradients.push_back(&(storage->get_eb_grad()));
    cdr_gradients.push_back(&(storage->get_gradient())); // Should already have been computed
  }

  // Extrapolate states to the EB and floor them so we cannot get negative values on the boundary. This
  // won't hurt mass conservation because the mass hasn't been injected yet
  TimeStepper::extrapolate_to_eb(eb_states, m_cdr->get_phase(), cdr_states);
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    data_ops::floor(*eb_states[idx], 0.0);
  }

  // We should already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIVData eb_gradient;
  m_amr->allocate(eb_gradient, m_cdr->get_phase(), SpaceDim);
  for (int i = 0; i < cdr_states.size(); i++){
    TimeStepper::extrapolate_to_eb(eb_gradient, m_cdr->get_phase(), *cdr_gradients[i]);
    TimeStepper::project_flux(*eb_gradients[i], eb_gradient);
  }
}

void euler_maruyama::compute_cdr_eb_fluxes(){
  CH_TIME("euler_maruyama::compute_cdr_eb_fluxes()");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::compute_cdr_eb_fluxes()";
  }

  Vector<EBAMRCellData*> states = m_cdr->getPhis();

  Vector<EBAMRIVData*> cdr_fluxes;
  Vector<EBAMRIVData*> extrap_cdr_fluxes;
  Vector<EBAMRIVData*> extrap_cdr_densities;
  Vector<EBAMRIVData*> extrap_cdr_velocities;
  Vector<EBAMRIVData*> extrap_cdr_gradients;
  Vector<EBAMRIVData*> extrap_rte_fluxes;

  cdr_fluxes = m_cdr->getEbFlux();

  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
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
  TimeStepper::compute_extrapolated_fluxes(extrap_cdr_fluxes, states, cdr_velocities, m_cdr->get_phase());
  TimeStepper::compute_extrapolated_velocities(extrap_cdr_velocities, cdr_velocities, m_cdr->get_phase());

  // Compute RTE flux on the boundary
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);

    EBAMRIVData& flux_eb = storage->get_eb_flux();
    solver->compute_boundary_flux(flux_eb, solver->getPhi());
    extrap_rte_fluxes.push_back(&flux_eb);
  }

  const EBAMRIVData& E = m_fieldSolver_scratch->get_E_eb();
  TimeStepper::compute_cdr_fluxes(cdr_fluxes,
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
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    RefCountedPtr<cdr_storage>& storage = euler_maruyama::get_cdr_storage(solver_it);

    cdr_states.push_back(&(solver->getPhi()));
    domain_states.push_back(&(storage->get_domain_state()));
    domain_gradients.push_back(&(storage->get_domain_grad()));
    cdr_gradients.push_back(&(storage->get_gradient())); // Should already be computed
  }

  // Extrapolate states to the domain faces
  TimeStepper::extrapolate_to_domain_faces(domain_states, m_cdr->get_phase(), cdr_states);

  // We already have the cell-centered gradients, extrapolate them to the EB and project the flux.
  EBAMRIFData grad;
  m_amr->allocate(grad, m_cdr->get_phase(), SpaceDim);
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const int idx = solver_it.index();
    if(solver->isMobile()){
      TimeStepper::extrapolate_to_domain_faces(grad, m_cdr->get_phase(), *cdr_gradients[idx]);
      TimeStepper::project_domain(*domain_gradients[idx], grad);
    }
    else{
      data_ops::set_value(*domain_gradients[idx], 0.0);
    }
  }
}

void euler_maruyama::compute_cdr_domain_fluxes(){
  CH_TIME("euler_maruyama::compute_cdr_domain_fluxes()");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::compute_cdr_domain_fluxes()" << endl;
  }

  Vector<EBAMRCellData*> states = m_cdr->getPhis();

  Vector<EBAMRIFData*>   cdr_fluxes;
  Vector<EBAMRIFData*>   extrap_cdr_fluxes;
  Vector<EBAMRIFData*>   extrap_cdr_densities;
  Vector<EBAMRIFData*>   extrap_cdr_velocities;
  Vector<EBAMRIFData*>   extrap_cdr_gradients;
  Vector<EBAMRIFData*>   extrap_rte_fluxes;

  Vector<EBAMRCellData*> cdr_velocities;
  Vector<EBAMRCellData*> cdr_gradients;

  cdr_fluxes = m_cdr->getDomainFlux();
  cdr_velocities = m_cdr->get_velocities();
  for (CdrIterator solver_it(*m_cdr); solver_it.ok(); ++solver_it){
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
  this->extrapolate_velo_to_domain_faces(extrap_cdr_velocities,   m_cdr->get_phase(), cdr_velocities);
  this->compute_extrapolated_domain_fluxes(extrap_cdr_fluxes,     states,             cdr_velocities, m_cdr->get_phase());
  this->extrapolate_vector_to_domain_faces(extrap_cdr_gradients,  m_cdr->get_phase(), cdr_gradients);

  // Compute RTE flux on domain faces
  for (rte_iterator solver_it(*m_rte); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver   = solver_it();
    RefCountedPtr<rte_storage>& storage = this->get_rte_storage(solver_it);

    EBAMRIFData& domain_flux = storage->get_domain_flux();
    solver->compute_domain_flux(domain_flux, solver->getPhi());
    extrap_rte_fluxes.push_back(&domain_flux);
  }

  const EBAMRIFData& E = m_fieldSolver_scratch->get_E_domain();

  // This fills the solvers' domain fluxes
  TimeStepper::compute_cdr_domain_fluxes(cdr_fluxes,
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
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const RefCountedPtr<species>& spec      = solver_it.getSpecies();
    const EBAMRIVData& solver_flux          = solver->getEbFlux();

    data_ops::incr(flux, solver_flux, spec->getChargeNumber()*units::s_Qe);
  }

  m_sigma->reset_cells(flux);
}

void euler_maruyama::compute_reaction_network(const Real a_dt){
  CH_TIME("euler_maruyama::compute_reaction_network");
  if(m_verbosity > 5){
    pout() << "euler_maruaya::compute_reaction_network" << endl;
  }

  // We have already computed E and the gradients of the CDR equations, so we will call the
  // TimeStepper version where all that crap is inputs. Saves memory and flops. 

  Vector<EBAMRCellData*> cdr_src = m_cdr->getSources();
  Vector<EBAMRCellData*> cdr_phi = m_cdr->getPhis();
  Vector<EBAMRCellData*> rte_src = m_rte->getSources();
  Vector<EBAMRCellData*> rte_phi = m_rte->getPhis();
  const EBAMRCellData& E = m_fieldSolver_scratch->get_E_cell();

  Vector<EBAMRCellData*> cdr_grad;
  for (CdrIterator solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_storage>& storage = get_cdr_storage(solver_it);

    EBAMRCellData& gradient = storage->get_gradient();
    cdr_grad.push_back(&gradient);
  }

  //  TimeStepper::advance_reaction_network(m_time, a_dt);
  TimeStepper::advance_reaction_network(cdr_src, rte_src, cdr_phi, cdr_grad, rte_phi, E, m_time, a_dt);
}

void euler_maruyama::advance_cdr(const Real a_dt){
  CH_TIME("euler_maruyama::advance_cdr");
  if(m_verbosity > 5){
    pout() << "euler_maruaya::advance_cdr" << endl;
  }

  for (auto solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<cdr_storage>& storage = euler_maruyama::get_cdr_storage(solver_it);

    EBAMRCellData& phi = solver->getPhi();
    EBAMRCellData& src = solver->getSource();
    
    EBAMRCellData& scratch  = storage->get_scratch();
    EBAMRCellData& scratch2 = storage->get_scratch2();

    // Compute hyperbolic term into scratch. Also include diffusion term if and only if we're using explicit diffusion
    const Real extrap_dt = m_extrap_advect ? a_dt : 0.0;
    if(!m_implicit_diffusion){
      solver->computeDivJ(scratch, phi, extrap_dt);
    }
    else{
      solver->computeDivF(scratch, phi, extrap_dt);
    }
    data_ops::scale(scratch, -1.0);

    // Increment with source term
    data_ops::incr(scratch, src, 1.0);  // scratch = [-div(F/J) + R]
    data_ops::scale(scratch, a_dt);     // scratch = [-div(F/J) + R]*dt
    data_ops::incr(phi, scratch, 1.0);  // Make phi = phi^k - dt*div(F/J) + dt*R

    solver->make_non_negative(phi);

    if(m_floor){ // Should we floor or not? Usually a good idea, and you can monitor the (hopefully negligible) injected mass
      if(m_debug){
	const Real mass_before = solver->computeMass();
	data_ops::floor(phi, 0.0);
	const Real mass_after = solver->computeMass();
	const Real rel_mass = (mass_after-mass_before)/mass_before;
	pout() << "euler_maruayma::injecting relative " << solver->getName() << " mass = " << rel_mass << endl;
      }
      else{
	data_ops::floor(phi, 0.0);
      }
    }
    

    // This is the implicit diffusion code. If we enter this routine then phi = phi^k - dt*div(F) + dt*R
    if(m_implicit_diffusion){
      // Solve implicit diffusion equation. This looks weird but we're solving
      //
      // phi^(k+1) = phi^k - dt*div(F) + dt*R + dt*div(D*div(phi^k+1))
      //
      // This discretization is equivalent to a diffusion-only discretization with phi^k -dt*div(F) + dt*R as initial solution
      // so we just use that for simplicity
      if(solver->isDiffusive()){
	data_ops::copy(scratch, phi); // Weird-ass initial solution, as explained above
	data_ops::set_value(scratch2, 0.0); // No source, those are a part of the initial solution
	solver->advanceEuler(phi, scratch, scratch2, a_dt);

	solver->make_non_negative(phi);

	if(m_floor){ // Should we floor or not? Usually a good idea, and you can monitor the (hopefully negligible) injected mass
	  if(m_debug){
	    const Real mass_before = solver->computeMass();
	    data_ops::floor(phi, 0.0);
	    const Real mass_after = solver->computeMass();
	    const Real rel_mass = (mass_after-mass_before)/mass_before;
	    pout() << "euler_maruayma::injecting relative " << solver->getName() << " mass = " << rel_mass << endl;
	  }
	  else{
	    data_ops::floor(phi, 0.0);
	  }
	}
      }
    }

    m_amr->averageDown(phi, m_cdr->get_phase());
    m_amr->interpGhost(phi, m_cdr->get_phase());
  }
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
  EBAMRIVData& sigma = m_sigma->getPhi();
  const EBAMRIVData& rhs = m_sigma->get_flux();
  data_ops::incr(sigma, rhs, a_dt);
}

void euler_maruyama::compute_cdr_velo(const Real a_time){
  CH_TIME("euler_maruyama::compute_cdr_velo");
  if(m_verbosity > 5){
    pout() << "euler_maruaya::compute_cdr_velo" << endl;
  }

  Vector<EBAMRCellData*> velocities = m_cdr->get_velocities();
  TimeStepper::compute_cdr_velocities(velocities, m_cdr->getPhis(), m_fieldSolver_scratch->get_E_cell(), a_time);
}

void euler_maruyama::compute_cdr_diffco(const Real a_time){
  CH_TIME("euler_maruyama::compute_cdr_diffco");
  if(m_verbosity > 5){
    pout() << "euler_maruaya::compute_cdr_diffco" << endl;
  }

  TimeStepper::compute_cdr_diffusion(m_fieldSolver_scratch->get_E_cell(), m_fieldSolver_scratch->get_E_eb());
}

void euler_maruyama::computeDt(Real& a_dt, TimeCode::which_code& a_timeCode){
  CH_TIME("euler_maruyama::computeDt");
  if(m_verbosity > 5){
    pout() << "euler_maruyama::computeDt" << endl;
  }

  Real dt = 1.E99;

  m_dt_cfl          = m_cdr->compute_cfl_dt();
  const Real dt_cfl = m_cfl*m_dt_cfl;
  if(dt_cfl < dt){
    dt = dt_cfl;
    a_timeCode = TimeCode::cfl;
  }

  const Real dt_relax = m_relax_time*this->compute_relaxation_time();
  if(dt_relax < dt){
    dt = dt_relax;
    a_timeCode = TimeCode::RelaxationTime;
  }

  if(dt < m_min_dt){
    dt = m_min_dt;
    a_timeCode = TimeCode::Hardcap;
  }

  if(dt > m_max_dt){
    dt = m_max_dt;
    a_timeCode = TimeCode::Hardcap;
  }

  // Diffusion step step constraint. If diffusion dt is the shortest scale, 
  if(m_whichDiffusion == whichDiffusion::Explicit){ // Have to accept time step constraint
    const Real dt_diffusion = m_cdr->compute_diffusive_dt();
    if(dt_diffusion < dt){
      dt = dt_diffusion;
      a_timeCode = TimeCode::Diffusion;
    }
  }
  else if(m_whichDiffusion == whichDiffusion::Automatic){ // If explicit diffusion dt is the shortest, go implicit
    const Real dt_diffusion = m_cdr->compute_diffusive_dt();
    if(dt_diffusion < dt){ // Use implicit diffusion
      m_implicit_diffusion = true;
    }
    else{ // Use explicit diffusion
      m_implicit_diffusion = false;
    }
  }

  a_dt = dt;
#include "CD_NamespaceFooter.H"

