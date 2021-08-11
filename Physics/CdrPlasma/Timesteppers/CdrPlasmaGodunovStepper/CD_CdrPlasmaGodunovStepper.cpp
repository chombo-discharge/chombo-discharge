/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaGodunovStepper.cpp
  @brief  Implementation of CD_CdrPlasmaGodunovStepper.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

#include <CD_CdrPlasmaGodunovStepper.H>
#include <CD_CdrPlasmaGodunovStorage.H>
#include <CD_DataOps.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

typedef CdrPlasmaGodunovStepper::CdrStorage     CdrStorage;
typedef CdrPlasmaGodunovStepper::FieldStorage   FieldStorage;
typedef CdrPlasmaGodunovStepper::RtStorage      RtStorage;
typedef CdrPlasmaGodunovStepper::SigmaStorage   SigmaStorage;

CdrPlasmaGodunovStepper::CdrPlasmaGodunovStepper(){
  m_className = "CdrPlasmaGodunovStepper";
  m_extrap_advect = true;
}

CdrPlasmaGodunovStepper::CdrPlasmaGodunovStepper(RefCountedPtr<CdrPlasmaPhysics>& a_physics){
  m_className    = "CdrPlasmaGodunovStepper";
  m_physics       = a_physics;
  m_extrap_advect = true;
}

CdrPlasmaGodunovStepper::~CdrPlasmaGodunovStepper(){

}

void CdrPlasmaGodunovStepper::parseOptions(){
  CH_TIME("CdrPlasmaGodunovStepper::parseOptions");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::parseOptions" << endl;
  }

  parseVerbosity();
  parseSolverVerbosity();
  parseFastPoisson();
  parseFastRadiativeTransfer();
  parseCFL();
  parseRelaxationTime();
  parseMinDt();
  parseMaxDt();
  parseSourceComputation();
  parseDiffusion();
  parseTransport();
  parseAdvection();
  parseFloor();
  parseDebug();
  parseFHD();
}

void CdrPlasmaGodunovStepper::parseRuntimeOptions(){
  CH_TIME("CdrPlasmaGodunovStepper::parseRuntimeOptions");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::parseRuntimeOptions" << endl;
  }

  parseVerbosity();
  parseSolverVerbosity();
  parseFastPoisson();
  parseFastRadiativeTransfer();
  parseCFL();
  parseRelaxationTime();
  parseMinDt();
  parseMaxDt();
  parseSourceComputation();
  parseDiffusion();
  parseTransport();
  parseAdvection();
  parseFloor();
  parseDebug();
  parseFHD();

  m_cdr->parseRuntimeOptions();
  m_rte->parseRuntimeOptions();
  m_fieldSolver->parseRuntimeOptions();
}

void CdrPlasmaGodunovStepper::parseDiffusion(){
  CH_TIME("CdrPlasmaGodunovStepper::parseDiffusion");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::parseDiffusion" << endl;
  }

  ParmParse pp("CdrPlasmaGodunovStepper");

  std::string str;
  pp.get("diffusion", str);
  if(str == "explicit"){
    m_implicit_diffusion = false;
    m_WhichDiffusionAlgorithm = WhichDiffusionAlgorithm::Explicit;
  }
  else if(str == "implicit"){
    m_implicit_diffusion = true;
    m_WhichDiffusionAlgorithm = WhichDiffusionAlgorithm::Implicit;
  }
  else if(str == "auto"){
    m_implicit_diffusion = true;
    m_WhichDiffusionAlgorithm = WhichDiffusionAlgorithm::Automatic;
  }
  else{
    MayDay::Abort("CdrPlasmaGodunovStepper_maruyama::parseDiffusion - unknown diffusion type requested");
  }
}

void CdrPlasmaGodunovStepper::parseTransport(){
  CH_TIME("CdrPlasmaGodunovStepper::parseTransport");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::parseTransport" << endl;
  }

  ParmParse pp("CdrPlasmaGodunovStepper");

  std::string str;
  pp.get("transport", str);
  if(str == "euler"){
    m_WhichTransportAlgorithm = WhichTransportAlgorithm::Euler;
  }
  else if(str == "rk2"){
    m_WhichTransportAlgorithm = WhichTransportAlgorithm::RK2;
  }
  else{
    MayDay::Abort("CdrPlasmaGodunovStepper::parseTransport - unknown transport algorithm requested");
  }
}

void CdrPlasmaGodunovStepper::parseAdvection(){
  CH_TIME("CdrPlasmaGodunovStepper::parseAdvection");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::parseAdvection" << endl;
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
    MayDay::Abort("CdrPlasmaGodunovStepper::parseAdvection - unknown argument");
  }
}

void CdrPlasmaGodunovStepper::parseFloor(){
  CH_TIME("CdrPlasmaGodunovStepper::parseFloor");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::parseFloor" << endl;
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
    MayDay::Abort("CdrPlasmaGodunovStepper::parseFloor - unknown argument requested.");
  }
}

void CdrPlasmaGodunovStepper::parseDebug(){
  CH_TIME("CdrPlasmaGodunovStepper::parseDebug");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::parseDebug" << endl;
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
    MayDay::Abort("CdrPlasmaGodunovStepper::parseDebug - unknown argument requested.");
  }
}

void CdrPlasmaGodunovStepper::parseFHD(){
  CH_TIME("CdrPlasmaGodunovStepper::parseFHD");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::parseFHD" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("fhd", m_fhd);
}

bool CdrPlasmaGodunovStepper::needToRegrid(){
  CH_TIME("CdrPlasmaGodunovStepper::deallocateInternals");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::needToRegrid" << endl;
  }

  return false;
}

RefCountedPtr<CdrStorage>& CdrPlasmaGodunovStepper::getCdrStorage(const CdrIterator<CdrSolver>& a_solverit){
  return m_cdr_scratch[a_solverit.index()];
}

RefCountedPtr<RtStorage>& CdrPlasmaGodunovStepper::getRtStorage(const RtIterator<RtSolver>& a_solverit){
  return m_rte_scratch[a_solverit.index()];
}

Real CdrPlasmaGodunovStepper::restrictDt(){
  CH_TIME("CdrPlasmaGodunovStepper::CdrPlasmaGodunovStepper");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::CdrPlasmaGodunovStepper" << endl;
  }

  return 1.E99;
}

Real CdrPlasmaGodunovStepper::advance(const Real a_dt){
  CH_TIME("CdrPlasmaGodunovStepper::advance");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::advance" << endl;
  }

  // INFO: When we enter here, the CdrSolver have updated ghost cells (either linear or quadratic) and should have been
  // filled with velocities and diffusion coefficients.
  // 
  // We must still do:
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
  Real t_tran = 0.0;
  Real t_rte  = 0.0;
  Real t_pois = 0.0;
  Real t_filE = 0.0;
  Real t_filV = 0.0;
  Real t_filD = 0.0;
  Real t_post = 0.0;
  Real t_tot  = 0.0;

  Real t0, t1;
  t0 = MPI_Wtime();
  t_tot  = -t0;

  // These calls are responsible for filling CDR and sigma solver boundary conditions
  // on the EB and on the domain walls
  t0 = MPI_Wtime();
  CdrPlasmaGodunovStepper::computeElectricFieldIntoScratch();       // Compute the electric field
  CdrPlasmaGodunovStepper::computeCdrGradients();        // Compute cdr gradients
  CdrPlasmaGodunovStepper::extrapolateSourceTerm(a_dt);            // If we used advective extrapolation, BCs are more work. 
  CdrPlasmaGodunovStepper::computeCdrEbStates();        // Extrapolate cell-centered stuff to EB centroids
  CdrPlasmaGodunovStepper::computeCdrEbFluxes();        // Extrapolate cell-centered fluxes to EB centroids
  CdrPlasmaGodunovStepper::computeCdrDomainStates();    // Extrapolate cell-centered states to domain edges
  CdrPlasmaGodunovStepper::computeCdrDomainFluxes();    // Extrapolate cell-centered fluxes to domain edges
  CdrPlasmaGodunovStepper::computeSigmaFlux();           // Update charge flux for sigma solver
  t1 = MPI_Wtime();
  t_filBC = t1 - t0;
  
  t0 = MPI_Wtime();
  CdrPlasmaGodunovStepper::advanceTransport(a_dt);
  t1 = MPI_Wtime();
  t_tran = t1 - t0;
  
  t0 = MPI_Wtime();
  if((m_timeStep +1) % m_fast_poisson == 0){
    CdrPlasmaStepper::solvePoisson();         // Update the Poisson equation
  }
  t1 = MPI_Wtime();
  t_pois = t1 - t0;

  t0 = MPI_Wtime();
  CdrPlasmaGodunovStepper::computeElectricFieldIntoScratch();       // Update electric fields too
  t1 = MPI_Wtime();
  t_filE = t1-t0;

  t0 = MPI_Wtime();
  //  CdrPlasmaGodunovStepper::computeCdrGradients();        // I'm not doing this. Gradients don't change that much. 
  CdrPlasmaGodunovStepper::computeReactionNetwork(a_dt); // Advance the reaction network. Put the result in solvers
  t1 = MPI_Wtime();
  t_reac = t1-t0;

  // Move Photons
  t0 = MPI_Wtime();
  if((m_timeStep +1) % m_fast_rte == 0){
    CdrPlasmaGodunovStepper::advanceRadiativeTransfer(a_dt);              // Update RTE equations
  }
  t1 = MPI_Wtime();
  t_rte = t1-t0;

  // Post step
  t0 = MPI_Wtime();
  CdrPlasmaGodunovStepper::postStep();
  t1 = MPI_Wtime();
  t_post = t1 - t0;
  
  // Update velocities and diffusion coefficients. We don't do sources here.
  t0 = MPI_Wtime();
  CdrPlasmaGodunovStepper::computeCdrDriftVelocities(m_time + a_dt);
  t1 = MPI_Wtime();
  t_filV = t1 - t0;
  t0 = MPI_Wtime();
  CdrPlasmaGodunovStepper::computeCdrDiffusionCoefficients(m_time + a_dt);
  t1 = MPI_Wtime();
  t_filD = t1 - t0;
  t_tot += t1;

  if(m_debug){
    pout() << endl;
    pout() << "CdrPlasmaGodunovStepper::advance breakdown:" << endl
	   << "BC fill   = " << 100.0*t_filBC/t_tot << "%" << endl
	   << "Reactions = " << 100.*t_reac/t_tot << "%" << endl
	   << "Transport = " << 100.*t_tran/t_tot << "%" << endl
	   << "RTE adv.  = " << 100.*t_rte/t_tot << "%" << endl
	   << "Poisson   = " << 100.*t_pois/t_tot << "%" << endl
	   << "Ecomp     = " << 100.*t_filE/t_tot << "%" << endl
	   << "post      = " << 100.*t_post/t_tot << "%" << endl
	   << "Vel       = " << 100.*t_filV/t_tot << "%" << endl
	   << "Dco       = " << 100.*t_filD/t_tot << "%" << endl
	   << "TOTAL = " << t_tot << "seconds" << endl;
    pout() << endl;
  }
  
  return a_dt;
}

void CdrPlasmaGodunovStepper::init(){
  CH_TIME("CdrPlasmaGodunovStepper::init");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::init" << endl;
  }

  // No need to do anything in this routine yet
}

void CdrPlasmaGodunovStepper::regridInternals(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){
  CH_TIME("CdrPlasmaGodunovStepper::regridInternals");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::regridInternals" << endl;
  }

  // Nothing to see here
}

void CdrPlasmaGodunovStepper::allocateInternals(){
  CH_TIME("CdrPlasmaGodunovStepper::allocateInternals");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::allocateInternals" << endl;
  }

  const int ncomp       = 1;
  const int num_species = m_physics->getNumCdrSpecies();
  const int num_Photons = m_physics->getNumRtSpecies();

  // Allocate cdr storage
  m_cdr_scratch.resize(num_species);
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    m_cdr_scratch[idx] = RefCountedPtr<CdrStorage> (new CdrStorage(m_amr, m_realm, m_cdr->getPhase(), ncomp));
    m_cdr_scratch[idx]->allocateStorage();
  }

  // Allocate RTE storage
  m_rte_scratch.resize(num_Photons);
  for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    m_rte_scratch[idx] = RefCountedPtr<RtStorage> (new RtStorage(m_amr, m_realm, m_rte->getPhase(), ncomp));
    m_rte_scratch[idx]->allocateStorage();
  }

  // Allocate Poisson storage
  m_fieldSolver_scratch = RefCountedPtr<FieldStorage> (new FieldStorage(m_amr, m_realm, m_cdr->getPhase(), ncomp));
  m_fieldSolver_scratch->allocateStorage();
  
  // Allocate sigma storage
  m_sigma_scratch = RefCountedPtr<SigmaStorage> (new SigmaStorage(m_amr, m_realm, m_cdr->getPhase(), ncomp));
  m_sigma_scratch->allocateStorage();
}

void CdrPlasmaGodunovStepper::deallocateInternals(){
  CH_TIME("CdrPlasmaGodunovStepper::deallocateInternals");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::deallocateInternals" << endl;
  }

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    m_cdr_scratch[idx]->deallocateStorage();
    m_cdr_scratch[idx] = RefCountedPtr<CdrStorage>(0);
  }

  for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    m_rte_scratch[idx]->deallocateStorage();
    m_rte_scratch[idx] = RefCountedPtr<RtStorage>(0);
  }

  m_cdr_scratch.resize(0);
  m_rte_scratch.resize(0);

  m_fieldSolver_scratch->deallocateStorage();
  m_fieldSolver_scratch = RefCountedPtr<FieldStorage>(0);
  
  m_sigma_scratch->deallocateStorage();
  m_sigma_scratch = RefCountedPtr<SigmaStorage>(0);
}

void CdrPlasmaGodunovStepper::computeElectricFieldIntoScratch(){
  CH_TIME("CdrPlasmaGodunovStepper::computeElectricFieldIntoScratch");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::computeElectricFieldIntoScratch" << endl;
  }
  
  EBAMRCellData& E_cell = m_fieldSolver_scratch->getElectricFieldCell();
  EBAMRIVData&   E_eb   = m_fieldSolver_scratch->getElectricFieldEb();
  EBAMRIFData&   E_dom  = m_fieldSolver_scratch->getElectricFieldDomain();

  const MFAMRCellData& phi = m_fieldSolver->getPotential();

  CdrPlasmaStepper::computeElectricField(E_cell, m_cdr->getPhase(), phi);     // Compute cell-centered field
  CdrPlasmaStepper::computeElectricField(E_eb,   m_cdr->getPhase(), E_cell);  // EB-centered field
  CdrPlasmaStepper::extrapolateToDomainFaces(E_dom, m_cdr->getPhase(), E_cell); // Domain centered field
}

void CdrPlasmaGodunovStepper::computeCdrGradients(){
  CH_TIME("CdrPlasmaGodunovStepper::computeCdrGradients");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::computeCdrGradients" << endl;
  }

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    RefCountedPtr<CdrSolver>& solver = solver_it();
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solver_it);

    EBAMRCellData& grad = storage->getGradient();
    m_amr->computeGradient(grad, solver->getPhi(), m_realm, phase::gas);
    m_amr->averageDown(grad, m_realm, m_cdr->getPhase());
    m_amr->interpGhost(grad, m_realm, m_cdr->getPhase());
  }
}

void CdrPlasmaGodunovStepper::computeCdrEbStates(){
  CH_TIME("CdrPlasmaGodunovStepper::computeCdrEbStates");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::computeCdrEbStates" << endl;
  }

  Vector<EBAMRCellData*> cdr_states;
  Vector<EBAMRIVData*>   eb_gradients;
  Vector<EBAMRIVData*>   eb_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solver_it);

    cdr_states.push_back(&(storage->getExtrap()));
    eb_states.push_back(&(storage->getEbState()));
    eb_gradients.push_back(&(storage->getEbGrad()));
    cdr_gradients.push_back(&(storage->getGradient())); // Should already have been computed
  }

  // Extrapolate states to the EB and floor them so we cannot get negative values on the boundary. This
  // won't hurt mass conservation because the mass hasn't been injected yet
  CdrPlasmaStepper::extrapolateToEb(eb_states, m_cdr->getPhase(), cdr_states);
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const int idx = solver_it.index();
    DataOps::floor(*eb_states[idx], 0.0);
  }

  // We should already have the cell-centered gradients, extrapolate them to the EB and project the flux. 
  EBAMRIVData eb_gradient;
  m_amr->allocate(eb_gradient, m_realm, m_cdr->getPhase(), SpaceDim);
  for (int i = 0; i < cdr_states.size(); i++){
    CdrPlasmaStepper::extrapolateToEb(eb_gradient, m_cdr->getPhase(), *cdr_gradients[i]);
    CdrPlasmaStepper::projectFlux(*eb_gradients[i], eb_gradient);
  }
}

void CdrPlasmaGodunovStepper::computeCdrEbFluxes(){
  CH_TIME("CdrPlasmaGodunovStepper::computeCdrEbFluxes()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::computeCdrEbFluxes()";
  }

  Vector<EBAMRCellData*> states;
  Vector<EBAMRIVData*> cdr_fluxes;
  Vector<EBAMRIVData*> extrap_cdr_fluxes;
  Vector<EBAMRIVData*> extrap_cdr_densities;
  Vector<EBAMRIVData*> extrap_cdr_velocities;
  Vector<EBAMRIVData*> extrap_cdr_gradients;
  Vector<EBAMRIVData*> extrap_rte_fluxes;

  cdr_fluxes = m_cdr->getEbFlux();

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrStorage>& storage = this->getCdrStorage(solver_it);

    EBAMRIVData& dens_eb = storage->getEbState();
    EBAMRIVData& velo_eb = storage->getEbVelo();
    EBAMRIVData& flux_eb = storage->getEbFlux();
    EBAMRIVData& grad_eb = storage->getEbGrad();
    EBAMRCellData& state = storage->getExtrap();

    states.push_back(&state);
    extrap_cdr_densities.push_back(&dens_eb);  // Computed in computeCdrEbStates
    extrap_cdr_velocities.push_back(&velo_eb); // Not yet computed
    extrap_cdr_fluxes.push_back(&flux_eb);     // Not yet computed
    extrap_cdr_gradients.push_back(&grad_eb);  // Computed in computeCdrEbStates
  }

  // Compute extrapolated fluxes and velocities at the EB
  Vector<EBAMRCellData*> cdr_velocities = m_cdr->getVelocities();
  CdrPlasmaStepper::computeExtrapolatedFluxes(extrap_cdr_fluxes, states, cdr_velocities, m_cdr->getPhase());
  CdrPlasmaStepper::computeExtrapolatedVelocities(extrap_cdr_velocities, cdr_velocities, m_cdr->getPhase());

  // Compute RTE flux on the boundary
  for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver   = solver_it();
    RefCountedPtr<RtStorage>& storage = this->getRtStorage(solver_it);

    EBAMRIVData& flux_eb = storage->getEbFlux();
    solver->computeBoundaryFlux(flux_eb, solver->getPhi());
    extrap_rte_fluxes.push_back(&flux_eb);
  }

  const EBAMRIVData& E = m_fieldSolver_scratch->getElectricFieldEb();
  CdrPlasmaStepper::computeCdrFluxes(cdr_fluxes,
				     extrap_cdr_fluxes,
				     extrap_cdr_densities,
				     extrap_cdr_velocities,
				     extrap_cdr_gradients,
				     extrap_rte_fluxes,
				     E,
				     m_time);
}

void CdrPlasmaGodunovStepper::computeCdrDomainStates(){
  CH_TIME("CdrPlasmaGodunovStepper::computeCdrDomainStates");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::computeCdrDomainStates" << endl;
  }

  Vector<EBAMRIFData*>   domain_gradients;
  Vector<EBAMRIFData*>   domain_states;
  Vector<EBAMRCellData*> cdr_states;
  Vector<EBAMRCellData*> cdr_gradients;
  
  for (auto solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solver_it);

    cdr_states.push_back(&(storage->getExtrap()));
    domain_states.push_back(&(storage->getDomainState()));
    domain_gradients.push_back(&(storage->getDomainGrad()));
    cdr_gradients.push_back(&(storage->getGradient())); // Should already be computed
  }

  // Extrapolate states to the domain faces
  CdrPlasmaStepper::extrapolateToDomainFaces(domain_states, m_cdr->getPhase(), cdr_states);

  // We already have the cell-centered gradients, extrapolate them to the EB and project the flux.
  EBAMRIFData grad;
  m_amr->allocate(grad, m_realm, m_cdr->getPhase(), SpaceDim);
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const int idx = solver_it.index();
    if(solver->isMobile()){
      CdrPlasmaStepper::extrapolateToDomainFaces(grad, m_cdr->getPhase(), *cdr_gradients[idx]);
      CdrPlasmaStepper::projectDomain(*domain_gradients[idx], grad);
    }
    else{
      DataOps::setValue(*domain_gradients[idx], 0.0);
    }
  }
}

void CdrPlasmaGodunovStepper::computeCdrDomainFluxes(){
  CH_TIME("CdrPlasmaGodunovStepper::computeCdrDomainFluxes()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::computeCdrDomainFluxes()" << endl;
  }

  Vector<EBAMRCellData*> states;
  Vector<EBAMRIFData*>   cdr_fluxes;
  Vector<EBAMRIFData*>   extrap_cdr_fluxes;
  Vector<EBAMRIFData*>   extrap_cdr_densities;
  Vector<EBAMRIFData*>   extrap_cdr_velocities;
  Vector<EBAMRIFData*>   extrap_cdr_gradients;
  Vector<EBAMRIFData*>   extrap_rte_fluxes;

  Vector<EBAMRCellData*> cdr_velocities;
  Vector<EBAMRCellData*> cdr_gradients;

  cdr_fluxes = m_cdr->getDomainFlux();
  cdr_velocities = m_cdr->getVelocities();
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrStorage>& storage = this->getCdrStorage(solver_it);

    EBAMRIFData& dens_domain = storage->getDomainState();
    EBAMRIFData& velo_domain = storage->getDomainVelo();
    EBAMRIFData& flux_domain = storage->getDomainFlux();
    EBAMRIFData& grad_domain = storage->getDomainGrad();
    EBAMRCellData& gradient  = storage->getGradient();
    EBAMRCellData& state     = storage->getExtrap();

    states.push_back(&state);
    extrap_cdr_densities.push_back(&dens_domain);  // Has not been computed
    extrap_cdr_velocities.push_back(&velo_domain); // Has not been computed
    extrap_cdr_fluxes.push_back(&flux_domain);     // Has not been computed
    extrap_cdr_gradients.push_back(&grad_domain);  // Has not been computed
    cdr_gradients.push_back(&gradient);
  }

  // Compute extrapolated velocities and fluxes at the domain faces
  this->extrapolateToDomainFaces(extrap_cdr_densities,         m_cdr->getPhase(), states);
  this->extrapolateVelocitiesVectorDomainFaces(extrap_cdr_velocities,   m_cdr->getPhase(), cdr_velocities);
  this->computeExtrapolatedDomainFluxes(extrap_cdr_fluxes,     states,             cdr_velocities, m_cdr->getPhase());
  this->extrapolateToVectorDomainFaces(extrap_cdr_gradients,  m_cdr->getPhase(), cdr_gradients);

  // Compute RTE flux on domain faces
  for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver   = solver_it();
    RefCountedPtr<RtStorage>& storage = this->getRtStorage(solver_it);

    EBAMRIFData& domain_flux = storage->getDomainFlux();
    solver->computeDomainFlux(domain_flux, solver->getPhi());
    extrap_rte_fluxes.push_back(&domain_flux);
  }

  const EBAMRIFData& E = m_fieldSolver_scratch->getElectricFieldDomain();

  // This fills the solvers' domain fluxes
  CdrPlasmaStepper::computeCdrDomainFluxes(cdr_fluxes,
					   extrap_cdr_fluxes,
					   extrap_cdr_densities,
					   extrap_cdr_velocities,
					   extrap_cdr_gradients,
					   extrap_rte_fluxes,
					   E,
					   m_time);
}

void CdrPlasmaGodunovStepper::computeSigmaFlux(){
  CH_TIME("CdrPlasmaGodunovStepper::computeSigmaFlux");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::computeSigmaFlux" << endl;
  }

  EBAMRIVData& flux = m_sigma->getFlux();
  DataOps::setValue(flux, 0.0);

  for (auto solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<CdrSolver>& solver = solver_it();
    const RefCountedPtr<CdrSpecies>& spec  = solver_it.getSpecies();
    const EBAMRIVData& solver_flux          = solver->getEbFlux();

    DataOps::incr(flux, solver_flux, spec->getChargeNumber()*Units::Qe);
  }

  m_sigma->resetCells(flux);
}

void CdrPlasmaGodunovStepper::computeReactionNetwork(const Real a_dt){
  CH_TIME("CdrPlasmaGodunovStepper::computeReactionNetwork");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::computeReactionNetwork" << endl;
  }

  // We have already computed E and the gradients of the CDR equations, so we will call the
  // CdrPlasmaStepper version where all that crap is inputs. Saves memory and flops. 

  Vector<EBAMRCellData*> cdr_src = m_cdr->getSources();
  Vector<EBAMRCellData*> cdr_phi = m_cdr->getPhis();
  Vector<EBAMRCellData*> rte_src = m_rte->getSources();
  Vector<EBAMRCellData*> rte_phi = m_rte->getPhis();
  const EBAMRCellData& E = m_fieldSolver_scratch->getElectricFieldCell();

  Vector<EBAMRCellData*> cdr_grad;
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrStorage>& storage = getCdrStorage(solver_it);

    EBAMRCellData& gradient = storage->getGradient();
    cdr_grad.push_back(&gradient);
  }

  // Compute all source terms
  CdrPlasmaStepper::advanceReactionNetwork(cdr_src, rte_src, cdr_phi, cdr_grad, rte_phi, E, m_time, a_dt);

  // Make phi = phi + S*dt
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();

    EBAMRCellData& phi = solver->getPhi();
    const EBAMRCellData& src = solver->getSource();

    DataOps::incr(phi, src, a_dt);
    if(m_floor){
      DataOps::floor(phi, 0.0);
    }

#if 0
    if(m_floor){ // Should we floor or not? Usually a good idea, and you can monitor the (hopefully negligible) injected mass
      if(m_debug){
	const Real mass_before = solver->computeMass();
	DataOps::floor(phi, 0.0);
	const Real mass_after = solver->computeMass();
	const Real rel_mass = (mass_after-mass_before)/mass_before;
	pout() << "CdrPlasmaGodunovStepper::computeReactionNetwork - injecting relative "
	       << solver->getName() << " mass = " << rel_mass << endl;
      }
      else{
	DataOps::floor(phi, 0.0);

      }
    }
#endif
  }
}

void CdrPlasmaGodunovStepper::advanceTransport(const Real a_dt){
  CH_TIME("CdrPlasmaGodunovStepper::advanceTransport");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::advanceTransport" << endl;
  }

  // Transport algorithm
  if(m_WhichTransportAlgorithm == WhichTransportAlgorithm::Euler){
    CdrPlasmaGodunovStepper::advanceTransportEuler(a_dt);
  }
  else if(m_WhichTransportAlgorithm == WhichTransportAlgorithm::RK2){
    CdrPlasmaGodunovStepper::advanceTransportRK2(a_dt);
  }
  
}

void CdrPlasmaGodunovStepper::advanceTransportEuler(const Real a_dt){
  CH_TIME("CdrPlasmaGodunovStepper::advanceTransportEuler");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::advanceTransportEuler" << endl;
  }

  for (auto solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();

    if(solver->isMobile() || solver->isDiffusive()){
      RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solver_it);

      EBAMRCellData& phi = solver->getPhi();
      EBAMRCellData& src = solver->getSource();
    
      EBAMRCellData& scratch  = storage->getScratch();
      EBAMRCellData& scratch2 = storage->getScratch2();

      // Compute hyperbolic term into scratch. Also include diffusion term if and only if we're using explicit diffusion
      const Real extrap_dt = (m_extrap_advect && solver->extrapolateSourceTerm()) ? a_dt : 0.0;
      if(!m_implicit_diffusion){
	solver->computeDivJ(scratch, phi, extrap_dt, true, true); // For explicit diffusion, scratch is computed as div(v*phi - D*grad(phi))
      }
      else{
	solver->computeDivF(scratch, phi, extrap_dt, true, true); // For implicit diffusion, sratch is computed as div(v*phi)
      }
      DataOps::scale(scratch, -1.0);     // scratch = -[div(F/J)]
      DataOps::scale(scratch, a_dt);     // scratch = [-div(F/J)]*dt
      DataOps::incr(phi, scratch, 1.0);  // Make phi = phi^k - dt*div(F/J)

      // Add random flux
      if(m_fhd && solver->isDiffusive()){
	solver->gwnDiffusionSource(scratch2, phi);
	DataOps::incr(phi, scratch2, a_dt);
      }

      if(m_floor){ // Should we floor or not? Usually a good idea, and you can monitor the (hopefully negligible) injected mass
	if(m_debug){
	  const Real mass_before = solver->computeMass();
	  DataOps::floor(phi, 0.0);
	  const Real mass_after = solver->computeMass();
	  const Real rel_mass = (mass_after-mass_before)/mass_before;
	  pout() << "CdrPlasmaGodunovStepper::advance_cdr - injecting relative " << solver->getName() << " mass = " << rel_mass << endl;
	}
	else{
	  DataOps::floor(phi, 0.0);
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
	  DataOps::copy(scratch, phi); // Weird-ass initial solution, as explained above
	  DataOps::setValue(scratch2, 0.0); // No source, those are a part of the initial solution
	  solver->advanceEuler(phi, scratch, scratch2, a_dt);

	  if(m_floor){ // Should we floor or not? Usually a good idea, and you can monitor the (hopefully negligible) injected mass
	    if(m_debug){
	      const Real mass_before = solver->computeMass();
	      DataOps::floor(phi, 0.0);
	      const Real mass_after = solver->computeMass();
	      const Real rel_mass = (mass_after-mass_before)/mass_before;
	      pout() << "CdrPlasmaGodunovStepper::advance_cdr (implicit diffusion) - inecting relative  "
		     << solver->getName() << " mass = " << rel_mass << endl;
	    }
	    else{
	      DataOps::floor(phi, 0.0);
	    }
	  }
	}
      }

      m_amr->averageDown(phi, m_realm, m_cdr->getPhase());
      m_amr->interpGhost(phi, m_realm, m_cdr->getPhase());
    }
  }

  // Advance the sigma equation
  EBAMRIVData& sigma = m_sigma->getPhi();
  const EBAMRIVData& rhs = m_sigma->getFlux();
  DataOps::incr(sigma, rhs, a_dt);
}

void CdrPlasmaGodunovStepper::advanceTransportRK2(const Real a_dt){
  CH_TIME("CdrPlasmaGodunovStepper::advanceTransportRK2");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::advanceTransportRK2" << endl;
  }

  // 1. Compute the "k1" coefficient. This is equivalent to using a semi-implicit
  //    euler method as the predictor for Heun's method
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solver_it);

    EBAMRCellData& phi     = solver->getPhi();
    EBAMRCellData& scratch = storage->getScratch();
    EBAMRCellData& k1      = storage->getScratch3();

    DataOps::copy(k1, phi);

    // Compute hyperbolic term into scratch. Also include diffusion term if and only if we're using explicit diffusion
    const Real extrap_dt = (m_extrap_advect && solver->extrapolateSourceTerm()) ? a_dt : 0.0;
    if(!m_implicit_diffusion){
      solver->computeDivJ(scratch, phi, extrap_dt, true, true); // For explicit diffusion, scratch is computed as div(v*phi - D*grad(phi))
    }
    else{
      solver->computeDivF(scratch, phi, extrap_dt, true, true); // For implicit diffusion, sratch is computed as div(v*phi)
    }
    DataOps::scale(scratch, -1.0);     // scratch = -[div(F/J)]
    DataOps::scale(scratch, a_dt);     // scratch = [-div(F/J)]*dt
    DataOps::incr(phi, scratch, 1.0);  // Make phi = phi^k - dt*div(F/J)

    // This is the implicit diffusion code. If we enter this routine then phi = phi^k - dt*div(F) + dt*R
    if(m_implicit_diffusion){
      DataOps::floor(phi, 0.0);
      
      // Solve implicit diffusion equation. This looks weird but we're solving
      //
      // phi^(k+1) = phi^k - dt*div(F) + dt*R + dt*div(D*div(phi^k+1))
      //
      // This discretization is equivalent to a diffusion-only discretization with phi^k -dt*div(F) + dt*R as initial solution
      // so we just use that for simplicity
      if(solver->isDiffusive()){
	EBAMRCellData& scratch2 = storage->getScratch2();
	
	DataOps::copy(scratch, phi);       // Make scratch = phiOld - dt*div(F/J)
	DataOps::setValue(scratch2, 0.0); // No source, those are a part of the initial solution
	solver->advanceEuler(phi, scratch, scratch2, a_dt);
      }
    }

    DataOps::floor(phi, 0.0);

    m_amr->averageDown(phi, m_realm, m_cdr->getPhase());
    m_amr->interpGhost(phi, m_realm, m_cdr->getPhase());

    // Compute k1 from phi^(k+1) = phi^k + k1
    DataOps::incr(k1, phi, -1.0);
    DataOps::scale(k1, -1.0);
  }

  // Update the sigma equation
  {
    EBAMRIVData& sigma = m_sigma->getPhi();
    EBAMRIVData& k1    = m_sigma_scratch->getScratch();
    DataOps::copy(k1, sigma);

    // Advance
    const EBAMRIVData& rhs = m_sigma->getFlux();
    DataOps::incr(sigma, rhs, a_dt);

    // Compute k1 from sigma^(k+1) = sigma^k + k1*dt
    DataOps::incr(k1, sigma, -1.0);
    DataOps::scale(k1, -1.0);
  }

  // 2. Compute the electric field and update boundary conditions and kinetic coefficients
  if((m_timeStep +1) % m_fast_poisson == 0){
    CdrPlasmaStepper::solvePoisson();         // Update the Poisson equation
    CdrPlasmaGodunovStepper::computeElectricFieldIntoScratch();       // Update electric fields too
  }
  CdrPlasmaGodunovStepper::computeCdrDriftVelocities(m_time + a_dt);
  CdrPlasmaGodunovStepper::computeCdrDiffusionCoefficients(m_time + a_dt);
  CdrPlasmaGodunovStepper::postStep();
  
  CdrPlasmaGodunovStepper::computeCdrGradients();        // Recompute cdr gradients after reaction advance (these might have changed)
  CdrPlasmaGodunovStepper::extrapolateSourceTerm(a_dt);            // BCs are more work with advective extrapolation
  CdrPlasmaGodunovStepper::computeCdrEbStates();        // Extrapolate cell-centered stuff to EB centroids
  CdrPlasmaGodunovStepper::computeCdrEbFluxes();        // Extrapolate cell-centered fluxes to EB centroids
  CdrPlasmaGodunovStepper::computeCdrDomainStates();    // Extrapolate cell-centered states to domain edges
  CdrPlasmaGodunovStepper::computeCdrDomainFluxes();    // Extrapolate cell-centered fluxes to domain edges
  CdrPlasmaGodunovStepper::computeSigmaFlux();           // Update charge flux for sigma solver

  // 3. Perform final advance, which will be the solution to 
  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver   = solver_it();
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solver_it);

    EBAMRCellData& phi      = solver->getPhi();
    EBAMRCellData& scratch  = storage->getScratch();
    const EBAMRCellData& k1 = storage->getScratch3();

    // Compute hyperbolic term into scratch. Also include diffusion term if and only if we're using explicit diffusion
    const Real extrap_dt = m_extrap_advect ? a_dt : 0.0;
    if(!m_implicit_diffusion){
      solver->computeDivJ(scratch, phi, extrap_dt, true, true); // For explicit diffusion, scratch is computed as div(v*phi - D*grad(phi))
      DataOps::scale(scratch, -1.0);
      DataOps::scale(scratch, a_dt);
      // Now make phiNew = phiOld + 0.5*(k1+k2)*dt but since phi = phiOld + k1, just do
      // phiNew = phi - 0.5*k1 + 0.5*scratch

      DataOps::incr(phi, k1,     -0.5);
      DataOps::incr(phi, scratch, 0.5);
    }
    else{ // Implicit diffusion is a bit more tricky
      solver->computeDivF(scratch, phi, extrap_dt, true, true); // For implicit diffusion, sratch is computed as div(v*phi)
      DataOps::scale(scratch, -a_dt);     // scratch = [-div(F)]*dt

      // Solve the stinking diffusion equation. This is weird but we want to solve phiNew = phiOld + 0.5*dt*(k1+f(phiNew)),
      // and phi holds phiOld + dt*k1 and f(phiNew) is semi-implicit. So we solve in the form
      //
      // phiNew = phiOld + 0.5*dt*(k1 - divF(phi) + div(D*grad(phiNew)))
      //
      DataOps::incr(phi, k1, -0.5);     // phi = phiOld + 0.5*k1
      DataOps::incr(phi, scratch, 0.5); // phi = phiOld + 0.5*(k1 - div(F))
      
      if(solver->isDiffusive()){
	EBAMRCellData& scratch2 = storage->getScratch2();
	
	DataOps::copy(scratch, phi);       // Weird-ass initial solution, as explained above
	DataOps::setValue(scratch2, 0.0); // No source, those are a part of the initial solution
	
	solver->advanceEuler(phi, scratch, scratch2, a_dt);
      }
    }

    DataOps::floor(phi, 0.0);

    m_amr->averageDown(phi, m_realm, m_cdr->getPhase());
    m_amr->interpGhost(phi, m_realm, m_cdr->getPhase());

    if(m_floor){ // Should we floor or not? Usually a good idea, and you can monitor the (hopefully negligible) injected mass
      if(m_debug){
	const Real mass_before = solver->computeMass();
	DataOps::floor(phi, 0.0);
	const Real mass_after = solver->computeMass();
	const Real rel_mass = (mass_after-mass_before)/mass_before;
	pout() << "CdrPlasmaGodunovStepper::advanceTransportRK2 - injecting relative " << solver->getName() << " mass = " << rel_mass << endl;
      }
      else{
	DataOps::floor(phi, 0.0);
      }
    }
  }

  { // Do the final sigma advance
    EBAMRIVData& sigma     = m_sigma->getPhi();
    const EBAMRIVData& rhs = m_sigma->getFlux();
    const EBAMRIVData& k1  = m_sigma_scratch->getScratch();
    DataOps::incr(sigma, k1, -0.5);
    DataOps::incr(sigma, rhs, 0.5*a_dt);
  }
}

void CdrPlasmaGodunovStepper::advanceRadiativeTransfer(const Real a_dt){
  CH_TIME("CdrPlasmaGodunovStepper::advanceRadiativeTransfer");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::advanceRadiativeTransfer" << endl;
  }

  // Source terms should already be in place so we can solve directly.
  for (RtIterator<RtSolver> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<RtSolver>& solver = solver_it();
    solver->advance(a_dt);
  }
}

void CdrPlasmaGodunovStepper::computeCdrDriftVelocities(const Real a_time){
  CH_TIME("CdrPlasmaGodunovStepper::computeCdrDriftVelocities");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::computeCdrDriftVelocities" << endl;
  }

  Vector<EBAMRCellData*> velocities = m_cdr->getVelocities();
  CdrPlasmaStepper::computeCdrDriftVelocities(velocities, m_cdr->getPhis(), m_fieldSolver_scratch->getElectricFieldCell(), a_time);
}

void CdrPlasmaGodunovStepper::computeCdrDiffusionCoefficients(const Real a_time){
  CH_TIME("CdrPlasmaGodunovStepper::computeCdrDiffusionCoefficients");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::computeCdrDiffusionCoefficients" << endl;
  }

  CdrPlasmaStepper::computeCdrDiffusion(m_fieldSolver_scratch->getElectricFieldCell(), m_fieldSolver_scratch->getElectricFieldEb());
}

void CdrPlasmaGodunovStepper::computeDt(Real& a_dt, TimeCode& a_timeCode){
  CH_TIME("CdrPlasmaGodunovStepper::computeDt");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::computeDt" << endl;
  }

  // Restrict by advection or advection-diffusion. 
  if(m_WhichDiffusionAlgorithm == WhichDiffusionAlgorithm::Explicit){
    m_dt_cfl   = m_cdr->computeAdvectionDiffusionDt();
    
    a_timeCode = TimeCode::AdvectionDiffusion;
    a_dt       = m_cfl*m_dt_cfl;
  }
  else if(m_WhichDiffusionAlgorithm == WhichDiffusionAlgorithm::Implicit){
    m_dt_cfl   = m_cdr->computeAdvectionDt();
    a_timeCode = TimeCode::Advection;

    a_dt = m_cfl*m_dt_cfl;
  }
  else if (m_WhichDiffusionAlgorithm == WhichDiffusionAlgorithm::Automatic){
    const Real advection_dt = m_cdr->computeAdvectionDt();
    const Real diffusion_dt = m_cdr->computeDiffusionDt();

    if(diffusion_dt < advection_dt){
      m_implicit_diffusion = true;

      m_dt_cfl   = advection_dt;
      a_dt       = m_cfl*m_dt_cfl;
      a_timeCode = TimeCode::Advection;
    }
    else {
      m_implicit_diffusion = false;

      m_dt_cfl   = m_cdr->computeAdvectionDiffusionDt();
      a_dt       = m_cfl*m_dt_cfl;
      a_timeCode = TimeCode::AdvectionDiffusion;
    }

    m_dt_cfl = a_dt/m_cfl;
  }

  // Below here we restrict by relaxation time and hardcaps. 
  const Real dt_relax = m_relax_time*this->computeRelaxationTime();
  if(dt_relax < a_dt){
    a_dt = dt_relax;
    a_timeCode = TimeCode::RelaxationTime;
  }

  if(a_dt < m_min_dt){
    a_dt = m_min_dt;
    a_timeCode = TimeCode::Hardcap;
  }

  if(a_dt > m_max_dt){
    a_dt = m_max_dt;
    a_timeCode = TimeCode::Hardcap;
  }

  m_timeCode = a_timeCode;
}

void CdrPlasmaGodunovStepper::postStep(){
  CH_TIME("CdrPlasmaGodunovStepper::postStep");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::postStep" << endl;
  }

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver> solver = solver_it();

    m_amr->averageDown(solver->getPhi(), m_realm, m_cdr->getPhase());
    m_amr->interpGhost(solver->getPhi(), m_realm, m_cdr->getPhase());
  }
}

void CdrPlasmaGodunovStepper::extrapolateSourceTerm(const Real a_dt){
  CH_TIME("CdrPlasmaGodunovStepper::extrapolateSourceTerm");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::extrapolateSourceTerm" << endl;
  }

  // TLDR: If we extrapolate the advective derivative and include the source term in the extrapolation,
  //       the boundary conditions should be computed from phi + 0.5*source*dt rather than just phi.

  for (CdrIterator<CdrSolver> solver_it = m_cdr->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<CdrSolver>& solver = solver_it();
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solver_it);

    const EBAMRCellData& state  = solver->getPhi();
    const EBAMRCellData& source = solver->getSource();

    EBAMRCellData& extrap = storage->getExtrap();

    DataOps::copy(extrap, state);
    if(m_extrap_advect && solver->extrapolateSourceTerm()) {
      DataOps::incr(extrap, source, 0.5*a_dt);
    }
  }
}

#include <CD_NamespaceFooter.H>
