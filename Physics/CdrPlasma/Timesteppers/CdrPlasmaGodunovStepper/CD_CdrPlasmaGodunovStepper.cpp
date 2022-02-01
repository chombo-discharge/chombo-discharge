/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaGodunovStepper.cpp
  @brief  Implementation of CD_CdrPlasmaGodunovStepper.H
  @author Robert Marskar
*/

// Std includes
#include <limits>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_CdrPlasmaGodunovStepper.H>
#include <CD_CdrPlasmaGodunovStorage.H>
#include <CD_Timer.H>
#include <CD_DataOps.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

typedef CdrPlasmaGodunovStepper::CdrStorage     CdrStorage;
typedef CdrPlasmaGodunovStepper::FieldStorage   FieldStorage;
typedef CdrPlasmaGodunovStepper::RtStorage      RtStorage;
typedef CdrPlasmaGodunovStepper::SigmaStorage   SigmaStorage;

CdrPlasmaGodunovStepper::CdrPlasmaGodunovStepper(RefCountedPtr<CdrPlasmaPhysics>& a_physics){
  CH_TIME("CdrPlasmaGodunovStepper::CdrPlasmaGodunovStepper()");

  // Default settings
  m_className    = "CdrPlasmaGodunovStepper";
  m_physics       = a_physics;
  m_extrapAdvect = true;
}

CdrPlasmaGodunovStepper::~CdrPlasmaGodunovStepper(){
  CH_TIME("CdrPlasmaGodunovStepper::~CdrPlasmaGodunovStepper()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::~CdrPlasmaGodunovStepper" << endl;
  }
}

void CdrPlasmaGodunovStepper::parseOptions(){
  CH_TIME("CdrPlasmaGodunovStepper::parseOptions()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::parseOptions()" << endl;
  }

  // Parse various options. 
  this->parseVerbosity();
  this->parseSolverVerbosity();
  this->parseFastPoisson();
  this->parseFastRadiativeTransfer();
  this->parseCFL();
  this->parseRelaxationTime();
  this->parseMinDt();
  this->parseMaxDt();
  this->parseSourceComputation();
  this->parseDiffusion();
  this->parseTransport();
  this->parseAdvection();
  this->parseFloor();
  this->parseDebug();
  this->parseFHD();
}

void CdrPlasmaGodunovStepper::parseRuntimeOptions(){
  CH_TIME("CdrPlasmaGodunovStepper::parseRuntimeOptions()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::parseRuntimeOptions()" << endl;
  }

  this->parseVerbosity();
  this->parseSolverVerbosity();
  this->parseFastPoisson();
  this->parseFastRadiativeTransfer();
  this->parseCFL();
  this->parseRelaxationTime();
  this->parseMinDt();
  this->parseMaxDt();
  this->parseSourceComputation();
  this->parseDiffusion();
  this->parseTransport();
  this->parseAdvection();
  this->parseFloor();
  this->parseDebug();
  this->parseFHD();

  // Solvers also parse their runtime options. 
  m_cdr        ->parseRuntimeOptions();
  m_rte        ->parseRuntimeOptions();
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
    m_implicitDiffusion = false;
    m_whichDiffusionAlgorithm = WhichDiffusionAlgorithm::Explicit;
  }
  else if(str == "implicit"){
    m_implicitDiffusion = true;
    m_whichDiffusionAlgorithm = WhichDiffusionAlgorithm::Implicit;
  }
  else if(str == "auto"){
    m_implicitDiffusion = true;
    m_whichDiffusionAlgorithm = WhichDiffusionAlgorithm::Automatic;
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
    m_transportAlgorithm = TransportAlgorithm::Euler;
  }
  else if(str == "rk2"){
    m_transportAlgorithm = TransportAlgorithm::RK2;
  }
  else if(str == "semi_implicit"){
    m_transportAlgorithm = TransportAlgorithm::SemiImplicit;
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
    m_extrapAdvect = true;
  }
  else if(str == "false"){
    m_extrapAdvect = false;
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
  CH_TIME("CdrPlasmaGodunovStepper::parseDebug()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::parseDebug()" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("debug", m_debug);
}

void CdrPlasmaGodunovStepper::parseFHD(){
  CH_TIME("CdrPlasmaGodunovStepper::parseFHD()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::parseFHD()" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("fhd", m_fhd);
}

RefCountedPtr<CdrStorage>& CdrPlasmaGodunovStepper::getCdrStorage(const CdrIterator<CdrSolver>& a_solverIt){
  return m_cdrScratch[a_solverIt.index()];
}

RefCountedPtr<RtStorage>& CdrPlasmaGodunovStepper::getRtStorage(const RtIterator<RtSolver>& a_solverIt){
  return m_rteScratch[a_solverIt.index()];
}

Real CdrPlasmaGodunovStepper::restrictDt(){
  CH_TIME("CdrPlasmaGodunovStepper::restrictDt()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::restrictDt()" << endl;
  }

  return std::numeric_limits<Real>::max();
}

Real CdrPlasmaGodunovStepper::advance(const Real a_dt){
  CH_TIME("CdrPlasmaGodunovStepper::advance(Real)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::advance(Real)" << endl;
  }

  // INFO: When we enter here, the CdrSolvers should have been filled with velocities and diffusion coefficients.
  //
  // The next steps are:
  //
  //    1. Advance transport
  //    2. Solve the Poisson equation (special option for semi-implicit transport)
  //    3. Solve the reactive problem.
  //    4. Solve the radiative transfer problem
  //    5. Put data back into solvers to prepare for the next time step. 

  Timer timer("GodunovStepper::advance");

  // 1. Solve the transport problem. Note that we call advanceTransport which holds the implementation. This differs for explicit and semi-implicit formulations.
  timer.startEvent("Transport");
  CdrPlasmaGodunovStepper::advanceTransport(a_dt);
  timer.stopEvent("Transport");  

  // 2. Solve the Poisson equation and compute the electric field. If we did a semi-implicit solve then the field has already been computed.
  if(m_transportAlgorithm != TransportAlgorithm::SemiImplicit){  
    timer.startEvent("Poisson");
    if((m_timeStep +1) % m_fastPoisson == 0){
      CdrPlasmaStepper::solvePoisson();
    }
    CdrPlasmaGodunovStepper::computeElectricFieldIntoScratch();
    timer.stopEvent("Poisson");
  }

  // 3. Solve the reactive problem. 
  timer.startEvent("Reactions");
  CdrPlasmaGodunovStepper::computeCdrGradients();        
  CdrPlasmaGodunovStepper::computeReactionNetwork(a_dt); 
  timer.stopEvent("Reactions");  

  // 4. Solve the radiative transfer problem. 
  timer.startEvent("Photons");
  if((m_timeStep +1) % m_fastRTE == 0){
    CdrPlasmaGodunovStepper::advanceRadiativeTransfer(a_dt);
  }
  timer.stopEvent("Photons");  

  // Do post step operations. 
  timer.startEvent("Post-step");
  CdrPlasmaGodunovStepper::postStep();
  timer.stopEvent("Post-step");  
  
  // 5. Update velocities and diffusion coefficients in order to prepare for the next time step. 
  timer.startEvent("Velocities");
  CdrPlasmaGodunovStepper::computeCdrDriftVelocities(m_time + a_dt);
  timer.stopEvent("Velocities");

  timer.startEvent("Diffu-coeffs");
  CdrPlasmaGodunovStepper::computeCdrDiffusionCoefficients(m_time + a_dt);
  timer.stopEvent("Diffu-coeffs");

  if(m_debug){
    timer.eventReport(pout(), false);
  }
  
  return a_dt;
}

void CdrPlasmaGodunovStepper::preRegrid(const int a_lbase, const int a_finestLevel) {
  CH_TIME("CdrPlasmaGodunovStepper::preRegrid(int, int)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::preRegrid(int, int)" << endl;
  }

  // Call the parent method -- this will put all solvers in pre-regrid mode. 
  CdrPlasmaStepper::preRegrid(a_lbase, a_finestLevel);

  // When we regrid we must compute the field on the new grid using the semi-implicit update. However, we only have
  // the conductivities and space charge on the old grids, so we must store them and interpolate them to the new grids.
  if(m_transportAlgorithm == TransportAlgorithm::SemiImplicit){
    m_amr->allocate(m_scratchConductivity,    m_realm, phase::gas, 1);
    m_amr->allocate(m_scratchSemiImplicitRho, m_realm, phase::gas, 1);

    DataOps::copy(m_scratchConductivity,    m_conductivityFactorCell);
    DataOps::copy(m_scratchSemiImplicitRho, m_semiImplicitRho       );    
  }
}

void CdrPlasmaGodunovStepper::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) {
  CH_TIME("CdrPlasmaGodunovStepper::regrid(int, int, int)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::regrid(int, int, int)" << endl;
  }

  // TLDR: If we are not using a semi-implicit scheme then we can just call the parent method. Modifications are needed for the
  //       semi-implicit scheme because we solve for the field using div((eps + dt*sigma^k/eps0)E^(k+1)) = -rho^(k+1)/eps0 but
  //       this means we need the conductivity and space charge at the previous time step for restoring the field on the new mesh.

  if(m_transportAlgorithm != TransportAlgorithm::SemiImplicit || m_timeStep == 0){ // Just use regular regrid when we start the simulation. 
    CdrPlasmaStepper::regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  }
  else{
    const Interval interv(0, 0);

    // Call parent method for setting up storage and regridding solvers. 
    this->allocateInternals();
    this->regridSolvers  (a_lmin, a_oldFinestLevel, a_newFinestLevel);
    this->regridInternals(a_lmin, a_oldFinestLevel, a_newFinestLevel);

    // When computing the electric field on the new mesh, we need the conductivity and space charge from the
    // previous step. CdrPlasmaGodunvoStepper::preRegrid will have been called before this routine, so we
    // have the conductivities (stored as sigma*dt/eps0) and space charge stored in scratch data holders. 
    for (int lvl = 0; lvl <= std::max(0, a_lmin-1); lvl++){
      m_scratchConductivity   [lvl]->copyTo(*m_conductivityFactorCell[lvl]);
      m_scratchSemiImplicitRho[lvl]->copyTo(*m_semiImplicitRho       [lvl]);      
    }

    // Regrid AMR levels. 
    for (int lvl = std::max(1, a_lmin); lvl <= a_newFinestLevel; lvl++){
      RefCountedPtr<EBPWLFineInterp>& interpolator = m_amr->getPwlInterpolator(m_realm, m_phase)[lvl];

      // Linearly interpolate (with limiters):
      interpolator->interpolate(*m_conductivityFactorCell[lvl], *m_conductivityFactorCell[lvl-1], interv);
      interpolator->interpolate(*m_semiImplicitRho       [lvl], *m_semiImplicitRho       [lvl-1], interv);

      // For parts of the new grid that overlap with the old grid, replace data with a copy.
      if(lvl <= std::min(a_oldFinestLevel, a_newFinestLevel)){
	m_scratchConductivity   [lvl]->copyTo(*m_conductivityFactorCell[lvl]);
	m_scratchSemiImplicitRho[lvl]->copyTo(*m_semiImplicitRho       [lvl]);      	
      }
    }

    // Coarsen the conductivity and space charge from the last step and update ghost cells. 
    m_amr->averageDown  (m_conductivityFactorCell, m_realm, m_phase);
    m_amr->interpGhostMG(m_conductivityFactorCell, m_realm, m_phase);

    m_amr->averageDown  (m_semiImplicitRho, m_realm, m_phase);
    m_amr->interpGhostMG(m_semiImplicitRho, m_realm, m_phase);

    // Set up the semi-implicit Poisson equation and solve it. 
    this->computeFaceConductivity (m_conductivityFactorFace, m_conductivityFactorEB, m_conductivityFactorCell);
    this->setupSemiImplicitPoisson(m_conductivityFactorFace, m_conductivityFactorEB, 1.0                     );

    // Now solve for the field on the new grids.
    const bool converged = this->solveSemiImplicitPoisson();

    // If we don't converge, try new Poisson solver settings
    if(!converged){ 
      pout() << "CdrPlasmaGodunovStepper::regrid - Poisson solver failed to converge during semi-implicit regrid." << endl;
    }

    // Now compute drift velocities and diffusion -- using the electric field from last time step on the new mesh. 
    CdrPlasmaStepper::computeCdrDriftVelocities();
    CdrPlasmaStepper::computeCdrDiffusion();

    // If we're doing a stationary RTE solve, we also have to recompute source terms
    if(this->stationaryRTE()){     // Solve RTE equations by using data that exists inside solvers
      const Real dummyDt = 0.0;

      // Need new source terms for RTE equations
      this->advanceReactionNetwork(m_time, dummyDt);
      this->solveRadiativeTransfer(dummyDt);    // Argument does not matter, it's a stationary solver.
    }
  }
}

void CdrPlasmaGodunovStepper::postRegrid() {
  CH_TIME("CdrPlasmaGodunovStepper::postRegrid()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::postRegrid()" << endl;
  }

  CdrPlasmaStepper::postRegrid();

  // Release memory. 
  if(m_transportAlgorithm == TransportAlgorithm::SemiImplicit){
    m_amr->deallocate(m_scratchConductivity   );
    m_amr->deallocate(m_scratchSemiImplicitRho);
  }
}

void CdrPlasmaGodunovStepper::postCheckpointSetup(){
  CH_TIME("CdrPlasmaGodunovStepper::postCheckpointSetup()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::postCheckpointSetup()" << endl;
  }

  // TLDR: Only the semi-implicit part overrides the parent method, and it is done because the field needs to be computed from
  //       a different equation. 

  if(m_transportAlgorithm == TransportAlgorithm::SemiImplicit){
    
    // When we enter this routine we will already have called read the checkpoint data into the conductivityFactor and semiimplicit space charge. We need
    // to set up the field solver with those quantities rather than the regular space charge. 
    this->computeFaceConductivity (m_conductivityFactorFace, m_conductivityFactorEB, m_conductivityFactorCell);
    this->setupSemiImplicitPoisson(m_conductivityFactorFace, m_conductivityFactorEB, 1.0                     );
  
    // Now solve for the field on the new grids.
    const bool converged = this->solveSemiImplicitPoisson();
    if(!converged){
      pout() << "CdrPlasmaGodunovStepper::postCheckpointSetup - Poisson solver failed to converge during restart." << endl;
    }
  
    // Now compute the drift and diffusion velocities and update the source terms. This will be OK because the
    // CdrPlasmaStepper functions fetch the electric field from the solver, but it already has the correct potential. 
    CdrPlasmaStepper::computeCdrDriftVelocities();
    CdrPlasmaStepper::computeCdrDiffusion();

    // If we are doing a stationary RTE then we probably have to update the elliptic equations.
    if(this->stationaryRTE()){
      const Real dummyDt = 0.0;

      // We also need new source terms for RTE equations.
      this->advanceReactionNetwork(m_time, dummyDt);
      this->solveRadiativeTransfer(dummyDt);    
    }
  }
  else{
    CdrPlasmaStepper::postCheckpointSetup();
  }
}

void CdrPlasmaGodunovStepper::regridInternals(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){
  CH_TIME("CdrPlasmaGodunovStepper::regridInternals(int, int, int)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::regridInternals(int, int, int)" << endl;
  }
}

bool CdrPlasmaGodunovStepper::solveSemiImplicitPoisson(){
  CH_TIME("CdrPlasmaGodunovStepper::solveSemiImplicitPoisson");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::solveSemiImplicitPoisson" << endl;
  }

  // Now set the semi-implicit space charge. 
  MFAMRCellData& rho      = m_fieldSolver->getRho();
  EBAMRCellData  rhoPhase = m_amr->alias(m_phase, rho);
  
  DataOps::setValue(rho, 0.0);
  DataOps::copy(rhoPhase, m_semiImplicitRho);

  m_amr->averageDown(rho, m_realm);
  m_amr->interpGhost(rho, m_realm);
  
  m_amr->interpToCentroids(rhoPhase, m_realm, m_phase);
  
  const bool converged = m_fieldSolver->solve(m_fieldSolver->getPotential(),
					      rho,
					      m_sigma->getPhi(),
					      false);

  return converged;
}

void CdrPlasmaGodunovStepper::allocateInternals(){
  CH_TIME("CdrPlasmaGodunovStepper::allocateInternals()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::allocateInternals()" << endl;
  }

  constexpr int nComp = 1;

  // TLDR: Since the advancement method requires some temporary variables (for holding various things at BCs, for example),
  //       we have organized the transient memory into classes. Here, we allocate those classes. 

  // Number of CDR solvers and RTE solvers that we have. 
  const int numCdrSpecies = m_physics->getNumCdrSpecies();
  const int numRteSpecies = m_physics->getNumRtSpecies();

  // Allocate storage for the CDR solvers. 
  m_cdrScratch.resize(numCdrSpecies);
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const int idx = solverIt.index();
    
    m_cdrScratch[idx] = RefCountedPtr<CdrStorage> (new CdrStorage(m_amr, m_realm, m_cdr->getPhase(), nComp));
    m_cdrScratch[idx]->allocateStorage();
  }

  // Allocate RTE storage.
  m_rteScratch.resize(numRteSpecies);
  for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
    const int idx = solverIt.index();
    
    m_rteScratch[idx] = RefCountedPtr<RtStorage> (new RtStorage(m_amr, m_realm, m_rte->getPhase(), nComp));
    m_rteScratch[idx]->allocateStorage();
  }

  // Allocate storage for field solver. 
  m_fieldScratch = RefCountedPtr<FieldStorage> (new FieldStorage(m_amr, m_realm, m_cdr->getPhase(), nComp));
  m_fieldScratch->allocateStorage();
  
  // Allocate storage for surface charge solver. 
  m_sigmaScratch = RefCountedPtr<SigmaStorage> (new SigmaStorage(m_amr, m_realm, m_cdr->getPhase(), nComp));
  m_sigmaScratch->allocateStorage();

  // Set up storage for semi-implicit field solves. 
  m_amr->allocate(m_semiImplicitRho,        m_realm, m_phase, nComp);
  m_amr->allocate(m_conductivityFactorCell, m_realm, m_phase, nComp);
  m_amr->allocate(m_conductivityFactorFace, m_realm, m_phase, nComp);
  m_amr->allocate(m_conductivityFactorEB,   m_realm, m_phase, nComp);

  // Initialize values. 
  DataOps::setValue(m_semiImplicitRho,        0.0);
  DataOps::setValue(m_conductivityFactorCell, 0.0);
  DataOps::setValue(m_conductivityFactorFace, 0.0);
  DataOps::setValue(m_conductivityFactorEB,   0.0);  
}

void CdrPlasmaGodunovStepper::deallocateInternals(){
  CH_TIME("CdrPlasmaGodunovStepper::deallocateInternals()");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::deallocateInternals()" << endl;
  }

  // TLDR: This routine simply deallocates the transient memory used by CdrPlasmaGodunovStepper.

  // Run through CDR solvers and deallocate the transient memory assocaited with them. 
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const int idx = solverIt.index();
    
    m_cdrScratch[idx]->deallocateStorage();
    m_cdrScratch[idx] = RefCountedPtr<CdrStorage>(0);
  }

  // Run through RTE solvers and deallocate the transient memory assocaited with them.   
  for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
    const int idx = solverIt.index();
    
    m_rteScratch[idx]->deallocateStorage();
    m_rteScratch[idx] = RefCountedPtr<RtStorage>(0);
  }

  m_cdrScratch.resize(0);
  m_rteScratch.resize(0);

  m_fieldScratch->deallocateStorage();
  m_fieldScratch = RefCountedPtr<FieldStorage>(0);
  
  m_sigmaScratch->deallocateStorage();
  m_sigmaScratch = RefCountedPtr<SigmaStorage>(0);

  // Deallocate storage for the semi-implicit solve. 
  m_amr->deallocate(m_semiImplicitRho       );
  m_amr->deallocate(m_conductivityFactorCell);
  m_amr->deallocate(m_conductivityFactorFace);
  m_amr->deallocate(m_conductivityFactorEB  );
}

void CdrPlasmaGodunovStepper::computeElectricFieldIntoScratch(){
  CH_TIME("CdrPlasmaGodunovStepper::computeElectricFieldIntoScratch");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::computeElectricFieldIntoScratch" << endl;
  }
  
  EBAMRCellData& E_cell = m_fieldScratch->getElectricFieldCell();
  EBAMRIVData&   E_eb   = m_fieldScratch->getElectricFieldEb();
  EBAMRIFData&   E_dom  = m_fieldScratch->getElectricFieldDomain();

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

  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const int idx = solverIt.index();
    RefCountedPtr<CdrSolver>& solver = solverIt();
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solverIt);

    EBAMRCellData& scratch = storage->getScratch();
    EBAMRCellData& grad    = storage->getGradient();

    scratch.copy(solver->getPhi());

    m_amr->averageDown  (scratch, m_realm, m_phase);
    m_amr->interpGhostMG(scratch, m_realm, m_phase);
    
    m_amr->computeGradient(grad, scratch, m_realm, phase::gas);
    
    m_amr->averageDown  (grad, m_realm, m_cdr->getPhase());
    m_amr->interpGhostMG(grad, m_realm, m_cdr->getPhase());
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
  
  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const RefCountedPtr<CdrSolver>& solver = solverIt();
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solverIt);

    cdr_states.push_back(&(storage->getExtrap()));
    eb_states.push_back(&(storage->getEbState()));
    eb_gradients.push_back(&(storage->getEbGrad()));
    cdr_gradients.push_back(&(storage->getGradient())); // Should already have been computed
  }

  // Extrapolate states to the EB and floor them so we cannot get negative values on the boundary. This
  // won't hurt mass conservation because the mass hasn't been injected yet
  CdrPlasmaStepper::extrapolateToEb(eb_states, m_cdr->getPhase(), cdr_states);
  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const int idx = solverIt.index();
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

  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<CdrStorage>& storage = this->getCdrStorage(solverIt);

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
  for (RtIterator<RtSolver> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<RtSolver>& solver   = solverIt();
    RefCountedPtr<RtStorage>& storage = this->getRtStorage(solverIt);

    EBAMRIVData& flux_eb = storage->getEbFlux();
    solver->computeBoundaryFlux(flux_eb, solver->getPhi());
    extrap_rte_fluxes.push_back(&flux_eb);
  }

  const EBAMRIVData& E = m_fieldScratch->getElectricFieldEb();
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
  
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const RefCountedPtr<CdrSolver>& solver = solverIt();
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solverIt);

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
  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const RefCountedPtr<CdrSolver>& solver = solverIt();
    const int idx = solverIt.index();
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
  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<CdrStorage>& storage = this->getCdrStorage(solverIt);

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
  this->extrapolateVelocitiesToDomainFaces(extrap_cdr_velocities,   m_cdr->getPhase(), cdr_velocities);
  this->computeExtrapolatedDomainFluxes(extrap_cdr_fluxes,     states,             cdr_velocities, m_cdr->getPhase());
  this->extrapolateVectorToDomainFaces(extrap_cdr_gradients,  m_cdr->getPhase(), cdr_gradients);

  // Compute RTE flux on domain faces
  for (RtIterator<RtSolver> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<RtSolver>& solver   = solverIt();
    RefCountedPtr<RtStorage>& storage = this->getRtStorage(solverIt);

    EBAMRIFData& domain_flux = storage->getDomainFlux();
    solver->computeDomainFlux(domain_flux, solver->getPhi());
    extrap_rte_fluxes.push_back(&domain_flux);
  }

  const EBAMRIFData& E = m_fieldScratch->getElectricFieldDomain();

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

  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const RefCountedPtr<CdrSolver>& solver = solverIt();
    const RefCountedPtr<CdrSpecies>& spec  = solverIt.getSpecies();
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
  const EBAMRCellData& E = m_fieldScratch->getElectricFieldCell();

  Vector<EBAMRCellData*> cdr_grad;
  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<CdrStorage>& storage = getCdrStorage(solverIt);

    EBAMRCellData& gradient = storage->getGradient();
    cdr_grad.push_back(&gradient);
  }

  // Compute all source terms
  CdrPlasmaStepper::advanceReactionNetwork(cdr_src, rte_src, cdr_phi, cdr_grad, rte_phi, E, m_time, a_dt);

  // Make phi = phi + S*dt
  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<CdrSolver>& solver = solverIt();

    EBAMRCellData& phi = solver->getPhi();
    const EBAMRCellData& src = solver->getSource();

    DataOps::incr(phi, src, a_dt);
    if(m_floor){
      DataOps::floor(phi, 0.0);
    }

#if 1
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

  switch(m_transportAlgorithm){
  case TransportAlgorithm::Euler:
    this->advanceTransportEuler(a_dt);
    break;
  case TransportAlgorithm::RK2:
    this->advanceTransportRK2(a_dt);
    break;
  case TransportAlgorithm::SemiImplicit:
    this->advanceTransportSemiImplicit(a_dt);
    break;
  default:
    MayDay::Error("CdrPlasmaGodunovStepper::advanceTransport - logic bust");
  }
}

void CdrPlasmaGodunovStepper::advanceTransportEuler(const Real a_dt){
  CH_TIME("CdrPlasmaGodunovStepper::advanceTransportEuler");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::advanceTransportEuler" << endl;
  }

  CdrPlasmaGodunovStepper::computeElectricFieldIntoScratch();   // Compute the electric field
  CdrPlasmaGodunovStepper::computeCdrGradients();               // Compute cdr gradients
  CdrPlasmaGodunovStepper::extrapolateSourceTerm(a_dt);         // If we used advective extrapolation, BCs are more work. 
  CdrPlasmaGodunovStepper::computeCdrEbStates();                // Extrapolate cell-centered stuff to EB centroids
  CdrPlasmaGodunovStepper::computeCdrEbFluxes();                // Extrapolate cell-centered fluxes to EB centroids
  CdrPlasmaGodunovStepper::computeCdrDomainStates();            // Extrapolate cell-centered states to domain edges
  CdrPlasmaGodunovStepper::computeCdrDomainFluxes();            // Extrapolate cell-centered fluxes to domain edges
  CdrPlasmaGodunovStepper::computeSigmaFlux();                  // Update charge flux for sigma solver    

  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<CdrSolver>& solver   = solverIt();

    if(solver->isMobile() || solver->isDiffusive()){
      RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solverIt);

      EBAMRCellData& phi = solver->getPhi();
      EBAMRCellData& src = solver->getSource();
    
      EBAMRCellData& scratch  = storage->getScratch();
      EBAMRCellData& scratch2 = storage->getScratch2();

      // Compute hyperbolic term into scratch. Also include diffusion term if and only if we're using explicit diffusion
      const Real extrap_dt = (m_extrapAdvect) ? a_dt : 0.0;
      if(!m_implicitDiffusion){
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
      if(m_implicitDiffusion){
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

void CdrPlasmaGodunovStepper::advanceTransportSemiImplicit(const Real a_dt){
  CH_TIME("CdrPlasmaGodunovStepper::advanceTransportSemiImplicit");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::advanceTransportSemiImplicit" << endl;
  }

  // Compute the conductivity first. We store it as sigma^k*a_dt/eps0
  this->computeCellConductivity(m_conductivityFactorCell);
  DataOps::scale(m_conductivityFactorCell, a_dt/Units::eps0);

  m_amr->averageDown  (m_conductivityFactorCell, m_realm, m_phase);
  m_amr->interpGhostMG(m_conductivityFactorCell, m_realm, m_phase);  

  // Average conductivity to faces and set up the semi-implicit poisson equation. 
  this->computeFaceConductivity (m_conductivityFactorFace, m_conductivityFactorEB, m_conductivityFactorCell);
  this->setupSemiImplicitPoisson(m_conductivityFactorFace, m_conductivityFactorEB, 1.0);

  // Compute the modified right-hand side. We store this as rho^k - dt*e * sum(Z * div(D*grad(phi))).
  DataOps::setValue(m_semiImplicitRho, 0.0);
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    const RefCountedPtr<CdrSolver>&  solver  = solverIt();
    const RefCountedPtr<CdrSpecies>& species = solverIt.getSpecies();

    const int Z = species->getChargeNumber();

    if(Z != 0){
      EBAMRCellData& phi = solver->getPhi();
      
      DataOps::incr(m_semiImplicitRho, phi, Z*Units::Qe);

      // If the solver is diffusive we must increment by the diffusion term as well. 
      if(solver->isDiffusive()){
	RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solverIt);
												 
	EBAMRCellData& divDgradPhi = storage->getScratch();

	solver->computeDivD(divDgradPhi, phi, false, false);

	DataOps::incr(m_semiImplicitRho, divDgradPhi, Z*a_dt*Units::Qe);
      }
    }
  }
  
  // Now solve the semi-implicit Poisson. Issue a warning if we didn't converge.
  const bool converged = this->solveSemiImplicitPoisson();
  if(!converged){
    pout() << "CdrPlasmaGodunovStepper::advanceTransportSemiImplicit -- semi-implicit Poisson solve did not converge!" << endl;
  }

  CdrPlasmaGodunovStepper::computeElectricFieldIntoScratch();   // Compute the electric field

  // Recompute velocities
  CdrPlasmaGodunovStepper::computeCdrDriftVelocities(m_time + a_dt);
  
  // Compute boundary conditions
  CdrPlasmaGodunovStepper::extrapolateSourceTerm(a_dt);         // If we used advective extrapolation, BCs are more work. 
  CdrPlasmaGodunovStepper::computeCdrEbStates();                // Extrapolate cell-centered stuff to EB centroids
  CdrPlasmaGodunovStepper::computeCdrEbFluxes();                // Extrapolate cell-centered fluxes to EB centroids
  CdrPlasmaGodunovStepper::computeCdrDomainStates();            // Extrapolate cell-centered states to domain edges
  CdrPlasmaGodunovStepper::computeCdrDomainFluxes();            // Extrapolate cell-centered fluxes to domain edges
  CdrPlasmaGodunovStepper::computeSigmaFlux();                  // Update charge flux for sigma solver    

  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<CdrSolver>& solver   = solverIt();

    if(solver->isMobile() || solver->isDiffusive()){
      RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solverIt);

      EBAMRCellData& phi = solver->getPhi();
    
      EBAMRCellData& scratch  = storage->getScratch();
      EBAMRCellData& scratch2 = storage->getScratch2();

      // Compute hyperbolic term into scratch. Also include diffusion term if and only if we're using explicit diffusion
      const Real extrap_dt = (m_extrapAdvect) ? a_dt : 0.0;
      solver->computeDivF(scratch, phi, extrap_dt, true, true); // For implicit diffusion, sratch is computed as div(v*phi)

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

  CdrPlasmaGodunovStepper::computeElectricFieldIntoScratch();   // Compute the electric field
  CdrPlasmaGodunovStepper::computeCdrGradients();               // Compute cdr gradients
  CdrPlasmaGodunovStepper::extrapolateSourceTerm(a_dt);         // If we used advective extrapolation, BCs are more work. 
  CdrPlasmaGodunovStepper::computeCdrEbStates();                // Extrapolate cell-centered stuff to EB centroids
  CdrPlasmaGodunovStepper::computeCdrEbFluxes();                // Extrapolate cell-centered fluxes to EB centroids
  CdrPlasmaGodunovStepper::computeCdrDomainStates();            // Extrapolate cell-centered states to domain edges
  CdrPlasmaGodunovStepper::computeCdrDomainFluxes();            // Extrapolate cell-centered fluxes to domain edges
  CdrPlasmaGodunovStepper::computeSigmaFlux();                  // Update charge flux for sigma solver    

  // 1. Compute the "k1" coefficient. This is equivalent to using a semi-implicit
  //    euler method as the predictor for Heun's method
  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<CdrSolver>& solver   = solverIt();
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solverIt);

    EBAMRCellData& phi     = solver->getPhi();
    EBAMRCellData& scratch = storage->getScratch();
    EBAMRCellData& k1      = storage->getScratch3();

    DataOps::copy(k1, phi);

    // Compute hyperbolic term into scratch. Also include diffusion term if and only if we're using explicit diffusion
    const Real extrap_dt = (m_extrapAdvect) ? a_dt : 0.0;
    if(!m_implicitDiffusion){
      solver->computeDivJ(scratch, phi, extrap_dt, true, true); // For explicit diffusion, scratch is computed as div(v*phi - D*grad(phi))
    }
    else{
      solver->computeDivF(scratch, phi, extrap_dt, true, true); // For implicit diffusion, sratch is computed as div(v*phi)
    }
    DataOps::scale(scratch, -1.0);     // scratch = -[div(F/J)]
    DataOps::scale(scratch, a_dt);     // scratch = [-div(F/J)]*dt
    DataOps::incr(phi, scratch, 1.0);  // Make phi = phi^k - dt*div(F/J)

    // This is the implicit diffusion code. If we enter this routine then phi = phi^k - dt*div(F) + dt*R
    if(m_implicitDiffusion){
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
    EBAMRIVData& k1    = m_sigmaScratch->getScratch();
    DataOps::copy(k1, sigma);

    // Advance
    const EBAMRIVData& rhs = m_sigma->getFlux();
    DataOps::incr(sigma, rhs, a_dt);

    // Compute k1 from sigma^(k+1) = sigma^k + k1*dt
    DataOps::incr(k1, sigma, -1.0);
    DataOps::scale(k1, -1.0);
  }

  // 2. Compute the electric field and update boundary conditions and kinetic coefficients
  if((m_timeStep +1) % m_fastPoisson == 0){
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
  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<CdrSolver>& solver   = solverIt();
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solverIt);

    EBAMRCellData& phi      = solver->getPhi();
    EBAMRCellData& scratch  = storage->getScratch();
    const EBAMRCellData& k1 = storage->getScratch3();

    // Compute hyperbolic term into scratch. Also include diffusion term if and only if we're using explicit diffusion
    const Real extrap_dt = m_extrapAdvect ? a_dt : 0.0;
    if(!m_implicitDiffusion){
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
    const EBAMRIVData& k1  = m_sigmaScratch->getScratch();
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
  for (RtIterator<RtSolver> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<RtSolver>& solver = solverIt();
    solver->advance(a_dt);
  }
}

void CdrPlasmaGodunovStepper::computeCdrDriftVelocities(const Real a_time){
  CH_TIME("CdrPlasmaGodunovStepper::computeCdrDriftVelocities");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::computeCdrDriftVelocities" << endl;
  }

  Vector<EBAMRCellData*> velocities = m_cdr->getVelocities();
  CdrPlasmaStepper::computeCdrDriftVelocities(velocities, m_cdr->getPhis(), m_fieldScratch->getElectricFieldCell(), a_time);
}

void CdrPlasmaGodunovStepper::computeCdrDiffusionCoefficients(const Real a_time){
  CH_TIME("CdrPlasmaGodunovStepper::computeCdrDiffusionCoefficients");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::computeCdrDiffusionCoefficients" << endl;
  }

  CdrPlasmaStepper::computeCdrDiffusion(m_fieldScratch->getElectricFieldCell(), m_fieldScratch->getElectricFieldEb());
}

void CdrPlasmaGodunovStepper::computeDt(Real& a_dt, TimeCode& a_timeCode){
  CH_TIME("CdrPlasmaGodunovStepper::computeDt(Real, TimeCode)");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::computeDt(Real, TimeCode)" << endl;
  }

  // TLDR: This routine really depends on what algorithms we use:
  //
  //       Explicit or semi-implicit -> restrict by advection, diffusion, and relaxation time
  //       Partially implicit        -> restrict by advection and relaxation time
  //
  // Note that the semi-implicit scheme does not require restriction by the relaxation time, but users will take
  // care of that through the input script. 
  

  // First, figure out what the transport time step must be for explicit and explicit-implicit methods. 
  if(m_whichDiffusionAlgorithm == WhichDiffusionAlgorithm::Explicit){
    m_dtCFL   = m_cdr->computeAdvectionDiffusionDt();
    
    a_timeCode = TimeCode::AdvectionDiffusion;
    a_dt       = m_cfl*m_dtCFL;
  }
  else if(m_whichDiffusionAlgorithm == WhichDiffusionAlgorithm::Implicit){
    m_dtCFL   = m_cdr->computeAdvectionDt();
    a_timeCode = TimeCode::Advection;

    a_dt = m_cfl*m_dtCFL;
  }
  else if (m_whichDiffusionAlgorithm == WhichDiffusionAlgorithm::Automatic){
    // If we are running with "automatic" diffusion, then we switch to implicit diffusion if the diffusive time step
    // is shorter than the advective time step. 
    const Real advectionDt = m_cdr->computeAdvectionDt();
    const Real diffusionDt = m_cdr->computeDiffusionDt();

    // These are the two cases we need: Either diffusion dt < advection dt and we switch to implicit diffusion and
    // only limit by advection. Or we run fully explicit and use both advection and diffusion as a limiting step. 
    if(diffusionDt < advectionDt){
      m_implicitDiffusion = true;

      m_dtCFL   = advectionDt;
      a_dt       = m_cfl*m_dtCFL;
      a_timeCode = TimeCode::Advection;
    }
    else {
      m_implicitDiffusion = false;

      m_dtCFL   = m_cdr->computeAdvectionDiffusionDt();
      a_dt       = m_cfl*m_dtCFL;
      a_timeCode = TimeCode::AdvectionDiffusion;
    }

    m_dtCFL = a_dt/m_cfl;
  }

  // Next, limit by the relaxation time. 
  const Real dtRelax = m_relaxTime*this->computeRelaxationTime();
  if(dtRelax < a_dt){
    a_dt       = dtRelax;
    a_timeCode = TimeCode::RelaxationTime;
  }

  // Limit by lower hardcap.
  if(a_dt < m_minDt){
    a_dt       = m_minDt;
    a_timeCode = TimeCode::Hardcap;
  }

  // Limit by upper hardcap.  
  if(a_dt > m_maxDt){
    a_dt       = m_maxDt;
    a_timeCode = TimeCode::Hardcap;
  }

  // Set time code also for the TimeStepper. 
  m_timeCode = a_timeCode;
}

void CdrPlasmaGodunovStepper::postStep(){
  CH_TIME("CdrPlasmaGodunovStepper::postStep");
  if(m_verbosity > 5){
    pout() << "CdrPlasmaGodunovStepper::postStep" << endl;
  }

  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<CdrSolver> solver = solverIt();

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

  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt){
    RefCountedPtr<CdrSolver>& solver = solverIt();
    RefCountedPtr<CdrStorage>& storage = CdrPlasmaGodunovStepper::getCdrStorage(solverIt);

    const EBAMRCellData& state  = solver->getPhi();
    const EBAMRCellData& source = solver->getSource();

    EBAMRCellData& extrap = storage->getExtrap();

    DataOps::copy(extrap, state);
    if(m_extrapAdvect) {
      DataOps::incr(extrap, source, 0.5*a_dt);
    }
  }
}

#ifdef CH_USE_HDF5
void CdrPlasmaGodunovStepper::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const {
  CH_TIME("CdrPlasmaGodunovStepper::writeCheckpointData(HDF5Handle, int)");
  if(m_verbosity > 3){
    pout() << "CdrPlasmaGodunovStepper::writeCheckpointData(HDF5Handle, int)" << endl;
  }

  // In addition to the checkpointed stuff from the parent class, the semi-implicit scheme also requires us to checkpoint
  // the factors used in the semi-implicit field solve. We write them here. 

  CdrPlasmaStepper::writeCheckpointData(a_handle, a_lvl);

  write(a_handle, *m_semiImplicitRho       [a_lvl], "semiImplicitRho"   );
  write(a_handle, *m_conductivityFactorCell[a_lvl], "conductivityFactor");  
}
#endif

#ifdef CH_USE_HDF5
void CdrPlasmaGodunovStepper::readCheckpointData(HDF5Handle& a_handle, const int a_lvl) {
  CH_TIME("CdrPlasmaGodunovStepper::readCheckpointData(HDF5Handle, int)");
  if(m_verbosity > 3){
    pout() << "CdrPlasmaGodunovStepper::readCheckpointData(HDF5Handle, int)" << endl;
  }

  // In addition to the checkpointed stuff from the parent class, the semi-implicit scheme also requires us to checkpoint
  // the factors used in the semi-implicit field solve. We read them here.   

  const Interval interv(0,0);

  CdrPlasmaStepper::readCheckpointData(a_handle, a_lvl);

  read(a_handle, *m_semiImplicitRho       [a_lvl], "semiImplicitRho",    m_amr->getGrids(m_realm)[a_lvl], interv, false);
  read(a_handle, *m_conductivityFactorCell[a_lvl], "conductivityFactor", m_amr->getGrids(m_realm)[a_lvl], interv, false);  
}
#endif      

#include <CD_NamespaceFooter.H>
