/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_AdvectionDiffusionStepper.cpp
  @brief  Implementation of CD_AdvectionDiffusionStepper.H
  @author Robert Marskar
  @data   March 2020
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_AdvectionDiffusionStepper.H>
#include <CD_AdvectionDiffusionSpecies.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>
  
using namespace Physics::AdvectionDiffusion;

AdvectionDiffusionStepper::AdvectionDiffusionStepper(){
  CH_TIME("AdvectionDiffusionStepper::AdvectionDiffusionStepper");
  
  ParmParse pp("AdvectionDiffusion");

  m_phase = phase::gas;

  pp.get("realm",      m_realm);
  pp.get("verbosity",  m_verbosity); 
  pp.get("fhd",        m_fhd); 
  pp.get("diffco",     m_diffCo);
  pp.get("omega",      m_omega);
  pp.get("cfl",        m_cfl);

  this->parseIntegrator();  
}

AdvectionDiffusionStepper::AdvectionDiffusionStepper(RefCountedPtr<CdrSolver>& a_solver) : AdvectionDiffusionStepper() {
  CH_TIME("AdvectionDiffusionStepper::AdvectionDiffusionStepper(full)");
  
  m_solver = a_solver;
}

AdvectionDiffusionStepper::~AdvectionDiffusionStepper(){
  CH_TIME("AdvectionDiffusionStepper::~AdvectionDiffusionStepper");  
}

void AdvectionDiffusionStepper::parseRuntimeOptions() {
  CH_TIME("AdvectionDiffusionStepper::parseRuntimeOptions");
  
  ParmParse pp("AdvectionDiffusion");
  
  pp.get("verbosity",  m_verbosity);
  pp.get("cfl",        m_cfl);

  this->parseIntegrator();
  
  m_solver->parseRuntimeOptions();
}

void AdvectionDiffusionStepper::parseIntegrator(){
  CH_TIME("AdvectionDiffusionStepper::parseIntegrator");
  
  ParmParse pp("AdvectionDiffusion");

  std::string str;  

  pp.get("integrator", str);

  if(str == "heun"){
    m_integrator = Integrator::Heun;
  }
  else if(str == "imex"){
    m_integrator = Integrator::EulerIMEX;
  }
  else{
    MayDay::Error("AdvectionDiffusionStepper::parseIntegrator -- logic bust");
  }  
}

void AdvectionDiffusionStepper::setupSolvers(){
  CH_TIME("AdvectionDiffusionStepper::setupSolvers");

  CH_assert(!m_solver.isNull());

  // Instantiate the species. 
  m_species = RefCountedPtr<CdrSpecies> (new AdvectionDiffusionSpecies());

  // Set up the solver. 
  m_solver->setVerbosity(m_verbosity);
  m_solver->setSpecies(m_species);
  m_solver->parseOptions();
  m_solver->setPhase(m_phase);
  m_solver->setAmr(m_amr);
  m_solver->setComputationalGeometry(m_computationalGeometry);
  m_solver->setRealm(m_realm);
  
  if(!m_solver->isMobile() && !m_solver->isDiffusive()){
    MayDay::Error("AdvectionDiffusionStepper::setupSolvers - can't turn off both advection AND diffusion");
  }
}

void AdvectionDiffusionStepper::postInitialize() {
  CH_TIME("AdvectionDiffusionStepper::postInitialize");  
}

void AdvectionDiffusionStepper::registerRealms() {
  CH_TIME("AdvectionDiffusionStepper::registerRealms");
  
  m_amr->registerRealm(m_realm);
}

void AdvectionDiffusionStepper::registerOperators(){
  CH_TIME("AdvectionDiffusionStepper::registerOperators");

  // Let the solver do this -- it knows what it needs. 
  m_solver->registerOperators();
}

void AdvectionDiffusionStepper::allocate() {
  CH_TIME("AdvectionDiffusionStepper::allocate");

  // Allocate storage for running the TimeStepper.
  
  m_solver->allocateInternals();

  m_amr->allocate(m_tmp, m_realm, m_phase, 1);
  m_amr->allocate(m_k1,  m_realm, m_phase, 1);
  m_amr->allocate(m_k2,  m_realm, m_phase, 1);
}

void AdvectionDiffusionStepper::initialData(){
  CH_TIME("AdvectionDiffusionStepper::initialData");

  // Fill the solver with initial data from the species.
  m_solver->initialData();

  // Set velocity, diffusion coefficient, and boundary conditions.
  m_solver->setSource(0.0);
  m_solver->setEbFlux(0.0);
  if(m_solver->isDiffusive()){
    m_solver->setDiffusionCoefficient(m_diffCo);
  }
  if(m_solver->isMobile()){
    this->setVelocity();
  }

  // Set flux functions
  auto fluxFunc = [](const RealVect a_pos, const Real a_time){
    return 0.0;
  };

  //  m_solver->setDomainFlux(fluxFunc);
}

void AdvectionDiffusionStepper::setVelocity(){
  CH_TIME("AdvectionDiffusionStepper::setVelocity");
  
  // Set the velocity.
  auto veloFunc = [omega=this->m_omega](const RealVect pos) -> RealVect {
    const Real r     = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
    const Real theta = atan2(pos[1], pos[0]);

    return RealVect(D_DECL(-r*omega*sin(theta), r*omega*cos(theta), 0.));
  };

  m_solver->setVelocity(veloFunc);  

  // Make sure we have updated ghost cells (they are needed for the flux computation)
  EBAMRCellData& vel = m_solver->getCellCenteredVelocity();
  m_amr->averageDown(vel, m_realm, m_phase);
  m_amr->interpGhost(vel, m_realm, m_phase);
}

#ifdef CH_USE_HDF5
void AdvectionDiffusionStepper::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const {
  CH_TIME("AdvectionDiffusionStepper::writeCheckpointData");
  
  m_solver->writeCheckpointLevel(a_handle, a_lvl);
}
#endif

#ifdef CH_USE_HDF5
void AdvectionDiffusionStepper::readCheckpointData(HDF5Handle& a_handle, const int a_lvl){
  CH_TIME("AdvectionDiffusionStepper::readCheckpointData");
  
  m_solver->readCheckpointLevel(a_handle, a_lvl);
}
#endif

void AdvectionDiffusionStepper::postCheckpointSetup(){
  CH_TIME("AdvectionDiffusionStepper::postCheckpointSetup");

  // Set velocity, diffusion coefficient, and boundary conditions. 
  m_solver->setSource(0.0);
  m_solver->setEbFlux(0.0);
  if(m_solver->isDiffusive()){
    m_solver->setDiffusionCoefficient(m_diffCo);
  }
  if(m_solver->isMobile()){
    this->setVelocity();
  }

  // Set flux functions
  auto fluxFunc = [](const RealVect a_pos, const Real a_time){
    return 0.0;
  };

  //m_solver->setDomainFlux(fluxFunc);  
}

int AdvectionDiffusionStepper::getNumberOfPlotVariables() const{
  CH_TIME("AdvectionDiffusionStepper::getNumberOfPlotVariables");

  // Not plotting anything of our own, so return whatever the solver wants to plot. 
  return m_solver->getNumberOfPlotVariables();
}

void AdvectionDiffusionStepper::writePlotData(EBAMRCellData&       a_output,
					      Vector<std::string>& a_plotVariableNames,
					      int&                 a_icomp) const {
  CH_TIME("AdvectionDiffusionStepper::writePlotData");

  // Append plot variable names
  a_plotVariableNames.append(m_solver->getPlotVariableNames());

  // Solver writes data to be plotted.
  m_solver->writePlotData(a_output, a_icomp);
}

void AdvectionDiffusionStepper::computeDt(Real& a_dt, TimeCode& a_timeCode){
  CH_TIME("AdvectionDiffusionStepper::computeDt");

  // TLDR: If we run explicit advection but implicit diffusion then we are only limited by the advective CFL. Otherwise,
  //       if diffusion is also explicit we need the advection-diffusion limited time step. 

  switch(m_integrator){
  case Integrator::Heun:
    {
      const Real dt = m_solver->computeAdvectionDiffusionDt();
    
      a_dt       = m_cfl*dt;
      a_timeCode = TimeCode::AdvectionDiffusion;
      break;
    }
  case Integrator::IMEX:
    {
      const Real dt = m_solver->computeAdvectionDt();

      a_dt       = m_cfl*dt;
      a_timeCode = TimeCode::Advection;
      break;
    }
  default:
    MayDay::Error("AdvectionDiffusionStepper::computeDt - logic bust");
    break;
  }
}

Real AdvectionDiffusionStepper::advance(const Real a_dt){
  CH_TIME("AdvectionDiffusionStepper::advance");

  // State to be advanced. 
  EBAMRCellData& state = m_solver->getPhi();

  switch(m_integrator){
  case Integrator::Heun:
    {
      const bool conservativeOnly = false;
      const bool addEbFlux        = true ;
      const bool addDomainFlux    = true ;        

      // Compute k1 coefficient
      m_solver->computeDivJ(m_k1, state, 0.0, conservativeOnly, addEbFlux, addDomainFlux);
      DataOps::copy(m_tmp, state);
      DataOps::incr(m_tmp, m_k1, -a_dt);

      // Compute k2 coefficient and final state
      m_solver->computeDivJ(m_k2, m_tmp, 0.0, conservativeOnly, addEbFlux, addDomainFlux);
      DataOps::incr(state, m_k1, -0.5*a_dt);
      DataOps::incr(state, m_k2, -0.5*a_dt); 

      // Add random diffusion flux. This is equivalent to a 1st order. Godunov splitting in an FHD context. 
      if(m_fhd){
	m_solver->gwnDiffusionSource(m_k1, state); // k1 holds random diffusion
	DataOps::incr(state, m_k1, a_dt);
      }
      
      break;
    }
  case Integrator::IMEX:
    {
      const bool addEbFlux        = true ;
      const bool addDomainFlux    = true ;

      if(m_solver->isDiffusive()){

	// Compute the finite volume approximation to kappa*div(F). The second "hook" is a debugging hook that includes redistribution when computing kappa*div(F). It
	// exists only for debugging/assurance reasons. 
	if(true){
	  const bool conservativeOnly = true;
	  
	  m_solver->computeDivF(m_k1, state, a_dt, conservativeOnly, addEbFlux, addDomainFlux);
	}
	else{
	  const bool conservativeOnly = false;
	  
	  m_solver->computeDivF(m_k1, state, a_dt, conservativeOnly, addEbFlux, addDomainFlux);
	  
	  DataOps::kappaScale(m_k1);
	}

	DataOps::scale(m_k1, -1.0);	
	
	// Use m_k1 as the old solution. 
	DataOps::copy(m_k2, state);

	// Do the Euler solve. 
	m_solver->advanceCrankNicholson(state, m_k2, m_k1, a_dt);
      }
      else{ // Purely inviscid advance. 
	m_solver->computeDivF(m_k1, state, a_dt, true, addEbFlux, addDomainFlux);

	DataOps::incr(state, m_k1, -a_dt);
      }
      
      break;
    }
  default:
    MayDay::Error("AdvectionDiffusionStepper - unknown integrator requested");
    break;
  }

  m_amr->averageDown(state, m_realm, m_phase);
  m_amr->interpGhost(state, m_realm, m_phase);  
  
  return a_dt;
}

void AdvectionDiffusionStepper::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt){
  CH_TIME("AdvectionDiffusionStepper::synchronizeSolverTimes");
  
  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;
  
  m_solver->setTime(a_step, a_time, a_dt);
}

void AdvectionDiffusionStepper::printStepReport(){
  CH_TIME("AdvectionDiffusionStepper::printStepReport");
}

bool AdvectionDiffusionStepper::needToRegrid(){
  CH_TIME("AdvectionDiffusionStepper::needToRegrid");
  
  return false;
}

void AdvectionDiffusionStepper::preRegrid(const int a_lbase, const int a_oldFinestLevel){
  CH_TIME("AdvectionDiffusionStepper::preRegrid");

  // TLDR: Call solver preRegrid function -- it knows what to do. 
  m_solver->preRegrid(a_lbase, a_oldFinestLevel);
}

void AdvectionDiffusionStepper::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){
  CH_TIME("AdvectionDiffusionStepper::regrid");
  
  // Regrid CDR solver
  m_solver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);

  // Set up the flow fields again. 
  m_solver->setSource(0.0);
  m_solver->setEbFlux(0.0);
  if(m_solver->isDiffusive()){
    m_solver->setDiffusionCoefficient(m_diffCo);
  }
  if(m_solver->isMobile()){
    this->setVelocity();
  }

  // Allocate memory for RK steps
  m_amr->allocate(m_tmp, m_realm, m_phase, 1);
  m_amr->allocate(m_k1,  m_realm, m_phase, 1);
  m_amr->allocate(m_k2,  m_realm, m_phase, 1);
}

void AdvectionDiffusionStepper::postRegrid() {
  CH_TIME("AdvectionDiffusionStepper::postRegrid");

  // Nothing to see here. 
}

#include <CD_NamespaceFooter.H>
