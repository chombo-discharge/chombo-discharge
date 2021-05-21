/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
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
  ParmParse pp("AdvectionDiffusion");

  m_phase = phase::gas;

  pp.get("realm",      m_realm);
  pp.get("verbosity",  m_verbosity); 
  pp.get("fhd",        m_fhd); 
  pp.get("diffco",     m_faceCenteredDiffusionCoefficient);
  pp.get("omega",      m_omega);
  pp.get("cfl",        m_cfl);
  pp.get("integrator", m_integrator);
}

AdvectionDiffusionStepper::AdvectionDiffusionStepper(RefCountedPtr<CdrSolver>& a_solver) : AdvectionDiffusionStepper() {
  m_solver = a_solver;
}

AdvectionDiffusionStepper::~AdvectionDiffusionStepper(){
}

void AdvectionDiffusionStepper::parseRuntimeOptions() {
  ParmParse pp("AdvectionDiffusion");

  pp.get("verbosity",  m_verbosity);
  pp.get("cfl",        m_cfl);
  pp.get("integrator", m_integrator);

  m_solver->parseRuntimeOptions();
}

void AdvectionDiffusionStepper::setupSolvers(){
  m_species = RefCountedPtr<CdrSpecies> (new AdvectionDiffusionSpecies());

  // Solver setup
  m_solver->setVerbosity(m_verbosity);
  m_solver->setSpecies(m_species);
  m_solver->parseOptions();
  m_solver->setPhase(m_phase);
  m_solver->setAmr(m_amr);
  m_solver->setComputationalGeometry(m_computationalGeometry);
  m_solver->sanityCheck();
  m_solver->setRealm(m_realm);
  
  if(!m_solver->isMobile() && !m_solver->isDiffusive()){
    MayDay::Abort("AdvectionDiffusionStepper::setupSolvers - can't turn off both advection AND diffusion");
  }
}

void AdvectionDiffusionStepper::postInitialize() {
}

void AdvectionDiffusionStepper::registerRealms() {
  m_amr->registerRealm(m_realm);
}

void AdvectionDiffusionStepper::registerOperators(){
  m_solver->registerOperators();
}

void AdvectionDiffusionStepper::allocate() {
  m_solver->allocateInternals();

  m_amr->allocate(m_tmp, m_realm, m_phase, 1);
  m_amr->allocate(m_k1,  m_realm, m_phase, 1);
  m_amr->allocate(m_k2,  m_realm, m_phase, 1);
}

void AdvectionDiffusionStepper::initialData(){
  m_solver->initialData();       // Fill initial through the cdr species

  m_solver->setSource(0.0);
  m_solver->setEbFlux(0.0);
  m_solver->setDomainFlux(0.0);
  if(m_solver->isDiffusive()){
    m_solver->setDiffusionCoefficient(m_faceCenteredDiffusionCoefficient);
  }
  if(m_solver->isMobile()){
    this->setVelocity();
  }

}

void AdvectionDiffusionStepper::setVelocity(){
  const int finest_level = m_amr->getFinestLevel();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    this->setVelocity(lvl);
  }

  EBAMRCellData& vel = m_solver->getCellCenteredVelocity();
  m_amr->averageDown(vel, m_realm, m_phase);
  m_amr->interpGhost(vel, m_realm, m_phase);
}

void AdvectionDiffusionStepper::setVelocity(const int a_level){
  // TLDR: This code goes down to each cell on grid level a_level and sets the velocity to omega*r
    
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_level];
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box = dbl.get(dit());

    EBCellFAB& vel = (*(m_solver->getCellCenteredVelocity())[a_level])[dit()];
    BaseFab<Real>& vel_reg = vel.getSingleValuedFAB();

    vel.setVal(0.0);
      
    // Regular cells
    for (BoxIterator bit(box); bit.ok(); ++bit){
      const IntVect iv = bit();
      const RealVect pos = m_amr->getProbLo() + (RealVect(iv) + 0.5*RealVect::Unit)*m_amr->getDx()[a_level];

      const Real r     = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
      const Real theta = atan2(pos[1],pos[0]);

      vel_reg(iv,0) = -r*m_omega*sin(theta);
      vel_reg(iv,1) =  r*m_omega*cos(theta);
    }

    // Irregular and multicells
    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_level])[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){

      const VolIndex vof = vofit();
      const IntVect iv   = vof.gridIndex();
      const RealVect pos = m_amr->getProbLo() + (RealVect(iv) + 0.5*RealVect::Unit)*m_amr->getDx()[a_level];

      const Real r     = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
      const Real theta = atan2(pos[1],pos[0]);

      vel(vof,0) = -r*m_omega*sin(theta);
      vel(vof,1) =  r*m_omega*cos(theta);
    }
  }
}

void AdvectionDiffusionStepper::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const {
  m_solver->writeCheckpointLevel(a_handle, a_lvl);
}

void AdvectionDiffusionStepper::readCheckpointData(HDF5Handle& a_handle, const int a_lvl){
  m_solver->readCheckpointLevel(a_handle, a_lvl);
}

void AdvectionDiffusionStepper::postCheckpointSetup(){
  m_solver->setSource(0.0);
  m_solver->setEbFlux(0.0);
  m_solver->setDomainFlux(0.0);
  if(m_solver->isDiffusive()){
    m_solver->setDiffusionCoefficient(m_faceCenteredDiffusionCoefficient);
  }
  if(m_solver->isMobile()){
    this->setVelocity();
  }
}

int AdvectionDiffusionStepper::getNumberOfPlotVariables() const{
  return m_solver->getNumberOfPlotVariables();
}

void AdvectionDiffusionStepper::writePlotData(EBAMRCellData&       a_output,
						  Vector<std::string>& a_plotVariableNames,
						  int&                 a_icomp) const {
  a_plotVariableNames.append(m_solver->getPlotVariableNames());
  m_solver->writePlotData(a_output, a_icomp);
}

void AdvectionDiffusionStepper::computeDt(Real& a_dt, TimeCode& a_timeCode){
  if(m_integrator == 0){      // Explicit diffusion. 
    const Real dt = m_solver->computeAdvectionDiffusionDt();
    
    a_dt       = m_cfl*dt;
    a_timeCode = TimeCode::AdvectionDiffusion;
  }
  else if(m_integrator == 1){ // Implicit diffusion
    const Real dt = m_solver->computeAdvectionDt();

    a_dt       = m_cfl*dt;
    a_timeCode = TimeCode::Advection;
  }
  else{
    MayDay::Abort("AdvectionDiffusionStepper::computeDt - logic bust");
  }
}

Real AdvectionDiffusionStepper::advance(const Real a_dt){
  EBAMRCellData& state = m_solver->getPhi();
  
  if(m_integrator == 0){ //   Use Heun's method
    m_solver->computeDivJ(m_k1, state, 0.0);
    DataOps::copy(m_tmp, state);
    DataOps::incr(m_tmp, m_k1, -a_dt); // m_tmp = phi - dt*div(J)

    m_solver->computeDivJ(m_k2, m_tmp, 0.0);
    DataOps::incr(state, m_k1, -0.5*a_dt);
    DataOps::incr(state, m_k2, -0.5*a_dt); // Done with deterministic update.

    //m_solver->make_non_negative(state);

    // Add random diffusion flux. This is equivalent to a 1st order. Godunov splitting
    if(m_fhd){
      m_solver->gwnDiffusionSource(m_k1, state); // k1 holds random diffusion
      DataOps::incr(state, m_k1, a_dt);
    }
  
    m_amr->averageDown(state, m_realm, m_phase);
    m_amr->interpGhost(state, m_realm, m_phase);
  }
  else if(m_integrator == 1){
    m_solver->computeDivF(m_k1, state, a_dt);
    DataOps::incr(state, m_k1, -a_dt);
    m_amr->averageDown(state, m_realm, m_phase);

    if(m_solver->isDiffusive()){
      DataOps::copy(m_k2, state); // Now holds phiOld - dt*div(F)
      if(m_fhd){
	m_solver->gwnDiffusionSource(m_k1, state); // k1 holds random diffusion
	DataOps::incr(m_k2, m_k1, a_dt);
      }
      DataOps::setValue(m_k1, 0.0);
      m_solver->advanceEuler(state, m_k2, m_k1, a_dt);
    }

    m_amr->averageDown(state, m_realm, m_phase);
    m_amr->averageDown(state, m_realm, m_phase);
  }
  else{
    MayDay::Abort("AdvectionDiffusionStepper - unknown integrator requested");
  }
  
  return a_dt;
}

void AdvectionDiffusionStepper::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt){
  m_timeStep = a_step;
  m_time = a_time;
  m_dt   = a_dt;
  
  m_solver->setTime(a_step, a_time, a_dt);
}

void AdvectionDiffusionStepper::printStepReport(){

}

bool AdvectionDiffusionStepper::needToRegrid(){
  return false;
}

void AdvectionDiffusionStepper::preRegrid(const int a_lbase, const int a_oldFinestLevel){
  m_solver->preRegrid(a_lbase, a_oldFinestLevel);
}

void AdvectionDiffusionStepper::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){

  // Regrid CDR solver
  m_solver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_solver->setSource(0.0);
  m_solver->setEbFlux(0.0);
  m_solver->setDomainFlux(0.0);
  if(m_solver->isDiffusive()){
    m_solver->setDiffusionCoefficient(m_faceCenteredDiffusionCoefficient);
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
}

#include <CD_NamespaceFooter.H>
