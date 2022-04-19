/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TracerParticleStepperImplem.H
  @brief  Implementation of CD_TracerParticleStepper.H
  @author Robert Marskar
*/

#ifndef CD_TracerParticleStepperImplem_H
#define CD_TracerParticleStepperImplem_H

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_TracerParticleStepper.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::TracerParticle;

template <typename P>
inline
TracerParticleStepper<P>::TracerParticleStepper() {
  CH_TIME("TracerParticleStepper::TracerParticleStepper");

  m_realm = Realm::Primal;
  m_phase = phase::gas;
  
  this->parseOptions();
}

template <typename P>
inline
TracerParticleStepper<P>::~TracerParticleStepper() {
  CH_TIME("TracerParticleStepper::~TracerParticleStepper");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::~TracerParticleStepper" << endl;
  }
}

template <typename P>
inline
void TracerParticleStepper<P>::setupSolvers() {
  CH_TIME("TracerParticleStepper::setupSolvers()");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::setupSolvers()" << endl;
  }

  m_solver = RefCountedPtr<TracerParticleSolver<P> > (new TracerParticleSolver<P>(m_amr, m_computationalGeometry));

  m_solver->setPhase(m_phase);  
  m_solver->setRealm(m_realm);
  m_solver->parseOptions();  
}

template <typename P>
inline
void TracerParticleStepper<P>::allocate() {
  CH_TIME("TracerParticleStepper::allocate()");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::allocate()" << endl;
  }

  m_amr->allocate(m_velocity, m_realm, m_phase, SpaceDim);
  m_solver->allocate();
}

template <typename P>
inline
void TracerParticleStepper<P>::initialData() {
  CH_TIME("TracerParticleStepper::initialData()");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::initialData()" << endl;
  }

  this->setVelocity();
  this->initialParticles();

  m_solver->setVelocity(m_velocity);
}

template <typename P>
inline
void TracerParticleStepper<P>::registerRealms() {
  CH_TIME("TracerParticleStepper::registerRealms()");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::registerRealms()" << endl;
  }

  m_amr->registerRealm(m_realm);
}

template <typename P>
inline
void TracerParticleStepper<P>::registerOperators() {
  CH_TIME("TracerParticleStepper::registerOperators()");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::registerOperators()" << endl;
  }

  m_solver->registerOperators();
}

template <typename P>
inline
void TracerParticleStepper<P>::parseOptions() {
  CH_TIME("TracerParticleStepper::parseOptions()");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::parseOptions()" << endl;
  }

  this->parseIntegrator();  
}

template <typename P>
inline
void TracerParticleStepper<P>::parseRuntimeOptions() {
  CH_TIME("TracerParticleStepper::parseRuntimeOptions()");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::parseRuntimeOptions()" << endl;
  }

  this->parseIntegrator();

  m_solver->parseRuntimeOptions();  
}

template <typename P>
inline
void TracerParticleStepper<P>::parseIntegrator() {
  CH_TIME("TracerParticleStepper::parseIntegrator()");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::parseIntegrator()" << endl;
  }

  ParmParse pp("TracerParticleStepper");

  std::string str;

  pp.get("verbosity",   m_verbosity);
  pp.get("cfl",         m_cfl      );
  pp.get("integration", str        );
  if(str == "euler"){
    m_algorithm = IntegrationAlgorithm::Euler;
  }
  else if(str == "rk2"){
    m_algorithm = IntegrationAlgorithm::RK2;
  }
  else if(str == "rk4"){
    m_algorithm = IntegrationAlgorithm::RK4;
  }
  else{
    MayDay::Error("TracerParticleStepper::parseRuntimeOptions -- logic bust");
  }  
}

#ifdef CH_USE_HDF5
template <typename P>
inline
void TracerParticleStepper<P>::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const {
  CH_TIME("TracerParticleStepper::writeCheckpointData(HDF5Handle, int)");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::writeCheckpointData(HDF5Handle, int)" << endl;
  }

  m_solver->writeCheckpointLevel(a_handle, a_lvl);
}
#endif

#ifdef CH_USE_HDF5
template <typename P>
inline
void TracerParticleStepper<P>::readCheckpointData(HDF5Handle& a_handle, const int a_lvl) {
  CH_TIME("TracerParticleStepper::readCheckpointData(HDF5Handle, int)");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::readCheckpointData(HDF5Handle, int)" << endl;
  }

  m_solver->readCheckpointLevel(a_handle, a_lvl);
}
#endif

template <typename P>
inline
void TracerParticleStepper<P>::writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const {
  CH_TIME("TracerParticleStepper::writePlotData(EBAMRCellData, Vector<std::string>, int)");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::writePlotData(EBAMRCellData, Vector<std::string>, int)" << endl;
  }

  m_solver->writePlotData(a_output, a_icomp);
  a_plotVariableNames.append(m_solver->getPlotVariableNames());
}

template <typename P>
inline
int TracerParticleStepper<P>::getNumberOfPlotVariables() const {
  CH_TIME("TracerParticleStepper::getNumberOfPlotVariables()");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::getNumberOfPlotVariables()" << endl;
  }

  return m_solver->getNumberOfPlotVariables();
}

template <typename P>
inline
void TracerParticleStepper<P>::computeDt(Real& a_dt, TimeCode& a_timeCode) {
  CH_TIME("TracerParticleStepper::computeDt()");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::computeDt()" << endl;
  }
}

template <typename P>
inline
Real TracerParticleStepper<P>::advance(const Real a_dt) {
  CH_TIME("TracerParticleStepper::advance(Real)");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::advance(Real)" << endl;
  }

  switch(m_algorithm){
  case IntegrationAlgorithm::Euler:
    {
      this->advanceParticlesEuler(a_dt);

      break;
    }
  case IntegrationAlgorithm::RK2:
    {
      this->advanceParticlesRK2(a_dt);

      break;
    }
  case IntegrationAlgorithm::RK4:
    {
      this->advanceParticlesRK4(a_dt);

      break;
    }
  default:
    {
      MayDay::Error("TracerParticleStepper::advance -- logic bust");
    }
  }  

  return a_dt;
}

template <typename P>
inline
void TracerParticleStepper<P>::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("TracerParticleStepper::synchronizeSolverTimes");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::synchronizeSolverTimes" << endl;
  }
  
  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;

  m_solver->setTime(a_step, a_time, a_dt);
}

template <typename P>
inline
void TracerParticleStepper<P>::preRegrid(const int a_lmin, const int a_oldFinestLevel) {
  CH_TIME("TracerParticleStepper::preRegrid(int, int)");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::preRegrid(int, int)" << endl;
  }

  m_solver->preRegrid(a_lmin, a_oldFinestLevel);
}

template <typename P>
inline
void TracerParticleStepper<P>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) {
  CH_TIME("TracerParticleStepper::regrid(int, int, int)");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::regrid(int, int, int)" << endl;
  }

  // Define velocity field on the new mesh. 
  m_amr->reallocate(m_velocity, m_phase, a_lmin);
  DataOps::setValue(m_velocity, 0.0);
  

  // Regrid tracer particles. 
  m_solver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
}

template <typename P>
inline
void TracerParticleStepper<P>::postRegrid() {
  CH_TIME("TracerParticleStepper::postRegrid()");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::postRegrid()" << endl;
  }

  this->setVelocity();
}

template <typename P>
inline
void TracerParticleStepper<P>::setVelocity() {
  CH_TIME("TracerParticleStepper::setVelocity()");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::setVelocity()" << endl;
  }

  // Various velocity field definitions that we can use. 
  auto diagonalVelocity = [](const RealVect& a_position) -> RealVect {
    return RealVect::Unit;
  };

  DataOps::setValue(m_velocity, diagonalVelocity, m_amr->getProbLo(), m_amr->getDx());

  m_amr->averageDown(m_velocity, m_realm, m_phase);
  m_amr->interpGhost(m_velocity, m_realm, m_phase);  
}

template <typename P>
inline
void TracerParticleStepper<P>::initialParticles() {
  CH_TIME("TracerParticleStepper::initialParticles()");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::initialParticles()" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const ProblemDomain&     domain = m_amr->getDomains()[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit){
      const Box      cellBox = dbl  [dit()];
      const EBISBox& ebisbox = ebisl[dit()];
    }
  }
}

template <typename P>
inline
void TracerParticleStepper<P>::printStepReport(){
  CH_TIME("TracerParticleStepper::printStepReport");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::printStepReport()" << endl;
  }  
}

template <typename P>
inline
void TracerParticleStepper<P>::advanceParticlesEuler(const Real a_dt) {
  CH_TIME("TracerParticleStepper::advanceParticlesEuler()");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::advanceParticlesEuler()" << endl;
  }

  MayDay::Error("TracerParticleStepper::advanceParticlesEuler -- not implemented");    
}

template <typename P>
inline
void TracerParticleStepper<P>::advanceParticlesRK2(const Real a_dt) {
  CH_TIME("TracerParticleStepper::advanceParticlesRK2()");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::advanceParticlesRK2()" << endl;
  }

  MayDay::Error("TracerParticleStepper::advanceParticlesRK2 -- not implemented");
}

template <typename P>
inline
void TracerParticleStepper<P>::advanceParticlesRK4(const Real a_dt) {
  CH_TIME("TracerParticleStepper::advanceParticlesRK4()");
  if(m_verbosity > 5){
    pout() << "TracerParticleStepper::advanceParticlesRK4()" << endl;
  }

  MayDay::Error("TracerParticleStepper::advanceParticlesRK4 -- not implemented");  
}

#include <CD_NamespaceFooter.H>

#endif
