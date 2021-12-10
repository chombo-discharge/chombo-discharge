/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RtPhysicsStepperImplem.H
  @brief  Implementation of CD_RtPhysicsStepper.H
  @author Robert Marskar
*/

// Std includes
#include <math.h>

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>
#include <PolyGeom.H>

// Our includes
#include <CD_RtPhysicsStepper.H>
#include <CD_RtPhysicsSpecies.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::RadiativeTransfer;

template <typename T>
RtPhysicsStepper<T>::RtPhysicsStepper(){
  CH_TIME("RtPhysicsStepper::RtPhysicsStepper");

  // Get shit for setting the source term
  ParmParse pp("RtPhysicsStepper");
  Vector<Real> v;

  pp.get("verbosity", m_verbosity);
  pp.get("realm",     m_realm);
  pp.get("dt",        m_dt);
}
    
template <typename T>
RtPhysicsStepper<T>::~RtPhysicsStepper(){
  CH_TIME("RtPhysicsStepper::~RtPhysicsStepper");  
}

template <typename T>
void RtPhysicsStepper<T>::setupSolvers(){
  CH_TIME("RtPhysicsStepper::setupSolvers");
  
  m_species = RefCountedPtr<RtSpecies> (new RtPhysicsSpecies());

  // Solver setup
  m_solver = RefCountedPtr<RtSolver> (new T());
  m_solver->setVerbosity(m_verbosity);
  m_solver->setRtSpecies(m_species);
  m_solver->parseOptions();
  m_solver->setPhase(phase::gas);
  m_solver->setAmr(m_amr);
  m_solver->setComputationalGeometry(m_computationalGeometry);
  m_solver->setRealm(m_realm);
}

template <typename T>
void RtPhysicsStepper<T>::parseRuntimeOptions() {
  CH_TIME("RtPhysicsStepper::parseRuntimeOptions");
  
  ParmParse pp("RtPhysicsStepper");

  pp.get("verbosity", m_verbosity);
  pp.get("dt",        m_dt);
      
  m_solver->parseRuntimeOptions();
}

template <typename T>
void RtPhysicsStepper<T>::allocate() {
  CH_TIME("RtPhysicsStepper::allocate");
  
  m_solver->allocateInternals();
}

template <typename T>
void RtPhysicsStepper<T>::initialData(){
  CH_TIME("RtPhysicsStepper::initialData");
  
  m_solver->initialData();
}

template <typename T>
void RtPhysicsStepper<T>::postInitialize(){
  CH_TIME("RtPhysicsStepper::postInitialize");
  
  this->setGaussianSource();  
}

#ifdef CH_USE_HDF5      
template <typename T>
void RtPhysicsStepper<T>::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const {
  CH_TIME("RtPhysicsStepper::writeCheckpointData");
  
  m_solver->writeCheckpointLevel(a_handle, a_lvl);
}
#endif

#ifdef CH_USE_HDF5      
template <typename T>
void RtPhysicsStepper<T>::readCheckpointData(HDF5Handle& a_handle, const int a_lvl){
  CH_TIME("RtPhysicsStepper::readCheckpointData");
  
  m_solver->readCheckpointLevel(a_handle, a_lvl);
}
#endif

template <typename T>
void RtPhysicsStepper<T>::postCheckpointSetup() {
  CH_TIME("RtPhysicsStepper::postCheckpointSetup");
  
  this->setGaussianSource();
}

template <typename T>
int RtPhysicsStepper<T>::getNumberOfPlotVariables() const {
  CH_TIME("RtPhysicsStepper::getNumberOfPlotVariables");
  
  return m_solver->getNumberOfPlotVariables();
}

template <typename T>
void RtPhysicsStepper<T>::writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const {
  CH_TIME("RtPhysicsStepper::writePlotData");
  
  a_plotVariableNames.append(m_solver->getPlotVariableNames());
  m_solver->writePlotData(a_output, a_icomp);
}

template <typename T>
void RtPhysicsStepper<T>::computeDt(Real& a_dt, TimeCode& a_timeCode) {
  CH_TIME("RtPhysicsStepper::computeDt");
  
  a_dt = m_dt;
}

template <typename T>
Real RtPhysicsStepper<T>::advance(const Real a_dt) {
  CH_TIME("RtPhysicsStepper::advance");
  
  m_solver->advance(a_dt);

  return a_dt;
}

template <typename T>
void RtPhysicsStepper<T>::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("RtPhysicsStepper::synchronizeSolverTimes");
  
  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;

  m_solver->setTime(m_timeStep, m_time, m_dt);
}

template <typename T>
void RtPhysicsStepper<T>::printStepReport() {
  CH_TIME("RtPhysicsStepper::printStepReport");  
}

template <typename T>
void RtPhysicsStepper<T>::registerRealms() {
  CH_TIME("RtPhysicsStepper::registerRealms");
  
  m_amr->registerRealm(m_realm);
}

template <typename T>
void RtPhysicsStepper<T>::registerOperators() {
  CH_TIME("RtPhysicsStepper::registerOperators");
  
  m_solver->registerOperators();
}

template <typename T>
void RtPhysicsStepper<T>::preRegrid(const int a_base, const int a_oldFinestLevel) {
  CH_TIME("RtPhysicsStepper::preRegrid");
  
  m_solver->preRegrid(a_base, a_oldFinestLevel);
}

template <typename T>
void RtPhysicsStepper<T>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) {
  CH_TIME("RtPhysicsStepper::regrid");
  
  m_solver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
}

template <typename T>
void RtPhysicsStepper<T>::postRegrid(){
  CH_TIME("RtPhysicsStepper::postRegrid");  

  this->setGaussianSource();  
}

template <typename T>
void RtPhysicsStepper<T>::setGaussianSource(){
  CH_TIME("RtPhysicsStepper::setGaussianSource");
  
  ParmParse pp("RtPhysicsStepper");

  Vector<Real> v(SpaceDim);
  RealVect blobCenter;
  Real     blobAmp;
  Real     blobRad;
  
  pp.get   ("blob_amplitude", blobAmp);
  pp.get   ("blob_radius",    blobRad);
  pp.getarr("blob_center",    v, 0, SpaceDim); blobCenter = RealVect(D_DECL(v[0], v[1], v[2]));

  auto sourceFunction = [blobAmp, blobRad, blobCenter](const RealVect a_position) -> Real {
    const RealVect dist = a_position-blobCenter;
    const Real dist2    = dist.dotProduct(dist);
    return blobAmp*exp(-dist2/(2*blobRad*blobRad));
  };

  m_solver->setSource(sourceFunction);
}

#include <CD_NamespaceFooter.H>
