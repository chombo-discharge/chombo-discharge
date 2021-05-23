/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RtPhysicsStepperImplem.H
  @brief  Implementation of CD_RtPhysicsStepper.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>
#include <PolyGeom.H>

// Our includes
#include <CD_RtPhysicsStepper.H>
#include <CD_RtPhysicsSpecies.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::RadiativeTransfer;

template <typename T>
RtPhysicsStepper<T>::RtPhysicsStepper(){

  // Get shit for setting the source term
  ParmParse pp("RtPhysicsStepper");
  Vector<Real> v;

  pp.get   ("verbosity",      m_verbosity);
  pp.get   ("realm",          m_realm);
  pp.get   ("blob_amplitude", m_blob_amplitude);
  pp.get   ("blob_radius",    m_blob_radius);
  pp.getarr("blob_center",    v, 0, SpaceDim); m_blob_center = RealVect(D_DECL(v[0], v[1], v[2]));
}
    
template <typename T>
RtPhysicsStepper<T>::~RtPhysicsStepper(){}

template <typename T>
void RtPhysicsStepper<T>::setupSolvers(){
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
  ParmParse pp("RtPhysicsStepper");

  pp.get   ("verbosity", m_verbosity);
      
  m_solver->parseRuntimeOptions();
}

template <typename T>
void RtPhysicsStepper<T>::allocate() {
  m_solver->allocateInternals();
}

template <typename T>
void RtPhysicsStepper<T>::initialData(){
  m_solver->initialData();

  this->set_gaussian_source();
}

template <typename T>
void RtPhysicsStepper<T>::postInitialize(){
}

template <typename T>
void RtPhysicsStepper<T>::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const {
  m_solver->writeCheckpointLevel(a_handle, a_lvl);
}

template <typename T>
void RtPhysicsStepper<T>::readCheckpointData(HDF5Handle& a_handle, const int a_lvl){
  m_solver->readCheckpointLevel(a_handle, a_lvl);
}

template <typename T>
void RtPhysicsStepper<T>::postCheckpointSetup() {
  this->set_gaussian_source();
}

template <typename T>
int RtPhysicsStepper<T>::getNumberOfPlotVariables() const {
  return m_solver->getNumberOfPlotVariables();
}

template <typename T>
void RtPhysicsStepper<T>::writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const {
  a_plotVariableNames.append(m_solver->getPlotVariableNames());
  m_solver->writePlotData(a_output, a_icomp);
}

template <typename T>
void RtPhysicsStepper<T>::computeDt(Real& a_dt, TimeCode& a_timeCode) {
  a_dt = 1.E-10;
}

template <typename T>
Real RtPhysicsStepper<T>::advance(const Real a_dt) {
  m_solver->advance(a_dt);

  return a_dt;
}

template <typename T>
void RtPhysicsStepper<T>::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) {
  m_timeStep = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

template <typename T>
void RtPhysicsStepper<T>::printStepReport() {
}

// Regrid routines
template <typename T>
bool RtPhysicsStepper<T>::needToRegrid() {
  return false;
}

template <typename T>
void RtPhysicsStepper<T>::registerRealms() {
  m_amr->registerRealm(m_realm);
}

template <typename T>
void RtPhysicsStepper<T>::registerOperators() {
  m_solver->registerOperators();
}

template <typename T>
void RtPhysicsStepper<T>::preRegrid(const int a_base, const int a_oldFinestLevel) {
  m_solver->preRegrid(a_base, a_oldFinestLevel);
}

template <typename T>
void RtPhysicsStepper<T>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) {
  m_solver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);

  // Do the rest of the setup
  this->set_gaussian_source();
}

template <typename T>
void RtPhysicsStepper<T>::postRegrid(){
}

template <typename T>
void RtPhysicsStepper<T>::set_gaussian_source(){

  EBAMRCellData& source = m_solver->getSource();
  const RealVect prob_lo = m_amr->getProbLo();

  //
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, phase::gas)[lvl];
    const Real dx                = m_amr->getDx()[lvl];


    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      const Box box = dbl.get(dit());

      FArrayBox& fab = (*source[lvl])[dit()].getFArrayBox();
	  
      for (BoxIterator bit(box); bit.ok(); ++bit){
	const IntVect iv   = bit();
	const RealVect pos = prob_lo + dx*(RealVect(iv) + 0.5*RealVect::Unit);

	const RealVect v = (pos - m_blob_center)/m_blob_radius;
	const Real factor = PolyGeom::dot(v,v);
	const Real value = m_blob_amplitude*exp(-0.5*factor);

	fab(iv,0) = int(value);
      }
    }
  }
}

#include <CD_NamespaceFooter.H>
