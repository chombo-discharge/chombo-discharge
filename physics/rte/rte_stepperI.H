/*!
  @file   rte_stepperI.H
  @brief  Implementation of rte_stepper.H
  @author Robert Marskar
  @data   May 2020
*/

#include "rte_stepper.H"
#include "rte_physics_species.H"

#include <ParmParse.H>
#include <PolyGeom.H>

#include "CD_NamespaceHeader.H"
using namespace physics::rte;

template <typename T>
rte_stepper<T>::rte_stepper(){

  // Get shit for setting the source term
  ParmParse pp("rte_stepper");
  Vector<Real> v;

  pp.get   ("verbosity",      m_verbosity);
  pp.get   ("Realm",          m_Realm);
  pp.get   ("blob_amplitude", m_blob_amplitude);
  pp.get   ("blob_radius",    m_blob_radius);
  pp.getarr("blob_center",    v, 0, SpaceDim); m_blob_center = RealVect(D_DECL(v[0], v[1], v[2]));
}
    
template <typename T>
rte_stepper<T>::~rte_stepper(){}

template <typename T>
void rte_stepper<T>::setup_solvers(){
  m_species = RefCountedPtr<rte_species> (new rte_physics_species());

  // Solver setup
  m_solver = RefCountedPtr<rte_solver> (new T());
  m_solver->setVerbosity(m_verbosity);
  m_solver->set_rte_species(m_species);
  m_solver->parseOptions();
  m_solver->setPhase(phase::gas);
  m_solver->setAmr(m_amr);
  m_solver->setComputationalGeometry(m_computationalGeometry);
  m_solver->setRealm(m_Realm);
}

template <typename T>
void rte_stepper<T>::parseRuntimeOptions() {
  ParmParse pp("rte_stepper");

  pp.get   ("verbosity", m_verbosity);
      
  m_solver->parseRuntimeOptions();
}

template <typename T>
void rte_stepper<T>::allocate() {
  m_solver->allocateInternals();
}

template <typename T>
void rte_stepper<T>::initialData(){
  m_solver->initialData();

  this->set_gaussian_source();
}

template <typename T>
void rte_stepper<T>::postInitialize(){
}

template <typename T>
void rte_stepper<T>::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const {
  m_solver->writeCheckpointLevel(a_handle, a_lvl);
}

template <typename T>
void rte_stepper<T>::readCheckpointData(HDF5Handle& a_handle, const int a_lvl){
  m_solver->readCheckpointLevel(a_handle, a_lvl);
}

template <typename T>
void rte_stepper<T>::postCheckpointSetup() {
  this->set_gaussian_source();
}

template <typename T>
int rte_stepper<T>::getNumberOfPlotVariables() const {
  return m_solver->getNumberOfPlotVariables();
}

template <typename T>
void rte_stepper<T>::writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const {
  a_plotVariableNames.append(m_solver->get_plotVariableNames());
  m_solver->writePlotData(a_output, a_icomp);
}

template <typename T>
void rte_stepper<T>::computeDt(Real& a_dt, TimeCode& a_timeCode) {
  a_dt = 1.E-10;
}

template <typename T>
Real rte_stepper<T>::advance(const Real a_dt) {
  m_solver->advance(a_dt);

  return a_dt;
}

template <typename T>
void rte_stepper<T>::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) {
  m_timeStep = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

template <typename T>
void rte_stepper<T>::printStepReport() {
}

// Regrid routines
template <typename T>
bool rte_stepper<T>::needToRegrid() {
  return false;
}

template <typename T>
void rte_stepper<T>::registerRealms() {
  m_amr->registerRealm(m_Realm);
}

template <typename T>
void rte_stepper<T>::registerOperators() {
  m_solver->registerOperators();
}

template <typename T>
void rte_stepper<T>::preRegrid(const int a_base, const int a_oldFinestLevel) {
  m_solver->preRegrid(a_base, a_oldFinestLevel);
}

template <typename T>
void rte_stepper<T>::deallocate() {
  m_solver->deallocateInternals();
}

template <typename T>
void rte_stepper<T>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) {
  m_solver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);

  // Do the rest of the setup
  this->set_gaussian_source();
}

template <typename T>
void rte_stepper<T>::postRegrid(){
}

template <typename T>
void rte_stepper<T>::set_gaussian_source(){

  EBAMRCellData& source = m_solver->getSource();
  const RealVect prob_lo = m_amr->getProbLo();

  //
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_Realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_Realm, phase::gas)[lvl];
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
#include "CD_NamespaceFooter.H"
