/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_ItoLayoutImplem.H
  @brief  Implementation of CD_ItoLayout.H
  @author Robert Marskar
*/

#ifndef CD_ItoLayoutImplem_H
#define CD_ItoLayoutImplem_H

// Our includes
#include <CD_ItoLayout.H>
#include "ito_iterator.H"
#include <CD_NamespaceHeader.H>

template <class T>
ItoLayout<T>::ItoLayout(){
  m_isDefined = false;
}

template <class T>
ItoLayout<T>::ItoLayout(const Vector<RefCountedPtr<ItoSpecies> >& a_species){
  this->define(a_species);
  m_solvers.resize(0);
}

template <class T>
ItoLayout<T>::~ItoLayout(){

}

template <class T>
ito_iterator<T> ItoLayout<T>::iterator(){
  return ito_iterator<T>(*this);
}

template <class T>
const std::string ItoLayout<T>::getRealm() const{
  return m_realm;
}

template <class T>
void ItoLayout<T>::setRealm(const std::string a_realm){
  m_realm = a_realm;

  for (auto solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->setRealm(m_realm);
  }
}

template <class T>
void ItoLayout<T>::define(const Vector<RefCountedPtr<ItoSpecies> >& a_species){
  m_isDefined = true;
}

template <class T>
void ItoLayout<T>::parseOptions(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->parseOptions();
  }
}

template <class T>
void ItoLayout<T>::parseRuntimeOptions(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->parseRuntimeOptions();
  }
}

template <class T>
void ItoLayout<T>::setPVRBuffer(const int a_buffer){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->setPVRBuffer(a_buffer);
  }
}

template <class T>
void ItoLayout<T>::setHalobuffer(const int a_buffer){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->setHalobuffer(a_buffer);
  }
}

template <class T>
void ItoLayout<T>::allocateInternals(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->allocateInternals();
  }
}

template <class T>
void ItoLayout<T>::addSolver(RefCountedPtr<T> a_solver){
  m_solvers.push_back(a_solver);
}

template <class T>
void ItoLayout<T>::interpolateVelocities(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->interpolateVelocities();
  }
}

template <class T>
void ItoLayout<T>::preRegrid(const int a_lbase, const int a_finestLevel){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->preRegrid(a_lbase, a_finestLevel);
  }
}

template <class T>
void ItoLayout<T>::initialData(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->initialData();
  }
}

template <class T>
void ItoLayout<T>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  }
}

template <class T>
void ItoLayout<T>::registerOperators(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->registerOperators();
  }
}

template <class T>
void ItoLayout<T>::setAmr(const RefCountedPtr<AmrMesh>& a_amr){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->setAmr(a_amr);
  }
}

template <class T>
void ItoLayout<T>::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->setComputationalGeometry(a_computationalGeometry);
  }
}

template <class T>
void ItoLayout<T>::setPhase(phase::which_phase a_phase){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->setPhase(a_phase);
  }
}

template <class T>
void ItoLayout<T>::setVerbosity(const int a_verbosity){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->setVerbosity(a_verbosity);
  }
}

template <class T>
void ItoLayout<T>::setTime(const int a_step, const Real a_time, const Real a_dt){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->setTime(a_step, a_time, a_dt);
  }
}

template <class T>
void ItoLayout<T>::sortParticlesByCell(const ItoSolver::WhichContainer a_container){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->sortParticlesByCell(a_container);
  }
}

template <class T>
void ItoLayout<T>::sortParticlesByPatch(const ItoSolver::WhichContainer a_container){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->sortParticlesByPatch(a_container);
  }
}

template <class T>
void ItoLayout<T>::makeSuperparticles(const ItoSolver::WhichContainer a_container, const int a_ppc){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->makeSuperparticles(a_container, a_ppc);
  }
}

template <class T>
void ItoLayout<T>::depositParticles(){
  this->depositParticles(ItoSolver::WhichContainer::bulk);
}

template <class T>
void ItoLayout<T>::depositParticles(const ItoSolver::WhichContainer a_container){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->depositParticles(a_container);
  }
}

template <class T>
void ItoLayout<T>::remap(){
  this->remap(ItoSolver::WhichContainer::bulk);
}

template <class T>
void ItoLayout<T>::remap(const ItoSolver::WhichContainer a_container){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->remap(a_container);
  }
}

template <class T>
size_t ItoLayout<T>::getNumParticles(const ItoSolver::WhichContainer a_container, const bool a_local) {
  size_t ret = 0;
  for (auto iter = this->iterator(); iter.ok(); ++iter){
    ret += iter()->getNumParticles(a_container, a_local);
  }

  return ret;
}

template <class T>
Real ItoLayout<T>::computeDt(){
  Real minDt = 1.E99;

  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    const Real thisDt = iter()->computeDt();
    minDt = Min(minDt, thisDt);
  }
  
  return minDt;
}

template <class T>
Real ItoLayout<T>::computeAdvectiveDt(){
  Real minDt = 1.E99;

  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    const Real thisDt = iter()->computeAdvectiveDt();
    minDt = Min(minDt, thisDt);
  }
  
  return minDt;
}

template <class T>
Real ItoLayout<T>::computeMinDt(const Real a_maxCellsToMove){
  Real minDt = 1.E99;

  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    const Real thisDt = iter()->computeMinDt(a_maxCellsToMove);
    minDt = Min(minDt, thisDt);
  }
  
  return minDt;
}

template <class T>
Real ItoLayout<T>::computeMinDriftDt(const Real a_maxCellsToMove){
  Real minDt = 1.E99;

  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    const Real thisDt = iter()->computeMinDriftDt(a_maxCellsToMove);
    minDt = Min(minDt, thisDt);
  }
  
  return minDt;
}

template <class T>
Real ItoLayout<T>::computeMinDiffusionDt(const Real a_maxCellsToMove){
  Real minDt = 1.E99;

  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    const Real thisDt = iter()->computeMinDiffusionDt(a_maxCellsToMove);
    minDt = Min(minDt, thisDt);
  }
  
  return minDt;
}

template <class T>
Real ItoLayout<T>::computeDiffusiveDt() {
  Real minDt = 1.E99;

  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    const Real thisDt = iter()->computeDiffusiveDt();
    minDt = Min(minDt, thisDt);
  }
  
  return minDt;
}

template <class T>
Vector<RefCountedPtr<T> >& ItoLayout<T>::getSolvers(){
  return m_solvers;
}

template <class T>
Vector<RefCountedPtr<ItoSpecies> >& ItoLayout<T>::getSpecies(){
  return m_species;
}

template <class T>
Vector<EBAMRCellData*> ItoLayout<T>::getVelocityFunctions(){
  Vector<EBAMRCellData*> ret(m_solvers.size(), nullptr);
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.index()] = &(iter()->getVelocityFunction());
  }

  return ret;
}

template <class T>
Vector<EBAMRCellData*> ItoLayout<T>::getDensities(){
  Vector<EBAMRCellData*> ret(m_solvers.size(), nullptr);
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.index()] = &(iter()->getPhi());
  }

  return ret;
}

template <class T>
Vector<EBAMRCellData*> ItoLayout<T>::getScratch(){
  Vector<EBAMRCellData*> ret(m_solvers.size(), nullptr);
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.index()] = &(iter()->getScratch());
  }

  return ret;
}

template <class T>
Vector<EBAMRCellData*> ItoLayout<T>::getDiffusionFunctions(){
  Vector<EBAMRCellData*> ret(m_solvers.size(), nullptr);
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.index()] = &(iter()->getDiffusionFunction());
  }

  return ret;
}

template <class T>
Vector<EBAMRCellData*> ItoLayout<T>::getMobilityFunctions(){
  Vector<EBAMRCellData*> ret(m_solvers.size(), nullptr);
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.index()] = &(iter()->getMobilityFunction());
  }

  return ret;
}

template <class T>
Vector<ParticleContainer<ItoParticle>* > ItoLayout<T>::getParticles(const ItoSolver::WhichContainer a_container){
  Vector<ParticleContainer<ItoParticle>* > ret(m_solvers.size(), nullptr);
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.index()] = &(iter()->getParticles(a_container));
  }
  
  return ret;
}

template <class T>
phase::which_phase ItoLayout<T>::getPhase() const {
  return m_phase;
}

template <class T, class S>
itoFactory<T, S>::itoFactory(){

}

template <class T, class S>
itoFactory<T, S>::~itoFactory(){

}

template <class T, class S>
RefCountedPtr<ItoLayout<T> > itoFactory<T, S>::newLayout(const Vector<RefCountedPtr<ItoSpecies> > a_species) const{

  auto ito = RefCountedPtr<ItoLayout<T> > (new ItoLayout<T>(a_species));

  // Cast solvers and instantiate them
  for (int i = 0; i < a_species.size(); i++){
    RefCountedPtr<T> solver = RefCountedPtr<T> (static_cast<T*> (new S()));
    solver->setSpecies(a_species[i]);
    solver->setVerbosity(-1);
    ito->addSolver(solver);
  }

  return ito;
}

#include <CD_NamespaceFooter.H>

#endif
