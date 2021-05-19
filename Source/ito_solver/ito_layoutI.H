/*!
  @file   ito_layout.cpp
  @brief  Implementation of ito_layout.H
  @author Robert Marskar
  @date   May 2020
*/

#ifndef _ITO_LAYOUTI_H
#define _ITO_LAYOUTI_H

#include "ito_layout.H"
#include "ito_iterator.H"

#include "CD_NamespaceHeader.H"

template <class T>
ito_layout<T>::ito_layout(){
  m_isDefined = false;
}

template <class T>
ito_layout<T>::ito_layout(const Vector<RefCountedPtr<ito_species> >& a_species){
  this->define(a_species);
  m_solvers.resize(0);
}

template <class T>
ito_layout<T>::~ito_layout(){

}

template <class T>
ito_iterator<T> ito_layout<T>::iterator(){
  return ito_iterator<T>(*this);
}

template <class T>
const std::string ito_layout<T>::getRealm() const{
  return m_realm;
}

template <class T>
void ito_layout<T>::setRealm(const std::string a_realm){
  m_realm = a_realm;

  for (auto solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->setRealm(m_realm);
  }
}

template <class T>
void ito_layout<T>::define(const Vector<RefCountedPtr<ito_species> >& a_species){
  m_isDefined = true;
}

template <class T>
void ito_layout<T>::parseOptions(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->parseOptions();
  }
}

template <class T>
void ito_layout<T>::parseRuntimeOptions(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->parseRuntimeOptions();
  }
}

template <class T>
void ito_layout<T>::setPVRBuffer(const int a_buffer){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->setPVRBuffer(a_buffer);
  }
}

template <class T>
void ito_layout<T>::setHalobuffer(const int a_buffer){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->setHalobuffer(a_buffer);
  }
}

template <class T>
void ito_layout<T>::allocateInternals(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->allocateInternals();
  }
}

template <class T>
void ito_layout<T>::add_solver(RefCountedPtr<T> a_solver){
  m_solvers.push_back(a_solver);
}

template <class T>
void ito_layout<T>::interpolate_velocities(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->interpolate_velocities();
  }
}

template <class T>
void ito_layout<T>::preRegrid(const int a_lbase, const int a_finestLevel){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->preRegrid(a_lbase, a_finestLevel);
  }
}

template <class T>
void ito_layout<T>::initialData(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->initialData();
  }
}

template <class T>
void ito_layout<T>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  }
}

template <class T>
void ito_layout<T>::registerOperators(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->registerOperators();
  }
}

template <class T>
void ito_layout<T>::setAmr(const RefCountedPtr<AmrMesh>& a_amr){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->setAmr(a_amr);
  }
}

template <class T>
void ito_layout<T>::setComputationalGeometry(const RefCountedPtr<computational_geometry>& a_computationalGeometry){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->setComputationalGeometry(a_computationalGeometry);
  }
}

template <class T>
void ito_layout<T>::setPhase(phase::which_phase a_phase){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->setPhase(a_phase);
  }
}

template <class T>
void ito_layout<T>::setVerbosity(const int a_verbosity){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->setVerbosity(a_verbosity);
  }
}

template <class T>
void ito_layout<T>::setTime(const int a_step, const Real a_time, const Real a_dt){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->setTime(a_step, a_time, a_dt);
  }
}

template <class T>
void ito_layout<T>::sortParticlesByCell(const ito_solver::which_container a_container){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->sortParticlesByCell(a_container);
  }
}

template <class T>
void ito_layout<T>::sortParticlesByPatch(const ito_solver::which_container a_container){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->sortParticlesByPatch(a_container);
  }
}

template <class T>
void ito_layout<T>::make_superparticles(const ito_solver::which_container a_container, const int a_ppc){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->make_superparticles(a_container, a_ppc);
  }
}

template <class T>
void ito_layout<T>::deposit_particles(){
  this->deposit_particles(ito_solver::which_container::bulk);
}

template <class T>
void ito_layout<T>::deposit_particles(const ito_solver::which_container a_container){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->deposit_particles(a_container);
  }
}

template <class T>
void ito_layout<T>::remap(){
  this->remap(ito_solver::which_container::bulk);
}

template <class T>
void ito_layout<T>::remap(const ito_solver::which_container a_container){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->remap(a_container);
  }
}

template <class T>
size_t ito_layout<T>::get_num_particles(const ito_solver::which_container a_container, const bool a_local) {
  size_t ret = 0;
  for (auto iter = this->iterator(); iter.ok(); ++iter){
    ret += iter()->get_num_particles(a_container, a_local);
  }

  return ret;
}

template <class T>
Real ito_layout<T>::computeDt(){
  Real minDt = 1.E99;

  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    const Real thisDt = iter()->computeDt();
    minDt = Min(minDt, thisDt);
  }
  
  return minDt;
}

template <class T>
Real ito_layout<T>::compute_advective_dt(){
  Real minDt = 1.E99;

  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    const Real thisDt = iter()->compute_advective_dt();
    minDt = Min(minDt, thisDt);
  }
  
  return minDt;
}

template <class T>
Real ito_layout<T>::compute_min_dt(const Real a_maxCellsToMove){
  Real minDt = 1.E99;

  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    const Real thisDt = iter()->compute_min_dt(a_maxCellsToMove);
    minDt = Min(minDt, thisDt);
  }
  
  return minDt;
}

template <class T>
Real ito_layout<T>::compute_min_drift_dt(const Real a_maxCellsToMove){
  Real minDt = 1.E99;

  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    const Real thisDt = iter()->compute_min_drift_dt(a_maxCellsToMove);
    minDt = Min(minDt, thisDt);
  }
  
  return minDt;
}

template <class T>
Real ito_layout<T>::compute_min_diffusion_dt(const Real a_maxCellsToMove){
  Real minDt = 1.E99;

  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    const Real thisDt = iter()->compute_min_diffusion_dt(a_maxCellsToMove);
    minDt = Min(minDt, thisDt);
  }
  
  return minDt;
}

template <class T>
Real ito_layout<T>::compute_diffusive_dt() {
  Real minDt = 1.E99;

  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    const Real thisDt = iter()->compute_diffusive_dt();
    minDt = Min(minDt, thisDt);
  }
  
  return minDt;
}

template <class T>
Vector<RefCountedPtr<T> >& ito_layout<T>::get_solvers(){
  return m_solvers;
}

template <class T>
Vector<RefCountedPtr<ito_species> >& ito_layout<T>::getSpecies(){
  return m_species;
}

template <class T>
Vector<EBAMRCellData*> ito_layout<T>::get_velo_funcs(){
  Vector<EBAMRCellData*> ret(m_solvers.size(), nullptr);
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.index()] = &(iter()->get_velo_func());
  }

  return ret;
}

template <class T>
Vector<EBAMRCellData*> ito_layout<T>::get_densities(){
  Vector<EBAMRCellData*> ret(m_solvers.size(), nullptr);
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.index()] = &(iter()->getPhi());
  }

  return ret;
}

template <class T>
Vector<EBAMRCellData*> ito_layout<T>::get_scratch(){
  Vector<EBAMRCellData*> ret(m_solvers.size(), nullptr);
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.index()] = &(iter()->get_scratch());
  }

  return ret;
}

template <class T>
Vector<EBAMRCellData*> ito_layout<T>::get_diffco_funcs(){
  Vector<EBAMRCellData*> ret(m_solvers.size(), nullptr);
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.index()] = &(iter()->get_diffco_func());
  }

  return ret;
}

template <class T>
Vector<EBAMRCellData*> ito_layout<T>::get_mobility_funcs(){
  Vector<EBAMRCellData*> ret(m_solvers.size(), nullptr);
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.index()] = &(iter()->get_mobility_func());
  }

  return ret;
}

template <class T>
Vector<ParticleContainer<ito_particle>* > ito_layout<T>::getParticles(const ito_solver::which_container a_container){
  Vector<ParticleContainer<ito_particle>* > ret(m_solvers.size(), nullptr);
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.index()] = &(iter()->getParticles(a_container));
  }
  
  return ret;
}

template <class T>
phase::which_phase ito_layout<T>::get_phase() const {
  return m_phase;
}

template <class T, class S>
ito_factory<T, S>::ito_factory(){

}

template <class T, class S>
ito_factory<T, S>::~ito_factory(){

}

template <class T, class S>
RefCountedPtr<ito_layout<T> > ito_factory<T, S>::new_layout(const Vector<RefCountedPtr<ito_species> > a_species) const{

  auto ito = RefCountedPtr<ito_layout<T> > (new ito_layout<T>(a_species));

  // Cast solvers and instantiate them
  for (int i = 0; i < a_species.size(); i++){
    RefCountedPtr<T> solver = RefCountedPtr<T> (static_cast<T*> (new S()));
    solver->setSpecies(a_species[i]);
    solver->setVerbosity(-1);
    ito->add_solver(solver);
  }

  return ito;
}
#include "CD_NamespaceFooter.H"

#endif
