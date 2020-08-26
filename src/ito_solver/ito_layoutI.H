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
const std::string ito_layout<T>::get_realm() const{
  return m_realm;
}

template <class T>
void ito_layout<T>::set_realm(const std::string a_realm){
  m_realm = a_realm;

  for (auto solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->set_realm(m_realm);
  }
}

template <class T>
void ito_layout<T>::define(const Vector<RefCountedPtr<ito_species> >& a_species){
  m_isDefined = true;
}

template <class T>
void ito_layout<T>::parse_options(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->parse_options();
  }
}

template <class T>
void ito_layout<T>::allocate_internals(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->allocate_internals();
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
void ito_layout<T>::pre_regrid(const int a_lbase, const int a_finest_level){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->pre_regrid(a_lbase, a_finest_level);
  }
}

template <class T>
void ito_layout<T>::initial_data(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->initial_data();
  }
}

template <class T>
void ito_layout<T>::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->regrid(a_lmin, a_old_finest_level, a_new_finest_level);
  }
}

template <class T>
void ito_layout<T>::register_operators(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->register_operators();
  }
}

template <class T>
void ito_layout<T>::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->set_amr(a_amr);
  }
}

template <class T>
void ito_layout<T>::set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->set_computational_geometry(a_compgeom);
  }
}

template <class T>
void ito_layout<T>::set_phase(phase::which_phase a_phase){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->set_phase(a_phase);
  }
}

template <class T>
void ito_layout<T>::set_verbosity(const int a_verbosity){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->set_verbosity(a_verbosity);
  }
}

template <class T>
void ito_layout<T>::set_time(const int a_step, const Real a_time, const Real a_dt){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->set_time(a_step, a_time, a_dt);
  }
}

template <class T>
void ito_layout<T>::sort_particles_by_cell(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->sort_particles_by_cell();
  }
}

template <class T>
void ito_layout<T>::sort_particles_by_patch(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->sort_particles_by_patch();
  }
}

template <class T>
void ito_layout<T>::make_superparticles(const int a_ppc){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->make_superparticles(a_ppc);
  }
}
template <class T>
void ito_layout<T>::update_invalid_particles(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->update_invalid_particles();
  }
}

template <class T>
void ito_layout<T>::copy_invalid_particles_to_particles(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->copy_invalid_particles_to_particles();
  }
}

template <class T>
void ito_layout<T>::move_invalid_particles_to_particles(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->move_invalid_particles_to_particles();
  }
}

template <class T>
void ito_layout<T>::deposit_particles(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->deposit_particles();
  }
}

template <class T>
void ito_layout<T>::remap(){
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    iter()->remap();
  }
}

template <class T>
size_t ito_layout<T>::get_num_particles(const bool a_local) {
  size_t ret = 0;
  for (auto iter = this->iterator(); iter.ok(); ++iter){
    ret += iter()->get_num_particles(a_local);
  }

  return ret;
}

template <class T>
size_t ito_layout<T>::get_num_eb_particles(const bool a_local) {
  size_t ret = 0;
  for (auto iter = this->iterator(); iter.ok(); ++iter){
    ret += iter()->get_num_eb_particles(a_local);
  }

  return ret;
}

template <class T>
size_t ito_layout<T>::get_num_domain_particles(const bool a_local) {
  size_t ret = 0;
  for (auto iter = this->iterator(); iter.ok(); ++iter){
    ret += iter()->get_num_domain_particles(a_local);
  }

  return ret;
}

template <class T>
size_t ito_layout<T>::get_num_source_particles(const bool a_local) {
  size_t ret = 0;
  for (auto iter = this->iterator(); iter.ok(); ++iter){
    ret += iter()->get_num_source_particles(a_local);
  }

  return ret;
}

template <class T>
Real ito_layout<T>::compute_dt(){
  Real minDt = 1.E99;

  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    const Real thisDt = iter()->compute_dt();
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
Vector<RefCountedPtr<T> >& ito_layout<T>::get_solvers(){
  return m_solvers;
}

template <class T>
Vector<RefCountedPtr<ito_species> >& ito_layout<T>::get_species(){
  return m_species;
}

template <class T>
Vector<EBAMRCellData*> ito_layout<T>::get_velocities(){
  Vector<EBAMRCellData*> ret(m_solvers.size(), nullptr);
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.get_solver()] = &(iter()->get_velo_cell());
  }

  return ret;
}

template <class T>
Vector<EBAMRCellData*> ito_layout<T>::get_densities(){
  Vector<EBAMRCellData*> ret(m_solvers.size(), nullptr);
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.get_solver()] = &(iter()->get_state());
  }

  return ret;
}

template <class T>
Vector<EBAMRCellData*> ito_layout<T>::get_diffusion(){
  Vector<EBAMRCellData*> ret(m_solvers.size(), nullptr);
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.get_solver()] = &(iter()->get_diffco_cell());
  }

  return ret;
}

template <class T>
Vector<particle_container<ito_particle>* > ito_layout<T>::get_particles(){
  Vector<particle_container<ito_particle>* > ret(m_solvers.size(), nullptr);
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.get_solver()] = &(iter()->get_particles());
  }
  
  return ret;
}

template <class T>
Vector<particle_container<ito_particle>* > ito_layout<T>::get_eb_particles(){
  
  Vector<particle_container<ito_particle>* > ret(m_solvers.size(), nullptr);
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.get_solver()] = &(iter()->get_eb_particles());
  }
  
  return ret;
}

template <class T>
Vector<particle_container<ito_particle>* > ito_layout<T>::get_domain_particles(){

  Vector<particle_container<ito_particle>* > ret(m_solvers.size(), nullptr);
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.get_solver()] = &(iter()->get_domain_particles());
  }
  
  return ret;
}

template <class T>
Vector<particle_container<ito_particle>* > ito_layout<T>::get_source_particles(){

  Vector<particle_container<ito_particle>* > ret(m_solvers.size(), nullptr);
  for (ito_iterator<T> iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.get_solver()] = &(iter()->get_source_particles());
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
    solver->set_species(a_species[i]);
    solver->set_verbosity(-1);
    ito->add_solver(solver);
  }

  return ito;
}

#endif
