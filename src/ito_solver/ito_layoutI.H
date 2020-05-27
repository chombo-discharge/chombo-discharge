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
ito_iterator<T> ito_layout<T>::iterator(const species_iteration::which_species a_mode){
  return ito_iterator<T>(*this, a_mode);
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
