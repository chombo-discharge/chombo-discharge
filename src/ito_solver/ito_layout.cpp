/*!
  @file   ito_layout.cpp
  @brief  Implementation of ito_layout.H
  @author Robert Marskar
  @date   May 2020
*/

#include "ito_layout.H"
#include "ito_iterator.H"

ito_layout::ito_layout(){
  m_isDefined = false;
}

ito_layout::ito_layout(const Vector<RefCountedPtr<ito_species> >& a_species){
  this->define(a_species);
  m_solvers.resize(0);
}

ito_layout::~ito_layout(){

}

ito_iterator ito_layout::iterator(const species_iteration::which_species a_mode){
  return ito_iterator(*this, a_mode);
}

void ito_layout::define(const Vector<RefCountedPtr<ito_species> >& a_species){
  m_isDefined = true;
}

void ito_layout::parse_options(){
  for (ito_iterator iter = this->iterator(); iter.ok(); ++iter){
    iter()->parse_options();
  }
}

void ito_layout::allocate_internals(){
  for (ito_iterator iter = this->iterator(); iter.ok(); ++iter){
    iter()->allocate_internals();
  }
}

void ito_layout::deallocate_internals(){
  for (ito_iterator iter = this->iterator(); iter.ok(); ++iter){
    //    iter()->deallocate_internals();
  }
}

void ito_layout::pre_regrid(const int a_lbase, const int a_finest_level){
  for (ito_iterator iter = this->iterator(); iter.ok(); ++iter){
    iter()->pre_regrid(a_lbase, a_finest_level);
  }
}

void ito_layout::initial_data(){
  for (ito_iterator iter = this->iterator(); iter.ok(); ++iter){
    iter()->initial_data();
  }
}

void ito_layout::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
  for (ito_iterator iter = this->iterator(); iter.ok(); ++iter){
    iter()->regrid(a_lmin, a_old_finest_level, a_new_finest_level);
  }
}

void ito_layout::register_operators(){
  for (ito_iterator iter = this->iterator(); iter.ok(); ++iter){
    iter()->register_operators();
  }
}

void ito_layout::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  for (ito_iterator iter = this->iterator(); iter.ok(); ++iter){
    iter()->set_amr(a_amr);
  }
}

void ito_layout::set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom){
  for (ito_iterator iter = this->iterator(); iter.ok(); ++iter){
    iter()->set_computational_geometry(a_compgeom);
  }
}

void ito_layout::set_phase(phase::which_phase a_phase = phase::gas){
  for (ito_iterator iter = this->iterator(); iter.ok(); ++iter){
    iter()->set_phase(a_phase);
  }
}

void ito_layout::set_verbosity(const int a_verbosity){
  for (ito_iterator iter = this->iterator(); iter.ok(); ++iter){
    iter()->set_verbosity(a_verbosity);
  }
}

void ito_layout::set_time(const int a_step, const Real a_time, const Real a_dt){
  for (ito_iterator iter = this->iterator(); iter.ok(); ++iter){
    iter()->set_time(a_step, a_time, a_dt);
  }
}

void ito_layout::write_plot_file(){
  for (ito_iterator iter = this->iterator(); iter.ok(); ++iter){
    //    iter()->write_plot_file();
  }
}

Real ito_layout::compute_min_drift_dt(const Real a_maxCellsToMove){
  Real minDt = 1.E99;

  for (ito_iterator iter = this->iterator(); iter.ok(); ++iter){
    const Real thisDt = iter()->compute_min_drift_dt(a_maxCellsToMove);
    minDt = Min(minDt, thisDt);
  }
  
  return minDt;
}

Real ito_layout::compute_min_diffusion_dt(const Real a_maxCellsToMove){
  Real minDt = 1.E99;

  for (ito_iterator iter = this->iterator(); iter.ok(); ++iter){
    const Real thisDt = iter()->compute_min_diffusion_dt(a_maxCellsToMove);
    minDt = Min(minDt, thisDt);
  }
  
  return minDt;
}

Vector<RefCountedPtr<ito_solver> >& ito_layout::get_solvers(){
  return m_solvers;
}

Vector<RefCountedPtr<ito_species> >& ito_layout::get_species(){
  return m_species;
}

Vector<particle_container<ito_particle>* > ito_layout::get_particles(){
  Vector<particle_container<ito_particle>* > ret(m_solvers.size(), nullptr);
  for (ito_iterator iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.get_solver()] = &(iter()->get_particles());
  }
  
  return ret;
}

Vector<particle_container<ito_particle>* > ito_layout::get_eb_particles(){
  
  Vector<particle_container<ito_particle>* > ret(m_solvers.size(), nullptr);
  for (ito_iterator iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.get_solver()] = &(iter()->get_eb_particles());
  }
  
  return ret;
}

Vector<particle_container<ito_particle>* > ito_layout::get_domain_particles(){

  Vector<particle_container<ito_particle>* > ret(m_solvers.size(), nullptr);
  for (ito_iterator iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.get_solver()] = &(iter()->get_domain_particles());
  }
  
  return ret;
}

Vector<particle_container<ito_particle>* > ito_layout::get_source_particles(){

  Vector<particle_container<ito_particle>* > ret(m_solvers.size(), nullptr);
  for (ito_iterator iter = this->iterator(); iter.ok(); ++iter){
    ret[iter.get_solver()] = &(iter()->get_source_particles());
  }
  
  return ret;
}

phase::which_phase ito_layout::get_phase() const {
  return m_phase;
}
