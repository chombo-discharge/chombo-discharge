/*!
  @file   pic_solver.cpp
  @brief  Implementation of pic_solver.H
  @author Robert Marskar
  @date   June 2019
*/

#include "pic_solver.H"

  
pic_solver::pic_solver(){

}
pic_solver::~pic_solver(){

}

void pic_solver::advance(const Real a_dt){

}

void pic_solver::deposit(EBAMRCellData& a_density) const {

}

void pic_solver::interpolate_force(const EBAMRCellData& a_E){

}

void pic_solver::set_deposition_type(){

}

void pic_solver::set_bisect_step(){

}

void pic_solver::insert_particles(EBAMRParticles& a_particles){

}

void pic_solver::allocate_internals(){

}

void pic_solver::cache_state(){

}

void pic_solver::set_time(const int a_step, const Real a_time, const Real a_dt){
  m_step = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

void pic_solver::deallocate_internals(){

}

void pic_solver::regrid(const int a_old_finest_level, const int a_new_finest_level){

}

void pic_solver::write_plot_file() const {

}

Real pic_solver::get_time() const{
  return m_time;
}
  
int pic_solver::get_step() const{
  return m_step;
}

EBAMRParticles& pic_solver::get_particles(){
  return m_particles;
  
}

EBAMRParticles& pic_solver::get_eb_particles(){
  MayDay::Abort("pic_solver::get_eb_particles - not returning the correct particles. Routine is not done.");
  return m_particles;
}

EBAMRParticles& pic_solver::get_domain_particles(){
  MayDay::Abort("pic_solver::get_domain_particles - not returning the correct particles. Routine is not done.");
  return m_particles;
}
EBAMRParticles& pic_solver::get_coll_particles(){
  MayDay::Abort("pic_solver::get_coll_particles - not returning the correct particles. Routine is not done.");
  return m_particles;
}
