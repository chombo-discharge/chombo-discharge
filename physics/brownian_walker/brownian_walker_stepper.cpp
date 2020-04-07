/*!
  @file   brownian_walker_stepper.cpp
  @brief  Implementation of brownian_walker_stepper.H
  @author Robert Marskar
  @data   March 2020
*/

#include "brownian_walker_stepper.H"
#include "brownian_walker_species.H"

#include <ParmParse.H>

using namespace physics::brownian_walker;

brownian_walker_stepper::brownian_walker_stepper(){
  ParmParse pp("advection_diffusion");

  pp.get("diffco",   m_diffco);
  pp.get("omega",    m_omega);
  pp.get("cfl",      m_cfl);
}

brownian_walker_stepper::brownian_walker_stepper(RefCountedPtr<ito_solver>& a_solver){
  m_solver = a_solver;
}

brownian_walker_stepper::~brownian_walker_stepper(){

}

void brownian_walker_stepper::initial_data(){
  CH_TIME("brownian_walker_stepper::initial_data");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::initial_data" << endl;
  }
  
  m_solver->initial_data();

  if(m_solver->is_diffusive()){
    MayDay::Warning("yes, is diffusive");
  }
  if(m_solver->is_mobile()){
    MayDay::Warning("yes, is mobile");
  }
}

void brownian_walker_stepper::write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const{
  CH_TIME("brownian_walker_stepper::write_checkpoint_data");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::write_checkpoint_data" << endl;
  }

  m_solver->write_checkpoint_level(a_handle, a_lvl);
}

void brownian_walker_stepper::read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) {
  CH_TIME("brownian_walker_stepper::read_checkpoint_data");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::read_checkpoint_data" << endl;
  }
  
  m_solver->read_checkpoint_level(a_handle, a_lvl);
}

void brownian_walker_stepper::post_checkpoint_setup() {
  CH_TIME("brownian_walker_stepper::post_checkpoint_setup");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::post_checkpoint_setup" << endl;
  }

  MayDay::Abort("brownian_walker_stepper::post_checkpoint_setup - not implemented");
}

int brownian_walker_stepper::get_num_plot_vars() const {
  CH_TIME("brownian_walker_stepper::get_num_plot_vars");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::get_num_plot_vars" << endl;
  }

  return m_solver->get_num_plotvars();
}

void brownian_walker_stepper::write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const {
  CH_TIME("brownian_walker_stepper::write_plot_data");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::write_plot_data" << endl;
  }
  a_plotvar_names.append(m_solver->get_plotvar_names());
  m_solver->write_plot_data(a_output, a_icomp);
}

void brownian_walker_stepper::compute_dt(Real& a_dt, time_code::which_code& a_timecode) {
  CH_TIME("brownian_walker_stepper::compute_dt");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::compute_dt" << endl;
  }
  
  MayDay::Abort("brownian_walker_stepper::compute_dt - not implemented yet");
}

void brownian_walker_stepper::synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("brownian_walker_stepper::synchronize_solver_times");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::synchronize_solver_times" << endl;
  }
  
  m_solver->set_time(a_step, a_time, a_dt);

  m_step = a_step;
  m_time = a_time;
  m_dt   = a_dt;
}

void brownian_walker_stepper::print_step_report() {
  CH_TIME("brownian_walker_stepper::print_step_report");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::print_step_report" << endl;
  }

  // Do nothing
}

bool brownian_walker_stepper::need_to_regrid() {
  return false;
}

void brownian_walker_stepper::cache() {

}

void brownian_walker_stepper::deallocate() {
  CH_TIME("brownian_walker_stepper::deallocate");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::deallocate" << endl;
  }
}

void brownian_walker_stepper::setup_solvers() {
  CH_TIME("brownian_walker_stepper::setup_solvers");
  if(m_verbosity > 5){
    pout() << "brownian_walker_stepper::setup_solvers" << endl;
  }

  m_species = RefCountedPtr<ito_species> (new brownian_walker_species());

  m_solver->set_verbosity(-1);
  m_solver->parse_options();
  m_solver->set_amr(m_amr);
  m_solver->set_species(m_species);
  m_solver->set_amr(m_amr);
  m_solver->set_phase(phase::gas);
  m_solver->set_computational_geometry(m_compgeom);
  m_solver->allocate_internals();
}

Real brownian_walker_stepper::advance(const Real a_dt) {

}

void brownian_walker_stepper::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level) {

}
