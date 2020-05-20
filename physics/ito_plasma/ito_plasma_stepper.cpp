/*!
  @file   ito_plasma_stepper.cpp
  @brief  Implementation of ito_plasma_stepper.H
  @author Robert Marskar
  @date   May 2020
*/

#include "ito_plasma_stepper.H"

using namespace physics::ito_plasma;

ito_plasma_stepper::ito_plasma_stepper(){
  m_name = "ito_plasma_stepper";
}

ito_plasma_stepper::~ito_plasma_stepper(){
}

void ito_plasma_stepper::setup_solvers(){
}

void ito_plasma_stepper::initial_data(){
}

void ito_plasma_stepper::post_checkpoint_setup(){
}

bool ito_plama_stepper::need_to_regrid(){
  CH_TIME("ito_plasma::need_to_regrid");
  if(m_verbosity > 5){
    pout() << m_name + "::need_to_regrid" << endl;
  }

  return false;
}

void ito_plasma_stepper::write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl){}
void ito_plasma_stepper::read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl){}
void ito_plasma_stepper::write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const {}
void ito_plasma_stepper::synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt){}
void ito_plasma_stepper::print_step_report(){}
void ito_plasma_stepper::compute_dt(Real& a_dt, time_code::which_code& a_timecode){}
void ito_plasma_stepper::register_operators(){
} 
void ito_plasma_stepper::pre_regrid(const int a_lmin, const int a_old_finest_level){
}
void ito_plasma_stepper::deallocate(){
}
void ito_plasma_stepper::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
}

int  ito_plasma_stepper::get_num_plot_vars() const {}

Real ito_plasma_stepper::advance(const Real a_dt){}



