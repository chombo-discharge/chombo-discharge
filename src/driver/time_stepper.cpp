/*!
  @file   time_stepper.cpp
  @brief  Implementation of time_stepper.H
  @author Robert Marskar
  @date   March 2020
*/

#include "time_stepper.H"

void time_stepper::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  CH_TIME("time_stepper::set_amr");
  if(m_verbosity > 5){
    pout() << "time_stepper::set_amr" << endl;
  }

  m_amr = a_amr;
}

void time_stepper::set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom){
  CH_TIME("time_stepper::set_computational_geometry");
  if(m_verbosity > 5){
    pout() << "time_stepper::set_computational_geometry" << endl;
  }

  m_compgeom = a_compgeom;
}

bool time_stepper::need_to_regrid(){
  CH_TIME("time_stepper::need_to_regrid");
  if(m_verbosity > 5){
    pout() << "time_stepper::need_to_regrid" << endl;
  }

  return false;
}

int time_stepper::get_redistribution_regsize() const {
  CH_TIME("time_stepper::get_redistribution_regsize");
  if(m_verbosity > 5){
    pout() << "time_stepper::get_redistribution_regsize" << endl;
  }

  return 1;
}
