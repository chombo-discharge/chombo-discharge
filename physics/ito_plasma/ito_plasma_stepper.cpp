/*!
  @file   ito_plasma_stepper.cpp
  @brief  Implementation of ito_plasma_stepper.H
  @author Robert Marskar
  @date   May 2020
*/

#include "ito_plasma_stepper.H"

using namespace physics::ito_plasma;

ito_plasma_stepper::ito_plasma_stepper(){
  m_verbosity = -1;
  m_name      = "ito_plasma_stepper";
  m_phase     = phase::gas;
}

ito_plasma_stepper::ito_plasma_stepper(RefCountedPtr<ito_plasma_physics>& a_physics){
  m_physics = a_physics;
}

ito_plasma_stepper::~ito_plasma_stepper(){
}

void ito_plasma_stepper::set_verbosity(const int a_verbosity){
  CH_TIME("ito_plasma_stepper::set_verbosity");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::set_verbosity" << endl;
  }
  m_verbosity = a_verbosity;
}

void ito_plasma_stepper::setup_solvers(){
  CH_TIME("ito_plasma_stepper::setup_solver");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setup_solvers" << endl;
  }
}

void ito_plasma_stepper::initial_data(){
  CH_TIME("ito_plasma_stepper::initial_data");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::initial_data" << endl;
  }
}

void ito_plasma_stepper::post_checkpoint_setup(){
  CH_TIME("ito_plasma_stepper::post_checkpoint_setup");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::post_checkpoint_setup" << endl;
  }
}

void ito_plasma_stepper::write_checkpoint_data(HDF5Handle& a_handle, const int a_lvl) const {
  CH_TIME("ito_plasma_stepper::write_checkpoint_data");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::write_checkpoint_data" << endl;
  }
}

void ito_plasma_stepper::read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl){
  CH_TIME("ito_plasma_stepper::read_checkpoint_data");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::read_checkpoint_data" << endl;
  }
}

void ito_plasma_stepper::write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const {
  CH_TIME("ito_plasma_stepper::write_plot_data");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::write_plot_data" << endl;
  }
}

void ito_plasma_stepper::synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt){
  CH_TIME("ito_plasma_stepper::synchronize_solver_times");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::synchronize_solver_times" << endl;
  }
}

void ito_plasma_stepper::print_step_report(){
  CH_TIME("ito_plasma_stepper::print_step_report");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::print_step_report" << endl;
  }
}

void ito_plasma_stepper::compute_dt(Real& a_dt, time_code::which_code& a_timecode){
  CH_TIME("ito_plasma_stepper::compute_dt");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_dt" << endl;
  }
}

void ito_plasma_stepper::register_operators(){
  CH_TIME("ito_plasma_stepper::register_operators");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::register_operators" << endl;
  }

  m_ito->register_operators();
  m_poisson->register_operators();
  m_rte->register_operators();
  m_sigma->register_operators();
}
  
void ito_plasma_stepper::pre_regrid(const int a_lmin, const int a_old_finest_level){
  CH_TIME("ito_plasma_stepper::pre_regrid");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::pre_regrid" << endl;
  }
}

void ito_plasma_stepper::deallocate(){
  CH_TIME("ito_plasma_stepper::deallocate");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::deallocate" << endl;
  }
}

void ito_plasma_stepper::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
  CH_TIME("ito_plasma_stepper::regrid");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::regrid" << endl;
  }
}

int  ito_plasma_stepper::get_num_plot_vars() const {
  CH_TIME("ito_plasma_stepper::get_num_plot_vars");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::get_num_plot_vars" << endl;
  }
}

void ito_plasma_stepper::set_ito(RefCountedPtr<ito_layout<ito_solver> >& a_ito){
  CH_TIME("ito_plasma_stepper::set_ito");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::set_ito" << endl;
  }

  m_ito = a_ito;
}

void ito_plasma_stepper::set_poisson(RefCountedPtr<poisson_solver>& a_poisson){
  CH_TIME("ito_plasma_stepper::set_poisson");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::set_poisson" << endl;
  }

  m_poisson = a_poisson;
}

void ito_plasma_stepper::set_rte(RefCountedPtr<rte_layout<mc_photo> >& a_rte){
  CH_TIME("ito_plasma_stepper::set_rte");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::set_rte" << endl;
  }
  
  m_rte = a_rte;
}

void ito_plasma_stepper::set_potential(Real (*a_potential)(const Real a_time)){
  CH_TIME("ito_plasma_stepper::set_potential");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::set_potential" << endl;
  }

  m_potential = a_potential;
}

Real ito_plasma_stepper::get_time() const{
  CH_TIME("ito_plasma_stepper::get_time");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::get_time" << endl;
  }

  return m_time;
}

void ito_plasma_stepper::compute_E(MFAMRCellData& a_E, const MFAMRCellData& a_potential){
  CH_TIME("ito_plasma_stepper::compute_E(mfamrcell,mfamrcell)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_E(mfamrcell, mfamrcell" << endl;
  }
}

void ito_plasma_stepper::compute_E(EBAMRCellData& a_E, const phase::which_phase a_phase){
  CH_TIME("ito_plasma_stepper::compute_E(ebamrcell, phase)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_E(ebamrcell, phase" << endl;
  }
}

void ito_plasma_stepper::compute_E(EBAMRCellData& a_E, const phase::which_phase a_phase, const MFAMRCellData& a_potential){
  CH_TIME("ito_plasma_stepper::compute_E(ebamrcell, phase, mfamrcell)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_E(ebamrcell, phase mfamrcell" << endl;
  }
}

void ito_plasma_stepper::compute_E(EBAMRFluxData& a_E_face, const phase::which_phase a_phase, const EBAMRCellData& a_E_cell){
  CH_TIME("ito_plasma_stepper::compute_E(ebamrflux, phase, mfamrcell)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_E(ebamrflux, phase mfamrcell" << endl;
  }
}

void ito_plasma_stepper::compute_E(EBAMRIVData& a_E_eb,  const phase::which_phase a_phase, const EBAMRCellData& a_E_cell){
  CH_TIME("ito_plasma_stepper::compute_E(ebamriv, phase, ebamrcell)");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::compute_E(ebamriv, phase ebamrcell " << endl;
  }
}
