/*!
  @file   ito_plasma_stepper.cpp
  @brief  Implementation of ito_plasma_stepper.H
  @author Robert Marskar
  @date   May 2020
*/

#include "ito_plasma_stepper.H"

using namespace physics::ito_plasma;

Real ito_plasma_stepper::s_constant_one(const RealVect a_pos){
  return 1.0;
}

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

  // Parse class options
  this->parse_options();

  // Set up solvers
  this->setup_ito();
  this->setup_poisson();
  this->setup_rte();
  this->setup_sigma();


  // Allocate internal stuff
  this->allocate_internals();
}

void ito_plasma_stepper::setup_ito(){
  CH_TIME("ito_plasma_stepper::setup_ito");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setup_ito" << endl;
  }

  m_ito->set_verbosity(m_verbosity);
  m_ito->parse_options();
  m_ito->set_amr(m_amr);
  m_ito->set_phase(m_phase);
  m_ito->set_computational_geometry(m_compgeom);
  m_ito->allocate_internals(); // Allocate some internal storage
}

void ito_plasma_stepper::setup_poisson(){
  CH_TIME("ito_plasma_stepper::setup_poisson");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setup_poisson" << endl;
  }

  m_poisson->set_verbosity(m_verbosity);
  m_poisson->parse_options();
  m_poisson->set_amr(m_amr);
  m_poisson->set_computational_geometry(m_compgeom);

  m_poisson->set_poisson_wall_func(0, Side::Lo, m_wall_func_x_lo); // Set function-based Poisson on xlo
  m_poisson->set_poisson_wall_func(0, Side::Hi, m_wall_func_x_hi); // Set function-based Poisson on xhi
  m_poisson->set_poisson_wall_func(1, Side::Lo, m_wall_func_y_lo); // Set function-based Poisson on ylo
  m_poisson->set_poisson_wall_func(1, Side::Hi, m_wall_func_y_hi); // Set function-based Poisson on yhi
#if CH_SPACEDIM==3
  m_poisson->set_poisson_wall_func(2, Side::Lo, m_wall_func_z_lo); // Set function-based Poisson on zlo
  m_poisson->set_poisson_wall_func(2, Side::Hi, m_wall_func_z_hi); // Set function-based Poisson on zhi
#endif
  m_poisson->set_potential(m_potential); // Needs to happen AFTER set_poisson_wall_func

  m_poisson->sanity_check();
  m_poisson->allocate_internals();
}

void ito_plasma_stepper::setup_rte(){
  CH_TIME("ito_plasma_stepper::setup_rte");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setup_rte" << endl;
  }

  m_rte->set_verbosity(m_verbosity);
  m_rte->parse_options();
  m_rte->set_phase(m_phase);
  m_rte->set_amr(m_amr);
  m_rte->set_computational_geometry(m_compgeom);
  m_rte->sanity_check();
  m_rte->allocate_internals();
}

void ito_plasma_stepper::setup_sigma(){
  CH_TIME("ito_plasma_stepper::setup_sigma");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::setup_sigma" << endl;
  }

  m_sigma = RefCountedPtr<sigma_solver> (new sigma_solver());
  m_sigma->set_amr(m_amr);
  m_sigma->set_verbosity(m_verbosity);
  m_sigma->set_computational_geometry(m_compgeom);
  m_sigma->allocate_internals();
}

void ito_plasma_stepper::set_poisson_wall_func(const int a_dir, const Side::LoHiSide a_side, Real (*a_func)(const RealVect a_pos)){
  CH_TIME("ito_plasma_stepper::set_poisson_wall_func(dir, side, func)");
  if(m_verbosity > 4){
    pout() << "ito_plasma_stepper::set_poisson_wall_func(dir, side, func)" << endl;
  }

  if(a_dir == 0){
    if(a_side == Side::Lo){
      m_wall_func_x_lo = a_func;
    }
    else if(a_side == Side::Hi){
      m_wall_func_x_hi = a_func;
    }
  }
  else if(a_dir == 1){
    if(a_side == Side::Lo){
      m_wall_func_y_lo = a_func;
    }
    else if(a_side == Side::Hi){
      m_wall_func_y_hi = a_func;
    }
  }
#if CH_SPACEDIM==3
  else if(a_dir == 2){
    if(a_side == Side::Lo){
      m_wall_func_z_lo = a_func;
    }
    else if(a_side == Side::Hi){
      m_wall_func_z_hi = a_func;
    }
  }
#endif
}

void ito_plasma_stepper::set_poisson_wall_func(Real (*a_func)(const RealVect a_pos)){
  CH_TIME("ito_plasma_stepper::set_poisson_wall_func(func)");
  if(m_verbosity > 4){
    pout() << "ito_plasma_stepper::set_poisson_wall_func(func)" << endl;
  }

  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      this->set_poisson_wall_func(dir, sit(), a_func);
    }
  }
}

void ito_plasma_stepper::initial_data(){
  CH_TIME("ito_plasma_stepper::initial_data");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::initial_data" << endl;
  }

  m_ito->initial_data();
  m_rte->initial_data();

  // Solve Poisson equation
  //  this->solve_poisson();
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

  for (ito_iterator<ito_solver> solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<ito_solver>& solver = solver_it();
    solver->write_checkpoint_level(a_handle, a_lvl);
  }

  for (rte_iterator<mc_photo> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    const RefCountedPtr<mc_photo>& solver = solver_it();
    solver->write_checkpoint_level(a_handle, a_lvl);
  }

  m_poisson->write_checkpoint_level(a_handle, a_lvl);
  m_sigma->write_checkpoint_level(a_handle, a_lvl);
}

void ito_plasma_stepper::read_checkpoint_data(HDF5Handle& a_handle, const int a_lvl){
  CH_TIME("ito_plasma_stepper::read_checkpoint_data");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::read_checkpoint_data" << endl;
  }

  for (ito_iterator<ito_solver> solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    solver->read_checkpoint_level(a_handle, a_lvl);
  }

  for (rte_iterator<mc_photo> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<mc_photo>& solver = solver_it();
    solver->read_checkpoint_level(a_handle, a_lvl);
  }

  m_poisson->read_checkpoint_level(a_handle, a_lvl);
  m_sigma->read_checkpoint_level(a_handle, a_lvl);
}

void ito_plasma_stepper::write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp) const {
  CH_TIME("ito_plasma_stepper::write_plot_data");
  if(m_verbosity > 5){
    pout() << "ito_plasma_stepper::write_plot_data" << endl;
  }

  // Poisson solver copies over its output data
  a_plotvar_names.append(m_poisson->get_plotvar_names());
  m_poisson->write_plot_data(a_output, a_icomp);

  // Surface charge solver writes
  a_plotvar_names.append(m_sigma->get_plotvar_names());
  m_sigma->write_plot_data(a_output, a_icomp);

  // Ito solvers copy their output data
  for (ito_iterator<ito_solver> solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    a_plotvar_names.append(solver->get_plotvar_names());
    solver->write_plot_data(a_output, a_icomp);
  }

  // RTE solvers copy their output data
  for (rte_iterator<mc_photo> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<mc_photo>& solver = solver_it();
    a_plotvar_names.append(solver->get_plotvar_names());
    solver->write_plot_data(a_output, a_icomp);
  }

  // Write the current to the output
#if 0 // Activate when implemented
  this->write_J(a_output, a_icomp);
#endif
  a_plotvar_names.push_back("x-J");
  a_plotvar_names.push_back("y-J");
  if(SpaceDim == 3){
    a_plotvar_names.push_back("z-J");
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

  int ncomp = 0;
  
  for (ito_iterator<ito_solver> solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<ito_solver>& solver = solver_it();
    ncomp += solver->get_num_plotvars();
  }
  
  for (rte_iterator<mc_photo> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<mc_photo>& solver = solver_it();
    ncomp += solver->get_num_plotvars();
  }

  ncomp += m_poisson->get_num_plotvars();
  ncomp += m_sigma->get_num_plotvars();
  ncomp += SpaceDim; // For plotting the current density

  return ncomp;
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
