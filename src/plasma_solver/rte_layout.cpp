/*!
  @file   rte_layout.cpp
  @brief  Implementation of rte_layout.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "rte_layout.H"
#include "rte_iterator.H"
#include "eddington_sp1.H"
#include "units.H"

rte_layout::rte_layout(const RefCountedPtr<plasma_kinetics> a_plaskin){
  m_photons = a_plaskin->get_photons();
  m_solvers.resize(m_photons.size());

  for (int i = 0; i < a_plaskin->get_num_photons(); i++){
    m_solvers[i] = RefCountedPtr<rte_solver> (new eddington_sp1());
    m_solvers[i]->set_photon_group(m_photons[i]);
  }

  this->set_verbosity(10);
  this->set_phase(phase::gas);
}

rte_layout::~rte_layout(){

}

void rte_layout::allocate_internals(){
  CH_TIME("rte_layout::allocate_internals");
  if(m_verbosity > 6){
    pout() << "rte_layout::allocate_internals" << endl;
  }

  for (rte_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    solver->allocate_internals();
  }
}

void rte_layout::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  CH_TIME("rte_layout::set_amr");
  if(m_verbosity > 5){
    pout() << "rte_layout::set_amr" << endl;
  }

  for (rte_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    solver->set_amr(a_amr);
  }
}

void rte_layout::set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom){
  CH_TIME("rte_layout::set_computational_geometry");
  if(m_verbosity > 5){
    pout() << "rte_layout::set_computational_geometry" << endl;
  }

  for (rte_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    solver->set_computational_geometry(a_compgeom);
  }
}

void rte_layout::set_phase(phase::which_phase a_phase){
  CH_TIME("rte_layout::set_phase");
  if(m_verbosity > 5){
    pout() << "rte_layout::set_phase" << endl;
  }

  for (rte_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    solver->set_phase(a_phase);
  }
}

void rte_layout::set_verbosity(const int a_verbosity){
  CH_TIME("rte_layout::set_verbosity");
  if(m_verbosity > 5){
    pout() << "rte_layout::set_verbosity" << endl;
  }

  m_verbosity = a_verbosity;

  for (rte_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    solver->set_verbosity(a_verbosity);
  }
}

void rte_layout::set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom){
  CH_TIME("rte_layout::set_physical_domain");
  if(m_verbosity > 5){
    pout() << "rte_layout::set_physical_domain" << endl;
  }

  for (rte_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    solver->set_physical_domain(a_physdom);
  }
}

void rte_layout::sanity_check(){
  CH_TIME("rte_layout::sanity_check");
  if(m_verbosity > 5){
    pout() << "rte_layout::sanity_check" << endl;
  }

  for (rte_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    solver->sanity_check();
  }
}

void rte_layout::set_time(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("rte_layout::set_time");
  if(m_verbosity > 5){
    pout() << "rte_layout::set_time" << endl;
  }

  for (rte_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    solver->set_time(a_step, a_time, a_dt);
  }
}

void rte_layout::regrid(){
  CH_TIME("rte_layout::regrid");
  if(m_verbosity > 5){
    pout() << "rte_layout::regrid" << endl;
  }

  for (rte_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    solver->regrid();
  }
}

void rte_layout::set_source(const EBAMRCellData& a_source){
  CH_TIME("rte_layout::set_source(ebamrcell)");
  if(m_verbosity > 5){
    pout() << "rte_layout::set_source(ebamrcell)" << endl;
  }

  for (rte_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    solver->set_source(a_source);
  }
}

void rte_layout::set_source(const Real a_source){
  CH_TIME("rte_layout::set_source(constant)");
  if(m_verbosity > 5){
    pout() << "rte_layout::set_source(constant)" << endl;
  }

  for (rte_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    solver->set_source(a_source);
  }
}

void rte_layout::write_plot_file(){
  CH_TIME("rte_layout::write_plot_file");
  if(m_verbosity > 5){
    pout() << "rte_layout::write_plot_file" << endl;
  }

  for (rte_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    solver->write_plot_file();
  }
}

void rte_layout::advance(const Real a_dt){
  CH_TIME("rte_layout::advance");
  if(m_verbosity > 6){
    pout() << "rte_layout::advance" << endl;
  }

  for (rte_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    solver->advance(a_dt);
  }
}

void rte_layout::initial_data(){
  CH_TIME("rte_layout::initial_data");
  if(m_verbosity > 6){
    pout() << "rte_layout::initial_data" << endl;
  }

  for (rte_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    solver->initial_data();
  }
}

bool rte_layout::is_stationary(){
  CH_TIME("rte_layout::is_stationary");
  if(m_verbosity > 5){
    pout() << "rte_layout::is_stationary" << endl;
  }

  bool stationary = true;

  for (rte_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();

    if(solver->is_stationary() == false){
      stationary = false;
    }
  }

  return stationary;
}


phase::which_phase rte_layout::get_phase(){
  CH_TIME("rte_layout::get_phase");
  if(m_verbosity > 5){
    pout() << "rte_layout::get_phase" << endl;
  }
  
 for (rte_iterator solver_it(*this); solver_it.ok(); ++solver_it){
   RefCountedPtr<rte_solver>& solver = solver_it();
   return solver->get_phase();
 }
}

Vector<RefCountedPtr<rte_solver> >& rte_layout::get_solvers(){
  return m_solvers;
}

Vector<RefCountedPtr<photon_group> >& rte_layout::get_photons(){
  return m_photons;
}

Vector<EBAMRCellData*> rte_layout::get_sources(){
  CH_TIME("rte_layout::get_sources");
  if(m_verbosity > 5){
    pout() << "rte_layout::get_sources" << endl;
  }

  Vector<EBAMRCellData*> sources;

  for (rte_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    sources.push_back(&(solver->get_source()));
  }

  return sources;
}


Vector<EBAMRCellData*> rte_layout::get_states(){
  CH_TIME("rte_layout::get_states");
  if(m_verbosity > 5){
    pout() << "rte_layout::get_states" << endl;
  }

  Vector<EBAMRCellData*> states;

  for (rte_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<rte_solver>& solver = solver_it();
    states.push_back(&(solver->get_state()));
  }

  return states;
}
