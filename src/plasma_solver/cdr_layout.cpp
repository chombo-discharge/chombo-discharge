/*!
  @file   cdr_layout.cpp
  @brief  Implementation of cdr_layout.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "cdr_layout.H"
#include "cdr_iterator.H"
#include "cdr_sg.H"
#include "cdr_gdnv.H"

cdr_layout::cdr_layout(const RefCountedPtr<plasma_kinetics> a_plaskin){

  m_species = a_plaskin->get_species();
  m_solvers.resize(m_species.size());
  for (int i = 0; i < a_plaskin->get_num_species(); i++){
    m_solvers[i] = RefCountedPtr<cdr_solver> (new cdr_sg());
    m_solvers[i]->set_species(m_species[i]);
  }

  this->set_verbosity(-1);
}

cdr_layout::~cdr_layout(){

}

void cdr_layout::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
  CH_TIME("cdr_layout::set_amr");
  if(m_verbosity > 5){
    pout() << "cdr_layout::set_amr" << endl;
  }

  for (cdr_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->set_amr(a_amr);
  }
}

void cdr_layout::set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom){
  CH_TIME("cdr_layout::set_computational_geometry");
  if(m_verbosity > 5){
    pout() << "cdr_layout::set_computational_geometry" << endl;
  }

  for (cdr_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->set_computational_geometry(a_compgeom);
  }
}

void cdr_layout::set_phase(phase::which_phase a_phase){
  CH_TIME("cdr_layout::set_verbosity");
  if(m_verbosity > 5){
    pout() << "cdr_layout::set_verbosity" << endl;
  }

  for (cdr_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->set_phase(a_phase);
  }
}

void cdr_layout::set_verbosity(const int a_verbosity){
  CH_TIME("cdr_layout::set_verbosity");
  if(m_verbosity > 5){
    pout() << "cdr_layout::set_verbosity" << endl;
  }

  for (cdr_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->set_verbosity(a_verbosity);
  }
}

void cdr_layout::set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom){
  CH_TIME("cdr_layout::set_physical_domain");
  if(m_verbosity > 5){
    pout() << "cdr_layout::set_physical_domain" << endl;
  }

  for (cdr_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->set_physical_domain(a_physdom);
  }
}

void cdr_layout::set_mass_redist(const bool a_mass_redist){
  CH_TIME("cdr_layout::set_mass_redist");
  if(m_verbosity > 5){
    pout() << "cdr_layout::set_mass_redist" << endl;
  }

  for (cdr_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->set_mass_redist(a_mass_redist);
  }
}

void cdr_layout::set_time(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("cdr_layout::set_time");
  if(m_verbosity > 5){
    pout() << "cdr_layout::set_time" << endl;
  }

  for (cdr_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->set_time(a_step, a_time, a_dt);
  }
}

void cdr_layout::regrid(){
  CH_TIME("cdr_layout::regrid");
  if(m_verbosity > 5){
    pout() << "cdr_layout::regrid" << endl;
  }

  for (cdr_iterator solver_it(*this); solver_it.ok(); ++solver_it){
    RefCountedPtr<cdr_solver>& solver = solver_it();
    solver->regrid();
  }
}

Vector<RefCountedPtr<cdr_solver> >& cdr_layout::get_solvers(){
  return m_solvers;
}

Vector<RefCountedPtr<species> >& cdr_layout::get_species(){
  return m_species;
}
