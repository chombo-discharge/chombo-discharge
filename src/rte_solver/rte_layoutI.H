/*!
  @file   rte_layoutI.H
  @brief  Implementation of rte_layout.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _RTE_LAYOUTI_H
#define _RTE_LAYOUTI_H

#include "rte_iterator.H"
#include "rte_layout.H"

namespace ChomboDischarge {

  template <class T>
  rte_layout<T>::rte_layout(const Vector<RefCountedPtr<rte_species> >& a_rte_species){
    m_species = a_rte_species;
    m_solvers.resize(0); // Solvers added through the add_solver function
  }

  template <class T>
  rte_layout<T>::~rte_layout(){

  }

  template <class T>
  rte_iterator<T> rte_layout<T>::iterator(){
    return rte_iterator<T>(*this);
  }

  template <class T>
  const std::string rte_layout<T>::get_realm() const{
    return m_realm;
  }

  template <class T>
  void rte_layout<T>::set_realm(const std::string a_realm){
    m_realm = a_realm;

    for (auto solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      solver_it()->set_realm(m_realm);
    }
  }

  template <class T>
  void rte_layout<T>::add_solver(RefCountedPtr<T> a_solver){
    m_solvers.push_back(a_solver);
  }

  template <class T>
  void rte_layout<T>::parse_options(){
    CH_TIME("rte_layout<T>::parse_options");
    if(m_verbosity > 6){
      pout() << "rte_layout<T>::parse_options" << endl;
    }

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      solver->parse_options();
    }
  }

  template <class T>
  void rte_layout<T>::parse_runtime_options(){
    CH_TIME("rte_layout<T>::parse_runtime_options");
    if(m_verbosity > 6){
      pout() << "rte_layout<T>::parse_runtime_options" << endl;
    }

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      solver->parse_runtime_options();
    }
  }

  template <class T>
  void rte_layout<T>::allocate_internals(){
    CH_TIME("rte_layout<T>::allocate_internals");
    if(m_verbosity > 6){
      pout() << "rte_layout<T>::allocate_internals" << endl;
    }

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      solver->allocate_internals();
    }
  }

  template <class T>
  void rte_layout<T>::pre_regrid(const int a_base, const int a_finest_level){
    CH_TIME("rte_layout<T>::pre_regrid");
    if(m_verbosity > 6){
      pout() << "rte_layout<T>::pre_regrid" << endl;
    }

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      solver->pre_regrid(a_base, a_finest_level);
    }
  }

  template <class T>
  void rte_layout<T>::deallocate_internals(){
    CH_TIME("rte_layout<T>::deallocate_internals");
    if(m_verbosity > 6){
      pout() << "rte_layout<T>::deallocate_internals" << endl;
    }

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      solver->deallocate_internals();
    }
  }

  template <class T>
  void rte_layout<T>::set_amr(const RefCountedPtr<amr_mesh>& a_amr){
    CH_TIME("rte_layout<T>::set_amr");
    if(m_verbosity > 5){
      pout() << "rte_layout<T>::set_amr" << endl;
    }

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      solver->set_amr(a_amr);
    }
  }

  template <class T>
  void rte_layout<T>::set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom){
    CH_TIME("rte_layout<T>::set_computational_geometry");
    if(m_verbosity > 5){
      pout() << "rte_layout<T>::set_computational_geometry" << endl;
    }

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      solver->set_computational_geometry(a_compgeom);
    }
  }

  template <class T>
  void rte_layout<T>::set_phase(phase::which_phase a_phase){
    CH_TIME("rte_layout<T>::set_phase");
    if(m_verbosity > 5){
      pout() << "rte_layout<T>::set_phase" << endl;
    }

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      solver->set_phase(a_phase);
    }
  }

  template <class T>
  void rte_layout<T>::set_verbosity(const int a_verbosity){
    CH_TIME("rte_layout<T>::set_verbosity");

    m_verbosity = a_verbosity;
    if(m_verbosity > 5){
      pout() << "rte_layout<T>::set_verbosity" << endl;
    }

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      solver->set_verbosity(a_verbosity);
    }
  }

  template <class T>
  void rte_layout<T>::sanity_check(){
    CH_TIME("rte_layout<T>::sanity_check");
    if(m_verbosity > 5){
      pout() << "rte_layout<T>::sanity_check" << endl;
    }

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      solver->sanity_check();
    }
  }

  template <class T>
  void rte_layout<T>::set_time(const int a_step, const Real a_time, const Real a_dt) {
    CH_TIME("rte_layout<T>::set_time");
    if(m_verbosity > 5){
      pout() << "rte_layout<T>::set_time" << endl;
    }

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      solver->set_time(a_step, a_time, a_dt);
    }
  }

  template <class T>
  void rte_layout<T>::regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level){
    CH_TIME("rte_layout<T>::regrid");
    if(m_verbosity > 5){
      pout() << "rte_layout<T>::regrid" << endl;
    }

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      solver->regrid(a_lmin, a_old_finest_level, a_new_finest_level);
    }
  }

  template <class T>
  void rte_layout<T>::register_operators(){
    CH_TIME("rte_layout<T>::register_operators");
    if(m_verbosity > 5){
      pout() << "rte_layout<T>::register_operators" << endl;
    }

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      solver->register_operators();
    }
  }

  template <class T>
  void rte_layout<T>::set_source(const EBAMRCellData& a_source){
    CH_TIME("rte_layout<T>::set_source(ebamrcell)");
    if(m_verbosity > 5){
      pout() << "rte_layout<T>::set_source(ebamrcell)" << endl;
    }

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      solver->set_source(a_source);
    }
  }

  template <class T>
  void rte_layout<T>::set_source(const Real a_source){
    CH_TIME("rte_layout<T>::set_source(constant)");
    if(m_verbosity > 5){
      pout() << "rte_layout<T>::set_source(constant)" << endl;
    }

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      solver->set_source(a_source);
    }
  }

  template <class T>
  void rte_layout<T>::set_stationary(const bool a_stationary){
    CH_TIME("rte_layout<T>::set_stationary");
    if(m_verbosity > 5){ 
      pout() << "rte_layout<T>::set_stationary" << endl;
    }

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      solver->set_stationary(a_stationary);
    }
  }

  template <class T>
  void rte_layout<T>::write_plot_file(){
    CH_TIME("rte_layout<T>::write_plot_file");
    if(m_verbosity > 5){
      pout() << "rte_layout<T>::write_plot_file" << endl;
    }

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      solver->write_plot_file();
    }
  }

  template <class T>
  void rte_layout<T>::advance(const Real a_dt){
    CH_TIME("rte_layout<T>::advance");
    if(m_verbosity > 6){
      pout() << "rte_layout<T>::advance" << endl;
    }

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      solver->advance(a_dt, solver->get_state(), solver->get_source());
    }
  }

  template <class T>
  void rte_layout<T>::initial_data(){
    CH_TIME("rte_layout<T>::initial_data");
    if(m_verbosity > 6){
      pout() << "rte_layout<T>::initial_data" << endl;
    }

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      solver->initial_data();
    }
  }

  template <class T>
  bool rte_layout<T>::is_stationary(){
    CH_TIME("rte_layout<T>::is_stationary");
    if(m_verbosity > 5){
      pout() << "rte_layout<T>::is_stationary" << endl;
    }

    bool stationary = true;

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();

      if(solver->is_stationary() == false){
	stationary = false;
      }
    }

    return stationary;
  }
  template <class T>

  phase::which_phase rte_layout<T>::get_phase(){
    CH_TIME("rte_layout<T>::get_phase");
    if(m_verbosity > 5){
      pout() << "rte_layout<T>::get_phase" << endl;
    }

    phase::which_phase p;
    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      p = solver->get_phase();
    }

    return p;
  }

  template <class T>
  Vector<RefCountedPtr<T> >& rte_layout<T>::get_solvers(){
    return m_solvers;
  }

  template <class T>
  Vector<RefCountedPtr<rte_species> >& rte_layout<T>::get_species(){
    return m_species;
  }

  template <class T>
  Vector<EBAMRCellData*> rte_layout<T>::get_sources(){
    CH_TIME("rte_layout<T>::get_sources");
    if(m_verbosity > 5){
      pout() << "rte_layout<T>::get_sources" << endl;
    }

    Vector<EBAMRCellData*> sources;

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      sources.push_back(&(solver->get_source()));
    }

    return sources;
  }

  template <class T>
  Vector<EBAMRCellData*> rte_layout<T>::get_states(){
    CH_TIME("rte_layout<T>::get_states");
    if(m_verbosity > 5){
      pout() << "rte_layout<T>::get_states" << endl;
    }

    Vector<EBAMRCellData*> states;

    for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
      RefCountedPtr<T>& solver = solver_it();
      states.push_back(&(solver->get_state()));
    }

    return states;
  }

  template <class T, class S>
  rte_factory<T, S>::rte_factory(){
  }

  template <class T, class S>
  rte_factory<T, S>::~rte_factory(){
  }


  template <class T, class S>
  RefCountedPtr<rte_layout<T> > rte_factory<T, S>::new_layout(const Vector<RefCountedPtr<rte_species> >& a_species) const{
    // Build rte layout
    auto rte = RefCountedPtr<rte_layout<T> > (new rte_layout<T> (a_species));
    auto spe = a_species;

    // Cast solvers and instantiate them
    for (int i = 0; i < a_species.size(); i++){
      auto solver = RefCountedPtr<T> (static_cast<T*> (new S()));
      solver->set_rte_species(spe[i]);
      solver->set_phase(phase::gas);
      solver->set_verbosity(-1);
      rte->add_solver(solver);
    }

    return rte;
  }
}

#endif
