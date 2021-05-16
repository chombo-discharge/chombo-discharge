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

#include "CD_NamespaceHeader.H"

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
const std::string rte_layout<T>::getRealm() const{
  return m_Realm;
}

template <class T>
void rte_layout<T>::setRealm(const std::string a_realm){
  m_Realm = a_realm;

  for (auto solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    solver_it()->setRealm(m_Realm);
  }
}

template <class T>
void rte_layout<T>::add_solver(RefCountedPtr<T> a_solver){
  m_solvers.push_back(a_solver);
}

template <class T>
void rte_layout<T>::parseOptions(){
  CH_TIME("rte_layout<T>::parseOptions");
  if(m_verbosity > 6){
    pout() << "rte_layout<T>::parseOptions" << endl;
  }

  for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<T>& solver = solver_it();
    solver->parseOptions();
  }
}

template <class T>
void rte_layout<T>::parseRuntimeOptions(){
  CH_TIME("rte_layout<T>::parseRuntimeOptions");
  if(m_verbosity > 6){
    pout() << "rte_layout<T>::parseRuntimeOptions" << endl;
  }

  for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<T>& solver = solver_it();
    solver->parseRuntimeOptions();
  }
}

template <class T>
void rte_layout<T>::allocateInternals(){
  CH_TIME("rte_layout<T>::allocateInternals");
  if(m_verbosity > 6){
    pout() << "rte_layout<T>::allocateInternals" << endl;
  }

  for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<T>& solver = solver_it();
    solver->allocateInternals();
  }
}

template <class T>
void rte_layout<T>::preRegrid(const int a_base, const int a_finestLevel){
  CH_TIME("rte_layout<T>::preRegrid");
  if(m_verbosity > 6){
    pout() << "rte_layout<T>::preRegrid" << endl;
  }

  for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<T>& solver = solver_it();
    solver->preRegrid(a_base, a_finestLevel);
  }
}

template <class T>
void rte_layout<T>::deallocateInternals(){
  CH_TIME("rte_layout<T>::deallocateInternals");
  if(m_verbosity > 6){
    pout() << "rte_layout<T>::deallocateInternals" << endl;
  }

  for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<T>& solver = solver_it();
    solver->deallocateInternals();
  }
}

template <class T>
void rte_layout<T>::setAmr(const RefCountedPtr<AmrMesh>& a_amr){
  CH_TIME("rte_layout<T>::setAmr");
  if(m_verbosity > 5){
    pout() << "rte_layout<T>::setAmr" << endl;
  }

  for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<T>& solver = solver_it();
    solver->setAmr(a_amr);
  }
}

template <class T>
void rte_layout<T>::setComputationalGeometry(const RefCountedPtr<computational_geometry>& a_computationalGeometry){
  CH_TIME("rte_layout<T>::setComputationalGeometry");
  if(m_verbosity > 5){
    pout() << "rte_layout<T>::setComputationalGeometry" << endl;
  }

  for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<T>& solver = solver_it();
    solver->setComputationalGeometry(a_computationalGeometry);
  }
}

template <class T>
void rte_layout<T>::setPhase(phase::which_phase a_phase){
  CH_TIME("rte_layout<T>::setPhase");
  if(m_verbosity > 5){
    pout() << "rte_layout<T>::setPhase" << endl;
  }

  for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<T>& solver = solver_it();
    solver->setPhase(a_phase);
  }
}

template <class T>
void rte_layout<T>::setVerbosity(const int a_verbosity){
  CH_TIME("rte_layout<T>::setVerbosity");

  m_verbosity = a_verbosity;
  if(m_verbosity > 5){
    pout() << "rte_layout<T>::setVerbosity" << endl;
  }

  for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<T>& solver = solver_it();
    solver->setVerbosity(a_verbosity);
  }
}

template <class T>
void rte_layout<T>::sanityCheck(){
  CH_TIME("rte_layout<T>::sanityCheck");
  if(m_verbosity > 5){
    pout() << "rte_layout<T>::sanityCheck" << endl;
  }

  for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<T>& solver = solver_it();
    solver->sanityCheck();
  }
}

template <class T>
void rte_layout<T>::setTime(const int a_step, const Real a_time, const Real a_dt) {
  CH_TIME("rte_layout<T>::setTime");
  if(m_verbosity > 5){
    pout() << "rte_layout<T>::setTime" << endl;
  }

  for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<T>& solver = solver_it();
    solver->setTime(a_step, a_time, a_dt);
  }
}

template <class T>
void rte_layout<T>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel){
  CH_TIME("rte_layout<T>::regrid");
  if(m_verbosity > 5){
    pout() << "rte_layout<T>::regrid" << endl;
  }

  for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<T>& solver = solver_it();
    solver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  }
}

template <class T>
void rte_layout<T>::registerOperators(){
  CH_TIME("rte_layout<T>::registerOperators");
  if(m_verbosity > 5){
    pout() << "rte_layout<T>::registerOperators" << endl;
  }

  for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<T>& solver = solver_it();
    solver->registerOperators();
  }
}

template <class T>
void rte_layout<T>::setSource(const EBAMRCellData& a_source){
  CH_TIME("rte_layout<T>::setSource(ebamrcell)");
  if(m_verbosity > 5){
    pout() << "rte_layout<T>::setSource(ebamrcell)" << endl;
  }

  for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<T>& solver = solver_it();
    solver->setSource(a_source);
  }
}

template <class T>
void rte_layout<T>::setSource(const Real a_source){
  CH_TIME("rte_layout<T>::setSource(constant)");
  if(m_verbosity > 5){
    pout() << "rte_layout<T>::setSource(constant)" << endl;
  }

  for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<T>& solver = solver_it();
    solver->setSource(a_source);
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
void rte_layout<T>::writePlotFile(){
  CH_TIME("rte_layout<T>::writePlotFile");
  if(m_verbosity > 5){
    pout() << "rte_layout<T>::writePlotFile" << endl;
  }

  for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<T>& solver = solver_it();
    solver->writePlotFile();
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
    solver->advance(a_dt, solver->getPhi(), solver->getSource());
  }
}

template <class T>
void rte_layout<T>::initialData(){
  CH_TIME("rte_layout<T>::initialData");
  if(m_verbosity > 6){
    pout() << "rte_layout<T>::initialData" << endl;
  }

  for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<T>& solver = solver_it();
    solver->initialData();
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
Vector<EBAMRCellData*> rte_layout<T>::getSources(){
  CH_TIME("rte_layout<T>::getSources");
  if(m_verbosity > 5){
    pout() << "rte_layout<T>::getSources" << endl;
  }

  Vector<EBAMRCellData*> sources;

  for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<T>& solver = solver_it();
    sources.push_back(&(solver->getSource()));
  }

  return sources;
}

template <class T>
Vector<EBAMRCellData*> rte_layout<T>::getPhis(){
  CH_TIME("rte_layout<T>::getPhis");
  if(m_verbosity > 5){
    pout() << "rte_layout<T>::getPhis" << endl;
  }

  Vector<EBAMRCellData*> states;

  for (rte_iterator<T> solver_it = this->iterator(); solver_it.ok(); ++solver_it){
    RefCountedPtr<T>& solver = solver_it();
    states.push_back(&(solver->getPhi()));
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
    solver->setPhase(phase::gas);
    solver->setVerbosity(-1);
    rte->add_solver(solver);
  }

  return rte;
}
#include "CD_NamespaceFooter.H"

#endif
