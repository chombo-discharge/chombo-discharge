/*!
  @file   rte_layoutI.H
  @brief  Implementation of rte_layout.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "rte_layout.H"

template <class T>
rte_factory<T>::rte_factory(){
}

template <class T>
rte_factory<T>::~rte_factory(){
}


template <class T>
RefCountedPtr<rte_layout> rte_factory<T>::new_layout(const RefCountedPtr<plasma_kinetics> a_plaskin) const{

  // Build rte layout
  auto rte = RefCountedPtr<rte_layout> (new rte_layout(a_plaskin));
  auto spe = a_plaskin->get_photons();

  // Cast solvers and instantiate them
  for (int i = 0; i < a_plaskin->get_num_photons(); i++){
    auto solver = RefCountedPtr<rte_solver> (new T());
    solver->set_photon_group(spe[i]);
    solver->set_phase(phase::gas);
    solver->set_verbosity(-1);
    rte->add_solver(solver);
  }

  return rte;
}
