/*!
  @file   cdr_layout.H
  @brief  Declaration of a class that holds a set of cdr_solvers
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "cdr_layout.H"

template <class T>
cdr_factory<T>::cdr_factory(){
}

template <class T>
cdr_factory<T>::~cdr_factory(){
}


template <class T>
RefCountedPtr<cdr_layout> cdr_factory<T>::new_layout(const RefCountedPtr<plasma_kinetics> a_plaskin) const{

  // Build cdr layout
  auto cdr = RefCountedPtr<cdr_layout> (new cdr_layout(a_plaskin));
  auto spe = a_plaskin->get_species();

  // Cast solvers and instantiate them
  for (int i = 0; i < a_plaskin->get_num_species(); i++){
    auto solver = RefCountedPtr<cdr_solver> (new T());
    solver->set_species(spe[i]);
    solver->set_phase(phase::gas);
    solver->set_verbosity(-1);
    cdr->add_solver(solver);
  }

  return cdr;
}
