/*!
  @file   cdr_layoutI.H
  @brief  Implementation of cdr_layout.H
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
RefCountedPtr<cdr_layout> cdr_factory<T>::new_layout(const Vector<RefCountedPtr<cdr_species> > a_species) const{

  // Build cdr layout
  auto cdr = RefCountedPtr<cdr_layout> (new cdr_layout(a_species));

  // Cast solvers and instantiate them
  for (int i = 0; i < a_species.size(); i++){
    auto solver = RefCountedPtr<cdr_solver> (new T());
    solver->set_species(a_species[i]);
    solver->set_phase(phase::gas);
    solver->set_verbosity(-1);
    cdr->add_solver(solver);
  }

  return cdr;
}
