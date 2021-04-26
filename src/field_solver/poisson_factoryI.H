/*!
  @file  poisson_factoryI.H
  @brief Implementation of poisson_factory.H
*/

#include "poisson_factory.H"

template <class T>
poisson_factory<T>::poisson_factory(){
}


template <class T>
poisson_factory<T>::~poisson_factory(){
}

template <class T>
RefCountedPtr<poisson_solver> poisson_factory<T>::new_solver() const {
  return RefCountedPtr<poisson_solver> (new T);
}
