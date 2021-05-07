/*!
  @file  field_solver_factoryI.H
  @brief Implementation of field_solver_factory.H
*/

#ifndef _POISSON_FACTORYI_H_
#define _POISSON_FACTORYI_H_

#include "field_solver_factory.H"

namespace ChomboDischarge {

  template <class T>
  field_solver_factory<T>::field_solver_factory(){
  }


  template <class T>
  field_solver_factory<T>::~field_solver_factory(){
  }

  template <class T>
  RefCountedPtr<field_solver> field_solver_factory<T>::new_solver() const {
    return RefCountedPtr<field_solver> (new T);
  }
}

#endif
