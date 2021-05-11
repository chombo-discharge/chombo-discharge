/*!
  @file  FieldSolverFactoryI.H
  @brief Implementation of FieldSolverFactory.H
*/

#ifndef _POISSON_FACTORYI_H_
#define _POISSON_FACTORYI_H_

#include "CD_FieldSolverFactory.H"

#include "CD_NamespaceHeader.H"

template <class T>
FieldSolverFactory<T>::FieldSolverFactory(){
}


template <class T>
FieldSolverFactory<T>::~FieldSolverFactory(){
}

template <class T>
RefCountedPtr<FieldSolver> FieldSolverFactory<T>::new_solver() const {
  return RefCountedPtr<FieldSolver> (new T);
}
#include "CD_NamespaceFooter.H"

#endif
