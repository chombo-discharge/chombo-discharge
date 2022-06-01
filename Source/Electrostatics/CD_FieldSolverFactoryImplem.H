/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_FieldSolverFactoryImplem.H
  @brief  Implementation of CD_FieldSolverFactory.H
  @author Robert Marskar
*/

#ifndef CD_FieldSolverFactoryImplem_H
#define CD_FieldSolverFactoryImplem_H

// Our includes
#include <CD_FieldSolverFactory.H>
#include <CD_NamespaceHeader.H>

template <class T>
FieldSolverFactory<T>::FieldSolverFactory()
{
  CH_TIME("FieldSolverFactory<T>::FieldSolverFactory()");
}

template <class T>
FieldSolverFactory<T>::~FieldSolverFactory()
{
  CH_TIME("FieldSolverFactory<T>::~FieldSolverFactory()");
}

template <class T>
RefCountedPtr<FieldSolver>
FieldSolverFactory<T>::newSolver() const
{
  CH_TIME("FieldSolverFactory<T>::newSolver()");

  return RefCountedPtr<FieldSolver>(new T);
}

#include <CD_NamespaceFooter.H>

#endif
