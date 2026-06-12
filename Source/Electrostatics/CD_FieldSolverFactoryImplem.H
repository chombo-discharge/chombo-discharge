/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_FieldSolverFactoryImplem.H
  @brief  Implementation of CD_FieldSolverFactory.H
  @author Robert Marskar
*/

#ifndef CD_FIELDSOLVERFACTORYIMPLEM_H
#define CD_FIELDSOLVERFACTORYIMPLEM_H

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
