/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_FieldSolverFactory.H
  @brief  Declaration of a template factory class for making FieldSolvers
  @author Robert Marskar
*/

#ifndef CD_FIELDSOLVERFACTORY_H
#define CD_FIELDSOLVERFACTORY_H

// Chombo includes
#include <RefCountedPtr.H>

// Our includes
#include <CD_FieldSolver.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Factory class for poisson solvers. Used to instantiate any FieldSolver-derived classes with automatic casting to the base class. 
  @details This class is mostly used for easily making FieldSolvers of any type into a FieldSolver (i.e. casting to smart pointer), but as a side-effect
  it also provides a factory method for generating any number of field solvers. 
*/
template <class T>
class FieldSolverFactory
{
public:
  /**
    @brief Constructor (does nothing)
  */
  FieldSolverFactory();

  /**
    @brief Constructor (does nothing)
  */
  ~FieldSolverFactory();

  /**
    @brief Factory method. 
    @details This instantiates a bare pointer of type T and casts it to RefCountedPtr<FieldSolver>. 
  */
  RefCountedPtr<FieldSolver>
  newSolver() const;
};

#include <CD_NamespaceFooter.H>

#include <CD_FieldSolverFactoryImplem.H>

#endif
