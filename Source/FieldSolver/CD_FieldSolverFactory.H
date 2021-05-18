/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_FieldSolverFactory.H
  @brief  Declaration of a template factory class for making FieldSolvers
  @author Robert Marskar
*/

#ifndef CD_FieldSolverFactory_H
#define CD_FieldSolverFactory_H

// Chombo includes
#include <RefCountedPtr.H>

// Our includes
#include <CD_FieldSolver.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Factory class for poisson solvers. Used to instantiate any FieldSolver-derived class. 
  @details This class is mostly used for easily making FieldSolvers of any type into a FieldSolver (i.e. casting to smart pointer), but also provides a factory
  method for generating any number of field solvers. 
*/
template <class T>
class FieldSolverFactory {
public:

  /*!
    @brief Constructor. Does nothing. 
  */
  FieldSolverFactory();

  /*!
    @brief Constructor. Does nothing. 
  */
  ~FieldSolverFactory();

  /*!
    @brief Factory method. 
    @details This instantiates a bare pointer of type T and casts it to RefCountedPtr<FieldSolver>. 
  */
  RefCountedPtr<FieldSolver> newSolver() const;
};

#include <CD_NamespaceFooter.H>

#include <CD_FieldSolverFactoryImplem.H>

#endif
