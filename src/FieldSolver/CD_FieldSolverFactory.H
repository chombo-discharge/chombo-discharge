/*!
  @file FieldSolverFactory.H
  @brief Factory class for Poisson solvers
*/

#ifndef _POISSON_FACTORY_
#define _POISSON_FACTORY_

#include "CD_FieldSolver.H"

#include "CD_NamespaceHeader.H"

/*!
  @brief Factory class for poisson solvers. Used to instantiate any FieldSolver-derived class
*/
template <class T>
class FieldSolverFactory {
public:
  FieldSolverFactory();
  ~FieldSolverFactory();
  RefCountedPtr<FieldSolver> new_solver() const;
};
#include "CD_NamespaceFooter.H"

#include "CD_FieldSolverFactoryImplem.H"

#endif
