/*!
  @file field_solver_factory.H
  @brief Factory class for Poisson solvers
*/

#ifndef _POISSON_FACTORY_
#define _POISSON_FACTORY_

#include "field_solver.H"

#include "CD_NamespaceHeader.H"

/*!
  @brief Factory class for poisson solvers. Used to instantiate any field_solver-derived class
*/
template <class T>
class field_solver_factory {
public:
  field_solver_factory();
  ~field_solver_factory();
  RefCountedPtr<field_solver> new_solver() const;
};
#include "CD_NamespaceFooter.H"

#include "field_solver_factoryI.H"

#endif
