/*!
  @file poisson_factory.H
  @brief Factory class for Poisson solvers
*/

#ifndef _POISSON_FACTORY_
#define _POISSON_FACTORY_

#include "poisson_solver.H"

/*!
  @brief Factory class for poisson solvers. Used to instantiate any poisson_solver-derived class
*/
template <class T>
class poisson_factory {
public:
  poisson_factory();
  ~poisson_factory();
  RefCountedPtr<poisson_solver> new_solver() const;
};

#endif
