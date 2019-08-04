/*!
  @file   euler_maruyama.H
  @brief  Declaration of the Euler-Maruyama method
  @author Robert Marskar
  @date   Aug. 2019
*/

#ifndef _EULER_MARUYAMA_
#define _EULER_MARUYAMA_

#include "type_definitions.H"
#include "amr_mesh.H"
#include "time_stepper.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"

/*!
  @brief Class for evolving plasma equations using the Euler-Maruyama scheme
*/
class euler_maruyama : public time_stepper {
public:
  
  // Forward declare nested classes. I don't want to clutter this file.
  // class cdr_storage;
  // class poisson_storage;
  // class rte_storage;
  // class sigma_storage;
  
public:

  euler_maruyama();
  ~euler_maruyama();

  // Pure functions that must be overwritten
  Real advance(const Real a_dt);
  void init_source_terms();
  void regrid_internals();
  void deallocate_internals();
  bool need_to_regrid();

protected:

  bool m_extrap_advect;

  Real restrict_dt();
};

#endif
