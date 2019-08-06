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
  class cdr_storage;
  class poisson_storage;
  class rte_storage;
  class sigma_storage;

  euler_maruyama();
  ~euler_maruyama();

  // Pure functions that must be overwritten
  Real advance(const Real a_dt);
  void init_source_terms();
  void regrid_internals();
  void deallocate_internals();
  bool need_to_regrid();

protected:

  // Scratch storage
  Vector<RefCountedPtr<cdr_storage> > m_cdr_scratch;
  Vector<RefCountedPtr<rte_storage> > m_rte_scratch;
  RefCountedPtr<poisson_storage>      m_poisson_scratch;
  RefCountedPtr<sigma_storage>        m_sigma_scratch;

  bool m_extrap_advect;

  // Get functions for storage
  RefCountedPtr<cdr_storage>& get_cdr_storage(const cdr_iterator& a_solverit);
  RefCountedPtr<rte_storage>& get_rte_storage(const rte_iterator& a_solverit);

  Real restrict_dt();

  // Aux functions
  void compute_E_into_scratch();
  void compute_cdr_gradients();
  void compute_cdr_eb_states();
  void compute_cdr_eb_fluxes();


};

#endif
