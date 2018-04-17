/*!
  @file   euler_f.H
  @brief  Declaration of a backward Euler method for advancing the plasma equations
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _euler_f_
#define _euler_f_

#include "type_definitions.H"
#include "amr_mesh.H"
#include "time_stepper.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"

/*!
  @brief Class for evolving plasma equations with a second order Runge-Kutta method
*/
class euler_f : public time_stepper {
public:

  class cdr_storage; 
  class poisson_storage;
  class rte_storage;
  class sigma_storage;

  /*!
    @brief Constructor
  */
  euler_f();

  /*!
    @brief Destructor
  */
  virtual ~euler_f();

  /*!
    @brief Advance method
  */
  virtual Real advance(const Real a_dt);

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocate_internals();

  /*!
    @brief Regrid internals
  */
  virtual void regrid_internals();
  
protected:

  /*!
    @brief Scratch storage
  */
  Vector<RefCountedPtr<cdr_storage> > m_cdr_scratch;

  /*!
    @brief Scratch storage
  */
  RefCountedPtr<poisson_storage> m_poisson_scratch;

  /*!
    @brief Scratch storage
  */
  Vector<RefCountedPtr<rte_storage> > m_rte_scratch;

  /*!
    @brief Scratch storage
  */
  RefCountedPtr<sigma_storage> m_sigma_scratch;

  /*!
    @brief RK2 method. 
  */
  Real m_alpha;

  /*!
    @brief Get storage
  */
  virtual RefCountedPtr<cdr_storage>& get_cdr_storage(const cdr_iterator& a_solverit);

  /*!
    @brief Get storage
  */
  virtual RefCountedPtr<rte_storage>& get_rte_storage(const rte_iterator& a_solverit);

  /*!
    @brief Allocate storage
  */
  virtual void allocate_cdr_storage();

  /*!
    @brief Allocate storage
  */
  virtual void allocate_poisson_storage();

  /*!
    @brief Allocate storage
  */
  virtual void allocate_rte_storage();

  /*!
    @brief Allocate storage
  */
  virtual void allocate_sigma_storage();

  /*!
    @brief Compute E at the start of the time step
  */
  virtual void compute_E_at_start_of_time_step(); 

  /*!
    @brief Compute the velocities at the start of the time step
  */
  virtual void compute_cdr_velo_at_start_of_time_step();

  /*!
    @brief Compute the cdr states at the EB at the start of the time step
  */
  virtual void compute_cdr_eb_states_at_start_of_time_step();

  /*!
    @brief Compute the diffusion coefficients at the start of the time step
  */
  virtual void compute_cdr_diffco_at_start_of_time_step();

  /*!
    @brief Compute the cdr source terms at the start of the time step
  */
  virtual void compute_cdr_sources_at_start_of_time_step();

  /*!
    @brief Compute the cdr fluxes the start of the time step
  */
  virtual void compute_cdr_fluxes_at_start_of_time_step();

  /*!
    @brief Compute the sigma flux the start of the time step
  */
  virtual void compute_sigma_flux_at_start_of_time_step();

  /*!
    @brief Compute k1 and intermediate state for cdr equations
  */
  virtual void advance_cdr(const Real a_dt);

  /*!
    @brief Compute k1 and intermediate state for sigma equations
  */
  virtual void advance_sigma(const Real a_dt);

  /*!
    @brief Solve Poisson equation by using the intermediate cdr and sigma states
  */
  virtual void advance_poisson();

  /*!
    @brief Compute the electric field by using the new Poisson solution
  */
  virtual void compute_E_after_poisson();

  /*!
    @brief Compute stationary solution of RTE equations by using intermediate states
  */
  virtual void advance_rte_stationary(const Real a_dt);

  /*!
    @brief Compute transient solution of RTE equations onto intermediate state
  */
  virtual void advance_rte_transient(const Real a_dt);

  /*!
    @brief Restrict dt
  */
  virtual Real restrict_dt();
  
};

#endif

