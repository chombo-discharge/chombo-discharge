/*!
  @file   splitstep_euler_f.H
  @brief  Declaration of a split step forward Euler method
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _splitstep_forward_euler_
#define _splitstep_forward_euler_

#include "type_definitions.H"
#include "amr_mesh.H"
#include "time_stepper.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"

/*!
  @brief Class for evolving plasma equations with a second order Runge-Kutta method
*/
class splitstep_euler_f : public time_stepper {
public:
  
  // Forward declare nested classes. I don't want to clutter this file.
  class cdr_storage;
  class poisson_storage;
  class rte_storage;
  class sigma_storage;

  /*!
    @brief Constructor
  */
  splitstep_euler_f();

  /*!
    @brief Destructor
  */
  virtual ~splitstep_euler_f();

  /*!
    @brief Advance method
  */
  virtual Real advance(const Real a_dt);

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
    @brief Get storage
  */
  virtual RefCountedPtr<cdr_storage>& get_cdr_storage(const cdr_iterator& a_solverit);

  /*!
    @brief Get storage
  */
  virtual RefCountedPtr<rte_storage>& get_rte_storage(const rte_iterator& a_solverit);

  /*!
    @brief Allocate storage for cdr equations
  */
  virtual void allocate_cdr_storage();

  /*!
    @brief Allocate storage for Poisson equation
  */
  virtual void allocate_poisson_storage();

  /*!
    @brief Allocate storage for RTE equation
  */
  virtual void allocate_rte_storage();

  /*!
    @brief Allocate storage for sigma
  */
  virtual void allocate_sigma_storage();

  /*!
    @brief Compute E at the start of the time step
  */
  virtual void compute_E_at_start_of_time_step();

  /*!
    @brief Compute the cdr velocities at the start of the times step
  */
  virtual void compute_cdr_velo_at_start_of_time_step();

  /*!
    @brief Compute the cdr velocities at the start of the times step
  */
  virtual void compute_cdr_diffco_at_start_of_time_step();

  /*!
    @brief Compute the cdr fluxes at the start of the time step
  */
  virtual void compute_cdr_fluxes_at_start_of_time_step();

  /*!
    @brief Compute the sigma flux at the start of the time step
  */
  virtual void compute_sigma_flux_at_start_of_time_step();

  /*!
    @brief Set source term to zero 
  */
  virtual void set_cdr_source_to_zero_at_start_of_time_step();

  /*!
    @brief Advance transport part of CDR equations
  */
  virtual void advance_cdr_transport(const Real a_dt);

  /*!
    @brief Advance sigma part of transport
  */
  virtual void advance_sigma_transport(const Real a_dt);

  /*!
    @brief Solve the Poisson equation in the transport step
  */
  virtual void solve_poisson_transport();

  /*!
    @brief Solve the RTE equations for the transport step
  */
  virtual void advance_rte_transport_stationary();

  /*!
    @brief Solve the transient RTE equations
  */
  virtual void advance_rte_transport_transient(const Real a_dt);

  /*!
    @brief Compute Electric field after transport step
  */
  virtual void compute_E_after_transport();

  /*!
    @brief Compute the cdr sources after the transport step
  */
  virtual void compute_cdr_source_after_transport();

  /*!
    @brief Source term advancement for cdr equations
  */
  virtual void advance_cdr_source(const Real a_dt);

  /*!
    @brief Source term advancement for sigma equations (this just copies back into the solver)
  */
  virtual void advance_sigma_source(const Real a_dt);

  /*!
    @brief Recompute Poisson after source term step
  */
  virtual void solve_poisson_source();

  /*!
    @brief Compute Electric field after source step
  */
  virtual void compute_E_after_source();
  
  /*!
    @brief Recompute Poisson after source term step
  */
  virtual void advance_rte_source_stationary();

  /*!
    @brief Recompute Poisson after source term step
  */
  virtual void advance_rte_source_transient(const Real a_dt);
  
  /*!
    @brief Restrict dt
  */
  virtual Real restrict_dt();
};

#endif
