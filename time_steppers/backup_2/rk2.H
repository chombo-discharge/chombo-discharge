/*!
  @file   rk2.H
  @brief  Declaration of a second order Runge-Kutta methods for advancing plasma equations
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _rk2_
#define _rk2_

#include "type_definitions.H"
#include "amr_mesh.H"
#include "time_stepper.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"

/*!
  @brief Class for evolving plasma equations with a second order Runge-Kutta method
  @details This class supports various types of Runge-Kutta through the RK2 tableu. 

  Class options
  -------------

      rk2.alpha = 0.5 # Set alpha. 0.5 = Heuns method, 1.0 = Midpoint method

*/
class rk2 : public time_stepper {
public:
  
  // Forward declare nested classes. I don't want to clutter this file.
  class cdr_storage;
  class poisson_storage;
  class rte_storage;
  class sigma_storage;


  
public:

  /*!
    @brief Constructor
  */
  rk2();

  /*!
    @brief Destructor
  */
  virtual ~rk2();

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
    @brief
  */
  virtual void allocate_cdr_storage();

  /*!
    @brief
  */
  virtual void allocate_poisson_storage();

  /*!
    @brief
  */
  virtual void allocate_rte_storage();

  /*!
    @brief
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
  virtual void advance_cdr_k1(const Real a_dt);

  /*!
    @brief Compute k1 and intermediate state for sigma equations
  */
  virtual void advance_sigma_k1(const Real a_dt);

  /*!
    @brief Solve Poisson equation by using the intermediate cdr and sigma states
  */
  virtual void solve_poisson_k1();

  /*!
    @brief Compute the electric field by using the intermediate poisson state
  */
  virtual void compute_E_after_k1();

  /*!
    @brief Compute stationary solution of RTE equations by using intermediate states
  */
  virtual void advance_rte_k1_stationary(const Real a_dt);

  /*!
    @brief Compute transient solution of RTE equations onto intermediate state
  */
  virtual void advance_rte_k1_transient(const Real a_dt);

  /*!
    @brief Compute the states at the EB after the k1 advance
  */
  virtual void compute_cdr_eb_states_after_k1();

  /*!
    @brief Compute the cdr velocities after doing the k1 advance
  */  
  virtual void compute_cdr_velo_after_k1(const Real a_dt);

  /*!
    @brief Compute the cdr diffusion coefficients after doing the k1 advance
  */  
  virtual void compute_cdr_diffco_after_k1(const Real a_dt);
  
  /*!
    @brief Compute the cdr sources after doing the k1 advance
  */  
  virtual void compute_cdr_sources_after_k1(const Real a_dt);

  /*!
    @brief Compute the cdr fluxes after doing the k1 advance
  */  
  virtual void compute_cdr_fluxes_after_k1(const Real a_dt);

  /*!
    @brief Compute the cdr fluxes after doing the k1 advance
  */  
  virtual void compute_sigma_flux_after_k1();

  /*!
    @brief Advance cdr equations to final state
  */  
  virtual void advance_cdr_k2(const Real a_dt);

  /*!
    @brief Advance cdr equations to final state
  */  
  virtual void advance_sigma_k2(const Real a_dt);

  /*!
    @brief Solve poisson equation by using final states
  */  
  virtual void solve_poisson_k2();

  /*!
    @brief Compute E by using final states
  */  
  virtual void compute_E_after_k2();

  /*!
    @brief Compute stationary RTE solutions by using final states
  */    
  virtual void advance_rte_k2_stationary(const Real a_dt);

  /*!
    @brief Compute transient RTE solutions onto final state
  */    
  virtual void advance_rte_k2_transient(const Real a_dt);

  /*!
    @brief Restrict dt
  */
  virtual Real restrict_dt();
  
};

#endif

