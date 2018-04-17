/*!
  @file   splitstep_tga.H
  @brief  Declaration of a split-step method for advancing plasma equations. 
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _SPLITSTEP_TGA_
#define _SPLITSTEP_TGA_

#include "type_definitions.H"
#include "amr_mesh.H"
#include "time_stepper.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"

/*!
  @brief Class for evolving plasma equations with a split-step RK2/TGA method. 
  @details This class uses RK2 advancements for the advection and source terms, and an implicit
  TGA method for the diffusion part. 

  Class options
  -------------

  splitstep_tga.rk2_alpha          = 1.0 # Set alpha. 1.0 = Heuns method, 0.5 = Midpoint method
  splitstep_tga.turn_off_advection = false
  splitstep_tga.turn_off_diffusion = false
  splitstep_tga.turn_off_source    = false
  splitstep_tga.turn_off_poisson   = false
  splitstep_tga.turn_off_rte       = false

*/
class splitstep_tga : public time_stepper {
public:

  // Forward declare nested classes. I don't want to clutter this file. 
  class cdr_storage;
  class poisson_storage;
  class rte_storage;
  class sigma_storage;

  /*!
    @brief Constructor
  */
  splitstep_tga();

  /*!
    @brief Destructor
  */
  virtual ~splitstep_tga();

  /*!
    @brief Advance method
  */
  virtual Real advance(const Real a_dt);

  /*!
    @brief Regrid internals
  */
  virtual void regrid_internals();

  /*!
    @brief Compute dt
  */
  virtual void compute_dt(Real& a_dt, time_code::which_code& a_timecode);

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
    @brief Semi-implicit or fully explicit coupling
  */
  bool m_simpi;

  /*!
    @brief Do advection advancement or not
  */
  bool m_do_advection;

  /*!
    @brief Do diffusion advancement or not
  */
  bool m_do_diffusion;

  /*!
    @brief Do source term advancement or not
  */
  bool m_do_source;

  /*!
    @brief Solve RTE equations
  */
  bool m_do_rte;

  /*!
    @brief Solve RTE equations
  */
  bool m_do_poisson;

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
    @brief Restrict dt
  */
  virtual Real restrict_dt();

  /*!
    @brief Backup solutions
  */
  virtual void cache_solutions();

  /*!
    @brief Advection advance
  */
  virtual void advance_advection(const Real a_dt);

  /*!
    @brief Advect diffusion
  */
  virtual void advance_diffusion(const Real a_dt);

  /*!
    @brief Advance source
  */
  virtual void advance_sources(const Real a_dt);

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
  virtual void advance_advection_cdr_k1(const Real a_dt);

  /*!
    @brief Compute k1 and intermediate state for sigma equations
  */
  virtual void advance_advection_sigma_k1(const Real a_dt);

  /*!
    @brief Solve Poisson equation by using the intermediate cdr and sigma states
  */
  virtual void solve_poisson_k1();

  /*!
    @brief Compute the E-field after k1
  */
  virtual void compute_E_after_k1();

  /*!
    @brief Compute the states at the EB after the k1 advance
  */
  virtual void compute_cdr_eb_states_after_k1();

  /*!
    @brief Compute the cdr velocities after doing the k1 advance
  */  
  virtual void compute_cdr_velo_after_k1(const Real a_dt);

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
  virtual void advance_advection_cdr_k2(const Real a_dt);

  /*!
    @brief Advance cdr equations to final state
  */  
  virtual void advance_advection_sigma_k2(const Real a_dt);

  /*!
    @brief Solve poisson equation by using final states
  */  
  virtual void solve_poisson_k2();

  /*!
    @brief Compute cdr sources using whatever is available in the solvers
  */
  virtual void compute_cdr_sources_at_start_of_time_step();

  /*!
    @brief Compute k1 and intermediate state for cdr equations
  */
  virtual void advance_source_cdr_k1(const Real a_dt);

  /*!
    @brief Stationary RTE advance
  */
  virtual void advance_rte_k1_stationary(const Real a_dt);

  /*!
    @brief Transient RTE advance
  */
  virtual void advance_rte_k1_transient(const Real a_dt);

  /*!
    @brief Compute cdr sources using whatever is available in the solvers
  */
  virtual void compute_cdr_sources_after_k1(const Real a_dt);

  /*!
    @brief Advance cdr equations to final state
  */
  virtual void advance_source_cdr_k2(const Real a_dt);

  /*!
    @brief Stationary RTE advance
  */
  virtual void advance_rte_k2_stationary(const Real a_dt);

  /*!
    @brief Transient RTE advance
  */
  virtual void advance_rte_k2_transient(const Real a_dt);
};

#endif
