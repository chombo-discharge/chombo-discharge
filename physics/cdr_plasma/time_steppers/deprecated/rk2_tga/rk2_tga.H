/*!
  @file   rk2_tga.H
  @brief  Declaration of a split-step method for advancing plasma equations. 
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _RK2_TGA_
#define _RK2_TGA_

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

  rk2_tga.rk2_alpha          = 1.0    # Set alpha. 1.0 = Heuns method, 0.5 = Midpoint method
  rk2_tga.turn_off_advection = false  # Turn off advection, debugging only
  rk2_tga.turn_off_diffusion = false  # Turn off diffusion, debugging only
  rk2_tga.turn_off_source    = false  # Turn off source term, debugging only
  rk2_tga.turn_off_poisson   = false  # Turn off Poisson, debugging only
  rk2_tga.turn_off_rte       = false  # Turn off RTE, debugging only

*/
class rk2_tga : public time_stepper {
public:

  // Forward declare nested classes. I don't want to clutter this file. 
  class cdr_storage;
  class poisson_storage;
  class rte_storage;
  class sigma_storage;

  /*!
    @brief Constructor
  */
  rk2_tga();

  /*!
    @brief Destructor
  */
  virtual ~rk2_tga();

  /*!
    @brief Advance method
  */
  virtual Real advance(const Real a_dt);

  /*!
    @brief Deallocate storage
  */
  virtual void deallocate_internals();

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
    @brief Do advection advancement or not
  */
  bool m_do_advec_src;

  /*!
    @brief Do diffusion advancement or not
  */
  bool m_do_diffusion;

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
    @brief Advection and source advance
  */
  virtual void advance_advection_source(const Real a_dt);

  /*!
    @brief Advect diffusion
  */
  virtual void advance_diffusion(const Real a_dt);

  /*!
    @brief Compute E at the start of the time step
  */
  virtual void compute_E_into_scratch();

  /*!
    @brief Compute the velocities
  */
  virtual void compute_cdr_velo(const Real a_time);

  /*!
    @brief Compute the cdr states at the EB at the start of the time step
  */
  virtual void compute_cdr_eb_states();

  /*!
    @brief Compute the cdr fluxes the start of the time step
  */
  virtual void compute_cdr_fluxes(const Real a_dt);

  /*!
    @brief Compute the sigma flux the start of the time step
  */
  virtual void compute_sigma_flux_into_scratch();

    /*!
    @brief Compute k1 and intermediate state for cdr equations
  */
  virtual void advance_advection_source_cdr_k1(const Real a_dt);

  /*!
    @brief Compute k1 and intermediate state for sigma equations
  */
  virtual void advance_advection_sigma_k1(const Real a_dt);

  /*!
    @brief Advance cdr equations to final state
  */  
  virtual void advance_advection_source_cdr_k2(const Real a_dt);

  /*!
    @brief Advance cdr equations to final state
  */  
  virtual void advance_advection_sigma_k2(const Real a_dt);

  /*!
    @brief Compute cdr sources using whatever is available in the solvers
  */
  virtual void compute_cdr_sources_into_scratch(const Real a_time);

  /*!
    @brief Stationary RTE advance
  */
  virtual void advance_rte_stationary(const Real a_time);
};

#endif
