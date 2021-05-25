/*!
  @file   rk2_tga.H
  @brief  Declaration of a split-step method for advancing plasma equations. 
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _RK2_TGA_
#define _RK2_TGA_

#include "CD_AmrMesh.H"
#include "CD_TimeStepper.H"
#include <CD_CdrIterator.H>
#include <CD_RtIterator.H>

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
class rk2_tga : public TimeStepper {
public:

  // Forward declare nested classes. I don't want to clutter this file. 
  class CdrStorage;
  class FieldStorage;
  class RtStorage;
  class SigmaStorage;

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
  virtual void deallocateInternals();

  /*!
    @brief Regrid internals
  */
  virtual void regridInternals();

  /*!
    @brief Compute dt
  */
  virtual void computeDt(Real& a_dt, TimeCode::which_code& a_timeCode);

protected:

  /*!
    @brief Scratch storage
  */
  Vector<RefCountedPtr<CdrStorage> > m_cdr_scratch;

  /*!
    @brief Scratch storage
  */
  RefCountedPtr<FieldStorage> m_fieldSolver_scratch;

  /*!
    @brief Scratch storage
  */
  Vector<RefCountedPtr<RtStorage> > m_rte_scratch;

  /*!
    @brief Scratch storage
  */
  RefCountedPtr<SigmaStorage> m_sigma_scratch;

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
  virtual RefCountedPtr<CdrStorage>& get_CdrStorage(const CdrIterator& a_solverit);

  /*!
    @brief Get storage
  */
  virtual RefCountedPtr<RtStorage>& get_RtStorage(const RtIterator& a_solverit);

  /*!
    @brief Allocate storage
  */
  virtual void allocateCdrStorage();

  /*!
    @brief Allocate storage
  */
  virtual void allocateFieldStorage();

  /*!
    @brief Allocate storage
  */
  virtual void allocateRtStorage();

  /*!
    @brief Allocate storage
  */
  virtual void allocateSigmaStorage();

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
  virtual void computeCdrVelo(const Real a_time);

  /*!
    @brief Compute the cdr states at the EB at the start of the time step
  */
  virtual void computeCdrEbStates();

  /*!
    @brief Compute the cdr fluxes the start of the time step
  */
  virtual void compute_cdr_fluxes(const Real a_dt);

  /*!
    @brief Compute the sigma flux the start of the time step
  */
  virtual void computeSigmaFlux_into_scratch();

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
#include "CD_NamespaceFooter.H"

#endif
