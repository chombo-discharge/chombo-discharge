/*!
  @file   rk2.H
  @brief  Declaration of a second order Runge-Kutta methods for advancing plasma equations
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _rk2_
#define _rk2_

#include "CD_AmrMesh.H"
#include "CD_TimeStepper.H"
#include <CD_CdrIterator.H>
#include <CD_RtIterator.H>

/*!
  @brief Class for evolving plasma equations with a second order Runge-Kutta method
  @details This class supports various types of Runge-Kutta through the RK2 tableu. 

  Class options
  -------------

  rk2.alpha = 0.5 # Set alpha. 0.5 = Heuns method, 1.0 = Midpoint method

*/
class rk2 : public TimeStepper {
public:
  
  // Forward declare nested classes. I don't want to clutter this file.
  class CdrStorage;
  class FieldStorage;
  class RtStorage;
  class SigmaStorage;


  
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
  virtual void regridInternals();

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocateInternals();
  
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
    @brief Get storage
  */
  virtual RefCountedPtr<CdrStorage>& get_CdrStorage(const CdrIterator& a_solverit);

  /*!
    @brief Get storage
  */
  virtual RefCountedPtr<RtStorage>& get_RtStorage(const RtIterator& a_solverit);

  /*!
    @brief
  */
  virtual void allocateCdrStorage();

  /*!
    @brief
  */
  virtual void allocateFieldStorage();

  /*!
    @brief
  */
  virtual void allocateRtStorage();

  /*!
    @brief
  */
  virtual void allocateSigmaStorage();

  /*!
    @brief Compute E at the start of the time step
  */
  virtual void compute_E_at_start_of_time_step();

  /*!
    @brief Compute the velocities at the start of the time step
  */
  virtual void computeCdrVelo_at_start_of_time_step();

  /*!
    @brief Compute the cdr states at the EB at the start of the time step
  */
  virtual void computeCdrEbStates_at_start_of_time_step();

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
  virtual void computeSigmaFlux_at_start_of_time_step();

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
  virtual void computeCdrEbStates_after_k1();

  /*!
    @brief Compute the cdr velocities after doing the k1 advance
  */  
  virtual void computeCdrVelo_after_k1(const Real a_dt);

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
  virtual void computeSigmaFlux_after_k1();

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
  
#include "CD_NamespaceFooter.H"

#endif

