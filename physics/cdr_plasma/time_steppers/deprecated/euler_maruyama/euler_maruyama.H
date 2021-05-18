/*!
  @file   euler_maruyama.H
  @brief  Declaration of the Euler-Maruyama method
  @author Robert Marskar
  @date   Aug. 2019
*/

#ifndef _EULER_MARUYAMA_
#define _EULER_MARUYAMA_

#include "CD_AmrMesh.H"
#include "CD_TimeStepper.H"
#include <CD_CdrIterator.H>
#include "rte_iterator.H"

/*!
  @brief Class for evolving plasma equations using the Euler-Maruyama scheme
*/
class euler_maruyama : public TimeStepper {
public:
  
  // Forward declare nested classes. I don't want to clutter this file.
  class cdr_storage;
  class poisson_storage;
  class rte_storage;
  class sigma_storage;

  euler_maruyama();
  ~euler_maruyama();

  // Pure functions that must be overwritten
  Real advance(const Real a_dt) override;
  void parseOptions() override;
  void parseRuntimeOptions() override;
  void init();
  void regridInternals(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;
  void allocateInternals() override;
  void deallocateInternals() override;
  void computeDt(Real& a_dt, TimeCode::which_code& a_timeCode) override;
  bool needToRegrid() override;

protected:

  enum whichDiffusion{
    Explicit  = 0,
    Implicit  = 1,
    Automatic = 2
  };

  // Scratch storage
  Vector<RefCountedPtr<cdr_storage> > m_cdr_scratch;
  Vector<RefCountedPtr<rte_storage> > m_rte_scratch;
  RefCountedPtr<poisson_storage>      m_fieldSolver_scratch;
  RefCountedPtr<sigma_storage>        m_sigma_scratch;

  bool m_extrap_advect;
  bool m_debug;
  bool m_floor;
  bool m_implicit_diffusion;

  whichDiffusion m_whichDiffusion;

  // Get functions for storage
  RefCountedPtr<cdr_storage>& get_cdr_storage(const CdrIterator& a_solverit);
  RefCountedPtr<rte_storage>& get_rte_storage(const rte_iterator& a_solverit);

  Real restrict_dt();

  // Aux functions
  void compute_E_into_scratch();
  void compute_cdr_gradients();
  void compute_cdr_eb_states();
  void compute_cdr_eb_fluxes();
  void compute_cdr_domain_states();
  void compute_cdr_domain_fluxes();
  void compute_sigma_flux();
  void compute_reaction_network(const Real a_dt);
  void advance_cdr(const Real a_dt);
  void advance_rte(const Real a_dt);
  void advance_sigma(const Real a_dt);

  void compute_cdr_velo(const Real a_time);
  void compute_cdr_diffco(const Real a_time);

  // Various option parsing
  void parse_diffusion();
  void parse_advection();
  void parse_floor();
  void parse_debug();
#include "CD_NamespaceFooter.H"

#endif
