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
#include <CD_RtIterator.H>

/*!
  @brief Class for evolving plasma equations using the Euler-Maruyama scheme
*/
class euler_maruyama : public TimeStepper {
public:
  
  // Forward declare nested classes. I don't want to clutter this file.
  class CdrStorage;
  class FieldStorage;
  class RtStorage;
  class SigmaStorage;

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
  Vector<RefCountedPtr<CdrStorage> > m_cdr_scratch;
  Vector<RefCountedPtr<RtStorage> > m_rte_scratch;
  RefCountedPtr<FieldStorage>      m_fieldSolver_scratch;
  RefCountedPtr<SigmaStorage>        m_sigma_scratch;

  bool m_extrap_advect;
  bool m_debug;
  bool m_floor;
  bool m_implicit_diffusion;

  whichDiffusion m_whichDiffusion;

  // Get functions for storage
  RefCountedPtr<CdrStorage>& get_CdrStorage(const CdrIterator& a_solverit);
  RefCountedPtr<RtStorage>& get_RtStorage(const RtIterator& a_solverit);

  Real restrict_dt();

  // Aux functions
  void compute_E_into_scratch();
  void computeCdrGradients();
  void computeCdrEbStates();
  void computeCdrEbFluxes();
  void computeCdrDomainStates();
  void computeCdrDomainFluxes();
  void computeSigmaFlux();
  void computeReactionNetwork(const Real a_dt);
  void advance_cdr(const Real a_dt);
  void advance_rte(const Real a_dt);
  void advance_sigma(const Real a_dt);

  void computeCdrVelo(const Real a_time);
  void compute_cdr_diffco(const Real a_time);

  // Various option parsing
  void parse_diffusion();
  void parse_advection();
  void parse_floor();
  void parse_debug();
#include "CD_NamespaceFooter.H"

#endif
