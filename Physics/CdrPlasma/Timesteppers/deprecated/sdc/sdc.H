/*!
  @file   sdc.H
  @brief  Declaration of a semi-implicit spectral deferred correction class
  @author Robert Marskar
  @date   Feb. 2019
*/

#ifndef _SDC_
#define _SDC_

#include "CD_AmrMesh.H"
#include "CD_TimeStepper.H"
#include <CD_CdrIterator.H>
#include <CD_RtIterator.H>

extern "C" void dgesv_(int*    N,
		       int*    NRHS,
		       double* A,
		       int*    LDA,
		       int*    IPIV,
		       double* B,
		       int*    LDB,
		       int*    INFO);

/*!
  @brief Class for evolving plasma equations with a multirate method using semi-implicit spectral deferred corrections
*/
class sdc : public TimeStepper {
public:
  
  // Forward declare nested classes. I don't want to clutter this file.
  class CdrStorage;
  class FieldStorage;
  class RtStorage;
  class SigmaStorage;
  
public:

  sdc();
  ~sdc();
  
  Real advance(const Real a_dt);

  void init();
  void cache_internals();
  void regridInternals(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);
  void allocateInternals();
  void deallocateInternals();
  void parseOptions();
  bool needToRegrid();

  Vector<EBAMRCellData*> get_cdr_errors();
  
protected:

  Vector<RefCountedPtr<CdrStorage> > m_cdr_scratch;
  Vector<RefCountedPtr<RtStorage> > m_rte_scratch;
  RefCountedPtr<FieldStorage>      m_fieldSolver_scratch;
  RefCountedPtr<SigmaStorage>        m_sigma_scratch;

  // Extra storage outside of solver storage
  EBAMRCellData m_scratch1;
  EBAMRCellData m_scratchD;

  Vector<EBAMRCellData*> m_dummy_rte;
  Vector<EBAMRCellData>  m_cache_FR0;

  int m_max_retries;
  int m_p;
  int m_k;
  int m_error_norm;
  int m_min_corr;
  int m_error_idx;

  Real m_accum_cfl;
  Real m_regrid_cfl;
  Real m_max_growth;
  Real m_minCFL;
  Real m_maxCFL;
  Real m_err_thresh;
  Real m_safety;
  Real m_decrease_safe;
  Real m_new_dt; // dt for next time step
  Real m_extrap_dt;

  bool m_semi_implicit;
  bool m_freeze_photoi;
  bool m_extrap_advect;
  bool m_adaptive_dt;
  bool m_print_report;
  bool m_have_err; 

  // Basically only debugging stuff
  bool m_do_advec_src;  
  bool m_do_diffusion;  
  bool m_do_poisson;
  bool m_do_rte;
  bool m_compute_D;
  bool m_compute_v;
  bool m_compute_S;
  bool m_consistent_E;
  bool m_consistent_rte;
  bool m_have_dt_err;
  bool m_profile_steps;

  // Error things
  Vector<Real> m_cdr_error;
  Real m_sigma_error;
  Real m_max_error;
  Real m_pre_error;

  // Stuff for Lobatto intervals, Vandermonde matrix, integration
  // weights, Lobatto nodes etc. 
  Vector<Vector<Real> > m_qmj;
  Vector<Vector<Real> > m_vandermonde;
  Vector<Real> m_nodes;
  Vector<Real> m_tm;
  Vector<Real> m_dtm;

  TimeCode::which_code m_timeCode;

  std::string m_which_nodes;

  // Get functions for storage
  RefCountedPtr<CdrStorage>& get_CdrStorage(const CdrIterator& a_solverit);
  RefCountedPtr<RtStorage>& get_RtStorage(const RtIterator& a_solverit);

  // Allocation functions
  void allocateCdrStorage();
  void allocateFieldStorage();
  void allocateRtStorage();
  void allocateSigmaStorage();

  // Get functions for substep data. 
  Vector<EBAMRCellData*> get_cdr_phik(const int a_m);
  Vector<EBAMRCellData*> get_rte_phik(const int a_m);
  EBAMRIVData& get_sigmak(const int a_m);

  // Routines for intervals and quadrature
  void setupQuadratureNodes(const int a_p);
  void setupUniformNodes(const int a_p);
  void setupLobattoNodes(const int a_p);
  void setupChebyshevNodes(const int a_p);
  void setupQmj(const int a_p);
  void setupSubintervals(const Real a_time, const Real a_dt);
  void quad(EBAMRCellData& a_quad, const Vector<EBAMRCellData>& a_integrand, const int a_m);
  void quad(EBAMRIVData& a_quad, const Vector<EBAMRIVData>& a_integrand, const int a_m);

  // Routines for copying data at tm[0] and tm[p]
  void copyCdrToPhiM0();
  void copySigmaToM0();
  void copy_rte_to_phi_m0();

  void copyPhiPToCdr();
  void copy_phi_p_to_rte();
  void copySigmaPToSigma();
  void copy_FRp_to_FR0();

  void integrateRtTransient(const Real a_dt);
  void integrateRtStationary();

  // These are the standard SDC integration routines
  void sweep(const Real a_dt, const Real a_time, const bool a_lagged_terms);
  void sweep_explicit(const Real a_dt, const Real a_time, const bool a_lagged_terms);
  void sweep_semi_implicit(const Real a_dt, const Real a_time, const bool a_lagged_terms);

  // Explicit divergence routines
  void computeDivD(const int a_m, const bool a_corrector);
  void computeDivF(const int a_m, const bool a_corrector);

  // Advance routine for SDC step
  void substep_rte(const int a_m, const bool a_corrector);
  void substep_cdr(const int a_m, const bool a_corrector);

  // Semi-implicit routines for Poisson equation
  void compute_semi_implicit_mobilities(const int a_m, const bool a_corrector);
  void compute_semi_implicit_rho(const int a_m, const bool a_corrector);
  void set_semi_implicit_permittivities();
  void solve_semi_implicit_poisson(const int a_m);
  void compute_semi_implicit_cdr_velocities(const int a_m, const Real a_time);

  // Advance routines
  void advance_sigma(const int a_m, const bool a_corrector);
  
  // Make integrands quadrature-ready
  void reconcileIntegrands();

  // SDC corrector routines
  void initializeErrors();
  void finalizeErrors();

  // Step size control routines
  void computeNewDt(bool& a_accept_step, const Real a_dt, const int a_num_corrections);
  void adaptiveReport(const Real a_first_dt, const Real a_dt, const Real a_new_dt, const int a_corr, const int a_rej, const Real a_max_err);

  // Poisson, RTE, and diffusion coefficient updates
  void updateField(); // This uses the solver states
  void updateField(const Vector<EBAMRCellData*>& a_densities, const EBAMRIVData& a_sigma);
  void updateDiffusionCoefficients();

  // Store & restore functions
  void storeSolvers();
  void restoreSolvers();

  // Aux functions
  void computeReactionNetwork(const int a_m, const Real a_time, const Real a_dt, const bool a_corrector);
  void compute_E_into_scratch();
  void computeCdrEbStates();
  void computeCdrEbStates(const Vector<EBAMRCellData*>& a_phis);
  void computeCdrDomainStates();
  void computeCdrDomainStates(const Vector<EBAMRCellData*>& a_phis);
  void computeCdrGradients();
  void computeCdrGradients(const Vector<EBAMRCellData*>& a_phis);
  void compute_cdr_fluxes(const Real a_time);
  void compute_cdr_fluxes(const Vector<EBAMRCellData*>& a_phis, const Real a_time);
  void computeCdrDomainFluxes(const Real a_time);
  void computeCdrDomainFluxes(const Vector<EBAMRCellData*>& a_phis, const Real a_time);
  void computeCdrVelo(const Real a_time);
  void computeCdrVelo(const Vector<EBAMRCellData*>& a_phis, const Real a_time);
  void computeDt(Real& a_dt, TimeCode::which_code& a_timeCode);
  void computeSigmaFlux();

  //
  Real getMaxError();
  Real restrict_dt() ;
  Real getMaxNodeDistance();

  // Profiling
  void writeStepProfile(const Real a_dt,
			  const Real a_error,
			  const int  a_substeps,
			  const int  a_corrections,
			  const int  a_rejections);

  // Parse optiosn that are not really used
  void parseNodes();
  void parse_semi_implicit();
  void parseAdaptiveOptions();
  void parseDebugOptions();
  void parseAdvectionOptions();
  void parse_photoi_options();
#include "CD_NamespaceFooter.H"

#endif
