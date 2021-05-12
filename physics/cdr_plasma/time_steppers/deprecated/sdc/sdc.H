/*!
  @file   sdc.H
  @brief  Declaration of a semi-implicit spectral deferred correction class
  @author Robert Marskar
  @date   Feb. 2019
*/

#ifndef _SDC_
#define _SDC_

#include "CD_AmrMesh.H"
#include "time_stepper.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"

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
class sdc : public time_stepper {
public:
  
  // Forward declare nested classes. I don't want to clutter this file.
  class cdr_storage;
  class poisson_storage;
  class rte_storage;
  class sigma_storage;
  
public:

  sdc();
  ~sdc();
  
  Real advance(const Real a_dt);

  void init();
  void cache_internals();
  void regrid_internals(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level);
  void allocate_internals();
  void deallocate_internals();
  void parseOptions();
  bool need_to_regrid();

  Vector<EBAMRCellData*> get_cdr_errors();
  
protected:

  Vector<RefCountedPtr<cdr_storage> > m_cdr_scratch;
  Vector<RefCountedPtr<rte_storage> > m_rte_scratch;
  RefCountedPtr<poisson_storage>      m_fieldSolver_scratch;
  RefCountedPtr<sigma_storage>        m_sigma_scratch;

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

  time_code::which_code m_timecode;

  std::string m_which_nodes;

  // Get functions for storage
  RefCountedPtr<cdr_storage>& get_cdr_storage(const cdr_iterator& a_solverit);
  RefCountedPtr<rte_storage>& get_rte_storage(const rte_iterator& a_solverit);

  // Allocation functions
  void allocate_cdr_storage();
  void allocate_poisson_storage();
  void allocate_rte_storage();
  void allocate_sigma_storage();

  // Get functions for substep data. 
  Vector<EBAMRCellData*> get_cdr_phik(const int a_m);
  Vector<EBAMRCellData*> get_rte_phik(const int a_m);
  EBAMRIVData& get_sigmak(const int a_m);

  // Routines for intervals and quadrature
  void setup_quadrature_nodes(const int a_p);
  void setup_uniform_nodes(const int a_p);
  void setup_lobatto_nodes(const int a_p);
  void setup_chebyshev_nodes(const int a_p);
  void setup_qmj(const int a_p);
  void setup_subintervals(const Real a_time, const Real a_dt);
  void quad(EBAMRCellData& a_quad, const Vector<EBAMRCellData>& a_integrand, const int a_m);
  void quad(EBAMRIVData& a_quad, const Vector<EBAMRIVData>& a_integrand, const int a_m);

  // Routines for copying data at tm[0] and tm[p]
  void copy_cdr_to_phi_m0();
  void copy_sigma_to_sigma_m0();
  void copy_rte_to_phi_m0();

  void copy_phi_p_to_cdr();
  void copy_phi_p_to_rte();
  void copy_sigma_p_to_sigma();
  void copy_FRp_to_FR0();

  void integrate_rte_transient(const Real a_dt);
  void integrate_rte_stationary();

  // These are the standard SDC integration routines
  void sweep(const Real a_dt, const Real a_time, const bool a_lagged_terms);
  void sweep_explicit(const Real a_dt, const Real a_time, const bool a_lagged_terms);
  void sweep_semi_implicit(const Real a_dt, const Real a_time, const bool a_lagged_terms);

  // Explicit divergence routines
  void compute_divD(const int a_m, const bool a_corrector);
  void compute_divF(const int a_m, const bool a_corrector);

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
  void reconcile_integrands();

  // SDC corrector routines
  void initialize_errors();
  void finalize_errors();

  // Step size control routines
  void compute_new_dt(bool& a_accept_step, const Real a_dt, const int a_num_corrections);
  void adaptive_report(const Real a_first_dt, const Real a_dt, const Real a_new_dt, const int a_corr, const int a_rej, const Real a_max_err);

  // Poisson, RTE, and diffusion coefficient updates
  void update_poisson(); // This uses the solver states
  void update_poisson(const Vector<EBAMRCellData*>& a_densities, const EBAMRIVData& a_sigma);
  void update_diffusion_coefficients();

  // Store & restore functions
  void store_solvers();
  void restore_solvers();

  // Aux functions
  void compute_reaction_network(const int a_m, const Real a_time, const Real a_dt, const bool a_corrector);
  void compute_E_into_scratch();
  void compute_cdr_eb_states();
  void compute_cdr_eb_states(const Vector<EBAMRCellData*>& a_states);
  void compute_cdr_domain_states();
  void compute_cdr_domain_states(const Vector<EBAMRCellData*>& a_states);
  void compute_cdr_gradients();
  void compute_cdr_gradients(const Vector<EBAMRCellData*>& a_states);
  void compute_cdr_fluxes(const Real a_time);
  void compute_cdr_fluxes(const Vector<EBAMRCellData*>& a_states, const Real a_time);
  void compute_cdr_domain_fluxes(const Real a_time);
  void compute_cdr_domain_fluxes(const Vector<EBAMRCellData*>& a_states, const Real a_time);
  void compute_cdr_velo(const Real a_time);
  void compute_cdr_velo(const Vector<EBAMRCellData*>& a_states, const Real a_time);
  void compute_dt(Real& a_dt, time_code::which_code& a_timecode);
  void compute_sigma_flux();

  //
  Real get_max_error();
  Real restrict_dt() ;
  Real get_max_node_distance();

  // Profiling
  void write_step_profile(const Real a_dt,
			  const Real a_error,
			  const int  a_substeps,
			  const int  a_corrections,
			  const int  a_rejections);

  // Parse optiosn that are not really used
  void parse_nodes();
  void parse_semi_implicit();
  void parse_adaptive_options();
  void parse_debug_options();
  void parse_advection_options();
  void parse_photoi_options();
#include "CD_NamespaceFooter.H"

#endif
