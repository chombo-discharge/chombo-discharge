/*!
  @file   sisdc.H
  @brief  Declaration of a semi-implicit spectral deferred correction class
  @author Robert Marskar
  @date   Feb. 2019
*/

#ifndef _SISDC_
#define _SISDC_

#include "CD_AmrMesh.H"
#include "CD_TimeStepper.H"
#include <CD_CdrIterator.H>
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
class sisdc : public TimeStepper {
public:
  
  // Forward declare nested classes. I don't want to clutter this file.
  class cdr_storage;
  class poisson_storage;
  class rte_storage;
  class sigma_storage;

  
public:

  sisdc();
  ~sisdc();
  
  Real advance(const Real a_dt);

  void init_source_terms();
  void regridInternals();
  void deallocateInternals();
  void parseOptions();
  bool needToRegrid();

  Vector<EBAMRCellData*> get_cdr_errors();
  
protected:

  Vector<RefCountedPtr<cdr_storage> > m_cdr_scratch;
  Vector<RefCountedPtr<rte_storage> > m_rte_scratch;
  RefCountedPtr<poisson_storage>      m_fieldSolver_scratch;
  RefCountedPtr<sigma_storage>        m_sigma_scratch;

  // Extra storage outside of solver storage
  EBAMRCellData scratch1;
  EBAMRCellData scratchD;

  int m_max_retries;
  int m_p;
  int m_k;
  int m_error_norm;
  int m_num_diff_corr;
  int m_min_corr;
  int m_error_idx;

  Real m_accum_cfl;
  Real m_regrid_cfl;
  Real m_max_growth;
  Real m_minCFL;
  Real m_maxCFL;
  Real m_cycleCFL;
  Real m_err_thresh;
  Real m_safety;
  Real m_decrease_safe;
  Real m_new_dt; // dt for next time step
  Real m_extrap_dt;
  Real m_min_cycle_cfl; // Smallest allowed global CFL for subcycling advances
  Real m_max_cycle_cfl; // Largest allowed global CFL for subcycling advances

  bool m_extrap_advect;
  bool m_cycle_sources; // Subcycle with source terms
  bool m_adaptive_dt;
  bool m_print_report;
  bool m_have_err; 
  bool m_useTGA; // Second order diffusion
  bool m_use_AP; // Use the AP scheme
  

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
  bool m_strong_diffu;
  bool m_have_dt_err;
  bool m_profile_steps;
  bool m_optimal_subcycling;
  bool m_multistep;

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
  RefCountedPtr<cdr_storage>& get_cdr_storage(const CdrIterator& a_solverit);
  RefCountedPtr<rte_storage>& get_rte_storage(const rte_iterator& a_solverit);

  // Allocation functions
  void allocate_cdr_storage();
  void allocate_poisson_storage();
  void allocate_rte_storage();
  void allocate_sigma_storage();

  // Get functions for substep data. 
  Vector<EBAMRCellData*> get_cdr_phik(const int a_m);
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

  void copy_phi_p_to_cdr();
  void copy_sigma_p_to_sigma();

  // Integrate RTE (only time dependent)
  void integrate_rte(const Real a_dt, const int a_m, const bool a_lagged_terms);

  // These are the standard SISDC integration routines
  void integrate(const Real a_dt, const Real a_time, const bool a_lagged_terms);
  void integrate_advection_reaction(const Real a_dt, const int a_m, const bool a_lagged_terms);
  void integrate_advection_nosubcycle(const Real a_dt, const int a_m, const bool a_lagged_terms);
  void integrate_advection_multistep(const Real a_dt, const int a_m, const bool a_lagged_terms);
  void integrate_advection_subcycle(const Real a_dt, const int a_m, const bool a_lagged_terms);
  void integrate_diffusion(const Real a_dt, const int a_m, const bool a_lagged_terms);

  // Here are the ''asymptotic-preserving'' AP methods
  void integrate_AP(const Real a_dt, const Real a_time, const bool a_lagged_terms);
  void integrate_AP_advection_reaction(const Real a_dt, const int a_m, const bool a_lagged_terms, const bool a_predictor);

  // Make integrands quadrature-ready
  void reconcile_integrands();

  // SISDC predictor routines
  void compute_FD_0();

  // SISDC corrector routines
  void initialize_errors();
#if 0 // Deprecated?
  void corrector_reconcile_gl_integrands();
  void corrector(const Real a_time, const Real a_dt);
  void corrector_advection_reaction(const int a_m, const Real a_dt);
  void corrector_advection_reaction_subcycle(const int a_m, const Real a_dt);
  void corrector_advection_reaction_nosubcycle(const int a_m, const Real a_dt);
  void corrector_diffusion(const int a_m);
  void corrector_diffusion_onestep(const int a_m);
  void corrector_diffusion_build_FD(const int a_m);
#endif
  void finalize_errors();

  // Step size control routines
  void compute_new_dt(bool& a_accept_step, const Real a_dt, const int a_num_corrections);
  void adaptive_report(const Real a_first_dt, const Real a_dt, const Real a_new_dt, const int a_corr, const int a_rej, const Real a_max_err);

  // Poisson, RTE, and diffusion coefficient updates
  void update_poisson(); // This uses the solver states
  void update_poisson(const Vector<EBAMRCellData*>& a_densities, const EBAMRIVData& a_sigma);
  void update_stationary_rte(const Real a_time); // This uses the solver states, make sure that is what you want
  void update_stationary_rte(const Vector<EBAMRCellData*>& a_densities, const Real a_time);
  void update_diffusion_coefficients();

  // Store & restore functions
  void store_solvers();
  void restore_solvers();

  // Aux functions
  void compute_E_into_scratch();
  void compute_cdr_eb_states();
  void compute_cdr_eb_states(const Vector<EBAMRCellData*>& a_phis);
  void compute_cdr_domain_states();
  void compute_cdr_domain_states(const Vector<EBAMRCellData*>& a_phis);
  void compute_cdr_gradients();
  void compute_cdr_gradients(const Vector<EBAMRCellData*>& a_phis);
  void compute_cdr_fluxes(const Real a_time);
  void compute_cdr_fluxes(const Vector<EBAMRCellData*>& a_phis, const Real a_time);
  void compute_cdr_domain_fluxes(const Real a_time);
  void compute_cdr_domain_fluxes(const Vector<EBAMRCellData*>& a_phis, const Real a_time);
  void compute_cdr_velo(const Real a_time);
  void compute_cdr_velo(const Vector<EBAMRCellData*>& a_phis, const Real a_time);
  void compute_cdr_sources(const Real a_time);
  void compute_cdr_sources(const Vector<EBAMRCellData*>& a_phis, const Real a_time);
  void computeDt(Real& a_dt, TimeCode::which_code& a_timeCode);
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

  // All the routines required for the subcycling advection routines
  void reset_finer_flux_registers_level(const int a_lvl, const int a_coarsest_level, const int a_finestLevel);
  void reset_redist_registers_level(const int a_lvl, const int a_coarsest_level, const int a_finestLevel);
  void update_redist_register(const LevelData<BaseIVFAB<Real> >& a_massDifference, const int a_solver, const int a_lvl);
  void update_flux_registers(LevelData<EBFluxFAB>& a_flux,
			     const int             a_solver,
			     const int             a_lvl,
			     const int             a_coarsest_level,
			     const int             a_finestLevel,
			     const Real            a_dt);
  void update_coarse_fine_register(const LevelData<BaseIVFAB<Real> >& a_massDifference,
				   const int a_solver,
				   const int a_lvl,
				   const int a_coarsest_level,
				   const int a_finestLevel);
  void reflux_level(EBAMRCellData& a_phi,
		    const int      a_solver,
		    const int      a_lvl,
		    const int      a_coarsest_level,
		    const int      a_finestLevel,
		    const Real     a_scale);
  void redist_level(LevelData<EBCellFAB>&       a_phi,
		    const int                   a_solver,   
		    const LevelData<EBCellFAB>& a_weights,
		    const int                   a_lvl);


  
  void subcycle_advect_amr(EBAMRFluxData& a_flux,
			   EBAMRFluxData& a_facePhi,
			   EBAMRCellData& a_divF_c,
			   EBAMRCellData& a_weights,
			   EBAMRIVData&   a_nonConservativeDivergence,
			   EBAMRIVData&   a_massDifference,
			   Vector<Real>&  a_tnew,
			   Vector<Real>&  a_told,
			   const int      a_m,
			   const int      a_lvl,
			   const int      a_coarsest_lvl,
			   const int      a_finestLevel,
			   const Real     a_dt);

  void subcycle_compute_advection_velocities();
  void subcycle_copy_states(const int a_m);
  void subcycle_copy_current_to_old_states(const int a_m, const int a_lvl);
  void subcycle_update_transport_bc(const int a_m, const int a_lvl, const Real a_time);
  void subcycle_update_sources(const int a_m, const int a_lvl, const Real a_time); 
  void subcycle_sync_levels(const int a_m, const int a_lvl, const int a_coarsest_level, const int a_finestLevel);
  
  void subcycle_integrate_level(LevelData<EBFluxFAB>&        a_flux,
				LevelData<EBFluxFAB>&        a_facePhi,
				LevelData<EBCellFAB>&        a_divF_c,
				LevelData<EBCellFAB>&        a_weights,
				LevelData<BaseIVFAB<Real> >& a_massDifference,
				LevelData<BaseIVFAB<Real> >& a_nonConservativeDivergence,
				const int                    a_m,
				const int                    a_lvl,
				const int                    a_coarsest_level,
				const int                    a_finestLevel,
				const Real                   a_coar_time_old,
				const Real                   a_coar_time_new,
				const Real                   a_time,
				const Real                   a_dt);

  // Parse optiosn that are not really used
  void parse_nodes();
  void parse_AP();
  void parse_diffusion_coupling();
  void parse_adaptive_options();
  void parse_subcycle_options();
  void parse_debug_options();
  void parse_advection_options();
#include "CD_NamespaceFooter.H"

#endif
