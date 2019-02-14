/*!
  @file   sisdc.H
  @brief  Declaration of a semi-implicit spectral deferred correction class
  @author Robert Marskar
  @date   Feb. 2019
*/

#ifndef _SISDC_
#define _SISDC_

#include "type_definitions.H"
#include "amr_mesh.H"
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
  @brief Class for evolving plasma equations with a multirate method using strongly stability preserving
  Runge-Kutta methods. 
*/
class sisdc : public time_stepper {
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
  
  void regrid_internals();
  void deallocate_internals();
  
protected:

  Vector<RefCountedPtr<cdr_storage> > m_cdr_scratch;
  Vector<RefCountedPtr<rte_storage> > m_rte_scratch;
  RefCountedPtr<poisson_storage>      m_poisson_scratch;
  RefCountedPtr<sigma_storage>        m_sigma_scratch;

  int m_p;
  int m_k;
  int m_error_norm;
  int m_num_diff_corr;

  Real m_minCFL;
  Real m_maxCFL;
  Real m_err_thresh;
  Real m_safety;
  Real m_dt_adapt;
  Real m_dt_cfl;

  bool m_adaptive_dt;

  // Basically6 only debugging stuff
  bool m_do_advec_src;  
  bool m_do_diffusion;  
  bool m_do_poisson;
  bool m_do_rte;
  bool m_compute_D;
  bool m_compute_v;
  bool m_compute_S;
  bool m_have_dtf;
  bool m_consistent_E;
  bool m_consistent_rte;
  bool m_print_diagno;
  bool m_write_diagno;
  bool m_strong_diffu;

  // Error things
  Vector<Real> m_cdr_error;
  Real m_sigma_error;
  Real m_max_error;

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

  // Store and reset functions
  void store_previous_solutions();
  void restore_previous_solutions();

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
  void gl_quad(EBAMRCellData& a_quad, const Vector<EBAMRCellData>& a_integrand, const int a_m);
  void gl_quad(EBAMRIVData& a_quad, const Vector<EBAMRIVData>& a_integrand, const int a_m);

  // Routines for copying data at tm[0] and tm[p]
  void copy_cdr_to_phi_m0();
  void copy_sigma_to_sigma_m0();
  void copy_phi_p_to_cdr();
  void copy_sigma_p_to_sigma();

  // SISDC predictor routines
  void predictor_compute_FD_0();
  void predictor(const Real a_time);
  void predictor_advection_reaction(const int a_m);
  void predictor_sigma(const int a_m);
  void predictor_diffusion(const int a_m);
  void predictor_diffusion_onestep(const int a_m);
  void predictor_diffusion_build_FD(const int a_m);

  // SISDC corrector routines
  void corrector_initialize_errors();
  void corrector_compute_FAR_p();
  void corrector_reconcile_gl_integrands();
  void corrector(const Real a_time, const Real a_dt);
  void corrector_advection_reaction(const int a_m, const Real a_dt);
  void corrector_sigma(const int a_m, const Real a_dt);
  void corrector_diffusion(const int a_m);
  void corrector_diffusion_onestep(const int a_m);
  void corrector_diffusion_build_FD(const int a_m);
  void corrector_finalize_errors();

  // Step size control routines
  void compute_new_dt(const Real a_dt, const int a_num_corrections);

  // Poisson, RTE, and diffusion coefficient updates
  void update_poisson(); // This uses the solver states
  void update_poisson(const Vector<EBAMRCellData*>& a_densities, const EBAMRIVData& a_sigma);
  void update_rte(const Real a_time); // This uses the solver states, make sure that is what you want
  void update_rte(const Vector<EBAMRCellData*>& a_densities, const Real a_time);
  void update_diffusion_coefficients();

  // Aux functions
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
  void compute_cdr_sources(const Real a_time);
  void compute_cdr_sources(const Vector<EBAMRCellData*>& a_states, const Real a_time);
  void compute_dt(Real& a_dt, time_code::which_code& a_timecode);
  void compute_sigma_flux();

  //
  Real get_max_error();
  Real restrict_dt() ;
  Real get_max_lobatto_distance();
};

#endif
