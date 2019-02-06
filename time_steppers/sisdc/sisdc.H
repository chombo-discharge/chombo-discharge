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

  int m_order;
  int m_min_order;
  int m_max_order;
  int m_error_norm;
  int m_num_diff_corr;

  Real m_minCFL;
  Real m_maxCFL;
  Real m_err_thresh;
  Real m_safety;
  Real m_dt_adapt;
  Real m_dt_cfl;

  bool m_adaptive_dt;
  bool m_adaptive_order;

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

  Vector<Real> m_cdr_error;
  Real m_sigma_error;
  Real m_max_error;

  Vector<Real> m_gl_weights;
  Vector<Real> m_gl_nodes;
  Vector<Real> m_tm;
  Vector<Real> m_dtm;

  time_code::which_code m_timecode;

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

  // Routines for intervals etc. 
  void setup_gauss_lobatto(const int a_order);
  void setup_subintervals(const int a_order, const Real a_time, const Real a_dt);

  // SISDC predictor routines
  void predictor(const Real a_time);
  void predictor_advection_reaction(const int a_m);
  void predictor_sigma(const int a_m);
  void predictor_diffusion(const int a_m);

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

  
  Real get_max_error();
  Real restrict_dt() ;
};

#endif
