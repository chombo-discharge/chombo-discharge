/*!
  @file   strang2.H
  @brief  Declaration of a time-adaptive Runge-Kutta methods for advancing plasma equations
  @author Robert Marskar
  @date   Sept. 2018
*/

#ifndef _STRANG2_
#define _STRANG2_

#include "type_definitions.H"
#include "amr_mesh.H"
#include "time_stepper.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"

/*!
  @brief Class for evolving plasma equations with a multirate method using strongly stability preserving
  Runge-Kutta methods. 
*/
class strang2 : public time_stepper {
public:
  
  // Forward declare nested classes. I don't want to clutter this file.
  class cdr_storage;
  class poisson_storage;
  class rte_storage;
  class sigma_storage;

  
public:

  strang2();
  ~strang2();
  
  Real advance(const Real a_dt);
  Real advance_single(const Real a_dt);
  
  void regrid_internals();
  void deallocate_internals();
  
protected:

  Vector<RefCountedPtr<cdr_storage> > m_cdr_scratch;
  Vector<RefCountedPtr<rte_storage> > m_rte_scratch;
  RefCountedPtr<poisson_storage>      m_poisson_scratch;
  RefCountedPtr<sigma_storage>        m_sigma_scratch;

  std::string m_splitting;

  int m_rk_order;
  int m_rk_stages;
  int m_error_norm;

  Real m_minCFL;
  Real m_maxCFL;
  Real m_dt_cfl;
  Real m_err_thresh;
  Real m_safety;
  Real m_dt_adapt;
  Real m_alpha;
  Real m_min_alpha;
  Real m_max_alpha;

  Vector<Real> m_cdr_error;
  Real m_sigma_error;
  Real m_max_error;

  bool m_fixed_order;
  bool m_print_diagno;
  bool m_write_diagno;
  bool m_use_embedded;
  bool m_compute_error;
  bool m_write_errf;
  bool m_do_advec_src;  
  bool m_do_diffusion;  
  bool m_do_poisson;
  bool m_do_rte;
  bool m_compute_D;
  bool m_compute_v;
  bool m_compute_S;
  bool m_adaptive_dt;
  bool m_have_dtf;
  bool m_multirate;
  bool m_consistent_E;
  bool m_consistent_rte;

  time_code::which_code m_timecode;

  RefCountedPtr<cdr_storage>& get_cdr_storage(const cdr_iterator& a_solverit);
  RefCountedPtr<rte_storage>& get_rte_storage(const rte_iterator& a_solverit);

  void allocate_cdr_storage();
  void allocate_poisson_storage();
  void allocate_rte_storage();
  void allocate_sigma_storage();
  
  void cache_solutions();
  void uncache_solutions();
  void copy_error_to_solvers();
  
  Real advance_adaptive(int& substeps, Real& a_dt, const Real a_time, const Real a_dtc);
  Real advance_fixed(const int a_substeps, const Real a_dt);
  Real advance_one_step(const Real a_time, const Real a_dt);

  void advance_rk(const Real a_time, const Real a_dt, const bool a_compute_err);
  void advance_rk22(const Real a_time, const Real a_dt, const bool a_compute_err);
  void advance_rk33(const Real a_time, const Real a_dt, const bool a_compute_err);
  
  void advance_diffusion(const Real a_time, const Real a_dt);
  void advance_euler_diffusion(const Real a_time, const Real a_dt);
  void advance_tga_diffusion(const Real a_time, const Real a_dt);

  void update_poisson();
  void update_rte(const Real a_time);
  void advance_rte_stationary(const Real a_time);

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
  void store_solvers();
  void restore_solvers();

  void compute_errors();
  
  Real get_max_error();
  Real restrict_dt() ;
};

#endif
