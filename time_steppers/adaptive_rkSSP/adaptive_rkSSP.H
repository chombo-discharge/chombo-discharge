/*!
  @file   adaptive_rkSSP.H
  @brief  Declaration of a time-adaptive Runge-Kutta methods for advancing plasma equations
  @author Robert Marskar
  @date   Sept. 2018
*/

#ifndef _MULTIRATE_RKSSP_
#define _MULTIRATE_RKSSP_

#include "type_definitions.H"
#include "amr_mesh.H"
#include "time_stepper.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"

/*!
  @brief Class for evolving plasma equations with a multirate method using strongly stability preserving
  Runge-Kutta methods. 
*/
class adaptive_rkSSP : public time_stepper {
public:
  
  // Forward declare nested classes. I don't want to clutter this file.
  class cdr_storage;
  class poisson_storage;
  class rte_storage;
  class sigma_storage;

  
public:

  adaptive_rkSSP();
  ~adaptive_rkSSP();
  
  Real advance(const Real a_dt);
  Real advance_single(const Real a_dt);
  
  void regrid_internals();
  void deallocate_internals();
  
protected:

  Vector<RefCountedPtr<cdr_storage> > m_cdr_scratch;
  Vector<RefCountedPtr<rte_storage> > m_rte_scratch;
  RefCountedPtr<poisson_storage>      m_poisson_scratch;
  RefCountedPtr<sigma_storage>        m_sigma_scratch;

  int m_order;
  int m_error_norm;

  Real m_minCFL;
  Real m_maxCFL;
  Real m_dt_cfl;
  Real m_err_thresh;
  Real m_safety;
  Real m_dt_adapt;

  Vector<Real> m_cdr_error;
  Real m_sigma_error;
  Real m_max_error;

  bool m_print_diagno;
  bool m_write_diagno;
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
  bool m_consistent_rte;

  time_code::which_code m_timecode;

  RefCountedPtr<cdr_storage>& get_cdr_storage(const cdr_iterator& a_solverit);
  RefCountedPtr<rte_storage>& get_rte_storage(const rte_iterator& a_solverit);

  void advance_fixed(const int a_substeps, const Real a_dt);
  void advance_adaptive(int& substeps, Real& a_dt, const Real a_time, const Real a_dtc);
  void advance_rk1(const Real a_time, const Real a_dt);
  void advance_rk2(const Real a_time, const Real a_dt);
  void advance_rk3(const Real a_time, const Real a_dt);
  void advance_diffusion(const Real a_time, const Real a_dt);
  void advance_euler_diffusion(const Real a_time, const Real a_dt);
  void advance_tga_diffusion(const Real a_time, const Real a_dt);
  void advance_rte_stationary(const Real a_time);
  void allocate_cdr_storage();
  void allocate_poisson_storage();
  void allocate_rte_storage();
  void allocate_sigma_storage();
  void cache_solutions();
  void uncache_solutions();
  void compute_cdr_gradients();
  void compute_errors();
  void compute_E_into_scratch();
  void compute_cdr_eb_states();
  void compute_cdr_domain_states();
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
  
  Real get_max_error();
  Real restrict_dt() ;
};

#endif
