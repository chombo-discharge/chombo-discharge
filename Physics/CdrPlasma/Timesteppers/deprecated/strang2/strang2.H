/*!
  @file   strang2.H
  @brief  Declaration of a time-adaptive Runge-Kutta methods for advancing plasma equations
  @author Robert Marskar
  @date   Sept. 2018
*/

#ifndef _STRANG2_
#define _STRANG2_

#include "CD_AmrMesh.H"
#include "CD_TimeStepper.H"
#include <CD_CdrIterator.H>
#include <CD_RtIterator.H>

/*!
  @brief Class for evolving plasma equations with a multirate method using strongly stability preserving
  Runge-Kutta methods. 
*/
class strang2 : public TimeStepper {
public:
  
  // Forward declare nested classes. I don't want to clutter this file.
  class CdrStorage;
  class FieldStorage;
  class RtStorage;
  class SigmaStorage;

  
public:

  strang2();
  ~strang2();
  
  Real advance(const Real a_dt);
  Real advance_single(const Real a_dt);
  
  void regridInternals();
  void deallocateInternals();
  
protected:

  Vector<RefCountedPtr<CdrStorage> > m_cdr_scratch;
  Vector<RefCountedPtr<RtStorage> > m_rte_scratch;
  RefCountedPtr<FieldStorage>      m_fieldSolver_scratch;
  RefCountedPtr<SigmaStorage>        m_sigma_scratch;

  std::string m_splitting;

  int m_rk_order;
  int m_rk_stages;
  int m_error_norm;
  int m_auto_min_stages;
  int m_auto_max_stages;

  Real m_minCFL;
  Real m_maxCFL;
  Real m_auto_safety_cfl;
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

  bool m_auto_stages;
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

  TimeCode::which_code m_timeCode;

  RefCountedPtr<CdrStorage>& get_CdrStorage(const CdrIterator& a_solverit);
  RefCountedPtr<RtStorage>& get_RtStorage(const RtIterator& a_solverit);

  void allocateCdrStorage();
  void allocateFieldStorage();
  void allocateRtStorage();
  void allocateSigmaStorage();
  
  void backup_solutions();
  void revert_backup();
  void copy_solvers_to_error();
  void copy_error_to_solvers();
  void copy_solvers_to_cache();
  void copy_cache_to_solvers();
  
  Real advance_adaptive(int& substeps, Real& a_dt, const Real a_time, const Real a_dtc);
  Real advance_fixed(const int a_substeps, const Real a_dt);
  Real advance_one_step(const Real a_time, const Real a_dt);

  void advance_rk(const Real a_time,   const Real a_dt);
  void advance_rkN2(const Real a_time, const Real a_dt, const int a_stages);
  void advance_rk33(const Real a_time, const Real a_dt);
  void advance_rk43(const Real a_time, const Real a_dt);
  void advance_rk53(const Real a_time, const Real a_dt);
  void advance_rk54(const Real a_time, const Real a_dt);
  
  void advance_diffusion(const Real a_time, const Real a_dt);
  void advanceEuler_diffusion(const Real a_time, const Real a_dt);
  void advanceTGA_diffusion(const Real a_time, const Real a_dt);

  void updateField();
  void update_rte(const Real a_time);
  void advance_rte_stationary(const Real a_time);


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
  void compute_cdr_sources(const Real a_time);
  void compute_cdr_sources(const Vector<EBAMRCellData*>& a_phis, const Real a_time);
  void computeDt(Real& a_dt, TimeCode::which_code& a_timeCode);
  void computeSigmaFlux();
  void storeSolvers();
  void restoreSolvers();

  void compute_errors();
  
  Real getMaxError();
  Real restrict_dt() ;
#include "CD_NamespaceFooter.H"

#endif
