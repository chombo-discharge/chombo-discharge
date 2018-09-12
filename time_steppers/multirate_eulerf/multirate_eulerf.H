/*!
  @file   multirate_eulerf.H
  @brief  Declaration of a second order Runge-Kutta methods for advancing plasma equations
  @author Robert Marskar
  @date   Sept. 2018
*/

#ifndef _multirate_eulerf_
#define _multirate_eulerf_

#include "type_definitions.H"
#include "amr_mesh.H"
#include "time_stepper.H"
#include "cdr_iterator.H"
#include "rte_iterator.H"

/*!
  @brief Class for evolving plasma equations with a multirate forward Euler method. 
*/
class multirate_eulerf : public time_stepper {
public:
  
  // Forward declare nested classes. I don't want to clutter this file.
  class cdr_storage;
  class poisson_storage;
  class rte_storage;
  class sigma_storage;

  
public:

  multirate_eulerf();
  ~multirate_eulerf();
  
  Real advance(const Real a_dt);
  
  void regrid_internals();
  void deallocate_internals();
  
protected:

  Vector<RefCountedPtr<cdr_storage> > m_cdr_scratch;
  Vector<RefCountedPtr<rte_storage> > m_rte_scratch;
  RefCountedPtr<poisson_storage>      m_poisson_scratch;
  RefCountedPtr<sigma_storage>        m_sigma_scratch;

  Real m_maxCFL;
  Real m_dt_cfl;

  Vector<Real> diag_times; 

  bool m_diagnostics;
  bool m_do_advec_src;  
  bool m_do_diffusion;  
  bool m_do_poisson;
  bool m_do_rte;


  RefCountedPtr<cdr_storage>& get_cdr_storage(const cdr_iterator& a_solverit);
  RefCountedPtr<rte_storage>& get_rte_storage(const rte_iterator& a_solverit);

  void advance_diffusion(const Real a_dt);
  void advance_multirate_advec_src(const int a_substeps, const Real a_dt);
  void advance_rte_stationary(const Real a_time);
  void allocate_cdr_storage();
  void allocate_poisson_storage();
  void allocate_rte_storage();
  void allocate_sigma_storage();
  void cache_solutions();
  void compute_E_into_scratch();
  void compute_cdr_eb_states();
  void compute_cdr_fluxes(const Real a_time);
  void compute_cdr_velo(const Real a_time);
  void compute_cdr_sources(const Real a_time);
  void compute_dt(Real& a_dt, time_code::which_code& a_timecode);
  void compute_sigma_flux();

  
  Real restrict_dt() ;
  
};

#endif
