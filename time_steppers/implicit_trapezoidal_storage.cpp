/*!
  @file   implicit_trapezoidal_storage.H
  @brief  Declaration of scratch storage for implicit_trapezoidal
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "implicit_trapezoidal_storage.H"

implicit_trapezoidal::cdr_storage::cdr_storage(){

}

implicit_trapezoidal::cdr_storage::cdr_storage(const RefCountedPtr<amr_mesh>& a_amr,
					       const phase::which_phase       a_phase,
					       const int                      a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

implicit_trapezoidal::cdr_storage::~cdr_storage(){
  
}

implicit_trapezoidal::rte_storage::rte_storage(){

}

implicit_trapezoidal::rte_storage::rte_storage(const RefCountedPtr<amr_mesh>& a_amr,
					       const phase::which_phase       a_phase,
					       const int                      a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

implicit_trapezoidal::rte_storage::~rte_storage(){
  
}

implicit_trapezoidal::poisson_storage::poisson_storage(){

}

implicit_trapezoidal::poisson_storage::poisson_storage(const RefCountedPtr<amr_mesh>& a_amr,
						       const phase::which_phase       a_phase,
						       const int                      a_ncomp){
  m_amr   = a_amr;
  m_ncomp = a_ncomp;
  m_phase = a_phase;
}

implicit_trapezoidal::poisson_storage::~poisson_storage(){
  
}

implicit_trapezoidal::sigma_storage::sigma_storage(){

}

implicit_trapezoidal::sigma_storage::sigma_storage(const RefCountedPtr<amr_mesh>& a_amr,
						   const phase::which_phase       a_phase,
						   const int                      a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

implicit_trapezoidal::sigma_storage::~sigma_storage(){
  
}

void implicit_trapezoidal::cdr_storage::allocate_storage(){
  m_amr->allocate(m_cache, m_phase, m_ncomp);
}

void implicit_trapezoidal::cdr_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
}

void implicit_trapezoidal::rte_storage::allocate_storage(){
  m_amr->allocate(m_cache, m_phase, m_ncomp);
}

void implicit_trapezoidal::rte_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
}

void implicit_trapezoidal::poisson_storage::allocate_storage(){
  m_amr->allocate(m_cache, m_ncomp);
}

void implicit_trapezoidal::poisson_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
}

void implicit_trapezoidal::sigma_storage::allocate_storage(){
  m_amr->allocate(m_cache, m_phase, m_ncomp);  
}

void implicit_trapezoidal::sigma_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);  
}
