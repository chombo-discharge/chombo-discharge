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
  m_amr->allocate(m_cache,       m_phase, m_ncomp);
  m_amr->allocate(m_scratch1,    m_phase, m_ncomp);
  m_amr->allocate(m_scratch2,    m_phase, m_ncomp);
  m_amr->allocate(m_scratch3,    m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV1,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV2,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV3,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV4,  m_phase, m_ncomp);
}

void implicit_trapezoidal::cdr_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
  m_amr->deallocate(m_scratch1);
  m_amr->deallocate(m_scratch2);
  m_amr->deallocate(m_scratch3);
  m_amr->deallocate(m_scratchIV1);
  m_amr->deallocate(m_scratchIV2);
  m_amr->deallocate(m_scratchIV3);
  m_amr->deallocate(m_scratchIV4);
}

void implicit_trapezoidal::rte_storage::allocate_storage(){
  m_amr->allocate(m_cache,      m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV,  m_phase, m_ncomp);
}

void implicit_trapezoidal::rte_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
  m_amr->deallocate(m_scratchIV);
}

void implicit_trapezoidal::poisson_storage::allocate_storage(){
  m_amr->allocate(m_cache,  m_ncomp);
  m_amr->allocate(m_E_cell, m_phase, SpaceDim);
  m_amr->allocate(m_E_face, m_phase, SpaceDim);
  m_amr->allocate(m_E_eb,   m_phase, SpaceDim);
}

void implicit_trapezoidal::poisson_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
  m_amr->deallocate(m_E_cell);
  m_amr->deallocate(m_E_face);
  m_amr->deallocate(m_E_eb);
}

void implicit_trapezoidal::sigma_storage::allocate_storage(){
  m_amr->allocate(m_cache,    m_phase, m_ncomp);
  m_amr->allocate(m_scratch1, m_phase, m_ncomp);  
}

void implicit_trapezoidal::sigma_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
  m_amr->deallocate(m_scratch1);  
}
