/*!
  @file   euler_f_storage.cpp
  @brief  Implementation of euler_f_storage.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "euler_f.H"
#include "euler_f_storage.H"

euler_f::cdr_storage::cdr_storage(){

}

euler_f::cdr_storage::cdr_storage(const RefCountedPtr<amr_mesh>& a_amr, const phase::which_phase a_phase, const int a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

euler_f::cdr_storage::~cdr_storage(){
  
}

void euler_f::cdr_storage::allocate_storage(){
  m_amr->allocate(m_scratchIV1,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV2,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV3,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV4,  m_phase, m_ncomp);
}

void euler_f::cdr_storage::deallocate_storage(){
  m_amr->deallocate(m_scratchIV1);
  m_amr->deallocate(m_scratchIV2);
  m_amr->deallocate(m_scratchIV3);
  m_amr->deallocate(m_scratchIV4);
}

euler_f::poisson_storage::poisson_storage(){

}

euler_f::poisson_storage::poisson_storage(const RefCountedPtr<amr_mesh>& a_amr,
					  const phase::which_phase       a_phase,
					  const int                      a_ncomp){
  m_amr   = a_amr;
  m_ncomp = a_ncomp;
  m_phase = a_phase;
}

euler_f::poisson_storage::~poisson_storage(){
  
}

void euler_f::poisson_storage::allocate_storage(){
  m_amr->allocate(m_E_cell, m_phase, SpaceDim);
  m_amr->allocate(m_E_face, m_phase, SpaceDim);
  m_amr->allocate(m_E_eb,   m_phase, SpaceDim);
}

void euler_f::poisson_storage::deallocate_storage(){
  m_amr->deallocate(m_E_cell);
  m_amr->deallocate(m_E_face);
  m_amr->deallocate(m_E_eb);
}

euler_f::rte_storage::rte_storage(){

}

euler_f::rte_storage::rte_storage(const RefCountedPtr<amr_mesh>& a_amr, const phase::which_phase a_phase, const int a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

euler_f::rte_storage::~rte_storage(){
  
}

void euler_f::rte_storage::allocate_storage(){
  m_amr->allocate(m_scratchIV,  m_phase, m_ncomp);
}

void euler_f::rte_storage::deallocate_storage(){
  m_amr->deallocate(m_scratchIV);
}

euler_f::sigma_storage::sigma_storage(){

}

euler_f::sigma_storage::sigma_storage(const RefCountedPtr<amr_mesh>& a_amr, const phase::which_phase a_phase, const int a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

euler_f::sigma_storage::~sigma_storage(){
  
}

void euler_f::sigma_storage::allocate_storage(){
  m_amr->allocate(m_rhs,  m_phase, m_ncomp);
}

void euler_f::sigma_storage::deallocate_storage(){
  m_amr->deallocate(m_rhs);
}
