/*!
  @file   rk3_storage.cpp
  @brief  Implementation of rk3_storage.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "rk3.H"
#include "rk3_storage.H"

rk3::cdr_storage::cdr_storage(){

}

rk3::cdr_storage::cdr_storage(const RefCountedPtr<amr_mesh>& a_amr, const phase::which_phase a_phase, const int a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

rk3::cdr_storage::~cdr_storage(){
  
}

void rk3::cdr_storage::allocate_storage(){
  m_amr->allocate(m_phi, m_phase, m_ncomp);
  m_amr->allocate(m_k1,  m_phase, m_ncomp);
  m_amr->allocate(m_k2,  m_phase, m_ncomp);
  m_amr->allocate(m_k3,  m_phase, m_ncomp);

  m_amr->allocate(m_scratchIV1,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV2,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV3,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV4,  m_phase, m_ncomp);

  m_amr->allocate(m_scratch, m_phase, m_ncomp); 
}

rk3::poisson_storage::poisson_storage(){

}

rk3::poisson_storage::poisson_storage(const RefCountedPtr<amr_mesh>& a_amr, const phase::which_phase a_phase, const int a_ncomp){
  m_amr   = a_amr;
  m_ncomp = a_ncomp;
  m_phase = a_phase;
}

rk3::poisson_storage::~poisson_storage(){
  
}

void rk3::poisson_storage::allocate_storage(){
  m_amr->allocate(m_phi,    m_ncomp);
  
  m_amr->allocate(m_E_cell, m_phase, SpaceDim);
  m_amr->allocate(m_E_face, m_phase, SpaceDim);
  m_amr->allocate(m_E_eb,   m_phase, SpaceDim);

  m_amr->allocate(m_scratch_phi, m_ncomp);
  m_amr->allocate(m_scratch_E,   m_phase, SpaceDim);
}

rk3::rte_storage::rte_storage(){

}

rk3::rte_storage::rte_storage(const RefCountedPtr<amr_mesh>& a_amr, const phase::which_phase a_phase, const int a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

rk3::rte_storage::~rte_storage(){
  
}

void rk3::rte_storage::allocate_storage(){
  m_amr->allocate(m_phi,        m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV,  m_phase, m_ncomp);
}

rk3::sigma_storage::sigma_storage(){

}

rk3::sigma_storage::sigma_storage(const RefCountedPtr<amr_mesh>& a_amr, const phase::which_phase a_phase, const int a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

rk3::sigma_storage::~sigma_storage(){
}

void rk3::sigma_storage::allocate_storage(){
  m_amr->allocate(m_phi, m_phase, m_ncomp);
  m_amr->allocate(m_k1,  m_phase, m_ncomp);
  m_amr->allocate(m_k2,  m_phase, m_ncomp);
  m_amr->allocate(m_k3,  m_phase, m_ncomp);
}
