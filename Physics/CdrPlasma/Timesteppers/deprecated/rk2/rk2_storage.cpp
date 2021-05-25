/*!
  @file   rk2_storage.cpp
  @brief  Implementation of rk2_storage.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "rk2.H"
#include "rk2_storage.H"

rk2::CdrStorage::CdrStorage(){

}

rk2::CdrStorage::CdrStorage(const RefCountedPtr<AmrMesh>& a_amr, const phase::which_phase a_phase, const int a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

rk2::CdrStorage::~CdrStorage(){
  
}

void rk2::CdrStorage::allocate_storage(){
  m_amr->allocate(m_phi, m_phase, m_ncomp);
  m_amr->allocate(m_k1,  m_phase, m_ncomp);
  m_amr->allocate(m_k2,  m_phase, m_ncomp);

  m_amr->allocate(m_scratchIV1,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV2,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV3,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV4,  m_phase, m_ncomp);

  m_amr->allocate(m_scratch, m_phase, m_ncomp); 
}

void rk2::CdrStorage::deallocate_storage(){
  m_amr->deallocate(m_phi);
  m_amr->deallocate(m_k1);
  m_amr->deallocate(m_k2);

  m_amr->deallocate(m_scratchIV1);
  m_amr->deallocate(m_scratchIV2);
  m_amr->deallocate(m_scratchIV3);
  m_amr->deallocate(m_scratchIV4);

  m_amr->deallocate(m_scratch);
}

rk2::FieldStorage::FieldStorage(){

}

rk2::FieldStorage::FieldStorage(const RefCountedPtr<AmrMesh>& a_amr, const phase::which_phase a_phase, const int a_ncomp){
  m_amr   = a_amr;
  m_ncomp = a_ncomp;
  m_phase = a_phase;
}

rk2::FieldStorage::~FieldStorage(){
  
}

void rk2::FieldStorage::allocate_storage(){
  m_amr->allocate(m_phi,    m_ncomp);
  
  m_amr->allocate(m_E_cell, m_phase, SpaceDim);
  m_amr->allocate(m_E_face, m_phase, SpaceDim);
  m_amr->allocate(m_E_eb,   m_phase, SpaceDim);

  m_amr->allocate(m_scratch_phi, m_ncomp);
  m_amr->allocate(m_scratch_E,   m_phase, SpaceDim);
}

void rk2::FieldStorage::deallocate_storage(){
  m_amr->deallocate(m_phi);
  
  m_amr->deallocate(m_E_cell);
  m_amr->deallocate(m_E_face);
  m_amr->deallocate(m_E_eb);

  m_amr->deallocate(m_scratch_phi);
  m_amr->deallocate(m_scratch_E);
}

rk2::RtStorage::RtStorage(){

}

rk2::RtStorage::RtStorage(const RefCountedPtr<AmrMesh>& a_amr, const phase::which_phase a_phase, const int a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

rk2::RtStorage::~RtStorage(){
  
}

void rk2::RtStorage::allocate_storage(){
  m_amr->allocate(m_phi,        m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV,  m_phase, m_ncomp);
}

void rk2::RtStorage::deallocate_storage(){
  m_amr->deallocate(m_phi);
  m_amr->deallocate(m_scratchIV);
}

rk2::SigmaStorage::SigmaStorage(){

}

rk2::SigmaStorage::SigmaStorage(const RefCountedPtr<AmrMesh>& a_amr, const phase::which_phase a_phase, const int a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

rk2::SigmaStorage::~SigmaStorage(){
}

void rk2::SigmaStorage::allocate_storage(){
  m_amr->allocate(m_phi, m_phase, m_ncomp);
  m_amr->allocate(m_k1,  m_phase, m_ncomp);
  m_amr->allocate(m_k2,  m_phase, m_ncomp);
}

void rk2::SigmaStorage::deallocate_storage(){
  m_amr->deallocate(m_phi);
  m_amr->deallocate(m_k1);
  m_amr->deallocate(m_k2);
#include "CD_NamespaceFooter.H"
