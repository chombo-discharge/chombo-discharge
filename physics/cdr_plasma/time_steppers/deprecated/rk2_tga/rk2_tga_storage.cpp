/*!
  @file   rk2_tga_storage.cpp
  @brief  Implementation of rk2_tga_storage.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "rk2_tga.H"
#include "rk2_tga_storage.H"

rk2_tga::cdr_storage::cdr_storage(){

}

rk2_tga::cdr_storage::cdr_storage(const RefCountedPtr<amr_mesh>& a_amr,
				  const phase::which_phase       a_phase,
				  const int                      a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

rk2_tga::cdr_storage::~cdr_storage(){
  
}

void rk2_tga::cdr_storage::allocate_storage(){
  m_amr->allocate(m_cache, m_phase, m_ncomp);
  m_amr->allocate(m_k1,    m_phase, m_ncomp);
  m_amr->allocate(m_k2,    m_phase, m_ncomp);

  m_amr->allocate(m_scratchIV1,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV2,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV3,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV4,  m_phase, m_ncomp);
}


void rk2_tga::cdr_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
  m_amr->deallocate(m_k1);
  m_amr->deallocate(m_k2);

  m_amr->deallocate(m_scratchIV1);
  m_amr->deallocate(m_scratchIV2);
  m_amr->deallocate(m_scratchIV3);
  m_amr->deallocate(m_scratchIV4);
}

rk2_tga::poisson_storage::poisson_storage(){

}

rk2_tga::poisson_storage::poisson_storage(const RefCountedPtr<amr_mesh>& a_amr,
					  const phase::which_phase       a_phase,
					  const int                      a_ncomp){
  m_amr   = a_amr;
  m_ncomp = a_ncomp;
  m_phase = a_phase;
}

rk2_tga::poisson_storage::~poisson_storage(){
  
}

void rk2_tga::poisson_storage::allocate_storage(){
  m_amr->allocate(m_cache,  m_ncomp);
  m_amr->allocate(m_E_cell, m_phase, SpaceDim);
  m_amr->allocate(m_E_face, m_phase, SpaceDim);
  m_amr->allocate(m_E_eb,   m_phase, SpaceDim);
}

void rk2_tga::poisson_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
  m_amr->deallocate(m_E_cell);
  m_amr->deallocate(m_E_face);
  m_amr->deallocate(m_E_eb);
}

rk2_tga::rte_storage::rte_storage(){

}

rk2_tga::rte_storage::rte_storage(const RefCountedPtr<amr_mesh>& a_amr,
				  const phase::which_phase       a_phase,
				  const int                      a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

rk2_tga::rte_storage::~rte_storage(){
  
}

void rk2_tga::rte_storage::allocate_storage(){
  m_amr->allocate(m_cache,      m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV,  m_phase, m_ncomp);
}

void rk2_tga::rte_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
  m_amr->deallocate(m_scratchIV);
}

rk2_tga::sigma_storage::sigma_storage(){

}

rk2_tga::sigma_storage::sigma_storage(const RefCountedPtr<amr_mesh>& a_amr,
				      const phase::which_phase       a_phase,
				      const int                      a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

rk2_tga::sigma_storage::~sigma_storage(){
}

void rk2_tga::sigma_storage::allocate_storage(){
  m_amr->allocate(m_cache, m_phase, m_ncomp);
  m_amr->allocate(m_k1,    m_phase, m_ncomp);
  m_amr->allocate(m_k2,    m_phase, m_ncomp);
}

void rk2_tga::sigma_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
  m_amr->deallocate(m_k1);
  m_amr->deallocate(m_k2);
#include "CD_NamespaceFooter.H"
