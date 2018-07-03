/*!
  @file   splitstep_rk2_tga_trapz_storage.cpp
  @brief  Implementation of splitstep_rk2_tga_trapz_storage.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "splitstep_rk2_tga_trapz.H"
#include "splitstep_rk2_tga_trapz_storage.H"

splitstep_rk2_tga_trapz::cdr_storage::cdr_storage(){

}

splitstep_rk2_tga_trapz::cdr_storage::cdr_storage(const RefCountedPtr<amr_mesh>& a_amr,
					const phase::which_phase       a_phase,
					const int                      a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

splitstep_rk2_tga_trapz::cdr_storage::~cdr_storage(){
  
}

void splitstep_rk2_tga_trapz::cdr_storage::allocate_storage(){
  m_amr->allocate(m_cache, m_phase, m_ncomp);
  m_amr->allocate(m_phi,   m_phase, m_ncomp);
  m_amr->allocate(m_k1,    m_phase, m_ncomp);
  m_amr->allocate(m_k2,    m_phase, m_ncomp);

  m_amr->allocate(m_scratchIV1,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV2,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV3,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV4,  m_phase, m_ncomp);
}


void splitstep_rk2_tga_trapz::cdr_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
  m_amr->deallocate(m_phi);
  m_amr->deallocate(m_k1);
  m_amr->deallocate(m_k2);

  m_amr->deallocate(m_scratchIV1);
  m_amr->deallocate(m_scratchIV2);
  m_amr->deallocate(m_scratchIV3);
  m_amr->deallocate(m_scratchIV4);
}

splitstep_rk2_tga_trapz::poisson_storage::poisson_storage(){

}

splitstep_rk2_tga_trapz::poisson_storage::poisson_storage(const RefCountedPtr<amr_mesh>& a_amr,
						const phase::which_phase       a_phase,
						const int                      a_ncomp){
  m_amr   = a_amr;
  m_ncomp = a_ncomp;
  m_phase = a_phase;
}

splitstep_rk2_tga_trapz::poisson_storage::~poisson_storage(){
  
}

void splitstep_rk2_tga_trapz::poisson_storage::allocate_storage(){
  m_amr->allocate(m_cache,  m_ncomp);
  m_amr->allocate(m_phi,    m_ncomp);
  m_amr->allocate(m_E_cell, m_phase, SpaceDim);
  m_amr->allocate(m_E_face, m_phase, SpaceDim);
  m_amr->allocate(m_E_eb,   m_phase, SpaceDim);
}

void splitstep_rk2_tga_trapz::poisson_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
  m_amr->deallocate(m_phi);
  m_amr->deallocate(m_E_cell);
  m_amr->deallocate(m_E_face);
  m_amr->deallocate(m_E_eb);
}

splitstep_rk2_tga_trapz::rte_storage::rte_storage(){

}

splitstep_rk2_tga_trapz::rte_storage::rte_storage(const RefCountedPtr<amr_mesh>& a_amr,
					const phase::which_phase       a_phase,
					const int                      a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

splitstep_rk2_tga_trapz::rte_storage::~rte_storage(){
  
}

void splitstep_rk2_tga_trapz::rte_storage::allocate_storage(){
  m_amr->allocate(m_cache,      m_phase, m_ncomp);
  m_amr->allocate(m_phi,        m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV,  m_phase, m_ncomp);
}

void splitstep_rk2_tga_trapz::rte_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
  m_amr->deallocate(m_phi);
  m_amr->deallocate(m_scratchIV);
}

splitstep_rk2_tga_trapz::sigma_storage::sigma_storage(){

}

splitstep_rk2_tga_trapz::sigma_storage::sigma_storage(const RefCountedPtr<amr_mesh>& a_amr,
					    const phase::which_phase       a_phase,
					    const int                      a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

splitstep_rk2_tga_trapz::sigma_storage::~sigma_storage(){
}

void splitstep_rk2_tga_trapz::sigma_storage::allocate_storage(){
  m_amr->allocate(m_cache, m_phase, m_ncomp);
  m_amr->allocate(m_phi,   m_phase, m_ncomp);
  m_amr->allocate(m_k1,    m_phase, m_ncomp);
  m_amr->allocate(m_k2,    m_phase, m_ncomp);
}

void splitstep_rk2_tga_trapz::sigma_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
  m_amr->deallocate(m_phi);
  m_amr->deallocate(m_k1);
  m_amr->deallocate(m_k2);
}
