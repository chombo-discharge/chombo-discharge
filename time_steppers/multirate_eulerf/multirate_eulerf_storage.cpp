/*!
  @file   multirate_eulerf_storage.cpp
  @brief  Implementation of multirate_eulerf_storage.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "multirate_eulerf.H"
#include "multirate_eulerf_storage.H"

multirate_eulerf::cdr_storage::cdr_storage(){

}

multirate_eulerf::cdr_storage::cdr_storage(const RefCountedPtr<amr_mesh>& a_amr,
					   const phase::which_phase       a_phase,
					   const int                      a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

multirate_eulerf::cdr_storage::~cdr_storage(){
  
}

void multirate_eulerf::cdr_storage::allocate_storage(){
  m_amr->allocate(m_cache,   m_phase, m_ncomp);
  m_amr->allocate(m_scratch, m_phase, m_ncomp);

  m_amr->allocate(m_scratchIV1,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV2,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV3,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV4,  m_phase, m_ncomp);
}


void multirate_eulerf::cdr_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
  m_amr->deallocate(m_scratch);

  m_amr->deallocate(m_scratchIV1);
  m_amr->deallocate(m_scratchIV2);
  m_amr->deallocate(m_scratchIV3);
  m_amr->deallocate(m_scratchIV4);
}

multirate_eulerf::poisson_storage::poisson_storage(){

}

multirate_eulerf::poisson_storage::poisson_storage(const RefCountedPtr<amr_mesh>& a_amr,
						   const phase::which_phase       a_phase,
						   const int                      a_ncomp){
  m_amr   = a_amr;
  m_ncomp = a_ncomp;
  m_phase = a_phase;
}

multirate_eulerf::poisson_storage::~poisson_storage(){
  
}

void multirate_eulerf::poisson_storage::allocate_storage(){
  m_amr->allocate(m_cache,  m_ncomp);
  m_amr->allocate(m_E_cell, m_phase, SpaceDim);
  m_amr->allocate(m_E_face, m_phase, SpaceDim);
  m_amr->allocate(m_E_eb,   m_phase, SpaceDim);
}

void multirate_eulerf::poisson_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
  m_amr->deallocate(m_E_cell);
  m_amr->deallocate(m_E_face);
  m_amr->deallocate(m_E_eb);
}

multirate_eulerf::rte_storage::rte_storage(){

}

multirate_eulerf::rte_storage::rte_storage(const RefCountedPtr<amr_mesh>& a_amr,
					   const phase::which_phase       a_phase,
					   const int                      a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

multirate_eulerf::rte_storage::~rte_storage(){
  
}

void multirate_eulerf::rte_storage::allocate_storage(){
  m_amr->allocate(m_cache,      m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV,  m_phase, m_ncomp);
}

void multirate_eulerf::rte_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
  m_amr->deallocate(m_scratchIV);
}

multirate_eulerf::sigma_storage::sigma_storage(){

}

multirate_eulerf::sigma_storage::sigma_storage(const RefCountedPtr<amr_mesh>& a_amr,
					       const phase::which_phase       a_phase,
					       const int                      a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

multirate_eulerf::sigma_storage::~sigma_storage(){
}

void multirate_eulerf::sigma_storage::allocate_storage(){
  m_amr->allocate(m_cache,   m_phase, m_ncomp);
  m_amr->allocate(m_scratch, m_phase, m_ncomp);
}

void multirate_eulerf::sigma_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
  m_amr->deallocate(m_scratch);
}
