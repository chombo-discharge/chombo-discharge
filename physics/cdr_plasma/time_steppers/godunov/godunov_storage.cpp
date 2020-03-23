/*!
  @file   godunov_storage.cpp
  @brief  Implementation of godunov_storage.H
  @author Robert Marskar
  @date   Aug. 2019
*/

#include "godunov.H"
#include "godunov_storage.H"

godunov::cdr_storage::cdr_storage(){

}

godunov::cdr_storage::cdr_storage(const RefCountedPtr<amr_mesh>& a_amr,
				  const phase::which_phase       a_phase,
				  const int                      a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

godunov::cdr_storage::~cdr_storage(){
  deallocate_storage();
}

void godunov::cdr_storage::allocate_storage(){
  m_amr->allocate(m_scratch,  m_phase, m_ncomp);
  m_amr->allocate(m_scratch2, m_phase, m_ncomp);
  m_amr->allocate(m_scratch3, m_phase, m_ncomp);
  m_amr->allocate(m_gradient, m_phase, SpaceDim);

  m_amr->allocate(m_scratchIVs,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIVD,  m_phase, SpaceDim);

  m_amr->allocate(m_scratchIV1,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV2,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV3,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV4,  m_phase, m_ncomp);

  m_amr->allocate(m_scratchIF1,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIF2,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIF3,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIF4,  m_phase, m_ncomp);
}

void godunov::cdr_storage::deallocate_storage(){
  m_amr->deallocate(m_scratch);
  m_amr->deallocate(m_scratch2);
  m_amr->deallocate(m_scratch3);
  m_amr->deallocate(m_gradient);

  m_amr->deallocate(m_scratchIVs);
  m_amr->deallocate(m_scratchIVD);

  m_amr->deallocate(m_scratchIV1);
  m_amr->deallocate(m_scratchIV2);
  m_amr->deallocate(m_scratchIV3);
  m_amr->deallocate(m_scratchIV4);

  m_amr->deallocate(m_scratchIF1);
  m_amr->deallocate(m_scratchIF2);
  m_amr->deallocate(m_scratchIF3);
  m_amr->deallocate(m_scratchIF4);
}

godunov::poisson_storage::poisson_storage(){

}

godunov::poisson_storage::poisson_storage(const RefCountedPtr<amr_mesh>& a_amr,
					  const phase::which_phase a_phase,
					  const int a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

godunov::poisson_storage::~poisson_storage(){
  deallocate_storage();
}

void godunov::poisson_storage::allocate_storage(){
  m_amr->allocate(m_E_cell,   m_phase, SpaceDim);
  m_amr->allocate(m_E_eb,     m_phase, SpaceDim);
  m_amr->allocate(m_E_dom,    m_phase, SpaceDim);
}

void godunov::poisson_storage::deallocate_storage(){
  m_amr->deallocate(m_E_cell);
  m_amr->deallocate(m_E_eb);
  m_amr->deallocate(m_E_dom);
}

godunov::rte_storage::rte_storage(){

}

godunov::rte_storage::rte_storage(const RefCountedPtr<amr_mesh>& a_amr,
				  const phase::which_phase a_phase,
				  const int a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

godunov::rte_storage::~rte_storage(){

}

void godunov::rte_storage::allocate_storage(){
  m_amr->allocate(m_scratchIV,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIF,  m_phase, m_ncomp);
}

void godunov::rte_storage::deallocate_storage(){
  m_amr->deallocate(m_scratchIV);
  m_amr->deallocate(m_scratchIF);
}

godunov::sigma_storage::sigma_storage(){

}

godunov::sigma_storage::sigma_storage(const RefCountedPtr<amr_mesh>& a_amr,
				      const phase::which_phase a_phase,
				      const int a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

godunov::sigma_storage::~sigma_storage(){

}

void godunov::sigma_storage::allocate_storage(){
  m_amr->allocate(m_scratch,  m_phase, m_ncomp);
}

void godunov::sigma_storage::deallocate_storage(){
  m_amr->deallocate(m_scratch);
}
