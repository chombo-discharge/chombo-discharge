/*!
  @file   euler_maruyama_storage.cpp
  @brief  Implementation of euler_maruyama_storage.H
  @author Robert Marskar
  @date   Aug. 2019
*/

#include "euler_maruyama.H"
#include "euler_maruyama_storage.H"

euler_maruyama::cdr_storage::cdr_storage(){

}

euler_maruyama::cdr_storage::cdr_storage(const RefCountedPtr<amr_mesh>& a_amr,
					 const phase::which_phase       a_phase,
					 const int                      a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

euler_maruyama::cdr_storage::~cdr_storage(){
  deallocate_storage();
}

void euler_maruyama::cdr_storage::allocate_storage(){
  m_amr->allocate(m_scratch,  m_phase, m_ncomp);
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

void euler_maruyama::cdr_storage::deallocate_storage(){
  m_amr->deallocate(m_scratch);
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

euler_maruyama::poisson_storage::poisson_storage(){

}

euler_maruyama::poisson_storage::poisson_storage(const RefCountedPtr<amr_mesh>& a_amr,
						 const phase::which_phase a_phase,
						 const int a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

euler_maruyama::poisson_storage::~poisson_storage(){
  deallocate_storage();
}

void euler_maruyama::poisson_storage::allocate_storage(){
  m_amr->allocate(m_E_cell,   m_phase, SpaceDim);
  m_amr->allocate(m_E_face,   m_phase, SpaceDim);
  m_amr->allocate(m_E_eb,     m_phase, SpaceDim);
  m_amr->allocate(m_E_dom,    m_phase, SpaceDim);
}

void euler_maruyama::poisson_storage::deallocate_storage(){
  m_amr->deallocate(m_E_cell);
  m_amr->deallocate(m_E_face);
  m_amr->deallocate(m_E_eb);
  m_amr->deallocate(m_E_dom);
}

euler_maruyama::rte_storage::rte_storage(){

}

euler_maruyama::rte_storage::rte_storage(const RefCountedPtr<amr_mesh>& a_amr,
					 const phase::which_phase a_phase,
					 const int a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

euler_maruyama::rte_storage::~rte_storage(){

}

void euler_maruyama::rte_storage::allocate_storage(){
  m_amr->allocate(m_scratchIV,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIF,  m_phase, m_ncomp);
}

void euler_maruyama::rte_storage::deallocate_storage(){
  m_amr->deallocate(m_scratchIV);
  m_amr->deallocate(m_scratchIF);
}

euler_maruyama::sigma_storage::sigma_storage(){

}

euler_maruyama::sigma_storage::sigma_storage(const RefCountedPtr<amr_mesh>& a_amr,
					   const phase::which_phase a_phase,
					   const int a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

euler_maruyama::sigma_storage::~sigma_storage(){

}
  
void euler_maruyama::sigma_storage::allocate_storage(){

}

void euler_maruyama::sigma_storage::deallocate_storage(){

}
