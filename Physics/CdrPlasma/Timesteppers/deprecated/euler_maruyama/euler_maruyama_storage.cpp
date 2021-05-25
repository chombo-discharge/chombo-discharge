/*!
  @file   euler_maruyama_storage.cpp
  @brief  Implementation of euler_maruyama_storage.H
  @author Robert Marskar
  @date   Aug. 2019
*/

#include "euler_maruyama.H"
#include "euler_maruyama_storage.H"

euler_maruyama::CdrStorage::CdrStorage(){

}

euler_maruyama::CdrStorage::CdrStorage(const RefCountedPtr<AmrMesh>& a_amr,
					 const phase::which_phase       a_phase,
					 const int                      a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

euler_maruyama::CdrStorage::~CdrStorage(){
  deallocate_storage();
}

void euler_maruyama::CdrStorage::allocate_storage(){
  m_amr->allocate(m_scratch,  m_phase, m_ncomp);
  m_amr->allocate(m_scratch2, m_phase, m_ncomp);
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

void euler_maruyama::CdrStorage::deallocate_storage(){
  m_amr->deallocate(m_scratch);
  m_amr->deallocate(m_scratch2);
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

euler_maruyama::FieldStorage::FieldStorage(){

}

euler_maruyama::FieldStorage::FieldStorage(const RefCountedPtr<AmrMesh>& a_amr,
						 const phase::which_phase a_phase,
						 const int a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

euler_maruyama::FieldStorage::~FieldStorage(){
  deallocate_storage();
}

void euler_maruyama::FieldStorage::allocate_storage(){
  m_amr->allocate(m_E_cell,   m_phase, SpaceDim);
  m_amr->allocate(m_E_eb,     m_phase, SpaceDim);
  m_amr->allocate(m_E_dom,    m_phase, SpaceDim);
}

void euler_maruyama::FieldStorage::deallocate_storage(){
  m_amr->deallocate(m_E_cell);
  m_amr->deallocate(m_E_eb);
  m_amr->deallocate(m_E_dom);
}

euler_maruyama::RtStorage::RtStorage(){

}

euler_maruyama::RtStorage::RtStorage(const RefCountedPtr<AmrMesh>& a_amr,
					 const phase::which_phase a_phase,
					 const int a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

euler_maruyama::RtStorage::~RtStorage(){

}

void euler_maruyama::RtStorage::allocate_storage(){
  m_amr->allocate(m_scratchIV,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIF,  m_phase, m_ncomp);
}

void euler_maruyama::RtStorage::deallocate_storage(){
  m_amr->deallocate(m_scratchIV);
  m_amr->deallocate(m_scratchIF);
}

euler_maruyama::SigmaStorage::SigmaStorage(){

}

euler_maruyama::SigmaStorage::SigmaStorage(const RefCountedPtr<AmrMesh>& a_amr,
					     const phase::which_phase a_phase,
					     const int a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

euler_maruyama::SigmaStorage::~SigmaStorage(){

}
  
void euler_maruyama::SigmaStorage::allocate_storage(){

}

void euler_maruyama::SigmaStorage::deallocate_storage(){

#include "CD_NamespaceFooter.H"
