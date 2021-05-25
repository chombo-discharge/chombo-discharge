/*!
  @file   sdc_storage.cpp
  @brief  Implementation of sdc_storage.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "sdc.H"
#include "sdc_storage.H"

sdc::CdrStorage::CdrStorage(){

}

sdc::CdrStorage::CdrStorage(const RefCountedPtr<AmrMesh>& a_amr,
			      const phase::which_phase       a_phase,
			      const int                      a_ncomp){
  m_amr       = a_amr;
  m_phase     = a_phase;
  m_ncomp     = a_ncomp;
}

sdc::CdrStorage::~CdrStorage(){
  deallocate_storage();
}

void sdc::CdrStorage::allocate_storage(const int a_p){
  m_p = a_p;

  m_amr->allocate(m_scratch,  m_phase, m_ncomp);
  m_amr->allocate(m_scratch2, m_phase, m_ncomp);
  m_amr->allocate(m_error,    m_phase, m_ncomp);
  m_amr->allocate(m_old,      m_phase, m_ncomp);
  m_amr->allocate(m_divF,     m_phase, m_ncomp);
  m_amr->allocate(m_cell_mob, m_phase, m_ncomp);
  m_amr->allocate(m_face_mob, m_phase, m_ncomp);
  m_amr->allocate(m_divD,     m_phase, m_ncomp);
  m_amr->allocate(m_divF,     m_phase, m_ncomp);
  m_amr->allocate(m_gradient, m_phase, SpaceDim);
  m_amr->allocate(m_scratchD, m_phase, SpaceDim);

  m_amr->allocate(m_eb_mob,      m_phase, m_ncomp);
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

  m_phi.resize(1+m_p);
  m_FA.resize(1+m_p);
  m_FD.resize(1+m_p);
  m_FR.resize(1+m_p);
  m_F.resize(1+m_p);

  for (int m = 0; m <= m_p; m++){
    m_amr->allocate(m_phi[m], m_phase, m_ncomp);
    m_amr->allocate(m_FA[m],  m_phase, m_ncomp);
    m_amr->allocate(m_FD[m],  m_phase, m_ncomp);
    m_amr->allocate(m_FR[m],  m_phase, m_ncomp);
    m_amr->allocate(m_F[m],   m_phase, m_ncomp);
  }
}

void sdc::CdrStorage::deallocate_storage(){

  m_amr->deallocate(m_scratch);
  m_amr->deallocate(m_scratch2);
  m_amr->deallocate(m_error);
  m_amr->deallocate(m_old);
  m_amr->deallocate(m_divF);
  m_amr->deallocate(m_divD);
  m_amr->deallocate(m_cell_mob);
  m_amr->deallocate(m_face_mob);
  m_amr->deallocate(m_gradient);
  m_amr->deallocate(m_scratchD);

  m_amr->deallocate(m_eb_mob);
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

  for (int m = 0; m <= m_p; m++){
    m_amr->deallocate(m_phi[m]);
    m_amr->deallocate(m_FA[m]);
    m_amr->deallocate(m_FD[m]);
    m_amr->deallocate(m_FR[m]);
    m_amr->deallocate(m_F[m]);
  }
}

sdc::FieldStorage::FieldStorage(){

}

sdc::FieldStorage::FieldStorage(const RefCountedPtr<AmrMesh>& a_amr,
				      const phase::which_phase       a_phase,
				      const int                      a_ncomp){
  m_amr    = a_amr;
  m_ncomp  = a_ncomp;
  m_phase  = a_phase;
}

sdc::FieldStorage::~FieldStorage(){
  deallocate_storage();
}

void sdc::FieldStorage::allocate_storage(const int a_p){
  m_p = a_p;
  
  m_amr->allocate(m_previous, m_ncomp);
  m_amr->allocate(m_E_cell,   m_phase, SpaceDim);
  m_amr->allocate(m_E_face,   m_phase, SpaceDim);
  m_amr->allocate(m_E_eb,     m_phase, SpaceDim);
  m_amr->allocate(m_E_dom,    m_phase, SpaceDim);
}

void sdc::FieldStorage::deallocate_storage(){
  m_amr->deallocate(m_previous);
  m_amr->deallocate(m_E_cell);
  m_amr->deallocate(m_E_face);
  m_amr->deallocate(m_E_eb);
  m_amr->deallocate(m_E_dom);
}

sdc::RtStorage::RtStorage(){

}

sdc::RtStorage::RtStorage(const RefCountedPtr<AmrMesh>& a_amr,
			      const phase::which_phase       a_phase,
			      const int                      a_ncomp){
  m_amr    = a_amr;
  m_phase  = a_phase;
  m_ncomp  = a_ncomp;


}

sdc::RtStorage::~RtStorage(){
  deallocate_storage();
}

void sdc::RtStorage::allocate_storage(const int a_p){
  m_p = a_p;
  
  m_amr->allocate(m_previous,   m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIF,  m_phase, m_ncomp);

  m_phi.resize(1+m_p);

  for (int m = 0; m <= m_p; m++){
    m_amr->allocate(m_phi[m], m_phase, m_ncomp);
  }
}

void sdc::RtStorage::deallocate_storage(){
  m_amr->deallocate(m_previous);
  m_amr->deallocate(m_scratchIV);
  m_amr->deallocate(m_scratchIF);

  for (int m = 0; m <= m_p; m++){
    m_amr->deallocate(m_phi[m]);
  }
}

sdc::SigmaStorage::SigmaStorage(){

}

sdc::SigmaStorage::SigmaStorage(const RefCountedPtr<AmrMesh>& a_amr,
				  const phase::which_phase       a_phase,
				  const int                      a_ncomp){
  m_amr    = a_amr;
  m_phase  = a_phase;
  m_ncomp  = a_ncomp;
}

sdc::SigmaStorage::~SigmaStorage(){
  deallocate_storage();
}

void sdc::SigmaStorage::allocate_storage(const int a_p){
  m_p = a_p;
  
  m_amr->allocate(m_scratch,  m_phase, m_ncomp);
  m_amr->allocate(m_error,    m_phase, m_ncomp);

  m_sigma.resize(1+m_p);
  m_Fnew.resize(1+m_p);
  m_Fold.resize(1+m_p);

  for (int m = 0; m <= m_p; m++){
    m_amr->allocate(m_sigma[m], m_phase, m_ncomp);
    m_amr->allocate(m_Fnew[m],  m_phase, m_ncomp);
    m_amr->allocate(m_Fold[m],  m_phase, m_ncomp);
  }
}

void sdc::SigmaStorage::deallocate_storage(){
  m_amr->deallocate(m_scratch);
  m_amr->deallocate(m_error);

  for (int m = 0; m <= m_p; m++){
    m_amr->deallocate(m_sigma[m]);
    m_amr->deallocate(m_Fnew[m]);
    m_amr->deallocate(m_Fold[m]);
  }
#include "CD_NamespaceFooter.H"
