/*!
  @file   sisdc_storage.cpp
  @brief  Implementation of sisdc_storage.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "sisdc.H"
#include "sisdc_storage.H"

sisdc::cdr_storage::cdr_storage(){

}

sisdc::cdr_storage::cdr_storage(const RefCountedPtr<amr_mesh>& a_amr,
				const phase::which_phase       a_phase,
				const int                      a_ncomp){
  m_amr       = a_amr;
  m_phase     = a_phase;
  m_ncomp     = a_ncomp;
}

sisdc::cdr_storage::~cdr_storage(){
  deallocate_storage();
}

void sisdc::cdr_storage::allocate_storage(const int a_order){
  m_order = Max(2, a_order); // Order = 1 is a special case
  
  m_amr->allocate(m_previous, m_phase, m_ncomp);
  m_amr->allocate(m_scratch,  m_phase, m_ncomp);
  m_amr->allocate(m_error,    m_phase, m_ncomp);
  m_amr->allocate(m_gradient, m_phase, SpaceDim);

  m_amr->allocate(m_scratchIV1,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV2,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV3,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV4,  m_phase, m_ncomp);

  m_amr->allocate(m_scratchIF1,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIF2,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIF3,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIF4,  m_phase, m_ncomp);

  m_phi.resize(m_order);
  m_FAR.resize(m_order);
  m_FD.resize(m_order);
  m_F.resize(m_order);

  for (int m = 0; m < m_order; m++){
    m_amr->allocate(m_phi[m], m_phase, m_ncomp);
    m_amr->allocate(m_FAR[m], m_phase, m_ncomp);
    m_amr->allocate(m_FD[m],  m_phase, m_ncomp);
    m_amr->allocate(m_F[m],   m_phase, m_ncomp);
  }
}


void sisdc::cdr_storage::deallocate_storage(){
  m_amr->deallocate(m_previous);
  m_amr->deallocate(m_scratch);
  m_amr->deallocate(m_error);
  m_amr->deallocate(m_gradient);

  m_amr->deallocate(m_scratchIV1);
  m_amr->deallocate(m_scratchIV2);
  m_amr->deallocate(m_scratchIV3);
  m_amr->deallocate(m_scratchIV4);

  m_amr->deallocate(m_scratchIF1);
  m_amr->deallocate(m_scratchIF2);
  m_amr->deallocate(m_scratchIF3);
  m_amr->deallocate(m_scratchIF4);

  for (int m = 0; m < m_order; m++){
    m_amr->deallocate(m_phi[m]);
    m_amr->deallocate(m_FAR[m]);
    m_amr->deallocate(m_FD[m]);
    m_amr->deallocate(m_F[m]);
  }

  m_phi.resize(0);
  m_FAR.resize(0);
  m_FD.resize(0);
  m_F.resize(0);
}

sisdc::poisson_storage::poisson_storage(){

}

sisdc::poisson_storage::poisson_storage(const RefCountedPtr<amr_mesh>& a_amr,
					const phase::which_phase       a_phase,
					const int                      a_ncomp){
  m_amr    = a_amr;
  m_ncomp  = a_ncomp;
  m_phase  = a_phase;
}

sisdc::poisson_storage::~poisson_storage(){
  deallocate_storage();
}

void sisdc::poisson_storage::allocate_storage(const int a_order){
  m_order = Max(2, a_order); // Order = 1 is a special case
  
  m_amr->allocate(m_previous, m_ncomp);
  m_amr->allocate(m_E_cell,   m_phase, SpaceDim);
  m_amr->allocate(m_E_face,   m_phase, SpaceDim);
  m_amr->allocate(m_E_eb,     m_phase, SpaceDim);
  m_amr->allocate(m_E_dom,    m_phase, SpaceDim);
}

void sisdc::poisson_storage::deallocate_storage(){
  m_amr->deallocate(m_previous);
  m_amr->deallocate(m_E_cell);
  m_amr->deallocate(m_E_face);
  m_amr->deallocate(m_E_eb);
  m_amr->deallocate(m_E_dom);
}

sisdc::rte_storage::rte_storage(){

}

sisdc::rte_storage::rte_storage(const RefCountedPtr<amr_mesh>& a_amr,
				const phase::which_phase       a_phase,
				const int                      a_ncomp){
  m_amr    = a_amr;
  m_phase  = a_phase;
  m_ncomp  = a_ncomp;
}

sisdc::rte_storage::~rte_storage(){
  deallocate_storage();
}

void sisdc::rte_storage::allocate_storage(const int a_order){
  m_order = Max(2, a_order); // Order = 1 is a special case
  
  m_amr->allocate(m_previous,   m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIF,  m_phase, m_ncomp);
}

void sisdc::rte_storage::deallocate_storage(){
  m_amr->deallocate(m_previous);
  m_amr->deallocate(m_scratchIV);
  m_amr->deallocate(m_scratchIF);
}

sisdc::sigma_storage::sigma_storage(){

}

sisdc::sigma_storage::sigma_storage(const RefCountedPtr<amr_mesh>& a_amr,
				    const phase::which_phase       a_phase,
				    const int                      a_ncomp){
  m_amr    = a_amr;
  m_phase  = a_phase;
  m_ncomp  = a_ncomp;
}

sisdc::sigma_storage::~sigma_storage(){
  deallocate_storage();
}

void sisdc::sigma_storage::allocate_storage(const int a_order){
  m_order = Max(2, a_order); // Order = 1 is a special case
  
  m_amr->allocate(m_previous, m_phase, m_ncomp);
  m_amr->allocate(m_scratch,  m_phase, m_ncomp);
  m_amr->allocate(m_error,    m_phase, m_ncomp);

  m_sigma.resize(m_order);
  m_Fsig.resize(m_order);

  for (int m = 0; m < m_order; m++){
    m_amr->allocate(m_sigma[m], m_phase, m_ncomp);
    m_amr->allocate(m_Fsig[m],  m_phase, m_ncomp);
  }
}

void sisdc::sigma_storage::deallocate_storage(){
  m_amr->deallocate(m_previous);
  m_amr->deallocate(m_scratch);
  m_amr->deallocate(m_error);

  for (int m = 0; m < m_order; m++){
    m_amr->deallocate(m_sigma[m]);
    m_amr->deallocate(m_Fsig[m]);
  }

  m_sigma.resize(0);
  m_Fsig.resize(0);
}
