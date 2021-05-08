/*!
  @file   strang2_storage.cpp
  @brief  Implementation of strang2_storage.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "strang2.H"
#include "strang2_storage.H"

strang2::cdr_storage::cdr_storage(){

}

strang2::cdr_storage::cdr_storage(const int a_stages,
				  const RefCountedPtr<amr_mesh>& a_amr,
				  const phase::which_phase       a_phase,
				  const int                      a_ncomp){
  m_stages    = a_stages;
  m_amr       = a_amr;
  m_phase     = a_phase;
  m_ncomp     = a_ncomp;
  m_has_extra = false;
}

strang2::cdr_storage::~cdr_storage(){
  deallocate_storage();
}

void strang2::cdr_storage::allocate_storage(){
  m_amr->allocate(m_cache,    m_phase, m_ncomp);
  m_amr->allocate(m_scratch,  m_phase, m_ncomp);
  m_amr->allocate(m_backup,   m_phase, m_ncomp);
  m_amr->allocate(m_previous, m_phase, m_ncomp);
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
}


void strang2::cdr_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
  m_amr->deallocate(m_scratch);
  m_amr->deallocate(m_backup);
  m_amr->deallocate(m_previous);
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
}

void strang2::cdr_storage::allocate_extra_storage(const int a_num_extra){
  if(m_has_extra == true){
    MayDay::Abort("strang2::cdr_storage::allocate_extra_storage - already allocated. Did you remember to deallocate first?");
  }

  if(a_num_extra > 0){
    m_extra_storage.resize(a_num_extra);

    const int ncomp = 1;

    for (int i = 0; i < m_extra_storage.size(); i++){
      m_extra_storage[i] = new EBAMRCellData();
      m_amr->allocate(*m_extra_storage[i], m_phase, ncomp);
    }
  }
  
  m_has_extra = true;
}

void strang2::cdr_storage::deallocate_extra_storage(){

  for (int i = 0; i < m_extra_storage.size(); i++){
    m_amr->deallocate(*m_extra_storage[i]);
  }
  m_extra_storage.resize(0);
  m_has_extra = false;
}

strang2::poisson_storage::poisson_storage(){

}

strang2::poisson_storage::poisson_storage(const int a_stages,
					  const RefCountedPtr<amr_mesh>& a_amr,
					  const phase::which_phase       a_phase,
					  const int                      a_ncomp){
  m_stages = a_stages;
  m_amr    = a_amr;
  m_ncomp  = a_ncomp;
  m_phase  = a_phase;
}

strang2::poisson_storage::~poisson_storage(){
  deallocate_storage();
}

void strang2::poisson_storage::allocate_storage(){
  m_amr->allocate(m_cache,   m_ncomp);
  m_amr->allocate(m_backup,  m_ncomp);
  m_amr->allocate(m_E_cell,  m_phase, SpaceDim);
  m_amr->allocate(m_E_face,  m_phase, SpaceDim);
  m_amr->allocate(m_E_eb,    m_phase, SpaceDim);
  m_amr->allocate(m_E_dom,   m_phase, SpaceDim);
}

void strang2::poisson_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
  m_amr->deallocate(m_backup);
  m_amr->deallocate(m_E_cell);
  m_amr->deallocate(m_E_face);
  m_amr->deallocate(m_E_eb);
  m_amr->deallocate(m_E_dom);
}

strang2::rte_storage::rte_storage(){

}

strang2::rte_storage::rte_storage(const int a_stages,
				  const RefCountedPtr<amr_mesh>& a_amr,
				  const phase::which_phase       a_phase,
				  const int                      a_ncomp){
  m_stages = a_stages;
  m_amr    = a_amr;
  m_phase  = a_phase;
  m_ncomp  = a_ncomp;
}

strang2::rte_storage::~rte_storage(){
  deallocate_storage();
}

void strang2::rte_storage::allocate_storage(){
  m_amr->allocate(m_cache,      m_phase, m_ncomp);
  m_amr->allocate(m_backup,     m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIF,  m_phase, m_ncomp);
}

void strang2::rte_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
  m_amr->deallocate(m_backup);
  m_amr->deallocate(m_scratchIV);
  m_amr->deallocate(m_scratchIF);
}

strang2::sigma_storage::sigma_storage(){

}

strang2::sigma_storage::sigma_storage(const int a_stages,
				      const RefCountedPtr<amr_mesh>& a_amr,
				      const phase::which_phase       a_phase,
				      const int                      a_ncomp){
  m_stages = a_stages;
  m_amr    = a_amr;
  m_phase  = a_phase;
  m_ncomp  = a_ncomp;
  m_has_extra = false;
}

strang2::sigma_storage::~sigma_storage(){
  deallocate_storage();
}

void strang2::sigma_storage::allocate_storage(){
  m_amr->allocate(m_cache,    m_phase, m_ncomp);
  m_amr->allocate(m_backup,   m_phase, m_ncomp);
  m_amr->allocate(m_scratch,  m_phase, m_ncomp);
  m_amr->allocate(m_previous, m_phase, m_ncomp);
  m_amr->allocate(m_error,    m_phase, m_ncomp);
}

void strang2::sigma_storage::deallocate_storage(){
  m_amr->deallocate(m_cache);
  m_amr->deallocate(m_backup);
  m_amr->deallocate(m_scratch);
  m_amr->deallocate(m_previous);
  m_amr->deallocate(m_error);
}

void strang2::sigma_storage::allocate_extra_storage(const int a_num_extra){
  if(m_has_extra == true){
    MayDay::Abort("strang2::sigma_storage::allocate_extra_storage - already allocated. Did you remember to deallocate first?");
  }

  if(a_num_extra > 0){
    m_extra_storage.resize(a_num_extra);

    const int ncomp = 1;

    for (int i = 0; i < m_extra_storage.size(); i++){
      m_extra_storage[i] = new EBAMRIVData();
      m_amr->allocate(*m_extra_storage[i], m_phase, ncomp);
    }
  }
  
  m_has_extra = true;
}

void strang2::sigma_storage::deallocate_extra_storage(){

  for (int i = 0; i < m_extra_storage.size(); i++){
    m_amr->deallocate(*m_extra_storage[i]);
  }
  m_extra_storage.resize(0);
  m_has_extra = false;
#include "CD_NamespaceFooter.H"
