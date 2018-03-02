/*!
  @file   splitstep_euler_f_storage.cpp
  @brief  Implementation of splitstep_euler_f_storage.H
  @author Robert Marskar
  @date   Feb. 2018
*/

#include "splitstep_euler_f.H"
#include "splitstep_euler_f_storage.H"

splitstep_euler_f::cdr_storage::cdr_storage(){

}

splitstep_euler_f::cdr_storage::cdr_storage(const RefCountedPtr<amr_mesh>& a_amr,
					    const phase::which_phase       a_phase,
					    const int                      a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

splitstep_euler_f::cdr_storage::~cdr_storage(){
  
}

void splitstep_euler_f::cdr_storage::allocate_storage(){
  m_amr->allocate(m_phi, m_phase, m_ncomp);

  m_amr->allocate(m_scratchIV1,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV2,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV3,  m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV4,  m_phase, m_ncomp);

  m_amr->allocate(m_scratch, m_phase, m_ncomp); 
}

splitstep_euler_f::poisson_storage::poisson_storage(){

}

splitstep_euler_f::poisson_storage::poisson_storage(const RefCountedPtr<amr_mesh>& a_amr,
						    const phase::which_phase       a_phase,
						    const int                      a_ncomp){
  m_amr   = a_amr;
  m_ncomp = a_ncomp;
  m_phase = a_phase;
}

splitstep_euler_f::poisson_storage::~poisson_storage(){

}

void splitstep_euler_f::poisson_storage::allocate_storage(){
  m_amr->allocate(m_phi,         m_ncomp);
  m_amr->allocate(m_scratch_phi, m_ncomp);
  
  m_amr->allocate(m_E_cell,    m_phase, SpaceDim);
  m_amr->allocate(m_E_face,    m_phase, SpaceDim);
  m_amr->allocate(m_E_eb,      m_phase, SpaceDim);
  m_amr->allocate(m_scratch_E, m_phase, SpaceDim);
}

splitstep_euler_f::rte_storage::rte_storage(){

}

splitstep_euler_f::rte_storage::rte_storage(const RefCountedPtr<amr_mesh>& a_amr,
					    const phase::which_phase       a_phase,
					    const int                      a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

splitstep_euler_f::rte_storage::~rte_storage(){
  
}

void splitstep_euler_f::rte_storage::allocate_storage(){
  m_amr->allocate(m_phi,        m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV,  m_phase, m_ncomp);
}

splitstep_euler_f::sigma_storage::sigma_storage(){

}

splitstep_euler_f::sigma_storage::sigma_storage(const RefCountedPtr<amr_mesh>& a_amr,
						const phase::which_phase       a_phase,
						const int                      a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

splitstep_euler_f::sigma_storage::~sigma_storage(){
}

void splitstep_euler_f::sigma_storage::allocate_storage(){
  m_amr->allocate(m_phi,     m_phase, m_ncomp);
  m_amr->allocate(m_scratch, m_phase, m_ncomp);
}
