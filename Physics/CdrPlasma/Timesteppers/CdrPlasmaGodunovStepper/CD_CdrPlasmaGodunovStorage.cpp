/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaGodunovStorage.cpp
  @brief  Implementation of CdrPlasmaGodunovStepper_storage.H
  @author Robert Marskar
  @date   Aug. 2019
*/

// Our includes
#include <CD_CdrPlasmaGodunovStepper.H>
#include <CD_CdrPlasmaGodunovStorage.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::CdrPlasma;

CdrPlasmaGodunovStepper::CdrStorage::CdrStorage(){

}

CdrPlasmaGodunovStepper::CdrStorage::CdrStorage(const RefCountedPtr<AmrMesh>& a_amr,
						const std::string              a_realm,
						const phase::which_phase       a_phase,
						const int                      a_ncomp){
  m_amr   = a_amr;
  m_realm = a_realm;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

CdrPlasmaGodunovStepper::CdrStorage::~CdrStorage(){
  deallocateStorage();
}

void CdrPlasmaGodunovStepper::CdrStorage::allocateStorage(){
  m_amr->allocate(m_scratch,     m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_scratch2,    m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_scratch3,    m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_cellExtr,    m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_gradient,    m_realm, m_phase, SpaceDim);
  m_amr->allocate(m_divDgradPhi, m_realm, m_phase, SpaceDim);  

  m_amr->allocate(m_scratchIVs, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_scratchIVD, m_realm, m_phase, SpaceDim);

  m_amr->allocate(m_scratchIV1, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV2, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV3, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_scratchIV4, m_realm, m_phase, m_ncomp);

  m_amr->allocate(m_scratchIF1, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_scratchIF2, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_scratchIF3, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_scratchIF4, m_realm, m_phase, m_ncomp);
}

void CdrPlasmaGodunovStepper::CdrStorage::deallocateStorage(){
  m_amr->deallocate(m_scratch    );
  m_amr->deallocate(m_scratch2   );
  m_amr->deallocate(m_scratch3   );
  m_amr->deallocate(m_cellExtr   );
  m_amr->deallocate(m_gradient   );
  m_amr->deallocate(m_divDgradPhi);  

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

CdrPlasmaGodunovStepper::FieldStorage::FieldStorage(){

}

CdrPlasmaGodunovStepper::FieldStorage::FieldStorage(const RefCountedPtr<AmrMesh>& a_amr,
						    const std::string              a_realm,
						    const phase::which_phase       a_phase,
						    const int                      a_ncomp){
  m_amr   = a_amr;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
  m_realm = a_realm;
}

CdrPlasmaGodunovStepper::FieldStorage::~FieldStorage(){
  deallocateStorage();
}

void CdrPlasmaGodunovStepper::FieldStorage::allocateStorage(){
  m_amr->allocate(m_E_cell, m_realm, m_phase, SpaceDim);
  m_amr->allocate(m_E_eb,   m_realm, m_phase, SpaceDim);
  m_amr->allocate(m_E_dom,  m_realm, m_phase, SpaceDim);
}

void CdrPlasmaGodunovStepper::FieldStorage::deallocateStorage(){
  m_amr->deallocate(m_E_cell);
  m_amr->deallocate(m_E_eb);
  m_amr->deallocate(m_E_dom);
}

CdrPlasmaGodunovStepper::RtStorage::RtStorage(){

}

CdrPlasmaGodunovStepper::RtStorage::RtStorage(const RefCountedPtr<AmrMesh>& a_amr,
					      const std::string              a_realm,
					      const phase::which_phase       a_phase,
					      const int                      a_ncomp){
  m_amr   = a_amr;
  m_realm = a_realm;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

CdrPlasmaGodunovStepper::RtStorage::~RtStorage(){

}

void CdrPlasmaGodunovStepper::RtStorage::allocateStorage(){
  m_amr->allocate(m_scratchIV, m_realm, m_phase, m_ncomp);
  m_amr->allocate(m_scratchIF, m_realm, m_phase, m_ncomp);
}

void CdrPlasmaGodunovStepper::RtStorage::deallocateStorage(){
  m_amr->deallocate(m_scratchIV);
  m_amr->deallocate(m_scratchIF);
}

CdrPlasmaGodunovStepper::SigmaStorage::SigmaStorage(){

}

CdrPlasmaGodunovStepper::SigmaStorage::SigmaStorage(const RefCountedPtr<AmrMesh>& a_amr,
						    const std::string              a_realm,
						    const phase::which_phase       a_phase,
						    const int                      a_ncomp){
  m_amr   = a_amr;
  m_realm = a_realm;
  m_phase = a_phase;
  m_ncomp = a_ncomp;
}

CdrPlasmaGodunovStepper::SigmaStorage::~SigmaStorage(){

}

void CdrPlasmaGodunovStepper::SigmaStorage::allocateStorage(){
  m_amr->allocate(m_scratch, m_realm, m_phase, m_ncomp);
}

void CdrPlasmaGodunovStepper::SigmaStorage::deallocateStorage(){
  m_amr->deallocate(m_scratch);
}

#include <CD_NamespaceFooter.H>
