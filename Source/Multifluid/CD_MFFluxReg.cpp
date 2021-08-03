/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFFluxReg.cpp
  @brief  Implementation of CD_MFFluxReg.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFFluxReg.H>
#include <CD_NamespaceHeader.H>

MFFluxReg::MFFluxReg(){
}

MFFluxReg::MFFluxReg(const Vector<RefCountedPtr<EBFluxRegister> >& a_fastfr){
  this->define(a_fastfr);
}

MFFluxReg::~MFFluxReg(){

}

void MFFluxReg::define(const Vector<RefCountedPtr<EBFluxRegister> >& a_fastfr){
  m_fastfr = a_fastfr;
}

const RefCountedPtr<EBFluxRegister>& MFFluxReg::getFluxRegPointer(const int a_phase) const {
  return m_fastfr[a_phase];
}

EBFluxRegister& MFFluxReg::getFluxReg(const int a_phase) {
  return *m_fastfr[a_phase];
}

const EBFluxRegister& MFFluxReg::getFluxReg(const int a_phase) const {
  return *m_fastfr[a_phase];
}

#include <CD_NamespaceFooter.H>
