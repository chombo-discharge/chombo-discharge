/*!
  @file   MFQuadCFInterp.cpp
  @brief  Implementation of MFQuadCFInterp.H
  @author Robert Marskar
  @date   Dec. 2017
*/

#include "MFQuadCFInterp.H"

#include "CD_NamespaceHeader.H"

MFQuadCFInterp::MFQuadCFInterp(){
}

MFQuadCFInterp::~MFQuadCFInterp(){
}

MFQuadCFInterp::MFQuadCFInterp(const Vector<RefCountedPtr<EBQuadCFInterp> >& a_quadcfi){
  this->define(a_quadcfi);
}

void MFQuadCFInterp::define(const Vector<RefCountedPtr<EBQuadCFInterp> >& a_quadcfi){
  m_quadcfi = a_quadcfi;
}

const RefCountedPtr<EBQuadCFInterp>& MFQuadCFInterp::getNWOEBQuadCFInterp_ptr(const int a_phase) const {
  return m_quadcfi[a_phase];
}

EBQuadCFInterp& MFQuadCFInterp::getNWOEBQuadCFInterp(const int a_phase) {
  return *m_quadcfi[a_phase];
}

const EBQuadCFInterp& MFQuadCFInterp::getNWOEBQuadCFInterp(const int a_phase) const {
  return *m_quadcfi[a_phase];
}

#include "CD_NamespaceFooter.H"
