/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFQuadCFInterp.cpp
  @brief  Implementation of CD_MFQuadCFInterp.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFQuadCFInterp.H>
#include <CD_NamespaceHeader.H>

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

const RefCountedPtr<EBQuadCFInterp>& MFQuadCFInterp::getEBQuadCFInterpPointer(const int a_phase) const {
  return m_quadcfi[a_phase];
}

EBQuadCFInterp& MFQuadCFInterp::getEBQuadCFInterp(const int a_phase) {
  return *m_quadcfi[a_phase];
}

const EBQuadCFInterp& MFQuadCFInterp::getEBQuadCFInterp(const int a_phase) const {
  return *m_quadcfi[a_phase];
}

#include <CD_NamespaceFooter.H>
