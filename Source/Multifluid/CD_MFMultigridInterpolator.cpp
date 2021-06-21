/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFMultigridInterpolator.cpp
  @brief  Implmementation of CD_MFMultigridInterpolator.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFMultigridInterpolator.H>
#include <CD_NamespaceHeader.H>

MFMultigridInterpolator::MFMultigridInterpolator(){

}

MFMultigridInterpolator::MFMultigridInterpolator(const Vector<RefCountedPtr<EBMultigridInterpolator> >& a_interpolators){
  this->define(a_interpolators);
}

MFMultigridInterpolator::~MFMultigridInterpolator(){

}

void MFMultigridInterpolator::define(const Vector<RefCountedPtr<EBMultigridInterpolator> >& a_interpolators){
  m_interpolators = a_interpolators;
}

const RefCountedPtr<EBMultigridInterpolator>& MFMultigridInterpolator::getInterpolator(const int a_phase) const {
  return m_interpolators[a_phase];
}

#include <CD_NamespaceFooter.H>
