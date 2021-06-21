/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFCoarAve.cpp
  @brief  Implementation of CD_MFCoarAve.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFCoarAve.H>
#include <CD_NamespaceHeader.H>

MFCoarAve::MFCoarAve(){
}

MFCoarAve::~MFCoarAve(){
}

MFCoarAve::MFCoarAve(const std::map<Phase, RefCountedPtr<EbCoarAve> >& a_aveOps){
  this->define(a_aveOps);
}

void MFCoarAve::define(const std::map<Phase, RefCountedPtr<EbCoarAve> >& a_aveOps){
  m_aveOps = a_aveOps;
}

const RefCountedPtr<EbCoarAve>& MFCoarAve::getAveOp(const Phase a_phase) const {
  return m_aveOps.at(a_phase);
}

#include <CD_NamespaceFooter.H>
