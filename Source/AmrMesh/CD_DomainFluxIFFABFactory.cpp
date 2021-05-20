/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_DomainFluxIFFABFactory.cpp
  @author Robert Marskar
  @brief  Implementation of CD_DomainFluxIFFABFactory.H
*/

// Our includes
#include <CD_DomainFluxIFFABFactory.H>
#include <CD_NamespaceHeader.H>

DomainFluxIFFABFactory::DomainFluxIFFABFactory(const EBISLayout& a_ebisl, const ProblemDomain& a_domain){
  m_ebisl  = a_ebisl;
  m_domain = a_domain;
}

DomainFluxIFFABFactory::~DomainFluxIFFABFactory(){

}

DomainFluxIFFAB* DomainFluxIFFABFactory::create(const Box& a_box, int a_ncomps, const DataIndex& a_dit) const {
  return new DomainFluxIFFAB(m_domain, m_ebisl[a_dit], a_box, a_ncomps);
}

#include <CD_NamespaceFooter.H>
