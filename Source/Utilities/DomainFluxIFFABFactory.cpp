/*!
  @file   DomainFluxIFFABFactory.cpp
  @author Robert Marskar
  @brief  Implementation of DomainFluxIFFABFactory.H
*/

#include "DomainFluxIFFABFactory.H"

#include "CD_NamespaceHeader.H"

DomainFluxIFFABFactory::DomainFluxIFFABFactory(const EBISLayout& a_ebisl, const ProblemDomain& a_domain){
  m_ebisl  = a_ebisl;
  m_domain = a_domain;
}

DomainFluxIFFABFactory::~DomainFluxIFFABFactory(){

}

DomainFluxIFFAB* DomainFluxIFFABFactory::create(const Box& a_box, int a_ncomps, const DataIndex& a_dit) const {
  return new DomainFluxIFFAB(m_domain, m_ebisl[a_dit], a_box, a_ncomps);
}
#include "CD_NamespaceFooter.H"
