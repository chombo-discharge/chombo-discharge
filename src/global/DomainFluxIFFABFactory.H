/*!
  @file   DomainFluxIFFABFactory.H
  @author Robert Marskar
  @brief  Declaration of factory class for DomainFluxIFFAB
*/

#ifndef _DOMAINFLUXIFFABFACTORY_
#define _DOMAINFLUXIFFABFACTORY_

#include "DomainFluxIFFAB.H"

#include <Vector.H>
#include <ProblemDomain.H>
#include <EBISLayout.H>
#include <DataIndex.H>

#include "CD_NamespaceHeader.H"
  
/*!
  @brief Factory class for DomainFluxIFFAB
*/
class DomainFluxIFFABFactory : public DataFactory<DomainFluxIFFAB> {
public:

  DomainFluxIFFABFactory(const EBISLayout& a_ebisl, const ProblemDomain& a_domain);
  ~DomainFluxIFFABFactory();

  virtual DomainFluxIFFAB* create(const Box& a_box, int a_ncomps, const DataIndex& a_dit) const;

private:
  EBISLayout    m_ebisl;
  ProblemDomain m_domain;
};
  
#include "CD_NamespaceFooter.H"

#endif
