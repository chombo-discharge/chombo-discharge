/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzEddingtonSP1DomainBCFactory.cpp
  @brief  Implementation of CD_EBHelmholtzEddingtonSP1DomainBCFactory.H
  @author Robert Marskar
*/

// Our includes
#include <CD_EBHelmholtzEddingtonSP1DomainBCFactory.H>
#include <CD_EBHelmholtzEddingtonSP1DomainBC.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzEddingtonSP1DomainBCFactory::EBHelmholtzEddingtonSP1DomainBCFactory(
  const EddingtonSP1DomainBc&     a_eddingtonBCs,
  const RefCountedPtr<RtSpecies>& a_species,
  const Real                      a_r1,
  const Real                      a_r2)
{
  m_eddingtonBCs = a_eddingtonBCs;
  m_species      = a_species;
  m_r1           = a_r1;
  m_r2           = a_r2;
}

EBHelmholtzEddingtonSP1DomainBCFactory::~EBHelmholtzEddingtonSP1DomainBCFactory()
{}

RefCountedPtr<EBHelmholtzDomainBC>
EBHelmholtzEddingtonSP1DomainBCFactory::create() const
{
  return RefCountedPtr<EBHelmholtzDomainBC>(new EBHelmholtzEddingtonSP1DomainBC(m_eddingtonBCs, m_species, m_r1, m_r2));
}

#include <CD_NamespaceFooter.H>
