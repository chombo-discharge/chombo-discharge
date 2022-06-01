/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_MFHelmholtzElectrostaticDomainBCFactory.cpp
  @brief  Implementation of CD_MFHelmholtzElectrostaticDomainBCFactory.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFHelmholtzElectrostaticDomainBCFactory.H>
#include <CD_EBHelmholtzElectrostaticDomainBC.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzElectrostaticDomainBCFactory::MFHelmholtzElectrostaticDomainBCFactory(
  const ElectrostaticDomainBc& a_electrostaticBCs)
{
  CH_TIME("MFHelmholtzElectrostaticDomainBCFactory::MFHelmholtzElectrostaticDomainBCFactory()");

  m_electrostaticBCs = a_electrostaticBCs;
}

MFHelmholtzElectrostaticDomainBCFactory::~MFHelmholtzElectrostaticDomainBCFactory()
{
  CH_TIME("MFHelmholtzElectrostaticDomainBCFactory::~MFHelmholtzElectrostaticDomainBCFactory()");
}

RefCountedPtr<EBHelmholtzDomainBC>
MFHelmholtzElectrostaticDomainBCFactory::create(const int a_iphase) const
{
  CH_TIME("MFHelmholtzElectrostaticDomainBCFactory::create(int)");

  return RefCountedPtr<EBHelmholtzDomainBC>(new EBHelmholtzElectrostaticDomainBC(m_electrostaticBCs));
}

#include <CD_NamespaceFooter.H>
