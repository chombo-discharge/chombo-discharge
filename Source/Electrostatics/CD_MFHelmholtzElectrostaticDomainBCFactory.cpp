/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
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
  : m_electrostaticBCs(a_electrostaticBCs)
{
  CH_TIME("MFHelmholtzElectrostaticDomainBCFactory::MFHelmholtzElectrostaticDomainBCFactory()");
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
