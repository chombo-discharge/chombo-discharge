/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_MFHelmholtzDirichletDomainBCFactory.cpp
  @brief  Implementation of CD_MFHelmholtzDirichletDomainBCFactory.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBHelmholtzDirichletDomainBC.H>
#include <CD_MFHelmholtzDirichletDomainBCFactory.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzDirichletDomainBCFactory::MFHelmholtzDirichletDomainBCFactory()
{
  CH_TIME("MFHelmholtzDirichletDomainBCFactory::MFHelmholtzDirichletDomainBCFactory()");

  m_useConstant = false;
  m_useFunction = false;
}

MFHelmholtzDirichletDomainBCFactory::MFHelmholtzDirichletDomainBCFactory(const Real a_value)
{
  CH_TIME("MFHelmholtzDirichletDomainBCFactory::MFHelmholtzDirichletDomainBCFactory(Real)");

  this->setValue(a_value);
}

MFHelmholtzDirichletDomainBCFactory::MFHelmholtzDirichletDomainBCFactory(
  const std::function<Real(const RealVect& a_pos)>& a_value)
{
  CH_TIME("MFHelmholtzDirichletDomainBCFactory::MFHelmholtzDirichletDomainBCFactory(std::function<Real(RealVect)>)");

  this->setValue(a_value);
}

MFHelmholtzDirichletDomainBCFactory::~MFHelmholtzDirichletDomainBCFactory()
{
  CH_TIME("MFHelmholtzDirichletDomainBCFactory::~MFHelmholtzDirichletDomainBCFactory()");
}

void
MFHelmholtzDirichletDomainBCFactory::setValue(const Real a_value)
{
  CH_TIME("MFHelmholtzDirichletDomainBCFactory::setValue(Real)");

  m_useConstant = true;
  m_useFunction = false;

  m_constantValue = a_value;
}

void
MFHelmholtzDirichletDomainBCFactory::setValue(const std::function<Real(const RealVect& a_pos)>& a_value)
{
  CH_TIME("MFHelmholtzDirichletDomainBCFactory::setValue(std::function<Real(RealVect)>)");

  m_useConstant = false;
  m_useFunction = true;

  m_functionValue = a_value;
}

RefCountedPtr<EBHelmholtzDomainBC>
MFHelmholtzDirichletDomainBCFactory::create(const int a_iphase) const
{
  CH_TIME("MFHelmholtzDirichletDomainBCFactory::create(int)");

  CH_assert(m_useConstant || m_useFunction);

  // Also issue an error because this error will break everything (although I don't see how you could end up there).
  if (!(m_useConstant || m_useFunction)) {
    MayDay::Error("MFHelmholtzDirichletDomainBCFactory::create - logic bust, not using function or constant!");
  }

  auto bc = new EBHelmholtzDirichletDomainBC();

  if (m_useConstant) {
    bc->setValue(m_constantValue);
  }
  else if (m_useFunction) {
    bc->setValue(m_functionValue);
  }

  return RefCountedPtr<EBHelmholtzDomainBC>(bc);
}

#include <CD_NamespaceFooter.H>
