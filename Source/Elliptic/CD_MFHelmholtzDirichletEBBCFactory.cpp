/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzDirichletEBBCFactory.cpp
  @brief  Implementation of CD_MFHelmholtzDirichletEBBCFactory.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_MFHelmholtzDirichletEBBCFactory.H>
#include <CD_MFHelmholtzDirichletEBBC.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzDirichletEBBCFactory::MFHelmholtzDirichletEBBCFactory(const int  a_order,
                                                                 const int  a_weight,
                                                                 const Real a_value)
{
  CH_TIME("MFHelmholtzDirichletEBBCFactory::MFHelmholtzDirichletEBBCFactory(int, int, Real)");

  CH_assert(a_order > 0);
  CH_assert(a_weight >= 0);

  m_order  = a_order;
  m_weight = a_weight;
  m_domainDropOrder = 0;  

  this->setValue(a_value);
}

MFHelmholtzDirichletEBBCFactory::MFHelmholtzDirichletEBBCFactory(
  const int                                         a_order,
  const int                                         a_weight,
  const std::function<Real(const RealVect& a_pos)>& a_value)
{
  CH_TIME("MFHelmholtzDirichletEBBCFactory::MFHelmholtzDirichletEBBCFactory(int, int, std::function<Real(RealVect))");

  CH_assert(a_order > 0);
  CH_assert(a_weight >= 0);

  m_order  = a_order;
  m_weight = a_weight;

  this->setValue(a_value);
}

MFHelmholtzDirichletEBBCFactory::~MFHelmholtzDirichletEBBCFactory()
{
  CH_TIME("MFHelmholtzDirichletEBBCFactory::~MFHelmholtzDirichletEBBCFactory()");
}

void
MFHelmholtzDirichletEBBCFactory::setValue(const Real a_value)
{
  CH_TIME("MFHelmholtzDirichletEBBCFactory::setValue(Real)");

  m_useConstant = true;
  m_useFunction = false;

  m_constantValue = a_value;
}

void
MFHelmholtzDirichletEBBCFactory::setValue(const std::function<Real(const RealVect& a_pos)>& a_value)
{
  CH_TIME("MFHelmholtzDirichletEBBCFactory::setValue(std::function<Real(RealVect)>)");

  m_useConstant = false;
  m_useFunction = true;

  m_functionValue = a_value;
}

void
MFHelmholtzDirichletEBBCFactory::setDomainDropOrder(const int a_domainSize) {
  CH_TIME("MFHelmholtzDirichletEBBCFactory::setDomainDropOrder()");

  m_domainDropOrder = a_domainSize;
}

RefCountedPtr<EBHelmholtzEBBC>
MFHelmholtzDirichletEBBCFactory::create(const int a_iphase, const RefCountedPtr<MFHelmholtzJumpBC>& a_jumpBC) const
{
  CH_TIME("MFHelmholtzDirichletEBBCFactory::create(int, RefCountedPtr<MFHelmholtzJumpBC>)");

  CH_assert(m_order > 0);
  CH_assert(m_weight >= 0);
  CH_assert(m_useFunction || m_useConstant);

  auto bc = new MFHelmholtzDirichletEBBC(a_iphase, a_jumpBC);

  bc->setOrder(m_order);
  bc->setWeight(m_weight);
  bc->setDomainDropOrder(m_domainDropOrder);

  if (m_useConstant) {
    bc->setValue(m_constantValue);
  }
  else if (m_useFunction) {
    bc->setValue(m_functionValue);
  }
  else {
    MayDay::Error("MFHelmholtzDirichletEBBCFactory::create(...) - logic bust. Not using constant or function");
  }

  return RefCountedPtr<EBHelmholtzEBBC>(bc);
}

#include <CD_NamespaceFooter.H>
