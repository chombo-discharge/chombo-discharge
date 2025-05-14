/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDirichletEBBCFactory.cpp
  @brief  Implementation of CD_EBHelmholtzDirichletEBBCFactory.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBHelmholtzDirichletEBBC.H>
#include <CD_EBHelmholtzDirichletEBBCFactory.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzDirichletEBBCFactory::EBHelmholtzDirichletEBBCFactory()
{
  CH_TIME("EBHelmholtzDirichletEBBCFactory::EBHelmholtzDirichletEBBCFactory()");

  m_order           = -1;
  m_weight          = -1;
  m_domainDropOrder = 0;
  m_useConstant     = false;
  m_useFunction     = false;
}

EBHelmholtzDirichletEBBCFactory::EBHelmholtzDirichletEBBCFactory(const int  a_order,
                                                                 const int  a_weight,
                                                                 const Real a_value)
{
  CH_TIME("EBHelmholtzDirichletEBBCFactory::EBHelmholtzDirichletEBBCFactory(int, int, Real)");

  CH_assert(a_order > 0);
  CH_assert(a_weight >= 0);

  this->setOrder(a_order);
  this->setWeight(a_weight);
  this->setValue(a_value);
}

EBHelmholtzDirichletEBBCFactory::EBHelmholtzDirichletEBBCFactory(
  const int                                         a_order,
  const int                                         a_weight,
  const std::function<Real(const RealVect& a_pos)>& a_value)
{
  CH_TIME("EBHelmholtzDirichletEBBCFactory::EBHelmholtzDirichletEBBCFactory(int, int, std::function<Real(RealVect)>)");

  CH_assert(a_order > 0);
  CH_assert(a_weight >= 0);

  this->setOrder(a_order);
  this->setWeight(a_weight);
  this->setValue(a_value);
}

EBHelmholtzDirichletEBBCFactory::~EBHelmholtzDirichletEBBCFactory()
{
  CH_TIME("EBHelmholtzDirichletEBBCFactory::~EBHelmholtzDirichletEBBCFactory()");
}

void
EBHelmholtzDirichletEBBCFactory::setOrder(const int a_order)
{
  CH_TIME("EBHelmholtzDirichletEBBCFactory::setOrder(int)");

  CH_assert(a_order > 0);

  m_order = a_order;
}

void
EBHelmholtzDirichletEBBCFactory::setWeight(const int a_weight)
{
  CH_TIME("EBHelmholtzDirichletEBBCFactory::setWeight(int)");

  CH_assert(a_weight >= 0);

  m_weight = a_weight;
}

void
EBHelmholtzDirichletEBBCFactory::setValue(const Real a_value)
{
  CH_TIME("EBHelmholtzDirichletEBBCFactory::setValue(Real)");

  m_useConstant = true;
  m_useFunction = false;

  m_constantValue = a_value;
}

void
EBHelmholtzDirichletEBBCFactory::setValue(const std::function<Real(const RealVect& a_pos)>& a_value)
{
  CH_TIME("EBHelmholtzDirichletEBBCFactory::setValue(std::function<Real(RealVect)>)");

  m_useConstant = false;
  m_useFunction = true;

  m_functionValue = a_value;
}

void
EBHelmholtzDirichletEBBCFactory::setDomainDropOrder(const int a_domainSize)
{
  CH_TIME("EBHelmholtzDirichletEBBCFactory::setDomainDropOrder()");

  m_domainDropOrder = a_domainSize;
}

void
EBHelmholtzDirichletEBBCFactory::setCoarseGridDropOrder(const bool a_dropOrder)
{
  m_dropOrder = a_dropOrder;
}

RefCountedPtr<EBHelmholtzEBBC>
EBHelmholtzDirichletEBBCFactory::create()
{
  CH_TIME("EBHelmholtzDirichletEBBCFactory::create()");

  CH_assert(m_useConstant || m_useFunction);
  CH_assert(m_order > 0);
  CH_assert(m_weight >= 0);

  // Also issue run-time error
  if (!(m_order > 0 && m_weight >= 0)) {
    MayDay::Error("EBHelmholtzDirichletEBBCFactory::create() - logic bust, must have m_order > 0 && m_weight >= 0");
  }

  auto bc = new EBHelmholtzDirichletEBBC();

  bc->setOrder(m_order);
  bc->setWeight(m_weight);
  bc->setDomainDropOrder(m_domainDropOrder);
  bc->setCoarseGridDropOrder(m_domainDropOrder);
  if (m_useConstant) {
    bc->setValue(m_constantValue);
  }
  else if (m_useFunction) {
    bc->setValue(m_functionValue);
  }
  else {
    MayDay::Error("EBHelmholtzDirichletEBBCFactory::create() - logic bust. Not using constant or function");
  }

  return RefCountedPtr<EBHelmholtzEBBC>(bc);
}

#include <CD_NamespaceFooter.H>
