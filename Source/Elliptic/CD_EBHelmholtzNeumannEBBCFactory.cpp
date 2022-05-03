/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzNeumannEBBCFactory.cpp
  @brief  Implementation of CD_EBHelmholtzNeumannEBBCFactory.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBHelmholtzNeumannEBBC.H>
#include <CD_EBHelmholtzNeumannEBBCFactory.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzNeumannEBBCFactory::EBHelmholtzNeumannEBBCFactory()
{
  CH_TIME("EBHelmholtzNeumannEBBCFactory::EBHelmholtzNeumannEBBCFactory()");

  m_multByBco = true;

  m_useConstant = false;
  m_useFunction = false;
}

EBHelmholtzNeumannEBBCFactory::EBHelmholtzNeumannEBBCFactory(const Real a_DphiDn)
{
  CH_TIME("EBHelmholtzNeumannEBBCFactory::EBHelmholtzNeumannEBBCFactory(Real)");

  this->setDphiDn(a_DphiDn);
}

EBHelmholtzNeumannEBBCFactory::EBHelmholtzNeumannEBBCFactory(const std::function<Real(const RealVect& a_pos)>& a_DphiDn)
{
  CH_TIME("EBHelmholtzNeumannEBBCFactory::EBHelmholtzNeumannEBBCFactory(std::function<Real(RealVect)>)");

  this->setDphiDn(a_DphiDn);
}

EBHelmholtzNeumannEBBCFactory::~EBHelmholtzNeumannEBBCFactory()
{
  CH_TIME("EBHelmholtzNeumannEBBCFactory::~EBHelmholtzNeumannEBBCFactory()");
}

void
EBHelmholtzNeumannEBBCFactory::setDphiDn(const Real a_DphiDn)
{
  CH_TIME("EBHelmholtzNeumannEBBCFactory::setDphiDn(Real)");

  m_multByBco = true;

  m_useConstant = true;
  m_useFunction = false;

  m_constantDphiDn = a_DphiDn;
}

void
EBHelmholtzNeumannEBBCFactory::setDphiDn(const std::function<Real(const RealVect& a_pos)>& a_DphiDn)
{
  CH_TIME("EBHelmholtzNeumannEBBCFactory::setDphiDn(std::function<Real(RealVect)>)");

  m_multByBco = true;

  m_useConstant = false;
  m_useFunction = true;

  m_functionDphiDn = a_DphiDn;
}

void
EBHelmholtzNeumannEBBCFactory::setBxDphiDn(const Real a_BxDphiDn)
{
  CH_TIME("EBHelmholtzNeumannEBBCFactory::setBxDphiDn(Real)");

  this->setDphiDn(a_BxDphiDn);

  m_multByBco = false;
}

void
EBHelmholtzNeumannEBBCFactory::setBxDphiDn(const std::function<Real(const RealVect& a_pos)>& a_BxDphiDn)
{
  CH_TIME("EBHelmholtzNeumannEBBCFactory::setBxDphiDn(std::function<Real(RealVect)>)");

  this->setDphiDn(a_BxDphiDn);

  m_multByBco = false;
}

RefCountedPtr<EBHelmholtzEBBC>
EBHelmholtzNeumannEBBCFactory::create()
{
  CH_TIME("EBHelmholtzNeumannEBBCFactory::create()");

  CH_assert(m_useConstant || m_useFunction);

  auto bc = new EBHelmholtzNeumannEBBC();

  if (m_multByBco) {
    if (m_useConstant) {
      bc->setDphiDn(m_constantDphiDn);
    }
    else if (m_useFunction) {
      bc->setDphiDn(m_functionDphiDn);
    }
    else {
      MayDay::Error("EBHelmholtzNeumannEBBCFactory::create - logic bust");
    }
  }
  else {
    if (m_useConstant) {
      bc->setBxDphiDn(m_constantDphiDn);
    }
    else if (m_useFunction) {
      bc->setBxDphiDn(m_functionDphiDn);
    }
    else {
      MayDay::Error("EBHelmholtzNeumannEBBCFactory::create() - logic bust");
    }
  }

  return RefCountedPtr<EBHelmholtzEBBC>(bc);
}

#include <CD_NamespaceFooter.H>
