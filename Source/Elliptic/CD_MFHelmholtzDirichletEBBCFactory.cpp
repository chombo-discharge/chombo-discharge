/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzDirichletEBBCFactory.cpp
  @brief  Implementation of CD_MFHelmholtzDirichletEBBCFactory.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFHelmholtzDirichletEBBCFactory.H>
#include <CD_MFHelmholtzDirichletEBBC.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzDirichletEBBCFactory::MFHelmholtzDirichletEBBCFactory()
  : EBHelmholtzDirichletEBBCFactory(){
}

MFHelmholtzDirichletEBBCFactory::MFHelmholtzDirichletEBBCFactory(const int a_order, const int a_weight, const Real a_value)
  : EBHelmholtzDirichletEBBCFactory(a_order, a_weight, a_value){
}

MFHelmholtzDirichletEBBCFactory::MFHelmholtzDirichletEBBCFactory(const int a_order, const int a_weight, const std::function<Real(const RealVect& a_pos)>& a_value)
  : EBHelmholtzDirichletEBBCFactory(a_order, a_weight, a_value){
}

MFHelmholtzDirichletEBBCFactory::~MFHelmholtzDirichletEBBCFactory(){
}

RefCountedPtr<EBHelmholtzEBBC> MFHelmholtzDirichletEBBCFactory::create() {
  auto bc = new MFHelmholtzDirichletEBBC();

  bc->setOrder(m_order);
  bc->setWeight(m_weight);
  if(m_useConstant){
    bc->setValue(m_constantValue);
  }
  else if(m_useFunction){
    bc->setValue(m_functionValue);
  }
  else{
    MayDay::Error("MFHelmholtzDirichletEBBCFactory::create() - logic bust. Not using constant or function");
  }

  return RefCountedPtr<EBHelmholtzEBBC>(bc);
}

#include <CD_NamespaceFooter.H>
