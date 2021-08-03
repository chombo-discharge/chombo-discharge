/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzEBBCFactory.cpp
  @brief  Implementation of CD_MFHelmholtzEBBCFactory.H
  @author Robert Marskar
*/

// Our includes
#include <CD_MFHelmholtzEBBCFactory.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzEBBCFactory::MFHelmholtzEBBCFactory(){

}

MFHelmholtzEBBCFactory::~MFHelmholtzEBBCFactory(){

}

void MFHelmholtzEBBCFactory::setOrder(const int a_order){
  CH_assert(a_order > 0);
  m_order = a_order;
}

void MFHelmholtzEBBCFactory::setWeight(const int a_weight){
  CH_assert(a_weight > 0);
  m_weight = a_weight;
}

#include <CD_NamespaceFooter.H>
