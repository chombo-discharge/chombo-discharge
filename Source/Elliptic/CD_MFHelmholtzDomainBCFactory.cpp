/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_MFHelmholtzDomainBCFactory.cpp
  @brief  Implementation of CD_MFHelmholtzDomainBCFactory.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_MFHelmholtzDomainBCFactory.H>
#include <CD_NamespaceHeader.H>

MFHelmholtzDomainBCFactory::MFHelmholtzDomainBCFactory(){
  CH_TIME("MFHelmholtzDomainBCFactory::MFHelmholtzDomainBCFactory()");
}

MFHelmholtzDomainBCFactory::~MFHelmholtzDomainBCFactory(){
  CH_TIME("MFHelmholtzDomainBCFactory::~MFHelmholtzDomainBCFactory()");
}

#include <CD_NamespaceFooter.H>
