/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RadiativeTransferSpecies.cpp
  @brief  Implementation of CD_RadiativeTransferSpecies.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>

// Our includes
#include <CD_RadiativeTransferSpecies.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::RadiativeTransfer;

RadiativeTransferSpecies::RadiativeTransferSpecies(){
  CH_TIME("RadiativeTransferSpecies::RadiativeTransferSpecies");
  
  // This stuff sets the name and a constnat kappa taken from the input script
  m_name     = "RadiativeTransferSpecies";
  m_constant = true;

  ParmParse pp("RadiativeTransferStepper");
  pp.get("kappa", m_kappa);
}

RadiativeTransferSpecies::~RadiativeTransferSpecies(){
  CH_TIME("RadiativeTransferSpecies::~RadiativeTransferSpecies");
}

#include <CD_NamespaceFooter.H>
