/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_RtPhysicsSpecies.cpp
  @brief  Implementation of CD_RtPhysicsSpecies.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_RtPhysicsSpecies.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::RadiativeTransfer;

RtPhysicsSpecies::RtPhysicsSpecies(){

  // This stuff sets the name and a constnat kappa taken from the input script
  m_name = "RtPhysicsSpecies";
  m_constant = true;

  ParmParse pp("RtPhysicsStepper");
  pp.get("kappa", m_kappa);
}

RtPhysicsSpecies::~RtPhysicsSpecies(){

}

#include <CD_NamespaceFooter.H>
