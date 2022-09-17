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

RadiativeTransferSpecies::RadiativeTransferSpecies()
{
  CH_TIME("RadiativeTransferSpecies::RadiativeTransferSpecies");

  // This stuff sets the name and a constnat kappa taken from the input script
  m_name = "RadiativeTransferSpecies";

  // Set constant kappa.
  ParmParse pp("RadiativeTransferStepper");

  Real kappa;
  pp.get("kappa", kappa);

  // Set the spatially varying absorption coefficient to a constant.
  m_kappa = [kappa](const RealVect a_pos) -> Real {
    return kappa;
  };
}

RadiativeTransferSpecies::~RadiativeTransferSpecies()
{
  CH_TIME("RadiativeTransferSpecies::~RadiativeTransferSpecies");
}

Real
RadiativeTransferSpecies::getAbsorptionCoefficient(const RealVect a_pos) const
{
  return m_kappa(a_pos);
}

#include <CD_NamespaceFooter.H>
