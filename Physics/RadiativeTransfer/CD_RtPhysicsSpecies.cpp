/*!
  @file   RtPhysicsSpecies.cpp
  @brief  Implementation of RtPhysicsSpecies.H
  @author Robert Marskar
  @date   May 2020
*/

#include <CD_RtPhysicsSpecies.H>

#include <ParmParse.H>

#include "CD_NamespaceHeader.H"
using namespace physics::rte;

RtPhysicsSpecies::RtPhysicsSpecies(){

  // This stuff sets the name and a constnat kappa taken from the input script
  m_name = "RtPhysicsSpecies";
  m_constant = true;

  ParmParse pp("RtPhysicsStepper");
  pp.get("kappa", m_kappa);
}

RtPhysicsSpecies::~RtPhysicsSpecies(){

}
#include "CD_NamespaceFooter.H"
