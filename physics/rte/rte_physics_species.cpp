/*!
  @file   rte_physics_species.cpp
  @brief  Implementation of rte_physics_species.H
  @author Robert Marskar
  @date   May 2020
*/

#include "rte_physics_species.H"

#include <ParmParse.H>

namespace ChomboDischarge {
  using namespace physics::rte;

  rte_physics_species::rte_physics_species(){

    // This stuff sets the name and a constnat kappa taken from the input script
    m_name = "rte_physics_species";
    m_constant = true;

    ParmParse pp("rte_stepper");
    pp.get("kappa", m_kappa);
  }

  rte_physics_species::~rte_physics_species(){

  }
}
