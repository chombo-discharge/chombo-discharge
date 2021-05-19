/*!
  @file   RtPhysicsSpecies.H
  @brief  Simple species for RTE physics
  @author Robert Marskar
  @date   March 2020
*/

#ifndef _RTE_PHYSICS_SPECIES_
#define _RTE_PHYSICS_SPECIES_

#include <CD_RtSpecies.H>

#include "CD_NamespaceHeader.H"
namespace physics {
  namespace rte {

    class RtPhysicsSpecies : public RtSpecies {
    public:
      RtPhysicsSpecies();
      ~RtPhysicsSpecies();
    };
  }
}
#include "CD_NamespaceFooter.H"

#endif
