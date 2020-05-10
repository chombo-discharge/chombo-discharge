/*!
  @file   rte_physics_species.H
  @brief  Simple species for RTE physics
  @author Robert Marskar
  @date   March 2020
*/

#ifndef _RTE_PHYSICS_SPECIES_
#define _RTE_PHYSICS_SPECIES_

#include "rte_species.H"

namespace physics {
  namespace rte {

    class rte_physics_species : public rte_species {
    public:
      rte_physics_species();
      ~rte_physics_species();
    };
  }
}

#endif
