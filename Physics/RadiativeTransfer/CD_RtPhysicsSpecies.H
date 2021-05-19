/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_RtPhysicsSpecies.H
  @brief  Simple species for Physics/RadiativeTransfer
  @author Robert Marskar
  @date   March 2020
*/

#ifndef CD_RtPhysicsSpecies_H
#define CD_RtPhysicsSpecies_H

// Our includes
#include <CD_RtSpecies.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace RadiativeTransfer {

    class RtPhysicsSpecies : public RtSpecies {
    public:
      RtPhysicsSpecies();
      ~RtPhysicsSpecies();
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
