/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
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

    /*!
      @brief Implementation of RtSpecies for usage in RadiativeTransfer module.
      @details This reads a constant kappa from the input file which it uses as Beer's length.
    */
    class RtPhysicsSpecies : public RtSpecies {
    public:

      /*!
	@brief Constructor (reads from input script)
      */
      RtPhysicsSpecies();

      /*!
	@brief Destructor (does nothing)
      */      
      ~RtPhysicsSpecies();
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
