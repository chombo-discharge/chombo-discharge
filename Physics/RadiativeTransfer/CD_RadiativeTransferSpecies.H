/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RadiativeTransferSpecies.H
  @brief  Simple species for Physics/RadiativeTransfer
  @author Robert Marskar
  @date   March 2020
*/

#ifndef CD_RadiativeTransferSpecies_H
#define CD_RadiativeTransferSpecies_H

// Our includes
#include <CD_RtSpecies.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace RadiativeTransfer {

    /*!
      @brief Implementation of RtSpecies for usage in RadiativeTransfer module.
      @details This reads a constant kappa from the input file which it uses as Beer's length.
    */
    class RadiativeTransferSpecies : public RtSpecies {
    public:

      /*!
	@brief Constructor (reads from input script)
      */
      RadiativeTransferSpecies();

      /*!
	@brief Destructor (does nothing)
      */      
      ~RadiativeTransferSpecies();
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
