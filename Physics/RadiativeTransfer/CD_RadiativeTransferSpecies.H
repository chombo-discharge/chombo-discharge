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
      virtual ~RadiativeTransferSpecies();

      /*!
	@brief Get absorption coefficient (i.e. the inverse absorption length) at physical coordinates. 
	@param[in] a_pos Physical coordinates. 
      */
      virtual Real getAbsorptionCoefficient(const RealVect a_pos) const override;

    protected:

      /*!
	@brief Absorption coefficient.
	@details Storing this as a std::function to show to easily use lambdas for returning spatially varying coefficients. 
      */
      std::function<Real(const RealVect a_pos)> m_kappa;
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
