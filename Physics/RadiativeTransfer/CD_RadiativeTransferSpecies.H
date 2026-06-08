/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_RadiativeTransferSpecies.H
  @brief  Declaration of the Physics::RadiativeTransfer::RadiativeTransferSpecies RtSpecies
  @author Robert Marskar
*/

#ifndef CD_RADIATIVETRANSFERSPECIES_H
#define CD_RADIATIVETRANSFERSPECIES_H

// Our includes
#include <CD_RtSpecies.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace RadiativeTransfer {

    /**
      @brief Implementation of RtSpecies for usage in RadiativeTransfer module.
      @details This reads a constant kappa from the input file which it uses as Beer's length.
    */
    class RadiativeTransferSpecies : public RtSpecies
    {
    public:
      /**
	@brief Constructor (reads from input script)
      */
      RadiativeTransferSpecies();

      /**
	@brief Destructor (does nothing)
      */
      virtual ~RadiativeTransferSpecies();

      /**
	@brief Get absorption coefficient (i.e. the inverse absorption length) at physical coordinates.
	@param[in] a_pos Physical coordinates.
	@return Absorption coefficient (inverse absorption length) at a_pos.
      */
      virtual Real
      getAbsorptionCoefficient(const RealVect a_pos) const override;

    protected:
      /**
	@brief Absorption coefficient.
	@details Storing this as a std::function to show to easily use lambdas for returning spatially varying coefficients. 
      */
      std::function<Real(const RealVect a_pos)> m_kappa;
    };
  } // namespace RadiativeTransfer
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
