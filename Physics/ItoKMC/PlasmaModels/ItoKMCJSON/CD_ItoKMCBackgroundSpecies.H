/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCBackgroundSpecies.H
  @brief  Declaration of a background species class for usage with ItoKMCJSON
  @author Robert Marskar
*/

#ifndef CD_ItoKMCBackgroundSpecies_H
#define CD_ItoKMCBackgroundSpecies_H

// Std includes
#include <string>

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoKMC {

    /*!
      @brief Encapsulation of a background species for usage with ItoKMCJSON
    */
    class ItoKMCBackgroundSpecies
    {
    public:
      /*!
	@brief Alias for position-dependent molar fraction
      */
      using MolarFraction = std::function<Real(const RealVect a_pos)>;

      /*!
	@brief Default constructor. Must subsequently call the define function. 
      */
      ItoKMCBackgroundSpecies() noexcept;

      /*!
	@brief Full constructor. Calls the define function. 
	@param[in] a_name          Species name
	@param[in] a_molarFraction Molar fraction
      */
      ItoKMCBackgroundSpecies(const std::string& a_name, const MolarFraction& a_molarFraction) noexcept;

      /*!
	@brief Destructor
      */
      virtual ~ItoKMCBackgroundSpecies() noexcept;

      /*!
	@brief Define function. Puts object in usable state.
	@param[in] a_name          Species name
	@param[in] a_molarFraction Molar fraction
      */
      void
      define(const std::string& a_name, const MolarFraction& a_molarFraction) noexcept;

      /*!
	@brief Get molar fraction at the physical coordinates.
	@param[in] a_pos Physical position
      */
      Real
      molarFraction(const RealVect a_pos) const noexcept;

      /*!
	@brief Get the species name
	@return m_name
      */
      std::string
      getName() const noexcept;

    protected:
      /*!
	@brief Make sure model has initial data
      */
      bool m_isDefined;

      /*!
	@brief Neutral species name
      */
      std::string m_name;

      /*!
	@brief Initial data function. 
      */
      MolarFraction m_molarFraction;
    };
  } // namespace ItoKMC
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
