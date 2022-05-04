/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_NeutralSpeciesJSON.H
  @brief  Declaration of a class that encapsulates a neutral species for use with CdrPlasmaJSON.
  @author Robert Marskar
*/

#ifndef CD_NeutralSpeciesJSON_H
#define CD_NeutralSpeciesJSON_H

// Std includes
#include <list>

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    /*!
      @brief Encapsulation of species settings and initial conditiosn for use with CdrPlasma. This is for tracked species.
    */
    class NeutralSpeciesJSON
    {
    public:
      /*!
	@brief Alias for initial data function
      */
      using NumberDensityFunction = std::function<Real(const RealVect a_pos)>;

      /*!
	@brief Default constructor. Must subsequently call define. 
      */
      NeutralSpeciesJSON();

      /*!
	@brief Full constructor. Calls the define function. 
	@param[in] a_name          Species name
	@param[in] a_molarFraction Molar fraction
	@param[in] a_function      Number density function
      */
      NeutralSpeciesJSON(const std::string a_name, const Real a_molarFraction, const NumberDensityFunction a_function);

      /*!
	@brief Destructor
      */
      virtual ~NeutralSpeciesJSON();

      /*!
	@brief Full constructor. Calls the define function. 
	@param[in] a_name          Species name
	@param[in] a_molarFraction Molar fraction
	@param[in] a_function      Number density function
      */
      void
      define(const std::string a_name, const Real a_molarFraction, const NumberDensityFunction a_function);

      /*!
	@brief Get (number) density at the physical coordinates.
	@param[in] a_pos Physical position
      */
      Real
      operator()(const RealVect a_pos) const;

      /*!
	@brief Get the molar fraction
	@return m_molarFraction
      */
      Real
      getMolarFraction() const;

      /*!
	@brief Get the species name
	@return m_name
      */
      std::string
      getName() const;

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
	@brief Molar fraction
      */
      Real m_molarFraction;

      /*!
	@brief Initial data function. Set by CdrPlasmaGenericModel
      */
      NumberDensityFunction m_function;
    };
  } // namespace CdrPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
