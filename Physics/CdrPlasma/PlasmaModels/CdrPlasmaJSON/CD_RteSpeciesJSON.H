/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RteSpeciesJSON.H
  @brief  Declaration of a class that encapsulates a RTE species for use with CdrPlasmaJSON
  @author Robert Marskar
*/

#ifndef CD_RteSpeciesJSON_H
#define CD_RteSpeciesJSON_H

// Std includes
#include <list>

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_RtSpecies.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    /*!
      @brief Encapsulation of radiative transfer species for usage with CdrPlasmaJSON. 
    */
    class RteSpeciesJSON : public RtSpecies {
    public:

      /*!
	@brief Alias for absorption length function. 
	@details Users must provide this in order to set the absorption coefficient (i.e., kappa) for the radiative transfer equation. 
      */
      using KappaFunction = std::function<Real(const RealVect a_pos)>;

      /*!
	@brief Default constructor. Must subsequently call define. 
      */
      RteSpeciesJSON() = default;

      /*!
	@brief Full constructor. Calls the define function. 
	@param[in] a_name        Species name
	@param[in] a_kappaFunction Initial data function
      */
      RteSpeciesJSON(const std::string&   a_name,
		     const KappaFunction& a_kappaFunction);      

      /*!
	@brief Disallowed constructor.
	@param[in] a_other Other species
      */
      RteSpeciesJSON(const RteSpeciesJSON& a_other) = delete;

      /*!
	@brief Disallowed move constructor.
	@param[in] a_other Other species
      */
      RteSpeciesJSON(const RteSpeciesJSON&& a_other) = delete;

      /*!
	@brief Disallowed copy assignment. 
	@param[in] a_other Other species
      */
      RteSpeciesJSON& operator=(const RteSpeciesJSON& a_other) = delete;

      /*!
	@brief Disallowed move assignment
	@param[in] a_other Other species
      */
      RteSpeciesJSON& operator=(const RteSpeciesJSON&& a_other) = delete;

      /*!
	@brief Destructor
      */
      virtual ~RteSpeciesJSON();

      /*!
	@brief Define function.
	@param[in] a_name          Species name
	@param[in] a_kappaFunction Function for setting absorption coefficient
      */
      void define(const std::string& a_name, const KappaFunction& a_kappaFunction);                  

      /*!
	@brief Initial data function. 
	@param[in] a_pos  Physical coordinates 
	@param[in] a_time Time
      */
      Real getAbsorptionCoefficient(const RealVect a_pos) const override final;
      
    protected:

      /*!
	@brief Is defined or not. 
      */
      bool m_isDefined;

      /*!
	@brief Absorption coefficient function. 
      */
      KappaFunction m_absorptionFunction;
    };
  }
}
  
#include <CD_NamespaceFooter.H>

#endif
