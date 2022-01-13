/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrSpeciesJSON.H
  @brief  Declaration of a class that encapsulates a CDR species for use with CdrPlasmaJSON
  @author Robert Marskar
*/

#ifndef CD_CdrSpeciesJSON_H
#define CD_CdrSpeciesJSON_H

// Std includes
#include <list>

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_CdrSpecies.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    /*!
      @brief Encapsulation of species settings and initial conditions for CDR species in CdrPlasma. 
    */
    class CdrSpeciesJSON : public CdrSpecies {
    public:

      /*!
	@brief Alias for initial data function.
      */
      using InitFunction = std::function<Real(const RealVect a_pos, const Real a_time)>;

      /*!
	@brief Default constructor. Must subsequently call define.
      */
      CdrSpeciesJSON() = default;

      /*!
	@brief Full constructor. Calls the define function. 
	@param[in] a_name        Species name
	@param[in] a_Z           Charge number
	@param[in] a_mobile      Mobile species or not
	@param[in] a_diffusive   Diffusive species or not
	@param[in] a_initialData Initial data function
      */
      CdrSpeciesJSON(const std::string  a_name,
		     const int          a_Z,
		     const bool         a_diffusive,
		     const bool         a_mobile,
		     const InitFunction a_initialData);      

      /*!
	@brief Disallowed constructor.
	@param[in] a_other Other species
      */
      CdrSpeciesJSON(const CdrSpeciesJSON& a_other) = delete;

      /*!
	@brief Disallowed move constructor.
	@param[in] a_other Other species
      */
      CdrSpeciesJSON(const CdrSpeciesJSON&& a_other) = delete;

      /*!
	@brief Disallowed copy assignment. 
	@param[in] a_other Other species
      */
      CdrSpeciesJSON& operator=(const CdrSpeciesJSON& a_other) = delete;

      /*!
	@brief Disallowed move assignment
	@param[in] a_other Other species
      */
      CdrSpeciesJSON& operator=(const CdrSpeciesJSON&& a_other) = delete;

      /*!
	@brief Destructor
      */
      virtual ~CdrSpeciesJSON();

      /*!
	@brief Define function.
	@param[in] a_name        Species name
	@param[in] a_Z           Charge number
	@param[in] a_mobile      Mobile species or not
	@param[in] a_diffusive   Diffusive species or not
	@param[in] a_initialData Initial data function
      */
      void define(const std::string  a_name,
		  const int          a_Z,
		  const bool         a_diffusive,
		  const bool         a_mobile,
		  const InitFunction a_initialData);                  

      /*!
	@brief Initial data function. 
	@param[in] a_pos  Physical coordinates 
	@param[in] a_time Time
      */
      Real initialData(const RealVect a_pos, const Real a_time) const override final;
      
    protected:

      /*!
	@brief Make sure model has initial data
      */
      bool m_isDefined;

      /*!
	@brief Initial data function. 
      */
      InitFunction m_initFunction;
    };
  }
}
  
#include <CD_NamespaceFooter.H>

#endif
