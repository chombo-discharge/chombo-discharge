/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaGenericModel.H
  @brief  Declaration of a generic CdrPlasma model.
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaGenericModel_H
#define CD_CdrPlasmaGenericModel_H

// Third-party includes
#include <nlohmann/json.hpp>

// Our includes
#include <CD_CdrPlasmaPhysics.H>
#include <CD_CdrPlasmaSpecies.H>
#include <CD_LookupTable.H>
#include <CD_NamespaceHeader.H>

using json = nlohmann::json;

namespace Physics {
  namespace CdrPlasma {

    class CdrPlasmaGenericModel {
    public:

      /*!
	@brief Default constructor.
      */
      CdrPlasmaGenericModel();

      /*!
	@brief Destructor
      */
      virtual ~CdrPlasmaGenericModel();

    protected:

      /*!
	@brief Verbose or not
      */
      bool m_verbose;

      /*!
	@brief JSON definition. This is populated when calling parseJSON.
      */
      json m_json;      

      /*!
	@brief Input file name
      */
      std::string m_jsonFile;

      /*!
	@brief int-string encoding of the tracked species. 
	@details This is needed because we sometimes need to use the species name for indexing in the vector. 
      */
      std::map<int, std::string> m_trackedCdrSpeciesMap;

      /*!
	@brief int-string encoding of the untracked species. 
	@details This is needed because we sometimes need to use the species name for indexing in the vector. 
      */
      std::map<int, std::string> m_untrackedCdrSpeciesMap;      

      /*!
	@brief Tracked CDR species. Note that this is the one that is needed by CdrPlasmaPhysics, i.e. they become CdrSolvers. 
      */
      std::vector<RefCountedPtr<CdrSpecies> > m_trackedCdrSpecies;

      /*!
	@brief Untracked species. These do not become CdrSolvers. 
      */
      std::vector<RefCountedPtr<CdrSpecies> > m_untrackedCdrSpecies;      

      /*!
	@brief Parse class options
      */
      virtual void parseOptions();

      /*!
	@brief Parse the JSON file
      */
      virtual void parseJSON();
      
      /*!
	@brief Instantiate species
      */
      virtual void instantiateCdrSpecies();
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
