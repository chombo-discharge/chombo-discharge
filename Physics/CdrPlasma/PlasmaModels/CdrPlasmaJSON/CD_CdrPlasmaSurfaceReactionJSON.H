/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaSurfaceReactionJSON.H
  @brief  Declaration of a class that encapsulates reactions for usage with CdrPlasmaJSON
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaSurfaceReactionJSON_H
#define CD_CdrPlasmaSurfaceReactionJSON_H

// Std includes
#include <map>
#include <vector>
#include <list>

// Our includes
#include <CD_CdrPlasmaPhysics.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    /*!
      @brief Class for encapsulating reactive process of type A + B + Y -> X3 +  for CdrPlasmaPhysics. This class just stores a list
      of product target species. 
      @details This class is specifically for surface reactions. When including more species on the left hand side, e.g. A + B -> C, it is the same
      as specifying A -> C, B -> C with equal reaction rates. 
    */
    class CdrPlasmaSurfaceReactionJSON {
    public:

      /*!
	@brief Disallowed weak construction.
      */
      CdrPlasmaSurfaceReactionJSON();

      /*!
	@brief Full constructor. 
	@param[in] a_plasmaReactants  Plasma species on the left-hand side of the reaction
	@param[in] a_photonReactants  Photon species on the left hand side of the reaction.
	@param[in] a_plasmaProducts   Plasma species on the right-hand side of the reaction. 
      */
      CdrPlasmaSurfaceReactionJSON(const std::list<int> a_plasmaReactants,
				   const std::list<int> a_photonReactants,
				   const std::list<int> a_plasmaProducts);

      /*!
	@brief Destructor
      */
      virtual ~CdrPlasmaSurfaceReactionJSON();      

      /*!
	@brief Get the plasma reactants
	@return Returns m_plasmaReactants
      */
      const std::list<int>& getPlasmaReactants() const;

      /*!
	@brief Get the photo-reactants
	@return Returns m_photonReactants
      */
      const std::list<int>& getPhotonReactants() const;      

      /*!
	@brief Get the plasma products
	@return Returns m_plasmaProducts
      */
      const std::list<int>& getPlasmaProducts() const;

    protected:

      /*!
	@brief Plasma species on the left-hand side of the reaction
      */
      std::list<int> m_plasmaReactants;

      /*!
	@brief Photon species on the left-hand side of the reaction
      */
      std::list<int> m_photonReactants;            

      /*!
	@brief Plasma species on the right-hand side of the reaction
      */      
      std::list<int> m_plasmaProducts;
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
