/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaReactionJSON.H
  @brief  Declaration of a class that encapsulates reactions for usage with CdrPlasmaJSON
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaReactionJSON_H
#define CD_CdrPlasmaReactionJSON_H

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
      @brief Class for encapsulating reactive process of type X1 + X2 -> X3 + Y for CdrPlasmaPhysics. This class just stores a list
      of product and target species.
    */
    class CdrPlasmaReactionJSON {
    public:

      /*!
	@brief Disallowed weak construction.
      */
      CdrPlasmaReactionJSON() = delete;

      /*!
	@brief Full constructor. 
	@param[in] a_plasmaReactants  Plasma species on the left-hand side of the reaction
	@param[in] a_neutralReactants Neutral species on the left-hand side of the reaction
	@param[in] a_plasmaProducts   Plasma species on the right-hand side of the reaction. 
	@param[in] a_photonProducts   Photons on the right-hand side of the reaction.
      */
      CdrPlasmaReactionJSON(const std::list<int> a_plasmaReactants,
			    const std::list<int> a_neutralReactants,			    
			    const std::list<int> a_plasmaProducts,
			    const std::list<int> a_photonProducts);

      /*!
	@brief Destructor
      */
      virtual ~CdrPlasmaReactionJSON();      

      /*!
	@brief Get the plasma reactants
	@return Returns m_plasmaReactants
      */
      const std::list<int>& getPlasmaReactants() const;

      /*!
	@brief Get the neutral reactants
	@return Returns m_neutralReactants
      */
      const std::list<int>& getNeutralReactants() const;

      /*!
	@brief Get the plasma products
	@return Returns m_plasmaProducts
      */
      const std::list<int>& getPlasmaProducts() const;

      /*!
	@brief Get the photon products
	@return Returns m_photonProducts
      */
      const std::list<int>& getPhotonProducts() const;            

    protected:

      /*!
	@brief Plasma species on the left-hand side of the reaction
      */
      std::list<int> m_plasmaReactants;

      /*!
	@brief Neutral species on the left-hand side of the reaction
      */
      std::list<int> m_neutralReactants;      

      /*!
	@brief Particles on the right-hand side of the reaction
      */      
      std::list<int> m_plasmaProducts;

      /*!
	@brief Photons on the right-hand side of the reaction
      */                  
      std::list<int> m_photonProducts;
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
