/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaPhotoReactionJSON.H
  @brief  Declaration of a class that encapsulates reactions for usage with CdrPlasmaJSON
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaPhotoReactionJSON_H
#define CD_CdrPlasmaPhotoReactionJSON_H

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
      of target species, but it can also compute the result of stoichiometric reactions. The general format is reactions of the type [Xin] -> [Xout] + [Yout] where
      [Xin] are left hand side particle reactions, and [Xout] and [Yout] are right-hand side particle and photon products. 

      The reaction rate is computed as r = k * n[X1]n[X2]. Note that if X1 or X2 is a neutral species, they should be absorbed into the reaction rate prior to
      evaluating the reaction. 
    */
    class CdrPlasmaPhotoReactionJSON {
    public:

      /*!
	@brief Default constructor. Don't use. 
      */
      CdrPlasmaPhotoReactionJSON();

      /*!
	@brief Full constructor. 
	@param[in] a_plasmaReactants  Plasma species on the left-hand side of the reaction
	@param[in] a_neutralReactants Neutral species on the left-hand side of the reaction
	@param[in] a_photonReactants  Photon species on the left hand side of the reaction.
	@param[in] a_plasmaProducts   Plasma species on the right-hand side of the reaction. 
	@param[in] a_photonProducts   Photons on the right-hand side of the reaction.
      */
      CdrPlasmaPhotoReactionJSON(const std::list<int> a_plasmaReactants,
				 const std::list<int> a_neutralReactants,
				 const std::list<int> a_photonReactants,
				 const std::list<int> a_plasmaProducts,
				 const std::list<int> a_neutralProducts);

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
	@brief Get the photo-reactants
	@return Returns m_photonReactants
      */
      const std::list<int>& getPhotonReactants() const;      

      /*!
	@brief Get the plasma products
	@return Returns m_plasmaProducts
      */
      const std::list<int>& getPlasmaProducts() const;

      /*!
	@brief Get the neutral products
	@return Returns m_neutralProducts
      */
      const std::list<int>& getNeutralProducts() const;            
      
      /*!
	@brief Destructor
      */
      virtual ~CdrPlasmaPhotoReactionJSON();

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
	@brief Photon species on the left-hand side of the reaction
      */
      std::list<int> m_photonReactants;            

      /*!
	@brief Plasma species on the right-hand side of the reaction
      */      
      std::list<int> m_plasmaProducts;

      /*!
	@brief Neutral species on the right-hand side of the reaction
      */                  
      std::list<int> m_neutralProducts;
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
