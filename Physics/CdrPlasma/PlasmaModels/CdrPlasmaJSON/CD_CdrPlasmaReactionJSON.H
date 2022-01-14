/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
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
      of target species, but it can also compute the result of stoichiometric reactions. The general format is reactions of the type [Xin] -> [Xout] + [Yout] where
      [Xin] are left hand side particle reactions, and [Xout] and [Yout] are right-hand side particle and photon products. 

      The reaction rate is computed as r = k * n[X1]n[X2]. Note that if X1 or X2 is a neutral species, they should be absorbed into the reaction rate prior to
      evaluating the reaction. 
    */
    class CdrPlasmaReactionJSON {
    public:

      /*!
	@brief Default constructor. Must subsequently call the define method
      */
      CdrPlasmaReactionJSON();

      /*!
	@brief Full constructor. 
	@param[in] a_plasmaReactants Plasma species on the left-hand side of the reaction
	@param[in] a_plasmaReactants Neutral species on the left-hand side of the reaction
	@param[in] a_plasmaProducts  Plasma species on the right-hand side of the reaction. 
	@param[in] a_photonProducts  Photons on the right-hand side of the reaction.
      */
      CdrPlasmaReactionJSON(const std::list<int> a_plasmaReactants,
			    const std::list<int> a_neutralReactants,			    
			    const std::list<int> a_plasmaProducts,
			    const std::list<int> a_photonProducts);

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
      
      /*!
	@brief Add the result of this reaction to the product vectors. This will ONLY iterate over the plasma species so if you had 
	neutral species involved, their contribution should be incorporated in a_reactionRate. 
	@param[inout] a_plasmaSources   On output, entries will be incremented by the number of consumed/produced particle products. 
	@param[inout] a_photonSources   On output, entries will be incremented by the number of produced photon products. 
	@param[inout] a_plasmaDensities Input particle densities. Note, these are evaluated per unit volume. 
	@param[in]    a_reactionRate    Basic reaction rate for reaction X1 + X2 -> X3 + Y etc. 
      */
      void addProducts(std::vector<Real>&       a_plasmaProducts,
		       std::vector<Real>&       a_photonProducts,
		       const std::vector<Real>& a_plasmaDensities,
		       const Real               a_reactionRate) const;
      
      /*!
	@brief Destructor
      */
      virtual ~CdrPlasmaReactionJSON();

    protected:

      /*!
	@brief Is defined or not
      */
      bool m_isDefined;

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
