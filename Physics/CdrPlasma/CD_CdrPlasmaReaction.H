/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrPlasmaReaction.H
  @brief  Declaration of a class that encapsulates reactions
  @author Robert Marskar
*/

#ifndef CD_CdrPlasmaReaction_H
#define CD_CdrPlasmaReaction_H

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
      @brief Class for encapsulating reactive process of type X1 + X2 -> X3 + Y for CdrPlasmaPhysics. 
      @details This only processes stoichiometric reactions. The general format is reactions of the type [Xin] -> [Xout] + [Yout] where
      [Xin] are left hand side particle reactions, and [Xout] and [Yout] are right-hand side particle and photon products. 

      The reaction rate is computed as r = k * n[X1]n[X2]. Note that if X1 or X2 is a neutral species, they should be absorbed into the reaction rate prior to
      evaluating the reaction. 

      When using this class, objects can be easily constructed as CdrPlasmaReaction myReaction({0,1}, {0,0,2}, {0}) which is equivalent to 
      the reaction X0 + X1 -> X0 + X0 + X2 and no photon reactants or products. 

      Note that, for generality, this class using integer indexing of the involved species. 
    */
    class CdrPlasmaReaction {
    public:

      /*!
	@brief Default constructor. Must subsequently call the define method
      */
      CdrPlasmaReaction();

      /*!
	@brief Full constructor. 
	@param[in] a_particleReactants Particles on the left-hand side of the reaction. 
	@param[in] a_particleProducts  Particles on the right-hand side of the reaction. 
	@param[in] a_photonProducts    Photons on the right-hand side of the reaction.
      */
      CdrPlasmaReaction(const std::list<int> a_particleReactants,
			const std::list<int> a_particleProducts,
			const std::list<int> a_photonProducts);            

      /*!
	@brief Full constructor. Permits construction as CdrPlasmaReaction({0,0,1}, {0,2}) etc. 
	@param[in] a_particleReactants Particles on the left-hand side of the reaction. 
	@param[in] a_particleProducts  Particles on the right-hand side of the reaction. 
	@param[in] a_photonProducts    Photons on the right-hand side of the reaction.
      */
      CdrPlasmaReaction(const std::initializer_list<int> a_particleReactants,
			const std::initializer_list<int> a_particleProducts,
			const std::initializer_list<int> a_photonProducts);

      /*!
	@brief Full constructor, assuming no photon products. 
	@param[in] a_particleReactants Particles on the left-hand side of the reaction. 
	@param[in] a_particleProducts  Particles on the right-hand side of the reaction. 
      */
      CdrPlasmaReaction(const std::list<int> a_particleReactants,
			const std::list<int> a_particleProducts);      

      /*!
	@brief Full constructor, assuming no photon products. Permits construction as CdrPlasmaReaction({0,0,1}, {0,2}) etc. 
	@param[in] a_particleReactants Particles on the left-hand side of the reaction. 
	@param[in] a_particleProducts  Particles on the right-hand side of the reaction. 
      */
      CdrPlasmaReaction(const std::initializer_list<int> a_particleReactants,
			const std::initializer_list<int> a_particleProducts);


      /*!
	@brief Get reaction rate
	@return Returns m_reactionRate.
      */
      Real& rate() const;

      /*!
	@brief Fire a reaction. This will fire the reaction using m_reactionRate. 
      */
      void fire(std::vector<Real>&       a_particleSources,
		std::vector<Real>&       a_photonSources,
		const std::vector<Real>& a_particleDensities) const;
      
      /*!
	@brief Fire a reaction. This will fire the reaction using the input rate. 
	@details Beware of the reaction rate units: If two-body reactions like X1 + X2 -> null are involved, the source 
	@param[inout] a_particleSources   On output, entries will be incremented by the number of consumed/produced particle products. 
	@param[inout] a_photonSources     On output, entries will be incremented by the number of produced photon products. 
	@param[inout] a_particleDensities Input particle densities. Note, these are evaluated per unit volume. 
	@param[in]    a_reactionRate      Basic reaction rate for reaction X1 + X2 -> X3 + Y etc. 
      */
      void fire(std::vector<Real>&       a_particleSources,
		std::vector<Real>&       a_photonSources,
		const std::vector<Real>& a_particleDensities,
		const Real               a_reactionRate) const;

      /*!
	@brief Destructor
      */
      virtual ~CdrPlasmaReaction();

    protected:

      /*!
	@brief Is defined or not
      */
      bool m_isDefined;

      /*!
	@brief Reaction rate
      */
      mutable Real m_reactionRate;

      /*!
	@brief Particles on the left-hand side of the reaction
      */
      std::list<int> m_particleReactants;

      /*!
	@brief Particles on the right-hand side of the reaction
      */      
      std::list<int> m_particleProducts;

      /*!
	@brief Photons on the right-hand side of the reaction
      */                  
      std::list<int> m_photonProducts;

      /*!
	@brief Change in the particle numbers due to firing of exactly one reaction
      */
      std::map<int, int> m_particleJump;

      /*!
	@brief Compute state changes
      */
      virtual void computeStateChanges();
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
