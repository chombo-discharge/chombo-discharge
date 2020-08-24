/*!
  @file   ito_reaction.H
  @brief  Declaration of a class for holding reaction types in ito_plasma physics
  @author Robert Marskar
  @date   Aug. 2020
*/

#ifndef _ITO_REACTION_H_
#define _ITO_REACTION_H_

#include <Vector.H>
#include <RealVect.H>

namespace physics {
  namespace ito_plasma {
    
    /*!
      @brief Reaction class for ito_plasma_physics. 
      @details To use this class you should initialize a reaction e.g. ito_reaction({0,1,2}, {0,4}, {0}) which is equivalent
      to the reaction X_0 + X_1 + X_2 => X_1 + X_4 + Y_0. 
    */
    class ito_reaction {
    public:

      /*!
	@brief Allow, but you need full ctor to do anything. 
      */
      ito_reaction() = default;

      /*!
	@brief Full Ctor, includes possible photons. 
      */
      ito_reaction(const std::initializer_list<int> a_reactants,
		   const std::initializer_list<int> a_particle_products,
		   const std::initializer_list<int> a_photon_products) : m_reactants(a_reactants),
									 m_particle_products(a_particle_products),
									 m_photon_products(a_photon_products){

      }

      /*!
	@brief Partial ctor. No photons generated. 
      */
      ito_reaction(const std::initializer_list<int> a_reactants,
		   const std::initializer_list<int> a_particle_products) : m_reactants(a_reactants),
									   m_particle_products(a_particle_products),
									   m_photon_products(0) {

      }

      /*!
	@brief Dtor. Not much to see here.
      */
      ~ito_reaction(){

      }

      /*!
	@brief Allows a state to jump with N reactions
      */
      inline void jump_state(Vector<int>& a_particles, const int a_num_reactions) const;

      /*!
	@brief Allows a state to jump with N reactions. Includes generation of photons. 
      */
      inline void jump_state(Vector<int>& a_particles, Vector<int>& a_photons, const int a_num_reactions) const;

      /*!
	@brief Return a modifiable version of the reaction rate.
      */
      inline Real& rate() const;

      /*!
	@brief Get propensity function.
      */
      inline Real propensity(const Vector<int>& a_particles) const;
      
    protected:

      // Reaction rate, mutable because I want ito_reaction to be a const object, but reaction rate may need to change
      mutable Real m_rate;

      std::vector<int> m_reactants;         // Particle reactants
      std::vector<int> m_particle_products; // Particle products
      std::vector<int> m_photon_products;   // Photon products
    };
  }
}

#include "ito_reactionI.H"

#endif
