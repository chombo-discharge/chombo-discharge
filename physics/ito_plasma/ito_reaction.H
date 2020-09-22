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
	@brief Full ctor. Includes photons and energies. LEA fodder. 
      */
      ito_reaction(const std::initializer_list<int>                   a_reactants,
		   const std::initializer_list<int>                   a_particle_products,
		   const std::initializer_list<int>                   a_photon_products,
		   const std::initializer_list<std::pair<int, Real> > a_energy_jumps) : m_reactants(a_reactants),
											m_particle_products(a_particle_products),
											m_photon_products(a_photon_products),
											m_energy_jumps(a_energy_jumps){
	this->compute_state_change();
	this->compute_energy_change();
      }

      /*!
	@brief Full ctor. Includes energies but not photons. LEA fodder. 
      */
      ito_reaction(const std::initializer_list<int>                   a_reactants,
		   const std::initializer_list<int>                   a_particle_products,
		   const std::initializer_list<std::pair<int, Real> > a_energy_jumps) : m_reactants(a_reactants),
											m_particle_products(a_particle_products),
											m_photon_products(0),
											m_energy_jumps(a_energy_jumps){
	this->compute_state_change();
	this->compute_energy_change();
      }

      /*!
	@brief Full Ctor, includes possible photons but disregards energy. LFA fodder. 
      */
      ito_reaction(const std::initializer_list<int> a_reactants,
		   const std::initializer_list<int> a_particle_products,
		   const std::initializer_list<int> a_photon_products) : m_reactants(a_reactants),
									 m_particle_products(a_particle_products),
									 m_photon_products(a_photon_products),
									 m_energy_jumps(0) {
	this->compute_state_change();
	this->compute_energy_change();
      }

      /*!
	@brief Partial ctor. No photons generated and no energy. LFA fodder. 
      */
      ito_reaction(const std::initializer_list<int> a_reactants,
		   const std::initializer_list<int> a_particle_products) : m_reactants(a_reactants),
									   m_particle_products(a_particle_products),
									   m_photon_products(0),
									   m_energy_jumps(0){
	this->compute_state_change();
	this->compute_energy_change();
      }

      /*!
	@brief Dtor. Not much to see here.
      */
      ~ito_reaction(){

      }

      /*!
	@brief Allows a state to jump with N reactions
      */
      inline
      void jump_state(Vector<unsigned long long>& a_particles, const unsigned long long a_num_reactions) const;

      /*!
	@brief Allows a state to jump with N reactions. Includes generation of photons. 
      */
      inline
      void jump_state(Vector<unsigned long long>& a_particles,
		      Vector<unsigned long long>& a_photons,
		      const unsigned long long    a_num_reactions) const;

      /*!
	@brief Allows the energies to jump with N reactions
      */
      inline
      void jump_energy(Vector<Real>& a_energies, const Vector<Real>& a_mean_energies, const unsigned long long a_num_reactions) const;

      /*!
	@brief Return a modifiable version of the reaction rate.
      */
      inline
      Real& rate() const;

      /*!
	@brief Get propensity function.
      */
      inline
      Real propensity(const Vector<unsigned long long>& a_particles) const;

      /*!
	@brief Get gi
      */
      inline
      int get_gi(const Vector<unsigned long long>& a_particles) const;

      /*!
	@brief Get reactants
      */
      inline
      const std::list<int>& get_reactants() const;

      /*!
	@brief Get reactants
      */
      inline
      const std::list<int>& get_particle_products() const;

      /*!
	@brief Get reactants
      */
      inline
      const std::list<int>& get_photon_products() const;

      /*!
	@brief Get the state change. You need to iterate through this map. 
      */
      inline
      const std::map<int, int>& get_state_change() const;

      /*!
	@brief Get the state change for a particular species
      */
      inline
      const int get_state_change(const int a_idx) const;

      /*!
	@brief Get the energy change. You need to iterate through this map. 
      */
      inline
      const std::map<int, Real>& get_energy_change() const;

      /*!
	@brief Get the energy change for a particular species
      */
      inline
      const Real get_energy_change(const int a_idx) const;
      
    protected:

      // Reaction rate, mutable because I want ito_reaction to be a const object, but reaction rate may need to change
      mutable Real m_rate;

      std::list<int> m_reactants;                      // Particle reactants
      std::list<int> m_particle_products;              // Particle products
      std::list<int> m_photon_products;                // Photon products
      std::list<std::pair<int, Real> > m_energy_jumps; // Energy changes

      std::map<int, int>  m_stateChange;  // Change in species due to firing of one reaction. First index is species, second is change. 
      std::map<int, Real> m_energyChange; // Change in energyk due to firing of one reaction. First index is species, second is change. 

      /*!
	@brief Compute the state change
      */
      inline void compute_state_change();

      /*!
	@brief Compute the state change
      */
      inline void compute_energy_change();
      
    };
  }
}

#include "ito_reactionI.H"

#endif
