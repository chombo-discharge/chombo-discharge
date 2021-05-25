/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaReaction.H
  @brief  Declaration of a class for holding reaction types in ito_plasma physics
  @author Robert Marskar
*/

#ifndef CD_ItoPlasmaReaction_H
#define CD_ItoPlasmaReaction_H

// Chombo includes
#include <Vector.H>
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoPlasma {
    
    /*!
      @brief Reaction class for ItoPlasmaPhysics. 
      @details To use this class you should initialize a reaction e.g. ItoPlasmaReaction({0,1,2}, {0,4}, {0}) which is equivalent
      to the reaction X_0 + X_1 + X_2 => X_1 + X_4 + Y_0. 
    */
    class ItoPlasmaReaction {
    public:

      /*!
	@brief Allow, but you need full ctor to do anything. 
      */
      ItoPlasmaReaction() = default;

      /*!
	@brief Full ctor. Includes Photons and energies. LEA fodder. 
      */
      ItoPlasmaReaction(const std::initializer_list<int>                   a_reactants,
			const std::initializer_list<int>                   a_particle_products,
			const std::initializer_list<int>                   a_photon_products,
			const std::initializer_list<std::pair<int, Real> > a_energy_jumps) : m_reactants(a_reactants),
											     m_particle_products(a_particle_products),
											     m_photon_products(a_photon_products),
											     m_energy_jumps(a_energy_jumps){
	this->computeStateChange();
	this->computeEnergyChange();
      }

      /*!
	@brief Full ctor. Includes energies but not Photons. LEA fodder. 
      */
      ItoPlasmaReaction(const std::initializer_list<int>                   a_reactants,
			const std::initializer_list<int>                   a_particle_products,
			const std::initializer_list<std::pair<int, Real> > a_energy_jumps) : m_reactants(a_reactants),
											     m_particle_products(a_particle_products),
											     m_photon_products(0),
											     m_energy_jumps(a_energy_jumps){
	this->computeStateChange();
	this->computeEnergyChange();
      }

      /*!
	@brief Full Ctor, includes possible Photons but disregards energy. LFA fodder. 
      */
      ItoPlasmaReaction(const std::initializer_list<int> a_reactants,
			const std::initializer_list<int> a_particle_products,
			const std::initializer_list<int> a_photon_products) : m_reactants(a_reactants),
									      m_particle_products(a_particle_products),
									      m_photon_products(a_photon_products),
									      m_energy_jumps(0) {
	this->computeStateChange();
	this->computeEnergyChange();
      }

      /*!
	@brief Partial ctor. No Photons generated and no energy. LFA fodder. 
      */
      ItoPlasmaReaction(const std::initializer_list<int> a_reactants,
			const std::initializer_list<int> a_particle_products) : m_reactants(a_reactants),
										m_particle_products(a_particle_products),
										m_photon_products(0),
										m_energy_jumps(0){
	this->computeStateChange();
	this->computeEnergyChange();
      }

      /*!
	@brief Dtor. Not much to see here.
      */
      ~ItoPlasmaReaction(){

      }

      /*!
	@brief Allows a state to jump with N reactions
      */
      inline
      void jumpState(Vector<long long>& a_particles, const long long a_num_reactions) const;

      /*!
	@brief Allows a state to jump with N reactions. Includes generation of Photons. 
      */
      inline
      void jumpState(Vector<long long>& a_particles,
		     Vector<long long>& a_Photons,
		     const long long    a_num_reactions) const;

      /*!
	@brief Allows the energies to jump with N reactions. Assumes no energy sources. 
      */
      inline
      void jumpEnergy(Vector<Real>& a_energies, const Vector<Real>& a_mean_energies, const long long a_num_reactions) const;

      /*!
	@brief Allows the energies to jump with N reactions. Allows energy sources. 
      */
      inline
      void jumpEnergy(Vector<Real>& a_energies, const Vector<Real>& a_mean_energies, const Vector<Real>& a_sources, const long long a_num_reactions, const Real a_dt) const;

      /*!
	@brief Return a modifiable version of the reaction rate.
      */
      inline
      Real& rate() const;

      /*!
	@brief Get propensity function.
      */
      inline
      Real propensity(const Vector<long long>& a_particles) const;

      /*!
	@brief Get gi
      */
      inline
      int getGi(const Vector<long long>& a_particles) const;

      /*!
	@brief Get reactants
      */
      inline
      const std::list<int>& getReactants() const;

      /*!
	@brief Get reactants
      */
      inline
      const std::list<int>& getParticleProdcuts() const;

      /*!
	@brief Get reactants
      */
      inline
      const std::list<int>& getPhotonProducts() const;

      /*!
	@brief Get the state change. You need to iterate through this map. 
      */
      inline
      const std::map<int, int>& getStateChange() const;

      /*!
	@brief Get the state change for a particular species
      */
      inline
      const int getStateChange(const int a_idx) const;

      /*!
	@brief Get the energy change. You need to iterate through this map. 
      */
      inline
      const std::map<int, Real>& getEnergychange() const;

      /*!
	@brief Get the energy change for a particular species
      */
      inline
      const Real getEnergychange(const int a_idx) const;
      
    protected:

      // Reaction rate, mutable because I want ItoPlasmaReaction to be a const object, but reaction rate may need to change
      mutable Real m_rate;

      std::list<int> m_reactants;                      // Particle reactants
      std::list<int> m_particle_products;              // Particle products
      std::list<int> m_photon_products;                // Photon products
      std::list<std::pair<int, Real> > m_energy_jumps; // Energy changes

      std::map<int, int>  m_phiChange;  // Change in species due to firing of one reaction. First index is species, second is change. 
      std::map<int, Real> m_energyChange; // Change in energyk due to firing of one reaction. First index is species, second is change. 

      /*!
	@brief Compute the state change
      */
      inline void computeStateChange();

      /*!
	@brief Compute the state change
      */
      inline void computeEnergyChange();
      
    };
  }
}

#include <CD_NamespaceFooter.H>

#include <CD_ItoPlasmaReactionImplem.H>

#endif
