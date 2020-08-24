/*!
  @file   ito_plasma_physics.H
  @brief  Main file for describing Ito-based plasma physics
  @author Robert Marskar
  @date   June 2020
*/

#ifndef _ITO_PLASMA_PHYSICS_
#define _ITO_PLASMA_PHYSICS_

#include <RealVect.H>
#include <RefCountedPtr.H>
#include <LoHiSide.H>
#include <List.H>

#include <chrono>
#include <random>

#include "ito_species.H"
#include "rte_species.H"
#include "photon.H"
#include "ito_particle.H"

namespace physics {
  namespace ito_plasma {

    /*!
      @brief Abstract clas for specifying plasma kinetics
    */
    class ito_plasma_physics {
    public:
      
      /*!
	@brief Ctor. Does nothing.
      */
      ito_plasma_physics();

      /*!
	@brief Dtor. Does nothing. 
      */
      ~ito_plasma_physics();

      /*!
	@brief Compute Townsend ionization length
      */
      virtual Real compute_alpha(const RealVect a_E) const = 0;

      /*!
	@brief Get all species
      */
      const Vector<RefCountedPtr<ito_species> >& get_ito_species() const;

      /*!
	@brief Get all photons
      */
      const Vector<RefCountedPtr<rte_species> >& get_rte_species() const;

      /*!
	@brief Return number of ion equations
      */
      int get_num_ito_species() const;

      /*!
	@brief Return number of RTE equations
      */
      int get_num_rte_species() const;

      /*!
	@brief Set initial surface charge. Default is 0, override if you want.
      */
      virtual Real initial_sigma(const Real a_time, const RealVect a_pos) const;

      /*!
	@brief Compute ito velocities
      */
      virtual Vector<RealVect> compute_ito_velocities(const Real         a_time,
						      const RealVect     a_pos,
						      const RealVect     a_E,
						      const Vector<Real> a_cdr_densities) const = 0;

      /*!
	@brief Compute diffusion coefficients
      */
      virtual Vector<Real> compute_ito_diffusion(const Real         a_time,
						 const RealVect     a_pos,
						 const RealVect     a_E,
						 const Vector<Real> a_cdr_densities) const = 0;

      /*!
	@brief Reaction network advance
      */
      virtual void advance_reaction_network(Vector<List<ito_particle>* >& a_particles,
					    Vector<List<photon>* >&       a_photons,
					    Vector<List<photon>* >&       a_newPhotons,
					    const RealVect                a_E,           
					    const RealVect                a_cellPos,
					    const RealVect                a_centroidPos,
					    const RealVect                a_bndryCentroid,
					    const RealVect                a_bndryNormal,
					    const RealVect                a_lo,
					    const RealVect                a_hi,
					    const Real                    a_dx,
					    const Real                    a_kappa, 
					    const Real                    a_dt) const = 0;

      
    protected:

      int m_seed;
      mutable std::mt19937_64 m_rng;
      mutable std::uniform_real_distribution<Real> m_udist11;

      Vector<RefCountedPtr<ito_species> > m_ito_species;
      Vector<RefCountedPtr<rte_species> > m_rte_species;

      int m_num_ito_species;
      int m_num_rte_species;
      int m_ppc;
      int m_poisson_switch;
      int m_init_weights;

      /*!
	@brief Draw a valid random position somewhere in a cell
      */
      inline RealVect random_position(const RealVect a_pos,
				      const RealVect a_lo,
				      const RealVect a_hi,
				      const RealVect a_bndryCentroid,
				      const RealVect a_normal,
				      const Real     a_dx,
				      const Real     a_kappa) const;

      /*!
	@brief Draw a random position somewhere in the cube (a_lo, a_hi) such that the position is
	(pos-bndryCentroid).dotProduct.a_normal > 0, i.e. inside the valid region. 
      */
      inline RealVect random_position(const RealVect a_lo,
				      const RealVect a_hi,
				      const RealVect a_bndryCentroid,
				      const RealVect a_normal) const;

      /*!
	@brief Return a random position in the cube (a_lo, a_hi);
      */
      inline RealVect random_position(const RealVect a_lo, const RealVect a_hi) const;

      /*!
	@brief Draw a random direction in space
      */
      inline RealVect random_direction() const;

      /*!
	@brief Check if advance is valid. This returns false if there is a negative number of particles in the cell
      */
      inline bool valid_advance(const Vector<List<ito_particle>* >& a_particles,
				const Vector<int>&                  a_newParticles) const;

      /*!
	@brief Add particles
      */
      inline void generate_particles(List<ito_particle>& a_particles,
				     const int           a_numNewParticles,
				     const RealVect      a_pos,
				     const RealVect      a_lo,
				     const RealVect      a_hi,
				     const RealVect      a_bndryCentroid,
				     const RealVect      a_normal,
				     const Real          a_dx,
				     const Real          a_kappa) const;

      /*!
	@brief Add particles
      */
      inline void add_particles(List<ito_particle>& a_particles,
				const int           a_numNewParticles,
				const RealVect      a_pos,
				const RealVect      a_lo,
				const RealVect      a_hi,
				const RealVect      a_bndryCentroid,
				const RealVect      a_bndryNormal,
				const Real          a_dx,
				const Real          a_kappa) const;

      /*!
	@brief Remove particles
      */
      inline void remove_particles(List<ito_particle>& a_particles,
				   const int           a_numParticlesToRemove,
				   const RealVect      a_pos,
				   const RealVect      a_lo,
				   const RealVect      a_hi,
				   const RealVect      a_bndryCentroid,
				   const RealVect      a_normal,
				   const Real          a_dx,
				   const Real          a_kappa) const;

      /*!
	@brief Add photoionization for specific transition
      */
      inline void add_photoionization(List<ito_particle>& a_particles1,
				      List<ito_particle>& a_particles2,
				      const List<photon>& a_photons) const;

      /*!
	@brief Add photons
      */
      inline void add_photons(List<photon>&      a_photons,
			      const rte_species& a_species,
			      const int          a_num_photons,
			      const RealVect     a_pos,
			      const RealVect     a_lo,
			      const RealVect     a_hi,
			      const RealVect     a_bndryCentroid,
			      const RealVect     a_bndryNormal,
			      const Real         a_dx,
			      const Real         a_kappa) const;

      /*!
	@brief Compute particle weights for superparticles
      */
      inline void compute_particle_weights(int& a_weight, int& a_num, int& remainder, const int a_numNewParticles) const;

      /*!
	@brief Poisson reaction
      */
      inline int poisson_reaction(const Real a_propensity, const Real a_dt) const;

      /*!
	@brief Draw Gaussian distributed particles
      */
      inline void draw_gaussian_particles(List<ito_particle>& a_particles1,
					  List<ito_particle>& a_particles2,
					  const int           a_numParticles,
					  const RealVect      a_blobCenter,
					  const Real          a_blobRadius,
					  const Real          a_weights);
    };

    /*!
      @brief Reaction class for ito_plasma_physics
    */
    class reaction {
    public:

      reaction () = delete;

      /*!
	@brief Full Ctor. This is the only one you will need. 
      */
      reaction(const Vector<int>& a_reactants,
	       const Vector<int>& a_particle_products,
	       const Vector<int>& a_photon_products){
	this->define_reaction(a_reactants, a_particle_products, a_photon_products);
      }
      
      ~reaction (){}

      inline Real& rate() const {
	return m_rate;
      }

      /*!
	@brief Note that 
      */
      inline Real propensity(const Vector<int>& a_particles) const {
	Real a = m_rate;

	// Make a local copy of the particles
	Vector<int> particles = a_particles;

	// Iterate
	for (int i = 0; i < m_reactants.size(); i++){
	  const int reactant = m_reactants[i];
	  
	  a *= particles[reactant];
	  particles[reactant]--;
	}

	return Max(0.0, a);
      }

      /*!
	@brief Get the reactants
      */
      inline const Vector<int>& get_reactants() const {
	return m_reactants;
      }

      /*!
	@brief Get the particle products
      */
      inline const Vector<int>& get_particle_products() const {
	return m_particle_products;
      }

      /*!
	@brief Get the photon products
      */
      inline const Vector<int>& get_photon_products() const {
	return m_photon_products;
      }

      /*!
	@brief Define function
      */
      inline void define_reaction(const Vector<int>& a_reactants,
				  const Vector<int>& a_particle_products,
				  const Vector<int>& a_photon_products){
	this->set_reactants(a_reactants);
	this->set_particle_products(a_particle_products);
	this->set_photon_products(a_photon_products);
      }
      
    protected:

      mutable Real m_rate;

      Vector<int> m_reactants;
      Vector<int> m_particle_products;
      Vector<int> m_photon_products;

      inline void set_reactants(const Vector<int>& a_reactants){
	m_reactants = a_reactants;
      }

      inline void set_particle_products(const Vector<int>& a_particle_products){
	m_particle_products = a_particle_products;
      }

      inline void set_photon_products(const Vector<int>& a_photon_products){
	m_photon_products = a_photon_products;
      }
    };
  }
}

#include "ito_plasma_physicsI.H"

#endif
