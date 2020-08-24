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

      /*!
	@brief Draw a valid random position somewhere in a cell
      */
      RealVect random_position(const RealVect a_pos,
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
      RealVect random_position(const RealVect a_lo,
			       const RealVect a_hi,
			       const RealVect a_bndryCentroid,
			       const RealVect a_normal) const;

      /*!
	@brief Return a random position in the cube (a_lo, a_hi);
      */
      RealVect random_position(const RealVect a_lo, const RealVect a_hi) const;

      /*!
	@brief Draw a random direction in space
      */
      RealVect random_direction() const;

      /*!
	@brief Check if advance is valid. This returns false if there is a negative number of particles in the cell
      */
      bool valid_advance(const Vector<List<ito_particle>* >& a_particles,
			 const Vector<int>&                  a_newParticles) const;

      /*!
	@brief Add particles
      */
      void generate_particles(List<ito_particle>& a_particles,
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
      void add_particles(List<ito_particle>& a_particles,
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
      void remove_particles(List<ito_particle>& a_particles,
			    const int           a_numParticlesToRemove,
			    const RealVect      a_pos,
			    const RealVect      a_lo,
			    const RealVect      a_hi,
			    const RealVect      a_bndryCentroid,
			    const RealVect      a_normal,
			    const Real          a_dx,
			    const Real          a_kappa) const;

      /*!
	@brief Compute particle weights for superparticles
      */
      void compute_particle_weights(int& a_weight, int& a_num, int& remainder, const int a_numNewParticles) const;

    };
  }
}
#endif
