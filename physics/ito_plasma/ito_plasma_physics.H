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
#include <memory>
#include <deque>

#include "ito_species.H"
#include "rte_species.H"
#include "photon.H"
#include "ito_particle.H"
#include "ito_reaction.H"
#include "photo_reaction.H"

namespace physics {
  namespace ito_plasma {

    // Ease up on the typing.
    using ReactionPtrs  = std::deque<std::shared_ptr<ito_reaction> >;
    using ReactionTypes = std::pair<ReactionPtrs, ReactionPtrs>;

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
	@brief Update reaction rates
      */
      virtual void update_reaction_rates(const RealVect a_E, const Real a_dx, const Real a_kappa) const = 0;
      
      /*!
	@brief Reaction network advance. Trying to inline for performance....
      */
      inline
      void advance_reaction_network(Vector<List<ito_particle>* >& a_particles,
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
				    const Real                    a_dt);

      
    protected:

      int m_seed;
      mutable std::mt19937_64 m_rng;
      mutable std::uniform_real_distribution<Real> m_udist01;
      mutable std::uniform_real_distribution<Real> m_udist11;

      Vector<RefCountedPtr<ito_species> > m_ito_species;
      Vector<RefCountedPtr<rte_species> > m_rte_species;

      int m_num_ito_species;
      int m_num_rte_species;
      int m_ppc;
      int m_poisson_switch;
      int m_init_weights;
      int m_Ncrit;

      Real m_eps;

      // These are all the reactions
      std::map<std::string, ito_reaction>   m_reactions; // SHould it be mutable...????
      mutable std::map<std::string, photo_reaction> m_photo_reactions;

      /*!
	@brief Draw a valid random position somewhere in a cell
      */
      inline
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
      inline
      RealVect random_position(const RealVect a_lo,
			       const RealVect a_hi,
			       const RealVect a_bndryCentroid,
			       const RealVect a_normal) const;

      /*!
	@brief Return a random position in the cube (a_lo, a_hi);
      */
      inline
      RealVect random_position(const RealVect a_lo, const RealVect a_hi) const;

      /*!
	@brief Draw a random direction in space
      */
      inline
      RealVect random_direction() const;

      /*!
	@brief Check if advance is valid. This returns false if there is a negative number of particles in the cell
      */
      inline
      bool valid_advance(const Vector<List<ito_particle>* >& a_particles,
			 const Vector<int>&                  a_newParticles) const;

      /*!
	@brief Reconcile the number of particles
      */
      inline
      void reconcile_particles(Vector<List<ito_particle>* >& a_particles,
			       const Vector<int>&            a_newNumParticles,
			       const Vector<int>&            a_oldNumParticles,
			       const RealVect                a_pos,
			       const RealVect                a_lo,
			       const RealVect                a_hi,
			       const RealVect                a_bndryCentroid,
			       const RealVect                a_bndryNormal,
			       const Real                    a_dx,
			       const Real                    a_kappa) const;
      /*!
	@brief Reconcile number of particles
      */
      inline
      void reconcile_particles(List<ito_particle>& a_particles,
			       const int           a_numNewParticles,
			       const int           a_oldNewParticles,
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
      inline
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
      inline
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
	@brief Reconcile the number of particles
      */
      inline
      void reconcile_photons(Vector<List<photon>* >& a_newPhotons,
			     const Vector<int>&      a_numNewPhotons,
			     const RealVect          a_pos,
			     const RealVect          a_lo,
			     const RealVect          a_hi,
			     const RealVect          a_bndryCentroid,
			     const RealVect          a_bndryNormal,
			     const Real              a_dx,
			     const Real              a_kappa) const;
      /*!
	@brief Reconcile number of particles
      */
      inline
      void reconcile_photons(List<photon>&  a_newPhotons,
			     const int      a_numNewPhotons,
			     const RealVect a_pos,
			     const RealVect a_lo,
			     const RealVect a_hi,
			     const RealVect a_bndryCentroid,
			     const RealVect a_normal,
			     const Real     a_dx,
			     const Real     a_kappa) const;

      /*!
	@brief Add photons
      */
      inline
      void reconcile_photons(List<photon>&      a_photons,
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
	@brief Reconcile photoionization mechanism
      */
      inline
      void reconcile_photoionization(Vector<List<ito_particle>* >& a_particles,
				     const Vector<List<photon>* >& a_photons) const;

      /*!
	@brief Compute particle weights for superparticles
      */
      inline
      void compute_particle_weights(int& a_weight, int& a_num, int& remainder, const int a_numNewParticles) const;

      /*!
	@brief Poisson reaction
      */
      inline
      int poisson_reaction(const Real a_propensity, const Real a_dt) const;

      /*!
	@brief Draw Gaussian distributed particles
      */
      inline
      void draw_gaussian_particles(List<ito_particle>& a_particles1,
				   List<ito_particle>& a_particles2,
				   const int           a_numParticles,
				   const RealVect      a_blobCenter,
				   const Real          a_blobRadius,
				   const Real          a_weights);

      /*!
	@brief Get the particle count
      */
      inline
      Vector<int> get_particle_count(const Vector<List<ito_particle>* >& a_particles) const;

      /*!
	@brief Get the particle count
      */
      inline
      int get_particle_count(const List<ito_particle>& a_particles) const;

      /*!
	@brief Performing a tau-leaping step. Mostly while we wait for hybrid...
      */
      inline
      void tau_leap(Vector<int>& a_particles, Vector<int>& a_newPhotons, const Real a_dt) const;

      /*!
	@brief Perform hybrid advance
      */
      inline
      void hybrid_advance(Vector<int>& a_particles, Vector<int>& a_newPhotons, const Real a_dt);

      /*!
	@brief Get critical and non-critical reactions
      */
      inline
      ReactionTypes get_critical_reactions(const Vector<int>& a_particles) const;

      /*!
	@brief Get the time to the next critical reaction. 
      */
      inline
      Real get_critical_step(const Vector<int>&  a_particles,
			     const ReactionPtrs& a_critical_reactions) const;

      /*!
	@brief Get the time to the next critical reaction. 
      */
      inline
      Real get_non_critical_step(const Vector<int>&  a_particles,
				 const ReactionPtrs& a_non_critical_reactions) const;

      /*!
	@brief Compute total propensity
      */
      inline
      Real total_propensity(const Vector<int>& a_particles) const;


      /*!
	@brief Critical propensity function
      */
      inline
      Real critical_propensity(const Vector<int>& a_particles, const ReactionPtrs& a_critical_reactions) const;

      /*!
	@brief Compute critical propensities
      */
      inline
      std::deque<Real> critical_propensities(const Vector<int>& a_particles, const ReactionPtrs& a_critical_reactions) const;


      /*!
	@brief Performing a single SSA step using a subset of reactions. 
      */
      inline
      void ssa_step(Vector<int>&             a_particles,
		    Vector<int>&             a_newPhotons,
		    const ReactionPtrs&      a_critical_reactions,
		    const std::deque<Real>&  a_propensities,
		    const Real               a_A) const;



      /*!
	@brief Performing a single SSA step using all reactions. 
      */
      inline
      void ssa_step(Vector<int>&             a_particles,
		    Vector<int>&             a_newPhotons,
		    const std::deque<Real>&  a_propensities,
		    const Real               a_A) const;


      /*!
	@brief Perform tau-leaping with a subsect of reactions
      */
      inline void tau_step(Vector<int>&             a_particles,
			   Vector<int>&             a_newPhotons,
			   const ReactionPtrs&      a_non_critical_reactions,
			   const std::deque<Real>&  a_propensities,
			   const Real               a_dt) const;


      /*!
	@brief Perform tau-leaping with a subsect of reactions
      */
      inline void tau_step(Vector<int>&             a_particles,
			   Vector<int>&             a_newPhotons,
			   const std::deque<Real>&  a_propensities,
			   const Real               a_dt) const;
    };
  }
}

#include "ito_plasma_physicsI.H"

#endif
