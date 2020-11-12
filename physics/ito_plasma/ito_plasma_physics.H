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
#include <fstream>
#include <sstream>

#include "ito_species.H"
#include "rte_species.H"
#include "photon.H"
#include "ito_particle.H"
#include "ito_reaction.H"
#include "photo_reaction.H"
#include "lookup_table.H"

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

      enum class coupling {LFA, LEA };
      
      /*!
	@brief Ctor. Does nothing.
      */
      ito_plasma_physics();

      /*!
	@brief Dtor. Does nothing. 
      */
      ~ito_plasma_physics();

      /*!
	@brief Compute a time step
      */
      virtual Real compute_dt(const RealVect a_E, const RealVect a_pos, const Vector<Real> a_densities) const = 0;

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
	@brief Get the coupling, LFA/LEA
      */
      coupling get_coupling() const;

      /*!
	@brief Set initial surface charge. Default is 0, override if you want.
      */
      virtual Real initial_sigma(const Real a_time, const RealVect a_pos) const;

      /*!
	@brief Compute the ito mobilities in the lfa approximation
      */
      virtual Vector<Real> compute_ito_mobilities_lfa(const Real a_time, const RealVect a_pos, const RealVect a_E) const = 0;

      /*!
	@brief Compute diffusion coefficients
      */
      virtual Vector<Real> compute_ito_diffusion_lfa(const Real         a_time,
						     const RealVect     a_pos,
						     const RealVect     a_E,
						     const Vector<Real> a_cdr_densities) const = 0;

      /*!
	@brief Update reaction rates for LFA updates. Since the rates don't change during the advance, they can
	be set to a constant. 
      */
      virtual void update_reaction_rates_lfa(const RealVect a_E, const Real a_dx, const Real a_kappa) const = 0;

      /*!
	@brief Update reaction rates for LFA updates. Since the rates don't change during the advance, they can
	be set to a constant. 
      */
      virtual void update_reaction_rates_lea(const RealVect a_E, const Vector<Real> a_mean_energies, const Real a_dx, const Real a_kappa) const = 0;

      /*!
	@brief Update reaction rates
      */
      inline
      void update_reaction_rates(const RealVect a_E, const Vector<Real> a_mean_energies, const Real a_dx, const Real a_kappa) const;


      /*!
	@brief Reaction network advance in the LEA.
      */
      inline
      void advance_reaction_network(Vector<List<ito_particle>* >& a_particles,
				    Vector<List<photon>* >&       a_photons,
				    Vector<List<photon>* >&       a_newPhotons,
				    const Vector<Real>&           a_sources, 
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

      enum class algorithm { ssa, tau, hybrid };



      int m_seed;
      int m_table_entries;
      
      algorithm m_algorithm;
      coupling  m_coupling;
      
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
      int m_NSSA;

      Real m_SSAlim;
      Real m_eps;

      // These are all the reactions
      std::map<std::string, ito_reaction> m_reactions; 
      mutable std::map<std::string, photo_reaction> m_photo_reactions;

      // Lookup tables if you need them. 
      std::map<std::string, lookup_table> m_tables;

      /*!
	@brief Add a table
      */
      void add_table(const std::string a_table_name, const std::string a_file);

      /*!
	@brief Read a file and put it in a lookup table
      */
      void read_file(lookup_table& a_table, const std::string a_file);

      /*!
	@brief Draw a valid random position somewhere in a cell
      */
      inline
      RealVect random_position(const RealVect a_cellPos,
			       const RealVect a_centroidPos,
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
	@brief Returns true if all particle numbers are >= 0
      */
      inline
      bool valid_state(const Vector<long long>& a_particles) const;

      /*!
	@brief Check charge conservation in a cell
      */
      inline
      void check_charge_conservation(const Vector<long long>& a_oldNumParticles,
				     const Vector<long long>& a_newNumParticles) const;

      /*!
	@brief Compute total charge in a cell
      */
      inline
      long long compute_charge(const Vector<List<ito_particle>* >& a_particles) const;

      /*!
	@brief Reconcile the number of particles
      */
      inline
      void reconcile_particles(Vector<List<ito_particle>* >& a_particles,
			       const Vector<long long>&      a_newNumParticles,
			       const Vector<long long>&      a_oldNumParticles,
			       const RealVect                a_cellPos,
			       const RealVect                a_centroidPos,
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
			       const long long     a_numNewParticles,
			       const long long     a_oldNewParticles,
			       const RealVect      a_cellPos,
			       const RealVect      a_centroidPos,
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
			 const long long     a_numParticlesToAdd,
			 const RealVect      a_cellPos,
			 const RealVect      a_centroidPos,
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
      void remove_particles(List<ito_particle>& a_particles, const long long a_numParticlesToRemove) const;

      /*!
	@brief Set mean particle energy
      */
      inline
      void set_mean_particle_energy(Vector<List<ito_particle>* >& a_particles, const Vector<Real> a_mean_energy) const;

      /*!
	@brief Set mean particle energy
      */
      inline
      void set_mean_particle_energy(List<ito_particle>& a_particles, const Real a_mean_energy) const;

      /*!
	@brief Check that all weights are > 0
      */
      inline
      void sanity_check(const List<ito_particle>& a_lit) const;

      /*!
	@brief Check that all weights are > 0
      */
      inline
      long long get_total_weight(const List<ito_particle>& a_lit) const;
      
      /*!
	@brief Compute particle weights
      */
      inline
      void compute_weights(ListIterator<ito_particle>& a_lit,
			   Vector<long long>&          a_particleWeights,
			   long long&                  a_minWeight) const;

      /*!
	@brief Trim particle weights down to a_minWeight
      */
      inline
      void trim_weights(long long&         a_numRemoved,
			Vector<long long>& a_particleWeights,
			const long long    a_numParticlesToRemove,
			const long long    a_minWeight) const;

      /*!
	@brief Remove particle weight down to a_numParticlesToRemove
      */
      inline
      void remove_weights(Vector<long long>& a_particleWeights, const long long a_numParticlesToRemove) const;

      /*!
	@brief Remove particle weight down to a_numParticlesToRemove
      */
      inline
      void reconcile_weights(List<ito_particle>& a_particles, const Vector<long long>& a_particleWeights) const;

      /*!
	@brief Reconcile the number of particles
      */
      inline
      void reconcile_photons(Vector<List<photon>* >&  a_newPhotons,
			     const Vector<long long>& a_numNewPhotons,
			     const RealVect           a_cellPos,
			     const RealVect           a_centroidPos,
			     const RealVect           a_lo,
			     const RealVect           a_hi,
			     const RealVect           a_bndryCentroid,
			     const RealVect           a_bndryNormal,
			     const Real               a_dx,
			     const Real               a_kappa) const;
      /*!
	@brief Reconcile number of photons
      */
      inline
      void reconcile_photons(List<photon>&      a_newPhotons,
			     const rte_species& a_species,
			     const long long    a_numNewPhotons,
			     const RealVect     a_cellPos,
			     const RealVect     a_centroidPos,
			     const RealVect     a_lo,
			     const RealVect     a_hi,
			     const RealVect     a_bndryCentroid,
			     const RealVect     a_normal,
			     const Real         a_dx,
			     const Real         a_kappa) const;

      /*!
	@brief Reconcile photoionization mechanism
      */
      inline
      void reconcile_photoionization(Vector<List<ito_particle>* >& a_particles,
				     Vector<Real>&                 a_mean_energies,
				     const Vector<long long>&      a_numParticles,
				     const Vector<List<photon>* >& a_photons) const;

      /*!
	@brief Compute particle weights for superparticles
      */
      inline
      void compute_particle_weights(long long&      a_weight,
				    long long&      a_num,
				    long long&      a_remainder,
				    const long long a_numNewParticles) const;

      /*!
	@brief Poisson reaction
      */
      inline
      long long poisson_reaction(const Real a_propensity, const Real a_dt) const;

      /*!
	@brief Compute particles to draw per rank
      */
      inline
      int compute_particles_per_rank(const int a_numParticles) const;

      /*!
	@brief Draw Gaussian distributed particles
      */
      inline
      void draw_gaussian_particles(List<ito_particle>& a_particles1,
				   List<ito_particle>& a_particles2,
				   const int           a_numParticles,
				   const RealVect      a_blobCenter,
				   const Real          a_blobRadius,
				   const Real          a_weight,
				   const Real          a_energy1,
				   const Real          a_energy2);

      /*!
	@brief Draw sphere particles
      */
      inline
      void draw_sphere_particles(List<ito_particle>& a_particles1,
				 List<ito_particle>& a_particles2,
				 const int           a_numParticles,
				 const RealVect      a_blobCenter,
				 const Real          a_blobRadius,
				 const Real          a_weight,
				 const Real          a_energy1,
				 const Real          a_energy2);

      /*!
	@brief Draw particles in a box
      */
      inline
      void draw_box_particles(List<ito_particle>& a_particles1,
			      List<ito_particle>& a_particles2,
			      const int           a_numParticles,
			      const RealVect      a_boxLo,
			      const RealVect      a_boxHi,
			      const Real          a_weight,
			      const Real          a_energy1,
			      const Real          a_energy2);

      /*!
	@brief Get the particle count
      */
      inline
      Vector<long long> get_particle_count(const Vector<List<ito_particle>* >& a_particles) const;
      
      /*!
	@brief Get the particle count
      */
      inline
      long long get_particle_count(const List<ito_particle>& a_particles) const;

      /*!
	@brief Get the particle energies
      */
      inline
      Vector<Real> get_particle_energy(const Vector<List<ito_particle>* >& a_particles) const;

      /*!
	@brief Get the particle energies
      */
      inline
      Real get_particle_energy(const List<ito_particle>& a_particles) const;

      /*!
	@brief Get the mean particle energies
      */
      inline
      Vector<Real> get_particle_mean_energy(const Vector<Real>& a_energies, const Vector<long long>& a_particle_count) const;

      /*!
	@brief Get the particle energies
      */
      inline
      Real get_particle_mean_energy(const Real& a_energy, const long long& a_particle_count) const;


      /*!
	@brief Get critical and non-critical reactions
      */
      inline
      ReactionTypes get_critical_reactions(const Vector<long long>& a_particles) const;

      /*!
	@brief Get the time to the next critical reaction. 
      */
      inline
      Real get_critical_step(const Vector<long long>& a_particles,
			     const ReactionPtrs&               a_critical_reactions) const;

      /*!
	@brief Get the time to the next critical reaction. 
      */
      inline
      Real get_non_critical_step(const Vector<long long>& a_particles,
				 const ReactionPtrs&      a_non_critical_reactions) const;

      /*!
	@brief Compute the total propensity for all reactions. 
      */
      inline
      Real propensity(const Vector<long long>& a_particles) const;

      /*!
	@brief Compute total propensity for a subset of reactions. 
      */
      inline
      Real propensity(const Vector<long long>& a_particles, const ReactionPtrs& a_critical_reactions) const;

      /*!
	@brief Computes all propensities
      */
      inline
      std::deque<Real> propensities(const Vector<long long>& a_particles) const;

      /*!
	@brief Compute critical propensities
      */
      inline
      std::deque<Real> propensities(const Vector<long long>& a_particles, const ReactionPtrs& a_critical_reactions) const;

      /*!
	@brief Do a tau-step with all reactions. You must update the reaction rates before entering this routine! This calls the other version. 
      */
      inline
      void tau_step(Vector<long long>&  a_particles,
		    Vector<long long>&  a_newPhotons,
		    Vector<Real>&       a_mean_energies,
		    const Vector<Real>& a_sources,
		    const Real          a_dx,
		    const Real          a_kappa,
		    const Real          a_dt) const;

      /*!
	@brief Do a tau-step with a subset of reactions. You must update the reaction rates before entering this routine!
      */
      inline
      void tau_step(Vector<long long>&  a_particles,
		    Vector<long long>&  a_newPhotons,
		    Vector<Real>&       a_mean_energies,
		    const Vector<Real>& a_sources,
		    const ReactionPtrs& a_reactions,
		    const Real          a_dx,
		    const Real          a_kappa,
		    const Real          a_dt) const;

      /*!
	@brief Do a tau-step with a subset of reactions. You must update the reaction rates before entering this routine! Propensities are not
	       recomputed with this version. 
      */
      inline
      void tau_step(Vector<long long>&      a_particles,
		    Vector<long long>&      a_newPhotons,
		    Vector<Real>&           a_mean_energies,
		    const Vector<Real>&     a_sources,
		    const ReactionPtrs&     a_reactions,
		    const std::deque<Real>& a_propensities,
		    const Real              a_dx,
		    const Real              a_kappa,
		    const Real              a_dt) const;

      /*!
	@brief New tau leaping advance. If using the LFA then a_mean_energies and a_sources are dummy arguments. 
      */
      inline
      void advance_tau(Vector<long long>&  a_particles,
		       Vector<long long>&  a_newPhotons,
		       Vector<Real>&       a_mean_energies,
		       const Vector<Real>& a_sources,
		       const Real          a_dt,
		       const RealVect      a_E,
		       const Real          a_dx,
		       const Real          a_kappa) const;

      /*!
	@brief Advance particles. If usign the LFA then a_mean_energies and a_sources are dummy arugments. 
      */
      inline
      void advance_particles(Vector<long long>&  a_particles,
			     Vector<long long>&  a_newPhotons,
			     Vector<Real>&       a_mean_energies,
			     const Vector<Real>& a_sources,
			     const Real          a_dt,
			     const RealVect      a_E,
			     const Real          a_dx,
			     const Real          a_kappa) const;

      /*!
	@brief Do an SSA-step with a subset of reactions. You must update the reaction rates before entering this routine!
      */
      inline
      void ssa_step(Vector<long long>&      a_particles,
		    Vector<long long>&      a_newPhotons,
		    Vector<Real>&           a_mean_energies,
		    const Vector<Real>&     a_sources,
		    const std::deque<Real>& a_propensities,
		    const ReactionPtrs&     a_reactions,
		    const Real              a_A,
		    const Real              a_dx,
		    const Real              a_kappa,
		    const Real              a_dt) const;

      /*!
	@brief This is a version that steps with all reactions, but where the propensities have been computed in advance. 
      */
      inline
      void ssa_step(Vector<long long>&      a_particles,
		    Vector<long long>&      a_newPhotons,
		    Vector<Real>&           a_mean_energies,
		    const Vector<Real>&     a_sources,
		    const std::deque<Real>& a_propensities,
		    const Real              a_A,
		    const Real              a_dx,
		    const Real              a_kappa,
		    const Real              a_dt) const;

      /*!
	@brief Do an SSA-step with all reactions. You must update the reaction rates before entering this routine! This computes all propensities and
	       then calls the other version. 
      */
      inline
      void ssa_step(Vector<long long>&  a_particles,
		    Vector<long long>&  a_newPhotons,
		    Vector<Real>&       a_mean_energies,
		    const Vector<Real>& a_sources,
		    const Real          a_dx,
		    const Real          a_kappa,
		    const Real          a_dt) const;

      /*!
	@brief Increase/decrease the mean energies with a volumetric source term
      */
      inline
      void energy_step(Vector<Real>&            a_mean_energies,
		       const Vector<long long>& a_particles,
		       const Vector<Real>&      a_sources,
		       const Real               a_dx,
		       const Real               a_kappa,
		       const Real               a_dt) const;

      /*!
	@brief New SSA advance. If using the LFA then a_mean_energies and a_sources are dummy arguments. 
      */
      inline
      void advance_ssa(Vector<long long>&  a_particles,
		       Vector<long long>&  a_newPhotons,
		       Vector<Real>&       a_mean_energies,
		       const Vector<Real>& a_sources,
		       const Real          a_dt,
		       const RealVect      a_E,
		       const Real          a_dx,
		       const Real          a_kappa) const;

      /*!
	@brief New hybrid advance. If using the LFA then a_mean_energies and a_sources are dummy arguments. 
      */
      inline
      void advance_hybrid(Vector<long long>&  a_particles,
			  Vector<long long>&  a_newPhotons,
			  Vector<Real>&       a_mean_energies,
			  const Vector<Real>& a_sources,
			  const Real          a_dt,
			  const RealVect      a_E,
			  const Real          a_dx,
			  const Real          a_kappa) const;

    private:
      
      // trim from start (in place)
      inline
      void ltrim(std::string &s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
      }
      
      // trim from end (in place)
      inline
      void rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
      }
      
      // trim from both ends (in place)
      inline
      void trim(std::string &s) {
	ltrim(s);
	rtrim(s);
      }
    };
  }
}

#include "ito_plasma_physicsI.H"

#endif
