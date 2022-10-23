/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaPhysics.H
  @brief  Main file for describing Ito-based plasma physics
  @author Robert Marskar
*/

#ifndef CD_ItoPlasmaPhysics_H
#define CD_ItoPlasmaPhysics_H

// Std includes
#include <chrono>
#include <random>
#include <memory>
#include <deque>
#include <fstream>
#include <sstream>

// Chombo includes
#include <RealVect.H>
#include <RefCountedPtr.H>
#include <LoHiSide.H>
#include <List.H>

// Our includes
#include <CD_ItoSpecies.H>
#include <CD_RtSpecies.H>
#include <CD_Photon.H>
#include <CD_ItoParticle.H>
#include <CD_ItoPlasmaReaction.H>
#include <CD_ItoPlasmaPhotoReaction.H>
#include <CD_LookupTable1D.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoPlasma {

    // Ease up on the typing.
    using ReactionPtrs  = std::deque<std::shared_ptr<ItoPlasmaReaction>>;
    using ReactionTypes = std::pair<ReactionPtrs, ReactionPtrs>;

    /*!
      @brief Abstract clas for specifying plasma kinetics using Kinetic Monte Carlo. 
    */
    class ItoPlasmaPhysics
    {
    public:
      enum class coupling
      {
        LFA,
        LEA
      };

      /*!
	@brief Constructor. Does nothing.
      */
      ItoPlasmaPhysics();

      /*!
	@brief Destructor. Does nothing. 
      */
      virtual ~ItoPlasmaPhysics();

      /*!
	@brief Compute a time step
      */
      virtual Real
      computeDt(const RealVect a_E, const RealVect a_pos, const Vector<Real> a_densities) const = 0;

      /*!
	@brief Compute Townsend ionization length
      */
      virtual Real
      computeAlpha(const RealVect a_E) const = 0;

      /*!
	@brief Get all species
      */
      const Vector<RefCountedPtr<ItoSpecies>>&
      getItoSpecies() const;

      /*!
	@brief Get all Photons
      */
      const Vector<RefCountedPtr<RtSpecies>>&
      getRtSpecies() const;

      /*!
	@brief Return number of ion equations
      */
      int
      getNumPlasmaSpecies() const;

      /*!
	@brief Return number of RTE equations
      */
      int
      getNumPhotonSpecies() const;

      /*!
	@brief Get the coupling, LFA/LEA
      */
      coupling
      getCoupling() const;

      /*!
	@brief Set initial surface charge. Default is 0, override if you want.
      */
      virtual Real
      initialSigma(const Real a_time, const RealVect a_pos) const;

      /*!
	@brief Compute the ito mobilities in the lfa approximation
      */
      virtual Vector<Real>
      computeItoMobilitiesLFA(const Real a_time, const RealVect a_pos, const RealVect a_E) const = 0;

      /*!
	@brief Compute diffusion coefficients
      */
      virtual Vector<Real>
      computeItoDiffusionLFA(const Real         a_time,
                             const RealVect     a_pos,
                             const RealVect     a_E,
                             const Vector<Real> a_cdr_densities) const = 0;

      /*!
	@brief Update reaction rates for LFA updates. Since the rates don't change during the advance, they can
	be set to a constant. 
      */
      virtual void
      updateReactionRatesLFA(const RealVect a_E, const Real a_dx, const Real a_kappa) const = 0;

      /*!
	@brief Update reaction rates for LFA updates. Since the rates don't change during the advance, they can
	be set to a constant. 
      */
      virtual void
      updateReactionRatesLEA(const RealVect     a_E,
                             const Vector<Real> a_mean_energies,
                             const Real         a_dx,
                             const Real         a_kappa) const = 0;

      /*!
	@brief Update reaction rates
      */
      inline void
      updateReactionRates(const RealVect     a_E,
                          const Vector<Real> a_mean_energies,
                          const Real         a_dx,
                          const Real         a_kappa) const;

      /*!
	@brief Advance particles. If usign the LFA then a_mean_energies and a_sources are dummy arugments. 
      */
      inline void
      advanceParticles(Vector<long long>&  a_particles,
                       Vector<long long>&  a_newPhotons,
                       Vector<Real>&       a_mean_energies,
                       const Vector<Real>& a_sources,
                       const Real          a_dt,
                       const RealVect      a_E,
                       const Real          a_dx,
                       const Real          a_kappa) const;

      /*!
	@brief Reaction network advance in the LEA.
      */
      inline void
      advanceReactionNetwork(Vector<List<ItoParticle>*>& a_particles,
                             Vector<List<Photon>*>&      a_Photons,
                             Vector<List<Photon>*>&      a_newPhotons,
                             const Vector<Real>&         a_sources,
                             const RealVect              a_E,
                             const RealVect              a_cellPos,
                             const RealVect              a_centroidPos,
                             const RealVect              a_bndryCentroid,
                             const RealVect              a_bndryNormal,
                             const RealVect              a_lo,
                             const RealVect              a_hi,
                             const Real                  a_dx,
                             const Real                  a_kappa,
                             const Real                  a_dt);

      /*!
	@brief Reconcile the number of particles
      */
      inline void
      reconcilePhotons(Vector<List<Photon>*>&   a_newPhotons,
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
	@brief Reconcile the number of particles
      */
      inline void
      reconcileParticles(Vector<List<ItoParticle>*>& a_particles,
                         const Vector<long long>&    a_newNumParticles,
                         const Vector<long long>&    a_oldNumParticles,
                         const RealVect              a_cellPos,
                         const RealVect              a_centroidPos,
                         const RealVect              a_lo,
                         const RealVect              a_hi,
                         const RealVect              a_bndryCentroid,
                         const RealVect              a_bndryNormal,
                         const Real                  a_dx,
                         const Real                  a_kappa) const;

      /*!
	@brief Reconcile number of Photons
      */
      inline void
      reconcilePhotons(List<Photon>&    a_newPhotons,
                       const RtSpecies& a_species,
                       const long long  a_numNewPhotons,
                       const RealVect   a_cellPos,
                       const RealVect   a_centroidPos,
                       const RealVect   a_lo,
                       const RealVect   a_hi,
                       const RealVect   a_bndryCentroid,
                       const RealVect   a_normal,
                       const Real       a_dx,
                       const Real       a_kappa) const;

      /*!
	@brief Reconcile photoionization mechanism
      */
      inline void
      reconcilePhotoionization(Vector<List<ItoParticle>*>&  a_particles,
                               Vector<Real>&                a_mean_energies,
                               const Vector<long long>&     a_numParticles,
                               const Vector<List<Photon>*>& a_Photons) const;

      /*!
	@brief Set mean particle energy
      */
      inline void
      setMeanParticleEnergy(Vector<List<ItoParticle>*>& a_particles, const Vector<Real> a_mean_energy) const;

    protected:
      enum class algorithm
      {
        ssa,
        tau,
        hybrid
      };

      int m_seed;
      int m_table_entries;

      algorithm m_algorithm;
      coupling  m_coupling;

      mutable std::mt19937_64                      m_rng;
      mutable std::uniform_real_distribution<Real> m_udist01;
      mutable std::uniform_real_distribution<Real> m_udist11;

      Vector<RefCountedPtr<ItoSpecies>> m_ItoSpecies;
      Vector<RefCountedPtr<RtSpecies>>  m_rtSpecies;

      int m_num_ItoSpecies;
      int m_numRtSpecies;
      int m_ppc;
      int m_fieldSolver_switch;
      int m_init_weights;
      int m_Ncrit;
      int m_NSSA;

      Real m_SSAlim;
      Real m_eps;

      // These are all the reactions
      std::map<std::string, ItoPlasmaReaction>              m_reactions;
      mutable std::map<std::string, ItoPlasmaPhotoReaction> m_photoReactions;

      // Lookup tables if you need them.
      std::map<std::string, LookupTable1D<2>> m_tables;

      /*!
	@brief Add a table
      */
      void
      addTable(const std::string a_table_name, const std::string a_file);

      /*!
	@brief Read a file and put it in a lookup table
      */
      void
      readFile(LookupTable1D<2>& a_table, const std::string a_file);

      /*!
	@brief Draw a valid random position somewhere in a cell
      */
      inline RealVect
      randomPosition(const RealVect a_cellPos,
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
      inline RealVect
      randomPosition(const RealVect a_lo,
                     const RealVect a_hi,
                     const RealVect a_bndryCentroid,
                     const RealVect a_normal) const;

      /*!
	@brief Return a random position in the cube (a_lo, a_hi);
      */
      inline RealVect
      randomPosition(const RealVect a_lo, const RealVect a_hi) const;

      /*!
	@brief Draw a random direction in space
      */
      inline RealVect
      randomDirection() const;

      /*!
	@brief Returns true if all particle numbers are >= 0
      */
      inline bool
      isStateValid(const Vector<long long>& a_particles) const;

      /*!
	@brief Check charge conservation in a cell
      */
      inline void
      checkChargeConservation(const Vector<long long>& a_oldNumParticles,
                              const Vector<long long>& a_newNumParticles) const;

      /*!
	@brief Compute total charge in a cell
      */
      inline long long
      computeCharge(const Vector<List<ItoParticle>*>& a_particles) const;

      /*!
	@brief Reconcile number of particles
      */
      inline void
      reconcileParticles(List<ItoParticle>& a_particles,
                         const long long    a_numNewParticles,
                         const long long    a_oldNewParticles,
                         const RealVect     a_cellPos,
                         const RealVect     a_centroidPos,
                         const RealVect     a_lo,
                         const RealVect     a_hi,
                         const RealVect     a_bndryCentroid,
                         const RealVect     a_normal,
                         const Real         a_dx,
                         const Real         a_kappa) const;

      /*!
	@brief Add particles
      */
      inline void
      addParticles(List<ItoParticle>& a_particles,
                   const long long    a_numParticlesToAdd,
                   const RealVect     a_cellPos,
                   const RealVect     a_centroidPos,
                   const RealVect     a_lo,
                   const RealVect     a_hi,
                   const RealVect     a_bndryCentroid,
                   const RealVect     a_bndryNormal,
                   const Real         a_dx,
                   const Real         a_kappa) const;

      /*!
	@brief Remove particles
      */
      inline void
      remove_particles(List<ItoParticle>& a_particles, const long long a_numParticlesToRemove) const;

      /*!
	@brief Set mean particle energy
      */
      inline void
      setMeanParticleEnergy(List<ItoParticle>& a_particles, const Real a_mean_energy) const;

      /*!
	@brief Check that all weights are > 0
      */
      inline void
      sanityCheck(const List<ItoParticle>& a_lit) const;

      /*!
	@brief Check that all weights are > 0
      */
      inline long long
      getTotalWeight(const List<ItoParticle>& a_lit) const;

      /*!
	@brief Compute particle weights
      */
      inline void
      computeWeights(ListIterator<ItoParticle>& a_lit,
                     Vector<long long>&         a_particleWeights,
                     long long&                 a_minWeight) const;

      /*!
	@brief Trim particle weights down to a_minWeight
      */
      inline void
      trimWeights(long long&         a_numRemoved,
                  Vector<long long>& a_particleWeights,
                  const long long    a_numParticlesToRemove,
                  const long long    a_minWeight) const;

      /*!
	@brief Remove particle weight down to a_numParticlesToRemove
      */
      inline void
      removeWeights(Vector<long long>& a_particleWeights, const long long a_numParticlesToRemove) const;

      /*!
	@brief Remove particle weight down to a_numParticlesToRemove
      */
      inline void
      reconcileWeights(List<ItoParticle>& a_particles, const Vector<long long>& a_particleWeights) const;

      /*!
	@brief Compute particle weights for superparticles
      */
      inline void
      computeParticleWeights(long long&      a_weight,
                             long long&      a_num,
                             long long&      a_remainder,
                             const long long a_numNewParticles) const;

      /*!
	@brief Poisson reaction
      */
      inline long long
      poissonReaction(const Real a_propensity, const Real a_dt) const;

      /*!
	@brief Compute particles to draw per rank
      */
      inline int
      computeParticlesPerRank(const int a_numParticles) const;

      /*!
	@brief Draw Gaussian distributed particles
      */
      inline void
      drawGaussianParticles(List<ItoParticle>& a_particles1,
                            List<ItoParticle>& a_particles2,
                            const int          a_numParticles,
                            const RealVect     a_blobCenter,
                            const Real         a_blobRadius,
                            const Real         a_weight,
                            const Real         a_energy1,
                            const Real         a_energy2);

      /*!
	@brief Draw sphere particles
      */
      inline void
      drawSphereParticles(List<ItoParticle>& a_particles1,
                          List<ItoParticle>& a_particles2,
                          const int          a_numParticles,
                          const RealVect     a_blobCenter,
                          const Real         a_blobRadius,
                          const Real         a_weight,
                          const Real         a_energy1,
                          const Real         a_energy2);

      /*!
	@brief Draw particles in a box
      */
      inline void
      drawBoxParticles(List<ItoParticle>& a_particles1,
                       List<ItoParticle>& a_particles2,
                       const int          a_numParticles,
                       const RealVect     a_boxLo,
                       const RealVect     a_boxHi,
                       const Real         a_weight,
                       const Real         a_energy1,
                       const Real         a_energy2);

      /*!
	@brief Get the particle count
      */
      inline Vector<long long>
      getParticleCount(const Vector<List<ItoParticle>*>& a_particles) const;

      /*!
	@brief Get the particle count
      */
      inline long long
      getParticleCount(const List<ItoParticle>& a_particles) const;

      /*!
	@brief Get the particle energies
      */
      inline Vector<Real>
      getParticleEnergy(const Vector<List<ItoParticle>*>& a_particles) const;

      /*!
	@brief Get the particle energies
      */
      inline Real
      getParticleEnergy(const List<ItoParticle>& a_particles) const;

      /*!
	@brief Get the mean particle energies
      */
      inline Vector<Real>
      getParticleMeanEnergy(const Vector<Real>& a_energies, const Vector<long long>& a_particle_count) const;

      /*!
	@brief Get the particle energies
      */
      inline Real
      getParticleMeanEnergy(const Real& a_energy, const long long& a_particle_count) const;

      /*!
	@brief Get critical and non-critical reactions
      */
      inline ReactionTypes
      getCriticalReactions(const Vector<long long>& a_particles) const;

      /*!
	@brief Get the time to the next critical reaction. 
      */
      inline Real
      getCriticalTimeStep(const Vector<long long>& a_particles, const ReactionPtrs& a_critical_reactions) const;

      /*!
	@brief Get the time to the next critical reaction. 
      */
      inline Real
      getNonCriticalTimeStep(const Vector<long long>& a_particles, const ReactionPtrs& a_non_critical_reactions) const;

      /*!
	@brief Compute the total propensity for all reactions. 
      */
      inline Real
      propensity(const Vector<long long>& a_particles) const;

      /*!
	@brief Compute total propensity for a subset of reactions. 
      */
      inline Real
      propensity(const Vector<long long>& a_particles, const ReactionPtrs& a_critical_reactions) const;

      /*!
	@brief Computes all propensities
      */
      inline std::deque<Real>
      propensities(const Vector<long long>& a_particles) const;

      /*!
	@brief Compute critical propensities
      */
      inline std::deque<Real>
      propensities(const Vector<long long>& a_particles, const ReactionPtrs& a_critical_reactions) const;

      /*!
	@brief Do a tau-step with all reactions. You must update the reaction rates before entering this routine! This calls the other version. 
      */
      inline void
      stepTau(Vector<long long>&  a_particles,
              Vector<long long>&  a_newPhotons,
              Vector<Real>&       a_mean_energies,
              const Vector<Real>& a_sources,
              const Real          a_dx,
              const Real          a_kappa,
              const Real          a_dt) const;

      /*!
	@brief Do a tau-step with a subset of reactions. You must update the reaction rates before entering this routine!
      */
      inline void
      stepTau(Vector<long long>&  a_particles,
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
      inline void
      stepTau(Vector<long long>&      a_particles,
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
      inline void
      advanceTau(Vector<long long>&  a_particles,
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
      inline void
      stepSSA(Vector<long long>&      a_particles,
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
      inline void
      stepSSA(Vector<long long>&      a_particles,
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
      inline void
      stepSSA(Vector<long long>&  a_particles,
              Vector<long long>&  a_newPhotons,
              Vector<Real>&       a_mean_energies,
              const Vector<Real>& a_sources,
              const Real          a_dx,
              const Real          a_kappa,
              const Real          a_dt) const;

      /*!
	@brief Increase/decrease the mean energies with a volumetric source term
      */
      inline void
      stepEnergy(Vector<Real>&            a_mean_energies,
                 const Vector<long long>& a_particles,
                 const Vector<Real>&      a_sources,
                 const Real               a_dx,
                 const Real               a_kappa,
                 const Real               a_dt) const;

      /*!
	@brief New SSA advance. If using the LFA then a_mean_energies and a_sources are dummy arguments. 
      */
      inline void
      advanceSSA(Vector<long long>&  a_particles,
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
      inline void
      advanceHybrid(Vector<long long>&  a_particles,
                    Vector<long long>&  a_newPhotons,
                    Vector<Real>&       a_mean_energies,
                    const Vector<Real>& a_sources,
                    const Real          a_dt,
                    const RealVect      a_E,
                    const Real          a_dx,
                    const Real          a_kappa) const;

    private:
      // trim from start (in place)
      inline void
      ltrim(std::string& s)
      {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
      }

      // trim from end (in place)
      inline void
      rtrim(std::string& s)
      {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
      }

      // trim from both ends (in place)
      inline void
      trim(std::string& s)
      {
        ltrim(s);
        rtrim(s);
      }
    };
  } // namespace ItoPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_ItoPlasmaPhysicsImplem.H>

#endif
