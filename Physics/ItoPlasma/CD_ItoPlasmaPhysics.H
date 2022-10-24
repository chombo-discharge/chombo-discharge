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
#include <CD_KMCDualState.H>
#include <CD_KMCDualStateReaction.H>
#include <CD_KMCSolver.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoPlasma {

    using KMCState      = KMCDualState<>;
    using KMCReaction   = KMCDualStateReaction<>;
    using ReactionPtrs  = std::deque<std::shared_ptr<ItoPlasmaReaction>>;
    using ReactionTypes = std::pair<ReactionPtrs, ReactionPtrs>;

    /*!
      @brief Base-class for interaction between Kinetic Monte Carlo and Ito-based plasma solvers. 
    */
    class ItoPlasmaPhysics
    {
    public:
      /*!
	@brief Enum for switching between LFA and LEA-based models. Soon deprecated
      */
      enum class coupling
      {
        LFA,
        LEA
      };

      /*!
	@brief Constructor. Does nothing.
      */
      inline ItoPlasmaPhysics();

      /*!
	@brief Destructor. Does nothing. 
      */
      inline virtual ~ItoPlasmaPhysics();

      /*!
	@brief Compute a time step to use. 
	@param[in] a_E         Electric field. 
	@param[in] a_pos       Position
	@param[in] a_densities Mesh-based plasma densities. 
      */
      virtual Real
      computeDt(const RealVect a_E, const RealVect a_pos, const Vector<Real> a_densities) const = 0;

      /*!
	@brief Compute Townsend effective ionization coefficient
	@param[in] a_E Electric field. 
      */
      virtual Real
      computeAlpha(const RealVect a_E) const = 0;

      /*!
	@brief Get all particle species
	@return m_plasmaSpecies
      */
      const Vector<RefCountedPtr<ItoSpecies>>&
      getItoSpecies() const;

      /*!
	@brief Get all photon species
	@return m_rtSpecies
      */
      const Vector<RefCountedPtr<RtSpecies>>&
      getRtSpecies() const;

      /*!
	@brief Return number of ion equations
      */
      inline int
      getNumPlasmaSpecies() const;

      /*!
	@brief Return number of RTE equations
      */
      inline int
      getNumPhotonSpecies() const;

      /*!
	@brief Get the coupling, LFA/LEA
      */
      inline coupling
      getCoupling() const;

      /*!
	@brief Set initial surface charge. Default is 0, override if you want.
      */
      inline virtual Real
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
	@brief Reconcile the number of particles.
	@details This will add/remove particles and potentially also adjust the particle weights.
	@param[inout] a_particles Computational particles
	@param[in]    a_newNumParticles New number of particles (i.e., after the KMC advance)
	@param[in]    a_oldNumParticles Previous number of particles (i.e., before the KMC advance)
	@param[in]    a_cellPos         Cell center position
	@param[in]    a_centroidPos     Cell centroid position
	@param[in]    a_loCorner        Low corner of minimum box enclosing the cut-cell
	@param[in]    a_hiCorner        High corner of minimum box enclosing the cut-cell
	@param[in]    a_bndryCentroid   Cut-cell boundary centroid
	@param[in]    a_bndryNormal     Cut-cell normal (pointing into the domain)
	@param[in]    a_dx              Grid resolution
	@param[in]    a_kappa           Cut-cell volume fraction. 
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
                         const Real                  a_kappa) const noexcept;

      /*!
	@brief Generate new photons. 
	@details This will add photons
	@param[in] a_newPhotons      New photons
	@param[in] a_numNewPhotons   Number of physical photons to be generated. 
	@param[in] a_cellPos         Cell center position
	@param[in] a_centroidPos     Cell centroid position
	@param[in] a_loCorner        Low corner of minimum box enclosing the cut-cell
	@param[in] a_hiCorner        High corner of minimum box enclosing the cut-cell
	@param[in] a_bndryCentroid   Cut-cell boundary centroid
	@param[in] a_bndryNormal     Cut-cell normal (pointing into the domain)
	@param[in] a_dx              Grid resolution
	@param[in] a_kappa           Cut-cell volume fraction. 
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
                       const Real               a_kappa) const noexcept;

      /*!
	@brief Reconcile photoionization reactions.
	@param[inout] a_particles       Particle products. 
	@param[in]    a_absorbedPhotons Photons absorbed on the mesh. 
	@details This runs through the photo-reactions and associates photo-ionization products. 
      */
      inline void
      reconcilePhotoionization(Vector<List<ItoParticle>*>&  a_particles,
                               const Vector<List<Photon>*>& a_absorbedPhotons) const noexcept;

    protected:
      /*!
	@brief Enum for switching between KMC algorithms
	@details 'SSA' is the Gillespie algorithm, 'Tau' is tau-leaping and 'Hybrid' is the Cao et. al. algorithm. 
      */
      enum class Algorithm
      {
        SSA,
        Tau,
        Hybrid
      };

      /*!
	@brief Class name. Used for options parsing
      */
      std::string m_className;

      /*!
	@brief Turn on/off debugging
      */
      bool m_debug;

      /*!
	@brief Algorithm to use for KMC advance
      */
      Algorithm m_algorithm;

      /*!
	@brief Kinetic Monte Carlo state
      */
      KMCState m_kmcState;

      /*!
	@brief Kinetic Monte Carlo solver. 
      */
      KMCSolver<KMCReaction, KMCState> m_kmcSolver;

      /*!
	@brief List of solver-tracked particle species. 
      */
      Vector<RefCountedPtr<ItoSpecies>> m_plasmaSpecies;

      /*!
	@brief List of solver-tracked photon species. 
      */
      Vector<RefCountedPtr<RtSpecies>> m_rtSpecies;

      // Soon to be deprecated.
      coupling m_coupling;
      int      m_numPlasmaSpecies;
      int      m_numRtSpecies;

      /*!
	@brief Maximum new number of particles 
      */
      int m_ppc;

      /*!
	@brief Solver setting for the Cao et. al algorithm. 
	@details Determines critical reactions. A reaction is critical if it is m_Ncrit firings away from depleting a reactant. 
      */
      int m_Ncrit;

      /*!
	@brief Solver setting for the Cao et. al algorithm. 
	@details Maximum number of SSA steps to run when switching into SSA-based advancement for non-critical reactions. 
      */
      int m_NSSA;

      /*!
	@brief Solver setting for the Cao et. al. algorithm. 
	@details Equal to the maximum permitted change in the relative propensity for non-critical reactions
      */
      Real m_SSAlim;

      /*!
	@brief Solver setting for the Cao et. al. algorithm. 
	@details Equal to the maximum permitted change in the relative propensity for non-critical reactions
      */
      Real m_eps;

      /*!
	@brief List of particle reactions
      */
      std::map<std::string, ItoPlasmaReaction> m_reactions;

      /*!
	@brief List of photoionization reactions
      */
      mutable std::map<std::string, ItoPlasmaPhotoReaction> m_photoReactions;

      /*!
	@brief Parse the maximum number of particles generated per cell
      */
      inline void
      parsePPC() noexcept;

      /*!
	@brief Parse the maximum number of particles generated per cell
      */
      inline void
      parseDebug() noexcept;

      /*!
	@brief Parse reaction algorithm
      */
      inline void
      parseAlgorithm() noexcept;

      /*!
	@brief Returns true if all particle numbers are >= 0
      */
      inline bool
      isStateValid(const Vector<long long>& a_particles) const;

      /*!
	@brief Remove particles from the input list.
	@details This will remove weight from the input particles if we can. Otherwise we remove full particles. 
	@param[inout] a_particles   List of (super-)particles to remove from. 
	@param[in]    a_numToRemove Number of physical particles to remove from the input list
      */
      inline void
      removeParticles(List<ItoParticle>& a_particles, const long long a_numToRemove) const;

      /*!
	@brief Get the number of physical particles. 
	@param[in] a_particles List of plasma species particles. 
      */
      inline Vector<long long>
      getNumPhysicalParticles(const Vector<List<ItoParticle>*>& a_particles) const;

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
    };
  } // namespace ItoPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_ItoPlasmaPhysicsImplem.H>

#endif
