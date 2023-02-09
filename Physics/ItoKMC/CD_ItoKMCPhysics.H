/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCPhysics.H
  @brief  Main file for describing Ito-based plasma physics
  @author Robert Marskar
*/

#ifndef CD_ItoKMCPhysics_H
#define CD_ItoKMCPhysics_H

// Std includes
#include <memory>
#include <vector>

// Chombo includes
#include <RealVect.H>
#include <RefCountedPtr.H>
#include <List.H>

// Our includes
#include <CD_ItoSpecies.H>
#include <CD_CdrSpecies.H>
#include <CD_RtSpecies.H>
#include <CD_Photon.H>
#include <CD_ItoParticle.H>
#include <CD_PointParticle.H>
#include <CD_ItoKMCPhotoReaction.H>
#include <CD_KMCDualState.H>
#include <CD_KMCDualStateReaction.H>
#include <CD_KMCSolver.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoKMC {

    /*!
      @brief Numerical representation of the KMC state. Can be floating-point or integer type
    */
    using FPR = Real;

    /*!
      @brief KMC state used in the Kinetic Monte Carlo advancement
    */
    using KMCState = KMCDualState<FPR>;

    /*!
      @brief KMC reaction used in the Kinetic Monte Carlo advancement
    */
    using KMCReaction = KMCDualStateReaction<KMCState, FPR>;

    /*!
      @brief KMC solverl type used in the Kinetic Monte carlo advancement
    */
    using KMCSolverType = KMCSolver<KMCReaction, KMCState, FPR>;

    /*!
      @brief Map to species type
      @details This is just for distinguishing between species that are treated with an Ito or CDR formalism
    */
    enum SpeciesType
    {
      Ito,
      CDR
    };

    /*!
      @brief Base class for interaction between Kinetic Monte Carlo and Ito-based plasma solvers. 
      @details When using this class, the user should populate m_kmcReactions with the appropriate reactions AND implement routines
      for computing the mobility/diffusion coefficeints. 
    */
    class ItoKMCPhysics
    {
    public:
      /*!
	@brief Constructor. Does nothing.
      */
      inline ItoKMCPhysics() noexcept;

      /*!
	@brief Destructor. Does nothing. 
      */
      inline virtual ~ItoKMCPhysics() noexcept;

      /*!
	@brief Compute a time step to use. 
	@param[in] a_E            Electric field. 
	@param[in] a_pos          Position
	@param[in] a_numParticles Number of reactive particles of each species
	@return Default returns infinity but user can implement a physics based time step. 
      */
      inline virtual Real
      computeDt(const RealVect a_E, const RealVect a_pos, const Vector<FPR> a_numParticles) const noexcept;

      /*!
	@brief Compute Townsend ionization coefficient
	@param[in] a_E Electric field. 
	@return Should return the Townsend ionization coefficient. 
      */
      virtual Real
      computeAlpha(const RealVect a_E) const = 0;

      /*!
	@brief Get all particle drift-diffusion species
	@return m_itoSpecies
      */
      const Vector<RefCountedPtr<ItoSpecies>>&
      getItoSpecies() const;

      /*!
	@brief Get all fluid drift-diffusion species
	@return m_cdrSpecies
      */
      const Vector<RefCountedPtr<CdrSpecies>>&
      getCdrSpecies() const;

      /*!
	@brief Get all photon species
	@return m_rtSpecies
      */
      const Vector<RefCountedPtr<RtSpecies>>&
      getRtSpecies() const;

      /*!
	@brief Return number of Ito solvers.
      */
      inline int
      getNumItoSpecies() const;

      /*!
	@brief Return number of CDR solvers.
      */
      inline int
      getNumCdrSpecies() const;

      /*!
	@brief Return total number of plasma species.
      */
      inline int
      getNumPlasmaSpecies() const;

      /*!
	@brief Return number of RTE solvers.
      */
      inline int
      getNumPhotonSpecies() const;

      /*!
	@brief Get the internal mapping between plasma species and Ito solvers
      */
      inline const std::map<int, std::pair<SpeciesType, int>>&
      getSpeciesMap() const noexcept;

      /*!
	@brief Parse run-time options
      */
      inline virtual void
      parseRuntimeOptions() noexcept;

      /*!
	@brief Set initial surface charge. Default is 0, override if you want.
	@param[in] a_time Simulation time
	@param[in] a_pos  Physical coordinate
      */
      inline virtual Real
      initialSigma(const Real a_time, const RealVect a_pos) const;

      /*!
	@brief Compute the Ito solver mobilities.
	@param[in] a_time Time
	@param[in] a_pos  Position
	@param[in] a_E    Electric field
	@return Must return a vector of non-negative mobility coefficients for the plasma species
      */
      virtual Vector<Real>
      computeMobilities(const Real a_time, const RealVect a_pos, const RealVect a_E) const noexcept = 0;

      /*!
	@brief Compute the Ito solver diffusion coefficients
	@param[in] a_time Time
	@param[in] a_pos  Position
	@param[in] a_E    Electric field
	@return Must return a vector of non-negative diffusion coefficients for the plasma species
      */
      virtual Vector<Real>
      computeDiffusionCoefficients(const Real a_time, const RealVect a_pos, const RealVect a_E) const noexcept = 0;

      /*!
	@brief Resolve particle and photon injection at the EB.
	@details Routine is here to handle charge injection, secondary emission etc.
	@param[out] a_outgoingParticles Outgoing plasma species particles. 
	@param[out] a_outgoingPhotons   Photons injected through the EB
	@param[in]  a_incomingParticles Particles that left the computational domain through the EB
	@param[in]  a_incomingPhotons   Photons that left the computational domain through the EB
	@param[in]  a_newNumParticles   Total number of particles in the cut-cell AFTER the transport step
	@param[in]  a_oldNumParticles   Total number of particles in the cut-cell BEFORE the transport step
	@param[in]  a_electricField     Electric field                  
	@param[in]  a_cellCenter        Physical position of the cell center. 
	@param[in]  a_cellCentroid      Cell centroid relative to the cell center (not multiplied by dx)
	@param[in]  a_bndryCentroid     EB face centroid relative to the cell center (not multiplied by dx)
	@param[in]  a_bndryNormal       Cut-cell normal vector. 
	@param[in]  a_bndryArea         Cut-cell boundary area - not multiplied by dx (2D) or dx^2 (3D)
	@param[in]  a_dx                Grid resolution on this level. 
	@param[in]  a_dt                Time step
	@param[in]  a_isDielectric      Dielectric or electrode. 
	@param[in]  a_matIndex          Material index (taken from computationalGeometry)
      */
      virtual void
      injectParticlesEB(Vector<List<ItoParticle>>&       a_outgoingParticles,
                        Vector<List<Photon>>&            a_outgoingPhotons,
                        const Vector<List<ItoParticle>>& a_incomingParticles,
                        const Vector<List<Photon>>&      a_incomingPhotons,
                        const Vector<FPR>&               a_newNumParticles,
                        const Vector<FPR>&               a_oldNumParticles,
                        const RealVect&                  a_E,
                        const RealVect&                  a_cellCenter,
                        const RealVect&                  a_cellCentroid,
                        const RealVect&                  a_bndryCentroid,
                        const RealVect&                  a_bndryNormal,
                        const Real                       a_bndryArea,
                        const Real                       a_dx,
                        const Real                       a_dt,
                        const bool                       a_isDielectric,
                        const int                        a_matIndex) const noexcept = 0;

      /*!
	@brief Advance particles. 
	@param[inout] a_numParticles  Number of physical particles
	@param[out]   a_numNewPhotons Number of new physical photons to generate (of each type)
	@param[in]    a_dt            Time step
	@param[in]    a_E             Electric field
	@param[in]    a_dx            Grid resolution
	@param[in]    a_kappa         Cut-cell volume fraction. 
      */
      inline void
      advanceKMC(Vector<FPR>&   a_numParticles,
                 Vector<FPR>&   a_numNewPhotons,
                 const Real     a_dt,
                 const RealVect a_E,
                 const Real     a_dx,
                 const Real     a_kappa) const;

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
	@note Public because this is called by ItoKMCStepper
      */
      inline void
      reconcileParticles(Vector<List<ItoParticle>*>& a_particles,
                         const Vector<FPR>&          a_newNumParticles,
                         const Vector<FPR>&          a_oldNumParticles,
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
	@note Public because this is called by ItoKMCStepper
      */
      inline void
      reconcilePhotons(Vector<List<Photon>*>& a_newPhotons,
                       const Vector<FPR>&     a_numNewPhotons,
                       const RealVect         a_cellPos,
                       const RealVect         a_centroidPos,
                       const RealVect         a_lo,
                       const RealVect         a_hi,
                       const RealVect         a_bndryCentroid,
                       const RealVect         a_bndryNormal,
                       const Real             a_dx,
                       const Real             a_kappa) const noexcept;

      /*!
	@brief Reconcile photoionization reactions.
	@param[inout] a_itoPrticles       Particle products placed in Ito solvers. 
	@param[inout] a_cdrParticles      Particle products placed in CDR solvers. 
	@param[in]    a_absorbedPhotons Photons absorbed on the mesh. 
	@details This runs through the photo-reactions and associates photo-ionization products. 
      */
      inline void
      reconcilePhotoionization(Vector<List<ItoParticle>*>&   a_itoParticles,
                               Vector<List<PointParticle>*>& a_cdrParticles,
                               const Vector<List<Photon>*>&  a_absorbedPhotons) const noexcept;

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
	@brief Enum for switching between various particle placement algorithms
      */
      enum class ParticlePlacement
      {
        Random,
        Centroid
      };

      /*!
	@brief Algorithm to use for KMC advance
      */
      Algorithm m_algorithm;

      /*!
	@brief Particle placement algorithm
      */
      ParticlePlacement m_particlePlacement;

      /*!
	@brief Map for associating a plasma species with an Ito solver. 
      */
      std::map<int, std::pair<SpeciesType, int>> m_speciesMap;

      /*!
	@brief Class name. Used for options parsing
      */
      std::string m_className;

      /*!
	@brief Turn on/off debugging
      */
      bool m_debug;

      /*!
	@brief Is defined or not
      */
      bool m_isDefined;

      /*!
	@brief Kinetic Monte Carlo state
      */
      mutable KMCState m_kmcState;

      /*!
	@brief Kinetic Monte Carlo solver. 
      */
      mutable KMCSolverType m_kmcSolver;

      /*!
	@brief List of reactions for the KMC solver
      */
      std::vector<std::shared_ptr<const KMCReaction>> m_kmcReactions;

      /*!
	@brief List of photoionization reactions
      */
      std::vector<std::shared_ptr<const ItoKMCPhotoReaction>> m_photoReactions;

      /*!
	@brief Random number generators for photoionization pathways. 
	@details The first index is the photon index, i.e. entry in m_photonSpecies. The second index in the
	map is an RNG generator for selecting photo-reactions. 
      */
      std::map<int, std::pair<std::discrete_distribution<int>, std::map<int, int>>> m_photoPathways;

      /*!
	@brief List of solver-tracked particle drift-diffusion species. 
      */
      Vector<RefCountedPtr<ItoSpecies>> m_itoSpecies;

      /*!
	@brief List of solver-tracked fluid drift-diffusion species. 
      */
      Vector<RefCountedPtr<CdrSpecies>> m_cdrSpecies;

      /*!
	@brief List of solver-tracked photon species. 
      */
      Vector<RefCountedPtr<RtSpecies>> m_rtSpecies;

      /*!
	@brief Maximum new number of particles generated by the chemistry advance
      */
      int m_maxNewParticles;

      /*!
	@brief Maximum new number of photons generated by the chemistry advance
      */
      int m_maxNewPhotons;

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
	@brief Define method -- defines all the internal machinery
      */
      inline void
      define() noexcept;

      /*!
	@brief Build internal representation of how we distinguish the Ito and CDR solvers
	@note This should ALWAYS be called after initializing the species since ItoKMCStepper will rely on it. 
      */
      inline void
      defineSpeciesMap() noexcept;

      /*!
	@brief Define the KMC solver and state
      */
      inline void
      defineKMC() noexcept;

      /*!
	@brief Define pathways for photo-reactions
      */
      inline void
      definePhotoPathways() noexcept;

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
	@brief Update reaction rates
	@param[in] a_E     Electric field
	@param[in] a_dx    Grid resolution
	@param[in] a_kappa Cut-cell volume fraction
	@note Must be implemented by the user.
      */
      virtual void
      updateReactionRates(const RealVect a_E, const Real a_dx, const Real a_kappa) const noexcept = 0;

      /*!
	@brief Remove particles from the input list.
	@details This will remove weight from the input particles if we can. Otherwise we remove full particles. 
	@param[inout] a_particles   List of (super-)particles to remove from. 
	@param[in]    a_numToRemove Number of physical particles to remove from the input list
      */
      inline void
      removeParticles(List<ItoParticle>& a_particles, const long long a_numToRemove) const;
    };
  } // namespace ItoKMC
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_ItoKMCPhysicsImplem.H>

#endif
