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
#include <memory>
#include <vector>

// Chombo includes
#include <RealVect.H>
#include <RefCountedPtr.H>
#include <List.H>

// Our includes
#include <CD_ItoSpecies.H>
#include <CD_RtSpecies.H>
#include <CD_Photon.H>
#include <CD_ItoParticle.H>
#include <CD_ItoPlasmaPhotoReaction.H>
#include <CD_KMCDualState.H>
#include <CD_KMCDualStateReaction.H>
#include <CD_KMCSolver.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoPlasma {

    using FPR           = long long;
    using KMCState      = KMCDualState<FPR>;
    using KMCReaction   = KMCDualStateReaction<KMCState, FPR>;
    using KMCSolverType = KMCSolver<KMCReaction, KMCState, FPR>;

    /*!
      @brief Base class for interaction between Kinetic Monte Carlo and Ito-based plasma solvers. 
      @details When using this class, the user should populate m_kmcReactions with the appropriate reactions AND implement routines
      for computing the mobility/diffusion coefficeints. 
    */
    class ItoPlasmaPhysics
    {
    public:
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
	@return Default returns infinity but user can implement a physics based time step. 
      */
      inline virtual Real
      computeDt(const RealVect a_E, const RealVect a_pos, const Vector<Real> a_densities) const noexcept;

      /*!
	@brief Compute Townsend ionization coefficient
	@param[in] a_E Electric field. 
	@return Should return the Townsend ionization coefficient. 
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
      computeItoMobilities(const Real a_time, const RealVect a_pos, const RealVect a_E) const noexcept = 0;

      /*!
	@brief Compute the Ito solver diffusion coefficients
	@param[in] a_time Time
	@param[in] a_pos  Position
	@param[in] a_E    Electric field
	@param[in] a_densities Mesh-based densities for the plasma species
	@return Must return a vector of non-negative diffusion coefficients for the plasma species
      */
      virtual Vector<Real>
      computeItoDiffusion(const Real         a_time,
                          const RealVect     a_pos,
                          const RealVect     a_E,
                          const Vector<Real> a_densities) const noexcept = 0;

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
      advanceKMC(Vector<long long>& a_numParticles,
                 Vector<long long>& a_numNewPhotons,
                 const Real         a_dt,
                 const RealVect     a_E,
                 const Real         a_dx,
                 const Real         a_kappa) const;

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
	@note Public because this is called by ItoPlasmaStepper
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
	@note Public because this is called by ItoPlasmaStepper
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
      std::vector<std::shared_ptr<const ItoPlasmaPhotoReaction>> m_photoReactions;

      /*!
	@brief List of solver-tracked particle species. 
      */
      Vector<RefCountedPtr<ItoSpecies>> m_plasmaSpecies;

      /*!
	@brief List of solver-tracked photon species. 
      */
      Vector<RefCountedPtr<RtSpecies>> m_rtSpecies;

      /*!
	@brief Maximum new number of particles 
      */
      long long m_ppc;

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
	@brief Define the KMC solver and state
      */
      inline void
      defineKMC() noexcept;

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
  } // namespace ItoPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_ItoPlasmaPhysicsImplem.H>

#endif
