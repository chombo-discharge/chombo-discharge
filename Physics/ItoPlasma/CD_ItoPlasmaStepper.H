/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaStepper.H
  @author Robert Marskar
  @brief  Abstract class for integrating the Ito plasma equations
*/

#ifndef CD_ItoPlasmaStepper_H
#define CD_ItoPlasmaStepper_H

// Std includes
#include <functional>

// Our includes
#include <CD_TimeStepper.H>
#include <CD_ItoPlasmaPhysics.H>
#include <CD_ItoLayout.H>
#include <CD_PointParticle.H>
#include <CD_RtLayout.H>
#include <CD_McPhoto.H>
#include "CD_FieldSolver.H"
#include "CD_SurfaceODESolver.H"
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoPlasma {

    /*!
      @brief Enum for differentiating between types of particles
    */
    enum class WhichParticles
    {
      All,
      AllMobile,
      AllDiffusive,
      ChargedMobile,
      ChargedDiffusive,
      AllMobileOrDiffusive,
      ChargedAndMobileOrDiffusive,
      Stationary,
    };

    /*!
      @brief Time stepper that advances the Ito-Poisson system of equations
    */
    class ItoPlasmaStepper : public TimeStepper
    {
    public:

      /*!
	@brief Default constructor. Sets default options. 
      */
      ItoPlasmaStepper();

      /*!
	@brief Full constructor. Calls the other constructor
	@param[in] a_physics Implementation of the Ito plasma physics model.
      */
      ItoPlasmaStepper(RefCountedPtr<ItoPlasmaPhysics>& a_physics);

      /*!
	@brief Destructor
      */
      virtual
      ~ItoPlasmaStepper();

      // Setup routines
      virtual void
      parseOptions() = 0;
      virtual void
      parseRuntimeOptions() = 0;
      virtual void
      allocateInternals() = 0;
      virtual void
      setupSolvers() override;
      virtual void
      allocate() override;
      virtual void
      initialData() override;
      virtual void
      initialSigma();
      virtual void
      postCheckpointSetup() override;
      virtual void
      postCheckpointPoisson();
      virtual void
      postInitialize() override;

      void
      setupIto();
      void
      setupPoisson();
      void
      setupRadiativeTransfer();
      void
      setupSigma();

      // IO routines
#ifdef CH_USE_HDF5
      void
      writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
#endif
#ifdef CH_USE_HDF5
      void
      readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
#endif
      virtual void
      writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const override;
      void
      writeJ(EBAMRCellData& a_output, int& a_icomp) const;
      void
      writeNumParticlesPerPatch(EBAMRCellData& a_output, int& a_icomp) const;
      int
      getNumberOfPlotVariables() const override;
      Vector<long int>
      getCheckpointLoads(const std::string a_realm, const int a_level) const override;

      // Advance routines
      virtual Real
      advance(const Real a_dt) = 0;
      virtual Real
      computeDt() override;
      void
      synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;
      void
      printStepReport() override;

      // Registration routines
      void
      registerRealms() override;
      void
      registerOperators() override;

      // Regrid routines.
      virtual void
      preRegrid(const int a_lmin, const int a_oldFinestLevel) override;
      virtual void
      regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;
      void
      postRegrid() override;
      bool
      loadBalanceThisRealm(const std::string a_realm) const override;
      void
      loadBalanceBoxes(Vector<Vector<int>>&             a_procs,
                       Vector<Vector<Box>>&             a_boxes,
                       const std::string                a_realm,
                       const Vector<DisjointBoxLayout>& a_grids,
                       const int                        a_lmin,
                       const int                        a_finestLevel) override;

      // Particle handling routines.
      void
      intersectParticles(const WhichParticles a_which, const EbRepresentation a_representation, const bool a_delete);
      void
      intersectParticles(const WhichParticles            a_which,
                         const ItoSolver::WhichContainer a_particles,
                         const ItoSolver::WhichContainer a_eb_particles,
                         const ItoSolver::WhichContainer a_domain_particles,
                         const EbRepresentation          a_representation,
                         const bool                      a_delete);

      void
      removeCoveredParticles(const WhichParticles   a_which,
                             const EbRepresentation a_representation,
                             const Real             a_tolerance);
      void
      removeCoveredParticles(const WhichParticles            a_which,
                             const ItoSolver::WhichContainer a_container,
                             const EbRepresentation          a_representation,
                             const Real                      a_tolerance);

      void
      transferCoveredParticles(const WhichParticles   a_which,
                               const EbRepresentation a_representation,
                               const Real             a_tolerance);
      void
      transferCoveredParticles(const WhichParticles            a_which,
                               const ItoSolver::WhichContainer a_containerFrom,
                               const ItoSolver::WhichContainer a_containerTo,
                               const EbRepresentation          a_representation,
                               const Real                      a_tolerance);

      // Internal routines
      void
      setVerbosity(const int a_verbosity);

      // Set routines
      void
      setIto(RefCountedPtr<ItoLayout<ItoSolver>>& a_cdr);
      void
      setFieldSolver(RefCountedPtr<FieldSolver>& a_poisson);
      void
      setRadiativeTransferSolvers(RefCountedPtr<RtLayout<McPhoto>>& a_rte);
      void
      setVoltage(const std::function<Real(const Real a_time)>& a_potential);

      // Get routines
      Real
      getTime() const;
      Real
      computeMaxElectricField(const phase::which_phase a_phase);

      // E-field computations
      void
      computeElectricField(MFAMRCellData& a_E, const MFAMRCellData& a_potential);
      void
      computeElectricField(EBAMRCellData& a_E, const phase::which_phase a_phase);
      void
      computeElectricField(EBAMRCellData& a_E, const phase::which_phase a_phase, const MFAMRCellData& a_potential);
      void
      computeElectricField(EBAMRFluxData& a_E_face, const phase::which_phase a_phase, const EBAMRCellData& a_E_cell);
      void
      computeElectricField(EBAMRIVData& a_E_eb, const phase::which_phase a_phase, const EBAMRCellData& a_E_cell);

      // Compute rho routines
      void
      computeSpaceChargeDensity();
      void
      computeSpaceChargeDensity(EBAMRCellData& a_rho, const phase::which_phase a_phase);
      void
      computeSpaceChargeDensity(MFAMRCellData& a_rho, const Vector<EBAMRCellData*>& a_densities);

      // Compute conductivity routines
      void
      computeConductivity(EBAMRCellData& a_conductivity);
      void
      computeConductivity(EBAMRCellData& a_conductivity, const Vector<ParticleContainer<ItoParticle>*>& a_particles);

      // Compute J routines
      void
      computeJ(EBAMRCellData& a_J, const Real a_dt);

      // Relaxation time
      Real
      computeRelaxationTime();
      Real
      computeRelaxationTime(const int a_level);
      Real
      computeRelaxationTime(const int a_level, const DataIndex a_dit);

      // Poisson-solving routines
      bool
      solvePoisson();
      bool
      solvePoisson(MFAMRCellData&                a_potential,
                   MFAMRCellData&                a_rhs,
                   const Vector<EBAMRCellData*>& a_densities,
                   const EBAMRIVData&            a_sigma);

      // Remap and deposition routines for solver subsets.
      void
      depositParticles(const WhichParticles a_WhichParticles);
      void
      depositParticles(const WhichParticles a_WhichParticles, const ItoSolver::WhichContainer a_container);
      void
      remapParticles(const WhichParticles a_WhichParticles);
      void
      remapParticles(const WhichParticles a_WhichParticles, const ItoSolver::WhichContainer a_container);

      // This sets the velocity function for each species. This is sgn(charge)*E
      void
      setItoVelocityFunctions();

      // Routines for computing velocities and diffusion coefficients
      void
      computeItoVelocities();
      void
      computeItoDiffusion();

      // Routines for computing mobilities in the LFA. This computes mesh-defined mobilities and interpolates them to the particle positions
      void
      computeItoMobilitiesLFA();
      void
      computeItoMobilitiesLFA(Vector<EBAMRCellData*>& a_meshMobilities, const EBAMRCellData& a_E, const Real a_time);
      void
      computeItoMobilitiesLFA(Vector<LevelData<EBCellFAB>*>& a_meshMobilities,
                              const LevelData<EBCellFAB>&    a_E,
                              const int                      a_level,
                              const Real                     a_time);
      void
      computeItoMobilitiesLFA(Vector<EBCellFAB*>& a_meshMobilities,
                              const EBCellFAB&    a_E,
                              const int           a_level,
                              const DataIndex     a_dit,
                              const Box           a_box,
                              const Real          a_time);

      // Routines for computing mobilities in the LEA. This runs through all the particles and calls a mobility computation function,
      // which you must implement in ItoSpecies
      void
      computeItoMobilitiesLEA();
      void
      computeItoMobilitiesLEA(Vector<ParticleContainer<ItoParticle>*>& a_particles);
      void
      computeItoMobilitiesLEA(Vector<ParticleData<ItoParticle>*>& a_particles);

      // Routines for computing diffusion coefficients.
      // These routines compute the diffusion coefficients by computing D(E)
      // as a mesh variable, and then interpolating that to the particle position.
      void
      computeItoDiffusionLFA();
      void
      computeItoDiffusionLFA(Vector<EBAMRCellData*>&       a_diffusionCoefficient,
                             const Vector<EBAMRCellData*>& a_densities,
                             const EBAMRCellData&          a_E,
                             const Real                    a_time);
      void
      computeItoDiffusionLFA(Vector<LevelData<EBCellFAB>*>&       a_diffusionCoefficient,
                             const Vector<LevelData<EBCellFAB>*>& a_densities,
                             const LevelData<EBCellFAB>&          a_E,
                             const int                            a_level,
                             const Real                           a_time);
      void
      computeItoDiffusionLFA(Vector<EBCellFAB*>&       a_diffusionCoefficient,
                             const Vector<EBCellFAB*>& a_densities,
                             const EBCellFAB&          a_E,
                             const int                 a_level,
                             const DataIndex           a_dit,
                             const Box                 a_box,
                             const Real                a_time);

      // Routines for computing diffusion coefficients in the LEA
      void
      computeItoDiffusion_lea();

      // Routines for computing the number of particles per cell. This is done on the particle Realm.
      void
      computeReactiveParticlesPerCell(EBAMRCellData& a_ppc);
      void
      computeReactiveParticlesPerCell(LevelData<EBCellFAB>& a_ppc, const int a_level);
      void
      computeReactiveParticlesPerCell(EBCellFAB&      a_ppc,
                                      const int       a_level,
                                      const DataIndex a_dit,
                                      const Box       a_box,
                                      const EBISBox&  a_ebisbox);

      // Routines for computing the mean energy of the particles in grid cells
      void
      computeReactiveMeanEnergiesPerCell(EBAMRCellData& a_mean_energies);
      void
      computeReactiveMeanEnergiesPerCell(LevelData<EBCellFAB>& a_mean_energies, const int a_level);
      void
      computeReactiveMeanEnergiesPerCell(EBCellFAB&      a_mean_energies,
                                         const int       a_level,
                                         const DataIndex a_dit,
                                         const Box       a_box,
                                         const EBISBox&  a_ebisbox);

      // New routines for chemistry advances.
      void
      advanceReactionNetworkNWO(const Real a_dt);
      void
      advanceReactionNetworkNWO(const EBAMRCellData& a_E, const EBAMRCellData& a_EdotJ, const Real a_dt);
      void
      advanceReactionNetworkNWO(LevelData<EBCellFAB>&       a_ppc,
                                LevelData<EBCellFAB>&       a_ypc,
                                LevelData<EBCellFAB>&       a_mean_energies,
                                const LevelData<EBCellFAB>& a_E,
                                const LevelData<EBCellFAB>& a_EdotJ,
                                const int                   a_level,
                                const Real                  a_dt);
      void
      advanceReactionNetworkNWO(EBCellFAB&       a_particlesPerCell,
                                EBCellFAB&       a_newPhotonsPerCell,
                                EBCellFAB&       a_meanParticleEnergies,
                                const EBCellFAB& a_E,
                                const EBCellFAB& a_EdotJ,
                                const int        a_level,
                                const DataIndex  a_dit,
                                const Box        a_box,
                                const Real       a_dx,
                                const Real       a_dt);

      // Routines for reconciling particles.
      void
      reconcileParticles(const EBAMRCellData& a_newParticlesPerCell,
                         const EBAMRCellData& a_oldParticlesPerCell,
                         const EBAMRCellData& a_meanParticleEnergies,
                         const EBAMRCellData& a_newPhotonsPerCell);
      void
      reconcileParticles(const LevelData<EBCellFAB>& a_newParticlesPerCell,
                         const LevelData<EBCellFAB>& a_oldParticlesPerCell,
                         const LevelData<EBCellFAB>& a_meanParticleEnergies,
                         const LevelData<EBCellFAB>& a_newPhotonsPerCell,
                         const int                   a_level);
      void
      reconcileParticles(const EBCellFAB& a_newParticlesPerCell,
                         const EBCellFAB& a_oldParticlesPerCell,
                         const EBCellFAB& a_meanParticleEnergies,
                         const EBCellFAB& a_newPhotonsPerCell,
                         const int        a_level,
                         const DataIndex  a_dit,
                         const Box        a_box,
                         const Real       a_dx);

      // Old routines for chemistry advances. The physics module will switch between LFA and LEA.
      void
      advanceReactionNetwork(const Real a_dt);
      void
      advanceReactionNetwork(Vector<ParticleContainer<ItoParticle>*>& a_particles,
                             Vector<ParticleContainer<Photon>*>&      a_Photons,
                             Vector<ParticleContainer<Photon>*>&      a_newPhotons,
                             const Vector<EBAMRCellData>&             a_sources,
                             const EBAMRCellData&                     a_E,
                             const Real                               a_dt);
      void
      advanceReactionNetwork(Vector<AMRCellParticles<ItoParticle>*>& a_particles,
                             Vector<AMRCellParticles<Photon>*>&      a_Photons,
                             Vector<AMRCellParticles<Photon>*>&      a_newPhotons,
                             const Vector<EBAMRCellData>&            a_sources,
                             const EBAMRCellData&                    a_E,
                             const Real                              a_dt);
      void
      advanceReactionNetwork(Vector<LayoutData<BinFab<ItoParticle>>*>& a_particles,
                             Vector<LayoutData<BinFab<Photon>>*>&      a_Photons,
                             Vector<LayoutData<BinFab<Photon>>*>&      a_newPhotons,
                             const Vector<LevelData<EBCellFAB>*>&      a_sources,
                             const LevelData<EBCellFAB>&               a_E,
                             const int                                 a_lvl,
                             const Real                                a_dt);
      void
      advanceReactionNetwork(Vector<BinFab<ItoParticle>*>& a_particles,
                             Vector<BinFab<Photon>*>&      a_Photons,
                             Vector<BinFab<Photon>*>&      a_newPhotons,
                             const Vector<EBCellFAB*>&     a_sources,
                             const EBCellFAB&              a_E,
                             const int                     a_lvl,
                             const DataIndex               a_dit,
                             const Box                     a_box,
                             const Real                    a_dx,
                             const Real                    a_dt);

      // Compute dt routines
      Real
      computePhysicsDt() const;
      Real
      computePhysicsDt(const EBAMRCellData& a_E, const Vector<EBAMRCellData*> a_densities) const;
      Real
      computePhysicsDt(const LevelData<EBCellFAB>&         a_E,
                       const Vector<LevelData<EBCellFAB>*> a_densities,
                       const int                           a_level) const;
      Real
      computePhysicsDt(const EBCellFAB&         a_E,
                       const Vector<EBCellFAB*> a_densities,
                       const int                a_level,
                       const DataIndex          a_dit,
                       const Box                a_box) const;

      // Routines for advancing Photons
      void
      advancePhotons(const Real a_dt);

      // Routines for sorting particles. These routines are here because RtLayout is general for RtSolver, and I'm too
      // lazy to write an McPhoto_layout. So here you go.
      void
      sortPhotonsByCell(McPhoto::WhichContainer a_which);
      void
      sortPhotonsByPatch(McPhoto::WhichContainer a_which);

      // Routines for computing EdotJ per species, including the diffusive contribution. a_source is defined on the particle Realm.
      void
      computeEdotJSource(const Real a_dt);
      void
      computeEdotJSource();
      void
      computeEdotJSourceNWO();
      void
      computeEdotJSourceNWO2(const Real a_dt);

    protected:
      /*!
	@brief An enum for encapsulating how time steps were restricted. 
      */
      enum class TimeCode
      {
        Advection,
        Diffusion,
        AdvectionDiffusion,
        Source,
        RelaxationTime,
        Restricted,
        Hardcap,
        Error,
        Physics
      };

      std::string        m_fluidRealm;   // Realm for Poisson and sigma solver. Defaults to Realm::Primal
      std::string        m_particleRealm; // Realm for particles. Defaults to Realm::Primal
      std::string        m_name;
      phase::which_phase m_phase;

      BoxSorting m_boxSort;

      RefCountedPtr<ItoPlasmaPhysics>     m_physics;
      RefCountedPtr<ItoLayout<ItoSolver>> m_ito;
      RefCountedPtr<RtLayout<McPhoto>>    m_rte;
      RefCountedPtr<FieldSolver>          m_fieldSolver;
      RefCountedPtr<SurfaceODESolver<1>>  m_sigma;

      TimeCode m_timeCode;

      // Filters that can be applied.
      std::deque<std::tuple<Real, int, int>> m_filters;

      std::function<Real(const Real a_time)> m_potential;

      bool m_regridSuperparticles;
      bool m_LoadBalancing;
      bool m_nwo_reactions; // For switching between the old and new reaction algorithms.
      bool m_filter_rho;

      int m_ppc;
      int m_LoadBalancing_idx;

      Real m_load_ppc;
      Real m_max_cells_hop;
      Real m_dt_relax;
      Real m_min_dt;
      Real m_max_dt;

      // Storage for current density. This is defined over the particle Realm.
      EBAMRCellData         m_J;
      mutable EBAMRCellData m_fluid_scratch1; // Scratch storage on the fluid Realm with 1 component
      mutable EBAMRCellData m_fluid_scratchD; // Scratch storage on the fluid Realm with SpaceDim components

      mutable EBAMRCellData m_particle_scratch1; // Scratch storage on the particle Realm with 1 component
      mutable EBAMRCellData m_particle_scratchD; // Scratch storage on the particle Realm with SpaceDim components

      mutable Vector<EBAMRCellData>
        m_fscratch1; // Scratch storage with a single component for all species. On the fluid Realm.
      mutable Vector<EBAMRCellData>
        m_fscratch2; // Scratch storage with a single component for all species. On the fluid Realm.

      EBAMRCellData
                    m_particle_ppc; // For holding the number of particles per cell for all species. Contains compoents = number of species
      EBAMRCellData m_particle_ypc; // For holding the number of generated Photons per cell on the particle Realm.
      EBAMRCellData m_particle_eps; // For holding the mean energy on the particle Realm
      EBAMRCellData m_particle_old; // For holding the mean energy on the particle Realm

      EBAMRCellData
                    m_fluid_ppc; // For holding the number of particles per cell for all species. Contains compoents = number of species
      EBAMRCellData m_fluid_ypc; // For holding the number of generated Photons per cell on the fluid Realm.
      EBAMRCellData m_fluid_eps; // For holding the mean energy on the particle Realm

      // Storage for the electric field on
      EBAMRCellData m_fluid_E;
      EBAMRCellData m_particle_E;

      // Storage for energy source terms
      Vector<EBAMRCellData> m_energy_sources;
      EBAMRCellData         m_EdotJ;

      // For parsing filters
      void
      parseFilters();

      // Set the particle buffers
      void
      setParticleBuffers();

      // Diagnostics tool
      void
      getParticleStatistics(Real&   a_avg,
                            Real&   a_sigma,
                            size_t& a_minParticles,
                            size_t& a_maxParticles,
                            int&    a_minRank,
                            int&    a_maxRank);
      void
      printTimerHead();
      void
      printTimerDiagnostics(Real& a_timer, const std::string a_prefix); // Timer should be in percentage
      void
      printTimerTail();

      // Do load balancing of the particle Realm.
      void
      loadBalanceParticleRealm(Vector<Vector<int>>&             a_procs,
                               Vector<Vector<Box>>&             a_boxes,
                               const std::string                a_realm,
                               const Vector<DisjointBoxLayout>& a_grids,
                               const int                        a_lmin,
                               const int                        a_finestLevel);

      Vector<RefCountedPtr<ItoSolver>>
      getLoadBalanceSolvers() const;
    };
  } // namespace ItoPlasma
};  // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
