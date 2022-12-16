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
    enum class SpeciesSubset
    {
      All,
      AllMobile,
      AllDiffusive,
      AllMobileOrDiffusive,
      AllMobileAndDiffusive,
      Charged,
      ChargedMobile,
      ChargedDiffusive,
      ChargedMobileOrDiffusive,
      ChargedMobileAndDiffusive,
      Stationary,
    };

    /*!
      @brief Base time stepper class that advances the Ito-Poisson system of equations
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
      virtual ~ItoPlasmaStepper();

      /*!
	@brief Parse options
      */
      virtual void
      parseOptions();

      /*!
	@brief Parse runtime configurable options
      */
      virtual void
      parseRuntimeOptions() override;

      /*!
	@brief Allocate storage for solvers.
      */
      virtual void
      allocate() override;

      /*!
	@brief Allocate "internal" storage.
      */
      virtual void
      allocateInternals();

      /*!
	@brief Set up solvers.
      */
      virtual void
      setupSolvers() override;

      /*!
	@brief Fill solvers with initial data
      */
      virtual void
      initialData() override;

      /*!
	@brief Perform post-checkpoint operations
	@details This allocates the internal storage, recomputes the electric field and recomputes the velocity and diffusion coefficients. 
      */
      virtual void
      postCheckpointSetup() override;

      /*!
	@brief Do some post-checkpoint operations for the electrostatic part. 
	@details This computes the electric field from the checkpointed potential and populates m_electricFieldFluid and m_electricFieldParticle.
      */
      virtual void
      postCheckpointPoisson();

      /*!
	@brief Post-initialization operations. Default does nothing
      */
      virtual void
      postInitialize() override;

#ifdef CH_USE_HDF5
      /*!
	@brief Read header data from checkpoint file.
	@param[inout] a_header HDF5 header.
      */
      virtual void
      readCheckpointHeader(HDF5HeaderData& a_header);
#endif

#ifdef CH_USE_HDF5
      /*!
	@brief Write header data to checkpoint file.
	@param[inout] a_header HDF5 header.
      */
      virtual void
      writeCheckpointHeader(HDF5HeaderData& a_header) const override;
#endif

#ifdef CH_USE_HDF5
      /*!  
	@brief Write checkpoint data to file
	@param[inout] a_handle HDF5 fil
	@param[in]    a_lvl    Grid level
      */
      virtual void
      writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
#endif

#ifdef CH_USE_HDF5
      /*!
	@brief Read checkpoint data from file
	@param[inout] a_handle HDF5 file
	@param[in]    a_lvl    Grid level
      */
      virtual void
      readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
#endif

      /*!
	@brief Write plot data to output holder. 
	@param[inout] a_output            Output data holder.
	@param[in]    a_plotVariableNames Plot variable names.
	@param[inout] a_icomp             Starting component in a_output to begin at. 
      */
      virtual void
      writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const override;

      /*!
	@brief Get number of plot variables for the output file
      */
      virtual int
      getNumberOfPlotVariables() const override;

      /*!
	@brief Get computational loads to be checkpointed. 
	@details This is used by Driver both for setting up load-balanced restarts AND for plotting the computational loads to a file. This routine is
	disjoint from loadBalanceBoxes because this routine is not part of a regrid. This means that we are not operating with temporarily load balanced
	grids, and all operators etc. are already in place.
	@param[in] a_realm Realm
	@param[in] a_level Grid level
	@return Returns computational loads for each box on grid level a_level. 
      */
      virtual Vector<long int>
      getCheckpointLoads(const std::string a_realm, const int a_level) const override;

      /*!
	@brief Advancement method. Needs to be implemented by subclasses.
	@param[in] a_dt Time step to be used for advancement
	@return    Returns the time step that was used. 
      */
      virtual Real
      advance(const Real a_dt) = 0;

      /*!
	@brief Compute a time step used for the advance method
      */
      virtual Real
      computeDt() override;

      /*!
	@brief Synchronize solver times for all the solvers
	@param[in] a_step Time step
	@param[in] a_time Simulation time
	@param[in] a_dt   Last time step used. 
      */
      virtual void
      synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;

      /*!
	@brief Print a step report. Used by Driver for user monitoring of simulation.
      */
      virtual void
      printStepReport() override;

      /*!
	@brief Register realms used for the simulation
      */
      virtual void
      registerRealms() override;

      /*!
	@brief Register operators used for the simulation
      */
      virtual void
      registerOperators() override;

      /*!
	@brief Perform pre-regrid operations.
	@details This stores the particles and mesh data on the old grids. 
	@param[in] a_lmin           The coarsest level that changes
	@param[in] a_oldFinestLevel The finest level before the regrid. 
      */
      virtual void
      preRegrid(const int a_lmin, const int a_oldFinestLevel) override;

      /*!
	@brief Regrid methods -- puts all data on the new mesh. 
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      virtual void
      regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /*!
	@brief Perform post-regrid operations -- default implementation does nothing. 
      */
      virtual void
      postRegrid() override;

      /*!
	@brief Load balancing query for a specified realm. If this returns true for a_realm, load balancing routines will be called during regrids for the input realm. 
	@param[in] a_realm Realm name
      */
      virtual bool
      loadBalanceThisRealm(const std::string a_realm) const override;

      /*!
	@brief Load balance grid boxes for a specified realm. 
	@param[out] a_procs       MPI ranks owning the various grid boxes. 
	@param[out] a_boxes       Grid boxes on every level (obtain them with a_grids[lvl].boxArray())
	@param[in]  a_realm       Realm identifier
	@param[in]  a_grids       Original grids
	@param[in]  a_lmin        Coarsest grid level that changed
	@param[in]  a_finestLevel New finest grid level
	@details This is only called by Driver if TimeStepper::loadBalanceThisRealm(a_realm) returned true. The default implementation
	uses volume-based loads for the grid patches. If the user wants to load balance boxes on a realm, this routine must be overwritten and
	he should compute loads for the various patches in a_grids and call LoadBalancing::makeBalance(Vector<int>, Vector<T>, Vector<Box>).
      */
      virtual void
      loadBalanceBoxes(Vector<Vector<int>>&             a_procs,
                       Vector<Vector<Box>>&             a_boxes,
                       const std::string                a_realm,
                       const Vector<DisjointBoxLayout>& a_grids,
                       const int                        a_lmin,
                       const int                        a_finestLevel) override;

      /*!
	@brief Set the Ito layout 
	@param[in] a_itoLayout Set of ItoSolver solvers. 
      */
      virtual void
      setIto(RefCountedPtr<ItoLayout<ItoSolver>>& a_itoLayout) noexcept;

      /*!
	@brief Set the field solver
	@param[in] a_fieldSolver Field solver
      */
      virtual void
      setFieldSolver(RefCountedPtr<FieldSolver>& a_fieldSolver) noexcept;

      /*!
	@brief Set the radiative transfer solvers
	@param[in] a_rteLayout Radiative transfer solvers
      */
      virtual void
      setRadiativeTransferSolvers(RefCountedPtr<RtLayout<McPhoto>>& a_rteLayout) noexcept;

      /*!
	@brief Set voltage used for the simulation
	@param[in] a_voltage Voltage
	@note This should be done BEFORE calling setupSolvers. 
      */
      virtual void
      setVoltage(const std::function<Real(const Real a_time)>& a_voltage) noexcept;

      /*!
	@brief Get current simulation time
	@return m_time
      */
      virtual Real
      getTime() const noexcept;

      /*!
	@brief Recompute the electric field onto the specified data holder
	@param[out] a_electricField Electric field on cell centroid (defined on the fluid realm). 
	@param[in]  a_phase Phase where we compute the electric field. 
      */
      virtual void
      computeElectricField(EBAMRCellData& a_electricField, const phase::which_phase a_phase) const noexcept;

    protected:
      /*!
	@brief An enum for encapsulating how time steps were restricted. 
      */
      enum class TimeCode
      {
        Advection,
        Diffusion,
        AdvectionDiffusion,
        RelaxationTime,
        Hardcap,
        Physics
      };

      /*!
	@brief Box sorting method when using dual-grid with particle load balancing.
      */
      BoxSorting m_boxSort;

      /*!
	@brief Time code for understanding how the time step was restricted. 
      */
      TimeCode m_timeCode;

      /*!
	@brief Realm used for the fluid part (i.e., electrostatic) part of the simulation
	@note This is hardcoded to Realm::Primal
      */
      std::string m_fluidRealm;

      /*!
	@brief Realm used for the particle part of the simulation. 
      */
      std::string m_particleRealm;

      /*!
	@brief Time stepper name
      */
      std::string m_name;

      /*!
	@brief Phase where we solve for the plasma. 
	@note This is hardcoded to phase::gas
      */
      phase::which_phase m_plasmaPhase;

      /*!
	@brief Implementation of ItoPlasmaPhysics. 
	@details All solver coupling occurs through the implementation of this class
      */
      RefCountedPtr<ItoPlasmaPhysics> m_physics;

      /*!
	@brief Ito solvers
      */
      RefCountedPtr<ItoLayout<ItoSolver>> m_ito;

      /*!
	@brief Radiative transfer solvers. 
      */
      RefCountedPtr<RtLayout<McPhoto>> m_rte;

      /*!
	@brief Field solver. 
      */
      RefCountedPtr<FieldSolver> m_fieldSolver;

      /*!
	@brief Surface charge solver. 
      */
      RefCountedPtr<SurfaceODESolver<1>> m_sigmaSolver;

      /*!
	@brief Voltage curve on the electrodes used in the simulation.
      */
      std::function<Real(const Real a_time)> m_voltage;

      /*!
	@brief Flag for abandoning simulation of Poisson solver fails
      */
      bool m_abortOnFailure;

      /*!
	@brief Using dual grid or not
      */
      bool m_dualGrid;

      /*!
	@brief Profile kernels or not
      */
      bool m_profile;

      /*!
	@brief Do or do not superparticle merging/splitting on regrids.
      */
      bool m_regridSuperparticles;

      /*!
	@brief Load balance particle realm or not
      */
      bool m_loadBalanceParticles;

      /*!
	@brief Load balance fluid realm or not
      */
      bool m_loadBalanceFluid;

      /*!
	@brief For switching between the old and new reaction algorithms
      */
      bool m_useNewReactionAlgorithm;

      /*!
	@brief Plot conductivity or not
      */
      bool m_plotConductivity;

      /*!
	@brief Plot conductivity or not
      */
      bool m_plotCurrentDensity;

      /*!
	@brief Plot number of particles per patch or not
      */
      bool m_plotParticlesPerPatch;

      /*!
	@brief Target number of particles per cell when squishing ItoParticle's into superparticles
      */
      int m_particlesPerCell;

      /*!
	@brief How often to merge superparticles
	@details If <= 1 we merge every time step. 
      */
      int m_mergeInterval;

      /*!
	@brief The "background" load per cell when using particle load balancing.
	@note This is used when computing the loads per patch on the particle realm such that the total load 
	is m_loadPerCell * numGridPoints + numParticles. 
      */
      Real m_loadPerCell;

      /*!
	@brief Accepted tolerance (relative to dx) for EB intersection
      */
      Real m_toleranceEB;

      /*!
	@brief CFL-like time step for advection
      */
      Real m_advectionCFL;

      /*!
	@brief CFL-like time step for diffusion
      */
      Real m_diffusionCFL;

      /*!
	@brief CFL-like time step for advection-diffusion
      */
      Real m_advectionDiffusionCFL;

      /*!
	@brief Factor proportional to the dielectric relaxation time dtRelax = eps0/sigma
	@details Used for restricting dt < m_relaxTimeFactor * eps0/sigma
      */
      Real m_relaxTimeFactor;

      /*!
	@brief Minimum permitted time step
      */
      Real m_minDt;

      /*!
	@brief Maximum permitted time step
      */
      Real m_maxDt;

      /*!
	@brief The advective time step (not multiplied by m_advectionCFL)
      */
      Real m_advectionDt;

      /*!
	@brief The diffusive time step (not multiplied by m_diffusionCFL)
      */
      Real m_diffusionDt;

      /*!
	@brief The advection-diffusion time step (not multiplied by m_advectionDiffusionCFL)
      */
      Real m_advectionDiffusionDt;

      /*!
	@brief The relaxation time eps0/sigma
      */
      Real m_relaxationTime;

      /*!
	@brief The physics-based time step
      */
      Real m_physicsDt;

      /*!
	@brief For holding the number of computational particles per cell when load balancing
      */
      Vector<EBAMRCellData> m_loadBalancePPC;

      /*!
	@brief Solver indices used when load-balancing the particle solvers.
      */
      Vector<int> m_loadBalanceIndices;

      /*!
	@brief Storage for current density.
	@note Defined over the fluid realm with SpaceDim components. 
      */
      EBAMRCellData m_currentDensity;

      /*!
	@brief Cell-centered conductivity.
	@note Defined on the fluid realm with 1 component
      */
      EBAMRCellData m_conductivityCell;

      /*!
	@brief Face-centered conductivity.
	@note Defined on the fluid realm with 1 component
      */
      EBAMRFluxData m_conductivityFace;

      /*!
	@brief EB-centered conductivity.
	@note Defined on the fluid realm with 1 component
      */
      EBAMRIVData m_conductivityEB;

      /*!
	@brief For holding the number of physical particles per cell
	@note Defined over the particle realm. 
      */
      EBAMRCellData m_oldPPC;

      /*!
	@brief For holding the number of physical particles per cell
	@note Defined over the particle realm. 
      */
      EBAMRCellData m_newPPC;

      /*!
	@brief For holding the number of physical particles per cell for all species.
	@note Defined on the particle realm with components = number of plasma species
      */
      EBAMRCellData m_particlePPC;

      /*!
	@brief For holding the number of generated photons per cell.
	@note Defined on the particle realm with components = number of photon species
      */
      EBAMRCellData m_particleYPC;

      /*!
	@brief For holding the mean particle energy
	@note Defined on the particle realm with components = number of plasma species
      */
      EBAMRCellData m_particleEPS;

      /*!
	@brief For holding the previous number of particles per cell for all species
	@note Defined on the particle realm with components = number of plasma species
      */
      EBAMRCellData m_particleOldPPC;

      /*!
	@brief For holding the number of particles per cell for all species.
	@note Defined on the fluid realm with components = number of plasma species
      */
      EBAMRCellData m_fluidPPC;

      /*!
	@brief For holding the number of generated photons per cell.
	@note Defined on the fluid realm with components = number of photon species
      */
      EBAMRCellData m_fluidYPC;

      /*!
	@brief For holding the mean particle energy
	@note Defined on the fluid realm with components = number of plasma species
      */
      EBAMRCellData m_fluidEPS;

      /*!
	@brief Storage for holding the plasma phase electric field on the fluid realm
      */
      EBAMRCellData m_electricFieldFluid;

      /*!
	@brief Storage for holding the plasma phase electric field on the particle realm
      */
      EBAMRCellData m_electricFieldParticle;

      /*!
	@brief Storage for holding the energy sources for each species.
	@note Defined on the particle realm but about to get deprecated when we switch to the new reaction algorithm. 
      */
      Vector<EBAMRCellData> m_energySources;

      /*!
	@brief Storage for holdnig E*J on the fluid realm. 1 component
      */
      EBAMRCellData m_EdotJ;

      /*!
	@brief Scratch storage on the fluid realm having 1 component
      */
      mutable EBAMRCellData m_fluidScratch1;

      /*!
	@brief Scratch storage on the fluid realm with SpaceDim components
      */
      mutable EBAMRCellData m_fluidScratchD;

      /*!
	@brief Scratch storage on the particle realm with 1 component
      */
      mutable EBAMRCellData m_particleScratch1;

      /*!
	@brief Scratch storage on the particle realm with SpaceDim components
      */
      mutable EBAMRCellData m_particleScratchD;

      /*!
	@brief Set up the Ito particle solvers
      */
      virtual void
      setupIto();

      /*!
	@brief Set up the electrostatic field solver
      */
      virtual void
      setupPoisson();

      /*!
	@brief Set up the radiative transfer solver
      */
      virtual void
      setupRadiativeTransfer();

      /*!
	@brief Set up the surface charge solver
      */
      virtual void
      setupSigma();

      /*!
	@brief Intersect a subset of the particles with the domain and embedded boundary. 
	@details This calls the other version, intersecting the Ito bulk particles with the EB and placing the collided particles
	in the appropriate particle containers (ItoSolver::WhichContainer::EB, for example). 
	@param[in] a_speciesSubset   Subset of particles to intersect. Can be charged, mobile, diffusive, etc. 
	@param[in] a_intersectionAlg Representation of the EB -- basically an algorithmic selection. 
	@param[in] a_delete          Delete original particles. If this is left to false, the particles are copied into the other containers. 
      */
      virtual void
      intersectParticles(const SpeciesSubset  a_speciesSubset,
                         const EBIntersection a_interactionAlg,
                         const bool           a_delete) noexcept;

      /*!
	@brief Intersect a subset of the particles with the domain and embedded boundary. 
	@param[in] a_speciesSubset   Subset of particles to intersect. Can be charged, mobile, diffusive, etc. 
	@param[in] a_containerBulk   ItoSolver container for "bulk" particles, i.e. particles that will be intersected. 
	@param[in] a_containerEB     ItoSolver container for EB particles. 
	@param[in] a_containerDomain ItoSolver container for domain particles. 
	@param[in] a_intersectionAlg EB intersection algorithm
	@param[in] a_delete          Delete original particles. If this is left to false, the particles are copied into the other containers. 
      */
      virtual void
      intersectParticles(const SpeciesSubset             a_speciesSubset,
                         const ItoSolver::WhichContainer a_containerBulk,
                         const ItoSolver::WhichContainer a_containerEB,
                         const ItoSolver::WhichContainer a_containerDomain,
                         const EBIntersection            a_intersectionAlg,
                         const bool                      a_delete) noexcept;

      /*!
	@brief Remove covered particles (i.e., particles inside the EB)
	@param[in] a_speciesSubset   Subset of particles to intersect. Can be charged, mobile, diffusive, etc. 
	@param[in] a_representation  Representation of the EB.
	@param[in] a_tolerance       Tolerance for removal. Removes up to distance a_tolerance from the EB.
	@note a_representation is basically an algorithmic selection that propagates into ItoSolver
      */
      virtual void
      removeCoveredParticles(const SpeciesSubset    a_which,
                             const EBRepresentation a_representation,
                             const Real             a_tolerance) noexcept;

      /*!
	@brief Remove covered particles (i.e., particles inside the EB)
	@param[in] a_speciesSubset   Subset of particles to intersect. Can be charged, mobile, diffusive, etc. 
	@param[in] a_container       Which of the containers in ItoSolver to remove from. 
	@param[in] a_representation  Representation of the EB.
	@param[in] a_tolerance       Tolerance for removal. Removes up to distance a_tolerance from the EB.
	@note a_representation is basically an algorithmic selection that propagates into ItoSolver
      */
      virtual void
      removeCoveredParticles(const SpeciesSubset             a_which,
                             const ItoSolver::WhichContainer a_container,
                             const EBRepresentation          a_representation,
                             const Real                      a_tolerance) noexcept;

      /*!
	@brief Transfer covered particles (i.e., particles inside the EB) from the ItoSolver bulk container to EB container. 
	@param[in] a_speciesSubset   Subset of particles to intersect. Can be charged, mobile, diffusive, etc. 
	@param[in] a_representation  Representation of the EB.
	@param[in] a_tolerance       Tolerance for transfer. Transfers up to distance a_tolerance from the EB.
	@note a_representation is basically an algorithmic selection that propagates into ItoSolver
      */
      virtual void
      transferCoveredParticles(const SpeciesSubset    a_speciesSubset,
                               const EBRepresentation a_representation,
                               const Real             a_tolerance) noexcept;

      /*!
	@brief Transfer covered particles (i.e., particles inside the EB) from the ItoSolver bulk container to EB container. 
	@param[in] a_speciesSubset   Subset of particles to intersect. Can be charged, mobile, diffusive, etc. 
	@param[in] a_containerFrom   Particle container (in ItoSolver) to transfer from
	@param[in] a_containerTo     Particle container (in ItoSolver) to transfer to
	@param[in] a_representation  Representation of the EB.
	@param[in] a_tolerance       Tolerance for transfer. Transfers up to distance a_tolerance from the EB.
	@note a_representation is basically an algorithmic selection that propagates into ItoSolver
      */
      virtual void
      transferCoveredParticles(const SpeciesSubset             a_speciesSubset,
                               const ItoSolver::WhichContainer a_containerFrom,
                               const ItoSolver::WhichContainer a_containerTo,
                               const EBRepresentation          a_representation,
                               const Real                      a_tolerance) noexcept;

      /*!
	@brief Write number of particles per patch to output holder. 
	@param[inout] a_output Output data holder.
	@param[inout] a_icomp  Starting component in a_output to begin at. 
      */
      virtual void
      writeNumParticlesPerPatch(EBAMRCellData& a_output, int& a_icomp) const;

      /*!
	@brief Get maximum density of the Ito species
	@param[out] a_maxDensity Maximum mesh density
	@param[out] a_minDensity Minium mesh density
	@param[out] a_maxSolver  Solver with highest density
	@param[out] a_minSolver  Solver with lowest density
      */
      virtual void
      getMaxMinItoDensity(Real&        a_maxDensity,
                          Real&        a_minDensity,
                          std::string& a_maxSolver,
                          std::string& a_minSolver) const noexcept;

      /*!
	@brief Compute some particle statistics
	@param[out] a_avgParticles Average numer of particles
	@param[out] a_sigma        Particle standard deviation (across MPI ranks)
	@param[out] a_minParticles Minimum number of particles
	@param[out] a_maxParticles Maximum number of particles
	@param[out] a_minRank      MPI rank with lowest number of particles
	@param[out] a_maxRank      MPI rank with highest number of particles
      */
      virtual void
      getParticleStatistics(Real& a_avgParticles,
                            Real& a_sigma,
                            Real& a_minParticles,
                            Real& a_maxParticles,
                            int&  a_minRank,
                            int&  a_maxRank);

      /*!
	@brief Routine called by loadBalanceBoxes and used for particle-based load balancing.
	@param[out] a_procs       MPI ranks owning the various grid boxes. 
	@param[out] a_boxes       Grid boxes on every level (obtain them with a_grids[lvl].boxArray())
	@param[in]  a_realm       Realm identifier
	@param[in]  a_grids       Original grids
	@param[in]  a_lmin        Coarsest grid level to load balance. 
	@param[in]  a_finestLevel Finest grid level
      */
      virtual void
      loadBalanceParticleRealm(Vector<Vector<int>>&             a_procs,
                               Vector<Vector<Box>>&             a_boxes,
                               const std::string                a_realm,
                               const Vector<DisjointBoxLayout>& a_grids,
                               const int                        a_lmin,
                               const int                        a_finestLevel) noexcept;

      /*!
	@brief Routine called by loadBalanceBoxes and used for particle-based load balancing.
	@param[out] a_procs       MPI ranks owning the various grid boxes. 
	@param[out] a_boxes       Grid boxes on every level (obtain them with a_grids[lvl].boxArray())
	@param[in]  a_realm       Realm identifier
	@param[in]  a_grids       Original grids
	@param[in]  a_lmin        Coarsest grid level to load balance. 
	@param[in]  a_finestLevel Finest grid level
      */
      virtual void
      loadBalanceFluidRealm(Vector<Vector<int>>&             a_procs,
                            Vector<Vector<Box>>&             a_boxes,
                            const std::string                a_realm,
                            const Vector<DisjointBoxLayout>& a_grids,
                            const int                        a_lmin,
                            const int                        a_finestLevel) noexcept;

      /*!
	@brief Get the solvers used for load balancing
	@details This depends on m_loadBalanceIndex -- if it's < 0 we use ALL plasma species for load balancing.
      */
      virtual Vector<RefCountedPtr<ItoSolver>>
      getLoadBalanceSolvers() const noexcept;

      /*!
	@brief Compute the maximum electric field (norm)
	@param[in] a_phase Phase where we compute the field. 
      */
      virtual Real
      computeMaxElectricField(const phase::which_phase a_phase) noexcept;

      /*!
	@brief Compute the space charge. Calls the other version. 
	@note This runs through the ItoSolver densities -- user should remember to deposit before calling this routine. 
      */
      virtual void
      computeSpaceChargeDensity() noexcept;

      /*!
	@brief Compute the space charge density (both phases)
	@param[out] a_rho   Space cahrge density
	@param[in]  a_phase Phase where we compute the space charge density      
      */
      virtual void
      computeSpaceChargeDensity(MFAMRCellData& a_rho, const Vector<EBAMRCellData*>& a_densities) noexcept;

      /*!
	@brief Compute the cell-centered conductiivty
	@param[out] a_conductivity Cell-centered conductivity
	@note Calls the other version. 
      */
      virtual void
      computeConductivityCell(EBAMRCellData& a_conductivity) noexcept;

      /*!
	@brief Compute the cell-centered conductiivty
	@param[out] a_conductivity Cell-centered conductivity
	@param[in]  a_particles    Ito particles.

      */
      virtual void
      computeConductivityCell(EBAMRCellData&                                 a_conductivity,
                              const Vector<ParticleContainer<ItoParticle>*>& a_particles) noexcept;

      /*!
	@brief Compute the current density. 
	@note This calls computeConductivityCell and uses m_electricFieldFluid. Make sure they're updated.
      */
      virtual void
      computeCurrentDensity(EBAMRCellData& a_J) noexcept;

      /*!
	@brief Compute the dielectric relaxation time
      */
      virtual Real
      computeRelaxationTime() noexcept;

      /*!
	@brief Solve the electrostatic problem
	@details This calls the other version with all lhs/rhs taken from the field solver. Space charge density is computed
	in this routine. 
	@return Returns true if the solver converged. 
      */
      virtual bool
      solvePoisson() noexcept;

      /*!
	@brief Deposit a subset of the ItoSolver particles on the mesh
	@note Calls the other version with ItoSolver::WhichContainer::Bulk
      */
      virtual void
      depositParticles(const SpeciesSubset a_speciesSubset) noexcept;

      /*!
	@brief Deposit a subset of the ItoSolver particles on the mesh
	@param[in] a_speciesSubset Subset of species to deposit
	@param[in] a_container     Which container in ItoSolver to deposit. 
      */
      virtual void
      depositParticles(const SpeciesSubset a_speciesSubset, const ItoSolver::WhichContainer a_container) noexcept;

      /*!
	@brief Remap a subset of ItoSolver particles
	@param[in] a_speciesSubset Subset of species to remap
	@note Calls the other version with ItoSolver::WhichContainer::Bulk
      */
      virtual void
      remapParticles(const SpeciesSubset a_speciesSubset) noexcept;

      /*!
	@brief Remap a subset of ItoSolver particles
	@param[in] a_speciesSubset Subset of species to remap
	@param[in] a_container     Which container in ItoSolver to remap.       
      */
      virtual void
      remapParticles(const SpeciesSubset a_speciesSubset, const ItoSolver::WhichContainer a_container) noexcept;

      /*!
	@brief Set the Ito velocity functions. This is sgn(charge) * E
	@note This should be set before ItoSolver computes velocities. 
      */
      virtual void
      setItoVelocityFunctions() noexcept;

      /*!
	@brief Compute ItoSolver velocities.
	@details This will compute the mobilities and then interpolate mu*E to particle positions.  
      */
      virtual void
      computeItoVelocities() noexcept;

      /*!
	@brief Compute ItoSolver velocities.
	@details This will compute the diffusion coefficients (in the specified coupling regime)
      */
      virtual void
      computeItoDiffusion() noexcept;

      /*!
	@brief Compute mesh-based mobilities for LFA coupling
	@note Calls the other version with mesh mobilities taken from ItoSolver solvers. 
      */
      virtual void
      computeItoMobilitiesLFA() noexcept;

      /*!
	@brief Compute mesh-based mobilities for LFA coupling
	@param[out] a_meshMobilities Mesh-based mobilities. Must be defined on the particle realm. 
	@param[in]  a_electricField  Electric field. Must be defined on the fluid realm. 
	@param[in]  a_time           Time
      */
      virtual void
      computeItoMobilitiesLFA(Vector<EBAMRCellData*>& a_meshMobilities,
                              const EBAMRCellData&    a_electricField,
                              const Real              a_time) noexcept;

      /*!
	@brief Level-routine when computing mesh-based mobilities for LFA coupling
	@param[out] a_meshMobilities Mesh-based mobilities. Must be defined on the fluid realm. 
	@param[in]  a_electricField  Electric field. Must be defined on the fluid realm. 
	@param[in]  a_level          Grid level
	@param[in]  a_time           Time
      */
      virtual void
      computeItoMobilitiesLFA(Vector<LevelData<EBCellFAB>*>& a_meshMobilities,
                              const LevelData<EBCellFAB>&    a_E,
                              const int                      a_level,
                              const Real                     a_time) noexcept;

      /*!
	@brief Patch-routine when computing mesh-based mobilities for LFA coupling
	@param[out] a_meshMobilities Mesh-based mobilities. Must be defined on the fluid realm. 
	@param[in]  a_electricField  Electric field. Must be defined on the fluid realm. 
	@param[in]  a_level          Grid level
	@param[in]  a_dit            Grid index
	@param[in]  a_cellBox        Computation box
	@param[in]  a_time           Time      
      */
      virtual void
      computeItoMobilitiesLFA(Vector<EBCellFAB*>& a_meshMobilities,
                              const EBCellFAB&    a_electricField,
                              const int           a_level,
                              const DataIndex     a_dit,
                              const Box           a_cellBox,
                              const Real          a_time) noexcept;

      /*!
	@brief Compute the IonSolver particle mobilities in the LEA coupling regime
	@note This one is very simple because the solvers just update the mobilities from appropriate coupling functions. 
      */
      virtual void
      computeItoMobilitiesLEA() noexcept;

      /*!
	@brief Compute mesh-based diffusion coefficients for LFA coupling
	@note Calls the other version with diffusion coefficients taken from ItoSolver solvers. 
      */
      virtual void
      computeItoDiffusionLFA() noexcept;

      /*!
	@brief Compute mesh-based diffusion coefficients for LFA coupling
	@param[out] a_diffusionCoefficients Mesh-based diffusion coefficients. Must be defined on the particle realm. 
	@param[in]  a_electricField         Electric field. Must be defined on the fluid realm. 
	@param[in]  a_time                  Time
      */
      virtual void
      computeItoDiffusionLFA(Vector<EBAMRCellData*>& a_diffusionCoefficient,
                             const EBAMRCellData&    a_electricField,
                             const Real              a_time) noexcept;

      /*!
	@brief Level-routine when computing mesh-based mobilities for LFA coupling
	@param[out] a_diffusionCoefficients Diffusion coefficients. Must be defined on the fluid realm. 
	@param[in]  a_electricField         Electric field. Must be defined on the fluid realm. 
	@param[in]  a_level                 Grid level
	@param[in]  a_time                  Time
      */
      virtual void
      computeItoDiffusionLFA(Vector<LevelData<EBCellFAB>*>& a_diffusionCoefficient,
                             const LevelData<EBCellFAB>&    a_electricField,
                             const int                      a_level,
                             const Real                     a_time) noexcept;

      /*!
	@brief Patch-routine when computing mesh-based diffusion coefficients in the LFA coupling
	@param[out] a_diffusionCoefficients Diffusion coefficients. Must be defined on the fluid realm. 
	@param[in]  a_electricField         Electric field. Must be defined on the fluid realm. 
	@param[in]  a_level                 Grid level
	@param[in]  a_time                  Time
	@param[in]  a_level                 Grid level
	@param[in]  a_dit                   Grid index
	@param[in]  a_cellBox               Computation box
	@param[in]  a_time                  Time      
      */
      virtual void
      computeItoDiffusionLFA(Vector<EBCellFAB*>& a_diffusionCoefficient,
                             const EBCellFAB&    a_E,
                             const int           a_level,
                             const DataIndex     a_dit,
                             const Box           a_box,
                             const Real          a_time) noexcept;

      /*!
	@brief Compute the IonSolver particle diffusion coefficients in the LEA coupling regime
	@note This one is very simple because the solvers just update the diffusion coefficients from appropriate coupling functions. 
      */
      virtual void
      computeItoDiffusionLEA() noexcept;

      /*!
	@brief Get the physical number of particles per cell.
	@param[out] a_ppc Number of physical particles per cell. Must be defined over the fluid realm.
       */
      virtual void
      getPhysicalParticlesPerCell(EBAMRCellData& a_ppc) const noexcept;

      /*!
	@brief Compute the number of reactive particles per cell
	@param[out] a_ppc Number of physical particles per grid cell for all species
	@note a_ppc Must be defined over the particle realm and have the same number of components as the number of plasma species involved in the simulation. 
      */
      virtual void
      computeReactiveParticlesPerCell(EBAMRCellData& a_ppc) noexcept;

      /*!
	@brief Compute the number of reactive particles per cell
	@param[out] a_ppc  Number of physical particles per grid cell for all species
	@param[in]  a_level Grid level
	@note a_ppc Must be defined over the particle realm and have the same number of components as the number of plasma species involved in the simulation. 
      */
      virtual void
      computeReactiveParticlesPerCell(LevelData<EBCellFAB>& a_ppc, const int a_level) noexcept;

      /*!
	@brief Compute the number of reactive particles per cell
	@param[out] a_ppc     Number of physical particles per grid cell for all species
	@param[in]  a_level   Grid level
	@param[in]  a_dit     Grid index
	@param[in]  a_box     Grid box
	@param[in]  a_ebisbox EB grid box
	@note a_ppc Must be defined over the particle realm and have the same number of components as the number of plasma species involved in the simulation. 
      */
      virtual void
      computeReactiveParticlesPerCell(EBCellFAB&      a_ppc,
                                      const int       a_level,
                                      const DataIndex a_dit,
                                      const Box       a_box,
                                      const EBISBox&  a_ebisbox) noexcept;

      /*!
	@brief Compute the mean particle energy in all grid cells.
	@param[out] a_meanEnergies Mean particle energies = sum(w * E)/sum(w) for all species. 
	@note a_meanEnergies must be defined over the particle realm and have the same number of components as the number of palsma species involved in the simulation. 
      */
      virtual void
      computeReactiveMeanEnergiesPerCell(EBAMRCellData& a_meanEnergies) noexcept;

      /*!
	@brief Compute the mean particle energy in all grid cells. Level version. 
	@param[out] a_meanEnergies Mean particle energies = sum(w * E)/sum(w) for all species. 
	@param[in]  a_level        Grid level
	@note a_meanEnergies must be defined over the particle realm and have the same number of components as the number of palsma species involved in the simulation. 
      */
      virtual void
      computeReactiveMeanEnergiesPerCell(LevelData<EBCellFAB>& a_meanEnergies, const int a_level) noexcept;

      /*!
	@brief Compute the mean particle energy in all grid cells. Patch version. 
	@param[out] a_meanEnergies Number of physical particles per grid cell for all species
	@param[in]  a_level        Grid level
	@param[in]  a_dit          Grid index
	@param[in]  a_box          Grid box
	@param[in]  a_ebisbox      EB grid box
	@note a_meanEnergies must be defined over the particle realm and have the same number of components as the number of palsma species involved in the simulation. 
      */
      virtual void
      computeReactiveMeanEnergiesPerCell(EBCellFAB&      a_meanEnergies,
                                         const int       a_level,
                                         const DataIndex a_dit,
                                         const Box       a_box,
                                         const EBISBox&  a_ebisbox) noexcept;

      /*!
	@brief Chemistry advance over time a_dt
	@param[in] a_dt Time step
	@note This will call the other version. The kernel will only run over valid grid cells (not covered by a finer grid. 
      */
      virtual void
      advanceReactionNetwork(const Real a_dt) noexcept;

      /*!
	@brief Chemistry advance over time a_dt. AMR version. 
	@param[in] a_electricField Electric field. Must be defined on the fluid realm. 
	@param[in] a_EdotJ         Source term EdotJ for each species. Must have same number of components as Ito layout and be defined on the fluid realm. 
	@param[in] a_dt            Time step
	@note This will call the level versions AND reconcileParticles. The kernel will only run over valid grid cells (not covered by a finer grid. 
      */
      virtual void
      advanceReactionNetwork(const EBAMRCellData& a_E, const EBAMRCellData& a_EdotJ, const Real a_dt) noexcept;

      /*!
	@brief Chemistry advance over time a_dt. Level version. 
	@param[inout] a_particlesPerCell      Number of particles per cell for each plasma species
	@param[inout] a_newPhotonsPerCell     Number of physical photons to be generated for each photon species
	@param[inout] a_meanParticleEnergies  Mean particle energies for each plasma species
	@param[in]    a_electricField         Electric field. 
	@param[in]    a_EdotJ                 EdotJ for each plasma species.
	@param[in]    a_level                 Grid level
	@param[in]    a_dt                    Time increment
	@note This should be called through the AMR signature. All data should be defined on the fluid realm. The kernel will
	only run over valid grid cells (not covered by a finer grid. 
      */
      virtual void
      advanceReactionNetwork(LevelData<EBCellFAB>&       a_particlesPerCell,
                             LevelData<EBCellFAB>&       a_newPhotonsPerCell,
                             LevelData<EBCellFAB>&       a_meanParticleEnergies,
                             const LevelData<EBCellFAB>& a_electricField,
                             const LevelData<EBCellFAB>& a_EdotJ,
                             const int                   a_level,
                             const Real                  a_dt) noexcept;

      /*!
	@brief Chemistry advance over time a_dt. Patch version. 
	@param[inout] a_particlesPerCell      Number of particles per cell for each plasma species
	@param[inout] a_newPhotonsPerCell     Number of physical photons to be generated for each photon species
	@param[inout] a_meanParticleEnergies  Mean particle energies for each plasma species
	@param[in]    a_electricField         Electric field. 
	@param[in]    a_EdotJ                 EdotJ for each plasma species.
	@param[in]    a_level                 Grid level
	@param[in]    a_dit                   Grid index.
	@param[in]    a_box                   Grid box
	@param[in]    a_dx                    Grid resolution.
	@param[in]    a_dt                    Time increment
	@note This should be called through the AMR signature. All data should be defined on the fluid realm. The kernel will
	only run over valid grid cells (not covered by a finer grid. 
      */
      virtual void
      advanceReactionNetwork(EBCellFAB&       a_particlesPerCell,
                             EBCellFAB&       a_newPhotonsPerCell,
                             EBCellFAB&       a_meanParticleEnergies,
                             const EBCellFAB& a_electricField,
                             const EBCellFAB& a_EdotJ,
                             const int        a_level,
                             const DataIndex  a_dit,
                             const Box        a_box,
                             const Real       a_dx,
                             const Real       a_dt) noexcept;

      /*!
	@brief Reconcile particles. At the bottom, this will call the physics interface for particle reconciliation
	@param[in] a_newParticlesPerCell  Particles per cell after the chemistry advance. 
	@param[in] a_oldParticlesPerCell  Particles per cell before the chemistry advance. 
	@param[in] a_meanParticleEnergies Mean particle energies after the chemistry advance. 
	@param[in] a_newPhotonsPerCell    Newly generated photons per cell after the chemistry avance. 
	@note All data should be defined on the particle realm since this is where the particle reconciliation takes place. 
      */
      virtual void
      reconcileParticles(const EBAMRCellData& a_newParticlesPerCell,
                         const EBAMRCellData& a_oldParticlesPerCell,
                         const EBAMRCellData& a_meanParticleEnergies,
                         const EBAMRCellData& a_newPhotonsPerCell) noexcept;

      /*!
	@brief Reconcile particles. At the bottom, this will call the physics interface for particle reconciliation
	@param[in] a_newParticlesPerCell  Particles per cell after the chemistry advance. 
	@param[in] a_oldParticlesPerCell  Particles per cell before the chemistry advance. 
	@param[in] a_meanParticleEnergies Mean particle energies after the chemistry advance. 
	@param[in] a_newPhotonsPerCell    Newly generated photons per cell after the chemistry avance. 
	@param[in] a_level                Grid level
	@note All data should be defined on the particle realm since this is where the particle reconciliation takes place. 
      */
      virtual void
      reconcileParticles(const LevelData<EBCellFAB>& a_newParticlesPerCell,
                         const LevelData<EBCellFAB>& a_oldParticlesPerCell,
                         const LevelData<EBCellFAB>& a_meanParticleEnergies,
                         const LevelData<EBCellFAB>& a_newPhotonsPerCell,
                         const int                   a_level) noexcept;

      /*!
	@brief Reconcile particles. At the bottom, this will call the physics interface for particle reconciliation
	@param[in] a_newParticlesPerCell  Particles per cell after the chemistry advance. 
	@param[in] a_oldParticlesPerCell  Particles per cell before the chemistry advance. 
	@param[in] a_meanParticleEnergies Mean particle energies after the chemistry advance. 
	@param[in] a_newPhotonsPerCell    Newly generated photons per cell after the chemistry avance. 
	@param[in] a_level                Grid level
	@param[in] a_dit                  Grid index
	@param[in] a_box                  Grid box
	@param[in] a_dx                   Grid resolution. 
	@note All data should be defined on the particle realm since this is where the particle reconciliation takes place. 
      */
      virtual void
      reconcileParticles(const EBCellFAB& a_newParticlesPerCell,
                         const EBCellFAB& a_oldParticlesPerCell,
                         const EBCellFAB& a_meanParticleEnergies,
                         const EBCellFAB& a_newPhotonsPerCell,
                         const int        a_level,
                         const DataIndex  a_dit,
                         const Box        a_box,
                         const Real       a_dx) noexcept;

      /*!
	@brief Resolve particle injection at EBs.
	@param[in] a_dt Time step.
	@note Calls the other version.
      */
      virtual void
      resolveParticlesEB(const Real a_dt) noexcept;

      /*!
	@brief Resolve particle injection at EBs.
	@param[in] a_dt Time step
      */
      virtual void
      resolveParticlesEB(Vector<ParticleContainer<ItoParticle>*>& a_inOutParticles,
                         Vector<ParticleContainer<Photon>*>&      a_inOutPhotons,
                         const EBAMRCellData&                     a_electricField,
                         const Real                               a_dt) noexcept;

      /*!
	@brief Inject particles through the EB
	@note Calls the other version with solver containers. 
      */
      virtual void
      injectParticlesEB() noexcept;

      /*!
	@brief Inject particles through the EB
      */
      virtual void
      injectParticlesEB(Vector<ParticleContainer<ItoParticle>*> a_bulkParticles,
                        Vector<ParticleContainer<ItoParticle>*> a_injectedParticles) noexcept;

      /*!
	@brief Compute a maximum time step from the physics interface.
	@note Calls the other version with m_electricFieldParticle and the ItoSolver densities. The kernels are only
	called on the valid region.
      */
      virtual Real
      computePhysicsDt() noexcept;

      /*!
	@brief Compute a maximum time step from the physics interface.
	@param[in] a_electricField Electric field. 
	@param[in] a_densities     ItoSolver mesh densities. 
	@note Input data must be defined on the fluid realms. Kernels are only called on the valid region. 
      */
      virtual Real
      computePhysicsDt(const EBAMRCellData& a_electricField) noexcept;

      /*!
	@brief Compute a maximum time step from the physics interface.
	@param[in] a_electricField    Electric field. 
	@param[in] a_particlesPerCell Number of particles per grid cell
	@param[in] a_level            Grid level
	@note Input data must be defined on the fluid realms. Kernels are only called on the valid region. 
      */
      virtual Real
      computePhysicsDt(const LevelData<EBCellFAB>& a_electricField,
                       const LevelData<EBCellFAB>& a_particlesPerCell,
                       const int                   a_level) noexcept;

      /*!
	@brief Compute a maximum time step from the physics interface.
	@param[in] a_electricField    Electric field. 
	@param[in] a_particlesPerCell Number of particles per grid cell
	@param[in] a_level            Grid level
	@param[in] a_dit              Grid index
	@param[in] a_box              Grid box
	@note Input data must be defined on the fluid realms. Kernels are only called on the valid region. 
      */
      virtual Real
      computePhysicsDt(const EBCellFAB& a_electricField,
                       const EBCellFAB& a_particlesPercCell,
                       const int        a_level,
                       const DataIndex  a_dit,
                       const Box        a_box) noexcept;

      /*!
	@brief Photon advancement routine
	@param[in] a_dt Time step
      */
      virtual void
      advancePhotons(const Real a_dt) noexcept;

      /*!
	@brief Sort photons by cells
	@param[in] a_which Which container to sort by cell
      */
      virtual void
      sortPhotonsByCell(const McPhoto::WhichContainer a_which) noexcept;

      /*!
	@brief Sort photons by patch
	@param[in] a_which Which container to sort by patch
      */
      virtual void
      sortPhotonsByPatch(const McPhoto::WhichContainer a_which) noexcept;

      /*!
	@brief Compute the energy source term for the various plasma species
	@param[in] a_dt Time step during which the particles move. 
	@note This fills m_EdotJ
      */
      virtual void
      computeEdotJSource(const Real a_dt) noexcept;

      /*!
	@brief Fill surface charge solver with initial data
      */
      virtual void
      initialSigma();

      /*!
	@brief Parse chattiness
	@details Sets m_verbosity from input script. 
      */
      virtual void
      parseVerbosity() noexcept;

      /*!
	@brief Parse abort on failure
      */
      virtual void
      parseAbort() noexcept;

      /*!
	@brief Parse the desired number of particles per cell. 
	@details Sets m_particlesPerCell from input script. 
      */
      virtual void
      parseSuperParticles() noexcept;

      /*!
	@brief Parse dual or single realm calculations
	@details Sets m_particleRealm = "ParticleRealm" if 'dual_grid' is true.
      */
      virtual void
      parseDualGrid() noexcept;

      /*!
	@brief Parse load balancing
      */
      virtual void
      parseLoadBalance() noexcept;

      /*!
	@brief Parse time step restrictions
      */
      virtual void
      parseTimeStepRestrictions() noexcept;

      /*!
	@brief Parse things related to the EB
      */
      virtual void
      parseParametersEB() noexcept;

      /*!
	@brief Parse plot variables
      */
      virtual void
      parsePlotVariables() noexcept;
    };
  } // namespace ItoPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
