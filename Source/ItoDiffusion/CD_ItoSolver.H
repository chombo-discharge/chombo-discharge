/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ItoSolver.H
 * @brief  Declaration of solver class for Ito diffusion
 * @author Robert Marskar
 */

#ifndef CD_ITOSOLVER_H
#define CD_ITOSOLVER_H

// Std includes
#include <optional>

// Our includes
#include <CD_AmrMesh.H>
#include <CD_ComputationalGeometry.H>
#include <CD_ItoSpecies.H>
#include <CD_ItoParticle.H>
#include <CD_ItoMergeParticle.H>
#include <CD_ParticleContainer.H>
#include <CD_EBRepresentation.H>
#include <CD_EBIntersection.H>
#include <CD_CellInfo.H>
#include <CD_ParticleManagement.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Base class for Ito diffusion particle models.
 * @details This is a particle class for solving Ito diffusion problems dX = v*dt + sqrt(2*D)*dW_t (dW_t is a Wiener
 * process over dt) over an AMR hierarchy. The diffusion coefficient is related to the hydrodynamic diffusion
 * coefficient by D(fluid) = sqrt(2*D(Ito)).
 * @note To use this class, the user must add a species and set the realm and phase where the solver will live.
 */
class ItoSolver
{
public:
  /**
   * @brief Enum class for distinguishing various types of particle containers.
   * @details This exists because the ItoSolver can partition particles into various containers, which is very useful
   * when one wants to add particles from a source term, remove particles that fall inside the EB, or parse boundary
   * conditions on the EB and domain faces. Here, Bulk    = "Active" particles, EB      = Particles on the EBs , Domain
   * = Particles on the domain sides, Source  = Source particles -- used to e.g. add new particles to the bulk
   * particles, Covered = Particles inside the EB, Scratch = Scratch particles
   */
  enum class WhichContainer
  {
    Bulk,
    EB,
    Domain,
    Source,
    Covered,
    Scratch
  };

  /**
   * @brief Constructor -- user must subsequently set the realm and, parse class options, set the species etc.
   */
  ItoSolver();

  /**
   * @brief Destructor (does nothing).
   */
  virtual ~ItoSolver();

  /**
   * @brief Deleted copy constructor.
   */
  ItoSolver(const ItoSolver&) = delete;

  /**
   * @brief Deleted copy assignment operator.
   */
  ItoSolver&
  operator=(const ItoSolver&) = delete;

  /**
   * @brief Default move constructor.
   */
  ItoSolver(ItoSolver&&) = default;

  /**
   * @brief Default move assignment operator.
   * @return Reference to this object
   */
  ItoSolver&
  operator=(ItoSolver&&) = default;

  /**
   * @brief Set the user-supplied per-cell particle merger used by merge_algorithm = external.
   * @details The built-in merge methods build their own per-cell merger internally; this one is used
   * only when merge_algorithm = external, applied cell-by-cell by makeSuperparticles().
   * @note The merger runs on ItoMergeParticle, not ItoParticle: every merge in this solver operates on
   * the reduced particle, and a per-cell merger sees exactly the columns a merge is entitled to change
   * (position, weight, energy). Out-of-tree mergers written against ItoParticle break at compile time.
   * @param[in] a_particleCellMerger Per-cell particle merger.
   */
  virtual void
  setParticleCellMerger(const ParticleManagement::ParticleMerger<ItoMergeParticle>& a_particleCellMerger) noexcept;

  /**
   * @brief Get this solver's name
   * @return Returns m_name
   */
  virtual std::string
  getName() const;

  /**
   * @brief Get the realm where this solver is registered.
   * @return m_realm
   */
  virtual std::string
  getRealm() const;

  /**
   * @brief Set the realm where this solver will live
   * @param[in] a_realm Realm name.
   */
  virtual void
  setRealm(const std::string& a_realm);

  /**
   * @brief Parse class options
   */
  virtual void
  parseOptions();

  /**
   * @brief Parse class runtime configurable options.
   */
  virtual void
  parseRuntimeOptions();

  /**
   * @brief Get the EB-particle intersection algorithm
   * @return Intersection algorithm
   */
  virtual EBIntersection
  getIntersectionAlgorithm() const noexcept;

  /**
   * @brief Set the initial data
   * @details This will add the initial particles and deposit them.
   */
  virtual void
  initialData();

  /**
   * @brief Fill a particle container randomly with particles such that we obtain the target uniform
   * density, with a per-level cap on the computational particle count.
   * @details The per-level form of the overload below. Used by initialData(), which seeds at
   * ItoSolver.particles_per_cell -- a per-level vector -- so that a cold start and the first merge agree
   * on the target by construction instead of seeding high and immediately merging back down.
   * @param[in,out] a_particles         Particle container to fill. Particles already in the container will be kept.
   * @param[in]     a_densityFunc       Density function.
   * @param[in]     a_maxParticlesPerCell Maximum computational particles generated per cell, per level;
   *                                    the last entry applies to any finer level.
   */
  virtual void
  generateParticlesFromDensity(ParticleContainer<ItoParticle>&              a_particles,
                               const std::function<Real(const RealVect x)>& a_densityFunc,
                               const Vector<int>&                           a_maxParticlesPerCell) const noexcept;

  /**
   * @brief Fill a particle container randomly with particles such that we obtain the target uniform density.
   * @param[in,out] a_particles Particle container to fill. Particles already in the container will be kept.
   * @param[in] a_densityFunc Density function.
   * @param[in] a_maxParticlesPerCell Maximum number of computational particles generated per cell.
   */
  virtual void
  generateParticlesFromDensity(ParticleContainer<ItoParticle>&              a_particles,
                               const std::function<Real(const RealVect x)>& a_densityFunc,
                               int                                          a_maxParticlesPerCell) const noexcept;

  /**
   * @brief Regrid this solver.
   * @param[in] a_lmin           Coarsest level where grids did not change.
   * @param[in] a_oldFinestLevel Finest AMR level before the regrid.
   * @param[in] a_newFinestLevel Finest AMR level after the regrid.
   */
  virtual void
  regrid(int a_lmin, int a_oldFinestLevel, int a_newFinestLevel);

  /**
   * @brief Allocate internal storage
   * @details This will allocate the required mesh data as well as all the particle data holders.
   * @note Storage for diffusion and advection fields are only allocate if the species is diffusive/mobile.
   */
  virtual void
  allocate();

  /**
   * @brief Clear a particle container -- this will delete all the particles in the input container.
   * @param[in] a_container Container
   */
  virtual void
  clear(WhichContainer a_container);

  /**
   * @brief Clear a particle container -- this will delete all the particles in the input container.
   * @param[in] a_particles Particle container to be emptied.
   */
  virtual void
  clear(ParticleContainer<ItoParticle>& a_particles) const;

  /**
   * @brief Deposit conductivities (i.e. mass*mobility / volume)
   * @details This deposits mass*mobility (not multiplied by charge)
   * @param[out] a_phi       Mesh data
   * @param[in]  a_particles Particle data
   */
  virtual void
  depositConductivity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const;

  /**
   * @brief Deposit conductivities (i.e. mass*mobility / volume)
   * @details This deposits mass*mobility (not multiplied by charge)
   * @param[out] a_phi                  Mesh data
   * @param[in]  a_particles            Particle data
   * @param[in]  a_deposition           Deposition method
   * @param[in]  a_coarseFineDeposition Coarse-fine deposition method.
   */
  virtual void
  depositConductivity(EBAMRCellData&                  a_phi,
                      ParticleContainer<ItoParticle>& a_particles,
                      DepositionType                  a_deposition,
                      CoarseFineDeposition            a_coarseFineDeposition) const;

  /**
   * @brief Deposit diffusivity (i.e. mass*D/volume)
   * @details This deposits mass*diffusion (not multiplied by charge)
   * @param[out] a_phi       Mesh data
   * @param[in]  a_particles Particle data
   * @note Calls the other versions with a_deposition = m_deposition
   */
  virtual void
  depositDiffusivity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const;

  /**
   * @brief Deposit diffusivity (i.e. mass*D/volume)
   * @details This deposits mass*mobility (not multiplied by charge)
   * @param[out] a_phi                  Mesh data
   * @param[in]  a_particles            Particle data
   * @param[in]  a_deposition           Deposition method
   * @param[in]  a_coarseFineDeposition Coarse-fine deposition method.
   */
  virtual void
  depositDiffusivity(EBAMRCellData&                  a_phi,
                     ParticleContainer<ItoParticle>& a_particles,
                     DepositionType                  a_deposition,
                     CoarseFineDeposition            a_coarseFineDeposition) const;

  /**
   * @brief Deposit energy densities (i.e. mass*energy/volume => total energy per unit volume)
   * @details This deposits mass*energy
   * @param[out] a_phi       Mesh data
   * @param[in]  a_particles Particle data
   * @note Calls the other versions with a_deposition = m_deposition
   */
  virtual void
  depositEnergyDensity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const;

  /**
   * @brief Deposit energy densities (i.e. mass*energy/volume => total energy per unit volume)
   * @details This deposits mass*energy
   * @param[out] a_phi                  Mesh data
   * @param[in]  a_particles            Particle data
   * @param[in]  a_deposition           Deposition method
   * @param[in]  a_coarseFineDeposition Coarse-fine deposition method.
   */
  virtual void
  depositEnergyDensity(EBAMRCellData&                  a_phi,
                       ParticleContainer<ItoParticle>& a_particles,
                       DepositionType                  a_deposition,
                       CoarseFineDeposition            a_coarseFineDeposition) const;

  /**
   * @brief Compute average mobility.
   * @details This computes the average mobility as mu_avg = sum(mass*mu)/sum(mass) by depositing on the mesh.
   * @param[out] a_phi       Mesh data
   * @param[in]  a_particles Particle data
   */
  virtual void
  computeAverageMobility(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const;

  /**
   * @brief Compute average diffusion coefficient.
   * @details This computes the average diffusion coefficient as D_avg = sum(mass*D)/sum(mass) by depositing on the
   * mesh.
   * @param[out] a_phi       Mesh data
   * @param[in]  a_particles Particle data
   */
  virtual void
  computeAverageDiffusion(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const;

  /**
   * @brief Compute average particle energy..
   * @details This computes the average energy as E_avg = sum(mass*energy)/sum(mass) by depositing on the mesh.
   * @param[out] a_phi       Mesh data
   * @param[in]  a_particles Particle data
   */
  virtual void
  computeAverageEnergy(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const;

  /**
   * @brief Deposit particles onto mesh.
   * @details This will deposit the mass (i.e., computational weight) "bulk" particles into m_phi.
   * @note Calls the other version with a_container = WhichContainer::Bulk
   */
  virtual void
  depositParticles();

  /**
   * @brief Deposit particles on to mesh.
   * @param[in] a_container Which container to deposit.
   * @details This will deposit mass (i.e., computational weight) of the the input particle container particles onto the
   * classes member 'm_phi'.
   * @note Calls the general version with arguments: m_phi, m_particles.at(a_container), m_deposition.
   */
  virtual void
  depositParticles(WhichContainer a_container);

  /**
   * @brief Do an NGP deposit of the SoA weight column on a specific grid level. Used for IO.
   * @param[out] a_output    Contains NGP deposition of particle weights. Ignores cut-cells.
   * @param[in]  a_particles SoA particles.
   * @param[in]  a_level     Grid level.
   */
  void
  depositWeightNGP(LevelData<EBCellFAB>&                 a_output,
                   const ParticleContainer<ItoParticle>& a_particles,
                   int                                   a_level) const noexcept;

  /**
   * @brief Do an NGP deposit of a gathered per-particle quantity on a specific grid level. Used for IO.
   * @param[out] a_output    Contains NGP deposition of the gathered quantity. Ignores cut-cells.
   * @param[in]  a_particles SoA particles.
   * @param[in]  a_level     Grid level.
   * @param[in]  a_gather    Per-particle value gatherer (leaf, index) -> Real.
   * @tparam     Gather      Callable (const ParticleSoA<ItoParticle>&, std::size_t) -> Real.
   */
  template <typename Gather>
  void
  depositGatheredNGP(LevelData<EBCellFAB>&                 a_output,
                     const ParticleContainer<ItoParticle>& a_particles,
                     int                                   a_level,
                     Gather                                a_gather) const noexcept;

  /**
   * @brief Coarsen the input data and interpolate its ghost cells.
   * @details This is the AMR synchronization that a deposit does *not* do for itself. Whether it is wanted depends on
   * what the caller does next: a deposit that is the final state of a field needs it, whereas a deposit that is one
   * term in a sum does not -- there, only the assembled sum needs synchronizing, and doing it per term produces values
   * that are immediately overwritten. Deposits therefore leave it to the caller.
   * @param[in,out] a_phi Cell-centered mesh data (one component).
   */
  void
  coarsenAndFillGhosts(EBAMRCellData& a_phi) const;

  /**
   * @brief Deposit the SoA weight column on the mesh (kappa-conservative + redistribution).
   * @details The container payload is a template parameter so that auxiliary point-particle containers -- e.g. the
   * rho-dagger and conductivity particles carried by the semi-implicit time steppers -- deposit through exactly the
   * same pipeline as the solver's own particles. This includes the cut-cell deposition flag and the redistribution,
   * both of which must be identical if the deposited densities are to be in the same normalization.
   * @note This leaves coarse levels un-averaged and ghost cells stale. Call coarsenAndFillGhosts() afterwards unless
   * the result is a term in a sum that the caller synchronizes itself.
   * @param[out] a_phi                  Mesh data -- must have exactly one component.
   * @param[in]  a_particles            SoA particles to be deposited.
   * @param[in]  a_deposition           Deposition method.
   * @param[in]  a_coarseFineDeposition Coarse-fine deposition strategy.
   * @tparam     P                      Payload type of the container.
   * @tparam     Traits                 Column descriptor for the payload.
   */
  template <typename P, typename Traits>
  void
  depositWeight(EBAMRCellData&                      a_phi,
                const ParticleContainer<P, Traits>& a_particles,
                DepositionType                      a_deposition,
                CoarseFineDeposition                a_coarseFineDeposition) const;

  /**
   * @brief Deposit a gathered per-particle quantity on the mesh (kappa-conservative + redistribution).
   * @param[out] a_phi                  Mesh data -- must have exactly one component.
   * @param[in]  a_particles            SoA particles to be deposited.
   * @param[in]  a_deposition           Deposition method.
   * @param[in]  a_coarseFineDeposition Coarse-fine deposition strategy.
   * @param[in]  a_gather               Per-particle value gatherer (leaf, index) -> Real.
   * @tparam     Gather                 Callable (const ParticleSoA<ItoParticle>&, std::size_t) -> Real.
   * @note This leaves coarse levels un-averaged and ghost cells stale. Call coarsenAndFillGhosts() afterwards unless
   * the result is a term in a sum that the caller synchronizes itself.
   */
  template <typename Gather>
  void
  depositGathered(EBAMRCellData&                        a_phi,
                  const ParticleContainer<ItoParticle>& a_particles,
                  DepositionType                        a_deposition,
                  CoarseFineDeposition                  a_coarseFineDeposition,
                  Gather                                a_gather) const;

  /**
   * @brief Remove particles that are inside the EB.
   * @param[in] a_tol              Tolerance -- removes particles if they are less than a_tol*dx away from the boundary.
   * @details Calls the other version with m_particles.
   * @param[in] a_representation Representation
   */
  virtual void
  removeCoveredParticles(EBRepresentation a_representation, Real a_tol);

  /**
   * @brief Remove particles that are inside the EB.
   * @param[in] a_container        Which particle container to remove particles from
   * @param[in] a_tol              Tolerance -- removes particles if they are less than a_tol*dx away from the boundary.
   * @note Calls the other version with ParticleContainer = m_particles.at(a_container)
   * @param[in] a_representation Representation
   */
  virtual void
  removeCoveredParticles(WhichContainer a_container, EBRepresentation a_representation, Real a_tol);

  /**
   * @brief Remove particles that are inside the EB
   * @param[in] a_particles        Particles to remove.
   * @param[in] a_tol              Tolerance -- removes particles if they are less than a_tol*dx away from the boundary.
   * @param[in] a_representation Representation
   */
  virtual void
  removeCoveredParticles(ParticleContainer<ItoParticle>& a_particles,
                         EBRepresentation                a_representation,
                         Real                            a_tol) const;

  /**
   * @brief Transfer particles that are covered (within a_tol distance from EB) to another container.
   * @details Calls the other version with m_particles and m_coveredParticles
   * @param[in] a_tol              Tolerance -- removes particles if they are less than a_tol*dx away from the boundary.
   * @param[in] a_representation Representation
   */
  virtual void
  transferCoveredParticles(EBRepresentation a_representation, Real a_tol);

  /**
   * @brief Transfer particles that are covered (within a_tol distance from EB) to another container.
   * @param[in] a_containerFrom    Which particle container to transfer particles from
   * @param[in] a_tol              Tolerance -- removes particles if they are less than a_tol*dx away from the boundary.
   * @param[in] a_containerTo Container to
   * @param[in] a_representation Representation
   */
  virtual void
  transferCoveredParticles(WhichContainer   a_containerFrom,
                           WhichContainer   a_containerTo,
                           EBRepresentation a_representation,
                           Real             a_tol);

  /**
   * @brief Transfer particles that are inside the EB to another container.
   * @param[in] a_particlesFrom    Which particle container to transfer particles from
   * @param[in] a_particlesTo      Which particle container to move particles into
   * @param[in] a_tol              Tolerance -- removes particles if they are less than a_tol*dx away from the boundary.
   * @param[in] a_representation Representation
   */
  virtual void
  transferCoveredParticles(ParticleContainer<ItoParticle>& a_particlesFrom,
                           ParticleContainer<ItoParticle>& a_particlesTo,
                           EBRepresentation                a_representation,
                           Real                            a_tol) const;

  /**
   * @brief Do boundary intersection tests.
   * @details This will intersect the particles in the "bulk" particles data holder with the domain faces and EBs. If a
   * particle crossed the EB it will be put into the "EB" particle data holder and likewise for the particles that
   * crossed the domain side.
   * @param[in] a_ebIntersection      Enum for switching between various types of intersection tests.
   * @param[in] a_deleteParticles     If true, the origin particle will also be removed from the bulk particle data
   * holder.
   * @param[in] a_nonDeletionModifier Optional input argument for letting the user manipulate particles that were
   * intersected but not deleted
   * @note This will call the other version with WhichContainer::Bulk, WhichContainer::EB, and WhichContainer::Domain.
   */
  virtual void
  intersectParticles(
    const EBIntersection                                               a_ebIntersection,
    const bool                                                         a_deleteParticles,
    const std::function<void(ParticleSoA<ItoParticle>&, std::size_t)>& a_nonDeletionModifier =
      [](ParticleSoA<ItoParticle>&, std::size_t) -> void {
      return;
    });

  /**
   * @brief Do boundary intersection tests.
   * @details This will intersect the particles in the a_particles data holder with the domain faces and EBs. If a
   * particle crossed the EB it will be put into the a_eb_particles particle data holder and likewise for the particles
   * that crossed the domain side.
   * @param[in] a_particles           Which particle container to use for the "active" particles
   * @param[in] a_ebParticles         Which particle container to put the EB particles
   * @param[in] a_domainParticles     Which particle container to put the domain particles
   * @param[in] a_ebIntersection      Enum for switching between various types of intersection tests.
   * @param[in] a_deleteParticles     If true, the origin particle will also be removed from the bulk particle data
   * holder.
   * @param[in] a_nonDeletionModifier Optional input argument for letting the user manipulate particles that were
   * intersected but not deleted
   * @note This will look up the ParticleContainer and call the other version.
   */
  virtual void
  intersectParticles(
    const WhichContainer                                               a_particles,
    const WhichContainer                                               a_ebParticles,
    const WhichContainer                                               a_domainParticles,
    const EBIntersection                                               a_ebIntersection,
    const bool                                                         a_deleteParticles,
    const std::function<void(ParticleSoA<ItoParticle>&, std::size_t)>& a_nonDeletionModifier =
      [](ParticleSoA<ItoParticle>&, std::size_t) -> void {
      return;
    });

  /**
   * @brief Do boundary intersection tests.
   * @details This will intersect the particles in the a_particles data holder with the domain faces and EBs. If a
   * particle crossed the EB it will be put into the a_eb_particles particle data holder and likewise for the particles
   * that crossed the domain side.
   * @param[in] a_particles           Which particle container to use for the "active" particles
   * @param[in] a_ebParticles         Which particle container to put the EB particles
   * @param[in] a_domainParticles     Which particle container to put the domain particles
   * @param[in] a_ebIntersection      Enum for switching between various types of intersection tests.
   * @param[in] a_deleteParticles     If true, the origin particle will also be removed from the bulk particle data
   * holder.
   * @param[in] a_nonDeletionModifier Optional input argument for letting the user manipulate particles that were
   * intersected but not deleted
   */
  virtual void
  intersectParticles(
    ParticleContainer<ItoParticle>&                                    a_particles,
    ParticleContainer<ItoParticle>&                                    a_ebParticles,
    ParticleContainer<ItoParticle>&                                    a_domainParticles,
    const EBIntersection                                               a_ebIntersection,
    const bool                                                         a_deleteParticles,
    const std::function<void(ParticleSoA<ItoParticle>&, std::size_t)>& a_nonDeletionModifier =
      [](ParticleSoA<ItoParticle>&, std::size_t) -> void {
      return;
    });

  /**
   * @brief Compute particle load on a specific grid level -- this will compute the number of particles in each box in
   * the input grids.
   * @param[out] a_loads Loads on each grid patch. This follows the order from a_dbl.boxArray()
   * @param[in]  a_dbl   Grids
   * @param[in]  a_level Grid level
   */
  virtual void
  computeLoads(Vector<long int>& a_loads, const DisjointBoxLayout& a_dbl, int a_level);

  /**
   * @brief Get number of plot variables
   * @return Returns the number of plot variables (Note: Scalars-> one variable and vectors -> SpaceDim variables).
   */
  virtual int
  getNumberOfPlotVariables() const;

  /**
   * @brief Get output plot names
   * @return Returns a string of all variables that this solver can add to plot files.
   */
  virtual Vector<std::string>
  getPlotVariableNames() const;

  /**
   * @brief Write plot data
   * @param[in,out] a_output Output data holder.
   * @param[in]    a_outputRealm Realm where a_output belongs
   * @param[in]    a_level  Grid level
   * @note This will write the plot data in this solver to a_output, starting on a_comp
   * @param[in] a_comp Comp
   */
  virtual void
  writePlotData(LevelData<EBCellFAB>& a_output,
                int&                  a_comp,
                const std::string&    a_outputRealm,
                int                   a_level) const noexcept;

  /**
   * @brief Get number of particles in a specified particle container
   * @param[in] a_whichContainer Container
   * @param[in] a_localOnly      Only report local particles (i.e. not reduced over MPI ranks)
   * @return Num particles
   */
  virtual unsigned long long
  getNumParticles(WhichContainer a_whichContainer, bool a_localOnly) const;

  /**
   * @brief Register operators
   * @details This will register the required operators for running this class.
   */
  virtual void
  registerOperators() const;

  /**
   * @brief Set computational geometry
   * @param[in] a_computationalGeometry Computational geometry.
   */
  virtual void
  setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

  /**
   * @brief Set the AmrMesh object
   * @param[in] a_amr AmrMesh object.
   */
  virtual void
  setAmr(const RefCountedPtr<AmrMesh>& a_amr);

  /**
   * @brief Set phase where the particles will live.
   * @param[in] a_phase Phase
   */
  virtual void
  setPhase(phase::which_phase a_phase);

  /**
   * @brief Set verbosity level for this solver.
   * @param[in] a_verbosity Verbosity level.
   */
  virtual void
  setVerbosity(int a_verbosity);

  /**
   * @brief Set the time for this solver.
   * @param[in] a_step Time step number
   * @param[in] a_time Time (in seconds)
   * @param[in] a_dt   Time step increment
   * @note This sets m_step=a_step, m_time=a_time, m_dt=a_dt
   */
  virtual void
  setTime(int a_step, Real a_time, Real a_dt);

#ifdef CH_USE_HDF5
  /**
   * @brief Write checkpoint data into HDF5 file.
   * @details This will call one of the other two functions writeCheckpointLevelParticles or writeCheckpointLevelFluid,
   * depending on how the user wants to checkpoint the data. In most cases the user will (and should) checkpoint the
   * particles.
   * @param[out] a_handle HDF5 file.
   * @param[in]  a_level Grid level
   */
  virtual void
  writeCheckpointLevel(HDF5Handle& a_handle, int a_level) const;
#endif

#ifdef CH_USE_HDF5
  /**
   * @brief Read checkpointed data from an HDF5 file.
   * @details This will call one of the functions readCheckPointLevelParticle or readCheckponitLevelFluid, depending on
   * how the checkpointing was made.
   * @param[out] a_handle HDF5 file.
   * @param[in]  a_level Grid level
   */
  virtual void
  readCheckpointLevel(HDF5Handle& a_handle, int a_level);
#endif

  /**
   * @brief Set the species
   * @param[in] a_species Species to be solved for.
   */
  virtual void
  setSpecies(const RefCountedPtr<ItoSpecies>& a_species);

  /**
   * @brief Perform pre-regrid operations.
   * @param[in] a_lbase          Coarsest level that changed during regrid.
   * @param[in] a_oldFinestLevel Finest grid level before the regrid operation.
   * @note This caches the bulk particles.
   */
  virtual void
  preRegrid(int a_lbase, int a_oldFinestLevel);

  /**
   * @brief Check if solver is mobile
   * @return Returns m_isMobile
   */
  virtual bool
  isMobile() const;

  /**
   * @brief Check if solver is diffusive
   * @return Returns m_isDiffusive
   */
  virtual bool
  isDiffusive() const;

  /**
   * @brief Get a general particle container
   * @param[in] a_container Which container to fetch.
   * @return Particles
   */
  virtual ParticleContainer<ItoParticle>&
  getParticles(WhichContainer a_container);

  /**
   * @brief Get a general particle container
   * @param[in] a_container Which container to fetch.
   * @return Particles
   */
  virtual const ParticleContainer<ItoParticle>&
  getParticles(WhichContainer a_container) const;

  /**
   * @brief Get the species
   * @return Returns m_species
   */
  virtual const RefCountedPtr<ItoSpecies>&
  getSpecies() const;

  /**
   * @brief Get the mesh data.
   * @return Returns m_phi
   */
  virtual EBAMRCellData&
  getPhi();

  /**
   * @brief Get cell-centered velocity mesh data
   * @return m_velocityFunction
   */
  virtual EBAMRCellData&
  getVelocityFunction();

  /**
   * @brief Get cell-centered velocity mesh data
   * @return m_velocityFunction
   */
  virtual const EBAMRCellData&
  getVelocityFunction() const;

  /**
   * @brief Get the diffusion function
   * @return m_diffusionFunction
   */
  virtual EBAMRCellData&
  getDiffusionFunction();

  /**
   * @brief Get the diffusion function
   * @return m_diffusionFunction
   */
  virtual const EBAMRCellData&
  getDiffusionFunction() const;

  /**
   * @brief Get mobility function
   * @return m_mobilityFunction
   */
  virtual EBAMRCellData&
  getMobilityFunction();

  /**
   * @brief Get mobility function
   * @return m_mobilityFunction
   */
  virtual const EBAMRCellData&
  getMobilityFunction() const;

  /**
   * @brief Utility function -- set a constant diffusion coefficient
   * @param[in] a_diffusionCoefficient Diffusion coefficient
   */
  virtual void
  setDiffusionFunction(Real a_diffusionCoefficient);

  /**
   * @brief Utility function -- set a constant velocity
   * @param[in] a_velocity Velocity
   */
  virtual void
  setVelocityFunction(const RealVect& a_velocity);

  /**
   * @brief Sets mobility coefficient for all particles.
   * @param[in] a_mobility Particle mobility.
   */
  virtual void
  setParticleMobility(Real a_mobility);

  /**
   * @brief Sets diffusion coefficient for all particles.
   * @param[in] a_diffCo Particle diffusion value.
   */
  virtual void
  setParticleDiffusion(Real a_diffCo);

  /**
   * @brief Interpolate the particle velocities.
   * @details This will compute the particle velocities as v = mu * V(Xp) where mu is the particle mobility and V(Xp) is
   * the interpolation of m_velocityFunction to the particle position.
   */
  virtual void
  interpolateVelocities();

  /**
   * @brief Interpolate the particle velocities.
   * @details This will compute the particle velocities as v = mu * V(Xp) where mu is the particle mobility and V(Xp) is
   * the interpolation of m_velocityFunction to the particle position.
   * @param[in] a_lvl Grid level
   * @param[in] a_dit Grid index
   */
  virtual void
  interpolateVelocities(int a_lvl, const DataIndex& a_dit);

  /**
   * @brief Interpolate mobilities
   * @details This will switch between the two ways of computing the particle mobility.
   */
  virtual void
  interpolateMobilities();

  /**
   * @brief Update mobilities parametrically from the particle energy.
   * @details This calls the mobility function implemented by ItoSpecies. The particle mobility is set mu =
   * f(p.energy()) where f is the mobility function in ItoSpecies
   */
  virtual void
  updateMobilities();

  /**
   * @brief Interpolate the diffusion field to the particle positions.
   * @details This computes D_p = Df(X_p) where Df is the diffusion field on the mesh.
   */
  virtual void
  interpolateDiffusion();

  /**
   * @brief Compute the cell-centered gradient of the mesh diffusion field.
   * @details Fills m_diffusionGradient with grad(Df), where Df is the mesh diffusion field, and then coarsens the
   * result and updates its ghost cells. The caller must have brought m_diffusionFunction's own ghost cells up to
   * date first -- the gradient stencil reaches into them and across refinement boundaries. Does nothing when the
   * solver is not diffusive or when ItoSolver.diffusion_grad_drift is false, in which case m_diffusionGradient
   * is not allocated.
   */
  virtual void
  computeDiffusionGradient();

  /**
   * @brief Interpolate the mesh diffusion gradient onto the particles' scratch vector storage.
   * @details This computes (grad D)(X_p) and writes it to the ItoParticle scratch vector columns, from which the
   * caller adds dt*(grad D) to the particle displacement. The scratch columns are set to zero when the solver is
   * not diffusive or the drift correction is disabled, so callers may add the term unconditionally. Requires a
   * preceding call to computeDiffusionGradient().
   */
  virtual void
  interpolateDiffusionGradient();

  /**
   * @brief Interpolate the mesh diffusion gradient onto the particles' scratch vector storage, on one patch.
   * @details This computes (grad D)(X_p) and writes it to the ItoParticle scratch vector columns. Requires a
   * preceding call to computeDiffusionGradient(). Public because the time steppers call it per patch, from
   * inside their own OpenMP loops over boxes.
   * @param[in] a_lvl Grid level
   * @param[in] a_dit Grid index
   */
  virtual void
  interpolateDiffusionGradient(int a_lvl, const DataIndex& a_dit);

  /**
   * @brief Check if the grad(D) drift correction is enabled.
   * @return Returns m_diffusionGradientDrift.
   */
  virtual bool
  isDiffusionGradientDrift() const noexcept;

  /**
   * @brief Update mobilities parametrically from the particle energy.
   * @details This calls the diffusion function implemented by ItoSpecies. The particle diffusion is set D =
   * f(p.energy()) where f is the diffusion function in ItoSpecies
   */
  virtual void
  updateDiffusion();

  /**
   * @brief Sort the input particle container by cell.
   * @param[in] a_container Container to be sorted.
   */
  virtual void
  organizeParticlesByCell(WhichContainer a_container);

  /**
   * @brief Sort the input particle container by patch.
   * @param[in] a_container Container to be sorted.
   */
  virtual void
  organizeParticlesByPatch(WhichContainer a_container);

  /**
   * @brief Make superparticles for a full container, merging to the solver's configured target.
   * @details The usual entry point: the target is ItoSolver.particles_per_cell and the method is
   * ItoSolver.merge_algorithm, so a caller that has no opinion about either does not have to carry one.
   * @param[in] a_container Which container to repartition into new superparticles.
   */
  virtual void
  makeSuperparticles(WhichContainer a_container);

  /**
   * @brief Make superparticles for a full container -- this is the AMR version that users will usually call.
   * @param[in] a_container        Which container to repartition into new superparticles
   * @param[in] a_particlesPerCell Target number of particles per cell
   */
  virtual void
  makeSuperparticles(WhichContainer a_container, int a_particlesPerCell);

  /**
   * @brief Make superparticles for a full container -- this is the AMR version that users will usually call.
   * @param[in] a_container        Which container to repartition into new superparticles
   * @param[in] a_particlesPerCell Target number of particles per cell. Vector indicates grid levels.
   */
  virtual void
  makeSuperparticles(WhichContainer a_container, const Vector<int>& a_particlesPerCell);

  /**
   * @brief Make superparticles with an explicitly chosen method, ignoring the configured
   * merge_algorithm.
   * @details Lets a caller run a different method than the solver's own (e.g. a stepper choosing a
   * regrid-time merge). Routes through mergeLite(); ParticleMergeMethod::None is a
   * no-op.
   * @param[in] a_container        Which container to repartition.
   * @param[in] a_particlesPerCell Target particle count per cell (per-level; the last entry applies to
   *                               any finer level).
   * @param[in] a_method           The merge method to run.
   */
  virtual void
  makeSuperparticles(WhichContainer                                a_container,
                     const Vector<int>&                            a_particlesPerCell,
                     const ParticleManagement::ParticleMergeMethod a_method);

  /**
   * @brief Get the configured target number of computational particles per cell, per AMR level.
   * @details ItoSolver.particles_per_cell. Exposed because a load estimator has to assume the same
   * target the merge will actually reduce to -- BrownianWalkerStepper's mesh balancer clamps its
   * per-cell physical count with it. Index it with the same clamp the merge uses: entry @c lvl if the
   * vector is long enough, otherwise the last entry.
   * @return The per-level target particle count.
   */
  virtual const Vector<int>&
  getParticlesPerCell() const noexcept;

  /**
   * @brief Get the reduced holder the bulk particles ride across a regrid in.
   * @details Only meaningful between preRegrid() and regrid(); the bulk ItoParticle container is empty
   * over that window, and this is where its particles actually are. Exposed for load balancers that
   * need to count or trial-regrid the population between those two calls -- see
   * BrownianWalkerStepper::loadBalanceBoxesParticles(). Reading getParticles(Bulk) there instead
   * yields zero loads and a valid but unbalanced layout, with nothing to indicate anything went wrong.
   * @return The reduced-particle regrid holder.
   */
  virtual ParticleContainer<ItoMergeParticle>&
  getRegridParticles() noexcept;

  /**
   * @brief Remap the bulk particle container.
   */
  virtual void
  remap();

  /**
   * @brief Remap all particles in the input container
   * @param[in] a_container Particle container
   */
  virtual void
  remap(WhichContainer a_container);

  /**
   * @brief Remap all particle containers
   */
  virtual void
  remapAll();

  /**
   * @brief Draw a random N-dimensional Gaussian number from a normal distribution with zero with and unit standard
   * deviation.
   * @note The distribution is truncated at m_normalDistributionTruncation -- values above that threshold will be
   * replaced by m_normalDistributionTruncation.
   * @return Return value
   */
  inline RealVect
  randomGaussian() const;

  /**
   * @brief Draw a random direction in N-dimensional space.
   * @details We use the algorithm by Marsaglia (1972).
   * @return Return value
   */
  inline RealVect
  randomDirection() const;

  /**
   * @brief Compute a time step for the advance -- this calls the level function.
   * @details This computes the time step differently whether or not diffusion and advection are active. The Ito
   * particle model does not have a fundamental time step limitation, so these limits "replicate" the time step
   * selections in a 1D fluid model. If we only use advection advection the time step is computed as dt = dx/sum(|V_i|)
   * = dtA. If only diffusion is active the time step is computed as dt = (dx*dx)/(2*SpaceDim*D) = dtD. If both
   * advection and diffusion are active the time step is computed as dt = 1/(1/dtA + 1/dtD).
   * @return Computed dt
   * @note The grad(D) drift correction (ItoSolver.diffusion_grad_drift) is not accounted for here. These
   * estimates are built from the particle velocity and diffusion columns, and grad(D) lives on neither -- it is
   * interpolated at the point of use and added to the displacement by the time stepper. The omission is
   * deliberate rather than an oversight: the correction displaces a particle by dt*grad(D), which for a streamer
   * in atmospheric air (D ~ 0.1 m^2/s varying over ~10 um, so |grad D| ~ 1e4 m/s) is about two percent of the
   * diffusion hop sqrt(2*D*dt) at dt = 1 ps, and the hop is what the diffusive limit already bounds. A model
   * that made grad(D) comparable to the drift velocity would make these estimates optimistic.
   */
  virtual Real
  computeDt() const;

  /**
   * @brief Compute a time step for the advance -- this returns the maximum permitted time step on the input grid level.
   * @details This computes the time step differently whether or not diffusion and advection are active. The Ito
   * particle model does not have a fundamental time step limitation, so these limits "replicate" the time step
   * selections in a 1D fluid model. If we only use advection advection the time step is computed as dt = dx/Vmax = dtA,
   * where Vmax is the largest velocity component along any of the directions. If only diffusion is active the time step
   * is computed as dt = (dx*dx)/(2*SpaceDim*D) = dtD. If both advection and diffusion are active the time step is
   * computed as dt = 1/(1/dtA + 1/dtD).
   * @param[in] a_lvl Grid level
   * @return Computed dt
   * @note Does not account for the grad(D) drift correction -- see computeDt().
   */
  virtual Real
  computeDt(int a_lvl) const;

  /**
   * @brief Compute a time step for the advance -- this returns the maximum permitted time step on the input grid patch.
   * @details This computes the time step differently whether or not diffusion and advection are active. The Ito
   * particle model does not have a fundamental time step limitation, so these limits "replicate" the time step
   * selections in a 1D fluid model. If we only use advection advection the time step is computed as dt = dx/Vmax = dtA,
   * where Vmax is the largest velocity component along any of the directions. If only diffusion is active the time step
   * is computed as dt = (dx*dx)/(2*SpaceDim*D) = dtD. If both advection and diffusion are active the time step is
   * computed as dt = 1/(1/dtA + 1/dtD).
   * @param[in] a_lvl Grid level
   * @param[in] a_dit Grid index
   * @return Computed dt
   * @note Does not account for the grad(D) drift correction -- see computeDt().
   */
  virtual Real
  computeDt(int a_lvl, const DataIndex& a_dit) const;

  /**
   * @brief Compute the largest possible time step such that the particles does not move more than a specified number of
   * grid cells on any grid level.
   * @details This computes the time step differently whether or not diffusion and advection are active. The Ito
   * particle model does not a fundamental time step limitation, so this function does its best at computing a time step
   * such that the particles do not move more than the specified amount of grid cells. If only advection is active the
   * time step is computed as dt = a_maxCellsToMove * dx/vMax where vMax is the largest velocity component. If only
   * diffusion is active the usual kernel is dX = dWt which yields a "hop" of length L_d = sqrt(2*D*dt)*N0 where N0 is
   * the diffusion hop in one of the coordinates. We do not need to worry about SpaceDim-related corrections because the
   * hops in each direction are independent, so solving for dt we find dtD = L_d^2/(2*D*N0^2). If both advectino and
   * diffusion are active then we can move up to L_d = v_d*dt + sqrt(2*D*dt)*N0 in any one coordinate direction (this is
   * the usual Ito kernel). This requires a solution to a quadratic equation. Fortunately, this is easy to solve for.
   * @param[in] a_maxCellsToMove
   * @return Computed hop dt
   */
  virtual Real
  computeHopDt(Real a_maxCellsToMove) const;

  /**
   * @brief Compute the largest possible time step such that the particles does not move more than a specified number of
   * grid cells on the input grid level.
   * @details This computes the time step differently whether or not diffusion and advection are active. The Ito
   * particle model does not a fundamental time step limitation, so this function does its best at computing a time step
   * such that the particles do not move more than the specified amount of grid cells. If only advection is active the
   * time step is computed as dt = a_maxCellsToMove * dx/vMax where vMax is the largest velocity component. If only
   * diffusion is active the usual kernel is dX = dWt which yields a "hop" of length L_d = sqrt(2*D*dt)*N0 where N0 is
   * the diffusion hop in one of the coordinates. We do not need to worry about SpaceDim-related corrections because the
   * hops in each direction are independent, so solving for dt we find dtD = L_d^2/(2*D*N0^2). If both advectino and
   * diffusion are active then we can move up to L_d = v_d*dt + sqrt(2*D*dt)*N0 in any one coordinate direction (this is
   * the usual Ito kernel). This requires a solution to a quadratic equation. Fortunately, this is easy to solve for.
   * @param[in] a_maxCellsToMove Maximum number of cells to move with a standard Ito kernel dX = v*dt + sqrt(2*D*dt)*N
   * @param[in] a_lvl            Grid level
   * @return Computed hop dt
   */
  virtual Real
  computeHopDt(Real a_maxCellsToMove, int a_lvl) const;

  /**
   * @brief Compute the largest possible time step such that the particles does not move more than a specified number of
   * grid cells on the input grid level.
   * @details This computes the time step differently whether or not diffusion and advection are active. The Ito
   * particle model does not a fundamental time step limitation, so this function does its best at computing a time step
   * such that the particles do not move more than the specified amount of grid cells. If only advection is active the
   * time step is computed as dt = a_maxCellsToMove * dx/vMax where vMax is the largest velocity component. If only
   * diffusion is active the usual kernel is dX = dWt which yields a "hop" of length L_d = sqrt(2*D*dt)*N0 where N0 is
   * the diffusion hop in one of the coordinates. We do not need to worry about SpaceDim-related corrections because the
   * hops in each direction are independent, so solving for dt we find dtD = L_d^2/(2*D*N0^2). If both advectino and
   * diffusion are active then we can move up to L_d = v_d*dt + sqrt(2*D*dt)*N0 in any one coordinate direction (this is
   * the usual Ito kernel). This requires a solution to a quadratic equation. Fortunately, this is easy to solve for.
   * @param[in] a_maxCellsToMove Maximum number of cells to move with a standard Ito kernel dX = v*dt + sqrt(2*D*dt)*N
   * @param[in] a_lvl            Grid level
   * @param[in] a_dit            Grid index
   * @return Computed hop dt
   */
  virtual Real
  computeHopDt(Real a_maxCellsToMove, int a_lvl, const DataIndex& a_dit) const;

  /**
   * @brief Compute advection time step dt = dx/vMax where vMax is the largest velocity component of the particle.
   * @return Computed advective dt
   * @note The grad(D) drift correction (ItoSolver.diffusion_grad_drift) is not accounted for here. These
   * estimates are built from the particle velocity and diffusion columns, and grad(D) lives on neither -- it is
   * interpolated at the point of use and added to the displacement by the time stepper. The omission is
   * deliberate rather than an oversight: the correction displaces a particle by dt*grad(D), which for a streamer
   * in atmospheric air (D ~ 0.1 m^2/s varying over ~10 um, so |grad D| ~ 1e4 m/s) is about two percent of the
   * diffusion hop sqrt(2*D*dt) at dt = 1 ps, and the hop is what the diffusive limit already bounds. A model
   * that made grad(D) comparable to the drift velocity would make these estimates optimistic.
   */
  virtual Real
  computeAdvectiveDt() const;

  /**
   * @brief Compute the drift dt. This computes the minimum dt = dx/vMax on the input level.
   * @param[in] a_lvl Grid level
   * @return Computed advective dt
   * @note Does not account for the grad(D) drift correction -- see computeDt().
   */
  virtual Real
  computeAdvectiveDt(int a_lvl) const;

  /**
   * @brief Compute the drift dt. This computes the minimum dt = dx/vMax on one level and one box.
   * @param[in] a_lvl Grid level
   * @param[in] a_dit Grid index
   * @return Computed advective dt
   * @note Does not account for the grad(D) drift correction -- see computeDt().
   */
  virtual Real
  computeAdvectiveDt(int a_lvl, const DataIndex& a_dit) const;

  /**
   * @brief Compute the diffusive dt. This computes dt = dx*dx/(2*SpaceDim*D) for all particles
   * @return Computed diffusive dt
   * @note The grad(D) drift correction (ItoSolver.diffusion_grad_drift) is not accounted for here. These
   * estimates are built from the particle velocity and diffusion columns, and grad(D) lives on neither -- it is
   * interpolated at the point of use and added to the displacement by the time stepper. The omission is
   * deliberate rather than an oversight: the correction displaces a particle by dt*grad(D), which for a streamer
   * in atmospheric air (D ~ 0.1 m^2/s varying over ~10 um, so |grad D| ~ 1e4 m/s) is about two percent of the
   * diffusion hop sqrt(2*D*dt) at dt = 1 ps, and the hop is what the diffusive limit already bounds. A model
   * that made grad(D) comparable to the drift velocity would make these estimates optimistic.
   */
  virtual Real
  computeDiffusiveDt() const;

  /**
   * @brief Compute the diffusive dt. This computes dt = dx*dx/(2*SpaceDim*D) for all particles on the input level
   * @param[in] a_lvl Grid level
   * @return Computed diffusive dt
   * @note Does not account for the grad(D) drift correction -- see computeDt().
   */
  virtual Real
  computeDiffusiveDt(int a_lvl) const;

  /**
   * @brief Compute the diffusive dt. This computes dt = dx*dx/(2*SpaceDim*D) for all particles on the input grid patch
   * @param[in] a_lvl Grid level
   * @param[in] a_dit Grid index
   * @return Computed diffusive dt
   * @note Does not account for the grad(D) drift correction -- see computeDt().
   */
  virtual Real
  computeDiffusiveDt(int a_lvl, const DataIndex& a_dit) const;

  /**
   * @brief Get deposition method
   * @return Returns m_deposition
   */
  DepositionType
  getDeposition() const;

  /**
   * @brief Get the coarse-fine deposition strategy
   * @return Coarse fine deposition
   */
  CoarseFineDeposition
  getCoarseFineDeposition() const;

  /**
   * @brief Return phase
   * @return Returns m_phase
   */
  phase::which_phase
  getPhase() const;

protected:
  /**
   * @brief Default component in data holders (do not touch).
   */
  static constexpr int m_comp = 0;

  /**
   * @brief Default number of component in data holders (do not touch).
   */
  static constexpr int m_nComp = 1;

  /**
   * @brief How to checkpoint files
   * @details Particles => Write particles to HDF5. Numbers => Write particle numbers to HDF5 (and lose information)
   */
  enum class WhichCheckpoint
  {
    Particles,
    Numbers
  };

  /**
   * @brief Enum for deciding how to interpolate particle mobilities.
   * @details Direct => Interpolate particle mobilities from the mesh. Velocity => Set particle mobilitys by computing
   * mu = |v|/|V| where v is the particle velocity and V is the "velocity field".
   */
  enum class WhichMobilityInterpolation
  {
    Direct,
    Velocity
  };

  /**
   * @brief How to checkpoint files. particles => write particles to HDF5. numbers => write numbers to HDF5
   */
  WhichCheckpoint m_checkpointing;

  /**
   * @brief Switch for deciding how to interpolate mobilities, i.e. interpolating either mu*E or just mu (to the
   * particle position)
   */
  WhichMobilityInterpolation m_mobilityInterp;

  /**
   * @brief User-supplied per-cell merger for merge_algorithm = external (set via
   * setParticleCellMerger()). The built-in methods build their own merger locally in the corresponding
   * mergeXxx().
   */
  ParticleManagement::ParticleMerger<ItoMergeParticle> m_particleCellMerger;

  /**
   * @brief The configured merge method (which specific algorithm).
   * @details Set by parseParticleMerger() from ItoSolver.merge_algorithm. makeSuperparticles()
   * dispatches on this to run the right algorithm. Defaults to EqualWeightKD (the default
   * merge_algorithm).
   */
  ParticleManagement::ParticleMergeMethod m_mergeMethod = ParticleManagement::ParticleMergeMethod::EqualWeightKD;

  /**
   * @brief Target number of computational particles per cell, per AMR level.
   * @details Set by parseParticleMerger() from ItoSolver.particles_per_cell. The solver owns this
   * because the solver is where merging happens -- it is the target every merge reduces to, the count
   * initialData() seeds the bulk container at, and (through getParticlesPerCell()) the number a
   * particle-count load estimate has to assume. A time stepper still owns *when* to merge during a run
   * (ItoKMCStepper.merge_interval); it no longer owns *what to merge to*.
   *
   * Per-level, for compatibility with the ItoKMCStepper key it replaces. Only the cell-based methods
   * honour it per level; the whole-container methods use one crowding threshold for the whole
   * hierarchy and read entry 0, so parseParticleMerger() aborts on a non-uniform vector rather than
   * silently truncating it.
   */
  Vector<int> m_particlesPerCell;

  /**
   * @brief Which merge method to run inside ItoSolver::regrid(), or nullopt for "whatever
   * merge_algorithm says".
   * @details Set by parseParticleMerger() from ItoSolver.regrid_superparticles. A separate knob from
   * merge_algorithm because the two are genuinely independent: the regrid merge runs once per regrid
   * on a population that has just been redistributed, and a run may reasonably want a cheaper (or no)
   * method there than it uses per step. ParticleMergeMethod::None turns the regrid merge off; the
   * extract/rebuild round-trip still happens, since preRegrid() has already reduced the particles by
   * then.
   */
  std::optional<ParticleManagement::ParticleMergeMethod> m_regridMergeMethod;

  /**
   * @brief nn_pair_tree/nn_pair_onecell/nn_pair_hash: iterate the local trivial tier to convergence within each
   * round.
   * @details Forwarded to ParticleManagement::mergeNearestNeighborsTree()/mergeNearestNeighborsOneCell()/
   * mergeNearestNeighborsHash(). Only used for nn_pair_tree/nn_pair_onecell/nn_pair_hash.
   */
  bool m_nnPairIterate = true;

  /**
   * @brief nn_pair_tree/nn_pair_onecell/nn_pair_hash: number of fallback candidates per query. The k-nearest
   * search is done with k = 1 + this. Only used for nn_pair_tree/nn_pair_onecell/nn_pair_hash.
   */
  int m_nnPairFallback = 1;

  /**
   * @brief nn_pair_tree/nn_pair_hash: maximum merge distance in whole cells (Chebyshev). A negative
   * value means unbounded (std::nullopt). Only used for nn_pair_tree/nn_pair_hash -- nn_pair_onecell's merge
   * distance is structurally fixed at 1 and does not read this.
   */
  int m_nnPairMaxCellDistance = -1;

  /**
   * @brief nn_pair_tree/nn_pair_onecell/nn_pair_hash: maximum number of drain rounds run per makeSuperparticles()
   * call.
   * @details Each round merges nearest-neighbor pairs and roughly halves an over-full cell's surplus; the
   * loop also stops early once a round merges nothing. This caps the round count so the cost stays bounded
   * and predictable rather than running to full convergence. Only used for nn_pair_tree/nn_pair_onecell/nn_pair_hash.
   */
  int m_nnPairMaxRounds = 3;

  /**
   * @brief Leaf size, in cell widths, at or below which the tree build switches from count-median
   * to weight-median splitting.
   * @details Forwarded to the a_splitWeightLeafDx argument of ParticleManagement::mergeKDCarve()
   * and mergeKDPatch() -- see those functions' own docs. Larger values apply the weight median to bigger leaves, giving
   * more uniform super-particle weights at some cost to how evenly leaves are apportioned between cells; 0 disables it
   * entirely. Governs the shared tree build, so it applies to both kd_carve and kd_patch.
   */
  Real m_kdSplitWeightLeafDx = 1.0;

  /**
   * @brief Number of particles used when restarting a simulation -- this is relevant only when restarting from a
   * "fluid" checkpoint file.
   */
  int m_restartPPC;

  /**
   * @brief Realm where this solve lives.
   */
  std::string m_realm;

  /**
   * @brief Computational geometry.
   */
  RefCountedPtr<ComputationalGeometry> m_computationalGeometry;

  /**
   * @brief AMR; needed for grid stuff
   */
  RefCountedPtr<AmrMesh> m_amr;

  /**
   * @brief Species that this solver solves for.
   */
  RefCountedPtr<ItoSpecies> m_species;

  /**
   * @brief Phase where this solver lives.
   */
  phase::which_phase m_phase;

  /**
   * @brief Solver name
   */
  std::string m_name;

  /**
   * @brief Class name
   * @details ItoSolver for parent class -- derived classes might be named something else.
   */
  std::string m_className;

  /**
   * @brief Truncation value for normal distribution.
   * @details This is used when drawing numbers from the Gaussian distribution, replacing all values above
   * m_normalDistributionTruncation by m_normalDistributionTruncation.
   */
  Real m_normalDistributionTruncation;

  /**
   * @brief Bisection step size to use for particle intersection tests with EBs.
   */
  Real m_bisectionStep;

  /**
   * @brief Verbosity level for this solver.
   */
  int m_verbosity;

  /**
   * @brief Time step
   */
  int m_timeStep;

  /**
   * @brief Current time
   */
  Real m_time;

  /**
   * @brief Time step
   */
  Real m_dt;

  /**
   * @brief NGP deposition in cut cells or not
   */
  bool m_forceIrregDepositionNGP;

  /**
   * @brief NGP interpolation in cut cells or not
   */
  bool m_forceIrregInterpolationNGP;

  /**
   * @brief Force usage of NGP when depositing "halo" particles.
   */
  bool m_forceHaloNGP;

  /**
   * @brief Use mass redistribution or not.
   */
  bool m_useRedistribution;

  /**
   * @brief Flag for blending the deposition clouds with the "non-conservative" divergence.
   */
  bool m_blendConservation;

  /**
   * @brief If true, solver is diffusive.
   */
  bool m_isDiffusive;

  /**
   * @brief If true, solver is mobile, i.e. advection is turned on.
   */
  bool m_isMobile;

  /**
   * @brief If true, grad(D) is added to the particle drift.
   * @details Set from ItoSolver.diffusion_grad_drift. An Ito-interpretation Euler-Maruyama update with a drift
   * of v alone transports the flux v*n - grad(D*n); adding grad(D) to the drift makes it transport v*n - D*grad(n)
   * instead, which is the flux CdrSolver integrates. The two coincide only where D is constant in space.
   */
  bool m_diffusionGradientDrift;

  /**
   * @brief Flag for outputting m_phi to plot files.
   */
  bool m_plotPhi;

  /**
   * @brief Flag for outputting m_velocityFunction to plot files.
   */
  bool m_plotVelocity;

  /**
   * @brief Flag for outputting m_diffusionFunction to plot files.
   */
  bool m_plotDiffCo;

  /**
   * @brief Flag for outputting m_diffusionGradient to plot files.
   */
  bool m_plotDiffusionGradient;

  /**
   * @brief Flag for depositing and plotting the bulk particles on the mesh.
   */
  bool m_plotParticles;

  /**
   * @brief Flag for depositing and plotting the EB particles on the mesh.
   */
  bool m_plotParticlesEB;

  /**
   * @brief Flag for depositing and plotting the domain particles on the mesh.
   */
  bool m_plotParticlesDomain;

  /**
   * @brief Flag for depositing and plotting the source particles on the mesh.
   */
  bool m_plotParticlesSource;

  /**
   * @brief Flag for depositing and plotting the covered particles on the mesh.
   */
  bool m_plotParticlesCovered;

  /**
   * @brief Flag for plotting the energy density on the mesh
   */
  bool m_plotEnergyDensity;

  /**
   * @brief Flag for plotting the average particle energy on the mesh
   */
  bool m_plotAverageEnergy;

  /**
   * @brief Deposition method when depositing particles to the mesh.
   */
  DepositionType m_deposition;

  /**
   * @brief Coarse-fine deposition strategy
   */
  CoarseFineDeposition m_coarseFineDeposition;

  /**
   * @brief Deposition method when depositing particles to mesh during plotting.
   */
  DepositionType m_plotDeposition;

  /**
   * @brief Storage for particle mesh data.
   */
  EBAMRCellData m_phi;

  /**
   * @brief Mobility function -- used when interpolating particle mobilities.
   */
  EBAMRCellData m_mobilityFunction;

  /**
   * @brief Cell-centered field used for interpolating velocities
   */
  EBAMRCellData m_velocityFunction;

  /**
   * @brief Diffusion-centerer field used for interpolating diffusion coefficients
   */
  EBAMRCellData m_diffusionFunction;

  /**
   * @brief Cell-centered gradient of m_diffusionFunction, used for the grad(D) drift correction.
   * @details Allocated only when the solver is diffusive and m_diffusionGradientDrift is true; otherwise this is a
   * bare pointer holder. Recomputed from m_diffusionFunction by computeDiffusionGradient(), so it carries no state
   * across a regrid and is never checkpointed.
   */
  EBAMRCellData m_diffusionGradient;

  /**
   * @brief Per-cell occupancy scratch for the kd merges (see ParticleManagement::mergeKDCarve).
   * @details Owned here, not by the merge, so the allocation happens once per regrid rather than once
   * per merge call. Allocated and released with the other mesh fields, in allocate()/regrid() and
   * preRegrid() respectively.
   */
  EBAMRFAB m_kdMergeCellHistogram;

  /**
   * @brief Per-cell merge-quota scratch for the kd merges. Same lifetime as m_kdMergeCellHistogram.
   */
  EBAMRFAB m_kdMergeLeafQuota;

  /**
   * @brief Scratch storage for holding the non-conservative deposition
   */
  mutable EBAMRIVData m_depositionNC;

  /**
   * @brief Scratch storage for holding the mass difference when using hybrid deposition
   */
  mutable EBAMRIVData m_massDiff;

  /**
   * @brief Various particle containers with identifiers.
   */
  std::map<WhichContainer, ParticleContainer<ItoParticle>> m_particleContainers;

  /**
   * @brief The bulk particles, in reduced form, between preRegrid() and regrid().
   * @details At a de-refining regrid 2^D fine cells collapse into one coarse cell, so between the
   * redistribution and the merge each rank transiently holds up to 2^D times the target particles per
   * cell. This container is what that spike is held in: preRegrid() extracts the bulk container into
   * it, the redistribution and load balancing happen on it, and regrid() merges it and only then
   * rebuilds ItoParticles. 53 B per particle in 3-D instead of 149 B, and the merge runs before the
   * population is ever materialised at the larger size.
   *
   * Bulk only. The EB/Domain/Source/Covered/Scratch containers keep the plain
   * ParticleContainer::preRegrid() path -- they are small, none of them is merged, and Scratch is
   * never populated.
   *
   * Live only across a regrid: regrid() releases it before returning. It is a member rather than a
   * local because preRegrid() and regrid() are separate calls with the load-balancing phase between
   * them, and because BrownianWalkerStepper's particle load balancer has to count it.
   */
  ParticleContainer<ItoMergeParticle> m_regridParticles;

  /**
   * @brief Algorithm for EB intersection
   */
  EBIntersection m_intersectionAlg;

  /**
   * @brief Parse RNG options -- this parses the RNG seed and instantiates the distributions.
   */
  void
  parseRNG();

  /**
   * @brief Parse the normal distribution truncation level
   */
  void
  parseTruncation();

  /**
   * @brief Parse particle deposition methods.
   */
  void
  parseDeposition();

  /**
   * @brief Parse the super-particle merger
   */
  void
  parseParticleMerger();

  /**
   * @brief Run one merge method over an already-extracted reduced-particle container.
   * @details The single place that knows the *shape* of each merge method, so that no caller has to.
   * Two shapes exist and they are not interchangeable:
   *   - cell-based (equal_weight_kd, reinitialize, reinitialize_bvh, nn_sfc, external): reduce each
   *     cell to at most its per-level target, refilling only within that same per-cell budget. These
   *     must NOT be followed by the container-wide split.
   *   - whole-container (kd_carve, kd_patch, kd_skin_nn, nn_pair_{tree,onecell,hash}): drain over-full
   *     cells against one crowding threshold for the whole hierarchy, then split under-full cells back
   *     up to it. These MUST split, and they mint scratch ids that are invalidated afterwards.
   * ParticleMergeMethod::None does neither -- not even the split, which would otherwise refill every
   * under-full cell in the domain at every regrid. Getting that wrong is silent: nothing fails, the
   * population simply grows.
   * @param[in,out] a_particles        Reduced-particle container to merge in place.
   * @param[in]     a_particlesPerCell Target particle count per cell (per-level for the cell-based
   *                                   methods; the whole-container methods use the first entry only).
   * @param[in]     a_method           Which merge method to run.
   */
  void
  mergeLite(ParticleContainer<ItoMergeParticle>&    a_particles,
            const Vector<int>&                      a_particlesPerCell,
            ParticleManagement::ParticleMergeMethod a_method);

  /**
   * @brief Apply a per-cell merger over an entire reduced-particle container, cell by cell.
   * @details Shared machinery for every cell-based merge: the container is cell-sorted here
   * (organizeParticlesByCell()); for each non-empty cell the particles are extracted, passed to
   * a_merger with the cell's CellInfo and target count, and the result rebuilt into the leaf.
   * @param[in,out] a_particles        Reduced-particle container to repartition.
   * @param[in]     a_particlesPerCell Target particle count per cell (per-level; the last entry
   *                                   applies to any finer level).
   * @param[in]     a_merger           The per-cell merger to apply.
   */
  void
  applyCellMerger(ParticleContainer<ItoMergeParticle>&                        a_particles,
                  const Vector<int>&                                          a_particlesPerCell,
                  const ParticleManagement::ParticleMerger<ItoMergeParticle>& a_merger);

  /**
   * @brief Set every particle id in a reduced-particle container to ParticleSoA::s_invalidID.
   * @details Run by mergeLite() after the whole-container methods, whose ids come from a rank-local
   * counter that restarts every drain round. Those are merge scratch, not particle identity, and are
   * indistinguishable from real ids downstream -- where they would reach H5Part output. "No id" is the
   * honest answer for a particle a merge produced, and is what an append()-built particle carries.
   * @param[in,out] a_particles Container whose id column is reset.
   */
  void
  invalidateParticleIDs(ParticleContainer<ItoMergeParticle>& a_particles) const;

  /**
   * @brief Extract a container's ItoParticles into a minimal ItoMergeParticle container on the
   * same realm/patch, then clear the ItoParticles.
   * @details Copies position + weight + particleID + energy only; velocity/mobility/diffusion are
   * dropped (they are recomputed from the field after the merge). No remap is needed since
   * positions/ownership are unchanged. Ghost particles are cleared first -- the copy loop has no
   * isGhost() test and would promote a ghost to an owned particle. Each source leaf's arena is
   * released as soon as its reduced copy exists, so the two representations overlap by one patch
   * rather than one level.
   * @param[in]  a_container Which container to extract from.
   * @param[out] a_merge     The resulting minimal-payload container, allocated on m_realm.
   */
  void
  extractIntoMergeContainer(WhichContainer a_container, ParticleContainer<ItoMergeParticle>& a_merge);

  /**
   * @brief Split every under-full cell up to the target count.
   * @details The whole-container merge/drain phase only reduces over-full cells; this brings cells
   * below the target back up to it via greedy heaviest-particle splitting (same rule as nn_sfc).
   * Called only from mergeLite(), and only for the methods whose shape calls for it -- see there.
   * @param[in,out] a_merge                     Reduced-particle container to split (cell-sorted internally).
   * @param[in]     a_numParticlesPerCellThresh Target particle count per cell.
   */
  void
  splitFromMergeContainer(ParticleContainer<ItoMergeParticle>& a_merge, int a_numParticlesPerCellThresh);

  /**
   * @brief Rebuild a container's ItoParticles from the reduced particles.
   * @details Rebuilds into the same patch (no remap); velocity/mobility/diffusion default to zero and
   * are recomputed from the field before use. particleID is carried across. Each reduced leaf's arena
   * is released as its ItoParticle copy is built. Leaves a_container patch-organized.
   * @param[in]     a_container Which container to rebuild into.
   * @param[in,out] a_merge     The reduced-particle container to consume.
   */
  void
  rebuildFromMergeContainer(WhichContainer a_container, ParticleContainer<ItoMergeParticle>& a_merge);

  /**
   * @brief Which EBGeometry point-cloud index backs mergeNnPairSearch()'s per-patch
   * spatial search -- see ParticleManagement::mergeNearestNeighborsTree()/mergeNearestNeighborsHash().
   */
  enum class NnPairSearchBackend
  {
    Tree, ///< EBGeometry::PointCloudBVH -- ParticleManagement::mergeNearestNeighborsTree().
    Hash  ///< EBGeometry::PointCloudHashGrid -- ParticleManagement::mergeNearestNeighborsHash().
  };

  /**
   * @brief Cell-based merge: equal-weight KD partition, one weighted-centroid particle per leaf.
   * @details Recursively partitions each cell into at most a_particlesPerCell equal-weight KD leaves
   * and reduces each to one weighted-centroid particle. Requires particle weights >= 1 to split.
   * Honours the per-level target.
   * @param[in,out] a_particles        Reduced-particle container to repartition.
   * @param[in]     a_particlesPerCell Target particle count per cell (per-level).
   */
  void
  mergeEqualWeightKD(ParticleContainer<ItoMergeParticle>& a_particles, const Vector<int>& a_particlesPerCell);

  /**
   * @brief Cell-based merge: redistribute each cell's physical particles into near-equal integer
   * weights at random cell positions.
   * @details All output particles carry the same weight-averaged energy. Requires integer weights.
   * Honours the per-level target.
   * @param[in,out] a_particles        Reduced-particle container to repartition.
   * @param[in]     a_particlesPerCell Target particle count per cell (per-level).
   */
  void
  mergeReinitialize(ParticleContainer<ItoMergeParticle>& a_particles, const Vector<int>& a_particlesPerCell);

  /**
   * @brief Cell-based merge: as mergeEqualWeightKD(), with reinitialized leaf positions.
   * @details Cut-cells use the weighted centroid (to stay inside the EB); full cells draw a random
   * point in the leaf bounding box. Honours the per-level target.
   * @param[in,out] a_particles        Reduced-particle container to repartition.
   * @param[in]     a_particlesPerCell Target particle count per cell (per-level).
   */
  void
  mergeReinitializeBVH(ParticleContainer<ItoMergeParticle>& a_particles, const Vector<int>& a_particlesPerCell);

  /**
   * @brief Cell-based merge: Hilbert-curve sort, then merge adjacent pairs down to the target.
   * @details Splits the heaviest particle when a cell is below target, so unlike the other cell-based
   * methods it does refill -- but per cell against the per-level target, which is why mergeLite() must
   * still not run the container-wide split after it. Does not require integer weights.
   * @param[in,out] a_particles        Reduced-particle container to repartition.
   * @param[in]     a_particlesPerCell Target particle count per cell (per-level).
   */
  void
  mergeNnSfc(ParticleContainer<ItoMergeParticle>& a_particles, const Vector<int>& a_particlesPerCell);

  /**
   * @brief Whole-container merge: distributed nearest-neighbour pair merge -- shared implementation
   * behind the nn_pair_tree and nn_pair_hash selectors.
   * @details Whole-container, cross-patch/cross-rank AMR merge. Identical in every respect between the
   * two backends (callbacks, id allocation, drain loop); only the
   * ParticleManagement::mergeNearestNeighborsTree()/mergeNearestNeighborsHash() call inside the drain
   * loop differs, selected at compile time via @c Backend. Requires the particle ghost mask registered
   * on this realm (see registerOperators()); it does its own ghost fill and needs no cell sorting. Uses
   * a single crowding threshold -- the coarsest level's a_particlesPerCell entry. Drains over-full
   * cells only; mergeLite() splits afterwards.
   * @tparam Backend Which point-cloud index to search with (see NnPairSearchBackend).
   * @param[in,out] merge              Reduced-particle container to repartition.
   * @param[in]     a_particlesPerCell Target particle count (only the first entry is used).
   */
  template <NnPairSearchBackend Backend>
  void
  mergeNnPairSearch(ParticleContainer<ItoMergeParticle>& merge, const Vector<int>& a_particlesPerCell);

  /**
   * @brief Whole-container merge: distributed nearest-neighbour pair merge, per-cell PointCloudBVH search.
   * @details As mergeNnPairSearch(), but the merge distance is structurally fixed at Chebyshev cell
   * distance 1 (own cell and Moore-adjacent neighbours only), so nn_pair_max_cell_dist is not honoured
   * by this method. Requires a width-1 particle ghost mask registered on this realm. Drains over-full
   * cells only; mergeLite() splits afterwards.
   * @param[in,out] merge              Reduced-particle container to repartition.
   * @param[in]     a_particlesPerCell Target particle count (only the first entry is used).
   */
  void
  mergeNnPairOneCell(ParticleContainer<ItoMergeParticle>& merge, const Vector<int>& a_particlesPerCell);

  /**
   * @brief Whole-container merge: kd-tree with a nearest-neighbour skin (the kd_skin_nn selector).
   * @details Same tree build and per-cell quota as mergeKDImpl(), but the boundary tier is the
   * nearest-neighbour pair merge instead of the carve arbitration. The particles are split across two
   * containers: ParticleManagement::mergeKDInterior() commits every leaf holding no ghost into one, and
   * what is left -- the contested skin -- is drained from the other by
   * ParticleManagement::mergeNearestNeighborsOneCell(), whose search radius is fixed at Chebyshev cell
   * distance 1 and so matches the width-1 halo exactly. The skin tier drains each cell to what remains
   * of its target after the interior results are counted against it, so the two tiers share one budget.
   * Fills a width-1 particle ghost mask itself, which must be registered on this realm. Tuned by
   * nn_pair_iterate/nn_pair_fallback/nn_pair_max_rounds, shared with the nn_pair family;
   * nn_pair_max_cell_dist does not apply. Drains over-full cells only; mergeLite() splits afterwards.
   * @param[in,out] merge              Reduced-particle container to repartition.
   * @param[in]     a_particlesPerCell Target particle count (only the first entry is used).
   */
  void
  mergeKDSkinNn(ParticleContainer<ItoMergeParticle>& merge, const Vector<int>& a_particlesPerCell);

  /**
   * @brief Whole-container merge: whole-patch kd-tree build -- shared implementation behind the
   * kd_carve and kd_patch selectors.
   * @details A single, non-iterative pass: no drain loop, and no nn_pair_iterate/nn_pair_fallback/
   * nn_pair_max_cell_dist equivalents, since there is no configurable search radius or convergence
   * loop. With a_enableBoundaryCarve the merge fills a width-1 particle ghost mask (which must be
   * registered on this realm) and arbitrates leaves straddling a patch boundary; without it the merge
   * is patch-local -- no ghosts, no communication, and a cell straddled by two patches is reduced by
   * each patch independently. Needs no cell sorting. Uses a single crowding threshold -- the coarsest
   * level's a_particlesPerCell entry. Drains over-full cells only; mergeLite() splits afterwards.
   * @param[in,out] merge                 Reduced-particle container to repartition.
   * @param[in]     a_particlesPerCell    Target particle count (only the first entry is used).
   * @param[in]     a_enableBoundaryCarve True to fill ghosts and run the carve tier, false for the
   *                                      patch-local merge.
   */
  void
  mergeKDImpl(ParticleContainer<ItoMergeParticle>& merge,
              const Vector<int>&                   a_particlesPerCell,
              bool                                 a_enableBoundaryCarve);

  /**
   * @brief Parse EB intersection algorithms
   */
  void
  parseIntersectionEB();

  /**
   * @brief Parse class verbosity
   */
  void
  parseVerbosity();

  /**
   * @brief Parse plot variables
   */
  void
  parsePlotVariables();

  /**
   * @brief Parse diffusion hop
   */
  void
  parseDiffusionHop();

  /**
   * @brief Parse whether or not to add the grad(D) drift correction.
   */
  void
  parseDiffusionGradient();

  /**
   * @brief Parse whether or not to use redistribution.
   */
  void
  parseRedistribution();

  /**
   * @brief Parse whether or not to compute a "non-conservative" divergence when redistributing mass.
   */
  void
  parseDivergenceComputation();

  /**
   * @brief Parse checkpointing method
   */
  void
  parseCheckpointing();

#ifdef CH_USE_HDF5
  /**
   * @brief Write checkpoint data into HDF5 file -- this version dumps that particles to the HDF5 file..
   * @details This will not write the full Ito partilces -- in order to trim memory and checkpoint file sizes this
   * function writse the particle data to HDF5 as SimpleItoParticle which only stores position, mass, and energy. Things
   * like mobility, diffusion coefficient, velocity etc. are parametrically dependent and are thus not stored.
   * @param[out] a_handle HDF5 file.
   * @param[in]  a_level Grid level
   */
  virtual void
  writeCheckPointLevelParticles(HDF5Handle& a_handle, int a_level) const;
#endif

#ifdef CH_USE_HDF5
  /**
   * @brief Write checkpoint data into HDF5 file -- this version writes the number of particles per cell rather than the
   * particles.
   * @details This will write the number of particles in each grid cell -- this yields very-low-memory HDF5 files but
   * when the simulation is restarted from this state then the particle distribution as well as the particle energy in
   * each cell is lost. Use with caution!
   * @param[out] a_handle HDF5 file.
   * @param[in]  a_level Grid level
   */
  virtual void
  writeCheckPointLevelFluid(HDF5Handle& a_handle, int a_level) const;
#endif

#ifdef CH_USE_HDF5
  /**
   * @brief Read checkpointed particles from  an HDF5 file.
   * @details This will read particles from the HDF5 file ala writeCheckPointLevelParticles.
   * @param[out] a_handle HDF5 file.
   * @param[in]  a_level Grid level
   */
  virtual void
  readCheckpointLevelParticles(HDF5Handle& a_handle, int a_level);
#endif

#ifdef CH_USE_HDF5
  /**
   * @brief Read checkpointed particle numberse from  an HDF5 file and instantiate the particles from that.
   * @details This will read particle numbers in each cell ala writeCheckPointLevelFluid and instantiate the particles
   * from that.
   * @param[out] a_handle HDF5 file.
   * @param[in]  a_level Grid level
   */
  virtual void
  readCheckpointLevelFluid(HDF5Handle& a_handle, int a_level);
#endif

  /**
   * @brief Restart particles from a specified number of particles in the grid cell.
   * @details This will instantiate the bulk particles by randomly drawing new particles in each grid cell.
   * @param[in] a_particlesPerCell Mesh data showing how many particles are in each cell.
   * @param[in] a_level            Grid level.
   * @param[in] a_newPPC           Desired number of computational particles per cell
   */
  void
  drawNewParticles(const LevelData<EBCellFAB>& a_particlesPerCell, int a_level, int a_newPPC);

  /**
   * @brief Redistribute mass in an AMR context.
   * @details We will have deposited particles into each cell, i.e. phi_i = m_i/dx^3. To obtain the true density in an
   * EB context we need to divide by kappa such that phi_i = m_i/(kappa_i*dx^3). Unfortunately, this is numerically
   * unstable because kappa_i can be arbitrarily small. We can set the density as phi_i = m_i/dx^3 but we will be
   * missing a mass m_i*(1-kappa_i) from the cell. This mass can be smooshed into neighboring cells such as to make the
   * total scheme conservative. As an option, we can also use a "non-conservative" divergence but the scheme is not
   * guaranteed to be non-negative. All of this is just to say that we take the mass that fell inside the EB and put it
   * back into the domain.
   * @param[in,out] a_phi Quantity to be redistributed.
   */
  virtual void
  redistributeAMR(EBAMRCellData& a_phi) const;

  /**
   * @brief Make the "non-conservative" kappa deposition -- computing depositionNC = sum(kappa*depositionC)/sum(kappa)
   * in a neighborhood around each vof.
   * @param[out] a_depositionNC     Non-conservative deposition
   * @param[in]  a_depositionKappaC Conserved deposition
   */
  virtual void
  depositNonConservative(EBAMRIVData& a_depositionNC, const EBAMRCellData& a_depositionKappaC) const;

  /**
   * @brief Make the "hybrid" deposition phiH = kappa*phiC + (1-kappa)*phiNC. On input, a_depositionH should contain
   * phiC
   * @param[in,out] a_depositionH    On input, contains phiC. On output, contain phiH
   * @param[out]   a_massDifference On output, contains mass loss in each cut-cell.
   * @param[in]    a_depositionNC   The "non-conservative" deposited variable.
   */
  virtual void
  depositHybrid(EBAMRCellData& a_depositionH, EBAMRIVData& a_massDifference, const EBAMRIVData& a_depositionNC) const;

  /**
   * @brief Interpolate mobilities -- this will switch between the two ways of computing the particle mobility.
   * @param[in] a_lvl               Grid level.
   * @param[in] a_dit               Grid index.
   * @param[in] a_velocityMagnitude Velocity magnitude.
   */
  virtual void
  interpolateMobilities(int a_lvl, const DataIndex& a_dit, const EBCellFAB& a_velocityMagnitude) noexcept;

  /**
   * @brief Directly interpolate mobilities. Interpolates for all particles in the specified grid patch.
   * @param[in] a_lvl Grid level.
   * @param[in] a_dit Grid index.
   */
  virtual void
  interpolateMobilitiesDirect(int a_lvl, const DataIndex& a_dit) noexcept;

  /**
   * @brief Interpolate mobilities through the velocity. Interpolates for all particles in the specified grid patch.
   * @param[in] a_lvl               Grid level.
   * @param[in] a_dit               Grid index.
   * @param[in] a_velocityMagnitude Velocity magnitude.
   */
  virtual void
  interpolateMobilitiesVelocity(int a_lvl, const DataIndex& a_dit, const EBCellFAB& a_velocityMagnitude) noexcept;

  /**
   * @brief Update mobilities parametrically from the particle energy.
   * @details This calls the mobility function implemented by ItoSpecies. The particle mobility is set mu =
   * f(p.energy()) where f is the mobility function in ItoSpecies
   * @param[in] a_level Grid level
   * @param[in] a_dit   Grid index
   */
  virtual void
  updateMobilities(int a_level, const DataIndex& a_dit);

  /**
   * @brief Interpolate the diffusion field to the particle positions.
   * @details This computes D_p = Df(X_p) where Df is the diffusion field on the mesh.
   * @param[in] a_lvl Grid level
   * @param[in] a_dit Grid index
   */
  virtual void
  interpolateDiffusion(int a_lvl, const DataIndex& a_dit);

  /**
   * @brief Update mobilities parametrically from the particle energy.
   * @details This calls the diffusion function implemented by ItoSpecies. The particle diffusion is set D =
   * f(p.energy()) where f is the diffusion function in ItoSpecies
   * @param[in] a_level Grid level
   * @param[in] a_dit   Grid index
   */
  virtual void
  updateDiffusion(int a_level, const DataIndex& a_dit);

  /**
   * @brief Write data to output. Convenience function.
   * @param[in,out] a_output Output data holder.
   * @param[in]    a_data   Data to write.
   * @param[in]    a_outputRealm Realm where a_output belongs
   * @param[in]    a_level  Grid level
   * @param[in]    a_interpToCentroids If true, a_data will be interpolated to cell centroids before writing to
   * a_output.
   * @param[in]    a_interpGhost If true, interpolate ghost cells
   * @param[in] a_comp Comp
   */
  virtual void
  writeData(LevelData<EBCellFAB>& a_output,
            int&                  a_comp,
            const EBAMRCellData&  a_data,
            const std::string&    a_outputRealm,
            int                   a_level,
            bool                  a_interpToCentroids,
            bool                  a_interpGhost) const noexcept;
};

#include <CD_NamespaceFooter.H>

#include <CD_ItoSolverImplem.H>

#endif
