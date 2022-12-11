/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoSolver.H
  @brief  Declaration of solver class for Ito diffusion
  @author Robert Marskar
*/

#ifndef CD_ItoSolver_H
#define CD_ItoSolver_H

// Our includes
#include <CD_AmrMesh.H>
#include <CD_ComputationalGeometry.H>
#include <CD_ItoSpecies.H>
#include <CD_ItoParticle.H>
#include <CD_EBParticleMesh.H>
#include <CD_ParticleContainer.H>
#include <CD_EBRepresentation.H>
#include <CD_EBIntersection.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for Ito diffusion particle models.
  @details This is a particle class for solving Ito diffusion problems dX = v*dt + sqrt(2*D)*dW_t (dW_t is a Wiener process over dt) over an AMR hierarchy. The diffusion
  coefficient is related to the hydrodynamic diffusion coefficient by D(fluid) = sqrt(2*D(Ito)). 
  @note To use this class, the user must add a species and set the realm and phase where the solver will live. 
*/
class ItoSolver
{
public:
  /*!
    @brief Enum class for distinguishing various types of particle containers.
    @details This exists because the ItoSolver can partition particles into various containers, which is very useful when one wants to add particles from
    a source term, remove particles that fall inside the EB, or parse boundary conditions on the EB and domain faces. Here,
    Bulk    = "Active" particles,
    EB      = Particles on the EBs ,
    Domain  = Particles on the domain sides,
    Source  = Source particles -- used to e.g. add new particles to the bulk particles,
    Covered = Particles inside the EB,
    Scratch = Scratch particles
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

  /*!
    @brief Constructor -- user must subsequently set the realm and, parse class options, set the species etc. 
  */
  ItoSolver();

  /*!
    @brief Destructor (does nothing).
  */
  virtual ~ItoSolver();

  /*!
    @brief Get this solver's name 
    @return Returns m_name
  */
  virtual std::string
  getName() const;

  /*!
    @brief Get the realm where this solver is registered.
    @return m_realm
  */
  virtual const std::string
  getRealm() const;

  /*!
    @brief Set the realm where this solver will live
    @param[in] a_realm Realm name. 
  */
  virtual void
  setRealm(const std::string a_realm);

  /*!
    @brief Parse class options
  */
  virtual void
  parseOptions();

  /*!
    @brief Parse class runtime configurable options.
  */
  virtual void
  parseRuntimeOptions();

  /*!
    @brief Get output plot names
    @return Returns a string of all variables that this solver can add to plot files. 
  */
  virtual Vector<std::string>
  getPlotVariableNames() const;

  /*!
    @brief Set the initial data
    @details This will add the initial particles and deposit them. 
  */
  virtual void
  initialData();

  /*!
    @brief Regrid this solver. 
    @param[in] a_lmin           Coarsest level where grids did not change. 
    @param[in] a_oldFinestLevel Finest AMR level before the regrid. 
    @param[in] a_newFinestLevel Finest AMR level after the regrid. 
  */
  virtual void
  regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Allocate internal storage
    @details This will allocate the required mesh data as well as all the particle data holders. 
    @note Storage for diffusion and advection fields are only allocate if the species is diffusive/mobile. 
  */
  virtual void
  allocateInternals();

  /*!
    @brief Clear a particle container -- this will delete all the particles in the input container. 
    @param[in] a_whichContainer Particle container to be cleared.
  */
  virtual void
  clear(const WhichContainer a_container);

  /*!
    @brief Clear a particle container -- this will delete all the particles in the input container. 
    @param[in] a_particles Particle container to be emptied. 
  */
  virtual void
  clear(ParticleContainer<ItoParticle>& a_particles) const;

  /*!
    @brief Clear particles -- this will deleted all the particles in the input container. 
    @param[in] a_particles Particle container to be emptied. 
  */
  virtual void
  clear(AMRParticles<ItoParticle>& a_particles) const;

  /*!
    @brief Deposit conductivities (i.e. mass*mobility / volume)
    @details This deposits mass*mobility (not multiplied by charge)
    @param[out] a_phi       Mesh data
    @param[in]  a_particles Particle data
  */
  virtual void
  depositConductivity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const;

  /*!
    @brief Deposit conductivities (i.e. mass*mobility / volume)
    @details This deposits mass*mobility (not multiplied by charge)
    @param[out] a_phi                  Mesh data
    @param[in]  a_particles            Particle data
    @param[in]  a_deposition           Deposition method
    @param[in]  a_coarseFineDeposition Coarse-fine deposition method. 
  */
  virtual void
  depositConductivity(EBAMRCellData&                  a_phi,
                      ParticleContainer<ItoParticle>& a_particles,
                      const DepositionType            a_deposition,
                      const CoarseFineDeposition      a_coarseFineDeposition) const;

  /*!
    @brief Deposit diffusivity (i.e. mass*D/volume)
    @details This deposits mass*diffusion (not multiplied by charge)
    @param[out] a_phi       Mesh data
    @param[in]  a_particles Particle data
    @note Calls the other versions with a_deposition = m_deposition
  */
  virtual void
  depositDiffusivity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const;

  /*!
    @brief Deposit diffusivity (i.e. mass*D/volume)
    @details This deposits mass*mobility (not multiplied by charge)
    @param[out] a_phi                  Mesh data
    @param[in]  a_particles            Particle data
    @param[in]  a_deposition           Deposition method
    @param[in]  a_coarseFineDeposition Coarse-fine deposition method. 
  */
  virtual void
  depositDiffusivity(EBAMRCellData&                  a_phi,
                     ParticleContainer<ItoParticle>& a_particles,
                     const DepositionType            a_deposition,
                     const CoarseFineDeposition      a_coarseFineDeposition) const;

  /*!
    @brief Deposit energy densities (i.e. mass*energy/volume => total energy per unit volume)
    @details This deposits mass*energy
    @param[out] a_phi       Mesh data
    @param[in]  a_particles Particle data
    @note Calls the other versions with a_deposition = m_deposition
  */
  virtual void
  depositEnergyDensity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const;

  /*!
    @brief Deposit energy densities (i.e. mass*energy/volume => total energy per unit volume)
    @details This deposits mass*energy
    @param[out] a_phi                  Mesh data
    @param[in]  a_particles            Particle data
    @param[in]  a_deposition           Deposition method
    @param[in]  a_coarseFineDeposition Coarse-fine deposition method. 
  */
  virtual void
  depositEnergyDensity(EBAMRCellData&                  a_phi,
                       ParticleContainer<ItoParticle>& a_particles,
                       const DepositionType            a_deposition,
                       const CoarseFineDeposition      a_coarseFineDeposition) const;

  /*!
    @brief Compute average mobility. 
    @details This computes the average mobility as mu_avg = sum(mass*mu)/sum(mass) by depositing on the mesh. 
    @param[out] a_phi       Mesh data
    @param[in]  a_particles Particle data
  */
  virtual void
  computeAverageMobility(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const;

  /*!
    @brief Compute average diffusion coefficient. 
    @details This computes the average diffusion coefficient as D_avg = sum(mass*D)/sum(mass) by depositing on the mesh. 
    @param[out] a_phi       Mesh data
    @param[in]  a_particles Particle data
  */
  virtual void
  computeAverageDiffusion(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const;

  /*!
    @brief Compute average particle energy.. 
    @details This computes the average energy as E_avg = sum(mass*energy)/sum(mass) by depositing on the mesh. 
    @param[out] a_phi       Mesh data
    @param[in]  a_particles Particle data
  */
  virtual void
  computeAverageEnergy(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const;

  /*!
    @brief Deposit particles onto mesh. 
    @details This will deposit the mass (i.e., computational weight) "bulk" particles into m_phi. 
    @note Cells the other version with a_container = WhichContainer::Bulk
  */
  virtual void
  depositParticles();

  /*!
    @brief Deposit particles on to mesh. 
    @param[in] a_container Which container to deposit.
    @details This will deposit mass (i.e., computational weight) of the the input particle container particles onto the classes member 'm_phi'. 
    @note Calls the general version with arguments: m_phi, m_particles.at(a_container), m_deposition.
  */
  virtual void
  depositParticles(const WhichContainer a_container);

  /*!
    @brief Generic deposition function which deposits a particle field onto the mesh using a specified deposition method. 
    @details The template parameters indicate the particle type and quantity to be deposited. The second template parameter must be a pointer to a member function in 
    the particle class with signature 'const Real& P::function() const'. E.g. 'const Real& P::mass() const' which is the default quantity to be deposited. 
    @param[out] a_phi                  Mesh data -- must have exactly one component. 
    @param[in]  a_particles            Particles to be deposited
  */
  template <class P, const Real& (P::*particleScalarField)() const>
  void
  depositParticles(EBAMRCellData& a_phi, ParticleContainer<P>& a_particles) const;

  /*!
    @brief Generic deposition function which deposits a particle field onto the mesh using a specified deposition method. 
    @details The template parameters indicate the particle type and quantity to be deposited. The second template parameter must be a pointer to a member function in 
    the particle class with signature 'const Real& P::function() const'. E.g. 'const Real& P::mass() const' which is the default quantity to be deposited. 
    @param[out] a_phi                  Mesh data -- must have exactly one component. 
    @param[in]  a_particles            Particles to be deposited
    @param[in]  a_deposition           Deposition method
    @param[in]  a_coarseFineDeposition Coarse-fine deposition strategy
  */
  template <class P, const Real& (P::*particleScalarField)() const>
  void
  depositParticles(EBAMRCellData&             a_phi,
                   ParticleContainer<P>&      a_particles,
                   const DepositionType       a_deposition,
                   const CoarseFineDeposition a_coarseFineDeposition) const;

  /*!
    @brief Generic deposition function which deposits a particle field onto the mesh using a specified deposition method. 
    @details Just like the version above, but with a different function signature.
    @param[out] a_phi                  Mesh data -- must have exactly one component. 
    @param[in]  a_particles            Particles to be deposited
    @param[in]  a_deposition           Deposition method
    @param[in]  a_coarseFineDeposition Coarse-fine deposition strategy
  */
  template <class P, Real (P::*particleScalarField)() const>
  void
  depositParticles(EBAMRCellData&             a_phi,
                   ParticleContainer<P>&      a_particles,
                   const DepositionType       a_deposition,
                   const CoarseFineDeposition a_coarseFineDeposition) const;

  /*!
    @brief Remove particles that are inside the EB. 
    @param[in] a_ebRepresentation EB representation (e.g. implicit function, discrete, voxel)
    @param[in] a_tol              Tolerance -- removes particles if they are less than a_tol*dx away from the boundary. 
    @details Calls the other version with m_particles.
  */
  virtual void
  removeCoveredParticles(const EBRepresentation a_representation, const Real a_tol);

  /*!
    @brief Remove particles that are inside the EB. 
    @param[in] a_container        Which particle container to remove particles from
    @param[in] a_ebRepresentation EB representation (e.g. implicit function, discrete, voxel)
    @param[in] a_tol              Tolerance -- removes particles if they are less than a_tol*dx away from the boundary. 
    @note Calls the other version with ParticleContainer = m_particles.at(a_container)
  */
  virtual void
  removeCoveredParticles(const WhichContainer a_container, const EBRepresentation a_representation, const Real a_tol);

  /*!
    @brief Remove particles that are inside the EB
    @param[in] a_particles        Particles to remove. 
    @param[in] a_ebRepresentation EB representation (e.g. implicit function, discrete, voxel)
    @param[in] a_tol              Tolerance -- removes particles if they are less than a_tol*dx away from the boundary. 
  */
  virtual void
  removeCoveredParticles(ParticleContainer<ItoParticle>& a_particles,
                         const EBRepresentation          a_representation,
                         const Real                      a_tol) const;

  /*!
    @brief Transfer particles that are covered (within a_tol distance from EB) to another container. 
    @details Calls the other version with m_particles and m_coveredParticles
    @param[in] a_ebRepresentation EB representation (e.g. implicit function, discrete, voxel)
    @param[in] a_tol              Tolerance -- removes particles if they are less than a_tol*dx away from the boundary. 
  */
  virtual void
  transferCoveredParticles(const EBRepresentation a_representation, const Real a_tol);

  /*!
    @brief Transfer particles that are covered (within a_tol distance from EB) to another container. 
    @param[in] a_containerFrom    Which particle container to transfer particles from
    @param[in] a_containerFrom    Which particle container to move particles into
    @param[in] a_ebRepresentation EB representation (e.g. implicit function, discrete, voxel)
    @param[in] a_tol              Tolerance -- removes particles if they are less than a_tol*dx away from the boundary. 
  */
  virtual void
  transferCoveredParticles(const WhichContainer   a_containerFrom,
                           const WhichContainer   a_containerTo,
                           const EBRepresentation a_representation,
                           const Real             a_tol);

  /*!
    @brief Transfer particles that are inside the EB to another container. 
    @param[in] a_particlesFrom    Which particle container to transfer particles from
    @param[in] a_particlesTo      Which particle container to move particles into
    @param[in] a_ebRepresentation EB representation (e.g. implicit function, discrete, voxel)
    @param[in] a_tol              Tolerance -- removes particles if they are less than a_tol*dx away from the boundary. 
  */
  virtual void
  transferCoveredParticles(ParticleContainer<ItoParticle>& a_particlesFrom,
                           ParticleContainer<ItoParticle>& a_particlesTo,
                           const EBRepresentation          a_representation,
                           const Real                      a_tol) const;

  /*!
    @brief Do boundary intersection tests.
    @details This will intersect the particles in the "bulk" particles data holder with the domain faces and EBs. If a particle crossed the EB it will
    be put into the "EB" particle data holder and likewise for the particles that crossed the domain side. 
    @param[in] a_ebIntersection   Enum for switching between various types of intersection tests. 
    @param[in] a_deleteParticles  If true, the origin particle will also be removed from the bulk particle data holder.
    @note This will call the other version with WhichContainer::Bulk, WhichContainer::EB, and WhichContainer::Domain.
  */
  virtual void
  intersectParticles(const EBIntersection a_ebIntersection, const bool a_deleteParticles);

  /*!
    @brief Do boundary intersection tests.
    @details This will intersect the particles in the a_particles data holder with the domain faces and EBs. If a particle crossed the EB it will
    be put into the a_eb_particles particle data holder and likewise for the particles that crossed the domain side. 
    @param[in] a_particles        Which particle container to use for the "active" particles
    @param[in] a_ebParticles      Which particle container to put the EB particles
    @param[in] a_domainParticles  Which particle container to put the domain particles
    @param[in] a_ebIntersection   Enum for switching between various types of intersection tests. 
    @param[in] a_deleteParticles  If true, the origin particle will also be removed from the bulk particle data holder.
    @note This will look up the ParticleContainer and call the other version.
  */
  virtual void
  intersectParticles(const WhichContainer a_particles,
                     const WhichContainer a_ebParticles,
                     const WhichContainer a_domainParticles,
                     const EBIntersection a_ebIntersection,
                     const bool           a_deleteParticles);

  /*!
    @brief Do boundary intersection tests.
    @details This will intersect the particles in the a_particles data holder with the domain faces and EBs. If a particle crossed the EB it will
    be put into the a_eb_particles particle data holder and likewise for the particles that crossed the domain side. 
    @param[in] a_particles        Which particle container to use for the "active" particles
    @param[in] a_ebParticles      Which particle container to put the EB particles
    @param[in] a_domainParticles  Which particle container to put the domain particles
    @param[in] a_ebIntersection   Enum for switching between various types of intersection tests. 
    @param[in] a_deleteParticles  If true, the origin particle will also be removed from the bulk particle data holder. 
  */
  virtual void
  intersectParticles(ParticleContainer<ItoParticle>& a_particles,
                     ParticleContainer<ItoParticle>& a_ebParticles,
                     ParticleContainer<ItoParticle>& a_domainParticles,
                     const EBIntersection            a_ebIntersection,
                     const bool                      a_deleteParticles);

  /*!
    @brief Compute particle load on a specific grid level -- this will compute the number of particles in each box in the input grids. 
    @param[out] a_loads Loads on each grid patch. This follows the order from a_dbl.boxArray()
    @param[in]  a_dbl   Grids
    @param[in]  a_level Grid level
  */
  virtual void
  computeLoads(Vector<long int>& a_loads, const DisjointBoxLayout& a_dbl, const int a_level);

  /*!
    @brief Get number of plot variables
    @return Returns the number of plot variables (Note: Scalars-> one variable and vectors -> SpaceDim variables).
  */
  virtual int
  getNumberOfPlotVariables() const;

  /*!
    @brief Get number of particles in a specified particle container
    @param[in] a_whichContainer Which container to count particles in. 
    @param[in] a_localOnly      Only report local particles (i.e. not reduced over MPI ranks)
  */
  virtual unsigned long long
  getNumParticles(const WhichContainer a_container, const bool a_localOnly) const;

  /*!
    @brief Register operators
    @details This will register the required operators for running this class. 
  */
  virtual void
  registerOperators() const;

  /*!
    @brief Set computational geometry
    @param[in] a_computationalGeometry Computational geometry. 
  */
  virtual void
  setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

  /*!
    @brief Set the AmrMesh object
    @param[in] a_amr AmrMesh object. 
  */
  virtual void
  setAmr(const RefCountedPtr<AmrMesh>& a_amr);

  /*!
    @brief Set phase where the particles will live.
    @param[in] a_phase Phase
  */
  virtual void
  setPhase(const phase::which_phase a_phase);

  /*!
    @brief Set verbosity level for this solver. 
    @param[in] a_verbosity Verbosity level.
  */
  virtual void
  setVerbosity(const int a_verbosity);

  /*!
    @brief Set the time for this solver. 
    @param[in] a_step Time step number
    @param[in] a_time Time (in seconds)
    @param[in] a_dt   Time step increment
    @note This sets m_step=a_step, m_time=a_time, m_dt=a_dt
  */
  virtual void
  setTime(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Write plot data
    @param[inout] a_output Output data holder.
    @param[inout] a_icomp  Starting component where this solver begins writing the output. 
    @note This will write the plot data in this solver to a_output, starting on a_comp
  */
  virtual void
  writePlotData(EBAMRCellData& a_output, int& a_comp);

#ifdef CH_USE_HDF5
  /*!
    @brief Write checkpoint data into HDF5 file. 
    @details This will call one of the other two functions writeCheckpointLevelParticles or writeCheckpointLevelFluid, depending on how the user wants to 
    checkpoint the data. In most cases the user will (and should) checkpoint the particles.
    @param[out] a_handle HDF5 file. 
    @param[in]  a_level Grid level
  */
  virtual void
  writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const;
#endif

#ifdef CH_USE_HDF5
  /*!
    @brief Read checkpointed data from an HDF5 file.
    @details This will call one of the functions readCheckPointLevelParticle or readCheckponitLevelFluid, depending on how the checkpointing was made. 
    @param[out] a_handle HDF5 file. 
    @param[in]  a_level Grid level
  */
  virtual void
  readCheckpointLevel(HDF5Handle& a_handle, const int a_level);
#endif

  /*!
    @brief Set the species
    @param[in] a_species Species to be solved for. 
  */
  virtual void
  setSpecies(const RefCountedPtr<ItoSpecies>& a_species);

  /*!
    @brief Perform pre-regrid operations. 
    @param[in] a_lbase          Coarsest level that changed during regrid. 
    @param[in] a_oldFinestLevel Finest grid level before the regrid operation. 
    @note This caches the bulk particles. 
  */
  virtual void
  preRegrid(const int a_lbase, const int a_oldFinestLevel);

  /*!
    @brief Add particles to a contain. This adds into a specific grid level and patch. The user can delete the input particles. 
    @param[inout] a_inputParticles Input particles -- can be destroyed if a_destructive=true
    @param[in]    a_lvl            Grid level
    @param[in]    a_dit            Grid index
    @param[in]    a_destructive    If true, delete the input particles. 
  */
  template <ItoSolver::WhichContainer C = WhichContainer::Bulk>
  void
  addParticles(ListBox<ItoParticle>& a_inputParticles,
               const int             a_lvl,
               const DataIndex       a_dit,
               const bool            a_destructive);

  /*!
    @brief Check if solver is mobile
    @return Returns m_isMobile
  */
  virtual bool
  isMobile() const;

  /*!
    @brief Check if solver is diffusive
    @return Returns m_isDiffusive
  */
  virtual bool
  isDiffusive() const;

  /*!
    @brief Get a general particle container
    @param[in] a_container Which container to fetch.
  */
  virtual ParticleContainer<ItoParticle>&
  getParticles(const WhichContainer a_container);

  /*!
    @brief Get a general particle container
    @param[in] a_container Which container to fetch.
  */
  virtual const ParticleContainer<ItoParticle>&
  getParticles(const WhichContainer a_container) const;

  /*!
    @brief Get the species
    @return Returns m_species
  */
  virtual const RefCountedPtr<ItoSpecies>&
  getSpecies() const;

  /*!
    @brief Get the mesh data.
    @return Returns m_phi
  */
  virtual EBAMRCellData&
  getPhi();

  /*!
    @brief Get cell-centered velocity mesh data
    @return m_velocityFunction
  */
  virtual EBAMRCellData&
  getVelocityFunction();

  /*!
    @brief Get cell-centered velocity mesh data
    @return m_velocityFunction
  */
  virtual const EBAMRCellData&
  getVelocityFunction() const;

  /*!
    @brief Get the diffusion function
    @return m_diffusionFunction
  */
  virtual EBAMRCellData&
  getDiffusionFunction();

  /*!
    @brief Get the diffusion function
    @return m_diffusionFunction
  */
  virtual const EBAMRCellData&
  getDiffusionFunction() const;

  /*!
    @brief Get scratch storage
    @return Returns m_scratch
  */
  virtual EBAMRCellData&
  getScratch();

  /*!
    @brief Get scratch storage
    @return Returns m_scratch
  */
  virtual const EBAMRCellData&
  getScratch() const;

  /*!
    @brief Get mobility function
    @return m_mobilityFunction
  */
  virtual EBAMRCellData&
  getMobilityFunction();

  /*!
    @brief Get mobility function
    @return m_mobilityFunction
  */
  virtual const EBAMRCellData&
  getMobilityFunction() const;

  /*!
    @brief Utility function -- set a constant diffusion coefficient
    @param[in] a_diffusionCoefficient Diffusion coefficient
  */
  virtual void
  setDiffusionFunction(const Real a_diffusionCoefficient);

  /*!
    @brief Utility function -- set a constant velocity
    @param[in] a_velocity Velocity
  */
  virtual void
  setVelocityFunction(const RealVect a_velocity);

  /*!
    @brief Sets mobility coefficient for all particles.
    @param[in] a_mobility Particle mobility.
  */
  virtual void
  setParticleMobility(const Real a_mobility);

  /*!
    @brief Sets diffusion coefficient for all particles.
    @param[in] a_diffusion Particle diffusion value.
  */
  virtual void
  setParticleDiffusion(const Real a_diffusion);

  /*!
    @brief Interpolate the particle velocities.
    @details This will compute the particle velocities as v = mu * V(Xp) where mu is the particle mobility and V(Xp) is the interpolation of m_velocityFunction
    to the particle position. 
  */
  virtual void
  interpolateVelocities();

  /*!
    @brief Interpolate the particle velocities.
    @details This will compute the particle velocities as v = mu * V(Xp) where mu is the particle mobility and V(Xp) is the interpolation of m_velocityFunction
    to the particle position. 
    @param[in] a_level Grid level
    @param[in] a_dit   Grid index
  */
  virtual void
  interpolateVelocities(const int a_level, const DataIndex& a_dit);

  /*!
    @brief Interpolate mobilities
    @details This will switch between the two ways of computing the particle mobility.
  */
  virtual void
  interpolateMobilities();

  /*!
    @brief Interpolate mobilities -- this will switch between the two ways of computing the particle mobility.
    @param[in] a_level Grid level
    @param[in] a_dit   Grid index
  */
  virtual void
  interpolateMobilities(const int a_level, const DataIndex& a_dit);

  /*!
    @brief Update mobilities parametrically from the particle energy. 
    @details This calls the mobility function implemented by ItoSpecies. The particle mobility is set mu = f(p.energy()) where f is the mobility function in ItoSpecies
  */
  virtual void
  updateMobilities();

  /*!
    @brief Update mobilities parametrically from the particle energy.
    @details This calls the mobility function implemented by ItoSpecies. The particle mobility is set mu = f(p.energy()) where f is the mobility function in ItoSpecies
    @param[in] a_level Grid level
    @param[in] a_dit   Grid index
  */
  virtual void
  updateMobilities(const int a_level, const DataIndex a_dit);

  /*!
    @brief Interpolate the diffusion field to the particle positions. 
    @details This computes D_p = Df(X_p) where Df is the diffusion field on the mesh. 
  */
  virtual void
  interpolateDiffusion();

  /*!
    @brief Interpolate the diffusion field to the particle positions. 
    @details This computes D_p = Df(X_p) where Df is the diffusion field on the mesh. 
    @param[in] a_level Grid level
    @param[in] a_dit   Grid index
  */
  virtual void
  interpolateDiffusion(const int a_level, const DataIndex& a_dit);

  /*!
    @brief Update mobilities parametrically from the particle energy. 
    @details This calls the diffusion function implemented by ItoSpecies. The particle diffusion is set D = f(p.energy()) where f is the diffusion function in ItoSpecies
  */
  virtual void
  updateDiffusion();

  /*!
    @brief Update mobilities parametrically from the particle energy. 
    @details This calls the diffusion function implemented by ItoSpecies. The particle diffusion is set D = f(p.energy()) where f is the diffusion function in ItoSpecies
    @param[in] a_level Grid level
    @param[in] a_dit   Grid index    
  */
  virtual void
  updateDiffusion(const int a_level, const DataIndex a_dit);

  /*!
    @brief Sort the input particle container by cell.
    @param[in] a_container Container to be sorted. 
  */
  virtual void
  organizeParticlesByCell(const WhichContainer a_container);

  /*!
    @brief Sort the input particle container by patch.
    @param[in] a_container Container to be sorted. 
  */
  virtual void
  organizeParticlesByPatch(const WhichContainer a_container);

  /*!
    @brief Make superparticles for a full container -- this is the AMR version that users will usually call. 
    @param[in] a_container        Which container to repartition into new superparticles
    @param[in] a_particlesPerCell Target number of particles per cell
  */
  virtual void
  makeSuperparticles(const WhichContainer a_container, const int a_particlesPerCell);

  /*!
    @brief Make superparticles for a particle container -- this is used for a specified grid level. 
    @param[in] a_container        Which container to repartition into new superparticles
    @param[in] a_particlesPerCell Target number of particles per cell
    @param[in] a_level            Grid level
  */
  virtual void
  makeSuperparticles(const WhichContainer a_container, const int a_particlesPerCell, const int a_level);

  /*!
    @brief Make superparticles for a particle container -- this is used for a specified grid level and patch
    @param[in] a_container        Which container to repartition into new superparticles
    @param[in] a_particlesPerCell Target number of particles per cell
    @param[in] a_level            Grid level
    @param[in] a_dit              Grid index
  */
  virtual void
  makeSuperparticles(const WhichContainer a_container,
                     const int            a_particlesPerPatch,
                     const int            a_level,
                     const DataIndex      a_dit);

  /*!
    @brief Remap the bulk particle container.
  */
  virtual void
  remap();

  /*!
    @brief Remap all particles in the input container
    @param[in] a_container Particle container
  */
  virtual void
  remap(const WhichContainer a_container);

  /*!
    @brief Remap all particle containers
  */
  virtual void
  remapAll();

  /*!
    @brief Draw a random N-dimensional Gaussian number from a normal distribution with zero with and unit standard deviation.
    @note The distribution is truncated at m_normalDistributionTruncation -- values above that threshold will be replaced by m_normalDistributionTruncation.
  */
  inline RealVect
  randomGaussian() const;

  /*!
    @brief Draw a random direction in N-dimensional space. 
    @details We use the algorithm by Marsaglia (1972).
  */
  inline RealVect
  randomDirection() const;

  /*!
    @brief Compute a time step for the advance -- this calls the level function. 
    @details This computes the time step differently whether or not diffusion and advection are active. The Ito particle model does not have a fundamental 
    time step limitation, so these limits "replicate" the time step selections in a 1D fluid model. 
    If we only use advection advection the time step is computed as dt = dx/Vmax = dtA, where Vmax is the largest velocity component along any of the directions. 
    If only diffusion is active the time step is computed as dt = (dx*dx)/(2*D) = dtD. 
    If both advection and diffusion are active the time step is computed as dt = 1/(1/dtA + 1/dtD). 
  */
  virtual Real
  computeDt() const;

  /*!
    @brief Compute a time step for the advance -- this returns the maximum permitted time step on the input grid level.
    @details This computes the time step differently whether or not diffusion and advection are active. The Ito particle model does not have a fundamental 
    time step limitation, so these limits "replicate" the time step selections in a 1D fluid model. 
    If we only use advection advection the time step is computed as dt = dx/Vmax = dtA, where Vmax is the largest velocity component along any of the directions. 
    If only diffusion is active the time step is computed as dt = (dx*dx)/(2*D) = dtD. 
    If both advection and diffusion are active the time step is computed as dt = 1/(1/dtA + 1/dtD). 
    @param[in] a_lvl Grid level
  */
  virtual Real
  computeDt(const int a_lvl) const;

  /*!
    @brief Compute a time step for the advance -- this returns the maximum permitted time step on the input grid patch.
    @details This computes the time step differently whether or not diffusion and advection are active. The Ito particle model does not have a fundamental 
    time step limitation, so these limits "replicate" the time step selections in a 1D fluid model. 
    If we only use advection advection the time step is computed as dt = dx/Vmax = dtA, where Vmax is the largest velocity component along any of the directions. 
    If only diffusion is active the time step is computed as dt = (dx*dx)/(2*D) = dtD. 
    If both advection and diffusion are active the time step is computed as dt = 1/(1/dtA + 1/dtD). 
    @param[in] a_lvl Grid level
    @param[in] a_dit Grid index
  */
  virtual Real
  computeDt(const int a_lvl, const DataIndex& a_dit) const;

  /*!
    @brief Compute the largest possible time step such that the particles does not move more than a specified number of grid cells on any grid level.
    @details This computes the time step differently whether or not diffusion and advection are active. The Ito particle model does not a fundamental
    time step limitation, so this function does its best at computing a time step such that the particles do not move more than the specified amount of grid cells.
    If only advection is active the time step is computed as dt = a_maxCellsToMove * dx/vMax where vMax is the largest velocity component. 
    If only diffusion is active the usual kernel is dX = dWt which yields a "hop" of length L_d = sqrt(2*D*dt)*N0 where N0 is the diffusion hop in one of the coordinates.
    We do not need to worry about SpaceDim-related corrections because the hops in each direction are independent, so solving for dt we find dtD = L_d^2/(2*D*N0^2). 
    If both advectino and diffusion are active then we can move up to L_d = v_d*dt + sqrt(2*D*dt)*N0 in any one coordinate direction (this is the usual Ito kernel). This
    requires a solution to a quadratic equation. Fortunately, this is easy to solve for. 
    @param[in] a_maxCellsToMove
  */
  virtual Real
  computeHopDt(const Real a_maxCellsToMove) const;

  /*!
    @brief Compute the largest possible time step such that the particles does not move more than a specified number of grid cells on the input grid level.
    @details This computes the time step differently whether or not diffusion and advection are active. The Ito particle model does not a fundamental
    time step limitation, so this function does its best at computing a time step such that the particles do not move more than the specified amount of grid cells.
    If only advection is active the time step is computed as dt = a_maxCellsToMove * dx/vMax where vMax is the largest velocity component. 
    If only diffusion is active the usual kernel is dX = dWt which yields a "hop" of length L_d = sqrt(2*D*dt)*N0 where N0 is the diffusion hop in one of the coordinates.
    We do not need to worry about SpaceDim-related corrections because the hops in each direction are independent, so solving for dt we find dtD = L_d^2/(2*D*N0^2). 
    If both advectino and diffusion are active then we can move up to L_d = v_d*dt + sqrt(2*D*dt)*N0 in any one coordinate direction (this is the usual Ito kernel). This
    requires a solution to a quadratic equation. Fortunately, this is easy to solve for. 
    @param[in] a_maxCellsToMove Maximum number of cells to move with a standard Ito kernel dX = v*dt + sqrt(2*D*dt)*N
    @param[in] a_lvl            Grid level
  */
  virtual Real
  computeHopDt(const Real a_maxCellsToMove, const int a_lvl) const;

  /*!
    @brief Compute the largest possible time step such that the particles does not move more than a specified number of grid cells on the input grid level.
    @details This computes the time step differently whether or not diffusion and advection are active. The Ito particle model does not a fundamental
    time step limitation, so this function does its best at computing a time step such that the particles do not move more than the specified amount of grid cells.
    If only advection is active the time step is computed as dt = a_maxCellsToMove * dx/vMax where vMax is the largest velocity component. 
    If only diffusion is active the usual kernel is dX = dWt which yields a "hop" of length L_d = sqrt(2*D*dt)*N0 where N0 is the diffusion hop in one of the coordinates.
    We do not need to worry about SpaceDim-related corrections because the hops in each direction are independent, so solving for dt we find dtD = L_d^2/(2*D*N0^2). 
    If both advectino and diffusion are active then we can move up to L_d = v_d*dt + sqrt(2*D*dt)*N0 in any one coordinate direction (this is the usual Ito kernel). This
    requires a solution to a quadratic equation. Fortunately, this is easy to solve for. 
    @param[in] a_maxCellsToMove Maximum number of cells to move with a standard Ito kernel dX = v*dt + sqrt(2*D*dt)*N
    @param[in] a_lvl            Grid level
    @param[in] a_dit            Grid index
  */
  virtual Real
  computeHopDt(const Real a_maxCellsToMove, const int a_lvl, const DataIndex& a_dit) const;

  /*!
    @brief Compute advection time step dt = dx/vMax where vMax is the largest velocity component of the particle. 
  */
  virtual Real
  computeAdvectiveDt() const;

  /*!
    @brief Compute the drift dt. This computes the minimum dt = dx/vMax on the input level. 
    @param[in] a_lvl Grid level
  */
  virtual Real
  computeAdvectiveDt(const int a_lvl) const;

  /*!
    @brief Compute the drift dt. This computes the minimum dt = dx/vMax on one level and one box. 
    @param[in] a_lvl Grid level
    @param[in] a_dit Grid index
  */
  virtual Real
  computeAdvectiveDt(const int a_lvl, const DataIndex& a_dit) const;

  /*!
    @brief Compute the diffusive dt. This computes dt = dx*dx/(2*D) for all particles
  */
  virtual Real
  computeDiffusiveDt() const;

  /*!
    @brief Compute the diffusive dt. This computes dt = dx*dx/(2*D) for all particles on the input level
    @param[in] a_lvl Grid level
  */
  virtual Real
  computeDiffusiveDt(const int a_lvl) const;

  /*!
    @brief Compute the diffusive dt. This computes dt = dx*dx/(2*D) for all particles on the input grid patch
    @param[in] a_lvl Grid level
    @param[in] a_dit Grid index
  */
  virtual Real
  computeDiffusiveDt(const int a_lvl, const DataIndex& a_dit) const;

  /*!
    @brief Get deposition method
    @return Returns m_deposition
  */
  DepositionType
  getDeposition() const;

  /*!
    @brief Get the coarse-fine deposition strategy
  */
  CoarseFineDeposition
  getCoarseFineDeposition() const;

  /*!
    @brief Return phase
    @return Returns m_phase
  */
  phase::which_phase
  getPhase() const;

protected:
  /*!
    @brief Default component in data holders (do not touch). 
  */
  static constexpr int m_comp = 0;

  /*!
    @brief Default number of component in data holders (do not touch).
  */
  static constexpr int m_nComp = 1;

  /*! 
    @brief How to checkpoint files
    @details Particles => Write particles to HDF5. Numbers => Write particle numbers to HDF5 (and lose information)
  */
  enum class WhichCheckpoint
  {
    Particles,
    Numbers
  };

  /*!
    @brief Enum for deciding how to interpolate particle mobilities. 
    @details Direct => Interpolate particle mobilities from the mesh. Velocity => Set particle mobilitys by computing mu = |v|/|V| where v is the particle velocity and
    V is the "velocity field". 
  */
  enum class WhichMobilityInterpolation
  {
    Direct,
    Velocity
  };

  /*!
    @brief How to checkpoint files. particles => write particles to HDF5. numbers => write numbers to HDF5
  */
  WhichCheckpoint m_checkpointing;

  /*!
    @brief Switch for deciding how to interpolate mobilities, i.e. interpolating either mu*E or just mu (to the particle position)
  */
  WhichMobilityInterpolation m_mobilityInterp;

  /*!
    @brief Number of particles used when restarting a simulation -- this is relevant only when restarting from a "fluid" checkpoint file. 
  */
  int m_restartPPC;

  /*!
    @brief Realm where this solve lives. 
  */
  std::string m_realm;

  /*!
    @brief Computational geometry. 
  */
  RefCountedPtr<ComputationalGeometry> m_computationalGeometry;

  /*!
    @brief AMR; needed for grid stuff
  */
  RefCountedPtr<AmrMesh> m_amr;

  /*!
    @brief Species that this solver solves for. 
  */
  RefCountedPtr<ItoSpecies> m_species;

  /*!
    @brief Phase where this solver lives. 
  */
  phase::which_phase m_phase;

  /*!
    @brief Solver name
  */
  std::string m_name;

  /*!
    @brief Class name
    @details ItoSolver for parent class -- derived classes might be named something else.
  */
  std::string m_className;

  /*!
    @brief Uniform integer distribution between [0,SpaceDim-1]
  */
  mutable std::uniform_int_distribution<int> m_uniformDistribution0d;

  /*!
    @brief Truncation value for normal distribution. 
    @details This is used when drawing numbers from the Gaussian distribution, replacing all values above m_normalDistributionTruncation by m_normalDistributionTruncation.
  */
  Real m_normalDistributionTruncation;

  /*!
    @brief Bisection step size to use for particle intersection tests with EBs.
  */
  Real m_bisectionStep;

  /*!
    @brief Verbosity level for this solver.
  */
  int m_verbosity;

  /*!
    @brief Time step
  */
  int m_timeStep;

  /*!
    @brief Size of refinement boundary halo
  */
  int m_haloBuffer;

  /*!
    @brief Current time
  */
  Real m_time;

  /*!
    @brief Time step
  */
  Real m_dt;

  /*!
    @brief NGP deposition in cut cells or not
  */
  bool m_forceIrregDepositionNGP;

  /*!
    @brief NGP interpolation in cut cells or not
  */
  bool m_forceIrregInterpolationNGP;

  /*!
    @brief Force usage of NGP when depositing "halo" particles. 
  */
  bool m_forceHaloNGP;

  /*!
    @brief Use mass redistribution or not. 
  */
  bool m_useRedistribution;

  /*!
    @brief Flag for blending the deposition clouds with the "non-conservative" divergence. 
  */
  bool m_blendConservation;

  /*!
    @brief If true, solver is diffusive. 
  */
  bool m_isDiffusive;

  /*!
    @brief If true, solver is mobile, i.e. advection is turned on. 
  */
  bool m_isMobile;

  /*!
    @brief Flag for outputting m_phi to plot files. 
  */
  bool m_plotPhi;

  /*!
    @brief Flag for outputting m_velocityFunction to plot files. 
  */
  bool m_plotVelocity;

  /*!
    @brief Flag for outputting m_diffusionFunction to plot files. 
  */
  bool m_plotDiffCo;

  /*!
    @brief Flag for depositing and plotting the bulk particles on the mesh. 
  */
  bool m_plotParticles;

  /*!
    @brief Flag for depositing and plotting the EB particles on the mesh. 
  */
  bool m_plotParticlesEB;

  /*!
    @brief Flag for depositing and plotting the domain particles on the mesh. 
  */
  bool m_plotParticlesDomain;

  /*!
    @brief Flag for depositing and plotting the source particles on the mesh. 
  */
  bool m_plotParticlesSource;

  /*!
    @brief Flag for depositing and plotting the covered particles on the mesh. 
  */
  bool m_plotParticlesCovered;

  /*!
    @brief Flag for plotting the energy density on the mesh
  */
  bool m_plotEnergyDensity;

  /*!
    @brief Flag for plotting the average particle energy on the mesh
  */
  bool m_plotAverageEnergy;

  /*!
    @brief Deposition method when depositing particles to the mesh. 
  */
  DepositionType m_deposition;

  /*!
    @brief Coarse-fine deposition strategy
  */
  CoarseFineDeposition m_coarseFineDeposition;

  /*!
    @brief Deposition method when depositing particles to mesh during plotting. 
  */
  DepositionType m_plotDeposition;

  /*!
    @brief Storage for particle mesh data. 
  */
  EBAMRCellData m_phi;

  /*!
    @brief Mobility function -- used when interpolating particle mobilities. 
  */
  EBAMRCellData m_mobilityFunction;

  /*!
    @brief Cell-centered field used for interpolating velocities
  */
  EBAMRCellData m_velocityFunction;

  /*!
    @brief Diffusion-centerer field used for interpolating diffusion coefficients
  */
  EBAMRCellData m_diffusionFunction;

  /*!
    @brief Scratch data
  */
  mutable EBAMRCellData m_scratch;

  /*!
    @brief Scratch storage for holding the non-conservative deposition
  */
  mutable EBAMRIVData m_depositionNC;

  /*!
    @brief Scratch storage for holding the mass difference when using hybrid deposition
  */
  mutable EBAMRIVData m_massDiff;

  /*!
    @brief Various particle containers with identifiers. 
  */
  std::map<WhichContainer, ParticleContainer<ItoParticle>> m_particleContainers;

  /*!
    @brief Parse RNG options -- this parses the RNG seed and instantiates the distributions. 
  */
  void
  parseRNG();

  /*!
    @brief Parse the normal distribution truncation level
  */
  void
  parseTruncation();

  /*!
    @brief Parse particle deposition methods.
  */
  void
  parseDeposition();

  /*!
    @brief Parse bisection step
  */
  void
  parseBisectStep();

  /*!
    @brief Parse class verbosity
  */
  void
  parseVerbosity();

  /*!
    @brief Parse plot variables
  */
  void
  parsePlotVariables();

  /*!
    @brief Parse diffusion hop
  */
  void
  parseDiffusionHop();

  /*!
    @brief Parse whether or not to use redistribution. 
  */
  void
  parseRedistribution();

  /*!
    @brief Parse whether or not to compute a "non-conservative" divergence when redistributing mass.
  */
  void
  parseDivergenceComputation();

  /*!
    @brief Parse checkpointing method
  */
  void
  parseCheckpointing();

#ifdef CH_USE_HDF5
  /*!
    @brief Write checkpoint data into HDF5 file -- this version dumps that particles to the HDF5 file.. 
    @details This will not write the full Ito partilces -- in order to trim memory and checkpoint file sizes this function writse
    the particle data to HDF5 as SimpleItoParticle which only stores position, mass, and energy. Things like mobility, diffusion coefficient, velocity etc.
    are parametrically dependent and are thus not stored. 
    @param[out] a_handle HDF5 file. 
    @param[in]  a_level Grid level
  */
  virtual void
  writeCheckPointLevelParticles(HDF5Handle& a_handle, const int a_level) const;
#endif

#ifdef CH_USE_HDF5
  /*!
    @brief Write checkpoint data into HDF5 file -- this version writes the number of particles per cell rather than the particles.
    @details This will write the number of particles in each grid cell -- this yields very-low-memory HDF5 files but when the simulation is restarted
    from this state then the particle distribution as well as the particle energy in each cell is lost. Use with caution!
    @param[out] a_handle HDF5 file. 
    @param[in]  a_level Grid level
  */
  virtual void
  writeCheckPointLevelFluid(HDF5Handle& a_handle, const int a_level) const;
#endif

#ifdef CH_USE_HDF5
  /*!
    @brief Read checkpointed particles from  an HDF5 file.
    @details This will read particles from the HDF5 file ala writeCheckPointLevelParticles. 
    @param[out] a_handle HDF5 file. 
    @param[in]  a_level Grid level
  */
  virtual void
  readCheckpointLevelParticles(HDF5Handle& a_handle, const int a_level);
#endif

#ifdef CH_USE_HDF5
  /*!
    @brief Read checkpointed particle numberse from  an HDF5 file and instantiate the particles from that. 
    @details This will read particle numbers in each cell ala writeCheckPointLevelFluid and instantiate the particles from that. 
    @param[out] a_handle HDF5 file. 
    @param[in]  a_level Grid level
  */
  virtual void
  readCheckpointLevelFluid(HDF5Handle& a_handle, const int a_level);
#endif

  /*!
    @brief Restart particles from a specified number of particles in the grid cell. 
    @details This will instantiate the bulk particles by randomly drawing new particles in each grid cell. 
    @param[in] a_particlesPerCell Mesh data showing how many particles are in each cell. 
    @param[in] a_level            Grid level. 
    @param[in] a_newPPC           Desired number of computational particles per cell
  */
  void
  drawNewParticles(const LevelData<EBCellFAB>& a_particlesPerCell, const int a_level, const int a_newPPC);

  /*!
    @brief Compute the cell-centered deposition -- this is the main deposition function.
    @details The template parameters indicate the particle type and quantity to be deposited. The second template parameter must be a pointer to a member function in 
    the particle class with signature 'const Real& P::function() const'. E.g. 'const Real& P::mass() const' which is the default quantity to be deposited. 
    @param[out] a_phi                  Cell-centered mesh data. Must have exactly one component. 
    @param[in]  a_particles            Particles to be deposited. 
    @param[in]  a_deposition           Deposition type
    @param[in]  a_coarseFineDeposition Coarse-fine deposition strategy
  */
  template <class P, const Real& (P::*particleScalarField)() const>
  void
  depositKappaConservative(EBAMRCellData&             a_phi,
                           ParticleContainer<P>&      a_particles,
                           const DepositionType       a_deposition,
                           const CoarseFineDeposition a_coarseFineDeposition) const;

  /*!
    @brief Compute the cell-centered deposition -- this is the main deposition function.
    @details Just like the version above except that the function signature is Real P::particleScalarField() const.
    @param[out] a_phi                  Cell-centered mesh data. Must have exactly one component. 
    @param[in]  a_particles            Particles to be deposited. 
    @param[in]  a_deposition           Deposition type
    @param[in]  a_coarseFineDeposition Coarse-fine deposition strategy
  */
  template <class P, Real (P::*particleScalarField)() const>
  void
  depositKappaConservative(EBAMRCellData&             a_phi,
                           ParticleContainer<P>&      a_particles,
                           const DepositionType       a_deposition,
                           const CoarseFineDeposition a_coarseFineDeposition) const;

  /*!
    @brief Redistribute mass in an AMR context. 
    @details We will have deposited particles into each cell, i.e. phi_i = m_i/dx^3. To obtain the true density in an EB context we need to divide by kappa such that 
    phi_i = m_i/(kappa_i*dx^3). Unfortunately, this is numerically unstable because kappa_i can be arbitrarily small. We can set the density as phi_i = m_i/dx^3 but
    we will be missing a mass m_i*(1-kappa_i) from the cell. This mass can be smooshed into neighboring cells such as to make the total scheme conservative. As an option,
    we can also use a "non-conservative" divergence but the scheme is not guaranteed to be non-negative. All of this is just to say that we take the mass that fell inside
    the EB and put it back into the domain. 
    @param[inout] a_phi Quantity to be redistributed. 
  */
  void
  redistributeAMR(EBAMRCellData& a_phi) const;

  /*!
    @brief Make the "non-conservative" kappa deposition -- computing depositionNC = sum(kappa*depositionC)/sum(kappa) in a neighborhood around each vof.
    @param[out] a_depositionNC     Non-conservative deposition
    @param[in]  a_depositionKappaC Conserved deposition
  */
  void
  depositNonConservative(EBAMRIVData& a_depositionNC, const EBAMRCellData& a_depositionKappaC) const;

  /*!
    @brief Make the "hybrid" deposition phiH = kappa*phiC + (1-kappa)*phiNC. On input, a_depositionH should contain phiC
    @param[inout] a_depositionH    On input, contains phiC. On output, contain phiH
    @param[out]   a_massDifference On output, contains mass loss in each cut-cell.
    @param[in]    a_depositionNC   The "non-conservative" deposited variable.
  */
  void
  depositHybrid(EBAMRCellData& a_depositionH, EBAMRIVData& a_massDifference, const EBAMRIVData& a_depositionNC) const;

  /*!
    @brief Increment redistribution objects with mass to be redistributed
    @param[in] a_massDifference mass to be redistributed. 
  */
  void
  incrementRedist(const EBAMRIVData& a_massDifference) const;

  /*!
    @brief Do redistribution on all levels, no coarse-fine stuff
    @param[inout] a_phi Data to be redistributed into. 
  */
  void
  levelRedist(EBAMRCellData& a_phi) const;

  /*!
    @brief Increment the coarse-fine redist with mass difference (only for EB x CF)
    @param[in] a_massDifference mass to be redistributed. 
  */
  void
  coarseFineIncrement(const EBAMRIVData& m_massDifference) const;

  /*!
    @brief Do coarse-fine mass redistribution
    @param[inout] a_phi Data to be redistributed into. 
  */
  void
  coarseFineRedistribution(EBAMRCellData& a_phi) const;

  /*!
    @brief Superparticle merging with BVH trees
    @param[inout] a_particles        Particles to be merged/split
    @param[in]    a_particlesPerCell Target number of particles per cell
  */
  void
  makeSuperparticles(List<ItoParticle>& a_particles, const int a_particlesPerCell);

  /*!
    @brief Write data to output. Convenience function which just copies data from one data holder to the output data holder. 
    @param[inout] a_output            Output data holder.
    @param[inout] a_icomp             Starting component where this solver begins writing the output. 
    @param[in]    a_data              Data to write. 
    @param[in]    a_interpToCentroids If true, a_data will be interpolated to cell centroids before writing to a_output. 
  */
  virtual void
  writeData(EBAMRCellData& a_output, int& a_comp, const EBAMRCellData& a_data, const bool a_interpToCentroids) const;

  /*!
    @brief Directly interpolate mobilities. Interpolates for all particles in the specified grid patch. 
    @param[in] a_level Grid level.
    @param[in] a_dit   Grid index.
  */
  virtual void
  interpolateMobilitiesDirect(const int a_level, const DataIndex& a_dit);

  /*!
    @brief Interpolate mobilities through the velocity. Interpolates for all particles in the specified grid patch. 
    @param[in] a_level Grid level.
    @param[in] a_dit   Grid index.
  */
  virtual void
  interpolateMobilitiesVelocity(const int a_level, const DataIndex& a_dit);
};

#include <CD_NamespaceFooter.H>

#include <CD_ItoSolverImplem.H>

#endif
