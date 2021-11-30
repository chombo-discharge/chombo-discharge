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

// Std includes
#include <random>

// Chombo includes
#include <Particle.H>
#include <ParticleData.H>
#include <BinItem.H>

// Our includes
#include <CD_AmrMesh.H>
#include <CD_ComputationalGeometry.H>
#include <CD_ItoSpecies.H>
#include <CD_ItoParticle.H>
#include <CD_EBParticleMesh.H>
#include <CD_ParticleContainer.H>
#include <CD_ItoMerge.H>
#include <CD_PointMass.H>
#include <CD_EbRepresentation.H>
#include <CD_NamespaceHeader.H>

  
/*!
  @brief Base class for Ito diffusion particle models.
  @details This is a particle class for solving Ito diffusion problems dX = v*dt + sqrt(2*D)*dW_t (dW_t is a Wiener process over dt) over an AMR hierarchy. The diffusion
  coefficient is related to the hydrodynamic diffusion coefficient by D(fluid) = sqrt(2*D(Ito)). 
  @note To use this class, the user must add a species and set the realm and phase where the solver will live. 
*/
class ItoSolver {
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
  enum class WhichContainer {
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
  virtual std::string getName() const;

  /*!
    @brief Get the realm where this solver is registered.
    @return m_realm
  */
  virtual const std::string getRealm() const;

  /*!
    @brief Set the realm where this solver will live
    @param[in] a_realm Realm name. 
  */
  virtual void setRealm(const std::string a_realm);

  /*!
    @brief Parse class options
  */
  virtual void parseOptions();

  /*!
    @brief Parse class runtime configurable options.
  */
  virtual void parseRuntimeOptions();

  /*!
    @brief Get output plot names
    @return Returns a string of all variables that this solver can add to plot files. 
  */
  virtual Vector<std::string> getPlotVariableNames() const;

  /*!
    @brief Set the initial data
    @details This will add the initial particles and deposit them. 
  */
  virtual void initialData();

  /*!
    @brief Regrid this solver. 
    @param[in] a_lmin           Coarsest level where grids did not change. 
    @param[in] a_oldFinestLevel Finest AMR level before the regrid. 
    @param[in] a_newFinestLevel Finest AMR level after the regrid. 
  */
  virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Allocate internal storage
    @details This will allocate the required mesh data as well as all the particle data holders. 
    @note Storage for diffusion and advection fields are only allocate if the species is diffusive/mobile. 
  */
  virtual void allocateInternals();

  /*!
    @brief Clear a particle container -- this will delete all the particles in the input container. 
    @param[in] a_whichContainer Particle container to be cleared.
  */
  virtual void clear(const WhichContainer a_container);

  /*!
    @brief Clear a particle container -- this will delete all the particles in the input container. 
    @param[in] a_particles Particle container to be emptied. 
  */
  virtual void clear(ParticleContainer<ItoParticle>& a_particles) const;

  /*!
    @brief Clear particles -- this will deleted all the particles in the input container. 
    @param[in] a_particles Particle container to be emptied. 
  */
  virtual void clear(AMRParticles<ItoParticle>& a_particles) const;

  /*!
    @brief Deposit conductivities (i.e. mass*mobility / volume)
    @details This deposits mass*mobility (not multiplied by charge)
    @param[out] a_phi       Mesh data
    @param[in]  a_particles Particle data
    @note Calls the other versions with a_deposition = m_deposition
  */
  virtual void depositConductivity(EBAMRCellData& a_phi, const ParticleContainer<ItoParticle>& a_particles) const;

  /*!
    @brief Deposit conductivities (i.e. mass*mobility / volume)
    @details This deposits mass*mobility (not multiplied by charge)
    @param[out] a_phi        Mesh data
    @param[in]  a_particles  Particle data
    @param[in]  a_deposition Deposition method
  */
  virtual void depositConductivity(EBAMRCellData& a_phi, const ParticleContainer<ItoParticle>& a_particles, const DepositionType a_deposition) const;

  /*!
    @brief Deposit diffusivity (i.e. mass*D/volume)
    @details This deposits mass*diffusion (not multiplied by charge)
    @param[out] a_phi       Mesh data
    @param[in]  a_particles Particle data
    @note Calls the other versions with a_deposition = m_deposition
  */
  virtual void depositDiffusivity(EBAMRCellData& a_phi, const ParticleContainer<ItoParticle>& a_particles) const;

  /*!
    @brief Deposit diffusivity (i.e. mass*D/volume)
    @details This deposits mass*mobility (not multiplied by charge)
    @param[out] a_phi        Mesh data
    @param[in]  a_particles  Particle data
    @param[in]  a_deposition Deposition method
  */
  virtual void depositDiffusivity(EBAMRCellData& a_phi, const ParticleContainer<ItoParticle>& a_particles, const DepositionType a_deposition) const;

  /*!
    @brief Deposit energy densities (i.e. mass*energy/volume => total energy per unit volume)
    @details This deposits mass*energy
    @param[out] a_phi       Mesh data
    @param[in]  a_particles Particle data
    @note Calls the other versions with a_deposition = m_deposition
  */
  virtual void depositEnergyDensity(EBAMRCellData& a_phi, const ParticleContainer<ItoParticle>& a_particles) const;

  /*!
    @brief Deposit energy densities (i.e. mass*energy/volume => total energy per unit volume)
    @details This deposits mass*energy
    @param[out] a_phi        Mesh data
    @param[in]  a_particles  Particle data
    @param[in]  a_deposition Deposition method
  */
  virtual void depositEnergyDensity(EBAMRCellData& a_phi, const ParticleContainer<ItoParticle>& a_particles, const DepositionType a_deposition) const;

  /*!
    @brief Compute average mobility. 
    @details This computes the average mobility as mu_avg = sum(mass*mu)/sum(mass) by depositing on the mesh. 
    @param[out] a_phi       Mesh data
    @param[in]  a_particles Particle data
  */
  virtual void computeAverageMobility(EBAMRCellData& a_phi, const ParticleContainer<ItoParticle>& a_particles) const;

  /*!
    @brief Compute average diffusion coefficient. 
    @details This computes the average diffusion coefficient as D_avg = sum(mass*D)/sum(mass) by depositing on the mesh. 
    @param[out] a_phi       Mesh data
    @param[in]  a_particles Particle data
  */
  virtual void computeAverageDiffusion(EBAMRCellData& a_phi, const ParticleContainer<ItoParticle>& a_particles) const;

  /*!
    @brief Compute average particle energy.. 
    @details This computes the average energy as E_avg = sum(mass*energy)/sum(mass) by depositing on the mesh. 
    @param[out] a_phi       Mesh data
    @param[in]  a_particles Particle data
  */
  virtual void computeAverageEnergy(EBAMRCellData& a_phi, const ParticleContainer<ItoParticle>& a_particles) const;

  /*!
    @brief Deposit particles onto mesh. 
    @details This will deposit the mass (i.e., computational weight) "bulk" particles into m_phi. 
    @note Cells the other version with a_container = WhichContainer::Bulk
  */
  virtual void depositParticles();

  /*!
    @brief Deposit particles on to mesh. 
    @param[in] a_container Which container to deposit.
    @details This will deposit mass (i.e., computational weight) of the the input particle container particles onto the classes member 'm_phi'. 
    @note Calls the general version with arguments: m_phi, m_particles.at(a_container), m_deposition.
  */
  virtual void depositParticles(const WhichContainer a_container);

  /*!
    @brief Generic deposition function which deposits a particle field onto the mesh using a specified deposition method. 
    @details The template parameters indicate the particle type and quantity to be deposited. The second template parameter must be a pointer to a member function in 
    the particle class with signature 'const Real& P::function() const'. E.g. 'const Real& P::mass() const' which is the default quantity to be deposited. 
    @param[out] a_phi        Mesh data -- must have exactly one component. 
    @param[in]  a_particles  Particles to be deposited
    @param[in]  a_deposition Deposition method
    @note The coarse-fine deposition handling is defined by the input parameters in m_pvrBuffer, m_haloBuffer, and m_forceIrregHaloNGP (see the class documentation).
  */
  template <class P, const Real&(P::*particleScalarField)() const>    
  void depositParticles(EBAMRCellData& a_phi, const ParticleContainer<P>& a_particles, const DepositionType a_deposition) const;

  /*!
    @brief Generic deposition function which deposits a particle field onto the mesh using a specified deposition method. 
    @details Just like the version above, but with a different function signature.
    @param[out] a_phi        Mesh data -- must have exactly one component. 
    @param[in]  a_particles  Particles to be deposited
    @param[in]  a_deposition Deposition method
    @note The coarse-fine deposition handling is defined by the input parameters in m_pvrBuffer, m_haloBuffer, and m_forceIrregHaloNGP (see the class documentation).
  */
  template <class P, Real(P::*particleScalarField)() const>
  void depositParticles(EBAMRCellData& a_phi, const ParticleContainer<P>& a_particles, const DepositionType a_deposition) const;

  /*!
    @brief Remove particles that are inside the EB. 
    @param[in] a_ebRepresentation EB representation (e.g. implicit function, discrete, voxel)
    @param[in] a_tol              Tolerance -- removes particles if they are less than a_tol*dx away from the boundary. 
    @details Calls the other version with m_particles.
  */
  virtual void removeCoveredParticles(const EbRepresentation a_representation, const Real a_tol);

  /*!
    @brief Remove particles that are inside the EB. 
    @param[in] a_container        Which particle container to remove particles from
    @param[in] a_ebRepresentation EB representation (e.g. implicit function, discrete, voxel)
    @param[in] a_tol              Tolerance -- removes particles if they are less than a_tol*dx away from the boundary. 
    @note Calls the other version with ParticleContainer = m_particles.at(a_container)
  */
  virtual void removeCoveredParticles(const WhichContainer a_container, const EbRepresentation a_representation, const Real a_tol);

  /*!
    @brief Remove particles that are inside the EB
    @param[in] a_particles        Particles to remove. 
    @param[in] a_ebRepresentation EB representation (e.g. implicit function, discrete, voxel)
    @param[in] a_tol              Tolerance -- removes particles if they are less than a_tol*dx away from the boundary. 
  */
  virtual void removeCoveredParticles(ParticleContainer<ItoParticle>& a_particles, const EbRepresentation a_representation, const Real a_tol) const;

  /*!
    @brief Transfer particles that are covered (within a_tol distance from EB) to another container. 
    @details This is the default version which moves particles from m_particles to m_covered_particles
  */
  virtual void transferCoveredParticles(const EbRepresentation a_representation, const Real a_tol);

  /*!
    @brief Which particles to remove
  */
  virtual void transferCoveredParticles(const WhichContainer   a_containerFrom,
					const WhichContainer   a_containerTo,
					const EbRepresentation a_representation,
					const Real              a_tol);

  /*!
    @brief Transfer particles that are inside the EB to another container. 
  */
  virtual void transferCoveredParticles(ParticleContainer<ItoParticle>& a_containerFrom,
					ParticleContainer<ItoParticle>& a_containerTo,
					const EbRepresentation           a_representation,
					const Real                        a_tol);
					    

  /*!
    @brief Do boundary intersection tests
  */
  virtual void intersectParticles(const EbRepresentation a_representation, const bool a_delete);

  /*!
    @brief Do boundary intersection tests.
  */
  virtual void intersectParticles(const WhichContainer       a_particles,
				  const WhichContainer       a_eb_particles,
				  const WhichContainer       a_dom_particles,
				  const EbRepresentation a_representation,				   
				  const bool              a_delete);

  /*!
    @brief Do boundary intersection tests.
  */
  virtual void intersectParticles(ParticleContainer<ItoParticle>& a_particles,
				  ParticleContainer<ItoParticle>& a_eb_particles,
				  ParticleContainer<ItoParticle>& a_dom_particles,
				  const EbRepresentation           a_representation,
				  const bool                        a_delete);

  /*!
    @brief Compute particle load on specific levle
  */
  virtual void computeLoads(Vector<long int>& a_loads, const DisjointBoxLayout& a_dbl, const int a_level);

  /*!
    @brief Get number of plot variables
    @return Returns the number of plot variables (Note: Scalars-> one variable and vectors -> SpaceDim variables).
  */
  virtual int getNumberOfPlotVariables() const;

  /*!
    @brief Get PVR buffer
  */
  virtual int getPVRBuffer() const;

  /*!
    @brief Get the halo buffer
  */
  virtual int getHaloBuffer() const;

  /*!
    @brief Set the PVR buffer
  */
  virtual void setPVRBuffer(const int a_buffer);

  /*!
    @brief Set the halo buffer
  */
  virtual void setHalobuffer(const int a_buffer);

  /*!
    @brief Get number of particles in a container
  */
  virtual size_t getNumParticles(const WhichContainer a_container, const bool a_local) const;

  /*!
    @brief Register operators
  */
  virtual void registerOperators();

  /*!
    @brief Set computational geometry
  */
  virtual void setComputationalGeometry(const RefCountedPtr<ComputationalGeometry> a_computationalGeometry);

  /*!
    @brief Set the amr object
  */
  virtual void setAmr(const RefCountedPtr<AmrMesh>& a_amr);

  /*!
    @brief Set phase
  */
  virtual void setPhase(phase::which_phase a_phase);

  /*!
    @brief Set verbosity.
  */
  virtual void setVerbosity(const int a_verbosity);

  /*!
    @brief Get current time
  */
  virtual void setTime(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Write plot data
  */
  virtual void writePlotData(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Write data to output. Convenience function. 
  */
  virtual void writeData(EBAMRCellData& a_output, int& a_comp, const EBAMRCellData& a_data, const bool a_interp);

  /*!
    @brief Write checkpoint data into handle
  */
#ifdef CH_USE_HDF5
  virtual void writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const;
#endif

  /*!
    @brief Checkpoint particles
  */
#ifdef CH_USE_HDF5  
  virtual void writeCheckPointLevelParticles(HDF5Handle& a_handle, const int a_level) const;
#endif

  /*!
    @brief Checkpoint particles
  */
#ifdef CH_USE_HDF5  
  virtual void writeCheckPointLevelFluid(HDF5Handle& a_handle, const int a_level) const;
#endif

  /*!
    @brief Read checkpoint data from handle
  */
#ifdef CH_USE_HDF5  
  virtual void readCheckpointLevel(HDF5Handle& a_handle, const int a_level);
#endif

  /*!
    @brief Set the species
  */
  virtual void setSpecies(RefCountedPtr<ItoSpecies> a_species);

  /*!
    @brief Cache the particles
  */
  virtual void preRegrid(const int a_base, const int a_oldFinestLevel);

  /*!
    @brief Add particles. Level and patch version. 
  */
  virtual void addParticles(ListBox<ItoParticle>& a_part, const int a_lvl, const DataIndex a_dit, const bool a_destructive);

  /*!
    @brief Check if solver is mobile
  */
  virtual bool isMobile() const;
  
  /*!
    @brief Check if solver is diffusive
  */
  virtual bool isDiffusive() const;

  /*!
    @brief Get a general particle container
  */
  virtual ParticleContainer<ItoParticle>& getParticles(const WhichContainer a_container);

  /*!
    @brief Get a general particle container
  */
  virtual const ParticleContainer<ItoParticle>& getParticles(const WhichContainer a_container) const;

  /*!
    @brief Get the species
  */
  virtual RefCountedPtr<ItoSpecies>& getSpecies();

  /*!
    @brief Get state
  */
  virtual EBAMRCellData& getPhi();

  /*!
    @brief Get cell-centered velocity
  */
  virtual EBAMRCellData& getVelocityFunction();

  /*!
    @brief Get cell-centered velocity
  */
  virtual EBAMRCellData& getDiffusionFunction();

  /*!
    @brief Get scratch storage
  */
  virtual EBAMRCellData& getScratch();

  /*!
    @brief Get mobility function
  */
  virtual EBAMRCellData& getMobilityFunction();

  /*!
    @brief Set diffusion coefficient
  */
  virtual void setDiffusionFunction(const Real a_diffusionCoefficient);

  /*!
    @brief Set velocity
  */
  virtual void setVelocityFunction(const RealVect a_vel);

  /*!
    @brief Sets constant mobility for all particles
  */
  virtual void set_mobility(const Real a_mobility);

  /*!
    @brief Interpolate velocities
  */
  virtual void interpolateVelocities();

  /*!
    @brief Interpolate velocities
  */
  virtual void interpolateVelocities(const int a_level, const DataIndex& a_dit);

  /*!
    @brief Interpolate mobilities
  */
  virtual void interpolateMobilities();

  /*!
    @brief Interpolate mobilities
  */
  virtual void interpolateMobilities(const int a_level, const DataIndex& a_dit);

  /*!
    @brief Update mobilities
    @details This calls the mobility function implemented by ItoSpecies
  */
  virtual void updateMobilities();

  /*!
    @brief Update mobilities
    @details This calls the mobility function implemented by ItoSpecies
  */
  virtual void updateMobilities(const int a_level, const DataIndex a_dit);

  /*!
    @brief Interpolate diffusion
  */
  virtual void interpolateDiffusion();

  /*!
    @brief Interpolate velocities
  */
  virtual void interpolateDiffusion(const int a_level, const DataIndex& a_dit);

  /*!
    @brief Update diffusion coefficients
    @details This calls the diffusion function implemented by ItoSpecies
  */
  virtual void updateDiffusion();

  /*!
    @brief Update diffusion coefficients
    @details This calls the diffusion function implemented by ItoSpecies
  */
  virtual void updateDiffusion(const int a_level, const DataIndex a_dit);

  /*!
    @brief Sort the m_particles container by cell
  */
  virtual void sortParticlesByCell(const WhichContainer a_container);

  /*!
    @brief Sort the m_particles container by cell
  */
  virtual void sortParticlesByPatch(const WhichContainer a_container);
  
  /*!
    @brief Make superparticles for a full container -- this is the AMR version that users will usually call. 
    @param[in] a_container        Which container to repartition into new superparticles
    @param[in] a_particlesPerCell Target number of particles per cell
  */
  virtual void makeSuperparticles(const WhichContainer a_container, const int a_particlesPerCell);

  /*!
    @brief Make superparticles for a particle container -- this is used for a specified grid level. 
    @param[in] a_container        Which container to repartition into new superparticles
    @param[in] a_particlesPerCell Target number of particles per cell
    @param[in] a_level            Grid level
  */
  virtual void makeSuperparticles(const WhichContainer a_container, const int a_particlesPerCell, const int a_level);

  /*!
    @brief Make superparticles for a particle container -- this is used for a specified grid level and patch
    @param[in] a_container        Which container to repartition into new superparticles
    @param[in] a_particlesPerCell Target number of particles per cell
    @param[in] a_level            Grid level
    @param[in] a_dit              Grid index
  */
  virtual void makeSuperparticles(const WhichContainer a_container, const int a_particlesPerPatch, const int a_level, const DataIndex a_dit);

  /*!
    @brief Remap all particles
  */
  void remap();

  /*!
    @brief Remap all particles
  */
  void remap(const WhichContainer a_container);

  /*!
    @brief Random Gaussian
  */
  inline
  RealVect randomGaussian();

  /*!
    @brief Random direction
  */
  inline
  RealVect randomDirection();

  /*!
    @brief Compute min dt
  */
  Real computeDt() const;

  /*!
    @brief Compute drift dt. This returns the minimum of (dx/dt)*a_maxCellsToMove over all levels. 
  */
  Real computeDt(const int a_lvl) const;

  /*!
    @brief Compute drift dt. This returns the minimum of (dx/dt)*a_maxCellsToMove over all levels. 
  */
  Real computeDt(const int a_lvl, const DataIndex a_dit, const Real a_dx) const;

  /*!
    @brief Compute drift dt. This returns the minimum of (dx/dt)*a_maxCellsToMove over all levels. 
  */
  Real computeMinDt(const Real a_maxCellsToMove) const;

  /*!
    @brief Compute drift dt. This returns the minimum of (dx/dt)*a_maxCellsToMove over all levels. 
  */
  Real computeMinDt(const Real a_maxCellsToMove, const int a_lvl) const;

  /*!
    @brief Compute drift dt. This returns the minimum of (dx/dt)*a_maxCellsToMove over all levels. 
  */
  Real computeMinDt(const Real a_maxCellsToMove, const int a_lvl, const DataIndex a_dit, const Real a_dx) const;

  /*!
    @brief Compute drift dt. This returns the minimum of (dx/dt)*a_maxCellsToMove over all levels. 
  */
  Real computeMinDriftDt(const Real a_maxCellsToMove) const;

  /*!
    @brief Compute advection time step
  */
  Real computeAdvectiveDt() const;

  /*!
    @brief Compute the drift dt. This returns (dx/dt)*a_maxCellsToMove
  */
  Vector<Real> computeDriftDt(const Real a_maxCellsToMove) const;

  /*!
    @brief Compute the drift dt. This computes the minimum dt = dx/v on every level
  */
  Vector<Real> computeDriftDt() const;

  /*!
    @brief Compute the drift dt. This computes the minimum dt = dx/v on one level
  */
  Real computeDriftDt(const int a_lvl) const;

  /*!
    @brief Compute the drift dt. This computes the minimum dt = dx/v on one level and one box. 
  */
  Real computeDriftDt(const int a_lvl, const DataIndex& a_dit, const RealVect a_dx) const;

  /*!
    @brief Compute diffusion dt
  */
  Real computeMinDiffusionDt(const Real a_maxCellsToMove) const;

  /*!
    @brief Compute diffusion dt as dx*maxCellsToMove/(sqrt(2D)*eps) on all levels
  */
  Vector<Real> computeDiffusionDt(const Real a_maxCellsToMove) const;

  /*!
    @brief Compute the diffusive dt
  */
  Real computeDiffusiveDt() const;
  
  /*!
    @brief Compute diffusion dt as dx/(sqrt(2D)*eps) on all levels
  */
  Vector<Real> computeDiffusionDt() const;

  /*!
    @brief Compute diffusion dt as dx/(sqrt(2D)*eps) on one level
  */
  Real computeDiffusionDt(const int a_lvl) const;

  /*!
    @brief Compute diffusion dt as dx/(sqrt(2D)*eps) on one level and one patch
  */
  Real computeDiffusionDt(const int a_lvl, const DataIndex& a_dit, const RealVect a_dx) const;

  /*!
    @brief Get deposition method
  */
  DepositionType getDeposition() const;

  /*!
    @brief Return phase
  */
  phase::which_phase getPhase() const;
  
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
  */
  enum class WhichCheckpoint {
    particles,
    numbers
  };

  /*!
    @brief Enum for deciding how to interpolate mobilities. 
  */
  enum class WhichMobilityInterpolation {
    mobility,
    velocity
  };

  /*!
    @brief How to checkpoint files. particles => write particles to HDF5. numbers => write numbers to HDF5
  */
  WhichCheckpoint m_checkpointing;

  /*!
    @brief Switch for deciding how to interpolate mobilities, i.e. interpolating either mu*E or just mu (to the particle position)
  */
  WhichMobilityInterpolation m_WhichMobilityInterpolation;

  /*!
    @brief Get number of particles used for restarting
  */
  int m_restartPPC;

  /*!
    @brief Splitting direction
  */
  int m_directionKD;

  /*!
    @brief Realm
  */
  std::string m_realm;

  /*!
    @brief Boundary volume hierarchy merging tree
  */
  ItoMerge::Tree<PointMass> m_mergeTree;

  /*!
    @brief Computational geometry. 
  */
  RefCountedPtr<ComputationalGeometry> m_computationalGeometry;

  /*!
    @brief AMR; needed for grid stuff
  */
  RefCountedPtr<AmrMesh> m_amr;

  /*!
    @brief Species
  */
  RefCountedPtr<ItoSpecies> m_species;

  /*!
    @brief Phase
  */
  phase::which_phase m_phase;

  /*!
    @brief Name
  */
  std::string m_name;

  /*!
    @brief Class name
  */
  std::string m_className;

  /*!
    @brief Random number generator engine.
  */
  mutable std::mt19937_64 m_rng;

  /*!
    @brief Uniform real distribution on [0,1]
  */
  mutable std::uniform_real_distribution<Real> m_uniformDistribution01; 

  /*!
    @brief Uniform real distribution on [-1,1]
  */
  mutable std::uniform_real_distribution<Real> m_uniformDistribution11;

  /*!
    @brief Gaussian distribution centered on 0, std=1
  */
  mutable std::normal_distribution<Real> m_normalDistribution01;

  /*!
    @brief Uniform integer distribution between [0,SpaceDim-1]
  */
  mutable std::uniform_int_distribution<int> m_uniformDistribution0d;

  /*!
    @brief Maximum exponential
  */
  Real m_normalDistributionTruncation;

  /*!
    @brief Random number generator seed.
  */
  Real m_rngSeed;

  /*!
    @brief Bisection step for particle trajectories
  */
  Real m_bisectionStep;

  /*!
    @brief Verbosity
  */
  int m_verbosity;

  /*!
    @brief Time step
  */
  int m_timeStep;

  /*!
    @brief Buffer for PVR
  */
  int m_pvrBuffer;

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
    @brief Length of maximum diffusion hop relative to the cell size
  */
  Real m_maxDiffusionHop;

  /*!
    @brief NGP deposition in cut cells or not
  */
  bool m_forceIrregDepositionNGP;

  /*!
    @brief NGP interpolation in cut cells or not
  */
  bool m_forceIrregInterpolationNGP;

  /*!
    @brief Use ngp when depositing halo particles or not
  */
  bool m_forceHaloNGP;

  /*!
    @brief Redistribute or not
  */
  bool m_useRedistribution;

  /*!
    @brief Flag for blending the deposition clouds with the nonconservative divergence
  */
  bool m_blendConservation;

  /*!
    @brief Solver is diffusive
  */
  bool m_isDiffusive;

  /*!
    @brief Also solve for advection/convection
  */
  bool m_isMobile;

  /*!
    @brief Output state
  */
  bool m_plotPhi;

  /*!
    @brief Output state
  */
  bool m_plotVelocity;

  /*!
    @brief Output state
  */
  bool m_plotDiffCo;

  /*!
    @brief Plot particles
  */
  bool m_plotParticles;

  /*!
    @brief Plot EB particles
  */
  bool m_plotParticlesEB;

  /*!
    @brief Plot domain particles
  */
  bool m_plotParticlesDomain;

  /*!
    @brief Plot source particles
  */
  bool m_plotParticlesSource;

  /*!
    @brief Plot covered particles
  */
  bool m_plotParticlesCovered;

  /*!
    @brief Plot energy density
  */
  bool m_plotEnergyDensity;

  /*!
    @brief Plot average energy or not. 
  */
  bool m_plotAverageEnergy;

  /*!
    @brief Mass/charge deposition type
  */
  DepositionType m_deposition; 

  /*!
    @brief Plot deposition type
  */
  DepositionType m_plotDeposition; 

  /*!
    @brief Deposited particles
  */
  EBAMRCellData m_phi;

  /*!
    @brief Mobility function
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
  std::map<WhichContainer, ParticleContainer<ItoParticle> > m_particleContainers;

  /*!
    @brief Parse superparticles
  */
  void parseSuperParticles();
  
  /*!
    @brief Parse RNG options
  */
  void parseRng();

  /*!
    @brief Parse deposition options
  */
  void parseDeposition();

  /*!
    @brief Parse bisection step
  */
  void parseBisectStep();

  /*!
    @brief Parse PVR buffer
  */
  void parsePvrBuffer();

  /*!
    @brief Parse plot variables
  */
  void parsePlotVariables();

  /*!
    @brief Parse diffusion hop
  */
  void parseDiffusionHop();

  /*!
    @brief Parse the conservation
  */
  void parseRedistribution();

  /*!
    @brief Parse the conservation
  */
  void parseDivergenceComputation();

  /*!
    @brief Parse checkpointing method
  */
  void parseCheckpointing();

  /*!
    @brief Compute the cell-centered deposition -- this is the main deposition function.
    @details The template parameters indicate the particle type and quantity to be deposited. The second template parameter must be a pointer to a member function in 
    the particle class with signature 'const Real& P::function() const'. E.g. 'const Real& P::mass() const' which is the default quantity to be deposited. 
    @param[out] a_phi        Cell-centered mesh data. Must have exactly one component. 
    @param[in]  a_particles  Particles to be deposited. 
    @param[in]  a_deposition Deposition type
    @note The coarse-fine deposition is handled internally in this routine (through parameters m_pvrBuffer and m_haloBuffer)
  */
  template <class P,  const Real&(P::*particleScalarField)() const>
  void depositKappaConservative(EBAMRCellData&              a_phi,
				const ParticleContainer<P>& a_particles,
				const DepositionType        a_deposition) const;

  /*!
    @brief Compute the cell-centered deposition -- this is the main deposition function.
    @details Just like the version above except that the function signature is Real P::particleScalarField() const.
    @param[out] a_phi        Cell-centered mesh data. Must have exactly one component. 
    @param[in]  a_particles  Particles to be deposited. 
    @param[in]  a_deposition Deposition type
    @note The coarse-fine deposition is handled internally in this routine (through parameters m_pvrBuffer and m_haloBuffer)
  */
  template <class P,  Real(P::*particleScalarField)() const>
  void depositKappaConservative(EBAMRCellData&              a_phi,
				const ParticleContainer<P>& a_particles,
				const DepositionType        a_deposition) const;  

  /*!
    @brief Redistribution method
  */
  void redistributeAMR(EBAMRCellData& a_phi) const;

  /*!
    @brief Make the "non-conservative" kappa deposition
  */
  void depositNonConservative(EBAMRIVData& a_depositionNC, const EBAMRCellData& a_depositionKappaC) const;

  /*!
    @brief Make the hybrid deposition. Also compute the mass difference
  */
  void depositHybrid(EBAMRCellData& a_depositionH, EBAMRIVData& a_massDifference, const EBAMRIVData& a_depositionNC) const;

  /*!
    @brief Increment redistribution objects with mass to be redistributed
  */
  void incrementRedist(const EBAMRIVData& a_massDifference) const;

  /*!
    @brief Do redistribution on all levels, no coarse-fine stuff
  */
  void levelRedist(EBAMRCellData& a_phi) const;

  /*!
    @brief Increment the coarse-fine redist with mass difference (only for EB x CF)
  */
  void coarseFineIncrement(const EBAMRIVData& m_massDifference) const;

  /*!
    @brief Do coarse-fine redistribution
  */
  void coarseFineRedistribution(EBAMRCellData& a_phi) const;

  /*!
    @brief Superparticle merging with BVH trees
  */
  void mergeBVH(List<ItoParticle>& a_particles, const int a_particlesPerCell);

  /*!
    @brief Restart particles
  */
  void restartParticles(LevelData<EBCellFAB>& a_num_particles, const int a_level);

  /*!
    @brief Transfer particles that are inside the EB
  */
  void transferCoveredParticles_if(ParticleContainer<ItoParticle>& a_src, ParticleContainer<ItoParticle>& a_dst, const Real a_tol);

  /*!
    @brief Do boundary intersection tests with an implicit function
  */
  void intersectParticlesIF(ParticleContainer<ItoParticle>& a_particles,
			    ParticleContainer<ItoParticle>& a_eb_particles,
			    ParticleContainer<ItoParticle>& a_domain_particles,
			    const bool                        a_delete);

  /*!
    @brief Interpolate mobilities
  */
  virtual void interpolateMobilitiesMu(const int a_level, const DataIndex& a_dit);

  /*!
    @brief Interpolate mobilities
  */
  virtual void interpolateMobilitiesVel(const int a_level, const DataIndex& a_dit);

  /*!
    @brief Sign function
  */
  inline
  Real sign(const Real& a) const;

  /*!
    @brief Draw a valid random position somewhere in a cell
  */
  inline
  RealVect randomPosition(const RealVect a_pos,
			  const RealVect a_lo,
			  const RealVect a_hi,
			  const RealVect a_bndryCentroid,
			  const RealVect a_normal,
			  const Real     a_dx,
			  const Real     a_kappa) const;

  /*!
    @brief Draw a random position somewhere in a cell
  */
  inline
  RealVect randomPosition(const RealVect a_lo,
			  const RealVect a_hi,
			  const RealVect a_bndryCentroid,
			  const RealVect a_normal) const;

  /*!
    @brief Return a random position in the cube (a_lo, a_hi);
  */
  inline
  RealVect randomPosition(const RealVect a_lo, const RealVect a_hi) const;
};

#include <CD_NamespaceFooter.H>

#include <CD_ItoSolverImplem.H>

#endif
