/*!
  @file   ito_solver.H
  @brief  Declaration of an abstract class for Ito diffusion
  @author Robert Marskar
  @date   April 2020
*/

#ifndef _ITO_SOLVER_
#define _ITO_SOLVER_

#include "CD_AmrMesh.H"
#include <CD_ComputationalGeometry.H>
#include <CD_ItoSpecies.H>
#include <CD_ItoParticle.H>
#include <CD_EbParticleInterp.H>
#include <CD_ParticleContainer.H>
#include "bvh.H"
#include <CD_PointMass.H>
#include <CD_EbRepresentation.H>

#include <Particle.H>
#include <ParticleData.H>
#include <BinItem.H>

#include <random>

#define BOOST_RANDOM 0

#if BOOST_RANDOM
#include <boost/random.hpp>
namespace RAN = boost::random;
#else
namespace RAN = std;
#endif

#include "CD_NamespaceHeader.H"
  
/*!
  @brief Base class for Ito diffusion particle models
  @details This is a particle class for solving Ito diffusion problems dX = vdt + Z*D*t over an AMR hierarchy. The diffusion
  coefficient is related to the hydrodynamic diffusion coefficient by D(fluid = sqrt(2*D(Ito))
*/
class ito_solver {
public:

  /*!
    @brief Enum class for interfacing into individual containers. 
  */
  enum class which_container {
    bulk,
    eb,
    domain,
    source,
    covered,
    scratch
  };

  /*!
    @brief Default ctor
  */
  ito_solver();

  /*!
    @Default dtor
  */
  virtual ~ito_solver();

  /*!
    @brief Get solver name
  */
  virtual std::string getName();

  /*!
    @brief Get Realm
  */
  virtual const std::string getRealm() const;

  /*!
    @brief Set Realm
  */
  virtual void setRealm(const std::string a_realm);

  /*!
    @brief Parse options
  */
  virtual void parseOptions();

  /*!
    @brief Parse options
  */
  virtual void parseRuntimeOptions();

  /*!
    @brief Get output plot names
  */
  virtual Vector<std::string> getPlotVariableNames() const;

  /*!
    @brief Set the initial data
  */
  virtual void initialData();

  /*!
    @brief Regrid solver
  */
  virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Allocate internal storage
  */
  virtual void allocateInternals();

  /*!
    @brief Clear a particle container
  */
  virtual void clear(const which_container a_container);

  /*!
    @brief Empty particle container
  */
  virtual void clear(ParticleContainer<ItoParticle>& a_particles);

  /*!
    @brief Empty particle container
  */
  virtual void clear(AMRParticles<ItoParticle>& a_particles);

  /*!
    @brief Set mass to conductivity
    @details This sets mass = mass*mobility
  */
  virtual void set_mass_to_conductivity(ParticleContainer<ItoParticle>& a_particles);

  /*!
    @brief Unset mass to conductivity
    @details This sets mass = mass/mobility
  */
  virtual void unset_mass_to_conductivity(ParticleContainer<ItoParticle>& a_particles);

  /*!
    @brief Set mass to diffusivity
    @details This sets mass = mass*diffusion_coefficient
  */
  virtual void set_mass_to_diffusivity(ParticleContainer<ItoParticle>& a_particles);

  /*!
    @brief Unset mass to diffusivity
    @details This sets mass = mass/diffusion_coefficient
  */
  virtual void unset_mass_to_diffusivity(ParticleContainer<ItoParticle>& a_particles);

  /*!
    @brief Set mass to energy
    @details This sets mass = mass*particle_energy
  */
  virtual void set_mass_to_energy(ParticleContainer<ItoParticle>& a_particles);

  /*!
    @brief Unset mass to energy
    @details This sets mass = mass/particle_energy
  */
  virtual void unset_mass_to_energy(ParticleContainer<ItoParticle>& a_particles);

  /*!
    @brief Deposit conductivities
    @details This deposits weight*mobility (not multiplied by charge)
  */
  virtual void deposit_conductivity();

  /*!
    @brief Deposit conductivities
    @details This deposits weight*mobility (not multiplied by charge)
  */
  virtual void deposit_conductivity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles);

  /*!
    @brief Deposit conductivities
    @details This deposits weight*mobility (not multiplied by charge)
  */
  virtual void deposit_conductivity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles, const DepositionType::Which a_deposition);

  /*!
    @brief Deposit diffusivity
    @details This deposits weight*D (not multiplied by charge)
  */
  virtual void deposit_diffusivity();

  /*!
    @brief Deposit diffusivities
    @details This deposits weight*D (not multiplied by charge)
  */
  virtual void deposit_diffusivity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles);

  /*!
    @brief Deposit diffusivity
    @details This deposits weight*D (not multiplied by charge)
  */
  virtual void deposit_diffusivity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles, const DepositionType::Which a_deposition);

  /*!
    @brief Deposit diffusivity
    @details This deposits weight*energy
  */
  virtual void deposit_energy_density();

  /*!
    @brief Deposit energy
    @details This deposits weight*energy
  */
  virtual void deposit_energy_density(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles);

  /*!
    @brief Deposit energy
    @details This deposits weight*energy
  */
  virtual void deposit_energy_density(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles, const DepositionType::Which a_deposition);

  /*!
    @brief Compute average mobility. Must sort the particles by cell first. 
    
  */
  virtual void compute_average_mobility(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles);

  /*!
    @brief Compute average diffusion. Must sort the particles by cell first. 
  */
  virtual void compute_average_diffusion(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles);

  /*!
    @brief Compute average energy. Must sort the particles by cell first. 
  */
  virtual void compute_average_energy(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles);

  /*!
    @brief Deposit particles
  */
  virtual void deposit_particles();

  /*!
    @brief Deposit particles
  */
  virtual void deposit_particles(const which_container a_container);

  /*!
    @brief Deposit the particles onto a_phi. The result is put on the a_comp compoennt
  */
  template <class T>
  void deposit_particles(EBAMRCellData& a_phi, const ParticleContainer<T>& a_particles);

  /*!
    @brief Deposit the particles onto a_phi. The result is put on the a_comp compoennt
  */
  template <class T>
  void deposit_particles(EBAMRCellData& a_phi, const ParticleContainer<T>& a_particles, const DepositionType::Which  a_deposition);

  /*!
    @brief Deposit the total weight of the particles in the grid cells
  */
  virtual void deposit_weights(EBAMRCellData& a_weights, const ParticleContainer<ItoParticle>& a_particles);

  /*!
    @brief Remove particles that are inside the EB
    @details Calls the other version with m_particles
  */
  virtual void remove_covered_particles(const EbRepresentation a_representation, const Real a_tol);

  /*!
    @brief Which particles to remove
  */
  virtual void remove_covered_particles(const which_container a_container, const EbRepresentation a_representation, const Real a_tol);

  /*!
    @brief Remove particles that are inside the EB
  */
  virtual void remove_covered_particles(ParticleContainer<ItoParticle>& a_particles, const EbRepresentation a_representation, const Real a_tol);

  /*!
    @brief Transfer particles that are covered (within a_tol distance from EB) to another container. 
    @details This is the default version which moves particles from m_particles to m_covered_particles
  */
  virtual void transfer_covered_particles(const EbRepresentation a_representation, const Real a_tol);

  /*!
    @brief Which particles to remove
  */
  virtual void transfer_covered_particles(const which_container   a_containerFrom,
					  const which_container   a_containerTo,
					  const EbRepresentation a_representation,
					  const Real              a_tol);

  /*!
    @brief Transfer particles that are inside the EB to another container. 
  */
  virtual void transfer_covered_particles(ParticleContainer<ItoParticle>& a_containerFrom,
					  ParticleContainer<ItoParticle>& a_containerTo,
					  const EbRepresentation           a_representation,
					  const Real                        a_tol);
					    

  /*!
    @brief Do boundary intersection tests
  */
  virtual void intersect_particles(const EbRepresentation a_representation, const bool a_delete);

  /*!
    @brief Do boundary intersection tests.
  */
  virtual void intersect_particles(const which_container       a_particles,
				   const which_container       a_eb_particles,
				   const which_container       a_dom_particles,
				   const EbRepresentation a_representation,				   
				   const bool              a_delete);

  /*!
    @brief Do boundary intersection tests.
  */
  virtual void intersect_particles(ParticleContainer<ItoParticle>& a_particles,
				   ParticleContainer<ItoParticle>& a_eb_particles,
				   ParticleContainer<ItoParticle>& a_dom_particles,
				   const EbRepresentation           a_representation,
				   const bool                        a_delete);

  /*!
    @brief Compute particle load on specific levle
  */
  virtual void compute_loads(Vector<long int>& a_loads, const DisjointBoxLayout& a_dbl, const int a_level);

  /*!
    @brief Get number of plot variables
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
  virtual size_t get_num_particles(const which_container a_container, const bool a_local) const;

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
  virtual void writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const;

  /*!
    @brief Checkpoint particles
  */
  virtual void writeCheckpointLevel_particles(HDF5Handle& a_handle, const int a_level) const;

  /*!
    @brief Checkpoint particles
  */
  virtual void writeCheckpointLevel_fluid(HDF5Handle& a_handle, const int a_level) const;

  /*!
    @brief Read checkpoint data from handle
  */
  virtual void readCheckpointLevel(HDF5Handle& a_handle, const int a_level);

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
  virtual ParticleContainer<ItoParticle>& getParticles(const which_container a_container);

  /*!
    @brief Get a general particle container
  */
  virtual const ParticleContainer<ItoParticle>& getParticles(const which_container a_container) const;

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
  virtual EBAMRCellData& get_velo_func();

  /*!
    @brief Get cell-centered velocity
  */
  virtual EBAMRCellData& get_diffco_func();

  /*!
    @brief Get scratch storage
  */
  virtual EBAMRCellData& get_scratch();

  /*!
    @brief Get mobility function
  */
  virtual EBAMRCellData& get_mobility_func();

  /*!
    @brief Set diffusion coefficient
  */
  virtual void setDiffusionCoefficient_func(const Real a_diffusionCoefficient);

  /*!
    @brief Set velocity
  */
  virtual void setVelocity_func(const RealVect a_vel);

  /*!
    @brief Sets constant mobility for all particles
  */
  virtual void set_mobility(const Real a_mobility);

  /*!
    @brief Interpolate velocities
  */
  virtual void interpolate_velocities();

  /*!
    @brief Interpolate velocities
  */
  virtual void interpolate_velocities(const int a_level, const DataIndex& a_dit);

  /*!
    @brief Interpolate mobilities
  */
  virtual void interpolate_mobilities();

  /*!
    @brief Interpolate mobilities
  */
  virtual void interpolate_mobilities(const int a_level, const DataIndex& a_dit);

  /*!
    @brief Update mobilities
    @details This calls the mobility function implemented by ItoSpecies
  */
  virtual void update_mobilities();

  /*!
    @brief Update mobilities
    @details This calls the mobility function implemented by ItoSpecies
  */
  virtual void update_mobilities(const int a_level, const DataIndex a_dit);

  /*!
    @brief Interpolate diffusion
  */
  virtual void interpolate_diffusion();

  /*!
    @brief Interpolate velocities
  */
  virtual void interpolate_diffusion(const int a_level, const DataIndex& a_dit);

  /*!
    @brief Update diffusion coefficients
    @details This calls the diffusion function implemented by ItoSpecies
  */
  virtual void update_diffusion();

  /*!
    @brief Update diffusion coefficients
    @details This calls the diffusion function implemented by ItoSpecies
  */
  virtual void update_diffusion(const int a_level, const DataIndex a_dit);

  /*!
    @brief Sort the m_particles container by cell
  */
  virtual void sortParticlesByCell(const which_container a_container);

  /*!
    @brief Sort the m_particles container by cell
  */
  virtual void sortParticlesByPatch(const which_container a_container);
  
  /*!
    @brief Make superparticles
  */
  virtual void make_superparticles(const which_container a_container, const int a_particlesPerPatch);

  /*!
    @brief Make superparticles on specified level
  */
  virtual void make_superparticles(const which_container a_container, const int a_particlesPerPatch, const int a_level);

  /*!
    @brief Make superparticles on specified level and patch
  */
  virtual void make_superparticles(const which_container a_container, const int a_particlesPerPatch, const int a_level, const DataIndex a_dit);

  /*!
    @brief Remap all particles
  */
  void remap();

  /*!
    @brief Remap all particles
  */
  void remap(const which_container a_container);

  /*!
    @brief Random Gaussian
  */
  inline
  RealVect random_gaussian();

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
  Real compute_min_dt(const Real a_maxCellsToMove) const;

  /*!
    @brief Compute drift dt. This returns the minimum of (dx/dt)*a_maxCellsToMove over all levels. 
  */
  Real compute_min_dt(const Real a_maxCellsToMove, const int a_lvl) const;

  /*!
    @brief Compute drift dt. This returns the minimum of (dx/dt)*a_maxCellsToMove over all levels. 
  */
  Real compute_min_dt(const Real a_maxCellsToMove, const int a_lvl, const DataIndex a_dit, const Real a_dx) const;

  /*!
    @brief Compute drift dt. This returns the minimum of (dx/dt)*a_maxCellsToMove over all levels. 
  */
  Real compute_min_drift_dt(const Real a_maxCellsToMove) const;

  /*!
    @brief Compute advection time step
  */
  Real compute_advective_dt() const;

  /*!
    @brief Compute the drift dt. This returns (dx/dt)*a_maxCellsToMove
  */
  Vector<Real> compute_drift_dt(const Real a_maxCellsToMove) const;

  /*!
    @brief Compute the drift dt. This computes the minimum dt = dx/v on every level
  */
  Vector<Real> compute_drift_dt() const;

  /*!
    @brief Compute the drift dt. This computes the minimum dt = dx/v on one level
  */
  Real compute_drift_dt(const int a_lvl) const;

  /*!
    @brief Compute the drift dt. This computes the minimum dt = dx/v on one level and one box. 
  */
  Real compute_drift_dt(const int a_lvl, const DataIndex& a_dit, const RealVect a_dx) const;

  /*!
    @brief Compute diffusion dt
  */
  Real compute_min_diffusion_dt(const Real a_maxCellsToMove) const;

  /*!
    @brief Compute diffusion dt as dx*maxCellsToMove/(sqrt(2D)*eps) on all levels
  */
  Vector<Real> computeDiffusionDt(const Real a_maxCellsToMove) const;

  /*!
    @brief Compute the diffusive dt
  */
  Real compute_diffusive_dt() const;
  
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
  DepositionType::Which getDeposition() const;

  /*!
    @brief Return phase
  */
  phase::which_phase getPhase() const;
  
protected:

  /*! 
    @brief How to checkpoint files
  */
  enum class which_checkpoint {
    particles,
    numbers
  };

  /*!
    @brief Enum for deciding how to interpolate mobilities. 
  */
  enum class mobility_interp {
    mobility,
    velocity
  };

  /*!
    @brief How to checkpoint files. particles => write particles to HDF5. numbers => write numbers to HDF5
  */
  which_checkpoint m_checkpointing;

  /*!
    @brief Switch for deciding how to interpolate mobilities, i.e. interpolating either mu*E or just mu (to the particle position)
  */
  mobility_interp m_mobility_interp;

  /*!
    @brief Get number of particles used for restarting
  */
  int m_ppc_restart;

  /*!
    @brief Splitting direction
  */
  int m_kd_direction;

  /*!
    @brief Realm
  */
  std::string m_realm;

  /*!
    @brief Fluid Realm
    @note This is for later, in case we want to move redistribution and averaging/ghost cell interpolating onto a different Realm. 
  */
  std::string m_fluid_Realm;

  /*!
    @brief Boundary volume hierarchy merging tree
  */
  bvh_tree<PointMass> m_tree;

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
    @brief RNG engine
  */
  //  std::mt19937_64 m_rng;
  RAN::mt19937_64 m_rng;

  /*!
    @brief Uniform real distribution on [0,1]
  */
  RAN::uniform_real_distribution<Real> m_udist01; 

  /*!
    @brief Uniform real distribution on [-1,1]
  */
  RAN::uniform_real_distribution<Real> m_udist11;

  /*!
    @brief Gaussian distribution centered on 0, std=1
  */
  RAN::normal_distribution<Real> m_gauss01;

  /*!
    @brief Uniform integer distribution between [0,SpaceDim-1]
  */
  RAN::uniform_int_distribution<int> m_udist0d;

  /*!
    @brief Maximum exponential
  */
  Real m_normal_max;

  /*!
    @brief Initial seed for RNG
  */
  Real m_seed_rng;

  /*!
    @brief Bisection step for particle trajectories
  */
  Real m_bisect_step;

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
  int m_pvr_buffer;

  /*!
    @brief Size of refinement boundary halo
  */
  int m_halo_buffer;

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
  Real m_max_diffusion_hop;

  /*!
    @brief Maximum probability that a diffusion hop is m_max_diffusion_hop long
  */
  Real m_max_hop_probability;

  /*!
    @brief Epsilon-parameter for diffusion hopping
  */
  Real m_hop_eps;

  /*!
    @brief NGP deposition in cut cells or not
  */
  bool m_irreg_ngp_deposition;

  /*!
    @brief NGP interpolation in cut cells or not
  */
  bool m_irreg_ngp_interpolation;

  /*!
    @brief Use ngp when depositing halo particles or not
  */
  bool m_ngp_halo;

  /*!
    @brief Redistribute or not
  */
  bool m_redistribute;

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
  bool m_plotDiffusionCoefficient;

  /*!
    @brief Plot particles
  */
  bool m_plot_particles;

  /*!
    @brief Plot EB particles
  */
  bool m_plot_eb_particles;

  /*!
    @brief Plot domain particles
  */
  bool m_plot_domain_particles;

  /*!
    @brief Plot source particles
  */
  bool m_plot_source_particles;

  /*!
    @brief Plot covered particles
  */
  bool m_plot_covered_particles;

  /*!
    @brief Plot energy density
  */
  bool m_plot_energy_density;

  /*!
    @brief Plot average energy or not. 
  */
  bool m_plot_average_energy;

  /*!
    @brief Mass/charge deposition type
  */
  DepositionType::Which m_deposition; 

  /*!
    @brief Plot deposition type
  */
  DepositionType::Which m_plot_deposition; 

  /*!
    @brief Deposited particles
  */
  EBAMRCellData m_phi;

  /*!
    @brief Scratch data
  */
  EBAMRCellData m_scratch;

  /*!
    @brief Mobility function
  */
  EBAMRCellData m_mobility_func;

  /*!
    @brief Cell-centered field used for interpolating velocities
  */
  EBAMRCellData m_velo_func;

  /*!
    @brief Diffusion-centerer field used for interpolating diffusion coefficients
  */
  EBAMRCellData m_faceCenteredDiffusionCoefficient_cell;

  /*!
    @brief Scratch storage for holding the non-conservative deposition
  */
  EBAMRIVData m_depositionNC;

  /*!
    @brief Scratch storage for holding the mass difference when using hybrid deposition
  */
  EBAMRIVData m_massDiff;

  /*!
    @brief Various particle containers with identifiers. 
  */
  std::map<which_container, ParticleContainer<ItoParticle> > m_ParticleContainers;

  /*!
    @brief Parse superparticles
  */
  void parse_superparticles();
  
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
  void parse_diffusion_hop();

  /*!
    @brief Parse the conservation
  */
  void parse_redistribution();

  /*!
    @brief Parse the conservation
  */
  void parseDivergenceComputation();

  /*!
    @brief Parse checkpointing method
  */
  void parse_checkpointing();

  /*!
    @brief This computes the "conservative" deposition, multiplied by kappa
  */
  template <class T>
  void depositKappaConservative(EBAMRCellData&               a_phi,
				 const ParticleContainer<T>& a_particles,
				 const DepositionType::Which  a_deposition);

  /*!
    @brief Just as depositKappaConservative, but using a PVR buffer where particle clouds on coarse levels are
    interpolated to fine levels. 
  */
  template <class T>
  void depositKappaConservativeWithPVR(EBAMRCellData&               a_phi,
					const ParticleContainer<T>& a_particles,
					const DepositionType::Which  a_deposition);

  /*!
    @brief Just as depositKappaConservative, but using halo particles when depositing particles from coarse levels to fine
    levels
  */
  template <class T>
  void depositKappaConservativeNoPVR_native(EBAMRCellData&               a_phi,
					     const ParticleContainer<T>& a_particles,
					     const DepositionType::Which  a_deposition);

  /*!
    @brief Just as depositKappaConservative, but using halo particles when depositing particles from coarse levels to fine
    levels
  */
  template <class T>
  void depositKappaConservativeNoPVR_NGP(EBAMRCellData&               a_phi,
					  const ParticleContainer<T>& a_particles,
					  const DepositionType::Which  a_deposition);

  /*!
    @brief Make the "non-conservative" kappa deposition
  */
  void depositNonConservative(EBAMRIVData& a_depositionNC, const EBAMRCellData& a_depositionKappaC);

  /*!
    @brief Make the hybrid deposition. Also compute the mass difference
  */
  void depositHybrid(EBAMRCellData& a_depositionH, EBAMRIVData& a_massDifference, const EBAMRIVData& a_depositionNC);

  /*!
    @brief Increment redistribution objects with mass to be redistributed
  */
  void incrementRedist(const EBAMRIVData& a_massDifference);

  /*!
    @brief Do redistribution on all levels, no coarse-fine stuff
  */
  void levelRedist(EBAMRCellData& a_phi);

  /*!
    @brief Increment the coarse-fine redist with mass difference (only for EB x CF)
  */
  void coarseFineIncrement(const EBAMRIVData& m_massDifference);

  /*!
    @brief Do coarse-fine redistribution
  */
  void coarseFineRedistribution(EBAMRCellData& a_phi);

  /*!
    @brief Superparticle merging with BVH trees
  */
  void bvh_merge(List<ItoParticle>& a_particles, const int a_particlesPerCell);

  /*!
    @brief Restart particles
  */
  void restart_particles(LevelData<EBCellFAB>& a_num_particles, const int a_level);

  /*!
    @brief Remove particles that are inside the EB using the implicit function
  */
  void remove_covered_particles_if(ParticleContainer<ItoParticle>& a_particles, const Real a_tol);

  /*!
    @brief Remove particles that are inside the EB using the discrete information
  */
  void remove_covered_particles_discrete(ParticleContainer<ItoParticle>& a_particles);
  
  /*!
    @brief Remove particles only in covered cells. 
  */
  void remove_covered_particles_voxels(ParticleContainer<ItoParticle>& a_particles);

  /*!
    @brief Transfer particles that are inside the EB
  */
  void transfer_covered_particles_if(ParticleContainer<ItoParticle>& a_src, ParticleContainer<ItoParticle>& a_dst, const Real a_tol);

  /*!
    @brief Do boundary intersection tests with an implicit function
  */
  void intersect_particles_if(ParticleContainer<ItoParticle>& a_particles,
			      ParticleContainer<ItoParticle>& a_eb_particles,
			      ParticleContainer<ItoParticle>& a_domain_particles,
			      const bool                        a_delete);

  /*!
    @brief Interpolate mobilities
  */
  virtual void interpolate_mobilities_mu(const int a_level, const DataIndex& a_dit);

  /*!
    @brief Interpolate mobilities
  */
  virtual void interpolate_mobilities_vel(const int a_level, const DataIndex& a_dit);

  /*!
    @brief Sign function
  */
  inline
  Real sign(const Real& a) const;

  /*!
    @brief Draw a valid random position somewhere in a cell
  */
  inline
  RealVect random_position(const RealVect a_pos,
			   const RealVect a_lo,
			   const RealVect a_hi,
			   const RealVect a_bndryCentroid,
			   const RealVect a_normal,
			   const Real     a_dx,
			   const Real     a_kappa);

  /*!
    @brief Draw a random position somewhere in a cell
  */
  RealVect random_position(const RealVect a_lo,
			   const RealVect a_hi,
			   const RealVect a_bndryCentroid,
			   const RealVect a_normal);

  /*!
    @brief Return a random position in the cube (a_lo, a_hi);
  */
  RealVect random_position(const RealVect a_lo, const RealVect a_hi);
};
#include "CD_NamespaceFooter.H"

#include "ito_solverI.H"

#endif
