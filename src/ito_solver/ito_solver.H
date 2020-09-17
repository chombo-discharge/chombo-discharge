/*!
  @file   ito_solver.H
  @brief  Declaration of an abstract class for Ito diffusion
  @author Robert Marskar
  @date   April 2020
*/

#ifndef _ITO_SOLVER_
#define _ITO_SOLVER_

#include "amr_mesh.H"
#include "computational_geometry.H"
#include "ito_species.H"
#include "ito_particle.H"
#include "EBParticleInterp.H"
#include "particle_container.H"
#include "bvh.H"
#include "point_mass.H"

#include <Particle.H>
#include <ParticleData.H>
#include <ParticleValidRegion.H>
#include <BinItem.H>

#include <random>

/*!
  @brief Base class for Ito diffusion particle models
  @details This is a particle class for solving Ito diffusion problems dX = vdt + Z*D*t over an AMR hierarchy. The diffusion
  coefficient is related to the hydrodynamic diffusion coefficient by D(fluid = sqrt(2*D(Ito))
*/
class ito_solver {
public:

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
  virtual std::string get_name();

  /*!
    @brief Get realm
  */
  virtual const std::string get_realm() const;

  /*!
    @brief Set realm
  */
  virtual void set_realm(const std::string a_realm);

  /*!
    @brief Parse options
  */
  virtual void parse_options();

  /*!
    @brief Get output plot names
  */
  virtual Vector<std::string> get_plotvar_names() const;

  /*!
    @brief Set the initial data
  */
  virtual void initial_data();

  /*!
    @brief Regrid solver
  */
  virtual void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level);

  /*!
    @brief Allocate internal storage
  */
  virtual void allocate_internals();

  /*!
    @brief Empty particle container
  */
  virtual void clear(particle_container<ito_particle>& a_particles);

  /*!
    @brief Empty particle container
  */
  virtual void clear(AMRParticles<ito_particle>& a_particles);

  /*!
    @brief Set mass to conductivity
    @details This sets mass = mass*mobility
  */
  virtual void set_mass_to_conductivity(particle_container<ito_particle>& a_particles);

  /*!
    @brief Unset mass to conductivity
    @details This sets mass = mass/mobility
  */
  virtual void unset_mass_to_conductivity(particle_container<ito_particle>& a_particles);

  /*!
    @brief Deposit conductivities
    @details This deposits weight*mobility (not multiplied by charge)
  */
  virtual void deposit_conductivity();

  /*!
    @brief Deposit conductivities
    @details This deposits weight*mobility (not multiplied by charge)
  */
  virtual void deposit_conductivity(EBAMRCellData& a_state, particle_container<ito_particle>& a_particles);

    /*!
    @brief Deposit conductivities
    @details This deposits weight*mobility (not multiplied by charge)
  */
  virtual void deposit_conductivity(EBAMRCellData& a_state, particle_container<ito_particle>& a_particles, const DepositionType::Which a_deposition);

  /*!
    @brief Deposit particles
  */
  virtual void deposit_particles();

  /*!
    @brief Deposit the particles onto a_state. The result is put on the a_comp compoennt
  */
  virtual void deposit_particles(EBAMRCellData& a_state, const AMRParticles<ito_particle>& a_particles);

  /*!
    @brief Deposit the particles onto a_state. The result is put on the a_comp compoennt
  */
  virtual void deposit_particles(EBAMRCellData&                    a_state,
				 const AMRParticles<ito_particle>& a_particles,
				 const DepositionType::Which       a_deposition);

  /*!
    @brief Deposit the total weight of the particles in the grid cells
  */
  virtual void deposit_weights(EBAMRCellData& a_weights, const AMRParticles<ito_particle>& a_particles);

  /*!
    @brief Remove particles that are inside the EB
    @details Calls the other version with m_particles
  */
  virtual void remove_eb_particles();

  /*!
    @brief Remove particles that are inside the EB
  */
  virtual void remove_eb_particles(particle_container<ito_particle>& a_particles);

  /*!
    @brief Do boundary intersection tests
  */
  virtual void intersect_particles();

  /*!
    @brief Compute particle load on specific levle
  */
  virtual void compute_loads(Vector<long int>& a_loads, const DisjointBoxLayout& a_dbl, const int a_level);

  /*!
    @brief Get number of plot variables
  */
  virtual int get_num_plotvars() const;

  /*!
    @brief Get number of particles in m_particles
  */
  virtual size_t get_num_particles(const bool a_local) const;

  /*!
    @brief Get number of particles in m_particles
  */
  virtual size_t get_num_eb_particles(const bool a_local) const;

  /*!
    @brief Get number of particles in m_particles
  */
  virtual size_t get_num_domain_particles(const bool a_local) const;

  /*!
    @brief Get number of particles in m_particles
  */
  virtual size_t get_num_source_particles(const bool a_local) const;

  /*!
    @brief Register operators
  */
  virtual void register_operators();

  /*!
    @brief Set computational geometry
  */
  virtual void set_computational_geometry(const RefCountedPtr<computational_geometry> a_compgeom);

  /*!
    @brief Set the amr object
  */
  virtual void set_amr(const RefCountedPtr<amr_mesh>& a_amr);

  /*!
    @brief Set phase
  */
  virtual void set_phase(phase::which_phase a_phase);

  /*!
    @brief Set verbosity.
  */
  virtual void set_verbosity(const int a_verbosity);

  /*!
    @brief Get current time
  */
  virtual void set_time(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Write plot data
  */
  virtual void write_plot_data(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Write data to output. Convenience function. 
  */
  virtual void write_data(EBAMRCellData& a_output, int& a_comp, const EBAMRCellData& a_data, const bool a_interp);

  /*!
    @brief Write checkpoint data into handle
  */
  virtual void write_checkpoint_level(HDF5Handle& a_handle, const int a_level) const;

  /*!
    @brief Read checkpoint data from handle
  */
  virtual void read_checkpoint_level(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Set the species
  */
  virtual void set_species(RefCountedPtr<ito_species> a_species);

  /*!
    @brief Cache the particles
  */
  virtual void pre_regrid(const int a_base, const int a_old_finest_level);

  /*!
    @brief Add particles. Container version. 
  */
  virtual void add_particles(particle_container<ito_particle>& a_part, const bool a_destructive);

  /*!
    @brief Add particles. AMR version. 
  */
  virtual void add_particles(AMRParticles<ito_particle>& a_part, const bool a_destructive);

  /*!
    @brief Add particles. Level version. 
  */
  virtual void add_particles(ParticleData<ito_particle>& a_part, const int a_level, const bool a_destructive);

  /*!
    @brief Add particles. Level and patch version. 
  */
  virtual void add_particles(ListBox<ito_particle>& a_part, const int a_lvl, const DataIndex a_dit, const bool a_destructive);

  /*!
    @brief Check if solver is mobile
  */
  virtual bool is_mobile() const;
  
  /*!
    @brief Check if solver is diffusive
  */
  virtual bool is_diffusive() const;

  /*!
    @brief Get the particle container
  */
  virtual particle_container<ito_particle>& get_particles();

  /*!
    @brief Get the particle container
  */
  virtual particle_container<ito_particle>& get_eb_particles();

  /*!
    @brief Get the particle container
  */
  virtual particle_container<ito_particle>& get_domain_particles();

  /*!
    @brief Get the particle container
  */
  virtual particle_container<ito_particle>& get_source_particles();

  /*!
    @brief Get the particle container
  */
  virtual particle_container<ito_particle>& get_scratch_particles();

  /*!
    @brief Get the species
  */
  virtual RefCountedPtr<ito_species>& get_species();

  /*!
    @brief Get state
  */
  virtual EBAMRCellData& get_state();

  /*!
    @brief Get cell-centered velocity
  */
  virtual EBAMRCellData& get_velo_func();

  /*!
    @brief Get cell-centered velocity
  */
  virtual EBAMRCellData& get_diffco_func();

  /*!
    @brief Set diffusion coefficient
  */
  virtual void set_diffco_func(const Real a_diffco);

  /*!
    @brief Set velocity
  */
  virtual void set_velocity_func(const RealVect a_vel);

  /*!
    @brief Interpolate velocities
  */
  virtual void interpolate_velocities();

  /*!
    @brief Interpolate velocities
  */
  virtual void interpolate_velocities(const int a_level, const DataIndex& a_dit);

  /*!
    @brief Interpolate diffusion
  */
  virtual void interpolate_diffusion();

  /*!
    @brief Interpolate velocities
  */
  virtual void interpolate_diffusion(const int a_level, const DataIndex& a_dit);

  /*!
    @brief Sort the m_particles container by cell
  */
  virtual void sort_particles_by_cell();

  /*!
    @brief Sort the m_particles container by cell
  */
  virtual void sort_particles_by_patch();
  
  /*!
    @brief Make superparticles
  */
  virtual void make_superparticles(const int a_particlesPerPatch);

  /*!
    @brief Make superparticles on specified level
  */
  virtual void make_superparticles(const int a_particlesPerPatch, const int a_level);

  /*!
    @brief Make superparticles on specified level and patch
  */
  virtual void make_superparticles(const int a_particlesPerPatch, const int a_level, const DataIndex a_dit);

  /*!
    @brief Remap all particles
  */
  void remap();

  /*!
    @brief Random Gaussian
  */
  RealVect random_gaussian();

  /*!
    @brief Random direction
  */
  RealVect random_direction();

  /*!
    @brief Compute min dt
  */
  Real compute_dt() const;

    /*!
    @brief Compute drift dt. This returns the minimum of (dx/dt)*a_maxCellsToMove over all levels. 
  */
  Real compute_dt(const int a_lvl) const;

  /*!
    @brief Compute drift dt. This returns the minimum of (dx/dt)*a_maxCellsToMove over all levels. 
  */
  Real compute_dt(const int a_lvl, const DataIndex a_dit, const Real a_dx) const;

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
  Vector<Real> compute_diffusion_dt(const Real a_maxCellsToMove) const;

  /*!
    @brief Compute the diffusive dt
  */
  Real compute_diffusive_dt() const;
  
  /*!
    @brief Compute diffusion dt as dx/(sqrt(2D)*eps) on all levels
  */
  Vector<Real> compute_diffusion_dt() const;

  /*!
    @brief Compute diffusion dt as dx/(sqrt(2D)*eps) on one level
  */
  Real compute_diffusion_dt(const int a_lvl) const;

  /*!
    @brief Compute diffusion dt as dx/(sqrt(2D)*eps) on one level and one patch
  */
  Real compute_diffusion_dt(const int a_lvl, const DataIndex& a_dit, const RealVect a_dx) const;

  /*!
    @brief Return phase
  */
  phase::which_phase get_phase() const;
  
protected:

  /*! 
    @brief How to checkpoint files
  */
  enum class which_checkpoint {
    particles,
    numbers
  };

  /*!
    @brief How to checkpoint files. particles => write particles to HDF5. numbers => write numbers to HDF5
  */
  which_checkpoint m_checkpointing;

  /*!
    @brief Get number of particles used for restarting
  */
  int m_ppc_restart;

  /*!
    @brief Realm
  */
  std::string m_realm;

  /*!
    @brief Boundary volume hierarchy merging tree
  */
  bvh_tree<point_mass> m_tree;

  /*!
    @brief Computational geometry. 
  */
  RefCountedPtr<computational_geometry> m_compgeom;

  /*!
    @brief AMR; needed for grid stuff
  */
  RefCountedPtr<amr_mesh> m_amr;

  /*!
    @brief Species
  */
  RefCountedPtr<ito_species> m_species;

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
  std::string m_class_name;

  /*!
    @brief RNG engine
  */
  std::mt19937_64 m_rng; 

  /*!
    @brief Uniform real distribution on [0,1]
  */
  std::uniform_real_distribution<Real> m_udist01; 

  /*!
    @brief Uniform real distribution on [-1,1]
  */
  std::uniform_real_distribution<Real> m_udist11;

  /*!
    @brief Gaussian distribution centered on 0, std=1
  */
  std::normal_distribution<Real> m_gauss01;

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
  int m_step;

  /*!
    @brief Buffer for PVR
  */
  int m_pvr_buffer;

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
    @brief Redistribute or not
  */
  bool m_redistribute;

  /*!
    @brief Flag for blending the deposition clouds with the nonconservative divergence
  */
  bool m_blend_conservation;

  /*!
    @brief Solver is diffusive
  */
  bool m_diffusive;

  /*!
    @brief Also solve for advection/convection
  */
  bool m_mobile;

  /*!
    @brief Output state
  */
  bool m_plot_phi;

  /*!
    @brief Output state
  */
  bool m_plot_vel;

  /*!
    @brief Output state
  */
  bool m_plot_dco;

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
    @brief Plot domain particles
  */
  bool m_plot_source_particles;

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
  EBAMRCellData m_state;

  /*!
    @brief Scratch data
  */
  EBAMRCellData m_scratch;

  /*!
    @brief Cell-centered field used for interpolating velocities
  */
  EBAMRCellData m_velo_func;

  /*!
    @brief Diffusion-centerer field used for interpolating diffusion coefficients
  */
  EBAMRCellData m_diffco_cell;

  /*!
    @brief Scratch storage for holding the non-conservative deposition
  */
  EBAMRIVData m_depositionNC;

  /*!
    @brief Scratch storage for holding the mass difference when using hybrid deposition
  */
  EBAMRIVData m_massDiff;

  /*!
    @brief All particles
  */
  particle_container<ito_particle> m_particles;

  /*!
    @brief This is a particle container which acts as scratch space when you want to add particles through the EB
  */
  particle_container<ito_particle> m_eb_particles;

  /*!
    @brief This is a particle container which acts as scratch space when you want to add particles through the domain side
  */
  particle_container<ito_particle> m_domain_particles;

  /*!
    @brief Source particles in case you want it. 
  */
  particle_container<ito_particle> m_source_particles;

  /*!
    @brief Data holder for copying particles that are INSIDE the EB
  */
  particle_container<ito_particle> m_scratch_particles;

  /*!
    @brief Parse superparticles
  */
  void parse_superparticles();
  
  /*!
    @brief Parse RNG options
  */
  void parse_rng();

  /*!
    @brief Parse deposition options
  */
  void parse_deposition();

  /*!
    @brief Parse bisection step
  */
  void parse_bisect_step();

  /*!
    @brief Parse PVR buffer
  */
  void parse_pvr_buffer();

  /*!
    @brief Parse plot variables
  */
  void parse_plot_vars();

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
  void parse_conservation();

  /*!
    @brief Parse checkpointing method
  */
  void parse_checkpointing();

  /*!
    @brief This computes the "conservative" deposition, multiplied by kappa
  */
  void deposit_kappaConservative(EBAMRCellData&                    a_state,
				 const AMRParticles<ito_particle>& a_particles,
				 const DepositionType::Which       a_deposition);

  /*!
    @brief Make the "non-conservative" kappa deposition
  */
  void deposit_nonConservative(EBAMRIVData& a_depositionNC, const EBAMRCellData& a_depositionKappaC);

  /*!
    @brief Make the hybrid deposition. Also compute the mass difference
  */
  void deposit_hybrid(EBAMRCellData& a_depositionH, EBAMRIVData& a_mass_diff, const EBAMRIVData& a_depositionNC);

  /*!
    @brief Increment redistribution objects with mass to be redistributed
  */
  void increment_redist(const EBAMRIVData& a_mass_diff);

  /*!
    @brief Do redistribution on all levels, no coarse-fine stuff
  */
  void level_redistribution(EBAMRCellData& a_state);

  /*!
    @brief Increment the coarse-fine redist with mass difference (only for EB x CF)
  */
  void coarse_fine_increment(const EBAMRIVData& m_mass_diff);

  /*!
    @brief Do coarse-fine redistribution
  */
  void coarse_fine_redistribution(EBAMRCellData& a_state);

  /*!
    @brief Superparticle merging with BVH trees
  */
  void bvh_merge(List<ito_particle>& a_particles, const int a_particlesPerCell);
  
  /*!
    @brief Sign function
  */
  inline Real sign(const Real& a) const;



  /*!
    @brief Restart particles
  */
  void restart_particles(LevelData<EBCellFAB>& a_num_particles, const int a_level);

  /*!
    @brief Draw a valid random position somewhere in a cell
  */
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
#endif
