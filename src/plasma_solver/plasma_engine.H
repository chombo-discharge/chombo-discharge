/*!
  @file plasma_engine.H
  @brief Declaration of the main class for time/space advancement of streamer equations
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _PLASMA_ENGINE_
#define _PLASMA_ENGINE_

#include <RefCountedPtr.H>

#include "computational_geometry.H"
#include "plasma_kinetics.H"
#include "time_stepper.H"
#include "physical_domain.H"
#include "amr_mesh.H"
#include "cell_tagger.H"
#include "mfis.H"
#include "tags.H"
#include "geo_coarsener.H"

#include <iostream>
#include <fstream>

/*!
  @brief Main class for time/space advancement of streamer equations. This class runs a timestepper and controls input/output, 
  regrid intervals etc. 

  Class options
  -------------

      plasma_engine.verbosity                       = 1             # Engine verbosity
      plasma_engine.plot_interval                   = 10            # Plot interval
      plasma_engine.regrid_interval                 = 10            # Regrid interval
      plasma_engine.checkpoint_interval             = 10            # Checkpoint interval
      plasma_engine.init_regrids                    = 0             # Number of initial regrids
      plasma_engine.start_time                      = 0             # Start time (fresh simulations only)
      plasma_engine.stop_time                       = 1.0           # Stop time
      plasma_engine.max_steps                       = 100           # Maximum number of steps
      plasma_engine.geometry_only                   = false         # Special option that ONLY plots the geometry
      plasma_engine.ebis_memory_load_balance        = false         # Use memory as loads for EBIS generation
      plasma_engine.num_plot_ghost                  = 1             # Number of ghost cells to include in plots
      plasma_engine.output_mode                     = full          # Output mode
      plasma_engine.output_directory                = ./            # Output directory
      plasma_engine.output_names                    = simulation    # Simulation output names
      plasma_engine.restart                         = false         # Do restart or not
      plasma_engine.restart_file                    = default_file  # Restart file
      plasma_engine.restart_mode                    = full          # Restart mode. 
      plasma_engine.refine_geometry                 = -1            # Refine geometry, -1 => Refine all the way down
      plasma_engine.refine_electrodes               = -1            # Refine electrode surfaces. -1 => equal to refine_geometry
      plasma_engine.refine_dielectrics              = -1            # Refine dielectric surfaces. -1 => equal to refine_geometry
      plasma_engine.refine_electrode_gas_interface  = -1            # Refine electrode-gas interfaces. -1 => ----"-----
      plasma_engine.refine_dielectric_gas_interface = -1            # Refine dielectric-gas interfaces. -1 => ----"-----
      plasma_engine.refine_solid_gas_interface      = -1            # Refine solid-gas interfaces. -1 => ----"-----
      plasma_engine.refine_solid_solid_interface    = -1            # Refine solid-solid interfaces. -1 => ----"-----

*/
class plasma_engine {
public:


  
  /*!
    @brief Full constructor. 
  */
  plasma_engine(const RefCountedPtr<physical_domain>&        a_physdom,
		const RefCountedPtr<computational_geometry>& a_compgeom,
		const RefCountedPtr<plasma_kinetics>&        a_plaskin,
		const RefCountedPtr<time_stepper>&           a_timestepper,
		const RefCountedPtr<amr_mesh>&               a_amr,
		const RefCountedPtr<cell_tagger>&            a_celltagger = RefCountedPtr<cell_tagger>(NULL),
		const RefCountedPtr<geo_coarsener>&          a_geocoarsen = RefCountedPtr<geo_coarsener>(NULL));

  /*!
    @brief Constructor
  */
  ~plasma_engine();

  /*!
    @brief Setup and run simulations
  */
  virtual void setup_and_run();

  /*!
    @brief Only set up and solve the Poisson equation
  */
  virtual void setup_poisson_only();

  /*!
    @brief This function will read two checkpoint files, and compute the error norms on the finest level in a_chk_coar data 
  */
  virtual void compute_norm(std::string a_chk_coarse, std::string a_chk_fine);

  /*!
    @brief Run the solvers untils a_time or a_max_steps has been reached
  */
  virtual void run(const Real a_start_time, const Real a_end_time, const int a_max_steps);

  /*!
    @brief Setup function
  */
  virtual void setup(const int a_init_regrids, const bool a_restart, const std::string a_restart_file);

  /*!
    @brief Verbosity
  */
  virtual void set_verbosity(const int a_verbosity);

  /*!
    @brief Set amr
  */
  virtual void set_amr(const RefCountedPtr<amr_mesh>& a_amr);


  /*!
    @brief Set geometric refinement depth
  */
  virtual void set_geom_refinement_depth(const int a_depth = -1);

  /*!
    @brief Set geometrical refinement depths. This provides a little more flexibility than the basic call. 
    @param[in] a_depth1 Conductor tag depth
    @param[in] a_depth2 Dielectric tag depth
    @param[in] a_depth3 Conductor-gas interface tag depth
    @param[in] a_depth4 Dielectric-gas interface tag depth
    @param[in] a_depth5 Solid-gas interface tag depth
    @param[in] a_depth6 Solid-solid interface tag depth
  */
  virtual void set_geom_refinement_depth(const int a_depth1,
					 const int a_depth2,
					 const int a_depth3,
					 const int a_depth4,
					 const int a_depth5,
					 const int a_depth6);

  /*!
    @brief Set output mode
  */
  virtual void set_output_mode(const output_mode::which_mode a_mode);

  /*!
    @brief Set restart mode
  */
  virtual void set_restart_mode(const restart_mode::which_mode a_mode);

  /*!
    @brief Set memory report mode
  */
  virtual void set_memory_report_mode(const memory_report_mode::which_mode a_memory_mode);

  /*!
    @brief Set output directory
  */
  virtual void set_output_directory(const std::string a_output_dir);

  /*!
    @brief Set output file names
  */
  virtual void set_output_file_names(const std::string a_output_names);
  
  /*!
    @brief Set potential function
  */
  virtual void set_potential(Real (*a_potential)(const Real a_time));

  /*!
    @brief Set wall potential distribution on domain faces
  */
  virtual void set_poisson_wall_func(const int a_dir, const Side::LoHiSide a_side, Real (*a_func)(const RealVect a_pos));

  /*!
    @brief Set wall function
  */
  virtual void set_poisson_wall_func(Real (*a_func)(const RealVect a_pos));
  
  /*!
    @brief Set the regrid interval
  */
  virtual void set_regrid_interval(const int a_regrid_interval);

  /*!
    @brief Set the plot interval
  */
  virtual void set_plot_interval(const int a_plot_interval);

  /*!
    @brief Set the checkpoint interval
  */
  virtual void set_checkpoint_interval(const int a_chk_interval);

  /*!
    @brief Set number of ghost cells to plot
  */
  virtual void set_num_plot_ghost(const int a_num_plot_ghost);

  /*!
    @brief Set restart
  */
  virtual void set_restart(const bool a_restart);

  /*!
    @brief Set load balance for memory
  */
  virtual void set_ebis_memory_load_balance(const bool a_balance);

  /*!
    @brief Write EBIS to HDF5
  */
  virtual void set_write_ebis(const bool a_write_ebis);

    /*!
    @brief Read EBIS to HDF5
  */
  virtual void set_read_ebis(const bool a_read_ebis);
  
  /*!
    @brief Only geometry
  */
  virtual void set_geom_only(const bool a_geom_only);

  /*!
    @brief Restart step file
  */
  virtual void set_restart_step(const int a_restart_step);

  /*!
    @brief Set simulation start time
  */
  virtual void set_start_time(const Real a_start_time);

  /*!
    @brief Set simulation stop time
  */
  virtual void set_stop_time(const Real a_stop_time);

  /*!
    @brief Set maximum number of steps
  */
  virtual void set_max_steps(const int a_max_steps);
  
  /*!
    @brief Set the maximum plot depth
  */
  virtual void set_max_plot_depth(const int a_max_plot_depth);

  /*!
    @brief Set the maximum plot depth
  */
  virtual void set_max_chk_depth(const int a_max_chk_depth);

  /*!
    @brief Set number of initial regrids
  */
  virtual void set_init_regrids(const int a_init_regrids);

  /*!
    @brief Set tag growth
  */
  virtual void set_grow_tags(const int a_grow_tags);

  /*!
    @brief Dump mass for cdr solvers
  */
  virtual void set_dump_mass(const bool a_dump_mass);

    /*!
    @brief Dump charge for cdr and sigma solvers
  */
  virtual void set_dump_charge(const bool a_dump_charge);

  /*!
    @brief Sanity check
  */
  virtual void sanity_check();

protected:

  /*!
    @brief Constructor
  */
  plasma_engine();

  /*!
    @brief Default function that is passed to mapped boundary conditions
  */
  static Real s_constant_one(const RealVect a_pos);

  /*!
    @brief Index space
  */
  RefCountedPtr<mfis> m_mfis;

  /*!
    @brief AMR
  */
  RefCountedPtr<amr_mesh> m_amr;

  /*!
    @brief Geometry
  */
  RefCountedPtr<computational_geometry> m_compgeom;

  /*!
    @brief Plasma kinetics
  */
  RefCountedPtr<plasma_kinetics> m_plaskin;

  /*!
    @brief Physical domain
  */
  RefCountedPtr<physical_domain> m_physdom;

  /*!
    @brief Time stepper
  */
  RefCountedPtr<time_stepper> m_timestepper;

  /*!
    @brief Cell tagger
  */
  RefCountedPtr<cell_tagger> m_celltagger;

    /*!
    @brief Geometric coarsener
  */
  RefCountedPtr<geo_coarsener> m_geocoarsen;

  /*!
    @brief Time code
  */
  time_code::which_code m_timecode;

  /*!
    @brief Output mode
  */
  output_mode::which_mode m_output_mode;

  /*!
    @brief Restart mode
  */
  restart_mode::which_mode m_restart_mode;

  /*!
    @brief Memory report mode
  */
  memory_report_mode::which_mode m_memory_mode;

  /*!
    @brief Output directory
  */
  std::string m_output_dir;

  /*!
    @brief Output filenames
  */
  std::string m_output_names;

  /*!
    @brief Name of the file holding the finest EBIS info for the gas phase
  */
  std::string m_ebis_gas_file;

    /*!
    @brief Name of the file holding the finest EBIS info for the gas phase
  */
  std::string m_ebis_sol_file;

  /*!
    @brief Verbosity
  */
  int m_verbosity;

  /*!
    @brief Time step
  */
  int m_step;

  /*!
    @brief Regrid interval
  */
  int m_regrid_interval;

  /*!
    @brief Checkpoint interval
  */
  int m_chk_interval;

  /*!
    @brief Plot interval
  */
  int m_plot_interval;

  /*!
    @brief Geometric tag depth
  */
  int m_geom_tag_depth;

  /*!
    @brief Geometric tag depth for dielectrics
  */
  int m_dielectric_tag_depth;

  /*!
    @brief Geometric tag depth for conductors
  */
  int m_conductor_tag_depth;

  /*!
    @brief Geometric tag depth for gas-solid interfaces
  */
  int m_gas_solid_interface_tag_depth;

  /*!
    @brief Geometry tag depth for solid-solid interfaces
  */
  int m_solid_solid_interface_tag_depth;

  /*!
    @brief Geometric tag depth for gas-dielectric interfaces
  */
  int m_gas_dielectric_interface_tag_depth;

  /*!
    @brief Geometric tag depth for gas-conductor interfaces
  */
  int m_gas_conductor_interface_tag_depth;

  /*!
    @brief Maximum number of steps
  */
  int m_max_steps;

  /*!
    @brief Maximum plot depth
  */
  int m_max_plot_depth;

  /*!
    @brief Maximum plot depth
  */
  int m_max_chk_depth;

  /*!
    @brief Number of initial regrids
  */
  int m_init_regrids;

  /*!
    @brief Number of ghost cells to plot
  */
  int m_num_plot_ghost;

  /*!
    @brief Restart step
  */
  int m_restart_step;

  /*!
    @brief Grow tagged cells
  */
  int m_grow_tags;

  /*!
    @brief Tags
  */
  Vector<IntVectSet> m_geom_tags;

  /*!
    @brief Tags
  */
  EBAMRTags m_tags;

  /*!
    @brief Cached tags. Needed when we regrid and change box ownership.
  */
  EBAMRTags m_cached_tags;

  /*!
    @brief Potential
  */
  Real (*m_potential)(const Real a_time);

  /*!
    @brief Function pointer for Poisson BC on domain walls
  */
  Real (*m_wall_func_x_lo)(const RealVect a_pos);

    /*!
    @brief Function pointer for Poisson BC on domain walls
  */
  Real (*m_wall_func_x_hi)(const RealVect a_pos);

    /*!
    @brief Function pointer for Poisson BC on domain walls
  */
  Real (*m_wall_func_y_lo)(const RealVect a_pos);

  /*!
    @brief Function pointer for Poisson BC on domain walls
  */
  Real (*m_wall_func_y_hi)(const RealVect a_pos);

#if CH_SPACEDIM==3
  /*!
    @brief Function pointer for Poisson BC on domain walls
  */
  Real (*m_wall_func_z_lo)(const RealVect a_pos);

  /*!
    @brief Function pointer for Poisson BC on domain walls
  */
  Real (*m_wall_func_z_hi)(const RealVect a_pos);
#endif

  /*!
    @brief Simulation start time
  */
  Real m_start_time;

  /*!
    @brief Simulation stop time
  */
  Real m_stop_time;

  /*!
    @brief Time increment
  */
  Real m_dt;

  /*!
    @brief Time
  */
  Real m_time;

  /*!
    @brief 
  */
  Real m_capacitance;

  /*!
    @brief MPI timer
  */
  Real m_wallclock_start;
  
  /*!
    @brief MPI timer
  */
  Real m_wallclock1;

  /*!
    @brief MPI timer
  */
  Real m_wallclock2;

  /*!
    @brief Potential has been set
  */
  bool m_potential_set;

  /*!
    @brief Restart or not
  */
  bool m_restart;

  /*!
    @brief Geometry only
  */
  bool m_geometry_only;

  /*!
    @brief Use memory as load balance for EBIS
  */
  bool m_ebis_memory_load_balance;

  /*!
    @brief Write EBIS to HDF5. Might be useful for restarting simulations
  */
  bool m_write_ebis;

  /*!
    @brief Read EBIS file from HDF5. Useful for restarting simulation.
  */
  bool m_read_ebis;

  /*!
    @brief Dump total mass (or charge) for every cdr solver at every time_step
  */
  bool m_dump_mass;

  /*!
    @brief Dump charge for cdr and surface charge solvers at every time step
  */
  bool m_dump_charge;

  /*!
    @brief Use the old or new IO functions
  */
  bool m_new_io;

  // Output parameters that go into checkpoint files
  Vector<Real> m_time_history;        // Simulation time
  Vector<Real> m_voltage_history;     // Applied voltage
  Vector<Real> m_dt_history;          // Time steps
  Vector<Real> m_exec_history;        // Execution time per time step
  Vector<Real> m_charge_history;      // Total charge in domain
  Vector<Real> m_ohmic_induced_I;     // Ohmic induced current
  Vector<Real> m_capacitive_induced;  // Capacitively induced current
  Vector<Real> m_domain_I;            // Current through domain faces
  Vector<Real> m_electrode_I;         // Current through electrode EBs


  /*!
    @brief Allocate internals
  */
  virtual void allocate_internals();

  /*!
    @brief Add potential to the output data holder
  */
  virtual void add_potential_to_output(EBAMRCellData& a_output, const int a_curVar);

  /*!
    @brief Add the electric field to the output data holder
  */
  virtual void add_electric_field_to_output(EBAMRCellData& a_output, const int a_curVar);

  /*!
    @brief Add the electric field to the output data holder
  */
  virtual void add_space_charge_to_output(EBAMRCellData& a_output, const int a_curVar);

  /*!
    @brief Add the surface charge to the output data holder
  */
  virtual void add_surface_charge_to_output(EBAMRCellData& a_output, const int a_curVar);

  /*!
    @brief Add the current density to output
  */
  virtual void add_current_density_to_output(EBAMRCellData& a_output, const int a_curVar);

  /*!
    @brief Add the ion densities to the output data holder
  */
  virtual void add_cdr_densities_to_output(EBAMRCellData& a_output, const int a_curVar);

  /*!
    @brief Add the ion velocities to the output data holder
  */
  virtual void add_cdr_velocities_to_output(EBAMRCellData& a_output, const int a_curVar);

  /*!
    @brief Add the ion source terms to the output data holder
  */
  virtual void add_cdr_source_to_output(EBAMRCellData& a_output, const int a_curVar);

  /*!
    @brief Add the photon densities to the output data holder
  */
  virtual void add_rte_densities_to_output(EBAMRCellData& a_output, const int a_curVar);

  /*!
    @brief Add the photon source terms to the output data holder
  */
  virtual void add_rte_source_to_output(EBAMRCellData& a_output, const int a_curVar);

  /*!
    @brief Add tracer fields to output data holder
  */
  virtual void add_tracer_fields_to_output(EBAMRCellData& a_output, const int a_curVar);

  /*!
    @brief Cache tags
  */
  virtual void cache_tags(const EBAMRTags& a_tags);

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocate_internals();

  /*!
    @brief Write the geometry to file
  */
  virtual void write_ebis();

  /*!
    @brief Get geometric tags
  */
  virtual void get_geom_tags();

  /*!
    @brief Get boxes and points for everything
  */
  virtual void get_loads_and_boxes(long long& a_myPoints,
				   long long& a_myPointsGhosts,
				   long long& a_myBoxes,
				   long long& a_totalPoints,
				   long long& a_totalPointsGhosts,
				   long long& a_totalBoxes,
				   Vector<long long>& a_my_level_boxes,
				   Vector<long long>& a_total_level_boxes,
				   Vector<long long>& a_my_level_points,
				   Vector<long long>& a_total_level_points,
				   const int& a_finestLevel,
				   const Vector<DisjointBoxLayout>& a_grids);

  /*!
    @brief Write a grid report
  */
  virtual void grid_report();

  /*!
    @brief Memory report
  */
  virtual void memory_report(const memory_report_mode::which_mode a_mode);

  /*!
    @brief Read checkpoint file
  */
  virtual void read_checkpoint_file(const std::string& a_restart_file);

  /*!
    @brief Read checkpoint file
  */
  virtual void old_read_checkpoint_file(const std::string& a_restart_file);
  
  /*!
    @brief Read checkpoint file
  */
  virtual void new_read_checkpoint_file(const std::string& a_restart_file);

  /*!
    @brief Read vector data up to and including a_step
  */
  virtual void read_vector_data(HDF5HeaderData& a_header, Vector<Real>& a_data, const std::string a_name, const int a_elements);

  /*!
    @brief Do a regrid step. This regrid method regrids data inside each solver. The exception is the Poisson solver
    whose solution is recomputed. The same goes for stationary RTE methods. 
  */
  virtual void regrid(const bool a_use_initial_data);

  /*!
    @brief Regrid internal storage for this class
  */
  virtual void regrid_internals(const int a_old_finest_level, const int a_new_finest_level);

  /*!
    @brief Write a regrid report
  */
  virtual void regrid_report(const Real a_total_time,
			     const Real a_tag_time,
			     const Real a_base_regrid_time,
			     const Real a_solver_regrid_time,
			     const Real a_elliptic_solve_time,
			     const Real a_solver_filling_time);
  
  /*!
    @brief Set the computational geometry
  */
  virtual void set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom);

  /*!
    @brief Set the plasma kinetics
  */
  virtual void set_plasma_kinetics(const RefCountedPtr<plasma_kinetics>& a_plaskin);

  /*!
    @brief Set the time stepper
  */
  virtual void set_time_stepper(const RefCountedPtr<time_stepper>& a_timestepper);
  
    /*!
    @brief Set physical domain
  */
  virtual void set_physical_domain(const RefCountedPtr<physical_domain>& a_physdom);

  /*!
    @brief Set the cell tagger
  */
  virtual void set_cell_tagger(const RefCountedPtr<cell_tagger>& a_celltagger);

  /*!
    @brief Set the geo coarsener
  */
  virtual void set_geo_coarsen(const RefCountedPtr<geo_coarsener>& a_geocoarsen);

  /*!
    @brief Setup function
  */
  virtual void setup_fresh(const int a_init_regrids);

  /*!
    @brief Set for restart
  */
  virtual void setup_for_restart(const int a_init_regrids, const std::string a_restart_file);

  /*!
    @brief Set up for geometry only
  */
  virtual void setup_geometry_only();

  /*!
    @brief Write a time step report
  */
  virtual void step_report(const Real a_start_time, const Real a_end_time, const int a_max_steps);

  /*!
    @brief Write a plot file
  */
  virtual void write_plot_file();

  /*!
    @brief Write a checkpoint file
  */
  virtual void write_checkpoint_file();
  
  /*!
    @brief Write a checkpoint file
  */
  virtual void old_write_checkpoint_file();

  /*!
    @brief Write a checkpoint file
  */
  virtual void new_write_checkpoint_file();

  /*!
    @brief Write vector data to a HDF5Header. This writes up to a_step
  */
  virtual void write_vector_data(HDF5HeaderData& a_header, const Vector<Real>& a_data, const std::string a_name, const int a_elem);

  /*!
    @brief Write the geometry only
  */
  virtual void write_geometry();
  
  /*!
    @brief Open file for mass dump
  */
  virtual void open_mass_dump_file(ofstream& a_file);

  /*!
    @brief Open file for mass dump
  */
  virtual void open_charge_dump_file(ofstream& a_file);

  /*!
    @brief Dump mass
  */
  virtual void dump_mass(ofstream& a_file);

  /*!
    @brief Dump charge
  */
  virtual void dump_charge(ofstream& a_file);

  /*!
    @brief Close file for mass dump
  */
  virtual void close_mass_dump_file(ofstream& a_file);

    /*!
    @brief Close file for mass dump
  */
  virtual void close_charge_dump_file(ofstream& a_file);

  /*!
    @brief Tag cells for refinement. This computes cell tags and global tags (union of cell tags with geometric tags);
  */
  virtual bool tag_cells(Vector<IntVectSet>& a_all_tags, EBAMRTags& a_cell_tags);

  /*!
    @brief Get output variable names for a given output mode
  */
  virtual Vector<string> get_output_variable_names();
};
#endif
