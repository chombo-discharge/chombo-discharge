/*!
  @file driver.H
  @brief Declaration of the main class for time/space advancement of streamer equations
  @author Robert Marskar
  @date Nov. 2017
*/

#ifndef _DRIVER_
#define _DRIVER_

#include <RefCountedPtr.H>

#include "computational_geometry.H"
#include "time_stepper.H"
#include "amr_mesh.H"
#include "cell_tagger.H"
#include "mfis.H"
#include "geo_coarsener.H"

#include <iostream>
#include <fstream>

/*!
  @brief Main class for time/space advancement of streamer equations. This class runs a timestepper and controls input/output, 
  regrid intervals etc. 
*/
class driver {
public:
  
  /*!
    @brief Full constructor. 
  */
  driver(const RefCountedPtr<computational_geometry>& a_compgeom,
	 const RefCountedPtr<time_stepper>&           a_timestepper,
	 const RefCountedPtr<amr_mesh>&               a_amr,
	 const RefCountedPtr<cell_tagger>&            a_celltagger = RefCountedPtr<cell_tagger>(NULL),
	 const RefCountedPtr<geo_coarsener>&          a_geocoarsen = RefCountedPtr<geo_coarsener>(NULL));

  /*!
    @brief Constructor
  */
  ~driver();

  /*!
    @brief Setup and run simulations
  */
  virtual void setup_and_run();

  // /*!
  //   @brief This function will read two checkpoint files, and compute the error norms on the finest level in a_chk_coar data 
  // */
  // virtual void compute_norm(std::string a_chk_coarse, std::string a_chk_fine);

  // /*!
  //   @brief Compute the error norm on the coarsest level 
  // */
  // virtual void compute_coarse_norm(const std::string a_chk1, const std::string a_chk2, const int a_species = -1);

  virtual void parse_regrid();
  virtual void parse_restart();
  virtual void parse_memrep();
  virtual void parse_output_directory();
  virtual void parse_geo_refinement();
  virtual void parse_coarsen();
  virtual void parse_verbosity();
  virtual void parse_output_file_names();
  virtual void parse_output_intervals();
  virtual void parse_num_plot_ghost();
  virtual void parse_grow_tags();
  virtual void parse_geom_only();
  virtual void parse_ebis_memory_load_balance();
  virtual void parse_write_ebis();
  virtual void parse_read_ebis();
  virtual void parse_simulation_time();
  virtual void parse_file_depth();
  virtual void parse_plot_vars();
  virtual void parse_geometry_generation();
  
  /*!
    @brief Run the solvers untils a_time or a_max_steps has been reached
  */
  virtual void run(const Real a_start_time, const Real a_end_time, const int a_max_steps);

  /*!
    @brief Setup function
  */
  virtual void setup(const int a_init_regrids, const bool a_restart, const std::string a_restart_file);

  /*!
    @brief Set amr
  */
  virtual void set_amr(const RefCountedPtr<amr_mesh>& a_amr);

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

  // /*!
  //   @brief Dump mass for cdr solvers
  // */
  // virtual void set_dump_mass(const bool a_dump_mass);

  // /*!
  //   @brief Dump charge for cdr and sigma solvers
  // */
  // virtual void set_dump_charge(const bool a_dump_charge);

  /*!
    @brief Output on centroids or cell centers
  */
  virtual void set_output_centroids(const bool a_output_centroids);

  /*!
    @brief Sanity check
  */
  virtual void sanity_check();

protected:

  /*!
    @brief Constructor
  */
  driver();

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
    @brief Geometry generation method
  */
  std::string m_geometry_generation;

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
    @brief Geometry scan level
  */
  int m_geo_scan_level;

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
  EBAMRBool m_cached_tags;

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
    @brief Allow grid coarsening
  */
  bool m_allow_coarsen;

  /*!
    @brief Write memory report
  */
  bool m_write_memory;

  /*!
    @brief Use recursive regrid (true) or regrid all levels at once (false)
  */
  bool m_recursive_regrid;

  /*!
    @brief Output on centroids or not
  */
  bool m_output_centroids;

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
    @brief Plot tagged cells
  */
  bool m_plot_tags;

  /*!
    @brief Plot MPI ranks
  */
  bool m_plot_ranks;

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
    @brief Read vector data up to and including a_step
  */
  virtual void read_vector_data(HDF5HeaderData& a_header, Vector<Real>& a_data, const std::string a_name, const int a_elements);

  /*!
    @brief Do a regrid step. This regrid method regrids data inside each solver. The exception is the Poisson solver
    whose solution is recomputed. The same goes for stationary RTE methods. 
  */
  virtual void regrid(const int a_lmin, const int a_lmax, const bool a_use_initial_data);

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
			     const Real a_solver_regrid_time);
  
  /*!
    @brief Set the computational geometry
  */
  virtual void set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom);

  /*!
    @brief Set the time stepper
  */
  virtual void set_time_stepper(const RefCountedPtr<time_stepper>& a_timestepper);

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
    @brief Write the memory usage to file
  */
  virtual void write_memory_usage();

  /*!
    @brief Write a checkpoint file
  */
  virtual void write_checkpoint_file();
  
  /*!
    @brief Write checkpoint data
  */
  virtual void write_checkpoint_level(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Write checkpoint data
  */
  virtual void read_checkpoint_level(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Write vector data to a HDF5Header. This writes up to a_step
  */
  virtual void write_vector_data(HDF5HeaderData& a_header, const Vector<Real>& a_data, const std::string a_name, const int a_elem);

  /*!
    @brief Write the geometry only
  */
  virtual void write_geometry();
  
  // /*!
  //   @brief Open file for mass dump
  // */
  // virtual void open_mass_dump_file(ofstream& a_file);

  // /*!
  //   @brief Open file for mass dump
  // */
  // virtual void open_charge_dump_file(ofstream& a_file);

  // /*!
  //   @brief Dump mass
  // */
  // virtual void dump_mass(ofstream& a_file);

  // /*!
  //   @brief Dump charge
  // */
  // virtual void dump_charge(ofstream& a_file);

  // /*!
  //   @brief Close file for mass dump
  // */
  // virtual void close_mass_dump_file(ofstream& a_file);

  // /*!
  //   @brief Close file for mass dump
  // */
  // virtual void close_charge_dump_file(ofstream& a_file);

  /*!
    @brief Tag cells for refinement. This computes cell tags and global tags (union of cell tags with geometric tags);
  */
  virtual bool tag_cells(Vector<IntVectSet>& a_all_tags, EBAMRTags& a_cell_tags);

  /*!
    @brief Get number of plot variables
  */
  virtual int get_num_plot_vars() const;

  /*!
    @brief Get names for the plot variables
  */
  virtual Vector<std::string> get_plotvar_names() const;

  /*!
    @brief Write internal plot data
  */
  virtual void write_plot_data(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Write cell tags to file
  */
  virtual void write_tags(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Write MPI ranks to output
  */
  virtual void write_ranks(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Return the finest tag level
  */
  virtual int get_finest_tag_level(const EBAMRTags& a_cell_tags) const;
};
#endif
