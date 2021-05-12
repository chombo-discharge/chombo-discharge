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
#include "CD_AmrMesh.H"
#include "cell_tagger.H"
#include "mfis.H"
#include "geo_coarsener.H"

#include <iostream>
#include <fstream>

#include "CD_NamespaceHeader.H"

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
	 const RefCountedPtr<AmrMesh>&               a_amr,
	 const RefCountedPtr<cell_tagger>&            a_celltagger = RefCountedPtr<cell_tagger>(NULL),
	 const RefCountedPtr<geo_coarsener>&          a_geocoarsen = RefCountedPtr<geo_coarsener>(NULL));

  /*!
    @brief Constructor
  */
  ~driver();

  /*!
    @brief Setup and run simulations
  */
  void setup_and_run(const std::string a_input_file);
  
  /*!
    @brief Run the solvers untils a_time or a_max_steps has been reached
  */
  void run(const Real a_start_time, const Real a_end_time, const int a_max_steps);

  /*!
    @brief Setup function
  */
  void setup(const std::string a_input_file, const int a_init_regrids, const bool a_restart, const std::string a_restart_file);

  /*!
    @brief Parse options
  */
  void parseOptions();

  /*!
    @brief Parse runtime options
  */
  void parseRuntimeOptions();

protected:

  /*!
    @brief Memory report mode enum
  */
  enum class memory_report_mode {
    overall,
    unfreed,
    allocated
  };

  /*!
    @brief Constructor
  */
  driver() = delete;

  /*!
    @brief Index space
  */
  RefCountedPtr<mfis> m_multifluidIndexSpace;

  /*!
    @brief AMR
  */
  RefCountedPtr<AmrMesh> m_amr;

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
  time_code m_timecode;

  /*!
    @brief Memory report mode
  */
  memory_report_mode m_memory_mode;

  /*!
    @brief Input file
  */
  std::string m_input_file;

  /*!
    @brief
  */
  std::string m_Realm;

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
    @brief Write regrid files or not
  */
  bool m_write_regrid_files;

  /*!
    @brief Write restart files or not
  */
  bool m_write_restart_files;

  /*!
    @brief Allow grid coarsening
  */
  bool m_allow_coarsen;

  /*!
    @brief Write memory report
  */
  bool m_write_memory;

  /*!
    @brief Write loads to file or not
  */
  bool m_write_loads;

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
  bool m_ebis_memory_LoadBalancing;

  /*!
    @brief Plot tagged cells
  */
  bool m_plot_tags;

  /*!
    @brief Plot MPI ranks
  */
  bool m_plot_ranks;

  /*!
    @brief Plot level sets
  */
  bool m_plot_levelset;

  /*!
    @brief Set amr
  */
  void set_amr(const RefCountedPtr<AmrMesh>& a_amr);

  /*!
    @brief Sanity check
  */
  void sanityCheck();

  /*!
    @brief Parse geometrical refinement
  */
  void parse_geo_refinement();

  /*!
    @brief Parse plot variables
  */
  void parse_plot_vars();

  /*!
    @brief Parse option for geometry generation
  */
  void parse_geometry_generation();

  /*!
    @brief Create output directories
  */
  void create_output_directories();

  /*!
    @brief Allocate internals
  */
  void allocate_internals();
  
  /*!
    @brief Cache tags
  */
  void cache_tags(const EBAMRTags& a_tags);

  /*!
    @brief Deallocate internal storage
  */
  void deallocate_internals();

  /*!
    @brief Write the geometry to file
  */
  void write_ebis();

  /*!
    @brief Get geometric tags
  */
  void get_geom_tags();

  /*!
    @brief Get boxes and points for everything
  */
  void get_loads_and_boxes(long long& a_myPoints,
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
  void grid_report();

  /*!
    @brief Number formatting game. This just prints the number 123456 as 123,456
  */
  const std::string number_fmt(const long long a_number, char a_sep = ',') const;

  /*!
    @brief Number formatting game. This just prints the number 123456 as 123,456
  */
  const Vector<std::string> number_fmt(const Vector<long long> a_number, char a_sep = ',') const;

  /*!
    @brief Memory report
  */
  void memory_report(const memory_report_mode a_mode);

  /*!
    @brief Read checkpoint file
  */
  void read_checkpoint_file(const std::string& a_restart_file);

  /*!
    @brief Read vector data up to and including a_step
  */
  void read_vector_data(HDF5HeaderData& a_header, Vector<Real>& a_data, const std::string a_name, const int a_elements);

  /*!
    @brief Do a regrid step. This regrid method regrids data inside each solver. The exception is the Poisson solver
    whose solution is recomputed. The same goes for stationary RTE methods. 
  */
  void regrid(const int a_lmin, const int a_lmax, const bool a_use_initial_data);

  /*!
    @brief Regrid internal storage for this class
  */
  void regrid_internals(const int a_old_finest_level, const int a_new_finest_level);

  /*!
    @brief Write a regrid report
  */
  void regrid_report(const Real a_total_time,
		     const Real a_tag_time,
		     const Real a_base_regrid_time,
		     const Real a_solver_regrid_time);
  
  /*!
    @brief Set the computational geometry
  */
  void set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom);

  /*!
    @brief Set the time stepper
  */
  void set_time_stepper(const RefCountedPtr<time_stepper>& a_timestepper);

  /*!
    @brief Set the cell tagger
  */
  void set_cell_tagger(const RefCountedPtr<cell_tagger>& a_celltagger);

  /*!
    @brief Set the geo coarsener
  */
  void set_geo_coarsen(const RefCountedPtr<geo_coarsener>& a_geocoarsen);

  /*!
    @brief Setup function
  */
  void setup_fresh(const int a_init_regrids);

  /*!
    @brief Set for restart
  */
  void setup_for_restart(const int a_init_regrids, const std::string a_restart_file);

  /*!
    @brief Set up for geometry only
  */
  void setup_geometry_only();

  /*!
    @brief Check if restart file exists
  */
  void check_restart_file(const std::string a_restart_file) const;

  /*!
    @brief Write a time step report
  */
  void step_report(const Real a_start_time, const Real a_end_time, const int a_max_steps);

  /*!
    @brief Write a plot file
  */
  void write_plot_file();

  /*!
    @brief Write a regrid file
  */
  void write_regrid_file();

  /*!
    @brief Write a restart file
  */
  void write_restart_file();

  /*!
    @brief Write a crash file
  */
  void write_crash_file();

  /*!
    @brief Write a plot file
  */
  void write_plot_file(const std::string a_filename);

  /*!
    @brief Write the memory usage to file
  */
  void write_memory_usage();

  /*!
    @brief Write computational loads
  */
  void write_computational_loads();

  /*!
    @brief Write a checkpoint file
  */
  void write_checkpoint_file();
  
  /*!
    @brief Write checkpoint data
  */
  void write_checkpoint_level(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Write tags
  */
  void write_checkpoint_tags(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Write checkpoint Realm loads
  */
  void write_checkpoint_Realm_loads(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Write checkpoint data
  */
  void read_checkpoint_level(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Read Realm loads. 
  */
  void read_checkpoint_Realm_loads(Vector<long int>& a_loads, HDF5Handle& a_handle, const std::string a_Realm, const int a_level);

  /*!
    @brief Write vector data to a HDF5Header. This writes up to a_step
  */
  void write_vector_data(HDF5HeaderData& a_header, const Vector<Real>& a_data, const std::string a_name, const int a_elem);

  /*!
    @brief Write the geometry only
  */
  void write_geometry();

  /*!
    @brief Tag cells for refinement. This computes cell tags and global tags (union of cell tags with geometric tags);
  */
  bool tag_cells(Vector<IntVectSet>& a_all_tags, EBAMRTags& a_cell_tags);

  /*!
    @brief Get number of plot variables
  */
  int get_num_plot_vars() const;

  /*!
    @brief Get names for the plot variables
  */
  Vector<std::string> get_plotvar_names() const;

  /*!
    @brief Write internal plot data
  */
  void write_plot_data(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Write cell tags to file
  */
  void write_tags(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Write MPI ranks to output
  */
  void write_ranks(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Write level sets
  */
  void write_levelset(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Return the finest tag level
  */
  int get_finest_tag_level(const EBAMRTags& a_cell_tags) const;

  /*!
    @brief Rebuild Parmparse
  */
  void rebuildParmParse() const;
};
#include "CD_NamespaceFooter.H"
#endif
