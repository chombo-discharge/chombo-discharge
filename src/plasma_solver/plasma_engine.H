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

/*!
  @brief Main class for time/space advancement of streamer equations
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
		const RefCountedPtr<cell_tagger>&            a_celltagger = RefCountedPtr<cell_tagger>(NULL));

  /*!
    @brief Constructor
  */
  ~plasma_engine();

  /*!
    @brief Do initial regrids. This regrid method uses the initial data inside each solver. 
  */
  virtual void initial_regrids(const int a_init_regrids);

  /*!
    @brief Run the solvers untils a_time or a_max_steps has been reached
  */
  virtual void run(const Real a_start_time, const Real a_end_time, const int a_max_steps);

  /*!
    @brief Setup function
  */
  virtual void setup_fresh(const int a_init_regrids);

  /*!
    @brief Set for restart
  */
  virtual void setup_for_restart(const std::string a_restart_file);

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
    @brief Sanity check
  */
  virtual void sanity_check();

protected:

  /*!
    @brief Constructor
  */
  plasma_engine();

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
    @brief Time stpeper
  */
  RefCountedPtr<time_stepper> m_timestepper;

  /*!
    @brief Set the cell tagger
  */
  RefCountedPtr<cell_tagger> m_celltagger;

  /*!
    @brief Time code
  */
  time_code::which_code m_timecode;

  /*!
    @brief Output mode
  */
  output_mode::which_mode m_output_mode;

  /*!
    @brief Output directory
  */
  std::string m_output_dir;

  /*!
    @brief Output filenames
  */
  std::string m_output_names;

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
    @brief Tags
  */
  Vector<IntVectSet> m_cell_tags;

  /*!
    @brief Tags
  */
  Vector<IntVectSet> m_geom_tags;

  /*!
    @brief Tag layouts
  */
  Vector<RefCountedPtr<LayoutData<IntVectSet> > > m_layout_tags;

  /*!
    @brief Potential
  */
  Real (*m_potential)(const Real a_time);

  /*!
    @brief Time increment
  */
  Real m_dt;

  /*!
    @brief Time
  */
  Real m_time;

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
    @brief Write a plot file
  */
  virtual void write_plot_file();

  /*!
    @brief Write a checkpoint file
  */
  virtual void write_checkpoint_file();

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
    @brief Read checkpoint file
  */
  virtual void read_checkpoint_file();

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
				   const int& a_finestLevel,
				   const Vector<DisjointBoxLayout>& a_grids);

  /*!
    @brief Write a grid report
  */
  virtual void grid_report();

  /*!
    @brief Do a regrid step. This regrid method regrids data inside each solver. The exception is the Poisson solver
    whose solution is recomputed. The same goes for stationary RTE methods. 
  */
  virtual void regrid();
  
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
    @brief Write a time step report
  */
  virtual void step_report(const Real a_start_time, const Real a_end_time, const int a_max_steps);

  /*!
    @brief Tag cells for refinement. This includes cells from geometries and the cell tagger
  */
  virtual void tag_cells(Vector<IntVectSet>& a_ivs);

  /*!
    @brief Get output variable names for a given output mode
  */
  virtual Vector<string> get_output_variable_names();
};
#endif
