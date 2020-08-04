/*!
  @file sigma_solver.H
  @brief  Declaration of a surface charge solver
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _SIGMA_SOLVER_
#define _SIGMA_SOLVER_

#include "computational_geometry.H"
#include "amr_mesh.H"

/*!
  @brief Surface charge solver
*/
class sigma_solver {
public:

  /*!
    @brief Constructor
  */
  sigma_solver();

  /*!
    @brief Destructor
  */
  ~sigma_solver();

  /*!
    @brief Get the realm
  */
  const std::string get_realm() const;

  /*!
    @brief Set the realm
  */
  virtual void set_realm(const std::string a_realm); 

  /*!
    @brief Allocate internal storage
  */
  virtual void allocate_internals();

  /*!
    @brief Cache state for regridding
  */
  virtual void pre_regrid(const int a_lbase, const int a_old_finest_level);

  /*!
    @brief Compute the right-hand side
  */
  virtual void compute_rhs(EBAMRIVData& a_rhs);

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocate_internals();

  /*!
    @brief Regrid method
  */
  virtual void regrid(const int a_lmin, const int a_old_finest_level, const int a_new_finest_level);

  /*!
    @brief Register operators
  */
  virtual void register_operators();
  
  /*!
    @brief Reset data holder on cells
  */
  virtual void reset_cells(EBAMRIVData& a_data);

  /*!
    @brief Set the amr object
  */
  virtual void set_amr(const RefCountedPtr<amr_mesh>& a_amr);

  /*!
    @brief Set computational geometry
  */
  virtual void set_computational_geometry(const RefCountedPtr<computational_geometry>& a_compgeom);

  /*!
    @brief Set phase
    @details This must be done BEFORE callilng set_computational_geometry
  */
  virtual void set_phase(phase::which_phase a_phase = phase::gas);

  /*!
    @brief Convenience function. Set the surface charge
  */
  virtual void set_sigma(const EBAMRIVData& a_sigma);

  /*!
    @brief Convenience function. Set the surface charge. 
  */
  virtual void set_sigma(const Real a_sigma);

  /*!
    @brief Set verbosity.
  */
  virtual void set_verbosity(const int a_verbosity);

  /*!
    @brief Get current time
  */
  virtual void set_time(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Write checkpoint data into handle
  */
  virtual void write_checkpoint_level(HDF5Handle& a_handle, const int a_level) const;

  /*!
    @brief Read checkpoint data from handle
  */
  virtual void read_checkpoint_level(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Write output data to a_output
  */
  virtual void write_plot_data(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Set variables to be plotted
  */
  virtual void set_plot_variables();

  /*!
    @brief Number of variables to be plotted
  */
  virtual int get_num_plotvars();

  /*!
    @brief Get output plot names
  */
  virtual Vector<std::string> get_plotvar_names() const;

  /*!
    @brief Compute total charge
  */
  virtual Real compute_charge();

  /*!
    @brief Get internal state
  */
  virtual EBAMRIVData& get_state();

  /*!
    @brief Get internal state
  */
  virtual EBAMRIVData& get_flux();
  
protected:

  /*!
    @brief Realm
  */
  std::string m_realm;

  /*!
    @brief Index spaces
  */
  RefCountedPtr<mfis> m_mfis;

  /*!
    @brief Computational geometry
  */
  RefCountedPtr<computational_geometry> m_compgeom;

  /*!
    @brief AMR; needed for grid stuff
  */
  RefCountedPtr<amr_mesh> m_amr;

  /*!
    @brief State
  */
  EBAMRIVData m_state;

  /*!
    @brief Cached state for regridding. DO NOT TOUCH!
  */
  EBAMRIVData m_cache;

  /*!
    @brief Flux
  */
  EBAMRIVData m_flux;

  /*!
    @brief Phase
  */
  phase::which_phase m_phase;

  /*!
    @brief Solver verbosity
  */
  int m_verbosity;

  /*!
    @brief Time step
  */
  int m_step;

  /*!
    @brief Current time
  */
  Real m_time;

  /*!
    @brief Time increment
  */
  Real m_dt;

  /*!
    @brief Output sigma
  */
  bool m_plot_phi;
  
};
#endif

