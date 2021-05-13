/*!
  @file sigma_solver.H
  @brief  Declaration of a surface charge solver
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _SIGMA_SOLVER_
#define _SIGMA_SOLVER_

#include "computational_geometry.H"
#include "CD_AmrMesh.H"

#include "CD_NamespaceHeader.H"
  
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
    @brief Get the Realm
  */
  const std::string get_Realm() const;

  /*!
    @brief Set the Realm
  */
  virtual void set_Realm(const std::string a_Realm); 

  /*!
    @brief Allocate internal storage
  */
  virtual void allocateInternals();

  /*!
    @brief Cache state for regridding
  */
  virtual void pre_regrid(const int a_lbase, const int a_oldFinestLevel);

  /*!
    @brief Compute the right-hand side
  */
  virtual void compute_rhs(EBAMRIVData& a_rhs);

  /*!
    @brief Deallocate internal storage
  */
  virtual void deallocateInternals();

  /*!
    @brief Regrid method
  */
  virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Register operators
  */
  virtual void registerOperators();
  
  /*!
    @brief Reset data holder on cells
  */
  virtual void reset_cells(EBAMRIVData& a_data);

  /*!
    @brief Set the amr object
  */
  virtual void setAmr(const RefCountedPtr<AmrMesh>& a_amr);

  /*!
    @brief Set computational geometry
  */
  virtual void setComputationalGeometry(const RefCountedPtr<computational_geometry>& a_computationalGeometry);

  /*!
    @brief Set phase
    @details This must be done BEFORE callilng setComputationalGeometry
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
  virtual void writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const;

  /*!
    @brief Read checkpoint data from handle
  */
  virtual void readCheckpointLevel(HDF5Handle& a_handle, const int a_level);

  /*!
    @brief Write output data to a_output
  */
  virtual void writePlotData(EBAMRCellData& a_output, int& a_comp);

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
  std::string m_Realm;

  /*!
    @brief Index spaces
  */
  RefCountedPtr<mfis> m_multifluidIndexSpace;

  /*!
    @brief Computational geometry
  */
  RefCountedPtr<computational_geometry> m_computationalGeometry;

  /*!
    @brief AMR; needed for grid stuff
  */
  RefCountedPtr<AmrMesh> m_amr;

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
#include "CD_NamespaceFooter.H"
#endif

