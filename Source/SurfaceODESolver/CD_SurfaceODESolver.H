/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_SurfaceODESolver.H
  @brief  Declaration of a cut-cell ODE solver. 
  @author Robert Marskar
*/

#ifndef CD_SurfaceODESolver_H
#define CD_SurfaceODESolver_H

// Our includes
#include <CD_AmrMesh.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Surface ODE solver
  @details This is a basic solver that acts as an ODE solver on cut-cells. 
*/
template <size_t N = 1>
class SurfaceODESolver
{
public:
  /*!
    @brief Default constructor. Must subsequently set AmrMesh.
    @details Sets realm to primal and phase to phase::gas.
  */
  SurfaceODESolver();

  /*!
    @brief Full constructor. 
    @details Sets AmrMesh to input and sets realm to primal and phase to phase::gas.
    @param[in] a_amr      AmrMesh reference
  */
  SurfaceODESolver(const RefCountedPtr<AmrMesh>& a_amr);

  /*!
    @brief Disallowed copy constructor
  */
  SurfaceODESolver(const SurfaceODESolver&) = delete;

  /*!
    @brief Disallowed move constructor
  */
  SurfaceODESolver(const SurfaceODESolver&&) = delete;

  /*!
    @brief Destructor (does nothing).
  */
  virtual ~SurfaceODESolver();

  /*!
    @brief Disallowed assignment operator
  */
  SurfaceODESolver&
  operator=(const SurfaceODESolver&) = delete;

  /*!
    @brief Disallowed assignment operator
  */
  SurfaceODESolver&
  operator=(const SurfaceODESolver&&) = delete;

  /*!
    @brief Parse solver options
  */
  virtual void
  parseOptions() noexcept;

  /*!
    @brief Parse runtime solver options
  */
  virtual void
  parseRuntimeOptions() noexcept;

  /*!
    @brief Set the amr object
    @param[in] a_amrMesh AmrMesh object. 
  */
  virtual void
  setAmr(const RefCountedPtr<AmrMesh>& a_amrMesh) noexcept;

  /*!
    @brief Set the realm.
    @param[in] a_realm Realm name. 
  */
  virtual void
  setRealm(const std::string a_realm) noexcept;

  /*!
    @brief Set the solver name
    @param[in] a_name Solver name
  */
  virtual void
  setName(const std::string a_name) noexcept;

  /*!
    @brief Get the realm where the solver is registered.
    @return Returns m_realm
  */
  std::string
  getRealm() const noexcept;

  /*!
    @brief Set phase.
    @param[in] a_phase Phase.
  */
  virtual void
  setPhase(const phase::which_phase a_phase) noexcept;

  /*!
    @brief Get phase
    @return m_phase
  */
  virtual phase::which_phase
  getPhase() const noexcept;

  /*!
    @brief Set verbosity.
    @param[in] a_verbosity Verbosity level. 
  */
  virtual void
  setVerbosity(const int a_verbosity) noexcept;

  /*!
    @brief Get verbosity
    @return m_verbosity
  */
  virtual int
  getVerbosity() const noexcept;

  /*!
    @brief Set the time. 
    @param[in] a_step Time step number
    @param[in] a_time Time (in seconds)
    @param[in] a_dt   Time step increment
    @note This sets m_step=a_step, m_time=a_time, m_dt=a_dt
  */
  virtual void
  setTime(const int a_step, const Real a_time, const Real a_dt) noexcept;

  /*!
    @brief Get current time step
    @return m_step
  */
  virtual Real
  getStep() const noexcept;

  /*!
    @brief Get current time
    @return m_time
  */
  virtual Real
  getTime() const noexcept;

  /*!
    @brief Get last time step
    @return m_dt
  */
  virtual Real
  getDt() const noexcept;

  /*!
    @brief Convenience function for setting m_phi
    @param[in] a_phi Values for all components
  */
  virtual void
  setPhi(const Real a_phi);

  /*!
    @brief Convenience function for setting m_phi
    @param[in] a_phi Component-wise values. 
  */
  virtual void
  setPhi(const std::array<Real, N>& a_phi);

  /*!
    @brief Convenience function for setting m_phi
    @param[in] a_phi Values
    @note a_phi must have N components. 
  */
  virtual void
  setPhi(const EBAMRIVData& a_phi);

  /*!
    @brief Get internal state
    @return Returns m_phi
  */
  virtual EBAMRIVData&
  getPhi() noexcept;

  /*!
    @brief Get internal state
    @return Returns m_phi
  */
  virtual const EBAMRIVData&
  getPhi() const noexcept;

  /*!
    @brief Convenience function for setting m_rhs
    @param[in] a_rhs Values for all components
  */
  virtual void
  setRHS(const Real a_rhs);

  /*!
    @brief Convenience function for setting m_rhs
    @param[in] a_rhs Component-wise values. 
  */
  virtual void
  setRHS(const std::array<Real, N>& a_rhs);

  /*!
    @brief Convenience function for setting m_rhs
    @param[in] a_rhs Values per cell and component. 
    @note a_rhs must have N components. 
  */
  virtual void
  setRHS(const EBAMRIVData& a_rhs);

  /*!
    @brief Get internal state
    @return Returns m_rhs
  */
  virtual EBAMRIVData&
  getRHS();

  /*!
    @brief Get internal state
    @return Returns m_rhs
  */
  virtual const EBAMRIVData&
  getRHS() const noexcept;

  /*!
    @brief Allocate internal storage for this class. 
    @details Allocates m_phi and m_rhs
  */
  virtual void
  allocate() noexcept;

  /*!
    @brief Deallocate internal storage.
  */
  virtual void
  deallocate() noexcept;

  /*!
    @brief Register operators. 
  */
  virtual void
  registerOperators() noexcept;

  /*!
    @brief Pre-regrid function. 
    @details This stores the data on the old mesh so it can be regridded later.
    @param[in] a_lbase          Coarsest level which will change during regrids. 
    @param[in] a_oldFinestLevel Finest level before the regrid operation. 
  */
  virtual void
  preRegrid(const int a_lbase, const int a_oldFinestLevel) noexcept;

  /*!
    @brief Regrid function.
    @param[in] a_lmin           Coarsest level where grids did not change. 
    @param[in] a_oldFinestLevel Finest AMR level before the regrid. 
    @param[in] a_newFinestLevel Finest AMR level after the regrid. 
    @details This interpolates or coarsens conservatively, e.g. sigma_c = sum(A_f * sigma_f)/A_c if we coarsen.
  */
  virtual void
  regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) noexcept;

  /*!
    @brief Reset m_phi on electrode cells. 
    @details This will set a_data to the specified value (in electrode cells)
    @param[in] a_val  Value
    @note Does not include ghost cells. 
  */
  virtual void
  resetElectrodes(const Real a_value) noexcept;

  /*!
    @brief Reset the input data holder on electrode cells. 
    @details This will set a_data to the specified value (in electrode cells)
    @param[in] a_data Input data to be set.
    @param[in] a_val  Value
    @note Does not include ghost cells. 
  */
  virtual void
  resetElectrodes(EBAMRIVData& a_phi, const Real a_value) const noexcept;

  /*!
    @brief Reset m_phi on dielectric cells. 
    @details This will set a_data to the specified value (in dielectric cells)
    @param[in] a_val  Value
    @note Does not include ghost cells. 
  */
  virtual void
  resetDielectrics(const Real a_value) noexcept;

  /*!
    @brief Reset the input data holder on dielectric cells. 
    @details This will set a_data to the specified value (in dielectric cells)
    @param[in] a_data Input data to be set.
    @param[in] a_val  Value
    @note Does not include ghost cells. 
  */
  virtual void
  resetDielectrics(EBAMRIVData& a_phi, const Real a_value) const noexcept;

#ifdef CH_USE_HDF5
  /*!
    @brief Write checkpoint data into handle
    @paramo[out] a_handle HDF5 file. 
    @param[in]   a_level Grid level
  */
  virtual void
  writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const noexcept;
#endif

#ifdef CH_USE_HDF5
  /*!
    @brief Read checkpoint data from HDF5 file. 
    @param[in] a_handle HDF5 handle.
    @param[in] a_level  Grid level
  */
  virtual void
  readCheckpointLevel(HDF5Handle& a_handle, const int a_level) noexcept;
#endif

#ifdef CH_USE_HDF5
  /*!
    @brief Write a plot file
    @details The name of the plot file is m_name.stepXXXXXXX.DIM.hdf5
  */
  virtual void
  writePlotFile() const noexcept;
#endif

  /*!
    @brief Write output data to a_output
    @param[inout] a_output Output data holder.
    @param[inout] a_icomp  Starting component where this solver begins writing the output. 
    @note This will write the plot data in this solver to a_output, starting on a_comp
    @note This routine writes m_phi and m_rhs into a_ouput. 
  */
  virtual void
  writePlotData(EBAMRCellData& a_output, int& a_comp) const noexcept;

  /*!
    @brief Get number of variables to be plotted. 
  */
  virtual int
  getNumberOfPlotVariables() const noexcept;

  /*!
    @brief Get output plot names
    @return Return a list of plot variable names. 
  */
  virtual Vector<std::string>
  getPlotVariableNames() const noexcept;

protected:
  /*!
    @brief Realm where this solver is registered. 
  */
  std::string m_realm;

  /*!
    @brief Class name
  */
  std::string m_className;

  /*!
    @brief Solver name
  */
  std::string m_name;

  /*!
    @brief Phase where the data lives. 
  */
  phase::which_phase m_phase;

  /*!
    @brief Conservative regrid or not
  */
  bool m_conservativeRegrid;

  /*!
    @brief Plot phi or not
  */
  bool m_plotPhi;

  /*!
    @brief Plot rhs or not
  */
  bool m_plotRHS;

  /*!
    @brief Index spaces
  */
  RefCountedPtr<MultiFluidIndexSpace> m_mfis;

  /*!
    @brief AMR; needed for grid stuff
  */
  RefCountedPtr<AmrMesh> m_amr;

  /*!
    @brief Cut-cells that are electrode cells. This does not include ghost cells. 

  */
  Vector<RefCountedPtr<LayoutData<VoFIterator>>> m_electrodeVoFs;

  /*!
    @brief Cut-cells that are also dielectric cells. This does not include ghost cells. 
  */
  Vector<RefCountedPtr<LayoutData<VoFIterator>>> m_dielectricVoFs;

  /*!
    @brief State. Stored on cut-cells. 
  */
  EBAMRIVData m_phi;

  /*!
    @brief Cached state used during regridding.
  */
  EBAMRIVData m_cache;

  /*!
    @brief Right-hand side.
  */
  EBAMRIVData m_rhs;

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
    @brief Parse verbosity
  */
  void
  parseVerbosity() noexcept;

  /*!
    @brief Parse regrid method
  */
  void
  parseRegrid() noexcept;

  /*!
    @brief Parse plot variables
  */
  void
  parsePlotVariables() noexcept;

  /*!
    @brief Define iterators for iterating over cell subsets
  */
  void
  defineVoFIterators() noexcept;
};

#include <CD_NamespaceFooter.H>

#include <CD_SurfaceODESolverImplem.H>

#endif
