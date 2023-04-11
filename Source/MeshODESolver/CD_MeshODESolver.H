/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MeshODESolver.H
  @brief  Encapsulation of an ODE solver on the mesh. 
  @author Robert Marskar
*/

#ifndef CD_MeshODESolver_H
#define CD_MeshODESolver_H

// Std includes
#include <array>

// Our includes
#include <CD_AmrMesh.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for solving dy/dt = f on an AMR hierarchy. 
  @details The template parameter is the number of variables in y and f.
*/
template <size_t N = 1>
class MeshODESolver
{
public:
  /*!
    @brief Alias for right-hand side
  */
  using RHSFunction = std::function<std::array<Real, N>(const std::array<Real, N>&, const Real)>;

  /*!
    @brief Default constructor. Must subsequently set everything through public member functions.
  */
  MeshODESolver();

  /*!
    @brief Full constructor. 
    @param[in] a_amr AMR core.
  */
  MeshODESolver(const RefCountedPtr<AmrMesh>& a_amr) noexcept;

  /*!
    @brief Disallowed copy constructor
  */
  MeshODESolver(const MeshODESolver&) = delete;

  /*!
    @brief Disallowed move constructor
  */
  MeshODESolver(const MeshODESolver&&) = delete;

  /*!
    @brief Destructor
  */
  virtual ~MeshODESolver();

  /*!
    @brief Disallowed assignment operator
  */
  MeshODESolver&
  operator=(const MeshODESolver&) = delete;

  /*!
    @brief Disallowed move assignement operator
  */
  MeshODESolver&
  operator=(const MeshODESolver&&) = delete;

  /*!
    @brief Set AmrMesh.
    @param[in] a_amrMesh AmrMesh
  */
  virtual void
  setAmr(const RefCountedPtr<AmrMesh>& a_amrMesh) noexcept;

  /*!
    @brief Parse class options.
  */
  virtual void
  parseOptions() noexcept;

  /*!
    @brief Parse run-time configurable class options.
  */
  virtual void
  parseRuntimeOptions() noexcept;

  /*!
    @brief Allocate internal storage
  */
  virtual void
  allocate() noexcept;

  /*!
    @brief Set phi for a specific component. 
    @param[in] a_phiFunc Phi function
    @param[in] a_comp    Component
  */
  virtual void
  setPhi(const std::function<Real(const RealVect& a_pos)>& a_phiFunc, const size_t a_comp) noexcept;

  /*!
    @brief Set phi everywhere. 
    @param[in] a_phiFunc Phi function. 
  */
  virtual void
  setPhi(const std::function<std::array<Real, N>(const RealVect& a_pos)>& a_phiFunc) noexcept;

  /*!
    @brief Set right-hand side for specified component. 
    @param[in] a_srcFunc Source term function.
    @param[in] a_comp    Component
  */
  virtual void
  setRHS(const std::function<Real(const RealVect& a_pos)>& a_srcFunc, const size_t a_comp) noexcept;

  /*!
    @brief Compute right-hand side from left-hand side. I.e. compute f = f(y,t). 
    @param[in] a_rhsFunction Function for computing the right-hand side. 
  */
  virtual void
  computeRHS(const RHSFunction& a_rhsFunction) noexcept;

  /*!
    @brief Compute right-hand side from left-hand side. I.e. compute f = f(y,t). 
    @param[out] a_rhs         Right-hand side on the mesh. 
    @param[in]  a_rhsFunction Function for computing the right-hand side. 
  */
  virtual void
  computeRHS(EBAMRCellData& rhs, const RHSFunction& a_rhsFunction) const noexcept;

  /*!
    @brief Get the solution vector (left-hand side of equation).
  */
  EBAMRCellData&
  getPhi() noexcept;

  /*!
    @brief Get the solution vector (left-hand side of equation).
  */
  const EBAMRCellData&
  getPhi() const noexcept;

  /*!
    @brief Get the solution vector (left-hand side of equation).
  */
  EBAMRCellData&
  getRHS() noexcept;

  /*!
    @brief Get the solution vector (left-hand side of equation).
  */
  const EBAMRCellData&
  getRHS() const noexcept;

  /*!
    @brief Perform pre-regrid operations. 
    @param[in] a_lbase          Coarsest level that changed during regrid. 
    @param[in] a_oldFinestLevel Finest grid level before the regrid operation. 
    @note This copies m_phi onto m_cache
  */
  virtual void
  preRegrid(const int a_lbase, const int a_oldFinestLevel) noexcept;

  /*!
    @brief Regrid this solver. 
    @param[in] a_lmin           Coarsest level where grids did not change. 
    @param[in] a_oldFinestLevel Finest AMR level before the regrid. 
    @param[in] a_newFinestLevel Finest AMR level after the regrid. 
    @details This linearly interpolates (potentially with limiters) m_phi to the new grids. 
  */
  virtual void
  regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) noexcept;

  /*!
    @brief Get the realm where this solver is registered. 
    @return Returns realm name. 
  */
  virtual std::string
  getRealm() const noexcept;

  /*!
    @brief Get solver name
    @note Not necessarily equal to class name (we can have many CdrSolvers instantiated). 
  */
  virtual std::string
  getName() const noexcept;

  /*!
    @brief Set solver name
    @param[in] a_name Solver name
  */
  virtual void
  setName(const std::string& a_name) noexcept;

  /*!
    @brief Get output plot names
    @return Return a list of plot variable names. 
  */
  virtual Vector<std::string>
  getPlotVariableNames() const noexcept;

  /*!
    @brief Get number of output fields
    @return Returns number of plot variables include in writePlotFile() and writePlotData()
  */
  virtual int
  getNumberOfPlotVariables() const noexcept;

  /*!
    @brief Register operators for AMR operations
    @note This includes operators for redistribution, flux registers, regridding, ghost cell interpolation, and conservative coarsening. 
  */
  virtual void
  registerOperators() const noexcept;

  /*!
    @brief Set the realm for this solver.
    @param[in] a_realm Realm identifier
  */
  virtual void
  setRealm(const std::string a_realm) noexcept;

  /*!
    @brief Set phase.
    @param[in] a_phase Phase.
  */
  virtual void
  setPhase(phase::which_phase a_phase) noexcept;

  /*!
    @brief Set verbosity.
    @param[in] a_verbosity Verbosity level. 
  */
  virtual void
  setVerbosity(const int a_verbosity) noexcept;

  /*!
    @brief Set the time for this solver. 
    @param[in] a_step Time step number
    @param[in] a_time Time (in seconds)
    @param[in] a_dt   Time step increment
    @note This sets m_step=a_step, m_time=a_time, m_dt=a_dt
  */
  virtual void
  setTime(const int a_step, const Real a_time, const Real a_dt) noexcept;

  /*!
    @brief Write plot file
    @details The name of the plot file is m_name.stepXXXXX.DIM.hdf5
    @note This calls writePlotData(...)
  */
  virtual void
  writePlotFile() const noexcept;

  /*!
    @brief Write output data to a_output
    @param[inout] a_output Output data holder.
    @param[inout] a_icomp  Starting component where this solver begins writing the output. 
    @note This will write the plot data in this solver to a_output, starting on a_comp
    @note This routine writes m_phi on centroids (and not cell-centers). 
  */
  virtual void
  writePlotData(EBAMRCellData& a_output, int& a_icomp) const noexcept;

  /*!
    @brief Write checkpoint data into HDF5 file. 
    @paramo[out] a_handle HDF5 file. 
    @param[in]   a_level Grid level
  */
#ifdef CH_USE_HDF5
  virtual void
  writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const noexcept;
#endif

  /*!
    @brief Read checkpoint data from HDF5 file. 
    @param[in] a_handle HDF5 handle.
    @param[in] const int a_level Grid level
  */
#ifdef CH_USE_HDF5
  virtual void
  readCheckpointLevel(HDF5Handle& a_handle, const int a_level) noexcept;
#endif

protected:
  /*!
    @brief AMR mesh reference
  */
  RefCountedPtr<AmrMesh> m_amr;

  /*!
    @brief Chattiness
  */
  int m_verbosity;

  /*!
    @brief Time step
  */
  int m_timeStep;

  /*!
    @brief Current time
  */
  Real m_time;

  /*!
    @brief Last time step increment. 
  */
  Real m_dt;

  /*!
    @brief Plot solution vector or not. 
  */
  bool m_plotPhi;

  /*!
    @brief Plot source term or not.
  */
  bool m_plotRHS;

  /*!
    @brief Use slopes when regridding
  */
  bool m_regridSlopes;

  /*!
    @brief Solver name
  */
  std::string m_name;

  /*!
    @brief Class name. 
    @details This will be different for different implementations of CdrSolver
  */
  std::string m_className;

  /*!
    @brief Realm where solver is registered
  */
  std::string m_realm;

  /*!
    @brief Phase where solver data lives. 
  */
  phase::which_phase m_phase;

  /*!
    @brief Solution vector (left hand side of equation). 
  */
  EBAMRCellData m_phi;

  /*!
    @brief Source term (right-hand side of equation). 
  */
  EBAMRCellData m_rhs;

  /*!
    @brief Cache storage when doing regrids.
  */
  EBAMRCellData m_cache;

  /*!
    @brief Write data to output. Convenience function. 
    @param[inout] a_output Output data holder.
    @param[inout] a_icomp  Starting component where this solver begins writing the output. 
    @param[in]    a_data   Data to write. 
    @param[in]    a_interp If true, a_data will be interpolated to cell centroids before writing to a_output. 
  */
  virtual void
  writeData(EBAMRCellData& a_output, int& a_comp, const EBAMRCellData& a_data, const bool a_interp) const noexcept;

  /*!
    @brief Parse plot variables
  */
  virtual void
  parsePlotVariables() noexcept;
};

#include <CD_NamespaceFooter.H>

#include <CD_MeshODESolverImplem.H>

#endif
