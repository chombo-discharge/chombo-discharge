/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MeshODESolverImplem.H
  @brief  Implementation of CD_MeshODESolver.H
  @author Robert Marskar
*/

#ifndef CD_MeshODESolverImplem_H
#define CD_MeshODESolverImplem_H

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_MeshODESolver.H>

template <size_t N>
MeshODESolver<T>::MeshODESolver()
{
  CH_TIME("MeshODESolver::MeshODESolver");

  // Default settings. 
  m_verbosity = -1;
  m_name      = "MeshODESolver";
  m_className = "MeshODESolver";
  m_phase     = phase::gas;
  m_realm     = Realm::Primal;

  this->parseOptions();
}

template <size_t N>
MeshODESolver<T>::~MeshODESolver()
{
  CH_TIME("MeshODESolver::~MeshODESolver");
}

template <size_t N>
void
MeshODESolver<T>::parseOptions() noexcept
{
  CH_TIME("MeshODESolver::parseOptions");

  ParmParse pp(m_className.c_str());
  
  pp.get("verbosity", m_verbosity);
  pp.get("use_regrid_slopes", m_regridSlopes);
  
  if(m_verbosity > 5) {
    pout() << m_name + "::parseOptions()" << endl;
  }
}

template <size_t N>
void
MeshODESolver<T>::parseRuntimeOptions() noexcept
{
  CH_TIME("MeshODESolver::parseRuntimeOptions()");
  if(m_verbosity > 5) {
    pout() << m_name + "::parseRuntimeOptions()" << endl;
  }
  
  ParmParse pp(m_className.c_str());
  
  pp.get("verbosity", m_verbosity);
  pp.get("use_regrid_slopes", m_regridSlopes);  
}

template <size_t N>
EBAMRCellData&
MeshODESolver<T>::getPhi() noexcept
{
  CH_TIME("MeshODESolver::getPhi()");
  if(m_verbosity > 5) {
    pout() << m_name + "::getPhi()" << endl;
  }
  
  return m_phi;
}

template <size_t N>
const EBAMRCellData&
MeshODESolver<T>::getPhi() const noexcept
{
  CH_TIME("MeshODESolver::getPhi()");
  if(m_verbosity > 5) {
    pout() << m_name + "::getPhi()" << endl;
  }
  
  return m_phi;
}

template <size_t N>
EBAMRCellData&
MeshODESolver<T>::getRHS() noexcept
{
  CH_TIME("MeshODESolver::getRHS()");
  if(m_verbosity > 5) {
    pout() << m_name + "::getRHS()" << endl;
  }
  
  return m_RHS;
}

template <size_t N>
const EBAMRCellData&
MeshODESolver<T>::getRHS() const noexcept
{
  CH_TIME("MeshODESolver::getRHS()");
  if(m_verbosity > 5) {
    pout() << m_name + "::getRHS()" << endl;
  }
  
  return m_RHS;
}

template <size_t N>
void
MeshODESolver<T>::preRegrid(const int a_lbase, const int a_oldFinestLevel) noexcept {
  CH_TIME("MeshODESolver::preRegrid(int, int)");
  if(m_verbosity > 5) {
    pout() << m_name + "::preRegrid(int, int)" << endl;
  }

  m_amr->allocate(m_cache, m_realm, m_phase, N);

  m_cache.copy(m_phi);
}

template <size_t N>
void
MeshODESolver<T>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) noexcept
{
  CH_TIME("MeshODESolver::regrid(int, int, int)");
  if(m_verbosity > 5) {
    pout() << m_name + "::regrid(int, int, int)" << endl;
  }
}

template <size_t N>
std::string
MeshODESolver<T>::getRealm() const noexcept
{
  CH_TIME("MeshODESolver::getRealm()");
  if(m_verbosity > 5) {
    pout() << m_name + "::getRealm()" << endl;
  }

  return m_realm;
}

template <size t N>
std::string
MeshODESolver<T>::getName() const noexcept
{
  CH_TIME("MeshODESolver::getName()");
  if(m_verbosity > 5) {
    pout() << m_name + "::getName()" << endl;
  }

  return m_name;
}

  /*!
    @brief Get output plot names
    @return Return a list of plot variable names. 
  */
  virtual Vector<std::string>
  getPlotVariableNames() const;

  /*!
    @brief Get number of output fields
    @return Returns number of plot variables include in writePlotFile() and writePlotData()
  */
  virtual int
  getNumberOfPlotVariables() const;

  /*!
    @brief Register operators for AMR operations
    @note This includes operators for redistribution, flux registers, regridding, ghost cell interpolation, and conservative coarsening. 
  */
  virtual void
  registerOperators();

  /*!
    @brief Set the realm for this solver.
    @param[in] a_realm Realm identifier
  */
  virtual void
  setRealm(const std::string a_realm);

  /*!
    @brief Set phase.
    @param[in] a_phase Phase.
  */
  virtual void
  setPhase(phase::which_phase a_phase);

  /*!
    @brief Set verbosity.
    @param[in] a_verbosity Verbosity level. 
  */
  virtual void
  setVerbosity(const int a_verbosity);

  /*!
    @brief Set the time for this solver. 
    @param[in] a_step Time step number
    @param[in] a_time Time (in seconds)
    @param[in] a_dt   Time step increment
    @note This sets m_step=a_step, m_time=a_time, m_dt=a_dt
  */
  virtual void
  setTime(const int a_step, const Real a_time, const Real a_dt);

  /*!
    @brief Write plot file
    @details The name of the plot file is m_name.stepXXXXX.DIM.hdf5
    @note This calls writePlotData(...)
  */
  virtual void
  writePlotFile();

  /*!
    @brief Write output data to a_output
    @param[inout] a_output Output data holder.
    @param[inout] a_icomp  Starting component where this solver begins writing the output. 
    @note This will write the plot data in this solver to a_output, starting on a_comp
    @note This routine writes m_phi on centroids (and not cell-centers). 
  */
  virtual void
  writePlotData(EBAMRCellData& a_output, int& a_icomp);

  /*!
    @brief Write checkpoint data into HDF5 file. 
    @paramo[out] a_handle HDF5 file. 
    @param[in]   a_level Grid level
  */
#ifdef CH_USE_HDF5
  virtual void
  writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const;
#endif

  /*!
    @brief Read checkpoint data from HDF5 file. 
    @param[in] a_handle HDF5 handle.
    @param[in] const int a_level Grid level
  */
#ifdef CH_USE_HDF5
  virtual void
  readCheckpointLevel(HDF5Handle& a_handle, const int a_level);
#endif  

#include <CD_NamespaceFooter.H>

#endif
