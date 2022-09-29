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
#include <CD_ComputationalGeometry.H>
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
    @brief Default constructor.
  */
  SurfaceODESolver();

  /*!
    @brief Destructor (does nothing).
  */
  virtual ~SurfaceODESolver();

  /*!
    @brief Get the realm where the solver is registered.
    @return Returns m_realm
  */
  const std::string
  getRealm() const;

  /*!
    @brief Set the realm.
    @param[in] a_realm Realm name. 
  */
  virtual void
  setRealm(const std::string a_realm);

  /*!
    @brief Set phase.
    @param[in] a_phase Phase.
    @details This sets the phase where all the data will live. Note that the dielectric cells can be on the "inside" or "outside" of the dielectric and
    this function thus controls where the data is allocated. That might seem silly, but note that this essentially determines from which side the flux comes. We
    have this design because we consider the surface charge solver be a single-phase solver (in the future, this might change). 
    @note This must be done BEFORE calling setComputationalGeometry
  */
  virtual void
  setPhase(phase::which_phase a_phase);

  /*!
    @brief Allocate internal storage for this class. 
  */
  virtual void
  allocateInternals();

  /*!
    @brief Compute the right-hand side (i.e. the flux).
    @details This just copies the input flux into the the internal data holder for the flux. 
    @param[in] a_flux Incoming charge flux
  */
  virtual void
  computeRHS(EBAMRIVData& a_flux);

  /*!
    @brief Deallocate internal storage.
  */
  virtual void
  deallocateInternals();

  /*!
    @brief Register operators. 
  */
  virtual void
  registerOperators();  

  /*!
    @brief Pre-regrid function. 
    @details This stores the data on the old mesh so it can be regridded later.
    @param[in] a_lbase          Coarsest level which will change during regrids. 
    @param[in] a_oldFinestLevel Finest level before the regrid operation. 
  */
  virtual void
  preRegrid(const int a_lbase, const int a_oldFinestLevel);  

  /*!
    @brief Regrid function.
    @param[in] a_lmin           Coarsest level where grids did not change. 
    @param[in] a_oldFinestLevel Finest AMR level before the regrid. 
    @param[in] a_newFinestLevel Finest AMR level after the regrid. 
    @details This interpolates or coarsens conservatively, e.g. sigma_c = sum(A_f * sigma_f)/A_c if we coarsen.
  */
  virtual void
  regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel);

  /*!
    @brief Reset data holder on non-dielectric cells. 
    @details This will set a_data to zero except in cells that are dielectric surface cells. 
    @param[in] a_data Input data to be reset. 
  */
  virtual void
  resetCells(EBAMRIVData& a_data);

  /*!
    @brief Set the amr object
    @param[in] a_amrMesh AmrMesh object. 
  */
  virtual void
  setAmr(const RefCountedPtr<AmrMesh>& a_amrMesh);

  /*!
    @brief Set the computational geometry
    @param[in] a_computationalGeometry The computational geometry.
  */
  virtual void
  setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);


  /*!
    @brief Convenience function. Set the surface charge equal to the input argument. 
    @param[in] a_sigma Surface charge.
  */
  virtual void
  setPhi(const EBAMRIVData& a_sigma);

  /*!
    @brief Convenience function. Set the surface charge equal to the input argument. 
    @param[in] a_sigma Surface charge.
  */
  virtual void
  setPhi(const Real a_sigma);

  /*!
    @brief Set verbosity.
    @param[in] a_verbosity Verbosity level. 
  */
  virtual void
  setVerbosity(const int a_verbosity);

  /*!
    @brief Set the time. 
    @param[in] a_step Time step number
    @param[in] a_time Time (in seconds)
    @param[in] a_dt   Time step increment
    @note This sets m_step=a_step, m_time=a_time, m_dt=a_dt
  */
  virtual void
  setTime(const int a_step, const Real a_time, const Real a_dt);

#ifdef CH_USE_HDF5
  /*!
    @brief Write checkpoint data into handle
    @paramo[out] a_handle HDF5 file. 
    @param[in]   a_level Grid level
  */
  virtual void
  writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const;
#endif

#ifdef CH_USE_HDF5
  /*!
    @brief Read checkpoint data from HDF5 file. 
    @param[in] a_handle HDF5 handle.
    @param[in] const int a_level Grid level
  */
  virtual void
  readCheckpointLevel(HDF5Handle& a_handle, const int a_level);
#endif

  /*!
    @brief Write output data to a_output
    @param[inout] a_output Output data holder.
    @param[inout] a_icomp  Starting component where this solver begins writing the output. 
    @note This will write the plot data in this solver to a_output, starting on a_comp
    @note This routine writes m_phi and m_rhs into a_ouput. 
  */
  virtual void
  writePlotData(EBAMRCellData& a_output, int& a_comp);

  /*!
    @brief Number of variables to be plotted. 
    @details Returns 2 - the surface charge and the incoming flux. 
  */
  virtual int
  getNumberOfPlotVariables();

  /*!
    @brief Get output plot names
    @return Return a list of plot variable names. 
  */
  virtual Vector<std::string>
  getPlotVariableNames() const;

  /*!
    @brief Compute total charge
    @return Computes the total charge on the dielectric surface. 
  */
  virtual Real
  computeCharge();

  /*!
    @brief Get internal state
    @return Returns m_phi
  */
  virtual EBAMRIVData&
  getPhi();

  /*!
    @brief Get internal state
    @return Returns m_rhs
  */
  virtual EBAMRIVData&
  getFlux();

protected:
  /*!
    @brief Component number in data holder
  */
  static constexpr int m_comp = 0;

  /*!
    @brief Number of components that this solver solves for
  */
  static constexpr int m_nComp = 1;

  /*!
    @brief Realm where this solver is registered. 
  */
  std::string m_realm;

  /*!
    @brief Index spaces
  */
  RefCountedPtr<MultiFluidIndexSpace> m_multifluidIndexSpace;

  /*!
    @brief Computational geometry
  */
  RefCountedPtr<ComputationalGeometry> m_computationalGeometry;

  /*!
    @brief AMR; needed for grid stuff
  */
  RefCountedPtr<AmrMesh> m_amr;

  /*!
    @brief State -- this holds the surface charge density. 
  */
  EBAMRIVData m_phi;

  /*!
    @brief Cached state for regridding. DO NOT TOUCH!
  */
  EBAMRIVData m_cache;

  /*!
    @brief Surface charge flux. 
  */
  EBAMRIVData m_rhs;

  /*!
    @brief Phase where the data lives. 
  */
  phase::which_phase m_phase;

  /*!
    @brief Solver verbosity
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
    @brief Time increment
  */
  Real m_dt;
};

#include <CD_NamespaceFooter.H>

#include <CD_SurfaceODESolverImplem.H>

#endif
