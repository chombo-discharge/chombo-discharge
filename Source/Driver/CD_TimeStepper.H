/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TimeStepper.H
  @brief  Declaration of main (abstract) time stepper class. 
  @author Robert Marskar
*/

#ifndef CD_TimeStepper_H
#define CD_TimeStepper_H

// Chombo includes
#include <CH_HDF5.H>

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_MultiFluidIndexSpace.H>
#include <CD_AmrMesh.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for advancing equations. 
  @details This class is used by Driver for advancing sets of equations. In short, this class should implement a time-stepping routine for a
  set of solvers. One should also implement routines for setting up the solvers, allocating necessary memory, regridding etc. 
*/
class TimeStepper
{
public:
  /*!
    @brief Default constructor (does nothing)
  */
  TimeStepper();

  /*!
    @brief Default destructor (does nothing)
  */
  virtual ~TimeStepper();

  /*!
    @brief Set AmrMesh
    @param[in] a_amr AmrMesh
  */
  void
  setAmr(const RefCountedPtr<AmrMesh>& a_amr);

  /*!
    @brief Set the computational geometry
    @param[in] a_computationalGeometry The computational geometry.
  */
  void
  setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

  /*!
    @brief Set up solvers
  */
  virtual void
  setupSolvers() = 0;

  /*!
    @brief Allocate data for the time stepper and solvers. 
  */
  virtual void
  allocate() = 0;

  /*!
    @brief Fill solvers with initial data
  */
  virtual void
  initialData() = 0;

  /*!
    @brief Post-initialize operations to be performed at end of setup stage. 
  */
  virtual void
  postInitialize() = 0;

  /*!
    @brief Post-initialize operations to be performed after filling solvers with data read from checkpoint files. 
  */
  virtual void
  postCheckpointSetup() = 0;

  /*!
    @brief Register realms to be used for the simulation.
  */
  virtual void
  registerRealms() = 0;

  /*!
    @brief Register operators to be used for the simulation
  */
  virtual void
  registerOperators() = 0;

  /*!
    @brief Parse runtime options
    @details Override this routine if your time stepper can use run-time configuration of the solvers that it advances. This can e.g.
    be the CFL condition.
  */
  virtual void
  parseRuntimeOptions();

#ifdef CH_USE_HDF5
  /*!
    @brief Read header data from checkpoint file.
    @param[inout] a_header HDF5 header.
  */
  virtual void
  readCheckpointHeader(HDF5HeaderData& a_header);
#endif

#ifdef CH_USE_HDF5
  /*!
    @brief Write header data to checkpoint file.
    @param[inout] a_header HDF5 header.
  */
  virtual void
  writeCheckpointHeader(HDF5HeaderData& a_header) const;
#endif

#ifdef CH_USE_HDF5
  /*!  
    @brief Write checkpoint data to file
    @param[inout] a_handle HDF5 fil
    @param[in]    a_lvl    Grid level
    @details Implement this routine for checkpointing data for restarts. This will typically call the solvers' checkpointing routine. 
  */
  virtual void
  writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const = 0;
#endif

#ifdef CH_USE_HDF5
  /*!
    @brief Read checkpoint data from file
    @param[inout] a_handle HDF5 file
    @param[in]    a_lvl    Grid level
    @details Implement this routine for reading data for restarts. This will typically call the solvers' checkponiting routine. 
  */
  virtual void
  readCheckpointData(HDF5Handle& a_handle, const int a_lvl) = 0;
#endif

  /*!
    @brief Get the number of plot variables for this time stepper. 
    @details This is necessary because Driver, not TimeStepper, is responsible for allocating the necessary memory. 
    @return Returns number of plot variables that will be written during writePlotData
  */
  virtual int
  getNumberOfPlotVariables() const = 0;

  /*!
    @brief Get plot variable names
  */
  virtual Vector<std::string>
  getPlotVariableNames() const = 0;

  /*!
    @brief Write plot data to output holder. 
    @param[inout] a_output Output data holder.
    @param[inout] a_icomp  Starting component in a_output to begin at. 
    @param[in]    a_level  Grid level
  */
  virtual void
  writePlotData(LevelData<EBCellFAB>& a_output, int& a_icomp, const int a_level) const = 0;

  /*!
    @brief An option for calling special functions prior to plotting data. 
    Called by Driver in the IMMEDIATELY before writing the plot file. 
  */
  virtual void
  prePlot();

  /*!
    @brief Get computational loads to be checkpointed. 
    @details This is used by Driver both for setting up load-balanced restarts AND for plotting the computational loads to a file. This routine is
    disjoint from loadBalanceBoxes because this routine is not part of a regrid. This means that we are not operating with temporarily load balanced
    grids where the but the final ones. 
    @note The default implementation uses the box volume as a proxy for the load. You should overwrite this if you load balance your application, and also
    make sure that the loads returned from this routine are consistent with what you put in loadBalanceBoxes. 

    Also note that the return vector has the same ordering as the DisjointBoxLayout's boxes on the input grid level. See the implementation for further details. 
    @param[in] a_realm Realm
    @param[in] a_level Grid level
    @return Returns computational loads for each box on grid level a_level. 
    
  */
  virtual Vector<long int>
  getCheckpointLoads(const std::string a_realm, const int a_level) const;

  /*!
    @brief Compute a time step to be used by Driver. 
  */
  virtual Real
  computeDt() = 0;

  /*!
    @brief Advancement method. The implementation of this method should advance all equations of motion
    @param[in] a_dt Time step to be used for advancement
    @return    Returns the time step that was used. 
    @note The return value does not need to equal a_dt. Adaptive time stepping methods will generally return != a_dt.
  */
  virtual Real
  advance(const Real a_dt) = 0;

  /*!
    @brief Synchronzie solver times and time steps
    @param[in] a_step Time step
    @param[in] a_time Time (in seconds)
    @param[in] a_dt   Time step that was used. 
  */
  virtual void
  synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) = 0;

  /*!
    @brief Print a step report. 
    @details This is called by Driver after time step. The routine can be used to display use information about the simulation progress. 
  */
  virtual void
  printStepReport() = 0;

  /*!
    @brief Perform pre-regrid operations.
    @details This should include all copying all data which should be interpolated to the new grids. It can also include deallocating memory in case the regrid
    operation takes a lot of memory. 
    @param[in] a_lmin           The coarsest level that changes
    @param[in] a_oldFinestLevel The finest level before the regrid. 
  */
  virtual void
  preRegrid(const int a_lmin, const int a_oldFinestLevel) = 0;

  /*!
    @brief Time stepper regrid method. 
    @param[in] a_lmin           The coarsest level that changed. 
    @param[in] a_oldFinestLevel The finest level before the regrid. 
    @param[in] a_newFinestLevel The finest level after the regrid. 
  */
  virtual void
  regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) = 0;

  /*!
    @brief Perform post-regrid operations. 
    @details This includes all operations to be done AFTER interpolating data to new grids. 
  */
  virtual void
  postRegrid() = 0;

  /*!
    @brief Fuction which can have Driver do regrids at arbitrary points in the simulation. 
    @details This is called by Driver at every time step and if it returns true, Driver will regrid, regardless of whether or not it is 
    on the regular "regrid intervals". 
  */
  virtual bool
  needToRegrid();

  /*!
    @brief Load balancing query for a specified realm. If this returns true for a_realm, load balancing routines will be called during regrids. 
    @param[in] a_realm Realm name
  */
  virtual bool
  loadBalanceThisRealm(const std::string a_realm) const;

  /*!
    @brief Load balance grid boxes for a specified realm. 
    @param[out] a_procs       MPI ranks owning the various grid boxes. 
    @param[out] a_boxes       Grid boxes on every level (obtain them with a_grids[lvl].boxArray())
    @param[in]  a_realm       Realm identifier
    @param[in]  a_grids       Original grids
    @param[in]  a_lmin        Coarsest grid level that changed
    @param[in]  a_finestLevel New finest grid level
    @details This is only called by Driver if TimeStepper::loadBalanceThisRealm(a_realm) returned true. The default implementation
    uses volume-based loads for the grid patches. If the user wants to load balance boxes on a realm, this routine must be overwritten and
    he should compute loads for the various patches in a_grids and call LoadBalancing::makeBalance(Vector<int>, Vector<T>, Vector<Box>).
  */
  virtual void
  loadBalanceBoxes(Vector<Vector<int>>&             a_procs,
                   Vector<Vector<Box>>&             a_boxes,
                   const std::string                a_realm,
                   const Vector<DisjointBoxLayout>& a_grids,
                   const int                        a_lmin,
                   const int                        a_finestLevel);

protected:
  /*!
    @brief Class verbosity
  */
  int m_verbosity;

  /*!
    @brief Time step
  */
  int m_timeStep;

  /*!
    @brief TIme
  */
  Real m_time;

  /*!
    @brief Previous time step size
  */
  Real m_dt;

  /*!
    @brief AmrMesh. 
  */
  RefCountedPtr<AmrMesh> m_amr;

  /*!
    @brief Computational geometry. 
  */
  RefCountedPtr<ComputationalGeometry> m_computationalGeometry;
};

#include <CD_NamespaceFooter.H>

#endif
