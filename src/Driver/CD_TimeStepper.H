/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_TimeStepper.H
  @brief  Declaration of main (abstract) time stepper class. 
  @author Robert Marskar
*/

#ifndef CD_TimeStepper_H
#define CD_TimeStepper_H

// Our includes
#include <computational_geometry.H>
#include <mfis.H>
#include <CD_AmrMesh.H>
#include <CD_TimeCode.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Abstract class for performing time stepping
*/
class TimeStepper {
public:

  /*!
    @brief Default constructor -- does nothing
  */
  TimeStepper();

  /*!
    @brief Default destructor -- does nothing
  */
  virtual ~TimeStepper();

  // Set routines
  void setAmr(const RefCountedPtr<AmrMesh>& a_amr);
  void setComputationalGeometry(const RefCountedPtr<computational_geometry>& a_computationalGeometry);

  // Setup routines
  virtual void setup_solvers() = 0;
  virtual void allocate() = 0;
  virtual void initialData() = 0;
  virtual void postInitialize() = 0;
  virtual void postCheckpointSetup() = 0;
  virtual int getRedistributionRegSize() const;

  // Registration routines
  virtual void registerRealms() = 0;
  virtual void registerOperators() = 0;

  // Runtime options
  virtual void parseRuntimeOptions();

  // IO routines
  virtual void writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const = 0;
  virtual void readCheckpointData(HDF5Handle& a_handle, const int a_lvl) = 0;
  virtual void writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const = 0;
  virtual int  getNumberOfPlotVariables() const = 0;
  virtual Vector<long int> getCheckpointLoads(const std::string a_realm, const int a_level) const;

  // ADVANCE ROUTINES
  /*!
    @brief Compute a time step
    @param[out] a_dt       Time step to be used.
    @param[out] a_timeCode Identifier for how the time step was computed.
  */
  virtual void computeDt(Real& a_dt, TimeCode& a_timeCode) = 0;

  /*!
    @brief Advancement method. The implementation of this method should advance all equations of motion
    @param[in] a_dt Time step to be used for advancement
    @return    Returns the time step that was used. 
    @note The return value does not need to equal a_dt. Adaptive time stepping methods will generally return != a_dt.
  */
  virtual Real advance(const Real a_dt) = 0;

  /*!
    @brief Synchronzie solver times and time steps
    @param[in] a_step Time step
    @param[in] a_time Time (in seconds)
    @param[in] a_dt   Time step that was used. 
  */
  virtual void synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt) = 0;

  /*!
    @brief Print a step report. This is called by Driver at every time step and you can output any information you like.
  */
  virtual void print_stepReport() = 0;

  // REGRID ROUTINES
  /*!
    @brief Perform pre-regrid operations.
    @details This should include all copying all data which should be interpolated to the new grids. It can also include deallocating memory in case the regrid
    operation takes a lot of memory. 
    @param[in] a_lmin           The coarsest level that changes
    @param[in] a_oldFinestLevel The previous finest level. 
  */
  virtual void preRegrid(const int a_lmin, const int a_oldFinestLevel) = 0;

  /*!
    @brief Perform post-regrid operations. 
    @details This includes all operations to be done after interpolating data to new grids (such as updating electrostatics, for example). 
    @param[in] a_lmin           The coarsest level that changed. 
    @param[in] a_oldFinestLevel The previous finest level. 
    @param[in] a_newFinestLevel The new finest level. 
  */
  virtual void postRegrid() = 0;

  /*!
    @brief Time stepper regrid method. 
    @param[in] a_lmin           The coarsest level that changed. 
    @param[in] a_oldFinestLevel The previous finest level. 
    @param[in] a_newFinestLevel The new finest level. 
  */
  virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) = 0;

  /*!
    @brief Fuction which can have Driver do regrids at arbitrary points in the simulation. This is called at every time step and if it returns true, Driver will regrid. 
  */
  virtual bool needToRegrid();

  /*!
    @brief Load balancing query for a specific realm. If this returns true for a_realm, a_realm will be load balanced. 
    @param[in] a_realm Realm name
  */
  virtual bool loadBalanceThisRealm(const std::string a_realm) const;

  /*!
    @brief Load balance grid boxes for a specific realm. 
    @param[out] a_procs       MPI ranks owning the various grid boxes. 
    @param[out] a_boxes       Grid boxes on every level
    @param[in]  a_realm       Realm identifier
    @param[in]  a_grids       Original grids
    @param[in]  a_lmin        Coarsest grid level that changed
    @param[in]  a_finestLevel New finest grid level
    @details This is only called by Driver if TimeStepper::loadBalanceThisRealm(a_realm) returned true
  */
  virtual void loadBalanceBoxes(Vector<Vector<int> >&            a_procs,
				Vector<Vector<Box> >&            a_boxes,
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
  int m_step;
  
  /*!
    @brief TIme
  */
  Real m_time;

  /*!
    @brief Previous time step size
  */
  Real m_dt;

  /*!
    @brief Time code for
  */
  TimeCode m_timeCode;

  /*!
    @brief AMR
  */
  RefCountedPtr<AmrMesh> m_amr;

  /*!
    @brief Physical domain
  */
  RefCountedPtr<computational_geometry> m_computationalGeometry;
};

#include <CD_NamespaceFooter.H>

#endif
