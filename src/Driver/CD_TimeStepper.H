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
#include <CD_NamespaceHeader.H>
  
enum class TimeCode {
  Advection,
  Diffusion,
  AdvectionDiffusion,  
  Source,
  RelaxationTime,
  Restricted,
  Hardcap,
  Error,
  Physics
};

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
  virtual Vector<long int> getCheckpointLoads(const std::string a_Realm, const int a_level) const;

  // Advance routines
  virtual void computeDt(Real& a_dt, TimeCode& a_timeCode) = 0;
  virtual Real advance(const Real a_dt) = 0;
  virtual void synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt) = 0;
  virtual void print_stepReport() = 0;

  // Regrid routines
  virtual void deallocate() = 0;
  virtual void preRegrid(const int a_lmin, const int a_oldFinestLevel) = 0;
  virtual void postRegrid() = 0;
  virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) = 0;
  virtual bool needToRegrid();
  virtual bool loadBalanceThisRealm(const std::string a_Realm) const;
  virtual void loadBalanceBoxes(Vector<Vector<int> >&            a_procs,
				Vector<Vector<Box> >&            a_boxes,
				const std::string                a_Realm,
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
