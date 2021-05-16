/*!
  @file   geometry_stepper.H
  @brief  Declaration of a TimeStepper class without any solver functionality, i.e. for geometries only. 
  @author Robert Marskar
  @date   March 2020
*/

#ifndef _GEOMETRY_STEPPER_
#define _GEOMETRY_STEPPER_

#include "CD_TimeStepper.H"

#include "CD_NamespaceHeader.H"
namespace physics {
  namespace geometry {
    
    /*!
      @brief Abstract class for performing time stepping
    */
    class geometry_stepper: public TimeStepper {
    public:

      geometry_stepper();
      ~geometry_stepper();

      // Setup routines
      virtual void setup_solvers() override;
      virtual void allocate() override;
      virtual void initialData() override;
      virtual void postInitialize() override;
      virtual void postCheckpointSetup() override;

      // Registration routines
      virtual void registerRealms() override;
      virtual void registerOperators() override;

      // IO routines
      virtual void writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
      virtual void readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
      virtual void writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const override;
      virtual int  getNumberOfPlotVariables() const override;

      // Advance routines
      virtual void computeDt(Real& a_dt, TimeCode& a_timeCode) override;
      virtual Real advance(const Real a_dt) override;
      virtual void synchronize_solver_times(const int a_step, const Real a_time, const Real a_dt) override;
      virtual void print_stepReport() override;

      // Regrid routines
      virtual void deallocate() override;
      virtual void preRegrid(const int a_lmin, const int a_oldFinestLevel) override;
      virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;
      virtual void postRegrid() override;

    };
  }
}
#include "CD_NamespaceFooter.H"

#endif
