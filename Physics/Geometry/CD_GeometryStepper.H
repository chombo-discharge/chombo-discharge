/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_GeometryStepper.H
  @brief  Declaration of a TimeStepper class without any solver functionality, i.e. for geometries only. 
  @author Robert Marskar
*/

#ifndef _GEOMETRY_STEPPER_
#define _GEOMETRY_STEPPER_

// Our includes
#include <CD_TimeStepper.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace Geometry {
    
    /*!
      @brief Abstract class for performing time stepping
    */
    class GeometryStepper: public TimeStepper {
    public:

      GeometryStepper();
      ~GeometryStepper();

      // Setup routines
      virtual void setupSolvers() override;
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
      virtual void synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;
      virtual void printStepReport() override;

      // Regrid routines
      virtual void preRegrid(const int a_lmin, const int a_oldFinestLevel) override;
      virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;
      virtual void postRegrid() override;

    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
