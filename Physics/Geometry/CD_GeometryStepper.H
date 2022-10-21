/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
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

  /*!
    @brief Namespace for encapsulating the Geometry module -- i.e. functionality for just plotting various geometries. 
  */
  namespace Geometry {

    /*!
      @brief This implementation does absolutely nothing -- it instantiates a TimeStepper with empty functionality. 
      @details The usage of this module is just to make it easy to develop/implement geometries. The user will simply set up an "application" with this module
      and use Driver.geometry_only=true when running his application. 
    */
    class GeometryStepper : public TimeStepper
    {
    public:
      /*!
	@brief Destructor (does nothing)
      */
      GeometryStepper();

      /*!
	@brief Destructor (does nothing)
      */
      ~GeometryStepper();

      /*!
	@brief Set up solvers (does nothing)
      */
      virtual void
      setupSolvers() override;

      /*!
	@brief Allocate solver internals (does nothing)
      */
      virtual void
      allocate() override;

      /*!
	@brief Fill solver initial data (does nothing)
      */
      virtual void
      initialData() override;

      /*!
	@brief Post initialization function (does nothing)
      */
      virtual void
      postInitialize() override;

      /*!
	@brief Post-restart function (does nothing)
      */
      virtual void
      postCheckpointSetup() override;

      /*!
	@brief Realm register -- does nothing
      */
      virtual void
      registerRealms() override;

      /*!
	@brief operator registration (does nothing)
      */
      virtual void
      registerOperators() override;

#ifdef CH_USE_HDF5
      /*!  
	@brief Write checkpoint data to file (does nothing)
	@param[inout] a_handle HDF5 fil
	@param[in]    a_lvl    Grid level
      */
      virtual void
      writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
#endif
#ifdef CH_USE_HDF5
      /*!
	@brief Read checkpoint data from file (does nothing)
	@param[inout] a_handle HDF5 file
	@param[in]    a_lvl    Grid level
      */
      virtual void
      readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
#endif

      /*!
	@brief Write plot data to file (does nothing)
	@param[inout] a_output            Output data holder.
	@param[in]    a_plotVariableNames Plot variable names.
	@param[inout] a_icomp             Starting component in a_output to begin at. 
      */
      virtual void
      writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const override;

      /*!
	@brief Get the number of plot variables for this time stepper. 
	@return Returns 0
      */
      virtual int
      getNumberOfPlotVariables() const override;

      /*!
	@brief Compute a time step to be used by Driver. 
	@return Returns 0
      */
      virtual Real
      computeDt() override;

      /*!
	@brief Advancement routine. Does nothing. 
	@param[in] a_dt Time step to be used for advancement
	@return    Returns the time step that was used. 
	@note This returns the largest time step possible to ensure end-of-simulation.
      */
      virtual Real
      advance(const Real a_dt) override;

      /*!
	@brief Synchronzie solver times and time steps (does nothing)
	@param[in] a_step Time step
	@param[in] a_time Time (in seconds)
	@param[in] a_dt   Time step that was used. 
      */
      virtual void
      synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;

      /*!
	@brief Print a step report (does nothing)
	@details This is called by Driver after time step. The routine can be used to display use information about the simulation progress. 
      */
      virtual void
      printStepReport() override;

      /*!
	@brief Perform pre-regrid operations (does nothing)
	@param[in] a_lmin           The coarsest level that changes
	@param[in] a_oldFinestLevel The finest level before the regrid. 
      */
      virtual void
      preRegrid(const int a_lmin, const int a_oldFinestLevel) override;

      /*!
	@brief Time stepper regrid method (does nothing)
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      virtual void
      regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /*!
	@brief Perform post-regrid operations (does nothing)
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      virtual void
      postRegrid() override;
    };
  } // namespace Geometry
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
