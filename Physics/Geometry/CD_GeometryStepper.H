/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_GeometryStepper.H
  @brief  Declaration of the Physics::Geometry::GeometryStepper TimeStepper
  @author Robert Marskar
*/

#ifndef CD_GEOMETRYSTEPPER_H
#define CD_GEOMETRYSTEPPER_H

// Our includes
#include <CD_TimeStepper.H>
#include <CD_NamespaceHeader.H>

namespace Physics {

  /**
    @brief Namespace for encapsulating the Geometry module — functionality for visualising geometries without any physics.
  */
  namespace Geometry {

    /**
      @brief Minimal TimeStepper that performs no physics, intended for geometry development and visualisation.
      @details All lifecycle methods are no-ops. The intended use is to set `Driver.geometry_only=true`
      in the input script, which causes Driver to write a single plot file of the computational geometry
      and then exit. This makes it easy to inspect and tune complex EB geometries without building a
      full physics application.
    */
    class GeometryStepper : public TimeStepper
    {
    public:
      /**
        @brief Constructor.
      */
      GeometryStepper();

      /**
        @brief Destructor.
      */
      virtual ~GeometryStepper();

      /**
        @brief Set up solvers (no-op).
      */
      virtual void
      setupSolvers() override;

      /**
        @brief Allocate solver internals (no-op).
      */
      virtual void
      allocate() override;

      /**
        @brief Fill solver initial data (no-op).
      */
      virtual void
      initialData() override;

      /**
        @brief Post-initialization hook (no-op).
      */
      virtual void
      postInitialize() override;

      /**
        @brief Post-restart hook (no-op).
      */
      virtual void
      postCheckpointSetup() override;

      /**
        @brief Register realms (no-op).
      */
      virtual void
      registerRealms() override;

      /**
        @brief Register operators (no-op).
      */
      virtual void
      registerOperators() override;

#ifdef CH_USE_HDF5
      /**
        @brief Write checkpoint data to file (no-op).
        @param[in,out] a_handle HDF5 file handle.
        @param[in]    a_lvl    Grid level.
      */
      virtual void
      writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
#endif

#ifdef CH_USE_HDF5
      /**
        @brief Read checkpoint data from file (no-op).
        @param[in,out] a_handle HDF5 file handle.
        @param[in]    a_lvl    Grid level.
      */
      virtual void
      readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
#endif

      /**
        @brief Get the number of plot variables.
        @return Always returns 0.
      */
      virtual int
      getNumberOfPlotVariables() const override;

      /**
        @brief Get plot variable names.
        @return Always returns an empty vector.
      */
      virtual Vector<std::string>
      getPlotVariableNames() const override;

      /**
        @brief Write plot data to output holder (no-op).
        @param[in,out] a_output       Output data holder.
        @param[in,out] a_icomp        Starting component in a_output.
        @param[in]    a_realm        Output realm name.
        @param[in]    a_level        Grid level.
      */
      virtual void
      writePlotData(LevelData<EBCellFAB>& a_output,
                    int&                  a_icomp,
                    const std::string     a_realm,
                    const int             a_level) const override;

      /**
        @brief Compute a time step.
        @return Always returns 0.
      */
      virtual Real
      computeDt() override;

      /**
        @brief Advance by one time step (no-op).
        @param[in] a_dt Time step size.
        @return Returns std::numeric_limits<Real>::max() to drive Driver to terminate after one step.
      */
      virtual Real
      advance(const Real a_dt) override;

      /**
        @brief Synchronize solver times and time steps (no-op).
        @param[in] a_step Current time step number.
        @param[in] a_time Current simulation time in seconds.
        @param[in] a_dt   Time step that was used.
      */
      virtual void
      synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;

      /**
        @brief Print a step report (no-op).
      */
      virtual void
      printStepReport() override;

      /**
        @brief Perform pre-regrid operations (no-op).
        @param[in] a_lmin           The coarsest level that changes.
        @param[in] a_oldFinestLevel The finest level before the regrid.
      */
      virtual void
      preRegrid(const int a_lmin, const int a_oldFinestLevel) override;

      /**
        @brief Regrid method (no-op).
        @param[in] a_lmin           The coarsest level that changed.
        @param[in] a_oldFinestLevel The finest level before the regrid.
        @param[in] a_newFinestLevel The finest level after the regrid.
      */
      virtual void
      regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /**
        @brief Post-regrid operations (no-op).
      */
      virtual void
      postRegrid() override;
    };
  } // namespace Geometry
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
