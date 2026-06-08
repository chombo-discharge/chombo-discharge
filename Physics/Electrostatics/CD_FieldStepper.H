/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_FieldStepper.H
  @brief  Declaration of the Physics::Electrostatics::FieldStepper TimeStepper
  @author Robert Marskar
*/

#ifndef CD_FIELDSTEPPER_H
#define CD_FIELDSTEPPER_H

// Our includes
#include <CD_TimeStepper.H>
#include <CD_SurfaceODESolver.H>
#include <CD_FieldSolver.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace Electrostatics {

    /**
      @brief TimeStepper for solving the electrostatic Poisson equation, optionally with surface charge.
      @details This module solves the electrostatic problem on a static or adaptively refined AMR hierarchy.
      The template parameter @p T is the concrete FieldSolver type (e.g. FieldSolverGMG).

      The solver decomposes the potential into contributions from the applied voltage and from space/surface
      charges. On live electrodes the voltage defaults to 1 V; use setRho() and setSigma() before calling
      postInitialize() to override the default zero charges.

      Adaptive mesh refinement is supported when a CellTagger is attached via Driver — the stepper
      re-solves the Poisson equation after every regrid in postRegrid().
    */
    template <class T>
    class FieldStepper : public TimeStepper
    {
    public:
      /**
	@brief Constructor -- parses some input options
      */
      FieldStepper();

      /**
	@brief Destructor (does nothing)
      */
      virtual ~FieldStepper();

      /**
	@brief Solver setup routine. This instantiates the FieldSolver and parses input options
      */
      void
      setupSolvers() override;

      /**
	@brief Set initial data -- this sets the space and surface charges
      */
      void
      initialData() override;

      /**
	@brief Post-initialization routine. This solves the Poisson equation. 
      */
      void
      postInitialize() override;

      /**
	@brief Allocation method -- allocates memory and internal data for solvers
      */
      void
      allocate() override;

      /**
	@brief Get the electrostatic potential.
	@return Reference to the multifluid AMR cell data holding the potential.
      */
      MFAMRCellData&
      getPotential();

#ifdef CH_USE_HDF5
      /**
	@brief Write checkpoint data to file.
	@note This module does not support checkpoint-restart; calling this is an error.
	@param[inout] a_handle HDF5 file handle.
	@param[in]    a_lvl    Grid level.
      */
      void
      writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
#endif

#ifdef CH_USE_HDF5
      /**
	@brief Read checkpoint data from file -- this physics does not checkpoint data.
	@param[inout] a_handle HDF5 file
	@param[in]    a_lvl    Grid level
      */
      void
      readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
#endif

      /**
	@brief Post-initialize operations to be performed -- does nothing. 
      */
      void
      postCheckpointSetup() override
      {}

      /**
	@brief Get number of plot variables contributed by this TimeStepper.
	@return Number of plot variables written by writePlotData().
      */
      int
      getNumberOfPlotVariables() const override;

      /**
	@brief Get plot variable names.
	@return Vector of plot variable name strings.
      */
      Vector<std::string>
      getPlotVariableNames() const override;

      /**
	@brief Write plot data to output holder.
	@param[inout] a_output       Output data holder.
	@param[inout] a_icomp        Starting component in a_output; incremented by the number of components written.
	@param[in]    a_outputRealm  Realm where a_output lives.
	@param[in]    a_level        Grid level.
      */
      void
      writePlotData(LevelData<EBCellFAB>& a_output,
                    int&                  a_icomp,
                    const std::string     a_outputRealm,
                    const int             a_level) const override;

      /**
	@brief Compute time step
	@details This is a static solver so it returns 0.0
      */
      Real
      computeDt() override
      {
        return 0.0;
      }

      /**
	@brief Perform a single time step with step a_dt.
	@param[in] a_dt Time step
	@return a_dt
      */
      Real
      advance(const Real a_dt) override;

      /**
	@brief Synchronize solver times and time steps.
	@param[in] a_step Current time step number.
	@param[in] a_time Current simulation time in seconds.
	@param[in] a_dt   Time step that was used.
      */
      void
      synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;

      /**
	@brief Print step report -- this does nothing. 
      */
      void
      printStepReport() override {};

      /**
	@brief Register simulation realms to be used for this simulation module. 
      */
      void
      registerRealms() override;

      /**
	@brief Register operators for this simulation module
      */
      void
      registerOperators() override;

      /**
	@brief Perform pre-regrid operations.
	@param[in] a_lbase        The coarsest level that changes.
	@param[in] a_finestLevel  The finest level before the regrid.
      */
      void
      preRegrid(const int a_lbase, const int a_finestLevel) override;

      /**
	@brief Load balancing query for a specified realm. If this returns true for a_realm, load balancing routines will be called during regrids. 
	@param[in] a_realm Realm name
      */
      bool
      loadBalanceThisRealm(const std::string a_realm) const override;

      /**
	@brief Load balance grid boxes for a specified realm. 
	@param[out] a_procs       MPI ranks owning the various grid boxes. 
	@param[out] a_boxes       Grid boxes on every level (obtain them with a_grids[lvl].boxArray())
	@param[in]  a_realm       Realm identifier
	@param[in]  a_grids       Original grids
	@param[in]  a_lmin        Coarsest grid level that changed
	@param[in]  a_finestLevel New finest grid level
	@details This is only called by Driver if TimeStepper::loadBalanceThisRealm(a_realm) returned true. The default implementation
	uses volume-based loads for the grid patches. If the user wants to load balance boxes on a realm, this routine must be overwritten and
	he should compute loads for the various patches in a_grids and call LoadBalancing::makeBalance on each level. It is up to the user/programmer
	to decide if load balancing should be done independently on each level, or if loads per MPI rank are accumulated across levels.
      */
      void
      loadBalanceBoxes(Vector<Vector<int>>&             a_procs,
                       Vector<Vector<Box>>&             a_boxes,
                       const std::string                a_realm,
                       const Vector<DisjointBoxLayout>& a_grids,
                       const int                        a_lmin,
                       const int                        a_finestLevel) override;

      /**
	@brief Regrid method -- regrids the potential distribution in FieldSolver (but does not solve the Poisson equation, that is done in postRegrid()).
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      void
      regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /**
	@brief Perform post-regrid operations -- this will resolve the Poisson equation. 
	@details This includes all operations to be done AFTER interpolating data to new grids. 
      */
      void
      postRegrid() override;

      /**
        @brief Set the space charge distribution.
        @param[in] a_rho   Callable returning the charge density (C/m^3) at a given position.
        @param[in] a_phase Phase (gas or solid) in which to apply the space charge.
      */
      void
      setRho(const std::function<Real(const RealVect& a_pos)>& a_rho, const phase::which_phase a_phase) noexcept;

      /**
        @brief Set the surface charge distribution.
        @param[in] a_sigma Callable returning the surface charge density (C/m^2) at a given position.
      */
      void
      setSigma(const std::function<Real(const RealVect& a_pos)>& a_sigma) noexcept;

    protected:
      /**
	@brief If true, the simulation will be introspectively load balanced. 
      */
      bool m_loadBalance;

      /**
	@brief If using load balancing, this specifies how boxes will be sorted (in space) before load balancing. 
      */
      BoxSorting m_boxSort;

      /**
	@brief Reference to the FieldSolver
      */
      RefCountedPtr<FieldSolver> m_fieldSolver;

      /**
	@brief Reference to the surface charge solver.
      */
      RefCountedPtr<SurfaceODESolver<1>> m_sigma;

      /**
	@brief Realm where the solver lives
      */
      std::string m_realm;

      /**
	@brief Space charge inside gas phase.
      */
      std::function<Real(const RealVect& a_pos)> m_rhoGas;

      /**
	@brief Space charge inside dielectric.
      */
      std::function<Real(const RealVect& a_pos)> m_rhoDielectric;

      /**
	@brief Surface charge distribution
      */
      std::function<Real(const RealVect& a_pos)> m_surfaceChargeDensity;

      /**
	@brief Solve the Poisson equation and compute the electric field.
	@details Called from postInitialize() and postRegrid().
      */
      void
      solvePoisson();
    };
  } // namespace Electrostatics
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_FieldStepperImplem.H>

#endif
