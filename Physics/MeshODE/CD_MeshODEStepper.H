/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_MeshODEStepper.H
  @brief  Declaration of the Physics::MeshODE::MeshODEStepper TimeStepper
  @author Robert Marskar
*/

#ifndef CD_MESHODESTEPPER_H
#define CD_MESHODESTEPPER_H

// Std includes
#include <array>

// Our includes
#include <CD_TimeStepper.H>
#include <CD_MeshODESolver.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace MeshODE {

    /**
      @brief Time integration algorithms available for advancing the ODE.
    */
    enum class IntegrationAlgorithm
    {
      Euler, ///< Explicit forward Euler (first order).
      RK2,   ///< Explicit second-order Runge-Kutta (Heun's method).
      RK4    ///< Explicit classical fourth-order Runge-Kutta.
    };

    /**
      @brief TimeStepper for advancing a system of N coupled ODEs on an AMR mesh.
      @details At each time step the right-hand side @f$ \mathbf{f}(\mathbf{y}, t) @f$ is evaluated
      point-wise on the mesh and the ODE @f$ \dot{\mathbf{y}} = \mathbf{f}(\mathbf{y},t) @f$ is
      integrated with the chosen algorithm (Euler, RK2, or RK4). The template parameter @p N sets
      the number of solution components stored per cell.

      The RHS function and initial data function must be set by the derived application before
      Driver calls setupSolvers().
    */
    template <size_t N>
    class MeshODEStepper : public TimeStepper
    {
    public:
      /**
	@brief Constructor. Does nothing. 
      */
      MeshODEStepper();

      /**
	@brief Copy constructor. Disallowed.
      */
      MeshODEStepper(const MeshODEStepper<N>&) = delete;

      /**
	@brief Move constructor. Disallowed.
      */
      MeshODEStepper(const MeshODEStepper<N>&&) = delete;

      /**
	@brief Destructor. 
      */
      virtual ~MeshODEStepper();

      /**
	@brief Copy assignment operator. Disallowed.
      */
      MeshODEStepper&
      operator=(const MeshODEStepper<N>&) = delete;

      /**
	@brief Move assignment operator. Disallowed.
      */
      MeshODEStepper&
      operator=(const MeshODEStepper<N>&&) = delete;

      /**
	@brief Instantiate the ODE solver. 
      */
      void
      setupSolvers() override;

      /**
	@brief Allocate storage for solvers and time stepper. 
      */
      void
      allocate() override;

      /**
	@brief Fill problem with initial data. 
      */
      void
      initialData() override;

      /**
	@brief Perform any post-initialization steps. 
      */
      void
      postInitialize() override;

      /**
	@brief Post checkpoint operations. 
      */
      void
      postCheckpointSetup() override;

      /**
	@brief Register realms. Primal is the only realm we need. 
      */
      void
      registerRealms() override;

      /**
	@brief Register operators. 
	@details Only need the solver's operators.
      */
      void
      registerOperators() override;

      /**
	@brief Parse options
      */
      void
      parseOptions();

      /**
	@brief Parse runtime options
      */
      void
      parseRuntimeOptions() override;

#ifdef CH_USE_HDF5
      /**  
	@brief Write checkpoint data to file
	@param[inout] a_handle HDF5 file
	@param[in]    a_lvl    Grid level
      */
      virtual void
      writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
#endif

#ifdef CH_USE_HDF5
      /**
	@brief Read checkpoint data from file
	@param[inout] a_handle HDF5 file
	@param[in]    a_lvl    Grid level
      */
      virtual void
      readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
#endif
      /**
	@brief Get the number of plot variables for this time stepper. 
	@details This is necessary because Driver, not TimeStepper, is responsible for allocating the necessary memory. 
	@return Returns number of plot variables that will be written during writePlotData
      */
      virtual int
      getNumberOfPlotVariables() const override;

      /**
	@brief Get plot variable names.
	@return Vector of variable name strings contributed by this TimeStepper.
      */
      virtual Vector<std::string>
      getPlotVariableNames() const override;

      /**
	@brief Write plot data to output holder.
	@param[inout] a_output       Output data holder.
	@param[inout] a_icomp        Starting component in a_output; incremented on exit.
	@param[in]    a_outputRealm  Realm where a_output lives.
	@param[in]    a_level        Grid level.
      */
      virtual void
      writePlotData(LevelData<EBCellFAB>& a_output,
                    int&                  a_icomp,
                    const std::string     a_outputRealm,
                    const int             a_level) const override;

      /**
	@brief Compute a time step.
	@return The fixed time step m_dt read from input options.
      */
      virtual Real
      computeDt() override;

      /**
	@brief Advancement method. Swaps between various kernels. 
	@param[in] a_dt Time step to be used for advancement
	@return    Returns the time step that was used (equal to a_dt in this case). 
      */
      virtual Real
      advance(const Real a_dt) override;

      /**
	@brief Synchronize solver times and time steps.
	@param[in] a_step Current time step number.
	@param[in] a_time Current simulation time in seconds.
	@param[in] a_dt   Time step that was used.
      */
      virtual void
      synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;

      /**
	@brief Print a step report. 
	@details This is called by Driver after time step. The routine can be used to display use information about the simulation progress. Not used by this class. 
      */
      virtual void
      printStepReport() override
      {}

      /**
	@brief Perform pre-regrid operations.
	@param[in] a_lmin           The coarsest level that changes
	@param[in] a_oldFinestLevel The finest level before the regrid. 
      */
      virtual void
      preRegrid(const int a_lmin, const int a_oldFinestLevel) override;

      /**
	@brief Time stepper regrid method. 
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      virtual void
      regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /**
	@brief Perform post-regrid operations. 
	@details This includes all operations to be done AFTER interpolating data to new grids. 
      */
      virtual void
      postRegrid() override;

    protected:
      /**
	@brief Integration algorithm
      */
      IntegrationAlgorithm m_algorithm;

      /**
	@brief Time step to use
      */
      Real m_dt;

      /**
	@brief Realm where solver and m_velocity lives.
      */
      std::string m_realm;

      /**
	@brief Phase (gas/solid) where the solver lives.
      */
      phase::which_phase m_phase;

      /**
	@brief MeshODE solver holding the solution on the AMR hierarchy.
      */
      RefCountedPtr<MeshODESolver<N>> m_solver;

      /**
	@brief Right-hand side function @f$ \mathbf{f}(\mathbf{y}, t) @f$ evaluated point-wise on the mesh.
      */
      std::function<std::array<Real, N>(const std::array<Real, N> y, Real t)> m_rhsFunction;

      /**
	@brief Callable that returns the initial state @f$ \mathbf{y}_0(\mathbf{x}) @f$ at a given position.
      */
      std::function<std::array<Real, N>(const RealVect& a_position)> m_initialData;

      /**
	@brief Parse integrator
      */
      virtual void
      parseIntegrator();

      /**
	@brief Parse chattiness
      */
      virtual void
      parseVerbosity();

      /**
	@brief Parse the problem type. 
      */
      virtual void
      parseProblem();

      /**
	@brief Advance using the explicit Euler rule
	@param[in] a_dt Advanced time step
      */
      virtual void
      advanceEuler(const Real a_dt);

      /**
	@brief Advance using a second-order Runge-Kutta method. 
	@param[in] a_dt Advanced time step
      */
      virtual void
      advanceRK2(const Real a_dt);

      /**
	@brief Advance using a fourth order Runge-Kutta method. 
	@param[in] a_dt Advanced time step
      */
      virtual void
      advanceRK4(const Real a_dt);
    };
  } // namespace MeshODE
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_MeshODEStepperImplem.H>

#endif
