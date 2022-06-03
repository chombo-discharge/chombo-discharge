/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MeshODEStepper.H
  @brief  Declaration of a TimeStepper implementation for advancing an ODE on a mesh. 
  @author Robert Marskar
*/

#ifndef CD_MeshODEStepper_H
#define CD_MeshODEStepper_H

// Our includes
#include <CD_TimeStepper.H>
#include <CD_MeshODESolver.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace MeshODE {

    /*!
      @brief Various supported integration algorithms. 
    */
    enum class IntegrationAlgorithm
    {
      Euler,
      RK2,
      RK4
    };

    /*!
      @brief Implementation of TimeStepper for solving an ODE on a mesh. N is the number of variables. 
    */
    template <size_t N>
    class MeshODEStepper : public TimeStepper
    {
    public:
      /*!
	@brief Constructor. Does nothing. 
      */
      MeshODEStepper();

      /*!
	@brief Copy constructor. Disallowed.
      */
      MeshODEStepper(const MeshODEStepper<N>&) = delete;

      /*!
	@brief Move constructor. Disallowed.
      */
      MeshODEStepper(const MeshODEStepper<N>&&) = delete;

      /*!
	@brief Destructor. 
      */
      virtual ~MeshODEStepper();

      /*!
	@brief Copy assignment operator. Disallowed.
      */
      MeshODEStepper&
      operator=(const MeshODEStepper<N>&) = delete;

      /*!
	@brief Move assignment operator. Disallowed.
      */
      MeshODEStepper&
      operator=(const MeshODEStepper<N>&&) = delete;

      /*!
	@brief Instantiate the ODE solver. 
      */
      void
      setupSolvers() override;

      /*!
	@brief Allocate storage for solvers and time stepper. 
      */
      void
      allocate() override;

      /*!
	@brief Fill problem with initial data. 
      */
      void
      initialData() override;

      /*!
	@brief Perform any post-initialization steps. 
      */
      void
      postInitialize() override;

      /*!
	@brief Post checkpoint operations. 
      */
      void
      postCheckpointSetup() override;

      /*!
	@brief Register realms. Primal is the only realm we need. 
      */
      void
      registerRealms() override;

      /*!
	@brief Register operators. 
	@details Only need the solver's operators.
      */
      void
      registerOperators() override;

      /*!
	@brief Parse options
      */
      void
      parseOptions();

      /*!
	@brief Parse runtime options
      */
      void
      parseRuntimeOptions() override;

#ifdef CH_USE_HDF5
      /*!  
	@brief Write checkpoint data to file
	@param[inout] a_handle HDF5 file
	@param[in]    a_lvl    Grid level
      */
      virtual void
      writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
#endif

#ifdef CH_USE_HDF5
      /*!
	@brief Read checkpoint data from file
	@param[inout] a_handle HDF5 file
	@param[in]    a_lvl    Grid level
      */
      virtual void
      readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
#endif

      /*!
	@brief Write plot data to output holder. 
	@param[inout] a_output            Output data holder.
	@param[in]    a_plotVariableNames Plot variable names.
	@param[inout] a_icomp             Starting component in a_output to begin at. 
      */
      virtual void
      writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const override;

      /*!
	@brief Get the number of plot variables for this time stepper. 
	@details This is necessary because Driver, not TimeStepper, is responsible for allocating the necessary memory. 
	@return Returns number of plot variables that will be written during writePlotData
      */
      virtual int
      getNumberOfPlotVariables() const override;

      /*!
	@brief Compute a time step to be used by Driver. 
      */
      virtual Real
      computeDt() override;

      /*!
	@brief Advancement method. Swaps between various kernels. 
	@param[in] a_dt Time step to be used for advancement
	@return    Returns the time step that was used (equal to a_dt in this case). 
      */
      virtual Real
      advance(const Real a_dt) override;

      /*!
	@brief Synchronize solver times and time steps
	@param[in] a_step Time step
	@param[in] a_time Time (in seconds)
	@param[in] a_dt   Time step that was used. 
      */
      virtual void
      synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;

      /*!
	@brief Print a step report. 
	@details This is called by Driver after time step. The routine can be used to display use information about the simulation progress. Not used by this class. 
      */
      virtual void
      printStepReport() override
      {}

      /*!
	@brief Perform pre-regrid operations.
	@param[in] a_lmin           The coarsest level that changes
	@param[in] a_oldFinestLevel The finest level before the regrid. 
      */
      virtual void
      preRegrid(const int a_lmin, const int a_oldFinestLevel) override;

      /*!
	@brief Time stepper regrid method. 
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      virtual void
      regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /*!
	@brief Perform post-regrid operations. 
	@details This includes all operations to be done AFTER interpolating data to new grids. 
      */
      virtual void
      postRegrid() override;

    protected:
      /*!
	@brief Integration algorithm
      */
      IntegrationAlgorithm m_algorithm;

      /*!
	@brief Time step to use
      */
      Real m_dt;

      /*!
	@brief Realm where solver and m_velocity lives.
      */
      std::string m_realm;

      /*!
	@brief Phase (gas/solid) where the solver lives.
      */
      phase::which_phase m_phase;

      /*!
	@brief Tracer particle solvers
      */
      RefCountedPtr<MeshODESolver<N>> m_solver;

      /*!
	@brief Encapsulation of f = f(y,t)
      */
      std::function<std::array<Real, N>(const std::array<Real, N> y, Real t)> m_rhsFunction;

      /*!
	@brief Initial data function. 
      */
      std::function<std::array<Real, N>(const RealVect& a_position)> m_initialData;

      /*!
	@brief Parse integrator
      */
      virtual void
      parseIntegrator();

      /*!
	@brief Parse chattiness
      */
      virtual void
      parseVerbosity();

      /*!
	@brief Parse the problem type. 
      */
      virtual void
      parseProblem();

      /*!
	@brief Advance using the explicit Euler rule
	@param[in] a_dt Advanced time step
      */
      virtual void
      advanceEuler(const Real a_dt);

      /*!
	@brief Advance using a second-order Runge-Kutta method. 
	@param[in] a_dt Advanced time step
      */
      virtual void
      advanceRK2(const Real a_dt);

      /*!
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
