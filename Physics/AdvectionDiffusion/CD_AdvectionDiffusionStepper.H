/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_AdvectionDiffusionStepper.H
  @brief  TimeStepper class for only solving the the advection-diffusion equation
  @author Robert Marskar
*/

#ifndef CD_AdvectionDiffusionStepper
#define CD_AdvectionDiffusionStepper

// Our includes
#include <CD_TimeStepper.H>
#include <CD_CdrSolver.H>
#include <CD_AdvectionDiffusionSpecies.H>
#include <CD_NamespaceHeader.H>

namespace Physics {

  /*!
    @brief Physics functionality for running scalar advection-diffusion code
  */
  namespace AdvectionDiffusion {

    /*!
      @brief Implementation of TimeStepper which runs a single advection-diffusion-reaction solver.
    */
    class AdvectionDiffusionStepper : public TimeStepper
    {
    public:
      /*!
	@brief Default constructor -- only parses options. 
      */
      AdvectionDiffusionStepper();

      /*!
	@brief Full constructor
	@param[in] a_solver CDR solver. 
      */
      AdvectionDiffusionStepper(RefCountedPtr<CdrSolver>& a_solver);

      /*!
	@brief Destructor (does nothing)
      */
      virtual ~AdvectionDiffusionStepper();

      /*!
	@brief Fill solver with initial data
      */
      void
      initialData() override;

      /*!
	@brief Allocate memory
      */
      void
      allocate() override;

      /*!
	@brief Perform post-initialization routines (there are none)
      */
      void
      postInitialize() override
      {}

      // IO routines
#ifdef CH_USE_HDF5
      /*!  
	@brief Write checkpoint data to file
	@param[inout] a_handle HDF5 fil
	@param[in]    a_lvl    Grid level
      */
      void
      writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
#endif

#ifdef CH_USE_HDF5
      /*!
	@brief Read checkpoint data from file
	@param[inout] a_handle HDF5 file
	@param[in]    a_lvl    Grid level
      */
      void
      readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
#endif

      /*!
	@brief Perform post-checkpoint setup routines
	@details This will set source terms, velocities, and diffusion coefficients. 
      */
      void
      postCheckpointSetup() override;

      /*!
	@brief Get the number of plot variables
      */
      int
      getNumberOfPlotVariables() const override;

      /*!
	@brief Get plot variable name
      */
      Vector<std::string>
      getPlotVariableNames() const override;

      /*!
	@brief Write plot data to output holder. 
	@param[inout] a_output            Output data holder.
	@param[inout] a_icomp             Starting component in a_output to begin at. 
	@param[in]    a_outputRealm       Realm where a_output belongs
	@param[in]    a_level             Grid level
      */
      void
      writePlotData(LevelData<EBCellFAB>& a_output, int& a_icomp, const std::string a_outputRealm, const int a_level) const override;

      /*!
	@brief Compute a time step to be used by Driver. 
      */
      Real
      computeDt() override;

      /*!
	@brief Synchronize solver times
	@param[in] a_step Time step
	@param[in] a_time Time (in seconds)
	@param[in] a_dt   Time step that was used. 
      */
      void
      synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;

      /*!
	@brief Print a step report if you want
	@note Implementation does nothing. 
      */
      void
      printStepReport() override
      {}

      /*!
	@brief Parse runtime options
	@details Reads in verbosity level, integrator, and CFL number to use. 
      */
      void
      parseRuntimeOptions() override;

      /*!
	@brief Register all realms that will be used
      */
      void
      registerRealms() override;

      /*!
	@brief Register all operators that will be used.
      */
      void
      registerOperators() override;

      /*!
	@brief Perform pre-regrid operations. Calls CdrSolver::preRegrid
	@param[in] a_lmin           The coarsest level that changes
	@param[in] a_oldFinestLevel The finest level before the regrid. 
      */
      void
      preRegrid(const int a_lbase, const int a_finestLevel) override;

      /*!
	@brief Instantiate the solver. 
      */
      void
      setupSolvers() override;

      /*!
	@brief Advancement function -- advances one time step.
	@param[in] a_dt Time step to be used for advancement
	@return a_dt
      */
      Real
      advance(const Real a_dt) override;

      /*!
	@brief Regrid method -- calls regrid method for solver.
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      void
      regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /*!
	@brief Perform post-regrid routines. This does nothing
      */
      void
      postRegrid() override
      {}

      /*!
	@brief Set the CFL number
	@param[in] a_cfl;
      */
      void
      setCFL(const Real a_cfl);

      /*!
	@brief Set the initial data.. 
	@param[in] a_initData Initial data function. 
      */
      void
      setInitialData(const std::function<Real(const RealVect& a_position)>& a_initData) noexcept;

      /*!
	@brief Set the velocity field. 
	@param[in] a_velocity Velocity field. 
      */
      void
      setVelocity(const std::function<RealVect(const RealVect& a_position)>& a_velocity) noexcept;

      /*!
	@brief Set the diffusion coeffieint. 
	@param[in] a_diffusion Diffusion coefficient
      */
      void
      setDiffusionCoefficient(const std::function<Real(const RealVect& a_position)>& a_diffusion) noexcept;

    protected:
      /*!
	@brief Switch for different integrators
      */
      enum Integrator
      {
        Heun,
        IMEX
      };

      /*!
	@brief Reference to solver
      */
      RefCountedPtr<CdrSolver> m_solver;

      /*!
	@brief Reference to advected species
      */
      RefCountedPtr<AdvectionDiffusionSpecies> m_species;

      /*!
	@brief Realm where the solver lives
      */
      std::string m_realm;

      /*!
	@brief Phase where the solver lives
      */
      phase::which_phase m_phase;

      /*!
	@brief Velocity field.
      */
      std::function<RealVect(const RealVect& a_position)> m_velocity;

      /*!
	@brief Diffusion coefficient. 
      */
      std::function<Real(const RealVect& a_position)> m_diffCo;

      /*!
	@brief Initial data
      */
      std::function<Real(const RealVect& a_initialData)> m_initialData;

      /*!
	@brief Smallest acceptable time step
      */
      Real m_minDt;

      /*!
	@brief Largest acceptable time step
      */
      Real m_maxDt;

      /*!
	@brief CFL number to use
      */
      Real m_cfl;

      /*!
	@brief Forced CFL (overrides input options)
      */
      Real m_forceCFL;

      /*!
	@brief For activating debugging. This is an option.
      */
      bool m_debug;

      /*!
	@brief Mobile or not
      */
      bool m_mobile;

      /*!
	@brief Diffusive or not
      */
      bool m_diffusive;

      /*!
	@brief Integrator
      */
      Integrator m_integrator;

      /*!
	@brief Parse the integration method
      */
      void
      parseIntegrator();
    };
  } // namespace AdvectionDiffusion
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
