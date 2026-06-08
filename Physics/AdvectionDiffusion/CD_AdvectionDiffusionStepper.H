/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_AdvectionDiffusionStepper.H
  @brief  Declaration of the Physics::AdvectionDiffusion::AdvectionDiffusionStepper class
  @author Robert Marskar
*/

#ifndef CD_ADVECTIONDIFFUSIONSTEPPER_H
#define CD_ADVECTIONDIFFUSIONSTEPPER_H

// Our includes
#include <CD_TimeStepper.H>
#include <CD_CdrSolver.H>
#include <CD_AdvectionDiffusionSpecies.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace AdvectionDiffusion {

    /**
      @brief TimeStepper implementation for advancing a single advection-diffusion equation.
      @details This class drives a single CdrSolver through a sequence of time steps. The velocity
      field, diffusion coefficient, and initial data are specified as functors and can be set either
      through ParmParse input or programmatically via the setter methods. Two explicit integrators
      are supported: a second-order Heun (explicit trapezoidal) method and an IMEX method that
      treats advection explicitly and diffusion implicitly via Crank-Nicholson.
    */
    class AdvectionDiffusionStepper : public TimeStepper
    {
    public:
      /**
        @brief Default constructor. Parses all options from ParmParse and sets default velocity,
        diffusion coefficient, and initial data functors.
      */
      AdvectionDiffusionStepper();

      /**
        @brief Full constructor. Delegates to the default constructor and then stores the provided solver.
        @param[in] a_solver CDR solver to use for the advection-diffusion advance.
      */
      AdvectionDiffusionStepper(RefCountedPtr<CdrSolver>& a_solver);

      /**
        @brief Destructor.
      */
      virtual ~AdvectionDiffusionStepper();

      /**
        @brief Fill the solver with initial data and set velocity, diffusion coefficient, and source term.
      */
      void
      initialData() override;

      /**
        @brief Allocate internal storage in the solver.
      */
      void
      allocate() override;

      /**
        @brief Post-initialization hook. No operations are performed.
      */
      void
      postInitialize() override
      {}

#ifdef CH_USE_HDF5
      /**
        @brief Write the solver state to an HDF5 checkpoint file.
        @param[inout] a_handle HDF5 file handle.
        @param[in]    a_lvl    AMR level to write.
      */
      void
      writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
#endif

#ifdef CH_USE_HDF5
      /**
        @brief Read the solver state from an HDF5 checkpoint file.
        @param[inout] a_handle HDF5 file handle.
        @param[in]    a_lvl    AMR level to read.
      */
      void
      readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
#endif

      /**
        @brief Post-checkpoint setup. Re-applies velocity, diffusion coefficient, source, and EB flux
        after a restart so the solver fields are consistent with the loaded state.
      */
      void
      postCheckpointSetup() override;

      /**
        @brief Return the number of output variables contributed by this stepper.
        @return Number of plot variables (delegated to the solver).
      */
      int
      getNumberOfPlotVariables() const override;

      /**
        @brief Return the names of the output variables contributed by this stepper.
        @return Variable name list (delegated to the solver).
      */
      Vector<std::string>
      getPlotVariableNames() const override;

      /**
        @brief Write plot data into the output holder.
        @param[inout] a_output      Output data holder.
        @param[inout] a_icomp       Starting component index in a_output; incremented on return.
        @param[in]    a_outputRealm Realm that owns a_output.
        @param[in]    a_level       AMR level to write.
      */
      void
      writePlotData(LevelData<EBCellFAB>& a_output,
                    int&                  a_icomp,
                    const std::string     a_outputRealm,
                    const int             a_level) const override;

      /**
        @brief Compute the next time step.
        @details For Heun integration the advection-diffusion CFL limit is used; for IMEX only the
        advection CFL limit is used since diffusion is treated implicitly. The result is clamped to
        [m_minDt, m_maxDt].
        @return Time step size in seconds.
      */
      Real
      computeDt() override;

      /**
        @brief Update internal time/step counters and forward them to the solver.
        @param[in] a_step Current time step number.
        @param[in] a_time Current simulation time in seconds.
        @param[in] a_dt   Time step that was used.
      */
      void
      synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;

      /**
        @brief Print a step report. No output is produced.
      */
      void
      printStepReport() override
      {}

      /**
        @brief Parse runtime options from ParmParse. Re-reads verbosity, CFL bounds, and integrator.
      */
      void
      parseRuntimeOptions() override;

      /**
        @brief Register the Primal realm with the AMR mesh.
      */
      void
      registerRealms() override;

      /**
        @brief Register the operators required by the solver.
      */
      void
      registerOperators() override;

      /**
        @brief Pre-regrid hook. Notifies the solver to store data before the mesh changes.
        @param[in] a_lbase          Coarsest level that will change.
        @param[in] a_oldFinestLevel Finest level before regridding.
      */
      void
      preRegrid(const int a_lbase, const int a_finestLevel) override;

      /**
        @brief Instantiate the species and configure the solver with realm, phase, and AMR data.
      */
      void
      setupSolvers() override;

      /**
        @brief Advance the solution by one time step using the selected integrator.
        @param[in] a_dt Time step size in seconds.
        @return The time step actually used (equal to a_dt).
      */
      Real
      advance(const Real a_dt) override;

      /**
        @brief Regrid the solver and re-apply the flow fields on the new mesh.
        @param[in] a_lmin           Coarsest level that changed.
        @param[in] a_oldFinestLevel Finest level before regridding.
        @param[in] a_newFinestLevel Finest level after regridding.
      */
      void
      regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /**
        @brief Post-regrid hook. No operations are performed.
      */
      void
      postRegrid() override
      {}

      /**
        @brief Force a specific CFL number, overriding the value read from input.
        @param[in] a_cfl CFL number to use.
      */
      void
      setCFL(const Real a_cfl);

      /**
        @brief Override the initial data functor.
        @param[in] a_initData Functor returning the initial scalar concentration at a given position.
      */
      void
      setInitialData(const std::function<Real(const RealVect& a_position)>& a_initData) noexcept;

      /**
        @brief Override the velocity field functor.
        @param[in] a_velocity Functor returning the advection velocity at a given position.
      */
      void
      setVelocity(const std::function<RealVect(const RealVect& a_position)>& a_velocity) noexcept;

      /**
        @brief Override the diffusion coefficient functor.
        @param[in] a_diffusion Functor returning the diffusion coefficient at a given position.
      */
      void
      setDiffusionCoefficient(const std::function<Real(const RealVect& a_position)>& a_diffusion) noexcept;

    protected:
      /**
        @brief Supported time integration methods.
      */
      enum Integrator
      {
        Heun, ///< Second-order explicit Heun (trapezoidal) method.
        IMEX  ///< IMEX: explicit advection, implicit (Crank-Nicholson) diffusion.
      };

      /**
        @brief The CDR solver used for advancing the scalar field.
      */
      RefCountedPtr<CdrSolver> m_solver;

      /**
        @brief Species definition passed to the CDR solver.
      */
      RefCountedPtr<AdvectionDiffusionSpecies> m_species;

      /**
        @brief Name of the realm on which the solver operates.
      */
      std::string m_realm;

      /**
        @brief Computational phase (gas or solid) on which the solver operates.
      */
      phase::which_phase m_phase;

      /**
        @brief Functor providing the advection velocity at a given position.
      */
      std::function<RealVect(const RealVect& a_position)> m_velocity;

      /**
        @brief Functor providing the diffusion coefficient at a given position.
      */
      std::function<Real(const RealVect& a_position)> m_diffCo;

      /**
        @brief Functor providing the initial scalar concentration at a given position.
      */
      std::function<Real(const RealVect& a_initialData)> m_initialData;

      /**
        @brief Minimum permissible time step in seconds.
      */
      Real m_minDt;

      /**
        @brief Maximum permissible time step in seconds.
      */
      Real m_maxDt;

      /**
        @brief CFL number read from input.
      */
      Real m_cfl;

      /**
        @brief Programmatically forced CFL number. Overrides m_cfl when positive.
      */
      Real m_forceCFL;

      /**
        @brief If true, print mass-conservation diagnostics to stdout after each step.
      */
      bool m_debug;

      /**
        @brief If true, the scalar is advected by the velocity field.
      */
      bool m_mobile;

      /**
        @brief If true, the scalar undergoes diffusion.
      */
      bool m_diffusive;

      /**
        @brief Selected time integration method.
      */
      Integrator m_integrator;

      /**
        @brief Parse the integration method from ParmParse input.
      */
      void
      parseIntegrator();
    };
  } // namespace AdvectionDiffusion
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
