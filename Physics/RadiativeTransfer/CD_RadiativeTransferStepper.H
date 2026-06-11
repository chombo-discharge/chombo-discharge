/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_RadiativeTransferStepper.H
  @brief  Declaration of the Physics::RadiativeTransfer::RadiativeTransferStepper TimeStepper
  @author Robert Marskar
*/

#ifndef CD_RADIATIVETRANSFERSTEPPER_H
#define CD_RADIATIVETRANSFERSTEPPER_H

// Our includes
#include <CD_TimeStepper.H>
#include <CD_RtSolver.H>
#include <CD_NamespaceHeader.H>

namespace Physics {

  /**
    @brief Namespace for encapsulating the radiative transfer physics module
  */
  namespace RadiativeTransfer {

    /**
      @brief TimeStepper for advancing a single radiative transfer species on an AMR mesh.
      @details Supports both Monte Carlo photon transport and Eddington/Helmholtz continuum
      approximations through the common RtSolver interface. The template parameter @p T is the
      concrete RtSolver type. A source term modelled as a Gaussian blob is applied to the
      solver at initialization and after every regrid; the shape is set via the input script.
    */
    template <class T>
    class RadiativeTransferStepper : public TimeStepper
    {
    public:
      /**
	@brief Constructor. Reads a couple of input options.
      */
      RadiativeTransferStepper();

      /**
	@brief Destructor (does nothing)
      */
      virtual ~RadiativeTransferStepper();

      /**
	@brief Set up the solver
      */
      void
      setupSolvers() override;

      /**
	@brief Fill simulation with initial data
      */
      void
      initialData() override;

      /**
	@brief Allocate necessary memory for solvers.
      */
      void
      allocate() override;

      /**
	@brief Post-initialization functionality -- sets the source term. 
      */
      void
      postInitialize() override;

#ifdef CH_USE_HDF5
      /**  
	@brief Write checkpoint data to file
	@param[in,out] a_handle HDF5 fil
	@param[in]    a_lvl    Grid level
      */
      void
      writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
#endif

#ifdef CH_USE_HDF5
      /**
	@brief Read checkpoint data from file
	@param[in,out] a_handle HDF5 file
	@param[in]    a_lvl    Grid level
      */
      void
      readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
#endif

      /**
	@brief Perform post-checkpoint setup routines (sets the source in the solver)
      */
      void
      postCheckpointSetup() override;

      /**
	@brief Get number of plot variables for this physics module.
	@return Number of plot variables contributed by the RtSolver.
      */
      int
      getNumberOfPlotVariables() const override;

      /**
	@brief Get plot variable names.
	@return Vector of variable name strings.
      */
      Vector<std::string>
      getPlotVariableNames() const override;

      /**
	@brief Write plot data to output holder. 
	@param[in,out] a_output Output data holder.
	@param[in,out] a_icomp  Starting component in a_output to begin at. 
	@param[in]    a_outputRealm Realm where a_output belongs
	@param[in]    a_level  Grid level
      */
      void
      writePlotData(LevelData<EBCellFAB>& a_output,
                    int&                  a_icomp,
                    const std::string&    a_outputRealm,
                    const int             a_level) const override;

      /**
	@brief Compute a time step to be used by Driver.
	@return The fixed time step m_dt read from input options.
      */
      Real
      computeDt() override;

      /**
	@brief Advancement method. Calls the solver function. 
	@param[in] a_dt Time step to be used for advancement
	@return    a_dt
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
	@brief Print a step report (does nothing)
	@details This is called by Driver after time step. 
      */
      void
      printStepReport() override;

      /**
	@brief Parse runtime options.
      */
      void
      parseRuntimeOptions() override;

      /**
	@brief Register realms
      */
      void
      registerRealms() override;

      /**
	@brief Register operators -- calls the solver registration routine. 
      */
      void
      registerOperators() override;

      /**
	@brief Perform pre-regrid operations — calls the solver function.
	@param[in] a_lmin           The coarsest level that changes.
	@param[in] a_oldFinestLevel The finest level before the regrid.
      */
      void
      preRegrid(const int a_lmin, const int a_oldFinestLevel) override;

      /**
	@brief Regrid function. Calls the solver function. 
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      void
      regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /**
	@brief Post-regrid function. Sets a Gaussian source. 
      */
      void
      postRegrid() override;

      /**
	@brief Force usage of a time step
	@details To remove the forced time step restriction, set a_dt to something negative
	@param[in] a_dt Forced time step. 
      */
      void
      forceDt(const Real a_dt);

      /**
	@brief Get the time step.
	@return The current time step m_dt.
      */
      Real
      getDt() const;

      /**
	@brief Get the solver solution.
	@return Const reference to the RtSolver's photon density field on the AMR hierarchy.
      */
      const EBAMRCellData&
      getPhi() const;

    protected:
      /**
	@brief Realm where solver lives.
      */
      std::string m_realm;

      /**
	@brief Time step to use
      */
      Real m_dt;

      /**
	@brief Forced time step. 
	@details This exists because the solver uses m_dt as a runtime input but there are cases when we want to force a time step.
      */
      Real m_forceDt;

      /**
	@brief Reference to solver
      */
      RefCountedPtr<RtSolver> m_solver;

      /**
	@brief Reference to species
      */
      RefCountedPtr<RtSpecies> m_species;

      /**
	@brief Set a Gaussian source term in the RtSolver.
	@details Reads the blob center, amplitude, and width from the input script and fills the
	solver's source field with a Gaussian distribution. Called after initialization and each regrid.
      */
      void
      setGaussianSource();
    };
  } // namespace RadiativeTransfer
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_RadiativeTransferStepperImplem.H>

#endif
