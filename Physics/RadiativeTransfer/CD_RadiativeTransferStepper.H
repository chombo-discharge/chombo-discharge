/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RadiativeTransferStepper.H
  @brief  TimeStepper class for only evolving radiative transfer modules
  @author Robert Marskar
*/

#ifndef CD_RadiativeTransferStepper_H
#define CD_RadiativeTransferStepper_H

// Our includes
#include <CD_TimeStepper.H>
#include <CD_RtSolver.H>
#include <CD_NamespaceHeader.H>

namespace Physics {

  /*!
    @brief Namespace for encapsulating the radiative transfer physics module
  */
  namespace RadiativeTransfer {

    /*!
      @brief Implementation of TimeStepper for solving for a single radiative transfer species. This supports both Monte Carlo photons and Helmholtz approximations
      through the common RtSolver interface. 
    */
    template <class T>
    class RadiativeTransferStepper : public TimeStepper
    {
    public:
      /*!
	@brief Constructor. Reads a couple of input options.
      */
      RadiativeTransferStepper();

      /*!
	@brief Destructor (does nothing)
      */
      virtual ~RadiativeTransferStepper();

      /*!
	@brief Set up the solver
      */
      void
      setupSolvers() override;

      /*!
	@brief Fill simulation with initial data
      */
      void
      initialData() override;

      /*!
	@brief Allocate necessary memory for solvers.
      */
      void
      allocate() override;

      /*!
	@brief Post-initialization functionality -- sets the source term. 
      */
      void
      postInitialize() override;

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
	@brief Perform post-checkpoint setup routines (sets the source in the solver)
      */
      void
      postCheckpointSetup() override;

      /*!
	@brief Get number of plot variables for this physics module
      */
      int
      getNumberOfPlotVariables() const override;

      /*!
	@brief Get plot variable names
      */
      Vector<std::string>
      getPlotVariableNames() const override;

      /*!
	@brief Write plot data to output holder. 
	@param[inout] a_output Output data holder.
	@param[inout] a_icomp  Starting component in a_output to begin at. 
	@param[in]    a_outputRealm Realm where a_output belongs
	@param[in]    a_level  Grid level
      */
      void
      writePlotData(LevelData<EBCellFAB>& a_output,
                    int&                  a_icomp,
                    const std::string     a_outputRealm,
                    const int             a_level) const override;

      /*!
	@brief Compute a time step to be used by Driver. 
      */
      Real
      computeDt() override;

      /*!
	@brief Advancement method. Calls the solver function. 
	@param[in] a_dt Time step to be used for advancement
	@return    a_dt
      */
      Real
      advance(const Real a_dt) override;

      /*!
	@brief Synchronzie solver times and time steps
	@param[in] a_step Time step
	@param[in] a_time Time (in seconds)
	@param[in] a_dt   Time step that was used. 
      */
      void
      synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;

      /*!
	@brief Print a step report (does nothing)
	@details This is called by Driver after time step. 
      */
      void
      printStepReport() override;

      /*!
	@brief Parse runtime options.
      */
      void
      parseRuntimeOptions() override;

      /*!
	@brief Register realms
      */
      void
      registerRealms() override;

      /*!
	@brief Register operators -- calls the solver registration routine. 
      */
      void
      registerOperators() override;

      /*!
	@brief Perform pre-regrid operations -- calls the solver function. 
	operation takes a lot of memory. 
	@param[in] a_lmin           The coarsest level that changes
	@param[in] a_oldFinestLevel The finest level before the regrid. 
      */
      void
      preRegrid(const int a_lmin, const int a_oldFinestLevel) override;

      /*!
	@brief Regrid function. Calls the solver function. 
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      void
      regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /*!
	@brief Post-regrid function. Sets a Gaussian source. 
      */
      void
      postRegrid() override;

      /*!
	@brief Force usage of a time step
	@details To remove the forced time step restriction, set a_dt to something negative
	@param[in] a_dt Forced time step. 
      */
      void
      forceDt(const Real a_dt);

      /*!
	@brief Get the time step
      */
      Real
      getDt() const;

      /*!
	@brief Get the solver solution
      */
      const EBAMRCellData&
      getPhi() const;

    protected:
      /*!
	@brief Realm where solver lives.
      */
      std::string m_realm;

      /*!
	@brief Time step to use
      */
      Real m_dt;

      /*!
	@brief Forced time step. 
	@details This exists because the solver uses m_dt as a runtime input but there are cases when we want to force a time step.
      */
      Real m_forceDt;

      /*!
	@brief Reference to solver
      */
      RefCountedPtr<RtSolver> m_solver;

      /*!
	@brief Reference to species
      */
      RefCountedPtr<RtSpecies> m_species;

      /*!
	@brief For setting a Gaussian source in the radiative transfer equation solver
      */
      void
      setGaussianSource();
    };
  } // namespace RadiativeTransfer
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_RadiativeTransferStepperImplem.H>

#endif
