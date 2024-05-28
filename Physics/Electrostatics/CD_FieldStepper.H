/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_FieldStepper.H
  @brief  TimeStepper class for only solving the Poisson equation (with surface charge)
  @author Robert Marskar
*/

#ifndef CD_FieldStepper_H
#define CD_FieldStepper_H

// Our includes
#include <CD_TimeStepper.H>
#include <CD_SurfaceODESolver.H>
#include <CD_FieldSolver.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace Electrostatics {

    /*!
      @brief Implementation of TimeStepper for solving electrostatic problems. The template parameter is the FieldSolver implementation type. 
      @details This is just a one-step method -- it currently does not feature adaptive mesh refinement but uses the static grid set up by Driver. If the
      user wants to add mesh refinement he/she can add a CellTagger which flags cells for refinement. 

      The voltage on electrodes/domain faces is set to one. If the user wants to set a different voltage he will need to fetch the FieldSolver from the FieldStepper
      and call FieldStepper::setVoltage(...).

      The initial space charge is set to a Gaussian blob, specified in the input script. The initial surface charge is set to a constant value specified in the input script. 
    */
    template <class T>
    class FieldStepper : public TimeStepper
    {
    public:
      /*!
	@brief Constructor -- parses some input options
      */
      FieldStepper();

      /*!
	@brief Destructor (does nothing)
      */
      virtual ~FieldStepper();

      /*!
	@brief Solver setup routine. This instantiates the FieldSolver and parses input options
      */
      void
      setupSolvers() override;

      /*!
	@brief Set initial data -- this sets the space and surface charges
      */
      void
      initialData() override;

      /*!
	@brief Post-initialization routine. This solves the Poisson equation. 
      */
      void
      postInitialize() override;

      /*!
	@brief Allocation method -- allocates memory and internal data for solvers
      */
      void
      allocate() override;

      /*!
	@brief Get the potential.
      */
      MFAMRCellData&
      getPotential();

#ifdef CH_USE_HDF5
      /*!  
	@brief Write checkpoint data to file -- this physics module does not checkpoint data. 
	@param[inout] a_handle HDF5 fil
	@param[in]    a_lvl    Grid level
      */
      void
      writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
#endif

#ifdef CH_USE_HDF5
      /*!
	@brief Read checkpoint data from file -- this physics does not checkpoint data.
	@param[inout] a_handle HDF5 file
	@param[in]    a_lvl    Grid level
      */
      void
      readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
#endif

      /*!
	@brief Post-initialize operations to be performed -- does nothing. 
      */
      void
      postCheckpointSetup() override
      {}

      /*!
	@brief Get number of plot variables
	@return Returns plot variables for FieldSolver
      */
      int
      getNumberOfPlotVariables() const override;

      /*!
	@brief Get plot variable names
      */
      Vector<std::string>
      getPlotVariableNames() const override;

      /*!
	@brief Write plot data to file
	@param[inout] a_output            Output data holder.
	@param[inout] a_icomp             Starting component in a_output to begin at. 
	@param{in]    a_outputRealm       Realm where a_output belongs
	@param[in]    a_level             Grid level
      */
      void
      writePlotData(LevelData<EBCellFAB>& a_output,
                    int&                  a_icomp,
                    const std::string     a_outputRealm,
                    const int             a_level) const override;

      /*!
	@brief Compute time step
	@details This is a static solver so it returns 0.0
      */
      Real
      computeDt() override
      {
        return 0.0;
      }

      /*!
	@brief Perform a single time step with step a_dt.
	@param[in] a_dt Time step
	@return a_dt
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
	@brief Print step report -- this does nothing. 
      */
      void
      printStepReport() override{};

      /*!
	@brief Register simulation realms to be used for this simulation module. 
      */
      void
      registerRealms() override;

      /*!
	@brief Register operators for this simulation module
      */
      void
      registerOperators() override;

      /*!
	@brief Perform pre-regrid operations.
	@param[in] a_lmin           The coarsest level that changes
	@param[in] a_oldFinestLevel The finest level before the regrid. 
      */
      void
      preRegrid(const int a_lbase, const int a_finestLevel) override;

      /*!
	@brief Load balancing query for a specified realm. If this returns true for a_realm, load balancing routines will be called during regrids. 
	@param[in] a_realm Realm name
      */
      bool
      loadBalanceThisRealm(const std::string a_realm) const override;

      /*!
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

      /*!
	@brief Regrid method -- regrids the potential distribution in FieldSolver (but does not solve the Poisson equation, that is done in postRegrid()).
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      void
      regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /*!
	@brief Perform post-regrid operations -- this will resolve the Poisson equation. 
	@details This includes all operations to be done AFTER interpolating data to new grids. 
      */
      void
      postRegrid() override;

      /*!
	@brief Set space charge
	@param[in] a_rho Space charge distribution
      */
      void
      setRho(const std::function<Real(const RealVect& a_pos)>& a_rho, const phase::which_phase a_phase) noexcept;

      /*!
	@brief Set surface charge
	@param[in] a_sigma Surface charge distribution
      */
      void
      setSigma(const std::function<Real(const RealVect& a_pos)>& a_sigma) noexcept;

    protected:
      /*!
	@brief If true, the simulation will be introspectively load balanced. 
      */
      bool m_loadBalance;

      /*!
	@brief If using load balancing, this specifies how boxes will be sorted (in space) before load balancing. 
      */
      BoxSorting m_boxSort;

      /*!
	@brief Reference to the FieldSolver
      */
      RefCountedPtr<FieldSolver> m_fieldSolver;

      /*!
	@brief Reference to the surface charge solver.
      */
      RefCountedPtr<SurfaceODESolver<1>> m_sigma;

      /*!
	@brief Realm where the solver lives
      */
      std::string m_realm;

      /*!
	@brief Space charge inside gas phase.
      */
      std::function<Real(const RealVect& a_pos)> m_rhoGas;

      /*!
	@brief Space charge inside dielectric.
      */
      std::function<Real(const RealVect& a_pos)> m_rhoDielectric;

      /*!
	@brief Surface charge distribution
      */
      std::function<Real(const RealVect& a_pos)> m_surfaceChargeDensity;

      /*!
	@brief Internal routine for solving the Poisson equation
      */
      void
      solvePoisson();
    };
  } // namespace Electrostatics
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_FieldStepperImplem.H>

#endif
