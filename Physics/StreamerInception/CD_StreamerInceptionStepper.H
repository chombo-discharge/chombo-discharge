/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_StreamerInceptionStepper.H
  @brief  TimeStepper class for evaluating the streamer inception criterion.
  @author Robert Marskar
  @author Rasmus Hoholm
*/

#ifndef CD_StreamerInceptionStepper_H
#define CD_StreamerInceptionStepper_H

// Our includes
#include <CD_TimeStepper.H>
#include <CD_TracerParticleSolver.H>
#include <CD_TracerParticle.H>
#include <CD_FieldSolver.H>
#include <CD_FieldSolverMultigrid.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace StreamerInception {

    /*!
      @brief Integration algorithm for the particles
    */
    enum class IntegrationAlgorithm
    {
      Euler,
      Trapezoidal
    };

    /*!
      @brief Class for streamer inception integral evaluations. 
      @details P is the tracer particle type and F is the field solver type. 
      Note that the integration algorithm requires some transient storage on each particle (such as the previous position). 
      This is defined in the TracerParticle template arguments. 
    */
    template <typename P = TracerParticle<4, 2>, typename F = FieldSolverMultigrid>
    class StreamerInceptionStepper : public TimeStepper
    {
    public:
      /*!
	@brief Default constructor
      */
      StreamerInceptionStepper();

      /*!
	@brief Disallowed constructor
      */
      StreamerInceptionStepper(const StreamerInceptionStepper&) = delete;

      /*!
	@brief Disallowed constructor
      */
      StreamerInceptionStepper(const StreamerInceptionStepper&&) = delete;

      /*!
	@brief Destructor
      */
      virtual ~StreamerInceptionStepper();

      /*!
	@brief Disallowed assignment
      */
      StreamerInceptionStepper&
      operator=(const StreamerInceptionStepper&) = delete;

      /*!
	@brief Disallowed assignment
      */
      StreamerInceptionStepper&
      operator=(const StreamerInceptionStepper&&) = delete;

      /*!
	@brief Instantiate the tracer particle solver.
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
	@note Not needed for this time stepper beacuse checkpoint-restart is not supported. 
      */
      void
      postCheckpointSetup() override
      {}

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
	@note Module does not support checkpoint-restart and calling this is an error. 
      */
      virtual void
      writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
#endif

#ifdef CH_USE_HDF5
      /*!
	@brief Read checkpoint data from file
	@param[inout] a_handle HDF5 file
	@param[in]    a_lvl    Grid level.
	@note Module does not support checkpoint-restart and calling this is an error. 
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
	@brief Get the electric field
      */
      virtual const EBAMRCellData* const
      getElectricField() const noexcept;

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
	@brief Print a step report. Not used by this class. 
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

      /*!
	@brief Set the effective ionization coefficient.
	@param[in] a_alpha Effective ionization coefficient. E is the field in SI units. 
      */
      virtual void
      setAlpha(const std::function<Real(const Real& E)>& a_alpha) noexcept;

      /*!
	@brief Set the background ionization rate (e.g. from detachment, cosmic radiation etc). 
	@param[in] a_backgroundRate Background ionization rate (units of 1/m^3 s). 
      */
      virtual void
      setBackgroundRate(const std::function<Real(const Real& E)>& a_backgroundRate) noexcept;

      /*!
	@brief Set the field emission current
	@param[in] a_fieldEmission Field emission current density. 
      */
      virtual void
      setFieldEmission(const std::function<Real(const Real& E)>& a_currentDensity) noexcept;

    protected:
      /*!
	@brief Integration algorithm for K-value. 
      */
      IntegrationAlgorithm m_inceptionAlgorithm;

      /*!
	@brief Integration algorithm for drift integration. 
      */
      IntegrationAlgorithm m_driftAlgorithm;

      /*!
	@brief Tracer particle solver.
      */
      RefCountedPtr<TracerParticleSolver<P>> m_particleSolver;

      /*!
	@brief Field solver.
      */
      RefCountedPtr<FieldSolver> m_fieldSolver;

      /*!
	@brief Electric field for particles (with voltage = 1 V)
      */
      EBAMRCellData m_electricField;

      /*!
	@brief Inception integral values (positive polarity)
      */
      EBAMRCellData m_inceptionIntegralPlus;

      /*!
	@brief Inception integral values (negative polarity)
      */
      EBAMRCellData m_inceptionIntegralMinu;

      /*!
	@brief Drift times for positive polarity
      */
      EBAMRCellData m_driftTimesPlus;

      /*!
	@brief Drift times for negative polarity
      */
      EBAMRCellData m_driftTimesMinu;

      /*!
	@brief Background ionization rates
      */
      EBAMRCellData m_backgroundIonization;

      /*!
	@brief Field emission rates for positive voltage
      */
      EBAMRCellData m_emissionRatesPlus;

      /*!
	@brief Field emission rates for negative voltage
      */
      EBAMRCellData m_emissionRatesMinu;

      /*!
	@brief Computed inception voltage for positive polarity. 
      */
      EBAMRCellData m_inceptionVoltagePlus;

      /*!
	@brief Computed inception voltage for negative polarity. 
      */
      EBAMRCellData m_inceptionVoltageMinu;

      /*!
	@brief Variable step size for K-integration algorithm (factor of dx)
      */
      Real m_relativeInceptionStep;

      /*!
	@brief Variable step size for K-integration algorithm (factor of dx)
      */
      Real m_relativeDriftStep;

      /*!
	@brief Fixed step size (physical dimensions)
      */
      Real m_fixedInceptionStep;

      /*!
	@brief Fixed step size (physical dimensions)
      */
      Real m_fixedDriftStep;

      /*!
	@brief Inception criteria (read from input)
      */
      Real m_inceptionK;

      /*!
	@brief If true, used fixed integration step for inception integration
      */
      bool m_useFixedInceptionStep;

      /*!
	@brief If true, used fixed integration step for drift integration
      */
      bool m_useFixedDriftStep;

      /*!
	@brief If true, print a report at the end
      */
      bool m_printReport;

      /*!
	@brief If true, plot the K-values
      */
      bool m_plotInceptionIntegral;

      /*!
	@brief If true, plot the inception voltage
      */
      bool m_plotInceptionVoltage;

      /*!
	@brief If true, plot the drift time
      */
      bool m_plotDriftTime;

      /*!
	@brief If true, plot the background ionization rates
      */
      bool m_plotBackgroundIonization;

      /*!
	@brief If true, plot the field emission.
      */
      bool m_plotFieldEmission;

      /*! 
	@brief Max K value for every voltage (positive polarity)
      */
      std::vector<Real> m_maxKPlus;

      /*! 
	@brief Max K value for every voltage (negative polarity)
      */
      std::vector<Real> m_maxKMinu;

      /*!
	@brief Critical volumes of K values for each voltage (positive polarity)
      */
      std::vector<Real> m_criticalVolumePlus;

      /*!
	@brief Critical volumes of K values for each voltage (negative polarity)
      */
      std::vector<Real> m_criticalVolumeMinu;

      /*!
	@brief (Average) time to appearance of first electron (in critical volume)
	@note For positive voltage. 
      */
      std::vector<Real> m_firstElectronPlus;

      /*!
	@brief (Average) time to appearance of first electron (in critical volume)
	@note For negative voltage. 
      */
      std::vector<Real> m_firstElectronMinu;

      /*!
	@brief Voltage sweeps
      */
      std::vector<Real> m_voltageSweeps;

      /*!
	@brief Realm where the solver lives
      */
      std::string m_realm;

      /*!
	@brief Phase (gas/solid) where the solver lives
      */
      phase::which_phase m_phase;

      /*!
	@brief Effective ionization coefficient
      */
      std::function<Real(const Real& E)> m_alpha;

      /*!
	@brief Background ionization rate expression.
      */
      std::function<Real(const Real& E)> m_backgroundRate;

      /*!
	@brief Field emission current
      */
      std::function<Real(const Real& E)> m_fieldEmission;

      /*!
	@brief Distribute particles in every grid cell
      */
      virtual void
      initialParticles() noexcept;

      /*!
	@brief Solve streamer inception integral for each particle in each voltage and store K values in m_inceptionIntegral. 
	@details This is called in postInitialize().
      */
      virtual void
      computeInceptionIntegral() noexcept;

      /*!
	@brief Integrate the inception integral using the Euler rule.
	@param[in] a_voltage Voltage multiplier. 
      */
      virtual void
      inceptionIntegrateEuler(const Real& a_voltage) noexcept;

      /*!
	@brief K integral: Add integration parts after particles move
	@param[in] a_voltage   Voltage multiplier
      */
      virtual void
      inceptionIntegrateTrapezoidal(const Real& a_voltage) noexcept;

      /*!
	@brief Calculate particles' drift times until inside electrode or outside grid. 
      */
      virtual void
      computeDriftTimes() noexcept;

      /*!
	@brief Compute ion drift time using Euler rule
	@param[in] a_voltage Voltage multiplier
      */
      virtual void
      driftIntegrateEuler(const EBAMRCellData& a_K, const Real& a_voltage) noexcept;

      /*!
	@brief Compute ion drift time using Euler rule
	@param[in] a_voltage Voltage multiplier
      */
      virtual void
      driftIntegrateTrapezoidal(const EBAMRCellData& a_K, const Real& a_voltage) noexcept;

      /*!
	@brief Move particles back to their original position.
      */
      virtual void
      rewindParticles() noexcept;

      /*!
	@brief Reset particles. 
      */
      virtual void
      resetParticles() noexcept;

      /*!
	@brief Compute the detachment rate for all voltages
      */
      virtual void
      computeBackgroundIonization() noexcept;

      /*!
	@brief Compute field emission rates
      */
      virtual void
      computeFieldEmission() noexcept;

      /*!
	@brief Interpolate between K values to find voltage giving K_inception and store values in m_inceptionVoltage. 
	@details This is the volumetric version.
      */
      virtual void
      computeInceptionVoltageVolume() noexcept;

      /*!
	@brief Estimate the inception voltage
	@details This returns the inception voltage interpolated from the K-integral for positive and negative polarity. These
	SHOULD have the same value if the calculation is grid converged. 
      */
      virtual std::pair<Real, Real>
      computeInceptionVoltages() const noexcept;

      /*!
	@brief Compute the critical volume of the K values for each voltage
      */
      virtual void
      computeCriticalVolume() noexcept;

      /*!
	@brief Compute time to first critical electron.
      */
      virtual void
      computeFirstElectron() noexcept;

      /*!
	@brief Print report to the terminal. 
      */
      virtual void
      printReport() const noexcept;

      /*!
	@brief Check if particle is outside grid boundaries
	@param[in] a_pos  Particle position
	@param[in] a_probLo Lower left boundary of grid
	@param[in] a_probHi Upper right boundary of grid
	@return true if particle outside boundary, false elsewise
      */
      inline bool
      particleOutsideGrid(const RealVect& a_pos, const RealVect& a_probLo, const RealVect& a_probHi) const noexcept;

      /*!
	@brief Check if particle is inside electrode
	@param[in] a_pos Particle position
	@return true if particle inside electrode, false elsewise
      */
      inline bool
      particleInsideEB(const RealVect a_pos) const noexcept;
    };
  } // namespace StreamerInception
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_StreamerInceptionStepperImplem.H>

#endif
