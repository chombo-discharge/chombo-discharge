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
      @brief Mode for specifying whether the module is run statically or dynamically
    */
    enum class Mode
    {
      None,
      Static,
      Dynamic
    };

    /*!
      @brief Class for streamer inception integral evaluations. 
      @details P is the tracer particle type and F is the field solver type. 
      Note that the integration algorithm requires some transient storage on each particle (such as the previous position). 
      This is defined in the TracerParticle template arguments. 
    */
    template <typename P = TracerParticle<2, 2>, typename F = FieldSolverMultigrid, typename C = CdrCTU>
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
	@note Returns the CELL-CENTROID centered electric field. 
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
	@brief Print a step report. Used in dynamic simulations. 
      */
      virtual void
      printStepReport() override;

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
	@brief Set the voltage curve (used for dynamic mode)
	@param[in] a_voltageCurve Voltage curve
      */
      virtual void
      setVoltageCurve(const std::function<Real(const Real& a_time)>& a_voltageCurve) noexcept;

      /*!
	@brief Set space charge distribution
      */
      virtual void
      setRho(const std::function<Real(const RealVect& x)>& a_rho) noexcept;

      /*!
	@brief Set the negative ion density
	@param[in] a_density Negative ion density
      */
      virtual void
      setNegativeIonDensity(const std::function<Real(const RealVect x)>& a_density) noexcept;

      /*!
	@brief Set the negative ion mobility (field-dependent)
	@param[in] a_mobility Negative ion  mobility
      */
      virtual void
      setNegativeIonMobility(const std::function<Real(const Real E)>& a_mobility) noexcept;

      /*!
	@brief Set the ionization coefficient.
	@param[in] a_alpha Townsend ionization coefficient. E is the field in SI units. 
      */
      virtual void
      setAlpha(const std::function<Real(const Real& E)>& a_alpha) noexcept;

      /*!
	@brief Set the attachment coefficient.
	@param[in] a_eta Townsend attachment coefficient. E is the field in SI units. 
      */
      virtual void
      setEta(const std::function<Real(const Real& E)>& a_eta) noexcept;

      /*!
	@brief Get ionization coefficient
      */
      virtual const std::function<Real(const Real& E)>&
      getAlpha() const noexcept;

      /*!
	@brief Get attachment coefficient
      */
      virtual const std::function<Real(const Real& E)>&
      getEta() const noexcept;

      /*!
	@brief Set the background ionization rate (e.g. from cosmic radiation etc). 
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

      /*!
	@brief Get the solver mode
      */
      virtual Mode
      getMode() const noexcept;

    protected:
      /*!
	@brief Mode
      */
      Mode m_mode;

      /*!
	@brief Integration algorithm for K-value. 
      */
      IntegrationAlgorithm m_inceptionAlgorithm;

      /*!
	@brief Tracer particle solver.
      */
      RefCountedPtr<TracerParticleSolver<P>> m_tracerParticleSolver;

      /*!
	@brief Negative ion solver
      */
      RefCountedPtr<TracerParticleSolver<P>> m_negativeIonSolver;

      /*!
	@brief Field solver.
      */
      RefCountedPtr<FieldSolver> m_fieldSolver;

      /*!
	@brief Electric field.
      */
      EBAMRCellData m_electricFieldPlus;

      /*!
	@brief Negated electric field (= -E)
      */
      EBAMRCellData m_electricFieldMinu;

      /*!
	@brief Inception integral values. 
	@note Dynamic mode. 
      */
      EBAMRCellData m_inceptionIntegral;

      /*!
	@brief Inception integral values.
	@note Static mode, positive polarity
      */
      EBAMRCellData m_inceptionIntegralPlus;

      /*!
	@brief Inception integral values.
	@note Static mode, negative polarity
      */
      EBAMRCellData m_inceptionIntegralMinu;

      /*!
	@brief Background ionization rates
	@note Both modes. 
      */
      EBAMRCellData m_backgroundIonization;

      /*!
	@brief Field emission rates
	@note Static mode, positive polarity
      */
      EBAMRCellData m_emissionRatesPlus;

      /*!
	@brief Field emission rates
	@note Static mode, negative polarity
      */
      EBAMRCellData m_emissionRatesMinu;

      /*!
	@brief Field emission rates
	@note Used for dynamic mode. 
      */
      EBAMRCellData m_emissionRate;

      /*!
	@brief Computed inception voltage
	@note Static mode, positive polarity
      */
      EBAMRCellData m_inceptionVoltagePlus;

      /*!
	@brief Computed inception voltage
	@note Static mode, negative polarity
      */
      EBAMRCellData m_inceptionVoltageMinu;

      /*! 
	@brief Max K value for every voltage
	@note Static mode, positive polarity
      */
      std::vector<Real> m_maxKPlus;

      /*! 
	@brief Max K value for every voltage
	@note Static mode, negative polarity
      */
      std::vector<Real> m_maxKMinu;

      /*!
	@brief Max K value for every time step
	@note For dynamic runs
      */
      std::vector<std::pair<Real, Real>> m_maxK;

      /*!
	@brief Critical volumes of K values for each voltage
	@note Static mode, positive polarity
      */
      std::vector<Real> m_criticalVolumePlus;

      /*!
	@brief Critical volumes of K values for each voltage
	@note Static mode, positive polarity
      */
      std::vector<Real> m_criticalVolumeMinu;

      /*!
	@brief Rate of appearance of first electron (in critical volume)
	@note Static mode, positive polarity. 
      */
      std::vector<Real> m_RdotPlus;

      /*!
	@brief Rate of appearance of first electron (in critical volume)
	@note Static mode, negative polarity. 
      */
      std::vector<Real> m_RdotMinu;

      /*!
	@brief Voltage sweeps
	@note Static mode
      */
      std::vector<Real> m_voltageSweeps;

      /*!
	@brief Integral factor for cumulative integration.
	@details This is the integral_Vcr dne/dt * (1 - eta/alpha) dV + integral(j/e*dS)  stored at each time step. 
	@note Only for dynamic mode. 
      */
      std::vector<std::pair<Real, Real>> m_Rdot;

      /*!
	@brief Cumulative probability of breakdown by time t.
	@details First parameter is time, second parameter is breakdown probability. 
      */
      std::vector<std::pair<Real, Real>> m_inceptionProbability;

      /*!
	@brief Critical volume at various time instances. 
	@details First parameter is time, second parameter is the critical volume. 
      */
      std::vector<std::pair<Real, Real>> m_criticalVolume;

      /*!
	@brief Time step for dynamic mode
      */
      Real m_fixedDt;

      /*!
	@brief CFL-like particle time step for negative ions. 
      */
      Real m_cfl;

      /*!
	@brief Maximum permitted time step
      */
      Real m_maxDt;

      /*!
	@brief Minimum permitted time step
      */
      Real m_minDt;

      /*!
	@brief Variable step size for K-integration algorithm (factor of dx)
      */
      Real m_relativeInceptionStep;

      /*!
	@brief Fixed step size (physical dimensions)
      */
      Real m_fixedInceptionStep;

      /*!
	@brief Inception criteria (read from input)
      */
      Real m_inceptionK;

      /*!
	@brief If true, used fixed integration step for inception integration
      */
      bool m_useFixedInceptionStep;

      /*!
	@brief Use fixed dt or not
      */
      bool m_useFixedDt;

      /*!
	@brief If true, print a report at the end
      */
      bool m_writeReport;

      /*!
	@brief Plot Poisson solver or not.
      */
      bool m_plotPoisson;

      /*!
	@brief Plot tracer particle solver or not. 
      */
      bool m_plotTracer;

      /*!
	@brief Plot negative ions or not
      */
      bool m_plotNegativeIons;

      /*!
	@brief If true, plot the K-values
      */
      bool m_plotInceptionIntegral;

      /*!
	@brief If true, plot the inception voltage
      */
      bool m_plotInceptionVoltage;

      /*!
	@brief If true, plot the background ionization rates
      */
      bool m_plotBackgroundIonization;

      /*!
	@brief If true, plot the field emission.
      */
      bool m_plotFieldEmission;

      /*!
	@brief If true, plot alpha
      */
      bool m_plotAlpha;

      /*!
	@brief If true, plot eta
      */
      bool m_plotEta;

      /*!
	@brief Time-dependent voltage
      */
      std::function<Real(const Real a_time)> m_voltageCurve;

      /*!
	@brief Realm where the solver lives
      */
      std::string m_realm;

      /*!
	@brief Output file
      */
      std::string m_outputFile;

      /*!
	@brief Phase (gas/solid) where the solver lives
      */
      phase::which_phase m_phase;

      /*!
	@brief Space charge distribution
      */
      std::function<Real(const RealVect& x)> m_rho;

      /*!
	@brief Negative ion density
      */
      std::function<Real(const RealVect& x)> m_negativeIonDensity;

      /*!
	@brief Negative ion mobility
      */
      std::function<Real(const Real& E)> m_negativeIonMobility;

      /*!
	@brief Detachment probability for negative ion.s
      */
      std::function<Real(const Real& E)> m_detachmentProbability;

      /*!
	@brief Ionization coefficient. 
      */
      std::function<Real(const Real& E)> m_alpha;

      /*!
	@brief Attachment coefficient.
      */
      std::function<Real(const Real& E)> m_eta;

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
      seedUniformParticles() noexcept;

      /*!
	@brief Add particles to every cell where alpha - eta > 0.0
      */
      virtual void
      seedIonizationParticles(const Real& a_voltage) noexcept;

      /*!
	@brief Create negative ions
      */
      virtual void
      seedNegativeIons() noexcept;

      /*!
	@brief Solve streamer inception integral for each particle in each voltage and store K values in m_inceptionIntegral. 
	@details This is called in postInitialize() only. 
	@note For static mode only. 
      */
      virtual void
      computeInceptionIntegralStatic() noexcept;

      /*!
	@brief Solve streamer inception integral 
	@details Called in postInitialize and the advance method
	@note Dynamic mode only. 
	@param[in] a_voltage Voltage multiplier
      */
      virtual void
      computeInceptionIntegralDynamic(const Real& a_voltage) noexcept;

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
	@brief Move negative ions using
	@param[in] a_voltage Voltage multiplier
	@param[in] a_dt      Time step. 
	@note Only for dynamic runs. 
      */
      virtual void
      advanceNegativeIons(const Real& a_voltage, const Real& a_dt) noexcept;

      /*!
	@brief Move negative ions using an Euler scheme. 
	@param[in] a_voltage Voltage multiplier
	@param[in] a_dt      Time step. 
	@note Only for dynamic runs. 
      */
      virtual void
      advanceNegativeIonsEuler(const Real& a_voltage, const Real& a_dt) noexcept;

      /*!
	@brief Move negative ions using Heun's method. 
	@param[in] a_voltage Voltage multiplier
	@param[in] a_dt      Time step. 
	@note Only for dynamic runs. 
      */
      virtual void
      advanceNegativeIonsTrapezoidal(const Real& a_voltage, const Real& a_dt) noexcept;

      /*!
	@brief Compute integral_Vcr(dne/dt * (1 - eta/alpha) dV)
	@param[in] a_voltage Voltage multiplier. 
      */
      virtual Real
      computeRdot(const Real& a_voltage) const noexcept;

      /*!
	@brief Move particles back to their original position.
      */
      virtual void
      rewindTracerParticles() noexcept;

      /*!
	@brief Reset particles. 
      */
      virtual void
      resetTracerParticles() noexcept;

      /*!
	@brief Compute the background ionization rate for all voltages
	@note For static simulations. 
      */
      virtual void
      computeBackgroundIonizationStatic() noexcept;

      /*!
	@brief Compute field emission rates
	@note For static simulations. 
      */
      virtual void
      computeFieldEmissionStatic() noexcept;

      /*!
	@brief Compute field emission rates
	@param[out] a_emissionRate Field emission rate. 
	@param[in]  a_voltage Voltage multiplier
	@note For plotting the emission rate in dynamic simulations. 
      */
      virtual void
      computeFieldEmission(EBAMRCellData& a_emissionRate, const Real& a_voltage) const noexcept;

      /*!
	@brief Evaluate a function f = f(E) in a volume. 
	@param[out] a_data    Mesh data
	@param[in]  a_voltage Voltage multiplier
	@param[im]  a_func    Function to evaluate
	@note For plotting the emission rate in dynamic simulations. 
      */
      virtual void
      evaluateFunction(EBAMRCellData&                           a_emissionRate,
                       const Real&                              a_voltage,
                       const std::function<Real(const Real E)>& a_func) const noexcept;

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
	@note For the static method
      */
      virtual void
      computeCriticalVolumeStatic() noexcept;

      /*!
	@brief Compute the critical volume of the K values for each voltage
	@note Dynamic approach. 
      */
      virtual Real
      computeCriticalVolumeDynamic() const noexcept;

      /*!
	@brief Compute time to first critical electron.
      */
      virtual void
      computeRdotStatic() noexcept;

      /*!
	@brief Print report to the terminal. 
	@note For static runs
      */
      virtual void
      writeReportStatic() const noexcept;

      /*!
	@brief Print report to the terminal. 
	@note For dynamic runs. Called in the constructor. 
      */
      virtual void
      writeReportDynamic() const noexcept;

      /*!
	@brief Write plot data for the 'static' mode
	@param[inout] a_output            Output data holder.
	@param[in]    a_plotVariableNames Plot variable names.
	@param[inout] a_icomp             Starting component in a_output to begin at. 
      */
      virtual void
      writePlotDataStatic(EBAMRCellData&       a_output,
                          Vector<std::string>& a_plotVariableNames,
                          int&                 a_icomp) const noexcept;

      /*!
	@brief Write plot data for the 'static' mode
	@param[inout] a_output            Output data holder.
	@param[in]    a_plotVariableNames Plot variable names.
	@param[inout] a_icomp             Starting component in a_output to begin at. 
      */
      virtual void
      writePlotDataDynamic(EBAMRCellData&       a_output,
                           Vector<std::string>& a_plotVariableNames,
                           int&                 a_icomp) const noexcept;

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

      /*!
	@brief Compute time step such that no ions move more than grid cell (on their current level)
	@param[in] a_voltage Voltage multiplier
	@return Smallest time step. 
	@note Particle velocities should be interpolated before calling this routine. 
      */
      inline Real
      computeNegativeIonCFL(const Real& a_voltage) const noexcept;
    };
  } // namespace StreamerInception
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_StreamerInceptionStepperImplem.H>

#endif
