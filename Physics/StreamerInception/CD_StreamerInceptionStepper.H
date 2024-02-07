/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Copyright © 2022 NTNU.       
 * Copyright © 2022 Fanny Skirbekk.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_StreamerInceptionStepper.H
  @brief  TimeStepper class for evaluating the streamer inception criterion.
  @author Robert Marskar (SINTEF)
  @author Rasmus Hoholm (SINTEF)
  @author Fanny Skirbekk (NTNU)
*/

#ifndef CD_StreamerInceptionStepper_H
#define CD_StreamerInceptionStepper_H

// Our includes
#include <CD_TimeStepper.H>
#include <CD_TracerParticleSolver.H>
#include <CD_TracerParticle.H>
#include <CD_FieldSolver.H>
#include <CD_FieldSolverMultigrid.H>
#include <CD_CdrSolver.H>
#include <CD_CdrCTU.H>
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
      @brief Step size selection method
    */
    enum class StepSizeMethod
    {
      Fixed,
      Dx,
      Alpha
    };

    /*!
      @brief Algorithm for negative ion advancement
    */
    enum class TransportAlgorithm
    {
      Euler,
      Heun,
      ImExCTU
    };

    /*!
      @brief For specifying whether the module is run in stationary or transient mode
    */
    enum class Mode
    {
      Stationary,
      Transient
    };

    /*!
      @brief For specifying how the time step was restricted
    */
    enum class TimeStepRestriction
    {
      Unknown,
      CDR,
      VoltageCurve,
      MinHardcap,
      MaxHardcap
    };

    /*!
      @brief Class for streamer inception integral evaluations. 
      @details P is the tracer particle type
               F is the field solver type. 
               C is the convection-diffusion-reaction solver type. 
    */
    template <typename P = TracerParticle<2, 2>, typename F = FieldSolverMultigrid, typename C = CdrCTU>
    class StreamerInceptionStepper : public TimeStepper
    {
    public:
      EBAMRCellData m_electricFieldPlus;
      EBAMRCellData m_electricFieldMinu;
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
	@brief Get the number of plot variables for this time stepper. 
	@details This is necessary because Driver, not TimeStepper, is responsible for allocating the necessary memory. 
	@return Returns number of plot variables that will be written during writePlotData
      */
      virtual int
      getNumberOfPlotVariables() const override;

      /*!
	@brief Get plot variable names
      */
      virtual Vector<std::string>
      getPlotVariableNames() const override;

      /*!
	@brief Write plot data to output holder. 
	@param[inout] a_output Output data holder.
	@param[inout] a_icomp  Starting component in a_output to begin at. 
	@param[in]    a_outputRealm Realm where a_output belongs
	@param[in]    a_level  Grid level
      */
      virtual void
      writePlotData(LevelData<EBCellFAB>& a_output,
                    int&                  a_icomp,
                    const std::string     a_outputRealm,
                    const int             a_level) const override;

      /*!
	@brief Get plot variable names for stationary mode
      */
      virtual Vector<std::string>
      getStationaryPlotVariableNames() const noexcept;

      /*!
	@brief Get plot variable names for transient mode
      */
      virtual Vector<std::string>
      getTransientPlotVariableNames() const noexcept;

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
	@brief Advance negative ions
      */
      virtual void
      advanceIons(const Real a_dt) noexcept;

      /*!
	@brief Synchronize solver times and time steps
	@param[in] a_step Time step
	@param[in] a_time Time (in seconds)
	@param[in] a_dt   Time step that was used. 
      */
      virtual void
      synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;

      /*!
	@brief Print a step report. Used in transient simulations. 
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
	@brief Set the voltage curve (used for transient mode)
	@param[in] a_voltageCurve Voltage curve
      */
      virtual void
      setVoltageCurve(const std::function<Real(const Real& a_time)>& a_voltageCurve) noexcept;

      /*!
	@brief Set space charge distribution
	@param[in] a_rho Space charge distribution
      */
      virtual void
      setRho(const std::function<Real(const RealVect& x)>& a_rho) noexcept;

      /*!
	@brief Set surface charge distribution
	@param[in] a_sigma Surface charge distribution
      */
      virtual void
      setSigma(const std::function<Real(const RealVect& x)>& a_sigma) noexcept;

      /*!
	@brief Set the negative ion density
	@param[in] a_density Negative ion density
      */
      virtual void
      setIonDensity(const std::function<Real(const RealVect x)>& a_density) noexcept;

      /*!
	@brief Set the negative ion mobility (field-dependent)
	@param[in] a_mobility Negative ion  mobility
      */
      virtual void
      setIonMobility(const std::function<Real(const Real E)>& a_mobility) noexcept;

      /*!
	@brief Set the negative ion diffusion coefficient (field-dependent)
	@param[in] a_diffCo Negative ion diffusion coefficient
      */
      virtual void
      setIonDiffusion(const std::function<Real(const Real E)>& a_diffCo) noexcept;

      /*!
	@brief Set the ionization coefficient.
	@param[in] a_alpha Townsend ionization coefficient. E is the field in SI units. 
      */
      virtual void
      setAlpha(const std::function<Real(const Real& E, const RealVect& x)>& a_alpha) noexcept;

      /*!
	@brief Set the attachment coefficient.
	@param[in] a_eta Townsend attachment coefficient. E is the field in SI units. 
      */
      virtual void
      setEta(const std::function<Real(const Real& E, const RealVect& x)>& a_eta) noexcept;

      /*!
	@brief Get ionization coefficient
      */
      virtual const std::function<Real(const Real& E, const RealVect& x)>&
      getAlpha() const noexcept;

      /*!
	@brief Get attachment coefficient
      */
      virtual const std::function<Real(const Real& E, const RealVect& x)>&
      getEta() const noexcept;

      /*!
	@brief Set the background ionization rate (e.g. from cosmic radiation etc). 
	@param[in] a_backgroundRate Background ionization rate (units of 1/m^3 s). 
      */
      virtual void
      setBackgroundRate(const std::function<Real(const Real& E, const RealVect& x)>& a_backgroundRate) noexcept;

      /*!
	@brief Set the detachment rate for negative ions.
	@details The total detachment rate is dn_e/dt = k*n_ion; this sets the constant k. 
	@param[in] a_detachmentRate Detachment rate (in units of 1/s). 
      */
      virtual void
      setDetachmentRate(const std::function<Real(const Real& E, const RealVect& x)>& a_detachmentRate) noexcept;

      /*!
	@brief Set the field emission current
	@param[in] a_fieldEmission Field emission current density. 
      */
      virtual void
      setFieldEmission(const std::function<Real(const Real& E, const RealVect& x)>& a_currentDensity) noexcept;

      /*!
	@brief Set the secondary emission coefficient
	@param[in] a_coeff Secondary emission coefficient. 
      */
      virtual void
      setSecondaryEmission(const std::function<Real(const Real& E, const RealVect& x)>& a_coeff) noexcept;

      /*!
	@brief Get the solver mode
      */
      virtual Mode
      getMode() const noexcept;

      /*!
	@brief Get the electric field
      */
      virtual const EBAMRCellData*
      getElectricField() const noexcept
      {
        return &m_electricFieldHomo;
      }

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
	@brief Step size selection method
      */
      StepSizeMethod m_stepSizeMethod;

      /*!
	@brief Transport algorithm for advancing negative ions in the transient mode
      */
      TransportAlgorithm m_transportAlgorithm;

      /*!
	@brief Time step restriction
      */
      TimeStepRestriction m_timeStepRestriction;

      /*!
	@brief Tracer particle solver.
      */
      RefCountedPtr<TracerParticleSolver<P>> m_tracerParticleSolver;

      /*!
	@brief Field solver.
      */
      RefCountedPtr<FieldSolver> m_fieldSolver;

      /*!
	@brief Negative ion solver.
      */
      RefCountedPtr<CdrSolver> m_ionSolver;

      /*!
	@brief Storage for electric field. Obtained via superposition
      */
      EBAMRCellData m_electricField;

      /*!
	@brief Electric field without space/surface charge and V = 1 on live electrodes
      */
      EBAMRCellData m_electricFieldHomo;

      /*!
	@brief Electric field with space/surface charge and V = 0 on all electrodes
      */
      EBAMRCellData m_electricFieldInho;

      /*!
	@brief Inception integral values. 
	@note Transient mode. 
      */
      EBAMRCellData m_inceptionIntegral;

      /*!
	@brief Inception integral values.
	@note Stationary mode, positive polarity
      */
      EBAMRCellData m_inceptionIntegralPlus;

      /*!
	@brief Inception integral values.
	@note Stationary mode, negative polarity
      */
      EBAMRCellData m_inceptionIntegralMinu;

      /*!
	@brief Background ionization rates
	@note Both modes. 
      */
      EBAMRCellData m_backgroundIonization;

      /*!
	@brief Detachment rates.
	@note Both modes. 
      */
      EBAMRCellData m_detachment;

      /*!
	@brief Field emission rates
	@note Stationary mode, positive polarity
      */
      EBAMRCellData m_emissionRatesPlus;

      /*!
	@brief Field emission rates
	@note Stationary mode, negative polarity
      */
      EBAMRCellData m_emissionRatesMinu;

      /*!
	@brief Field emission rates
	@note Used for transient mode. 
      */
      EBAMRCellData m_emissionRate;

      /*!
	@brief Computed inception voltage
	@note Stationary mode, positive polarity
      */
      EBAMRCellData m_inceptionVoltagePlus;

      /*!
	@brief Computed inception voltage
	@note Stationary mode, negative polarity
      */
      EBAMRCellData m_inceptionVoltageMinu;

      /*!
	@brief Computed streamer inception voltage
	@note Stationary mode, positive polarity
      */
      EBAMRCellData m_streamerInceptionVoltagePlus;

      /*!
	@brief Computed streamer inception voltage
	@note Stationary mode, negative polarity
      */
      EBAMRCellData m_streamerInceptionVoltageMinu;

      /*!
	@brief Computed Townsend inception voltage
	@note Stationary mode, positive polarity
      */
      EBAMRCellData m_townsendInceptionVoltagePlus;

      /*!
	@brief Computed Townsend inception voltage
	@note Stationary mode, negative polarity
      */
      EBAMRCellData m_townsendInceptionVoltageMinu;

      /*!
	@brief Townsend criterion
	@note Transient mode. 
      */
      EBAMRCellData m_townsendCriterion;

      /*!
	@brief Townsend criterion for starting positive ions
      */
      EBAMRCellData m_townsendCriterionPlus;

      /*!
	@brief Townsend criterion for starting positive ions
      */
      EBAMRCellData m_townsendCriterionMinu;

      /*! 
	@brief Max K value for every voltage
	@note Stationary mode, positive polarity
      */
      std::vector<Real> m_maxKPlus;

      /*! 
	@brief Max K value for every voltage
	@note Stationary mode, negative polarity
      */
      std::vector<Real> m_maxKMinu;

      /*! 
	@brief Max Townsend value (exp(K)/gamma) for every voltage
	@note Stationary mode, positive polarity
      */
      std::vector<Real> m_maxTPlus;

      /*! 
	@brief Max Townsend value (exp(K)/gamma) for every voltage
	@note Stationary mode, negative polarity
      */
      std::vector<Real> m_maxTMinu;

      /*!
	@brief Max K value for every time step
	@note For transient runs
      */
      std::vector<std::pair<Real, Real>> m_maxK;

      /*!
	@brief Max Townsend value for every time step
	@note For transient runs
      */
      std::vector<std::pair<Real, Real>> m_maxT;

      /*!
	@brief Critical volumes of K values for each voltage
	@note Stationary mode, positive polarity
      */
      std::vector<Real> m_criticalVolumePlus;

      /*!
	@brief Critical volumes of K values for each voltage
	@note Stationary mode, negative polarity
      */
      std::vector<Real> m_criticalVolumeMinu;

      /*!
	@brief Critical area of K values for each voltage
	@note Stationary mode, positive polarity
      */
      std::vector<Real> m_criticalAreaPlus;

      /*!
	@brief Critical area of K values for each voltage
	@note Stationary mode, negative polarity
      */
      std::vector<Real> m_criticalAreaMinu;

      /*!
	@brief Ionization volumes for each voltage
	@note Stationary mode
      */
      std::vector<Real> m_ionizationVolume;

      /*!
	@brief Rate of appearance of first electron (in critical volume)
	@note Stationary mode, positive polarity. 
      */
      std::vector<Real> m_RdotPlus;

      /*!
	@brief Rate of appearance of first electron (in critical volume)
	@note Stationary mode, negative polarity. 
      */
      std::vector<Real> m_RdotMinu;

      /*!
	@brief Voltage sweeps
	@note Stationary mode
      */
      std::vector<Real> m_voltageSweeps;

      /*!
	@brief Integral factor for cumulative integration.
	@details This is the integral_Vcr dne/dt * (1 - eta/alpha) dV + integral(j/e*dS)  stored at each time step. 
	@note Only for transient mode. 
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
	@brief Critical area at various time instances. 
	@details First parameter is time, second parameter is the critical area. 
      */
      std::vector<std::pair<Real, Real>> m_criticalArea;

      /*!
	@brief Ionization volume at various time instances.
	@details First parameter is time, second parameter is the ionization volume.
      */
      std::vector<std::pair<Real, Real>> m_ionizationVolumeTransient;

      /*!
	@brief CFL-step for negative ions.
      */
      Real m_cfl;

      /*!
	@brief First time step
      */
      Real m_firstDt;

      /*!
	@brief Maximum permitted time step
      */
      Real m_maxDt;

      /*!
	@brief Minimum permitted time step.
      */
      Real m_minDt;

      /*!
	@brief Permitted relative deviation in V(t) when computing voltage-curve based time steps.
      */
      Real m_epsVoltage;

      /*!
	@brief Maximum permitted time step growth 
      */
      Real m_maxDtGrowth;

      /*!
	@brief Fixed step size (physical dimensions)
	@brief Step size factor
      */
      Real m_stepSizeFactor;

      /*!
	@brief Inception criteria (read from input)
      */
      Real m_inceptionK;

      /*!
	@brief Profile calculation or not
      */
      bool m_profile;

      /*!
	@brief Debug or not
      */
      bool m_debug;

      /*!
	@brief Evaluate secondary emission
      */
      bool m_evaluateTownsend;

      /*!
	@brief Use full reconstruction of critical volume or not
      */
      bool m_fullIntegration;

      /*!
	@brief Ion transport on/off
      */
      bool m_ionTransport;

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
	@brief If true, plot electron detachment rate, i.e. dne/dt due to detachment from negative ions. 
      */
      bool m_plotDetachment;

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
	@brief Plot the Townsend-criterion for initiatory ions
      */
      bool m_plotTownsend;

      /*!
	@brief Time-dependent voltage
      */
      std::function<Real(const Real a_time)> m_voltageCurve;

      /*!
	@brief Realm where the solver lives
      */
      std::string m_realm;

      /*!
	@brief output data file
      */
      std::string m_dataOutputFile;

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
	@brief Surface charge distribution
      */
      std::function<Real(const RealVect& x)> m_sigma;

      /*!
	@brief Negative ion density
      */
      std::function<Real(const RealVect& x)> m_initialIonDensity;

      /*!
	@brief Negative ion mobility
      */
      std::function<Real(const Real& E)> m_ionMobility;

      /*!
	@brief Negative ion diffusion coefficient
      */
      std::function<Real(const Real& E)> m_ionDiffusion;

      /*!
	@brief Ionization coefficient. 
      */
      std::function<Real(const Real& E, const RealVect& x)> m_alpha;

      /*!
	@brief Attachment coefficient.
      */
      std::function<Real(const Real& E, const RealVect& x)> m_eta;

      /*!
	@brief Detachment probability for negative ion.
      */
      std::function<Real(const Real& E, const RealVect& x)> m_detachmentRate;

      /*!
	@brief Background ionization rate expression.
      */
      std::function<Real(const Real& E, const RealVect& x)> m_backgroundRate;

      /*!
	@brief Field emission current
      */
      std::function<Real(const Real& E, const RealVect& x)> m_fieldEmission;

      /*!
	@brief Secondary emission coefficient
      */
      std::function<Real(const Real& E, const RealVect& x)> m_secondaryEmission;

      /*!
	@brief Distribute particles in every grid cell
      */
      virtual void
      seedUniformParticles() noexcept;

      /*!
	@brief Add particles to every cell where alpha - eta > 0.0
	@param[in] a_voltage Input voltage
      */
      virtual void
      seedIonizationParticles(const Real a_voltage) noexcept;

      /*!
	@brief Solve streamer inception integral for each particle in each voltage and store K values in m_inceptionIntegral. 
	@details This is called in postInitialize() only. 
	@note For stationary mode only. 
      */
      virtual void
      computeInceptionIntegralStationary() noexcept;

      /*!
	@brief Solve streamer inception integral 
	@details Called in postInitialize and the advance method
	@note Transient mode only. 
	@param[in] a_voltage Voltage multiplier
      */
      virtual void
      computeInceptionIntegralTransient(const Real& a_voltage) noexcept;

      /*!
	@brief Integrate the inception integral using the Euler rule.
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
	@brief Solve for the Townsend criterion for each particle in each voltage. 
	@details This is called in postInitialize() only. 
	@note For stationary mode only. 
      */
      virtual void
      computeTownsendCriterionStationary() noexcept;

      /*!
	@brief Solve for the Townsend criterion for each particle in each voltage. 
	@param[in] a_voltage Voltage multiplier
	@note For transient mode
      */
      virtual void
      computeTownsendCriterionTransient(const Real& a_voltage) noexcept;

      /*!
	@brief Track particles (positive ions) using an Euler rule and check if the collide with a cathode
	@param[in] a_voltage Voltage multiplier
      */
      virtual void
      townsendTrackEuler(const Real& a_voltage) noexcept;

      /*!
	@brief Track particles (positive ions) using a trapezoidal rule and check if the collide with a cathode
	@param[in] a_voltage Voltage multiplier
      */
      virtual void
      townsendTrackTrapezoidal(const Real& a_voltage) noexcept;

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
	@note For stationary simulations. 
      */
      virtual void
      computeBackgroundIonizationStationary() noexcept;
      /*!
	@brief Compute the detachment ionization rate for all voltages
	@note For stationary simulations. 
      */
      virtual void
      computeDetachmentStationary() noexcept;

      /*!
	@brief Compute field emission rates
	@note For stationary simulations. 
      */
      virtual void
      computeFieldEmissionStationary() noexcept;

      /*!
	@brief Compute field emission rates
	@param[out] a_emissionRate Field emission rate. 
	@param[in]  a_voltage Voltage multiplier
	@note For plotting the emission rate in transient simulations. 
      */
      virtual void
      computeFieldEmission(EBAMRCellData& a_emissionRate, const Real& a_voltage) const noexcept;

      /*!
	@brief Evaluate a function f = f(E, x) in a volume. 
	@param[out] a_data    Mesh data
	@param[in]  a_voltage Voltage multiplier
	@param[in]  a_func    Function to evaluate
	@note For plotting the emission rate in transient simulations. 
      */
      virtual void
      evaluateFunction(EBAMRCellData&                                             a_data,
                       const Real&                                                a_voltage,
                       const std::function<Real(const Real E, const RealVect x)>& a_func) const noexcept;

      /*!
	@brief Evaluate a function f = f(E,x) in a volume. 
	@param[out] a_data    Mesh data
	@param[in]  a_voltage Voltage multiplier
	@param[in]  a_func    Function to evaluate
	@param[in]  a_level   Grid level
	@note For plotting the emission rate in transient simulations. 
      */
      virtual void
      evaluateFunction(LevelData<EBCellFAB>&                                      a_data,
                       const Real&                                                a_voltage,
                       const std::function<Real(const Real E, const RealVect x)>& a_func,
                       const int                                                  a_level) const noexcept;

      /*!
	@brief Interpolate between K values to find voltage giving K_inception and store values in m_inceptionVoltage. 
	@details This is the volumetric version.
      */
      virtual void
      computeInceptionVoltageVolume() noexcept;

      /*!
	@brief Compute the minimum inception voltage and the starting electron position.
	@param[in] a_Uinc Inception voltage for all starting positions. 
      */
      virtual std::pair<Real, RealVect>
      computeMinimumInceptionVoltage(const EBAMRCellData& a_Uinc) const noexcept;

      /*!
	@brief Compute the critical volume of the K values for each voltage
	@note For the stationary method
      */
      virtual void
      computeCriticalVolumeStationary() noexcept;

      /*!
	@brief Compute the critical volume of the K values for each voltage
	@note Transient approach. 
      */
      virtual Real
      computeCriticalVolumeTransient() const noexcept;

      /*!
	@brief Compute the critical area of the K values for each voltage
	@note For the stationary method
      */
      virtual void
      computeCriticalAreaStationary() noexcept;

      /*!
	@brief Compute the critical area of the K values for each voltage
	@note For the transient method
      */
      virtual Real
      computeCriticalAreaTransient() const noexcept;

      /*!
	@brief Compute the ionization volume for each voltage
	@note For the stationary method
      */
      virtual void
      computeIonizationVolumeStationary() noexcept;

      /*!
	@brief Compute the ionization volume for each voltage
	@note For the transient method
      */
      virtual Real
      computeIonizationVolumeTransient(const Real& a_voltage) const noexcept;

      /*!
	@brief Compute time to first critical electron.
      */
      virtual void
      computeRdotStationary() noexcept;

      /*!
	@brief Print data report
	@note For transient runs
      */
      void
      writeDataReportTransient(const Real& a_voltage) const noexcept;

      /*!
	@brief Print report to the terminal. 
	@note For stationary runs
      */
      virtual void
      writeReportStationary() const noexcept;

      /*!
	@brief Print report to the terminal. 
	@note For transient runs. Called in the constructor. 
      */
      virtual void
      writeReportTransient() const noexcept;

      /*!
	@brief Write plot data for the 'stationary' mode
	@param[inout] a_output Output data holder.
	@param[inout] a_icomp  Starting component in a_output to begin at. 
	@param[in]    a_outputRealm Realm where a_output belongs
	@param[in]    a_level  Grid level
      */
      virtual void
      writePlotDataStationary(LevelData<EBCellFAB>& a_output,
                              int&                  a_icomp,
                              const std::string     a_outputRealm,
                              const int             a_level) const noexcept;

      /*!
	@brief Write plot data for the 'transient' mode
	@param[inout] a_output            Output data holder.
	@param[in]    a_plotVariableNames Plot variable names.
	@param[in]    a_outputRealm Realm where a_output belongs
	@param[inout] a_icomp             Starting component in a_output to begin at. 
      */
      virtual void
      writePlotDataTransient(LevelData<EBCellFAB>& a_output,
                             int&                  a_icomp,
                             const std::string     a_outputRealm,
                             const int             a_level) const noexcept;

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
	@brief Set the negative ion velocity. Note
	@param[in] a_voltage Current voltage
      */
      inline void
      computeIonVelocity(const Real& a_voltage) noexcept;

      /*!
	@brief Set the negative ion diffusion coefficient
	@param[in] a_voltage Current voltage
      */
      inline void
      computeIonDiffusion(const Real& a_voltage) noexcept;

      /*!
	@brief Parse simulation mode
      */
      inline void
      parseMode() noexcept;

      /*!
	@brief Parse voltage levels
      */
      inline void
      parseVoltages() noexcept;

      /*!
	@brief Parse output settings
      */
      inline void
      parseOutput() noexcept;

      /*!
	@brief Parse class verbosity
      */
      inline void
      parseVerbosity() noexcept;

      /*!
	@brief Parse the inception algorithm
      */
      inline void
      parseInceptionAlgorithm() noexcept;

      /*!
	@brief Parse the transport algorithm
      */
      inline void
      parseTransportAlgorithm() noexcept;

      /*!
	@brief Parse plot variables
      */
      inline void
      parsePlotVariables() noexcept;

      /*!
	@brief Calculate the total electric field = inhomogeneous + V * homogeneous
	@param[in] a_sumField Sum of contributions
	@param[in] a_inhomogeneousField Inhomogeneous field contribution
	@param[in] a_homogeneous Homogeneous field contribution
	@param[in] a_voltage Electrode voltage
      */
      inline void
      superposition(EBAMRCellData&       a_sumField,
                    const EBAMRCellData& a_inhomogeneousField,
                    const EBAMRCellData& a_homogeneousField,
                    const Real           a_voltage) const noexcept;

      /*!
	@brief Update m_electricField with the input voltage
	@param[in] a_voltage Electrode voltage
      */
      inline void
      superposition(const Real a_voltage) noexcept;

      /*!
	@brief Update the input field with the scaled voltage
	@param[out] a_sumField Sum of fields
	@param[in] a_voltage Electrode voltage
      */
      inline void
      superposition(EBAMRCellData& a_sumField, const Real a_voltage) const noexcept;

      /*!
	@brief Write data to output. Convenience function used for IO.
	@param[inout] a_output Output data holder.
	@param[inout] a_icomp  Starting component where this solver begins writing the output. 
	@param[in]    a_data   Data to write. 
	@param[in]    a_outputRealm Realm where a_output belongs
	@param[in]    a_level  Grid level
	@param[in]    a_interpToCentroids If true, a_data will be interpolated to cell centroids before writing to a_output. 
	@param[in]    a_interpGhost If true, interpolate ghost cells
      */
      virtual void
      writeData(LevelData<EBCellFAB>& a_output,
                int&                  a_comp,
                const EBAMRCellData&  a_data,
                const std::string     a_outputRealm,
                const int             a_level,
                const bool            a_interpToCentroids,
                const bool            a_interpGhost) const noexcept;
    };
  } // namespace StreamerInception
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_StreamerInceptionStepperImplem.H>

#endif
