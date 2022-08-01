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
      Trapezoidal,
      Simpsons
    };

    /*!
      @brief Class for streamer inception integral evaluations. 
      @details P is the tracer particle type and F is the field solver type. Note that the integration algorithm requires some transient storage
      on each particle (such as the previous position). This is defined in the TracerParticle template arguments. 
    */
    template <typename P = TracerParticle<4, 1>, typename F = FieldSolverMultigrid>
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
	@details This performs the inception integral evaluation and all that. 
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
	@brief Set the electron detachment rate
	@param[in] a_detachmentRate Detachment rate (units of 1/s). 
      */
      virtual void
      setDetachmentRate(const std::function<Real(const Real& E)>& a_detachmentRate) noexcept;

    protected:
      /*!
	@brief Integration algorithm.
      */
      IntegrationAlgorithm m_integrationAlgorithm;

      /*!
	@brief Tracer particle solver.
      */
      RefCountedPtr<TracerParticleSolver<P>> m_particleSolver;

      /*!
	@brief Field solver.
      */
      RefCountedPtr<FieldSolver> m_fieldSolver;

      /*!
	@brief Velocity field for particles. 
      */
      EBAMRCellData m_velocity;

      /*!
	@brief Inception integral values (K-values)
      */
      EBAMRCellData m_inceptionIntegral;

      /*!
	@brief Drift times (t*µV)
      */
      EBAMRCellData m_driftTimes;

      /*!
	@brief Computed inception voltages. 
      */
      EBAMRCellData m_inceptionVoltage;

      /*!
	@brief Inception criteria (read from input)
      */
      Real m_inceptionK;

      /*! 
	@brief Max K value for every voltage
      */
      std::vector<Real> m_maxKValues;

      /*!
	@brief Critical volumes of K values for each voltage
      */
      std::vector<Real> m_criticalVolumes;

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
	@brief Electron detachment rate
      */
      std::function<Real(const Real& E)> m_detachmentRate;

      /*!
	@brief Voltage sweeps
      */
      std::vector<Real> m_voltageSweeps;

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
	@brief K integral: Move particles as long as integration isn't stopped
	@param[in] a_voltage Voltage multiplier
      */
      virtual void
      K_moveParticlesTrapezoidal(const Real& a_voltage) noexcept;

      /*!
	@brief K integral: Add integration parts after particles move
	@param[in] a_voltage Voltage multiplier
      */
      virtual void
      K_addIntegrationTrapezoidal(const Real& a_voltage) noexcept;

      /*!
	@brief K integral: Move particles as long as integration isn't stopped
	@param[in] a_particles Particles of m_particleSolver
	@param[in] a_voltage   Voltage multiplier
      */
      virtual void
      K_moveParticlesSimpsons(ParticleContainer<P>& a_particles, const Real& a_voltage) noexcept;

      /*!
	@brief K integral: Add integration parts after particles move
	@param[in] a_particles Particles of m_particleSolver
	@param[in] a_voltage   Voltage multiplier
      */
      virtual void
      K_addIntegrationSimpsons(ParticleContainer<P>& a_particles, const Real& voltage) noexcept;

      /*!
	@brief K integral: Move particles one step according to the integration algorithm
	@param[in] a_particles Particles to move. 
	@param[in] a_voltage   Voltage multiplier
	@param[in] a_step      Step number for Simpsons' (doesn't affect Trapezoidal)
      */
      virtual void
      K_moveParticlesOneStep(ParticleContainer<P>& a_particles, const Real& a_voltage, const int& a_step) noexcept;

      /*!
	@brief 
	@param[in] a_particles Particles to evaluate. 
	@return true if integration shall continue, false if not
      */
      virtual bool
      keepIntegrating(ParticleContainer<P>& a_particles) const noexcept;

      /*!
	@brief Move particles back to their original position.
	@param[in] a_particles Particles of m_particleSolver
      */
      virtual void
      rewindParticles(ParticleContainer<P>& a_particles) noexcept;

      /*!
	@brief Set all particles' masses to 0
	@param[in] a_particles  Particles of m_particleSolver
      */
      virtual void
      resetParticleMass(ParticleContainer<P>& a_particles) noexcept;

      /*!
	@brief Calculate particles' drift times until inside electrode or outside grid. 
      */
      virtual void
      computeDriftTimes() noexcept;

      /*!
	@brief Drift time integral (trapezoidal method): Add time step to particles' drift time tracker
	@param[in] a_particles  Particles of m_particleSolver
	@param[in] a_dxSign     Sign of dx (0=positive, 1=negative)
      */
      virtual void
      driftTime_moveParticlesTrapezoidal(ParticleContainer<P>& a_particles, const int& a_dxSign) noexcept;

      /*!
	@brief Drift time integral (trapezoidal method): Add integration parts after particles move
	@param[in] a_particles Particles of m_particleSolver
      */
      virtual void
      driftTime_addIntegrationTrapezoidal(ParticleContainer<P>& a_particles) noexcept;

      /*!
	@brief Drift time integral (Simpson's method): Add time step to particles' drift time tracker
	@param[in] a_particles Particles of m_particleSolver
	@param[in] a_dxSign    Sign of dx (0=positive, 1=negative)
      */
      virtual void
      driftTime_moveParticlesSimpsons(ParticleContainer<P>& a_particles, const int& dxSign) noexcept;

      /*!
	@brief Drift time integral (Simpson's method): Add integration parts after particles move
	@param[in] a_particles  Particles of m_particleSolver
      */
      virtual void
      driftTime_addIntegrationSimpsons(ParticleContainer<P>& a_particles) noexcept;

      /*!
	@brief Drift time integral: Move particles one step according to the integration algorithm
	@param[in] a_particles Particles of m_particleSolver
	@param[in] a_dxSign    Sign of dx (0=positive, 1=negative)
	@param[in] a_step      Step number for Simpsons' (doesn't affect Trapezoidal)
      */
      virtual void
      driftTime_moveParticlesOneStep(ParticleContainer<P>& a_particles, const int& a_dxSign, const int& a_step) noexcept;

      /*!
	@brief Interpolate between K values to find voltage giving K_inception and store values in m_inceptionVoltage
      */
      virtual void
      computeInceptionVoltage() noexcept;

      /*!
	@brief Compute the critical volume of the K values for each voltage
      */
      virtual void
      computeCriticalVolume() noexcept;

      /*!
	@brief Print m_voltageSweeps, m_maxKValues, and m_criticalVolumes to terminal
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
      particleInsideElectrode(const RealVect& a_pos) const noexcept;
    };
  } // namespace StreamerInception
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_StreamerInceptionStepperImplem.H>

#endif
