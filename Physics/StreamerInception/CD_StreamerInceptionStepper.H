/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_StreamerInceptionStepper.H
  @brief  TimeStepper class for evaluating the streamer inception criterion.
  @author Robert Marskar
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
      Simpsons,
      Simpsons38,
      Booles
    };

    /*!
      @brief Class for volume evaluation of the streamer inception integral. 
      @details P is the tracer particle type and F is the field solver type. 
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
	@brief Distribute particles in every grid cell
      */
      void
      initialParticles();

      /*!
	@brief Perform any post-initialization steps. 
      */
      void
      postInitialize() override;

      /*!
	@brief Solve streamer inception integral for each particle in each voltage and store K values in m_Kvalues. Called in postInitialize().
      */
      void
      solveKValues();
      
      /*!
	@brief K integral: Move particles as long as integration isn't stopped
	@param[in] amrParticles  Particles of m_particleSolver
      */
      void
      moveParticlesTrapezoidal(ParticleContainer<P>& amrParticles, const Real& voltage); 

      /*!
	@brief K integral: Add integration parts after particles move
	@param[in] amrParticles  Particlexs of m_particleSolver
      */
      void
      addIntegrationTrapezoidal(ParticleContainer<P>& amrParticles, const Real& voltage);

      /*!
	@brief K integral: Move particles as long as integration isn't stopped
	@param[in] amrParticles  Particles of m_particleSolver
      */
      void
      moveParticlesSimpsons(ParticleContainer<P>& amrParticles, const Real& voltage);

      /*!
	@brief K integral: Add integration parts after particles move
	@param[in] amrParticles  Particles of m_particleSolver
      */
      void
      addIntegrationSimpsons(ParticleContainer<P>& amrParticles, const Real& voltage);

      /*!
	@brief 
	@param[in] amrParticles  Particles of m_particleSolver
	@param[out] true if integration shall continue, false if not
      */
      bool
      keepGoing(ParticleContainer<P>& amrParticles) const;

      /*!
	@brief Move particles back to their original position after integration finished
	@param[in] amrParticles  Particles of m_particleSolver
      */
      void
      moveParticlesInitialPosition(ParticleContainer<P>& amrParticles);

      /*!
	@brief Set all particles' masses to 0
	@param[in] amrParticles  Particles of m_particleSolver
      */
      void
      resetParticleMass(ParticleContainer<P>& amrParticles);

      /*!
	@brief Calculate particles' drift times until inside electrode or outside grid
      */
      void
      solveDriftTimes();

      /*!
	@brief Drift time integral (trapezoidal method): Add time step to particles' drift time tracker
	@param[in] amrParticles  Particles of m_particleSolver
	@param[in] dxSign  Sign of dx (0=positive, 1=negative)
     */
      void
      moveParticlesDriftTimeTrapezoidal(ParticleContainer<P>& amrParticles, const int dxSign);

      /*!
	@brief Drift time integral (trapezoidal method): Add integration parts after particles move
	@param[in] amrParticles  Particles of m_particleSolver
      */
      void
      addIntegrationDriftTimeTrapezoidal(ParticleContainer<P>& amrParticles);

      /*!
	@brief Drift time integral (Simpson's method): Add time step to particles' drift time tracker
	@param[in] amrParticles  Particles of m_particleSolver
	@param[in] dxSign  Sign of dx (0=positive, 1=negative)
      */
      void
      moveParticlesDriftTimeSimpsons(ParticleContainer<P>& amrParticles, const int dxSign);

      /*!
	@brief Drift time integral (Simpson's method): Add integration parts after particles move
	@param[in] amrParticles  Particles of m_particleSolver
      */
      void
      addIntegrationDriftTimeSimpsons(ParticleContainer<P>& amrParticles);      

      /*!
	@brief Post checkpoint operations. Not needed for this time stepper. 
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
	@brief Set the effective ionization coefficient.
      */
      virtual void
      setAlpha(const std::function<Real(const Real& E)>& a_alpha) noexcept;

      /*!
	@brief Set voltage for field solver. 
	@param[in] a_voltage Input voltage. 
      */
      virtual void
      setVoltage(const Real a_voltage) noexcept;

      /*!
	@brief 

      /*!
	@brief Check if particle is outside grid boundaries
	@param[in] pos  Particle position
	@param[in] probLo Lower left boundary of grid
	@param[in] probHi Upper right boundary of grid
	@param[out]  true if particle outside boundary, false elsewise
      */
      bool
      particleOutsideGrid(const RealVect pos, const RealVect probLo, const RealVect probHi);

      /*!
	@brief Check if particle is inside electrode
	@param[in] pos  Particle position
	@param[out] true if particle inside electrode, false elsewise
      */
      bool
      particleInsideElectrode(const RealVect pos);

    protected:
      /*!
	@brief Integration algorithm
      */
      IntegrationAlgorithm m_integrationAlgorithm;

      /*!
	@brief Tracer particle solver
      */
      RefCountedPtr<TracerParticleSolver<P>> m_particleSolver;

      /*!
	@brief Field solver
      */
      RefCountedPtr<FieldSolver> m_fieldSolver;

      /*!
	@brief Velocity field for particles. 
      */
      EBAMRCellData m_velocity;

      /*!
	@brief K-values -- > should have input-specified number of components
      */
      EBAMRCellData m_Kvalues;

      /*!
	@brief Drift time
      */
      EBAMRCellData m_driftTimes;

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
	@brief Voltage sweeps
      */
      std::vector<Real> m_voltageSweeps;

      /*!
	@brief Scale electrostatic solution.
	@param[in] a_scale Scaling factor for potential and electric field. 
      */
      virtual void
      scaleElectrostatics(const Real& a_scale) noexcept;
    };
  } // namespace StreamerInception
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_StreamerInceptionStepperImplem.H>

#endif
