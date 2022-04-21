/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TracerParticleStepper.H
  @brief  Declaration of a TimeStepper implementation for advancing a tracer particle solver. 
  @author Robert Marskar
*/

#ifndef CD_TracerParticleStepper_H
#define CD_TracerParticleStepper_H

// Our includes
#include <CD_TracerParticle.H>
#include <CD_TracerParticleSolver.H>
#include <CD_TracerParticlePhysics.H>
#include <CD_TimeStepper.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace TracerParticle {

    /*!
      @brief Various supported integration algorithms. 
    */
    enum class IntegrationAlgorithm {
      Euler,
      RK2,
      RK4
    };

    /*!
      @brief Enum for handling various velocity fields.
    */
    enum class VelocityField {
      Diagonal,
      Rotational
    };

    /*!
      @brief Implementation of TimeStepper for advancing tracer particles in a fixed velocity field. 
      @details The template requirements on the particle type P are the same as for TracerParticleSolver with the addition 
      of templated functions of the type RealVect& P::vector<size_t>() which returns a RealVect. This is used for the Runge-Kutta
      particle advection kernels. This class will work fine with TracerParticle<0, M> where M >= 4. 
    */
    template <typename P>
    class TracerParticleStepper : public TimeStepper {
    public:

      /*!
	@brief Constructor. Does nothing. 
      */
      TracerParticleStepper();

      /*!
	@brief Copy constructor. Disallowed.
	@param[in] a_other Other
      */
      TracerParticleStepper(const TracerParticleStepper<P>& a_other) = delete;

      /*!
	@brief Move constructor. Disallowed.
	@param[in] a_other Other
      */
      TracerParticleStepper(const TracerParticleStepper<P>&& a_other) = delete;      

      /*!
	@brief Destructor. 
      */
      virtual ~TracerParticleStepper();

      /*!
	@brief Copy assignment operator. Disallowed.
	@param[in] a_other Other
      */
      TracerParticleStepper& operator=(const TracerParticleStepper<P>& a_other) = delete;

      /*!
	@brief Move assignment operator. Disallowed.
	@param[in] a_other Other
      */
      TracerParticleStepper& operator=(const TracerParticleStepper<P>&& a_other) = delete;      

      /*!
	@brief Instantiate the tracer particle solver.
      */
      void setupSolvers() override;

      /*!
	@brief Allocate storage for solvers and time stepper. 
      */
      void allocate() override;

      /*!
	@brief Fill problem with initial data. 
      */
      void initialData() override;

      /*!
	@brief Perform any post-initialization steps. Not needed for this time stepper. 
      */
      void postInitialize() override { }

      /*!
	@brief Post checkpoint operations. Not needed for this time stepper. 
      */
      void postCheckpointSetup() override { }

      /*!
	@brief Register realms. Primal is the only realm we need. 
      */
      void registerRealms() override;

      /*!
	@brief Register operators. 
	@details Only need the solver's operators.
      */
      void registerOperators() override;

      /*!
	@brief Parse options
      */
      void parseOptions();

      /*!
	@brief Parse runtime options
      */
      void parseRuntimeOptions() override;

#ifdef CH_USE_HDF5
      /*!  
	@brief Write checkpoint data to file
	@param[inout] a_handle HDF5 file
	@param[in]    a_lvl    Grid level
      */
      virtual void writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
#endif

#ifdef CH_USE_HDF5
      /*!
	@brief Read checkpoint data from file
	@param[inout] a_handle HDF5 file
	@param[in]    a_lvl    Grid level
      */
      virtual void readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
#endif

      /*!
	@brief Write plot data to output holder. 
	@param[inout] a_output            Output data holder.
	@param[in]    a_plotVariableNames Plot variable names.
	@param[inout] a_icomp             Starting component in a_output to begin at. 
      */
      virtual void writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const override;

      /*!
	@brief Get the number of plot variables for this time stepper. 
	@details This is necessary because Driver, not TimeStepper, is responsible for allocating the necessary memory. 
	@return Returns number of plot variables that will be written during writePlotData
      */
      virtual int getNumberOfPlotVariables() const override;

      /*!
	@brief Compute a time step to be used by Driver. 
	@param[out] a_dt       Time step to be used.
	@param[out] a_timeCode Identifier for how the time step was computed. 
      */
      virtual void computeDt(Real& a_dt, TimeCode& a_timeCode) override;

      /*!
	@brief Advancement method. The implementation of this method should advance all equations of motion
	@param[in] a_dt Time step to be used for advancement
	@return    Returns the time step that was used. 
	@note The return value does not need to equal a_dt. Adaptive time stepping methods will generally return != a_dt.
      */
      virtual Real advance(const Real a_dt) override;

      /*!
	@brief Synchronzie solver times and time steps
	@param[in] a_step Time step
	@param[in] a_time Time (in seconds)
	@param[in] a_dt   Time step that was used. 
      */
      virtual void synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;

      /*!
	@brief Print a step report. 
	@details This is called by Driver after time step. The routine can be used to display use information about the simulation progress. 
      */
      virtual void printStepReport() override;

      /*!
	@brief Perform pre-regrid operations.
	@details This should include all copying all data which should be interpolated to the new grids. It can also include deallocating memory in case the regrid
	operation takes a lot of memory. 
	@param[in] a_lmin           The coarsest level that changes
	@param[in] a_oldFinestLevel The finest level before the regrid. 
      */
      virtual void preRegrid(const int a_lmin, const int a_oldFinestLevel) override;

      /*!
	@brief Time stepper regrid method. 
	@param[in] a_lmin           The coarsest level that changed. 
	@param[in] a_oldFinestLevel The finest level before the regrid. 
	@param[in] a_newFinestLevel The finest level after the regrid. 
      */
      virtual void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

      /*!
	@brief Perform post-regrid operations. 
	@details This includes all operations to be done AFTER interpolating data to new grids. 
      */
      virtual void postRegrid() override;

    protected:

      /*!
	@brief Integration algorithm
      */
      IntegrationAlgorithm m_algorithm;

      /*!
	@brief Velocity field type
      */
      VelocityField m_velocityField;

      /*!
	@brief Realm where solver and m_velocity lives
      */
      std::string m_realm;

      /*!
	@brief Phase (gas/solid) where the solver lives
      */
      phase::which_phase m_phase;

      /*!
	@brief Tracer particle solvers
      */
      RefCountedPtr<TracerParticleSolver<P> > m_solver;

      /*!
	@brief Tracer particle velocity field
      */
      EBAMRCellData m_velocity;

      /*!
	@brief Particle "CFL" time step
      */
      Real m_cfl;

      /*!
	@brief Particle blob radius
      */
      Real m_blobRadius;

      /*!
	@brief Particle blob center
      */
      RealVect m_blobCenter;

      /*!
	@brief Number of particles
      */
      size_t m_numInitialParticles;

      /*!
	@brief Fill initial particles
      */
      virtual void initialParticles();

      /*!
	@brief Set the velocity on the mesh.
      */
      virtual void setVelocity();

      /*!
	@brief Parse reactive integrator.
      */
      virtual void parseIntegrator();

      /*!
	@brief Parse velocity field.
      */
      virtual void parseVelocityField();

      /*!
	@brief Parse initial conditions. 
      */
      virtual void parseInitialConditions();      

      /*!
	@brief Advance particles using explicit Euler rule
	@param[in] a_dt Advanced time step
      */
      virtual void advanceParticlesEuler(const Real a_dt);

      /*!
	@brief Advance particles using second order Runge-Kutta
	@param[in] a_dt Advanced time step
      */
      virtual void advanceParticlesRK2(const Real a_dt);

      /*!
	@brief Advance particles using fourth order Runge-Kutta
	@param[in] a_dt Advanced time step
      */
      virtual void advanceParticlesRK4(const Real a_dt);            
    };
  }
}

#include <CD_NamespaceFooter.H>

#include <CD_TracerParticleStepperImplem.H>

#endif
