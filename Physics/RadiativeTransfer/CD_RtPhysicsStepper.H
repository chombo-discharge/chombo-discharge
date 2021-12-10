/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RtPhysicsStepper.H
  @brief  TimeStepper class for only evolving radiative transfer modules
  @author Robert Marskar
*/

#ifndef CD_RtPhysicsStepper_H
#define CD_RtPhysicsStepper_H

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
    class RtPhysicsStepper : public TimeStepper {
    public:

      /*!
	@brief Constructor. Reads a couple of input options.
      */
      RtPhysicsStepper();

      /*!
	@brief Destructor (does nothing)
      */
      ~RtPhysicsStepper();

      // Stuff below here is derived from TimeStepper.
      // ----------------------------------------------
      void setupSolvers() override;
      void initialData() override;
      void allocate() override;
      void postInitialize() override;

      // IO routines
#ifdef CH_USE_HDF5
      void writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
#endif
#ifdef CH_USE_HDF5      
      void readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
#endif
      void postCheckpointSetup() override;
      int getNumberOfPlotVariables() const override;
      void writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const override;

      // Advance routines
      void computeDt(Real& a_dt, TimeCode& a_timeCode) override;
      Real advance(const Real a_dt) override;
      void synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;
      void printStepReport() override;

      // Runtime options
      void parseRuntimeOptions() override;

      // Registration routines
      void registerRealms() override;
      void registerOperators() override;

      // Regrid routines
      bool needToRegrid() override;
      void preRegrid(const int a_lmin, const int a_oldFinestLevel) override;
      void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;
      void postRegrid() override;

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
	@brief Reference to solver
      */
      RefCountedPtr<RtSolver>  m_solver;

      /*!
	@brief Reference to species
      */      
      RefCountedPtr<RtSpecies> m_species;

      /*!
	@brief For setting a Gaussian source in the radiative transfer equation solver
      */
      void setGaussianSource();
    };
  }
}

#include <CD_NamespaceFooter.H>

#include <CD_RtPhysicsStepperImplem.H>

#endif
