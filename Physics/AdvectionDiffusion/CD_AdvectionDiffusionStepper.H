/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_AdvectionDiffusionStepper.H
  @brief  TimeStepper class for only solving the the advection-diffusion equation
  @author Robert Marskar
*/

#ifndef CD_AdvectionDiffusionStepper
#define CD_AdvectionDiffusionStepper

// Our includes
#include <CD_TimeStepper.H>
#include <CD_CdrSolver.H>
#include <CD_NamespaceHeader.H>

namespace Physics {

  /*!
    @brief Physics functionality for running scalar advection-diffusion code
  */
  namespace AdvectionDiffusion {

    /*!
      @brief Implementation of TimeStepper which runs a single advection-diffusion-reaction solver.
    */
    class AdvectionDiffusionStepper : public TimeStepper {
    public:

      AdvectionDiffusionStepper();
      AdvectionDiffusionStepper(RefCountedPtr<CdrSolver>& a_solver);
      ~AdvectionDiffusionStepper();

      // Setup routines
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

      void synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;
      void printStepReport() override;

      //
      void parseRuntimeOptions() override;

      // Registration routines
      void registerRealms() override;
      void registerOperators() override;

      // Regrid routines
      bool needToRegrid() override;

      void preRegrid(const int a_lbase, const int a_finestLevel) override;

      // Pure routines for derived classes
      void setupSolvers() override;
      Real advance(const Real a_dt) override;
      void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;
      void postRegrid() override;
      

    protected:

      std::string m_realm;

      phase::which_phase m_phase;

      // CDR solver & species
      RefCountedPtr<CdrSolver>  m_solver;
      RefCountedPtr<CdrSpecies> m_species;

      // Diffusion coefficient and velocity
      Real m_cfl;
      Real m_faceCenteredDiffusionCoefficient;
      Real m_omega;

      // Use stochastic diffusion
      bool m_fhd;

      // Integrator (yes, I'm too lazy to write an enum for holding this)
      int m_integrator;

      void setVelocity();
      void setVelocity(const int a_level);

    private:

      EBAMRCellData m_tmp;
      EBAMRCellData m_k1;
      EBAMRCellData m_k2;
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
