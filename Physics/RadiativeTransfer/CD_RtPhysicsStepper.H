/*!
  @file   RtPhysicsStepper.H
  @brief  TimeStepper class for only solving the Rte equation (with surface charge)
  @author Robert Marskar
  @date   May 2020
*/

#ifndef _RTE_STEPPER_
#define _RTE_STEPPER_

#include "CD_TimeStepper.H"
#include <CD_RtSolver.H>

#include "CD_NamespaceHeader.H"
namespace physics {
  namespace rte {

    /*!
      @brief Class that just solves the rte equation using the RtSolver advance method
    */
    template <class T>
    class RtPhysicsStepper : public TimeStepper {
    public:
      RtPhysicsStepper();
      ~RtPhysicsStepper();

      // Stuff below here is derived from TimeStepper.
      // ----------------------------------------------
      void setup_solvers() override;
      void initialData() override;
      void allocate() override;
      void postInitialize() override;

      // IO routines
      void writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
      void readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
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

      std::string m_realm;

      // Solver and species
      RefCountedPtr<RtSolver>  m_solver;
      RefCountedPtr<RtSpecies> m_species;

      // Stuff for source term
      Real     m_blob_radius;    // Gaussian radius
      Real     m_blob_amplitude; // Gaussian amplitude
      RealVect m_blob_center;    // Center

      void set_gaussian_source();
    };
  }
}
#include "CD_NamespaceFooter.H"

#include <CD_RtPhysicsStepperImplem.H>

#endif
