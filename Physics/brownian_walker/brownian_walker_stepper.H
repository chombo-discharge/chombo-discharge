/*!
  @file   brownian_walker_stepper.H
  @brief  TimeStepper class for only solving the Brownian_Walker equation (with surface charge)
  @author Robert Marskar
  @data   March 2020
*/

#ifndef _BROWNIAN_WALKER_STEPPER_
#define _BROWNIAN_WALKER_STEPPER_

#include "CD_TimeStepper.H"
#include "ito_solver.H"
#include "ito_layout.H"

#include "CD_NamespaceHeader.H"
namespace physics {
  namespace brownian_walker {

    /*!
      @brief Class that just solves the brownian_walker equation
    */
    class brownian_walker_stepper : public TimeStepper {
    public:

      brownian_walker_stepper();
      brownian_walker_stepper(RefCountedPtr<ito_solver>& a_solver);
      ~brownian_walker_stepper();

      // Setup routines
      void initialData() override;
      void allocate() override;
      void postInitialize() override;

      // Load balancing
      bool loadBalanceThisRealm(const std::string a_realm) const override;
      void loadBalanceBoxes(Vector<Vector<int> >&            a_procs,
			      Vector<Vector<Box> >&            a_boxes,
			      const std::string                a_realm,
			      const Vector<DisjointBoxLayout>& a_grids,
			      const int                        a_lmin,
			      const int                        a_finestLevel) override;

      // IO routines
      void writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override;
      void readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
      void postCheckpointSetup() override;
      int getNumberOfPlotVariables() const override;
      void writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const override;

      // Advance routines
      void computeDt(Real& a_dt, TimeCode& a_timeCode) override;

      void synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;
      void printStepReport() override;
      
      // Registration routines
      void registerOperators() override;
      void registerRealms() override;

      // Runtime options
      void parseRuntimeOptions() override;

      // Regrid routines
      bool needToRegrid() override;

      void preRegrid(const int a_lbase, const int a_oldFinestLevel) override;
      void deallocate() override;

      // Pure routines for derived classes
      void setup_solvers() override;
      Real advance(const Real a_dt) override;
      void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;
      void postRegrid() override;

    protected:

      // Realm
      std::string m_realm; 

      // Phase
      phase::which_phase m_phase;

      // Load balance or not
      bool m_LoadBalancing;

      // Ito solver & species
      RefCountedPtr<ito_solver>  m_solver;
      RefCountedPtr<ito_species> m_species;

      // Diffusion coefficient and velocity
      Real m_faceCenteredDiffusionCoefficient;
      Real m_omega;

      // Particles per cell
      int m_ppc;

      // Maximum hops and stuff like that
      Real m_max_cells_hop;

      void setVelocity();
      void setVelocity(const int a_level);
    };
  }
}
#include "CD_NamespaceFooter.H"
#endif
