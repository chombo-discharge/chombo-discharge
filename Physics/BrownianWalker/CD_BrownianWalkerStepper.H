/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_BrownianWalkerStepper.H
  @brief  TimeStepper class for only solving the Brownian_Walker equation (with surface charge)
  @author Robert Marskar
*/

#ifndef CD_BrownianWalkerStepper_H
#define CD_BrownianWalkerStepper_H

// Our includes
#include "CD_TimeStepper.H"
#include <CD_ItoSolver.H>
#include <CD_ItoLayout.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace BrownianWalker {

    /*!
      @brief Class that just solves the BrownianWalker equation
    */
    class BrownianWalkerStepper : public TimeStepper {
    public:

      BrownianWalkerStepper();
      BrownianWalkerStepper(RefCountedPtr<ItoSolver>& a_solver);
      ~BrownianWalkerStepper();

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

      // Pure routines for derived classes
      void setupSolvers() override;
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
      RefCountedPtr<ItoSolver>  m_solver;
      RefCountedPtr<ItoSpecies> m_species;

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

#include <CD_NamespaceFooter.H>

#endif
