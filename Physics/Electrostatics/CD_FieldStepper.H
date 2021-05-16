/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_FieldStepper.H
  @brief  TimeStepper class for only solving the Poisson equation (with surface charge)
  @author Robert Marskar
*/

#ifndef CD_FieldStepper_H
#define CD_FieldStepper_H

// Our includes
#include <CD_TimeStepper.H>
#include <sigma_solver.H>
#include <CD_FieldSolver.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace Electrostatics {

    /*!
      @brief Class that just solves the poisson equation
    */
    template <class T>
    class FieldStepper : public TimeStepper {
    public:

      static Real s_constant_one(const RealVect a_pos);
      static Real s_potential_one(const Real a_time);

      FieldStepper();
      ~FieldStepper();

      // Stuff below here is derived from TimeStepper.
      // ----------------------------------------------
      void setup_solvers() override;
      void initialData() override;
      void postInitialize() override;
      void allocate() override;

      // IO routines
      void writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const override{}
      void readCheckpointData(HDF5Handle& a_handle, const int a_lvl) override;
      void postCheckpointSetup() override {}
      int getNumberOfPlotVariables() const override;
      void writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const override;

      // Advance routines
      void computeDt(Real& a_dt, TimeCode& a_timeCode) override {a_dt = 0.0;}
      Real advance(const Real a_dt) override;
      void synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) override;
      void printStepReport() override {};

      // Registration routines
      void registerRealms() override;
      void registerOperators() override;

      // Regrid routines
      bool needToRegrid() override {return false;}
      void preRegrid(const int a_lbase, const int a_finestLevel) override;
      void deallocate() override;
      bool loadBalanceThisRealm(const std::string a_Realm) const override;
      void loadBalanceBoxes(Vector<Vector<int> >&            a_procs,
			      Vector<Vector<Box> >&            a_boxes,
			      const std::string                a_Realm,
			      const Vector<DisjointBoxLayout>& a_grids,
			      const int                        a_lmin,
			      const int                        a_finestLevel) override;
      void regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;
      void postRegrid() override;

      // Solve routine
      void solve_poisson();

    protected:

      bool m_LoadBalancing;
      BoxSorting m_boxSort;

      RefCountedPtr<FieldSolver> m_fieldSolver;
      RefCountedPtr<sigma_solver>   m_sigma;

      std::string m_Realm;

      Real m_init_rho;
      Real m_init_sigma;
    };
  }
}

#include <CD_NamespaceFooter.H>

#include <CD_FieldStepperImplem.H>

#endif
