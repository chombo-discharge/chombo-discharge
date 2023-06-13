/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_FieldStepperImplem.H
  @brief  Implementation of CD_FieldStepper.H
  @author Robert Marskar
*/

#ifndef CD_FieldStepperImplem_H
#define CD_FieldStepperImplem_H

// Std includes
#include <math.h>

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>

// Our includes
#include <CD_FieldStepper.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace Electrostatics {

    template <class T>
    FieldStepper<T>::FieldStepper()
    {
      CH_TIME("FieldStepper::FieldStepper");

      m_verbosity = -1;

      ParmParse pp("FieldStepper");

      Real     rho0       = 0.0;
      Real     sigma0     = 0.0;
      Real     blobRadius = 1.0;
      RealVect blobCenter = RealVect::Zero;

      std::string  str;
      Vector<Real> vec(SpaceDim);

      // Read in parameters
      pp.get("init_rho", rho0);
      pp.get("init_sigma", sigma0);
      pp.get("rho_radius", blobRadius);
      pp.get("load_balance", m_loadBalance);
      pp.get("box_sorting", str);
      pp.get("realm", m_realm);
      pp.get("verbosity", m_verbosity);
      pp.getarr("rho_center", vec, 0, SpaceDim);

      blobCenter = RealVect(D_DECL(vec[0], vec[1], vec[2]));

      m_rhoGas = [a = rho0, c = blobCenter, r = blobRadius](const RealVect x) -> Real {
        const Real d = (x - c).dotProduct(x - c);
        return a * exp(-d / (2 * r * r));
      };

      m_rhoDielectric        = m_rhoGas;
      m_surfaceChargeDensity = [s = sigma0](const RealVect& x) -> Real {
        return s;
      };

      if (str == "none") {
        m_boxSort = BoxSorting::None;
      }
      if (str == "std") {
        m_boxSort = BoxSorting::Std;
      }
      else if (str == "shuffle") {
        m_boxSort = BoxSorting::Shuffle;
      }
      else if (str == "morton") {
        m_boxSort = BoxSorting::Morton;
      }
      else {
        MayDay::Error("FieldStepper::FieldStepper - unknown box sorting method requested for argument 'BoxSorting'");
      }
    }

    template <class T>
    FieldStepper<T>::~FieldStepper()
    {
      CH_TIME("FieldStepper::~FieldStepper");
    }

    template <class T>
    void
    FieldStepper<T>::setupSolvers()
    {
      CH_TIME("FieldStepper::setupSolvers");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::setupSolvers" << endl;
      }

      // Define a voltage function to be used in the simulation.
      auto voltage = [](const Real a_time) -> Real {
        return 1.0;
      };

      // Instantiate the FieldSolver and set it up.
      m_fieldSolver = RefCountedPtr<FieldSolver>(new T());
      m_fieldSolver->parseOptions();
      m_fieldSolver->setAmr(m_amr);
      m_fieldSolver->setComputationalGeometry(m_computationalGeometry);
      m_fieldSolver->setVoltage(voltage);
      m_fieldSolver->setRealm(m_realm);
      m_fieldSolver->setTime(0, 0.0, 0.0);
      m_fieldSolver->setVerbosity(m_verbosity);

      // Instantiate the surface charge solver and set it up.
      m_sigma = RefCountedPtr<SurfaceODESolver<1>>(new SurfaceODESolver<1>(m_amr));
      m_sigma->setPhase(phase::gas);
      m_sigma->setRealm(m_realm);
      m_sigma->setName("Surface charge");
      m_sigma->parseOptions();
      m_sigma->setTime(0, 0.0, 0.0);
      m_sigma->setVerbosity(m_verbosity);
    }

    template <class T>
    void
    FieldStepper<T>::registerOperators()
    {
      CH_TIME("FieldStepper::registerOperators");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::registerOperators" << endl;
      }

      m_fieldSolver->registerOperators();
      m_sigma->registerOperators();
    }

    template <class T>
    void
    FieldStepper<T>::registerRealms()
    {
      CH_TIME("FieldStepper::registerRealms");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::registerRealms" << endl;
      }

      m_amr->registerRealm(m_realm);
    }

    template <class T>
    void
    FieldStepper<T>::allocate()
    {
      CH_TIME("FieldStepper::allocate");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::allocate" << endl;
      }

      m_fieldSolver->allocate();
      m_sigma->allocate();
    }

    template <class T>
    MFAMRCellData&
    FieldStepper<T>::getPotential()
    {
      CH_TIME("FieldStepper::getPotential");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::getPotential" << endl;
      }

      return m_fieldSolver->getPotential();
    }

    template <class T>
    void
    FieldStepper<T>::initialData()
    {
      CH_TIME("FieldStepper::initialData");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::initialData" << endl;
      }
    }

    template <class T>
    void
    FieldStepper<T>::solvePoisson()
    {
      CH_TIME("FieldStepper::solvePoisson");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::solvePoisson" << endl;
      }

      // Solve using rho from fieldsolver and sigma from SurfaceODESolver.
      MFAMRCellData& phi   = m_fieldSolver->getPotential();
      MFAMRCellData& rho   = m_fieldSolver->getRho();
      EBAMRIVData&   sigma = m_sigma->getPhi();

      // Solve and issue error message if we did not converge.
      const bool converged = m_fieldSolver->solve(phi, rho, sigma);

      if (!converged) {
        MayDay::Warning("FieldStepper<T>::solvePoisson - did not converge");
      }

      // Compute the electric field.
      m_fieldSolver->computeElectricField();
    }

    template <class T>
    void
    FieldStepper<T>::postInitialize()
    {
      CH_TIME("FieldStepper::postInitialize");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::postInitialize" << endl;
      }

      // Fetch the potential and surface/space charges and initialize them. The potential is set to zero
      // and the space/surface charge set from m_rhoGas, m_rhoDielectric, and m_surfaceChargeDensity
      MFAMRCellData& state = m_fieldSolver->getPotential();
      EBAMRIVData&   sigma = m_sigma->getPhi();

      DataOps::setValue(state, 0.0);
      DataOps::setValue(sigma, m_surfaceChargeDensity, m_amr->getProbLo(), m_amr->getDx(), 0);
      m_sigma->resetElectrodes(sigma, 0.0);
      m_fieldSolver->setRho(m_rhoGas);

      // Solve Poisson equation.
      this->solvePoisson();
    }

    template <class T>
    Real
    FieldStepper<T>::advance(const Real a_dt)
    {
      CH_TIME("FieldStepper::advance");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::advance" << endl;
      }

      MayDay::Error("FieldStepper<T>::advance - callling this is an error. Please set Driver.max_steps = 0");

      return std::numeric_limits<Real>::max();
    }

#ifdef CH_USE_HDF5
    template <class T>
    void
    FieldStepper<T>::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const
    {
      CH_TIME("FieldStepper::writeCheckpointData");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::writeCheckpointData" << endl;
      }
    }
#endif

#ifdef CH_USE_HDF5
    template <class T>
    void
    FieldStepper<T>::readCheckpointData(HDF5Handle& a_handle, const int a_lvl)
    {
      CH_TIME("FieldStepper::readCheckpointData");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::readCheckpointData" << endl;
      }

      MayDay::Error("FieldStepper<T>::readCheckpointData - checkpointing not supported for this module");
    }
#endif

    template <class T>
    int
    FieldStepper<T>::getNumberOfPlotVariables() const
    {
      CH_TIME("FieldStepper::getNumberOfPlotVariables");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::getNumberOfPlotVariables" << endl;
      }

      int ncomp = 0;

      ncomp += m_fieldSolver->getNumberOfPlotVariables();
      ncomp += m_sigma->getNumberOfPlotVariables();

      return ncomp;
    }

    template <class T>
    Vector<std::string>
    FieldStepper<T>::getPlotVariableNames() const
    {
      CH_TIME("FieldStepper::getPlotVariableNames");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::getPlotVariableNames" << endl;
      }

      Vector<std::string> plotVars;

      plotVars.append(m_fieldSolver->getPlotVariableNames());
      plotVars.append(m_sigma->getPlotVariableNames());

      return plotVars;
    }

    template <class T>
    void
    FieldStepper<T>::writePlotData(LevelData<EBCellFAB>& a_output,
                                   int&                  a_icomp,
                                   const std::string     a_outputRealm,
                                   const int             a_level) const
    {
      CH_TIME("FieldStepper::writePlotData");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::writePlotData" << endl;
      }

      CH_assert(a_level >= 0);
      CH_assert(a_level <= m_amr->getFinestLevel());

      // Write into plot data holder memory.
      m_fieldSolver->writePlotData(a_output, a_icomp, a_outputRealm, a_level);
      m_sigma->writePlotData(a_output, a_icomp, a_outputRealm, a_level);
    }

    template <class T>
    void
    FieldStepper<T>::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt)
    {
      CH_TIME("FieldStepper::synchronizeSolverTimes");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::synchronizeSolverTimes" << endl;
      }

      m_timeStep = a_step;
      m_time     = a_time;
      m_dt       = a_dt;

      m_fieldSolver->setTime(a_step, a_time, a_dt);
    }

    template <class T>
    void
    FieldStepper<T>::preRegrid(const int a_lbase, const int a_oldFinestLevel)
    {
      CH_TIME("FieldStepper::preRegrid");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::preRegrid" << endl;
      }

      m_fieldSolver->preRegrid(a_lbase, a_oldFinestLevel);
      m_sigma->preRegrid(a_lbase, a_oldFinestLevel);
    }

    template <class T>
    void
    FieldStepper<T>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
    {
      CH_TIME("FieldStepper::regrid");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::regrid" << endl;
      }

      // TLDR: The FieldSolver regrid methods just regrids the data -- it does not re-solve the Poisson equation. Here,
      //       that is done in postRegrid.

      m_fieldSolver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
      m_sigma->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
    }

    template <class T>
    void
    FieldStepper<T>::postRegrid()
    {
      CH_TIME("FieldStepper::postRegrid");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::postRegrid" << endl;
      }

      this->solvePoisson();
    }

    template <class T>
    bool
    FieldStepper<T>::loadBalanceThisRealm(const std::string a_realm) const
    {
      CH_TIME("FieldStepper::loadBalanceThisRealm");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::loadBalanceThisRealm" << endl;
      }

      return (a_realm == m_realm) && m_loadBalance;
    }

    template <class T>
    void
    FieldStepper<T>::loadBalanceBoxes(Vector<Vector<int>>&             a_procs,
                                      Vector<Vector<Box>>&             a_boxes,
                                      const std::string                a_realm,
                                      const Vector<DisjointBoxLayout>& a_grids,
                                      const int                        a_lmin,
                                      const int                        a_finestLevel)
    {
      CH_TIME("FieldStepper::loadBalanceBoxes");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::loadBalanceBoxes" << endl;
      }

      CH_assert(m_loadBalance && a_realm == m_realm);

      // TLDR: This code tries to compute a load for each grid patch by applying a relaxation operator to each box. This means that the load
      //       should be a decent estimate that takes into account boundary conditions, coarse-fine interface arithmetic, and enlargened stencils
      //       near the embedded boundary.

      a_procs.resize(1 + a_finestLevel);
      a_boxes.resize(1 + a_finestLevel);

      // We need to make AmrMesh restore some operators that we need in order to create a multigrid object. Fortunately, FieldSolver has routines
      // for doing that but it will not know if AmrMesh has updated it's operators or not. So, we need to regrid them.
      m_amr->regridOperators(a_lmin);

      // Field solver implementation gets the responsibility of computing loads on each level.
      for (int lvl = 0; lvl <= a_finestLevel; lvl++) {
        Vector<long long> loads = m_fieldSolver->computeLoads(a_grids[lvl], lvl);

        // Do the desired sorting and load balancing
        a_boxes[lvl] = a_grids[lvl].boxArray();

        LoadBalancing::sort(a_boxes[lvl], loads, m_boxSort);
        LoadBalancing::makeBalance(a_procs[lvl], loads, a_boxes[lvl]);
      }
    }

    template <class T>
    void
    FieldStepper<T>::setRho(const std::function<Real(const RealVect& a_pos)>& a_rho,
                            const phase::which_phase                          a_phase) noexcept
    {
      CH_TIME("FieldStepper::setRho");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::setRho" << endl;
      }

      if (a_phase == phase::gas) {
        m_rhoGas = a_rho;
      }
      else if (a_phase == phase::solid) {
        m_rhoDielectric = a_rho;
      }
    }

    template <class T>
    void
    FieldStepper<T>::setSigma(const std::function<Real(const RealVect& a_pos)>& a_sigma) noexcept
    {
      CH_TIME("FieldStepper::setSigma");
      if (m_verbosity > 5) {
        pout() << "FieldStepper<T>::setSigma" << endl;
      }

      m_surfaceChargeDensity = a_sigma;
    }
  } // namespace Electrostatics
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
