/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_StreamerInceptionStepperImplem.H
  @brief  Implementation of CD_StreamerInceptionStepper.H
  @author Robert Marskar
  @author Rasmus Hoholm
*/

#ifndef CD_StreamerInceptionStepperImplem_H
#define CD_StreamerInceptionStepperImplem_H

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_StreamerInceptionStepper.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::StreamerInception;

template <typename P, typename F>
StreamerInceptionStepper<P, F>::StreamerInceptionStepper()
{
  CH_TIME("StreamerInceptionStepper::StreamerInceptionStepper");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::StreamerInceptionStepper" << endl;
  }

  m_realm = Realm::Primal;
  m_phase = phase::gas;

  m_verbosity = -1;

  this->parseOptions();

  // Default stuff.
  m_alpha = [](const Real E) {
    return 1.0;
  };

  m_backgroundRate = [](const Real E) {
    return 1. / E;
  };

  m_fieldEmission = [](const Real E) {
    return 1. / E;
  };

  m_ionMobility = [](const Real E) {
    return 1.0;
  };
}

template <typename P, typename F>
StreamerInceptionStepper<P, F>::~StreamerInceptionStepper()
{
  CH_TIME("StreamerInceptionStepper::~StreamerInceptionStepper");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::~StreamerInceptionStepper" << endl;
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::setupSolvers()
{
  CH_TIME("StreamerInceptionStepper::setupSolvers");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setupSolvers" << endl;
  }

  // Define a voltage function to be used in the simulation.
  auto voltage = [](const Real a_time) -> Real {
    return 1.0;
  };

  // Instantiate the field solver.
  m_fieldSolver = RefCountedPtr<FieldSolver>(new F());
  m_fieldSolver->setVerbosity(m_verbosity);
  m_fieldSolver->parseOptions();
  m_fieldSolver->setAmr(m_amr);
  m_fieldSolver->setComputationalGeometry(m_computationalGeometry);
  m_fieldSolver->setVoltage(voltage);
  m_fieldSolver->setRealm(m_realm);
  m_fieldSolver->setTime(0, 0.0, 0.0);

  // Instantiate the tracer particle solver.
  m_particleSolver = RefCountedPtr<TracerParticleSolver<P>>(new TracerParticleSolver<P>());
  m_particleSolver->parseOptions();
  m_particleSolver->setAmr(m_amr);
  m_particleSolver->setComputationalGeometry(m_computationalGeometry);
  m_particleSolver->setRealm(m_realm);
  m_particleSolver->setPhase(m_phase);
  m_particleSolver->setTime(0, 0.0, 0.0);
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::allocate()
{
  CH_TIME("StreamerInceptionStepper::allocate");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::allocate" << endl;
  }

  m_fieldSolver->allocateInternals();
  m_particleSolver->allocate();

  m_amr->allocate(m_electricField, m_realm, m_phase, SpaceDim);
  m_amr->allocate(m_inceptionIntegral, m_realm, m_phase, m_voltageSweeps.size());
  m_amr->allocate(m_backgroundIonization, m_realm, m_phase, m_voltageSweeps.size());
  m_amr->allocate(m_emissionRatesPlus, m_realm, m_phase, m_voltageSweeps.size());
  m_amr->allocate(m_emissionRatesMinu, m_realm, m_phase, m_voltageSweeps.size());
  m_amr->allocate(m_inceptionVoltage, m_realm, m_phase, 1);
  m_amr->allocate(m_driftTimesPlus, m_realm, m_phase, m_voltageSweeps.size());
  m_amr->allocate(m_driftTimesMinu, m_realm, m_phase, m_voltageSweeps.size());

  DataOps::setValue(m_inceptionIntegral, 0.0);
  DataOps::setValue(m_backgroundIonization, 0.0);
  DataOps::setValue(m_emissionRatesPlus, 0.0);
  DataOps::setValue(m_emissionRatesMinu, 0.0);
  DataOps::setValue(m_inceptionVoltage, 0.0);
  DataOps::setValue(m_driftTimesPlus, 0.0);
  DataOps::setValue(m_driftTimesMinu, 0.0);
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::initialData()
{
  CH_TIME("StreamerInceptionStepper::initialData");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::initialData" << endl;
  }

  // Set rho = sigma = 0 for field solver and then solve the Poisson equation.
  m_fieldSolver->setRho(0.0);
  m_fieldSolver->setSigma(0.0);

  const bool converged = m_fieldSolver->solve(m_fieldSolver->getPotential(),
                                              m_fieldSolver->getRho(),
                                              m_fieldSolver->getSigma());

  if (!converged) {
    MayDay::Error("StreamerInceptionStepper::initialData -- could not solve Poisson equation. Aborting application.");
  }

  // Get the gas-side electric field and give it to the particle solvers.
  DataOps::setValue(m_electricField, 0.0);
  DataOps::copy(m_electricField, m_amr->alias(phase::gas, m_fieldSolver->getElectricField()));
  m_amr->averageDown(m_electricField, m_realm, m_phase);
  m_amr->interpGhost(m_electricField, m_realm, m_phase);
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::registerRealms()
{
  CH_TIME("StreamerInceptionStepper::registerRealms");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::registerRealms" << endl;
  }

  m_amr->registerRealm(m_realm);
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::registerOperators()
{
  CH_TIME("StreamerInceptionStepper::registerOperators");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::registerOperators" << endl;
  }

  m_fieldSolver->registerOperators();
  m_particleSolver->registerOperators();
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::parseOptions()
{
  CH_TIME("StreamerInceptionStepper::parseOptions");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::parseOptions" << endl;
  }

  ParmParse pp("StreamerInceptionStepper");

  Real voltageLo    = 0.0;
  Real voltageHi    = 0.0;
  Real stepSize     = 0.0;
  int  voltageSteps = 0;

  std::string str;

  m_plotInceptionIntegral    = false;
  m_plotInceptionVoltage     = false;
  m_plotBackgroundIonization = false;
  m_plotFieldEmission        = false;
  m_plotDriftTime            = false;

  // Read in parameters
  pp.get("verbosity", m_verbosity);
  pp.get("voltage_lo", voltageLo);
  pp.get("voltage_hi", voltageHi);
  pp.get("voltage_steps", voltageSteps);
  pp.get("K_inception", m_inceptionK);
  pp.get("integration", str);
  pp.get("print_report", m_printReport);

  if (str == "euler") {
    m_integrationAlgorithm = IntegrationAlgorithm::Euler;
  }
  else if (str == "trapz") {
    m_integrationAlgorithm = IntegrationAlgorithm::Trapezoidal;
  }
  else if (str == "simp") {
    m_integrationAlgorithm = IntegrationAlgorithm::Simpson;
  }
  else {
    MayDay::Error("StreamerInceptionStepper::StreamerInceptionStepper -- logic bust");
  }

  // Get plot variables.
  const int num = pp.countval("plt_vars");
  if (num > 0) {
    Vector<std::string> plotVars(num);
    pp.getarr("plt_vars", plotVars, 0, num);

    for (int i = 0; i < num; i++) {
      if (plotVars[i] == "K") {
        m_plotInceptionIntegral = true;
      }
      else if (plotVars[i] == "Uinc") {
        m_plotInceptionVoltage = true;
      }
      else if (plotVars[i] == "bg_rate") {
        m_plotBackgroundIonization = true;
      }
      else if (plotVars[i] == "emission") {
        m_plotFieldEmission = true;
      }
      else if (plotVars[i] == "drift") {
        m_plotDriftTime = true;
      }
    }
  }

  // Get step size
  pp.get("step_size", str, 0);
  pp.get("step_size", stepSize, 1);
  if (str == "dx") {
    this->setRelativeStep(stepSize);
  }
  else if (str == "fixed") {
    this->setFixedStep(stepSize);
  }

  // Define voltages
  const Real dV = (voltageHi - voltageLo) / voltageSteps;

  m_voltageSweeps.push_back(voltageLo);

  for (size_t i = 0; i < voltageSteps; i++) {
    m_voltageSweeps.push_back(m_voltageSweeps[i] + dV);
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::parseRuntimeOptions()
{
  CH_TIME("StreamerInceptionStepper::parseRuntimeOptions");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::parseRuntimeOptions" << endl;
  }
}

#ifdef CH_USE_HDF5
template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const
{
  CH_TIME("StreamerInceptionStepper::writeCheckpointData");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::writeCheckpointData" << endl;
  }
}
#endif

#ifdef CH_USE_HDF5
template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::readCheckpointData(HDF5Handle& a_handle, const int a_lvl)
{
  CH_TIME("StreamerInceptionStepper::readCheckpointData");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::readCheckpointData" << endl;
  }

  MayDay::Error("StreamerInceptionStepper::readCheckpointData -- restart not supported. Use Driver.restart=0");
}
#endif

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::writePlotData(EBAMRCellData&       a_output,
                                              Vector<std::string>& a_plotVariableNames,
                                              int&                 a_icomp) const
{
  CH_TIME("StreamerInceptionStepper::writePlotData");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::writePlotData" << endl;
  }

  // Append plot variable names.
  a_plotVariableNames.append(m_fieldSolver->getPlotVariableNames());
  a_plotVariableNames.append(m_particleSolver->getPlotVariableNames());

  // Write into plot data holder memory.
  m_fieldSolver->writePlotData(a_output, a_icomp);
  m_particleSolver->writePlotData(a_output, a_icomp);

  // Write our K values.
  if (m_plotInceptionIntegral) {
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const Interval srcInterval(0, m_voltageSweeps.size() - 1);
      const Interval dstInterval(a_icomp, a_icomp + m_voltageSweeps.size() - 1);
      m_inceptionIntegral[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
    }
    a_icomp += m_voltageSweeps.size();

    std::string varName;
    for (const Real& V : m_voltageSweeps) {
      varName = "K-value, V = ";
      varName += std::to_string(V);
      a_plotVariableNames.push_back(varName);
    }
  }

  // Write our drift time values.
  if (m_plotDriftTime) {

    // Positive polarity
    {
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const Interval srcInterval(0, m_voltageSweeps.size() - 1);
        const Interval dstInterval(a_icomp, a_icomp + m_voltageSweeps.size() - 1);
        m_driftTimesPlus[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
      }
      a_icomp += m_voltageSweeps.size();

      std::string varName;
      for (const Real& V : m_voltageSweeps) {
        varName = "Drift time for V = +";
        varName += std::to_string(V);
        a_plotVariableNames.push_back(varName);
      }
    }

    // Negative polarity
    {
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const Interval srcInterval(0, m_voltageSweeps.size() - 1);
        const Interval dstInterval(a_icomp, a_icomp + m_voltageSweeps.size() - 1);
        m_driftTimesMinu[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
      }
      a_icomp += m_voltageSweeps.size();

      std::string varName;
      for (const Real& V : m_voltageSweeps) {
        varName = "Drift time for V = -";
        varName += std::to_string(V);
        a_plotVariableNames.push_back(varName);
      }
    }
  }

  // Write the inception voltage.
  if (m_plotInceptionVoltage) {
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const Interval srcInterval(0, 0);
      const Interval dstInterval(a_icomp, a_icomp);
      m_inceptionVoltage[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
    }
    a_icomp += 1;

    a_plotVariableNames.push_back("Inception voltage");
  }

  // Write the background ionization rates
  if (m_plotBackgroundIonization) {
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const Interval srcInterval(0, m_voltageSweeps.size() - 1);
      const Interval dstInterval(a_icomp, a_icomp + m_voltageSweeps.size() - 1);
      m_backgroundIonization[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
    }
    a_icomp += m_voltageSweeps.size();

    std::string varName;
    for (const Real& V : m_voltageSweeps) {
      varName = "Background ionization rate, V =  ";
      varName += std::to_string(V);
      a_plotVariableNames.push_back(varName);
    }
  }

  // Write field emission rates
  if (m_plotFieldEmission) {
    // Positive polarity
    {
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const Interval srcInterval(0, m_voltageSweeps.size() - 1);
        const Interval dstInterval(a_icomp, a_icomp + m_voltageSweeps.size() - 1);
        m_emissionRatesPlus[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
      }

      std::string varName;
      for (const Real& V : m_voltageSweeps) {
        varName = "Emission rate, V = +";
        varName += std::to_string(V);
        a_plotVariableNames.push_back(varName);
      }

      a_icomp += m_voltageSweeps.size();
    }

    // Negative polarity
    {
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const Interval srcInterval(0, m_voltageSweeps.size() - 1);
        const Interval dstInterval(a_icomp, a_icomp + m_voltageSweeps.size() - 1);
        m_emissionRatesMinu[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
      }

      std::string varName;
      for (const Real& V : m_voltageSweeps) {
        varName = "Emission rate, V = -";
        varName += std::to_string(V);
        a_plotVariableNames.push_back(varName);
      }

      a_icomp += m_voltageSweeps.size();
    }
  }
}

template <typename P, typename F>
int
StreamerInceptionStepper<P, F>::getNumberOfPlotVariables() const
{
  CH_TIME("StreamerInceptionStepper::getNumberOfPlotVariables");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::getNumberOfPlotVariables" << endl;
  }

  int ncomp = 0;

  ncomp += m_fieldSolver->getNumberOfPlotVariables();
  ncomp += m_particleSolver->getNumberOfPlotVariables();

  // K-values
  if (m_plotInceptionIntegral) {
    ncomp += m_voltageSweeps.size();
  }

  // Drift times (+/-) direction
  if (m_plotDriftTime) {
    ncomp += 2 * m_voltageSweeps.size();
  }

  // Inception voltage.
  if (m_plotInceptionVoltage) {
    ncomp += 1;
  }

  // Background ionization rates
  if (m_plotBackgroundIonization) {
    ncomp += m_voltageSweeps.size();
  }

  // Field emission rates
  if (m_plotFieldEmission) {
    ncomp += 2 * m_voltageSweeps.size();
  }

  return ncomp;
}

template <typename P, typename F>
Real
StreamerInceptionStepper<P, F>::computeDt()
{
  CH_TIME("StreamerInceptionStepper::computeDt");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeDt" << endl;
  }

  return std::numeric_limits<Real>::max();
}

template <typename P, typename F>
Real
StreamerInceptionStepper<P, F>::advance(const Real a_dt)
{
  CH_TIME("StreamerInceptionStepper::advance");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::advance" << endl;
  }

  MayDay::Error("StreamerInceptionStepper::advance -- logic bust (Use Driver.max_steps = 0)");

  return a_dt;
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt)
{
  CH_TIME("StreamerInceptionStepper::synchronizeSolverTimes");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::synchronizeSolverTimes" << endl;
  }

  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;

  m_fieldSolver->setTime(a_step, a_time, a_dt);
  m_particleSolver->setTime(a_step, a_time, a_dt);
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::preRegrid(const int a_lmin, const int a_oldFinestLevel)
{
  CH_TIME("StreamerInceptionStepper::preRegrid");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::preRegrid" << endl;
  }

  m_fieldSolver->preRegrid(a_lmin, a_oldFinestLevel);
  m_particleSolver->preRegrid(a_lmin, a_oldFinestLevel);
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("StreamerInceptionStepper::regrid");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::regrid" << endl;
  }

  // Regrid tracer particles and field
  m_fieldSolver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_particleSolver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);

  // Regrid velocity field
  m_amr->reallocate(m_electricField, m_phase, a_lmin);
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::postRegrid()
{
  CH_TIME("StreamerInceptionStepper::postRegrid");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::postRegrid" << endl;
  }

  // This can be left empty because initialData is called right after the regrid method and calling
  // the advance method is an error.
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::setFixedStep(const Real a_step) noexcept
{
  CH_TIME("StreamerInceptionStepper::setFixedStep");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setFixedStep" << endl;
  }

  m_useFixedStep = true;
  m_fixedStep    = a_step;
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::setRelativeStep(const Real a_step) noexcept
{
  CH_TIME("StreamerInceptionStepper::setRelativeStep");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setRelativeStep" << endl;
  }

  m_useFixedStep = false;
  m_relativeStep = a_step;
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::setAlpha(const std::function<Real(const Real& E)>& a_alpha) noexcept
{
  CH_TIME("StreamerInceptionStepper::setAlpha");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setAlpha" << endl;
  }

  m_alpha = a_alpha;
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::setIonMobility(const std::function<Real(const Real& E)>& a_mobility) noexcept
{
  CH_TIME("StreamerInceptionStepper::setIonMobility");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setIonMobility" << endl;
  }

  m_ionMobility = a_mobility;
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::setBackgroundIonization(
  const std::function<Real(const Real& E)>& a_backgroundIonization) noexcept
{
  CH_TIME("StreamerInceptionStepper::setBackgroundIonization");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setBackgroundIonization" << endl;
  }

  m_backgroundRate = a_backgroundIonization;
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::setFieldEmission(const std::function<Real(const Real& E)>& a_currentDensity) noexcept
{
  CH_TIME("StreamerInceptionStepper::setFieldEmission");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setFieldEmission" << endl;
  }

  m_fieldEmission = a_currentDensity;
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::initialParticles() noexcept
{
  CH_TIME("StreamerInceptionStepper::initialParticles");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::initialParticles" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();
  amrParticles.clearParticles();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

    const LevelData<BaseFab<bool>>& validCellsLD = *m_amr->getValidCells(m_realm)[lvl];

    const Real dx = m_amr->getDx()[lvl];

    ParticleData<P>& levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const EBISBox&       ebisbox    = ebisl[dit()];
      const BaseFab<bool>& validCells = validCellsLD[dit()];

      List<P>& particles = levelParticles[dit()].listItems();

      auto regularKernel = [&](const IntVect& iv) -> void {
        if (validCells(iv, 0) && ebisbox.isRegular(iv)) {
          const RealVect pos = m_amr->getProbLo() + (0.5 * RealVect::Unit + RealVect(iv)) * dx;
          particles.add(P(pos));
        }
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        if (validCells(vof.gridIndex())) {
          const RealVect pos = m_amr->getProbLo() + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);
          particles.add(P(pos));
        }
      };

      // Execute kernels over appropriate regions.
      const Box    cellBox = dbl[dit()];
      VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);

      // Define particles' masses (i.e. integration & drift time tracker) and integration flag
      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        // Set mass to zero and do a backup of the initial position (we need to rewind particles
        // back later on). The scalar<0> flag is the flag we use for start-stop criteria for the
        // particle integration.
        p.mass()               = 0.0;
        p.template vect<0>()   = p.position();
        p.template scalar<0>() = 1.0;
      }
    }
  }

  // Remove particles inside the EB.
  m_amr->removeCoveredParticlesIF(amrParticles, m_phase);

  amrParticles.remap();
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::postInitialize()
{
  CH_TIME("StreamerInceptionStepper::postInitialize");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::postInitialize" << endl;
  }

  // Initialize the particle solver.
  m_particleSolver->setVelocity(m_electricField);

  this->initialParticles();

  // Solve K integral for each particle in each voltage and store K values in m_inceptionIntegral
  this->computeInceptionIntegral();

  // Solve drift time integral in both directions and store values in m_driftTimes
  this->computeDriftTimes();

  // Compute detachment rates
  this->computeBackgroundIonization();

  // Compute field emission rates
  this->computeFieldEmission();

  // Interpolate between K values to find voltage giving K_inception and store values in m_inceptionVoltage
  this->computeInceptionVoltageVolume();

  // Compute critical volume of K values for each voltage
  this->computeCriticalVolume();

  // Compute time to first critical electron
  this->computeFirstElectron();

  // Print voltages, max K values, and critical volumes to terminal
  if (m_printReport) {
    this->printReport();
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::computeInceptionIntegral() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeInceptionIntegral");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeInceptionIntegral" << endl;
  }

  // TLDR: For this integration we compute the value of the inception integral
  //
  //          K = int alpha(E)dl
  //
  //       The integration runs from the start position of the particle and along a field line until either the
  //       particle leaves the domain or alpha < 0.0.
  //
  //       For the Euler rule we get
  //
  //          T += alpha(E(x)) * dx
  //
  //       For the trapezoidal rule we get
  //
  //          T += 0.5 * dx * [alpha(E(x)) + alpha(E(x+dx))]
  //
  //       For Simpsons rule we get
  //
  //          T += dx / 6 * [alpha(E(x)) + 4 * alpha(E(x+dx/2)) + alpha(E(x+dx))]

  // Transient storage we can deposit particles onto.
  EBAMRCellData K;
  EBAMRCellData Kplus;
  EBAMRCellData Kminu;

  m_amr->allocate(K, m_realm, m_phase, 1);
  m_amr->allocate(Kplus, m_realm, m_phase, 1);
  m_amr->allocate(Kminu, m_realm, m_phase, 1);

  DataOps::setValue(K, 0.0);
  DataOps::setValue(Kplus, 0.0);
  DataOps::setValue(Kminu, 0.0);

  // Integrate for each voltage
  for (int i = 0; i < m_voltageSweeps.size(); ++i) {

    // Polarities. Note that we are dealing with electrons ions so
    // for positive polarity the particles move opposite to the field.
    std::vector<Real> polarities{-1.0, 1.0};

    for (const auto& p : polarities) {

      // Draw a new set of particles and update the velocity.
      this->initialParticles();

      // Scale E-field by +/-
      DataOps::scale(m_electricField, p);

      // Switch between various integration algorithms. At the end of this
      // the K-value time should be stored on the particle mass.
      switch (m_integrationAlgorithm) {
      case IntegrationAlgorithm::Euler: {
	this->inceptionIntegrateEuler(m_voltageSweeps[i]);

        break;
      }
      case IntegrationAlgorithm::Trapezoidal: {
        this->inceptionIntegrateTrapezoidal(m_voltageSweeps[i]);

        break;
      }
      case IntegrationAlgorithm::Simpson: {
        this->inceptionIntegrateSimpson(m_voltageSweeps[i]);

        break;
      }
      default: {
        MayDay::Error("StreamerInceptionStepper::computeInceptionIntegral -- logic bust");

        break;
      }
      }

      // Move particles back to initial position and deposit them on the mesh.
      this->rewindParticles();

      if (p > 0.0) {
        m_particleSolver->deposit(Kplus);
      }
      else {
        m_particleSolver->deposit(Kminu);
      }

      // Scale field back again.
      DataOps::scale(m_electricField, p);
    }

    // Make K = max(Kplus, Kminu)
    DataOps::max(K, Kplus, Kminu);

    m_amr->averageDown(K, m_realm, m_phase);
    m_amr->interpGhost(K, m_realm, m_phase);

    // Get max K value.
    Real maxK = -std::numeric_limits<Real>::max();
    Real minK = +std::numeric_limits<Real>::max();

    DataOps::getMaxMin(maxK, minK, K, 0);

    m_maxKValues.push_back(maxK);

    // Copy K values to m_inceptionIntegral so we can plot them later.
    DataOps::copy(m_inceptionIntegral, K, Interval(i, i), Interval(0, 0));
  }

  this->initialParticles();
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::inceptionIntegrateEuler(const Real& a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::inceptionIntegrateEuler");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::inceptionIntegrateEuler" << endl;
  }

  const RealVect probLo = m_amr->getProbLo();
  const RealVect probHi = m_amr->getProbHi();

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  while (this->keepIntegrating()) {

    m_particleSolver->interpolateVelocities();

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const Real step = m_useFixedStep ? m_fixedStep : m_relativeStep * m_amr->getDx()[lvl];

      const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      
      for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
	const EBISBox& ebisbox = ebisl[dit()];
	
        for (ListIterator<P> lit(amrParticles[lvl][dit()].listItems()); lit.ok(); ++lit) {

          P& p = lit();

	  if(p.template scalar<0>() > 0.0) {
	    const RealVect v      = p.velocity();
	    const RealVect E      = a_voltage * p.velocity();
	    const RealVect unitV  = v / v.vectorLength();
	    const RealVect newPos = p.position() + step * unitV;
	    const Real     alpha  = m_alpha(E.vectorLength());

	    // Stop integration if particle ends up inside the EB or outside the domain, or
	    // alpha < 0.0.
	    const bool outsideGrid = this->particleOutsideGrid(newPos, probLo, probHi);
	    const bool insideEB    = this->particleInsideEB(newPos);

	    if (outsideGrid || insideEB || alpha < 0.0) {
	      p.template scalar<0>() = -1.0;
	      
	      p.position() = p.template vect<0>();
	    }
	    else {
	      p.position() = newPos;

	      p.mass() += step * alpha;
	    }
	  }
	}
      }

    }

    amrParticles.remap();
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::inceptionIntegrateTrapezoidal(const Real& a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::inceptionIntegrateTrapezoidal");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::inceptionIntegrateTrapezoidal" << endl;
  }

  MayDay::Abort("CD_StreamerInceptionStepperImplem.H - inceptionIntegrateTrapezoidal -- not implemented");
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::inceptionIntegrateSimpson(const Real& a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::inceptionIntegrateSimpson");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::inceptionIntegrateSimpson" << endl;
  }

  MayDay::Abort("CD_StreamerInceptionStepperImplem.H - inceptionIntegrateSimpson -- not implemented");
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::computeDriftTimes() noexcept
{

  CH_TIME("StreamerInceptionStepper::computeDriftTimes");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeDriftTimes" << endl;
  }

  // TLDR: For this integration we compute the time it takes for ions to drift into the critical volume. For each ion
  //       we compute
  //
  //          T = int dl/|v| where |v| = mu * E.
  //
  //       The integration runs from the start position of the particle and until it arrives inside the critical
  //       volume. To figure out the critical volume, we interpolate it onto the particle's "mass". If the particle
  //       starts out inside the critical volume then T = 0.0, and if it never reaches the critical volume then we
  //       set T = inf.
  //
  //       For the Euler rule we get
  //
  //          T += dx / |v(x)|
  //
  //       For the trapezoidal rule we get
  //
  //          T += 0.5 * dx /(|v(x)| + |v(x+dx)|)
  //
  //       For Simpsons rule we get
  //
  //          T += dx /  6 * [1/|v(x)| + 4 / v(x + dx/2) + 1 / v(x+dx)]

  // Transient storage we can deposit into.
  EBAMRCellData K;
  EBAMRCellData T;

  m_amr->allocate(T, m_realm, m_phase, 1);
  m_amr->allocate(K, m_realm, m_phase, 1);

  DataOps::setValue(T, 0.0);
  DataOps::setValue(K, 0.0);

  // Integrate for each voltage
  for (int i = 0; i < m_voltageSweeps.size(); ++i) {

    // Copy current K-value to the above data holder.
    DataOps::copy(K, m_inceptionIntegral, Interval(0, 0), Interval(i, i));
    m_amr->averageDown(K, m_realm, m_phase);
    m_amr->interpGhost(K, m_realm, m_phase);

    // Polarities. Note that we are dealing with negative ions so
    // for positive polarity the particles move opposite to the field.
    std::vector<Real> polarities{-1.0, 1.0};

    for (const auto& p : polarities) {

      // Draw a new set of particles (negative ions) and update the velocity.
      this->initialParticles();

      // Scale E-field by +/-
      DataOps::scale(m_electricField, p);

      // Switch between various integration algorithms. At the end of this
      // the drift time should be stored on the particle mass.
      switch (m_integrationAlgorithm) {
      case IntegrationAlgorithm::Euler: {

        this->driftIntegrateEuler(K, m_voltageSweeps[i]);

        break;
      }
      case IntegrationAlgorithm::Trapezoidal: {

        MayDay::Abort("CD_StreamerInceptionStepper::computeDriftTimes -- Trapezoidal not supported");

        break;
      }
      case IntegrationAlgorithm::Simpson: {

        MayDay::Abort("CD_StreamerInceptionStepper::computeDriftTimes -- Simpsons not supported");

        break;
      }
      default: {
        MayDay::Error("StreamerInceptionStepper::computeDriftTimes -- logic bust");

        break;
      }
      }

      // Move particles back to initial position and deposit them on the mesh.
      this->rewindParticles();
      m_particleSolver->deposit(T);

      // Copy drift time to solver data holder
      if (p < 0.0) {
        DataOps::copy(m_driftTimesPlus, T, Interval(i, i), Interval(0, 0));
      }
      else {
        DataOps::copy(m_driftTimesMinu, T, Interval(i, i), Interval(0, 0));
      }

      // Scale field back.
      DataOps::scale(m_electricField, p);
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::driftIntegrateEuler(const EBAMRCellData& a_K, const Real& a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::driftIntegrateEuler");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::driftIntegrateEuler" << endl;
  }

  const RealVect probLo = m_amr->getProbLo();
  const RealVect probHi = m_amr->getProbHi();

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  while (this->keepIntegrating()) {

    // Get K-value and velocity at particle position.
    m_particleSolver->interpolateMass(a_K);
    m_particleSolver->interpolateVelocities();

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const Real step = m_useFixedStep ? m_fixedStep : m_relativeStep * m_amr->getDx()[lvl];

      for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
        for (ListIterator<P> lit(amrParticles[lvl][dit()].listItems()); lit.ok(); ++lit) {

          P& p = lit();

	  if(p.template scalar<0>() > 0.0) {
	    const RealVect E    = a_voltage * p.velocity();
	    const Real     absE = E.vectorLength();
	    const Real     mu   = m_ionMobility(E.vectorLength());

	    const RealVect vel    = E;
	    const RealVect unitV  = vel/vel.vectorLength();
	    const RealVect newPos = p.position() + unitV * step;

	    // Stop integration if particle ends up inside the EB, outside the domain, or if
	    // the K-value is less than the inception criterion.
	    const bool outsideGrid = this->particleOutsideGrid(newPos, probLo, probHi);
	    const bool insideEB    = this->particleInsideEB(newPos);
	    const bool insideCrit  = p.mass() >= m_inceptionK;

	    if (insideCrit) {
	      p.template scalar<0>() = -1.0;
	    }
	    else if (outsideGrid || insideEB) {
	      // Stop integration and discard result because the ion
	      // did not enter the critical volume.
	      p.template scalar<0>() = -1.0;
	      p.template scalar<1>() = std::numeric_limits<Real>::quiet_NaN();
	    }
	    else {
	      p.position() = newPos;

	      p.template scalar<1>() += step / (mu * absE);
	    }
	  }
	}
      }
    }

    m_particleSolver->remap();
  }

  // Put the integration variable onto the particle "mass"
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
      for (ListIterator<P> lit(amrParticles[lvl][dit()].listItems()); lit.ok(); ++lit) {

        P& p = lit();

        p.mass() = p.template scalar<1>();
      }
    }
  }
}

template <typename P, typename F>
bool
StreamerInceptionStepper<P, F>::keepIntegrating() const noexcept
{
  CH_TIME("StreamerInceptionStepper::keepIntegrating");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::keepIntegrating" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  // If keepGoing > 0 then we keep integrating. Otherwise we stop.
  int keepGoing = -1;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
      for (ListIterator<P> lit(amrParticles[lvl][dit()].listItems()); lit.ok(); ++lit) {
        if (lit().template scalar<0>() > 0.0) {
          keepGoing = 1;
        }
      }
    }
  }

  return (ParallelOps::max(keepGoing) > 0) ? true : false;
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::rewindParticles() noexcept
{
  CH_TIME("StreamerInceptionStepper::rewindParticles");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::rewindParticles" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
      for (ListIterator<P> lit(amrParticles[lvl][dit()].listItems()); lit.ok(); ++lit) {
        P& p = lit();

	p.position() = p.template vect<0>();
      }
    }
  }

  amrParticles.remap();
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::computeBackgroundIonization() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeBackgroundIonization");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeBackgroundIonization" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {

      EBCellFAB& bgIonization    = (*m_backgroundIonization[lvl])[dit()];
      FArrayBox& bgIonizationReg = bgIonization.getFArrayBox();

      const EBCellFAB& electricField    = (*m_electricField[lvl])[dit()];
      const FArrayBox& electricFieldReg = electricField.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        const Real E = RealVect(D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)))
                         .vectorLength();

        for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
          bgIonizationReg(iv, i) = m_backgroundRate(m_voltageSweeps[i] * E);
        }
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const Real E = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)))
                         .vectorLength();

        for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
          bgIonization(vof, i) = m_backgroundRate(m_voltageSweeps[i] * E);
        }
      };

      // Execute kernels
      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      BoxLoops::loop(dbl[dit()], regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::computeFieldEmission() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeFieldEmission");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeFieldEmission" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const EBISBox& ebisbox = ebisl[dit()];

      // Here, emissionPlus is when we have the standard voltage (V = 1) and
      // emissionMinu is when the voltage is reverted.
      EBCellFAB& emissionPlus = (*m_emissionRatesPlus[lvl])[dit()];
      EBCellFAB& emissionMinu = (*m_emissionRatesMinu[lvl])[dit()];

      emissionPlus.setVal(0.0);
      emissionMinu.setVal(0.0);

      const EBCellFAB& electricField    = (*m_electricField[lvl])[dit()];
      const FArrayBox& electricFieldReg = electricField.getFArrayBox();

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const RealVect E = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));

        const Real normalE    = E.dotProduct(ebisbox.normal(vof));
        const Real magnitudeE = E.vectorLength();

        if (normalE > 0.0) {
          // Field points out of the electrode but for the opposite polarity it will
          // point into the electrode and the electrons can then escape.
          for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
            emissionPlus(vof, i) = 0.0;
            emissionMinu(vof, i) = m_fieldEmission(magnitudeE * m_voltageSweeps[i]);
          }
        }
        else if (normalE < 0.0) {
          // Opposite case of the above.
          for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
            emissionPlus(vof, i) = m_fieldEmission(magnitudeE * m_voltageSweeps[i]);
            emissionMinu(vof, i) = 0.0;
          }
        }
      };

      // Execute kernels
      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
      BoxLoops::loop(vofit, irregularKernel);
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::computeInceptionVoltageVolume() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeInceptionVoltageVolume");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeInceptionVoltageVolume" << endl;
  }

  // TLDR: This routine runs through all the K-values in each cell and estimates the
  //       inception voltage using linear interpolation.

  if (m_voltageSweeps.size() < 2) {
    DataOps::setValue(m_inceptionVoltage, std::numeric_limits<Real>::quiet_NaN());

    MayDay::Warning("StreamerInceptionStepper::computeInceptionVoltageVolume -- not enough voltages for estimating "
                    "inception voltage");
  }
  else {

    constexpr int comp = 0;

    // Function which interpolates the inception voltage if possible. Used in the kernels.
    auto interpVinc = [Kinc = this->m_inceptionK, &V = this->m_voltageSweeps](const std::vector<Real>& K) -> Real {
      Real Vinc = std::numeric_limits<Real>::quiet_NaN();
      for (size_t i = 0; i < K.size() - 1; i++) {

        if (K[i] <= Kinc && K[i + 1] > Kinc) {
          Vinc = V[i] + (Kinc - K[i]) * (V[i + 1] - V[i]) / (K[i + 1] - K[i]);

          break;
        }
      }

      return Vinc;
    };

    // Iterate through m_inceptionIntegral data and calculate inception voltage.
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); ++lvl) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        EBCellFAB& inceptionVoltage    = (*m_inceptionVoltage[lvl])[dit()];
        FArrayBox& inceptionVoltageReg = inceptionVoltage.getFArrayBox();

        const EBCellFAB& inceptionIntegral    = (*m_inceptionIntegral[lvl])[dit()];
        const FArrayBox& inceptionIntegralReg = inceptionIntegral.getFArrayBox();

        // Regular kernel.
        auto regularKernel = [&](const IntVect& iv) -> void {
          std::vector<Real> Kvalues;
          for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
            Kvalues.emplace_back(inceptionIntegralReg(iv, i));
          }

          inceptionVoltageReg(iv, comp) = interpVinc(Kvalues);
        };

        // Irregular kernel.
        auto irregularKernel = [&](const VolIndex& vof) -> void {
          std::vector<Real> Kvalues;
          for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
            Kvalues.emplace_back(inceptionIntegral(vof, i));
          }

          inceptionVoltage(vof, comp) = interpVinc(Kvalues);
        };

        // Kernel regions.
        const Box&   cellBox = dbl[dit()];
        VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

        BoxLoops::loop(cellBox, regularKernel);
        BoxLoops::loop(vofit, irregularKernel);
      }
    }

    // Coarsen data.
    m_amr->averageDown(m_inceptionVoltage, m_realm, m_phase);
    m_amr->interpGhost(m_inceptionVoltage, m_realm, m_phase);
  }
}

template <typename P, typename F>
Real
StreamerInceptionStepper<P, F>::computeInceptionVoltage() const noexcept
{
  CH_TIME("StreamerInceptionStepper::computeInceptionVoltage");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeInceptionVoltage" << endl;
  }

  Real Vinc = std::numeric_limits<Real>::quiet_NaN();

  if (m_voltageSweeps.size() < 2) {
    MayDay::Warning("StreamerInceptionStepper::computeInceptionVoltageVolume -- not enough voltages for estimating "
                    "inception voltage");
  }
  else {

    constexpr int comp = 0;

    auto interpVinc =
      [&Vinc, &K = this->m_maxKValues, &Kinc = this->m_inceptionK, &V = this->m_voltageSweeps]() -> void {
      for (size_t i = 0; i < K.size() - 1; i++) {

        if (K[i] <= Kinc && K[i + 1] > Kinc) {
          Vinc = V[i] + (Kinc - K[i]) * (V[i + 1] - V[i]) / (K[i + 1] - K[i]);

          break;
        }
      }
    };

    interpVinc();
  }

  return Vinc;
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::computeCriticalVolume() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeCriticalVolume");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeCriticalVolume" << endl;
  }

  // Solve critical volume of K values for each voltage
  for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
    Real criticalVolume = 0.0;

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); ++lvl) {
      const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

      const LevelData<BaseFab<bool>>& validCellsLD = *m_amr->getValidCells(m_realm)[lvl];

      const Real dx = m_amr->getDx()[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
        const EBISBox&       ebisbox    = ebisl[dit()];
        const BaseFab<bool>& validCells = validCellsLD[dit()];

        EBCellFAB&     cellData = (*m_inceptionIntegral[lvl])[dit()];
        BaseFab<Real>& regData  = cellData.getFArrayBox();

        auto regularKernel = [&](const IntVect& iv) -> void {
          if ((regData(iv, i) >= m_inceptionK) && validCells(iv, 0) && ebisbox.isRegular(iv)) {
            criticalVolume += std::pow(dx, SpaceDim);
          }
        };

        auto irregularKernel = [&](const VolIndex& vof) -> void {
          if ((cellData(vof, i) >= m_inceptionK) && validCells(vof.gridIndex())) {
            criticalVolume += ebisbox.volFrac(vof) * std::pow(dx, SpaceDim);
          }
        };

        // Kernel regions.
        const Box    cellBox = dbl[dit()];
        VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

        BoxLoops::loop(cellBox, regularKernel);
        BoxLoops::loop(vofit, irregularKernel);
      }
    }

    // Add critical volume to m_criticalVolumes
    m_criticalVolumes.push_back(ParallelOps::sum(criticalVolume));
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::computeFirstElectron() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeFirstElectron");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeFirstElectron" << endl;
  }

  // TLDR: This routine computes the average time to the appereance of an electron
  //       anywhere in the critical volume. We do this by assuming that the various
  //       processes involved are Poisson processes P(r1), P(r2), ..., where r1
  //       is the rate for the first Poisson process. The total rate is then given
  //       by R = r1 + r2 + ... and the average time is 1/R.

  // Total rate.
  std::vector<Real> Rplus(m_voltageSweeps.size(), std::numeric_limits<Real>::min());
  std::vector<Real> Rminu(m_voltageSweeps.size(), std::numeric_limits<Real>::min());

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real&              dx    = m_amr->getDx()[lvl];

    const Real vol = std::pow(dx, SpaceDim);

    for (DataIterator dit(dbl); dit.ok(); ++dit) {

      const EBISBox&       ebisbox    = ebisl[dit()];
      const BaseFab<bool>& validCells = (*m_amr->getValidCells(m_realm)[lvl])[dit()];

      const EBCellFAB& bgRate    = (*m_backgroundIonization[lvl])[dit()];
      const FArrayBox& bgRateReg = bgRate.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        if (validCells(iv) && ebisbox.isRegular(iv)) {
          for (size_t i = 0; i < m_voltageSweeps.size(); i++) {

            // Add background ionization rate.
            Rplus[i] += bgRateReg(iv, i) * vol;
            Rminu[i] += bgRateReg(iv, i) * vol;
          }
        }
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const IntVect iv = vof.gridIndex();
        if (validCells(iv) && ebisbox.isIrregular(iv)) {

          const Real volFrac = ebisbox.volFrac(vof);

          for (size_t i = 0; i < m_voltageSweeps.size(); i++) {

            // Add background ionization rate.
            Rplus[i] += bgRate(vof, i) * vol * volFrac;
            Rminu[i] += bgRate(vof, i) * vol * volFrac;
          }
        }
      };

      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      BoxLoops::loop(dbl[dit()], regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }

    for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
      m_firstElectronPlus.push_back(1. / ParallelOps::sum(Rplus[i]));
      m_firstElectronMinu.push_back(1. / ParallelOps::sum(Rminu[i]));
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::printReport() const noexcept
{
  CH_TIME("StreamerInceptionStepper::printReport");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::printReport" << endl;
  }

  // Compute the inception voltage.
#ifdef CH_MPI
  if (procID() == 0) {
#endif
    const Real Uinc = this->computeInceptionVoltage();

    std::cout << "# ===================================================================\n";
    if (std::isnan(Uinc)) {
      std::cout << "# Could not compute inception voltage\n";
    }
    else {
      std::cout << "# Inception voltage = " << Uinc << "\n";
    }
#ifdef CH_MPI
  }
#endif

#ifdef CH_MPI
  if (procID() == 0) {
#endif
    std::cout << "# ===================================================================\n";
    std::cout << left << setw(15) << setfill(' ') << "# Voltage";
    std::cout << left << setw(15) << setfill(' ') << "Max K";
    std::cout << left << setw(20) << setfill(' ') << "Crit. vol (m^3)";
    std::cout << left << setw(20) << setfill(' ') << "Time lag + (s)";
    std::cout << left << setw(20) << setfill(' ') << "Time lag - (s)" << std::endl;
    std::cout << "# ===================================================================\n";

    for (int i = 0; i < m_voltageSweeps.size(); i++) {
      std::cout << left << setw(15) << setfill(' ') << m_voltageSweeps[i];
      std::cout << left << setw(15) << setfill(' ') << m_maxKValues[i];
      std::cout << left << setw(20) << setfill(' ') << m_criticalVolumes[i];
      std::cout << left << setw(20) << setfill(' ') << m_firstElectronPlus[i];
      std::cout << left << setw(20) << setfill(' ') << m_firstElectronMinu[i];
      std::cout << endl;
    }
    std::cout << "# ===================================================================\n";
#ifdef CH_MPI
  }
#endif
}

template <typename P, typename F>
bool
StreamerInceptionStepper<P, F>::particleOutsideGrid(const RealVect& a_pos,
                                                    const RealVect& a_probLo,
                                                    const RealVect& a_probHi) const noexcept
{
#ifndef NDEBUG
  CH_TIME("StreamerInceptionStepper::particleOutsideGrid");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::particleOutsideGrid" << endl;
  }
#endif

  bool isOutside = false;

  for (int dir = 0; dir < SpaceDim; dir++) {
    if (a_pos[dir] < a_probLo[dir] || a_pos[dir] > a_probHi[dir]) {
      isOutside = true;
    }
  }

  return isOutside;
}

template <typename P, typename F>
bool
StreamerInceptionStepper<P, F>::particleInsideEB(const RealVect a_pos) const noexcept
{
#ifndef NDEBUG
  CH_TIME("StreamerInceptionStepper::particleInsideEB");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::particleInsideEB" << endl;
  }
#endif

  const RefCountedPtr<BaseIF>& implicitFunction = m_amr->getBaseImplicitFunction(m_phase);

  return (implicitFunction->value(a_pos) > 0.0) ? true : false;
}

#include <CD_NamespaceFooter.H>

#endif
