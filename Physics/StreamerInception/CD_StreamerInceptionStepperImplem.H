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
    return 0.0;
  };

  m_fieldEmission = [](const Real E) {
    return 0.0;
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
  m_amr->allocate(m_inceptionIntegralPlus, m_realm, m_phase, m_voltageSweeps.size());
  m_amr->allocate(m_inceptionIntegralMinu, m_realm, m_phase, m_voltageSweeps.size());
  m_amr->allocate(m_backgroundIonization, m_realm, m_phase, m_voltageSweeps.size());
  m_amr->allocate(m_emissionRatesPlus, m_realm, m_phase, m_voltageSweeps.size());
  m_amr->allocate(m_emissionRatesMinu, m_realm, m_phase, m_voltageSweeps.size());
  m_amr->allocate(m_inceptionVoltagePlus, m_realm, m_phase, 1);
  m_amr->allocate(m_inceptionVoltageMinu, m_realm, m_phase, 1);  
  m_amr->allocate(m_driftTimesPlus, m_realm, m_phase, m_voltageSweeps.size());
  m_amr->allocate(m_driftTimesMinu, m_realm, m_phase, m_voltageSweeps.size());

  DataOps::setValue(m_inceptionIntegralPlus, 0.0);
  DataOps::setValue(m_inceptionIntegralMinu, 0.0);
  DataOps::setValue(m_backgroundIonization, 0.0);
  DataOps::setValue(m_emissionRatesPlus, 0.0);
  DataOps::setValue(m_emissionRatesMinu, 0.0);
  DataOps::setValue(m_inceptionVoltagePlus, 0.0);
  DataOps::setValue(m_inceptionVoltageMinu, 0.0);  
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

  std::string str1;
  std::string str2;  

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
  pp.get("print_report", m_printReport);
  
  // Get the inception algorithm
  pp.get("inception_alg", str1, 0);
  pp.get("inception_alg", str2, 1);
  pp.get("inception_alg", stepSize, 2);    
  if (str1 == "euler") {
    m_inceptionAlgorithm = IntegrationAlgorithm::Euler;
  }
  else if (str1 == "trapz") {
    m_inceptionAlgorithm = IntegrationAlgorithm::Trapezoidal;
  }
  else {
    MayDay::Error("Expected 'euler' or 'trapz' for 'StreamerInceptionStepper.inception_alg'");    
  }
  
  if(str2 == "dx") {
    m_useFixedInceptionStep    = false;
    m_relativeInceptionStep = stepSize;
  }
  else if(str2 == "fixed") {
    m_useFixedInceptionStep    = true;
    m_fixedInceptionStep    = stepSize;
  }
  else {
    MayDay::Error("Expected 'dx' or 'fixed' for 'StreamerInceptionStepper.inception_alg'");
  }

  // Get the drift algorithm
  pp.get("drift_alg", str1, 0);
  pp.get("drift_alg", str2, 1);
  pp.get("drift_alg", stepSize, 2);    
  if (str1 == "euler") {
    m_driftAlgorithm = IntegrationAlgorithm::Euler;
  }
  else if (str1 == "trapz") {
    m_driftAlgorithm = IntegrationAlgorithm::Trapezoidal;
  }
  else {
    MayDay::Error("Expected 'euler' or 'trapz' for 'StreamerInceptionStepper.drift_alg'");    
  }
  
  if(str2 == "dx") {
    m_useFixedDriftStep    = false;
    m_relativeDriftStep = stepSize;
  }
  else if(str2 == "fixed") {
    m_useFixedDriftStep    = true;
    m_fixedDriftStep    = stepSize;
  }
  else {
    MayDay::Error("Expected 'dx' or 'fixed' for 'StreamerInceptionStepper.drift_alg'");
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

  const std::string prefix = "StreamerInceptionStepper/";

  // Append plot variable names.
  a_plotVariableNames.append(m_fieldSolver->getPlotVariableNames());
  a_plotVariableNames.append(m_particleSolver->getPlotVariableNames());

  // Write into plot data holder memory.
  m_fieldSolver->writePlotData(a_output, a_icomp);
  m_particleSolver->writePlotData(a_output, a_icomp);

  // Write the inception voltage.
  if (m_plotInceptionVoltage) {

    // Positive polarity
    DataOps::copy(a_output, m_inceptionVoltagePlus, Interval(a_icomp, a_icomp), Interval(0,0));
    a_icomp += 1;
    a_plotVariableNames.push_back(prefix + "Inception voltage +");

    // Negative polarity
    DataOps::copy(a_output, m_inceptionVoltageMinu, Interval(a_icomp, a_icomp), Interval(0,0));
    a_icomp += 1;
    a_plotVariableNames.push_back(prefix + "Inception voltage -");    
  }

  // Write our K values.
  if (m_plotInceptionIntegral) {

    // Positive polarity
    {
      const Interval srcInterval(0, m_voltageSweeps.size() - 1);
      const Interval dstInterval(a_icomp, a_icomp + m_voltageSweeps.size() - 1);      
      DataOps::copy(a_output, m_inceptionIntegralPlus, dstInterval, srcInterval);

      a_icomp += m_voltageSweeps.size();

      std::string varName;
      for (const Real& V : m_voltageSweeps) {
        varName = prefix + "K-value/ V = +";
        varName += std::to_string(V);
        a_plotVariableNames.push_back(varName);
      }
    }

    // Negative polarity
    {
      const Interval srcInterval(0, m_voltageSweeps.size() - 1);
      const Interval dstInterval(a_icomp, a_icomp + m_voltageSweeps.size() - 1);      
      DataOps::copy(a_output, m_inceptionIntegralPlus, dstInterval, srcInterval);      

      a_icomp += m_voltageSweeps.size();

      std::string varName;
      for (const Real& V : m_voltageSweeps) {
        varName = prefix + "K-value/ V = -";
        varName += std::to_string(V);
        a_plotVariableNames.push_back(varName);
      }
    }
  }

  // Write our drift time values.
  if (m_plotDriftTime) {

    // Positive polarity
    {
      const Interval srcInterval(0, m_voltageSweeps.size() - 1);
      const Interval dstInterval(a_icomp, a_icomp + m_voltageSweeps.size() - 1);      
      DataOps::copy(a_output, m_driftTimesPlus, dstInterval, srcInterval);
      
      a_icomp += m_voltageSweeps.size();

      std::string varName;
      for (const Real& V : m_voltageSweeps) {
        varName = prefix + "Drift time by mu/ V = +";
        varName += std::to_string(V);
        a_plotVariableNames.push_back(varName);
      }
    }

    // Negative polarity
    {
      const Interval srcInterval(0, m_voltageSweeps.size() - 1);
      const Interval dstInterval(a_icomp, a_icomp + m_voltageSweeps.size() - 1);      
      DataOps::copy(a_output, m_driftTimesMinu, dstInterval, srcInterval);      

      a_icomp += m_voltageSweeps.size();

      std::string varName;
      for (const Real& V : m_voltageSweeps) {
        varName = prefix + "Drift time by mu/ V = -";
        varName += std::to_string(V);
        a_plotVariableNames.push_back(varName);
      }
    }
  }

  // Write the background ionization rates
  if (m_plotBackgroundIonization) {
    const Interval srcInterval(0, m_voltageSweeps.size() - 1);
    const Interval dstInterval(a_icomp, a_icomp + m_voltageSweeps.size() - 1);
    DataOps::copy(a_output, m_backgroundIonization, dstInterval, srcInterval);

    a_icomp += m_voltageSweeps.size();

    std::string varName;
    for (const Real& V : m_voltageSweeps) {
      varName = prefix + "Background ionization rate/ V =  ";
      varName += std::to_string(V);
      a_plotVariableNames.push_back(varName);
    }
  }

  // Write field emission rates
  if (m_plotFieldEmission) {
    // Positive polarity
    {
      const Interval srcInterval(0, m_voltageSweeps.size() - 1);
      const Interval dstInterval(a_icomp, a_icomp + m_voltageSweeps.size() - 1);
      DataOps::copy(a_output, m_emissionRatesPlus, dstInterval, srcInterval);

      a_icomp += m_voltageSweeps.size();      

      std::string varName;
      for (const Real& V : m_voltageSweeps) {
        varName = prefix + "Field emission rate/ V = +";
        varName += std::to_string(V);
        a_plotVariableNames.push_back(varName);
      }
    }

    // Negative polarity
    {

      const Interval srcInterval(0, m_voltageSweeps.size() - 1);
      const Interval dstInterval(a_icomp, a_icomp + m_voltageSweeps.size() - 1);
      DataOps::copy(a_output, m_emissionRatesMinu, dstInterval, srcInterval);

      a_icomp += m_voltageSweeps.size();            

      std::string varName;
      for (const Real& V : m_voltageSweeps) {
        varName = prefix + "Field emission rate/ V = -";
        varName += std::to_string(V);
        a_plotVariableNames.push_back(varName);
      }
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
    ncomp += 2 * m_voltageSweeps.size();
  }

  // Drift times (+/-) direction
  if (m_plotDriftTime) {
    ncomp += 2 * m_voltageSweeps.size();
  }

  // Inception voltage.
  if (m_plotInceptionVoltage) {
    ncomp += 2;
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
const EBAMRCellData* const
StreamerInceptionStepper<P, F>::getElectricField() const noexcept
{
  CH_TIME("StreamerInceptionStepper::getElectricField");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::getElectricField" << endl;
  }
  return &m_electricField;
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

  m_amr->reallocate(m_inceptionIntegralPlus, m_phase, a_lmin);
  m_amr->reallocate(m_inceptionIntegralMinu, m_phase, a_lmin);
  m_amr->reallocate(m_backgroundIonization, m_phase, a_lmin);
  m_amr->reallocate(m_emissionRatesPlus, m_phase, a_lmin);
  m_amr->reallocate(m_emissionRatesMinu, m_phase, a_lmin);
  m_amr->reallocate(m_inceptionVoltagePlus, m_phase, a_lmin);
  m_amr->reallocate(m_inceptionVoltageMinu, m_phase, a_lmin);
  m_amr->reallocate(m_driftTimesPlus, m_phase, a_lmin);
  m_amr->reallocate(m_driftTimesMinu, m_phase, a_lmin);
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

      // Define particles' weights (i.e. integration & drift time tracker) and integration flag
      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        // Set mass to zero and do a backup of the initial position (we need to rewind particles
        // back later on). The scalar<0> flag is the flag we use for start-stop criteria for the
        // particle integration.
        p.weight()           = 0.0;
        p.template real<0>() = 1.0;
        p.template vect<0>() = p.position();
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
  if(m_plotDriftTime) {
    this->computeDriftTimes();
  }

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

  // Transient storage we can deposit particles onto.
  EBAMRCellData K;
  EBAMRCellData Kplus;
  EBAMRCellData Kminu;

  m_amr->allocate(K, m_realm, m_phase, 1);
  m_amr->allocate(Kplus, m_realm, m_phase, 1);
  m_amr->allocate(Kminu, m_realm, m_phase, 1);

  // Integrate for each voltage
  for (int i = 0; i < m_voltageSweeps.size(); ++i) {

    DataOps::setValue(K, 0.0);
    DataOps::setValue(Kplus, 0.0);
    DataOps::setValue(Kminu, 0.0);

    // Polarities. Note that we are dealing with electrons so
    // for positive polarity the particles move opposite to the field.
    std::vector<Real> polarities{-1.0, 1.0};

    for (const auto& p : polarities) {

      // Scale E-field by +/-
      DataOps::scale(m_electricField, p);

      this->resetParticles();

      // Switch between various integration algorithms. At the end of this
      // the K-value time should be stored on the particle weight and the particles
      // should be back in their original grid cells.
      switch (m_inceptionAlgorithm) {
      case IntegrationAlgorithm::Euler: {
        this->inceptionIntegrateEuler(m_voltageSweeps[i]);

        break;
      }
      case IntegrationAlgorithm::Trapezoidal: {
        this->inceptionIntegrateTrapezoidal(m_voltageSweeps[i]);

        break;
      }
      default: {
        MayDay::Error("StreamerInceptionStepper::computeInceptionIntegral -- logic bust");

        break;
      }
      }

      // Deposit particles on mesh and copy to relevant data holder.
      if (p < 0.0) {
        m_particleSolver->deposit(Kplus);

	m_amr->averageDown(Kplus, m_realm, m_phase);
	m_amr->interpGhost(Kplus, m_realm, m_phase);
	
        DataOps::copy(m_inceptionIntegralPlus, Kplus, Interval(i, i), Interval(0, 0));

	// Get max K-value
	Real maxK = -std::numeric_limits<Real>::max();
	Real minK = +std::numeric_limits<Real>::max();

	DataOps::getMaxMin(maxK, minK, Kplus, 0);

	m_maxKPlus.push_back(maxK);
      }
      else {
        m_particleSolver->deposit(Kminu);

	m_amr->averageDown(Kplus, m_realm, m_phase);
	m_amr->interpGhost(Kplus, m_realm, m_phase);
	
        DataOps::copy(m_inceptionIntegralMinu, Kminu, Interval(i, i), Interval(0, 0));

	// Get max K-value
	Real maxK = -std::numeric_limits<Real>::max();
	Real minK = +std::numeric_limits<Real>::max();

	DataOps::getMaxMin(maxK, minK, Kminu, 0);

	m_maxKMinu.push_back(maxK);	
      }

      // Scale field back again.
      DataOps::scale(m_electricField, p);
    }
  }
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

  // Allocate a data holder for holding the processed particles. This
  // will be faster because then we only have to iterate through the
  // particles that are actually still moving.
  ParticleContainer<P> amrProcessedParticles;
  m_amr->allocate(amrProcessedParticles, 0, m_realm);

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  m_particleSolver->remap();
  m_particleSolver->interpolateVelocities();

  while (amrParticles.getNumberOfValidParticlesGlobal() > 0) {

    // Euler integration.
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const Real dx   = m_amr->getDx()[lvl];
      const Real step = m_useFixedInceptionStep ? m_fixedInceptionStep : (m_relativeInceptionStep * dx);

      for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
        List<P>& solverParticles    = amrParticles[lvl][dit()].listItems();
        List<P>& processedParticles = amrProcessedParticles[lvl][dit()].listItems();

        for (ListIterator<P> lit(solverParticles); lit.ok();) {
          P& p = lit();

          const RealVect v     = p.velocity();
          const Real     dt    = step / v.vectorLength();
          const Real     alpha = m_alpha(a_voltage * v.vectorLength());

          // Update particle position
          const RealVect newPos = p.position() + dt * v;

          // Stop integration for particles that move into regions alpha < 0.0,
          // inside the EB or outside of the domain.
          const bool negativeAlpha = (alpha <= 0.0);
          const bool outsideDomain = this->particleOutsideGrid(newPos, probLo, probHi);
          const bool insideEB      = this->particleInsideEB(newPos);

          if (negativeAlpha || outsideDomain || insideEB) {
            processedParticles.transfer(lit);
          }
          else {
            p.position() = newPos;
            p.weight() += step * alpha;

            ++lit;
          }
        }
      }
    }

    // Update velocities.
    m_particleSolver->remap();
    m_particleSolver->interpolateVelocities();
  }

  ParticleOps::copyDestructive(amrParticles, amrProcessedParticles);

  this->rewindParticles();
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::inceptionIntegrateTrapezoidal(const Real& a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::inceptionIntegrateTrapezoidal");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::inceptionIntegrateTrapezoidal" << endl;
  }

  // TLDR: We move the particle using Heun's method.
  //
  //          p^(k+1) = p^k + 0.5 * dt * [v(p^k) + v(p^l)]
  //
  //       where p^l = p^k + dt * v(p^k). We will set dt = d/|v(p^k)|. Once we have the
  //       particle endpoints we can compute the contribution to the inception integral
  //       by using the trapezoidal (quadrature) rule
  //
  //          T += 0.5 * dx * [alpha(E(p^k)) + alpha(E(p^(k+1)))],
  //
  //       where dx = |p^(k+1) - p^k|.

  const RealVect probLo = m_amr->getProbLo();
  const RealVect probHi = m_amr->getProbHi();

  // Allocate a data holder for holding the processed particles. This
  // will be faster because then we only have to iterate through the
  // particles that are actually still moving.
  ParticleContainer<P> amrProcessedParticles;
  m_amr->allocate(amrProcessedParticles, 0, m_realm);

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  m_particleSolver->remap();
  m_particleSolver->interpolateVelocities();

  while (amrParticles.getNumberOfValidParticlesGlobal() > 0) {

    // Euler stage.
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const Real dx   = m_amr->getDx()[lvl];
      const Real step = m_useFixedInceptionStep ? m_fixedInceptionStep : (m_relativeInceptionStep * dx);      

      for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
        List<P>& solverParticles    = amrParticles[lvl][dit()].listItems();
        List<P>& processedParticles = amrProcessedParticles[lvl][dit()].listItems();

        for (ListIterator<P> lit(solverParticles); lit.ok();) {
          P& p = lit();

          const RealVect v     = p.velocity();
          const Real     dt    = step / v.vectorLength();
          const Real     alpha = m_alpha(a_voltage * v.vectorLength());

          // Storage alpha(p^k), v(p^k), and the step size.
          p.template real<0>() = alpha;
          p.template real<1>() = dt;
          p.template vect<1>() = v;

          // Update position.
          const RealVect newPos = p.position() + v * dt;

          // Stop integration for particles that move into regions alpha < 0.0,
          // inside the EB or outside of the domain.
          const bool negativeAlpha = (alpha <= 0.0);
          const bool outsideDomain = this->particleOutsideGrid(newPos, probLo, probHi);
          const bool insideEB      = this->particleInsideEB(newPos);

          if (negativeAlpha || outsideDomain || insideEB) {
            processedParticles.transfer(lit);
          }
          else {
            p.position() = newPos;

            ++lit;
          }
        }
      }
    }

    // Remap and update velocities.
    m_particleSolver->remap();
    m_particleSolver->interpolateVelocities();

    // Second stage.
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
        List<P>& solverParticles    = amrParticles[lvl][dit()].listItems();
        List<P>& processedParticles = amrProcessedParticles[lvl][dit()].listItems();

        for (ListIterator<P> lit(solverParticles); lit.ok();) {

          P& p = lit();

          const Real     dt  = p.template real<1>();
          const RealVect vk  = p.template vect<1>();
          const RealVect vk1 = p.velocity();

          // Note the weird subtraction since p.position() was updated
          // to p^k + dt * v^k.
          const RealVect newPos = p.position() + 0.5 * dt * (vk1 - vk);

          // Compute new alpha.
          const Real alpha1 = p.template real<0>();
          const Real alpha2 = m_alpha(a_voltage * vk1.vectorLength());
          const Real d      = 0.5 * dt * (vk + vk1).vectorLength();

          // Stop integration for particles that move into regions alpha < 0.0,
          // inside the EB or outside of the domain.
          const bool negativeAlpha = (alpha1 + alpha2) <= 0.0;
          const bool outsideDomain = this->particleOutsideGrid(newPos, probLo, probHi);
          const bool insideEB      = this->particleInsideEB(newPos);

          if (negativeAlpha || outsideDomain || insideEB) {
            processedParticles.transfer(lit);
          }
          else {
            p.position() = newPos;
            p.weight() += 0.5 * d * (alpha1 + alpha2);

            ++lit;
          }
        }
      }
    }

    // Remap and update velocities.
    m_particleSolver->remap();
    m_particleSolver->interpolateVelocities();
  }

  // Copy processed particles over to the solver particles.
  ParticleOps::copyDestructive(amrParticles, amrProcessedParticles);

  this->rewindParticles();
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
  //       volume. To figure out the critical volume, we interpolate it onto the particle's weight. If the particle
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

  // Transient storage we can deposit into.
  EBAMRCellData K;
  EBAMRCellData T;

  m_amr->allocate(T, m_realm, m_phase, 1);
  m_amr->allocate(K, m_realm, m_phase, 1);

  DataOps::setValue(T, 0.0);
  DataOps::setValue(K, 0.0);

  // Integrate for each voltage
  for (int i = 0; i < m_voltageSweeps.size(); ++i) {

    // Polarities. Note that we are dealing with negative ions so
    // for positive polarity the particles move opposite to the field.
    std::vector<Real> polarities{-1.0, 1.0};

    for (const auto& p : polarities) {

      this->resetParticles();

      // Scale E-field by +/-
      DataOps::scale(m_electricField, p);

      // Copy current K-value to the above data holder.
      if (p < 0.0) {
        DataOps::copy(K, m_inceptionIntegralPlus, Interval(0, 0), Interval(i, i));
      }
      else {
        DataOps::copy(K, m_inceptionIntegralMinu, Interval(0, 0), Interval(i, i));
      }

      // Switch between various integration algorithms. At the end of this
      // the drift time should be stored on the particle weight.
      switch (m_driftAlgorithm) {
      case IntegrationAlgorithm::Euler: {

        this->driftIntegrateEuler(K, m_voltageSweeps[i]);

        break;
      }
      case IntegrationAlgorithm::Trapezoidal: {

        this->driftIntegrateTrapezoidal(K, m_voltageSweeps[i]);

        break;
      }
      default: {
        MayDay::Error("StreamerInceptionStepper::computeDriftTimes -- logic bust");

        break;
      }
      }

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

  // Allocate a data holder for holding the processed particles. This
  // will be faster because then we only have to iterate through the
  // particles that are actually still moving.
  ParticleContainer<P> amrProcessedParticles;
  m_amr->allocate(amrProcessedParticles, 0, m_realm);

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  // Reset transient storage for particle integration. Need to use it because we interpolate
  // the K-value onto the particle weight, so the integration variable goes here.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
      for (ListIterator<P> lit(amrParticles[lvl][dit()].listItems()); lit.ok(); ++lit) {
        lit().template real<0>() = 0.0;
      }
    }
  }

  // Get K-value and velocity at particle position.
  m_particleSolver->remap();
  m_particleSolver->interpolateWeight(a_K);
  m_particleSolver->interpolateVelocities();

  while (amrParticles.getNumberOfValidParticlesGlobal() > 0) {

    // Euler kernel.
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const Real dx   = m_amr->getDx()[lvl];
      const Real step = m_useFixedDriftStep ? m_fixedDriftStep : (m_relativeDriftStep * dx);            

      for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
        List<P>& solverParticles    = amrParticles[lvl][dit()].listItems();
        List<P>& processedParticles = amrProcessedParticles[lvl][dit()].listItems();

        for (ListIterator<P> lit(solverParticles); lit.ok();) {
          P& p = lit();

          const RealVect v  = p.velocity();
          const Real     dt = step / v.vectorLength();

          // Update the particle position.
          const RealVect newPos = p.position() + dt * v;

          // Stop integration if particle ends up inside the EB, outside the domain, or if
          // the K-value is less than the inception criterion.
          const bool outsideGrid = this->particleOutsideGrid(newPos, probLo, probHi);
          const bool insideEB    = this->particleInsideEB(newPos);
          const bool insideCrit  = (p.weight() >= m_inceptionK);

          if (insideCrit) {
            p.weight() = p.template real<0>();

            processedParticles.transfer(lit);
          }
          else if (outsideGrid || insideEB) {
            p.weight() = std::numeric_limits<Real>::quiet_NaN();

            processedParticles.transfer(lit);
          }
          else {
            p.position() = newPos;

            p.template real<0>() += step / v.vectorLength();

            ++lit;
          }
        }
      }
    }

    amrProcessedParticles.remap();

    m_particleSolver->remap();
    m_particleSolver->interpolateWeight(a_K);
    m_particleSolver->interpolateVelocities();
  }

  ParticleOps::copyDestructive(amrParticles, amrProcessedParticles);

  this->rewindParticles();
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::driftIntegrateTrapezoidal(const EBAMRCellData& a_K, const Real& a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::driftIntegrateTrapezoidal");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::driftIntegrateTrapezoidal" << endl;
  }

  // TLDR: We move the particle using Heun's method.
  //
  //          p^(k+1) = p^k + 0.5 * dt * [v(p^k) + v(p^l)]
  //
  //       where p^l = p^k + dt * v(p^k). We will set dt = d/|v(p^k)|. Once we have the
  //       particle endpoints we can compute the contribution to the inception integral
  //       by using the trapezoidal (quadrature) rule
  //
  //          T += 0.5 * dx * [1/|v^k| + 1/|v^(k+)1|]
  //
  //       where dx = |p^(k+1) - p^k|.

  const RealVect probLo = m_amr->getProbLo();
  const RealVect probHi = m_amr->getProbHi();

  // Allocate a data holder for holding the processed particles. This
  // will be faster because then we only have to iterate through the
  // particles that are actually still moving.
  ParticleContainer<P> amrProcessedParticles;
  m_amr->allocate(amrProcessedParticles, 0, m_realm);

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  // Get K-value and velocity at particle position.
  m_particleSolver->remap();
  m_particleSolver->interpolateWeight(a_K);
  m_particleSolver->interpolateVelocities();

  // Reset transient storage for particle integration. Need to use it because we interpolate
  // the K-value onto the particle weight, so the integration variable goes here.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
      for (ListIterator<P> lit(amrParticles[lvl][dit()].listItems()); lit.ok(); ++lit) {
        lit().template real<0>() = 0.0;
      }
    }
  }

  while (amrParticles.getNumberOfValidParticlesGlobal() > 0) {

    // Euler stage.
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const Real dx   = m_amr->getDx()[lvl];
      const Real step = m_useFixedDriftStep ? m_fixedDriftStep : (m_relativeDriftStep * dx);                  

      for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
        List<P>& solverParticles    = amrParticles[lvl][dit()].listItems();
        List<P>& processedParticles = amrProcessedParticles[lvl][dit()].listItems();

        for (ListIterator<P> lit(solverParticles); lit.ok();) {
          P& p = lit();

          const RealVect v  = p.velocity();
          const Real     dt = step / v.vectorLength();

          // Store the above.
          p.template real<1>() = dt;
          p.template vect<1>() = v;

          // Update the particle position.
          const RealVect newPos = p.position() + dt * v;

          // Stop integration if particle ends up inside the EB, outside the domain, or if
          // the K-value is less than the inception criterion.
          const bool outsideGrid = this->particleOutsideGrid(newPos, probLo, probHi);
          const bool insideEB    = this->particleInsideEB(newPos);
          const bool insideCrit  = (p.weight() >= m_inceptionK);

          if (insideCrit) {
            p.weight() = p.template real<0>();

            processedParticles.transfer(lit);
          }
          else if (outsideGrid || insideEB) {
            p.weight() = std::numeric_limits<Real>::quiet_NaN();

            processedParticles.transfer(lit);
          }
          else {
            p.position() = newPos;

            ++lit;
          }
        }
      }
    }

    amrProcessedParticles.remap();

    m_particleSolver->remap();
    m_particleSolver->interpolateWeight(a_K);
    m_particleSolver->interpolateVelocities();

    // Second stage.
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {

      for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
        List<P>& solverParticles    = amrParticles[lvl][dit()].listItems();
        List<P>& processedParticles = amrProcessedParticles[lvl][dit()].listItems();

        for (ListIterator<P> lit(solverParticles); lit.ok();) {
          P& p = lit();

          const Real dt = p.template real<1>();

          const RealVect vk  = p.template vect<1>();
          const RealVect vk1 = p.velocity();

          // Update the particle position and compute distance moved.
          const RealVect newPos = p.position() + 0.5 * dt * (vk + vk1);
          const Real     d      = 0.5 * dt * (vk + vk1).vectorLength();

          // Stop integration if particle ends up inside the EB, outside the domain, or if
          // the K-value is less than the inception criterion.
          const bool outsideGrid = this->particleOutsideGrid(newPos, probLo, probHi);
          const bool insideEB    = this->particleInsideEB(newPos);
          const bool insideCrit  = (p.weight() >= m_inceptionK);

          if (insideCrit) {
            p.weight() = p.template real<0>();

            processedParticles.transfer(lit);
          }
          else if (outsideGrid || insideEB) {
            p.weight() = std::numeric_limits<Real>::quiet_NaN();

            processedParticles.transfer(lit);
          }
          else {
            p.position() = newPos;

            p.template real<0>() += 0.5 * d * (1. / vk.vectorLength() + 1. / vk1.vectorLength());

            ++lit;
          }
        }
      }
    }

    amrProcessedParticles.remap();    

    m_particleSolver->remap();
    m_particleSolver->interpolateWeight(a_K);
    m_particleSolver->interpolateVelocities();
  }

  ParticleOps::copyDestructive(amrParticles, amrProcessedParticles);

  this->rewindParticles();
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
StreamerInceptionStepper<P, F>::resetParticles() noexcept
{
  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
      for (ListIterator<P> lit(amrParticles[lvl][dit()].listItems()); lit.ok(); ++lit) {
        P& p = lit();

        p.weight()           = 0.0;
        p.template real<0>() = 1.0;
      }
    }
  }
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
    DataOps::setValue(m_inceptionVoltagePlus, std::numeric_limits<Real>::quiet_NaN());
    DataOps::setValue(m_inceptionVoltageMinu, std::numeric_limits<Real>::quiet_NaN());    

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
        EBCellFAB& inceptionVoltagePlus    = (*m_inceptionVoltagePlus[lvl])[dit()];
        EBCellFAB& inceptionVoltageMinu    = (*m_inceptionVoltageMinu[lvl])[dit()];
	
        FArrayBox& inceptionVoltagePlusReg = inceptionVoltagePlus.getFArrayBox();
        FArrayBox& inceptionVoltageMinuReg = inceptionVoltageMinu.getFArrayBox();	

        const EBCellFAB& inceptionIntegralPlus    = (*m_inceptionIntegralPlus[lvl])[dit()];
        const EBCellFAB& inceptionIntegralMinu    = (*m_inceptionIntegralMinu[lvl])[dit()];
	
        const FArrayBox& inceptionIntegralPlusReg = inceptionIntegralPlus.getFArrayBox();
        const FArrayBox& inceptionIntegralMinuReg = inceptionIntegralMinu.getFArrayBox();	

        // Regular kernel.
        auto regularKernel = [&](const IntVect& iv) -> void {
          std::vector<Real> Kplus;
          std::vector<Real> Kminu;
	  
          for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
            Kplus.emplace_back(inceptionIntegralPlusReg(iv, i));
            Kminu.emplace_back(inceptionIntegralMinuReg(iv, i));	    
          }

          inceptionVoltagePlusReg(iv, comp) = interpVinc(Kplus);
          inceptionVoltageMinuReg(iv, comp) = interpVinc(Kminu);	  
        };

        // Irregular kernel.
        auto irregularKernel = [&](const VolIndex& vof) -> void {
          std::vector<Real> Kplus;
          std::vector<Real> Kminu;
	  
          for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
            Kplus.emplace_back(inceptionIntegralPlus(vof, i));
            Kminu.emplace_back(inceptionIntegralMinu(vof, i));	    
          }

          inceptionVoltagePlus(vof, comp) = interpVinc(Kplus);
          inceptionVoltageMinu(vof, comp) = interpVinc(Kminu);	  
        };

        // Kernel regions.
        const Box&   cellBox = dbl[dit()];
        VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

        BoxLoops::loop(cellBox, regularKernel);
        BoxLoops::loop(vofit, irregularKernel);
      }
    }

    // Coarsen data.
    m_amr->averageDown(m_inceptionVoltagePlus, m_realm, m_phase);
    m_amr->interpGhost(m_inceptionVoltageMinu, m_realm, m_phase);
  }
}

template <typename P, typename F>
std::pair<Real, Real>
StreamerInceptionStepper<P, F>::computeInceptionVoltages() const noexcept
{
  CH_TIME("StreamerInceptionStepper::computeInceptionVoltages");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeInceptionVoltages" << endl;
  }

  Real Vplus = std::numeric_limits<Real>::quiet_NaN();
  Real Vminu = std::numeric_limits<Real>::quiet_NaN();  

  if (m_voltageSweeps.size() < 2) {
    MayDay::Warning("StreamerInceptionStepper::computeInceptionVoltageVolume -- not enough voltages for estimating "
                    "inception voltage");
  }
  else {

    constexpr int comp = 0;

    auto interpVinc = [&Kinc = this->m_inceptionK, &V = this->m_voltageSweeps](const std::vector<Real>& K) -> Real {
      Real Vinc = std::numeric_limits<Real>::quiet_NaN();
      
      for (size_t i = 0; i < K.size() - 1; i++) {

        if (K[i] <= Kinc && K[i + 1] > Kinc) {
          Vinc = V[i] + (Kinc - K[i]) * (V[i + 1] - V[i]) / (K[i + 1] - K[i]);

          break;
        }
      }

      return Vinc;
    };


    Vplus = interpVinc(m_maxKPlus);
    Vminu = interpVinc(m_maxKPlus);    
  }

  return std::make_pair(Vplus, Vminu);
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
    Real criticalVolumePlus = 0.0;
    Real criticalVolumeMinu = 0.0;    

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); ++lvl) {
      const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

      const LevelData<BaseFab<bool>>& validCellsLD = *m_amr->getValidCells(m_realm)[lvl];

      const Real dx = m_amr->getDx()[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
        const EBISBox&       ebisbox    = ebisl[dit()];
        const BaseFab<bool>& validCells = validCellsLD[dit()];

        const EBCellFAB& inceptionIntegralPlus = (*m_inceptionIntegralPlus[lvl])[dit()];
        const EBCellFAB& inceptionIntegralMinu = (*m_inceptionIntegralMinu[lvl])[dit()];
	
        const FArrayBox& inceptionIntegralPlusReg = inceptionIntegralPlus.getFArrayBox();
        const FArrayBox& inceptionIntegralMinuReg = inceptionIntegralMinu.getFArrayBox();	

        auto regularKernel = [&](const IntVect& iv) -> void {
	  if(validCells(iv, 0) && ebisbox.isRegular(iv)) {
	    if (inceptionIntegralPlusReg(iv, i) >= m_inceptionK) {
	      criticalVolumePlus += std::pow(dx, SpaceDim);
	    }
	    if (inceptionIntegralMinuReg(iv, i) >= m_inceptionK) {
	      criticalVolumeMinu += std::pow(dx, SpaceDim);
	    }
	  }
        };

        auto irregularKernel = [&](const VolIndex& vof) -> void {
	  if(validCells(vof.gridIndex())) {
	    const Real kappa = ebisbox.volFrac(vof);
	    
	    if (inceptionIntegralPlus(vof, i) >= m_inceptionK) {
	      criticalVolumePlus += kappa * std::pow(dx, SpaceDim);
	    }
	    if (inceptionIntegralMinu(vof, i) >= m_inceptionK) {
	      criticalVolumeMinu += kappa * std::pow(dx, SpaceDim);
	    }
	  }
        };

        // Kernel regions.
        const Box    cellBox = dbl[dit()];
        VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

        BoxLoops::loop(cellBox, regularKernel);
        BoxLoops::loop(vofit, irregularKernel);
      }
    }

    m_criticalVolumePlus.push_back(ParallelOps::sum(criticalVolumePlus));
    m_criticalVolumeMinu.push_back(ParallelOps::sum(criticalVolumeMinu));    
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
    const std::pair<Real, Real> Uinc = this->computeInceptionVoltages();

    std::cout << "# ===================================================================\n";
    if (std::isnan(Uinc.first) || std::isnan(Uinc.second)) {
      std::cout << "# Could not compute inception voltage\n";
    }
    else {
      std::cout << "# Inception voltage +/- = " << Uinc.first << "/" << Uinc.second << "\n";
    }
#ifdef CH_MPI
  }
#endif

#ifdef CH_MPI
  if (procID() == 0) {
#endif
    std::cout << "# ===================================================================\n";
    std::cout << left << setw(15) << setfill(' ') << "# +/-V";
    std::cout << left << setw(15) << setfill(' ') << "Max K(+)";
    std::cout << left << setw(15) << setfill(' ') << "Max K(-)";    
    std::cout << left << setw(20) << setfill(' ') << "Crit. vol(+) (m^3)";
    std::cout << left << setw(20) << setfill(' ') << "Crit. vol(-) (m^3)";    
    std::cout << left << setw(20) << setfill(' ') << "Time lag + (s)";
    std::cout << left << setw(20) << setfill(' ') << "Time lag - (s)" << std::endl;
    std::cout << "# ===================================================================\n";

    for (int i = 0; i < m_voltageSweeps.size(); i++) {
      std::cout << left << setw(15) << setfill(' ') << m_voltageSweeps[i];
      std::cout << left << setw(15) << setfill(' ') << m_maxKPlus[i];
      std::cout << left << setw(15) << setfill(' ') << m_maxKMinu[i];      
      std::cout << left << setw(20) << setfill(' ') << m_criticalVolumePlus[i];
      std::cout << left << setw(20) << setfill(' ') << m_criticalVolumeMinu[i];      
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
    if (a_pos[dir] <= a_probLo[dir] || a_pos[dir] >= a_probHi[dir]) {
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

  return (implicitFunction->value(a_pos) >= 0.0) ? true : false;
}

#include <CD_NamespaceFooter.H>

#endif
