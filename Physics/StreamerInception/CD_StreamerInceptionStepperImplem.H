/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Copyright © 2022 NTNU.
 * Copyright © 2022 Fanny Skirbekk. 
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_StreamerInceptionStepperImplem.H
  @brief  Implementation of CD_StreamerInceptionStepper.H
  @author Robert Marskar
  @author Rasmus Hoholm
  @author Fanny Skirbekk
*/

#ifndef CD_StreamerInceptionStepperImplem_H
#define CD_StreamerInceptionStepperImplem_H

// Std includes
#include <iostream>
#include <fstream>

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_Units.H>
#include <CD_StreamerInceptionSpecies.H>
#include <CD_StreamerInceptionStepper.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::StreamerInception;

template <typename P, typename F, typename C>
StreamerInceptionStepper<P, F, C>::StreamerInceptionStepper()
{
  CH_TIME("StreamerInceptionStepper::StreamerInceptionStepper");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::StreamerInceptionStepper" << endl;
  }

  // Default settings.
  m_verbosity          = -1;
  m_realm              = Realm::Primal;
  m_phase              = phase::gas;
  m_mode               = Mode::Stationary;
  m_transportAlgorithm = TransportAlgorithm::ImExCTU;

  this->parseOptions();

  m_rho = [](const RealVect x) -> Real {
    return 0.0;
  };

  m_alpha = [](const Real E) -> Real {
    return 0.0;
  };

  m_eta = [](const Real E) -> Real {
    return 0.0;
  };

  m_voltageCurve = [](const Real a_time) -> Real {
    return 1.0;
  };

  m_backgroundRate = [](const Real E, const RealVect x) -> Real {
    return 0.0;
  };

  m_detachmentRate = [](const Real E, const RealVect x) -> Real {
    return 0.0;
  };

  m_fieldEmission = [](const Real E, const RealVect x) -> Real {
    return 0.0;
  };

  m_initialIonDensity = [](const RealVect x) -> Real {
    return 0.0;
  };

  m_ionMobility = [](const Real E) -> Real {
    return 0.0;
  };

  m_ionDiffusion = [](const Real E) -> Real {
    return 0.0;
  };
}

template <typename P, typename F, typename C>
StreamerInceptionStepper<P, F, C>::~StreamerInceptionStepper()
{
  CH_TIME("StreamerInceptionStepper::~StreamerInceptionStepper");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::~StreamerInceptionStepper" << endl;
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::setupSolvers()
{
  CH_TIME("StreamerInceptionStepper::setupSolvers");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setupSolvers" << endl;
  }

  // Always solve for the field using a voltage of one (and then scale up/down later on)
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
  m_tracerParticleSolver = RefCountedPtr<TracerParticleSolver<P>>(new TracerParticleSolver<P>());
  m_tracerParticleSolver->parseOptions();
  m_tracerParticleSolver->setAmr(m_amr);
  m_tracerParticleSolver->setComputationalGeometry(m_computationalGeometry);
  m_tracerParticleSolver->setRealm(m_realm);
  m_tracerParticleSolver->setPhase(m_phase);
  m_tracerParticleSolver->setTime(0, 0.0, 0.0);
  m_tracerParticleSolver->setName("Tracer particle solver");
  m_tracerParticleSolver->setVolumeScale(true);

  // Instantiate the ion solver.
  m_ionSolver = RefCountedPtr<CdrSolver>(new C());

  m_ionSolver->parseOptions();
  m_ionSolver->setPhase(m_phase);
  m_ionSolver->setAmr(m_amr);
  m_ionSolver->setComputationalGeometry(m_computationalGeometry);
  m_ionSolver->setRealm(m_realm);

  auto species = RefCountedPtr<CdrSpecies>(new StreamerInceptionSpecies(m_initialIonDensity, true, true));

  m_ionSolver->setSpecies(species);
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::allocate()
{
  CH_TIME("StreamerInceptionStepper::allocate");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::allocate" << endl;
  }

  m_fieldSolver->allocate();
  m_tracerParticleSolver->allocate();
  m_ionSolver->allocate();

  // Needed for all modes.
  m_amr->allocate(m_electricFieldPlus, m_realm, m_phase, SpaceDim);
  m_amr->allocate(m_electricFieldMinu, m_realm, m_phase, SpaceDim);

  switch (m_mode) {
  case Mode::Stationary: {
    m_amr->allocate(m_inceptionIntegralPlus, m_realm, m_phase, m_voltageSweeps.size());
    m_amr->allocate(m_inceptionIntegralMinu, m_realm, m_phase, m_voltageSweeps.size());
    m_amr->allocate(m_inceptionVoltagePlus, m_realm, m_phase, 1);
    m_amr->allocate(m_inceptionVoltageMinu, m_realm, m_phase, 1);

    if (m_plotFieldEmission) {
      m_amr->allocate(m_emissionRatesPlus, m_realm, m_phase, m_voltageSweeps.size());
      m_amr->allocate(m_emissionRatesMinu, m_realm, m_phase, m_voltageSweeps.size());

      DataOps::setValue(m_emissionRatesPlus, 0.0);
      DataOps::setValue(m_emissionRatesMinu, 0.0);
    }

    if (m_plotBackgroundIonization) {
      m_amr->allocate(m_backgroundIonization, m_realm, m_phase, m_voltageSweeps.size());

      DataOps::setValue(m_backgroundIonization, 0.0);
    }

    if (m_plotDetachment) {
      m_amr->allocate(m_detachment, m_realm, m_phase, m_voltageSweeps.size());

      DataOps::setValue(m_detachment, 0.0);
    }

    DataOps::setValue(m_inceptionIntegralPlus, 0.0);
    DataOps::setValue(m_inceptionIntegralMinu, 0.0);
    DataOps::setValue(m_inceptionVoltagePlus, 0.0);
    DataOps::setValue(m_inceptionVoltageMinu, 0.0);

    break;
  }
  case Mode::Transient: {
    m_amr->allocate(m_inceptionIntegral, m_realm, m_phase, 1);
    m_amr->allocate(m_emissionRate, m_realm, m_phase, 1);
    m_amr->allocate(m_backgroundIonization, m_realm, m_phase, 1);
    m_amr->allocate(m_detachment, m_realm, m_phase, 1);

    DataOps::setValue(m_inceptionIntegral, 0.0);
    DataOps::setValue(m_emissionRate, 0.0);
    DataOps::setValue(m_backgroundIonization, 0.0);

    break;
  }
  default: {
    break;
  }
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::initialData()
{
  CH_TIME("StreamerInceptionStepper::initialData");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::initialData" << endl;
  }

  // Set rho = sigma = 0 for field solver and then solve the Poisson equation.
  m_fieldSolver->setRho(m_rho);
  m_fieldSolver->setSigma(0.0);

  const bool converged =
    m_fieldSolver->solve(m_fieldSolver->getPotential(), m_fieldSolver->getRho(), m_fieldSolver->getSigma(), false);

  if (!converged) {
    MayDay::Warning("StreamerInceptionStepper::initialData -- could not solve Poisson equation. ");
  }

  // Get the gas-side electric field and give it to the particle solvers.
  DataOps::copy(m_electricFieldPlus, m_amr->alias(phase::gas, m_fieldSolver->getElectricField()));

  m_amr->conservativeAverage(m_electricFieldPlus, m_realm, m_phase);
  m_amr->interpGhostMG(m_electricFieldPlus, m_realm, m_phase);

  m_amr->interpToCentroids(m_electricFieldPlus, m_realm, m_phase);

  // Also store -E
  DataOps::copy(m_electricFieldMinu, m_electricFieldPlus);
  DataOps::scale(m_electricFieldMinu, -1.0);
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::registerRealms()
{
  CH_TIME("StreamerInceptionStepper::registerRealms");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::registerRealms" << endl;
  }

  m_amr->registerRealm(m_realm);
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::registerOperators()
{
  CH_TIME("StreamerInceptionStepper::registerOperators");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::registerOperators" << endl;
  }

  m_fieldSolver->registerOperators();
  m_tracerParticleSolver->registerOperators();
  m_ionSolver->registerOperators();
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::parseOptions()
{
  CH_TIME("StreamerInceptionStepper::parseOptions");

  this->parseVerbosity();
  this->parseMode();
  this->parseVoltages();
  this->parseOutput();
  this->parsePlotVariables();
  this->parseInceptionAlgorithm();
  this->parseTransportAlgorithm();
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::parseRuntimeOptions()
{
  CH_TIME("StreamerInceptionStepper::parseRuntimeOptions");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::parseRuntimeOptions" << endl;
  }

  this->parseVerbosity();
  this->parsePlotVariables();
  this->parseInceptionAlgorithm();
  this->parseTransportAlgorithm();
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::parseVerbosity() noexcept
{
  CH_TIME("StreamerInceptionStepper::parseVerbosity");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::parseVerbosity" << endl;
  }

  ParmParse pp("StreamerInceptionStepper");
  pp.get("verbosity", m_verbosity);
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::parseMode() noexcept
{
  CH_TIME("StreamerInceptionStepper::parseMode");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::parseMode" << endl;
  }

  ParmParse pp("StreamerInceptionStepper");

  std::string str;

  pp.get("mode", str);
  if (str == "stationary") {
    m_mode = Mode::Stationary;
  }
  else if (str == "transient") {
    m_mode = Mode::Transient;
  }
  else {
    MayDay::Error("Expected 'none', 'stationary', or 'transient' for 'StreamerInceptionStepper.mode'");
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::parseVoltages() noexcept
{
  CH_TIME("StreamerInceptionStepper::parseVoltages");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::parseVoltages" << endl;
  }

  ParmParse pp("StreamerInceptionStepper");

  std::string str;

  Real voltageLo    = 0.0;
  Real voltageHi    = 0.0;
  Real stepSize     = 0.0;
  int  voltageSteps = 0;

  pp.get("voltage_lo", voltageLo);
  pp.get("voltage_hi", voltageHi);
  pp.get("voltage_steps", voltageSteps);
  pp.get("K_inception", m_inceptionK);

  // Define voltages
  const Real dV = (voltageHi - voltageLo) / (1 + voltageSteps);

  m_voltageSweeps.push_back(voltageLo);

  for (size_t i = 0; i < voltageSteps; i++) {
    m_voltageSweeps.push_back(m_voltageSweeps[i] + dV);
  }

  m_voltageSweeps.push_back(voltageHi);
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::parseOutput() noexcept
{
  CH_TIME("StreamerInceptionStepper::parseOutput");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::parseOutput" << endl;
  }

  ParmParse pp("StreamerInceptionStepper");

  pp.get("output_file", m_outputFile);
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::parseInceptionAlgorithm() noexcept
{
  CH_TIME("StreamerInceptionStepper::parseInceptionAlgorithm");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::parseInceptionAlgorithm" << endl;
  }

  ParmParse pp("StreamerInceptionStepper");

  std::string str;

  // Get the inception algorithm
  pp.get("inception_alg", str, 0);
  if (str == "euler") {
    m_inceptionAlgorithm = IntegrationAlgorithm::Euler;
  }
  else if (str == "trapz") {
    m_inceptionAlgorithm = IntegrationAlgorithm::Trapezoidal;
  }
  else {
    MayDay::Error("Expected 'euler' or 'trapz' for 'StreamerInceptionStepper.inception_alg'");
  }

  // Get the step size selection
  Real stepSize;
  pp.get("inception_alg", str, 1);
  pp.get("inception_alg", stepSize, 2);
  if (str == "dx") {
    m_useFixedInceptionStep = false;
    m_relativeInceptionStep = stepSize;
  }
  else if (str == "fixed") {
    m_useFixedInceptionStep = true;
    m_fixedInceptionStep    = stepSize;
  }
  else {
    MayDay::Error("Expected 'dx' or 'fixed' for 'StreamerInceptionStepper.inception_alg'");
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::parseTransportAlgorithm() noexcept
{
  CH_TIME("StreamerInceptionStepper::parseTransportAlgorithm");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::parseTransportAlgorithm" << endl;
  }

  ParmParse pp("StreamerInceptionStepper");

  std::string str;

  pp.get("transport_alg", str);
  pp.get("ion_transport", m_ionTransport);
  pp.get("cfl", m_cfl);
  pp.get("max_dt", m_maxDt);
  pp.get("min_dt", m_minDt);
  pp.get("delta_v_percent", m_deltaVPercent);
  pp.get("end_time", m_endTime);

  m_minDt = std::max(m_minDt, 0.0);

  if (str == "euler") {
    m_transportAlgorithm == TransportAlgorithm::Euler;
  }
  else if (str == "heun") {
    m_transportAlgorithm == TransportAlgorithm::Heun;
  }
  else if (str == "imex") {
    m_transportAlgorithm == TransportAlgorithm::ImExCTU;
  }
  else {
    MayDay::Error("Expected 'euler', 'heun', or 'imex' for 'StreamerInceptionStepper.transport_alg'");
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::parsePlotVariables() noexcept
{
  CH_TIME("StreamerInceptionStepper::parsePlotVariables");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::parsePlotVariables" << endl;
  }

  m_plotPoisson              = false;
  m_plotTracer               = false;
  m_plotNegativeIons         = false;
  m_plotInceptionIntegral    = false;
  m_plotInceptionVoltage     = false;
  m_plotBackgroundIonization = false;
  m_plotDetachment           = false;
  m_plotFieldEmission        = false;
  m_plotAlpha                = false;
  m_plotEta                  = false;

  ParmParse pp("StreamerInceptionStepper");

  // Get plot variables.
  const int num = pp.countval("plt_vars");
  if (num > 0) {
    Vector<std::string> plotVars(num);
    pp.getarr("plt_vars", plotVars, 0, num);

    for (int i = 0; i < num; i++) {
      if (plotVars[i] == "poisson") {
        m_plotPoisson = true;
      }
      else if (plotVars[i] == "tracer") {
        m_plotTracer = true;
      }
      else if (plotVars[i] == "ions") {
        m_plotNegativeIons = true;
      }
      else if (plotVars[i] == "K") {
        m_plotInceptionIntegral = true;
      }
      else if (plotVars[i] == "Uinc") {
        m_plotInceptionVoltage = true;
      }
      else if (plotVars[i] == "bg_rate") {
        m_plotBackgroundIonization = true;
      }
      else if (plotVars[i] == "detachment") {
        m_plotDetachment = true;
      }
      else if (plotVars[i] == "emission") {
        m_plotFieldEmission = true;
      }
      else if (plotVars[i] == "alpha") {
        m_plotAlpha = true;
      }
      else if (plotVars[i] == "eta") {
        m_plotEta = true;
      }
    }
  }
}

#ifdef CH_USE_HDF5
template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const
{
  CH_TIME("StreamerInceptionStepper::writeCheckpointData");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::writeCheckpointData" << endl;
  }
}
#endif

#ifdef CH_USE_HDF5
template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::readCheckpointData(HDF5Handle& a_handle, const int a_lvl)
{
  CH_TIME("StreamerInceptionStepper::readCheckpointData");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::readCheckpointData" << endl;
  }

  MayDay::Error("StreamerInceptionStepper::readCheckpointData -- restart not supported. Use Driver.restart=0");
}
#endif

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::writePlotData(EBAMRCellData&       a_output,
                                                 Vector<std::string>& a_plotVariableNames,
                                                 int&                 a_icomp) const
{
  CH_TIME("StreamerInceptionStepper::writePlotData");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::writePlotData" << endl;
  }

  if (m_plotPoisson) {
    m_fieldSolver->writePlotData(a_output, a_icomp);
    a_plotVariableNames.append(m_fieldSolver->getPlotVariableNames());
  }

  if (m_plotTracer) {
    m_tracerParticleSolver->writePlotData(a_output, a_icomp);
    a_plotVariableNames.append(m_tracerParticleSolver->getPlotVariableNames());
  }

  if (m_plotNegativeIons) {
    m_ionSolver->writePlotData(a_output, a_icomp);
    a_plotVariableNames.append(m_ionSolver->getPlotVariableNames());
  }

  // Write data -- this differs between the 'stationary' and 'transient' modes
  switch (m_mode) {
  case Mode::Stationary: {
    this->writePlotDataStationary(a_output, a_plotVariableNames, a_icomp);

    break;
  }
  case Mode::Transient:
    this->writePlotDataTransient(a_output, a_plotVariableNames, a_icomp);

    break;
  default: {
    break;
  }
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::writePlotDataStationary(EBAMRCellData&       a_output,
                                                           Vector<std::string>& a_plotVariableNames,
                                                           int&                 a_icomp) const noexcept
{
  CH_TIME("StreamerInceptionStepper::writePlotDataStationary");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::writePlotDataStationary" << endl;
  }

  const std::string prefix = "StreamerInceptionStepper/";

  // Write the inception voltage.
  if (m_plotInceptionVoltage) {

    // Positive polarity
    DataOps::copy(a_output, m_inceptionVoltagePlus, Interval(a_icomp, a_icomp), Interval(0, 0));
    a_icomp += 1;
    a_plotVariableNames.push_back(prefix + "Inception voltage +");

    // Negative polarity
    DataOps::copy(a_output, m_inceptionVoltageMinu, Interval(a_icomp, a_icomp), Interval(0, 0));
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
      DataOps::copy(a_output, m_inceptionIntegralMinu, dstInterval, srcInterval);

      a_icomp += m_voltageSweeps.size();

      std::string varName;
      for (const Real& V : m_voltageSweeps) {
        varName = prefix + "K-value/ V = -";
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

  // Write detachment rates
  if (m_plotDetachment) {
    const Interval srcInterval(0, m_voltageSweeps.size() - 1);
    const Interval dstInterval(a_icomp, a_icomp + m_voltageSweeps.size() - 1);
    DataOps::copy(a_output, m_detachment, dstInterval, srcInterval);

    a_icomp += m_voltageSweeps.size();

    std::string varName;
    for (const Real& V : m_voltageSweeps) {
      varName = prefix + "Detachment rate/ V =  ";
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

  // Alpha value
  if (m_plotAlpha) {
    EBAMRCellData alpha;
    m_amr->allocate(alpha, m_realm, m_phase, 1);

    for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
      this->evaluateFunction(alpha, m_voltageSweeps[i], m_alpha);

      m_amr->conservativeAverage(alpha, m_realm, m_phase);
      m_amr->interpGhost(alpha, m_realm, m_phase);

      DataOps::copy(a_output, alpha, Interval(a_icomp, a_icomp), Interval(0, 0));

      a_icomp++;

      const std::string varName = prefix + "Alpha coefficient/ V = " + std::to_string(m_voltageSweeps[i]);

      a_plotVariableNames.push_back(varName);
    }
  }

  // Eta
  if (m_plotEta) {
    EBAMRCellData eta;
    m_amr->allocate(eta, m_realm, m_phase, 1);

    for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
      this->evaluateFunction(eta, m_voltageSweeps[i], m_eta);

      m_amr->conservativeAverage(eta, m_realm, m_phase);
      m_amr->interpGhost(eta, m_realm, m_phase);

      DataOps::copy(a_output, eta, Interval(a_icomp, a_icomp), Interval(0, 0));

      a_icomp++;

      const std::string varName = prefix + "Eta coefficient/ V = " + std::to_string(m_voltageSweeps[i]);

      a_plotVariableNames.push_back(varName);
    }
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::writePlotDataTransient(EBAMRCellData&       a_output,
                                                          Vector<std::string>& a_plotVariableNames,
                                                          int&                 a_icomp) const noexcept
{
  CH_TIME("StreamerInceptionStepper::writePlotDataTransient");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::writePlotDataTransient" << endl;
  }

  // Add potential and electric field to output
  {
    EBAMRCellData potential;
    EBAMRCellData E;

    m_amr->allocate(potential, m_realm, m_phase, 1);
    m_amr->allocate(E, m_realm, m_phase, SpaceDim);

    int comp1 = 0;
    int comp2 = 0;
    m_fieldSolver->writeMultifluidData(potential, comp1, m_fieldSolver->getPotential(), phase::gas, true);
    m_fieldSolver->writeMultifluidData(E, comp2, m_fieldSolver->getElectricField(), phase::gas, true);

    DataOps::scale(potential, m_voltageCurve(m_time));
    DataOps::scale(E, m_voltageCurve(m_time));

    // Copy potential to output data.
    DataOps::copy(a_output, potential, Interval(a_icomp, a_icomp), Interval(0, 0));

    a_plotVariableNames.push_back("Electric potential");
    a_icomp += 1;

    // Copy electric field to output data
    const Interval srcInterv(0, SpaceDim - 1);
    const Interval dstInterv(a_icomp, a_icomp + SpaceDim - 1);
    DataOps::copy(a_output, E, dstInterv, srcInterv);

    a_plotVariableNames.push_back("x-E");
    a_plotVariableNames.push_back("y-E");
    if (SpaceDim == 3) {
      a_plotVariableNames.push_back("z-E");
    }

    a_icomp += SpaceDim;
  }

  // Plot inception integral
  if (m_plotInceptionIntegral) {
    DataOps::copy(a_output, m_inceptionIntegral, Interval(a_icomp, a_icomp), Interval(0, 0));

    a_icomp++;

    a_plotVariableNames.push_back("Inception integral");
  }

  // Plot background ionization rate
  if (m_plotBackgroundIonization) {
    EBAMRCellData bgIonization;

    m_amr->allocate(bgIonization, m_realm, m_phase, 1);

    this->evaluateFunction(bgIonization, m_voltageCurve(m_time), m_backgroundRate);

    m_amr->conservativeAverage(bgIonization, m_realm, m_phase);
    m_amr->interpGhost(bgIonization, m_realm, m_phase);

    DataOps::copy(a_output, bgIonization, Interval(a_icomp, a_icomp), Interval(0, 0));

    a_icomp++;

    a_plotVariableNames.push_back("Background rate");
  }

  // Plot electron detachment rate
  if (m_plotDetachment) {
    EBAMRCellData detachment;

    m_amr->allocate(detachment, m_realm, m_phase, 1);

    this->evaluateFunction(detachment, m_voltageCurve(m_time), m_detachmentRate);

    const EBAMRCellData& phi = m_ionSolver->getPhi();

    DataOps::multiply(detachment, phi);

    m_amr->conservativeAverage(detachment, m_realm, m_phase);
    m_amr->interpGhost(detachment, m_realm, m_phase);

    DataOps::copy(a_output, detachment, Interval(a_icomp, a_icomp), Interval(0, 0));

    a_icomp++;

    a_plotVariableNames.push_back("Detachment rate");
  }

  if (m_plotFieldEmission) {
    DataOps::copy(a_output, m_emissionRate, Interval(a_icomp, a_icomp), Interval(0, 0));

    a_icomp++;

    a_plotVariableNames.push_back("Field emission");
  }

  // Alpha value
  if (m_plotAlpha) {
    EBAMRCellData alpha;
    m_amr->allocate(alpha, m_realm, m_phase, 1);

    this->evaluateFunction(alpha, m_voltageCurve(m_time), m_alpha);

    m_amr->conservativeAverage(alpha, m_realm, m_phase);
    m_amr->interpGhost(alpha, m_realm, m_phase);

    DataOps::copy(a_output, alpha, Interval(a_icomp, a_icomp), Interval(0, 0));

    a_icomp++;

    a_plotVariableNames.push_back("Townsend alpha coefficient");
  }

  // Eta
  if (m_plotEta) {
    EBAMRCellData eta;
    m_amr->allocate(eta, m_realm, m_phase, 1);

    this->evaluateFunction(eta, m_voltageCurve(m_time), m_eta);

    m_amr->conservativeAverage(eta, m_realm, m_phase);
    m_amr->interpGhost(eta, m_realm, m_phase);

    DataOps::copy(a_output, eta, Interval(a_icomp, a_icomp), Interval(0, 0));

    a_icomp++;

    a_plotVariableNames.push_back("Townsend eta coefficient");
  }
}

template <typename P, typename F, typename C>
int
StreamerInceptionStepper<P, F, C>::getNumberOfPlotVariables() const
{
  CH_TIME("StreamerInceptionStepper::getNumberOfPlotVariables");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::getNumberOfPlotVariables" << endl;
  }

  int ncomp = 0;

  if (m_plotPoisson) {
    ncomp += m_fieldSolver->getNumberOfPlotVariables();
  }
  if (m_plotTracer) {
    ncomp += m_tracerParticleSolver->getNumberOfPlotVariables();
  }
  if (m_plotNegativeIons) {
    ncomp += m_ionSolver->getNumberOfPlotVariables();
  }

  switch (m_mode) {
  case Mode::Stationary: {
    // K-values
    if (m_plotInceptionIntegral) {
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

    // Detachment rates
    if (m_plotDetachment) {
      ncomp += m_voltageSweeps.size();
    }

    // Field emission rates
    if (m_plotFieldEmission) {
      ncomp += 2 * m_voltageSweeps.size();
    }

    // Plotting alpha
    if (m_plotAlpha) {
      ncomp += m_voltageSweeps.size();
    }

    // Plotting eta
    if (m_plotEta) {
      ncomp += m_voltageSweeps.size();
    }

    break;
  }
  case Mode::Transient: {

    // Always write the electric field and potential.
    ncomp += SpaceDim + 1;

    if (m_plotInceptionIntegral) {
      ncomp += 1;
    }
    if (m_plotBackgroundIonization) {
      ncomp += 1;
    }
    if (m_plotDetachment) {
      ncomp += 1;
    }
    if (m_plotFieldEmission) {
      ncomp += 1;
    }
    if (m_plotAlpha) {
      ncomp += 1;
    }
    if (m_plotEta) {
      ncomp += 1;
    }

    break;
  }
  default: {
    break;
  }
  }

  return ncomp;
}

template <typename P, typename F, typename C>
const EBAMRCellData* const
StreamerInceptionStepper<P, F, C>::getElectricField() const noexcept
{
  CH_TIME("StreamerInceptionStepper::getElectricField");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::getElectricField" << endl;
  }

  return &m_electricFieldPlus;
}

template <typename P, typename F, typename C>
Real
StreamerInceptionStepper<P, F, C>::computeDt()
{
  CH_TIME("StreamerInceptionStepper::computeDt");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeDt" << endl;
  }

  Real dt = std::numeric_limits<Real>::max();

  if (m_mode == Mode::Transient) {

    switch (m_transportAlgorithm) {
    case TransportAlgorithm::Euler: {
      dt = m_cfl * m_ionSolver->computeAdvectionDiffusionDt();

      break;
    }
    case TransportAlgorithm::Heun: {
      dt = m_cfl * m_ionSolver->computeAdvectionDiffusionDt();

      break;
    }
    case TransportAlgorithm::ImExCTU: {
      dt = m_cfl * m_ionSolver->computeAdvectionDt();

      break;
    }
    default: {
      MayDay::Error("StreamerInceptionStepper::computDt -- logic bust");

      break;
    }
    }
    auto computeGradient = [this](const Real & a_timeStep){
			     return (m_voltageCurve(m_time+a_timeStep)-m_voltageCurve(m_time))/ a_timeStep;
			   };
    Real prevDt = m_minDt;
    Real currDt = 2*m_minDt;
    Real currGradient = computeGradient(currDt);
    const Real currVoltage = m_voltageCurve(m_time);
    Real deltaV = currGradient*currDt;
    while(std::abs(deltaV) <= std::abs(m_deltaVPercent*currVoltage) && (m_time+currDt) < m_endTime){
      prevDt = currDt;
      currDt += m_minDt;
      currGradient = computeGradient(currDt);
      deltaV = currGradient*currDt;
    }
    dt = prevDt;
    
    dt = std::min(dt, m_maxDt);
    dt = std::max(dt, m_minDt);
  }

  return dt;
}

template <typename P, typename F, typename C>
Real
StreamerInceptionStepper<P, F, C>::advance(const Real a_dt)
{
  CH_TIME("StreamerInceptionStepper::advance");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::advance" << endl;
  }

  if (m_mode == Mode::Transient) {

    const Real curTime    = m_time + a_dt;
    const Real curVoltage = m_voltageCurve(curTime);

    // Move the negative ions -- still want this for the time step.
    this->computeIonVelocity(curVoltage);
    if (m_ionTransport) {
      this->computeIonDiffusion(curVoltage);
      this->advanceIons(a_dt);
    }

    // Seed new particles and compute the inception integral.
    m_tracerParticleSolver->setVelocity(m_electricFieldMinu);

    this->seedIonizationParticles(curVoltage);
    this->computeInceptionIntegralTransient(curVoltage);

    // Compute the critical volume and the electron appearance rate.
    const Real Vcr  = this->computeCriticalVolumeTransient();
    const Real Rdot = this->computeRdot(curVoltage);

    m_criticalVolume.emplace_back(curTime, Vcr);
    m_Rdot.emplace_back(curTime, Rdot);

    // Compute the inception probability using the trapezoidal rule.
    if (m_Rdot.size() >= 2) {
      Real p = 0.0;

      for (size_t i = 0; i < m_Rdot.size() - 1; i++) {
        const Real dt = m_Rdot[i + 1].first - m_Rdot[i].first;

        p += 0.5 * dt * (m_Rdot[i + 1].second + m_Rdot[i].second);
      }

      m_inceptionProbability.emplace_back(curTime, 1.0 - exp(-p));
    }

    // Get the maximum K-value
    Real maxK = -std::numeric_limits<Real>::max();
    Real minK = +std::numeric_limits<Real>::max();

    DataOps::getMaxMin(maxK, minK, m_inceptionIntegral, 0);

    m_maxK.emplace_back(m_time + a_dt, maxK);

    this->writeReportTransient();
  }
  else {
    MayDay::Error("StreamerInceptionStepper::advance -- must have 'StreamerInceptionStepper.mode = transient'");
  }

  return a_dt;
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::advanceIons(const Real a_dt) noexcept
{
  CH_TIME("StreamerInceptionStepper::advanceIons");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::advanceIons" << endl;
  }

  EBAMRCellData& phi = m_ionSolver->getPhi();

  // Use an outflow BC for transport. Note that extrapolateAdvectiveFluxToEB computes
  // -n.(phi * v). Since we compute phi^(k+1) = phi^k - dt*sum(fluxes) this means that a negative
  // flux is an incoming flux.
  EBAMRIVData& ebFlux = m_ionSolver->getEbFlux();
  m_ionSolver->extrapolateAdvectiveFluxToEB(ebFlux);
  DataOps::floor(ebFlux, 0.0);

  switch (m_transportAlgorithm) {
  case TransportAlgorithm::Euler: {
    EBAMRCellData divJ;
    m_amr->allocate(divJ, m_realm, m_phase, 1.0);

    m_ionSolver->computeDivJ(divJ, phi, 0.0, false, true, false);

    DataOps::incr(phi, divJ, -a_dt);

    break;
  }
  case TransportAlgorithm::Heun: {
    // Transient storage
    EBAMRCellData yp;
    EBAMRCellData k1;
    EBAMRCellData k2;

    m_amr->allocate(yp, m_realm, m_phase, 1);
    m_amr->allocate(k1, m_realm, m_phase, 1);
    m_amr->allocate(k2, m_realm, m_phase, 1);

    // Compute k1 coefficient
    m_ionSolver->computeDivJ(k1, phi, 0.0, false, true, false);
    DataOps::copy(yp, phi);
    DataOps::incr(yp, k1, -a_dt);

    // Compute k2 coefficient and final state
    m_ionSolver->computeDivJ(k2, yp, 0.0, false, true, false);
    DataOps::incr(phi, k1, -0.5 * a_dt);
    DataOps::incr(phi, k2, -0.5 * a_dt);

    break;
  }
  case TransportAlgorithm::ImExCTU: {
    const bool addEbFlux     = true;
    const bool addDomainFlux = true;

    // Transient storage
    EBAMRCellData k1;
    EBAMRCellData k2;

    m_amr->allocate(k1, m_realm, m_phase, 1);
    m_amr->allocate(k2, m_realm, m_phase, 1);

    m_ionSolver->computeDivF(k1, phi, a_dt, false, true, true);

    DataOps::kappaScale(k1);
    DataOps::scale(k1, -1.0);

    // Use k1 as the old solution.
    DataOps::copy(k2, phi);

    // Do the Euler solve.
    m_ionSolver->advanceCrankNicholson(phi, k2, k1, a_dt);

    break;
  }
  default: {
    MayDay::Error("StreamerInceptionStepper::advanceIons -- logic bust");
  }
  }

  DataOps::floor(phi, 0.0);

  m_amr->average(phi, m_realm, m_phase, Average::Conservative);
  m_amr->interpGhost(phi, m_realm, m_phase);
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt)
{
  CH_TIME("StreamerInceptionStepper::synchronizeSolverTimes");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::synchronizeSolverTimes" << endl;
  }

  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;

  m_fieldSolver->setTime(a_step, a_time, a_dt);
  m_tracerParticleSolver->setTime(a_step, a_time, a_dt);
  m_ionSolver->setTime(a_step, a_time, a_dt);
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::printStepReport()
{
  CH_TIME("StreamerInceptionStepper::printStepReport");
#ifndef NDEBUG
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::printStepReport" << endl;
  }
#endif

  pout() << "                                ** Voltage                = " << m_voltageCurve(m_time) << endl;
  pout() << "                                ** Crit. volume           = " << m_criticalVolume.back().second << endl;
  pout() << "                                ** Inception probability  = " << m_inceptionProbability.back().second
         << endl;
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::preRegrid(const int a_lmin, const int a_oldFinestLevel)
{
  CH_TIME("StreamerInceptionStepper::preRegrid");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::preRegrid" << endl;
  }

  m_fieldSolver->preRegrid(a_lmin, a_oldFinestLevel);
  m_tracerParticleSolver->preRegrid(a_lmin, a_oldFinestLevel);
  m_ionSolver->preRegrid(a_lmin, a_oldFinestLevel);
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("StreamerInceptionStepper::regrid");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::regrid" << endl;
  }

  // Regrid tracer particles and field
  m_fieldSolver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_tracerParticleSolver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_ionSolver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);

  // Regrid velocity field
  m_amr->reallocate(m_electricFieldPlus, m_phase, a_lmin);
  m_amr->reallocate(m_electricFieldMinu, m_phase, a_lmin);

  switch (m_mode) {
  case Mode::Stationary: {
    m_amr->reallocate(m_inceptionIntegralPlus, m_phase, a_lmin);
    m_amr->reallocate(m_inceptionIntegralMinu, m_phase, a_lmin);
    m_amr->reallocate(m_inceptionVoltagePlus, m_phase, a_lmin);
    m_amr->reallocate(m_inceptionVoltageMinu, m_phase, a_lmin);

    if (m_plotBackgroundIonization) {
      m_amr->reallocate(m_backgroundIonization, m_phase, a_lmin);
    }

    if (m_plotDetachment) {
      m_amr->reallocate(m_detachment, m_phase, a_lmin);
    }

    if (m_plotFieldEmission) {
      m_amr->reallocate(m_emissionRatesPlus, m_phase, a_lmin);
      m_amr->reallocate(m_emissionRatesMinu, m_phase, a_lmin);
    }

    break;
  }
  case Mode::Transient: {

    m_amr->reallocate(m_inceptionIntegral, m_phase, a_lmin);
    m_amr->reallocate(m_emissionRate, m_phase, a_lmin);
    m_amr->reallocate(m_backgroundIonization, m_phase, a_lmin);

    break;
  }
  default: {
    break;
  }
  }

  // Solve Poisson equation again.
  m_fieldSolver->setRho(m_rho);
  m_fieldSolver->setSigma(0.0);

  const bool converged = m_fieldSolver->solve(m_fieldSolver->getPotential(),
                                              m_fieldSolver->getRho(),
                                              m_fieldSolver->getSigma());

  if (!converged) {
    MayDay::Error("StreamerInceptionStepper::initialData -- could not solve Poisson equation. Aborting application.");
  }

  // Get the gas-side electric field
  DataOps::copy(m_electricFieldPlus, m_amr->alias(phase::gas, m_fieldSolver->getElectricField()));

  m_amr->average(m_electricFieldPlus, m_realm, m_phase, Average::Conservative);
  m_amr->interpGhostMG(m_electricFieldPlus, m_realm, m_phase);

  m_amr->interpToCentroids(m_electricFieldPlus, m_realm, m_phase);

  // Also store -E
  DataOps::copy(m_electricFieldMinu, m_electricFieldPlus);
  DataOps::scale(m_electricFieldMinu, -1.0);

  // Re-compute the ion velocity and diffusion coefficients
  this->computeIonVelocity(m_voltageCurve(m_time));
  this->computeIonDiffusion(m_voltageCurve(m_time));
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::postRegrid()
{
  CH_TIME("StreamerInceptionStepper::postRegrid");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::postRegrid" << endl;
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::setVoltageCurve(const std::function<Real(const Real& E)>& a_voltageCurve) noexcept
{
  CH_TIME("StreamerInceptionStepper::setVoltageCurve");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setVoltageCurve" << endl;
  }

  m_voltageCurve = a_voltageCurve;
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::setRho(const std::function<Real(const RealVect& x)>& a_rho) noexcept
{
  CH_TIME("StreamerInceptionStepper::setRho");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setRho" << endl;
  }

  m_rho = a_rho;
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::setIonDensity(const std::function<Real(const RealVect x)>& a_density) noexcept
{
  CH_TIME("StreamerInceptionStepper::setIonDensity");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setIonDensity" << endl;
  }

  m_initialIonDensity = a_density;
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::setIonMobility(const std::function<Real(const Real x)>& a_mobility) noexcept
{
  CH_TIME("StreamerInceptionStepper::setIonMobility");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setIonMobility" << endl;
  }

  m_ionMobility = a_mobility;
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::setIonDiffusion(const std::function<Real(const Real x)>& a_diffCo) noexcept
{
  CH_TIME("StreamerInceptionStepper::setIonDiffusion");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setIonDiffusion" << endl;
  }

  m_ionDiffusion = a_diffCo;
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::setAlpha(const std::function<Real(const Real& E)>& a_alpha) noexcept
{
  CH_TIME("StreamerInceptionStepper::setAlpha");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setAlpha" << endl;
  }

  m_alpha = a_alpha;
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::setEta(const std::function<Real(const Real& E)>& a_eta) noexcept
{
  CH_TIME("StreamerInceptionStepper::setEta");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setEta" << endl;
  }

  m_eta = a_eta;
}

template <typename P, typename F, typename C>
const std::function<Real(const Real& E)>&
StreamerInceptionStepper<P, F, C>::getAlpha() const noexcept
{
  return m_alpha;
}

template <typename P, typename F, typename C>
const std::function<Real(const Real& E)>&
StreamerInceptionStepper<P, F, C>::getEta() const noexcept
{
  return m_eta;
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::setBackgroundRate(
  const std::function<Real(const Real& E, const RealVect& x)>& a_backgroundRate) noexcept
{
  CH_TIME("StreamerInceptionStepper::setBackgroundRate");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setBackgroundRate" << endl;
  }

  m_backgroundRate = a_backgroundRate;
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::setDetachmentRate(
  const std::function<Real(const Real& E, const RealVect& x)>& a_detachmentRate) noexcept
{
  CH_TIME("StreamerInceptionStepper::setDetachmentRate");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setDetachmentRate" << endl;
  }

  m_detachmentRate = a_detachmentRate;
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::setFieldEmission(
  const std::function<Real(const Real& E, const RealVect& x)>& a_currentDensity) noexcept
{
  CH_TIME("StreamerInceptionStepper::setFieldEmission");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setFieldEmission" << endl;
  }

  m_fieldEmission = a_currentDensity;
}

template <typename P, typename F, typename C>
Mode
StreamerInceptionStepper<P, F, C>::getMode() const noexcept
{
  return m_mode;
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::seedUniformParticles() noexcept
{
  CH_TIME("StreamerInceptionStepper::seedUniformParticles");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::seedUniformParticles" << endl;
  }

  ParticleContainer<P>& amrParticles = m_tracerParticleSolver->getParticles();
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
          P p;
          p.position() = m_amr->getProbLo() + (0.5 * RealVect::Unit + RealVect(iv)) * dx;

          particles.add(p);
        }
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        if (validCells(vof.gridIndex())) {
          P p;
          p.position() = m_amr->getProbLo() + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);

          particles.add(p);
        }
      };

      // Execute kernels over appropriate regions.
      const Box    cellBox = dbl[dit()];
      VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        // Set mass to zero and do a backup of the initial position (we need to rewind particles
        // back later on). The scalar<0> flag is the flag we use for start-stop criteria for the
        // particle integration.
        p.weight()           = 0.0;
        p.template vect<0>() = p.position();
      }
    }
  }

  // Remove particles inside the EB.
  m_amr->removeCoveredParticlesIF(amrParticles, m_phase);

  amrParticles.remap();
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::seedIonizationParticles(const Real& a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::seedIonizationParticles");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::seedIonizationParticles" << endl;
  }

  ParticleContainer<P>& amrParticles = m_tracerParticleSolver->getParticles();
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

      const EBCellFAB& electricField    = (*m_electricFieldPlus[lvl])[dit()];
      const FArrayBox& electricFieldReg = electricField.getFArrayBox();

      if (!ebisbox.isAllCovered()) {

        auto regularKernel = [&](const IntVect& iv) -> void {
          if (validCells(iv, 0) && ebisbox.isRegular(iv)) {

            const RealVect EE = RealVect(
              D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));
            const Real E = std::abs(a_voltage * EE.vectorLength());

            const Real alpha = m_alpha(E);
            const Real eta   = m_eta(E);

            if (alpha >= eta) {
              P p;

              p.position() = m_amr->getProbLo() + (0.5 * RealVect::Unit + RealVect(iv)) * dx;

              particles.add(p);
            }
          }
        };

        auto irregularKernel = [&](const VolIndex& vof) -> void {
          const IntVect iv = vof.gridIndex();

          if (validCells(iv, 0) && ebisbox.isIrregular(iv)) {

            const RealVect EE = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));
            const Real     E  = std::abs(a_voltage * EE.vectorLength());

            const Real alpha = m_alpha(E);
            const Real eta   = m_eta(E);

            if (alpha > eta) {
              P p;

              p.position() = m_amr->getProbLo() + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);

              particles.add(p);
            }
          }
        };

        // Execute kernels over appropriate regions.
        const Box    cellBox = dbl[dit()];
        VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

        BoxLoops::loop(cellBox, regularKernel);
        BoxLoops::loop(vofit, irregularKernel);
      }

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        p.weight()           = 0.0;
        p.template vect<0>() = p.position();
      }
    }
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::postInitialize()
{
  CH_TIME("StreamerInceptionStepper::postInitialize");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::postInitialize" << endl;
  }

  // Initialize the negative ion solver
  m_ionSolver->initialData();
  this->computeIonVelocity(m_voltageCurve(m_time));
  this->computeIonDiffusion(m_voltageCurve(m_time));

  switch (m_mode) {
  case Mode::Stationary: {

    this->computeInceptionIntegralStationary();
    this->computeInceptionVoltageVolume();
    this->computeCriticalVolumeStationary();
    this->computeRdotStationary();

    // Compute background ionization and field emission on the mesh. Mostly for plotting.
    if (m_plotBackgroundIonization) {
      this->computeBackgroundIonizationStationary();
    }

    if (m_plotDetachment) {
      this->computeDetachmentStationary();
    }

    if (m_plotFieldEmission) {
      this->computeFieldEmissionStationary();
    }

    this->writeReportStationary();

    break;
  }
  case Mode::Transient: {

    // Seed tracer particles and set velocity
    this->seedIonizationParticles(m_voltageCurve(m_time));
    m_tracerParticleSolver->setVelocity(m_electricFieldMinu);

    // Compute the 'Rdot' value at current time (which is probably t = 0).
    this->computeInceptionIntegralTransient(m_voltageCurve(m_time));

    // Get the maximum K-value
    Real maxK = -std::numeric_limits<Real>::max();
    Real minK = +std::numeric_limits<Real>::max();

    DataOps::getMaxMin(maxK, minK, m_inceptionIntegral, 0);

    m_Rdot.emplace_back(m_time, this->computeRdot(m_voltageCurve(m_time)));
    m_criticalVolume.emplace_back(m_time, this->computeCriticalVolumeTransient());
    m_inceptionProbability.emplace_back(m_time, 0.0);
    m_maxK.emplace_back(m_time, maxK);

    break;
  }
  default: {
    break;
  }
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::computeInceptionIntegralStationary() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeInceptionIntegralStationary");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeInceptionIntegralStationary" << endl;
  }

  // TLDR: For this integration we compute the value of the inception integral
  //
  //          K = int alpha_eff(E)dl
  //
  //       The integration runs from the start position of the particle and along a field line until either the
  //       particle leaves the domain or alpha < 0.0.
  //
  //       For the Euler rule we get
  //
  //          T += alpha_eff(E(x)) * dx
  //
  //       For the trapezoidal rule we get
  //
  //          T += 0.5 * dx * [alpha_eff(E(x)) + alpha_eff(E(x+dx))]

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
    std::vector<Real> polarities{1.0, -1.0};

    this->seedIonizationParticles(m_voltageSweeps[i]);

    for (const auto& p : polarities) {

      if (p > 0.0) {
        m_tracerParticleSolver->setVelocity(m_electricFieldMinu);
      }
      else {
        m_tracerParticleSolver->setVelocity(m_electricFieldPlus);
      }

      this->resetTracerParticles();

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
        MayDay::Error("StreamerInceptionStepper::computeInceptionIntegralStationary -- logic bust");

        break;
      }
      }

      // Deposit particles on mesh and copy to relevant data holder.
      if (p > 0.0) {
        m_tracerParticleSolver->deposit(Kplus);

        m_amr->conservativeAverage(Kplus, m_realm, m_phase);
        m_amr->interpGhost(Kplus, m_realm, m_phase);

        DataOps::copy(m_inceptionIntegralPlus, Kplus, Interval(i, i), Interval(0, 0));

        // Get max K-value
        Real maxK = -std::numeric_limits<Real>::max();
        Real minK = +std::numeric_limits<Real>::max();

        DataOps::getMaxMin(maxK, minK, Kplus, 0);

        m_maxKPlus.push_back(maxK);
      }
      else {
        m_tracerParticleSolver->deposit(Kminu);

        m_amr->conservativeAverage(Kminu, m_realm, m_phase);
        m_amr->interpGhost(Kminu, m_realm, m_phase);

        DataOps::copy(m_inceptionIntegralMinu, Kminu, Interval(i, i), Interval(0, 0));

        // Get max K-value
        Real maxK = -std::numeric_limits<Real>::max();
        Real minK = +std::numeric_limits<Real>::max();

        DataOps::getMaxMin(maxK, minK, Kminu, 0);

        m_maxKMinu.push_back(maxK);
      }
    }
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::computeInceptionIntegralTransient(const Real& a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::computeInceptionIntegralTransient");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeInceptionIntegralTransient" << endl;
  }

  // Compute K-value on each particle using specified algorithm.
  switch (m_inceptionAlgorithm) {
  case IntegrationAlgorithm::Euler: {
    this->inceptionIntegrateEuler(a_voltage);

    break;
  }
  case IntegrationAlgorithm::Trapezoidal: {
    this->inceptionIntegrateTrapezoidal(a_voltage);

    break;
  }
  default: {
    MayDay::Error("StreamerInceptionStepper::computeInceptionIntegralTransient - logic bust");

    break;
  }
  }

  // Deposit the particles onto m_inceptionIntegral
  m_tracerParticleSolver->deposit(m_inceptionIntegral);

  m_amr->conservativeAverage(m_inceptionIntegral, m_realm, m_phase);
  m_amr->interpGhost(m_inceptionIntegral, m_realm, m_phase);
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::inceptionIntegrateEuler(const Real& a_voltage) noexcept
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
  m_amr->allocate(amrProcessedParticles, m_realm);

  ParticleContainer<P>& amrParticles = m_tracerParticleSolver->getParticles();

  m_tracerParticleSolver->remap();
  m_tracerParticleSolver->interpolateVelocities();

  const size_t outcastsBefore  = amrParticles.getNumberOfOutcastParticlesGlobal();
  const size_t particlesBefore = amrParticles.getNumberOfValidParticlesGlobal();

  const Real sgn = (a_voltage > 0.0) - (a_voltage < 0.0);

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

          const RealVect vel      = sgn * p.velocity();
          const Real     v        = vel.vectorLength();
          const Real     E        = std::abs(a_voltage * v);
          const Real     dt       = step / v;
          const Real     alpha    = m_alpha(E);
          const Real     eta      = m_eta(E);
          const Real     alphaEff = alpha - eta;

          // Update particle position
          const RealVect newPos = p.position() + dt * vel;

          // Stop integration for particles that move into regions alpha < 0.0,
          // inside the EB or outside of the domain.
          const bool negativeAlpha = alphaEff < 0.0;
          const bool outsideDomain = this->particleOutsideGrid(newPos, probLo, probHi);
          const bool insideEB      = this->particleInsideEB(newPos);

          if (negativeAlpha || outsideDomain || insideEB) {
            processedParticles.transfer(lit);
          }
          else {
            p.position() = newPos;
            p.weight() += step * alphaEff;

            ++lit;
          }
        }
      }
    }

    // Update velocities.
    m_tracerParticleSolver->remap();
    m_tracerParticleSolver->interpolateVelocities();
  }

  ParticleOps::copyDestructive(amrParticles, amrProcessedParticles);

  this->rewindTracerParticles();

  const size_t outcastsAfter  = amrParticles.getNumberOfOutcastParticlesGlobal();
  const size_t particlesAfter = amrParticles.getNumberOfValidParticlesGlobal();

  if (particlesBefore != particlesAfter) {
    MayDay::Error("StreamerInceptionStepper::inceptionIntegrateEuler - lost/gained particles!");
  }
  if (outcastsBefore != outcastsAfter) {
    MayDay::Error("StreamerInceptionStepper::inceptionIntegrateEuler - lost/gained outcast particles!");
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::inceptionIntegrateTrapezoidal(const Real& a_voltage) noexcept
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
  //          T += 0.5 * dx * [alpha_eff(E(p^k)) + alpha_eff(E(p^(k+1)))],
  //
  //       where dx = |p^(k+1) - p^k|.

  const RealVect probLo = m_amr->getProbLo();
  const RealVect probHi = m_amr->getProbHi();

  // Allocate a data holder for holding the processed particles. This
  // will be faster because then we only have to iterate through the
  // particles that are actually still moving.
  ParticleContainer<P> amrProcessedParticles;
  m_amr->allocate(amrProcessedParticles, m_realm);

  ParticleContainer<P>& amrParticles = m_tracerParticleSolver->getParticles();

  m_tracerParticleSolver->remap();
  m_tracerParticleSolver->interpolateVelocities();

  const Real sgn = (a_voltage > 0.0) - (a_voltage < 0.0);

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

          const RealVect v        = p.velocity();
          const Real     dt       = step / v.vectorLength();
          const Real     E        = std::abs(a_voltage * v.vectorLength());
          const Real     alpha    = m_alpha(E);
          const Real     eta      = m_eta(E);
          const Real     alphaEff = alpha - eta;

          // Storage alpha(p^k), v(p^k), and the step size.
          p.template real<0>() = alphaEff;
          p.template real<1>() = dt;
          p.template vect<1>() = v;

          // Update position.
          const RealVect newPos = p.position() + sgn * v * dt;

          // Stop integration for particles that move into regions alpha < 0.0,
          // inside the EB or outside of the domain.
          const bool negativeAlpha = (alphaEff <= 0.0);
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
    m_tracerParticleSolver->remap();
    m_tracerParticleSolver->interpolateVelocities();

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
          const Real     E   = std::abs(a_voltage * vk1.vectorLength());

          // Note the weird subtraction since p.position() was updated
          // to p^k + dt * v^k.
          const RealVect newPos = p.position() + 0.5 * sgn * dt * (vk1 - vk);

          // Compute new alpha.
          const Real alpha1 = p.template real<0>();
          const Real alpha  = m_alpha(E);
          const Real eta    = m_eta(E);
          const Real alpha2 = alpha - eta;
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
    m_tracerParticleSolver->remap();
    m_tracerParticleSolver->interpolateVelocities();
  }

  // Copy processed particles over to the solver particles.
  ParticleOps::copyDestructive(amrParticles, amrProcessedParticles);

  this->rewindTracerParticles();
}

template <typename P, typename F, typename C>
Real
StreamerInceptionStepper<P, F, C>::computeRdot(const Real& a_voltage) const noexcept
{
  CH_TIME("StreamerInceptionStepper::computeRdot");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeRdot" << endl;
  }

  Real Rdot = 0.0;

  // TLDR: We are computing the integral of dne/dt * (1-eta/alpha) over the critical volume and integral(j/e * dS) over the "critical surface"

  const RealVect probLo = m_amr->getProbLo();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real               dx    = m_amr->getDx()[lvl];

    const Real vol  = std::pow(dx, SpaceDim);
    const Real area = std::pow(dx, SpaceDim - 1);

    for (DataIterator dit(dbl); dit.ok(); ++dit) {

      const EBISBox&       ebisbox    = ebisl[dit()];
      const BaseFab<bool>& validCells = (*m_amr->getValidCells(m_realm)[lvl])[dit()];

      const EBCellFAB& electricField     = (*m_electricFieldPlus[lvl])[dit()];
      const EBCellFAB& inceptionIntegral = (*m_inceptionIntegral[lvl])[dit()];
      const EBCellFAB& fieldEmission     = (*m_emissionRate[lvl])[dit()];
      const EBCellFAB& ionDensity        = (*m_ionSolver->getPhi()[lvl])[dit()];

      const FArrayBox& electricFieldReg     = electricField.getFArrayBox();
      const FArrayBox& inceptionIntegralReg = inceptionIntegral.getFArrayBox();
      const FArrayBox& ionDensityReg        = ionDensity.getFArrayBox();

      // Add contribution from background ionization in regular cells.
      auto regularKernel = [&](const IntVect& iv) -> void {
        if (ebisbox.isRegular(iv) && validCells(iv, 0)) {
          if (inceptionIntegralReg(iv, 0) >= m_inceptionK) {
            const RealVect pos = probLo + (0.5 * RealVect::Unit + RealVect(iv)) * dx;
            const RealVect EE  = RealVect(
              D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));
            const Real E = std::abs(a_voltage * EE.vectorLength());

            const Real alpha = m_alpha(E);
            const Real eta   = m_eta(E);
            const Real k     = m_detachmentRate(E, pos);
            const Real dndt  = m_backgroundRate(E, pos) + k * ionDensityReg(iv, 0);

            CH_assert(alpha >= eta);

            Rdot += dndt * (1.0 - eta / alpha) * vol;
          }
        }
      };

      // Add contribution from background ionization and field emission in cut-cells.
      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const IntVect iv = vof.gridIndex();
        if (ebisbox.isIrregular(iv) && validCells(iv, 0)) {
          if (inceptionIntegral(vof, 0) >= m_inceptionK) {
            const RealVect pos = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);
            const RealVect EE  = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));
            const Real     E   = std::abs(a_voltage * EE.vectorLength());

            const Real kappa    = ebisbox.volFrac(vof);
            const Real areaFrac = ebisbox.bndryArea(vof);
            const Real alpha    = m_alpha(E);
            const Real eta      = m_eta(E);
            const Real k        = m_detachmentRate(E, pos);
            const Real j        = m_fieldEmission(E, pos);
            const Real dndt     = m_backgroundRate(E, pos) + k * ionDensity(vof, 0);

            Rdot += dndt * (1.0 - eta / alpha) * kappa * vol;
            Rdot += j / Units::Qe * areaFrac * area;
          }
        }
      };

      // Run the kernels.
      Box          cellBox = dbl[dit()];
      VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }

  return ParallelOps::sum(Rdot);
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::rewindTracerParticles() noexcept
{
  CH_TIME("StreamerInceptionStepper::rewindTracerParticles");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::rewindTracerParticles" << endl;
  }

  ParticleContainer<P>& amrParticles = m_tracerParticleSolver->getParticles();

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

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::resetTracerParticles() noexcept
{
  CH_TIME("StreamerInceptionStepper::resetTracerParticles");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::resetTracerParticles" << endl;
  }

  ParticleContainer<P>& amrParticles = m_tracerParticleSolver->getParticles();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
      for (ListIterator<P> lit(amrParticles[lvl][dit()].listItems()); lit.ok(); ++lit) {
        P& p = lit();

        p.weight() = 0.0;
      }
    }
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::computeBackgroundIonizationStationary() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeBackgroundIonizationStationary");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeBackgroundIonizationStationary" << endl;
  }

  const RealVect probLo = m_amr->getProbLo();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real               dx    = m_amr->getDx()[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const EBISBox& ebisbox = ebisl[dit()];

      EBCellFAB& bgIonization    = (*m_backgroundIonization[lvl])[dit()];
      FArrayBox& bgIonizationReg = bgIonization.getFArrayBox();

      const EBCellFAB& electricField    = (*m_electricFieldPlus[lvl])[dit()];
      const FArrayBox& electricFieldReg = electricField.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        const RealVect pos = probLo + (0.5 * RealVect::Unit + RealVect(iv)) * dx;
        const RealVect EE = RealVect(D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));
        const Real     E  = EE.vectorLength();

        for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
          bgIonizationReg(iv, i) = m_backgroundRate(std::abs(m_voltageSweeps[i] * E), pos);
        }
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const RealVect pos = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);
        const RealVect EE  = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));
        const Real     E   = EE.vectorLength();

        for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
          bgIonization(vof, i) = m_backgroundRate(std::abs(m_voltageSweeps[i] * E), pos);
        }
      };

      // Execute kernels
      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      BoxLoops::loop(dbl[dit()], regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::computeDetachmentStationary() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeDetachmentStationary");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeDetachmentStationary" << endl;
  }

  const RealVect probLo = m_amr->getProbLo();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real               dx    = m_amr->getDx()[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const EBISBox& ebisbox = ebisl[dit()];

      EBCellFAB& detachment    = (*m_detachment[lvl])[dit()];
      FArrayBox& detachmentReg = detachment.getFArrayBox();

      const EBCellFAB& ionDensity    = (*m_ionSolver->getPhi()[lvl])[dit()];
      const EBCellFAB& electricField = (*m_electricFieldPlus[lvl])[dit()];

      const FArrayBox& ionDensityReg    = ionDensity.getFArrayBox();
      const FArrayBox& electricFieldReg = electricField.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        const RealVect pos = probLo + (0.5 * RealVect::Unit + RealVect(iv)) * dx;
        const RealVect EE = RealVect(D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));
        const Real     E  = EE.vectorLength();
        const Real     phi = ionDensityReg(iv, 0);

        for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
          detachmentReg(iv, i) = m_detachmentRate(std::abs(m_voltageSweeps[i] * E), pos) * phi;
        }
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const RealVect pos = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);
        const RealVect EE  = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));
        const Real     E   = EE.vectorLength();
        const Real     phi = ionDensity(vof, 0);

        for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
          detachment(vof, i) = m_detachmentRate(std::abs(m_voltageSweeps[i] * E), pos) * phi;
        }
      };

      // Execute kernels
      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      BoxLoops::loop(dbl[dit()], regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }

  m_amr->conservativeAverage(m_detachment, m_realm, m_phase);
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::computeFieldEmissionStationary() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeFieldEmissionStationary");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeFieldEmissionStationary" << endl;
  }

  const RealVect probLo = m_amr->getProbLo();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real&              dx    = m_amr->getDx()[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const EBISBox& ebisbox = ebisl[dit()];

      // Here, emissionPlus is when we have the standard voltage (V = 1) and
      // emissionMinu is when the voltage is reverted.
      EBCellFAB& emissionPlus = (*m_emissionRatesPlus[lvl])[dit()];
      EBCellFAB& emissionMinu = (*m_emissionRatesMinu[lvl])[dit()];

      emissionPlus.setVal(0.0);
      emissionMinu.setVal(0.0);

      const EBCellFAB& electricField    = (*m_electricFieldPlus[lvl])[dit()];
      const FArrayBox& electricFieldReg = electricField.getFArrayBox();

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const RealVect pos = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);
        const RealVect E   = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));

        const Real normalE    = E.dotProduct(ebisbox.normal(vof));
        const Real magnitudeE = E.vectorLength();

        if (normalE > 0.0) {
          // Field points out of the electrode but for the opposite polarity it will
          // point into the electrode and the electrons can then escape.
          for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
            emissionPlus(vof, i) = 0.0;
            emissionMinu(vof, i) = m_fieldEmission(std::abs(magnitudeE * m_voltageSweeps[i]), pos);
          }
        }
        else if (normalE < 0.0) {
          // Opposite case of the above.
          for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
            emissionPlus(vof, i) = m_fieldEmission(std::abs(magnitudeE * m_voltageSweeps[i]), pos);
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

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::computeFieldEmission(EBAMRCellData& a_emissionRate,
                                                        const Real&    a_voltage) const noexcept
{
  CH_TIME("StreamerInceptionStepper::computeFieldEmission");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeFieldEmission" << endl;
  }

  CH_assert(a_emissionRate[0]->nComp() == 1);

  const RealVect probLo = m_amr->getProbLo();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real               dx    = m_amr->getDx()[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const EBISBox& ebisbox = ebisl[dit()];

      // Here, emissionPlus is when we have the standard voltage (V = 1) and
      // emissionMinu is when the voltage is reverted.
      EBCellFAB& emission = (*a_emissionRate[lvl])[dit()];

      emission.setVal(0.0);

      const EBCellFAB& electricField    = (*m_electricFieldPlus[lvl])[dit()];
      const FArrayBox& electricFieldReg = electricField.getFArrayBox();

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const RealVect pos = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);
        const RealVect EE  = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));
        const Real     E   = std::abs(a_voltage * EE.vectorLength());

        if (EE.dotProduct(ebisbox.normal(vof)) > 0.0) {
          emission(vof, 0) = m_fieldEmission(E, pos);
        }
      };

      // Execute kernels
      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
      BoxLoops::loop(vofit, irregularKernel);
    }
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::evaluateFunction(EBAMRCellData&                           a_data,
                                                    const Real&                              a_voltage,
                                                    const std::function<Real(const Real E)>& a_func) const noexcept
{
  CH_TIME("StreamerInceptionStepper::evaluateFunction");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::evaluateFunction" << endl;
  }

  CH_assert(a_data[0]->nComp() == 1);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const EBISBox& ebisbox = ebisl[dit()];

      // Here, emissionPlus is when we have the standard voltage (V = 1) and
      // emissionMinu is when the voltage is reverted.
      EBCellFAB& data    = (*a_data[lvl])[dit()];
      FArrayBox& dataReg = data.getFArrayBox();

      const EBCellFAB& electricField    = (*m_electricFieldPlus[lvl])[dit()];
      const FArrayBox& electricFieldReg = electricField.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        const RealVect EE = RealVect(D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));
        const Real     E  = std::abs(a_voltage * EE.vectorLength());

        dataReg(iv, 0) = a_func(E);
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const RealVect EE = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));
        const Real     E  = std::abs(a_voltage * EE.vectorLength());

        data(vof, 0) = a_func(E);
      };

      // Execute kernels
      Box          cellBox = dbl[dit()];
      VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::evaluateFunction(
  EBAMRCellData&                                             a_data,
  const Real&                                                a_voltage,
  const std::function<Real(const Real E, const RealVect x)>& a_func) const noexcept
{
  CH_TIME("StreamerInceptionStepper::evaluateFunction");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::evaluateFunction" << endl;
  }

  CH_assert(a_data[0]->nComp() == 1);

  const RealVect probLo = m_amr->getProbLo();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real               dx    = m_amr->getDx()[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const EBISBox& ebisbox = ebisl[dit()];

      // Here, emissionPlus is when we have the standard voltage (V = 1) and
      // emissionMinu is when the voltage is reverted.
      EBCellFAB& data    = (*a_data[lvl])[dit()];
      FArrayBox& dataReg = data.getFArrayBox();

      data.setVal(0.0);

      const EBCellFAB& electricField    = (*m_electricFieldPlus[lvl])[dit()];
      const FArrayBox& electricFieldReg = electricField.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        const RealVect pos = probLo + (0.5 * RealVect::Unit + RealVect(iv)) * dx;
        const RealVect EE = RealVect(D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));
        const Real     E  = std::abs(a_voltage * EE.vectorLength());

        dataReg(iv, 0) = a_func(E, pos);
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const RealVect pos = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);
        const RealVect EE  = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));
        const Real     E   = std::abs(a_voltage * EE.vectorLength());

        data(vof, 0) = a_func(E, pos);
      };

      // Execute kernels
      Box          cellBox = dbl[dit()];
      VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::computeInceptionVoltageVolume() noexcept
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
        EBCellFAB& inceptionVoltagePlus = (*m_inceptionVoltagePlus[lvl])[dit()];
        EBCellFAB& inceptionVoltageMinu = (*m_inceptionVoltageMinu[lvl])[dit()];

        FArrayBox& inceptionVoltagePlusReg = inceptionVoltagePlus.getFArrayBox();
        FArrayBox& inceptionVoltageMinuReg = inceptionVoltageMinu.getFArrayBox();

        const EBCellFAB& inceptionIntegralPlus = (*m_inceptionIntegralPlus[lvl])[dit()];
        const EBCellFAB& inceptionIntegralMinu = (*m_inceptionIntegralMinu[lvl])[dit()];

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
    m_amr->conservativeAverage(m_inceptionVoltagePlus, m_realm, m_phase);
    m_amr->interpGhost(m_inceptionVoltageMinu, m_realm, m_phase);
  }
}

template <typename P, typename F, typename C>
std::pair<Real, RealVect>
StreamerInceptionStepper<P, F, C>::computeMinimumInceptionVoltage(const EBAMRCellData& a_Uinc) const noexcept
{
  CH_TIME("StreamerInceptionStepper::computeMinimumInceptionVoltage");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeMinimumInceptionVoltage" << endl;
  }

  Real     UInc = std::numeric_limits<Real>::infinity();
  RealVect xInc = RealVect::Zero;

  const RealVect probLo = m_amr->getProbLo();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real&              dx    = m_amr->getDx()[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const BaseFab<bool>& validCells = (*m_amr->getValidCells(m_realm)[lvl])[dit()];
      const EBISBox&       ebisBox    = ebisl[dit()];

      const EBCellFAB& voltage    = (*a_Uinc[lvl])[dit()];
      const FArrayBox& voltageReg = voltage.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        if (validCells(iv, 0) && ebisBox.isRegular(iv)) {
          const Real& U = voltageReg(iv, 0);

          const RealVect pos = probLo + (0.5 * RealVect::Unit + RealVect(iv)) * dx;

          if (!(std::isnan(U))) {
            UInc = std::min(UInc, std::abs(U));
            xInc = pos;
          }
        }
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const IntVect iv = vof.gridIndex();
        if (validCells(iv, 0) && ebisBox.isIrregular(iv)) {
          const Real& U = voltage(vof, 0);

          const RealVect centroid = ebisBox.centroid(vof);
          const RealVect pos      = probLo + (0.5 * RealVect::Unit + RealVect(iv) + centroid) * dx;

          if (!(std::isnan(U))) {
            UInc = std::min(UInc, std::abs(U));
            xInc = pos;
          }
        }
      };

      Box          cellBox = dbl[dit()];
      VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }

  return ParallelOps::min(UInc, xInc);
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::computeCriticalVolumeStationary() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeCriticalVolumeStationary");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeCriticalVolumeStationary" << endl;
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
          if (validCells(iv, 0) && ebisbox.isRegular(iv)) {
            if (inceptionIntegralPlusReg(iv, i) >= m_inceptionK) {
              criticalVolumePlus += std::pow(dx, SpaceDim);
            }
            if (inceptionIntegralMinuReg(iv, i) >= m_inceptionK) {
              criticalVolumeMinu += std::pow(dx, SpaceDim);
            }
          }
        };

        auto irregularKernel = [&](const VolIndex& vof) -> void {
          if (validCells(vof.gridIndex())) {
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

template <typename P, typename F, typename C>
Real
StreamerInceptionStepper<P, F, C>::computeCriticalVolumeTransient() const noexcept
{
  CH_TIME("StreamerInceptionStepper::computeCriticalVolumeTransient");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeCriticalVolumeTransient" << endl;
  }

  Real Vcr = 0.0;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); ++lvl) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

    const LevelData<BaseFab<bool>>& validCellsLD = *m_amr->getValidCells(m_realm)[lvl];

    const Real dx  = m_amr->getDx()[lvl];
    const Real vol = std::pow(dx, SpaceDim);

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      const EBISBox&       ebisbox    = ebisl[dit()];
      const BaseFab<bool>& validCells = validCellsLD[dit()];

      const EBCellFAB& inceptionIntegral    = (*m_inceptionIntegral[lvl])[dit()];
      const FArrayBox& inceptionIntegralReg = inceptionIntegral.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        if (validCells(iv, 0) && ebisbox.isRegular(iv)) {
          if (inceptionIntegralReg(iv, 0) >= m_inceptionK) {
            Vcr += vol;
          }
        }
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        if (validCells(vof.gridIndex())) {
          if (inceptionIntegral(vof, 0) >= m_inceptionK) {
            Vcr += ebisbox.volFrac(vof) * vol;
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

  return ParallelOps::sum(Vcr);
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::computeRdotStationary() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeRdotStationary");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeRdotStationary" << endl;
  }

  // TLDR: This routine computes the average time to the appereance of an electron
  //       anywhere in the critical volume. We do this by assuming that we have a
  //       background ionization rate m_backgroundIonization(E).

  // Total rate.
  std::vector<Real> Rplus(m_voltageSweeps.size(), 0.0);
  std::vector<Real> Rminu(m_voltageSweeps.size(), 0.0);

  const RealVect probLo = m_amr->getProbLo();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real&              dx    = m_amr->getDx()[lvl];

    const Real vol  = std::pow(dx, SpaceDim);
    const Real area = std::pow(dx, SpaceDim - 1);

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const EBISBox&       ebisbox    = ebisl[dit()];
      const BaseFab<bool>& validCells = (*m_amr->getValidCells(m_realm)[lvl])[dit()];

      const EBCellFAB& electricField         = (*m_electricFieldPlus[lvl])[dit()];
      const EBCellFAB& inceptionIntegralPlus = (*m_inceptionIntegralPlus[lvl])[dit()];
      const EBCellFAB& inceptionIntegralMinu = (*m_inceptionIntegralMinu[lvl])[dit()];
      const EBCellFAB& ionDensity            = (*m_ionSolver->getPhi()[lvl])[dit()];

      const FArrayBox& electricFieldReg         = electricField.getFArrayBox();
      const FArrayBox& inceptionIntegralPlusReg = inceptionIntegralPlus.getFArrayBox();
      const FArrayBox& inceptionIntegralMinuReg = inceptionIntegralMinu.getFArrayBox();
      const FArrayBox& ionDensityReg            = ionDensity.getFArrayBox();

      // Regular kernel for background ionization plus detachment from negative ions.
      auto regularKernel = [&](const IntVect& iv) -> void {
        if (validCells(iv) && ebisbox.isRegular(iv)) {
          const RealVect pos = probLo + (0.5 * RealVect::Unit + RealVect(iv)) * dx;
          const RealVect EE  = RealVect(
            D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));

          for (size_t i = 0; i < m_voltageSweeps.size(); i++) {

            const Real E = std::abs(m_voltageSweeps[i] * EE.vectorLength());

            const Real alpha = m_alpha(E);
            const Real eta   = m_eta(E);
            const Real k     = m_detachmentRate(E, pos);
            const Real dndt  = m_backgroundRate(E, pos) + k * ionDensityReg(iv, 0);

            // Add background ionization rate.
            if (inceptionIntegralPlusReg(iv, i) >= m_inceptionK) {
              Rplus[i] += dndt * (1.0 - eta / alpha) * vol;
            }

            if (inceptionIntegralMinuReg(iv, i) >= m_inceptionK) {
              Rminu[i] += dndt * (1.0 - eta / alpha) * vol;
            }
          }
        }
      };

      // Irregular kernel for background ionization + detachment from negative ions.
      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const IntVect iv = vof.gridIndex();
        if (validCells(iv, 0)) {
          const Real     volFrac = ebisbox.volFrac(vof);
          const RealVect pos     = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);
          const RealVect EE = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));

          for (size_t i = 0; i < m_voltageSweeps.size(); i++) {

            const Real E = std::abs(m_voltageSweeps[i] * EE.vectorLength());

            const Real alpha = m_alpha(E);
            const Real eta   = m_eta(E);
            const Real k     = m_detachmentRate(E, pos);
            const Real dndt  = m_backgroundRate(E, pos) + k * ionDensity(vof, 0);

            // Add background ionization rate.
            if (inceptionIntegralPlus(vof, i) >= m_inceptionK) {
              Rplus[i] += dndt * (1.0 - eta / alpha) * vol * volFrac;
            }
            if (inceptionIntegralMinu(vof, i) >= m_inceptionK) {
              Rminu[i] += dndt * (1.0 - eta / alpha) * vol * volFrac;
            }
          }
        }
      };

      // Field emission irregular kernel.
      auto fieldEmissionKernel = [&](const VolIndex& vof) -> void {
        const IntVect iv = vof.gridIndex();

        if (validCells(iv, 0)) {
          const Real     areaFrac = ebisbox.bndryArea(vof);
          const RealVect pos      = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);
          const RealVect EE = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));

          for (size_t i = 0; i < m_voltageSweeps.size(); i++) {

            const Real E = std::abs(m_voltageSweeps[i] * EE.vectorLength());

            if (EE.dotProduct(ebisbox.normal(vof)) > 0.0) {
              if (inceptionIntegralMinu(vof, i) >= m_inceptionK) {
                Rminu[i] += m_fieldEmission(E, pos) * area * areaFrac;
              }
            }
            else if (EE.dotProduct(ebisbox.normal(vof)) < 0.0) {
              if (inceptionIntegralPlus(vof, i) >= m_inceptionK) {
                Rplus[i] += m_fieldEmission(E, pos) * area * areaFrac;
              }
            }
          }
        }
      };

      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      BoxLoops::loop(dbl[dit()], regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
      BoxLoops::loop(vofit, fieldEmissionKernel);
    }
  }

  for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
    m_RdotPlus.push_back(ParallelOps::sum(Rplus[i]));
    m_RdotMinu.push_back(ParallelOps::sum(Rminu[i]));
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::writeReportStationary() const noexcept
{
  CH_TIME("StreamerInceptionStepper::writeReportStationary");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::writeReportStationary" << endl;
  }

  const std::pair<Real, RealVect> UIncPlus = this->computeMinimumInceptionVoltage(m_inceptionVoltagePlus);
  const std::pair<Real, RealVect> UIncMinu = this->computeMinimumInceptionVoltage(m_inceptionVoltageMinu);

#ifdef CH_MPI
  if (procID() == 0) {
#endif

    std::ofstream output(m_outputFile, std::ofstream::out);

    // clang-format off
    output << "# ==================================================================================================================================" << "\n";
    if (std::isnan(UIncPlus.first) || std::isnan(UIncMinu.first)) {
      output << "# Could not compute inception voltage\n";
    }
    else {
      output << "# Inception voltage + = " << UIncPlus.first << ",\t x = " << UIncPlus.second << "\n";
      output << "# Inception voltage - = " << UIncMinu.first << ",\t x = " << UIncMinu.second << "\n";      
    }
    
    output << "# ==================================================================================================================================" << "\n";    
    output << left << setw(15) << setfill(' ') << "# +/- Voltage";
    output << left << setw(15) << setfill(' ') << "Max K(+)";
    output << left << setw(15) << setfill(' ') << "Max K(-)";
    output << left << setw(20) << setfill(' ') << "Crit. vol(+)";
    output << left << setw(20) << setfill(' ') << "Crit. vol(-)";
    output << left << setw(20) << setfill(' ') << "Rdot (+)";
    output << left << setw(20) << setfill(' ') << "Rdot (-)" << "\n";
    output << "# ==================================================================================================================================" << "\n";        

    for (int i = 0; i < m_voltageSweeps.size(); i++) {
      output << left << setw(15) << setfill(' ') << m_voltageSweeps[i];
      output << left << setw(15) << setfill(' ') << m_maxKPlus[i];
      output << left << setw(15) << setfill(' ') << m_maxKMinu[i];
      output << left << setw(20) << setfill(' ') << m_criticalVolumePlus[i];
      output << left << setw(20) << setfill(' ') << m_criticalVolumeMinu[i];
      output << left << setw(20) << setfill(' ') << m_RdotPlus[i];
      output << left << setw(20) << setfill(' ') << m_RdotMinu[i];
      output << endl;
    }
    output << "# ==================================================================================================================================" << "\n";
    // clang-format on
#ifdef CH_MPI
  }
#endif
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::writeReportTransient() const noexcept
{
  CH_TIME("StreamerInceptionStepper::writeReportTransient");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::writeReportTransient" << endl;
  }

#ifdef CH_MPI
  if (procID() == 0) {
#endif
    std::ofstream output(m_outputFile, std::ofstream::out);

    output << std::left << std::setw(15) << setfill(' ') << "# Time t";
    output << std::left << std::setw(15) << setfill(' ') << "V(t)";
    output << std::left << std::setw(15) << setfill(' ') << "max K(t)";
    output << std::left << std::setw(15) << setfill(' ') << "Vcr(t)";
    output << std::left << std::setw(15) << setfill(' ') << "Rdot(t)";
    output << std::left << std::setw(15) << setfill(' ') << "P(t)";
    output << std::left << std::setw(15) << setfill(' ') << "dP(t, t+dt)";
    output << "\n";

    // Compute dP(t, t+dt)
    std::vector<Real> dProb;
    for (size_t i = 0; i < m_Rdot.size() - 1; i++) {
      const Real dt   = m_Rdot[i + 1].first - m_Rdot[i].first;
      const Real prob = m_inceptionProbability[i].second;
      const Real Rdot = m_Rdot[i].second;

      dProb.emplace_back((1.0 - prob) * Rdot * dt);
    }
    dProb.emplace_back(0.0);

    for (size_t i = 0; i < m_Rdot.size(); i++) {
      const Real time = m_Rdot[i].first;
      output << std::left << std::setw(15) << time;
      output << std::left << std::setw(15) << m_voltageCurve(time);
      output << std::left << std::setw(15) << m_maxK[i].second;
      output << std::left << std::setw(15) << m_criticalVolume[i].second;
      output << std::left << std::setw(15) << m_Rdot[i].second;
      output << std::left << std::setw(15) << m_inceptionProbability[i].second;
      output << std::left << std::setw(15) << dProb[i];
      output << "\n";
    }

    output.close();

#ifdef CH_MPI
  }
#endif
}

template <typename P, typename F, typename C>
bool
StreamerInceptionStepper<P, F, C>::particleOutsideGrid(const RealVect& a_pos,
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

template <typename P, typename F, typename C>
bool
StreamerInceptionStepper<P, F, C>::particleInsideEB(const RealVect a_pos) const noexcept
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

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::computeIonVelocity(const Real& a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::computeIonVelocity");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeIonVelocity" << endl;
  }

  CH_assert(!(m_ionSolver.isNull()));
  CH_assert(m_ionSolver->isMobile());

  EBAMRCellData& vel = m_ionSolver->getCellCenteredVelocity();

  // Set v = -E
  DataOps::copy(vel, m_electricFieldPlus);
  DataOps::scale(vel, -a_voltage);

  // Allocate mesh data that holds mu and compute it on the mesh.
  EBAMRCellData mu;
  m_amr->allocate(mu, m_realm, m_phase, 1);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {

      const EBCellFAB& v    = (*vel[lvl])[dit()];
      const FArrayBox& vReg = v.getFArrayBox();

      EBCellFAB& MU    = (*mu[lvl])[dit()];
      FArrayBox& MUREG = MU.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        const RealVect EE = RealVect(D_DECL(vReg(iv, 0), vReg(iv, 1), vReg(iv, 2)));
        const Real     E  = EE.vectorLength();

        MUREG(iv, 0) = m_ionMobility(E);
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const RealVect EE = RealVect(D_DECL(v(vof, 0), v(vof, 1), v(vof, 2)));
        const Real     E  = EE.vectorLength();

        MU(vof, 0) = m_ionMobility(E);
      };

      VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
      Box          cellBox = dbl[dit()];

      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }

  DataOps::multiplyScalar(vel, mu);

  m_amr->average(vel, m_realm, m_phase, Average::Conservative);
  m_amr->interpGhost(vel, m_realm, m_phase);
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::computeIonDiffusion(const Real& a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::computeIonDiffusion");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeIonDiffusion" << endl;
  }

  CH_assert(!(m_ionSolver.isNull()));
  CH_assert(m_ionSolver->isMobile());

  // Compute the diffusion coefficient on cell centers.
  EBAMRCellData diffCoCell;
  m_amr->allocate(diffCoCell, m_realm, m_phase, 1);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {

      const EBCellFAB& electricField    = (*m_electricFieldPlus[lvl])[dit()];
      const FArrayBox& electricFieldReg = electricField.getFArrayBox();

      EBCellFAB& dCo    = (*diffCoCell[lvl])[dit()];
      FArrayBox& dCoReg = dCo.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        const RealVect EE = RealVect(D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));
        const Real     E  = EE.vectorLength();

        dCoReg(iv, 0) = m_ionDiffusion(a_voltage * E);
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const RealVect EE = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));
        const Real     E  = EE.vectorLength();

        dCo(vof, 0) = m_ionDiffusion(a_voltage * E);
      };

      VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
      Box          cellBox = dbl[dit()];

      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }

  // Coarsen and update ghost cells
  m_amr->average(diffCoCell, m_realm, m_phase, Average::Conservative);
  m_amr->interpGhost(diffCoCell, m_realm, m_phase);

  // Fetch from solver and average to grid faces.
  EBAMRFluxData& diffCoFace = m_ionSolver->getFaceCenteredDiffusionCoefficient();
  EBAMRIVData&   diffCoEB   = m_ionSolver->getEbCenteredDiffusionCoefficient();

  DataOps::setValue(diffCoFace, 0.0);
  DataOps::setValue(diffCoEB, 0.0); // Neumann BC.

  DataOps::averageCellToFace(diffCoFace,
                             diffCoCell,
                             m_amr->getDomains(),
                             1,
                             Interval(0, 0),
                             Interval(0, 0),
                             Average::Arithmetic);
}

#include <CD_NamespaceFooter.H>

#endif
