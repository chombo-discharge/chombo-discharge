/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Copyright © 2022 NTNU.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_StreamerInceptionStepperImplem.H
  @brief  Implementation of CD_StreamerInceptionStepper.H
  @author Robert Marskar (SINTEF)
  @author Rasmus Hoholm (SINTEF)
  @author Fanny Skirbekk (NTNU)
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

  MayDay::Warning("Put electric field magnitude in output -- use MFAMRCellData for homogeneous/inhomogeneous parts");

  // Default settings.
  m_verbosity           = -1;
  m_realm               = Realm::Primal;
  m_phase               = phase::gas;
  m_mode                = Mode::Stationary;
  m_transportAlgorithm  = TransportAlgorithm::ImExCTU;
  m_timeStepRestriction = TimeStepRestriction::Unknown;
  m_profile             = false;
  m_debug               = false;
  m_fullIntegration     = false;
  m_evaluateTownsend    = false;

  this->parseOptions();

  m_rho = [](const RealVect x) -> Real {
    return 0.0;
  };

  m_sigma = [](const RealVect x) -> Real {
    return 0.0;
  };

  m_alpha = [](const Real E, const RealVect x) -> Real {
    return 0.0;
  };

  m_eta = [](const Real E, const RealVect x) -> Real {
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

  m_secondaryEmission = [](const Real E, const RealVect x) -> Real {
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
  m_tracerParticleSolver->setDeposition(DepositionType::NGP);

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
  m_amr->allocate(m_potential, m_realm, 1);
  m_amr->allocate(m_potentialHomo, m_realm, 1);
  m_amr->allocate(m_potentialInho, m_realm, 1);

  m_amr->allocate(m_electricFieldHomo, m_realm, m_phase, SpaceDim);
  m_amr->allocate(m_electricFieldInho, m_realm, m_phase, SpaceDim);

  DataOps::setValue(m_potentialHomo, 0.0);
  DataOps::setValue(m_potentialInho, 0.0);

  switch (m_mode) {
  case Mode::Stationary: {
    m_amr->allocate(m_inceptionIntegralPlus, m_realm, m_phase, m_voltageSweeps.size());
    m_amr->allocate(m_inceptionIntegralMinu, m_realm, m_phase, m_voltageSweeps.size());
    m_amr->allocate(m_inceptionVoltagePlus, m_realm, m_phase, 1);
    m_amr->allocate(m_inceptionVoltageMinu, m_realm, m_phase, 1);
    m_amr->allocate(m_streamerInceptionVoltagePlus, m_realm, m_phase, 1);
    m_amr->allocate(m_streamerInceptionVoltageMinu, m_realm, m_phase, 1);
    m_amr->allocate(m_townsendInceptionVoltagePlus, m_realm, m_phase, 1);
    m_amr->allocate(m_townsendInceptionVoltageMinu, m_realm, m_phase, 1);

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

    m_amr->allocate(m_townsendCriterionPlus, m_realm, m_phase, m_voltageSweeps.size());
    m_amr->allocate(m_townsendCriterionMinu, m_realm, m_phase, m_voltageSweeps.size());

    DataOps::setValue(m_townsendCriterionPlus, 0.0);
    DataOps::setValue(m_townsendCriterionMinu, 0.0);

    DataOps::setValue(m_inceptionIntegralPlus, 0.0);
    DataOps::setValue(m_inceptionIntegralMinu, 0.0);
    DataOps::setValue(m_inceptionVoltagePlus, 0.0);
    DataOps::setValue(m_inceptionVoltageMinu, 0.0);

    break;
  }
  case Mode::Transient: {
    m_amr->allocate(m_inceptionIntegral, m_realm, m_phase, 1);
    m_amr->allocate(m_townsendCriterion, m_realm, m_phase, 1);
    m_amr->allocate(m_emissionRate, m_realm, m_phase, 1);
    m_amr->allocate(m_backgroundIonization, m_realm, m_phase, 1);
    m_amr->allocate(m_detachment, m_realm, m_phase, 1);

    DataOps::setValue(m_inceptionIntegral, 0.0);
    DataOps::setValue(m_townsendCriterion, 0.0);
    DataOps::setValue(m_emissionRate, 0.0);
    DataOps::setValue(m_backgroundIonization, 0.0);
    DataOps::setValue(m_detachment, 0.0);

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

  this->solvePoisson();
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::solvePoisson() noexcept
{
  CH_TIME("StreamerInceptionStepper::solvePoisson");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::solvePoisson" << endl;
  }

  // Solve for the inhomogeneous part
  m_fieldSolver->setRho(m_rho);
  m_fieldSolver->setSigma(m_sigma);
  m_fieldSolver->setVoltage([](const Real& a_time) {
    return 0.0;
  });

  // Solve using our m_potentialInho as an initial guess.
  DataOps::copy(m_fieldSolver->getPotential(), m_potentialInho);
  const bool convergedInho =
    m_fieldSolver->solve(m_fieldSolver->getPotential(), m_fieldSolver->getRho(), m_fieldSolver->getSigma(), false);

  if (!convergedInho) {
    MayDay::Warning("StreamerInceptionStepper::solvePoisson -- could not solve the inhomogeneous Poisson equation. ");
  }

  DataOps::copy(m_potentialInho, m_fieldSolver->getPotential());
  DataOps::copy(m_electricFieldInho, m_amr->alias(phase::gas, m_fieldSolver->getElectricField()));

  // Solve for the homogeneous part
  m_fieldSolver->setRho([](const RealVect& a_pos) {
    return 0.0;
  });
  m_fieldSolver->setSigma([](const RealVect& a_pos) {
    return 0.0;
  });
  m_fieldSolver->setVoltage([](const Real& a_time) {
    return 1.0;
  });

  DataOps::copy(m_fieldSolver->getPotential(), m_potentialHomo);
  const bool convergedHomo =
    m_fieldSolver->solve(m_fieldSolver->getPotential(), m_fieldSolver->getRho(), m_fieldSolver->getSigma(), false);

  if (!convergedHomo) {
    MayDay::Warning("StreamerInceptionStepper::solvePoisson -- could not solve the homogeneous Poisson equation. ");
  }

  DataOps::copy(m_potentialHomo, m_fieldSolver->getPotential());
  DataOps::copy(m_electricFieldHomo, m_amr->alias(phase::gas, m_fieldSolver->getElectricField()));

  if (m_debug) {

    // Do a dummy check if the homogeneous and inhomogeneous solutions with a potential of one, and check that the
    // difference is "sufficiently small"
    const Real testVoltage = 1.0;
    const Real errorThresh = 1.E-6;

    // This is the sum of the two potentials
    DataOps::setValue(m_potential, 0.0);
    DataOps::incr(m_potential, m_potentialHomo, testVoltage);
    DataOps::incr(m_potential, m_potentialInho, 1.0);

    // Do a true solution with a test voltage and space/surface charge.
    auto voltage = [V = testVoltage](const Real& a_time) {
      return V;
    };

    DataOps::copy(m_fieldSolver->getPotential(), m_potential);

    m_fieldSolver->setRho(m_rho);
    m_fieldSolver->setSigma(m_sigma);
    m_fieldSolver->setVoltage(voltage);
    m_fieldSolver->solve(m_fieldSolver->getPotential(), m_fieldSolver->getRho(), m_fieldSolver->getSigma(), false);

    // Do the difference between the two potentials and figure out the max error. It should be < 1.E-6 * testVoltage
    DataOps::incr(m_potential, m_fieldSolver->getPotential(), -1.0);

    EBAMRCellData phiGas = m_amr->alias(phase::gas, m_potential);

    Real max;
    Real min;

    DataOps::getMaxMin(max, min, phiGas, 0);

    if (std::max(std::abs(max), std::abs(min)) > errorThresh * testVoltage) {
      MayDay::Error("StreamerInceptionStepper::solvePoisson - debug test failed. Check your BCs!");
    }
  }
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
  pp.get("profile", m_profile);
  pp.query("debug", m_debug);
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
  pp.get("eval_townsend", m_evaluateTownsend);

  // Define voltages
  const Real deltaV = voltageHi - voltageLo;

  const Real dV = deltaV / (1 + voltageSteps);

  m_voltageSweeps.push_back(voltageLo);

  for (size_t i = 0; i < voltageSteps; i++) {
    m_voltageSweeps.push_back(m_voltageSweeps[i] + dV);
  }

  m_voltageSweeps.push_back(voltageHi);
  m_maxTPlus.resize(m_voltageSweeps.size(), 0.0);
  m_maxTMinu.resize(m_voltageSweeps.size(), 0.0);
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
  pp.get("full_integration", m_fullIntegration);
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
    m_stepSizeFactor = stepSize;

    m_stepSizeMethod = StepSizeMethod::Fixed;
  }
  else if (str == "fixed") {
    m_stepSizeFactor = stepSize;

    m_stepSizeMethod = StepSizeMethod::Dx;
  }
  else if (str == "alpha") {
    m_stepSizeFactor = stepSize;

    m_stepSizeMethod = StepSizeMethod::Alpha;
  }
  else {
    MayDay::Error("Expected 'fixed', 'dx', or 'alpha' for 'StreamerInceptionStepper.inception_alg'");
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
  pp.get("first_dt", m_firstDt);
  pp.get("min_dt", m_minDt);
  pp.get("max_dt", m_maxDt);
  pp.get("max_dt_growth", m_maxDtGrowth);
  pp.get("voltage_eps", m_epsVoltage);

  CH_assert(m_minDt >= 0.0);
  CH_assert(m_maxDt >= 0.0);
  CH_assert(m_firstDt > 0.0);

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
  m_plotTownsend             = false;

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
      else if (plotVars[i] == "T") {
        m_plotTownsend = true;
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
    // Electric potential for each voltage
    ncomp += 2 * m_voltageSweeps.size();

    // Electric field magnitude
    //    ncomp =+ 2 * m_voltageSweeps.size();
    
    // K-values
    if (m_plotInceptionIntegral) {
      ncomp += 2 * m_voltageSweeps.size();
    }

    // Inception voltage.
    if (m_plotInceptionVoltage) {
      ncomp += 6;
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

    // When plotting both, add the effective ionization coefficient
    if (m_plotAlpha && m_plotEta) {
      ncomp += m_voltageSweeps.size();
    }

    // Plotting gamma coefficient
    if (m_plotTownsend) {
      ncomp += 2 * m_voltageSweeps.size();
    }

    break;
  }
  case Mode::Transient: {

    // Always write the electric field and potential.
    ncomp += SpaceDim + 1;

    if (m_plotInceptionIntegral) {
      ncomp += 1;
    }
    if (m_plotTownsend) {
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
    if (m_plotAlpha && m_plotEta) {
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
Vector<std::string>
StreamerInceptionStepper<P, F, C>::getPlotVariableNames() const
{
  CH_TIME("StreamerInceptionStepper::getPlotVariableNames");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::getPlotVariableNames" << endl;
  }

  Vector<std::string> plotVars;

  if (m_plotPoisson) {
    Vector<std::string> poissonVars = m_fieldSolver->getPlotVariableNames();
    for (int i = 0; i < poissonVars.size(); i++) {
      poissonVars[i] = "Poisson/" + poissonVars[i];
    }
    plotVars.append(poissonVars);
  }

  if (m_plotTracer) {
    Vector<std::string> tracerVars = m_tracerParticleSolver->getPlotVariableNames();
    for (int i = 0; i < tracerVars.size(); i++) {
      tracerVars[i] = "Tracer/" + tracerVars[i];
    }    
    plotVars.append(tracerVars);
  }

  if (m_plotNegativeIons) {
    Vector<std::string> cdrVars = m_ionSolver->getPlotVariableNames();
    for (int i = 0; i < cdrVars.size(); i++) {
      cdrVars[i] = "CDR/" + cdrVars[i];
    }        
    plotVars.append(cdrVars);
  }

  switch (m_mode) {
  case Mode::Stationary: {
    plotVars.append(this->getStationaryPlotVariableNames());

    break;
  }
  case Mode::Transient: {
    plotVars.append(this->getTransientPlotVariableNames());

    break;
  }
  default: {
    MayDay::Error("StreamerInceptionStepper::getPlotVariableNames - logic bust");

    break;
  }
  }

  return plotVars;
}

template <typename P, typename F, typename C>
Vector<std::string>
StreamerInceptionStepper<P, F, C>::getStationaryPlotVariableNames() const noexcept
{
  CH_TIME("StreamerInceptionStepper::getStationaryPlotVariableNames");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::getStationaryPlotVariableNames" << endl;
  }

  Vector<std::string> plotVars;

  const std::string prefix = "StreamerInceptionStepper/";

  // Always plotted.
  for(const Real& V : m_voltageSweeps) {
    plotVars.push_back(prefix + "Potential/+/ V = +" + std::to_string(V));
    plotVars.push_back(prefix + "Potential/-/ V = -" + std::to_string(V));

    // plotVars.push_back(prefix + "E/+/ V = +" + std::to_string(V));
    // plotVars.push_back(prefix + "E/-/ V = -" + std::to_string(V));        
  }

  if (m_plotInceptionVoltage) {
    plotVars.push_back(prefix + "Minimum inception voltage +");
    plotVars.push_back(prefix + "Minimum inception voltage -");
    plotVars.push_back(prefix + "Streamer inception voltage +");
    plotVars.push_back(prefix + "Streamer inception voltage -");
    plotVars.push_back(prefix + "Townsend inception voltage +");
    plotVars.push_back(prefix + "Townsend inception voltage -");
  }

  if (m_plotInceptionIntegral) {
    std::string varName;

    for (const Real& V : m_voltageSweeps) {
      varName = prefix + "K-value/+/ V = +" + std::to_string(V);
      plotVars.push_back(varName);
    }

    for (const Real& V : m_voltageSweeps) {
      varName = prefix + "K-value/-/ V = -" + std::to_string(V);
      plotVars.push_back(varName);
    }
  }

  // Secondary emission coefficient for initiatory ions
  if (m_plotTownsend) {
    for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
      const std::string varName = prefix + "T-value/+/ V = +" + std::to_string(m_voltageSweeps[i]);

      plotVars.push_back(varName);
    }

    for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
      const std::string varName = prefix + "T-value/-/ V = -" + std::to_string(m_voltageSweeps[i]);

      plotVars.push_back(varName);
    }
  }

  // Write the background ionization rates
  if (m_plotBackgroundIonization) {
    std::string varName;

    for (const Real& V : m_voltageSweeps) {
      varName = prefix + "Background ionization rate/ V =  " + std::to_string(V);
      plotVars.push_back(varName);
    }
  }

  // Write detachment rates
  if (m_plotDetachment) {
    std::string varName;

    for (const Real& V : m_voltageSweeps) {
      varName = prefix + "Detachment rate/ V =  " + std::to_string(V);
      plotVars.push_back(varName);
    }
  }

  // Write field emission rates
  if (m_plotFieldEmission) {
    std::string varName;

    for (const Real& V : m_voltageSweeps) {
      varName = prefix + "Field emission rate/+/ V = +" + std::to_string(V);
      plotVars.push_back(varName);
    }

    for (const Real& V : m_voltageSweeps) {
      varName = prefix + "Field emission rate/-/ V = -" + std::to_string(V);
      plotVars.push_back(varName);
    }
  }

  // Alpha value
  if (m_plotAlpha) {
    for (size_t i = 0; i < m_voltageSweeps.size(); i++) {

      const std::string varName = prefix + "Alpha coefficient/ V = " + std::to_string(m_voltageSweeps[i]);

      plotVars.push_back(varName);
    }
  }

  // Eta
  if (m_plotEta) {
    for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
      const std::string varName = prefix + "Eta coefficient/ V = " + std::to_string(m_voltageSweeps[i]);

      plotVars.push_back(varName);
    }
  }

  // Effective alpha coefficient
  if (m_plotAlpha && m_plotEta) {
    for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
      const std::string varName = prefix + "Alpha_eff/ V = " + std::to_string(m_voltageSweeps[i]);

      plotVars.push_back(varName);
    }
  }

  return plotVars;
}

template <typename P, typename F, typename C>
Vector<std::string>
StreamerInceptionStepper<P, F, C>::getTransientPlotVariableNames() const noexcept
{
  CH_TIME("StreamerInceptionStepper::getTransientPlotVariableNames");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::getTransientPlotVariableNames" << endl;
  }

  Vector<std::string> plotVars;

  // Add potential and electric field to output
  plotVars.push_back("Electric potential");
  plotVars.push_back("x-E");
  plotVars.push_back("y-E");
  if (SpaceDim == 3) {
    plotVars.push_back("z-E");
  }

  if (m_plotInceptionIntegral) {
    plotVars.push_back("Inception integral");
  }

  if (m_plotTownsend) {
    plotVars.push_back("Townsend criterion");
  }

  if (m_plotBackgroundIonization) {
    plotVars.push_back("Background rate");
  }

  if (m_plotDetachment) {
    plotVars.push_back("Detachment rate");
  }

  if (m_plotFieldEmission) {
    plotVars.push_back("Field emission");
  }

  if (m_plotAlpha) {
    plotVars.push_back("Townsend alpha coefficient");
  }

  if (m_plotEta) {
    plotVars.push_back("Townsend eta coefficient");
  }

  if (m_plotAlpha && m_plotEta) {
    plotVars.push_back("Effective Townsend coefficient");
  }

  return plotVars;
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::writePlotData(LevelData<EBCellFAB>& a_output,
                                                 int&                  a_icomp,
                                                 const std::string     a_outputRealm,
                                                 const int             a_level) const
{
  CH_TIME("StreamerInceptionStepper::writePlotData");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::writePlotData" << endl;
  }

  if (m_plotPoisson) {
    m_fieldSolver->writePlotData(a_output, a_icomp, a_outputRealm, a_level);
  }

  if (m_plotTracer) {
    m_tracerParticleSolver->writePlotData(a_output, a_icomp, a_outputRealm, a_level);
  }

  if (m_plotNegativeIons) {
    m_ionSolver->writePlotData(a_output, a_icomp, a_outputRealm, a_level);
  }

  // Write internal data -- this differs between the 'stationary' and 'transient' modes
  switch (m_mode) {
  case Mode::Stationary: {
    this->writePlotDataStationary(a_output, a_icomp, a_outputRealm, a_level);

    break;
  }
  case Mode::Transient:
    this->writePlotDataTransient(a_output, a_icomp, a_outputRealm, a_level);

    break;
  default: {
    MayDay::Error("StreamerInceptionStepper::writePlotData - logic bust");

    break;
  }
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::writePlotDataStationary(LevelData<EBCellFAB>& a_output,
                                                           int&                  a_icomp,
                                                           const std::string     a_outputRealm,
                                                           const int             a_level) const noexcept
{
  CH_TIME("StreamerInceptionStepper::writePlotDataStationary");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::writePlotDataStationary" << endl;
  }

  const std::string prefix = "StreamerInceptionStepper/";

  for (const auto& V : m_voltageSweeps) {
    DataOps::setValue(m_potential, 0.0);
    DataOps::incr(m_potential, m_potentialHomo, Real(V));
    DataOps::incr(m_potential, m_potentialInho, 1.0);
    
    m_fieldSolver
    ->writeMultifluidData(a_output, a_icomp, m_potential, phase::gas, a_outputRealm, a_level, true);
    
    DataOps::setValue(m_potential, 0.0);
    DataOps::incr(m_potential, m_potentialHomo, -Real(V));
    DataOps::incr(m_potential, m_potentialInho, 1.0);
    
    m_fieldSolver
    ->writeMultifluidData(a_output, a_icomp, m_potential, phase::gas, a_outputRealm, a_level, true);
  }   

  if (m_plotInceptionVoltage) {
    this->writeData(a_output, a_icomp, m_inceptionVoltagePlus, a_outputRealm, a_level, false, true);
    this->writeData(a_output, a_icomp, m_inceptionVoltageMinu, a_outputRealm, a_level, false, true);
    this->writeData(a_output, a_icomp, m_streamerInceptionVoltagePlus, a_outputRealm, a_level, false, true);
    this->writeData(a_output, a_icomp, m_streamerInceptionVoltageMinu, a_outputRealm, a_level, false, true);
    this->writeData(a_output, a_icomp, m_townsendInceptionVoltagePlus, a_outputRealm, a_level, false, true);
    this->writeData(a_output, a_icomp, m_townsendInceptionVoltageMinu, a_outputRealm, a_level, false, true);
  }

  if (m_plotInceptionIntegral) {
    this->writeData(a_output, a_icomp, m_inceptionIntegralPlus, a_outputRealm, a_level, false, true);
    this->writeData(a_output, a_icomp, m_inceptionIntegralMinu, a_outputRealm, a_level, false, true);
  }

  if (m_plotTownsend) {
    this->writeData(a_output, a_icomp, m_townsendCriterionPlus, a_outputRealm, a_level, false, true);
    this->writeData(a_output, a_icomp, m_townsendCriterionMinu, a_outputRealm, a_level, false, true);
  }

  if (m_plotBackgroundIonization) {
    this->writeData(a_output, a_icomp, m_backgroundIonization, a_outputRealm, a_level, false, true);
  }

  if (m_plotDetachment) {
    this->writeData(a_output, a_icomp, m_detachment, a_outputRealm, a_level, false, true);
  }

  if (m_plotFieldEmission) {
    this->writeData(a_output, a_icomp, m_emissionRatesPlus, a_outputRealm, a_level, false, true);
    this->writeData(a_output, a_icomp, m_emissionRatesMinu, a_outputRealm, a_level, false, true);
  }

  if (m_plotAlpha) {
    LevelData<EBCellFAB> alpha;
    LevelData<EBCellFAB> alphaCoar;

    m_amr->allocate(alpha, m_realm, m_phase, a_level, 1);
    if (a_level > 0) {
      m_amr->allocate(alphaCoar, m_realm, m_phase, a_level - 1, 1);
    }

    for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
      this->evaluateFunction(alpha, m_voltageSweeps[i], m_alpha, a_level);
      if (a_level > 0) {
        this->evaluateFunction(alphaCoar, m_voltageSweeps[i], m_alpha, a_level - 1);
        m_amr->interpGhost(alpha, alphaCoar, a_level, m_realm, m_phase);
      }
      else {
        alpha.exchange();
      }

      m_amr->copyData(a_output, alpha, a_level, a_outputRealm, m_realm, Interval(a_icomp, a_icomp), Interval(0, 0));

      a_icomp++;
    }
  }

  if (m_plotEta) {
    LevelData<EBCellFAB> eta;
    LevelData<EBCellFAB> etaCoar;

    m_amr->allocate(eta, m_realm, m_phase, a_level, 1);
    if (a_level > 0) {
      m_amr->allocate(etaCoar, m_realm, m_phase, a_level - 1, 1);
    }

    for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
      this->evaluateFunction(eta, m_voltageSweeps[i], m_eta, a_level);
      if (a_level > 0) {
        this->evaluateFunction(etaCoar, m_voltageSweeps[i], m_eta, a_level - 1);
        m_amr->interpGhost(eta, etaCoar, a_level, m_realm, m_phase);
      }
      else {
        eta.exchange();
      }

      m_amr->copyData(a_output, eta, a_level, a_outputRealm, m_realm, Interval(a_icomp, a_icomp), Interval(0, 0));

      a_icomp++;
    }
  }

  if (m_plotAlpha && m_plotEta) {
    LevelData<EBCellFAB> alphaEff;
    LevelData<EBCellFAB> alphaEffCoar;

    m_amr->allocate(alphaEff, m_realm, m_phase, a_level, 1);
    if (a_level > 0) {
      m_amr->allocate(alphaEffCoar, m_realm, m_phase, a_level - 1, 1);
    }

    auto alphaFunc = [alpha = this->m_alpha, eta = this->m_eta](const Real E, const RealVect x) -> Real {
      return alpha(E, x) - eta(E, x);
    };

    for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
      this->evaluateFunction(alphaEff, m_voltageSweeps[i], alphaFunc, a_level);
      if (a_level > 0) {
        this->evaluateFunction(alphaEffCoar, m_voltageSweeps[i], alphaFunc, a_level - 1);
        m_amr->interpGhost(alphaEff, alphaEffCoar, a_level, m_realm, m_phase);
      }
      else {
        alphaEff.exchange();
      }

      // Do a copy, but note that this does not include ghost cells across the CF interface since
      // we are only computing for a single grid level.
      m_amr->copyData(a_output, alphaEff, a_level, a_outputRealm, m_realm, Interval(a_icomp, a_icomp), Interval(0, 0));

      a_icomp++;
    }
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::writePlotDataTransient(LevelData<EBCellFAB>& a_output,
                                                          int&                  a_icomp,
                                                          const std::string     a_outputRealm,
                                                          const int             a_level) const noexcept
{
  CH_TIME("StreamerInceptionStepper::writePlotDataTransient");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::writePlotDataTransient" << endl;
  }

  CH_assert(a_level >= 0);
  CH_assert(a_level <= m_amr->getFinestLevel());

  // Add potential and electric field to output. Note that we need to scale appropriately.
  LevelData<EBCellFAB> aliasData;
  aliasLevelData<EBCellFAB>(aliasData, &a_output, Interval(a_icomp, a_icomp + SpaceDim));

  m_fieldSolver
    ->writeMultifluidData(a_output, a_icomp, m_fieldSolver->getPotential(), phase::gas, a_outputRealm, a_level, true);
  m_fieldSolver->writeMultifluidData(a_output,
                                     a_icomp,
                                     m_fieldSolver->getElectricField(),
                                     phase::gas,
                                     a_outputRealm,
                                     a_level,
                                     true);

  DataOps::scale(aliasData, m_voltageCurve(m_time));

  if (m_plotInceptionIntegral) {
    this->writeData(a_output, a_icomp, m_inceptionIntegral, a_outputRealm, a_level, false, true);
  }

  if (m_plotTownsend) {
    this->writeData(a_output, a_icomp, m_townsendCriterion, a_outputRealm, a_level, false, true);
  }

  if (m_plotBackgroundIonization) {
    LevelData<EBCellFAB> bgIonization;
    m_amr->allocate(bgIonization, m_realm, m_phase, a_level, 1);

    this->evaluateFunction(bgIonization, m_voltageCurve(m_time), m_backgroundRate, a_level);

    // Do a copy, but note that this does not include ghost cells across the CF interface since
    // we are only computing for a single grid level.
    m_amr
      ->copyData(a_output, bgIonization, a_level, a_outputRealm, m_realm, Interval(a_icomp, a_icomp), Interval(0, 0));

    a_icomp++;
  }

  if (m_plotDetachment) {
    LevelData<EBCellFAB> detachRate;
    LevelData<EBCellFAB> detachRateCoar;

    m_amr->allocate(detachRate, m_realm, m_phase, a_level, 1);
    if (a_level > 0) {
      m_amr->allocate(detachRateCoar, m_realm, m_phase, a_level - 1, 1);
    }

    this->evaluateFunction(detachRate, m_voltageCurve(m_time), m_detachmentRate, a_level);
    if (a_level > 0) {
      this->evaluateFunction(detachRateCoar, m_voltageCurve(m_time), m_detachmentRate, a_level - 1);
      m_amr->interpGhost(detachRate, detachRateCoar, a_level, m_realm, m_phase);
    }

    m_amr->copyData(a_output, detachRate, a_level, a_outputRealm, m_realm, Interval(a_icomp, a_icomp), Interval(0, 0));

    a_icomp++;
  }

  if (m_plotFieldEmission) {
    this->writeData(a_output, a_icomp, m_emissionRate, a_outputRealm, a_level, false, false);
  }

  if (m_plotAlpha) {
    LevelData<EBCellFAB> alpha;
    LevelData<EBCellFAB> alphaCoar;

    m_amr->allocate(alpha, m_realm, m_phase, a_level, 1);
    if (a_level > 0) {
      m_amr->allocate(alphaCoar, m_realm, m_phase, a_level - 1, 1);
    }

    this->evaluateFunction(alpha, m_voltageCurve(m_time), m_alpha, a_level);
    if (a_level > 0) {
      this->evaluateFunction(alphaCoar, m_voltageCurve(m_time), m_alpha, a_level - 1);
      m_amr->interpGhost(alpha, alphaCoar, a_level, m_realm, m_phase);
    }
    else {
      alpha.exchange();
    }

    // Do a copy, but note that this does not include ghost cells across the CF interface since
    // we are only computing for a single grid level.
    m_amr->copyData(a_output, alpha, a_level, a_outputRealm, m_realm, Interval(a_icomp, a_icomp), Interval(0, 0));

    a_icomp++;
  }

  if (m_plotEta) {
    LevelData<EBCellFAB> eta;
    LevelData<EBCellFAB> etaCoar;

    m_amr->allocate(eta, m_realm, m_phase, a_level, 1);
    if (a_level > 0) {
      m_amr->allocate(etaCoar, m_realm, m_phase, a_level - 1, 1);
    }

    this->evaluateFunction(eta, m_voltageCurve(m_time), m_eta, a_level);
    if (a_level > 0) {
      this->evaluateFunction(etaCoar, m_voltageCurve(m_time), m_eta, a_level - 1);
      m_amr->interpGhost(eta, etaCoar, a_level, m_realm, m_phase);
    }
    else {
      eta.exchange();
    }

    // Do a copy, but note that this does not include ghost cells across the CF interface since
    // we are only computing for a single grid level.
    m_amr->copyData(a_output, eta, a_level, a_outputRealm, m_realm, Interval(a_icomp, a_icomp), Interval(0, 0));

    a_icomp++;
  }

  if (m_plotAlpha && m_plotEta) {
    LevelData<EBCellFAB> alphaEff;
    LevelData<EBCellFAB> alphaEffCoar;

    m_amr->allocate(alphaEff, m_realm, m_phase, a_level, 1);
    if (a_level > 0) {
      m_amr->allocate(alphaEffCoar, m_realm, m_phase, a_level - 1, 1);
    }

    auto alphaFunc = [alpha = this->m_alpha, eta = this->m_eta](const Real E, const RealVect x) -> Real {
      return alpha(E, x) - eta(E, x);
    };

    this->evaluateFunction(alphaEff, m_voltageCurve(m_time), alphaFunc, a_level);
    if (a_level > 0) {
      this->evaluateFunction(alphaEffCoar, m_voltageCurve(m_time), alphaFunc, a_level - 1);
      m_amr->interpGhost(alphaEff, alphaEffCoar, a_level, m_realm, m_phase);
    }
    else {
      alphaEff.exchange();
    }

    // Do a copy, but note that this does not include ghost cells across the CF interface since
    // we are only computing for a single grid level.
    m_amr->copyData(a_output, alphaEff, a_level, a_outputRealm, m_realm, Interval(a_icomp, a_icomp), Interval(0, 0));

    a_icomp++;
  }
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

    // Compute time step from ion solver.
    m_timeStepRestriction = TimeStepRestriction::Unknown;

    // Restrict by ion solver
    if (m_ionTransport) {
      Real ionDt = std::numeric_limits<Real>::infinity();
      switch (m_transportAlgorithm) {
      case TransportAlgorithm::Euler: {
        ionDt = m_cfl * m_ionSolver->computeAdvectionDiffusionDt();

        break;
      }
      case TransportAlgorithm::Heun: {
        ionDt = m_cfl * m_ionSolver->computeAdvectionDiffusionDt();

        break;
      }
      case TransportAlgorithm::ImExCTU: {
        ionDt = m_cfl * m_ionSolver->computeAdvectionDt();

        break;
      }
      default: {
        MayDay::Error("StreamerInceptionStepper::computDt -- logic bust");

        break;
      }
      }

      if (ionDt < dt) {
        dt                    = ionDt;
        m_timeStepRestriction = TimeStepRestriction::CDR;
      }
    }

    // Restrict by voltage curve.
    Real curveDt = std::numeric_limits<Real>::max();
    if (m_timeStep == 0) {
      curveDt = m_firstDt;
    }
    else {
      const Real preU = m_voltageCurve(m_time - m_dt);
      const Real curU = m_voltageCurve(m_time);
      const Real dVdt = std::abs(curU - preU) / dt;

      if (dVdt > std::numeric_limits<Real>::epsilon()) {
        curveDt = m_epsVoltage * std::abs(curU) / dVdt;
      }

      // Do not grow too fast.
      curveDt = std::min(m_dt * (1.0 + m_maxDtGrowth), curveDt);
      curveDt = std::max(m_dt * (1.0 - m_maxDtGrowth), curveDt);
    }

    if (curveDt < dt) {
      dt                    = curveDt;
      m_timeStepRestriction = TimeStepRestriction::VoltageCurve;
    }

    // Restrict by hardcaps
    if (dt > m_maxDt) {
      dt                    = m_maxDt;
      m_timeStepRestriction = TimeStepRestriction::MaxHardcap;
    }

    if (dt < m_minDt) {
      dt                    = m_minDt;
      m_timeStepRestriction = TimeStepRestriction::MinHardcap;
    }
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
    Timer timer("StreamerInceptionStepper::advance");

    const Real curTime    = m_time + a_dt;
    const Real curVoltage = m_voltageCurve(curTime);

    // Move the negative ions -- still want this for the time step.
    timer.startEvent("Ion advance");
    if (m_ionTransport) {
      this->computeIonVelocity(curVoltage);
      this->computeIonDiffusion(curVoltage);
      this->advanceIons(a_dt);
    }
    timer.stopEvent("Ion advance");

    // Seed new particles and compute the inception integral.
    timer.startEvent("Inception integral");
    this->seedIonizationParticles(curVoltage);
    this->computeInceptionIntegralTransient(curVoltage);
    timer.stopEvent("Inception integral");

    if (m_evaluateTownsend) {
      timer.startEvent("Townsend criterion");
      this->seedIonizationParticles(curVoltage);
      this->computeTownsendCriterionTransient(curVoltage);
      timer.stopEvent("Townsend criterion");
    }

    // Compute the critical volume, critical area, the ionization volume and the electron appearance rate.
    timer.startEvent("Compute Vcr");
    const Real Vcr  = this->computeCriticalVolumeTransient();
    const Real Acr  = this->computeCriticalAreaTransient();
    const Real Vion = this->computeIonizationVolumeTransient(curVoltage);
    const Real Rdot = this->computeRdot(curVoltage);
    timer.stopEvent("Compute Vcr");

    m_criticalVolume.emplace_back(curTime, Vcr);
    m_criticalArea.emplace_back(curTime, Acr);
    m_ionizationVolumeTransient.emplace_back(curTime, Vion);
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

    // Get the maximum K and T values
    timer.startEvent("Get max/min K");
    Real maxK = -std::numeric_limits<Real>::max();
    Real minK = +std::numeric_limits<Real>::max();

    Real maxT = -std::numeric_limits<Real>::max();
    Real minT = +std::numeric_limits<Real>::max();

    DataOps::getMaxMin(maxK, minK, m_inceptionIntegral, 0);
    DataOps::getMaxMin(maxT, minT, m_townsendCriterion, 0);

    if (!m_fullIntegration) {
      maxK = std::min(maxK, m_inceptionK);
      maxT = std::min(maxT, 1.0);
    }
    m_maxK.emplace_back(m_time + a_dt, maxK);
    m_maxT.emplace_back(m_time + a_dt, maxT);
    timer.stopEvent("Get max/min K");

    timer.startEvent("Write report");
    this->writeReportTransient();
    timer.stopEvent("Write report");

    if (m_profile) {
      timer.eventReport(pout(), false);
    }
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

  std::string timeStepMessage;
  switch (m_timeStepRestriction) {
  case TimeStepRestriction::Unknown: {
    timeStepMessage = "Unknown";

    break;
  }
  case TimeStepRestriction::CDR: {
    timeStepMessage = "CFL";

    break;
  }
  case TimeStepRestriction::VoltageCurve: {
    timeStepMessage = "Voltage curve";

    break;
  }
  case TimeStepRestriction::MinHardcap: {
    timeStepMessage = "Min hardcap";

    break;
  }
  case TimeStepRestriction::MaxHardcap: {
    timeStepMessage = "Max hardcap";

    break;
  }
  default: {
    MayDay::Warning("StreamerInceptionStepper::printStepReport - logic bust");

    break;
  }
  }

  // clang-format off
  pout() << "                                ** Voltage                = " << m_voltageCurve(m_time) << endl;
  pout() << "                                ** Crit. volume           = " << m_criticalVolume.back().second << endl;
  pout() << "                                ** Inception probability  = " << m_inceptionProbability.back().second << endl;
  pout() << "                                ** Time step restriction  = " << "'" << timeStepMessage << "'" << endl;
  // clang-format on
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

  m_amr->allocate(m_scratchHomo, m_realm, 1);
  m_amr->allocate(m_scratchInho, m_realm, 1);

  m_amr->copyData(m_scratchHomo, m_potentialHomo);
  m_amr->copyData(m_scratchInho, m_potentialInho);
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

  // Regrid electric field and potentials
  m_amr->reallocate(m_potential, a_lmin);
  m_amr->reallocate(m_potentialHomo, a_lmin);
  m_amr->reallocate(m_potentialInho, a_lmin);

  m_amr->reallocate(m_electricFieldHomo, m_phase, a_lmin);
  m_amr->reallocate(m_electricFieldInho, m_phase, a_lmin);

  const EBCoarseToFineInterp::Type interpType = EBCoarseToFineInterp::ConservativeMinMod;

  m_amr->interpToNewGrids(m_potentialHomo, m_scratchHomo, a_lmin, a_oldFinestLevel, a_newFinestLevel, interpType);
  m_amr->interpToNewGrids(m_potentialInho, m_scratchInho, a_lmin, a_oldFinestLevel, a_newFinestLevel, interpType);

  switch (m_mode) {
  case Mode::Stationary: {
    m_amr->reallocate(m_inceptionIntegralPlus, m_phase, a_lmin);
    m_amr->reallocate(m_inceptionIntegralMinu, m_phase, a_lmin);
    m_amr->reallocate(m_inceptionVoltagePlus, m_phase, a_lmin);
    m_amr->reallocate(m_inceptionVoltageMinu, m_phase, a_lmin);
    m_amr->reallocate(m_streamerInceptionVoltagePlus, m_phase, a_lmin);
    m_amr->reallocate(m_streamerInceptionVoltageMinu, m_phase, a_lmin);
    m_amr->reallocate(m_townsendInceptionVoltagePlus, m_phase, a_lmin);
    m_amr->reallocate(m_townsendInceptionVoltageMinu, m_phase, a_lmin);

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

    m_amr->reallocate(m_townsendCriterionPlus, m_phase, a_lmin);
    m_amr->reallocate(m_townsendCriterionMinu, m_phase, a_lmin);

    break;
  }
  case Mode::Transient: {

    m_amr->reallocate(m_inceptionIntegral, m_phase, a_lmin);
    m_amr->reallocate(m_townsendCriterion, m_phase, a_lmin);
    m_amr->reallocate(m_emissionRate, m_phase, a_lmin);
    m_amr->reallocate(m_backgroundIonization, m_phase, a_lmin);

    break;
  }
  default: {
    break;
  }
  }

  // Solve for the inhomogeneous and homogeneous Poisson fields again
  this->solvePoisson();

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

  m_scratchHomo.clear();
  m_scratchInho.clear();
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
StreamerInceptionStepper<P, F, C>::setSigma(const std::function<Real(const RealVect& x)>& a_sigma) noexcept
{
  CH_TIME("StreamerInceptionStepper::setSigma");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setSigma" << endl;
  }

  m_sigma = a_sigma;
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
StreamerInceptionStepper<P, F, C>::setAlpha(
  const std::function<Real(const Real& E, const RealVect& x)>& a_alpha) noexcept
{
  CH_TIME("StreamerInceptionStepper::setAlpha");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setAlpha" << endl;
  }

  m_alpha = a_alpha;
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::setEta(const std::function<Real(const Real& E, const RealVect& x)>& a_eta) noexcept
{
  CH_TIME("StreamerInceptionStepper::setEta");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setEta" << endl;
  }

  m_eta = a_eta;
}

template <typename P, typename F, typename C>
const std::function<Real(const Real& E, const RealVect& x)>&
StreamerInceptionStepper<P, F, C>::getAlpha() const noexcept
{
  return m_alpha;
}

template <typename P, typename F, typename C>
const std::function<Real(const Real& E, const RealVect& x)>&
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
void
StreamerInceptionStepper<P, F, C>::setSecondaryEmission(
  const std::function<Real(const Real& E, const RealVect& x)>& a_coefficient) noexcept
{
  CH_TIME("StreamerInceptionStepper::setSecondaryEmission");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setSecondaryEmission" << endl;
  }

  m_secondaryEmission = a_coefficient;
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
StreamerInceptionStepper<P, F, C>::seedIonizationParticles(const Real a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::seedIonizationParticles");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::seedIonizationParticles" << endl;
  }

  ParticleContainer<P>& amrParticles = m_tracerParticleSolver->getParticles();
  amrParticles.clearParticles();

  // Compute the electric field at the input voltage
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, SpaceDim);
  this->superposition(scratch, a_voltage);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

    const LevelData<BaseFab<bool>>& validCellsLD = *m_amr->getValidCells(m_realm)[lvl];

    const Real     dx     = m_amr->getDx()[lvl];
    const RealVect probLo = m_amr->getProbLo();

    ParticleData<P>& levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const EBISBox&       ebisbox    = ebisl[dit()];
      const BaseFab<bool>& validCells = validCellsLD[dit()];

      List<P>& particles = levelParticles[dit()].listItems();

      const EBCellFAB& electricField    = (*scratch[lvl])[dit()];
      const FArrayBox& electricFieldReg = electricField.getFArrayBox();

      if (!ebisbox.isAllCovered()) {

        auto regularKernel = [&](const IntVect& iv) -> void {
          if (validCells(iv, 0) && ebisbox.isRegular(iv)) {

            const RealVect x  = probLo + dx * (0.5 * RealVect::Unit + RealVect(iv));
            const RealVect EE = RealVect(
              D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));

            const Real E     = EE.vectorLength();
            const Real alpha = m_alpha(E, x);
            const Real eta   = m_eta(E, x);

            if (alpha > eta) {
              P p;

              p.position() = x;

              particles.add(p);
            }
          }
        };

        auto irregularKernel = [&](const VolIndex& vof) -> void {
          const IntVect iv = vof.gridIndex();

          if (validCells(iv, 0) && ebisbox.isIrregular(iv)) {

            const RealVect x  = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);
            const RealVect EE = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));
            const Real     E  = EE.vectorLength();

            const Real alpha = m_alpha(E, x);
            const Real eta   = m_eta(E, x);

            if (alpha > eta) {
              P p;

              p.position() = x;

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
    if (m_evaluateTownsend) {
      this->computeTownsendCriterionStationary();
    }
    this->computeCriticalVolumeStationary();
    this->computeCriticalAreaStationary();
    this->computeIonizationVolumeStationary();
    this->computeInceptionVoltageVolume();

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
    this->computeInceptionIntegralTransient(m_voltageCurve(m_time));

    // Evaluate the Townsend criterion if the user asks for it.
    if (m_evaluateTownsend) {
      this->seedIonizationParticles(m_voltageCurve(m_time));
      this->computeTownsendCriterionTransient(m_voltageCurve(m_time));
    }

    // Get the maximum K and T values
    Real maxK = -std::numeric_limits<Real>::max();
    Real minK = +std::numeric_limits<Real>::max();

    Real maxT = -std::numeric_limits<Real>::max();
    Real minT = +std::numeric_limits<Real>::max();

    DataOps::getMaxMin(maxK, minK, m_inceptionIntegral, 0);
    DataOps::getMaxMin(maxT, minT, m_townsendCriterion, 0);
    if (!m_fullIntegration) {
      maxK = std::min(maxK, m_inceptionK);
      maxT = std::min(maxT, 1.0);
    }

    m_Rdot.emplace_back(m_time, this->computeRdot(m_voltageCurve(m_time)));
    m_criticalVolume.emplace_back(m_time, this->computeCriticalVolumeTransient());
    m_criticalArea.emplace_back(m_time, this->computeCriticalAreaTransient());
    m_ionizationVolumeTransient.emplace_back(m_time, this->computeIonizationVolumeTransient(m_voltageCurve(m_time)));
    m_inceptionProbability.emplace_back(m_time, 0.0);
    m_maxK.emplace_back(m_time, maxK);
    m_maxT.emplace_back(m_time, maxT);

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
  EBAMRCellData Kplus;
  EBAMRCellData Kminu;

  m_amr->allocate(Kplus, m_realm, m_phase, 1);
  m_amr->allocate(Kminu, m_realm, m_phase, 1);

  // Integrate for each voltage
  for (int i = 0; i < m_voltageSweeps.size(); ++i) {

    DataOps::setValue(Kplus, 0.0);
    DataOps::setValue(Kminu, 0.0);

    // Polarities. Note that we are dealing with electrons so
    // for positive polarity the particles move opposite to the field.
    std::vector<Real> polarities{1.0, -1.0};

    this->seedIonizationParticles(m_voltageSweeps[i]);

    for (const auto& p : polarities) {

      this->resetTracerParticles();

      // Switch between various integration algorithms. At the end of this
      // the K-value time should be stored on the particle weight and the particles
      // should be back in their original grid cells.
      switch (m_inceptionAlgorithm) {
      case IntegrationAlgorithm::Euler: {
        this->inceptionIntegrateEuler(p * m_voltageSweeps[i]);

        break;
      }
      case IntegrationAlgorithm::Trapezoidal: {
        this->inceptionIntegrateTrapezoidal(p * m_voltageSweeps[i]);

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
        if (!m_fullIntegration) {
          maxK = std::min(maxK, m_inceptionK);
        }

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
        if (!m_fullIntegration) {
          maxK = std::min(maxK, m_inceptionK);
        }

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

  size_t outcastsBefore  = 0;
  size_t particlesBefore = 0;

  if (m_debug) {
    outcastsBefore  = amrParticles.getNumberOfOutcastParticlesGlobal();
    particlesBefore = amrParticles.getNumberOfValidParticlesGlobal();
  }

  // Allocate something that holds the velocity of the electrons
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, SpaceDim);
  this->superposition(scratch, a_voltage);
  DataOps::scale(scratch, -1.0);

  m_tracerParticleSolver->setVelocity(scratch);
  m_tracerParticleSolver->interpolateVelocities();

  while (amrParticles.getNumberOfValidParticlesGlobal() > 0) {

    // Euler integration.
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const Real dx = m_amr->getDx()[lvl];

      Real spaceStep = std::numeric_limits<Real>::max();
      Real alphaStep = std::numeric_limits<Real>::max();

      switch (m_stepSizeMethod) {
      case StepSizeMethod::Fixed: {
        spaceStep = m_stepSizeFactor;

        break;
      }
      case StepSizeMethod::Dx: {
        spaceStep = m_stepSizeFactor * dx;

        break;
      }
      case StepSizeMethod::Alpha: {
        spaceStep = m_stepSizeFactor * dx;
        alphaStep = m_stepSizeFactor;

        break;
      }
      default: {
        MayDay::Error("StreamerInceptionStepper::inceptionIntegrateEuler - logic bust");

        break;
      }
      }

      for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
        List<P>& solverParticles    = amrParticles[lvl][dit()].listItems();
        List<P>& processedParticles = amrProcessedParticles[lvl][dit()].listItems();

        for (ListIterator<P> lit(solverParticles); lit.ok();) {
          P& p = lit();

          const RealVect x        = p.position();
          const RealVect vel      = p.velocity();
          const Real     v        = vel.vectorLength();
          const Real     E        = v;
          const Real     alpha    = m_alpha(E, x);
          const Real     eta      = m_eta(E, x);
          const Real     alphaEff = alpha - eta;

          const Real     deltaX = std::min(spaceStep, alphaStep / alphaEff);
          const Real     dt     = deltaX / v;
          const RealVect newPos = p.position() + dt * vel;
          const Real     delta  = (newPos - x).vectorLength();

          const bool outsideDomain = this->particleOutsideGrid(newPos, probLo, probHi);
          const bool insideEB      = this->particleInsideEB(newPos);

          // If the particle struck the EB or domain we finish off the integration with a partial step
          Real s = 0.0;

          if (alphaEff < 0.0) {
            processedParticles.transfer(lit);
          }
          else if (insideEB) {
            const RefCountedPtr<BaseIF>& impFunc = m_amr->getBaseImplicitFunction(m_phase);

            if (ParticleOps::ebIntersectionBisect(impFunc, x, newPos, spaceStep, s)) {
              p.weight() = p.weight() + s * delta * alphaEff;
            }

            processedParticles.transfer(lit);
          }
          else if (outsideDomain) {
            if (ParticleOps::domainIntersection(x, newPos, m_amr->getProbLo(), m_amr->getProbHi(), s)) {
              p.weight() = p.weight() + s * delta * alphaEff;
            }

            processedParticles.transfer(lit);
          }

          else {
            p.position() = newPos;
            p.weight()   = p.weight() + deltaX * alphaEff;

            ++lit;
          }
        }

        // Transfer particles that completed their integration.
        if (!m_fullIntegration) {
          for (ListIterator<P> lit(solverParticles); lit.ok();) {
            P& p = lit();

            if (p.weight() >= m_inceptionK) {
              processedParticles.transfer(lit);
            }
            else {
              ++lit;
            }
          }
        }
      }
    }

    // Update velocities.
    m_tracerParticleSolver->remap();
    m_tracerParticleSolver->interpolateVelocities();
  }

  ParticleOps::copyDestructive(amrParticles, amrProcessedParticles);

  // Truncate the weights if we didn't run full integration
  if (!m_fullIntegration) {
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
        List<P>& solverParticles = amrParticles[lvl][dit()].listItems();

        for (ListIterator<P> lit(solverParticles); lit.ok(); ++lit) {
          lit().weight() = std::min(m_inceptionK, lit().weight());
        }
      }
    }
  }

  this->rewindTracerParticles();

  size_t outcastsAfter  = 0;
  size_t particlesAfter = 0;

  if (m_debug) {
    outcastsAfter  = amrParticles.getNumberOfOutcastParticlesGlobal();
    particlesAfter = amrParticles.getNumberOfValidParticlesGlobal();
  }

  if (particlesBefore != particlesAfter) {
    MayDay::Warning("StreamerInceptionStepper::inceptionIntegrateEuler - lost/gained particles!");
  }
  if (outcastsBefore != outcastsAfter) {
    MayDay::Warning("StreamerInceptionStepper::inceptionIntegrateEuler - lost/gained outcast particles!");
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

  size_t outcastsBefore  = 0;
  size_t particlesBefore = 0;

  if (m_debug) {
    outcastsBefore  = amrParticles.getNumberOfOutcastParticlesGlobal();
    particlesBefore = amrParticles.getNumberOfValidParticlesGlobal();
  }

  // Allocate something that holds the velocity of the electrons
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, SpaceDim);
  this->superposition(scratch, a_voltage);
  DataOps::scale(scratch, -1.0);

  m_tracerParticleSolver->setVelocity(scratch);
  m_tracerParticleSolver->interpolateVelocities();

  while (amrParticles.getNumberOfValidParticlesGlobal() > 0) {

    // Euler stage.
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const Real dx = m_amr->getDx()[lvl];

      Real spaceStep = std::numeric_limits<Real>::max();
      Real alphaStep = std::numeric_limits<Real>::max();

      switch (m_stepSizeMethod) {
      case StepSizeMethod::Fixed: {
        spaceStep = m_stepSizeFactor;

        break;
      }
      case StepSizeMethod::Dx: {
        spaceStep = m_stepSizeFactor * dx;

        break;
      }
      case StepSizeMethod::Alpha: {
        spaceStep = m_stepSizeFactor * dx;
        alphaStep = m_stepSizeFactor;

        break;
      }
      default: {
        MayDay::Error("StreamerInceptionStepper::inceptionIntegrateTrapezoidal - logic bust");

        break;
      }
      }

      for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
        List<P>& solverParticles    = amrParticles[lvl][dit()].listItems();
        List<P>& processedParticles = amrProcessedParticles[lvl][dit()].listItems();

        for (ListIterator<P> lit(solverParticles); lit.ok();) {
          P& p = lit();

          const RealVect x        = p.position();
          const RealVect vel      = p.velocity();
          const Real     v        = vel.vectorLength();
          const Real     E        = v;
          const Real     alpha    = m_alpha(E, x);
          const Real     eta      = m_eta(E, x);
          const Real     alphaEff = alpha - eta;

          const Real     deltaX = std::min(spaceStep, alphaStep / alphaEff);
          const Real     dt     = deltaX / v;
          const RealVect newPos = p.position() + dt * vel;
          const Real     delta  = (newPos - x).vectorLength();

          const bool outsideDomain = this->particleOutsideGrid(newPos, probLo, probHi);
          const bool insideEB      = this->particleInsideEB(newPos);

          // If particle hit the EB or domain we finish off with a partial Euler step
          Real s = 0.0;

          if (alphaEff < 0.0) {
            processedParticles.transfer(lit);
          }
          else if (insideEB) {
            const RefCountedPtr<BaseIF>& impFunc = m_amr->getBaseImplicitFunction(m_phase);

            if (ParticleOps::ebIntersectionBisect(impFunc, x, newPos, spaceStep, s)) {
              p.weight() = p.weight() + s * delta * alphaEff;
            }

            processedParticles.transfer(lit);
          }
          else if (outsideDomain) {
            if (ParticleOps::domainIntersection(x, newPos, m_amr->getProbLo(), m_amr->getProbHi(), s)) {
              p.weight() = p.weight() + s * delta * alphaEff;
            }

            processedParticles.transfer(lit);
          }
          else {
            // Do an Euler step, storaging alpha(p^k), v(p^k), and the time step size.
            p.template real<0>() = alphaEff;
            p.template real<1>() = dt;
            p.template vect<1>() = vel;

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
      const Real dx = m_amr->getDx()[lvl];

      Real spaceStep = std::numeric_limits<Real>::max();
      Real alphaStep = std::numeric_limits<Real>::max();

      switch (m_stepSizeMethod) {
      case StepSizeMethod::Fixed: {
        spaceStep = m_stepSizeFactor;

        break;
      }
      case StepSizeMethod::Dx: {
        spaceStep = m_stepSizeFactor * dx;

        break;
      }
      case StepSizeMethod::Alpha: {
        spaceStep = m_stepSizeFactor * dx;
        alphaStep = m_stepSizeFactor;

        break;
      }
      default: {
        MayDay::Error("StreamerInceptionStepper::inceptionIntegrateEuler - logic bust");

        break;
      }
      }

      for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
        List<P>& solverParticles    = amrParticles[lvl][dit()].listItems();
        List<P>& processedParticles = amrProcessedParticles[lvl][dit()].listItems();

        for (ListIterator<P> lit(solverParticles); lit.ok();) {

          P& p = lit();

          const Real     dt  = p.template real<1>();
          const RealVect vk  = p.template vect<1>();
          const RealVect vk1 = p.velocity();
          const RealVect x   = p.position();
          const Real     E   = vk1.vectorLength();

          // Note the weird subtraction since p.position() was updated
          // to p^k + dt * v^k.
          const RealVect oldPos = p.position() - dt * vk;
          const RealVect newPos = oldPos + 0.5 * dt * (vk + vk1);

          // Compute new alpha.
          const Real alphak  = p.template real<0>();
          const Real alphak1 = m_alpha(E, x) - m_eta(E, x);
          const Real delta   = (newPos - oldPos).vectorLength();

          // Stop integration for particles that move into regions alpha < 0.0,
          // inside the EB or outside of the domain.
          const bool negativeAlpha = (alphak + alphak1) < 0.0;
          const bool outsideDomain = this->particleOutsideGrid(newPos, probLo, probHi);
          const bool insideEB      = this->particleInsideEB(newPos);

          // If the particle wound up inside the EB we finish off the integration with a partial Euler step
          Real s = 0.0;

          if (negativeAlpha) {
            processedParticles.transfer(lit);
          }
          else if (insideEB) {
            const RefCountedPtr<BaseIF>& impFunc = m_amr->getBaseImplicitFunction(m_phase);

            if (ParticleOps::ebIntersectionBisect(impFunc, oldPos, newPos, spaceStep, s)) {
              p.weight() = p.weight() + s * delta * alphak;
            }

            processedParticles.transfer(lit);
          }
          else if (outsideDomain) {
            if (ParticleOps::domainIntersection(oldPos, newPos, m_amr->getProbLo(), m_amr->getProbHi(), s)) {
              p.weight() = p.weight() + s * delta * alphak;
            }

            processedParticles.transfer(lit);
          }
          else {
            p.position() = newPos;
            p.weight()   = p.weight() + 0.5 * delta * (alphak + alphak1);

            ++lit;
          }
        }

        // Transfer particles that completed their integration (if we're doing partial integration)
        if (!m_fullIntegration) {
          for (ListIterator<P> lit(solverParticles); lit.ok();) {
            P& p = lit();

            if (p.weight() >= m_inceptionK) {
              processedParticles.transfer(lit);
            }
            else {
              ++lit;
            }
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

  // Truncate the weights if we didn't run full integration
  if (!m_fullIntegration) {
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
        List<P>& solverParticles = amrParticles[lvl][dit()].listItems();

        for (ListIterator<P> lit(solverParticles); lit.ok(); ++lit) {
          lit().weight() = std::min(m_inceptionK, lit().weight());
        }
      }
    }
  }

  this->rewindTracerParticles();
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::computeTownsendCriterionStationary() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeTownsendCriterionStationary");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeTownsendCriterionStationary" << endl;
  }

  // Allocate storage for gamma-coefficient
  EBAMRCellData gamma;
  EBAMRCellData expK;

  m_amr->allocate(gamma, m_realm, m_phase, 1);
  m_amr->allocate(expK, m_realm, m_phase, 1);

  // Do all voltages and both polarities
  for (int i = 0; i < m_voltageSweeps.size(); ++i) {

    // Polarities. Here, we are dealing with positive ions so the particles move with the field.
    std::vector<Real> polarities{1.0, -1.0};

    this->seedIonizationParticles(m_voltageSweeps[i]);

    for (const auto& p : polarities) {

      DataOps::setValue(gamma, 0.0);

      this->resetTracerParticles();

      switch (m_inceptionAlgorithm) {
      case IntegrationAlgorithm::Euler: {
        this->townsendTrackEuler(p * m_voltageSweeps[i]);

        break;
      }
      case IntegrationAlgorithm::Trapezoidal: {
        this->townsendTrackTrapezoidal(p * m_voltageSweeps[i]);

        break;
      }
      default: {
        MayDay::Error("StreamerInceptionStepper::computeTownsendCriterionStationary -- logic bust");

        break;
      }
      }

      // Deposit the particles on the mesh and see if they struck a cathode surface.
      m_tracerParticleSolver->deposit(gamma);

      m_amr->conservativeAverage(gamma, m_realm, m_phase);
      m_amr->interpGhost(gamma, m_realm, m_phase);

      // For turning K into exp(K)-1
      auto exponentiate = [](const Real x) -> Real {
        return x > 0.0 ? exp(x) - 1 : 0.0;
      };

      if (p > 0.0) {
        EBAMRCellData Kplus = m_amr->slice(m_inceptionIntegralPlus, Interval(i, i));
        DataOps::copy(expK, Kplus);
        DataOps::compute(expK, exponentiate);
        DataOps::multiply(gamma, expK);

        // If we're running without full integration we truncate the Townsend value to 1.
        if (!m_fullIntegration) {
          DataOps::compute(gamma, [](const Real x) {
            return std::min(x, 1.0);
          });
        }

        DataOps::copy(m_townsendCriterionPlus, gamma, Interval(i, i), Interval(0, 0));

        Real minT = 0.0;
        Real maxT = 0.0;

        DataOps::getMaxMin(maxT, minT, gamma, 0);

        m_maxTPlus[i] = maxT;
      }
      else {
        EBAMRCellData Kminu = m_amr->slice(m_inceptionIntegralMinu, Interval(i, i));
        DataOps::copy(expK, Kminu);
        DataOps::compute(expK, exponentiate);
        DataOps::multiply(gamma, expK);

        // If we're running without full integration we truncate the Townsend value to 1.
        if (!m_fullIntegration) {
          DataOps::compute(gamma, [](const Real x) {
            return std::min(x, 1.0);
          });
        }

        DataOps::copy(m_townsendCriterionMinu, gamma, Interval(i, i), Interval(0, 0));

        Real minT = 0.0;
        Real maxT = 0.0;

        DataOps::getMaxMin(maxT, minT, gamma, 0);

        m_maxTMinu[i] = maxT;
      }
    }
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::computeTownsendCriterionTransient(const Real& a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::computeTownsendCriterionTransient");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeTownsendCriterionTransient" << endl;
  }

  // Compute T-value on each particle using specified algorithm.
  switch (m_inceptionAlgorithm) {
  case IntegrationAlgorithm::Euler: {
    this->townsendTrackEuler(a_voltage);

    break;
  }
  case IntegrationAlgorithm::Trapezoidal: {
    this->townsendTrackTrapezoidal(a_voltage);

    break;
  }
  default: {
    MayDay::Error("StreamerInceptionStepper::computeTownsendCriterionTransient - logic bust");

    break;
  }
  }

  // Compute gamma
  EBAMRCellData gamma;
  m_amr->allocate(gamma, m_realm, m_phase, 1);
  DataOps::setValue(gamma, 0.0);
  m_tracerParticleSolver->deposit(gamma);

  // Turn K into exp(K)-1
  auto exponentiate = [](const Real x) -> Real {
    return x > 0.0 ? exp(x) - 1 : 0.0;
  };
  DataOps::copy(m_townsendCriterion, m_inceptionIntegral);
  DataOps::compute(m_townsendCriterion, exponentiate);
  DataOps::multiply(m_townsendCriterion, gamma);

  // If we're running without full integration we truncate the Townsend value to 1.
  if (!m_fullIntegration) {
    DataOps::compute(m_townsendCriterion, [](const Real x) {
      return std::min(x, 1.0);
    });
  }

  m_amr->conservativeAverage(m_townsendCriterion, m_realm, m_phase);
  m_amr->interpGhost(m_townsendCriterion, m_realm, m_phase);
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::townsendTrackEuler(const Real& a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::townsendTrackEuler");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::townsendTrackEuler" << endl;
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

  size_t outcastsBefore  = 0;
  size_t particlesBefore = 0;

  if (m_debug) {
    outcastsBefore  = amrParticles.getNumberOfOutcastParticlesGlobal();
    particlesBefore = amrParticles.getNumberOfValidParticlesGlobal();
  }

  // Allocate something that holds the velocity of the ions
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, SpaceDim);
  this->superposition(scratch, m_electricFieldInho, m_electricFieldHomo, a_voltage);

  m_tracerParticleSolver->setVelocity(scratch);
  m_tracerParticleSolver->interpolateVelocities();

  while (amrParticles.getNumberOfValidParticlesGlobal() > 0) {

    // Euler integration.
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const Real dx = m_amr->getDx()[lvl];

      // Set the spatial step.
      Real spaceStep = std::numeric_limits<Real>::max();

      switch (m_stepSizeMethod) {
      case StepSizeMethod::Fixed: {
        spaceStep = m_stepSizeFactor;

        break;
      }
      case StepSizeMethod::Dx: {
        spaceStep = m_stepSizeFactor * dx;

        break;
      }
      case StepSizeMethod::Alpha: {
        spaceStep = m_stepSizeFactor * dx;

        break;
      }
      default: {
        MayDay::Error("StreamerInceptionStepper::townsendTrackEuler - logic bust");

        break;
      }
      }

      // Integrate particles until they leave alpha > 0 or strike a cathode surface.
      for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
        List<P>& solverParticles    = amrParticles[lvl][dit()].listItems();
        List<P>& processedParticles = amrProcessedParticles[lvl][dit()].listItems();

        for (ListIterator<P> lit(solverParticles); lit.ok();) {
          P& p = lit();

          const RealVect x      = p.position();
          const RealVect vel    = p.velocity();
          const Real     v      = vel.vectorLength();
          const Real     E      = v;
          const Real     deltaX = spaceStep;
          const Real     dt     = deltaX / v;
          const RealVect newPos = p.position() + dt * vel;

          const bool outsideDomain = this->particleOutsideGrid(newPos, probLo, probHi);
          const bool insideEB      = this->particleInsideEB(newPos);

          // Stop integration if avalanche phase is over
          const Real alpha         = m_alpha(E, x);
          const Real eta           = m_eta(E, x);
          const Real alphaEff      = alpha - eta;
          const bool negativeAlpha = alphaEff <= 0.0;

          p.weight() = 0.0;

          if (insideEB) {
            p.weight() = m_secondaryEmission(E, x);

            processedParticles.transfer(lit);
          }
          else if (outsideDomain || negativeAlpha) {
            processedParticles.transfer(lit);
          }
          else {
            p.position() = newPos;

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

  size_t outcastsAfter  = 0;
  size_t particlesAfter = 0;

  if (m_debug) {
    outcastsAfter  = amrParticles.getNumberOfOutcastParticlesGlobal();
    particlesAfter = amrParticles.getNumberOfValidParticlesGlobal();
  }

  if (particlesBefore != particlesAfter) {
    MayDay::Warning("StreamerInceptionStepper::townsendTrackEuler - lost/gained particles!");
  }
  if (outcastsBefore != outcastsAfter) {
    MayDay::Warning("StreamerInceptionStepper::townsendTrackEuler - lost/gained outcast particles!");
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::townsendTrackTrapezoidal(const Real& a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::townsendTrackTrapezoidal");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::townsendTrackTrapezoidal" << endl;
  }

  // TLDR: We move the particle using Heun's method.
  //
  //          p^(k+1) = p^k + 0.5 * dt * [v(p^k) + v(p^l)]
  //
  //       where p^l = p^k + dt * v(p^k). We will set dt = d/|v(p^k)|.

  const RealVect probLo = m_amr->getProbLo();
  const RealVect probHi = m_amr->getProbHi();

  // Allocate a data holder for holding the processed particles. This
  // will be faster because then we only have to iterate through the
  // particles that are actually still moving.
  ParticleContainer<P> amrProcessedParticles;
  m_amr->allocate(amrProcessedParticles, m_realm);

  ParticleContainer<P>& amrParticles = m_tracerParticleSolver->getParticles();

  m_tracerParticleSolver->remap();

  size_t outcastsBefore  = 0;
  size_t particlesBefore = 0;

  if (m_debug) {
    outcastsBefore  = amrParticles.getNumberOfOutcastParticlesGlobal();
    particlesBefore = amrParticles.getNumberOfValidParticlesGlobal();
  }

  // Allocate something that holds the velocity of the ions
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, SpaceDim);
  this->superposition(scratch, m_electricFieldInho, m_electricFieldHomo, a_voltage);

  m_tracerParticleSolver->setVelocity(scratch);
  m_tracerParticleSolver->interpolateVelocities();

  while (amrParticles.getNumberOfValidParticlesGlobal() > 0) {

    // Euler stage.
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const Real dx = m_amr->getDx()[lvl];

      Real spaceStep = std::numeric_limits<Real>::max();
      Real alphaStep = std::numeric_limits<Real>::max();

      switch (m_stepSizeMethod) {
      case StepSizeMethod::Fixed: {
        spaceStep = m_stepSizeFactor;

        break;
      }
      case StepSizeMethod::Dx: {
        spaceStep = m_stepSizeFactor * dx;

        break;
      }
      case StepSizeMethod::Alpha: {
        spaceStep = m_stepSizeFactor * dx;

        break;
      }
      default: {
        MayDay::Error("StreamerInceptionStepper::townsendTrackTrapezoidal - logic bust");

        break;
      }
      }

      // First stage
      for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
        List<P>& solverParticles    = amrParticles[lvl][dit()].listItems();
        List<P>& processedParticles = amrProcessedParticles[lvl][dit()].listItems();

        for (ListIterator<P> lit(solverParticles); lit.ok();) {
          P& p = lit();

          const RealVect x        = p.position();
          const RealVect vel      = p.velocity();
          const Real     v        = vel.vectorLength();
          const Real     E        = v;
          const Real     alpha    = m_alpha(E, x);
          const Real     eta      = m_eta(E, x);
          const Real     alphaEff = alpha - eta;

          // Compute a time step that we will use for the integration.
          const Real     deltaX = spaceStep;
          const Real     dt     = deltaX / v;
          const RealVect newPos = p.position() + vel * dt;

          const bool outsideDomain = this->particleOutsideGrid(newPos, probLo, probHi);
          const bool insideEB      = this->particleInsideEB(newPos);

          if (insideEB) {
            p.weight() = m_secondaryEmission(E, x);

            processedParticles.transfer(lit);
          }
          else if (alphaEff < 0.0) {
            processedParticles.transfer(lit);
          }
          else if (outsideDomain) {
            processedParticles.transfer(lit);
          }
          else {
            // Do an Euler step, storaging alpha(p^k), v(p^k), and the time step size.
            p.template real<0>() = alphaEff;
            p.template real<1>() = dt;
            p.template vect<1>() = vel;

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
          const RealVect x   = p.position();
          const Real     E   = vk1.vectorLength();

          // Note the weird subtraction since p.position() was updated
          // to p^k + dt * v^k.
          const RealVect oldPos = p.position() - dt * vk;
          const RealVect newPos = p.position() + 0.5 * dt * (vk1 - vk);

          // Compute new alpha.
          const Real alphak  = p.template real<0>();
          const Real alphak1 = m_alpha(E, x) - m_eta(E, x);

          // Stop integration for particles that move into regions alpha < 0.0,
          // inside the EB or outside of the domain.
          const bool negativeAlpha = (alphak + alphak1) < 0.0;
          const bool outsideDomain = this->particleOutsideGrid(newPos, probLo, probHi);
          const bool insideEB      = this->particleInsideEB(newPos);

          // If the particle wound up inside the EB we finish off the integration with a partial Euler step
          Real s = 0.0;

          if (insideEB) {
            p.weight() = m_secondaryEmission(E, oldPos);

            processedParticles.transfer(lit);
          }
          else if (outsideDomain) {
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
  }

  // Copy processed particles over to the solver particles.
  ParticleOps::copyDestructive(amrParticles, amrProcessedParticles);

  this->rewindTracerParticles();

  size_t outcastsAfter  = 0;
  size_t particlesAfter = 0;

  if (m_debug) {
    outcastsAfter  = amrParticles.getNumberOfOutcastParticlesGlobal();
    particlesAfter = amrParticles.getNumberOfValidParticlesGlobal();
  }

  if (particlesBefore != particlesAfter) {
    MayDay::Warning("StreamerInceptionStepper::townsendTrackTrapezoidal - lost/gained particles!");
  }
  if (outcastsBefore != outcastsAfter) {
    MayDay::Warning("StreamerInceptionStepper::townsendTrackTrapezoidal - lost/gained outcast particles!");
  }
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

  // Compute the electric field at the input voltage
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, SpaceDim);
  this->superposition(scratch, a_voltage);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real               dx    = m_amr->getDx()[lvl];

    const Real vol  = std::pow(dx, SpaceDim);
    const Real area = std::pow(dx, SpaceDim - 1);

    for (DataIterator dit(dbl); dit.ok(); ++dit) {

      const EBISBox&       ebisbox    = ebisl[dit()];
      const BaseFab<bool>& validCells = (*m_amr->getValidCells(m_realm)[lvl])[dit()];

      const EBCellFAB& electricField     = (*scratch[lvl])[dit()];
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
            const Real E = EE.vectorLength();

            const Real alpha = m_alpha(E, pos);
            const Real eta   = m_eta(E, pos);
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
            const Real     E   = EE.vectorLength();

            const Real kappa    = ebisbox.volFrac(vof);
            const Real areaFrac = ebisbox.bndryArea(vof);
            const Real alpha    = m_alpha(E, pos);
            const Real eta      = m_eta(E, pos);
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

  // Holds the electric field
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, SpaceDim);

  for (int i = 0; i < m_voltageSweeps.size(); i++) {

    // Compute the electric field into the scratch storage
    this->superposition(scratch, m_voltageSweeps[i]);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const Real               dx    = m_amr->getDx()[lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        const EBISBox& ebisbox = ebisl[dit()];

        EBCellFAB& bgIonization    = (*m_backgroundIonization[lvl])[dit()];
        FArrayBox& bgIonizationReg = bgIonization.getFArrayBox();

        const EBCellFAB& electricField    = (*scratch[lvl])[dit()];
        const FArrayBox& electricFieldReg = electricField.getFArrayBox();

        auto regularKernel = [&](const IntVect& iv) -> void {
          const RealVect pos = probLo + (0.5 * RealVect::Unit + RealVect(iv)) * dx;
          const RealVect EE  = RealVect(
            D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));
          const Real E = EE.vectorLength();

          bgIonizationReg(iv, i) = m_backgroundRate(E, pos);
        };

        auto irregularKernel = [&](const VolIndex& vof) -> void {
          const RealVect pos = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);
          const RealVect EE  = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));
          const Real     E   = EE.vectorLength();

          bgIonization(vof, i) = m_backgroundRate(E, pos);
        };

        // Execute kernels
        VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

        BoxLoops::loop(dbl[dit()], regularKernel);
        BoxLoops::loop(vofit, irregularKernel);
      }
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

  // Holds the electric field
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, SpaceDim);

  for (int i = 0; i < m_voltageSweeps.size(); i++) {

    // Compute the electric field into the scratch storage
    this->superposition(scratch, m_voltageSweeps[i]);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
      const Real               dx    = m_amr->getDx()[lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        const EBISBox& ebisbox = ebisl[dit()];

        EBCellFAB& detachment    = (*m_detachment[lvl])[dit()];
        FArrayBox& detachmentReg = detachment.getFArrayBox();

        const EBCellFAB& ionDensity    = (*m_ionSolver->getPhi()[lvl])[dit()];
        const EBCellFAB& electricField = (*scratch[lvl])[dit()];

        const FArrayBox& ionDensityReg    = ionDensity.getFArrayBox();
        const FArrayBox& electricFieldReg = electricField.getFArrayBox();

        auto regularKernel = [&](const IntVect& iv) -> void {
          const RealVect pos = probLo + (0.5 * RealVect::Unit + RealVect(iv)) * dx;
          const RealVect EE  = RealVect(
            D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));
          const Real E   = EE.vectorLength();
          const Real phi = ionDensityReg(iv, 0);

          detachmentReg(iv, i) = m_detachmentRate(E, pos) * phi;
        };

        auto irregularKernel = [&](const VolIndex& vof) -> void {
          const RealVect pos = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);
          const RealVect EE  = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));
          const Real     E   = EE.vectorLength();
          const Real     phi = ionDensity(vof, 0);

          detachment(vof, i) = m_detachmentRate(E, pos) * phi;
        };

        // Execute kernels
        VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

        BoxLoops::loop(dbl[dit()], regularKernel);
        BoxLoops::loop(vofit, irregularKernel);
      }
    }

    m_amr->conservativeAverage(m_detachment, m_realm, m_phase);
  }
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

  // Holds the electric field
  EBAMRCellData scratchPlus;
  EBAMRCellData scratchMinu;

  m_amr->allocate(scratchPlus, m_realm, m_phase, SpaceDim);
  m_amr->allocate(scratchMinu, m_realm, m_phase, SpaceDim);

  for (int i = 0; i < m_voltageSweeps.size(); i++) {

    // Compute the electric field into the scratch storage
    this->superposition(scratchPlus, m_voltageSweeps[i]);
    this->superposition(scratchMinu, -m_voltageSweeps[i]);

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

        const EBCellFAB& electricFieldPlus = (*scratchPlus[lvl])[dit()];
        const EBCellFAB& electricFieldMinu = (*scratchMinu[lvl])[dit()];

        auto irregularKernel = [&](const VolIndex& vof) -> void {
          const RealVect pos = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);

          const RealVect Eplus = RealVect(
            D_DECL(electricFieldPlus(vof, 0), electricFieldPlus(vof, 1), electricFieldPlus(vof, 2)));
          const RealVect Eminu = RealVect(
            D_DECL(electricFieldMinu(vof, 0), electricFieldMinu(vof, 1), electricFieldMinu(vof, 2)));

          const Real normalEplus = Eplus.dotProduct(ebisbox.normal(vof));
          const Real normalEminu = Eminu.dotProduct(ebisbox.normal(vof));

          emissionPlus(vof, i) = 0.0;
          emissionMinu(vof, i) = 0.0;

          if (normalEplus < 0.0) {
            emissionPlus(vof, i) = m_fieldEmission(Eplus.vectorLength(), pos);
          }
          if (normalEminu < 0.0) {
            emissionMinu(vof, i) = m_fieldEmission(Eminu.vectorLength(), pos);
          }
        };

        // Execute kernels
        VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
        BoxLoops::loop(vofit, irregularKernel);
      }
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

  // Compute the electric field into the scratch storage
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, SpaceDim);
  this->superposition(scratch, a_voltage);

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

      const EBCellFAB& electricField    = (*scratch[lvl])[dit()];
      const FArrayBox& electricFieldReg = electricField.getFArrayBox();

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const RealVect pos = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);
        const RealVect EE  = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));
        const Real     E   = EE.vectorLength();

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
StreamerInceptionStepper<P, F, C>::evaluateFunction(
  EBAMRCellData&                                             a_data,
  const Real&                                                a_voltage,
  const std::function<Real(const Real E, const RealVect x)>& a_func) const noexcept
{
  CH_TIME("StreamerInceptionStepper::evaluateFunction");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::evaluateFunction" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    this->evaluateFunction(*a_data[lvl], a_voltage, a_func, lvl);
  }
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::evaluateFunction(LevelData<EBCellFAB>& a_data,
                                                    const Real&           a_voltage,
                                                    const std::function<Real(const Real E, const RealVect x)>& a_func,
                                                    const int a_level) const noexcept
{
  CH_TIME("StreamerInceptionStepper::evaluateFunction(level)");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::evaluateFunction(level)" << endl;
  }

  CH_assert(a_level >= 0);
  CH_assert(a_level <= m_amr->getFinestLevel());

  // Compute the electric field
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, SpaceDim);
  this->superposition(scratch, a_voltage);  

  const RealVect probLo = m_amr->getProbLo();

  const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[a_level];
  const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[a_level];
  const Real               dx    = m_amr->getDx()[a_level];

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const EBISBox& ebisbox = ebisl[dit()];

    // Here, emissionPlus is when we have the standard voltage (V = 1) and
    // emissionMinu is when the voltage is reverted.
    EBCellFAB& data    = a_data[dit()];
    FArrayBox& dataReg = data.getFArrayBox();

    data.setVal(0.0);

    const EBCellFAB& electricField    = (*scratch[a_level])[dit()];
    const FArrayBox& electricFieldReg = electricField.getFArrayBox();

    auto regularKernel = [&](const IntVect& iv) -> void {
      const RealVect pos = probLo + (0.5 * RealVect::Unit + RealVect(iv)) * dx;
      const RealVect EE  = RealVect(D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));
      const Real     E   = EE.vectorLength();

      dataReg(iv, 0) = a_func(E, pos);
    };

    auto irregularKernel = [&](const VolIndex& vof) -> void {
      const RealVect pos = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);
      const RealVect EE  = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));
      const Real     E   = EE.vectorLength();

      data(vof, 0) = a_func(E, pos);
    };

    // Execute kernels
    Box          cellBox = dbl[dit()];
    VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[a_level])[dit()];

    BoxLoops::loop(cellBox, regularKernel);
    BoxLoops::loop(vofit, irregularKernel);
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

  CH_assert(m_mode == Mode::Stationary);

  if (m_voltageSweeps.size() < 2) {
    DataOps::setValue(m_inceptionVoltagePlus, std::numeric_limits<Real>::quiet_NaN());
    DataOps::setValue(m_inceptionVoltageMinu, std::numeric_limits<Real>::quiet_NaN());

    MayDay::Warning("StreamerInceptionStepper::computeInceptionVoltageVolume -- not enough voltages for estimating "
                    "inception voltage");
  }
  else {

    constexpr int comp = 0;

    // Function which interpolates the inception voltage if possible. Used in the kernels. This one does interpolation.
    auto calcUincInterp = [Kinc = this->m_inceptionK,
                           &V   = this->m_voltageSweeps](const std::vector<Real>& K,
                                                       const std::vector<Real>& T) -> std::array<Real, 3> {
      Real streamerInc = std::numeric_limits<Real>::quiet_NaN();
      Real townsendInc = std::numeric_limits<Real>::quiet_NaN();

      bool foundInception = false;

      // Streamer criterion
      for (size_t i = 0; i < K.size() - 1; i++) {
        if (K[i] <= Kinc && K[i + 1] > Kinc) {
          streamerInc = V[i] + (Kinc - K[i]) * (V[i + 1] - V[i]) / (K[i + 1] - K[i]);

          break;
        }
        else if (K[i] == Kinc) {
          streamerInc = V[i];

          break;
        }
      }

      // Townsend criterion
      for (size_t i = 0; i < T.size() - 1; i++) {
        if (T[i] <= 1.0 && T[i + 1] > 1) {
          townsendInc = V[i] + (1.0 - T[i]) * (V[i + 1] - V[i]) / (T[i + 1] - T[i]);

          break;
        }
        else if (T[i] == 1.0) {
          townsendInc = V[i];

          break;
        }
      }

      Real Uinc;

      if (std::isnan(streamerInc) && std::isnan(townsendInc)) {
        Uinc = std::numeric_limits<Real>::quiet_NaN();
      }
      else if (std::isnan(streamerInc) && !std::isnan(townsendInc)) {
        Uinc = townsendInc;
      }
      else if (!std::isnan(streamerInc) && std::isnan(townsendInc)) {
        Uinc = streamerInc;
      }
      else {
        Uinc = std::min(streamerInc, townsendInc);
      }

      return std::array<Real, 3>{Uinc, streamerInc, townsendInc};
    };

    // Same functionality as above, but without interpolation.
    auto calcUincNoInterp = [Kinc = this->m_inceptionK,
                             &V   = this->m_voltageSweeps](const std::vector<Real>& K,
                                                         const std::vector<Real>& T) -> std::array<Real, 3> {
      Real streamerInc = std::numeric_limits<Real>::quiet_NaN();
      Real townsendInc = std::numeric_limits<Real>::quiet_NaN();

      // Streamer criterion
      for (size_t i = 0; i < K.size() - 1; i++) {
        if (K[i] <= Kinc && K[i + 1] >= Kinc) {
          streamerInc = V[i];

          break;
        }
      }

      // Townsend criterion
      for (size_t i = 0; i < T.size() - 1; i++) {
        if (T[i] <= 1.0 && T[i + 1] >= 1) {
          townsendInc = V[i];

          break;
        }
      }

      Real Uinc;

      if (std::isnan(streamerInc) && std::isnan(townsendInc)) {
        Uinc = std::numeric_limits<Real>::quiet_NaN();
      }
      else if (std::isnan(streamerInc) && !std::isnan(townsendInc)) {
        Uinc = townsendInc;
      }
      else if (!std::isnan(streamerInc) && std::isnan(townsendInc)) {
        Uinc = streamerInc;
      }
      else {
        Uinc = std::min(streamerInc, townsendInc);
      }

      return std::array<Real, 3>{Uinc, streamerInc, townsendInc};
    };

    // Iterate through m_inceptionIntegral data and calculate inception voltage.
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); ++lvl) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        EBCellFAB& inceptionVoltagePlus = (*m_inceptionVoltagePlus[lvl])[dit()];
        EBCellFAB& inceptionVoltageMinu = (*m_inceptionVoltageMinu[lvl])[dit()];

        EBCellFAB& streamerInceptionVoltagePlus = (*m_streamerInceptionVoltagePlus[lvl])[dit()];
        EBCellFAB& streamerInceptionVoltageMinu = (*m_streamerInceptionVoltageMinu[lvl])[dit()];

        EBCellFAB& townsendInceptionVoltagePlus = (*m_townsendInceptionVoltagePlus[lvl])[dit()];
        EBCellFAB& townsendInceptionVoltageMinu = (*m_townsendInceptionVoltageMinu[lvl])[dit()];

        FArrayBox& inceptionVoltagePlusReg = inceptionVoltagePlus.getFArrayBox();
        FArrayBox& inceptionVoltageMinuReg = inceptionVoltageMinu.getFArrayBox();

        FArrayBox& streamerInceptionVoltagePlusReg = streamerInceptionVoltagePlus.getFArrayBox();
        FArrayBox& streamerInceptionVoltageMinuReg = streamerInceptionVoltageMinu.getFArrayBox();

        FArrayBox& townsendInceptionVoltagePlusReg = townsendInceptionVoltagePlus.getFArrayBox();
        FArrayBox& townsendInceptionVoltageMinuReg = townsendInceptionVoltageMinu.getFArrayBox();

        const EBCellFAB& inceptionIntegralPlus = (*m_inceptionIntegralPlus[lvl])[dit()];
        const EBCellFAB& inceptionIntegralMinu = (*m_inceptionIntegralMinu[lvl])[dit()];

        const FArrayBox& inceptionIntegralPlusReg = inceptionIntegralPlus.getFArrayBox();
        const FArrayBox& inceptionIntegralMinuReg = inceptionIntegralMinu.getFArrayBox();

        const EBCellFAB& townsendCriterionPlus = (*m_townsendCriterionPlus[lvl])[dit()];
        const EBCellFAB& townsendCriterionMinu = (*m_townsendCriterionMinu[lvl])[dit()];

        const FArrayBox& townsendCriterionPlusReg = townsendCriterionPlus.getFArrayBox();
        const FArrayBox& townsendCriterionMinuReg = townsendCriterionMinu.getFArrayBox();

        // Regular kernel.
        auto regularKernel = [&](const IntVect& iv) -> void {
          std::vector<Real> Kplus;
          std::vector<Real> Kminu;

          std::vector<Real> Tplus;
          std::vector<Real> Tminu;

          for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
            Kplus.emplace_back(inceptionIntegralPlusReg(iv, i));
            Kminu.emplace_back(inceptionIntegralMinuReg(iv, i));

            Tplus.emplace_back(townsendCriterionPlusReg(iv, i));
            Tminu.emplace_back(townsendCriterionMinuReg(iv, i));
          }

          std::array<Real, 3> inceptionVoltagesPlus;
          std::array<Real, 3> inceptionVoltagesMinu;

          if (m_fullIntegration) {
            inceptionVoltagesPlus = calcUincInterp(Kplus, Tplus);
            inceptionVoltagesMinu = calcUincInterp(Kminu, Tminu);
          }
          else {
            inceptionVoltagesPlus = calcUincNoInterp(Kplus, Tplus);
            inceptionVoltagesMinu = calcUincNoInterp(Kminu, Tminu);
          }

          inceptionVoltagePlusReg(iv, comp) = std::get<0>(inceptionVoltagesPlus);
          inceptionVoltageMinuReg(iv, comp) = std::get<0>(inceptionVoltagesMinu);

          streamerInceptionVoltagePlusReg(iv, comp) = std::get<1>(inceptionVoltagesPlus);
          streamerInceptionVoltageMinuReg(iv, comp) = std::get<1>(inceptionVoltagesMinu);

          townsendInceptionVoltagePlusReg(iv, comp) = std::get<2>(inceptionVoltagesPlus);
          townsendInceptionVoltageMinuReg(iv, comp) = std::get<2>(inceptionVoltagesMinu);
        };

        // Irregular kernel.
        auto irregularKernel = [&](const VolIndex& vof) -> void {
          std::vector<Real> Kplus;
          std::vector<Real> Kminu;

          std::vector<Real> Tplus;
          std::vector<Real> Tminu;

          for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
            Kplus.emplace_back(inceptionIntegralPlus(vof, i));
            Kminu.emplace_back(inceptionIntegralMinu(vof, i));

            Tplus.emplace_back(townsendCriterionPlus(vof, i));
            Tminu.emplace_back(townsendCriterionMinu(vof, i));
          }

          std::array<Real, 3> inceptionVoltagesPlus;
          std::array<Real, 3> inceptionVoltagesMinu;

          if (m_fullIntegration) {
            inceptionVoltagesPlus = calcUincInterp(Kplus, Tplus);
            inceptionVoltagesMinu = calcUincInterp(Kminu, Tminu);
          }
          else {
            inceptionVoltagesPlus = calcUincNoInterp(Kplus, Tplus);
            inceptionVoltagesMinu = calcUincNoInterp(Kminu, Tminu);
          }

          inceptionVoltagePlus(vof, comp) = std::get<0>(inceptionVoltagesPlus);
          inceptionVoltageMinu(vof, comp) = std::get<0>(inceptionVoltagesMinu);

          streamerInceptionVoltagePlus(vof, comp) = std::get<1>(inceptionVoltagesPlus);
          streamerInceptionVoltageMinu(vof, comp) = std::get<1>(inceptionVoltagesMinu);

          townsendInceptionVoltagePlus(vof, comp) = std::get<2>(inceptionVoltagesPlus);
          townsendInceptionVoltageMinu(vof, comp) = std::get<2>(inceptionVoltagesMinu);
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
            if (std::abs(U) < UInc) {
              UInc = std::abs(U);
              xInc = pos;
            }
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
            if (std::abs(U) < UInc) {
              UInc = std::abs(U);
              xInc = pos;
            }
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

        const EBCellFAB& townsendCritPlus = (*m_townsendCriterionPlus[lvl])[dit()];
        const EBCellFAB& townsendCritMinu = (*m_townsendCriterionMinu[lvl])[dit()];

        const FArrayBox& townsendCritPlusReg = townsendCritPlus.getFArrayBox();
        const FArrayBox& townsendCritMinuReg = townsendCritMinu.getFArrayBox();

        auto regularKernel = [&](const IntVect& iv) -> void {
          if (validCells(iv, 0) && ebisbox.isRegular(iv)) {
            if (inceptionIntegralPlusReg(iv, i) >= m_inceptionK || townsendCritPlusReg(iv, i) >= 1.0) {
              criticalVolumePlus += std::pow(dx, SpaceDim);
            }
            if (inceptionIntegralMinuReg(iv, i) >= m_inceptionK || townsendCritMinuReg(iv, i) >= 1.0) {
              criticalVolumeMinu += std::pow(dx, SpaceDim);
            }
          }
        };

        auto irregularKernel = [&](const VolIndex& vof) -> void {
          if (validCells(vof.gridIndex())) {
            const Real kappa = ebisbox.volFrac(vof);

            if (inceptionIntegralPlus(vof, i) >= m_inceptionK || townsendCritPlus(vof, i) >= 1.0) {
              criticalVolumePlus += kappa * std::pow(dx, SpaceDim);
            }
            if (inceptionIntegralMinu(vof, i) >= m_inceptionK || townsendCritMinu(vof, i) >= 1.0) {
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

      const EBCellFAB& townsendCriterion    = (*m_townsendCriterion[lvl])[dit()];
      const FArrayBox& townsendCriterionReg = townsendCriterion.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        if (validCells(iv, 0) && ebisbox.isRegular(iv)) {
          if (inceptionIntegralReg(iv, 0) >= m_inceptionK || townsendCriterionReg(iv, 0) >= 1.0) {
            Vcr += vol;
          }
        }
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        if (validCells(vof.gridIndex())) {
          if (inceptionIntegral(vof, 0) >= m_inceptionK || townsendCriterion(vof, 0) >= 1.0) {
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
StreamerInceptionStepper<P, F, C>::computeCriticalAreaStationary() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeCriticalAreaStationary");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeCriticalAreaStationary" << endl;
  }

  // Solve critical area of K values for each voltage
  for (size_t i = 0; i < m_voltageSweeps.size(); ++i) {
    Real criticalAreaPlus = 0.0;
    Real criticalAreaMinu = 0.0;

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

        const EBCellFAB& townsendCriterionPlus = (*m_townsendCriterionPlus[lvl])[dit()];
        const EBCellFAB& townsendCriterionMinu = (*m_townsendCriterionMinu[lvl])[dit()];

        auto irregularKernel = [&](const VolIndex& vof) -> void {
          if (validCells(vof.gridIndex())) {
            const Real boundaryArea = ebisbox.bndryArea(vof);

            if (inceptionIntegralPlus(vof, i) >= m_inceptionK || townsendCriterionPlus(vof, i) >= 1.0) {
              criticalAreaPlus += boundaryArea * std::pow(dx, SpaceDim - 1);
            }
            if (inceptionIntegralMinu(vof, i) >= m_inceptionK || townsendCriterionMinu(vof, i) >= 1.0) {
              criticalAreaMinu += boundaryArea * std::pow(dx, SpaceDim - 1);
            }
          }
        };

        // Kernel regions.
        VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

        BoxLoops::loop(vofit, irregularKernel);
      }
    }

    m_criticalAreaPlus.push_back(ParallelOps::sum(criticalAreaPlus));
    m_criticalAreaMinu.push_back(ParallelOps::sum(criticalAreaMinu));
  }
}

template <typename P, typename F, typename C>
Real
StreamerInceptionStepper<P, F, C>::computeCriticalAreaTransient() const noexcept
{
  CH_TIME("StreamerInceptionStepper::computeCriticalAreaTransient");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeCriticalAreaTransient" << endl;
  }

  Real Acr = 0.0;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); ++lvl) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

    const LevelData<BaseFab<bool>>& validCellsLD = *m_amr->getValidCells(m_realm)[lvl];

    const Real dx   = m_amr->getDx()[lvl];
    const Real area = std::pow(dx, SpaceDim - 1);

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      const EBISBox&       ebisbox    = ebisl[dit()];
      const BaseFab<bool>& validCells = validCellsLD[dit()];

      const EBCellFAB& inceptionIntegral = (*m_inceptionIntegral[lvl])[dit()];
      const EBCellFAB& townsendCriterion = (*m_townsendCriterion[lvl])[dit()];

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        if (validCells(vof.gridIndex())) {
          if (inceptionIntegral(vof, 0) >= m_inceptionK || townsendCriterion(vof, 0) >= 1.0) {
            Acr += ebisbox.bndryArea(vof) * area;
          }
        }
      };

      // Kernel regions.
      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      BoxLoops::loop(vofit, irregularKernel);
    }
  }

  return ParallelOps::sum(Acr);
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::computeIonizationVolumeStationary() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeIonizationVolumeStationary");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeIonizationVolumeStationary" << endl;
  }

  // Storage something that holds the electric filed
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, SpaceDim);

  // Solve ionization volume for each voltage
  for (size_t i = 0; i < m_voltageSweeps.size(); ++i) {
    Real ionizationVolume = 0.0;

    this->superposition(scratch, m_voltageSweeps[i]);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); ++lvl) {
      const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
      const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

      const LevelData<BaseFab<bool>>& validCellsLD = *m_amr->getValidCells(m_realm)[lvl];

      const Real     dx     = m_amr->getDx()[lvl];
      const RealVect probLo = m_amr->getProbLo();

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
        const EBISBox&       ebisbox    = ebisl[dit()];
        const BaseFab<bool>& validCells = validCellsLD[dit()];

        const EBCellFAB& electricField    = (*scratch[lvl])[dit()];
        const FArrayBox& electricFieldReg = electricField.getFArrayBox();

        auto regularKernel = [&](const IntVect& iv) -> void {
          if (validCells(iv, 0) && ebisbox.isRegular(iv)) {
            const RealVect x  = probLo + dx * (0.5 * RealVect::Unit + RealVect(iv));
            const RealVect EE = RealVect(
              D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));
            const Real E = EE.vectorLength();

            const Real alpha = m_alpha(E, x);
            const Real eta   = m_eta(E, x);

            if (alpha >= eta) {
              ionizationVolume += std::pow(dx, SpaceDim);
            }
          }
        };

        auto irregularKernel = [&](const VolIndex& vof) -> void {
          if (validCells(vof.gridIndex())) {

            const RealVect x  = probLo + Location::position(Location::Cell::Center, vof, ebisbox, dx);
            const RealVect EE = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));
            const Real     E  = EE.vectorLength();

            const Real alpha = m_alpha(E, x);
            const Real eta   = m_eta(E, x);

            const Real kappa = ebisbox.volFrac(vof);

            if (alpha >= eta) {
              ionizationVolume += kappa * std::pow(dx, SpaceDim);
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

    m_ionizationVolume.push_back(ParallelOps::sum(ionizationVolume));
  }
}

template <typename P, typename F, typename C>
Real
StreamerInceptionStepper<P, F, C>::computeIonizationVolumeTransient(const Real& a_voltage) const noexcept
{
  CH_TIME("StreamerInceptionStepper::computeIonizationVolumeTransient");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeIonizationVolumeTransient" << endl;
  }

  Real Vion = 0.0;

  // Calculate electric field
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, SpaceDim);
  this->superposition(scratch, m_electricFieldInho, m_electricFieldHomo, a_voltage);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); ++lvl) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

    const LevelData<BaseFab<bool>>& validCellsLD = *m_amr->getValidCells(m_realm)[lvl];

    const Real     dx     = m_amr->getDx()[lvl];
    const Real     vol    = std::pow(dx, SpaceDim);
    const RealVect probLo = m_amr->getProbLo();

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      const EBISBox&       ebisbox    = ebisl[dit()];
      const BaseFab<bool>& validCells = validCellsLD[dit()];

      const EBCellFAB& electricField    = (*scratch[lvl])[dit()];
      const FArrayBox& electricFieldReg = electricField.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        if (validCells(iv, 0) && ebisbox.isRegular(iv)) {

          const RealVect x  = probLo + dx * (0.5 * RealVect::Unit + RealVect(iv));
          const RealVect EE = RealVect(
            D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));
          const Real E = EE.vectorLength();

          const Real alpha = m_alpha(E, x);
          const Real eta   = m_eta(E, x);
          if (alpha >= eta) {
            Vion += vol;
          }
        }
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        if (validCells(vof.gridIndex())) {

          const RealVect x  = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);
          const RealVect EE = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));
          const Real     E  = EE.vectorLength();

          const Real alpha = m_alpha(E, x);
          const Real eta   = m_eta(E, x);
          if (alpha >= eta) {
            Vion = ebisbox.volFrac(vof) * vol;
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

  return ParallelOps::sum(Vion);
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

  const std::pair<Real, RealVect> streamerUIncPlus = this->computeMinimumInceptionVoltage(
    m_streamerInceptionVoltagePlus);
  const std::pair<Real, RealVect> streamerUIncMinu = this->computeMinimumInceptionVoltage(
    m_streamerInceptionVoltageMinu);

  const std::pair<Real, RealVect> townsendUIncPlus = this->computeMinimumInceptionVoltage(
    m_townsendInceptionVoltagePlus);
  const std::pair<Real, RealVect> townsendUIncMinu = this->computeMinimumInceptionVoltage(
    m_townsendInceptionVoltageMinu);

#ifdef CH_MPI
  if (procID() == 0) {
#endif

    std::ofstream output(m_outputFile, std::ofstream::out);

    const std::string lineBreak = "# " + std::string(178, '=');

    // clang-format off
    output << lineBreak << "\n";
    if (std::isnan(UIncPlus.first) || std::isnan(UIncMinu.first)) {
      output << "# Could not compute inception voltage\n";
    }
    else {
      if(m_fullIntegration){
	output << "# Minimum inception voltage(+)  = " << UIncPlus.first << ",\t x = " << UIncPlus.second << "\n";
	output << "# Minimum inception voltage(-)  = " << UIncMinu.first << ",\t x = " << UIncMinu.second << "\n";
	output << "# \n";
	output << "# Streamer inception voltage(+) = " << streamerUIncPlus.first << ",\t x = " << streamerUIncPlus.second << "\n";
	output << "# Streamer inception voltage(-) = " << streamerUIncMinu.first << ",\t x = " << streamerUIncMinu.second << "\n";
	output << "# \n";
	output << "# Townsend inception voltage(+) = " << townsendUIncPlus.first << ",\t x = " << townsendUIncPlus.second << "\n";
	output << "# Townsend inception voltage(-) = " << townsendUIncMinu.first << ",\t x = " << townsendUIncMinu.second << "\n";		
      }
      else{
	output << "# Minimum inception voltage(+)  >= " << UIncPlus.first << ",\t x = " << UIncPlus.second << "\n";
	output << "# Minimum inception voltage(-)  >= " << UIncMinu.first << ",\t x = " << UIncMinu.second << "\n";
	output << "# \n";
	output << "# Streamer inception voltage(+) >= " << streamerUIncPlus.first << ",\t x = " << streamerUIncPlus.second << "\n";
	output << "# Streamer inception voltage(-) >= " << streamerUIncMinu.first << ",\t x = " << streamerUIncMinu.second << "\n";
	output << "# \n";
	output << "# Townsend inception voltage(+) >= " << townsendUIncPlus.first << ",\t x = " << townsendUIncPlus.second << "\n";
	output << "# Townsend inception voltage(-) >= " << townsendUIncMinu.first << ",\t x = " << townsendUIncMinu.second << "\n";			
      }      
    }


    output << lineBreak << "\n";
    output << left << setw(15) << setfill(' ') << "# +/- Voltage";
    output << left << setw(15) << setfill(' ') << "Max K(+)";
    output << left << setw(15) << setfill(' ') << "Max K(-)";
    output << left << setw(20) << setfill(' ') << "Max T(+)";
    output << left << setw(20) << setfill(' ') << "Max T(-)";    
    output << left << setw(20) << setfill(' ') << "Crit. vol(+)";
    output << left << setw(20) << setfill(' ') << "Crit. vol(-)";
    output << left << setw(20) << setfill(' ') << "Crit. area(+)";
    output << left << setw(20) << setfill(' ') << "Crit. area(-)";
    output << left << setw(20) << setfill(' ') << "Ionization vol." << "\n";
    output << lineBreak << "\n";        

    for (int i = 0; i < m_voltageSweeps.size(); i++) {
      output << left << setw(15) << setfill(' ') << m_voltageSweeps[i];
      output << left << setw(15) << setfill(' ') << m_maxKPlus[i];
      output << left << setw(15) << setfill(' ') << m_maxKMinu[i];
      output << left << setw(20) << setfill(' ') << m_maxTPlus[i];
      output << left << setw(20) << setfill(' ') << m_maxTMinu[i];      
      output << left << setw(20) << setfill(' ') << m_criticalVolumePlus[i];
      output << left << setw(20) << setfill(' ') << m_criticalVolumeMinu[i];
      output << left << setw(20) << setfill(' ') << m_criticalAreaPlus[i];
      output << left << setw(20) << setfill(' ') << m_criticalAreaMinu[i];
      output << left << setw(20) << setfill(' ') << m_ionizationVolume[i];
      output << endl;
    }
    output << lineBreak << "\n";
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
    output << std::left << std::setw(15) << setfill(' ') << "max T(t)";
    output << std::left << std::setw(15) << setfill(' ') << "Vcr(t)";
    output << std::left << std::setw(15) << setfill(' ') << "Acr(t)";
    output << std::left << std::setw(15) << setfill(' ') << "Vion(t)";
    output << std::left << std::setw(15) << setfill(' ') << "lambda(t)";
    output << std::left << std::setw(15) << setfill(' ') << "P(t)";
    output << std::left << std::setw(15) << setfill(' ') << "dP(t, t+dt)";
    output << std::left << std::setw(15) << setfill(' ') << "Time lag";
    output << "\n";

    // Compute dP(t, t+dt)
    std::vector<Real> dProb;
    for (size_t i = 0; i < m_Rdot.size() - 1; i++) {
      const Real t    = m_Rdot[i].first;
      const Real dt   = m_Rdot[i + 1].first - m_Rdot[i].first;
      const Real prob = m_inceptionProbability[i].second;
      const Real Rdot = m_Rdot[i].second;

      dProb.emplace_back((1.0 - prob) * Rdot * dt);
    }
    dProb.emplace_back(0.0);

    // Compute the average waiting time.
    std::vector<Real> tau(m_Rdot.size(), 0.0);
    for (size_t i = 0; i < m_Rdot.size() - 1; i++) {
      const Real t1     = m_Rdot[i].first;
      const Real t2     = m_Rdot[i + 1].first;
      const Real dt     = t2 - t1;
      const Real prob   = m_inceptionProbability[i].second;
      const Real lambda = m_Rdot[i].second;

      tau[i + 1] = tau[i] + t1 * (1 - prob) * lambda * dt;
    }
    for (int i = 0; i < tau.size(); i++) {
      const Real prob = m_inceptionProbability[i].second;

      tau[i] = prob > 0.0 ? tau[i] / prob : std::numeric_limits<Real>::infinity();
    }

    for (size_t i = 0; i < m_Rdot.size(); i++) {
      const Real time = m_Rdot[i].first;

      output << std::left << std::setw(15) << time;
      output << std::left << std::setw(15) << m_voltageCurve(time);
      output << std::left << std::setw(15) << m_maxK[i].second;
      output << std::left << std::setw(15) << m_maxT[i].second;
      output << std::left << std::setw(15) << m_criticalVolume[i].second;
      output << std::left << std::setw(15) << m_criticalArea[i].second;
      output << std::left << std::setw(15) << m_ionizationVolumeTransient[i].second;
      output << std::left << std::setw(15) << m_Rdot[i].second;
      output << std::left << std::setw(15) << m_inceptionProbability[i].second;
      output << std::left << std::setw(15) << dProb[i];
      output << std::left << std::setw(15) << tau[i];
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

  // Compute electric field at the input voltage and set v = -E
  this->superposition(vel, a_voltage);
  DataOps::scale(vel, -1.0);

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

  m_amr->arithmeticAverage(vel, m_realm, m_phase);
  m_amr->interpGhostPwl(vel, m_realm, m_phase);
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

  // Compute the electric field at the input voltage.
  EBAMRCellData scratch;
  m_amr->allocate(scratch, m_realm, m_phase, SpaceDim);
  this->superposition(scratch, a_voltage);

  // Compute the diffusion coefficient on cell centers.
  EBAMRCellData diffCoCell;
  m_amr->allocate(diffCoCell, m_realm, m_phase, 1);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {

      const EBCellFAB& electricField    = (*scratch[lvl])[dit()];
      const FArrayBox& electricFieldReg = electricField.getFArrayBox();

      EBCellFAB& dCo    = (*diffCoCell[lvl])[dit()];
      FArrayBox& dCoReg = dCo.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        const RealVect EE = RealVect(D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));
        const Real     E  = EE.vectorLength();

        dCoReg(iv, 0) = m_ionDiffusion(E);
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const RealVect EE = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));
        const Real     E  = EE.vectorLength();

        dCo(vof, 0) = m_ionDiffusion(E);
      };

      VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];
      Box          cellBox = dbl[dit()];

      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }

  // Coarsen and update ghost cells
  m_amr->arithmeticAverage(diffCoCell, m_realm, m_phase);
  m_amr->interpGhostPwl(diffCoCell, m_realm, m_phase);

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

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::superposition(EBAMRCellData&       a_sumField,
                                                 const EBAMRCellData& a_inhomogeneousField,
                                                 const EBAMRCellData& a_homogeneousField,
                                                 const Real           a_voltage) const noexcept
{
  CH_TIME("StreamerInceptionStepper::superposition(full)");

  DataOps::setValue(a_sumField, 0.0);
  DataOps::incr(a_sumField, a_homogeneousField, a_voltage);
  DataOps::incr(a_sumField, a_inhomogeneousField, 1.0);

  m_amr->arithmeticAverage(a_sumField, m_realm, m_phase);
  m_amr->interpGhostPwl(a_sumField, m_realm, m_phase);
  m_amr->interpToCentroids(a_sumField, m_realm, m_phase);
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::superposition(EBAMRCellData& a_sumField, const Real a_voltage) const noexcept
{
  CH_TIME("StreamerInceptionStepper::superposition(short)");

  this->superposition(a_sumField, m_electricFieldInho, m_electricFieldHomo, a_voltage);
}

template <typename P, typename F, typename C>
void
StreamerInceptionStepper<P, F, C>::writeData(LevelData<EBCellFAB>& a_output,
                                             int&                  a_comp,
                                             const EBAMRCellData&  a_data,
                                             const std::string     a_outputRealm,
                                             const int             a_level,
                                             const bool            a_interpToCentroids,
                                             const bool            a_interpGhost) const noexcept
{

  CH_TIMERS("StreamerInceptionStepper::writeData");
  CH_TIMER("StreamerInceptionStepper::writeData::allocate", t1);
  CH_TIMER("StreamerInceptionStepper::writeData::local_copy", t2);
  CH_TIMER("StreamerInceptionStepper::writeData::interp_ghost", t3);
  CH_TIMER("StreamerInceptionStepper::writeData::interp_centroid", t4);
  CH_TIMER("StreamerInceptionStepper::writeData::final_copy", t5);
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::writeData" << endl;
  }

  // Number of components we are working with.
  const int numComp = a_data[a_level]->nComp();

  // Component ranges that we copy to/from.
  const Interval srcInterv(0, numComp - 1);
  const Interval dstInterv(a_comp, a_comp + numComp - 1);

  CH_START(t1);
  LevelData<EBCellFAB> scratch;
  m_amr->allocate(scratch, m_realm, m_phase, a_level, numComp);
  CH_STOP(t1);

  CH_START(t2);
  m_amr->copyData(scratch, *a_data[a_level], a_level, m_realm, m_realm);
  CH_START(t2);

  // Interpolate ghost cells
  CH_START(t3);
  if (a_level > 0 && a_interpGhost) {
    m_amr->interpGhost(scratch, *a_data[a_level - 1], a_level, m_realm, m_phase);
  }
  CH_STOP(t3);

  CH_START(t4);
  if (a_interpToCentroids) {
    m_amr->interpToCentroids(scratch, m_realm, m_phase, a_level);
  }
  CH_STOP(t4);

  DataOps::setCoveredValue(scratch, 0.0);

  CH_START(t5);
  m_amr->copyData(a_output, scratch, a_level, a_outputRealm, m_realm, dstInterv, srcInterv);
  CH_STOP(t5);

  a_comp += numComp;
}

template <typename P, typename F, typename C>
const EBAMRCellData*
StreamerInceptionStepper<P, F, C>::getElectricField() const noexcept {
  CH_TIMERS("StreamerInceptionStepper::getElectricField");

  return &m_electricFieldHomo;
}

#include <CD_NamespaceFooter.H>

#endif
