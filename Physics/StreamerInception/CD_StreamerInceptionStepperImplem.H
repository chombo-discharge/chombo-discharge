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

// Std includes
#include <iostream>
#include <fstream>

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_Units.H>
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

  // Default settings.
  m_verbosity = -1;
  m_realm     = Realm::Primal;
  m_phase     = phase::gas;
  m_mode      = Mode::Static;

  this->parseOptions();

  m_rho = [](const RealVect x) -> Real {
    return 0.0;
  };

  m_alpha = [](const Real E) -> Real {
    return 0.0; };

  m_eta = [](const Real E) -> Real { return 0.0; };

  m_voltageCurve = [](const Real a_time) -> Real { return 1.0; };

  m_backgroundRate = [](const Real E) -> Real { return 0.0; };

  m_fieldEmission = [](const Real E) -> Real { return 0.0; };

  m_negativeIonDensity = [](const RealVect x) -> Real { return 0.0;};

  m_negativeIonMobility = [](const Real E) -> Real { return 0.0; };

#if 1 // To-do list.
  std::cout << "StreamerInceptionStepperImplem.H -- add alpha-tagging to inception tagger\n";
  std::cout << "StreamerInceptionStepperImplem.H -- fix curvature tagger\n";
#endif  
}

template <typename P, typename F>
StreamerInceptionStepper<P, F>::~StreamerInceptionStepper()
{
  CH_TIME("StreamerInceptionStepper::~StreamerInceptionStepper");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::~StreamerInceptionStepper" << endl;
  }

  if (m_writeReport) {
    switch (m_mode) {
    case Mode::Static: {
      this->writeReportStatic();

      break;
    }
    case Mode::Dynamic: {
      this->writeReportDynamic();

      break;
    }
    default: {
      break;
    }
    }
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

  // Always solve for the field using a voltage of one (and then scale up/down later on)
  auto voltage = [](const Real a_time) -> Real { return 1.0; };

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

  // This one is needed for the dynamic runs which also track the ions.
  if(m_mode == Mode::Dynamic) {
    m_negativeIonSolver = RefCountedPtr<TracerParticleSolver<P>>(new TracerParticleSolver<P>());
    m_negativeIonSolver->parseOptions();
    m_negativeIonSolver->setAmr(m_amr);
    m_negativeIonSolver->setComputationalGeometry(m_computationalGeometry);
    m_negativeIonSolver->setRealm(m_realm);
    m_negativeIonSolver->setPhase(m_phase);
    m_negativeIonSolver->setTime(0, 0.0, 0.0);
    m_negativeIonSolver->setName("Negative ion solver");

    m_negativeIonSolver->setVolumeScale(false);
  }
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
  m_tracerParticleSolver->allocate();

  // Needed for all modes.
  m_amr->allocate(m_electricFieldPlus, m_realm, m_phase, SpaceDim);
  m_amr->allocate(m_electricFieldMinu, m_realm, m_phase, SpaceDim);

  switch (m_mode) {
  case Mode::Static: {
    m_amr->allocate(m_inceptionIntegralPlus, m_realm, m_phase, m_voltageSweeps.size());
    m_amr->allocate(m_inceptionIntegralMinu, m_realm, m_phase, m_voltageSweeps.size());
    m_amr->allocate(m_backgroundIonization, m_realm, m_phase, m_voltageSweeps.size());
    m_amr->allocate(m_emissionRatesPlus, m_realm, m_phase, m_voltageSweeps.size());
    m_amr->allocate(m_emissionRatesMinu, m_realm, m_phase, m_voltageSweeps.size());
    m_amr->allocate(m_inceptionVoltagePlus, m_realm, m_phase, 1);
    m_amr->allocate(m_inceptionVoltageMinu, m_realm, m_phase, 1);

    DataOps::setValue(m_inceptionIntegralPlus, 0.0);
    DataOps::setValue(m_inceptionIntegralMinu, 0.0);
    DataOps::setValue(m_backgroundIonization, 0.0);
    DataOps::setValue(m_emissionRatesPlus, 0.0);
    DataOps::setValue(m_emissionRatesMinu, 0.0);
    DataOps::setValue(m_inceptionVoltagePlus, 0.0);
    DataOps::setValue(m_inceptionVoltageMinu, 0.0);

    break;
  }
  case Mode::Dynamic: {
    m_amr->allocate(m_inceptionIntegral, m_realm, m_phase, 1);
    m_amr->allocate(m_emissionRate, m_realm, m_phase, 1);
    m_amr->allocate(m_backgroundIonization, m_realm, m_phase, 1);

    DataOps::setValue(m_inceptionIntegral, 0.0);
    DataOps::setValue(m_emissionRate, 0.0);
    DataOps::setValue(m_backgroundIonization, 0.0);

    m_negativeIonSolver->allocate();    

    break;
  }
  default: {
    break;
  }
  }
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
  m_fieldSolver->setRho(m_rho);
  m_fieldSolver->setSigma(0.0);

  const bool converged =
    m_fieldSolver->solve(m_fieldSolver->getPotential(), m_fieldSolver->getRho(), m_fieldSolver->getSigma(), false);

  if (!converged) {
    MayDay::Error("StreamerInceptionStepper::initialData -- could not solve Poisson equation. Aborting application.");
  }

  // Get the gas-side electric field and give it to the particle solvers.
  DataOps::copy(m_electricFieldPlus, m_amr->alias(phase::gas, m_fieldSolver->getElectricField()));

  m_amr->averageDown(m_electricFieldPlus, m_realm, m_phase);
  m_amr->averageDown(m_electricFieldPlus, m_realm, m_phase);

  m_amr->interpToCentroids(m_electricFieldPlus, m_realm, m_phase);

  // Also store -E
  DataOps::copy(m_electricFieldMinu, m_electricFieldPlus);
  DataOps::scale(m_electricFieldMinu, -1.0);
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
  m_tracerParticleSolver->registerOperators();

  if(m_mode == Mode::Dynamic) {
    m_negativeIonSolver->registerOperators();
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::parseOptions()
{
  CH_TIME("StreamerInceptionStepper::parseOptions");

  ParmParse pp("StreamerInceptionStepper");

  Real voltageLo    = 0.0;
  Real voltageHi    = 0.0;
  Real stepSize     = 0.0;
  int  voltageSteps = 0;

  std::string str1;
  std::string str2;

  m_plotPoisson              = false;
  m_plotTracer               = false;  
  m_plotNegativeIons         = false;    
  m_plotInceptionIntegral    = false;
  m_plotInceptionVoltage     = false;
  m_plotBackgroundIonization = false;
  m_plotFieldEmission        = false;
  m_plotAlpha                = false;
  m_plotEta                  = false;

  // Read in parameters
  pp.get("verbosity", m_verbosity);
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::parseOptions" << endl;
  }

  pp.get("voltage_lo", voltageLo);
  pp.get("voltage_hi", voltageHi);
  pp.get("voltage_steps", voltageSteps);
  pp.get("K_inception", m_inceptionK);
  pp.get("print_report", m_writeReport);
  pp.get("output_file", m_outputFile);
  pp.get("max_dt", m_maxDt);
  pp.get("min_dt", m_minDt);    

  pp.get("mode", str1);
  if (str1 == "none") {
    m_mode = Mode::None;
  }
  else if (str1 == "static") {
    m_mode = Mode::Static;
  }
  else if (str1 == "dynamic") {
    m_mode = Mode::Dynamic;
  }
  else {
    MayDay::Error("Expected 'none', 'static', or 'dynamic' for 'StreamerInceptionStepper.mode'");
  }

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

  if (str2 == "dx") {
    m_useFixedInceptionStep = false;
    m_relativeInceptionStep = stepSize;
  }
  else if (str2 == "fixed") {
    m_useFixedInceptionStep = true;
    m_fixedInceptionStep    = stepSize;
  }
  else {
    MayDay::Error("Expected 'dx' or 'fixed' for 'StreamerInceptionStepper.inception_alg'");
  }

  // Get the time step.
  pp.get("dt", str1, 0);
  if (str1 == "fixed") {
    m_useFixedDt = true;
    
    pp.get("dt", m_fixedDt, 1);
  }
  else if (str1 == "cfl") {
    m_useFixedDt = false;

    pp.get("dt", m_cfl, 1);
  }
  else {
    MayDay::Error("Expected 'fixed' or 'cfl' for 'StreamerInceptionStepper.dt'");
  }

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
      else if (plotVars[i] == "neg_ions") {
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

  // 
  if(m_plotPoisson) {
    m_fieldSolver->writePlotData(a_output, a_icomp);
    a_plotVariableNames.append(m_fieldSolver->getPlotVariableNames());    
  }

  if(m_plotTracer) {
    m_tracerParticleSolver->writePlotData(a_output, a_icomp);
    a_plotVariableNames.append(m_tracerParticleSolver->getPlotVariableNames());
  }



  // Write data -- this differs between the 'static' and 'dynamic' modes
  switch (m_mode) {
  case Mode::Static: {

  
    this->writePlotDataStatic(a_output, a_plotVariableNames, a_icomp);

    break;
  }
  case Mode::Dynamic:

    this->writePlotDataDynamic(a_output, a_plotVariableNames, a_icomp);

    break;
  default: {
    break;
  }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::writePlotDataStatic(EBAMRCellData&       a_output,
                                                    Vector<std::string>& a_plotVariableNames,
                                                    int&                 a_icomp) const noexcept
{
  CH_TIME("StreamerInceptionStepper::writePlotDataStatic");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::writePlotDataStatic" << endl;
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

      m_amr->averageDown(alpha, m_realm, m_phase);
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
      this->evaluateFunction(eta, m_voltageSweeps[i], m_alpha);

      m_amr->averageDown(eta, m_realm, m_phase);
      m_amr->interpGhost(eta, m_realm, m_phase);

      DataOps::copy(a_output, eta, Interval(a_icomp, a_icomp), Interval(0, 0));

      a_icomp++;

      const std::string varName = "Eta coefficient/ V = " + std::to_string(m_voltageSweeps[i]);

      a_plotVariableNames.push_back(varName);
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::writePlotDataDynamic(EBAMRCellData&       a_output,
                                                     Vector<std::string>& a_plotVariableNames,
                                                     int&                 a_icomp) const noexcept
{
  CH_TIME("StreamerInceptionStepper::writePlotDataDynamic");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::writePlotDataDynamic" << endl;
  }

  // Add potential and electric field to output
  {
    EBAMRCellData E;
    m_amr->allocate(E, m_realm, m_phase, SpaceDim);

    DataOps::copy(E, m_electricFieldPlus);
    DataOps::scale(E, m_voltageCurve(m_time));

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

  if(m_plotNegativeIons && m_mode == Mode::Dynamic) {
    m_negativeIonSolver->writePlotData(a_output, a_icomp);
    a_plotVariableNames.append(m_negativeIonSolver->getPlotVariableNames());
  }
      

  // Plot inception integral
  if (m_plotInceptionIntegral) {
    DataOps::copy(a_output, m_inceptionIntegral, Interval(a_icomp, a_icomp), Interval(0, 0));

    a_icomp++;

    a_plotVariableNames.push_back("Inception integral");
  }

  if (m_plotBackgroundIonization) {
    EBAMRCellData bgIonization;

    m_amr->allocate(bgIonization, m_realm, m_phase, 1);

    this->evaluateFunction(bgIonization, m_voltageCurve(m_time), m_backgroundRate);

    m_amr->averageDown(bgIonization, m_realm, m_phase);
    m_amr->interpGhost(bgIonization, m_realm, m_phase);

    DataOps::copy(a_output, bgIonization, Interval(a_icomp, a_icomp), Interval(0, 0));

    a_icomp++;

    a_plotVariableNames.push_back("Background rate");
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

    m_amr->averageDown(alpha, m_realm, m_phase);
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

    m_amr->averageDown(eta, m_realm, m_phase);
    m_amr->interpGhost(eta, m_realm, m_phase);

    DataOps::copy(a_output, eta, Interval(a_icomp, a_icomp), Interval(0, 0));

    a_icomp++;

    a_plotVariableNames.push_back("Townsend eta coefficient");
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

  if(m_plotPoisson) {
    ncomp += m_fieldSolver->getNumberOfPlotVariables();
  }
  if(m_plotTracer) {
    ncomp += m_tracerParticleSolver->getNumberOfPlotVariables();
  }


  switch (m_mode) {
  case Mode::Static: {
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
  case Mode::Dynamic: {

    // Always write the electric field.
    ncomp += SpaceDim;

    if (m_plotInceptionIntegral) {
      ncomp += 1;
    }
    if(m_plotNegativeIons) {
      ncomp += m_negativeIonSolver->getNumberOfPlotVariables();
    }          
    if (m_plotBackgroundIonization) {
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

template <typename P, typename F>
const EBAMRCellData* const
StreamerInceptionStepper<P, F>::getElectricField() const noexcept
{
  CH_TIME("StreamerInceptionStepper::getElectricField");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::getElectricField" << endl;
  }

  return &m_electricFieldPlus;
}

template <typename P, typename F>
Real
StreamerInceptionStepper<P, F>::computeDt()
{
  CH_TIME("StreamerInceptionStepper::computeDt");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeDt" << endl;
  }

  Real dt = std::numeric_limits<Real>::max();
  
  if(m_useFixedDt) {
    dt = m_fixedDt;
  }
  else {
    dt = m_cfl * this->computeNegativeIonCFL(m_voltageCurve(m_time));
  }

  dt = std::min(m_maxDt, dt);
  dt = std::max(m_minDt, dt);  

  return dt;
}

template <typename P, typename F>
Real
StreamerInceptionStepper<P, F>::advance(const Real a_dt)
{
  CH_TIME("StreamerInceptionStepper::advance");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::advance" << endl;
  }

  if (m_mode == Mode::Dynamic) {

    const Real curTime    = m_time + a_dt;
    const Real curVoltage = m_voltageCurve(curTime);
    
    m_tracerParticleSolver->setVelocity(m_electricFieldMinu);
    m_negativeIonSolver->setVelocity(m_electricFieldMinu);    

    // Seed new particles and compute the inception integral. 
    this->seedIonizationParticles(curVoltage);
    this->computeInceptionIntegralDynamic(curVoltage);
    
    // Move the negative ions. 
    this->advanceNegativeIons(curVoltage, a_dt);

    // Compute the critical volume and the electron appearance rate.
    const Real Vcr = this->computeCriticalVolumeDynamic();
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
  }
  else {
    MayDay::Error("StreamerInceptionStepper::advance -- must have 'StreamerInceptionStepper.mode = dynamic'");
  }

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
  m_tracerParticleSolver->setTime(a_step, a_time, a_dt);

  if(m_mode == Mode::Dynamic) {
    m_negativeIonSolver->setTime(a_step, a_time, a_dt);
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::printStepReport()
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

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::preRegrid(const int a_lmin, const int a_oldFinestLevel)
{
  CH_TIME("StreamerInceptionStepper::preRegrid");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::preRegrid" << endl;
  }

  m_fieldSolver->preRegrid(a_lmin, a_oldFinestLevel);
  m_tracerParticleSolver->preRegrid(a_lmin, a_oldFinestLevel);

  if(m_mode == Mode::Dynamic) {
    m_negativeIonSolver->preRegrid(a_lmin, a_oldFinestLevel);
  }
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
  m_tracerParticleSolver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  
  if(m_mode == Mode::Dynamic) {  
    m_negativeIonSolver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  }

  // Regrid velocity field
  m_amr->reallocate(m_electricFieldPlus, m_phase, a_lmin);
  m_amr->reallocate(m_electricFieldMinu, m_phase, a_lmin);

  switch (m_mode) {
  case Mode::Static: {
    m_amr->reallocate(m_inceptionIntegralPlus, m_phase, a_lmin);
    m_amr->reallocate(m_inceptionIntegralMinu, m_phase, a_lmin);
    m_amr->reallocate(m_backgroundIonization, m_phase, a_lmin);
    m_amr->reallocate(m_emissionRatesPlus, m_phase, a_lmin);
    m_amr->reallocate(m_emissionRatesMinu, m_phase, a_lmin);
    m_amr->reallocate(m_inceptionVoltagePlus, m_phase, a_lmin);
    m_amr->reallocate(m_inceptionVoltageMinu, m_phase, a_lmin);

    break;
  }
  case Mode::Dynamic: {

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

  const bool converged =
    m_fieldSolver->solve(m_fieldSolver->getPotential(), m_fieldSolver->getRho(), m_fieldSolver->getSigma());

  if (!converged) {
    MayDay::Error("StreamerInceptionStepper::initialData -- could not solve Poisson equation. Aborting application.");
  }

  // Get the gas-side electric field
  DataOps::copy(m_electricFieldPlus, m_amr->alias(phase::gas, m_fieldSolver->getElectricField()));

  m_amr->averageDown(m_electricFieldPlus, m_realm, m_phase);
  m_amr->interpGhostMG(m_electricFieldPlus, m_realm, m_phase);

  m_amr->interpToCentroids(m_electricFieldPlus, m_realm, m_phase);

  // Also store -E
  DataOps::copy(m_electricFieldMinu, m_electricFieldPlus);
  DataOps::scale(m_electricFieldMinu, -1.0);
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::postRegrid()
{
  CH_TIME("StreamerInceptionStepper::postRegrid");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::postRegrid" << endl;
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::setVoltageCurve(const std::function<Real(const Real& E)>& a_voltageCurve) noexcept
{
  CH_TIME("StreamerInceptionStepper::setVoltageCurve");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setVoltageCurve" << endl;
  }

  m_voltageCurve = a_voltageCurve;
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::setRho(const std::function<Real(const RealVect& x)>& a_rho) noexcept {
  CH_TIME("StreamerInceptionStepper::setRho");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setRho" << endl;
  }

  m_rho = a_rho;
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::setNegativeIonDensity(const std::function<Real(const RealVect x)>& a_density) noexcept {
  CH_TIME("StreamerInceptionStepper::setNegativeIonDensity");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setNegativeIonDensity" << endl;
  }

  m_negativeIonDensity = a_density;
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::setNegativeIonMobility(const std::function<Real(const Real x)>& a_mobility) noexcept {
  CH_TIME("StreamerInceptionStepper::setNegativeIonMobility");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setNegativeIonMobility" << endl;
  }

  m_negativeIonMobility = a_mobility;
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
StreamerInceptionStepper<P, F>::setEta(const std::function<Real(const Real& E)>& a_eta) noexcept
{
  CH_TIME("StreamerInceptionStepper::setEta");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setEta" << endl;
  }

  m_eta = a_eta;
}

template <typename P, typename F>
const std::function<Real(const Real& E)>&
StreamerInceptionStepper<P, F>::getAlpha() const noexcept
{
  return m_alpha;
}

template <typename P, typename F>
const std::function<Real(const Real& E)>&
StreamerInceptionStepper<P, F>::getEta() const noexcept
{
  return m_eta;
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::setBackgroundRate(const std::function<Real(const Real& E)>& a_backgroundRate) noexcept
{
  CH_TIME("StreamerInceptionStepper::setBackgroundRate");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setBackgroundRate" << endl;
  }

  m_backgroundRate = a_backgroundRate;
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
Mode
StreamerInceptionStepper<P, F>::getMode() const noexcept
{
  return m_mode;
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::seedUniformParticles() noexcept
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
StreamerInceptionStepper<P, F>::seedIonizationParticles(const Real& a_voltage) noexcept
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

            // clang-format off
				 const Real E = a_voltage * RealVect(D_DECL(electricFieldReg(iv, 0),
									    electricFieldReg(iv, 1),
									    electricFieldReg(iv, 2))).vectorLength();
            // clang-format on

            const Real alpha = m_alpha(E);
            const Real eta   = m_eta(E);

            if (alpha > eta) {
              const RealVect pos = m_amr->getProbLo() + (0.5 * RealVect::Unit + RealVect(iv)) * dx;
              particles.add(P(pos));
            }
          }
        };

        auto irregularKernel = [&](const VolIndex& vof) -> void {
          const IntVect iv = vof.gridIndex();

          if (validCells(iv, 0) && ebisbox.isIrregular(iv)) {

            // clang-format off
	    const Real E = a_voltage * RealVect(D_DECL(electricField(vof, 0),
						       electricField(vof, 1),
						       electricField(vof, 2))).vectorLength();
            // clang-format on

            const Real alpha = m_alpha(E);
            const Real eta   = m_eta(E);

            if (alpha > eta) {
              const RealVect pos = m_amr->getProbLo() + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);
              particles.add(P(pos));
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

        // Set mass to zero and do a backup of the initial position (we need to rewind particles
        // back later on). The scalar<0> flag is the flag we use for start-stop criteria for the
        // particle integration.
        p.weight()           = 0.0;
        p.template real<0>() = 1.0;
        p.template vect<0>() = p.position();
      }
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::seedNegativeIons() noexcept
{
  CH_TIME("StreamerInceptionStepper::seedNegativeIons");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::seedNegativeIons" << endl;
  }

  CH_assert(!(m_negativeIonSolver.isNull()));

  if(m_negativeIonSolver.isNull()) {
    MayDay::Abort("shit");
  }

  ParticleContainer<P>& amrParticles = m_negativeIonSolver->getParticles();
  amrParticles.clearParticles();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

    const LevelData<BaseFab<bool>>& validCellsLD = *m_amr->getValidCells(m_realm)[lvl];

    const Real dx  = m_amr->getDx()[lvl];
    const Real vol = std::pow(dx, SpaceDim);

    ParticleData<P>& levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const EBISBox&       ebisbox    = ebisl[dit()];
      const BaseFab<bool>& validCells = validCellsLD[dit()];

      List<P>& particles = levelParticles[dit()].listItems();

      if (!ebisbox.isAllCovered()) {
        auto regularKernel = [&](const IntVect& iv) -> void {
          if (validCells(iv, 0) && ebisbox.isRegular(iv)) {
	    const RealVect pos = m_amr->getProbLo() + (0.5 * RealVect::Unit + RealVect(iv)) * dx;

	    P p;
	    
	    p.position() = pos;
	    p.weight()   = m_negativeIonDensity(pos) * vol;

	    if(p.weight() > std::numeric_limits<Real>::min()) {
	      particles.add(p);
	    }
	  }
			     };

        auto irregularKernel = [&](const VolIndex& vof) -> void {
          const IntVect iv = vof.gridIndex();

          if (validCells(iv, 0) && ebisbox.isIrregular(iv)) {
	    const RealVect pos = m_amr->getProbLo() + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);

	    P p;

	    p.position() = pos;
	    p.weight()   = m_negativeIonDensity(pos) * vol * ebisbox.volFrac(vof);

	    if(p.weight() > std::numeric_limits<Real>::min()) {	    
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
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::postInitialize()
{
  CH_TIME("StreamerInceptionStepper::postInitialize");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::postInitialize" << endl;
  }

  switch (m_mode) {
  case Mode::Static: {

    this->computeInceptionIntegralStatic();
    this->computeBackgroundIonizationStatic();
    this->computeFieldEmissionStatic();
    this->computeInceptionVoltageVolume();
    this->computeCriticalVolumeStatic();
    this->computeFirstElectron();

    break;
  }
  case Mode::Dynamic: {

    // Seed particles and
    this->seedIonizationParticles(m_voltageCurve(m_time));

    m_tracerParticleSolver->setVelocity(m_electricFieldMinu);
    m_negativeIonSolver->setVelocity(m_electricFieldMinu);    

    // Compute the 'Rdot' value at current time (which is probably t = 0).
    this->computeInceptionIntegralDynamic(m_voltageCurve(m_time));

    // Get the maximum K-value
    Real maxK = -std::numeric_limits<Real>::max();
    Real minK = +std::numeric_limits<Real>::max();

    DataOps::getMaxMin(maxK, minK, m_inceptionIntegral, 0);    

    m_Rdot.emplace_back(m_time, this->computeRdot(m_voltageCurve(m_time)));
    m_criticalVolume.emplace_back(m_time, this->computeCriticalVolumeDynamic());
    m_inceptionProbability.emplace_back(m_time, 0.0);
    m_maxK.emplace_back(m_time, maxK);

    // Seed the negative ions
    this->seedNegativeIons();

    break;
  }
  default: {
    break;
  }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::computeInceptionIntegralStatic() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeInceptionIntegralStatic");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeInceptionIntegralStatic" << endl;
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
        MayDay::Error("StreamerInceptionStepper::computeInceptionIntegralStatic -- logic bust");

        break;
      }
      }

      // Deposit particles on mesh and copy to relevant data holder.
      if (p > 0.0) {
        m_tracerParticleSolver->deposit(Kplus);

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
        m_tracerParticleSolver->deposit(Kminu);

        m_amr->averageDown(Kminu, m_realm, m_phase);
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

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::computeInceptionIntegralDynamic(const Real& a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::computeInceptionIntegralDynamic");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeInceptionIntegralDynamic" << endl;
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
    MayDay::Error("StreamerInceptionStepper::computeInceptionIntegralDynamic - logic bust");

    break;
  }
  }

  // Deposit the particles onto m_inceptionIntegral
  m_tracerParticleSolver->deposit(m_inceptionIntegral);

  m_amr->averageDown(m_inceptionIntegral, m_realm, m_phase);
  m_amr->interpGhost(m_inceptionIntegral, m_realm, m_phase);
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

  ParticleContainer<P>& amrParticles = m_tracerParticleSolver->getParticles();

  m_tracerParticleSolver->remap();
  m_tracerParticleSolver->interpolateVelocities();

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

          const RealVect v        = p.velocity();
          const Real     dt       = step / v.vectorLength();
          const Real     alpha    = m_alpha(a_voltage * v.vectorLength());
          const Real     eta      = m_eta(a_voltage * v.vectorLength());
          const Real     alphaEff = alpha - eta;

          // Update particle position
          const RealVect newPos = p.position() + dt * v;

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
  //          T += 0.5 * dx * [alpha_eff(E(p^k)) + alpha_eff(E(p^(k+1)))],
  //
  //       where dx = |p^(k+1) - p^k|.

  const RealVect probLo = m_amr->getProbLo();
  const RealVect probHi = m_amr->getProbHi();

  // Allocate a data holder for holding the processed particles. This
  // will be faster because then we only have to iterate through the
  // particles that are actually still moving.
  ParticleContainer<P> amrProcessedParticles;
  m_amr->allocate(amrProcessedParticles, 0, m_realm);

  ParticleContainer<P>& amrParticles = m_tracerParticleSolver->getParticles();

  m_tracerParticleSolver->remap();
  m_tracerParticleSolver->interpolateVelocities();

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
          const Real     alpha    = m_alpha(a_voltage * v.vectorLength());
          const Real     eta      = m_eta(a_voltage * v.vectorLength());
          const Real     alphaEff = alpha - eta;

          // Storage alpha(p^k), v(p^k), and the step size.
          p.template real<0>() = alphaEff;
          p.template real<1>() = dt;
          p.template vect<1>() = v;

          // Update position.
          const RealVect newPos = p.position() + v * dt;

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

          // Note the weird subtraction since p.position() was updated
          // to p^k + dt * v^k.
          const RealVect newPos = p.position() + 0.5 * dt * (vk1 - vk);

          // Compute new alpha.
          const Real alpha1 = p.template real<0>();
          const Real alpha  = m_alpha(a_voltage * vk1.vectorLength());
          const Real eta    = m_eta(a_voltage * vk1.vectorLength());
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

template <typename P, typename F>
void StreamerInceptionStepper<P, F>::advanceNegativeIons(const Real& a_voltage, const Real& a_dt) noexcept
{
  CH_TIME("StreamerInceptionStepper::advanceNegativeIons");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::advanceNegativeIons" << endl;
  }

  switch(m_inceptionAlgorithm) {
  case IntegrationAlgorithm::Euler: {
    this->advanceNegativeIonsEuler(a_voltage, a_dt);

    break;
  };
  case IntegrationAlgorithm::Trapezoidal: {
    this->advanceNegativeIonsTrapezoidal(a_voltage, a_dt);

    break;
  };  
  default: {
    MayDay::Error("StreamerInceptionStepper::advanceNegativeIons - logic bust");

    break;
  }
  }
}

template <typename P, typename F>
void StreamerInceptionStepper<P, F>::advanceNegativeIonsEuler(const Real& a_voltage, const Real& a_dt) noexcept
{
  CH_TIME("StreamerInceptionStepper::advanceNegativeIonsEuler");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::advanceNegativeIonsEuler" << endl;
  }

  m_negativeIonSolver->remap();
  m_negativeIonSolver->interpolateVelocities();

  ParticleContainer<P>& amrParticles = m_negativeIonSolver->getParticles();  

  // Euler integration.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {

    for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
      List<P>& ionParticles = amrParticles[lvl][dit()].listItems();

      for (ListIterator<P> lit(ionParticles); lit.ok(); ++lit) {
	P& p = lit();

	// We have interpolated -E onto the particle position but the actual velocity
	// must be multiplied by the input voltage (because we solved with a voltage = 1V)
	// and by the ion mobility.
	const RealVect v        = a_voltage * p.velocity();
	const Real     Eminu    = v.vectorLength();
	const Real     mu       = m_negativeIonMobility(Eminu);

	p.template vect<0>() = p.position();
	p.position()         += mu * v * a_dt;
      }
    }
  }

  // Remove particles that wound up inside the EB. The remap function
  // removes those that left the domain.
  m_amr->removeCoveredParticlesIF(amrParticles, m_phase);  

  m_negativeIonSolver->remap();
  m_negativeIonSolver->interpolateVelocities();
}

template <typename P, typename F>
void StreamerInceptionStepper<P, F>::advanceNegativeIonsTrapezoidal(const Real& a_voltage, const Real& a_dt) noexcept
{
  CH_TIME("StreamerInceptionStepper::advanceNegativeIonsTrapezoidal");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::advanceNegativeIonsTrapezoidal" << endl;
  }

  m_negativeIonSolver->remap();
  m_negativeIonSolver->interpolateVelocities();

  ParticleContainer<P>& amrParticles = m_negativeIonSolver->getParticles();  

  // Euler stage.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {

    for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
      List<P>& ionParticles = amrParticles[lvl][dit()].listItems();

      for (ListIterator<P> lit(ionParticles); lit.ok(); ++lit) {
	P& p = lit();

	// We have interpolated -E onto the particle position but the actual velocity
	// must be multiplied by the input voltage (because we solved with a voltage = 1V)
	// and by the ion mobility.
	const RealVect v        = a_voltage * p.velocity();
	const Real     Eminu    = v.vectorLength();
	const Real     mu       = m_negativeIonMobility(Eminu);

	p.template vect<0>() = p.position();
	p.template vect<1>() = mu * v;
	
	p.position() += mu * v * a_dt;
      }
    }
  }

  m_negativeIonSolver->remap();
  m_amr->removeCoveredParticlesIF(amrParticles, m_phase);    
  m_negativeIonSolver->interpolateVelocities();

  // Second stage
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {

    for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
      List<P>& ionParticles = amrParticles[lvl][dit()].listItems();

      for (ListIterator<P> lit(ionParticles); lit.ok(); ++lit) {
	P& p = lit();

	const RealVect v        = a_voltage * p.velocity();
	const Real     Eminu    = v.vectorLength();
	const Real     mu       = m_negativeIonMobility(Eminu);

	const RealVect& xk = p.template vect<0>();	
	const RealVect& vk = p.template vect<1>();
	const RealVect  vk1 = mu * v;
	
	p.position() = xk + 0.5 * a_dt * (vk + vk1);
      }
    }
  }

  m_amr->removeCoveredParticlesIF(amrParticles, m_phase);  
  m_negativeIonSolver->remap();
  m_negativeIonSolver->interpolateVelocities();  
}

template <typename P, typename F>
Real
StreamerInceptionStepper<P, F>::computeRdot(const Real& a_voltage) const noexcept
{
  CH_TIME("StreamerInceptionStepper::computeRdot");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeRdot" << endl;
  }

  Real Rdot = 0.0;

  // TLDR: We are computing the integral of dne/dt * (1-eta/alpha) over the critical volume and integral(j/e * dS) over the "critical surface"

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

      const FArrayBox& electricFieldReg     = electricField.getFArrayBox();
      const FArrayBox& inceptionIntegralReg = inceptionIntegral.getFArrayBox();

      // Add contribution from background ionization in regular cells.
      auto regularKernel = [&](const IntVect& iv) -> void {
        if (ebisbox.isRegular(iv) && validCells(iv, 0)) {
          if (inceptionIntegralReg(iv, 0) >= m_inceptionK) {

            // clang-format off
	    const Real E = a_voltage * RealVect(D_DECL(electricFieldReg(iv, 0),
						       electricFieldReg(iv, 1),
						       electricFieldReg(iv, 2))).vectorLength();
            // clang-format on

            const Real alpha = m_alpha(E);
            const Real eta   = m_eta(E);
            const Real dndt  = m_backgroundRate(E);

            CH_assert(alpha >= eta);

            Rdot += dndt * (1.0 - eta / alpha) * vol;
          }
        }
      };

      // Add contribution from background ionization in cut-cells.
      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const IntVect iv = vof.gridIndex();
        if (ebisbox.isIrregular(iv) && validCells(iv, 0)) {
          if (inceptionIntegral(vof, 0) >= m_inceptionK) {

            // clang-format off
	    const Real E = a_voltage * RealVect(D_DECL(electricField(vof, 0),
						       electricField(vof, 1),
						       electricField(vof, 2))).vectorLength();
            // clang-format on

            const Real kappa    = ebisbox.volFrac(vof);
            const Real areaFrac = ebisbox.bndryArea(vof);
            const Real alpha    = m_alpha(E);
            const Real eta      = m_eta(E);
            const Real j        = m_fieldEmission(E);
            const Real dndt     = m_backgroundRate(E);

            Rdot += dndt * (1.0 - eta / alpha) * kappa * vol;
	    //            Rdot += j / Units::Qe * areaFrac * area;
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

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::rewindTracerParticles() noexcept
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

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::resetTracerParticles() noexcept
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

        p.weight()           = 0.0;
        p.template real<0>() = 1.0;
      }
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::computeBackgroundIonizationStatic() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeBackgroundIonizationStatic");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeBackgroundIonizationStatic" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {

      EBCellFAB& bgIonization    = (*m_backgroundIonization[lvl])[dit()];
      FArrayBox& bgIonizationReg = bgIonization.getFArrayBox();

      const EBCellFAB& electricField    = (*m_electricFieldPlus[lvl])[dit()];
      const FArrayBox& electricFieldReg = electricField.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        // clang-format off
        const Real E = RealVect(D_DECL(electricFieldReg(iv, 0),
				       electricFieldReg(iv, 1),
				       electricFieldReg(iv, 2))).vectorLength();
        // clang-format on

        for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
          bgIonizationReg(iv, i) = m_backgroundRate(m_voltageSweeps[i] * E);
        }
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        // clang-format off
        const Real E = RealVect(D_DECL(electricField(vof, 0),
				       electricField(vof, 1),
				       electricField(vof, 2))).vectorLength();
        // clang-format on

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
StreamerInceptionStepper<P, F>::computeFieldEmissionStatic() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeFieldEmissionStatic");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeFieldEmissionStatic" << endl;
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

      const EBCellFAB& electricField    = (*m_electricFieldPlus[lvl])[dit()];
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
StreamerInceptionStepper<P, F>::computeFieldEmission(EBAMRCellData& a_emissionRate,
                                                     const Real&    a_voltage) const noexcept
{
  CH_TIME("StreamerInceptionStepper::computeFieldEmission");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeFieldEmission" << endl;
  }

  CH_assert(a_emissionRate[0]->nComp() == 1);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const EBISBox& ebisbox = ebisl[dit()];

      // Here, emissionPlus is when we have the standard voltage (V = 1) and
      // emissionMinu is when the voltage is reverted.
      EBCellFAB& emission = (*a_emissionRate[lvl])[dit()];

      emission.setVal(0.0);

      const EBCellFAB& electricField    = (*m_electricFieldPlus[lvl])[dit()];
      const FArrayBox& electricFieldReg = electricField.getFArrayBox();

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        // clang-format off
        const RealVect E = a_voltage * RealVect(D_DECL(electricField(vof, 0),
						       electricField(vof, 1),
						       electricField(vof, 2)));
	// clang-format on. 

        if (E.dotProduct(ebisbox.normal(vof)) > 0.0) {
	  emission(vof, 0) = m_fieldEmission(E.vectorLength());
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
StreamerInceptionStepper<P, F>::evaluateFunction(EBAMRCellData& a_data, const Real& a_voltage, const std::function<Real(const Real E)>& a_func) const noexcept {
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
      EBCellFAB& data = (*a_data[lvl])[dit()];
      FArrayBox& dataReg = data.getFArrayBox();

      data.setVal(0.0);

      const EBCellFAB& electricField    = (*m_electricFieldPlus[lvl])[dit()];
      const FArrayBox& electricFieldReg = electricField.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {

	// clang-format off
        const RealVect E = a_voltage * RealVect(D_DECL(electricFieldReg(iv, 0),
						       electricFieldReg(iv, 1),
						       electricFieldReg(iv, 2)));
	// clang-format on. 

	dataReg(iv, 0) = a_func(E.vectorLength());
      };      

      auto irregularKernel = [&](const VolIndex& vof) -> void {

	// clang-format off
        const RealVect E = a_voltage * RealVect(D_DECL(electricField(vof, 0),
						       electricField(vof, 1),
						       electricField(vof, 2)));
	// clang-format on. 

	data(vof, 0) = a_func(E.vectorLength());
      };

      // Execute kernels
      Box cellBox = dbl[dit()];
      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      BoxLoops::loop(cellBox, regularKernel);      
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
    Vminu = interpVinc(m_maxKMinu);
  }

  return std::make_pair(Vplus, Vminu);
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::computeCriticalVolumeStatic() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeCriticalVolumeStatic");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeCriticalVolumeStatic" << endl;
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

template <typename P, typename F>
Real
StreamerInceptionStepper<P, F>::computeCriticalVolumeDynamic() const noexcept
{
  CH_TIME("StreamerInceptionStepper::computeCriticalVolumeDynamic");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeCriticalVolumeDynamic" << endl;
  }

  Real Vcr = 0.0;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); ++lvl) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    
    const LevelData<BaseFab<bool>>& validCellsLD = *m_amr->getValidCells(m_realm)[lvl];

    const Real dx = m_amr->getDx()[lvl];
    const Real vol = std::pow(dx, SpaceDim);

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      const EBISBox&       ebisbox    = ebisl[dit()];
      const BaseFab<bool>& validCells = validCellsLD[dit()];

      const EBCellFAB& inceptionIntegral = (*m_inceptionIntegral[lvl])[dit()];
      const FArrayBox& inceptionIntegralReg = inceptionIntegral.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
			     if (validCells(iv, 0) && ebisbox.isRegular(iv)) {
			       if(inceptionIntegralReg(iv, 0) >= m_inceptionK) {
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

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::computeFirstElectron() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeFirstElectron");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeFirstElectron" << endl;
  }

  // TLDR: This routine computes the average time to the appereance of an electron
  //       anywhere in the critical volume. We do this by assuming that we have a
  //       background ionization rate m_backgroundIonization(E).

  // Total rate.
  std::vector<Real> Rplus(m_voltageSweeps.size(), std::numeric_limits<Real>::min());
  std::vector<Real> Rminu(m_voltageSweeps.size(), std::numeric_limits<Real>::min());

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real&              dx    = m_amr->getDx()[lvl];

    const Real vol  = std::pow(dx, SpaceDim);
    const Real area = std::pow(dx, SpaceDim - 1);

    for (DataIterator dit(dbl); dit.ok(); ++dit) {

      const EBISBox&       ebisbox    = ebisl[dit()];
      const BaseFab<bool>& validCells = (*m_amr->getValidCells(m_realm)[lvl])[dit()];

      const EBCellFAB& backgroundRate        = (*m_backgroundIonization[lvl])[dit()];
      const EBCellFAB& inceptionIntegralPlus = (*m_inceptionIntegralPlus[lvl])[dit()];
      const EBCellFAB& inceptionIntegralMinu = (*m_inceptionIntegralMinu[lvl])[dit()];

      const FArrayBox& backgroundRateReg        = backgroundRate.getFArrayBox();
      const FArrayBox& inceptionIntegralPlusReg = inceptionIntegralPlus.getFArrayBox();
      const FArrayBox& inceptionIntegralMinuReg = inceptionIntegralMinu.getFArrayBox();

      const EBCellFAB& fieldEmissionPlus = (*m_emissionRatesPlus[lvl])[dit()];
      const EBCellFAB& fieldEmissionMinu = (*m_emissionRatesMinu[lvl])[dit()];

      // Background ionization kernel for regular cells.
      auto backgroundRegularKernel = [&](const IntVect& iv) -> void {
        if (validCells(iv) && ebisbox.isRegular(iv)) {
          for (size_t i = 0; i < m_voltageSweeps.size(); i++) {

            // Add background ionization rate.
            if (inceptionIntegralPlusReg(iv, i) >= m_inceptionK) {
              Rplus[i] += backgroundRateReg(iv, i) * vol;
            }

            if (inceptionIntegralMinuReg(iv, i) >= m_inceptionK) {
              Rminu[i] += backgroundRateReg(iv, i) * vol;
            }
          }
        }
      };

      // Background ionization kernel for irregular cells.
      auto backgroundIrregularKernel = [&](const VolIndex& vof) -> void {
        const IntVect iv = vof.gridIndex();
        if (validCells(iv, 0)) {
          const Real volFrac = ebisbox.volFrac(vof);

          for (size_t i = 0; i < m_voltageSweeps.size(); i++) {

            // Add background ionization rate.
            if (inceptionIntegralPlus(vof, i) >= m_inceptionK) {
              Rplus[i] += backgroundRate(vof, i) * vol * volFrac;
            }
            if (inceptionIntegralMinu(vof, i) >= m_inceptionK) {
              Rminu[i] += backgroundRate(vof, i) * vol * volFrac;
            }
          }
        }
      };

      // Field emission irregular kernel.
      auto fieldEmissionKernel = [&](const VolIndex& vof) -> void {
        const IntVect iv = vof.gridIndex();

        if (validCells(iv, 0)) {
          const Real areaFrac = ebisbox.bndryArea(vof);

          for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
            if (inceptionIntegralPlus(vof, i) >= m_inceptionK) {
              Rplus[i] += fieldEmissionPlus(vof, i) * area * areaFrac;
            }
            if (inceptionIntegralMinu(vof, i) >= m_inceptionK) {
              Rminu[i] += fieldEmissionMinu(vof, i) * area * areaFrac;
            }
          }
        }
      };

      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

      BoxLoops::loop(dbl[dit()], backgroundRegularKernel);
      BoxLoops::loop(vofit, backgroundIrregularKernel);
      BoxLoops::loop(vofit, fieldEmissionKernel);
    }
  }

  for (size_t i = 0; i < m_voltageSweeps.size(); i++) {
    m_firstElectronPlus.push_back(1. / ParallelOps::sum(Rplus[i]));
    m_firstElectronMinu.push_back(1. / ParallelOps::sum(Rminu[i]));
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::writeReportStatic() const noexcept
{
  CH_TIME("StreamerInceptionStepper::writeReportStatic");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::writeReportStatic" << endl;
  }

#ifdef CH_MPI
  if (procID() == 0) {
#endif

    std::ofstream output(m_outputFile, std::ofstream::out);
    
    const std::pair<Real, Real> Uinc = this->computeInceptionVoltages();

    // clang-format off
    output << "# ==================================================================================================================================" << "\n";
    if (std::isnan(Uinc.first) || std::isnan(Uinc.second)) {
      output << "# Could not compute inception voltage\n";
    }
    else {
      output << "# Inception voltage +/- = " << Uinc.first << "/" << Uinc.second << "\n";
    }
    
    output << "# ==================================================================================================================================" << "\n";    
    output << left << setw(15) << setfill(' ') << "# +/- Voltage";
    output << left << setw(15) << setfill(' ') << "Max K(+)";
    output << left << setw(15) << setfill(' ') << "Max K(-)";
    output << left << setw(20) << setfill(' ') << "Crit. vol(+)";
    output << left << setw(20) << setfill(' ') << "Crit. vol(-)";
    output << left << setw(20) << setfill(' ') << "First e(+)";
    output << left << setw(20) << setfill(' ') << "First e(-)" << "\n";
    output << "# ==================================================================================================================================" << "\n";        

    for (int i = 0; i < m_voltageSweeps.size(); i++) {
      output << left << setw(15) << setfill(' ') << m_voltageSweeps[i];
      output << left << setw(15) << setfill(' ') << m_maxKPlus[i];
      output << left << setw(15) << setfill(' ') << m_maxKMinu[i];
      output << left << setw(20) << setfill(' ') << m_criticalVolumePlus[i];
      output << left << setw(20) << setfill(' ') << m_criticalVolumeMinu[i];
      output << left << setw(20) << setfill(' ') << m_firstElectronPlus[i];
      output << left << setw(20) << setfill(' ') << m_firstElectronMinu[i];
      output << endl;
    }
    output << "# ==================================================================================================================================" << "\n";
    // clang-format on
#ifdef CH_MPI
  }
#endif
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::writeReportDynamic() const noexcept
{
  CH_TIME("StreamerInceptionStepper::writeReportDynamic");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::writeReportDynamic" << endl;
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

template <typename P, typename F>
Real
StreamerInceptionStepper<P, F>::computeNegativeIonCFL(const Real& a_voltage) const noexcept {
  CH_TIME("StreamerInceptionStepper::computeNegativeIonCFL");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeNegativeIonCFL" << endl;
  }

  Real dt = std::numeric_limits<Real>::max();

  m_negativeIonSolver->interpolateVelocities();

  ParticleContainer<P>& amrParticles = m_negativeIonSolver->getParticles();  

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Real dx = m_amr->getDx()[lvl];
    
    for (DataIterator dit(m_amr->getGrids(m_realm)[lvl]); dit.ok(); ++dit) {
      List<P>& ionParticles = amrParticles[lvl][dit()].listItems();

      for (ListIterator<P> lit(ionParticles); lit.ok(); ++lit) {
	P& p = lit();

	// We have interpolated -E onto the particle position but the actual velocity
	// must be multiplied by the input voltage (because we solved with a voltage = 1V)
	// and by the ion mobility.
	const RealVect v        = a_voltage * p.velocity();
	const Real     Eminu    = v.vectorLength();
	const Real     mu       = m_negativeIonMobility(Eminu);
	const Real     vMax     = mu * std::abs(v.maxDir(true));

	dt = std::min(dt, std::abs(dx/ (mu * Eminu)));
      }
    }
  }

  return ParallelOps::min(dt);
}

#include <CD_NamespaceFooter.H>

#endif
