/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_StreamerInceptionStepperImplem.H
  @brief  Implementation of CD_StreamerInceptionStepper.H
  @author Robert Marskar
*/

#ifndef CD_StreamerInceptionStepperImplem_H
#define CD_StreamerInceptionStepperImplem_H

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_StreamerInceptionStepper.H>
#include <CD_NamespaceHeader.H>
#include <iomanip> // for printing to terminal

// Tasks:
//
// 1. Compute inception integral with particle tracing
// 2. Compute critical volume
// 3. Compute negative ion drift times to electrode(s). Compute as time/(V * mobility)
// 4. Evaluate field emission rates.
// 5. Background ionization rates
// ** Add different integration methods: Simpson's, Simpson's 3/8, Boole's. Should create different moveParticles() for each function to compute alpha at different steps.
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
  
  // Default alpha
  m_alpha = [](const Real E) { return tanh(E); };
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
  auto voltage = [](const Real a_time) -> Real { return 1.0; };

  // Instantiate fieldSolver
  m_fieldSolver = RefCountedPtr<FieldSolver>(new F());
  m_fieldSolver->setVerbosity(m_verbosity);
  m_fieldSolver->parseOptions();
  m_fieldSolver->setAmr(m_amr);
  m_fieldSolver->setComputationalGeometry(m_computationalGeometry);
  m_fieldSolver->setVoltage(voltage);
  m_fieldSolver->setRealm(m_realm);
  m_fieldSolver->setTime(0, 0.0, 0.0);

  // Instantiate particleSolver
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

  m_amr->allocate(m_velocity, m_realm, m_phase, SpaceDim);
  m_amr->allocate(m_KValues, m_realm, m_phase, m_voltageSweeps.size());
  m_amr->allocate(m_VInception, m_realm, m_phase, 1);
  m_amr->allocate(m_driftTimes, m_realm, m_phase, 1 + m_voltageSweeps.size() );
  
  DataOps::setValue(m_KValues, 0.0);  
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
  const bool converged =
    m_fieldSolver->solve(m_fieldSolver->getPotential(), m_fieldSolver->getRho(), m_fieldSolver->getSigma());
  if (!converged) {
    MayDay::Error("StreamerInceptionStepper::initialData -- could not solve Poisson equation. Aborting application.");
  }

  // Get the gas-side electric field and give it to the particle solvers.
  DataOps::setValue(m_velocity, 0.0);
  DataOps::copy(m_velocity, m_amr->alias(phase::gas, m_fieldSolver->getElectricField()));
  m_amr->averageDown(m_velocity, m_realm, m_phase);
  m_amr->interpGhost(m_velocity, m_realm, m_phase);

  m_particleSolver->setVelocity(m_velocity);

  this->initialParticles();
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

  Real        voltageLo    = 0.0;
  Real        voltageHi    = 0.0;
  int         voltageSteps = 0;
  std::string str;

  // Read in parameters
  pp.get("verbosity", m_verbosity);
  pp.get("voltage_lo", voltageLo);
  pp.get("voltage_hi", voltageHi);
  pp.get("voltage_steps", voltageSteps);
  pp.get("K_inception", m_KInception);
  pp.get("integration", str);
  if (str == "trapezoidal") {
    m_integrationAlgorithm = IntegrationAlgorithm::Trapezoidal;
  }
  else if (str == "simpsons") {
    m_integrationAlgorithm = IntegrationAlgorithm::Simpsons;
  }
  else {
    MayDay::Error("StreamerInceptionStepper::StreamerInceptionStepper -- logic bust");
  }

  // Define voltages
  const Real dV = (voltageHi - voltageLo) / voltageSteps;
  m_voltageSweeps.push_back(voltageLo);
  for (int i = 0; i < voltageSteps; i++) {
    m_voltageSweeps.push_back(m_voltageSweeps[i] + dV);
  }  
}

// Not used in this program
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
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Interval srcInterval(0, m_voltageSweeps.size() - 1); 
    const Interval dstInterval(a_icomp, a_icomp + m_voltageSweeps.size() - 1); 
    m_KValues[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
  }
  a_icomp += m_voltageSweeps.size();

  std::string varName;
  for (const Real& V : m_voltageSweeps) {
    varName = "K at ";
    varName += std::to_string(V);
    varName += " V";
    a_plotVariableNames.push_back(varName);
  }
  
  // Write our drift time values.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Interval srcInterval(0, m_voltageSweeps.size());
    const Interval dstInterval(a_icomp, a_icomp + m_voltageSweeps.size());
    m_driftTimes[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
  }
  a_icomp += 1 + m_voltageSweeps.size();

  a_plotVariableNames.push_back("t*µV for positive charge");
  for (const Real& V : m_voltageSweeps) {
    varName = "t*µV for negative charge at ";
    varName += std::to_string(V);
    varName += " V";
    a_plotVariableNames.push_back(varName);
  }

  // Write our inception voltages
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Interval srcInterval(0, 0);
    const Interval dstInterval(a_icomp, a_icomp);
    m_VInception[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
  }
  a_icomp += 1;

  a_plotVariableNames.push_back("V_inception");
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

  ncomp += m_voltageSweeps.size(); // K values

  ncomp += 1 + m_voltageSweeps.size(); // Drift times

  ncomp += 1; // Inception voltages

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

// StreamerInceptionStepper is static. Any call to advance gives an error.
template <typename P, typename F>
Real
StreamerInceptionStepper<P, F>::advance(const Real a_dt)
{
  CH_TIME("StreamerInceptionStepper::advance");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::advance" << endl;
  }

  MayDay::Error(
    "StreamerInceptionStepper::advance -- this is a static solver. Calling the advance method is an error.");

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

// Combined implementation from FieldStepper and TracerParticleStepper
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
  m_amr->reallocate(m_velocity, m_phase, a_lmin);
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
StreamerInceptionStepper<P, F>::setVoltage(const Real a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::setVoltage");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setVoltage" << endl;
  }

  CH_assert(!m_fieldSolver.isNull());

  // Give field solver a constant voltage.
  m_fieldSolver->setVoltage([a_voltage](const Real& a_time) -> Real { return a_voltage; });
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::printStepReport()
{
  CH_TIME("StreamerInceptionStepper::printStepReport");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::printStepReport" << endl;
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::scaleElectrostatics(const Real& a_scale) noexcept
{
  CH_TIME("StreamerInceptionStepper::scaleElectrostatics");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::scaleElectrostatics" << endl;
  }

  CH_assert(!m_fieldSolver.isNull());

  MFAMRCellData& potential     = m_fieldSolver->getPotential();
  MFAMRCellData& electricField = m_fieldSolver->getElectricField();

  DataOps::scale(potential, a_scale);
  DataOps::scale(electricField, a_scale);
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::initialParticles()
{
  CH_TIME("StreamerInceptionStepper::initialParticles");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::initialParticles" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();
  amrParticles.clearParticles();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(amrParticles.getRealm())[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(amrParticles.getRealm(), m_phase)[lvl];

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
      
      const Box            cellBox    = dbl[dit()];      
      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];      

      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);

      // Define particles' masses (i.e. integration & drift time tracker) and integration flag
      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        p.mass()               = 0.0;
        p.template vect<0>()   = p.position();
        p.template scalar<0>() = 1.0;  // flag for keep integrating
      }
    }
  }

  // Remapping
  m_particleSolver->remap();

  // Remove particles inside the EB
  m_amr->removeCoveredParticlesIF(amrParticles, m_phase);
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::postInitialize()
{
  CH_TIME("StreamerInceptionStepper::postInitialize");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::postInitialize" << endl;
  }

  // Solve K integral for each particle in each voltage and store K values in m_KValues
  this->solveKValues();

  // Solve drift time integral in both directions and store values in m_driftTimes
  this->solveDriftTimes();
  
  // Interpolate between K values to find voltage giving K_inception and store values in m_VInception
  //this->solveVInception();

  // Compute critical volume of K values for each voltage
  //this->solveCriticalVolume();

  // Print voltages, max K values, and critical volumes to terminal
  //this->printToTerminal();
}


template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::solveKValues() {
  CH_TIME("StreamerInceptionStepper::solveKValues");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::solveKValues" << endl;
  }

  // Something transient storage we can deposit into. 
  EBAMRCellData K;
  m_amr->allocate(K, m_realm, m_phase, 1);

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();  
  
  // Integrate for each voltage
  for (int i = 0; i < m_voltageSweeps.size(); ++i) {

    // Update velocities and starting mass to zero.
    m_particleSolver->interpolateVelocities();     

    this->resetParticleMass();

    switch (m_integrationAlgorithm) {
    case IntegrationAlgorithm::Trapezoidal: {
      while (this->keepGoing()) {
        // Check integration flag and move particles if true
        this->K_moveParticlesTrapezoidal(m_voltageSweeps[i]);

        // Add integration parts to p.mass()
        this->K_addIntegrationTrapezoidal(m_voltageSweeps[i]);
      }
      
      break;
    }
    case IntegrationAlgorithm::Simpsons: {
      while (this->keepGoing()) {
        // Check integration flag and move particles if true
        this->K_moveParticlesSimpsons(m_voltageSweeps[i]);

        // Add integration parts to p.mass()
        this->K_addIntegrationSimpsons(m_voltageSweeps[i]);
      }

      break;
    }
    default: {
      MayDay::Error("StreamerInceptionStepper::solveKValues -- logic bust");
    }
    }

    // Move particles back to initial position. This includes remapping.
    this->moveParticlesInitialPosition();

    // Deposit
    m_particleSolver->deposit(K);

    // Set negative K values to 0
    this->zeroNegativeKValues(K);              
	      
    // Get max K value
    Real maxK;
    Real minK;
    DataOps::getMaxMin(maxK, minK, K, 0);
    m_maxKValues.push_back(maxK);

    // Copy K values to m_KValues
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      K[lvl]->copyTo(Interval(0,0), *m_KValues[lvl], Interval(i,i));
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::K_moveParticlesTrapezoidal(const Real& voltage)
{
  CH_TIME("StreamerInceptionStepper::K_moveParticlesTrapezoidal");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::K_moveParticlesTrapezoidal" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();  

  const size_t particlesBefore = amrParticles.getNumberOfValidParticesGlobal();

  // Move particles and remap
  this->K_moveParticlesOneStep(voltage, 0);
  
  const size_t particlesAfter = amrParticles.getNumberOfValidParticesGlobal();
  if (particlesBefore != particlesAfter) {
    MayDay::Warning("StreamerInceptionStepper::K_moveParticlesTrapezoidal -- logic bust, lost some particles");
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::K_addIntegrationTrapezoidal(const Real& voltage)
{
  CH_TIME("StreamerInceptionStepper::K_addIntegrationTrapezoidal");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::K_addIntegrationTrapezoidal" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(amrParticles.getRealm())[lvl];

    ParticleData<P>& levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        // Add integration parts if integration flag true
        if (p.template scalar<0>() > 0.0) {
          const RealVect  E               = p.velocity();
          const Real      Emagnitude      = voltage * E.vectorLength();
          const Real      endPointAlpha   = m_alpha(Emagnitude);
          const Real      dx              = p.template scalar<1>();	  
          const Real      startPointAlpha = p.template scalar<2>();
          p.mass() += 0.5 * (endPointAlpha + startPointAlpha) * dx;
        }
      }
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::K_moveParticlesSimpsons(const Real& voltage)
{
  CH_TIME("StreamerInceptionStepper::K_moveParticlesSimpsons");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::K_moveParticlesSimpsons" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();
  
  const size_t particlesBefore = amrParticles.getNumberOfValidParticesGlobal();

  // Step 1: Move particles dx/2 and remap
  this->K_moveParticlesOneStep(voltage, 1);

  // Step 2: Move particles dx/2 and remap
  this->K_moveParticlesOneStep(voltage, 2);
  
  const size_t particlesAfter = amrParticles.getNumberOfValidParticesGlobal();
  if (particlesBefore != particlesAfter) {
    MayDay::Warning("StreamerInceptionStepper::K_moveParticlesSimpsons -- logic bust, lost some particles");
  }
}


template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::K_addIntegrationSimpsons(const Real& voltage)
{
  CH_TIME("StreamerInceptionStepper::K_addIntegrationSimpsons");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::K_addIntegrationSimpsons" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(amrParticles.getRealm())[lvl];

    ParticleData<P>& levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        // Add integration parts if integration flag true
        if (p.template scalar<0>() > 0) {
          const RealVect  E               = p.velocity();
          const Real      Emagnitude      = voltage * E.vectorLength();
          const Real      endPointAlpha   = m_alpha(Emagnitude);
          const Real      dx   = p.template scalar<1>();
          const Real      startPointAlpha = p.template scalar<2>();
          const Real      midPointAlpha   = p.template scalar<3>();
          p.mass() += (dx / 6) * (startPointAlpha + 4 * midPointAlpha + endPointAlpha);
        }
      }
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::K_moveParticlesOneStep(const Real& voltage, const int& step)
{
  CH_TIME("StreamerInceptionStepper::K_moveParticlesOneStep");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::K_moveParticlesOneStep" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  // Get end points of grid
  const RealVect probLo = m_amr->getProbLo();
  const RealVect probHi = m_amr->getProbHi();
  
  // Move particles
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    // Particles move half a cell on this level. 
    const Real dx = 0.5*m_amr->getDx()[lvl];
    
    const DisjointBoxLayout& dbl            = m_amr->getGrids(amrParticles.getRealm())[lvl];
    
    ParticleData<P>&         levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

	// If integration flag true, move particle and check integration conditions
        if (p.template scalar<0>() > 0.0) {
          const RealVect  E                     = p.velocity();
          const RealVect  normalizedVelocity    = E / E.vectorLength();
          const Real      Emagnitude            = voltage * E.vectorLength();
          const Real      currentAlpha          = m_alpha(Emagnitude);

	  // Calculate new position based on integration algorithm
	  RealVect newPosition;
	  switch (m_integrationAlgorithm) {
	  case IntegrationAlgorithm::Trapezoidal: {
	    newPosition = p.position() + normalizedVelocity * dx;

	    // Store values to particle for subsequent integration
	    p.template scalar<1>() = dx;
	    p.template scalar<2>() = currentAlpha;
	    
	    break;
	  }
	  case IntegrationAlgorithm::Simpsons: {	      
	    newPosition = p.position() + normalizedVelocity * 0.5 * dx;

	    // Store values to particle for subsequent integration	    
	    if (step == 1) {
	      p.template scalar<1>() = dx;
	      p.template scalar<2>() = currentAlpha;
	    }
	    else if (step == 2) {
	      p.template scalar<3>() = currentAlpha;
	    }
	    else {
	      MayDay::Error("StreamerInceptionStepper::K_moveParticlesOneStep -- logic bust");
	    }
	    
	    break;
	  }
	  default: {
	    MayDay::Error("StreamerInceptionStepper::K_moveParticlesOneStep -- logic bust");
	  }
	  }

	  // Check keep integrating criteria. Move particle if fulfilled, flag if not.
	  if (currentAlpha <= 0.0) {
	    p.template scalar<0>() = -1.0;
	  }
	  else if (this->particleOutsideGrid(newPosition, probLo, probHi)) {
	    p.template scalar<0>() = -1.0;
	  }
	  else if (this->particleInsideElectrode(newPosition)) {
	    p.template scalar<0>() = -1.0;
	  }
	  else {
	    p.position() = newPosition;
	  }	  
        }
      }
    }
  }
  
  m_particleSolver->remap();
  m_particleSolver->interpolateVelocities();
}  
  

template <typename P, typename F>
bool
StreamerInceptionStepper<P, F>::keepGoing() const
{
  CH_TIME("StreamerInceptionStepper::keepGoing");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::keepGoing" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  // If keepGoing > 0 then we keep integrating. Otherwise we stop.
  int keepGoing = -1;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(amrParticles.getRealm())[lvl];

    ParticleData<P>& levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
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
StreamerInceptionStepper<P, F>::moveParticlesInitialPosition()
{
  CH_TIME("StreamerInceptionStepper::moveParticlesInitialPosition");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::moveParticlesInitialPosition" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();
  
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(amrParticles.getRealm())[lvl];

    ParticleData<P>& levelParticles = amrParticles[lvl];
    for (DataIterator dit(dbl); dit.ok(); ++dit) {

      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        p.position() = p.template vect<0>();
      }
    }
  }
  
  m_particleSolver->remap();
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::zeroNegativeKValues(EBAMRCellData& K)
{
  CH_TIME("StreamerInceptionStepper::zeroNegativeKValues");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::zeroNegativeKValues" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); ++lvl) {      
    LevelData<EBCellFAB>& levelData = *K[lvl];    
    const DisjointBoxLayout& dbl = levelData.disjointBoxLayout();
    
    const EBISLayout&        ebisl = m_amr->getEBISLayout(amrParticles.getRealm(), m_phase)[lvl];    

    const Real dx = m_amr->getDx()[lvl];
    
    ParticleData<P>& levelParticles = amrParticles[lvl]; 
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      const EBISBox&       ebisbox    = ebisl[dit()]; // used in irregularKernel
      
      EBCellFAB& patchData = levelData[dit()];  
      BaseFab<Real>& regularData = patchData.getSingleValuedFAB();

      List<P>& particles = levelParticles[dit()].listItems();

      auto regularKernel = [&](const IntVect& iv) -> void {
	if (regularData(iv, 0) < 0) {
	  regularData(iv, 0) = 0.0;
	}
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
	const IntVect iv = vof.gridIndex();
	if (regularData(iv, 0) < 0) {
	  regularData(iv, 0) = 0.0; 
	}
      };

      const Box cellBox = dbl[dit()];
      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];      

      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);      	  
    }
  }
}


template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::resetParticleMass()
{
  CH_TIME("StreamerInceptionStepper::resetParticleMass");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::resetParticleMass" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();
  
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl            = m_amr->getGrids(amrParticles.getRealm())[lvl];
    ParticleData<P>&         levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {

        // Reset mass and flag for keep integrating.
        P& p = lit();

        p.mass() = 0.0;
        p.template scalar<0>() = 1.0;
      }
    }
  }
}


template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::solveDriftTimes() {
  CH_TIME("StreamerInceptionStepper::solveDriftTimes");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::solveDriftTimes" << endl;
  }
  
  // Something transient storage we can deposit into. 
  EBAMRCellData t;
  m_amr->allocate(t, m_realm, m_phase, 1);

  // Solve drift time for positive ions
  solveDriftTime(t, 1, -1);

  // Solve drift time for negative ions for each voltage
  for (int voltInd = 0; voltInd < m_voltageSweeps.size(); voltInd++) {
    solveDriftTime(t, -1, voltInd);
  }
}


template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::solveDriftTime(EBAMRCellData& t, const int& dxSign, const int& voltInd)
{
  CH_TIME("StreamerInceptionStepper::solveDriftTime");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::solveDriftTime" << endl;
  }

  // Update velocities and starting mass to zero.
  m_particleSolver->interpolateVelocities();  
    
  this->resetParticleMass();  

  // Compute drift time integrals with input-based algorithm. 
  switch (m_integrationAlgorithm) {
  case IntegrationAlgorithm::Trapezoidal: {
    while (this->keepGoing()) {
      // Check integration flag and move particles if true
      this->driftTime_moveParticlesTrapezoidal(dxSign, voltInd);

      // Add integration parts to p.mass()
      this->driftTime_addIntegrationTrapezoidal();
    }

    break;
  }
  case IntegrationAlgorithm::Simpsons: {
    while (this->keepGoing()) {
      // Check integration flag and move particles if true	
      this->driftTime_moveParticlesSimpsons(dxSign, voltInd);
	
      // Add integration parts to p.mass()
      this->driftTime_addIntegrationSimpsons();
    }

    break;
  }
  default: {
    MayDay::Error("StreamerInceptionStepper::solveDriftTimes -- logic bust");
  }
  }
    
  // Move particles back to initial position. This includes remapping.
  this->moveParticlesInitialPosition();
  
  // Deposit onto t
  m_particleSolver->deposit(t);

  // Find storage index according to ion charge and voltage
  int index;
  if (dxSign == 1) {
    index = 0;
  }
  else if (dxSign == -1) {
    index = 1 + voltInd;
  }

  // Copy to m_driftTimes
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    t[lvl]->copyTo(Interval(0,0), *m_driftTimes[lvl], Interval(index, index));
  }       
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::driftTime_moveParticlesTrapezoidal(const int& dxSign, const int& voltInd) {
  CH_TIME("StreamerInceptionStepper::driftTime_moveParticlesTrapezoidal");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::driftTime_moveParticlesTrapezoidal" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  const size_t particlesBefore = amrParticles.getNumberOfValidParticesGlobal();  

  // Move particles and remap
  this->driftTime_moveParticlesOneStep(dxSign, 0, voltInd);

  const size_t particlesAfter = amrParticles.getNumberOfValidParticesGlobal();
  if (particlesBefore != particlesAfter) {
    MayDay::Warning("StreamerInceptionStepper::driftTime_moveParticlesTrapezoidal -- logic bust, lost some particles");
  }
}


template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::driftTime_addIntegrationTrapezoidal()
{
  CH_TIME("StreamerInceptionStepper:driftTime_addIntegrationTrapezoidal");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::driftTime_addIntegrationTrapezoidal" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(amrParticles.getRealm())[lvl];

    ParticleData<P>& levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        // Add integration parts if integration flag true
        if (p.template scalar<0>() > 0.0) {
          const RealVect  E               = p.velocity();
          const Real      Emagnitude      = E.vectorLength();
          const Real      dx   = p.template scalar<1>();	  
          const Real      startPointEmagnitude = p.template scalar<2>();
          p.mass() += (dx / 2) * ( (1/startPointEmagnitude) + (1/Emagnitude) );
        }
      }
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::driftTime_moveParticlesSimpsons(const int& dxSign, const int& voltInd) {
  CH_TIME("StreamerInceptionStepper::driftTime_moveparticlesSimpsons");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::driftTime_moveparticlesSimpsons" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  const size_t particlesBefore = amrParticles.getNumberOfValidParticesGlobal();  

  // Step 1: Move particles dx/2 and remap
  this->driftTime_moveParticlesOneStep(dxSign, 1, voltInd);

  // Step 2: Move particles dx/2 and remap
  this->driftTime_moveParticlesOneStep(dxSign, 2, voltInd);
  
  const size_t particlesAfter = amrParticles.getNumberOfValidParticesGlobal();
  if (particlesBefore != particlesAfter) {
    MayDay::Warning("StreamerInceptionStepper::driftTime_moveParticlesSimpsons -- logic bust, lost some particles");
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::driftTime_addIntegrationSimpsons()
{
  CH_TIME("StreamerInceptionStepper:driftTime_addIntegrationSimpsons");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::driftTime_addIntegrationSimpsons" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(amrParticles.getRealm())[lvl];

    ParticleData<P>& levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        // Add integration parts if integration flag true
        if (p.template scalar<0>() > 0.0) {
          const RealVect  E               = p.velocity();
          const Real      Emagnitude      = E.vectorLength();
          const Real      dx   = p.template scalar<1>();	  
          const Real      startPointEmagnitude = p.template scalar<2>();
	  const Real      midPointEmagnitude = p.template scalar<3>();
          p.mass() += (dx / 6) * ( (1/startPointEmagnitude) + 4 * (1/midPointEmagnitude) + (1/Emagnitude) );
        }
      }
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::driftTime_moveParticlesOneStep(const int& dxSign, const int& step, const int& voltInd)
{
  CH_TIME("StreamerInceptionStepper::driftTime_moveParticlesOneStep");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::driftTime_moveParticlesOneStep" << endl;    
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  const RealVect probLo = m_amr->getProbLo();
  const RealVect probHi = m_amr->getProbHi();  

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {    
    // Particles move half a cell on this level. 
    const Real dx = 0.5*m_amr->getDx()[lvl];

    const DisjointBoxLayout& dbl = m_amr->getGrids(amrParticles.getRealm())[lvl];

    ParticleData<P>& levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {      
      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();
	
	// If integration flag true, move particle and check integration conditions
	if (p.template scalar<0>() > 0.0) {
	  const RealVect E                     = p.velocity();
          const RealVect normalizedVelocity    = E / E.vectorLength();
          const Real     Emagnitude            = E.vectorLength();
	 
	  // Calculate new position based on integration algorithm and direction
	  RealVect newPosition;
	  switch (m_integrationAlgorithm) {	    
	  case IntegrationAlgorithm::Trapezoidal: {
	    newPosition = p.position() + normalizedVelocity * dxSign * dx;

	    // Store values to particle for subsequent integration
	    p.template scalar<1>() = dx;
	    p.template scalar<2>() = Emagnitude;	    
	    
	    break;
	  }
	  case IntegrationAlgorithm::Simpsons: {
	    newPosition = p.position() + normalizedVelocity * 0.5 * dxSign * dx;

	    // Store values to particle for subsequent integration
	    if (step == 1) {
	      p.template scalar<1>() = dx;
	      p.template scalar<2>() = Emagnitude;
	    }
	    else if (step == 2) {
	      p.template scalar<3>() = Emagnitude;
	    }
	    else {
	      MayDay::Error("StreamerInceptionStepper:driftTime_moveParticlesOneStep -- logic bust");
	    }
	    
	    break;
	  }
	  default: {
	    MayDay::Error("StreamerInceptionStepper::driftTime_moveParticlesOneStep -- logic bust");
	  }
	  }	  	  
	  
	  // Check keep integrating criteria. Move particle if fulfilled, flag if not.
	  if (this->particleOutsideGrid(newPosition, probLo, probHi)) {
            p.template scalar<0>() = -1.0;
          }
          else if (this->particleInsideElectrode(newPosition)) {
            p.template scalar<0>() = -1.0;
          }
          else {
            p.position() = newPosition;
          }	  	  
	}
      }
    }
  }
  
  m_particleSolver->remap();
  m_particleSolver->interpolateVelocities();

  // If negative ion, check if particles are inside critical volume
  if (dxSign == -1) {
    this->checkParticlesInsideCriticalVolume(voltInd);
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::checkParticlesInsideCriticalVolume(const int& voltInd)
{
  CH_TIME("StreamerInceptionStepper::checkParticlesInsideCriticalVolume");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::checkParticlesInsideCriticalVolume" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  const RealVect probLo = m_amr->getProbLo();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Real dx = 0.5*m_amr->getDx()[lvl];

    const DisjointBoxLayout& dbl            = m_amr->getGrids(amrParticles.getRealm())[lvl];
    
    ParticleData<P>&         levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      BaseFab<Real>& cellData = (*m_KValues[lvl])[dit()].getSingleValuedFAB();

      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
	P& p = lit();

	if (p.template scalar<0>() > 0.0) {
	  // Find integer spatial vector
	  RealVect rv = (p.position() - probLo)/dx - 0.5*RealVect::Unit;
	  IntVect iv;
	  for (int i = 0; i < SpaceDim; i++) {
	    iv[i] = round(rv[i]);
	  }

	  // Get K value of cell and check if greater than m_KInception. If so, flag particle.
	  Real K = cellData(iv, voltInd);
	  if (K >= m_KInception) {
	    p.template scalar<0>() = -1.0;
	  }
	}
      }
    }
  }
}



/* ISSUES: 
   - V_inception on edge of electrode (irregular). Remove? This is also true for K-values.
   - Visualization: logarithmic scale preferable, but not allowed by Visit because of 0 values (e.g. electrode)
*/
template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::solveVInception() {
  CH_TIME("StreamerInceptionStepper::solveVInception");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::solveVInception" << endl;
  }

  // Only solve if there are at least two voltages in m_voltageSweeps.
  if (m_voltageSweeps.size() < 2) {
    if (procID() == 0) {
      std::cout << "V_inception not solved -- need at least two voltages" << endl;
    }
    return;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();
  amrParticles.clearParticles();

  // Iterate through m_KValues data and calculate V_inception at each point
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); ++lvl) {      
    LevelData<EBCellFAB>& levelData = *m_KValues[lvl];    
    const DisjointBoxLayout& dbl = levelData.disjointBoxLayout();
    
    const EBISLayout&        ebisl = m_amr->getEBISLayout(amrParticles.getRealm(), m_phase)[lvl];

    const LevelData<BaseFab<bool>>& validCellsLD = *m_amr->getValidCells(m_realm)[lvl];    

    const Real dx = m_amr->getDx()[lvl];
    
    ParticleData<P>& levelParticles = amrParticles[lvl]; 
    
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      const EBISBox&       ebisbox    = ebisl[dit()]; 
      const BaseFab<bool>& validCells = validCellsLD[dit()];      
      
      const EBCellFAB& patchData = levelData[dit()];  
      const BaseFab<Real>& regularData = patchData.getSingleValuedFAB();

      List<P>& particles = levelParticles[dit()].listItems();

      auto regularKernel = [&](const IntVect& iv) -> void {       
	Real V_inception = 0.0;
	bool addP = false;
       
	// Calculate V_inception and draw particle only if largest K > 0.1
	if ( abs(regularData(iv, m_voltageSweeps.size()-1)) > 0.1 ) {
	  addP = true;
	  this->computeVInceptionCell(regularData, iv, V_inception);
	}
	
	// Add particle if addP
	if ( addP && validCells(iv, 0) && ebisbox.isRegular(iv) ) {
	  this->addParticleRegular(particles, iv, dx, V_inception);
	}	        
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
	Real V_inception = 0.0;
	const IntVect iv = vof.gridIndex();
	bool addP = true;

	// Calculate V_inception and draw particle only if largest K > 0.1
	if (abs(regularData(iv, m_voltageSweeps.size()-1)) > 0.1) {
	  addP = true;
	  this->computeVInceptionCell(regularData, iv, V_inception);
	}	

	// Add particle if addP
	if ( addP && validCells(iv) ) {
	  this->addParticleIrregular(particles, vof, ebisbox, dx, V_inception);
	}	
      };
      
      const Box cellBox = dbl[dit()];      
      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];      

      BoxLoops::loop(cellBox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);      
    }
  }

  // Remapping
  m_particleSolver->remap();

  // Remove particles inside the EB
  m_amr->removeCoveredParticlesIF(amrParticles, m_phase);

  // Deposit
  m_particleSolver->deposit(m_VInception);
}



template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::computeVInceptionCell(const BaseFab<Real>& regularData, const IntVect& iv, Real& V_inception)
{
  CH_TIME("StreamerInceptionStepper::computeVInceptionCell");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeVInceptionCell" << endl;
  } 

  // If inception smaller than smallest K value, interpolate with the two lowest K values
  if (regularData(iv, 0) >= m_KInception) {
    this->interpolate(m_voltageSweeps[0], m_voltageSweeps[1], regularData(iv, 0), regularData(iv, 1), m_KInception, V_inception);
  }

  // If inception greater than greatest K value, interpolate with the two largest K values
  else if (regularData(iv, m_voltageSweeps.size() - 1) <= m_KInception) {
    const int i1 = m_voltageSweeps.size() - 2;
    const int i2 = m_voltageSweeps.size() - 1;
    
    this->interpolate(m_voltageSweeps[i1], m_voltageSweeps[i2], regularData(iv, i1), regularData(iv, i2), m_KInception, V_inception);
  }

  // If inception within voltage interval, iterate through values until K value greater than m_KInception
  else {	  
    for (int val = 1; val < m_voltageSweeps.size(); val++) {
      if (regularData(iv, val) >= m_KInception) {
	this->interpolate(m_voltageSweeps[val-1], m_voltageSweeps[val], regularData(iv, val-1), regularData(iv, val), m_KInception, V_inception);
	break;
      }
    }
  }
  
  // Warning if V_inception<1 - there's probably something wrong, perhaps too few voltage_steps
  if (V_inception < 1) {
    MayDay::Error("StreamerInceptionStepper::computeVInceptionCell -- V_inception < 1. Something's probably wrong.");
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::interpolate(const Real& x1, const Real& x2, const Real& y1, const Real& y2, const Real& x, Real& y)
{ 
  CH_TIME("StreamerInceptionStepper::interpolate");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::interpolate" << endl;
  }

  const Real b  = (y2 - y1) / (x2 - x1);
  const Real y0 = y1 - b * x1;

  y = ((x - y0) / b);
}

 
template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::addParticleRegular(List<P>& particles, const IntVect iv, const Real dx, const Real mass)
{
  CH_TIME("StreamerInceptionStepper::addParticleRegular");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::addParticleRegular" << endl;
  }
  
  P p;
  p.position() = m_amr->getProbLo() + (0.5 * RealVect::Unit + RealVect(iv)) * dx;
  p.mass() = mass;
  particles.add(p);						   
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::addParticleIrregular(List<P>& particles, const VolIndex& vof, const EBISBox& ebisbox, const Real dx, const Real mass)
{
  CH_TIME("StreamerInceptionStepper::addParticleIrregular");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::addParticleIrregular" << endl;
  }

  P p;
  p.position() = m_amr->getProbLo() + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);
  p.mass() = mass;
  particles.add(p);  
}

  
// QUESTIONS TO ROBERT:
//   - validCells(iv, 0) eller validCells(iv, voltInd)? Har antatt førstnevnte.
template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::solveCriticalVolume() {
  CH_TIME("StreamerInceptionStepper::solveCriticalVolume");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::solveCriticalVolume" << endl;
  }
  
  // Solve critical volume of K values for each voltage
  for (int voltInd = 0; voltInd < m_voltageSweeps.size(); voltInd++) {
    Real criticalVolume = 0.0;
    
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); ++lvl) {
      LevelData<EBCellFAB>& levelData = *m_KValues[lvl];    
      const DisjointBoxLayout& dbl = levelData.disjointBoxLayout();
    
      const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

      const LevelData<BaseFab<bool>>& validCellsLD = *m_amr->getValidCells(m_realm)[lvl];    

      const Real dx = m_amr->getDx()[lvl];
       
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
	const EBISBox&       ebisbox    = ebisl[dit()]; 
	const BaseFab<bool>& validCells = validCellsLD[dit()];      
      
	EBCellFAB& patchData = levelData[dit()];  
	BaseFab<Real>& regularData = patchData.getSingleValuedFAB();

	auto regularKernel = [&](const IntVect& iv) -> void {
	  // Check if cell's K value >= m_KInception and if cell is valid
	  if ( (regularData(iv, voltInd) >= m_KInception) && validCells(iv, 0) && ebisbox.isRegular(iv)) {
	    // Add cell volume to critical volume
	    criticalVolume += pow(dx, SpaceDim);
	  }
	};

	auto irregularKernel = [&](const VolIndex& vof) -> void {
	  // Check if cell's K value >= m_KInception and if cell is valid
	  if ( (regularData(vof.gridIndex(), voltInd) >= m_KInception) && validCells(vof.gridIndex()) ) {
	    // Add cell volume to critical volume
	    criticalVolume += pow(dx, SpaceDim);
	  }
	};

	const Box cellBox = dbl[dit()];
	VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[dit()];

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
StreamerInceptionStepper<P, F>::printToTerminal() {
  CH_TIME("StreamerInceptionStepper::printToTerminal");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::printToTerminal" << endl;
  }

  auto printElement = [](const Real el) -> void {
    std::cout << left << setw(15) << setfill(' ') << el;
  };

  if (procID() == 0) {
    std::cout << left << setw(15) << setfill(' ') << "Voltage";
    std::cout << left << setw(15) << setfill(' ') << "Max K";
    std::cout << left << setw(15) << setfill(' ') << "Critical Volume";
    std::cout << "\n=============================================\n";

    for (int i = 0; i < m_voltageSweeps.size(); i++) {
      printElement(m_voltageSweeps[i]);
      printElement(m_maxKValues[i]);
      printElement(m_criticalVolumes[i]);
      std::cout << endl;
    }

    std::cout << "=============================================\n";
  }
}


template <typename P, typename F>
bool
StreamerInceptionStepper<P, F>::particleOutsideGrid(const RealVect pos, const RealVect probLo, const RealVect probHi)
{
  CH_TIME("StreamerInceptionStepper::particleOutsideGrid");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::particleOutsideGrid" << endl;
  }
  
  bool isOutside = false;
  
  for (int dir = 0; dir < SpaceDim; dir++) {
    if (pos[dir] <= probLo[dir] || pos[dir] >= probHi[dir]) {
      isOutside = true;
    }
  }

  return isOutside;
}

template <typename P, typename F>
bool
StreamerInceptionStepper<P, F>::particleInsideElectrode(const RealVect pos)
{
  CH_TIME("StreamerInceptionStepper::particleInsideElectrode");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::particleInsideElectrode" << endl;
  }
  
  const RefCountedPtr<BaseIF>& implicitFunction = m_amr->getBaseImplicitFunction(m_phase);

  return (implicitFunction->value(pos) >= 0.0) ? true : false;
}


#include <CD_NamespaceFooter.H>

#endif
