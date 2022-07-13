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
  m_amr->allocate(m_driftTimes, m_realm, m_phase, 2);
  
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
    MayDay::Error("StreamerInceptionStepper::initialize -- could not solve Poisson equation. Aborting application.");
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
    varName = "K for ";
    varName += std::to_string(V);
    varName += " V";
    a_plotVariableNames.push_back(varName);
  }
  
  // Write our drift time values.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Interval srcInterval(0, 1);
    const Interval dstInterval(a_icomp, a_icomp + 1);
    m_driftTimes[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
  }
  a_icomp += 2;

  a_plotVariableNames.push_back("t*µV for +E direction");
  a_plotVariableNames.push_back("t*µV for -E direction");

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

  ncomp += m_voltageSweeps.size();

  ncomp += 3;  // Drift times (+/- direction) and inception voltage

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
  //this->solveDriftTimes();
  
  // Interpolate between K values to find voltage giving K_inception and store values in m_VInception
  this->solveVInception();
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

    this->resetParticleMass(amrParticles);

    switch (m_integrationAlgorithm) {
    case IntegrationAlgorithm::Trapezoidal: {
      while (this->keepGoing(amrParticles)) {
        // Check integration flag and move particles if true
        this->K_moveParticlesTrapezoidal(amrParticles, m_voltageSweeps[i]);

        // Add integration parts to p.mass()
        this->K_addIntegrationTrapezoidal(amrParticles, m_voltageSweeps[i]);
      }
      
      break;
    }
    case IntegrationAlgorithm::Simpsons: {
      while (this->keepGoing(amrParticles)) {
        // Check integration flag and move particles if true
        this->K_moveParticlesSimpsons(amrParticles, m_voltageSweeps[i]);

        // Add integration parts to p.mass()
        this->K_addIntegrationSimpsons(amrParticles, m_voltageSweeps[i]);
      }

      break;
    }
    default: {
      MayDay::Error("StreamerInceptionStepper::solveKValues -- logic bust");
    }
    }

    // Move particles back to initial position. This includes remapping.
    this->moveParticlesInitialPosition(amrParticles);

    // Deposit
    m_particleSolver->deposit(K);

    // Set negative K values to 0
    this->zeroNegativeKValues(K, amrParticles);              
	      
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
StreamerInceptionStepper<P, F>::K_moveParticlesTrapezoidal(ParticleContainer<P>& amrParticles, const Real& voltage)
{
  CH_TIME("StreamerInceptionStepper::K_moveParticlesTrapezoidal");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::K_moveParticlesTrapezoidal" << endl;
  }

  const size_t particlesBefore = amrParticles.getNumberOfValidParticesGlobal();

  const RealVect probLo = m_amr->getProbLo();
  const RealVect probHi = m_amr->getProbHi(); 

  // Move particles
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
          const Real     Emagnitude            = voltage * E.vectorLength();
          const Real     startPointAlpha       = m_alpha(Emagnitude);

          const RealVect newPosition = p.position() + normalizedVelocity * dx;

          p.template scalar<1>() = dx;	  
          p.template scalar<2>() = startPointAlpha;

          // Check keep integrating criteria. Move particle if fulfilled, flag if not.
          if (startPointAlpha <= 0.0) {
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

  const size_t particlesAfter = amrParticles.getNumberOfValidParticesGlobal();
  if (particlesBefore != particlesAfter) {
    MayDay::Warning("StreamerInceptionStepper::K_moveParticlesTrapezoidal -- logic bust, lost some particles");
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::K_addIntegrationTrapezoidal(ParticleContainer<P>& amrParticles, const Real& voltage)
{
  CH_TIME("StreamerInceptionStepper::K_addIntegrationTrapezoidal");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::K_addIntegrationTrapezoidal" << endl;
  }

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
StreamerInceptionStepper<P, F>::K_moveParticlesSimpsons(ParticleContainer<P>& amrParticles, const Real& voltage)
{
  CH_TIME("StreamerInceptionStepper::K_moveParticlesSimpsons");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::K_moveParticlesSimpsons" << endl;
  }
  
  const size_t particlesBefore = amrParticles.getNumberOfValidParticesGlobal();

  const RealVect probLo = m_amr->getProbLo();
  const RealVect probHi = m_amr->getProbHi();
  
  // Move particles dx/2 (1/2)
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
          const Real      startPointAlpha       = m_alpha(Emagnitude);
	  
          const RealVect newPosition = p.position() + normalizedVelocity * 0.5 * dx;
	  
          p.template scalar<1>() = dx;
          p.template scalar<2>() = startPointAlpha;

	  // Check keep integrating criteria. Move particle if fulfilled, flag if not.
	  if (startPointAlpha <= 0.0) {
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

  // Move particles dx/2 (2/2)
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
        if (p.template scalar<0>() > 0) {
          const RealVect E                     = p.velocity();
          const RealVect  normalizedVelocity    = E / E.vectorLength();
          const Real      Emagnitude            = voltage * E.vectorLength();
          const Real      midPointAlpha         = m_alpha(Emagnitude);
	  
          p.template scalar<3>() = midPointAlpha;

          const RealVect newPosition = p.position() + normalizedVelocity * 0.5 * dx;

	  // Check keep integrating criteria. Move particle if fulfilled, flag if not.
	  if (midPointAlpha <= 0.0) {
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

  const size_t particlesAfter = amrParticles.getNumberOfValidParticesGlobal();
  if (particlesBefore != particlesAfter) {
    MayDay::Warning("StreamerInceptionStepper::K_moveParticlesSimpsons -- logic bust, lost some particles");
  }
}


template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::K_addIntegrationSimpsons(ParticleContainer<P>& amrParticles, const Real& voltage)
{
  CH_TIME("StreamerInceptionStepper::K_addIntegrationSimpsons");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::K_addIntegrationSimpsons" << endl;
  }

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
bool
StreamerInceptionStepper<P, F>::keepGoing(ParticleContainer<P>& amrParticles) const
{
  CH_TIME("StreamerInceptionStepper::keepGoing(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::keepGoing(ParticleContainer)" << endl;
  }

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
StreamerInceptionStepper<P, F>::moveParticlesInitialPosition(ParticleContainer<P>& amrParticles)
{
  CH_TIME("StreamerInceptionStepper::moveParticlesInitialPosition");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::moveParticlesInitialPosition" << endl;
  }

  // DELETE
  int countP = 0;
  
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(amrParticles.getRealm())[lvl];

    ParticleData<P>& levelParticles = amrParticles[lvl];
    for (DataIterator dit(dbl); dit.ok(); ++dit) {

      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        p.position() = p.template vect<0>();

	countP++; // DELETE
      }
    }
  }

  // DELETE
  std::cout << "Number of particles: " << countP << "\n\n";
  
  m_particleSolver->remap();
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::zeroNegativeKValues(EBAMRCellData& K, ParticleContainer<P>& amrParticles)
{
  CH_TIME("StreamerInceptionStepper::zeroNegativeKValues");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::zeroNegativeKValues" << endl;
  }

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
StreamerInceptionStepper<P, F>::resetParticleMass(ParticleContainer<P>& amrParticles)
{
  CH_TIME("StreamerInceptionStepper::resetParticleMass");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::resetParticleMass" << endl;
  }
  
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

// Add different integration methods after Simpson's fix.
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
  
  // Compute drift time integral twice - once for +dx and once for -dx direction.
  for (int i = 0; i < 2; i++) {    
    // Update velocities and starting mass to zero.
    m_particleSolver->interpolateVelocities();  
    
    ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();    
    
    this->resetParticleMass(amrParticles);  

    // Compute drift time integrals with input-based algorithm. 
    switch (m_integrationAlgorithm) {
    case IntegrationAlgorithm::Trapezoidal: {
      while (this->keepGoing(amrParticles)) {
	// Check integration flag and move particles if true
	this->driftTime_moveParticlesTrapezoidal(amrParticles, i);

	// Add integration parts to p.mass()
	this->driftTime_addIntegrationTrapezoidal(amrParticles);
      }

      break;
    }
    case IntegrationAlgorithm::Simpsons: {
      while (this->keepGoing(amrParticles)) {
	// Check integration flag and move particles if true	
	this->driftTime_moveParticlesSimpsons(amrParticles, i);
	
	// Add integration parts to p.mass()
	this->driftTime_addIntegrationSimpsons(amrParticles);
      }

      break;
    }
    default: {
      MayDay::Error("StreamerInceptionStepper::solveDriftTimes -- logic bust");
    }
    }
    
    // Move particles back to initial position. This includes remapping.
    this->moveParticlesInitialPosition(amrParticles);
  
    // Deposit and copy t values to m_driftTimes.
    m_particleSolver->deposit(t);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      t[lvl]->copyTo(Interval(0,0), *m_driftTimes[lvl], Interval(i,i));
    }     
  }    
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::driftTime_moveParticlesTrapezoidal(ParticleContainer<P>& amrParticles, const int dxSign) {
  CH_TIME("StreamerInceptionStepper::driftTime_moveParticlesTrapezoidal");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::driftTime_moveParticlesTrapezoidal" << endl;
  }

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
	  
	  // Find new position in direction decided by dxSign
	  RealVect newPosition;	  
	  if (dxSign == 0) {
	    newPosition = p.position() + normalizedVelocity * dx;
	  }
	  else if (dxSign == 1) {
	    newPosition = p.position() + normalizedVelocity * (-dx);
	  }	  
	  else {
	    MayDay::Error("StreamerInceptionStepper::driftTime_moveParticlesTrapezoidal -- logic bust");
	  }	 	  
	  	  
          p.template scalar<2>() = Emagnitude;
          p.template scalar<1>() = dx;	  
	  
	  // Check if particle outside grid or inside electrode. If not, move to new position.
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
}


template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::driftTime_addIntegrationTrapezoidal(ParticleContainer<P>& amrParticles){
  CH_TIME("StreamerInceptionStepper:driftTime_addIntegrationTrapezoidal");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::driftTime_addIntegrationTrapezoidal" << endl;
  }

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
StreamerInceptionStepper<P, F>::driftTime_moveParticlesSimpsons(ParticleContainer<P>& amrParticles, const int dxSign) {
  CH_TIME("StreamerInceptionStepper::driftTime_moveparticlesSimpsons");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::driftTime_moveparticlesSimpsons" << endl;
  }

  const RealVect probLo = m_amr->getProbLo();
  const RealVect probHi = m_amr->getProbHi();  

  // Move particles dx/2 (1/2)
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
	  
	  // Find new position in direction decided by dxSign
	  RealVect newPosition;	  
	  if (dxSign == 0) {
	    newPosition = p.position() + normalizedVelocity * 0.5 * dx;
	  }
	  else if (dxSign == 1) {
	    newPosition = p.position() + normalizedVelocity * 0.5 * (-dx);
	  }	  
	  else {
	    MayDay::Error("StreamerInceptionStepper::driftTime_moveParticlesTrapezoidal -- logic bust");
	  }	 	  

          p.template scalar<1>() = dx;	  	  
          p.template scalar<2>() = Emagnitude;
	  
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

  // Move particles dx/2 (2/2)
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
	  
	  // Find new position in direction decided by dxSign
	  RealVect newPosition;	  
	  if (dxSign == 0) {
	    newPosition = p.position() + normalizedVelocity * 0.5 * dx;
	  }
	  else if (dxSign == 1) {
	    newPosition = p.position() + normalizedVelocity * 0.5 * (-dx);
	  }	  
	  else {
	    MayDay::Error("StreamerInceptionStepper::driftTime_moveParticlesTrapezoidal -- logic bust");
	  }	 	  

          p.template scalar<3>() = Emagnitude;
	  
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

  // CHECK NUM PARTICLES BEFORE == PARTICLES AFTER?
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::driftTime_addIntegrationSimpsons(ParticleContainer<P>& amrParticles) {
  CH_TIME("StreamerInceptionStepper:driftTime_addIntegrationSimpsons");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::driftTime_addIntegrationSimpsons" << endl;
  }

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




/* ISSUES: 
   - V_inception on edge of electrode (irregular). Remove? This is also true for K-values.
   - Visualization: logarithmic scale preferable (if whole volume's values shall be shown), but not allowed by Visit because of 0 values (e.g. electrode)
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
      
      EBCellFAB& patchData = levelData[dit()];  
      BaseFab<Real>& regularData = patchData.getSingleValuedFAB();

      List<P>& particles = levelParticles[dit()].listItems();

      auto regularKernel = [&](const IntVect& iv) -> void {
	Real V_inception = 0.0;
	bool addP = true;
       
	// If largest K < 0.1, don't calculate V_inception and don't draw particle
	if (abs(regularData(iv, m_voltageSweeps.size()-1)) < 0.1) {
	  addP = false;
	}
	
	// If inception smaller than smallest K value, use slope decided by the two lowest K values
	else if (regularData(iv, 0) >= m_KInception) {
	  const Real slope = ( regularData(iv, 1) - regularData(iv, 0) ) / ( m_voltageSweeps[1] - m_voltageSweeps[0] );
	  const Real K_0 = regularData(iv, 0) - slope * m_voltageSweeps[0];
	  V_inception = ( m_KInception - K_0 ) / slope;
	}

	// If inception greater than greatest K value, use slope decided by the two largest K values
	else if (regularData(iv, m_voltageSweeps.size() - 1) <= m_KInception) {
	  const Real slope = ( regularData(iv, m_voltageSweeps.size() - 1) - regularData(iv, m_voltageSweeps.size() - 2) ) / ( m_voltageSweeps[m_voltageSweeps.size() - 1] - m_voltageSweeps[m_voltageSweeps.size() - 2] );
	  const Real K_0 = regularData(iv, m_voltageSweeps.size() - 1) - slope * m_voltageSweeps[m_voltageSweeps.size() - 1];
	  V_inception = ( m_KInception - K_0 ) / slope;

	  if (V_inception < 0) {
	    std::cout << "regularKernel greater.\t" << regularData(iv, m_voltageSweeps.size()-1) << "\t" << regularData(iv, m_voltageSweeps.size()-2) << "\t" << m_voltageSweeps[m_voltageSweeps.size()-1] << "\t" << m_voltageSweeps[m_voltageSweeps.size()-2] << "\n\n";
	  }	  
	}

	// If inception within voltage interval, iterate through values
	else {	  
	  for (int val = 1; val < m_voltageSweeps.size(); val++) {
	    if (regularData(iv, val) >= m_KInception) {
	      const Real slope = ( regularData(iv, val) - regularData(iv, val - 1) ) / ( m_voltageSweeps[val] - m_voltageSweeps[val - 1] );
	      const Real K_0 = regularData(iv, val - 1) - slope * m_voltageSweeps[val - 1];
	      V_inception = ( m_KInception - K_0 ) / slope;

	      if (V_inception < 0) {
		std::cout << "regularKernel within.\t" << regularData(iv, val) << "\t" << regularData(iv, val-1) << "\t" << m_voltageSweeps[val] << "\t" << m_voltageSweeps[val-1] << "\n\n";
	      }
	      break;
	    }
	  }
	}

	// Add particle if addP
	if (addP && validCells(iv, 0) && ebisbox.isRegular(iv)) {
	  P p;
	  p.position() = m_amr->getProbLo() + (0.5 * RealVect::Unit + RealVect(iv)) * dx;
	  p.mass() = V_inception;
	  particles.add(p);

	  // Warning if V_inception<1 - there's probably something wrong, maybe too few voltage_steps	  
	  if (V_inception < 1) {
	    std::cout << "WARNING: V_inception < 1" << endl;
	  }	 
	}	        
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
	Real V_inception = 0.0;
	const IntVect iv = vof.gridIndex(); 
	bool addP = true;
	
	// If largest K < 0.1, don't calculate V_inception and don't draw particle	
	if (abs(regularData(iv, m_voltageSweeps.size()-1)) < 0.1) {	  
	}
	
	// If inception smaller than smallest K value, use slope decided by the two lowest K values
	else if (regularData(iv, 0) >= m_KInception) {
	  const Real slope = ( regularData(iv, 1) - regularData(iv, 0) ) / ( m_voltageSweeps[1] - m_voltageSweeps[0] );
	  const Real K_0 = regularData(iv, 0) - slope * m_voltageSweeps[0];
	  V_inception = ( m_KInception - K_0 ) / slope;
	}

	// If inception greater than greatest K value, use slope decided by the two largest K values
	else if (regularData(iv, m_voltageSweeps.size() - 1) <= m_KInception) {
	  const Real slope = ( regularData(iv, m_voltageSweeps.size() - 1) - regularData(iv, m_voltageSweeps.size() - 2) ) / ( m_voltageSweeps[m_voltageSweeps.size() - 1] - m_voltageSweeps[m_voltageSweeps.size() - 2] );
	  const Real K_0 = regularData(iv, m_voltageSweeps.size() - 1) - slope * m_voltageSweeps[m_voltageSweeps.size() - 1];
	  V_inception = ( m_KInception - K_0 ) / slope;

	  if (V_inception < 0) {
	    std::cout << "irregularKernel greater.\t" << regularData(iv, m_voltageSweeps.size()-1) << "\t" << regularData(iv, m_voltageSweeps.size()-2) << "\t" << m_voltageSweeps[m_voltageSweeps.size()-1] << "\t" << m_voltageSweeps[m_voltageSweeps.size()-2] << "\n\n";
	  }	  
	}

	// If inception within voltage interval, iterate through values
	else {	  
	  for (int val = 1; val < m_voltageSweeps.size(); val++) {
	    if (regularData(iv, val) >= m_KInception) {
	      const Real slope = ( regularData(iv, val) - regularData(iv, val - 1) ) / ( m_voltageSweeps[val] - m_voltageSweeps[val - 1] );
	      const Real K_0 = regularData(iv, val - 1) - slope * m_voltageSweeps[val - 1];
	      V_inception = ( m_KInception - K_0 ) / slope;

	      if (V_inception < 0) {
		std::cout << "irregularKernel within.\t" << regularData(iv,val) << "\t" << regularData(iv, val-1) << "\t" << m_voltageSweeps[val] << "\t" << m_voltageSweeps[val-1] << "\n\n";
	      }
	      break;
	    }
	  }
	}

	// Add particle if addP
	if (addP && validCells(vof.gridIndex())) {	  
	  P p;
	  p.position() = m_amr->getProbLo() + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);
	  p.mass() = V_inception;
	  particles.add(p);

	  // Warning if V_inception<1 - there's probably something wrong, maybe too few voltage_steps
	  if (V_inception < 1) {
	    std::cout << "WARNING: V_inception < 1" << endl;
	  }	  
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
bool
StreamerInceptionStepper<P, F>::particleOutsideGrid(const RealVect pos, const RealVect probLo, const RealVect probHi)
{
  CH_TIME("StreamerInceptionStepper::particleOutsideGrid");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::particleOutsideGrid" << endl;
  }
  
  bool isOutside = false;

  for (int dir = 0; dir < SpaceDim; dir++) {
    if (pos[dir] < probLo[dir] || pos[dir] >= probHi[dir]) {
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
