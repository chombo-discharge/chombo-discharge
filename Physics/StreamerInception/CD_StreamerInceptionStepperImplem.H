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

// Tasks:
//
// 1. Compute inception integral with particle tracing.
// 2. Compute critical volume.
// 3. Compute negative ion drift times to electrode(s). Compute as time/(V * mobility).
// 4. Background ionization rates
// 5. Evaluate field emission rates.
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

  // Default alpha.
  m_alpha = [](const Real E) {
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

  // clang-format off
  m_amr->allocate(m_velocity,          m_realm, m_phase, SpaceDim);
  m_amr->allocate(m_inceptionIntegral, m_realm, m_phase, m_voltageSweeps.size());
  m_amr->allocate(m_inceptionVoltage,  m_realm, m_phase, 1);
  m_amr->allocate(m_driftTimes,        m_realm, m_phase, 2);
  // clang-format on

  DataOps::setValue(m_inceptionIntegral, 0.0);
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
  DataOps::setValue(m_velocity, 0.0);
  DataOps::copy(m_velocity, m_amr->alias(phase::gas, m_fieldSolver->getElectricField()));
  m_amr->averageDown(m_velocity, m_realm, m_phase);
  m_amr->interpGhost(m_velocity, m_realm, m_phase);
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
  int  voltageSteps = 0;

  std::string str;

  // Read in parameters
  pp.get("verbosity", m_verbosity);
  pp.get("voltage_lo", voltageLo);
  pp.get("voltage_hi", voltageHi);
  pp.get("voltage_steps", voltageSteps);
  pp.get("K_inception", m_inceptionK);
  pp.get("integration", str);
  pp.get("print_report", m_printReport);

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
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Interval srcInterval(0, m_voltageSweeps.size() - 1);
    const Interval dstInterval(a_icomp, a_icomp + m_voltageSweeps.size() - 1);
    m_inceptionIntegral[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
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

  // Write the inception voltage.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Interval srcInterval(0, 0);
    const Interval dstInterval(a_icomp, a_icomp);
    m_inceptionVoltage[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
  }
  a_icomp += 1;

  a_plotVariableNames.push_back("Inception voltage");
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
  ncomp += m_voltageSweeps.size();

  // Drift times (+/-) direction
  ncomp += 2;

  // Inception voltage.
  ncomp += 1;

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
StreamerInceptionStepper<P, F>::setDetachmentRate(const std::function<Real(const Real& E)>& a_detachmentRate) noexcept
{
  CH_TIME("StreamerInceptionStepper::setDetachmentRate");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::setDetachmentRate" << endl;
  }

  m_detachmentRate = a_detachmentRate;
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
  m_particleSolver->setVelocity(m_velocity);

  this->initialParticles();

  // Solve K integral for each particle in each voltage and store K values in m_inceptionIntegral
  this->computeInceptionIntegral();

  // Solve drift time integral in both directions and store values in m_driftTimes
  this->computeDriftTimes();

  // Interpolate between K values to find voltage giving K_inception and store values in m_inceptionVoltage
  this->computeInceptionVoltageVolume();

  // Compute critical volume of K values for each voltage
  this->computeCriticalVolume();

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

  // Transient storage we can deposit particles onto.
  EBAMRCellData K;
  m_amr->allocate(K, m_realm, m_phase, 1);

  // Integrate for each voltage
  for (int i = 0; i < m_voltageSweeps.size(); ++i) {

    // Update velocities and starting mass to zero.
    m_particleSolver->interpolateVelocities();

    ParticleContainer<P>& particles = m_particleSolver->getParticles();

    this->resetParticleMass();

    switch (m_integrationAlgorithm) {
    case IntegrationAlgorithm::Trapezoidal: {
      while (this->keepIntegrating()) {

        // Move particles and add integration contribution for the current step.
        this->K_moveParticlesTrapezoidal(m_voltageSweeps[i]);
        this->K_addIntegrationTrapezoidal(m_voltageSweeps[i]);
      }

      break;
    }
    case IntegrationAlgorithm::Simpsons: {
      while (this->keepIntegrating()) {

        // Move particles and add integration contribution for the current step.
        this->K_moveParticlesSimpsons(m_voltageSweeps[i]);
        this->K_addIntegrationSimpsons(m_voltageSweeps[i]);
      }

      break;
    }
    default: {
      MayDay::Error("StreamerInceptionStepper::computeInceptionIntegral -- logic bust");

      break;
    }
    }

    // Move particles back to initial position and deposit them on the mesh.
    this->rewindParticles();

    m_particleSolver->deposit(K);

    // Get max K value.
    Real maxK = -std::numeric_limits<Real>::max();
    Real minK = std::numeric_limits<Real>::max();

    DataOps::getMaxMin(maxK, minK, K, 0);

    m_maxKValues.push_back(maxK);

    // Copy K values to m_inceptionIntegral so we can plot them later.
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      K[lvl]->copyTo(Interval(0, 0), *m_inceptionIntegral[lvl], Interval(i, i));
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::K_moveParticlesTrapezoidal(const Real& a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::K_moveParticlesTrapezoidal");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::K_moveParticlesTrapezoidal" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

#ifndef NDEBUG
  const size_t particlesBefore = amrParticles.getNumberOfValidParticesGlobal();
#endif

  // Move particles
  this->K_moveParticlesOneStep(a_voltage, 0);

#ifndef NDEBUG
  const size_t particlesAfter = amrParticles.getNumberOfValidParticesGlobal();

  if (particlesBefore != particlesAfter) {
    MayDay::Warning("StreamerInceptionStepper::K_moveParticlesTrapezoidal -- logic bust, lost some particles");
  }
#endif
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::K_addIntegrationTrapezoidal(const Real& a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::K_addIntegrationTrapezoidal");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::K_addIntegrationTrapezoidal" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<P>& particles = amrParticles[lvl][dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        // Add integration parts if integration flag true
        if (p.template scalar<0>() > 0.0) {
          const RealVect E               = p.velocity();
          const Real     Emagnitude      = a_voltage * E.vectorLength();
          const Real     endPointAlpha   = m_alpha(Emagnitude);
          const Real     dx              = p.template scalar<1>();
          const Real     startPointAlpha = p.template scalar<2>();

          if (startPointAlpha >= 0.0 && endPointAlpha >= 0.0) {
            p.mass() += 0.5 * (endPointAlpha + startPointAlpha) * dx;
          }
          else {
            p.template scalar<0>() = -1.0;
          }
        }
      }
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::K_moveParticlesSimpsons(const Real& a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::K_moveParticlesSimpsons");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::K_moveParticlesSimpsons" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

#ifndef NDEBUG
  const size_t particlesBefore = amrParticles.getNumberOfValidParticesGlobal();
#endif

  // Step 1: Move particles dx/2 and remap
  this->K_moveParticlesOneStep(a_voltage, 1);

#ifndef NDEBUG
  const size_t particlesMiddle = amrParticles.getNumberOfValidParticesGlobal();
#endif

  // Step 2: Move particles dx/2 and remap
  this->K_moveParticlesOneStep(a_voltage, 2);

#ifndef NDEBUG
  const size_t particlesAfter = amrParticles.getNumberOfValidParticesGlobal();
  if (particlesBefore != particlesAfter || particlesBefore != particlesMiddle) {
    MayDay::Warning("StreamerInceptionStepper::K_moveParticlesSimpsons -- logic bust, lost some particles");
  }
#endif
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::K_addIntegrationSimpsons(const Real& a_voltage) noexcept
{
  CH_TIME("StreamerInceptionStepper::K_addIntegrationSimpsons");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::K_addIntegrationSimpsons" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    ParticleData<P>& levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        // Add integration parts if integration flag true
        if (p.template scalar<0>() > 0) {
          const RealVect E               = p.velocity();
          const Real     Emagnitude      = a_voltage * E.vectorLength();
          const Real     endPointAlpha   = m_alpha(Emagnitude);
          const Real     dx              = p.template scalar<1>();
          const Real     startPointAlpha = p.template scalar<2>();
          const Real     midPointAlpha   = p.template scalar<3>();

          if (startPointAlpha >= 0.0 && midPointAlpha >= 0.0 && endPointAlpha >= 0.0) {
            p.mass() += (dx / 6) * (startPointAlpha + 4 * midPointAlpha + endPointAlpha);
          }
          else {
            p.template scalar<0>() = -1.0;
          }
        }
      }
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::K_moveParticlesOneStep(const Real& a_voltage, const int& a_step) noexcept
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
    const Real dx = 0.5 * m_amr->getDx()[lvl];

    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    ParticleData<P>& levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        // If integration flag true, move particle and check integration conditions
        if (p.template scalar<0>() > 0.0) {
          const RealVect E                  = p.velocity();
          const RealVect normalizedVelocity = E / E.vectorLength();
          const Real     Emagnitude         = a_voltage * E.vectorLength();
          const Real     currentAlpha       = m_alpha(Emagnitude);

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
            if (a_step == 1) {
              p.template scalar<1>() = dx;
              p.template scalar<2>() = currentAlpha;
            }
            else if (a_step == 2) {
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
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

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
StreamerInceptionStepper<P, F>::rewindParticles() noexcept
{
  CH_TIME("StreamerInceptionStepper::rewindParticles");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::rewindParticles" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

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
StreamerInceptionStepper<P, F>::resetParticleMass() noexcept
{
  CH_TIME("StreamerInceptionStepper::resetParticleMass");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::resetParticleMass" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl            = m_amr->getGrids(m_realm)[lvl];
    ParticleData<P>&         levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {

        // Reset mass and flag for keep integrating.
        P& p = lit();

        p.mass()               = 0.0;
        p.template scalar<0>() = 1.0;
      }
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::computeDriftTimes() noexcept
{
  CH_TIME("StreamerInceptionStepper::computeDriftTimes");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::computeDriftTimes" << endl;
  }

  // Something transient storage we can deposit into.
  EBAMRCellData t;
  m_amr->allocate(t, m_realm, m_phase, 1);

  // Compute drift time integral twice - once for +dx and once for -dx direction.
  for (int i = 0; i < 2; i++) {
    // Update velocities and starting mass to zero.
    m_particleSolver->interpolateVelocities();

    ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

    this->resetParticleMass();

    // Compute drift time integrals with input-based algorithm.
    switch (m_integrationAlgorithm) {
    case IntegrationAlgorithm::Trapezoidal: {
      while (this->keepIntegrating()) {

        this->driftTime_moveParticlesTrapezoidal(i);
        this->driftTime_addIntegrationTrapezoidal();
      }

      break;
    }
    case IntegrationAlgorithm::Simpsons: {
      while (this->keepIntegrating()) {

        this->driftTime_moveParticlesSimpsons(i);
        this->driftTime_addIntegrationSimpsons();
      }

      break;
    }
    default: {
      MayDay::Error("StreamerInceptionStepper::computeDriftTimes -- logic bust");
    }
    }

    // Rewind particles and deposit drift times onto the mesh.
    this->rewindParticles();

    m_particleSolver->deposit(t);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      t[lvl]->copyTo(Interval(0, 0), *m_driftTimes[lvl], Interval(i, i));
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::driftTime_moveParticlesTrapezoidal(const int& a_dxSign) noexcept
{
  CH_TIME("StreamerInceptionStepper::driftTime_moveParticlesTrapezoidal");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::driftTime_moveParticlesTrapezoidal" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  const size_t particlesBefore = amrParticles.getNumberOfValidParticesGlobal();

  // Move particles and remap
  this->driftTime_moveParticlesOneStep(a_dxSign, 0);

  const size_t particlesAfter = amrParticles.getNumberOfValidParticesGlobal();
  if (particlesBefore != particlesAfter) {
    MayDay::Warning("StreamerInceptionStepper::driftTime_moveParticlesTrapezoidal -- logic bust, lost some particles");
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::driftTime_addIntegrationTrapezoidal() noexcept
{
  CH_TIME("StreamerInceptionStepper:driftTime_addIntegrationTrapezoidal");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::driftTime_addIntegrationTrapezoidal" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    ParticleData<P>& levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        // Add integration parts if integration flag true
        if (p.template scalar<0>() > 0.0) {
          const RealVect E                    = p.velocity();
          const Real     Emagnitude           = E.vectorLength();
          const Real     dx                   = p.template scalar<1>();
          const Real     startPointEmagnitude = p.template scalar<2>();
          p.mass() += (dx / 2) * ((1 / startPointEmagnitude) + (1 / Emagnitude));
        }
      }
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::driftTime_moveParticlesSimpsons(const int& a_dxSign) noexcept
{
  CH_TIME("StreamerInceptionStepper::driftTime_moveparticlesSimpsons");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::driftTime_moveparticlesSimpsons" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

#ifndef NDEBUG
  const size_t particlesBefore = amrParticles.getNumberOfValidParticesGlobal();
#endif

  // Step 1: Move particles dx/2 and remap
  this->driftTime_moveParticlesOneStep(a_dxSign, 1);

  // Step 2: Move particles dx/2 and remap
  this->driftTime_moveParticlesOneStep(a_dxSign, 2);

#ifndef NDEBUG  
  const size_t particlesAfter = amrParticles.getNumberOfValidParticesGlobal();
  if (particlesBefore != particlesAfter) {
    MayDay::Warning("StreamerInceptionStepper::driftTime_moveParticlesSimpsons -- logic bust, lost some particles");
  }
#endif
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::driftTime_addIntegrationSimpsons() noexcept
{
  CH_TIME("StreamerInceptionStepper:driftTime_addIntegrationSimpsons");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::driftTime_addIntegrationSimpsons" << endl;
  }

  ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    ParticleData<P>& levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        // Add integration parts if integration flag true
        if (p.template scalar<0>() > 0.0) {
          const RealVect E                    = p.velocity();
          const Real     Emagnitude           = E.vectorLength();
          const Real     dx                   = p.template scalar<1>();
          const Real     startPointEmagnitude = p.template scalar<2>();
          const Real     midPointEmagnitude   = p.template scalar<3>();
	  
          p.mass() += (dx / 6) * ((1 / startPointEmagnitude) + 4 * (1 / midPointEmagnitude) + (1 / Emagnitude));
        }
      }
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::driftTime_moveParticlesOneStep(const int& a_dxSign, const int& a_step) noexcept
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
    const Real dx = 0.5 * m_amr->getDx()[lvl];

    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];

    ParticleData<P>& levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        // If integration flag true, move particle and check integration conditions
        if (p.template scalar<0>() > 0.0) {
          const RealVect E                  = p.velocity();
          const RealVect normalizedVelocity = E / E.vectorLength();
          const Real     Emagnitude         = E.vectorLength();

          // Calculate new position based on integration algorithm and direction
          RealVect newPosition;
          switch (m_integrationAlgorithm) {
          case IntegrationAlgorithm::Trapezoidal: {
            if (a_dxSign == 0) {
              newPosition = p.position() + normalizedVelocity * dx;
            }
            else if (a_dxSign == 1) {
              newPosition = p.position() + normalizedVelocity * (-dx);
            }
            else {
              MayDay::Error("StreamerInceptionStepper::driftTime_moveParticlesOneStep -- logic bust");
            }

            // Store values to particle for subsequent integration
            p.template scalar<1>() = dx;
            p.template scalar<2>() = Emagnitude;

            break;
          }
          case IntegrationAlgorithm::Simpsons: {
            if (a_dxSign == 0) {
              newPosition = p.position() + normalizedVelocity * 0.5 * dx;
            }
            else if (a_dxSign == 1) {
              newPosition = p.position() + normalizedVelocity * 0.5 * (-dx);
            }
            else {
              MayDay::Error("StreamerInceptionStepper::driftTime_moveParticlesOneStep -- logic bust");
            }

            // Store values to particle for subsequent integration
            if (a_step == 1) {
              p.template scalar<1>() = dx;
              p.template scalar<2>() = Emagnitude;
            }
            else if (a_step == 2) {
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
StreamerInceptionStepper<P, F>::printReport() const noexcept
{
  CH_TIME("StreamerInceptionStepper::printReport");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::printReport" << endl;
  }

  auto printElement = [](const Real el) -> void {
    std::cout << left << setw(15) << setfill(' ') << el;
  };

  // Compute the inception voltage.
#ifdef CH_MPI
  if (procID() == 0) {
#endif
    const Real Uinc = this->computeInceptionVoltage();

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
    std::cout << "# =============================================\n";
    std::cout << left << setw(15) << setfill(' ') << "# Voltage";
    std::cout << left << setw(15) << setfill(' ') << "Max K";
    std::cout << left << setw(15) << setfill(' ') << "Critical volume\n";
    std::cout << "# =============================================\n";

    for (int i = 0; i < m_voltageSweeps.size(); i++) {
      printElement(m_voltageSweeps[i]);
      printElement(m_maxKValues[i]);
      printElement(m_criticalVolumes[i]);
      std::cout << endl;
    }

    std::cout << "# =============================================\n";
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
StreamerInceptionStepper<P, F>::particleInsideElectrode(const RealVect& a_pos) const noexcept
{
#ifndef NDEBUG
  CH_TIME("StreamerInceptionStepper::particleInsideElectrode");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::particleInsideElectrode" << endl;
  }
#endif

  const RefCountedPtr<BaseIF>& implicitFunction = m_amr->getBaseImplicitFunction(m_phase);

  return (implicitFunction->value(a_pos) >= 0.0) ? true : false;
}

#include <CD_NamespaceFooter.H>

#endif
