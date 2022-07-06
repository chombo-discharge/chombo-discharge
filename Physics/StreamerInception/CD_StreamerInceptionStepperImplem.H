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

  this->parseOptions();

  m_verbosity = -1;

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
  m_amr->allocate(m_Kvalues, m_realm, m_phase, m_voltageSweeps.size());
  
  DataOps::setValue(m_Kvalues, 0.0);  
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
        const RealVect pos = m_amr->getProbLo() + (0.5 * RealVect::Unit + RealVect(iv)) * dx;

        if (validCells(iv, 0) && ebisbox.isRegular(iv)) {
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

      // Define particles' masses (i.e. integration tracker)
      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        p.mass()               = 0.0;
        p.template vect<0>()   = p.position();
        p.template scalar<0>() = 1.0;
      }
    }
  }

  // Remapping
  m_particleSolver->remap();

  // // Remove particles inside the EB
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

  // Something transient storage we can deposit into. 
  EBAMRCellData K;
  m_amr->allocate(K, m_realm, m_phase, 1);

  // Integrate for each voltage
  for (int i = 0; i < m_voltageSweeps.size(); ++i) {

    // Update velocities and starting mass to zero.
    m_particleSolver->interpolateVelocities();
    
    ParticleContainer<P>& amrParticles = m_particleSolver->getParticles();    

    this->resetParticleMass(amrParticles);

    switch (m_integrationAlgorithm) {
    case IntegrationAlgorithm::Trapezoidal: {
      while (this->keepIntegrating(amrParticles)) {
        // Check contIntegral criteria and move particles if true
        this->moveParticlesTrapezoidal(amrParticles, m_voltageSweeps[i]);

        // Add integration parts to p.mass()
        this->addIntegrationTrapezoidal(amrParticles, m_voltageSweeps[i]);
      }

      break;
    }
    case IntegrationAlgorithm::Simpsons: {
      while (this->keepIntegrating(amrParticles)) {
        // Check contIntegral criteria and move particles if true
        this->moveParticlesSimpsons(amrParticles, m_voltageSweeps[i]);

        // Add integration parts to p.mass()
        this->addIntegrationSimpsons(amrParticles, m_voltageSweeps[i]);
      }

      break;
    }
    default: {
      MayDay::Error("StreamerInceptionStepper::postInitialize -- logic bust");
    }
    }

    // Move particles back to initial position. This includes remapping.
    this->moveParticlesInitialPosition(amrParticles);

    // Deposit and copy K values to m_Kvalues
    m_particleSolver->deposit(K);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      K[lvl]->copyTo(Interval(0,0), *m_Kvalues[lvl], Interval(i,i));
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::moveParticlesTrapezoidal(ParticleContainer<P>& amrParticles, const Real& voltage)
{
  CH_TIME("StreamerInceptionStepper::moveParticlesTrapezoidal");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::moveParticlesTrapezoidal" << endl;
  }

  const size_t particlesBefore = amrParticles.getNumberOfValidParticesGlobal();

  const RealVect probLo = m_amr->getProbLo();
  const RealVect probHi = m_amr->getProbHi();

  // Set distance moved per step dx
  //  const Real dx = 0.1;  //m_amr->getDx()[m_amr->getFinestLevel()];
  //  const Real dx = 0.1 * m_amr->getDx()[m_amr->getFinestLevel()];

  m_particleSolver->interpolateVelocities();

  // Move particles
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(amrParticles.getRealm())[lvl];

    ParticleData<P>& levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {

      // Particles move half a cell on this level. 
      const Real dx = 0.5*m_amr->getDx()[lvl];

      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        // Move particle and check if contIntegral if contIntegral not already false
        if (p.template scalar<0>() > 0.0) {

          const RealVect E                     = p.velocity();
          const Real     particleDistanceMoved = dx;
          const RealVect normalizedVelocity    = E / E.vectorLength();
          const Real     Emagnitude            = voltage * E.vectorLength();
          const Real     startPointAlpha       = m_alpha(Emagnitude);

          const RealVect newPosition = p.position() + normalizedVelocity * particleDistanceMoved;

          p.template scalar<2>() = startPointAlpha;
          p.template scalar<1>() = particleDistanceMoved;

          // Check if alpha<0
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
    MayDay::Warning("StreamerInceptionStepper::moveParticlesTrapezoidal -- logic bust, lost some particles");
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::moveParticlesSimpsons(ParticleContainer<P>& amrParticles, const Real& voltage)
{
  CH_TIME("StreamerInceptionStepper::moveParticlesSimpsons");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::moveParticlesSimpsons" << endl;
  }
  
  const size_t particlesBefore = amrParticles.getNumberOfValidParticesGlobal();

  const RealVect probLo = m_amr->getProbLo();
  const RealVect probHi = m_amr->getProbHi();

  // Set distance moved per step dx
  const Real dx = 0.1; //m_amr->getDx()[m_amr->getFinestLevel()];

  // Move particles dx/2 (1/2)
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const EBISLayout&        ebisl          = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const DisjointBoxLayout& dbl            = m_amr->getGrids(amrParticles.getRealm())[lvl];
    ParticleData<P>&         levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        // Move particle and check if contIntegral if contIntegral not already false
        if (p.template scalar<0>() > 0) {

          const RealVect& E                     = p.velocity();
          const Real      particleDistanceMoved = dx;
          const RealVect  normalizedVelocity    = E / E.vectorLength();
          const Real      Emagnitude            = voltage * E.vectorLength();
          const Real      startPointAlpha       = m_alpha(Emagnitude);

          // Check if alpha<0
          if (startPointAlpha < 0) {
            p.template scalar<0>() = -1;
            std::cout << "alpha<0" << endl;
            continue;
          }

          p.template scalar<1>() = particleDistanceMoved;
          p.template scalar<2>() = startPointAlpha;

          const RealVect newPosition = p.position() + normalizedVelocity * 0.5 * particleDistanceMoved;

          // Check if particle is outside of grid boundaries
          if (particleOutsideGrid(newPosition, probLo, probHi)) {
            p.template scalar<0>() = -1;
            continue;
          }

          // Check if particle is inside electrode
          if (particleInsideElectrode(newPosition)) {
            p.template scalar<0>() = -1;
            continue;
          }

          // Move particle if contIntegral
          p.position() = newPosition;
        }
      }
    }
  }
  m_particleSolver->remap();
  m_particleSolver->interpolateVelocities();

  // Move particles dx/2 (2/2)
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const EBISLayout&        ebisl          = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const DisjointBoxLayout& dbl            = m_amr->getGrids(amrParticles.getRealm())[lvl];
    ParticleData<P>&         levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        // Move particle and check if contIntegral if contIntegral not already false
        if (p.template scalar<0>() > 0) {

          const RealVect& E                     = p.velocity();
          const Real      particleDistanceMoved = dx;
          const RealVect  normalizedVelocity    = E / E.vectorLength();
          const Real      Emagnitude            = voltage * E.vectorLength();
          const Real      midPointAlpha         = m_alpha(Emagnitude);

          // Check if alpha<0
          if (midPointAlpha < 0) {
            p.template scalar<0>() = -1;
            std::cout << "alpha<0" << endl;
            continue;
          }

          p.template scalar<1>() = particleDistanceMoved;
          p.template scalar<3>() = midPointAlpha;

          const RealVect newPosition = p.position() + normalizedVelocity * 0.5 * particleDistanceMoved;

          // Check if particle is outside of grid boundaries
          if (particleOutsideGrid(newPosition, probLo, probHi)) {
            p.template scalar<0>() = -1;
            continue;
          }

          // Check if particle is inside electrode
          if (particleInsideElectrode(newPosition)) {
            p.template scalar<0>() = -1;
            continue;
          }

          // Move particle if contIntegral
          p.position() = newPosition;
        }
      }
    }
  }
  m_particleSolver->remap();
  m_particleSolver->interpolateVelocities();

  const size_t particlesAfter = amrParticles.getNumberOfValidParticesGlobal();
  if (particlesBefore != particlesAfter) {
    MayDay::Warning("StreamerInceptionStepper::moveParticlesSimpsons -- logic bust, lost some particles");
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::addIntegrationTrapezoidal(ParticleContainer<P>& amrParticles, const Real& voltage)
{
  CH_TIME("StreamerInceptionStepper::addIntegrationTrapezoidal");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::addIntegrationTrapezoidal" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(amrParticles.getRealm())[lvl];

    ParticleData<P>& levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        // Add integration parts if contIntegral true
        if (p.template scalar<0>() > 0.0) {
          const RealVect& E               = p.velocity();
          const Real      Emagnitude      = voltage * E.vectorLength();
          const Real      endPointAlpha   = m_alpha(Emagnitude);
          const Real      startPointAlpha = p.template scalar<2>();
          const Real      distanceMoved   = p.template scalar<1>();
          p.mass() += 0.5 * (endPointAlpha + startPointAlpha) * distanceMoved;
        }
      }
    }
  }
}

template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::addIntegrationSimpsons(ParticleContainer<P>& amrParticles, const Real& voltage)
{
  CH_TIME("StreamerInceptionStepper::addIntegrationSimpsons");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::addIntegrationSimpsons" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(amrParticles.getRealm())[lvl];

    ParticleData<P>& levelParticles = amrParticles[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<P>& particles = levelParticles[dit()].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        P& p = lit();

        // Add integration parts if contIntegral true
        if (p.template scalar<0>() > 0) {
          const RealVect& E               = p.velocity();
          const Real      Emagnitude      = voltage * E.vectorLength();
          const Real      endPointAlpha   = m_alpha(Emagnitude);
          const Real      distanceMoved   = p.template scalar<1>();
          const Real      startPointAlpha = p.template scalar<2>();
          const Real      midPointAlpha   = p.template scalar<3>();
          p.mass() += (distanceMoved / 6) * (startPointAlpha + 4 * midPointAlpha + endPointAlpha);
        }
      }
    }
  }
}

template <typename P, typename F>
bool
StreamerInceptionStepper<P, F>::keepIntegrating(ParticleContainer<P>& amrParticles) const
{
  CH_TIME("StreamerInceptionStepper::keepIntegrating(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::keepIntegrating(ParticleContainer)" << endl;
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

  amrParticles.remap();
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

// Done in the constructor
template <typename P, typename F>
void
StreamerInceptionStepper<P, F>::parseOptions()
{
  CH_TIME("StreamerInceptionStepper::parseOptions");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::parseOptions" << endl;
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

  // Write our K-values.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Interval srcInterval(0, m_voltageSweeps.size() - 1);
    const Interval dstInterval(a_icomp, a_icomp + m_voltageSweeps.size() - 1);
    m_Kvalues[lvl]->copyTo(srcInterval, *a_output[lvl], dstInterval);
  }
  a_icomp += m_voltageSweeps.size();

  std::string varName;
  for (const Real& V : m_voltageSweeps) {
    varName = "K for ";
    varName += std::to_string(V);
    varName += " V";
    a_plotVariableNames.push_back(varName);
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

  ncomp += m_voltageSweeps.size();

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
StreamerInceptionStepper<P, F>::printStepReport()
{
  CH_TIME("StreamerInceptionStepper::printStepReport");
  if (m_verbosity > 5) {
    pout() << "StreamerInceptionStepper::printStepReport" << endl;
  }
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

#include <CD_NamespaceFooter.H>

#endif
