/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_TracerParticleStepperImplem.H
  @brief  Implementation of CD_TracerParticleStepper.H
  @author Robert Marskar
*/

#ifndef CD_TRACERPARTICLESTEPPERIMPLEM_H
#define CD_TRACERPARTICLESTEPPERIMPLEM_H

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_ParticleManagement.H>
#include <CD_Random.H>
#include <CD_TracerParticleStepper.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::TracerParticle;

template <typename P>
inline TracerParticleStepper<P>::TracerParticleStepper()
{
  CH_TIME("TracerParticleStepper::TracerParticleStepper");

  m_realm = Realm::Primal;
  m_phase = phase::gas;

  this->parseOptions();
}

template <typename P>
inline TracerParticleStepper<P>::~TracerParticleStepper()
{
  CH_TIME("TracerParticleStepper::~TracerParticleStepper");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::~TracerParticleStepper" << endl;
  }
}

template <typename P>
inline void
TracerParticleStepper<P>::setupSolvers()
{
  CH_TIME("TracerParticleStepper::setupSolvers()");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::setupSolvers()" << endl;
  }

  m_solver = RefCountedPtr<TracerParticleSolver<P>>(new TracerParticleSolver<P>(m_amr, m_computationalGeometry));

  m_solver->setPhase(m_phase);
  m_solver->setRealm(m_realm);
  m_solver->parseOptions();
}

template <typename P>
inline void
TracerParticleStepper<P>::allocate()
{
  CH_TIME("TracerParticleStepper::allocate()");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::allocate()" << endl;
  }

  m_amr->allocate(m_velocity, m_realm, m_phase, SpaceDim);
  m_solver->allocate();
}

template <typename P>
inline void
TracerParticleStepper<P>::initialData()
{
  CH_TIME("TracerParticleStepper::initialData()");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::initialData()" << endl;
  }

  this->setVelocity();
  this->initialParticles();

  m_solver->setVelocity(m_velocity);
  m_solver->interpolateVelocities();
}

template <typename P>
inline void
TracerParticleStepper<P>::registerRealms()
{
  CH_TIME("TracerParticleStepper::registerRealms()");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::registerRealms()" << endl;
  }

  m_amr->registerRealm(m_realm);
}

template <typename P>
inline void
TracerParticleStepper<P>::registerOperators()
{
  CH_TIME("TracerParticleStepper::registerOperators()");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::registerOperators()" << endl;
  }

  m_solver->registerOperators();
}

template <typename P>
inline void
TracerParticleStepper<P>::parseOptions()
{
  CH_TIME("TracerParticleStepper::parseOptions()");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::parseOptions()" << endl;
  }

  this->parseIntegrator();
  this->parseVelocityField();
  this->parseInitialConditions();
}

template <typename P>
inline void
TracerParticleStepper<P>::parseRuntimeOptions()
{
  CH_TIME("TracerParticleStepper::parseRuntimeOptions()");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::parseRuntimeOptions()" << endl;
  }

  this->parseIntegrator();

  m_solver->parseRuntimeOptions();
}

template <typename P>
inline void
TracerParticleStepper<P>::parseIntegrator()
{
  CH_TIME("TracerParticleStepper::parseIntegrator()");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::parseIntegrator()" << endl;
  }

  ParmParse pp("TracerParticleStepper");

  std::string str;

  pp.get("verbosity", m_verbosity);
  pp.get("cfl", m_cfl);
  pp.get("integration", str);
  if (str == "euler") {
    m_algorithm = IntegrationAlgorithm::Euler;
  }
  else if (str == "rk2") {
    m_algorithm = IntegrationAlgorithm::RK2;
  }
  else if (str == "rk4") {
    m_algorithm = IntegrationAlgorithm::RK4;
  }
  else {
    MayDay::Error("TracerParticleStepper::parseIntegrator -- logic bust");
  }
}

template <typename P>
inline void
TracerParticleStepper<P>::parseVelocityField()
{
  CH_TIME("TracerParticleStepper::parseVelocityField()");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::parseVelocityField()" << endl;
  }

  ParmParse pp("TracerParticleStepper");

  int v;
  pp.get("velocity_field", v);

  if (v == 0) {
    m_velocityField = VelocityField::Diagonal;
  }
  else if (v == 1) {
    m_velocityField = VelocityField::Rotational;
  }
  else {
    MayDay::Error("TracerParticleStepper::parseVelocityField -- logic bust");
  }
}

template <typename P>
inline void
TracerParticleStepper<P>::parseInitialConditions()
{
  CH_TIME("TracerParticleStepper::parseInitialConditions()");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::parseInitialConditions()" << endl;
  }

  ParmParse pp("TracerParticleStepper");

  Real numParticles;
  pp.get("initial_particles", numParticles);

  m_numInitialParticles = size_t(std::max(0.0, numParticles));
}

#ifdef CH_USE_HDF5
template <typename P>
inline void
TracerParticleStepper<P>::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const
{
  CH_TIME("TracerParticleStepper::writeCheckpointData(HDF5Handle, int)");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::writeCheckpointData(HDF5Handle, int)" << endl;
  }

  m_solver->writeCheckpointLevel(a_handle, a_lvl);
}
#endif

#ifdef CH_USE_HDF5
template <typename P>
inline void
TracerParticleStepper<P>::readCheckpointData(HDF5Handle& a_handle, const int a_lvl)
{
  CH_TIME("TracerParticleStepper::readCheckpointData(HDF5Handle, int)");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::readCheckpointData(HDF5Handle, int)" << endl;
  }

  m_solver->readCheckpointLevel(a_handle, a_lvl);
}
#endif

template <typename P>
inline int
TracerParticleStepper<P>::getNumberOfPlotVariables() const
{
  CH_TIME("TracerParticleStepper::getNumberOfPlotVariables()");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::getNumberOfPlotVariables()" << endl;
  }

  return m_solver->getNumberOfPlotVariables();
}

template <typename P>
inline Vector<std::string>
TracerParticleStepper<P>::getPlotVariableNames() const
{
  CH_TIME("TracerParticleStepper::getPlotVariableNames()");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::getPlotVariableNames()" << endl;
  }

  return m_solver->getPlotVariableNames();
}

template <typename P>
inline void
TracerParticleStepper<P>::writePlotData(LevelData<EBCellFAB>& a_output,
                                        int&                  a_icomp,
                                        const std::string&    a_outputRealm,
                                        const int             a_level) const
{
  CH_TIME("TracerParticleStepper::writePlotData(EBAMRCellData, Vector<std::string>, int)");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::writePlotData(EBAMRCellData, Vector<std::string>, int)" << endl;
  }

  CH_assert(a_level >= 0);
  CH_assert(a_level <= m_amr->getFinestLevel());

  m_solver->writePlotData(a_output, a_icomp, a_outputRealm, a_level);
}

template <typename P>
inline Real
TracerParticleStepper<P>::computeDt()
{
  CH_TIME("TracerParticleStepper::computeDt()");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::computeDt()" << endl;
  }

  return m_cfl * m_solver->computeDt();
}

template <typename P>
inline Real
TracerParticleStepper<P>::advance(const Real a_dt)
{
  CH_TIME("TracerParticleStepper::advance(Real)");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::advance(Real)" << endl;
  }

  switch (m_algorithm) {
  case IntegrationAlgorithm::Euler: {
    this->advanceParticlesEuler(a_dt);

    break;
  }
  case IntegrationAlgorithm::RK2: {
    this->advanceParticlesRK2(a_dt);

    break;
  }
  case IntegrationAlgorithm::RK4: {
    this->advanceParticlesRK4(a_dt);

    break;
  }
  default: {
    MayDay::Error("TracerParticleStepper::advance -- logic bust");
  }
  }

  return a_dt;
}

template <typename P>
inline void
TracerParticleStepper<P>::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt)
{
  CH_TIME("TracerParticleStepper::synchronizeSolverTimes");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::synchronizeSolverTimes" << endl;
  }

  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;

  m_solver->setTime(a_step, a_time, a_dt);
}

template <typename P>
inline void
TracerParticleStepper<P>::preRegrid(const int a_lmin, const int a_oldFinestLevel)
{
  CH_TIME("TracerParticleStepper::preRegrid(int, int)");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::preRegrid(int, int)" << endl;
  }

  m_solver->preRegrid(a_lmin, a_oldFinestLevel);
}

template <typename P>
inline void
TracerParticleStepper<P>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("TracerParticleStepper::regrid(int, int, int)");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::regrid(int, int, int)" << endl;
  }

  // Define velocity field on the new mesh.
  m_amr->reallocate(m_velocity, m_phase, a_lmin);
  DataOps::setValue(m_velocity, 0.0);

  // Regrid tracer particles.
  m_solver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
}

template <typename P>
inline void
TracerParticleStepper<P>::postRegrid()
{
  CH_TIME("TracerParticleStepper::postRegrid()");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::postRegrid()" << endl;
  }

  // Update particle velocities.
  this->setVelocity();
  m_solver->interpolateVelocities();
}

template <typename P>
inline void
TracerParticleStepper<P>::setVelocity()
{
  CH_TIME("TracerParticleStepper::setVelocity()");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::setVelocity()" << endl;
  }

  std::function<RealVect(const RealVect a_position)> velFunc;

  switch (m_velocityField) {
  case VelocityField::Diagonal: {
    velFunc = [](const RealVect& a_position) -> RealVect {
      return RealVect::Unit;
    };

    break;
  }
  case VelocityField::Rotational: {
    velFunc = [](const RealVect pos) -> RealVect {
      const Real r     = pos.vectorLength();
      const Real theta = atan2(pos[1], pos[0]);

      return RealVect(D_DECL(-r * sin(theta), r * cos(theta), 0.));
    };

    break;
  }
  }

  DataOps::setValue(m_velocity,
                    velFunc,
                    m_amr->getProbLo(),
                    m_amr->getDx(),
                    m_amr->getMultiCutVofIterator(m_realm, m_phase));

  m_amr->conservativeAverage(m_velocity, m_realm, m_phase);
  m_amr->interpGhost(m_velocity, m_realm, m_phase);
}

template <typename P>
inline void
TracerParticleStepper<P>::initialParticles()
{
  CH_TIME("TracerParticleStepper::initialParticles()");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::initialParticles()" << endl;
  }

  // TLDR: This code draws distributed particles. Thanks to our nifty ParticleManagement and Random static classes this is just
  //       a matter of defining a distribution function. We define this as a lambda that draws particles uniformly distributed
  //       inside a box.
  //
  //       After creating the distribution we draw the particles (MPI rank distribution being taken care of under the hood) and
  //       set their weight to one. Finally, we add these particles to the ParticleContainer and remove the particles that lie inside
  //       the EB.

  // Create a random distribution which draw particles uniformly distributed in the domain.
  const RealVect probLo = m_amr->getProbLo();
  const RealVect probHi = m_amr->getProbHi();

  auto uniformDistribution = [probLo, probHi]() -> RealVect {
    RealVect ret = probLo;

    for (int dir = 0; dir < SpaceDim; dir++) {
      ret[dir] += (probHi - probLo)[dir] * Random::getUniformReal01();
    }

    return ret;
  };

  // Draw particles (unit weight) directly into the SoA container. drawRandomParticles partitions the
  // count across MPI ranks and remaps each particle to its owning patch/level/rank.
  ParticleContainerSoA<P>& solverParticles = m_solver->getParticles();

  ParticleManagement::drawRandomParticles(solverParticles, m_numInitialParticles, uniformDistribution);

  // Remove particles inside the EB.
  m_amr->removeCoveredParticlesIF(solverParticles, m_phase, 0.0);
}

template <typename P>
inline void
TracerParticleStepper<P>::advanceParticlesEuler(const Real a_dt)
{
  CH_TIME("TracerParticleStepper::advanceParticlesEuler()");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::advanceParticlesEuler()" << endl;
  }

  // TLDR: The new position is just x^(k+1) = x^k + dt*v^k

  ParticleContainerSoA<P>& amrParticles = m_solver->getParticles();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      ParticleSoA<P>&   leaf = amrParticles[lvl][din];
      const std::size_t n    = leaf.size();

      double* const pos[SpaceDim] = {D_DECL(leaf.positionColumn(0), leaf.positionColumn(1), leaf.positionColumn(2))};
      const ParticleReal* const vel[SpaceDim] = {
        D_DECL(leaf.template column<&P::vx>(), leaf.template column<&P::vy>(), leaf.template column<&P::vz>())};

      for (std::size_t i = 0; i < n; i++) {
        for (int dir = 0; dir < SpaceDim; dir++) {
          pos[dir][i] += vel[dir][i] * a_dt; // x^(k+1) = x^k + dt * v^k
        }
      }
    }
  }

  amrParticles.remap();
  m_amr->removeCoveredParticlesIF(amrParticles, m_phase, 0.0);

  m_solver->interpolateVelocities();
}

template <typename P>
inline void
TracerParticleStepper<P>::advanceParticlesRK2(const Real a_dt)
{
  CH_TIME("TracerParticleStepper::advanceParticlesRK2()");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::advanceParticlesRK2()" << endl;
  }

  // TLDR: The new positions are x^(k+1) = x^k + 0.5*dt*[ v(x^k) + v(x^*) ] where x^* = x^k + dt*v^k.

  ParticleContainerSoA<P>& amrParticles = m_solver->getParticles();

  // First step. Store old position (xk) and velocity (k1), then do the Euler advance.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      ParticleSoA<P>&   leaf = amrParticles[lvl][din];
      const std::size_t n    = leaf.size();

      double* const pos[SpaceDim] = {D_DECL(leaf.positionColumn(0), leaf.positionColumn(1), leaf.positionColumn(2))};
      const ParticleReal* const vel[SpaceDim] = {
        D_DECL(leaf.template column<&P::vx>(), leaf.template column<&P::vy>(), leaf.template column<&P::vz>())};
      ParticleReal* const xk[SpaceDim] = {
        D_DECL(leaf.template column<&P::xkx>(), leaf.template column<&P::xky>(), leaf.template column<&P::xkz>())};
      ParticleReal* const k1[SpaceDim] = {
        D_DECL(leaf.template column<&P::k1x>(), leaf.template column<&P::k1y>(), leaf.template column<&P::k1z>())};

      for (std::size_t i = 0; i < n; i++) {
        for (int dir = 0; dir < SpaceDim; dir++) {
          xk[dir][i] = pos[dir][i];          // store x^k
          k1[dir][i] = vel[dir][i];          // store v(x^k)
          pos[dir][i] += vel[dir][i] * a_dt; // x^* = x^k + dt * v(x^k)
        }
      }
    }
  }

  // Remap and interpolate the velocities again. This puts the velocity v = v(x^*) into the particles
  amrParticles.remap();
  m_amr->removeCoveredParticlesIF(amrParticles, m_phase, 0.0);
  m_solver->interpolateVelocities();

  // Do the second RK2 stage.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      ParticleSoA<P>&   leaf = amrParticles[lvl][din];
      const std::size_t n    = leaf.size();

      double* const pos[SpaceDim] = {D_DECL(leaf.positionColumn(0), leaf.positionColumn(1), leaf.positionColumn(2))};
      const ParticleReal* const vel[SpaceDim] = {
        D_DECL(leaf.template column<&P::vx>(), leaf.template column<&P::vy>(), leaf.template column<&P::vz>())};
      const ParticleReal* const xk[SpaceDim] = {
        D_DECL(leaf.template column<&P::xkx>(), leaf.template column<&P::xky>(), leaf.template column<&P::xkz>())};
      const ParticleReal* const k1[SpaceDim] = {
        D_DECL(leaf.template column<&P::k1x>(), leaf.template column<&P::k1y>(), leaf.template column<&P::k1z>())};

      for (std::size_t i = 0; i < n; i++) {
        for (int dir = 0; dir < SpaceDim; dir++) {
          // x^(k+1) = x^k + dt/2 * [ v(x^k) + v(x^*) ]
          pos[dir][i] = xk[dir][i] + 0.5 * a_dt * (k1[dir][i] + vel[dir][i]);
        }
      }
    }
  }

  // Remap and interpolate the velocities again
  amrParticles.remap();
  m_amr->removeCoveredParticlesIF(amrParticles, m_phase, 0.0);
  m_solver->interpolateVelocities();
}

template <typename P>
inline void
TracerParticleStepper<P>::advanceParticlesRK4(const Real a_dt)
{
  CH_TIME("TracerParticleStepper::advanceParticlesRK4()");
  if (m_verbosity > 5) {
    pout() << "TracerParticleStepper::advanceParticlesRK4()" << endl;
  }

  // TLDR: Just the standard RK4 method.

  ParticleContainerSoA<P>& amrParticles = m_solver->getParticles();

  const Real dtHalf  = a_dt / 2.0;
  const Real dtSixth = a_dt / 6.0;

  // k1 step.
  {
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        ParticleSoA<P>&   leaf = amrParticles[lvl][din];
        const std::size_t n    = leaf.size();

        double* const pos[SpaceDim] = {D_DECL(leaf.positionColumn(0), leaf.positionColumn(1), leaf.positionColumn(2))};
        const ParticleReal* const vel[SpaceDim] = {
          D_DECL(leaf.template column<&P::vx>(), leaf.template column<&P::vy>(), leaf.template column<&P::vz>())};
        ParticleReal* const xk[SpaceDim] = {
          D_DECL(leaf.template column<&P::xkx>(), leaf.template column<&P::xky>(), leaf.template column<&P::xkz>())};
        ParticleReal* const k1[SpaceDim] = {
          D_DECL(leaf.template column<&P::k1x>(), leaf.template column<&P::k1y>(), leaf.template column<&P::k1z>())};

        for (std::size_t i = 0; i < n; i++) {
          for (int dir = 0; dir < SpaceDim; dir++) {
            xk[dir][i]  = pos[dir][i];                       // store x^k
            k1[dir][i]  = vel[dir][i];                       // k1 = v(x^k)
            pos[dir][i] = xk[dir][i] + dtHalf * vel[dir][i]; // x = x^k + 0.5*dt*k1
          }
        }
      }
    }

    // Remap and compute v = v(x)
    amrParticles.remap();
    m_solver->interpolateVelocities();
  }

  // k2 step.
  {
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        ParticleSoA<P>&   leaf = amrParticles[lvl][din];
        const std::size_t n    = leaf.size();

        double* const pos[SpaceDim] = {D_DECL(leaf.positionColumn(0), leaf.positionColumn(1), leaf.positionColumn(2))};
        const ParticleReal* const vel[SpaceDim] = {
          D_DECL(leaf.template column<&P::vx>(), leaf.template column<&P::vy>(), leaf.template column<&P::vz>())};
        const ParticleReal* const xk[SpaceDim] = {
          D_DECL(leaf.template column<&P::xkx>(), leaf.template column<&P::xky>(), leaf.template column<&P::xkz>())};
        ParticleReal* const k2[SpaceDim] = {
          D_DECL(leaf.template column<&P::k2x>(), leaf.template column<&P::k2y>(), leaf.template column<&P::k2z>())};

        for (std::size_t i = 0; i < n; i++) {
          for (int dir = 0; dir < SpaceDim; dir++) {
            k2[dir][i]  = vel[dir][i];                       // k2 = v(x^k + 0.5*dt*k1)
            pos[dir][i] = xk[dir][i] + dtHalf * vel[dir][i]; // x = x^k + 0.5*dt*k2
          }
        }
      }
    }

    // Remap and compute v = v(x)
    amrParticles.remap();
    m_solver->interpolateVelocities();
  }

  // k3 step.
  {
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        ParticleSoA<P>&   leaf = amrParticles[lvl][din];
        const std::size_t n    = leaf.size();

        double* const pos[SpaceDim] = {D_DECL(leaf.positionColumn(0), leaf.positionColumn(1), leaf.positionColumn(2))};
        const ParticleReal* const vel[SpaceDim] = {
          D_DECL(leaf.template column<&P::vx>(), leaf.template column<&P::vy>(), leaf.template column<&P::vz>())};
        const ParticleReal* const xk[SpaceDim] = {
          D_DECL(leaf.template column<&P::xkx>(), leaf.template column<&P::xky>(), leaf.template column<&P::xkz>())};
        ParticleReal* const k3[SpaceDim] = {
          D_DECL(leaf.template column<&P::k3x>(), leaf.template column<&P::k3y>(), leaf.template column<&P::k3z>())};

        for (std::size_t i = 0; i < n; i++) {
          for (int dir = 0; dir < SpaceDim; dir++) {
            k3[dir][i]  = vel[dir][i];                     // k3 = v(x^k + 0.5*dt*k2)
            pos[dir][i] = xk[dir][i] + a_dt * vel[dir][i]; // x = x^k + dt*k3
          }
        }
      }
    }

    // Remap and compute v = v(x)
    amrParticles.remap();
    m_solver->interpolateVelocities();
  }

  // Final step.
  {
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        ParticleSoA<P>&   leaf = amrParticles[lvl][din];
        const std::size_t n    = leaf.size();

        double* const pos[SpaceDim] = {D_DECL(leaf.positionColumn(0), leaf.positionColumn(1), leaf.positionColumn(2))};
        const ParticleReal* const vel[SpaceDim] = {
          D_DECL(leaf.template column<&P::vx>(), leaf.template column<&P::vy>(), leaf.template column<&P::vz>())};
        const ParticleReal* const xk[SpaceDim] = {
          D_DECL(leaf.template column<&P::xkx>(), leaf.template column<&P::xky>(), leaf.template column<&P::xkz>())};
        const ParticleReal* const k1[SpaceDim] = {
          D_DECL(leaf.template column<&P::k1x>(), leaf.template column<&P::k1y>(), leaf.template column<&P::k1z>())};
        const ParticleReal* const k2[SpaceDim] = {
          D_DECL(leaf.template column<&P::k2x>(), leaf.template column<&P::k2y>(), leaf.template column<&P::k2z>())};
        const ParticleReal* const k3[SpaceDim] = {
          D_DECL(leaf.template column<&P::k3x>(), leaf.template column<&P::k3y>(), leaf.template column<&P::k3z>())};

        for (std::size_t i = 0; i < n; i++) {
          for (int dir = 0; dir < SpaceDim; dir++) {
            pos[dir][i] = xk[dir][i] + dtSixth * k1[dir][i] + dtHalf * k2[dir][i] + dtHalf * k3[dir][i] +
                          dtSixth * vel[dir][i];
          }
        }
      }
    }

    // Remap and compute v = v(x)
    amrParticles.remap();
    m_solver->interpolateVelocities();
  }

  m_amr->removeCoveredParticlesIF(amrParticles, m_phase, 0.0);
}

#include <CD_NamespaceFooter.H>

#endif
