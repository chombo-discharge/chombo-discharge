/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_BrownianWalkerStepper.cpp
 * @brief  Implementation of CD_BrownianWalkerStepper.H
 * @author Robert Marskar
 */

// Chombo includes
#include <ParmParse.H>
#include <CH_Timer.H>

// Our includes
#include <CD_BrownianWalkerStepper.H>
#include <CD_BrownianWalkerSpecies.H>
#include <CD_Random.H>
#include <CD_ParticleLoops.H>
#include <CD_ParallelOps.H>
#include <CD_EBCoarseToFineInterp.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::BrownianWalker;

BrownianWalkerStepper::BrownianWalkerStepper() : m_phase(phase::gas)
{
  CH_TIME("BrownianWalkerStepper::BrownianWalkerStepper");

  ParmParse pp("BrownianWalker");

  std::string str;
  pp.get("realm", m_realm);
  pp.get("diffco", m_diffCo);
  pp.get("mobility", m_mobility);
  pp.get("omega", m_omega);
  pp.get("verbosity", m_verbosity);
  pp.get("cfl", m_cfl);
  pp.get("load_balance", m_loadBalance);
  pp.get("which_balance", str);
  pp.get("verify_conservation", m_verifyConservation);

  if (str == "mesh") {
    m_whichLoadBalance = LoadBalancingMethod::Mesh;
  }
  else if (str == "particle") {
    m_whichLoadBalance = LoadBalancingMethod::Particle;
  }
  else {
    MayDay::Error(
      "BrownianWalkerStepper::BrownianWalkerStepper -- logic bust. Do not understand the load balancing argument 'which_balance'");
  }
}

BrownianWalkerStepper::BrownianWalkerStepper(RefCountedPtr<ItoSolver>& a_solver) : BrownianWalkerStepper()
{
  CH_TIME("BrownianWalkerStepper::BrownianWalkerStepper(full)");

  CH_assert(!a_solver.isNull());

  m_solver = a_solver;
}

BrownianWalkerStepper::~BrownianWalkerStepper()
{
  CH_TIME("BrownianWalkerStepper::~BrownianWalkerStepper");
}

void
BrownianWalkerStepper::parseRuntimeOptions()
{
  CH_TIME("BrownianWalkerStepper::parseRuntimeOptions");

  ParmParse pp("BrownianWalker");

  pp.get("verbosity", m_verbosity);
  pp.get("cfl", m_cfl);
  pp.get("load_balance", m_loadBalance);

  m_solver->parseRuntimeOptions();
}

void
BrownianWalkerStepper::initialData()
{
  CH_TIME("BrownianWalkerStepper::initialData");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::initialData" << endl;
  }

  // Fill initial particles and then make the desired number of superparticles.
  m_solver->initialData();
  this->makeSuperParticles();

  // Set advection and diffusion fields.
  this->setAdvectionDiffusion();

  // Set particle diffusion coefficient and mobility
  m_solver->setParticleDiffusion(m_diffCo);
  m_solver->setParticleMobility(m_mobility);

  m_solver->interpolateVelocities();
}

void
BrownianWalkerStepper::postInitialize()
{
  CH_TIME("BrownianWalkerStepper::postInitialize");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::postInitialize" << endl;
  }
}

void
BrownianWalkerStepper::setAdvectionDiffusion()
{
  CH_TIME("BrownianWalkerStepper::setAdvectionDiffusion");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::setAdvectionDiffusion" << endl;
  }

  if (m_solver->isDiffusive()) {
    this->setDiffusion();
  }

  if (m_solver->isMobile()) {
    this->setVelocity();
  }
}

void
BrownianWalkerStepper::setDiffusion()
{
  CH_TIME("BrownianWalkerStepper::setDiffusion");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::setDiffusion" << endl;
  }

  CH_assert(m_solver->isDiffusive());

  // The mesh diffusion field must hold the same coefficient the particles carry. The particle coefficients are
  // still assigned directly (setParticleDiffusion), so this field is not what moves the walkers -- but it is the
  // field that is plotted as 'dco', and it is the field a grad(D) drift correction would have to differentiate.
  // EBCellFAB::setVal covers ghost cells and the value is uniform, so no coarsening or ghost fill is needed here.
  m_solver->setDiffusionFunction(m_diffCo);
}

void
BrownianWalkerStepper::setVelocity()
{
  CH_TIME("BrownianWalkerStepper::setVelocity");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::setVelocity" << endl;
  }

  CH_assert(m_solver->isMobile());

  // TLDR: This just sets the velocity field everywhere.

  // Velocity field in solver
  EBAMRCellData& vel = m_solver->getVelocityFunction();

  // Nifty lambda describing the advective field
  auto veloFunc = [omega = this->m_omega](const RealVect& pos) -> RealVect {
    const Real r     = pos.vectorLength();
    const Real theta = atan2(pos[1], pos[0]);

    return RealVect(D_DECL(-r * omega * sin(theta), r * omega * cos(theta), 0.));
  };

  DataOps::setValue(vel, 0.0);
  DataOps::setValue(vel, veloFunc, m_amr->getProbLo(), m_amr->getDx(), m_amr->getMultiCutVofIterator(m_realm, m_phase));

  // Coarsen and update ghost cells.
  m_amr->conservativeAverage(vel, m_realm, m_phase);
  m_amr->interpGhost(vel, m_realm, m_phase);

  DataOps::setCoveredValue(vel, m_amr->getCoveredCells(m_realm, m_phase), 0, 0.0);
}

bool
BrownianWalkerStepper::loadBalanceThisRealm(const std::string& a_realm) const
{
  CH_TIME("BrownianWalkerStepper::loadBalanceThisRealm");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::loadBalanceThisRealm" << endl;
  }

  return m_loadBalance && (a_realm == m_realm);
}

void
BrownianWalkerStepper::loadBalanceBoxes(Vector<Vector<int>>&             a_procs,
                                        Vector<Vector<Box>>&             a_boxes,
                                        const std::string&               a_realm,
                                        const Vector<DisjointBoxLayout>& a_grids,
                                        const int                        a_lmin,
                                        const int                        a_finestLevel)
{
  CH_TIME("BrownianWalkerStepper::loadBalanceBoxes");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::loadBalanceBoxes" << endl;
  }

  CH_assert(m_loadBalance && a_realm == m_realm);

  switch (m_whichLoadBalance) {
  case LoadBalancingMethod::Mesh: {
    this->loadBalanceBoxesMesh(a_procs, a_boxes, a_realm, a_grids, a_lmin, a_finestLevel);

    break;
  }
  case LoadBalancingMethod::Particle: {
    this->loadBalanceBoxesParticles(a_procs, a_boxes, a_realm, a_grids, a_lmin, a_finestLevel);

    break;
  }
  default: {
    MayDay::Error("BrownianWalkerStepper::loadBalanceBoxes - logic bust");

    break;
  }
  }
}

#ifdef CH_USE_HDF5
void
BrownianWalkerStepper::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const
{
  CH_TIME("BrownianWalkerStepper::writeCheckpointData");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::writeCheckpointData" << endl;
  }

  m_solver->writeCheckpointLevel(a_handle, a_lvl);
}
#endif

#ifdef CH_USE_HDF5
void
BrownianWalkerStepper::readCheckpointData(HDF5Handle& a_handle, const int a_lvl)
{
  CH_TIME("BrownianWalkerStepper::readCheckpointData");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::readCheckpointData" << endl;
  }

  m_solver->readCheckpointLevel(a_handle, a_lvl);
}
#endif

void
BrownianWalkerStepper::postCheckpointSetup()
{
  CH_TIME("BrownianWalkerStepper::postCheckpointSetup");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::postCheckpointSetup" << endl;
  }

  // Update particle distribution
  m_solver->remap();
  this->makeSuperParticles();

  // Update advection-diffusion fields
  this->setAdvectionDiffusion();

  // Set particle diffusion coefficient and mobility
  m_solver->setParticleDiffusion(m_diffCo);
  m_solver->setParticleMobility(m_mobility);

  // Interpolate particle velocities.
  m_solver->interpolateVelocities();
}

int
BrownianWalkerStepper::getNumberOfPlotVariables() const
{
  CH_TIME("BrownianWalkerStepper::getNumberOfPlotVariables");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::getNumberOfPlotVariables" << endl;
  }

  const int numPlotVars = m_solver->getNumberOfPlotVariables();

  return numPlotVars;
}

Vector<std::string>
BrownianWalkerStepper::getPlotVariableNames() const
{
  CH_TIME("BrownianWalkerStepper::getPlotVariableNames");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::getPlotVariableNames" << endl;
  }

  return m_solver->getPlotVariableNames();
}

void
BrownianWalkerStepper::writePlotData(LevelData<EBCellFAB>& a_output,
                                     int&                  a_icomp,
                                     const std::string&    a_outputRealm,
                                     const int             a_level) const
{
  CH_TIME("BrownianWalkerStepper::writePlotData");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::writePlotData" << endl;
  }

  CH_assert(a_level >= 0);
  CH_assert(a_level <= m_amr->getFinestLevel());

  m_solver->writePlotData(a_output, a_icomp, a_outputRealm, a_level);
}

Real
BrownianWalkerStepper::computeDt()
{
  CH_TIME("BrownianWalkerStepper::computeDt");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::computeDt" << endl;
  }

  return m_cfl * m_solver->computeDt();
}

void
BrownianWalkerStepper::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt)
{
  CH_TIME("BrownianWalkerStepper::synchronizeSolverTimes");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::synchronizeSolverTimes" << endl;
  }

  // Solver needs to synchronize.
  m_solver->setTime(a_step, a_time, a_dt);

  // TimeStepper needs to synchronize.
  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;
}

void
BrownianWalkerStepper::printStepReport()
{
  CH_TIME("BrownianWalkerStepper::printStepReport");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::printStepReport" << endl;
  }

  // Do nothing
  const size_t localParticles  = m_solver->getNumParticles(ItoSolver::WhichContainer::Bulk, true);
  const size_t globalParticles = m_solver->getNumParticles(ItoSolver::WhichContainer::Bulk, false);

  pout() << "                                   #part = " << localParticles << " (" << globalParticles << ")" << endl;
}

bool
BrownianWalkerStepper::needToRegrid()
{
  CH_TIME("BrownianWalkerStepper::needToRegrid");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::needToRegrid" << endl;
  }

  return false;
}

void
BrownianWalkerStepper::preRegrid(const int a_lbase, const int a_oldFinestLevel)
{
  CH_TIME("BrownianWalkerStepper::preRegrid");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::preRegrid" << endl;
  }

  // clang-format off
  // TLDR: This does two things. The first is to deposit the number of particles per cell (ish) to the mesh. This can be used to load balance the application
  //       in the regrid step. The second this is that it puts all particle data holders in "regrid" mode.
  // clang-format on
  m_amr->allocate(m_regridPPC, m_realm, m_phase, 1);

  // ORDERING: this deposit must stay ahead of m_solver->preRegrid() below. That call moves the bulk
  // particles into the solver's reduced regrid holder and leaves the ItoParticle container empty, so
  // depositing afterwards would deposit nothing and loadBalanceBoxesMesh() would balance on zeros.
  // Deposit mass to scratch data holder. Then make sure the number of particles per cell
  m_solver->depositWeight(m_regridPPC,
                          m_solver->getParticles(ItoSolver::WhichContainer::Bulk),
                          DepositionType::NGP,
                          CoarseFineDeposition::Interp);

  m_solver->coarsenAndFillGhosts(m_regridPPC);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Real dx = m_amr->getDx()[lvl];
    const Real dV = std::pow(dx, SpaceDim);

    DataOps::scale(*m_regridPPC[lvl], dV);
  }

  // Solver knows what to do.
  m_solver->preRegrid(a_lbase, a_oldFinestLevel);
}

void
BrownianWalkerStepper::setupSolvers()
{
  CH_TIME("BrownianWalkerStepper::setupSolvers");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::setupSolvers" << endl;
  }

  m_species = RefCountedPtr<ItoSpecies>(new BrownianWalkerSpecies());

  m_solver->setVerbosity(m_verbosity);
  m_solver->parseOptions();
  m_solver->setAmr(m_amr);
  m_solver->setSpecies(m_species);
  m_solver->setPhase(m_phase);
  m_solver->setComputationalGeometry(m_computationalGeometry);
  m_solver->setRealm(m_realm);
}

void
BrownianWalkerStepper::registerRealms()
{
  CH_TIME("BrownianWalkerStepper::registerRealms");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::registerRealms" << endl;
  }

  m_amr->registerRealm(m_realm);
}

void
BrownianWalkerStepper::registerOperators()
{
  CH_TIME("BrownianWalkerStepper::registerOperators");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::registerOperators" << endl;
  }

  m_solver->registerOperators();
}

void
BrownianWalkerStepper::allocate()
{
  CH_TIME("BrownianWalkerStepper::allocate");

  // Allocate solver storage -- it knows what to do.
  m_solver->allocate();
}

Real
BrownianWalkerStepper::advance(const Real a_dt)
{
  CH_TIME("BrownianWalkerStepper::advance");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::advance" << endl;
  }

  CH_assert(m_solver->isMobile() || m_solver->isDiffusive());

  // TLDR: This function advances the particles using an Euler-Maruyama kernel. The steps are simply:
  //
  //          1. Compute Xnew = Xold + (V + grad(D))*dt + sqrt(2*D*dt)*N0 where N0 is a random number
  //          2. Remap the particles, assigning them to new grid boxes.
  //          3. Remove particles that struck the EB.
  //          4. Make super-particles.
  //          5. Update the particle velocities and diffusion coefficients.
  //          6. Deposit particles on mesh.
  //

  // The grad(D) drift correction, which makes the walkers transport v*n - D*grad(n) rather than the
  // v*n - grad(D*n) a plain Ito update gives. m_diffCo is uniform here, so the gradient is identically zero
  // and this changes nothing -- it is wired up so the code path is exercised and so that a spatially varying
  // diffusion coefficient would be handled correctly.
  const bool gradientDrift = m_solver->isDiffusive() && m_solver->isDiffusionGradientDrift();

  if (gradientDrift) {
    m_solver->computeDiffusionGradient();
  }

  // 1. Euler-Maruayma kernel on each patch.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    ParticleContainer<ItoParticle>& particles = m_solver->getParticles(ItoSolver::WhichContainer::Bulk);

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      // Particles that we iterate through.
      ParticleSoA<ItoParticle>& leaf = particles[lvl][din];

      const std::size_t n = leaf.size();

      // Euler step.
      if (m_solver->isMobile()) {
        double* oldPos[SpaceDim] = {D_DECL(leaf.column<&ItoParticle::old_x>(),
                                           leaf.column<&ItoParticle::old_y>(),
                                           leaf.column<&ItoParticle::old_z>())};

        const ParticleReal* v[SpaceDim] = {
          D_DECL(leaf.column<&ItoParticle::vx>(), leaf.column<&ItoParticle::vy>(), leaf.column<&ItoParticle::vz>())};

        double* const pos[SpaceDim] = {D_DECL(leaf.positionColumn(0), leaf.positionColumn(1), leaf.positionColumn(2))};

        ParticleLoops::loop(leaf, [&](const std::size_t i) {
          for (int dir = 0; dir < SpaceDim; dir++) {
            oldPos[dir][i] = pos[dir][i];
            pos[dir][i] += static_cast<Real>(v[dir][i]) * a_dt;
          }
        });
      }

      // grad(D) drift. This is a drift term like the one above, but it applies to any diffusive solver,
      // mobile or not. interpolateDiffusionGradient puts (grad D)(X_p) on the scratch vector columns.
      if (gradientDrift) {
        m_solver->interpolateDiffusionGradient(lvl, din);

        const ParticleReal* gradD[SpaceDim] = {D_DECL(leaf.column<&ItoParticle::scratch_x>(),
                                                      leaf.column<&ItoParticle::scratch_y>(),
                                                      leaf.column<&ItoParticle::scratch_z>())};

        double* const pos[SpaceDim] = {D_DECL(leaf.positionColumn(0), leaf.positionColumn(1), leaf.positionColumn(2))};

        ParticleLoops::loop(leaf, [&](const std::size_t i) {
          for (int dir = 0; dir < SpaceDim; dir++) {
            pos[dir][i] += static_cast<Real>(gradD[dir][i]) * a_dt;
          }
        });
      }

      // Diffusion hop
      if (m_solver->isDiffusive()) {
        const ParticleReal* D = leaf.column<&ItoParticle::diffusion>();

        for (std::size_t i = 0; i < n; i++) {
          const RealVect ran = Random::getNormal01() * Random::getDirection();
          const RealVect hop = ran * sqrt(2.0 * static_cast<Real>(D[i]) * a_dt);

          leaf.setPosition(i, leaf.position(i) + hop);
        }
      }
    }
  }

  // 2. Remap particles and assign them to correct patches. This discards particles outside the simulation domain.
  m_solver->remap();

  // 3. Particles that strike the EB are absorbed on it, and removed from the simulation.
  m_solver->removeCoveredParticles(EBRepresentation::ImplicitFunction, 0.0);

  // 4. Make new super-particles.
  this->makeSuperParticles();

  // 5. Update particle diffusion and velocities.
  m_solver->setParticleMobility(m_mobility);
  m_solver->setParticleDiffusion(m_diffCo);
  m_solver->interpolateVelocities();

  // Deposit onto mesh.
  m_solver->depositParticles();

  return a_dt;
}

void
BrownianWalkerStepper::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("BrownianWalkerStepper::regrid");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::regrid" << endl;
  }

  // Solver regrids. The super-particle merge happens inside that call now (see
  // ItoSolver.regrid_superparticles), on the reduced particles and before they are rebuilt as
  // ItoParticles -- merging again here would just re-merge an already-merged population.
  m_solver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
}

void
BrownianWalkerStepper::postRegrid()
{
  CH_TIME("BrownianWalkerStepper::postRegrid");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::postRegrid" << endl;
  }

  m_regridPPC.clear();

  // Update advection-diffusion fields
  this->setAdvectionDiffusion();

  // Set particle diffusion coefficient and mobility
  m_solver->setParticleDiffusion(m_diffCo);
  m_solver->setParticleMobility(m_mobility);

  // Interpolate particle velocities.
  m_solver->interpolateVelocities();
}

void
BrownianWalkerStepper::makeSuperParticles()
{
  CH_TIME("BrownianWalkerStepper::makeSuperParticles");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::makeSuperParticles" << endl;
  }

  if (m_solver->getParticlesPerCell()[0] > 0) {
    // A merge redistributes weight but must never create or destroy it. When requested, compute the
    // total particle weight on the container before and after the merge -- independent of the merge
    // scheme and of how the container is organized -- and abort if it drifts by more than round-off.
    // (makeSuperparticles() does any cell-sorting it needs internally and returns the container
    // patch-organized.)
    const Real weightBefore = m_verifyConservation ? this->computeTotalWeight() : 0.0;

    m_solver->makeSuperparticles(ItoSolver::WhichContainer::Bulk);

    if (m_verifyConservation) {
      const Real weightAfter = this->computeTotalWeight();
      const Real absDrift    = std::abs(weightAfter - weightBefore);

      if (absDrift > 1.E-9 * std::max(1.0, std::abs(weightBefore))) {
        MayDay::Abort("BrownianWalkerStepper::makeSuperParticles -- superparticle merge did not conserve "
                      "total particle weight");
      }
    }
  }
}

Real
BrownianWalkerStepper::computeTotalWeight()
{
  CH_TIME("BrownianWalkerStepper::computeTotalWeight");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::computeTotalWeight" << endl;
  }

  const ParticleContainer<ItoParticle>& particles = m_solver->getParticles(ItoSolver::WhichContainer::Bulk);

  Real weight = 0.0;

  for (int lvl = 0; lvl <= particles.getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = particles.getGrids()[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const auto&         leaf = particles[lvl][dit()];
      const double* const w    = leaf.weightColumn();

      weight = ParticleLoops::reduce(leaf, weight, [&](const Real a_acc, const std::size_t a_i) {
        return a_acc + w[a_i];
      });
    }
  }

  return ParallelOps::sum(weight);
}

void
BrownianWalkerStepper::loadBalanceBoxesMesh(Vector<Vector<int>>&             a_procs,
                                            Vector<Vector<Box>>&             a_boxes,
                                            const std::string&               a_realm,
                                            const Vector<DisjointBoxLayout>& a_grids,
                                            const int                        a_lmin,
                                            const int                        a_finestLevel)
{
  CH_TIME("BrownianWalkerStepper::loadBalanceBoxesMesh");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::loadBalanceBoxesMesh" << endl;
  }

  CH_assert(m_loadBalance && a_realm == m_realm);

  // clang-format off
  // TLDR: This routine is called AFTER AmrMesh::regridAMR which means that we have all EB-related information we need for building operators. We happen to
  //       know that ItoSolver computed the number of particles per cell in the preRegrid method and that these values are returned by a call to
  //       EBAMRCellData& ItoSolver::getScratch(). We take that data and regrid it onto the new grids. This requires us to manually build an operator which
  //       can do that interpolation.
  //
  //       Once we've put that data on the new mesh, we can simply compute the sum of all mesh data in each grid patch. That sum is equal to the number of particles
  //       in the patch, which we can use for load balancing.
  // clang-format on

  constexpr int comp = 0;

  // Make some storage for the number of particles per cell on the new grids.
  EBAMRCellData newParticlesPerCell;
  m_amr->allocate(newParticlesPerCell, a_realm, m_phase, 1);

  // Grid information.
  const Vector<RefCountedPtr<EBLevelGrid>>& eblg     = m_amr->getEBLevelGrid(a_realm, m_phase);
  const Vector<RefCountedPtr<EBLevelGrid>>& eblgCoFi = m_amr->getEBLevelGridCoFi(a_realm, m_phase);
  const Vector<int>&                        refRat   = m_amr->getRefinementRatios();

  // Copy old mesh data to new mesh data.
  for (int lvl = 0; lvl <= std::max(0, a_lmin - 1); lvl++) {
    m_regridPPC[lvl]->copyTo(*newParticlesPerCell[lvl]);
  }

  // Now regrid where we got new grids.
  for (int lvl = a_lmin; lvl <= a_finestLevel; lvl++) {
    if (lvl > 0) {
      EBCoarseToFineInterp fineInterp(*eblg[lvl], *eblgCoFi[lvl - 1], *eblg[lvl - 1], refRat[lvl - 1]);

      fineInterp.interpolate(*newParticlesPerCell[lvl],
                             *newParticlesPerCell[lvl - 1],
                             Interval(0, 0),
                             EBCoarseToFineInterp::Type::PWC);

      // Replace data where old region overlapped new region.
      if (lvl < std::min(newParticlesPerCell.size(), m_regridPPC.size())) {
        m_regridPPC[lvl]->copyTo(*newParticlesPerCell[lvl]);
      }
    }
  }

  // At this point we need to replace the data UNDERNEATH the fine grids and in the covered cells -- we don't want to
  // count data in those regions are valid data (because it does not represent particles). After that's done,
  // newParticlesPerCel should be a realistic representation of the number of particles per cell.
  constexpr Real zero = 0.0;

  DataOps::setInvalidValue(newParticlesPerCell, refRat, zero);
  DataOps::setCoveredValue(newParticlesPerCell, m_amr->getCoveredCells(a_realm, m_phase), zero);

  // We now have everything we need to start load balancing.
  a_procs.resize(1 + a_finestLevel);
  a_boxes.resize(1 + a_finestLevel);

  // These levels should not change.
  for (int lvl = 0; lvl < a_lmin; lvl++) {
    a_procs[lvl] = a_grids[lvl].procIDs();
    a_boxes[lvl] = a_grids[lvl].boxArray();
  }

  Loads rankLoads;
  rankLoads.resetLoads();

  // Load balance these levels.
  for (int lvl = a_lmin; lvl <= a_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = a_grids[lvl];
    const DataIterator&      dit = dbl.dataIterator();
    const Real               dx  = m_amr->getDx()[lvl];
    const Real               dV  = std::pow(dx, SpaceDim);

    // Boxes, loads, and ranks. This is stuff to be assigned. Note that the input boxes are lexicographically sorted
    // (this is what DisjointBoxLayout does).
    Vector<Box>      boxes = dbl.boxArray();
    Vector<long int> loads;
    Vector<int>      ranks;

    loads.resize(boxes.size(), 0L);
    ranks.resize(boxes.size());

    // Compute loads in each grid patch.
    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const Box            cellBox          = dbl[din];
      const EBCellFAB&     particlesPerCell = (*newParticlesPerCell[lvl])[din];
      const EBISBox&       ebisBox          = particlesPerCell.getEBISBox();
      const BaseFab<Real>& ppcFAB           = particlesPerCell.getSingleValuedFAB();

      // Compute the number of particles in the patch. Note that we try to estimate the number of computational
      // in each cell. This might not be the same as the contents in newParticlesPerCell*dV because the stepper will
      // later run with superparticle merging/splitting.
      Real sum = 0.0;

      // We will have at most this many particles per grid cell after the merge. This kernel does that.
      // The target is the solver's (ItoSolver.particles_per_cell), read with the same per-level clamp
      // the merge itself uses -- this loop is already per level, so a per-level target is honoured.
      const Vector<int>& ppcPerLevel = m_solver->getParticlesPerCell();
      const int          ppc         = (lvl < ppcPerLevel.size()) ? ppcPerLevel[lvl] : ppcPerLevel.back();

      auto kernel = [&](const IntVect& iv) -> void {
        if (!ebisBox.isCovered(iv)) {
          const Real numPhysParticles = std::abs(ppcFAB(iv, comp) * dV);
          const Real numCompParticles = std::min(numPhysParticles, Real(ppc));

          sum += numCompParticles;
        }
      };

      // Not vectorizable: local FP sum reduction over the box + out-of-line isCovered guard. One-time
      // load balancing (per regrid). sum is box-local and loads[din.intCode()] is per-box -> no race.
      BoxLoops::loop<D_DECL(1, 1, 1)>(cellBox, kernel);

      loads[din.intCode()] = lround(sum);
    }

    // If running with MPI, loads must be gathered on all ranks.
    ParallelOps::sum(loads);

    // Sort the boxes and loads using a Morton code. Then load balance the application.
    LoadBalancing::sort(boxes, loads, BoxSorting::Morton);
    LoadBalancing::makeBalance(ranks, rankLoads, loads, boxes);

    // Assign ranks to boxes
    a_boxes[lvl] = boxes;
    a_procs[lvl] = ranks;
  }
}

void
BrownianWalkerStepper::loadBalanceBoxesParticles(Vector<Vector<int>>&             a_procs,
                                                 Vector<Vector<Box>>&             a_boxes,
                                                 const std::string&               a_realm,
                                                 const Vector<DisjointBoxLayout>& a_grids,
                                                 const int                        a_lmin,
                                                 const int                        a_finestLevel)
{
  CH_TIME("BrownianWalkerStepper::loadBalanceBoxesParticles");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::loadBalanceBoxesParticles" << endl;
  }

  CH_assert(m_loadBalance && a_realm == m_realm);

  // clang-format off
  // TLDR: This load balancing method computes the number of particles in the new grids directly. It does so by remapping the particles
  //       to the new grids and then simply counting them. This is then used for load balancing. The downside of this method is that the
  //       particles needs to be remapped twice (once here, and once in BrownianWalkerStepper::regrid). So, this method is usually slower
  //       than the other one when the number of particles is large.
  // clang-format on

  // The bulk ItoParticle container is EMPTY at this point -- Driver::regrid() calls loadBalanceBoxes()
  // between the solver's preRegrid() and regrid(), and preRegrid() moved the particles into the
  // solver's reduced regrid holder. Count that instead. Reading getParticles(Bulk) here would produce
  // an all-zero load vector, which LoadBalancing happily turns into a valid, completely unbalanced
  // layout with nothing to indicate that anything went wrong.
  ParticleContainer<ItoMergeParticle>& particles = m_solver->getRegridParticles();

  // Regrid the particles onto the new mesh (SoA regrid rebuilds over the new layout from the preRegrid cache).
  // a_grids are the proxy grids (m_amr->getProxyGrids()) on which the Realm built its tile->box maps
  // during regridAmr, so we alias those rather than building our own.
  particles.regrid(a_grids,
                   m_amr->getDomains(),
                   m_amr->getDx(),
                   m_amr->getRefinementRatios(),
                   m_amr->getMinBlockSize(),
                   m_amr->getLevelTiles(a_realm),
                   a_finestLevel);

  a_procs.resize(1 + a_finestLevel);
  a_boxes.resize(1 + a_finestLevel);

  // These levels are not to be load balanced.
  for (int lvl = 0; lvl < a_lmin; lvl++) {
    a_procs[lvl] = a_grids[lvl].procIDs();
    a_boxes[lvl] = a_grids[lvl].boxArray();
  }

  Loads rankLoads;
  rankLoads.resetLoads();

  // Load balance these levels.
  for (int lvl = a_lmin; lvl <= a_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = a_grids[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    Vector<Box>      boxes;
    Vector<long int> loads;
    Vector<int>      ranks;

    // Boxes -- note that these are currently lexicographically ordered.
    boxes = dbl.boxArray();

    loads.resize(boxes.size(), 0L);
    ranks.resize(boxes.size());

    // Count the number of particles in each grid patch.
    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      loads[din.intCode()] = long(particles[lvl][din].size());
    }

    // If running with MPI, loads must be gathered on all ranks.
    ParallelOps::sum(loads);

    // Sort the boxes and loads using a Morton code. Then load balance the application.
    LoadBalancing::sort(boxes, loads, BoxSorting::Morton);
    LoadBalancing::makeBalance(ranks, rankLoads, loads, boxes);

    // Assign ranks to boxes
    a_boxes[lvl] = boxes;
    a_procs[lvl] = ranks;
  }

  // Put particles back (re-cache for the subsequent regrid).
  particles.preRegrid();
}

Vector<long int>
BrownianWalkerStepper::getCheckpointLoads(const std::string& a_realm, const int a_level) const
{
  CH_TIME("BrownianWalkerStepper::getCheckpointLoads(string, int)");
  if (m_verbosity > 5) {
    pout() << "BrownianWalkerStepper::getCheckpointLoads(string, int)" << endl;
  }

  Vector<long int> loads;

  const DisjointBoxLayout& dbl      = m_amr->getGrids(a_realm)[a_level];
  const DataIterator&      dit      = dbl.dataIterator();
  const Vector<Box>        boxArray = dbl.boxArray();

  loads.resize(boxArray.size(), 0L);

  // ORDERING: safe only because every caller runs outside the preRegrid()->regrid() window, where the
  // bulk container is empty. Driver reaches this through writeLoads() when writing a checkpoint, a plot
  // file, or the pre-regrid file -- and the pre-regrid file is written before Driver::regrid() is
  // entered at all, hence before preRegrid(). Anything that moves it inside that window reads zeros.
  const ParticleContainer<ItoParticle>& particles = m_solver->getParticles(ItoSolver::WhichContainer::Bulk);

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    loads[din.intCode()] = long(particles[a_level][din].size());
  }

  // If running with MPI, loads must be gathered on all ranks.
  ParallelOps::sum(loads);

  return loads;
}

#include <CD_NamespaceFooter.H>
