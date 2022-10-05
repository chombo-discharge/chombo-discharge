/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaStepper.cpp
  @brief  Implementation of CD_ItoPlasmaStepper.H
  @author Robert Marskar
*/

// Chombo includes
#include <EBArith.H>
#include <PolyGeom.H>
#include <ParmParse.H>

// Our includes
#include <CD_ItoPlasmaStepper.H>
#include <CD_DataOps.H>
#include <CD_Units.H>
#include <CD_FieldSolverMultigrid.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoPlasma;

ItoPlasmaStepper::ItoPlasmaStepper()
{
  CH_TIME("ItoPlasmaStepper::ItoPlasmaStepper");

  // Default settings.
  m_verbosity               = -1;
  m_name                    = "ItoPlasmaStepper";
  m_phase                   = phase::gas;
  m_dt                      = 0.0;
  m_time                    = 0.0;
  m_timeStep                = 0;
  m_loadPerCell             = 1.0;
  m_useNewReactionAlgorithm = true;
  m_regridSuperparticles    = true;
  m_fluidRealm              = Realm::Primal;
  m_particleRealm           = Realm::Primal;

  this->parseOptions();
}

ItoPlasmaStepper::ItoPlasmaStepper(RefCountedPtr<ItoPlasmaPhysics>& a_physics) : ItoPlasmaStepper()
{
  CH_TIME("ItoPlasmaStepper::ItoPlasmaStepper(RefCountrPtr<ItoPlasmaPhysics>)");

  m_physics = a_physics;
}

ItoPlasmaStepper::~ItoPlasmaStepper() {
  CH_TIME("ItoPlasmaStepper::~ItoPlasmaStepper");
}

void
ItoPlasmaStepper::parsOptions() {
  CH_TIME("ItoPlasmaStepper::parseOptions");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::parseOptions" << endl;
  }
}

void
ItoPlasmaStepper::parseRuntimeOptions() {
  CH_TIME("ItoPlasmaStepper::parseRuntimeOptions");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::parseRuntimeOptions" << endl;
  }
}

void
ItoPlasmaStepper::setVerbosity(const int a_verbosity)
{
  CH_TIME("ItoPlasmaStepper::setVerbosity");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::setVerbosity" << endl;
  }
  m_verbosity = a_verbosity;
}

void
ItoPlasmaStepper::setupSolvers()
{
  CH_TIME("ItoPlasmaStepper::setup_solver");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::setupSolvers" << endl;
  }

  // Parse class options
  this->parseOptions();

  // Set up solvers
  this->setupIto();
  this->setupPoisson();
  this->setupRadiativeTransfer();
  this->setupSigma();

  // Allocate internal stuff
  this->allocateInternals();

  // Set the particle buffers for the Ito solver
  this->setParticleBuffers();
}

void
ItoPlasmaStepper::setupIto()
{
  CH_TIME("ItoPlasmaStepper::setupIto");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::setupIto" << endl;
  }

  m_ito->setVerbosity(m_verbosity);
  m_ito->parseOptions();
  m_ito->setAmr(m_amr);
  m_ito->setPhase(m_phase);
  m_ito->setComputationalGeometry(m_computationalGeometry);
  m_ito->setRealm(m_particleRealm);
}

void
ItoPlasmaStepper::setupPoisson()
{
  CH_TIME("ItoPlasmaStepper::setupPoisson");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::setupPoisson" << endl;
  }

  m_fieldSolver->setVerbosity(m_verbosity);
  m_fieldSolver->parseOptions();
  m_fieldSolver->setAmr(m_amr);
  m_fieldSolver->setComputationalGeometry(m_computationalGeometry);
  m_fieldSolver->setVoltage(m_potential);
  m_fieldSolver->setRealm(m_fluidRealm);
}

void
ItoPlasmaStepper::setupRadiativeTransfer()
{
  CH_TIME("ItoPlasmaStepper::setupRadiativeTransfer");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::setupRadiativeTransfer" << endl;
  }

  m_rte->setVerbosity(m_verbosity);
  m_rte->parseOptions();
  m_rte->setPhase(m_phase);
  m_rte->setAmr(m_amr);
  m_rte->setComputationalGeometry(m_computationalGeometry);
  m_rte->setRealm(m_particleRealm);
  m_rte->sanityCheck();
}

void
ItoPlasmaStepper::setupSigma()
{
  CH_TIME("ItoPlasmaStepper::setupSigma");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::setupSigma" << endl;
  }

  m_sigma = RefCountedPtr<SurfaceODESolver<1>>(new SurfaceODESolver<1>(m_amr));
  m_sigma->setVerbosity(m_verbosity);
  m_sigma->setRealm(m_fluidRealm);
  m_sigma->setPhase(m_phase);
  m_sigma->setName("Surface charge");
}

void
ItoPlasmaStepper::setParticleBuffers()
{
  CH_TIME("ItoPlasmaStepper::setParticleBuffers");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::setParticleBuffers" << endl;
  }
}

void
ItoPlasmaStepper::allocate()
{
  CH_TIME("ItoPlasmaStepper::allocate");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::allocate" << endl;
  }

  m_ito->allocateInternals();
  m_rte->allocateInternals();
  m_fieldSolver->allocateInternals();
  m_sigma->allocate();
}

void
ItoPlasmaStepper::postInitialize()
{}

void
ItoPlasmaStepper::initialData()
{
  CH_TIME("ItoPlasmaStepper::initialData");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::initialData" << endl;
  }

  m_fieldSolver->setPermittivities(); // Set permittivities for Poisson operator
  m_ito->initialData();               // This deposits, of course.
  m_rte->initialData();
  this->initialSigma();

  m_ito->sortParticlesByCell(ItoSolver::WhichContainer::Bulk);
  m_ito->makeSuperparticles(ItoSolver::WhichContainer::Bulk, m_particlesPerCell);
  m_ito->sortParticlesByPatch(ItoSolver::WhichContainer::Bulk);

  // Solve Poisson equation and compute the E-field
  this->solvePoisson();

  // Fill solvers with velocities and diffusion coefficients
  this->computeItoVelocities();
  this->computeItoDiffusion();
}

void
ItoPlasmaStepper::initialSigma()
{
  CH_TIME("ItoPlasmaStepper::initialSigma");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::initialSigma" << endl;
  }

  const RealVect origin       = m_amr->getProbLo();
  const int      finest_level = m_amr->getFinestLevel();

  EBAMRIVData& sigma = m_sigma->getPhi();

  for (int lvl = 0; lvl <= finest_level; lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_fluidRealm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_fluidRealm, phase::gas)[lvl];
    const Real               dx    = m_amr->getDx()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      BaseIVFAB<Real>& state = (*sigma[lvl])[dit()];

      const EBISBox&    ebisbox = ebisl[dit()];
      const IntVectSet& ivs     = state.getIVS();
      const EBGraph&    ebgraph = state.getEBGraph();

      for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit) {
        const VolIndex& vof = vofit();
        const RealVect  pos = origin + vof.gridIndex() * dx + 0.5 * ebisbox.bndryCentroid(vof) * dx;

        for (int comp = 0; comp < state.nComp(); comp++) {
          state(vof, comp) = m_physics->initialSigma(m_time, pos);
        }
      }
    }
  }

  m_amr->conservativeAverage(sigma, m_fluidRealm, phase::gas);
  m_sigma->resetElectrodes(sigma, 0.0);
}

void
ItoPlasmaStepper::postCheckpointSetup()
{
  CH_TIME("ItoPlasmaStepper::postCheckpointSetup");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::postCheckpointSetup" << endl;
  }

  //this->solvePoisson();
  this->allocateInternals();

  m_ito->remap();

  this->postCheckpointPoisson();

  this->computeItoVelocities();
  this->computeItoDiffusion();
}

void
ItoPlasmaStepper::postCheckpointPoisson()
{
  CH_TIME("ItoPlasmaStepper::postCheckpointPoisson");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::postCheckpointPoisson" << endl;
  }

  // Do some post checkpointing stuff.
  m_fieldSolver->postCheckpoint();

  // Do ghost cells and then compute E
  MFAMRCellData& state = m_fieldSolver->getPotential();

  m_amr->conservativeAverage(state, m_fluidRealm);
  m_amr->interpGhost(state, m_fluidRealm);

  m_fieldSolver->computeElectricField(); // Solver checkpoints the potential. Now compute the field.

  // Interpolate the fields to centroids
  EBAMRCellData E;
  m_amr->allocatePointer(E);
  m_amr->alias(E, m_phase, m_fieldSolver->getElectricField());

  // Fluid Realm
  m_fluid_E.copy(E);
  m_amr->conservativeAverage(m_fluid_E, m_fluidRealm, m_phase);
  m_amr->interpGhostPwl(m_fluid_E, m_fluidRealm, m_phase);
  m_amr->interpToCentroids(m_fluid_E, m_fluidRealm, m_phase);

  // Particle Realm
  m_particle_E.copy(E);
  m_amr->conservativeAverage(m_particle_E, m_particleRealm, m_phase);
  m_amr->interpGhostPwl(m_particle_E, m_particleRealm, m_phase);
  m_amr->interpToCentroids(m_particle_E, m_particleRealm, m_phase);

  // Compute maximum E
  // const Real Emax = this->computeMaxElectricField(m_phase);
  // std::cout << Emax << std::endl;
}

#ifdef CH_USE_HDF5
void
ItoPlasmaStepper::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const
{
  CH_TIME("ItoPlasmaStepper::writeCheckpointData");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::writeCheckpointData" << endl;
  }

  for (ItoIterator<ItoSolver> solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    const RefCountedPtr<ItoSolver>& solver = solver_it();
    solver->writeCheckpointLevel(a_handle, a_lvl);
  }

  for (RtIterator<McPhoto> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it) {
    const RefCountedPtr<McPhoto>& solver = solver_it();
    solver->writeCheckpointLevel(a_handle, a_lvl);
  }

  m_fieldSolver->writeCheckpointLevel(a_handle, a_lvl);
  m_sigma->writeCheckpointLevel(a_handle, a_lvl);
}
#endif

#ifdef CH_USE_HDF5
void
ItoPlasmaStepper::readCheckpointData(HDF5Handle& a_handle, const int a_lvl)
{
  CH_TIME("ItoPlasmaStepper::readCheckpointData");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::readCheckpointData" << endl;
  }

  for (ItoIterator<ItoSolver> solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>& solver = solver_it();
    solver->readCheckpointLevel(a_handle, a_lvl);
  }

  for (RtIterator<McPhoto> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<McPhoto>& solver = solver_it();
    solver->readCheckpointLevel(a_handle, a_lvl);
  }

  m_fieldSolver->readCheckpointLevel(a_handle, a_lvl);
  m_sigma->readCheckpointLevel(a_handle, a_lvl);
}
#endif

void
ItoPlasmaStepper::writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const
{
  CH_TIME("ItoPlasmaStepper::writePlotData");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::writePlotData" << endl;
  }

  // Poisson solver copies over its output data
  a_plotVariableNames.append(m_fieldSolver->getPlotVariableNames());
  m_fieldSolver->writePlotData(a_output, a_icomp);

  // Surface charge solver writes
  a_plotVariableNames.append(m_sigma->getPlotVariableNames());
  m_sigma->writePlotData(a_output, a_icomp);

  // Ito solvers copy their output data
  for (ItoIterator<ItoSolver> solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>& solver = solver_it();
    a_plotVariableNames.append(solver->getPlotVariableNames());
    solver->writePlotData(a_output, a_icomp);
  }

  // RTE solvers copy their output data
  for (RtIterator<McPhoto> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<McPhoto>& solver = solver_it();
    a_plotVariableNames.append(solver->getPlotVariableNames());
    solver->writePlotData(a_output, a_icomp);
  }

  // Write the current to the output
  this->writeJ(a_output, a_icomp);
  a_plotVariableNames.push_back("x-J");
  a_plotVariableNames.push_back("y-J");
  if (SpaceDim == 3) {
    a_plotVariableNames.push_back("z-J");
  }

  // Write the number of particles per patch
  this->writeNumParticlesPerPatch(a_output, a_icomp);
  a_plotVariableNames.push_back("particles_per_patch");
}

void
ItoPlasmaStepper::writeJ(EBAMRCellData& a_output, int& a_icomp) const
{
  CH_TIME("ItoPlasmaStepper::writeJ");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::writeJ" << endl;
  }

  const Interval src_interv(0, SpaceDim - 1);
  const Interval dst_interv(a_icomp, a_icomp + SpaceDim - 1);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    if (m_J.getRealm() == a_output.getRealm()) {
      m_J[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
    }
    else {
      m_J[lvl]->copyTo(src_interv, *a_output[lvl], dst_interv);
    }
  }
  a_icomp += SpaceDim;
}

void
ItoPlasmaStepper::writeNumParticlesPerPatch(EBAMRCellData& a_output, int& a_icomp) const
{
  CH_TIME("ItoPlasmaStepper::writeNumParticlesPerPatch");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::writeNumParticlesPerPatch" << endl;
  }

  const Interval src_interv(0, 0);
  const Interval dst_interv(a_icomp, a_icomp);

  DataOps::setValue(m_particle_scratch1, 0.0);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    const ParticleContainer<ItoParticle>& particles = solver_it()->getParticles(ItoSolver::WhichContainer::Bulk);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[lvl];
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
        (*m_particle_scratch1[lvl])[dit()] += particles[lvl][dit].numItems();
      }
    }
  }

  // Copy to output holder
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    if (m_particle_scratch1.getRealm() == a_output.getRealm()) {
      m_particle_scratch1[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
    }
    else {
      m_particle_scratch1[lvl]->copyTo(src_interv, *a_output[lvl], dst_interv);
    }
  }

  a_icomp += 1;
}

void
ItoPlasmaStepper::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt)
{
  CH_TIME("ItoPlasmaStepper::synchronizeSolverTimes");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::synchronizeSolverTimes" << endl;
  }

  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;

  m_ito->setTime(a_step, a_time, a_dt);
  m_fieldSolver->setTime(a_step, a_time, a_dt);
  m_rte->setTime(a_step, a_time, a_dt);
  m_sigma->setTime(a_step, a_time, a_dt);
}

void
ItoPlasmaStepper::printStepReport()
{
  CH_TIME("ItoPlasmaStepper::printStepReport");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::printStepReport" << endl;
  }

  const Real Emax = this->computeMaxElectricField(m_phase);

  const size_t l_particles = m_ito->getNumParticles(ItoSolver::WhichContainer::Bulk, true);
  const size_t g_particles = m_ito->getNumParticles(ItoSolver::WhichContainer::Bulk, false);

  const size_t l_eb_particles = m_ito->getNumParticles(ItoSolver::WhichContainer::EB, true);
  const size_t g_eb_particles = m_ito->getNumParticles(ItoSolver::WhichContainer::EB, false);

  const size_t l_domain_particles = m_ito->getNumParticles(ItoSolver::WhichContainer::Domain, true);
  const size_t g_domain_particles = m_ito->getNumParticles(ItoSolver::WhichContainer::Domain, false);

  const size_t l_source_particles = m_ito->getNumParticles(ItoSolver::WhichContainer::Source, true);
  const size_t g_source_particles = m_ito->getNumParticles(ItoSolver::WhichContainer::Source, false);

  Real avg;
  Real sigma;

  int minRank;
  int maxRank;

  size_t minParticles;
  size_t maxParticles;

  // Compute some particle statistics
  this->getParticleStatistics(avg, sigma, minParticles, maxParticles, minRank, maxRank);

  // How was the time step restricted
  std::string str;
  switch (m_timeCode) {
  case TimeCode::Physics:
    str = "dt restricted by 'physics'";
    break;
  case TimeCode::Advection:
    str = "dt restricted by 'advection'";
    break;
  case TimeCode::RelaxationTime:
    str = "dt restricted by 'relaxation time'";
    break;
  case TimeCode::Hardcap:
    str = "dt restricted by 'hardcap'";
    break;
  default:
    str = "dt restricted by 'unspecified'";
    break;
  }
  pout() << "                                   " + str << endl;
  pout() << "                                   Emax      = " << Emax << endl
         << "                                   #part     = " << l_particles << " (" << g_particles << ")" << endl
         << "                                   #eb part  = " << l_eb_particles << " (" << g_eb_particles << ")" << endl
         << "                                   #dom part = " << l_domain_particles << " (" << g_domain_particles << ")"
         << endl
         << "                                   #src part = " << l_source_particles << " (" << g_source_particles << ")"
         << endl
         << "                                   #part min = " << minParticles << " (on rank = " << minRank << ")"
         << endl
         << "                                   #part max = " << maxParticles << " (on rank = " << maxRank << ")"
         << endl
         << "                                   #part avg = " << avg << endl
         << "                                   #part dev = " << sigma << " (" << 100. * sigma / avg << "%)" << endl;
}

void
ItoPlasmaStepper::getParticleStatistics(Real&   a_avg,
                                        Real&   a_sigma,
                                        size_t& a_minPart,
                                        size_t& a_maxPart,
                                        int&    a_minRank,
                                        int&    a_maxRank)
{
  CH_TIME("ItoPlasmaStepper::getParticleStatistics");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::getParticleStatistics" << endl;
  }

  const int srcProc = 0;
  const int nProc   = numProc();

  const size_t numLocal = m_ito->getNumParticles(ItoSolver::WhichContainer::Bulk, true);

  // Gather on source proc
  Vector<size_t> allCounts(nProc);
  gather(allCounts, numLocal, srcProc);

  // Compute and broadcast the average and the standard deviation
  if (procID() == srcProc) {

    a_avg = 0.0;
    for (int i = 0; i < nProc; i++) {
      a_avg += 1.0 * allCounts[i];
    }
    a_avg *= 1. / nProc;

    a_sigma = 0.0;
    for (int i = 0; i < nProc; i++) {
      a_sigma += std::pow(1.0 * allCounts[i] - a_avg, 2);
    }
    a_sigma = sqrt(a_sigma / nProc);
  }

  broadcast(a_avg, srcProc);
  broadcast(a_sigma, srcProc);

  // Get the minimum/maximum number of particles
  a_minRank = srcProc;
  a_maxRank = srcProc;

  a_minPart = std::numeric_limits<size_t>::max();
  a_maxPart = 0;

  if (procID() == srcProc) {
    for (int i = 0; i < nProc; i++) {
      if (allCounts[i] < a_minPart) {
        a_minPart = allCounts[i];
        a_minRank = i;
      }

      if (allCounts[i] > a_maxPart) {
        a_maxPart = allCounts[i];
        a_maxRank = i;
      }
    }
  }

  broadcast(a_minRank, srcProc);
  broadcast(a_maxRank, srcProc);

  broadcast(a_minPart, srcProc);
  broadcast(a_maxPart, srcProc);
}

void
ItoPlasmaStepper::printTimerDiagnostics(Real& a_timer, const std::string a_prefix)
{
  CH_TIME("ItoPlasmaStepper::printTimerDiagnostics");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::printTimerDiagnostics" << endl;
  }

  const int srcProc = 0;
  const int nProc   = numProc();

  const size_t numLocal = m_ito->getNumParticles(ItoSolver::WhichContainer::Bulk, true);

  // Gather all timers on source proc
  Vector<Real> allTimers(nProc);
  gather(allTimers, a_timer, srcProc);

  Real averageTime;
  Real sigmaTime;

  Real minTime;
  Real maxTime;

  int minRank;
  int maxRank;

  // Compute and broadcast the average and the standard deviation
  if (procID() == srcProc) {

    averageTime = 0.0;
    for (int i = 0; i < nProc; i++) {
      averageTime += allTimers[i];
    }
    averageTime *= 1. / nProc;

    sigmaTime = 0.0;
    for (int i = 0; i < nProc; i++) {
      sigmaTime += std::pow(1.0 * allTimers[i] - averageTime, 2);
    }
    sigmaTime = sqrt(sigmaTime / nProc);
  }

  broadcast(averageTime, srcProc);
  broadcast(sigmaTime, srcProc);

  // Get the minimum/maximum number of particles
  minRank = srcProc;
  maxRank = srcProc;

  minTime = std::numeric_limits<Real>::max();
  maxTime = 0;

  if (procID() == srcProc) {
    for (int i = 0; i < nProc; i++) {
      if (allTimers[i] < minTime) {
        minTime = allTimers[i];
        minRank = i;
      }

      if (allTimers[i] > maxTime) {
        maxTime = allTimers[i];
        maxRank = i;
      }
    }
  }

  broadcast(minRank, srcProc);
  broadcast(maxRank, srcProc);

  broadcast(minTime, srcProc);
  broadcast(maxTime, srcProc);

  // Fix formatting for the various fields
  std::stringstream ss_locTime;
  std::stringstream ss_minTime;
  std::stringstream ss_maxTime;
  std::stringstream ss_avgTime;
  std::stringstream ss_sigTime;

  ss_locTime << std::fixed << std::setprecision(2) << a_timer;
  ss_minTime << std::fixed << std::setprecision(2) << minTime;
  ss_maxTime << std::fixed << std::setprecision(2) << maxTime;
  ss_avgTime << std::fixed << std::setprecision(2) << averageTime;
  ss_sigTime << std::fixed << std::setprecision(2) << sigmaTime;

  pout() << std::left << std::setw(25) << a_prefix << " | " << std::right << std::setw(8) << ss_locTime.str() << " | "
         << std::right << std::setw(8) << ss_minTime.str() << " | " << std::right << std::setw(8) << ss_maxTime.str()
         << " | " << std::right << std::setw(8) << ss_avgTime.str() << " | " << std::right << std::setw(8)
         << ss_sigTime.str() << " | " << std::right << std::setw(8) << minRank << " | " << std::right << std::setw(8)
         << maxRank << " | " << endl;
}

void
ItoPlasmaStepper::printTimerHead()
{
  pout() << "--------------------------------------------------------------------------------------------------------"
         << endl
         << std::left << std::setw(25) << "Kernel"
         << " | " << std::right << std::setw(8) << "Loc."
         << " | " << std::right << std::setw(8) << "Min."
         << " | " << std::right << std::setw(8) << "Max."
         << " | " << std::right << std::setw(8) << "Avg."
         << " | " << std::right << std::setw(8) << "Dev."
         << " | " << std::right << std::setw(8) << "Min rank"
         << " | " << std::right << std::setw(8) << "Max rank"
         << " | " << endl
         << "-------------------------------------------------------------------------------------------------------|"
         << endl;
}

void
ItoPlasmaStepper::printTimerTail()
{
  pout()
    << "--------------------------------------------------------------------------------------------------------\n";
}

void
ItoPlasmaStepper::parseFilters()
{
  CH_TIME("ItoPlasmaStepper::computeDt");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::parseFilters" << endl;
  }

  ParmParse pp(m_name.c_str());

  // Build filters. Always uses a compensation step.
  for (int i = 0; i < 100; i++) {

    const int ndigits = round(log10(1.0 + 1.0 * i));
    char*     cstr    = new char[1 + ndigits];
    sprintf(cstr, "%d", 1 + i);

    const std::string str = "filter_" + std::string(cstr);

    if (pp.contains(str.c_str())) {
      Real alpha;
      int  stride;
      int  N;
      bool comp;

      pp.get(str.c_str(), alpha, 0);
      pp.get(str.c_str(), stride, 1);
      pp.get(str.c_str(), N, 2);
      pp.get(str.c_str(), comp, 3);

      m_filters.emplace_front(alpha, stride, N);
      if (comp) {
        const Real alphaC = (N + 1) - N * alpha;
        m_filters.emplace_front(alphaC, stride, 1);
      }
    }

    delete cstr;
  }
}

Real
ItoPlasmaStepper::computeDt()
{
  CH_TIME("ItoPlasmaStepper::computeDt");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeDt" << endl;
  }

  return m_max_cells_hop * m_ito->computeDt();
}

void
ItoPlasmaStepper::registerRealms()
{
  CH_TIME("ItoPlasmaStepper::registerRealms");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::registerRealms" << endl;
  }

  m_amr->registerRealm(m_fluidRealm);
  m_amr->registerRealm(m_particleRealm);
}

void
ItoPlasmaStepper::registerOperators()
{
  CH_TIME("ItoPlasmaStepper::registerOperators");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::registerOperators" << endl;
  }

  m_ito->registerOperators();
  m_fieldSolver->registerOperators();
  m_rte->registerOperators();
  m_sigma->registerOperators();
}

void
ItoPlasmaStepper::preRegrid(const int a_lmin, const int a_oldFinestLevel)
{
  CH_TIME("ItoPlasmaStepper::preRegrid");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::preRegrid" << endl;
  }

  m_ito->preRegrid(a_lmin, a_oldFinestLevel);
  m_fieldSolver->preRegrid(a_lmin, a_oldFinestLevel);
  m_rte->preRegrid(a_lmin, a_oldFinestLevel);
  m_sigma->preRegrid(a_lmin, a_oldFinestLevel);
}

void
ItoPlasmaStepper::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("ItoPlasmaStepper::regrid");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::regrid" << endl;
  }

  // Allocate new memory
  this->allocateInternals();

  // Regrid solvers
  m_ito->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_fieldSolver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_rte->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_sigma->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);

  if (m_regridSuperparticles) {
    m_ito->sortParticlesByCell(ItoSolver::WhichContainer::Bulk);
    m_ito->makeSuperparticles(ItoSolver::WhichContainer::Bulk, m_particlesPerCell);
    m_ito->sortParticlesByPatch(ItoSolver::WhichContainer::Bulk);
  }

  // Redeposit particles
  m_ito->depositParticles();

  // Recompute the electric field
  const bool converged = this->solvePoisson();
  if (!converged) {
    MayDay::Abort("ItoPlasmaStepper::regrid - Poisson solve did not converge after regrid!!!");
  }

  // Recompute new velocities and diffusion coefficients
  this->computeItoVelocities();
  this->computeItoDiffusion();
}

void
ItoPlasmaStepper::postRegrid()
{}

int
ItoPlasmaStepper::getNumberOfPlotVariables() const
{
  CH_TIME("ItoPlasmaStepper::getNumberOfPlotVariables");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::getNumberOfPlotVariables" << endl;
  }

  int ncomp = 0;

  for (ItoIterator<ItoSolver> solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>& solver = solver_it();
    ncomp += solver->getNumberOfPlotVariables();
  }

  for (RtIterator<McPhoto> solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<McPhoto>& solver = solver_it();
    ncomp += solver->getNumberOfPlotVariables();
  }

  ncomp += m_fieldSolver->getNumberOfPlotVariables();
  ncomp += m_sigma->getNumberOfPlotVariables();
  ncomp += SpaceDim; // For plotting the current density
  ncomp += 1;        // For plotting the number of particles per cell

  return ncomp;
}

void
ItoPlasmaStepper::setIto(RefCountedPtr<ItoLayout<ItoSolver>>& a_ito)
{
  CH_TIME("ItoPlasmaStepper::setIto");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::setIto" << endl;
  }

  m_ito = a_ito;
}

void
ItoPlasmaStepper::setFieldSolver(RefCountedPtr<FieldSolver>& a_poisson)
{
  CH_TIME("ItoPlasmaStepper::setFieldSolver");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::setFieldSolver" << endl;
  }

  m_fieldSolver = a_poisson;
}

void
ItoPlasmaStepper::setRadiativeTransferSolvers(RefCountedPtr<RtLayout<McPhoto>>& a_rte)
{
  CH_TIME("ItoPlasmaStepper::setRadiativeTransferSolvers");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::setRadiativeTransferSolvers" << endl;
  }

  m_rte = a_rte;
}

void
ItoPlasmaStepper::setVoltage(const std::function<Real(const Real a_time)>& a_potential)
{
  CH_TIME("ItoPlasmaStepper::setVoltage");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::setVoltage" << endl;
  }

  m_potential = a_potential;
}

Real
ItoPlasmaStepper::computeMaxElectricField(const phase::which_phase a_phase)
{
  CH_TIME("ItoPlasmaStepper::computeMaxElectricField");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeMaxElectricField" << endl;
  }

  // Get a handle to the E-field
  EBAMRCellData Ephase;
  m_amr->allocatePointer(Ephase);
  m_amr->alias(Ephase, m_phase, m_fieldSolver->getElectricField());

  // Interpolate to centroids
  EBAMRCellData E;
  m_amr->allocate(E, m_fluidRealm, m_phase, SpaceDim);
  DataOps::copy(E, Ephase);
  m_amr->interpToCentroids(E, m_fluidRealm, m_phase);

  Real max, min;
  DataOps::getMaxMinNorm(max, min, E);

  return max;
}

Real
ItoPlasmaStepper::getTime() const
{
  CH_TIME("ItoPlasmaStepper::getTime");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::getTime" << endl;
  }

  return m_time;
}

void
ItoPlasmaStepper::computeElectricField(MFAMRCellData& a_E, const MFAMRCellData& a_potential)
{
  CH_TIME("ItoPlasmaStepper::computeElectricField(mfamrcell,mfamrcell)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeElectricField(mfamrcell, mfamrcell" << endl;
  }

  m_fieldSolver->computeElectricField(a_E, a_potential);

  m_amr->conservativeAverage(a_E, m_fluidRealm);
  m_amr->interpGhost(a_E, m_fluidRealm);
}

void
ItoPlasmaStepper::computeElectricField(EBAMRCellData& a_E, const phase::which_phase a_phase)
{
  CH_TIME("ItoPlasmaStepper::computeElectricField(ebamrcell, phase)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeElectricField(ebamrcell, phase" << endl;
  }

  this->computeElectricField(a_E, a_phase, m_fieldSolver->getPotential());
}

void
ItoPlasmaStepper::computeElectricField(EBAMRCellData&           a_E,
                                       const phase::which_phase a_phase,
                                       const MFAMRCellData&     a_potential)
{
  CH_TIME("ItoPlasmaStepper::computeElectricField(ebamrcell, phase, mfamrcell)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeElectricField(ebamrcell, phase mfamrcell" << endl;
  }

  m_fieldSolver->computeElectricField(a_E, a_phase, a_potential);

  m_amr->conservativeAverage(a_E, m_fluidRealm, a_phase);
  m_amr->interpGhost(a_E, m_fluidRealm, a_phase);
}

void
ItoPlasmaStepper::computeElectricField(EBAMRFluxData&           a_E_face,
                                       const phase::which_phase a_phase,
                                       const EBAMRCellData&     a_E_cell)
{
  CH_TIME("ItoPlasmaStepper::computeElectricField(ebamrflux, phase, mfamrcell)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeElectricField(ebamrflux, phase mfamrcell" << endl;
  }

  CH_assert(a_E_face[0]->nComp() == SpaceDim);
  CH_assert(a_E_cell[0]->nComp() == SpaceDim);

  const int finest_level = m_amr->getFinestLevel();
  for (int lvl = 0; lvl <= finest_level; lvl++) {

    const DisjointBoxLayout& dbl    = m_amr->getGrids(m_fluidRealm)[lvl];
    const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_fluidRealm, a_phase)[lvl];
    const ProblemDomain&     domain = m_amr->getDomains()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      const EBCellFAB& E_cell  = (*a_E_cell[lvl])[dit()];
      const EBISBox&   ebisbox = ebisl[dit()];
      const EBGraph&   ebgraph = ebisbox.getEBGraph();
      const Box&       box     = dbl.get(dit());

      for (int dir = 0; dir < SpaceDim; dir++) {
        EBFaceFAB& E_face = (*a_E_face[lvl])[dit()][dir];
        E_face.setVal(0.0);

        EBLevelDataOps::averageCellToFace(E_face, E_cell, ebgraph, box, 0, dir, domain, dir, dir);
      }
    }
    a_E_face[lvl]->exchange();
  }
}

void
ItoPlasmaStepper::computeElectricField(EBAMRIVData&             a_E_eb,
                                       const phase::which_phase a_phase,
                                       const EBAMRCellData&     a_E_cell)
{
  CH_TIME("ItoPlasmaStepper::computeElectricField(ebamriv, phase, ebamrcell)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeElectricField(ebamriv, phase ebamrcell)" << endl;
  }

  CH_assert(a_E_eb[0]->nComp() == SpaceDim);
  CH_assert(a_E_cell[0]->nComp() == SpaceDim);

  const IrregAmrStencil<EbCentroidInterpolationStencil>& interp_stencil =
    m_amr->getEbCentroidInterpolationStencils(m_fluidRealm, a_phase);
  interp_stencil.apply(a_E_eb, a_E_cell);
}

void
ItoPlasmaStepper::computeSpaceChargeDensity()
{
  CH_TIME("ItoPlasmaStepper::computeSpaceChargeDensity()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeSpaceChargeDensity()" << endl;
  }

  this->computeSpaceChargeDensity(m_fieldSolver->getRho(), m_ito->getDensities());
}

void
ItoPlasmaStepper::computeSpaceChargeDensity(MFAMRCellData& a_rho, const Vector<EBAMRCellData*>& a_densities)
{
  CH_TIME("ItoPlasmaStepper::computeSpaceChargeDensity(rho, densities)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeSpaceChargeDensity(rho, densities)" << endl;
  }

  // TLDR: a_densities is from the ito solvers so it is defined over the particle Realm. But a_rho is defined over
  //       the fluid Realm so we need scratch storage we can copy into. We use m_fluid_scratch1 for that.

  // Reset
  DataOps::setValue(a_rho, 0.0);

  // Make alias
  EBAMRCellData rhoPhase;
  m_amr->allocatePointer(rhoPhase);
  m_amr->alias(rhoPhase, m_phase, a_rho);

  // Increment each solver
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    const RefCountedPtr<ItoSolver>&  solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();
    const int                        idx     = solver_it.index();
    const int                        q       = species->getChargeNumber();

    if (species->getChargeNumber() != 0) {
      m_fluid_scratch1.copy(*a_densities[idx]);
      DataOps::incr(rhoPhase, m_fluid_scratch1, q);
    }
  }

  DataOps::scale(a_rho, Units::Qe);

  m_amr->conservativeAverage(a_rho, m_fluidRealm);
  m_amr->interpGhost(a_rho, m_fluidRealm);

  // Interpolate to centroids
  m_amr->interpToCentroids(rhoPhase, m_fluidRealm, m_phase);
}

void
ItoPlasmaStepper::computeConductivity(EBAMRCellData& a_conductivity)
{
  CH_TIME("ItoPlasmaStepper::computeConductivity(conductivity)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeConductivity(conductivity)" << endl;
  }

  this->computeConductivity(a_conductivity, m_ito->getParticles(ItoSolver::WhichContainer::Bulk));
}

void
ItoPlasmaStepper::computeConductivity(EBAMRCellData&                                 a_conductivity,
                                      const Vector<ParticleContainer<ItoParticle>*>& a_particles)
{
  CH_TIME("ItoPlasmaStepper::computeConductivity(conductivity)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeConductivity(conductivity)" << endl;
  }

  DataOps::setValue(a_conductivity, 0.0);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();
    const int q   = species->getChargeNumber();

    if (Abs(q) > 0 && solver->isMobile()) {
      DataOps::setValue(m_particle_scratch1, 0.0);

      solver->depositConductivity(m_particle_scratch1, *a_particles[idx]);

      // Copy to fluid Realm and add to fluid stuff
      m_fluid_scratch1.copy(m_particle_scratch1);
      DataOps::incr(a_conductivity, m_fluid_scratch1, Abs(q));
    }
  }

  DataOps::scale(a_conductivity, Units::Qe);

  m_amr->conservativeAverage(a_conductivity, m_fluidRealm, m_phase);
  m_amr->interpGhostPwl(a_conductivity, m_fluidRealm, m_phase);

  // See if this helps....
  m_amr->interpToCentroids(a_conductivity, m_fluidRealm, m_phase);
}

void
ItoPlasmaStepper::computeJ(EBAMRCellData& a_J, const Real a_dt)
{
  CH_TIME("ItoPlasmaStepper::computeJ(J)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeJ(J)" << endl;
  }

  // TLDR: a_J is defined over the fluid Realm but the computation takes place on the particle Realm.
  //       If the Realms are different we compute on a scratch storage instead

  this->computeConductivity(m_fluid_scratch1);
  DataOps::copy(a_J, m_fluid_E);

  DataOps::multiplyScalar(a_J, m_fluid_scratch1);
}

Real
ItoPlasmaStepper::computeRelaxationTime()
{
  CH_TIME("ItoPlasmaStepper::computeRelaxationTime()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeRelaxationTime()" << endl;
  }

  Real dt = 1.E99;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Real thisDt = this->computeRelaxationTime(lvl);

    dt = Min(dt, thisDt);
  }

  return dt;
}

Real
ItoPlasmaStepper::computeRelaxationTime(const int a_level)
{
  CH_TIME("ItoPlasmaStepper::computeRelaxationTime(level)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeRelaxationTime(level)" << endl;
  }

  Real dt = 1.E99;

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_fluidRealm)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
    const Real thisDt = this->computeRelaxationTime(a_level, dit());

    dt = Min(dt, thisDt);
  }

#ifdef CH_MPI
  Real tmp    = dt;
  int  result = MPI_Allreduce(&dt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if (result != MPI_SUCCESS) {
    MayDay::Error("CdrSolver::compute_cfl_dt() - communication error on norm");
  }
  dt = tmp;
#endif

  return dt;
}

Real
ItoPlasmaStepper::computeRelaxationTime(const int a_level, const DataIndex a_dit)
{
  CH_TIME("ItoPlasmaStepper::computeRelaxationTime(level, dit)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeRelaxationTime(level, dit)" << endl;
  }

  const int  comp   = 0;
  const Real SAFETY = 1.E-10;

  const Box      box     = m_amr->getGrids(m_fluidRealm)[a_level].get(a_dit);
  const EBISBox& ebisbox = m_amr->getEBISLayout(m_fluidRealm, m_phase)[a_level][a_dit];

  // Get a handle to the E-field
  EBAMRCellData amrE;
  m_amr->allocatePointer(amrE);
  m_amr->alias(amrE, m_phase, m_fieldSolver->getElectricField());

  const EBCellFAB& E = (*amrE[a_level])[a_dit];
  const EBCellFAB& J = (*m_J[a_level])[a_dit];

  EBCellFAB dt(ebisbox, box, 1);
  EBCellFAB e_magnitude(ebisbox, box, 1);
  EBCellFAB j_magnitude(ebisbox, box, 1);

  e_magnitude.setVal(0.0);
  j_magnitude.setVal(0.0);

  DataOps::vectorLength(e_magnitude, E, box);
  DataOps::vectorLength(j_magnitude, J, box);
  j_magnitude += SAFETY;

  dt.setVal(Units::eps0);
  dt *= e_magnitude;
  dt /= j_magnitude;

  return dt.min(comp);
}

bool
ItoPlasmaStepper::solvePoisson()
{
  CH_TIME("ItoPlasmaStepper::solvePoisson()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::solvePoisson()" << endl;
  }

  // Computes rho
  this->computeSpaceChargeDensity();

  // Solves the Poisson equation
  const bool converged =
    m_fieldSolver->solve(m_fieldSolver->getPotential(), m_fieldSolver->getRho(), m_sigma->getPhi(), false);

  // Computes cell-centered E onto storage in the field solver
  m_fieldSolver->computeElectricField();

  // Code below here interpolates E to centroids on both Realms
  EBAMRCellData E;
  m_amr->allocatePointer(E);
  m_amr->alias(E, m_phase, m_fieldSolver->getElectricField());

  // Fluid Realm
  m_fluid_E.copy(E);
  m_amr->conservativeAverage(m_fluid_E, m_fluidRealm, m_phase);
  m_amr->interpGhostPwl(m_fluid_E, m_fluidRealm, m_phase);
  m_amr->interpToCentroids(m_fluid_E, m_fluidRealm, m_phase);

  // Particle Realm
  m_particle_E.copy(E);
  m_amr->conservativeAverage(m_particle_E, m_particleRealm, m_phase);
  m_amr->interpGhostPwl(m_particle_E, m_particleRealm, m_phase);
  m_amr->interpToCentroids(m_particle_E, m_particleRealm, m_phase);

  return converged;
}

bool
ItoPlasmaStepper::solvePoisson(MFAMRCellData&                a_potential,
                               MFAMRCellData&                a_rho,
                               const Vector<EBAMRCellData*>& a_densities,
                               const EBAMRIVData&            a_sigma)
{
  CH_TIME("ItoPlasmaStepper::solvePoisson(phi, rho, densities, sigma)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::solvePoisson(phi, rho, densities, sigma)" << endl;
  }

  this->computeSpaceChargeDensity(a_rho, a_densities);

  const bool converged = m_fieldSolver->solve(a_potential, a_rho, a_sigma, false);
  m_fieldSolver->computeElectricField();

  return converged;
}

void
ItoPlasmaStepper::intersectParticles(const WhichParticles   a_WhichParticles,
                                     const EbRepresentation a_representation,
                                     const bool             a_delete)
{
  CH_TIME("ItoPlasmaStepper::intersectParticles(WhichParticles, EbRepresentation)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::intersectParticles(WhichParticles, EbRepresentation)" << endl;
  }

  this->intersectParticles(a_WhichParticles,
                           ItoSolver::WhichContainer::Bulk,
                           ItoSolver::WhichContainer::EB,
                           ItoSolver::WhichContainer::Domain,
                           a_representation,
                           a_delete);
}

void
ItoPlasmaStepper::intersectParticles(const WhichParticles            a_WhichParticles,
                                     const ItoSolver::WhichContainer a_particles,
                                     const ItoSolver::WhichContainer a_eb_particles,
                                     const ItoSolver::WhichContainer a_domain_particles,
                                     const EbRepresentation          a_representation,
                                     const bool                      a_delete)
{
  CH_TIME("ItoPlasmaStepper::intersectParticles(WhichParticles, string, string, string, EbRepresentation, bool)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::intersectParticles(WhichParticles, string, string, string, EbRepresentation, bool)"
           << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = species->getChargeNumber() != 0;
#if 0
    switch (a_WhichParticles) {
    case WhichParticles::All:
      solver->intersectParticles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case WhichParticles::AllMobile:
      if (mobile)
        solver->intersectParticles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case WhichParticles::AllDiffusive:
      if (diffusive)
        solver->intersectParticles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case WhichParticles::ChargedMobile:
      if (charged && mobile)
        solver->intersectParticles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case WhichParticles::ChargedDiffusive:
      if (charged && diffusive)
        solver->intersectParticles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case WhichParticles::AllMobileOrDiffusive:
      if (mobile || diffusive)
        solver->intersectParticles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case WhichParticles::ChargedAndMobileOrDiffusive:
      if (charged && (mobile || diffusive))
        solver->intersectParticles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    case WhichParticles::Stationary:
      if (!mobile && !diffusive)
        solver->intersectParticles(a_particles, a_eb_particles, a_domain_particles, a_representation, a_delete);
      break;
    default:
      MayDay::Abort(
        "ItoPlasmaStepper::intersectParticles_particles(WhichParticles, string, string, string, EbRepresentation, bool) - logic bust");
    }
#endif
  }
}

void
ItoPlasmaStepper::removeCoveredParticles(const WhichParticles   a_WhichParticles,
                                         const EbRepresentation a_representation,
                                         const Real             a_tolerance)
{
  CH_TIME("ItoPlasmaStepper::removeCoveredParticles(WhichParticles, representation, tolerance)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::removeCoveredParticles(WhichParticles, representation, tolerance)" << endl;
  }

  this->removeCoveredParticles(a_WhichParticles, ItoSolver::WhichContainer::Bulk, a_representation, a_tolerance);
}

void
ItoPlasmaStepper::removeCoveredParticles(const WhichParticles            a_which,
                                         const ItoSolver::WhichContainer a_container,
                                         const EbRepresentation          a_representation,
                                         const Real                      a_tolerance)
{
  CH_TIME("ItoPlasmaStepper::removeCoveredParticles(WhichParticles, container, EbRepresentation, tolerance)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::removeCoveredParticles(WhichParticles, container, EbRepresentation, tolerance)"
           << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = species->getChargeNumber() != 0;

    switch (a_which) {
    case WhichParticles::All:
      solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      break;
    case WhichParticles::AllMobile:
      if (mobile)
        solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      break;
    case WhichParticles::AllDiffusive:
      if (diffusive)
        solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      break;
    case WhichParticles::ChargedMobile:
      if (charged && mobile)
        solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      break;
    case WhichParticles::ChargedDiffusive:
      if (charged && diffusive)
        solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      break;
    case WhichParticles::AllMobileOrDiffusive:
      if (mobile || diffusive)
        solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      break;
    case WhichParticles::ChargedAndMobileOrDiffusive:
      if (charged && (mobile || diffusive))
        solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      break;
    case WhichParticles::Stationary:
      if (!mobile && !diffusive)
        solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      break;
    default:
      MayDay::Abort("ItoPlasmaStepper::removeCoveredParticles_particles(which particles) - logic bust");
    }
  }
}

void
ItoPlasmaStepper::transferCoveredParticles(const WhichParticles   a_which,
                                           const EbRepresentation a_representation,
                                           const Real             a_tolerance)
{
  CH_TIME("ItoPlasmaStepper::transferCoveredParticles_particles(WhichParticles, EbRepresentation, tolerance)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::transferCoveredParticles_particles(WhichParticles, EbRepresentation, tolerance)"
           << endl;
  }

  this->transferCoveredParticles(a_which,
                                 ItoSolver::WhichContainer::Bulk,
                                 ItoSolver::WhichContainer::Covered,
                                 a_representation,
                                 a_tolerance);
}

void
ItoPlasmaStepper::transferCoveredParticles(const WhichParticles            a_which,
                                           const ItoSolver::WhichContainer a_containerFrom,
                                           const ItoSolver::WhichContainer a_containerTo,
                                           const EbRepresentation          a_representation,
                                           const Real                      a_tolerance)
{
  CH_TIME(
    "ItoPlasmaStepper::transferCoveredParticles_particles(WhichParticles, container, container, EbRepresentation, tolerance)");
  if (m_verbosity > 5) {
    pout()
      << "ItoPlasmaStepper::transferCoveredParticles_particles(WhichParticles, container, container, EbRepresentation, tolerance)"
      << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = species->getChargeNumber() != 0;

    switch (a_which) {
    case WhichParticles::All:
      solver->transferCoveredParticles(a_containerFrom, a_containerTo, a_representation, a_tolerance);
      break;
    case WhichParticles::AllMobile:
      if (mobile)
        solver->transferCoveredParticles(a_containerFrom, a_containerTo, a_representation, a_tolerance);
      break;
    case WhichParticles::AllDiffusive:
      if (diffusive)
        solver->transferCoveredParticles(a_containerFrom, a_containerTo, a_representation, a_tolerance);
      break;
    case WhichParticles::ChargedMobile:
      if (charged && mobile)
        solver->transferCoveredParticles(a_containerFrom, a_containerTo, a_representation, a_tolerance);
      break;
    case WhichParticles::ChargedDiffusive:
      if (charged && diffusive)
        solver->transferCoveredParticles(a_containerFrom, a_containerTo, a_representation, a_tolerance);
      break;
    case WhichParticles::AllMobileOrDiffusive:
      if (mobile || diffusive)
        solver->transferCoveredParticles(a_containerFrom, a_containerTo, a_representation, a_tolerance);
      break;
    case WhichParticles::ChargedAndMobileOrDiffusive:
      if (charged && (mobile || diffusive))
        solver->transferCoveredParticles(a_containerFrom, a_containerTo, a_representation, a_tolerance);
      break;
    case WhichParticles::Stationary:
      if (!mobile && !diffusive)
        solver->transferCoveredParticles(a_containerFrom, a_containerTo, a_representation, a_tolerance);
      break;
    default:
      MayDay::Abort("ItoPlasmaStepper::transferCoveredParticles_particles(...) - logic bust");
    }
  }
}

void
ItoPlasmaStepper::remapParticles(const WhichParticles a_WhichParticles)
{
  CH_TIME("ItoPlasmaStepper::remapParticles(WhichParticles)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::remapParticles(WhichParticles)" << endl;
  }

  this->remapParticles(a_WhichParticles, ItoSolver::WhichContainer::Bulk);
}

void
ItoPlasmaStepper::remapParticles(const WhichParticles a_WhichParticles, const ItoSolver::WhichContainer a_container)
{
  CH_TIME("ItoPlasmaStepper::remapParticles(WhichParticles)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::remapParticles(WhichParticles)" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = species->getChargeNumber() != 0;

    switch (a_WhichParticles) {
    case WhichParticles::All:
      solver->remap(a_container);
      break;
    case WhichParticles::AllMobile:
      if (mobile)
        solver->remap(a_container);
      break;
    case WhichParticles::AllDiffusive:
      if (diffusive)
        solver->remap(a_container);
      break;
    case WhichParticles::ChargedMobile:
      if (charged && mobile)
        solver->remap(a_container);
      break;
    case WhichParticles::ChargedDiffusive:
      if (charged && diffusive)
        solver->remap(a_container);
      break;
    case WhichParticles::AllMobileOrDiffusive:
      if (mobile || diffusive)
        solver->remap(a_container);
      break;
    case WhichParticles::ChargedAndMobileOrDiffusive:
      if (charged && (mobile || diffusive))
        solver->remap(a_container);
      break;
    case WhichParticles::Stationary:
      if (!mobile && !diffusive)
        solver->remap(a_container);
      break;
    default:
      MayDay::Abort("ItoPlasmaStepper::remapParticles(which particles) - logic bust");
    }
  }
}

void
ItoPlasmaStepper::depositParticles(const WhichParticles a_WhichParticles)
{
  CH_TIME("ItoPlasmaStepper::depositParticles(WhichParticles)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::depositParticles(WhichParticles)" << endl;
  }

  this->depositParticles(a_WhichParticles, ItoSolver::WhichContainer::Bulk);
}

void
ItoPlasmaStepper::depositParticles(const WhichParticles a_WhichParticles, const ItoSolver::WhichContainer a_container)
{
  CH_TIME("ItoPlasmaStepper::depositParticles(WhichParticles)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::depositParticles(WhichParticles)" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = species->getChargeNumber() != 0;

    switch (a_WhichParticles) {
    case WhichParticles::All:
      solver->depositParticles(a_container);
      break;
    case WhichParticles::AllMobile:
      if (mobile)
        solver->depositParticles(a_container);
      break;
    case WhichParticles::AllDiffusive:
      if (diffusive)
        solver->depositParticles(a_container);
      break;
    case WhichParticles::ChargedMobile:
      if (charged && mobile)
        solver->depositParticles(a_container);
      break;
    case WhichParticles::ChargedDiffusive:
      if (charged && diffusive)
        solver->depositParticles(a_container);
      break;
    case WhichParticles::AllMobileOrDiffusive:
      if (mobile || diffusive)
        solver->depositParticles(a_container);
      break;
    case WhichParticles::ChargedAndMobileOrDiffusive:
      if (charged && (mobile || diffusive))
        solver->depositParticles(a_container);
      break;
    case WhichParticles::Stationary:
      if (!mobile && !diffusive)
        solver->depositParticles(a_container);
      break;
    default:
      MayDay::Abort("ItoPlasmaStepper::depositParticles(WhichParticles) - logic bust");
    }
  }
}

void
ItoPlasmaStepper::setItoVelocityFunctions()
{
  CH_TIME("ItoPlasmaStepper::setItoVelocityFunctions");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::setItoVelocityFunctions" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    if (solver->isMobile()) {
      EBAMRCellData& velo_func = solver->getVelocityFunction();
      velo_func.copy(m_particle_E);

      const int q = species->getChargeNumber();
      const int s = (q > 0) - (q < 0);

      DataOps::scale(velo_func, s);
    }
  }
}

void
ItoPlasmaStepper::computeItoVelocities()
{
  CH_TIME("ItoPlasmaStepper::computeItoVelocities()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoVelocities()" << endl;
  }

  const ItoPlasmaPhysics::coupling which_coupling = m_physics->getCoupling();

  // Set velocity functions
  this->setItoVelocityFunctions();

  // Compute mobilities based on appropriate coupling
  switch (which_coupling) {
  case ItoPlasmaPhysics::coupling::LFA:
    this->computeItoMobilitiesLFA();
    break;
  case ItoPlasmaPhysics::coupling::LEA:
    this->computeItoMobilitiesLEA();
    break;
  }

  // Interpolate velocity function to particle position
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->interpolateVelocities(); // Interpolates v = +/- mu*E
  }
}

void
ItoPlasmaStepper::computeItoDiffusion()
{
  CH_TIME("ItoPlasmaStepper::computeItoDiffusion()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoDiffusion()" << endl;
  }

  const ItoPlasmaPhysics::coupling which_coupling = m_physics->getCoupling();

  // Compute mobilities based on appropriate coupling
  switch (which_coupling) {
  case ItoPlasmaPhysics::coupling::LFA:
    this->computeItoDiffusionLFA();
    break;
  case ItoPlasmaPhysics::coupling::LEA:
    this->computeItoDiffusion_lea();
    break;
  }
}

void
ItoPlasmaStepper::computeItoMobilitiesLFA()
{
  CH_TIME("ItoPlasmaStepper::computeItoMobilitiesLFA()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoMobilitiesLFA()" << endl;
  }

  Vector<EBAMRCellData*> meshMobilities = m_ito->getMobilityFunctions();
  this->computeItoMobilitiesLFA(meshMobilities, m_fluid_E, m_time);
}

void
ItoPlasmaStepper::computeItoMobilitiesLFA(Vector<EBAMRCellData*>& a_meshMobilities,
                                          const EBAMRCellData&    a_E,
                                          const Real              a_time)
{
  CH_TIME("ItoPlasmaStepper::computeItoMobilitiesLFA(mobilities, E, time)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoMobilitiesLFA(mobilities, E, time)" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {

    // Computation is done on the fluid Realm.
    Vector<LevelData<EBCellFAB>*> meshMobilities;
    for (int i = 0; i < a_meshMobilities.size(); i++) {
      meshMobilities.push_back(&(*(m_fscratch1[i])[lvl]));
    }

    this->computeItoMobilitiesLFA(meshMobilities, *a_E[lvl], lvl, a_time);
  }

  // Average down and interpolate ghost cells. Then interpolate mobilities to particle positions.
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    const int                 idx    = solver_it.index();
    RefCountedPtr<ItoSolver>& solver = solver_it();

    if (solver->isMobile()) {

#if 0 // In principle, we should be able to average down and interpolate on the fluid Realm and then copy directly to the particle Realm. \
      // But we need to make sure that EBAMRData::copy also gets ghost cells 
      m_amr->conservativeAverage(m_fscratch1[idx], m_fluidRealm, m_phase);
      m_amr->interpGhost(m_fscratch1[idx], m_fluidRealm, m_phase);

      a_meshMobilities[idx]->copy(m_fscratch1[idx]);
#else
      // Copy to particle Realm, build ghost cells and the interpolate the mobilities to particle positions.
      a_meshMobilities[idx]->copy(m_fscratch1[idx]);

      m_amr->conservativeAverage(*a_meshMobilities[idx], m_particleRealm, m_phase);
      m_amr->interpGhost(*a_meshMobilities[idx], m_particleRealm, m_phase);
#endif

      solver->interpolateMobilities();
    }
  }
}

void
ItoPlasmaStepper::computeItoMobilitiesLFA(Vector<LevelData<EBCellFAB>*>& a_meshMobilities,
                                          const LevelData<EBCellFAB>&    a_E,
                                          const int                      a_level,
                                          const Real                     a_time)
{
  CH_TIME("ItoPlasmaStepper::computeItoMobilitiesLFA(mobilities, E, level, time)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoMobilitiesLFA(mobilities, E, level, time)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_fluidRealm)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
    const EBCellFAB& E  = a_E[dit()];
    const Box        bx = dbl.get(dit());

    Vector<EBCellFAB*> meshMobilities;
    for (int i = 0; i < a_meshMobilities.size(); i++) {
      meshMobilities.push_back(&((*a_meshMobilities[i])[dit()]));
    }

    this->computeItoMobilitiesLFA(meshMobilities, E, a_level, dit(), bx, a_time);
  }
}

void
ItoPlasmaStepper::computeItoMobilitiesLFA(Vector<EBCellFAB*>& a_meshMobilities,
                                          const EBCellFAB&    a_E,
                                          const int           a_level,
                                          const DataIndex     a_dit,
                                          const Box           a_box,
                                          const Real          a_time)
{
  CH_TIME("ItoPlasmaStepper::computeItoMobilitiesLFA(meshMobilities, E, level, dit, box, time)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoMobilitiesLFA(meshMobilities, E, level, dit, box, time)" << endl;
  }

  const int            comp    = 0;
  const Real           dx      = m_amr->getDx()[a_level];
  const RealVect       prob_lo = m_amr->getProbLo();
  const BaseFab<Real>& E       = a_E.getSingleValuedFAB();

  // Do regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit) {
    const IntVect  iv  = bit();
    const RealVect pos = m_amr->getProbLo() + dx * (RealVect(iv) + 0.5 * RealVect::Unit);
    const RealVect e   = RealVect(D_DECL(E(iv, 0), E(iv, 1), E(iv, 2)));

    // Call ito_physics and compute diffusion for each particle species
    const Vector<Real> mobilities = m_physics->computeItoMobilitiesLFA(a_time, pos, e);

    // Put mobilities in data holder
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
      const int idx                                           = solver_it.index();
      (*a_meshMobilities[idx]).getSingleValuedFAB()(iv, comp) = mobilities[idx];
    }
  }

  // Do irregular cells
  VoFIterator& vofit = (*m_amr->getVofIterator(m_fluidRealm, m_phase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit) {
    const VolIndex& vof = vofit();
    const RealVect  e   = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));
    const RealVect  pos = EBArith::getVofLocation(vof, dx * RealVect::Unit, prob_lo);

    // Compute diffusion
    const Vector<Real> mobilities = m_physics->computeItoMobilitiesLFA(a_time, pos, e);

    // Put diffusion in the appropriate place.
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
      const int idx                       = solver_it.index();
      (*a_meshMobilities[idx])(vof, comp) = mobilities[idx];
    }
  }

  // Covered is bogus.
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    const int idx = solver_it.index();
    a_meshMobilities[idx]->setCoveredCellVal(0.0, comp);
  }
}

void
ItoPlasmaStepper::computeItoMobilitiesLEA()
{
  CH_TIME("ItoPlasmaStepper");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoMobilitiesLEA()" << endl;
  }

  // This is really simple because the solvers do this directly...
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->updateMobilities();
  }
}

void
ItoPlasmaStepper::computeItoDiffusionLFA()
{
  CH_TIME("ItoPlasmaStepper::computeItoDiffusionLFA()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoDiffusionLFA()" << endl;
  }

  Vector<EBAMRCellData*> diffco_funcs = m_ito->getDiffusionFunctions();
  Vector<EBAMRCellData*> densities    = m_ito->getDensities();

  this->computeItoDiffusionLFA(diffco_funcs, densities, m_particle_E, m_time);
}

void
ItoPlasmaStepper::computeItoDiffusionLFA(Vector<EBAMRCellData*>&       a_diffusionCoefficient_funcs,
                                         const Vector<EBAMRCellData*>& a_densities,
                                         const EBAMRCellData&          a_E,
                                         const Real                    a_time)
{
  CH_TIME("ItoPlasmaStepper::computeItoDiffusionLFA(velo, E, time)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoDiffusionLFA(velo, E, time)" << endl;
  }

  // TLDR: In this routine we make m_fscratch1 hold the diffusion coefficients on the fluid Realm and m_fscratch2 hold the particle densities on the fluid Realm.
  //       This requires a couple of copies.

  // 1. Copy particle Realm densities to fluid Realm scratch data.
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    const int idx = solver_it.index();

    m_fscratch2[idx].copy(*a_densities[idx]);
  }

  // 2. Compute on each level. On the fluid Realm.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {

    const int num_ItoSpecies = m_physics->getNumItoSpecies();

    Vector<LevelData<EBCellFAB>*> diffco_funcs(num_ItoSpecies);
    Vector<LevelData<EBCellFAB>*> densities(num_ItoSpecies);

    for (int idx = 0; idx < a_diffusionCoefficient_funcs.size(); idx++) {
      diffco_funcs[idx] = &(*(m_fscratch1[idx])[lvl]);
      densities[idx]    = &(*(m_fscratch2[idx])[lvl]);
    }

    this->computeItoDiffusionLFA(diffco_funcs, densities, *m_fluid_E[lvl], lvl, a_time);
  }

  // Average down, interpolate ghost cells, and then interpolate to particle positions
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    const int                 idx    = solver_it.index();
    RefCountedPtr<ItoSolver>& solver = solver_it();

    if (solver->isDiffusive()) {

#if 0 // In principle, we should be able to average down and interpolate ghost cells on the fluid Realm, and copy the entire result over to the particle Realm.
      m_amr->conservativeAverage(m_fscratch1[idx], m_fluidRealm, m_phase);
      m_amr->interpGhost(m_fscratch2[idx], m_fluidRealm, m_phase);
      a_diffusionCoefficient_funcs[idx]->copy(m_fluid_scratch1[idx]);
#else // Instead, we copy to the particle Realm and average down there, then interpolate.
      a_diffusionCoefficient_funcs[idx]->copy(m_fscratch1[idx]);

      m_amr->conservativeAverage(*a_diffusionCoefficient_funcs[idx], m_particleRealm, m_phase);
      m_amr->interpGhost(*a_diffusionCoefficient_funcs[idx], m_particleRealm, m_phase);
#endif

      solver->interpolateDiffusion();
    }
  }
}

void
ItoPlasmaStepper::computeItoDiffusionLFA(Vector<LevelData<EBCellFAB>*>&       a_diffusionCoefficient,
                                         const Vector<LevelData<EBCellFAB>*>& a_densities,
                                         const LevelData<EBCellFAB>&          a_E,
                                         const int                            a_level,
                                         const Real                           a_time)
{
  CH_TIME("ItoPlasmaStepper::computeItoDiffusionLFA(velo, E, level, time)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoDiffusionLFA(velo, E, level, time)" << endl;
  }

  const int num_ItoSpecies = m_physics->getNumItoSpecies();

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_fluidRealm)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
    const Box box = dbl.get(dit());

    Vector<EBCellFAB*> diffusion(num_ItoSpecies);
    Vector<EBCellFAB*> densities(num_ItoSpecies);
    ;

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
      const int idx = solver_it.index();

      if (solver_it()->isDiffusive()) {
        diffusion[idx] = &(*a_diffusionCoefficient[idx])[dit()];
      }
      densities[idx] = &(*a_densities[idx])[dit()];
    }

    this->computeItoDiffusionLFA(diffusion, densities, a_E[dit()], a_level, dit(), box, a_time);
  }
}

void
ItoPlasmaStepper::computeItoDiffusionLFA(Vector<EBCellFAB*>&       a_diffusionCoefficient,
                                         const Vector<EBCellFAB*>& a_densities,
                                         const EBCellFAB&          a_E,
                                         const int                 a_level,
                                         const DataIndex           a_dit,
                                         const Box                 a_box,
                                         const Real                a_time)
{
  CH_TIME("ItoPlasmaStepper::computeItoDiffusionLFA(velo, E, level, dit, time)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoDiffusionLFA(velo, E, level, dit, time)" << endl;
  }

  const int            comp    = 0;
  const Real           dx      = m_amr->getDx()[a_level];
  const RealVect       prob_lo = m_amr->getProbLo();
  const BaseFab<Real>& E       = a_E.getSingleValuedFAB();

  // Do regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit) {
    const IntVect  iv  = bit();
    const RealVect pos = m_amr->getProbLo() + dx * (RealVect(iv) + 0.5 * RealVect::Unit);
    const RealVect e   = RealVect(D_DECL(E(iv, 0), E(iv, 1), E(iv, 2)));

    // Make grid densities
    Vector<Real> densities;
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
      const int idx = solver_it.index();
      densities.push_back((*a_densities[idx]).getSingleValuedFAB()(iv, comp));
    }

    // Call ito_physics and compute diffusion for each particle species
    const Vector<Real> diffusion = m_physics->computeItoDiffusionLFA(a_time, pos, e, densities);

    // Put diffusion where they belong
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
      RefCountedPtr<ItoSolver>& solver = solver_it();
      if (solver->isDiffusive()) {
        const int idx                                                 = solver_it.index();
        (*a_diffusionCoefficient[idx]).getSingleValuedFAB()(iv, comp) = diffusion[idx];
      }
    }
  }

  // Do irregular cells
  VoFIterator& vofit = (*m_amr->getVofIterator(m_fluidRealm, m_phase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit) {
    const VolIndex& vof = vofit();
    const RealVect  e   = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));
    const RealVect  pos = EBArith::getVofLocation(vof, dx * RealVect::Unit, prob_lo);

    // Get densities
    Vector<Real> densities;
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
      const int idx = solver_it.index();
      densities.push_back((*a_densities[idx])(vof, comp));
    }

    // Compute diffusion
    const Vector<Real> diffusion = m_physics->computeItoDiffusionLFA(a_time, pos, e, densities);

    // Put diffusion in the appropriate place.
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
      if (solver_it()->isDiffusive()) {
        const int idx                             = solver_it.index();
        (*a_diffusionCoefficient[idx])(vof, comp) = diffusion[idx];
      }
    }
  }

  // Covered is bogus.
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    if (solver_it()->isDiffusive()) {
      const int idx = solver_it.index();
      a_diffusionCoefficient[idx]->setCoveredCellVal(0.0, comp);
    }
  }
}

void
ItoPlasmaStepper::computeItoDiffusion_lea()
{
  CH_TIME("ItoPlasmaStepper::computeItoDiffusion_lea()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoDiffusion_lea()" << endl;
  }

  // This is really simple because the solvers do this directly... No monkeying with interpolations or anything.
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->updateDiffusion();
  }
}

void
ItoPlasmaStepper::computeReactiveParticlesPerCell(EBAMRCellData& a_ppc)
{
  CH_TIME("ItoPlasmaStepper::computeReactiveParticlesPerCell(ppc)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeReactiveParticlesPerCell(ppc)" << endl;
  }

  DataOps::setValue(a_ppc, 0.0);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    this->computeReactiveParticlesPerCell(*a_ppc[lvl], lvl);
  }
}

void
ItoPlasmaStepper::computeReactiveParticlesPerCell(LevelData<EBCellFAB>& a_ppc, const int a_level)
{
  CH_TIME("ItoPlasmaStepper::computeReactiveParticlesPerCell(ppc, lvl)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeReactiveParticlesPerCell(ppc, lvl)" << endl;
  }

  const DisjointBoxLayout& dbl   = m_amr->getGrids(m_particleRealm)[a_level];
  const EBISLayout&        ebisl = m_amr->getEBISLayout(m_particleRealm, m_phase)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

    const Box      box     = dbl[dit()];
    const EBISBox& ebisbox = ebisl[dit()];

    this->computeReactiveParticlesPerCell(a_ppc[dit()], a_level, dit(), box, ebisbox);
  }
}

void
ItoPlasmaStepper::computeReactiveParticlesPerCell(EBCellFAB&      a_ppc,
                                                  const int       a_level,
                                                  const DataIndex a_dit,
                                                  const Box       a_box,
                                                  const EBISBox&  a_ebisbox)
{
  CH_TIME("ItoPlasmaStepper::computeReactiveParticlesPerCell(ppc, lvl, dit, box)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeReactiveParticlesPerCell(ppc, lvl, dit, box)" << endl;
  }

  BaseFab<Real>& numFab = a_ppc.getSingleValuedFAB();

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>& solver = solver_it();
    const int                 idx    = solver_it.index();

    const ParticleContainer<ItoParticle>& particles     = solver->getParticles(ItoSolver::WhichContainer::Bulk);
    const BinFab<ItoParticle>&            cellParticles = particles.getCellParticles(a_level, a_dit);

    // Regular cells.
    for (BoxIterator bit(a_box); bit.ok(); ++bit) {
      const IntVect iv = bit();

      Real num = 0.0;

      if (a_ebisbox.isRegular(iv)) {

        const List<ItoParticle>& listParticles = cellParticles(iv, 0);
        for (ListIterator<ItoParticle> lit(listParticles); lit.ok(); ++lit) {
          num += lit().weight();
        }
      }

      numFab(iv, idx) = num;
    }

    // Irregular cells.
    VoFIterator& vofit = (*m_amr->getVofIterator(m_particleRealm, m_phase)[a_level])[a_dit];
    for (vofit.reset(); vofit.ok(); ++vofit) {
      const VolIndex& vof        = vofit();
      const IntVect   iv         = vof.gridIndex();
      const RealVect  normal     = a_ebisbox.normal(vof);
      const RealVect  ebCentroid = a_ebisbox.bndryCentroid(vof);

      Real num = 0.0;

      const List<ItoParticle>& listParticles = cellParticles(iv, 0);

      for (ListIterator<ItoParticle> lit(listParticles); lit.ok(); ++lit) {
        const RealVect& pos = lit().position();

        if (PolyGeom::dot((pos - ebCentroid), normal) >= 0.0) {
          num += lit().weight();
        }
      }

      numFab(iv, idx) = num;
    }
  }
}

void
ItoPlasmaStepper::computeReactiveMeanEnergiesPerCell(EBAMRCellData& a_mean_energies)
{
  CH_TIME("ItoPlasmaStepper::compute_mean-energies_per_cell(EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeReactiveParticlesPerCell(EBAMRCellData)" << endl;
  }

  DataOps::setValue(a_mean_energies, 0.0);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    this->computeReactiveMeanEnergiesPerCell(*a_mean_energies[lvl], lvl);
  }
}

void
ItoPlasmaStepper::computeReactiveMeanEnergiesPerCell(LevelData<EBCellFAB>& a_mean_energies, const int a_level)
{
  CH_TIME("ItoPlasmaStepper::computeReactiveMeanEnergiesPerCell(ppc, lvl)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeReactiveMeanEnergiesPerCell(ppc, lvl)" << endl;
  }

  const DisjointBoxLayout& dbl   = m_amr->getGrids(m_particleRealm)[a_level];
  const EBISLayout&        ebisl = m_amr->getEBISLayout(m_particleRealm, m_phase)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

    const Box      box     = dbl[dit()];
    const EBISBox& ebisbox = ebisl[dit()];

    this->computeReactiveMeanEnergiesPerCell(a_mean_energies[dit()], a_level, dit(), box, ebisbox);
  }
}

void
ItoPlasmaStepper::computeReactiveMeanEnergiesPerCell(EBCellFAB&      a_mean_energies,
                                                     const int       a_level,
                                                     const DataIndex a_dit,
                                                     const Box       a_box,
                                                     const EBISBox&  a_ebisbox)
{
  CH_TIME("ItoPlasmaStepper::computeReactiveMeanEnergiesPerCell(ppc, lvl, dit, box)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeReactiveMeanEnergiesPerCell(ppc, lvl, dit, box)" << endl;
  }

  BaseFab<Real>& numFab = a_mean_energies.getSingleValuedFAB();

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>& solver = solver_it();
    const int                 idx    = solver_it.index();

    const ParticleContainer<ItoParticle>& particles     = solver->getParticles(ItoSolver::WhichContainer::Bulk);
    const BinFab<ItoParticle>&            cellParticles = particles.getCellParticles(a_level, a_dit);

    // Regular cells.
    for (BoxIterator bit(a_box); bit.ok(); ++bit) {
      const IntVect iv = bit();

      if (a_ebisbox.isRegular(iv)) {
        Real m = 0.0;
        Real E = 0.0;

        const List<ItoParticle>& listParticles = cellParticles(iv, 0);

        for (ListIterator<ItoParticle> lit(listParticles); lit.ok(); ++lit) {
          m += lit().weight();
          E += lit().weight() * lit().energy();
        }

        numFab(iv, idx) = E / m;
      }
    }

    // Irregular cells.
    VoFIterator& vofit = (*m_amr->getVofIterator(m_particleRealm, m_phase)[a_level])[a_dit];
    for (vofit.reset(); vofit.ok(); ++vofit) {
      const VolIndex& vof        = vofit();
      const IntVect   iv         = vof.gridIndex();
      const RealVect  normal     = a_ebisbox.normal(vof);
      const RealVect  ebCentroid = a_ebisbox.bndryCentroid(vof);

      Real m = 0.0;
      Real E = 0.0;

      const List<ItoParticle>& listParticles = cellParticles(iv, 0);

      for (ListIterator<ItoParticle> lit(listParticles); lit.ok(); ++lit) {
        const RealVect& pos = lit().position();

        if (PolyGeom::dot((pos - ebCentroid), normal) >= 0.0) {
          m += lit().weight();
          E += lit().weight() * lit().energy();
        }
      }

      numFab(iv, idx) = E / m;
    }
  }
}

void
ItoPlasmaStepper::advanceReactionNetworkNWO(const Real a_dt)
{
  CH_TIME("ItoPlasmaStepper::advanceReactionNetworkNWO(dt)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::advanceReactionNetworkNWO(dt)" << endl;
  }

  this->advanceReactionNetworkNWO(m_fluid_E, m_EdotJ, a_dt);
}

void
ItoPlasmaStepper::advanceReactionNetworkNWO(const EBAMRCellData& a_E, const EBAMRCellData& a_EdotJ, const Real a_dt)
{
  CH_TIME("ItoPlasmaStepper::advanceReactionNetwork(ppc, ypc, E, sources, dt)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::advanceReactionNetwork(ppc, ypc, E, sources, dt)" << endl;
  }

  // 1. Compute the number of particles per cell. Set the number of Photons to be generated per cell to zero.
  this->computeReactiveParticlesPerCell(m_particle_ppc);
  this->computeReactiveMeanEnergiesPerCell(m_particle_eps);

  m_fluid_ppc.copy(m_particle_ppc);
  m_fluid_eps.copy(m_particle_eps);

  DataOps::setValue(m_fluid_ypc, 0.0);
  DataOps::setValue(m_particle_ypc, 0.0);
  DataOps::copy(m_particle_old, m_particle_ppc);

  // 2. Solve for the new number of particles per cell. This also obtains the number of Photons to be generated in each cell.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    this->advanceReactionNetworkNWO(*m_fluid_ppc[lvl],
                                    *m_fluid_ypc[lvl],
                                    *m_fluid_eps[lvl],
                                    *a_E[lvl],
                                    *a_EdotJ[lvl],
                                    lvl,
                                    a_dt);
  }

  // 3. Copy the results to the particle Realm.
  m_particle_ppc.copy(m_fluid_ppc);
  m_particle_ypc.copy(m_fluid_ypc);
  m_particle_eps.copy(m_fluid_eps);

  // 4. Reconcile particles on the particle Realm. Not implemented (yet).
  this->reconcileParticles(m_particle_ppc, m_particle_old, m_particle_eps, m_particle_ypc);
}

void
ItoPlasmaStepper::advanceReactionNetworkNWO(LevelData<EBCellFAB>&       a_particlesPerCell,
                                            LevelData<EBCellFAB>&       a_newPhotonsPerCell,
                                            LevelData<EBCellFAB>&       a_meanParticleEnergies,
                                            const LevelData<EBCellFAB>& a_E,
                                            const LevelData<EBCellFAB>& a_EdotJ,
                                            const int                   a_level,
                                            const Real                  a_dt)
{
  CH_TIME("ItoPlasmaStepper::advanceReactionNetwork(ppc, ypc, energies, E, sources, level, dt)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::advanceReactionNetwork(ppc, ypc, energies, E, sources, level, dt)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_fluidRealm)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
    this->advanceReactionNetworkNWO(a_particlesPerCell[dit()],
                                    a_newPhotonsPerCell[dit()],
                                    a_meanParticleEnergies[dit()],
                                    a_E[dit()],
                                    a_EdotJ[dit()],
                                    a_level,
                                    dit(),
                                    dbl[dit()],
                                    m_amr->getDx()[a_level],
                                    a_dt);
  }
}

void
ItoPlasmaStepper::advanceReactionNetworkNWO(EBCellFAB&       a_particlesPerCell,
                                            EBCellFAB&       a_newPhotonsPerCell,
                                            EBCellFAB&       a_meanParticleEnergies,
                                            const EBCellFAB& a_E,
                                            const EBCellFAB& a_EdotJ,
                                            const int        a_level,
                                            const DataIndex  a_dit,
                                            const Box        a_box,
                                            const Real       a_dx,
                                            const Real       a_dt)
{
  CH_TIME("ItoPlasmaStepper::advanceReactionNetworkNWO(ppc, ypc, E, sources, level, dit, box, dx, dt)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::advanceReactionNetworkNWO(ppc, ypc, E, sources, level, dit, box, dx, dt)" << endl;
  }

  const int num_ItoSpecies = m_physics->getNumItoSpecies();
  const int num_rtSpecies  = m_physics->getNumRtSpecies();

  const RealVect prob_lo = m_amr->getProbLo();

  const EBISBox& ebisbox = m_amr->getEBISLayout(m_fluidRealm, m_phase)[a_level][a_dit];
  const EBISBox& ebgraph = m_amr->getEBISLayout(m_fluidRealm, m_phase)[a_level][a_dit];

  const BaseFab<Real>& Efab = a_E.getSingleValuedFAB();

  Vector<long long> particles(num_ItoSpecies);
  Vector<long long> newPhotons(num_rtSpecies);
  Vector<Real>      meanEnergies(num_ItoSpecies);
  Vector<Real>      energySources(num_ItoSpecies);

  // Regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit) {
    const IntVect iv = bit();

    if (ebisbox.isRegular(iv)) {
      const Real kappa = 1.0;
      const Real dV    = kappa * pow(a_dx, SpaceDim);

      const RealVect pos = prob_lo + a_dx * (RealVect(iv) + 0.5 * RealVect::Unit);
      const RealVect E   = RealVect(D_DECL(Efab(iv, 0), Efab(iv, 1), Efab(iv, 2)));

      // Initialize for this cell.
      for (int i = 0; i < num_ItoSpecies; i++) {
        particles[i]     = llround(a_particlesPerCell.getSingleValuedFAB()(iv, i));
        meanEnergies[i]  = a_meanParticleEnergies.getSingleValuedFAB()(iv, i);
        energySources[i] = a_EdotJ.getSingleValuedFAB()(iv, i) * dV / Units::Qe;
      }

      for (int i = 0; i < num_rtSpecies; i++) {
        newPhotons[i] = 0LL;
      }

      // Do the physics advance
      m_physics->advanceParticles(particles, newPhotons, meanEnergies, energySources, a_dt, E, a_dx, kappa);

      // Set result
      for (int i = 0; i < num_ItoSpecies; i++) {
        a_particlesPerCell.getSingleValuedFAB()(iv, i)     = 1.0 * particles[i];
        a_meanParticleEnergies.getSingleValuedFAB()(iv, i) = 1.0 * meanEnergies[i];
      }

      for (int i = 0; i < num_rtSpecies; i++) {
        a_newPhotonsPerCell.getSingleValuedFAB()(iv, i) = 1.0 * newPhotons[i];
      }
    }
  }

  // Irregular cells
  VoFIterator& vofit = (*m_amr->getVofIterator(m_fluidRealm, m_phase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit) {
    const VolIndex& vof   = vofit();
    const Real      kappa = ebisbox.volFrac(vof);
    const Real      dV    = kappa * pow(a_dx, SpaceDim);
    const RealVect  pos   = EBArith::getVofLocation(vof, a_dx * RealVect::Unit, prob_lo);
    const RealVect  E     = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));

    // Initialize for this cell.
    for (int i = 0; i < num_ItoSpecies; i++) {
      particles[i]     = a_particlesPerCell(vof, i);
      meanEnergies[i]  = a_meanParticleEnergies(vof, i);
      energySources[i] = a_EdotJ(vof, i) * dV / Units::Qe;
    }

    for (int i = 0; i < num_rtSpecies; i++) {
      newPhotons[i] = 0LL;
    }

    m_physics->advanceParticles(particles, newPhotons, meanEnergies, energySources, a_dt, E, a_dx, kappa);

    // Set result
    for (int i = 0; i < num_ItoSpecies; i++) {
      a_particlesPerCell(vof, i)     = 1.0 * particles[i];
      a_meanParticleEnergies(vof, i) = 1.0 * meanEnergies[i];
    }

    for (int i = 0; i < num_rtSpecies; i++) {
      a_newPhotonsPerCell(vof, i) = 1.0 * newPhotons[i];
    }
  }
}

void
ItoPlasmaStepper::reconcileParticles(const EBAMRCellData& a_newParticlesPerCell,
                                     const EBAMRCellData& a_oldParticlesPerCell,
                                     const EBAMRCellData& a_meanParticleEnergies,
                                     const EBAMRCellData& a_newPhotonsPerCell)
{
  CH_TIME("ItoPlasmaStepper::reconcileParticles(EBAMRCellData, EBAMRCellData, EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::reconcileParticles(EBAMRCellData, EBAMRCellData, EBAMRCellData)";
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {

    this->reconcileParticles(*a_newParticlesPerCell[lvl],
                             *a_oldParticlesPerCell[lvl],
                             *a_meanParticleEnergies[lvl],
                             *a_newPhotonsPerCell[lvl],
                             lvl);
  }
}

void
ItoPlasmaStepper::reconcileParticles(const LevelData<EBCellFAB>& a_newParticlesPerCell,
                                     const LevelData<EBCellFAB>& a_oldParticlesPerCell,
                                     const LevelData<EBCellFAB>& a_meanParticleEnergies,
                                     const LevelData<EBCellFAB>& a_newPhotonsPerCell,
                                     const int                   a_level)
{
  CH_TIME("ItoPlasmaStepper::reconcileParticles(LevelData<EBCellFAB>x3, int)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::reconcileParticles(LevelData<EBCellFAB>x3, int)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
    this->reconcileParticles(a_newParticlesPerCell[dit()],
                             a_oldParticlesPerCell[dit()],
                             a_meanParticleEnergies[dit()],
                             a_newPhotonsPerCell[dit()],
                             a_level,
                             dit(),
                             dbl[dit()],
                             m_amr->getDx()[a_level]);
  }
}

void
ItoPlasmaStepper::reconcileParticles(const EBCellFAB& a_newParticlesPerCell,
                                     const EBCellFAB& a_oldParticlesPerCell,
                                     const EBCellFAB& a_meanParticleEnergies,
                                     const EBCellFAB& a_newPhotonsPerCell,
                                     const int        a_level,
                                     const DataIndex  a_dit,
                                     const Box        a_box,
                                     const Real       a_dx)
{
  CH_TIME("ItoPlasmaStepper::reconcileParticles(EBCellFABx3, int, DataIndex, Box, Real)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::reconcileParticles(EBCellFABx3, int, DataIndex, Box, Real)" << endl;
  }

  const int num_ItoSpecies = m_physics->getNumItoSpecies();
  const int num_rtSpecies  = m_physics->getNumRtSpecies();

  const RealVect prob_lo = m_amr->getProbLo();

  const EBISBox& ebisbox = m_amr->getEBISLayout(m_particleRealm, m_phase)[a_level][a_dit];
  const EBISBox& ebgraph = m_amr->getEBISLayout(m_particleRealm, m_phase)[a_level][a_dit];

  Vector<BinFab<ItoParticle>*> particlesFAB(num_ItoSpecies);
  Vector<BinFab<Photon>*>      sourcePhotonsFAB(num_rtSpecies);
  Vector<BinFab<Photon>*>      bulkPhotonsFAB(num_rtSpecies);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>& solver = solver_it();
    const int                 idx    = solver_it.index();

    ParticleContainer<ItoParticle>& solverParticles = solver->getParticles(ItoSolver::WhichContainer::Bulk);

    particlesFAB[idx] = &(solverParticles.getCellParticles(a_level, a_dit));
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<McPhoto>& solver = solver_it();
    const int               idx    = solver_it.index();

    ParticleContainer<Photon>& solverBulkPhotons = solver->getBulkPhotons();
    ParticleContainer<Photon>& solverSourPhotons = solver->getSourcePhotons();

    bulkPhotonsFAB[idx]   = &(solverBulkPhotons.getCellParticles(a_level, a_dit));
    sourcePhotonsFAB[idx] = &(solverSourPhotons.getCellParticles(a_level, a_dit));
  }

  // Regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit) {
    const IntVect iv = bit();
    if (ebisbox.isRegular(iv)) {
      const RealVect cellPos       = prob_lo + a_dx * (RealVect(iv) + 0.5 * RealVect::Unit);
      const RealVect centroidPos   = cellPos;
      const RealVect lo            = -0.5 * RealVect::Unit;
      const RealVect hi            = 0.5 * RealVect::Unit;
      const RealVect bndryCentroid = RealVect::Zero;
      const RealVect bndryNormal   = RealVect::Zero;
      const Real     kappa         = 1.0;

      Vector<List<ItoParticle>*>       particles(num_ItoSpecies);
      Vector<List<Photon>*>            bulkPhotons(num_rtSpecies);
      Vector<List<Photon>*>            sourcePhotons(num_rtSpecies);
      Vector<RefCountedPtr<RtSpecies>> photoSpecies(num_rtSpecies);

      Vector<Real>      particleMeanEnergies(num_ItoSpecies);
      Vector<long long> numNewParticles(num_ItoSpecies);
      Vector<long long> numOldParticles(num_ItoSpecies);
      Vector<long long> numNewPhotons(num_rtSpecies);

      for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
        const int idx = solver_it.index();

        particles[idx]            = &((*particlesFAB[idx])(iv, 0));
        particleMeanEnergies[idx] = a_meanParticleEnergies.getSingleValuedFAB()(iv, idx);
        numNewParticles[idx]      = llround(a_newParticlesPerCell.getSingleValuedFAB()(iv, idx));
        numOldParticles[idx]      = llround(a_oldParticlesPerCell.getSingleValuedFAB()(iv, idx));
      }

      for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it) {
        const int idx = solver_it.index();

        bulkPhotons[idx]   = &((*bulkPhotonsFAB[idx])(iv, 0));
        sourcePhotons[idx] = &((*sourcePhotonsFAB[idx])(iv, 0));
        photoSpecies[idx]  = solver_it()->getSpecies();
        numNewPhotons[idx] = llround(a_newPhotonsPerCell.getSingleValuedFAB()(iv, idx));

        sourcePhotons[idx]->clear();
      }

      // Reconcile particles, Photons, and photoionization
      m_physics->reconcileParticles(particles,
                                    numNewParticles,
                                    numOldParticles,
                                    cellPos,
                                    centroidPos,
                                    lo,
                                    hi,
                                    bndryCentroid,
                                    bndryNormal,
                                    a_dx,
                                    kappa);
      m_physics->reconcilePhotons(sourcePhotons,
                                  numNewPhotons,
                                  cellPos,
                                  centroidPos,
                                  lo,
                                  hi,
                                  bndryCentroid,
                                  bndryNormal,
                                  a_dx,
                                  kappa);
      m_physics->reconcilePhotoionization(particles, particleMeanEnergies, numNewParticles, bulkPhotons);
      m_physics->setMeanParticleEnergy(particles, particleMeanEnergies);

      // Clear the bulk Photons - they have now been absorbed on the mesh.
      for (int i = 0; i < num_rtSpecies; i++) {
        //	bulkPhotons[i]->clear();
      }
    }
  }

  // Irregular cells
  VoFIterator& vofit = (*m_amr->getVofIterator(m_particleRealm, m_phase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit) {
    const VolIndex& vof           = vofit();
    const IntVect   iv            = vof.gridIndex();
    const Real      kappa         = ebisbox.volFrac(vof);
    const RealVect  cellPos       = EBArith::getVofLocation(vof, a_dx * RealVect::Unit, prob_lo);
    const RealVect  centroidPos   = ebisbox.centroid(vof);
    const RealVect  bndryCentroid = ebisbox.bndryCentroid(vof);
    const RealVect  bndryNormal   = ebisbox.normal(vof);

    RealVect lo = -0.5 * RealVect::Unit;
    RealVect hi = 0.5 * RealVect::Unit;
    if (kappa < 1.0) {
      DataOps::computeMinValidBox(lo, hi, bndryNormal, bndryCentroid);
    }

    Vector<List<ItoParticle>*>       particles(num_ItoSpecies);
    Vector<List<Photon>*>            bulkPhotons(num_rtSpecies);
    Vector<List<Photon>*>            sourcePhotons(num_rtSpecies);
    Vector<RefCountedPtr<RtSpecies>> photoSpecies(num_rtSpecies);

    Vector<Real>      particleMeanEnergies(num_ItoSpecies);
    Vector<long long> numNewParticles(num_ItoSpecies);
    Vector<long long> numOldParticles(num_ItoSpecies);
    Vector<long long> numNewPhotons(num_rtSpecies);

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
      const int idx = solver_it.index();

      particles[idx]            = &((*particlesFAB[idx])(iv, 0));
      particleMeanEnergies[idx] = a_meanParticleEnergies(vof, idx);
      numNewParticles[idx]      = llround(a_newParticlesPerCell(vof, idx));
      numOldParticles[idx]      = llround(a_oldParticlesPerCell(vof, idx));
    }

    for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it) {
      const int idx = solver_it.index();

      bulkPhotons[idx]   = &((*bulkPhotonsFAB[idx])(iv, 0));
      sourcePhotons[idx] = &((*sourcePhotonsFAB[idx])(iv, 0));
      photoSpecies[idx]  = solver_it()->getSpecies();
      numNewPhotons[idx] = llround(a_newPhotonsPerCell(vof, idx));

      sourcePhotons[idx]->clear();
    }

    // Reconcile particles, Photons, and photoionization
    m_physics->reconcileParticles(particles,
                                  numNewParticles,
                                  numOldParticles,
                                  cellPos,
                                  centroidPos,
                                  lo,
                                  hi,
                                  bndryCentroid,
                                  bndryNormal,
                                  a_dx,
                                  kappa);
    m_physics->reconcilePhotons(sourcePhotons,
                                numNewPhotons,
                                cellPos,
                                centroidPos,
                                lo,
                                hi,
                                bndryCentroid,
                                bndryNormal,
                                a_dx,
                                kappa);
    m_physics->reconcilePhotoionization(particles, particleMeanEnergies, numNewParticles, bulkPhotons);
    m_physics->setMeanParticleEnergy(particles, particleMeanEnergies);

    // Clear the bulk Photons - they have now been absorbed on the mesh.
    for (int i = 0; i < num_rtSpecies; i++) {
      //      bulkPhotons[i]->clear();
    }
  }
}

void
ItoPlasmaStepper::advanceReactionNetwork(const Real a_dt)
{
  CH_TIME("ItoPlasmaStepper::advanceReactionNetwork(a_dt)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::advanceReactionNetwork(a_dt)" << endl;
  }

  if (m_useNewReactionAlgorithm) {
    this->advanceReactionNetworkNWO(a_dt);
  }
  else {
    const int num_ItoSpecies = m_physics->getNumItoSpecies();
    const int num_rtSpecies  = m_physics->getNumRtSpecies();

    Vector<ParticleContainer<ItoParticle>*> particles(num_ItoSpecies);   // Current particles.
    Vector<ParticleContainer<Photon>*>      bulk_Photons(num_rtSpecies); // Photons absorbed on mesh
    Vector<ParticleContainer<Photon>*>      new_Photons(num_rtSpecies);  // Produced Photons go here.

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
      particles[solver_it.index()] = &(solver_it()->getParticles(ItoSolver::WhichContainer::Bulk));
    }

    for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it) {
      bulk_Photons[solver_it.index()] = &(solver_it()->getBulkPhotons());
      new_Photons[solver_it.index()]  = &(solver_it()->getSourcePhotons());
    }

    this->advanceReactionNetwork(particles, bulk_Photons, new_Photons, m_energy_sources, m_particle_E, a_dt);
  }
}

void
ItoPlasmaStepper::advanceReactionNetwork(Vector<ParticleContainer<ItoParticle>*>& a_particles,
                                         Vector<ParticleContainer<Photon>*>&      a_Photons,
                                         Vector<ParticleContainer<Photon>*>&      a_newPhotons,
                                         const Vector<EBAMRCellData>&             a_sources,
                                         const EBAMRCellData&                     a_E,
                                         const Real                               a_dt)
{
  CH_TIME(
    "ItoPlasmaStepper::advanceReactionNetwork(Vector<ParticleContainer*> x3, Vector<EBAMRCellData*>, EBAMRCellData, Real)");
  if (m_verbosity > 5) {
    pout()
      << "ItoPlasmaStepper::advanceReactionNetwork(Vector<ParticleContainer*> x3, Vector<EBAMRCellData*>, EBAMRCellData, Real)"
      << endl;
  }

  const int num_ItoSpecies = m_physics->getNumItoSpecies();
  const int num_rtSpecies  = m_physics->getNumRtSpecies();

  Vector<AMRCellParticles<ItoParticle>*> particles(num_ItoSpecies);
  Vector<AMRCellParticles<Photon>*>      Photons(num_ItoSpecies);
  Vector<AMRCellParticles<Photon>*>      newPhotons(num_ItoSpecies);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    const int idx  = solver_it.index();
    particles[idx] = &(a_particles[idx]->getCellParticles());
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it) {
    const int idx   = solver_it.index();
    Photons[idx]    = &(a_Photons[idx]->getCellParticles());
    newPhotons[idx] = &(a_newPhotons[idx]->getCellParticles());
  }

  //Advance reaction network
  this->advanceReactionNetwork(particles, Photons, newPhotons, a_sources, m_particle_E, a_dt);
}

void
ItoPlasmaStepper::advanceReactionNetwork(Vector<AMRCellParticles<ItoParticle>*>& a_particles,
                                         Vector<AMRCellParticles<Photon>*>&      a_Photons,
                                         Vector<AMRCellParticles<Photon>*>&      a_newPhotons,
                                         const Vector<EBAMRCellData>&            a_sources,
                                         const EBAMRCellData&                    a_E,
                                         const Real                              a_dt)
{
  CH_TIME(
    "ItoPlasmaStepper::advanceReactionNetwork(Vector<AMRCellParticles*> x 3, Vector<EBAMRCellData*>, EBAMRCellData, Real)");
  if (m_verbosity > 5) {
    pout()
      << "ItoPlasmaStepper::advanceReactionNetwork(Vector<AMRCellParticles*> x 3, Vector<EBAMRCellData*>, EBAMRCellData, Real)"
      << endl;
  }

  const int num_ItoSpecies = m_physics->getNumItoSpecies();
  const int num_rtSpecies  = m_physics->getNumRtSpecies();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    Vector<LayoutData<BinFab<ItoParticle>>*> particles(num_ItoSpecies);
    Vector<LayoutData<BinFab<Photon>>*>      Photons(num_rtSpecies);
    Vector<LayoutData<BinFab<Photon>>*>      newPhotons(num_rtSpecies);
    Vector<LevelData<EBCellFAB>*>            sources(num_ItoSpecies);

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
      const int idx  = solver_it.index();
      particles[idx] = &(*(*a_particles[idx])[lvl]);
      sources[idx]   = &(*(a_sources[idx])[lvl]);
    }

    for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it) {
      const int idx   = solver_it.index();
      Photons[idx]    = &(*(*a_Photons[idx])[lvl]);
      newPhotons[idx] = &(*(*a_newPhotons[idx])[lvl]);
    }

    this->advanceReactionNetwork(particles, Photons, newPhotons, sources, *a_E[lvl], lvl, a_dt);
  }
}

void
ItoPlasmaStepper::advanceReactionNetwork(Vector<LayoutData<BinFab<ItoParticle>>*>& a_particles,
                                         Vector<LayoutData<BinFab<Photon>>*>&      a_Photons,
                                         Vector<LayoutData<BinFab<Photon>>*>&      a_newPhotons,
                                         const Vector<LevelData<EBCellFAB>*>&      a_sources,
                                         const LevelData<EBCellFAB>&               a_E,
                                         const int                                 a_lvl,
                                         const Real                                a_dt)
{
  CH_TIME(
    "ItoPlasmaStepper::advanceReactionNetwork(Vector<LD<BinFab>* > x 3, Vector<LD<EBCellFAB>*>, EBAMRCellData, level, dt)");
  if (m_verbosity > 5) {
    pout()
      << "ItoPlasmaStepper::advanceReactionNetwork(Vector<LD<BinFab>* > x 3, Vector<LD<EBCellFAB>*>, EBAMRCellData, level, dt)"
      << endl;
  }

  const int num_ItoSpecies = m_physics->getNumItoSpecies();
  const int num_rtSpecies  = m_physics->getNumRtSpecies();

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[a_lvl];
  const Real               dx  = m_amr->getDx()[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
    const Box box = dbl.get(dit());

    Vector<BinFab<ItoParticle>*> particles(num_ItoSpecies);
    Vector<BinFab<Photon>*>      Photons(num_rtSpecies);
    ;
    Vector<BinFab<Photon>*> newPhotons(num_rtSpecies);
    Vector<EBCellFAB*>      sources(num_ItoSpecies);

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
      const int idx  = solver_it.index();
      particles[idx] = &((*a_particles[idx])[dit()]);
      sources[idx]   = &((*a_sources[idx])[dit()]);
    }

    for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it) {
      const int idx   = solver_it.index();
      Photons[idx]    = &((*a_Photons[idx])[dit()]);
      newPhotons[idx] = &((*a_newPhotons[idx])[dit()]);
    }

    this->advanceReactionNetwork(particles, Photons, newPhotons, sources, a_E[dit()], a_lvl, dit(), box, dx, a_dt);
  }
}

void
ItoPlasmaStepper::advanceReactionNetwork(Vector<BinFab<ItoParticle>*>& a_particles,
                                         Vector<BinFab<Photon>*>&      a_Photons,
                                         Vector<BinFab<Photon>*>&      a_newPhotons,
                                         const Vector<EBCellFAB*>&     a_sources,
                                         const EBCellFAB&              a_E,
                                         const int                     a_lvl,
                                         const DataIndex               a_dit,
                                         const Box                     a_box,
                                         const Real                    a_dx,
                                         const Real                    a_dt)
{
  CH_TIME(
    "ItoPlasmaStepper::advanceReactionNetwork(Vector<BinFab*> x 3, Vector<EBCellFAB*>, EBCellFAB, level, dit, box, dx, dt)");
  if (m_verbosity > 5) {
    pout()
      << "ItoPlasmaStepper::advanceReactionNetwork(Vector<BinFab*> x 3, Vector<EBCellFAB*>, EBCellFAB, level, dit, box, dx, dt)"
      << endl;
  }

  const int comp = 0;

  const int num_ItoSpecies = m_physics->getNumItoSpecies();
  const int num_rtSpecies  = m_physics->getNumRtSpecies();

  const RealVect prob_lo = m_amr->getProbLo();
  const RealVect dx      = a_dx * RealVect::Unit;

  const EBISBox& ebisbox = m_amr->getEBISLayout(m_particleRealm, m_phase)[a_lvl][a_dit];
  const EBISBox& ebgraph = m_amr->getEBISLayout(m_particleRealm, m_phase)[a_lvl][a_dit];

  const BaseFab<Real>& Efab = a_E.getSingleValuedFAB();

  // Regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit) {
    const IntVect iv = bit();

    if (ebisbox.isRegular(iv)) {
      const Real     kappa = 1.0;
      const RealVect pos   = prob_lo + a_dx * (RealVect(iv) + 0.5 * RealVect::Unit);
      const RealVect e     = RealVect(D_DECL(Efab(iv, 0), Efab(iv, 1), Efab(iv, 2)));

      Vector<List<ItoParticle>*> particles(num_ItoSpecies);
      Vector<List<Photon>*>      Photons(num_rtSpecies);
      Vector<List<Photon>*>      newPhotons(num_rtSpecies);
      Vector<Real>               sources(num_ItoSpecies);

      for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
        const int idx = solver_it.index();

        List<ItoParticle>& bp = (*a_particles[idx])(iv, comp);
        particles[idx]        = &bp;

        const BaseFab<Real>& sourcesFAB = a_sources[idx]->getSingleValuedFAB();
        sources[idx]                    = sourcesFAB(iv, comp);
      }

      for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it) {
        const int idx = solver_it.index();

        List<Photon>& bp    = (*a_Photons[idx])(iv, comp);
        List<Photon>& bpNew = (*a_newPhotons[idx])(iv, comp);

        Photons[idx]    = &bp;
        newPhotons[idx] = &bpNew;
      }

      // Dummy stuff for regular cells
      const RealVect lo = -0.5 * RealVect::Unit;
      const RealVect hi = 0.5 * RealVect::Unit;
      const RealVect n  = RealVect::Zero;
      const RealVect c  = RealVect::Zero;

      // Advance reactions
      m_physics
        ->advanceReactionNetwork(particles, Photons, newPhotons, sources, e, pos, c, c, n, lo, hi, a_dx, kappa, a_dt);
    }
  }

  // Now do the irregular cells
  VoFIterator& vofit = (*m_amr->getVofIterator(m_particleRealm, m_phase)[a_lvl])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit) {
    const VolIndex vof   = vofit();
    const IntVect  iv    = vof.gridIndex();
    const RealVect pos   = prob_lo + a_dx * (RealVect(iv) + 0.5 * RealVect::Unit);
    const RealVect cen   = ebisbox.centroid(vof);
    const Real     kappa = ebisbox.volFrac(vof);
    const RealVect e     = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));
    const RealVect n     = ebisbox.normal(vof);
    const RealVect ebc   = ebisbox.bndryCentroid(vof);

    // Compute a small box that encloses the cut-cell volume
    RealVect lo = -0.5 * RealVect::Unit;
    RealVect hi = 0.5 * RealVect::Unit;
    if (kappa < 1.0) {
      DataOps::computeMinValidBox(lo, hi, n, ebc);
    }

    Vector<List<ItoParticle>*> particles(num_ItoSpecies);
    Vector<List<Photon>*>      Photons(num_rtSpecies);
    Vector<List<Photon>*>      newPhotons(num_rtSpecies);
    Vector<Real>               sources(num_ItoSpecies);

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
      const int idx = solver_it.index();

      List<ItoParticle>& bp = (*a_particles[idx])(iv, comp);
      particles[idx]        = &bp;

      const BaseFab<Real>& sourcesFAB = a_sources[idx]->getSingleValuedFAB();
      sources[idx]                    = sourcesFAB(iv, comp);
    }

    for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it) {
      const int idx = solver_it.index();

      List<Photon>& bp    = (*a_Photons[idx])(iv, comp);
      List<Photon>& bpNew = (*a_newPhotons[idx])(iv, comp);

      Photons[idx]    = &bp;
      newPhotons[idx] = &bpNew;
    }

    // Advance reactions
    m_physics
      ->advanceReactionNetwork(particles, Photons, newPhotons, sources, e, pos, cen, ebc, n, lo, hi, a_dx, kappa, a_dt);
  }
}

Real
ItoPlasmaStepper::computePhysicsDt() const
{
  CH_TIME("ItoPlasmaStepper::computePhysicsDt()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computePhysicsDt()" << endl;
  }

  // TLDR: This is done on the particle Realm because of the densities (which are defined on the particle Realm).

  const Real dt = this->computePhysicsDt(m_particle_E, m_ito->getDensities());

  return dt;
}

Real
ItoPlasmaStepper::computePhysicsDt(const EBAMRCellData& a_E, const Vector<EBAMRCellData*> a_densities) const
{
  CH_TIME("ItoPlasmaStepper::computePhysicsDt(EBAMRCellFAB, Vector<EBAMRCellFAB*>)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computePhysicsDt(EBAMRCellFAB, Vector<EBAMRCellFAB*>)" << endl;
  }

  const int num_ItoSpecies = m_physics->getNumItoSpecies();

  Real minDt = 1.E99;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {

    Vector<LevelData<EBCellFAB>*> densities(num_ItoSpecies);

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
      const int idx = solver_it.index();

      densities[idx] = &(*(*a_densities[idx])[lvl]);
    }

    const Real levelDt = this->computePhysicsDt(*a_E[lvl], densities, lvl);

    minDt = Min(minDt, levelDt);
  }

  return minDt;
}

Real
ItoPlasmaStepper::computePhysicsDt(const LevelData<EBCellFAB>&         a_E,
                                   const Vector<LevelData<EBCellFAB>*> a_densities,
                                   const int                           a_level) const
{
  CH_TIME("ItoPlasmaStepper::computePhysicsDt(LevelData<EBCellFAB>, Vector<LevelData<EBCellFAB> *>, int)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computePhysicsDt(LevelData<EBCellFAB>, Vector<LevelData<EBCellFAB> *>, int)" << endl;
  }

  const int num_ItoSpecies = m_physics->getNumItoSpecies();

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[a_level];

  Real minDt = 1.E99;

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

    Vector<EBCellFAB*> densities(num_ItoSpecies);

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
      const int idx = solver_it.index();

      densities[idx] = &((*a_densities[idx])[dit()]);
    }

    const Real boxDt = this->computePhysicsDt(a_E[dit()], densities, a_level, dit(), dbl.get(dit()));

    minDt = Min(minDt, boxDt);
  }

  // MPI reduction....
#ifdef CH_MPI
  Real tmp    = 1.;
  int  result = MPI_Allreduce(&minDt, &tmp, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if (result != MPI_SUCCESS) {
    MayDay::Error("ItoPlasmaStepper::computePhysicsDt(lvl) - communication error on norm");
  }
  minDt = tmp;
#endif

  return minDt;
}

Real
ItoPlasmaStepper::computePhysicsDt(const EBCellFAB&         a_E,
                                   const Vector<EBCellFAB*> a_densities,
                                   const int                a_level,
                                   const DataIndex          a_dit,
                                   const Box                a_box) const
{
  CH_TIME("ItoPlasmaStepper::computePhysicsDt(EBCellFAB, Vector<EBCellFAB*>, lvl, dit, box)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computePhysicsDt(EBCellFAB, Vector<EBCellFAB*>, lvl, dit, box)" << endl;
  }

  Real minDt = 1.E99;

  const int num_ItoSpecies = m_physics->getNumItoSpecies();

  const int            comp    = 0;
  const Real           dx      = m_amr->getDx()[a_level];
  const RealVect       prob_lo = m_amr->getProbLo();
  const BaseFab<Real>& E       = a_E.getSingleValuedFAB();
  const EBISBox&       ebisbox = m_amr->getEBISLayout(m_particleRealm, m_phase)[a_level][a_dit];

  // Do regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit) {
    const IntVect  iv  = bit();
    const RealVect pos = m_amr->getProbLo() + dx * (RealVect(iv) + 0.5 * RealVect::Unit);
    const RealVect e   = RealVect(D_DECL(E(iv, 0), E(iv, 1), E(iv, 2)));

    if (ebisbox.isRegular(iv)) {
      Vector<Real> densities(num_ItoSpecies);
      for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
        const int idx = solver_it.index();

        const BaseFab<Real>& basefab = a_densities[idx]->getSingleValuedFAB();
        densities[idx]               = basefab(iv, comp);
      }

      const Real cellDt = m_physics->computeDt(e, pos, densities);

      minDt = Min(minDt, cellDt);
    }
  }

  // Do irregular cells
  VoFIterator& vofit = (*m_amr->getVofIterator(m_particleRealm, m_phase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit) {
    const VolIndex& vof = vofit();
    const RealVect  e   = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));
    const RealVect  pos = EBArith::getVofLocation(vof, dx * RealVect::Unit, prob_lo);

    Vector<Real> densities(num_ItoSpecies);

    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
      const int idx  = solver_it.index();
      densities[idx] = (*a_densities[idx])(vof, comp);
    }

    const Real cellDt = m_physics->computeDt(e, pos, densities);

    minDt = Min(minDt, cellDt);
  }

  return minDt;
}

void
ItoPlasmaStepper::advancePhotons(const Real a_dt)
{
  CH_TIME("ItoPlasmaStepper::advancePhotons(a_dt)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::advance_advancePhotons(a_dt)" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<McPhoto>& solver = solver_it();

    // Add source Photons and move the Photons
    ParticleContainer<Photon>& Photons       = solver->getPhotons();
    ParticleContainer<Photon>& bulkPhotons   = solver->getBulkPhotons();
    ParticleContainer<Photon>& ebPhotons     = solver->getEbPhotons();
    ParticleContainer<Photon>& domainPhotons = solver->getDomainPhotons();
    ParticleContainer<Photon>& sourcePhotons = solver->getSourcePhotons();

    if (solver->isInstantaneous()) {
      solver->clear(Photons);

      // Add source Photons
      Photons.addParticles(sourcePhotons);
      solver->clear(sourcePhotons);

      // Instantaneous advance
      solver->advancePhotonsInstantaneous(bulkPhotons, ebPhotons, domainPhotons, Photons);
    }
    else {
      // Add source Photons
      Photons.addParticles(sourcePhotons);
      solver->clear(sourcePhotons);

      // Stationary advance
      solver->advancePhotonsTransient(bulkPhotons, ebPhotons, domainPhotons, Photons, a_dt);
    }
  }
}

void
ItoPlasmaStepper::sortPhotonsByCell(McPhoto::WhichContainer a_which)
{
  CH_TIME("ItoPlasmaStepper::sortPhotonsByCell(McPhoto::WhichContainer)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::sortPhotonsByCell(McPhoto::WhichContainer)" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->sortPhotonsByCell(a_which);
  }
}

void
ItoPlasmaStepper::sortPhotonsByPatch(McPhoto::WhichContainer a_which)
{
  CH_TIME("ItoPlasmaStepper::sortPhotonsByPatch(McPhoto::WhichContainer)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::sortPhotonsByPatch(McPhoto::WhichContainer)" << endl;
  }

  for (auto solver_it = m_rte->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->sortPhotonsByPatch(a_which);
  }
}

bool
ItoPlasmaStepper::loadBalanceThisRealm(const std::string a_realm) const
{
  CH_TIME("TimeStepper::loadBalanceThisRealm");
  if (m_verbosity > 5) {
    pout() << "TimeStepper::loadBalanceThisRealm" << endl;
  }

  bool ret = false;

  if (a_realm == m_particleRealm && m_loadBalance) {
    ret = true;
  }

  return ret;
}

void
ItoPlasmaStepper::loadBalanceBoxes(Vector<Vector<int>>&             a_procs,
                                   Vector<Vector<Box>>&             a_boxes,
                                   const std::string                a_realm,
                                   const Vector<DisjointBoxLayout>& a_grids,
                                   const int                        a_lmin,
                                   const int                        a_finestLevel)
{
  CH_TIME("ItoPlasmaStepper::loadBalanceBoxes");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper_stepper::loadBalanceBoxes" << endl;
  }

  if (m_loadBalance && a_realm == m_particleRealm) {
    this->loadBalanceParticleRealm(a_procs, a_boxes, a_realm, a_grids, a_lmin, a_finestLevel);
  }
  else {
    MayDay::Abort("ItoPlasmaStepper::loadBalanceBoxes - shouldn't happen, how did you get here..?");
  }
}

Vector<long int>
ItoPlasmaStepper::getCheckpointLoads(const std::string a_realm, const int a_level) const
{
  CH_TIME("ItoPlasmaStepper::getCheckpointLoads(...)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper_stepper::getCheckpointLoads(...)" << endl;
  }

  const DisjointBoxLayout& dbl  = m_amr->getGrids(a_realm)[a_level];
  const int                nbox = dbl.size();

  Vector<long int> loads(nbox, 0L);
  if (m_loadBalance && a_realm == m_particleRealm) {
    Vector<RefCountedPtr<ItoSolver>> lb_solvers = this->getLoadBalanceSolvers();

    for (int isolver = 0; isolver < lb_solvers.size(); isolver++) {
      Vector<long int> solver_loads(nbox, 0L);
      lb_solvers[isolver]->computeLoads(solver_loads, dbl, a_level);

      for (int ibox = 0; ibox < nbox; ibox++) {
        loads[ibox] += solver_loads[ibox];
      }
    }

    // Now add the "constant" loads
    for (LayoutIterator lit = dbl.layoutIterator(); lit.ok(); ++lit) {
      const Box box  = dbl[lit()];
      const int ibox = lit().intCode();

      loads[ibox] += lround(m_loadPerCell * box.numPts());
    }
  }
  else {
    loads = TimeStepper::getCheckpointLoads(a_realm, a_level);
  }

  return loads;
}

void
ItoPlasmaStepper::loadBalanceParticleRealm(Vector<Vector<int>>&             a_procs,
                                           Vector<Vector<Box>>&             a_boxes,
                                           const std::string                a_realm,
                                           const Vector<DisjointBoxLayout>& a_grids,
                                           const int                        a_lmin,
                                           const int                        a_finestLevel)
{
  CH_TIME("ItoPlasmaStepper::loadBalanceParticleRealm(...)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper_stepper::loadBalanceParticleRealm(...)" << endl;
  }

  // Decompose the DisjointBoxLayout
  a_procs.resize(1 + a_finestLevel);
  a_boxes.resize(1 + a_finestLevel);

  for (int lvl = a_lmin; lvl <= a_finestLevel; lvl++) {
    a_procs[lvl] = a_grids[lvl].procIDs();
    a_boxes[lvl] = a_grids[lvl].boxArray();
  }

  // Get the particles that we will use for load balancing.
  Vector<RefCountedPtr<ItoSolver>> lb_solvers = this->getLoadBalanceSolvers();

  // Regrid particles onto the "dummy grids" a_grids
  for (int i = 0; i < lb_solvers.size(); i++) {
    ParticleContainer<ItoParticle>& particles = lb_solvers[i]->getParticles(ItoSolver::WhichContainer::Bulk);

    m_amr->remapToNewGrids(particles, a_lmin, a_finestLevel);

    // If we make superparticles during regrids, do it here so we can better estimate the computational loads for each patch. This way, if a grid is removed the realistic
    // load estimate of the underlying grid(s) is improved.
    if (m_regridSuperparticles) {
      particles.sortParticlesByCell();
      lb_solvers[i]->makeSuperparticles(ItoSolver::WhichContainer::Bulk, m_particlesPerCell);
      particles.sortParticlesByPatch();
    }
  }

  // Get loads on each level
  Vector<Vector<long int>> loads(1 + a_finestLevel);
  for (int lvl = 0; lvl <= a_finestLevel; lvl++) {
    loads[lvl] = this->getCheckpointLoads(a_realm, lvl);
  }

  // Do the actual load balancing
  LoadBalancing::sort(a_boxes, loads, m_boxSort);
  LoadBalancing::balanceLevelByLevel(a_procs, loads, a_boxes);
  //  LoadBalancing::hierarchy(a_procs, loads, a_boxes); If you want to try something crazy...

  // Go back to "pre-regrid" mode so we can get particles to the correct patches after load balancing.
  for (int i = 0; i < lb_solvers.size(); i++) {
    ParticleContainer<ItoParticle>& particles = lb_solvers[i]->getParticles(ItoSolver::WhichContainer::Bulk);
    particles.preRegrid(a_lmin);
  }
}

Vector<RefCountedPtr<ItoSolver>>
ItoPlasmaStepper::getLoadBalanceSolvers() const
{
  CH_TIME("ItoPlasmaStepper::getLoadBalanceSolvers()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::getLoadBalanceSolvers()" << endl;
  }

  Vector<RefCountedPtr<ItoSolver>> lb_solvers;

  if (m_loadBalance_idx < 0) {
    for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
      RefCountedPtr<ItoSolver>& solver = solver_it();

      lb_solvers.push_back(solver);
    }
  }
  else {
    RefCountedPtr<ItoSolver>& solver = m_ito->getSolvers()[m_loadBalance_idx];
    lb_solvers.push_back(solver);
  }

  return lb_solvers;
}

void
ItoPlasmaStepper::computeEdotJSource(const Real a_dt)
{
  CH_TIME("ItoPlasmaStepper::computeEdotJSource(a_dt)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeEdotJSource(a_dt)" << endl;
  }

  // Swap between these two.
  if (m_useNewReactionAlgorithm) {
    //    this->computeEdotJSourceNWO();
    this->computeEdotJSourceNWO2(a_dt);
  }
  else {
    this->computeEdotJSource();
  }
}

void
ItoPlasmaStepper::computeEdotJSource()
{
  CH_TIME("ItoPlasmaStepper::computeEdotJSource()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeEdotJSource()" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();
    const int q   = species->getChargeNumber();

    DataOps::setValue(m_energy_sources[idx], 0.0);

    // Do mobile contribution.
    if (q != 0 && solver->isMobile()) {

      // Drift contribution
      solver->depositConductivity(m_particle_scratch1,
                                  solver->getParticles(ItoSolver::WhichContainer::Bulk)); // Deposit mu*n
      DataOps::copy(
        m_particle_scratchD,
        m_particle_E); // Could use m_particle_E or solver's m_velo_func here, but m_velo_func = +/- E (depends on q)

      DataOps::multiplyScalar(m_particle_scratchD, m_particle_scratch1);           // m_particle_scratchD = mu*n*E
      DataOps::dotProduct(m_particle_scratch1, m_particle_E, m_particle_scratchD); // m_particle_scratch1 = mu*n*E*E
      DataOps::incr(m_energy_sources[idx], m_particle_scratch1, 1.0);              // a_source[idx] += mu*n*E*E
    }

    // Diffusive contribution
    if (q != 0 && solver->isDiffusive()) {

      // Compute the negative gradient of the diffusion term
      solver->depositDiffusivity(m_particle_scratch1, solver->getParticles(ItoSolver::WhichContainer::Bulk));
      m_amr->interpGhostMG(m_particle_scratch1, m_particleRealm, m_phase);
      m_amr->computeGradient(m_particle_scratchD, m_particle_scratch1, m_particleRealm, m_phase);
      DataOps::scale(m_particle_scratchD, -1.0); // scratchD = -grad(D*n)

      DataOps::dotProduct(m_particle_scratch1, m_particle_scratchD, m_particle_E); // m_particle_scratch1 = -E*grad(D*n)
      DataOps::incr(m_energy_sources[idx], m_particle_scratch1, 1.0);              // a_source[idx]
    }

    if (q != 0 && (solver->isMobile() || solver->isDiffusive())) {
      DataOps::scale(m_energy_sources[idx], Abs(q) * Units::Qe);
    }
  }
}

void
ItoPlasmaStepper::computeEdotJSourceNWO()
{
  CH_TIME("ItoPlasmaStepper::computeEdotJSourceNWO()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeEdotJSourceNWO()" << endl;
  }

  DataOps::setValue(m_EdotJ, 0.0);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();
    const int q   = species->getChargeNumber();

    // Do mobile contribution. Computes Z*e*E*mu*n*E*E
    if (q != 0 && solver->isMobile()) {
      solver->depositConductivity(m_particle_scratch1,
                                  solver->getParticles(ItoSolver::WhichContainer::Bulk)); // Deposit mu*n
      m_fluid_scratch1.copy(m_particle_scratch1);                                         // Copy mu*n to fluid Realm
      DataOps::copy(m_fluid_scratchD, m_fluid_E);                                         // m_fluid_scratchD = E
      DataOps::multiplyScalar(m_fluid_scratchD, m_fluid_scratch1);                        // m_fluid_scratchD = E*mu*n
      DataOps::dotProduct(m_fluid_scratch1, m_fluid_E, m_fluid_scratchD); // m_particle_scratch1 = E.dot.(E*mu*n)
      DataOps::scale(m_fluid_scratch1, Abs(q) * Units::Qe);               // m_particle_scratch1 = Z*e*mu*n*E*E

      m_amr->conservativeAverage(m_fluid_scratch1, m_fluidRealm, m_phase);
      m_amr->interpGhost(m_fluid_scratch1, m_fluidRealm, m_phase);
      DataOps::plus(m_EdotJ, m_fluid_scratch1, 0, idx, 1); // a_source[idx] += Z*e*mu*n*E*E
    }

    // Diffusive contribution. Computes -Z*e*E*grad(D*n)
    if (q != 0 && solver->isDiffusive()) {
      solver->depositDiffusivity(m_particle_scratch1,
                                 solver->getParticles(ItoSolver::WhichContainer::Bulk)); // Deposit D*n
      m_fluid_scratch1.copy(m_particle_scratch1);                                        // Copy D*n to fluid Realm
      m_amr->interpGhostMG(m_fluid_scratch1, m_fluidRealm, m_phase);
      m_amr->computeGradient(m_fluid_scratchD, m_fluid_scratch1, m_fluidRealm, m_phase); // scratchD = grad(D*n)
      DataOps::scale(m_fluid_scratchD, -1.0);                                            // scratchD = -grad(D*n)
      DataOps::dotProduct(m_fluid_scratch1, m_fluid_scratchD, m_fluid_E);                // scratch1 = -E.dot.grad(D*n)
      DataOps::scale(m_fluid_scratch1, Abs(q) * Units::Qe);                              // scratch1 = -Z*e*E*grad(D*n)

      m_amr->conservativeAverage(m_fluid_scratch1, m_fluidRealm, m_phase);
      m_amr->interpGhost(m_fluid_scratch1, m_fluidRealm, m_phase);

      DataOps::plus(m_EdotJ, m_fluid_scratch1, 0, idx, 1); // source  += -Z*e*E*grad(D*n)
    }
  }
}

void
ItoPlasmaStepper::computeEdotJSourceNWO2(const Real a_dt)
{
  CH_TIME("ItoPlasmaStepper::computeEdotJSourceNWO2(a_dt)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeEdotJSourceNWO2(a_dt)" << endl;
  }

  DataOps::setValue(m_EdotJ, 0.0);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();
    const int q   = species->getChargeNumber();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    ParticleContainer<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk);

    const DepositionType deposition = solver->getDeposition();

    if ((mobile || diffusive) && q != 0) {

      // We will interpolate m_particle_E onto particle velocity vectors.
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[lvl];

        for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
          const EBCellFAB& E       = (*m_particle_E[lvl])[dit()];
          const EBISBox&   ebisbox = E.getEBISBox();
          const FArrayBox& Efab    = E.getFArrayBox();
          const RealVect   dx      = m_amr->getDx()[lvl] * RealVect::Unit;
          const RealVect   origin  = m_amr->getProbLo();
          const Box        box     = dbl[dit()];

          List<ItoParticle>& particleList = particles[lvl][dit()].listItems();

          // This interpolates the velocity function on to the particle velocities
#if 1
          MayDay::Warning("EBParticleMesh should be replaced with call to AmrMesh as in ItoSolver");
#endif
          EBParticleMesh meshInterp(box, ebisbox, dx, origin);
          //	  meshInterp.interpolateVelocity(particleList, Efab, deposition);
          meshInterp.interpolate<ItoParticle, &ItoParticle::velocity>(particleList, E, deposition, true);

          // Go through the particles and set their mass to E.dot(X^new - X^old)
          for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
            ItoParticle& p = lit();

            const Real      m    = p.weight();
            const RealVect& v    = p.velocity(); // Actually = E(X^new)
            const RealVect& Xnew = p.position();
            const RealVect& Xold = p.oldPosition();

            p.tmp()    = m;
            p.weight() = m * PolyGeom::dot(v, Xnew - Xold);
          }
        }
      }

      // Deposit the result
      solver->depositParticles<ItoParticle, &ItoParticle::weight>(m_particle_scratch1, particles);
      m_fluid_scratch1.copy(m_particle_scratch1);

      // Scale by Qe/dt to make it Joule/dt. Then add to correct index
      DataOps::scale(m_fluid_scratch1, q * Units::Qe / a_dt);
      DataOps::plus(m_EdotJ, m_fluid_scratch1, 0, idx, 1);

      // Set p.weight() back to the original value
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[lvl];
        for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
          List<ItoParticle>& particleList = particles[lvl][dit()].listItems();

          for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
            ItoParticle& p = lit();
            p.weight()     = p.tmp();
          }
        }
      }
    }
  }
}

#include <CD_NamespaceFooter.H>
