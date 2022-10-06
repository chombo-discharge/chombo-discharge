/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaStepper.cpp
  @brief  Implementation of CD_ItoPlasmaStepper.H
  @author Robert Marskar
*/

// Std includes
#include <limits>

// Chombo includes
#include <EBArith.H>
#include <PolyGeom.H>
#include <ParmParse.H>

// Our includes
#include <CD_ItoPlasmaStepper.H>
#include <CD_ParticleOps.H>
#include <CD_DataOps.H>
#include <CD_ParallelOps.H>
#include <CD_Units.H>
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoPlasma;

ItoPlasmaStepper::ItoPlasmaStepper()
{
  CH_TIME("ItoPlasmaStepper::ItoPlasmaStepper");

  // Default settings.
  m_verbosity               = -1;
  m_profile                 = false;
  m_name                    = "ItoPlasmaStepper";
  m_plasmaPhase             = phase::gas;
  m_dt                      = 0.0;
  m_time                    = 0.0;
  m_timeStep                = 0;
  m_loadPerCell             = 1.0;
  m_useNewReactionAlgorithm = true;
  m_regridSuperparticles    = true;
  m_fluidRealm              = Realm::Primal;
  m_particleRealm           = Realm::Primal;
  m_advectionCFL            = 1.0;
  m_diffusionCFL            = std::numeric_limits<Real>::max();
  m_advectionDiffusionCFL   = std::numeric_limits<Real>::max();
  m_relaxTimeFactor         = std::numeric_limits<Real>::max();
  m_minDt                   = std::numeric_limits<Real>::min();
  m_maxDt                   = std::numeric_limits<Real>::max();

  this->parseOptions();
}

ItoPlasmaStepper::ItoPlasmaStepper(RefCountedPtr<ItoPlasmaPhysics>& a_physics) : ItoPlasmaStepper()
{
  CH_TIME("ItoPlasmaStepper::ItoPlasmaStepper(RefCountrPtr<ItoPlasmaPhysics>)");

  m_physics = a_physics;
}

ItoPlasmaStepper::~ItoPlasmaStepper() { CH_TIME("ItoPlasmaStepper::~ItoPlasmaStepper"); }

void
ItoPlasmaStepper::parseOptions()
{
  CH_TIME("ItoPlasmaStepper::parseOptions");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::parseOptions" << endl;
  }

  this->parseVerbosity();
  this->parseSuperparticles();
  this->parseDualGrid();
  this->parseLoadBalance();
  this->parseTimeStepRestrictions();
}

void
ItoPlasmaStepper::parseRuntimeOptions()
{
  CH_TIME("ItoPlasmaStepper::parseRuntimeOptions");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::parseRuntimeOptions" << endl;
  }

  this->parseVerbosity();
  this->parseSuperparticles();
  this->parseLoadBalance();
  this->parseTimeStepRestrictions();
}

void
ItoPlasmaStepper::parseVerbosity() noexcept
{
  CH_TIME("ItoPlasmaStepper::parseVerbosity");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::parseVerbosity" << endl;
  }

  ParmParse pp(m_name.c_str());

  pp.get("verbosity", m_verbosity);
  pp.get("profile", m_profile);
}

void
ItoPlasmaStepper::parseSuperparticles() noexcept
{
  CH_TIME("ItoPlasmaStepper::parseSuperparticles");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::parseSuperparticles" << endl;
  }

  ParmParse pp(m_name.c_str());

  pp.get("particles_per_cell", m_particlesPerCell);
  pp.get("merge_interval", m_mergeInterval);
  pp.get("regrid_superparticles", m_regridSuperparticles);

  if (m_particlesPerCell <= 0) {
    MayDay::Error("ItoPlasmaStepper::parseSuperparticles -- must have 'particles_per_cell' > 0");
  }

  if (m_mergeInterval <= 1) {
    m_mergeInterval = 1;
  }
}

void
ItoPlasmaStepper::parseDualGrid() noexcept
{
  CH_TIME("ItoPlasmaStepper::parseDualGrid");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::parseDualGrid" << endl;
  }

  ParmParse pp(m_name.c_str());

  bool dualGrid = false;
  pp.get("dual_grid", dualGrid);

  if (dualGrid) {
    m_particleRealm = "ParticleRealm";

    CH_assert(m_particleRealm != m_fluidRealm);
  }
  else {
    m_particleRealm = m_fluidRealm;
  }
}

void
ItoPlasmaStepper::parseLoadBalance() noexcept
{
  CH_TIME("ItoPlasmaStepper::parseLoadBalance");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::parseLoadBalance" << endl;
  }

  ParmParse pp(m_name.c_str());

  std::string str;

  pp.get("load_balance", m_loadBalance);
  pp.get("load_index", m_loadBalanceIndex);
  pp.get("load_per_cell", m_loadPerCell);

  // Box sorting for load balancing
  pp.get("box_sorting", str);
  if (str == "none") {
    m_boxSort = BoxSorting::None;
  }
  else if (str == "std") {
    m_boxSort = BoxSorting::Std;
  }
  else if (str == "shuffle") {
    m_boxSort = BoxSorting::Shuffle;
  }
  else if (str == "morton") {
    m_boxSort = BoxSorting::Morton;
  }
  else {
    const std::string err = "ItoPlasmaStepper::parseLoadBalance - 'box_sorting = " + str + "' not recognized";

    MayDay::Error(err.c_str());
  }
}

void
ItoPlasmaStepper::parseTimeStepRestrictions() noexcept
{
  CH_TIME("ItoPlasmaStepper::parseTimeStepRestrictions");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::parseTimeStepRestrictions" << endl;
  }

  ParmParse pp(m_name.c_str());

  pp.get("advection_cfl", m_advectionCFL);
  pp.get("diffusion_cfl", m_diffusionCFL);
  pp.get("advection_diffusion_cfl", m_advectionDiffusionCFL);
  pp.get("relax_dt", m_relaxTimeFactor);
  pp.get("min_dt", m_minDt);
  pp.get("max_dt", m_maxDt);

  if (m_relaxTimeFactor <= 0.0) {
    MayDay::Error("ItoPlasmaStepper::parseTimeStepRestrictions() - must have relax_dt > 0.0");
  }

  if (m_minDt < 0.0) {
    MayDay::Error("ItoPlasmaStepper::parseTimeStepRestrictions() - must have min_dt >= 0.0");
  }

  if (m_maxDt < 0.0) {
    MayDay::Error("ItoPlasmaStepper::parseTimeStepRestrictions() - must have max_dt >= 0.0");
  }

  if (m_advectionCFL <= 0.0) {
    MayDay::Error("ItoPlasmaStepper::parseTimeStepRestrictions() - must have advection_cfl >= 0.0");
  }

  if (m_diffusionCFL <= 0.0) {
    MayDay::Error("ItoPlasmaStepper::parseTimeStepRestrictions() - must have diffusion_cfl >= 0.0");
  }

  if (m_advectionDiffusionCFL <= 0.0) {
    MayDay::Error("ItoPlasmaStepper::parseTimeStepRestrictions() - must have advection_diffusion_cfl >= 0.0");
  }
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
}

void
ItoPlasmaStepper::setupIto()
{
  CH_TIME("ItoPlasmaStepper::setupIto");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::setupIto" << endl;
  }

  m_ito->parseOptions();
  m_ito->setAmr(m_amr);
  m_ito->setPhase(m_plasmaPhase);
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
  m_fieldSolver->setVoltage(m_voltage);
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
  m_rte->setPhase(m_plasmaPhase);
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
  m_sigma->setPhase(m_plasmaPhase);
  m_sigma->setName("Surface charge");
  m_sigma->setTime(0, 0.0, 0.0);
}

void
ItoPlasmaStepper::allocate()
{
  CH_TIME("ItoPlasmaStepper::allocate");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::allocate" << endl;
  }

  // Solver allocation.
  m_ito->allocateInternals();
  m_rte->allocateInternals();
  m_fieldSolver->allocateInternals();
  m_sigma->allocate();

  this->allocateInternals();
}

void
ItoPlasmaStepper::allocateInternals()
{
  CH_TIME("ItoPlasmaStepper::allocateInternals");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::allocateInternals" << endl;
  }

  const int numPlasmaSpecies = m_physics->getNumItoSpecies();
  const int numPhotonSpecies = m_physics->getNumRtSpecies();

  if (numPhotonSpecies <= 0) {
    MayDay::Warning("ItoPlasmaStepper::allocate -- how to handle case with no photon species?");
  }
  if (numPlasmaSpecies <= 0) {
    MayDay::Warning("ItoPlasmaStepper::allocate -- how to handle case with no plasma species?");
  }

  // Scratch data.
  m_amr->allocate(m_fluidScratch1, m_fluidRealm, m_plasmaPhase, 1);
  m_amr->allocate(m_fluidScratchD, m_fluidRealm, m_plasmaPhase, SpaceDim);

  m_amr->allocate(m_particleScratch1, m_particleRealm, m_plasmaPhase, 1);
  m_amr->allocate(m_particleScratchD, m_particleRealm, m_plasmaPhase, SpaceDim);

  // Allocate for energy sources
  m_energySources.resize(numPlasmaSpecies);
  for (int i = 0; i < m_energySources.size(); i++) {
    m_amr->allocate(m_energySources[i], m_particleRealm, m_plasmaPhase, 1);
  }

  // Conductivity things. Defined on the fluid realm.
  m_amr->allocate(m_conductivityCell, m_fluidRealm, m_plasmaPhase, 1);
  m_amr->allocate(m_conductivityFace, m_fluidRealm, m_plasmaPhase, 1);
  m_amr->allocate(m_conductivityEB, m_fluidRealm, m_plasmaPhase, 1);

  // Electric field data.
  m_amr->allocate(m_electricFieldParticle, m_particleRealm, m_plasmaPhase, SpaceDim);
  m_amr->allocate(m_electricFieldFluid, m_fluidRealm, m_plasmaPhase, SpaceDim);

  // Particles and photons per cell on the realms.
  m_amr->allocate(m_particlePPC, m_particleRealm, m_plasmaPhase, numPlasmaSpecies);
  m_amr->allocate(m_particleEPS, m_particleRealm, m_plasmaPhase, numPlasmaSpecies);
  m_amr->allocate(m_particleOldPPC, m_particleRealm, m_plasmaPhase, numPlasmaSpecies);
  m_amr->allocate(m_particleYPC, m_particleRealm, m_plasmaPhase, numPhotonSpecies);

  m_amr->allocate(m_fluidPPC, m_fluidRealm, m_plasmaPhase, numPlasmaSpecies);
  m_amr->allocate(m_fluidYPC, m_fluidRealm, m_plasmaPhase, numPhotonSpecies);
  m_amr->allocate(m_fluidEPS, m_fluidRealm, m_plasmaPhase, numPlasmaSpecies);

  // Aux things for current density, source energy source terms etc.
  m_amr->allocate(m_currentDensity, m_fluidRealm, m_plasmaPhase, SpaceDim);
  m_amr->allocate(m_EdotJ, m_fluidRealm, m_plasmaPhase, numPlasmaSpecies);
}

void
ItoPlasmaStepper::postInitialize()
{
  CH_TIME("ItoPlasmaStepper::postInitialize");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::postInitialize" << endl;
  }
}

void
ItoPlasmaStepper::initialData()
{
  CH_TIME("ItoPlasmaStepper::initialData");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::initialData" << endl;
  }

  m_fieldSolver->setPermittivities();
  m_ito->initialData();
  m_rte->initialData();
  this->initialSigma();

  // Make superparticles.
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

  const RealVect probLo = m_amr->getProbLo();

  EBAMRIVData& sigma = m_sigma->getPhi();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_sigma->getRealm())[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_sigma->getRealm(), m_sigma->getPhase())[lvl];
    const Real               dx    = m_amr->getDx()[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      BaseIVFAB<Real>& phi     = (*sigma[lvl])[dit()];
      const EBISBox&   ebisbox = ebisl[dit()];

      CH_assert(phi.nComp() == 1);

      auto kernel = [&](const VolIndex& vof) -> void {
        const RealVect pos = probLo + Location::position(Location::Cell::Boundary, vof, ebisbox, dx);

        phi(vof, 0) = m_physics->initialSigma(m_time, pos);
      };

      VoFIterator& vofit = (*m_amr->getVofIterator(m_sigma->getRealm(), m_sigma->getPhase())[lvl])[dit()];

      BoxLoops::loop(vofit, kernel);
    }
  }

  // Coarsen throughout the AMR hierarchy.
  m_amr->conservativeAverage(sigma, m_fluidRealm, m_sigma->getPhase());

  // Set surface charge to zero on electrode cut-cells.
  m_sigma->resetElectrodes(sigma, 0.0);
}

void
ItoPlasmaStepper::postCheckpointSetup()
{
  CH_TIME("ItoPlasmaStepper::postCheckpointSetup");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::postCheckpointSetup" << endl;
  }

  // Allocate internal storage.
  this->allocateInternals();

  MayDay::Warning("ItoPlasmaStepper::postCheckpointSetup -- check if remap is necessary");
  m_ito->remap();

  // Recompute the electric field.
  this->postCheckpointPoisson();

  // Compute velocities and diffusion coefficients so we're prepared for the next time step.
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

  // Update ghost cells and re-compute the electric field from the HDF5 data.
  MFAMRCellData& potential = m_fieldSolver->getPotential();

  m_amr->conservativeAverage(potential, m_fluidRealm);
  m_amr->interpGhostMG(potential, m_fluidRealm);

  m_fieldSolver->computeElectricField();

  // Fetch the electric field data on the plasma phase.
  const EBAMRCellData E = m_amr->alias(m_plasmaPhase, m_fieldSolver->getElectricField());

  // Copy onto the storage holding the electric field on the fluid realm. Then interpolate to centroids.
  m_electricFieldFluid.copy(E);
  m_amr->conservativeAverage(m_electricFieldFluid, m_fluidRealm, m_plasmaPhase);
  m_amr->interpGhostPwl(m_electricFieldFluid, m_fluidRealm, m_plasmaPhase);
  m_amr->interpToCentroids(m_electricFieldFluid, m_fluidRealm, m_plasmaPhase);

  // Copy onto the storage holding the electric field on the particle realm.
  m_electricFieldParticle.copy(E);
  m_amr->conservativeAverage(m_electricFieldParticle, m_particleRealm, m_plasmaPhase);
  m_amr->interpGhostPwl(m_electricFieldParticle, m_particleRealm, m_plasmaPhase);
  m_amr->interpToCentroids(m_electricFieldParticle, m_particleRealm, m_plasmaPhase);
}

#ifdef CH_USE_HDF5
void
ItoPlasmaStepper::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const
{
  CH_TIME("ItoPlasmaStepper::writeCheckpointData");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::writeCheckpointData" << endl;
  }

  for (ItoIterator<ItoSolver> solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->writeCheckpointLevel(a_handle, a_lvl);
  }

  for (RtIterator<McPhoto> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->writeCheckpointLevel(a_handle, a_lvl);
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

  for (ItoIterator<ItoSolver> solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->readCheckpointLevel(a_handle, a_lvl);
  }

  for (RtIterator<McPhoto> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->readCheckpointLevel(a_handle, a_lvl);
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
  for (ItoIterator<ItoSolver> solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>& solver = solverIt();

    a_plotVariableNames.append(solver->getPlotVariableNames());
    solver->writePlotData(a_output, a_icomp);
  }

  // RTE solvers copy their output data
  for (RtIterator<McPhoto> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<McPhoto>& solver = solverIt();

    a_plotVariableNames.append(solver->getPlotVariableNames());
    solver->writePlotData(a_output, a_icomp);
  }

  // Write the current to the output
  const Interval srcInterv(0, SpaceDim - 1);
  const Interval dstInterv(a_icomp, a_icomp + SpaceDim - 1);
  a_output.copy(srcInterv, m_currentDensity, dstInterv);

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
ItoPlasmaStepper::writeNumParticlesPerPatch(EBAMRCellData& a_output, int& a_icomp) const
{
  CH_TIME("ItoPlasmaStepper::writeNumParticlesPerPatch");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::writeNumParticlesPerPatch" << endl;
  }

  DataOps::setValue(m_particleScratch1, 0.0);

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    const ParticleContainer<ItoParticle>& particles = solverIt()->getParticles(ItoSolver::WhichContainer::Bulk);

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      for (DataIterator dit(m_amr->getGrids(m_particleRealm)[lvl]); dit.ok(); ++dit) {
        (*m_particleScratch1[lvl])[dit()] += particles[lvl][dit].numItems();
      }
    }
  }

  a_output.copy(Interval(0, 0), m_particleScratch1, Interval(a_icomp, a_icomp));

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

  const Real Emax = this->computeMaxElectricField(m_plasmaPhase);

  const size_t localParticlesBulk    = m_ito->getNumParticles(ItoSolver::WhichContainer::Bulk, true);
  const size_t globalParticlesBulk   = m_ito->getNumParticles(ItoSolver::WhichContainer::Bulk, false);
  const size_t localParticlesEB      = m_ito->getNumParticles(ItoSolver::WhichContainer::EB, true);
  const size_t globalParticlesEB     = m_ito->getNumParticles(ItoSolver::WhichContainer::EB, false);
  const size_t localParticlesDomain  = m_ito->getNumParticles(ItoSolver::WhichContainer::Domain, true);
  const size_t globalParticlesDomain = m_ito->getNumParticles(ItoSolver::WhichContainer::Domain, false);
  const size_t localParticlesSource  = m_ito->getNumParticles(ItoSolver::WhichContainer::Source, true);
  const size_t globalParticlesSource = m_ito->getNumParticles(ItoSolver::WhichContainer::Source, false);

  // Compute some global particle statistics
  Real avgParticles = 0.0;
  Real stdDev       = 0.0;

  Real minParticles = 0.0;
  Real maxParticles = 0.0;

  int minRank = 0;
  int maxRank = 0;

  this->getParticleStatistics(avgParticles, stdDev, minParticles, maxParticles, minRank, maxRank);

  // How was the time step restricted
  std::string str;
  switch (m_timeCode) {
  case TimeCode::Physics: {
    str = "dt restricted by 'physics'";

    break;
  }
  case TimeCode::Advection: {
    str = "dt restricted by 'advection'";

    break;
  }
  case TimeCode::RelaxationTime: {
    str = "dt restricted by 'relaxation time'";

    break;
  }
  case TimeCode::Hardcap: {
    str = "dt restricted by 'hardcap'";

    break;
  }
  default: {
    str = "dt restricted by 'unspecified'";

    break;
  }
  }

  // Print the step report.
  pout() << "                                   " + str << endl;
  pout() << "                                   Emax      = " << Emax << endl
         << "                                   #part     = " << localParticlesBulk << " (" << globalParticlesBulk
         << ")" << endl
         << "                                   #eb part  = " << localParticlesEB << " (" << globalParticlesEB << ")"
         << endl
         << "                                   #dom part = " << localParticlesDomain << " (" << globalParticlesDomain
         << ")" << endl
         << "                                   #src part = " << localParticlesSource << " (" << globalParticlesSource
         << ")" << endl
         << "                                   #part min = " << minParticles << " (on rank = " << minRank << ")"
         << endl
         << "                                   #part max = " << maxParticles << " (on rank = " << maxRank << ")"
         << endl
         << "                                   #part avg = " << avgParticles << endl
         << "                                   #part dev = " << stdDev << " (" << 100. * stdDev / avgParticles << "%)"
         << endl;
}

void
ItoPlasmaStepper::getParticleStatistics(Real& a_avgParticles,
                                        Real& a_sigma,
                                        Real& a_minParticles,
                                        Real& a_maxParticles,
                                        int&  a_minRank,
                                        int&  a_maxRank)
{
  CH_TIME("ItoPlasmaStepper::getParticleStatistics");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::getParticleStatistics" << endl;
  }

  // TLDR: We compute the number of particles, the standard deviation of the number of particles, as well
  //       as the ranks having the smallest/largest number of particles.

  const Real numParticles = 1.0 * m_ito->getNumParticles(ItoSolver::WhichContainer::Bulk, true);

  const std::pair<Real, int> minParticles = ParallelOps::minRank(numParticles);
  const std::pair<Real, int> maxParticles = ParallelOps::maxRank(numParticles);

  a_avgParticles = ParallelOps::average(numParticles);
  a_sigma        = ParallelOps::standardDeviation(numParticles);

  a_minParticles = minParticles.first;
  a_maxParticles = maxParticles.first;

  a_minRank = minParticles.second;
  a_maxRank = maxParticles.second;
}

Real
ItoPlasmaStepper::computeDt()
{
  CH_TIME("ItoPlasmaStepper::computeDt");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeDt" << endl;
  }

  Real dt = std::numeric_limits<Real>::max();

  // Compute various time steps.
  m_advectionDt          = m_ito->computeAdvectiveDt();
  m_diffusionDt          = m_ito->computeDiffusiveDt();
  m_advectionDiffusionDt = m_ito->computeDt();
  m_physicsDt            = this->computePhysicsDt();
  m_relaxationTime       = this->computeRelaxationTime();

  if (m_advectionCFL * m_advectionDt < dt) {
    dt         = m_advectionCFL * m_advectionDt;
    m_timeCode = TimeCode::Advection;
  }

  if (m_diffusionCFL * m_diffusionDt < dt) {
    dt         = m_diffusionCFL * m_diffusionDt;
    m_timeCode = TimeCode::Diffusion;
  }

  if (m_advectionDiffusionCFL * m_advectionDiffusionDt < dt) {
    dt         = m_advectionDiffusionCFL * m_advectionDiffusionDt;
    m_timeCode = TimeCode::AdvectionDiffusion;
  }

  if (m_physicsDt < dt) {
    dt         = m_physicsDt;
    m_timeCode = TimeCode::Physics;
  }

  if (m_minDt > dt) {
    dt         = m_minDt;
    m_timeCode = TimeCode::Hardcap;
  }

  if (m_maxDt < dt) {
    dt         = m_maxDt;
    m_timeCode = TimeCode::Hardcap;
  }

  return dt;
}

void
ItoPlasmaStepper::registerRealms()
{
  CH_TIME("ItoPlasmaStepper::registerRealms");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::registerRealms" << endl;
  }

  // TLDR: If using dual grid then m_particleRealm != m_fluidRealm and we'll have two realms.

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
{
  CH_TIME("ItoPlasmaStepper::postRegrid");
}

int
ItoPlasmaStepper::getNumberOfPlotVariables() const
{
  CH_TIME("ItoPlasmaStepper::getNumberOfPlotVariables");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::getNumberOfPlotVariables" << endl;
  }

  int ncomp = 0;

  // Ito solver variables.
  for (ItoIterator<ItoSolver> solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    ncomp += solverIt()->getNumberOfPlotVariables();
  }

  // RTE solver variables.
  for (RtIterator<McPhoto> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    ncomp += solverIt()->getNumberOfPlotVariables();
  }

  // Field solver variables.
  ncomp += m_fieldSolver->getNumberOfPlotVariables();

  // Surface charge solver variables.
  ncomp += m_sigma->getNumberOfPlotVariables();

  // Current density.
  ncomp += SpaceDim;

  // Number of particles per cell.
  ncomp += 1;

  return ncomp;
}

void
ItoPlasmaStepper::setIto(RefCountedPtr<ItoLayout<ItoSolver>>& a_ito) noexcept
{
  CH_TIME("ItoPlasmaStepper::setIto");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::setIto" << endl;
  }

  m_ito = a_ito;
}

void
ItoPlasmaStepper::setFieldSolver(RefCountedPtr<FieldSolver>& a_fieldSolver) noexcept
{
  CH_TIME("ItoPlasmaStepper::setFieldSolver");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::setFieldSolver" << endl;
  }

  m_fieldSolver = a_fieldSolver;
}

void
ItoPlasmaStepper::setRadiativeTransferSolvers(RefCountedPtr<RtLayout<McPhoto>>& a_rteLayout) noexcept
{
  CH_TIME("ItoPlasmaStepper::setRadiativeTransferSolvers");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::setRadiativeTransferSolvers" << endl;
  }

  m_rte = a_rteLayout;
}

void
ItoPlasmaStepper::setVoltage(const std::function<Real(const Real a_time)>& a_voltage) noexcept
{
  CH_TIME("ItoPlasmaStepper::setVoltage");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::setVoltage" << endl;
  }

  m_voltage = a_voltage;
}

Real
ItoPlasmaStepper::computeMaxElectricField(const phase::which_phase a_phase) noexcept
{
  CH_TIME("ItoPlasmaStepper::computeMaxElectricField");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeMaxElectricField" << endl;
  }

  // Get a handle to the E-field. Note that this is the cell-centered field!
  const EBAMRCellData cellCenteredE = m_amr->alias(a_phase, m_fieldSolver->getElectricField());

  // Interpolate to centroids
  EBAMRCellData centroidCenteredE;
  m_amr->allocate(centroidCenteredE, m_fluidRealm, a_phase, SpaceDim);

  DataOps::copy(centroidCenteredE, cellCenteredE);

  m_amr->interpToCentroids(centroidCenteredE, m_fluidRealm, m_plasmaPhase);

  Real max = 0.0;
  Real min = 0.0;

  DataOps::getMaxMinNorm(max, min, centroidCenteredE);

  return max;
}

void
ItoPlasmaStepper::computeElectricField(EBAMRCellData& a_electricField, const phase::which_phase a_phase) const noexcept
{
  CH_TIME("ItoPlasmaStepper::computeElectricField(EBAMRCellData, phase)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeElectricField(EBAMRCellData, phase)" << endl;
  }

  CH_assert(a_electricField.getRealm() == m_fluidRealm);

  m_fieldSolver->computeElectricField(a_electricField, a_phase, m_fieldSolver->getPotential());
}

Real
ItoPlasmaStepper::getTime() const noexcept
{
  CH_TIME("ItoPlasmaStepper::getTime");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::getTime" << endl;
  }

  return m_time;
}

void
ItoPlasmaStepper::computeSpaceChargeDensity() noexcept
{
  CH_TIME("ItoPlasmaStepper::computeSpaceChargeDensity()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeSpaceChargeDensity()" << endl;
  }

  this->computeSpaceChargeDensity(m_fieldSolver->getRho(), m_ito->getDensities());
}

void
ItoPlasmaStepper::computeSpaceChargeDensity(MFAMRCellData& a_rho, const Vector<EBAMRCellData*>& a_densities) noexcept
{
  CH_TIME("ItoPlasmaStepper::computeSpaceChargeDensity(rho, densities)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeSpaceChargeDensity(rho, densities)" << endl;
  }

  // TLDR: a_densities is from the Ito solvers so it is defined over the particle realm. This CAN be equal to
  //       the same realm on which a_rho is defined, but not necessarily. So, we use m_fluidScratch1 as temporary
  //       scratch storage.

  // Reset
  DataOps::setValue(a_rho, 0.0);

  // Alias for the plasma phase.
  EBAMRCellData rhoPhase = m_amr->alias(m_plasmaPhase, a_rho);

  // Increment each solver
  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    const RefCountedPtr<ItoSolver>&  solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();
    const int                        idx     = solverIt.index();
    const int                        q       = species->getChargeNumber();

    if (species->getChargeNumber() != 0) {
      m_fluidScratch1.copy(*a_densities[idx]);
      DataOps::incr(rhoPhase, m_fluidScratch1, q);
    }
  }

  DataOps::scale(a_rho, Units::Qe);

  m_amr->conservativeAverage(a_rho, m_fluidRealm);
  m_amr->interpGhost(a_rho, m_fluidRealm);

  // Interpolate to centroids.
  m_amr->interpToCentroids(rhoPhase, m_fluidRealm, m_plasmaPhase);
}

void
ItoPlasmaStepper::computeConductivityCell(EBAMRCellData& a_conductivity) noexcept
{
  CH_TIME("ItoPlasmaStepper::computeConductivityCell(EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeConductivityCell(EBAMRCellData)" << endl;
  }

  this->computeConductivityCell(a_conductivity, m_ito->getParticles(ItoSolver::WhichContainer::Bulk));
}

void
ItoPlasmaStepper::computeConductivityCell(EBAMRCellData&                                 a_conductivity,
                                          const Vector<ParticleContainer<ItoParticle>*>& a_particles) noexcept
{
  CH_TIME("ItoPlasmaStepper::computeConductivityCell(EBAMRCellData, Particles)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeConductivityCell(EBAMRCellData, Particles)" << endl;
  }

  // TLDR: This will deposit the particle conductivity on the mesh (onto m_particleScratch1) which is
  //       then added to the total conductivity.

  DataOps::setValue(a_conductivity, 0.0);

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();
    const int q   = species->getChargeNumber();

    if (q != 0 && solver->isMobile()) {
      solver->depositConductivity(m_particleScratch1, *a_particles[idx]);

      // Add to the fluid realm.
      m_fluidScratch1.copy(m_particleScratch1);
      DataOps::incr(a_conductivity, m_fluidScratch1, std::abs(q));
    }
  }

  DataOps::scale(a_conductivity, Units::Qe);

  m_amr->conservativeAverage(a_conductivity, m_fluidRealm, m_plasmaPhase);
  m_amr->interpGhostPwl(a_conductivity, m_fluidRealm, m_plasmaPhase);

  // Interpolate to centroids.
  m_amr->interpToCentroids(a_conductivity, m_fluidRealm, m_plasmaPhase);
}

void
ItoPlasmaStepper::computeCurrentDensity(EBAMRCellData& a_J) noexcept
{
  CH_TIME("ItoPlasmaStepper::computeCurrentDensity(EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeCurrentDensity(EBAMRCellData)" << endl;
  }

  CH_assert(a_J[0]->nComp() == SpaceDim);

  // TLDR: a_J is defined over the fluid Realm but the computation takes place on the particle Realm.
  //       If the Realms are different we compute on a scratch storage instead

  this->computeConductivityCell(m_fluidScratch1);
  DataOps::copy(a_J, m_electricFieldFluid);

  DataOps::multiplyScalar(a_J, m_fluidScratch1);
}

Real
ItoPlasmaStepper::computeRelaxationTime() noexcept
{
  CH_TIME("ItoPlasmaStepper::computeRelaxationTime()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeRelaxationTime()" << endl;
  }

  // TLDR: We compute eps0/conductivity directly.

  EBAMRCellData conductivity;
  EBAMRCellData relaxTime;

  m_amr->allocate(conductivity, m_fluidRealm, m_plasmaPhase, 1);
  m_amr->allocate(relaxTime, m_fluidRealm, m_plasmaPhase, 1);

  this->computeConductivityCell(conductivity);

  DataOps::setValue(relaxTime, Units::eps0);
  DataOps::divideFallback(relaxTime, conductivity, std::numeric_limits<Real>::max());

  m_amr->conservativeAverage(relaxTime, m_fluidRealm, m_plasmaPhase);

  Real min = 0.0;
  Real max = 0.0;

  DataOps::getMaxMinNorm(max, min, relaxTime);

  return min;
}

bool
ItoPlasmaStepper::solvePoisson() noexcept
{
  CH_TIME("ItoPlasmaStepper::solvePoisson()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::solvePoisson()" << endl;
  }

  // Compute the space charge density.
  this->computeSpaceChargeDensity();

  // Solve the Poisson equation and compute the cell-centered electric field.
  MFAMRCellData& phi   = m_fieldSolver->getPotential();
  MFAMRCellData& rho   = m_fieldSolver->getRho();
  EBAMRIVData&   sigma = m_sigma->getPhi();

  const bool converged = m_fieldSolver->solve(phi, rho, sigma, false);

  m_fieldSolver->computeElectricField();

  // Copy the electric field to appropriate data holders and perform center-to-centroid
  // interpolation.
  EBAMRCellData E;
  m_amr->allocatePointer(E);
  m_amr->alias(E, m_plasmaPhase, m_fieldSolver->getElectricField());

  // Fluid realm
  m_electricFieldFluid.copy(E);
  m_amr->conservativeAverage(m_electricFieldFluid, m_fluidRealm, m_plasmaPhase);
  m_amr->interpGhostPwl(m_electricFieldFluid, m_fluidRealm, m_plasmaPhase);
  m_amr->interpToCentroids(m_electricFieldFluid, m_fluidRealm, m_plasmaPhase);

  // Particle realm
  m_electricFieldParticle.copy(E);
  m_amr->conservativeAverage(m_electricFieldParticle, m_particleRealm, m_plasmaPhase);
  m_amr->interpGhostPwl(m_electricFieldParticle, m_particleRealm, m_plasmaPhase);
  m_amr->interpToCentroids(m_electricFieldParticle, m_particleRealm, m_plasmaPhase);

  return converged;
}

void
ItoPlasmaStepper::intersectParticles(const SpeciesSubset  a_speciesSubset,
                                     const EBIntersection a_intersectionAlg,
                                     const bool           a_delete) noexcept
{
  CH_TIME("ItoPlasmaStepper::intersectParticles(SpeciesSubset, EBRepresentation, bool)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::intersectParticles(SpeciesSubset, EBRepresentation, bool)" << endl;
  }

  this->intersectParticles(a_speciesSubset,
                           ItoSolver::WhichContainer::Bulk,
                           ItoSolver::WhichContainer::EB,
                           ItoSolver::WhichContainer::Domain,
                           a_intersectionAlg,
                           a_delete);
}

void
ItoPlasmaStepper::intersectParticles(const SpeciesSubset             a_speciesSubset,
                                     const ItoSolver::WhichContainer a_containerBulk,
                                     const ItoSolver::WhichContainer a_containerEB,
                                     const ItoSolver::WhichContainer a_containerDomain,
                                     const EBIntersection            a_intersectionAlg,
                                     const bool                      a_delete) noexcept
{
  CH_TIME("ItoPlasmaStepper::intersectParticles(SpeciesSubset, Containerx3, EBRepresentation, bool)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::intersectParticles(SpeciesSubset, Containerx3, EBRepresentation, bool)" << endl;
  }

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = (species->getChargeNumber() != 0);

    switch (a_speciesSubset) {
    case SpeciesSubset::All: {
      solver->intersectParticles(a_containerBulk, a_containerEB, a_containerDomain, a_intersectionAlg, a_delete);

      break;
    }
    case SpeciesSubset::AllMobile: {
      if (mobile) {
        solver->intersectParticles(a_containerBulk, a_containerEB, a_containerDomain, a_intersectionAlg, a_delete);
      }

      break;
    }
    case SpeciesSubset::AllDiffusive: {
      if (diffusive) {
        solver->intersectParticles(a_containerBulk, a_containerEB, a_containerDomain, a_intersectionAlg, a_delete);
      }

      break;
    }
    case SpeciesSubset::ChargedMobile: {
      if (charged && mobile) {
        solver->intersectParticles(a_containerBulk, a_containerEB, a_containerDomain, a_intersectionAlg, a_delete);
      }

      break;
    }
    case SpeciesSubset::ChargedDiffusive: {
      if (charged && diffusive) {
        solver->intersectParticles(a_containerBulk, a_containerEB, a_containerDomain, a_intersectionAlg, a_delete);
      }

      break;
    }
    case SpeciesSubset::AllMobileOrDiffusive: {
      if (mobile || diffusive) {
        solver->intersectParticles(a_containerBulk, a_containerEB, a_containerDomain, a_intersectionAlg, a_delete);
      }

      break;
    }
    case SpeciesSubset::ChargedAndMobileOrDiffusive: {
      if (charged && (mobile || diffusive)) {
        solver->intersectParticles(a_containerBulk, a_containerEB, a_containerDomain, a_intersectionAlg, a_delete);
      }

      break;
    }
    case SpeciesSubset::Stationary: {
      if (!mobile && !diffusive) {
        solver->intersectParticles(a_containerBulk, a_containerEB, a_containerDomain, a_intersectionAlg, a_delete);
      }

      break;
    }
    default: {
      MayDay::Abort("ItoPlasmaStepper::intersectParticles - logic bust");

      break;
    }
    }
  }
}

void
ItoPlasmaStepper::removeCoveredParticles(const SpeciesSubset    a_speciesSubset,
                                         const EBRepresentation a_representation,
                                         const Real             a_tolerance) noexcept
{
  CH_TIME("ItoPlasmaStepper::removeCoveredParticles(SpeciesSubset, EBRepresentation, Real)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::removeCoveredParticles(SpeciesSubset, EBRepresentation, Real)" << endl;
  }

  this->removeCoveredParticles(a_speciesSubset, ItoSolver::WhichContainer::Bulk, a_representation, a_tolerance);
}

void
ItoPlasmaStepper::removeCoveredParticles(const SpeciesSubset             a_which,
                                         const ItoSolver::WhichContainer a_container,
                                         const EBRepresentation          a_representation,
                                         const Real                      a_tolerance) noexcept
{
  CH_TIME("ItoPlasmaStepper::removeCoveredParticles(SpeciesSubset, container, EBRepresentation, tolerance)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::removeCoveredParticles(SpeciesSubset, container, EBRepresentation, tolerance)" << endl;
  }

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = (species->getChargeNumber() != 0);

    switch (a_which) {
    case SpeciesSubset::All: {
      solver->removeCoveredParticles(a_container, a_representation, a_tolerance);

      break;
    }
    case SpeciesSubset::AllMobile: {
      if (mobile) {
        solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      }

      break;
    }
    case SpeciesSubset::AllDiffusive: {
      if (diffusive) {
        solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      }

      break;
    }
    case SpeciesSubset::ChargedMobile: {
      if (charged && mobile) {
        solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      }

      break;
    }
    case SpeciesSubset::ChargedDiffusive: {
      if (charged && diffusive) {
        solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      }

      break;
    }
    case SpeciesSubset::AllMobileOrDiffusive: {
      if (mobile || diffusive) {
        solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      }

      break;
    }
    case SpeciesSubset::ChargedAndMobileOrDiffusive: {
      if (charged && (mobile || diffusive)) {
        solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      }

      break;
    }
    case SpeciesSubset::Stationary: {
      if (!mobile && !diffusive) {
        solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      }

      break;
    }
    default: {
      MayDay::Abort("ItoPlasmaStepper::removeCoveredParticles - logic bust");

      break;
    }
    }
  }
}

void
ItoPlasmaStepper::transferCoveredParticles(const SpeciesSubset    a_speciesSubset,
                                           const EBRepresentation a_representation,
                                           const Real             a_tolerance) noexcept
{
  CH_TIME("ItoPlasmaStepper::transferCoveredParticles(SpeciesSubset, EBRepresentation, Real)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::transferCoveredParticles(SpeciesSubset, EBRepresentation, Real)" << endl;
  }

  this->transferCoveredParticles(a_speciesSubset,
                                 ItoSolver::WhichContainer::Bulk,
                                 ItoSolver::WhichContainer::Covered,
                                 a_representation,
                                 a_tolerance);
}

void
ItoPlasmaStepper::transferCoveredParticles(const SpeciesSubset             a_speciesSubset,
                                           const ItoSolver::WhichContainer a_containerFrom,
                                           const ItoSolver::WhichContainer a_containerTo,
                                           const EBRepresentation          a_representation,
                                           const Real                      a_tolerance) noexcept
{
  CH_TIME("ItoPlasmaStepper::transferCoveredParticles(SpeciesSubset, Containerx2, EBRepresentation, Real)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::transferCoveredParticles(SpeciesSubset, Containerx2, EBRepresentation, Real)" << endl;
  }

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = (species->getChargeNumber() != 0);

    switch (a_speciesSubset) {
    case SpeciesSubset::All: {
      solver->transferCoveredParticles(a_containerFrom, a_containerTo, a_representation, a_tolerance);

      break;
    }
    case SpeciesSubset::AllMobile: {
      if (mobile) {
        solver->transferCoveredParticles(a_containerFrom, a_containerTo, a_representation, a_tolerance);
      }

      break;
    }
    case SpeciesSubset::AllDiffusive: {
      if (diffusive) {
        solver->transferCoveredParticles(a_containerFrom, a_containerTo, a_representation, a_tolerance);
      }

      break;
    }
    case SpeciesSubset::ChargedMobile: {
      if (charged && mobile) {
        solver->transferCoveredParticles(a_containerFrom, a_containerTo, a_representation, a_tolerance);
      }

      break;
    }
    case SpeciesSubset::ChargedDiffusive: {
      if (charged && diffusive) {
        solver->transferCoveredParticles(a_containerFrom, a_containerTo, a_representation, a_tolerance);
      }

      break;
    }
    case SpeciesSubset::AllMobileOrDiffusive: {
      if (mobile || diffusive) {
        solver->transferCoveredParticles(a_containerFrom, a_containerTo, a_representation, a_tolerance);
      }

      break;
    }
    case SpeciesSubset::ChargedAndMobileOrDiffusive: {
      if (charged && (mobile || diffusive)) {
        solver->transferCoveredParticles(a_containerFrom, a_containerTo, a_representation, a_tolerance);
      }

      break;
    }
    case SpeciesSubset::Stationary: {
      if (!mobile && !diffusive) {
        solver->transferCoveredParticles(a_containerFrom, a_containerTo, a_representation, a_tolerance);
      }

      break;
    }
    default: {
      MayDay::Abort("ItoPlasmaStepper::transferCoveredParticles - logic bust");

      break;
    }
    }
  }
}

void
ItoPlasmaStepper::remapParticles(const SpeciesSubset a_speciesSubset) noexcept
{
  CH_TIME("ItoPlasmaStepper::remapParticles(SpeciesSubset)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::remapParticles(SpeciesSubset)" << endl;
  }

  this->remapParticles(a_speciesSubset, ItoSolver::WhichContainer::Bulk);
}

void
ItoPlasmaStepper::remapParticles(const SpeciesSubset             a_speciesSubset,
                                 const ItoSolver::WhichContainer a_container) noexcept
{
  CH_TIME("ItoPlasmaStepper::remapParticles(SpeciesSubset)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::remapParticles(SpeciesSubset)" << endl;
  }

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = (species->getChargeNumber() != 0);

    switch (a_speciesSubset) {
    case SpeciesSubset::All: {
      solver->remap(a_container);

      break;
    }
    case SpeciesSubset::AllMobile: {
      if (mobile) {
        solver->remap(a_container);
      }

      break;
    }
    case SpeciesSubset::AllDiffusive: {
      if (diffusive) {
        solver->remap(a_container);
      }

      break;
    }
    case SpeciesSubset::ChargedMobile: {
      if (charged && mobile) {
        solver->remap(a_container);
      }

      break;
    }
    case SpeciesSubset::ChargedDiffusive: {
      if (charged && diffusive) {
        solver->remap(a_container);
      }

      break;
    }
    case SpeciesSubset::AllMobileOrDiffusive: {
      if (mobile || diffusive) {
        solver->remap(a_container);
      }

      break;
    }
    case SpeciesSubset::ChargedAndMobileOrDiffusive: {
      if (charged && (mobile || diffusive)) {
        solver->remap(a_container);
      }

      break;
    }
    case SpeciesSubset::Stationary: {
      if (!mobile && !diffusive) {
        solver->remap(a_container);
      }

      break;
    }
    default: {
      MayDay::Abort("ItoPlasmaStepper::remapParticles - logic bust");

      break;
    }
    }
  }
}

void
ItoPlasmaStepper::depositParticles(const SpeciesSubset a_speciesSubset) noexcept
{
  CH_TIME("ItoPlasmaStepper::depositParticles(SpeciesSubset)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::depositParticles(SpeciesSubset)" << endl;
  }

  this->depositParticles(a_speciesSubset, ItoSolver::WhichContainer::Bulk);
}

void
ItoPlasmaStepper::depositParticles(const SpeciesSubset             a_speciesSubset,
                                   const ItoSolver::WhichContainer a_container) noexcept
{
  CH_TIME("ItoPlasmaStepper::depositParticles(SpeciesSubset)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::depositParticles(SpeciesSubset)" << endl;
  }

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = (species->getChargeNumber() != 0);

    switch (a_speciesSubset) {
    case SpeciesSubset::All: {
      solver->depositParticles(a_container);

      break;
    }
    case SpeciesSubset::AllMobile: {
      if (mobile) {
        solver->depositParticles(a_container);
      }

      break;
    }
    case SpeciesSubset::AllDiffusive: {
      if (diffusive) {
        solver->depositParticles(a_container);
      }

      break;
    }
    case SpeciesSubset::ChargedMobile: {
      if (charged && mobile) {
        solver->depositParticles(a_container);
      }

      break;
    }
    case SpeciesSubset::ChargedDiffusive: {
      if (charged && diffusive) {
        solver->depositParticles(a_container);
      }

      break;
    }
    case SpeciesSubset::AllMobileOrDiffusive: {
      if (mobile || diffusive) {
        solver->depositParticles(a_container);
      }

      break;
    }
    case SpeciesSubset::ChargedAndMobileOrDiffusive: {
      if (charged && (mobile || diffusive)) {
        solver->depositParticles(a_container);
      }

      break;
    }
    case SpeciesSubset::Stationary: {
      if (!mobile && !diffusive) {
        solver->depositParticles(a_container);
      }

      break;
    }
    default: {
      MayDay::Abort("ItoPlasmaStepper::depositParticles - logic bust");

      break;
    }
    }
  }
}

void
ItoPlasmaStepper::setItoVelocityFunctions() noexcept
{
  CH_TIME("ItoPlasmaStepper::setItoVelocityFunctions");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::setItoVelocityFunctions" << endl;
  }

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    if (solver->isMobile()) {
      EBAMRCellData& velocityFunction = solver->getVelocityFunction();
      velocityFunction.copy(m_electricFieldParticle);

      const int Z = species->getChargeNumber();

      if (Z < 0) {
        DataOps::scale(velocityFunction, -1.0);
      }
      else if (Z == 0) {
        MayDay::Warning("ItoPlasmaStepper::setItoVelocityFunctions -- what to do about sign for neutral species?");
      }

      // Coarsen and update ghost cells.
      m_amr->conservativeAverage(velocityFunction, m_particleRealm, m_plasmaPhase);
      m_amr->interpGhost(velocityFunction, m_particleRealm, m_plasmaPhase);
    }
  }
}

void
ItoPlasmaStepper::computeItoVelocities() noexcept
{
  CH_TIME("ItoPlasmaStepper::computeItoVelocities()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoVelocities()" << endl;
  }

  const ItoPlasmaPhysics::coupling fieldCoupling = m_physics->getCoupling();

  // Set the ItoSolver velocity functions.
  this->setItoVelocityFunctions();

  // Compute mobilities based on appropriate coupling
  switch (fieldCoupling) {
  case ItoPlasmaPhysics::coupling::LFA: {
    this->computeItoMobilitiesLFA();

    break;
  }
  case ItoPlasmaPhysics::coupling::LEA: {
    this->computeItoMobilitiesLEA();

    break;
  }
  default: {
    MayDay::Error("ItoPlasmaStepper::computeItoVelocities - logic bust");

    break;
  }
  }

  // Interpolate velocity function to particle position so that particles get velocity v = +/- mu*E
  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->interpolateVelocities();
  }
}

void
ItoPlasmaStepper::computeItoDiffusion() noexcept
{
  CH_TIME("ItoPlasmaStepper::computeItoDiffusion()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoDiffusion()" << endl;
  }

  const ItoPlasmaPhysics::coupling fieldCoupling = m_physics->getCoupling();

  // Compute mobilities based on appropriate coupling
  switch (fieldCoupling) {
  case ItoPlasmaPhysics::coupling::LFA: {
    this->computeItoDiffusionLFA();

    break;
  }
  case ItoPlasmaPhysics::coupling::LEA: {
    this->computeItoDiffusionLEA();

    break;
  }
  default: {
    MayDay::Error("ItoPlasmaStepper::computeItoDiffusion - logic bust");

    break;
  }
  }
}

void
ItoPlasmaStepper::computeItoMobilitiesLFA() noexcept
{
  CH_TIME("ItoPlasmaStepper::computeItoMobilitiesLFA()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoMobilitiesLFA()" << endl;
  }

  Vector<EBAMRCellData*> meshMobilities = m_ito->getMobilityFunctions();

  this->computeItoMobilitiesLFA(meshMobilities, m_electricFieldFluid, m_time);
}

void
ItoPlasmaStepper::computeItoMobilitiesLFA(Vector<EBAMRCellData*>& a_meshMobilities,
                                          const EBAMRCellData&    a_electricField,
                                          const Real              a_time) noexcept
{
  CH_TIME("ItoPlasmaStepper::computeItoMobilitiesLFA(mobilities, E, time)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoMobilitiesLFA(mobilities, E, time)" << endl;
  }

  const int numPlasmaSpecies = m_physics->getNumItoSpecies();

  CH_assert(a_electricField.getRealm() == m_fluidRealm);
  CH_assert(a_meshMobilities.size() == numPlasmaSpecies);

  // The mesh mobilities belong on the particle realm (they are the ItoSolver mobilities) but we need to run
  // the computation on the fluid realm. So, create some transient storage for that.
  Vector<EBAMRCellData> fluidScratchMobilities(numPlasmaSpecies);
  for (int i = 0; i < numPlasmaSpecies; i++) {
    m_amr->allocate(fluidScratchMobilities[i], m_fluidRealm, m_plasmaPhase, 1);

    CH_assert(a_meshMobilities[i]->getRealm() == m_particleRealm);
  }

  // Now run the computation on the fluid realm, computing the mobilities into fluidScratchMobilities
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    Vector<LevelData<EBCellFAB>*> mobilities(numPlasmaSpecies);

    for (int i = 0; i < numPlasmaSpecies; i++) {
      mobilities[i] = &(*(fluidScratchMobilities[i])[lvl]);
    }

    // Run the level computation, which will fill mobilities aka fluidScratchMobilities.
    this->computeItoMobilitiesLFA(mobilities, *a_electricField[lvl], lvl, a_time);
  }

  // Copy the fluid realm data into the particle realm data and update ghost cells
  for (int i = 0; i < numPlasmaSpecies; i++) {
    a_meshMobilities[i]->copy(fluidScratchMobilities[i]);

    m_amr->conservativeAverage(*a_meshMobilities[i], m_particleRealm, m_plasmaPhase);
    m_amr->interpGhost(*a_meshMobilities[i], m_particleRealm, m_plasmaPhase);
  }

  // Interpolate mobilities to the particle position.
  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>& solver = solverIt();

    if (solver->isMobile()) {
      solver->interpolateMobilities();
    }
  }
}

void
ItoPlasmaStepper::computeItoMobilitiesLFA(Vector<LevelData<EBCellFAB>*>& a_meshMobilities,
                                          const LevelData<EBCellFAB>&    a_electricField,
                                          const int                      a_level,
                                          const Real                     a_time) noexcept
{
  CH_TIME("ItoPlasmaStepper::computeItoMobilitiesLFA(mobilities, E, level, time)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoMobilitiesLFA(mobilities, E, level, time)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_fluidRealm)[a_level];

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const EBCellFAB& E       = a_electricField[dit()];
    const Box        cellBox = dbl[dit()];

    Vector<EBCellFAB*> meshMobilities;
    for (int i = 0; i < a_meshMobilities.size(); i++) {
      meshMobilities.push_back(&((*a_meshMobilities[i])[dit()]));
    }

    this->computeItoMobilitiesLFA(meshMobilities, E, a_level, dit(), cellBox, a_time);
  }
}

void
ItoPlasmaStepper::computeItoMobilitiesLFA(Vector<EBCellFAB*>& a_meshMobilities,
                                          const EBCellFAB&    a_electricField,
                                          const int           a_level,
                                          const DataIndex     a_dit,
                                          const Box           a_box,
                                          const Real          a_time) noexcept
{
  CH_TIME("ItoPlasmaStepper::computeItoMobilitiesLFA(meshMobilities, E, level, dit, box, time)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoMobilitiesLFA(meshMobilities, E, level, dit, box, time)" << endl;
  }

  // TLDR: We go through each and every cell and call the physics interface. This includes cells covered by a finer grid
  //       but data is coarsened later anyways.

  const Real     dx      = m_amr->getDx()[a_level];
  const RealVect probLo  = m_amr->getProbLo();
  const EBISBox& ebisbox = m_amr->getEBISLayout(m_fluidRealm, m_plasmaPhase)[a_level][a_dit];

  // Handle to regular data.
  const FArrayBox&   electricFieldReg = a_electricField.getFArrayBox();
  Vector<FArrayBox*> meshMobilitiesReg;
  for (int i = 0; i < a_meshMobilities.size(); i++) {
    meshMobilitiesReg.push_back(&(a_meshMobilities[i]->getFArrayBox()));
  }

  // Regular kernel
  auto regularKernel = [&](const IntVect& iv) -> void {
    const RealVect pos = m_amr->getProbLo() + dx * (RealVect(iv) + 0.5 * RealVect::Unit);
    const RealVect E   = RealVect(D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));

    // Call ito_physics and compute mobilities for each particle species
    const Vector<Real> mobilities = m_physics->computeItoMobilitiesLFA(a_time, pos, E);

    // Put mobilities in appropriate data holder
    for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
      const int idx = solverIt.index();

      (*meshMobilitiesReg[idx])(iv, 0) = mobilities[idx];
    }
  };

  // Irregular kernel.
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    const RealVect e   = RealVect(D_DECL(a_electricField(vof, 0), a_electricField(vof, 1), a_electricField(vof, 2)));
    const RealVect pos = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);

    // Compute diffusion
    const Vector<Real> mobilities = m_physics->computeItoMobilitiesLFA(a_time, pos, e);

    // Put diffusion in the appropriate place.
    for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
      const int idx = solverIt.index();

      (*a_meshMobilities[idx])(vof, 0) = mobilities[idx];
    }
  };

  VoFIterator& vofit = (*m_amr->getVofIterator(m_fluidRealm, m_plasmaPhase)[a_level])[a_dit];

  // Run the kernels.
  BoxLoops::loop(a_box, regularKernel);
  BoxLoops::loop(vofit, irregularKernel);

  // Covered is bogus.
  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    const int idx = solverIt.index();
    a_meshMobilities[idx]->setCoveredCellVal(0.0, 0);
  }
}

void
ItoPlasmaStepper::computeItoMobilitiesLEA() noexcept
{
  CH_TIME("ItoPlasmaStepper::computeItoMobilitiesLEA()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoMobilitiesLEA()" << endl;
  }

  // This is really simple because the solvers do this directly...
  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->updateMobilities();
  }
}

void
ItoPlasmaStepper::computeItoDiffusionLFA() noexcept
{
  CH_TIME("ItoPlasmaStepper::computeItoDiffusionLFA()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoDiffusionLFA()" << endl;
  }

  Vector<EBAMRCellData*> diffusionCoefficients = m_ito->getDiffusionFunctions();
  Vector<EBAMRCellData*> densities             = m_ito->getDensities();

  this->computeItoDiffusionLFA(diffusionCoefficients, densities, m_electricFieldFluid, m_time);
}

void
ItoPlasmaStepper::computeItoDiffusionLFA(Vector<EBAMRCellData*>&       a_diffusionCoefficients,
                                         const Vector<EBAMRCellData*>& a_densities,
                                         const EBAMRCellData&          a_electricField,
                                         const Real                    a_time) noexcept
{
  CH_TIME("ItoPlasmaStepper::computeItoDiffusionLFA(Vector<EBAMRCellData>x2, EBAMRCellData, Real)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoDiffusionLFA(Vector<EBAMRCellData>x2, EBAMRCellData, Real)" << endl;
  }

  const int numPlasmaSpecies = m_physics->getNumItoSpecies();

  CH_assert(a_electricField.getRealm() == m_fluidRealm);
  CH_assert(a_diffusionCoefficients.size() == numPlasmaSpecies);
  CH_assert(a_densities.size() == numPlasmaSpecies);

  // The mesh diffusion coefficients belong on the particle realm (they are the ItoSolver diffusion coefficients) but we need to run
  // the computation on the fluid realm. So, create some transient storage for that.
  Vector<EBAMRCellData> fluidScratchDiffusion(numPlasmaSpecies);
  Vector<EBAMRCellData> fluidScratchDensities(numPlasmaSpecies);
  for (int i = 0; i < numPlasmaSpecies; i++) {
    m_amr->allocate(fluidScratchDiffusion[i], m_fluidRealm, m_plasmaPhase, 1);
    m_amr->allocate(fluidScratchDensities[i], m_fluidRealm, m_plasmaPhase, 1);

    // Copy particle realm data over to fluid realm
    fluidScratchDensities[i].copy(*(a_densities[i]));

    CH_assert(a_diffusionCoefficients[i].getRealm() == m_particleRealm);
    CH_assert(a_densities[i].getRealm() == m_particleRealm);
  }

  // Compute mesh-based diffusion coefficients on the fluid realm.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    Vector<LevelData<EBCellFAB>*> diffusionCoefficients(numPlasmaSpecies);
    Vector<LevelData<EBCellFAB>*> densities(numPlasmaSpecies);

    for (int i = 0; i < numPlasmaSpecies; i++) {
      diffusionCoefficients[i] = &(*(fluidScratchDiffusion[i])[lvl]);
      densities[i]             = &(*(fluidScratchDensities[i])[lvl]);
    }

    this->computeItoDiffusionLFA(diffusionCoefficients, densities, *a_electricField[lvl], lvl, a_time);
  }

  // Copy the fluid realm data over to the particle realm data and then coarsen and update ghost cells.
  for (int i = 0; i < numPlasmaSpecies; i++) {
    a_diffusionCoefficients[i]->copy(fluidScratchDiffusion[i]);

    m_amr->conservativeAverage(*a_diffusionCoefficients[i], m_particleRealm, m_plasmaPhase);
    m_amr->interpGhost(*a_diffusionCoefficients[i], m_particleRealm, m_plasmaPhase);
  }

  // Interpolate diffusion coefficients to particle positions.
  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>& solver = solverIt();

    if (solver->isDiffusive()) {
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

  const int numPlasmaSpecies = m_physics->getNumItoSpecies();

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_fluidRealm)[a_level];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
    const Box box = dbl.get(dit());

    Vector<EBCellFAB*> diffusion(numPlasmaSpecies);
    Vector<EBCellFAB*> densities(numPlasmaSpecies);
    ;

    for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
      const int idx = solverIt.index();

      if (solverIt()->isDiffusive()) {
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

  const int            comp   = 0;
  const Real           dx     = m_amr->getDx()[a_level];
  const RealVect       probLo = m_amr->getProbLo();
  const BaseFab<Real>& E      = a_E.getSingleValuedFAB();

  // Do regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit) {
    const IntVect  iv  = bit();
    const RealVect pos = m_amr->getProbLo() + dx * (RealVect(iv) + 0.5 * RealVect::Unit);
    const RealVect e   = RealVect(D_DECL(E(iv, 0), E(iv, 1), E(iv, 2)));

    // Make grid densities
    Vector<Real> densities;
    for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
      const int idx = solverIt.index();
      densities.push_back((*a_densities[idx]).getSingleValuedFAB()(iv, comp));
    }

    // Call ito_physics and compute diffusion for each particle species
    const Vector<Real> diffusion = m_physics->computeItoDiffusionLFA(a_time, pos, e, densities);

    // Put diffusion where they belong
    for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
      RefCountedPtr<ItoSolver>& solver = solverIt();
      if (solver->isDiffusive()) {
        const int idx                                                 = solverIt.index();
        (*a_diffusionCoefficient[idx]).getSingleValuedFAB()(iv, comp) = diffusion[idx];
      }
    }
  }

  // Do irregular cells
  VoFIterator& vofit = (*m_amr->getVofIterator(m_fluidRealm, m_plasmaPhase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit) {
    const VolIndex& vof = vofit();
    const RealVect  e   = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));
    const RealVect  pos = EBArith::getVofLocation(vof, dx * RealVect::Unit, probLo);

    // Get densities
    Vector<Real> densities;
    for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
      const int idx = solverIt.index();
      densities.push_back((*a_densities[idx])(vof, comp));
    }

    // Compute diffusion
    const Vector<Real> diffusion = m_physics->computeItoDiffusionLFA(a_time, pos, e, densities);

    // Put diffusion in the appropriate place.
    for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
      if (solverIt()->isDiffusive()) {
        const int idx                             = solverIt.index();
        (*a_diffusionCoefficient[idx])(vof, comp) = diffusion[idx];
      }
    }
  }

  // Covered is bogus.
  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    if (solverIt()->isDiffusive()) {
      const int idx = solverIt.index();
      a_diffusionCoefficient[idx]->setCoveredCellVal(0.0, comp);
    }
  }
}

void
ItoPlasmaStepper::computeItoDiffusionLEA()
{
  CH_TIME("ItoPlasmaStepper::computeItoDiffusionLEA()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computeItoDiffusionLEA()" << endl;
  }

  // This is really simple because the solvers do this directly... No monkeying with interpolations or anything.
  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->updateDiffusion();
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
  const EBISLayout&        ebisl = m_amr->getEBISLayout(m_particleRealm, m_plasmaPhase)[a_level];

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

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>& solver = solverIt();
    const int                 idx    = solverIt.index();

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
    VoFIterator& vofit = (*m_amr->getVofIterator(m_particleRealm, m_plasmaPhase)[a_level])[a_dit];
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
  const EBISLayout&        ebisl = m_amr->getEBISLayout(m_particleRealm, m_plasmaPhase)[a_level];

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

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>& solver = solverIt();
    const int                 idx    = solverIt.index();

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
    VoFIterator& vofit = (*m_amr->getVofIterator(m_particleRealm, m_plasmaPhase)[a_level])[a_dit];
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

  this->advanceReactionNetworkNWO(m_electricFieldFluid, m_EdotJ, a_dt);
}

void
ItoPlasmaStepper::advanceReactionNetworkNWO(const EBAMRCellData& a_E, const EBAMRCellData& a_EdotJ, const Real a_dt)
{
  CH_TIME("ItoPlasmaStepper::advanceReactionNetwork(ppc, ypc, E, sources, dt)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::advanceReactionNetwork(ppc, ypc, E, sources, dt)" << endl;
  }

  // 1. Compute the number of particles per cell. Set the number of Photons to be generated per cell to zero.
  this->computeReactiveParticlesPerCell(m_particlePPC);
  this->computeReactiveMeanEnergiesPerCell(m_particleEPS);

  m_fluidPPC.copy(m_particlePPC);
  m_fluidEPS.copy(m_particleEPS);

  DataOps::setValue(m_fluidYPC, 0.0);
  DataOps::setValue(m_particleYPC, 0.0);
  DataOps::copy(m_particleOldPPC, m_particlePPC);

  // 2. Solve for the new number of particles per cell. This also obtains the number of Photons to be generated in each cell.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    this->advanceReactionNetworkNWO(*m_fluidPPC[lvl],
                                    *m_fluidYPC[lvl],
                                    *m_fluidEPS[lvl],
                                    *a_E[lvl],
                                    *a_EdotJ[lvl],
                                    lvl,
                                    a_dt);
  }

  // 3. Copy the results to the particle Realm.
  m_particlePPC.copy(m_fluidPPC);
  m_particleYPC.copy(m_fluidYPC);
  m_particleEPS.copy(m_fluidEPS);

  // 4. Reconcile particles on the particle Realm. Not implemented (yet).
  this->reconcileParticles(m_particlePPC, m_particleOldPPC, m_particleEPS, m_particleYPC);
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

  const int numPlasmaSpecies = m_physics->getNumItoSpecies();
  const int numPhotonSpecies = m_physics->getNumRtSpecies();

  const RealVect probLo = m_amr->getProbLo();

  const EBISBox& ebisbox = m_amr->getEBISLayout(m_fluidRealm, m_plasmaPhase)[a_level][a_dit];
  const EBISBox& ebgraph = m_amr->getEBISLayout(m_fluidRealm, m_plasmaPhase)[a_level][a_dit];

  const BaseFab<Real>& Efab = a_E.getSingleValuedFAB();

  Vector<long long> particles(numPlasmaSpecies);
  Vector<long long> newPhotons(numPhotonSpecies);
  Vector<Real>      meanEnergies(numPlasmaSpecies);
  Vector<Real>      energySources(numPlasmaSpecies);

  // Regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit) {
    const IntVect iv = bit();

    if (ebisbox.isRegular(iv)) {
      const Real kappa = 1.0;
      const Real dV    = kappa * pow(a_dx, SpaceDim);

      const RealVect pos = probLo + a_dx * (RealVect(iv) + 0.5 * RealVect::Unit);
      const RealVect E   = RealVect(D_DECL(Efab(iv, 0), Efab(iv, 1), Efab(iv, 2)));

      // Initialize for this cell.
      for (int i = 0; i < numPlasmaSpecies; i++) {
        particles[i]     = llround(a_particlesPerCell.getSingleValuedFAB()(iv, i));
        meanEnergies[i]  = a_meanParticleEnergies.getSingleValuedFAB()(iv, i);
        energySources[i] = a_EdotJ.getSingleValuedFAB()(iv, i) * dV / Units::Qe;
      }

      for (int i = 0; i < numPhotonSpecies; i++) {
        newPhotons[i] = 0LL;
      }

      // Do the physics advance
      m_physics->advanceParticles(particles, newPhotons, meanEnergies, energySources, a_dt, E, a_dx, kappa);

      // Set result
      for (int i = 0; i < numPlasmaSpecies; i++) {
        a_particlesPerCell.getSingleValuedFAB()(iv, i)     = 1.0 * particles[i];
        a_meanParticleEnergies.getSingleValuedFAB()(iv, i) = 1.0 * meanEnergies[i];
      }

      for (int i = 0; i < numPhotonSpecies; i++) {
        a_newPhotonsPerCell.getSingleValuedFAB()(iv, i) = 1.0 * newPhotons[i];
      }
    }
  }

  // Irregular cells
  VoFIterator& vofit = (*m_amr->getVofIterator(m_fluidRealm, m_plasmaPhase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit) {
    const VolIndex& vof   = vofit();
    const Real      kappa = ebisbox.volFrac(vof);
    const Real      dV    = kappa * pow(a_dx, SpaceDim);
    const RealVect  pos   = EBArith::getVofLocation(vof, a_dx * RealVect::Unit, probLo);
    const RealVect  E     = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));

    // Initialize for this cell.
    for (int i = 0; i < numPlasmaSpecies; i++) {
      particles[i]     = a_particlesPerCell(vof, i);
      meanEnergies[i]  = a_meanParticleEnergies(vof, i);
      energySources[i] = a_EdotJ(vof, i) * dV / Units::Qe;
    }

    for (int i = 0; i < numPhotonSpecies; i++) {
      newPhotons[i] = 0LL;
    }

    m_physics->advanceParticles(particles, newPhotons, meanEnergies, energySources, a_dt, E, a_dx, kappa);

    // Set result
    for (int i = 0; i < numPlasmaSpecies; i++) {
      a_particlesPerCell(vof, i)     = 1.0 * particles[i];
      a_meanParticleEnergies(vof, i) = 1.0 * meanEnergies[i];
    }

    for (int i = 0; i < numPhotonSpecies; i++) {
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

  const int numPlasmaSpecies = m_physics->getNumItoSpecies();
  const int numPhotonSpecies = m_physics->getNumRtSpecies();

  const RealVect probLo = m_amr->getProbLo();

  const EBISBox& ebisbox = m_amr->getEBISLayout(m_particleRealm, m_plasmaPhase)[a_level][a_dit];
  const EBISBox& ebgraph = m_amr->getEBISLayout(m_particleRealm, m_plasmaPhase)[a_level][a_dit];

  Vector<BinFab<ItoParticle>*> particlesFAB(numPlasmaSpecies);
  Vector<BinFab<Photon>*>      sourcePhotonsFAB(numPhotonSpecies);
  Vector<BinFab<Photon>*>      bulkPhotonsFAB(numPhotonSpecies);

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>& solver = solverIt();
    const int                 idx    = solverIt.index();

    ParticleContainer<ItoParticle>& solverParticles = solver->getParticles(ItoSolver::WhichContainer::Bulk);

    particlesFAB[idx] = &(solverParticles.getCellParticles(a_level, a_dit));
  }

  for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<McPhoto>& solver = solverIt();
    const int               idx    = solverIt.index();

    ParticleContainer<Photon>& solverBulkPhotons = solver->getBulkPhotons();
    ParticleContainer<Photon>& solverSourPhotons = solver->getSourcePhotons();

    bulkPhotonsFAB[idx]   = &(solverBulkPhotons.getCellParticles(a_level, a_dit));
    sourcePhotonsFAB[idx] = &(solverSourPhotons.getCellParticles(a_level, a_dit));
  }

  // Regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit) {
    const IntVect iv = bit();
    if (ebisbox.isRegular(iv)) {
      const RealVect cellPos       = probLo + a_dx * (RealVect(iv) + 0.5 * RealVect::Unit);
      const RealVect centroidPos   = cellPos;
      const RealVect lo            = -0.5 * RealVect::Unit;
      const RealVect hi            = 0.5 * RealVect::Unit;
      const RealVect bndryCentroid = RealVect::Zero;
      const RealVect bndryNormal   = RealVect::Zero;
      const Real     kappa         = 1.0;

      Vector<List<ItoParticle>*>       particles(numPlasmaSpecies);
      Vector<List<Photon>*>            bulkPhotons(numPhotonSpecies);
      Vector<List<Photon>*>            sourcePhotons(numPhotonSpecies);
      Vector<RefCountedPtr<RtSpecies>> photoSpecies(numPhotonSpecies);

      Vector<Real>      particleMeanEnergies(numPlasmaSpecies);
      Vector<long long> numNewParticles(numPlasmaSpecies);
      Vector<long long> numOldParticles(numPlasmaSpecies);
      Vector<long long> numNewPhotons(numPhotonSpecies);

      for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
        const int idx = solverIt.index();

        particles[idx]            = &((*particlesFAB[idx])(iv, 0));
        particleMeanEnergies[idx] = a_meanParticleEnergies.getSingleValuedFAB()(iv, idx);
        numNewParticles[idx]      = llround(a_newParticlesPerCell.getSingleValuedFAB()(iv, idx));
        numOldParticles[idx]      = llround(a_oldParticlesPerCell.getSingleValuedFAB()(iv, idx));
      }

      for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
        const int idx = solverIt.index();

        bulkPhotons[idx]   = &((*bulkPhotonsFAB[idx])(iv, 0));
        sourcePhotons[idx] = &((*sourcePhotonsFAB[idx])(iv, 0));
        photoSpecies[idx]  = solverIt()->getSpecies();
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
      for (int i = 0; i < numPhotonSpecies; i++) {
        //	bulkPhotons[i]->clear();
      }
    }
  }

  // Irregular cells
  VoFIterator& vofit = (*m_amr->getVofIterator(m_particleRealm, m_plasmaPhase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit) {
    const VolIndex& vof           = vofit();
    const IntVect   iv            = vof.gridIndex();
    const Real      kappa         = ebisbox.volFrac(vof);
    const RealVect  cellPos       = EBArith::getVofLocation(vof, a_dx * RealVect::Unit, probLo);
    const RealVect  centroidPos   = ebisbox.centroid(vof);
    const RealVect  bndryCentroid = ebisbox.bndryCentroid(vof);
    const RealVect  bndryNormal   = ebisbox.normal(vof);

    RealVect lo = -0.5 * RealVect::Unit;
    RealVect hi = 0.5 * RealVect::Unit;
    if (kappa < 1.0) {
      DataOps::computeMinValidBox(lo, hi, bndryNormal, bndryCentroid);
    }

    Vector<List<ItoParticle>*>       particles(numPlasmaSpecies);
    Vector<List<Photon>*>            bulkPhotons(numPhotonSpecies);
    Vector<List<Photon>*>            sourcePhotons(numPhotonSpecies);
    Vector<RefCountedPtr<RtSpecies>> photoSpecies(numPhotonSpecies);

    Vector<Real>      particleMeanEnergies(numPlasmaSpecies);
    Vector<long long> numNewParticles(numPlasmaSpecies);
    Vector<long long> numOldParticles(numPlasmaSpecies);
    Vector<long long> numNewPhotons(numPhotonSpecies);

    for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
      const int idx = solverIt.index();

      particles[idx]            = &((*particlesFAB[idx])(iv, 0));
      particleMeanEnergies[idx] = a_meanParticleEnergies(vof, idx);
      numNewParticles[idx]      = llround(a_newParticlesPerCell(vof, idx));
      numOldParticles[idx]      = llround(a_oldParticlesPerCell(vof, idx));
    }

    for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
      const int idx = solverIt.index();

      bulkPhotons[idx]   = &((*bulkPhotonsFAB[idx])(iv, 0));
      sourcePhotons[idx] = &((*sourcePhotonsFAB[idx])(iv, 0));
      photoSpecies[idx]  = solverIt()->getSpecies();
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
    for (int i = 0; i < numPhotonSpecies; i++) {
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
    const int numPlasmaSpecies = m_physics->getNumItoSpecies();
    const int numPhotonSpecies = m_physics->getNumRtSpecies();

    Vector<ParticleContainer<ItoParticle>*> particles(numPlasmaSpecies);    // Current particles.
    Vector<ParticleContainer<Photon>*>      bulk_Photons(numPhotonSpecies); // Photons absorbed on mesh
    Vector<ParticleContainer<Photon>*>      new_Photons(numPhotonSpecies);  // Produced Photons go here.

    for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
      particles[solverIt.index()] = &(solverIt()->getParticles(ItoSolver::WhichContainer::Bulk));
    }

    for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
      bulk_Photons[solverIt.index()] = &(solverIt()->getBulkPhotons());
      new_Photons[solverIt.index()]  = &(solverIt()->getSourcePhotons());
    }

    this->advanceReactionNetwork(particles, bulk_Photons, new_Photons, m_energySources, m_electricFieldParticle, a_dt);
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

  const int numPlasmaSpecies = m_physics->getNumItoSpecies();
  const int numPhotonSpecies = m_physics->getNumRtSpecies();

  Vector<AMRCellParticles<ItoParticle>*> particles(numPlasmaSpecies);
  Vector<AMRCellParticles<Photon>*>      Photons(numPlasmaSpecies);
  Vector<AMRCellParticles<Photon>*>      newPhotons(numPlasmaSpecies);

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    const int idx  = solverIt.index();
    particles[idx] = &(a_particles[idx]->getCellParticles());
  }

  for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    const int idx   = solverIt.index();
    Photons[idx]    = &(a_Photons[idx]->getCellParticles());
    newPhotons[idx] = &(a_newPhotons[idx]->getCellParticles());
  }

  //Advance reaction network
  this->advanceReactionNetwork(particles, Photons, newPhotons, a_sources, m_electricFieldParticle, a_dt);
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

  const int numPlasmaSpecies = m_physics->getNumItoSpecies();
  const int numPhotonSpecies = m_physics->getNumRtSpecies();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    Vector<LayoutData<BinFab<ItoParticle>>*> particles(numPlasmaSpecies);
    Vector<LayoutData<BinFab<Photon>>*>      Photons(numPhotonSpecies);
    Vector<LayoutData<BinFab<Photon>>*>      newPhotons(numPhotonSpecies);
    Vector<LevelData<EBCellFAB>*>            sources(numPlasmaSpecies);

    for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
      const int idx  = solverIt.index();
      particles[idx] = &(*(*a_particles[idx])[lvl]);
      sources[idx]   = &(*(a_sources[idx])[lvl]);
    }

    for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
      const int idx   = solverIt.index();
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

  const int numPlasmaSpecies = m_physics->getNumItoSpecies();
  const int numPhotonSpecies = m_physics->getNumRtSpecies();

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[a_lvl];
  const Real               dx  = m_amr->getDx()[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
    const Box box = dbl.get(dit());

    Vector<BinFab<ItoParticle>*> particles(numPlasmaSpecies);
    Vector<BinFab<Photon>*>      Photons(numPhotonSpecies);
    ;
    Vector<BinFab<Photon>*> newPhotons(numPhotonSpecies);
    Vector<EBCellFAB*>      sources(numPlasmaSpecies);

    for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
      const int idx  = solverIt.index();
      particles[idx] = &((*a_particles[idx])[dit()]);
      sources[idx]   = &((*a_sources[idx])[dit()]);
    }

    for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
      const int idx   = solverIt.index();
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

  const int numPlasmaSpecies = m_physics->getNumItoSpecies();
  const int numPhotonSpecies = m_physics->getNumRtSpecies();

  const RealVect probLo = m_amr->getProbLo();
  const RealVect dx     = a_dx * RealVect::Unit;

  const EBISBox& ebisbox = m_amr->getEBISLayout(m_particleRealm, m_plasmaPhase)[a_lvl][a_dit];
  const EBISBox& ebgraph = m_amr->getEBISLayout(m_particleRealm, m_plasmaPhase)[a_lvl][a_dit];

  const BaseFab<Real>& Efab = a_E.getSingleValuedFAB();

  // Regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit) {
    const IntVect iv = bit();

    if (ebisbox.isRegular(iv)) {
      const Real     kappa = 1.0;
      const RealVect pos   = probLo + a_dx * (RealVect(iv) + 0.5 * RealVect::Unit);
      const RealVect e     = RealVect(D_DECL(Efab(iv, 0), Efab(iv, 1), Efab(iv, 2)));

      Vector<List<ItoParticle>*> particles(numPlasmaSpecies);
      Vector<List<Photon>*>      Photons(numPhotonSpecies);
      Vector<List<Photon>*>      newPhotons(numPhotonSpecies);
      Vector<Real>               sources(numPlasmaSpecies);

      for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
        const int idx = solverIt.index();

        List<ItoParticle>& bp = (*a_particles[idx])(iv, comp);
        particles[idx]        = &bp;

        const BaseFab<Real>& sourcesFAB = a_sources[idx]->getSingleValuedFAB();
        sources[idx]                    = sourcesFAB(iv, comp);
      }

      for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
        const int idx = solverIt.index();

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
  VoFIterator& vofit = (*m_amr->getVofIterator(m_particleRealm, m_plasmaPhase)[a_lvl])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit) {
    const VolIndex vof   = vofit();
    const IntVect  iv    = vof.gridIndex();
    const RealVect pos   = probLo + a_dx * (RealVect(iv) + 0.5 * RealVect::Unit);
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

    Vector<List<ItoParticle>*> particles(numPlasmaSpecies);
    Vector<List<Photon>*>      Photons(numPhotonSpecies);
    Vector<List<Photon>*>      newPhotons(numPhotonSpecies);
    Vector<Real>               sources(numPlasmaSpecies);

    for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
      const int idx = solverIt.index();

      List<ItoParticle>& bp = (*a_particles[idx])(iv, comp);
      particles[idx]        = &bp;

      const BaseFab<Real>& sourcesFAB = a_sources[idx]->getSingleValuedFAB();
      sources[idx]                    = sourcesFAB(iv, comp);
    }

    for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
      const int idx = solverIt.index();

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

  const Real dt = this->computePhysicsDt(m_electricFieldParticle, m_ito->getDensities());

  return dt;
}

Real
ItoPlasmaStepper::computePhysicsDt(const EBAMRCellData& a_E, const Vector<EBAMRCellData*> a_densities) const
{
  CH_TIME("ItoPlasmaStepper::computePhysicsDt(EBAMRCellFAB, Vector<EBAMRCellFAB*>)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::computePhysicsDt(EBAMRCellFAB, Vector<EBAMRCellFAB*>)" << endl;
  }

  const int numPlasmaSpecies = m_physics->getNumItoSpecies();

  Real minDt = 1.E99;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {

    Vector<LevelData<EBCellFAB>*> densities(numPlasmaSpecies);

    for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
      const int idx = solverIt.index();

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

  const int numPlasmaSpecies = m_physics->getNumItoSpecies();

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[a_level];

  Real minDt = 1.E99;

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

    Vector<EBCellFAB*> densities(numPlasmaSpecies);

    for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
      const int idx = solverIt.index();

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

  const int numPlasmaSpecies = m_physics->getNumItoSpecies();

  const int            comp    = 0;
  const Real           dx      = m_amr->getDx()[a_level];
  const RealVect       probLo  = m_amr->getProbLo();
  const BaseFab<Real>& E       = a_E.getSingleValuedFAB();
  const EBISBox&       ebisbox = m_amr->getEBISLayout(m_particleRealm, m_plasmaPhase)[a_level][a_dit];

  // Do regular cells
  for (BoxIterator bit(a_box); bit.ok(); ++bit) {
    const IntVect  iv  = bit();
    const RealVect pos = m_amr->getProbLo() + dx * (RealVect(iv) + 0.5 * RealVect::Unit);
    const RealVect e   = RealVect(D_DECL(E(iv, 0), E(iv, 1), E(iv, 2)));

    if (ebisbox.isRegular(iv)) {
      Vector<Real> densities(numPlasmaSpecies);
      for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
        const int idx = solverIt.index();

        const BaseFab<Real>& basefab = a_densities[idx]->getSingleValuedFAB();
        densities[idx]               = basefab(iv, comp);
      }

      const Real cellDt = m_physics->computeDt(e, pos, densities);

      minDt = Min(minDt, cellDt);
    }
  }

  // Do irregular cells
  VoFIterator& vofit = (*m_amr->getVofIterator(m_particleRealm, m_plasmaPhase)[a_level])[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit) {
    const VolIndex& vof = vofit();
    const RealVect  e   = RealVect(D_DECL(a_E(vof, 0), a_E(vof, 1), a_E(vof, 2)));
    const RealVect  pos = EBArith::getVofLocation(vof, dx * RealVect::Unit, probLo);

    Vector<Real> densities(numPlasmaSpecies);

    for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
      const int idx  = solverIt.index();
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

  for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<McPhoto>& solver = solverIt();

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

  for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->sortPhotonsByCell(a_which);
  }
}

void
ItoPlasmaStepper::sortPhotonsByPatch(McPhoto::WhichContainer a_which)
{
  CH_TIME("ItoPlasmaStepper::sortPhotonsByPatch(McPhoto::WhichContainer)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::sortPhotonsByPatch(McPhoto::WhichContainer)" << endl;
  }

  for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->sortPhotonsByPatch(a_which);
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

    // If we're load balancing with particles, get the number of particles per patch
    // from the relevant particle solvers. Since these are Ito solvers, the loads
    // are equal to the number of computational particles in the grid patches.
    Vector<RefCountedPtr<ItoSolver>> loadBalanceProxySolvers = this->getLoadBalanceSolvers();

    for (int isolver = 0; isolver < loadBalanceProxySolvers.size(); isolver++) {

      // This solver computes loads -- there's a parallel gather operation
      // under the hood here.
      Vector<long int> solverLoads(nbox, 0L);
      loadBalanceProxySolvers[isolver]->computeLoads(solverLoads, dbl, a_level);

      // Add to total loads.
      for (int ibox = 0; ibox < nbox; ibox++) {
        loads[ibox] += solverLoads[ibox];
      }
    }

    // Add the "constant" loads -- these are computational loads due to the "mesh" part. We use
    // a heuristic where we have m_loadPerCell "cost".
    for (LayoutIterator lit = dbl.layoutIterator(); lit.ok(); ++lit) {
      const Box box = dbl[lit()];

      loads[lit().intCode()] += lround(m_loadPerCell * box.numPts());
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
                                           const int                        a_finestLevel) noexcept
{
  CH_TIME("ItoPlasmaStepper::loadBalanceParticleRealm(...)");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::loadBalanceParticleRealm(...)" << endl;
  }

  // Decompose the DisjointBoxLayout
  a_procs.resize(1 + a_finestLevel);
  a_boxes.resize(1 + a_finestLevel);

  for (int lvl = a_lmin; lvl <= a_finestLevel; lvl++) {
    a_procs[lvl] = a_grids[lvl].procIDs();
    a_boxes[lvl] = a_grids[lvl].boxArray();
  }

  // Get the particles that we will use for load balancing.
  Vector<RefCountedPtr<ItoSolver>> loadBalanceProxySolvers = this->getLoadBalanceSolvers();

  // Regrid particles onto the "dummy grids" a_grids
  for (int i = 0; i < loadBalanceProxySolvers.size(); i++) {
    ParticleContainer<ItoParticle>& particles = loadBalanceProxySolvers[i]->getParticles(
      ItoSolver::WhichContainer::Bulk);

    m_amr->remapToNewGrids(particles, a_lmin, a_finestLevel);

    // If we make superparticles during regrids, do it here so we can better estimate the computational loads for each patch. This way, if a grid is removed the realistic
    // load estimate of the underlying grid(s) is improved.
    if (m_regridSuperparticles) {
      particles.sortParticlesByCell();
      loadBalanceProxySolvers[i]->makeSuperparticles(ItoSolver::WhichContainer::Bulk, m_particlesPerCell);
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
  for (int i = 0; i < loadBalanceProxySolvers.size(); i++) {
    ParticleContainer<ItoParticle>& particles = loadBalanceProxySolvers[i]->getParticles(
      ItoSolver::WhichContainer::Bulk);
    particles.preRegrid(a_lmin);
  }
}

Vector<RefCountedPtr<ItoSolver>>
ItoPlasmaStepper::getLoadBalanceSolvers() const noexcept
{
  CH_TIME("ItoPlasmaStepper::getLoadBalanceSolvers()");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::getLoadBalanceSolvers()" << endl;
  }

  Vector<RefCountedPtr<ItoSolver>> loadBalanceProxySolvers;

  if (m_loadBalanceIndex < 0) {
    for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
      loadBalanceProxySolvers.push_back(solverIt());
    }
  }
  else {
    RefCountedPtr<ItoSolver>& solver = m_ito->getSolvers()[m_loadBalanceIndex];
    loadBalanceProxySolvers.push_back(solver);
  }

  return loadBalanceProxySolvers;
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

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();
    const int q   = species->getChargeNumber();

    DataOps::setValue(m_energySources[idx], 0.0);

    // Do mobile contribution.
    if (q != 0 && solver->isMobile()) {

      // Drift contribution
      solver->depositConductivity(m_particleScratch1,
                                  solver->getParticles(ItoSolver::WhichContainer::Bulk)); // Deposit mu*n
      DataOps::copy(
        m_particleScratchD,
        m_electricFieldParticle); // Could use m_electricFieldParticle or solver's m_velo_func here, but m_velo_func = +/- E (depends on q)

      DataOps::multiplyScalar(m_particleScratchD, m_particleScratch1); // m_particleScratchD = mu*n*E
      DataOps::dotProduct(m_particleScratch1,
                          m_electricFieldParticle,
                          m_particleScratchD);                      // m_particleScratch1 = mu*n*E*E
      DataOps::incr(m_energySources[idx], m_particleScratch1, 1.0); // a_source[idx] += mu*n*E*E
    }

    // Diffusive contribution
    if (q != 0 && solver->isDiffusive()) {

      // Compute the negative gradient of the diffusion term
      solver->depositDiffusivity(m_particleScratch1, solver->getParticles(ItoSolver::WhichContainer::Bulk));
      m_amr->interpGhostMG(m_particleScratch1, m_particleRealm, m_plasmaPhase);
      m_amr->computeGradient(m_particleScratchD, m_particleScratch1, m_particleRealm, m_plasmaPhase);
      DataOps::scale(m_particleScratchD, -1.0); // scratchD = -grad(D*n)

      DataOps::dotProduct(m_particleScratch1,
                          m_particleScratchD,
                          m_electricFieldParticle);                 // m_particleScratch1 = -E*grad(D*n)
      DataOps::incr(m_energySources[idx], m_particleScratch1, 1.0); // a_source[idx]
    }

    if (q != 0 && (solver->isMobile() || solver->isDiffusive())) {
      DataOps::scale(m_energySources[idx], Abs(q) * Units::Qe);
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

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();
    const int q   = species->getChargeNumber();

    // Do mobile contribution. Computes Z*e*E*mu*n*E*E
    if (q != 0 && solver->isMobile()) {
      solver->depositConductivity(m_particleScratch1,
                                  solver->getParticles(ItoSolver::WhichContainer::Bulk)); // Deposit mu*n
      m_fluidScratch1.copy(m_particleScratch1);                                           // Copy mu*n to fluid Realm
      DataOps::copy(m_fluidScratchD, m_electricFieldFluid);                               // m_fluidScratchD = E
      DataOps::multiplyScalar(m_fluidScratchD, m_fluidScratch1);                          // m_fluidScratchD = E*mu*n
      DataOps::dotProduct(m_fluidScratch1,
                          m_electricFieldFluid,
                          m_fluidScratchD);                // m_particleScratch1 = E.dot.(E*mu*n)
      DataOps::scale(m_fluidScratch1, Abs(q) * Units::Qe); // m_particleScratch1 = Z*e*mu*n*E*E

      m_amr->conservativeAverage(m_fluidScratch1, m_fluidRealm, m_plasmaPhase);
      m_amr->interpGhost(m_fluidScratch1, m_fluidRealm, m_plasmaPhase);
      DataOps::plus(m_EdotJ, m_fluidScratch1, 0, idx, 1); // a_source[idx] += Z*e*mu*n*E*E
    }

    // Diffusive contribution. Computes -Z*e*E*grad(D*n)
    if (q != 0 && solver->isDiffusive()) {
      solver->depositDiffusivity(m_particleScratch1,
                                 solver->getParticles(ItoSolver::WhichContainer::Bulk)); // Deposit D*n
      m_fluidScratch1.copy(m_particleScratch1);                                          // Copy D*n to fluid Realm
      m_amr->interpGhostMG(m_fluidScratch1, m_fluidRealm, m_plasmaPhase);
      m_amr->computeGradient(m_fluidScratchD, m_fluidScratch1, m_fluidRealm, m_plasmaPhase); // scratchD = grad(D*n)
      DataOps::scale(m_fluidScratchD, -1.0);                                                 // scratchD = -grad(D*n)
      DataOps::dotProduct(m_fluidScratch1, m_fluidScratchD, m_electricFieldFluid); // scratch1 = -E.dot.grad(D*n)
      DataOps::scale(m_fluidScratch1, Abs(q) * Units::Qe);                         // scratch1 = -Z*e*E*grad(D*n)

      m_amr->conservativeAverage(m_fluidScratch1, m_fluidRealm, m_plasmaPhase);
      m_amr->interpGhost(m_fluidScratch1, m_fluidRealm, m_plasmaPhase);

      DataOps::plus(m_EdotJ, m_fluidScratch1, 0, idx, 1); // source  += -Z*e*E*grad(D*n)
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

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();
    const int q   = species->getChargeNumber();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    ParticleContainer<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk);

    const DepositionType deposition = solver->getDeposition();

    if ((mobile || diffusive) && q != 0) {

      // We will interpolate m_electricFieldParticle onto particle velocity vectors.
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[lvl];

        for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
          const EBCellFAB& E       = (*m_electricFieldParticle[lvl])[dit()];
          const EBISBox&   ebisbox = E.getEBISBox();
          const FArrayBox& Efab    = E.getFArrayBox();
          const RealVect   dx      = m_amr->getDx()[lvl] * RealVect::Unit;
          const RealVect   probLo  = m_amr->getProbLo();
          const Box        box     = dbl[dit()];

          List<ItoParticle>& particleList = particles[lvl][dit()].listItems();

          // This interpolates the velocity function on to the particle velocities
#if 1
          MayDay::Warning("EBParticleMesh should be replaced with call to AmrMesh as in ItoSolver");
#endif
          EBParticleMesh meshInterp(box, ebisbox, dx, probLo);
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
      solver->depositParticles<ItoParticle, &ItoParticle::weight>(m_particleScratch1, particles);
      m_fluidScratch1.copy(m_particleScratch1);

      // Scale by Qe/dt to make it Joule/dt. Then add to correct index
      DataOps::scale(m_fluidScratch1, q * Units::Qe / a_dt);
      DataOps::plus(m_EdotJ, m_fluidScratch1, 0, idx, 1);

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
