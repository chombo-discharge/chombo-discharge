/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoKMCStepperImplem.H
  @brief  Implementation of CD_ItoKMCStepper.H
  @author Robert Marskar
*/

#ifndef CD_ItoKMCStepperImplem_H
#define CD_ItoKMCStepperImplem_H

// Std includes
#include <limits>

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_EBFineInterp.H>
#include <CD_ItoKMCStepper.H>
#include <CD_ParticleOps.H>
#include <CD_DataOps.H>
#include <CD_ParallelOps.H>
#include <CD_Units.H>
#include <CD_Timer.H>
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

template <typename I, typename C, typename R, typename F>
ItoKMCStepper<I, C, R, F>::ItoKMCStepper() noexcept
{
  CH_TIME("ItoKMCStepper::ItoKMCStepper");

  m_verbosity                        = -1;
  m_profile                          = false;
  m_name                             = "ItoKMCStepper";
  m_plasmaPhase                      = phase::gas;
  m_dt                               = 0.0;
  m_time                             = 0.0;
  m_timeStep                         = 0;
  m_loadPerCell                      = 1.0;
  m_redistributeCDR                  = true;
  m_regridSuperparticles             = true;
  m_fluidRealm                       = Realm::Primal;
  m_particleRealm                    = Realm::Primal;
  m_minParticleAdvectionCFL          = 0.0;
  m_maxParticleAdvectionCFL          = 1.0;
  m_minParticleDiffusionCFL          = 0.0;
  m_maxParticleDiffusionCFL          = std::numeric_limits<Real>::max();
  m_minParticleAdvectionDiffusionCFL = std::numeric_limits<Real>::max();
  m_maxParticleAdvectionDiffusionCFL = std::numeric_limits<Real>::max();
  m_fluidAdvectionDiffusionCFL       = 0.5;
  m_relaxTimeFactor                  = std::numeric_limits<Real>::max();
  m_minDt                            = std::numeric_limits<Real>::min();
  m_maxDt                            = std::numeric_limits<Real>::max();
}

template <typename I, typename C, typename R, typename F>
ItoKMCStepper<I, C, R, F>::ItoKMCStepper(RefCountedPtr<ItoKMCPhysics>& a_physics) noexcept : ItoKMCStepper<I, C, R, F>()
{
  CH_TIME("ItoKMCStepper::ItoKMCStepper(RefCountrPtr<ItoKMCPhysics>)");

  m_physics = a_physics;

  if (m_physics->getNumPlasmaSpecies() == 0) {
    MayDay::Abort("ItoKMCStepper::ItoKMCStepper -- numPlasmaSpecies = 0, there's no problem to solve here!");
  }
}

template <typename I, typename C, typename R, typename F>
ItoKMCStepper<I, C, R, F>::~ItoKMCStepper() noexcept
{
  CH_TIME("ItoKMCStepper::~ItoKMCStepper");
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::parseOptions() noexcept
{
  CH_TIME("ItoKMCStepper::parseOptions");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseOptions" << endl;
  }

  this->parseVerbosity();
  this->parseExitOnFailure();
  this->parseRedistributeCDR();
  this->parsePlotVariables();
  this->parseSuperParticles();
  this->parseDualGrid();
  this->parseLoadBalance();
  this->parseTimeStepRestrictions();
  this->parseParametersEB();
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::parseRuntimeOptions() noexcept
{
  CH_TIME("ItoKMCStepper::parseRuntimeOptions");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRuntimeOptions" << endl;
  }

  this->parseVerbosity();
  this->parseExitOnFailure();
  this->parseRedistributeCDR();
  this->parsePlotVariables();
  this->parseSuperParticles();
  this->parseLoadBalance();
  this->parseTimeStepRestrictions();
  this->parseParametersEB();

  m_ito->parseRuntimeOptions();
  m_cdr->parseRuntimeOptions();
  m_fieldSolver->parseRuntimeOptions();
  m_rte->parseRuntimeOptions();
  m_sigmaSolver->parseRuntimeOptions();

  m_physics->parseRuntimeOptions();
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::parseVerbosity() noexcept
{
  CH_TIME("ItoKMCStepper::parseVerbosity");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseVerbosity" << endl;
  }

  ParmParse pp(m_name.c_str());

  pp.get("verbosity", m_verbosity);
  pp.get("profile", m_profile);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::parseExitOnFailure() noexcept
{
  CH_TIME("ItoKMCStepper::parseExitOnFailure");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseExitOnFailure" << endl;
  }

  ParmParse pp(m_name.c_str());

  pp.get("abort_on_failure", m_abortOnFailure);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::parseRedistributeCDR() noexcept
{
  CH_TIME("ItoKMCStepper::parseRedistributeCDR");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRedistributeCDR" << endl;
  }

  ParmParse pp(m_name.c_str());

  pp.get("redistribute_cdr", m_redistributeCDR);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::parsePlotVariables() noexcept
{
  CH_TIME("ItoKMCStepper::parsePlotVariables");
  if (m_verbosity > 5) {
    pout() << m_name + "::parsePlotVariables" << endl;
  }

  m_plotConductivity      = false;
  m_plotCurrentDensity    = false;
  m_plotParticlesPerPatch = false;

  // Read in plot variables.
  ParmParse pp(m_name.c_str());
  const int num = pp.countval("plt_vars");

  if (num > 0) {
    Vector<std::string> str(num);
    pp.getarr("plt_vars", str, 0, num);

    // Set plot variables
    for (int i = 0; i < num; i++) {
      if (str[i] == "conductivity") {
        m_plotConductivity = true;
      }
      else if (str[i] == "current_density") {
        m_plotCurrentDensity = true;
      }
      else if (str[i] == "particles_per_patch") {
        m_plotParticlesPerPatch = true;
      }
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::parseSuperParticles() noexcept
{
  CH_TIME("ItoKMCStepper::parseSuperParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseSuperParticles" << endl;
  }

  ParmParse pp(m_name.c_str());

  m_particlesPerCell.resize(pp.countval("particles_per_cell"));

  pp.get("merge_interval", m_mergeInterval);
  pp.get("regrid_superparticles", m_regridSuperparticles);
  pp.getarr("particles_per_cell", m_particlesPerCell, 0, m_particlesPerCell.size());

  for (int lvl = 0; lvl < m_particlesPerCell.size(); lvl++) {
    if (m_particlesPerCell[lvl] <= 0) {
      MayDay::Error("ItoKMCStepper::parseSuperParticles -- must have 'particles_per_cell' > 0");
    }
  }

  m_mergeInterval = std::max(m_mergeInterval, 1);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::parseDualGrid() noexcept
{
  CH_TIME("ItoKMCStepper::parseDualGrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseDualGrid" << endl;
  }

  ParmParse pp(m_name.c_str());

  pp.get("dual_grid", m_dualGrid);

  if (m_dualGrid) {
    m_particleRealm = "ParticleRealm";

    CH_assert(m_particleRealm != m_fluidRealm);
  }
  else {
    m_particleRealm = m_fluidRealm;
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::parseLoadBalance() noexcept
{
  CH_TIME("ItoKMCStepper::parseLoadBalance");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseLoadBalance" << endl;
  }

  ParmParse pp(m_name.c_str());

  std::string str;

  pp.get("load_balance_particles", m_loadBalanceParticles);
  pp.get("load_balance_fluid", m_loadBalanceFluid);
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
    const std::string err = "ItoKMCStepper::parseLoadBalance - 'box_sorting = " + str + "' not recognized";

    MayDay::Error(err.c_str());
  }

  // Get the load balancing index.
  const int numIndices = pp.countval("load_indices");

  if (numIndices > 0) {
    pp.getarr("load_indices", m_loadBalanceIndices, 0, numIndices);
  }
  else {
    const std::string err = "ItoKMCStepper::parseLoadBalance - 'load_indices' argument has zero entries";

    MayDay::Error(err.c_str());
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::parseTimeStepRestrictions() noexcept
{
  CH_TIME("ItoKMCStepper::parseTimeStepRestrictions");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseTimeStepRestrictions" << endl;
  }

  ParmParse pp(m_name.c_str());

  pp.get("min_particle_advection_cfl", m_minParticleAdvectionCFL);
  pp.get("max_particle_advection_cfl", m_maxParticleAdvectionCFL);
  pp.get("min_particle_diffusion_cfl", m_minParticleDiffusionCFL);
  pp.get("max_particle_diffusion_cfl", m_maxParticleDiffusionCFL);
  pp.get("min_particle_advection_diffusion_cfl", m_minParticleAdvectionDiffusionCFL);
  pp.get("max_particle_advection_diffusion_cfl", m_maxParticleAdvectionDiffusionCFL);
  pp.get("fluid_advection_diffusion_cfl", m_fluidAdvectionDiffusionCFL);
  pp.get("relax_dt_factor", m_relaxTimeFactor);
  pp.get("min_dt", m_minDt);
  pp.get("max_dt", m_maxDt);

  if (m_relaxTimeFactor <= 0.0) {
    MayDay::Error("ItoKMCStepper::parseTimeStepRestrictions() - must have relax_dt > 0.0");
  }

  if (m_minDt < 0.0) {
    MayDay::Error("ItoKMCStepper::parseTimeStepRestrictions() - must have min_dt >= 0.0");
  }

  if (m_maxDt < 0.0) {
    MayDay::Error("ItoKMCStepper::parseTimeStepRestrictions() - must have max_dt >= 0.0");
  }

  if (m_maxParticleAdvectionCFL <= 0.0) {
    MayDay::Error("ItoKMCStepper::parseTimeStepRestrictions() - must have max_particle_advection_cfl > 0.0");
  }

  if (m_minParticleAdvectionCFL < 0.0) {
    MayDay::Error("ItoKMCStepper::parseTimeStepRestrictions() - must have min_particle_advection_cfl >= 0.0");
  }

  if (m_maxParticleDiffusionCFL <= 0.0) {
    MayDay::Error("ItoKMCStepper::parseTimeStepRestrictions() - must have particle_diffusion_cfl > 0.0");
  }

  if (m_minParticleDiffusionCFL < 0.0) {
    MayDay::Error("ItoKMCStepper::parseTimeStepRestrictions() - must have particle_diffusion_cfl >= 0.0");
  }

  if (m_maxParticleAdvectionDiffusionCFL <= 0.0) {
    MayDay::Error("ItoKMCStepper::parseTimeStepRestrictions() - must have particle_advection_diffusion_cfl > 0.0");
  }

  if (m_minParticleAdvectionDiffusionCFL < 0.0) {
    MayDay::Error("ItoKMCStepper::parseTimeStepRestrictions() - must have particle_advection_diffusion_cfl >= 0.0");
  }

  if (m_fluidAdvectionDiffusionCFL <= 0.0) {
    MayDay::Error("ItoKMCStepper::parseTimeStepRestrictions() - must have fluid_advection_diffusion_cfl > 0.0");
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::parseParametersEB() noexcept
{
  CH_TIME("ItoKMCStepper::parseTimeStepRestrictions");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseTimeStepRestrictions" << endl;
  }

  ParmParse pp(m_name.c_str());

  std::string str;

  pp.get("eb_tolerance", m_toleranceEB);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::setupSolvers() noexcept
{
  CH_TIME("ItoKMCStepper::setupSolver");
  if (m_verbosity > 5) {
    pout() << m_name + "::setupSolvers" << endl;
  }

  this->setupIto();
  this->setupCdr();
  this->setupPoisson();
  this->setupRadiativeTransfer();
  this->setupSigma();
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::setupIto() noexcept
{
  CH_TIME("ItoKMCStepper::setupIto");
  if (m_verbosity > 5) {
    pout() << m_name + "::setupIto" << endl;
  }

  ItoFactory<ItoSolver, I> factory;
  m_ito = factory.newLayout(m_physics->getItoSpecies());

  m_ito->parseOptions();
  m_ito->setAmr(m_amr);
  m_ito->setPhase(m_plasmaPhase);
  m_ito->setComputationalGeometry(m_computationalGeometry);
  m_ito->setRealm(m_particleRealm);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::setupCdr() noexcept
{
  CH_TIME("ItoKMCStepper::setupCdr");
  if (m_verbosity > 5) {
    pout() << m_name + "::setupCdr" << endl;
  }

  CdrFactory<CdrSolver, C> factory;
  m_cdr = factory.newLayout(m_physics->getCdrSpecies());

  m_cdr->parseOptions();
  m_cdr->setAmr(m_amr);
  m_cdr->setPhase(m_plasmaPhase);
  m_cdr->setComputationalGeometry(m_computationalGeometry);
  m_cdr->setRealm(m_fluidRealm);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::setupRadiativeTransfer() noexcept
{
  CH_TIME("ItoKMCStepper::setupRadiativeTransfer");
  if (m_verbosity > 5) {
    pout() << m_name + "::setupRadiativeTransfer" << endl;
  }

  RtFactory<McPhoto, R> factory;
  m_rte = factory.newLayout(m_physics->getRtSpecies());

  m_rte->parseOptions();
  m_rte->setPhase(m_plasmaPhase);
  m_rte->setAmr(m_amr);
  m_rte->setComputationalGeometry(m_computationalGeometry);
  m_rte->setRealm(m_particleRealm);
  m_rte->sanityCheck();
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::setupPoisson() noexcept
{
  CH_TIME("ItoKMCStepper::setupPoisson");
  if (m_verbosity > 5) {
    pout() << m_name + "::setupPoisson" << endl;
  }

  m_fieldSolver = RefCountedPtr<FieldSolver>(new F());
  m_fieldSolver->parseOptions();
  m_fieldSolver->setAmr(m_amr);
  m_fieldSolver->setComputationalGeometry(m_computationalGeometry);
  m_fieldSolver->setVoltage(m_voltage);
  m_fieldSolver->setRealm(m_fluidRealm);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::setupSigma() noexcept
{
  CH_TIME("ItoKMCStepper::setupSigma");
  if (m_verbosity > 5) {
    pout() << m_name + "::setupSigma" << endl;
  }

  m_sigmaSolver = RefCountedPtr<SurfaceODESolver<1>>(new SurfaceODESolver<1>(m_amr));
  m_sigmaSolver->parseOptions();
  m_sigmaSolver->setRealm(m_fluidRealm);
  m_sigmaSolver->setPhase(m_plasmaPhase);
  m_sigmaSolver->setName("Surface charge");
  m_sigmaSolver->setTime(0, 0.0, 0.0);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::allocate() noexcept
{
  CH_TIME("ItoKMCStepper::allocate");
  if (m_verbosity > 5) {
    pout() << m_name + "::allocate" << endl;
  }

  m_ito->allocate();
  m_cdr->allocate();
  m_rte->allocate();
  m_fieldSolver->allocate();
  m_sigmaSolver->allocate();

  this->allocateInternals();
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::allocateInternals() noexcept
{
  CH_TIME("ItoKMCStepper::allocateInternals");
  if (m_verbosity > 5) {
    pout() << m_name + "::allocateInternals" << endl;
  }

  const int numItoSpecies    = m_physics->getNumItoSpecies();
  const int numCdrSpecies    = m_physics->getNumCdrSpecies();
  const int numPlasmaSpecies = m_physics->getNumPlasmaSpecies();
  const int numPhotonSpecies = m_physics->getNumPhotonSpecies();

  CH_assert(numPlasmaSpecies > 0);

  // Scratch data.
  m_amr->allocate(m_fluidScratch1, m_fluidRealm, m_plasmaPhase, 1);
  m_amr->allocate(m_fluidScratchD, m_fluidRealm, m_plasmaPhase, SpaceDim);
  m_amr->allocate(m_fluidScratchEB, m_fluidRealm, m_plasmaPhase, 1);

  m_amr->allocate(m_particleScratch1, m_particleRealm, m_plasmaPhase, 1);
  m_amr->allocate(m_particleScratchD, m_particleRealm, m_plasmaPhase, SpaceDim);
  m_amr->allocate(m_particleScratchEB, m_particleRealm, m_plasmaPhase, 1);

  // Storage for conductiviies on cell, cell faces, and EB faces.
  m_amr->allocate(m_conductivityCell, m_fluidRealm, m_plasmaPhase, 1);
  m_amr->allocate(m_conductivityFace, m_fluidRealm, m_plasmaPhase, 1);
  m_amr->allocate(m_conductivityEB, m_fluidRealm, m_plasmaPhase, 1);

  // Electric field data on both realms.
  m_amr->allocate(m_electricFieldParticle, m_particleRealm, m_plasmaPhase, SpaceDim);
  m_amr->allocate(m_electricFieldFluid, m_fluidRealm, m_plasmaPhase, SpaceDim);

  // Data for CDR-based solver mobilities
  m_cdrMobilities.resize(numCdrSpecies);
  m_cdrPhotoiProducts.resize(numCdrSpecies);
  for (int i = 0; i < numCdrSpecies; i++) {
    m_amr->allocate(m_cdrMobilities[i], m_fluidRealm, m_plasmaPhase, 1);
    m_amr->allocate(m_cdrPhotoiProducts[i], m_particleRealm);
  }

  // Storage for the density gradients
  m_fluidGradPhiCDR.resize(numCdrSpecies);
  m_fluidGradPhiIto.resize(numItoSpecies);
  for (int i = 0; i < numCdrSpecies; i++) {
    m_amr->allocate(m_fluidGradPhiCDR[i], m_fluidRealm, m_plasmaPhase, SpaceDim);
  }
  for (int i = 0; i < numItoSpecies; i++) {
    m_amr->allocate(m_fluidGradPhiIto[i], m_fluidRealm, m_plasmaPhase, SpaceDim);
  }  

  // Storage for secondary particle and photon emission
  m_secondaryParticles.resize(numItoSpecies);
  m_secondaryPhotons.resize(numPhotonSpecies);

  m_cdrFluxes.resize(numCdrSpecies);
  m_cdrFluxesExtrap.resize(numCdrSpecies);

  for (int i = 0; i < numItoSpecies; i++) {
    m_amr->allocate(m_secondaryParticles[i], m_particleRealm);
  }

  for (int i = 0; i < numPhotonSpecies; i++) {
    m_amr->allocate(m_secondaryPhotons[i], m_particleRealm);
  }

  for (int i = 0; i < numCdrSpecies; i++) {
    m_amr->allocate(m_cdrFluxes[i], m_particleRealm, m_plasmaPhase, 1);
    m_amr->allocate(m_cdrFluxesExtrap[i], m_particleRealm, m_plasmaPhase, 1);
  }

  // Current density.
  m_amr->allocate(m_currentDensity, m_fluidRealm, m_plasmaPhase, SpaceDim);

  // Storage required for the reaction network.
  m_amr->allocate(m_oldPPC, m_particleRealm, m_plasmaPhase, numPlasmaSpecies);
  m_amr->allocate(m_newPPC, m_particleRealm, m_plasmaPhase, numPlasmaSpecies);
  m_amr->allocate(m_fluidPPC, m_fluidRealm, m_plasmaPhase, numPlasmaSpecies);

  if (numItoSpecies > 0) {
    m_amr->allocate(m_particleItoPPC, m_particleRealm, m_plasmaPhase, numItoSpecies);
    m_amr->allocate(m_particleOldItoPPC, m_particleRealm, m_plasmaPhase, numItoSpecies);
  }
  else {
    // Allocate some dummy data -- makes it easier. Trust me.
    m_amr->allocate(m_particleItoPPC, m_particleRealm, m_plasmaPhase, 1);
    m_amr->allocate(m_particleOldItoPPC, m_particleRealm, m_plasmaPhase, 1);
  }

  if (numCdrSpecies > 0) {
    m_amr->allocate(m_fluidCdrPPC, m_fluidRealm, m_plasmaPhase, numCdrSpecies);
    m_amr->allocate(m_fluidOldCdrPPC, m_fluidRealm, m_plasmaPhase, numCdrSpecies);
  }
  else {
    m_amr->allocatePointer(m_fluidCdrPPC, m_fluidRealm);
    m_amr->allocatePointer(m_fluidOldCdrPPC, m_fluidRealm);
  }

  if (numPhotonSpecies > 0) {
    m_amr->allocate(m_particleYPC, m_particleRealm, m_plasmaPhase, numPhotonSpecies);
    m_amr->allocate(m_fluidYPC, m_fluidRealm, m_plasmaPhase, numPhotonSpecies);
  }
  else {
    // Allocate some dummy data -- makes it easier. Trust me.
    m_amr->allocate(m_particleYPC, m_particleRealm, m_plasmaPhase, 1);
    m_amr->allocate(m_fluidYPC, m_fluidRealm, m_plasmaPhase, 1);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::postInitialize() noexcept
{
  CH_TIME("ItoKMCStepper::postInitialize");
  if (m_verbosity > 5) {
    pout() << m_name + "::postInitialize" << endl;
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::initialData() noexcept
{
  CH_TIME("ItoKMCStepper::initialData");
  if (m_verbosity > 5) {
    pout() << m_name + "::initialData" << endl;
  }

  CH_assert(!(m_cdr.isNull()));
  CH_assert(!(m_ito.isNull()));
  CH_assert(!(m_rte.isNull()));
  CH_assert(!(m_sigmaSolver.isNull()));
  CH_assert(!(m_fieldSolver.isNull()));

  m_ito->initialData();
  m_cdr->initialData();
  m_rte->initialData();
  this->initialSigma();

  // Make superparticles.
  m_ito->organizeParticlesByCell(ItoSolver::WhichContainer::Bulk);
  m_ito->makeSuperparticles(ItoSolver::WhichContainer::Bulk, m_particlesPerCell);
  m_ito->organizeParticlesByPatch(ItoSolver::WhichContainer::Bulk);

  // Solve Poisson equation and compute the E-field
  m_fieldSolver->setPermittivities();
  this->solvePoisson();

  // Fill solvers with velocities and diffusion coefficients
  this->computeDriftVelocities();
  this->computeDiffusionCoefficients();
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::initialSigma() noexcept
{
  CH_TIME("ItoKMCStepper::initialSigma");
  if (m_verbosity > 5) {
    pout() << m_name + "::initialSigma" << endl;
  }

  const RealVect probLo = m_amr->getProbLo();

  EBAMRIVData& sigma = m_sigmaSolver->getPhi();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_sigmaSolver->getRealm())[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_sigmaSolver->getRealm(), m_sigmaSolver->getPhase())[lvl];
    const Real               dx    = m_amr->getDx()[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      BaseIVFAB<Real>& phi     = (*sigma[lvl])[dit()];
      const EBISBox&   ebisbox = ebisl[dit()];

      CH_assert(phi.nComp() == 1);

      auto kernel = [&](const VolIndex& vof) -> void {
        const RealVect pos = probLo + Location::position(Location::Cell::Boundary, vof, ebisbox, dx);

        phi(vof, 0) = m_physics->initialSigma(m_time, pos);
      };

      VoFIterator& vofit = (*m_amr->getVofIterator(m_sigmaSolver->getRealm(), m_sigmaSolver->getPhase())[lvl])[dit()];

      BoxLoops::loop(vofit, kernel);
    }
  }

  // Coarsen throughout the AMR hierarchy.
  m_amr->conservativeAverage(sigma, m_fluidRealm, m_sigmaSolver->getPhase());

  // Set surface charge to zero on electrode cut-cells.
  m_sigmaSolver->resetElectrodes(sigma, 0.0);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::postCheckpointSetup() noexcept
{
  CH_TIME("ItoKMCStepper::postCheckpointSetup");
  if (m_verbosity > 5) {
    pout() << m_name + "::postCheckpointSetup" << endl;
  }

  // Allocate internal storage.
  this->allocateInternals();

  m_ito->remap();

  // Recompute the electric field.
  this->postCheckpointPoisson();

  // Compute velocities and diffusion coefficients so we're prepared for the next time step.
  this->computeDriftVelocities();
  this->computeDiffusionCoefficients();
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::postCheckpointPoisson() noexcept
{
  CH_TIME("ItoKMCStepper::postCheckpointPoisson");
  if (m_verbosity > 5) {
    pout() << m_name + "::postCheckpointPoisson" << endl;
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

  // Set up the Poisson solver
  m_fieldSolver->setupSolver();
}

#ifdef CH_USE_HDF5
template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::writeCheckpointHeader(HDF5HeaderData& a_header) const noexcept
{
  CH_TIME("ItoKMCStepper::writeCheckpointHeader");
  if (m_verbosity > 5) {
    pout() << m_name + "::writeCheckpointHeader" << endl;
  }
}
#endif

#ifdef CH_USE_HDF5
template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::readCheckpointHeader(HDF5HeaderData& a_header) noexcept
{
  CH_TIME("ItoKMCStepper::readCheckpointHeader");
  if (m_verbosity > 5) {
    pout() << m_name + "::readCheckpointHeader" << endl;
  }
}
#endif

#ifdef CH_USE_HDF5
template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const noexcept
{
  CH_TIME("ItoKMCStepper::writeCheckpointData");
  if (m_verbosity > 5) {
    pout() << m_name + "::writeCheckpointData" << endl;
  }

  for (ItoIterator<ItoSolver> solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->writeCheckpointLevel(a_handle, a_lvl);
  }

  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->writeCheckpointLevel(a_handle, a_lvl);
  }

  for (RtIterator<McPhoto> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->writeCheckpointLevel(a_handle, a_lvl);
  }

  m_fieldSolver->writeCheckpointLevel(a_handle, a_lvl);
  m_sigmaSolver->writeCheckpointLevel(a_handle, a_lvl);
}
#endif

#ifdef CH_USE_HDF5
template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::readCheckpointData(HDF5Handle& a_handle, const int a_lvl) noexcept
{
  CH_TIME("ItoKMCStepper::readCheckpointData");
  if (m_verbosity > 5) {
    pout() << m_name + "::readCheckpointData" << endl;
  }

  for (ItoIterator<ItoSolver> solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->readCheckpointLevel(a_handle, a_lvl);
  }

  for (CdrIterator<CdrSolver> solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->readCheckpointLevel(a_handle, a_lvl);
  }

  for (RtIterator<McPhoto> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->readCheckpointLevel(a_handle, a_lvl);
  }

  m_fieldSolver->readCheckpointLevel(a_handle, a_lvl);
  m_sigmaSolver->readCheckpointLevel(a_handle, a_lvl);
}
#endif

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::writePlotData(EBAMRCellData&       a_output,
                                         Vector<std::string>& a_plotVariableNames,
                                         int&                 a_icomp) const noexcept
{
  CH_TIME("ItoKMCStepper::writePlotData");
  if (m_verbosity > 5) {
    pout() << m_name + "::writePlotData" << endl;
  }

  // Poisson solver copies over its output data
  a_plotVariableNames.append(m_fieldSolver->getPlotVariableNames());
  m_fieldSolver->writePlotData(a_output, a_icomp);

  // Surface charge solver writes
  a_plotVariableNames.append(m_sigmaSolver->getPlotVariableNames());
  m_sigmaSolver->writePlotData(a_output, a_icomp);

  // Ito solvers copy their output data
  m_ito->depositParticles();
  for (ItoIterator<ItoSolver> solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>& solver = solverIt();

    a_plotVariableNames.append(solver->getPlotVariableNames());
    solver->writePlotData(a_output, a_icomp);
  }

  // Cdr solvers output their data
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<CdrSolver>& solver = solverIt();

    a_plotVariableNames.append(solver->getPlotVariableNames());
    solver->writePlotData(a_output, a_icomp);
  }

  // RTE solvers copy their output data
  for (RtIterator<McPhoto> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<McPhoto>& solver = solverIt();

    a_plotVariableNames.append(solver->getPlotVariableNames());
    solver->writePlotData(a_output, a_icomp);
  }

  // Write the conductivity to the output
  if (m_plotConductivity) {
    a_output.copy(Interval(0, 0), m_conductivityCell, Interval(a_icomp, a_icomp));
    a_plotVariableNames.push_back("Conductivity");
    a_icomp++;
  }

  // Write the current to the output
  if (m_plotCurrentDensity) {
    const Interval srcInterv(0, SpaceDim - 1);
    const Interval dstInterv(a_icomp, a_icomp + SpaceDim - 1);
    a_output.copy(srcInterv, m_currentDensity, dstInterv);

    a_plotVariableNames.push_back("x-J");
    a_plotVariableNames.push_back("y-J");
    if (SpaceDim == 3) {
      a_plotVariableNames.push_back("z-J");
    }

    a_icomp += SpaceDim;
  }

  // Write the number of particles per patch
  if (m_plotParticlesPerPatch) {
    this->writeNumberOfParticlesPerPatch(a_output, a_icomp);
    a_plotVariableNames.push_back("Particles per patch");
    a_icomp++;
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::writeNumberOfParticlesPerPatch(EBAMRCellData& a_output, int& a_icomp) const noexcept
{
  CH_TIME("ItoKMCStepper::writeNumberOfParticlesPerPatch");
  if (m_verbosity > 5) {
    pout() << m_name + "::writeNumberOfParticlesPerPatch" << endl;
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

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt) noexcept
{
  CH_TIME("ItoKMCStepper::synchronizeSolverTimes");
  if (m_verbosity > 5) {
    pout() << m_name + "::synchronizeSolverTimes" << endl;
  }

  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;

  m_ito->setTime(a_step, a_time, a_dt);
  m_fieldSolver->setTime(a_step, a_time, a_dt);
  m_rte->setTime(a_step, a_time, a_dt);
  m_sigmaSolver->setTime(a_step, a_time, a_dt);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::printStepReport() noexcept
{
  CH_TIME("ItoKMCStepper::printStepReport");
  if (m_verbosity > 5) {
    pout() << m_name + "::printStepReport" << endl;
  }

  const Real Emax = this->computeMaxElectricField(m_plasmaPhase);

  const unsigned long long localParticlesBulk    = m_ito->getNumParticles(ItoSolver::WhichContainer::Bulk, true);
  const unsigned long long globalParticlesBulk   = m_ito->getNumParticles(ItoSolver::WhichContainer::Bulk, false);
  const unsigned long long localParticlesEB      = m_ito->getNumParticles(ItoSolver::WhichContainer::EB, true);
  const unsigned long long globalParticlesEB     = m_ito->getNumParticles(ItoSolver::WhichContainer::EB, false);
  const unsigned long long localParticlesDomain  = m_ito->getNumParticles(ItoSolver::WhichContainer::Domain, true);
  const unsigned long long globalParticlesDomain = m_ito->getNumParticles(ItoSolver::WhichContainer::Domain, false);
  const unsigned long long localParticlesSource  = m_ito->getNumParticles(ItoSolver::WhichContainer::Source, true);
  const unsigned long long globalParticlesSource = m_ito->getNumParticles(ItoSolver::WhichContainer::Source, false);

  Real avgParticles = 0.0;
  Real stdDev       = 0.0;

  Real minParticles = 0.0;
  Real maxParticles = 0.0;

  int minRank = 0;
  int maxRank = 0;

  this->getParticleStatistics(avgParticles, stdDev, minParticles, maxParticles, minRank, maxRank);

  Real maxDensity = -std::numeric_limits<Real>::max();
  Real minDensity = +std::numeric_limits<Real>::max();

  std::string maxSolver = "invalid solver";
  std::string minSolver = "invalid solver";

  this->getMaxMinItoDensity(maxDensity, minDensity, maxSolver, minSolver);
  this->getMaxMinCDRDensity(maxDensity, minDensity, maxSolver, minSolver);

  std::string str;
  switch (m_timeCode) {
  case TimeCode::Physics: {
    str = "dt restricted by 'Physics'";

    break;
  }
  case TimeCode::AdvectionIto: {
    str = "dt restricted by 'Advection (Ito)'";

    break;
  }
  case TimeCode::DiffusionIto: {
    str = "dt restricted by 'Diffusion (Ito)'";

    break;
  }
  case TimeCode::AdvectionDiffusionIto: {
    str = "dt restricted by 'AdvectionDiffusion (Ito)'";

    break;
  }
  case TimeCode::AdvectionDiffusionCDR: {
    str = "dt restricted by 'AdvectionDiffusion (CDR)'";

    break;
  }
  case TimeCode::RelaxationTime: {
    str = "dt restricted by 'Relaxation time'";

    break;
  }
  case TimeCode::Hardcap: {
    str = "dt restricted by 'Hardcap'";

    break;
  }
  default: {
    str = "dt restricted by 'Unspecified'";

    break;
  }
  }

  // Print the step report.

  //clang-format off
  const std::string whitespace = "                                   ";
  pout() << "                                   " + str << endl;
  pout() << whitespace + "Emax        = " << Emax << endl
         << whitespace + "Max density = " << maxDensity << " (" << maxSolver << ")" << endl
         << whitespace + "CFL (Ito)   = " << m_dt / m_particleAdvectionDiffusionDt << endl
         << whitespace + "CFL (CDR)   = " << m_dt / m_fluidAdvectionDiffusionDt << endl
         << whitespace + "dt/dt_relax = " << m_dt / m_relaxationTime << endl
         << whitespace + "#Particles  = " << DischargeIO::numberFmt(localParticlesBulk) << " ("
         << DischargeIO::numberFmt(globalParticlesBulk) << ")" << endl
         << whitespace + "#EB part.   = " << DischargeIO::numberFmt(localParticlesEB) << " ("
         << DischargeIO::numberFmt(globalParticlesEB) << ")" << endl
         << whitespace + "#Dom. part. = " << DischargeIO::numberFmt(localParticlesDomain) << " ("
         << DischargeIO::numberFmt(globalParticlesDomain) << ")" << endl
         << whitespace + "#Src. part. = " << DischargeIO::numberFmt(localParticlesSource) << " ("
         << DischargeIO::numberFmt(globalParticlesSource) << ")" << endl
         << whitespace + "#Min part.  = " << minParticles << " (on rank = " << minRank << ")" << endl
         << whitespace + "#Max part.  = " << maxParticles << " (on rank = " << maxRank << ")" << endl
         << whitespace + "#Avg. part. = " << avgParticles << endl
         << whitespace + "#Dev. part. = " << stdDev << " (" << 100. * stdDev / avgParticles << "%)" << endl;
  //clang-format on
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::getMaxMinItoDensity(Real&        a_maxDensity,
                                               Real&        a_minDensity,
                                               std::string& a_maxSolver,
                                               std::string& a_minSolver) const noexcept
{
  CH_TIME("ItoKMCStepper::getMaxMinDensity(Realx2, std::string2x)");
  if (m_verbosity > 5) {
    pout() << m_name + "::getMaxMinDensity(Realx2, std::string2x)" << endl;
  }

  // Go through each solver and find the max/min values.
  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    const RefCountedPtr<ItoSolver>& solver = solverIt();

    Real curMin = std::numeric_limits<Real>::max();
    Real curMax = -std::numeric_limits<Real>::max();

    DataOps::getMaxMin(curMax, curMin, solverIt()->getPhi(), 0);

    if (curMax > a_maxDensity) {
      a_maxDensity = curMax;
      a_maxSolver  = solver->getName();
    }

    if (curMin < a_minDensity) {
      a_minDensity = curMin;
      a_minSolver  = solver->getName();
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::getMaxMinCDRDensity(Real&        a_maxDensity,
                                               Real&        a_minDensity,
                                               std::string& a_maxSolver,
                                               std::string& a_minSolver) const noexcept
{
  CH_TIME("ItoKMCStepper::getMaxMinCDRDensity(Realx2, std::string2x)");
  if (m_verbosity > 5) {
    pout() << m_name + "::getMaxMinCDRDensity(Realx2, std::string2x)" << endl;
  }

  // Go through each solver and find the max/min values.
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    const RefCountedPtr<CdrSolver>& solver = solverIt();

    Real curMin = std::numeric_limits<Real>::max();
    Real curMax = -std::numeric_limits<Real>::max();

    DataOps::getMaxMin(curMax, curMin, solverIt()->getPhi(), 0);

    if (curMax > a_maxDensity) {
      a_maxDensity = curMax;
      a_maxSolver  = solver->getName();
    }

    if (curMin < a_minDensity) {
      a_minDensity = curMin;
      a_minSolver  = solver->getName();
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::getParticleStatistics(Real& a_avgParticles,
                                                 Real& a_sigma,
                                                 Real& a_minParticles,
                                                 Real& a_maxParticles,
                                                 int&  a_minRank,
                                                 int&  a_maxRank)
{
  CH_TIME("ItoKMCStepper::getParticleStatistics");
  if (m_verbosity > 5) {
    pout() << m_name + "::getParticleStatistics" << endl;
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

template <typename I, typename C, typename R, typename F>
Real
ItoKMCStepper<I, C, R, F>::computeDt()
{
  CH_TIME("ItoKMCStepper::computeDt");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDt" << endl;
  }

  Timer timer(m_name + "::computeDt");

  Real dt = std::numeric_limits<Real>::max();

  // Compute various time steps.
  timer.startEvent("Advection (Ito)");
  m_particleAdvectionDt = m_ito->computeAdvectiveDt();
  timer.stopEvent("Advection (Ito)");

  timer.startEvent("Diffusion (Ito)");
  m_particleDiffusionDt = m_ito->computeDiffusiveDt();
  timer.stopEvent("Diffusion (Ito)");

  timer.startEvent("AdvectionDiffusion (Ito)");
  m_particleAdvectionDiffusionDt = m_ito->computeDt();
  timer.stopEvent("AdvectionDiffusion (Ito)");

  timer.startEvent("AdvectionDiffusion (CDR)");
  m_fluidAdvectionDiffusionDt = m_cdr->computeAdvectionDiffusionDt();
  timer.stopEvent("AdvectionDiffusion (CDR)");

  timer.startEvent("Physics");
  m_physicsDt = this->computePhysicsDt();
  timer.stopEvent("Physics");

  timer.startEvent("Relaxation");
  m_relaxationTime = this->computeRelaxationTime();
  timer.stopEvent("Relaxation");

  if (m_maxParticleAdvectionCFL * m_particleAdvectionDt < dt) {
    dt         = m_maxParticleAdvectionCFL * m_particleAdvectionDt;
    m_timeCode = TimeCode::AdvectionIto;
  }

  if (m_maxParticleDiffusionCFL * m_particleDiffusionDt < dt) {
    dt         = m_maxParticleDiffusionCFL * m_particleDiffusionDt;
    m_timeCode = TimeCode::DiffusionIto;
  }

  if (m_maxParticleAdvectionDiffusionCFL * m_particleAdvectionDiffusionDt < dt) {
    dt         = m_maxParticleAdvectionDiffusionCFL * m_particleAdvectionDiffusionDt;
    m_timeCode = TimeCode::AdvectionDiffusionIto;
  }

  if (std::min(m_fluidAdvectionDiffusionCFL, 0.9) * m_fluidAdvectionDiffusionDt < dt) {
    dt         = std::min(m_fluidAdvectionDiffusionCFL, 0.9) * m_fluidAdvectionDiffusionDt;
    m_timeCode = TimeCode::AdvectionDiffusionCDR;
  }

  if (m_relaxTimeFactor * m_relaxationTime < dt) {
    dt         = m_relaxTimeFactor * m_relaxationTime;
    m_timeCode = TimeCode::RelaxationTime;
  }

  if (m_physicsDt < dt) {
    dt         = m_physicsDt;
    m_timeCode = TimeCode::Physics;
  }

  if (dt < m_minParticleAdvectionCFL * m_particleAdvectionDt) {
    dt = m_minParticleAdvectionCFL * m_particleAdvectionDt;

    m_timeCode = TimeCode::AdvectionIto;
  }

  if (dt < m_minParticleDiffusionCFL * m_particleDiffusionDt) {
    dt = m_minParticleDiffusionCFL * m_particleDiffusionDt;

    m_timeCode = TimeCode::DiffusionIto;
  }

  if (dt < m_minParticleAdvectionDiffusionCFL * m_particleAdvectionDiffusionDt) {
    dt = m_minParticleAdvectionDiffusionCFL * m_particleAdvectionDiffusionDt;

    m_timeCode = TimeCode::AdvectionDiffusionIto;
  }

  if (m_minDt > dt) {
    dt         = m_minDt;
    m_timeCode = TimeCode::Hardcap;
  }

  if (m_maxDt < dt) {
    dt         = m_maxDt;
    m_timeCode = TimeCode::Hardcap;
  }

  if (m_profile) {
    timer.eventReport(pout(), false);
  }

  return dt;
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::registerRealms() noexcept
{
  CH_TIME("ItoKMCStepper::registerRealms");
  if (m_verbosity > 5) {
    pout() << m_name + "::registerRealms" << endl;
  }

  // TLDR: If using dual grid then m_particleRealm != m_fluidRealm and we'll have two realms.
  m_amr->registerRealm(m_fluidRealm);
  m_amr->registerRealm(m_particleRealm);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::registerOperators() noexcept
{
  CH_TIME("ItoKMCStepper::registerOperators");
  if (m_verbosity > 5) {
    pout() << m_name + "::registerOperators" << endl;
  }

  m_ito->registerOperators();
  m_cdr->registerOperators();
  m_fieldSolver->registerOperators();
  m_rte->registerOperators();
  m_sigmaSolver->registerOperators();
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::preRegrid(const int a_lmin, const int a_oldFinestLevel) noexcept
{
  CH_TIME("ItoKMCStepper::preRegrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::preRegrid" << endl;
  }

  // If we are load balancing then we need to store the number of particles per cell on the old grids. This
  // will be used to estimate computational loads on the new grids.
  if (m_loadBalanceParticles) {
    Vector<RefCountedPtr<ItoSolver>> lbSolvers = this->getLoadBalanceSolvers();

    m_loadBalancePPC.resize(lbSolvers.size());

    // Allocate and compute number of computational particles per cell.
    for (int i = 0; i < lbSolvers.size(); i++) {
      m_amr->allocate(m_loadBalancePPC[i], m_particleRealm, m_plasmaPhase, 1);

      EBAMRCellData&                        compPPC   = m_loadBalancePPC[i];
      const ParticleContainer<ItoParticle>& particles = lbSolvers[i]->getParticles(ItoSolver::WhichContainer::Bulk);

      ParticleOps::getComputationalParticlesPerCell(compPPC, particles);
    }
  }

  m_ito->preRegrid(a_lmin, a_oldFinestLevel);
  m_cdr->preRegrid(a_lmin, a_oldFinestLevel);
  m_fieldSolver->preRegrid(a_lmin, a_oldFinestLevel);
  m_rte->preRegrid(a_lmin, a_oldFinestLevel);
  m_sigmaSolver->preRegrid(a_lmin, a_oldFinestLevel);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) noexcept
{
  CH_TIME("ItoKMCStepper::regrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::regrid" << endl;
  }

  this->allocateInternals();

  m_ito->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_cdr->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_fieldSolver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_rte->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_sigmaSolver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);

  if (m_regridSuperparticles) {
    m_ito->organizeParticlesByCell(ItoSolver::WhichContainer::Bulk);
    m_ito->makeSuperparticles(ItoSolver::WhichContainer::Bulk, m_particlesPerCell);
    m_ito->organizeParticlesByPatch(ItoSolver::WhichContainer::Bulk);
  }

  // Redeposit particles and update the electric field on the new mesh.
  m_ito->depositParticles();

  const bool converged = this->solvePoisson();
  if (!converged) {
    const std::string err = "ItoKMCStepper::regrid - Poisson solve did not converge after regrid!!!";

    if (m_abortOnFailure) {
      MayDay::Error(err.c_str());
    }
    else {
      MayDay::Warning(err.c_str());
    }
  }

  this->computeDriftVelocities();
  this->computeDiffusionCoefficients();
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::postRegrid() noexcept
{
  CH_TIME("ItoKMCStepper::postRegrid");

  if (m_loadBalanceParticles) {
    for (int i = 0; i < m_loadBalancePPC.size(); i++) {
      m_amr->deallocate(m_loadBalancePPC[i]);
    }
  }
}

template <typename I, typename C, typename R, typename F>
int
ItoKMCStepper<I, C, R, F>::getNumberOfPlotVariables() const noexcept
{
  CH_TIME("ItoKMCStepper::getNumberOfPlotVariables");
  if (m_verbosity > 5) {
    pout() << m_name + "::getNumberOfPlotVariables" << endl;
  }

  int numComp = 0;

  // Ito solver variables.
  for (ItoIterator<ItoSolver> solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    numComp += solverIt()->getNumberOfPlotVariables();
  }

  // Cdr solver variables
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    numComp += solverIt()->getNumberOfPlotVariables();
  }

  // RTE solver variables.
  for (RtIterator<McPhoto> solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    numComp += solverIt()->getNumberOfPlotVariables();
  }

  // Field solver variables.
  numComp += m_fieldSolver->getNumberOfPlotVariables();

  // Surface charge solver variables.
  numComp += m_sigmaSolver->getNumberOfPlotVariables();

  // Conductivity
  if (m_plotConductivity) {
    numComp += 1;
  }

  // Current density.
  if (m_plotCurrentDensity) {
    numComp += SpaceDim;
  }

  // Number of particles per patch
  if (m_plotParticlesPerPatch) {
    numComp += 1;
  }

  return numComp;
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::setVoltage(const std::function<Real(const Real a_time)>& a_voltage) noexcept
{
  CH_TIME("ItoKMCStepper::setVoltage");
  if (m_verbosity > 5) {
    pout() << m_name + "::setVoltage" << endl;
  }

  m_voltage = a_voltage;
}

template <typename I, typename C, typename R, typename F>
Real
ItoKMCStepper<I, C, R, F>::computeMaxElectricField(const phase::which_phase a_phase) const noexcept
{
  CH_TIME("ItoKMCStepper::computeMaxElectricField");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeMaxElectricField" << endl;
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

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeElectricField(EBAMRCellData&           a_electricField,
                                                const phase::which_phase a_phase) const noexcept
{
  CH_TIME("ItoKMCStepper::computeElectricField(EBAMRCellData, phase)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeElectricField(EBAMRCellData, phase)" << endl;
  }

  CH_assert(a_electricField.getRealm() == m_fluidRealm);

  m_fieldSolver->computeElectricField(a_electricField, a_phase, m_fieldSolver->getPotential());
}

template <typename I, typename C, typename R, typename F>
Real
ItoKMCStepper<I, C, R, F>::getTime() const noexcept
{
  CH_TIME("ItoKMCStepper::getTime");
  if (m_verbosity > 5) {
    pout() << m_name + "::getTime" << endl;
  }

  return m_time;
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeSpaceChargeDensity() noexcept
{
  CH_TIME("ItoKMCStepper::computeSpaceChargeDensity()");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeSpaceChargeDensity()" << endl;
  }

  this->computeSpaceChargeDensity(m_fieldSolver->getRho(), m_ito->getDensities(), m_cdr->getPhis());
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeSpaceChargeDensity(MFAMRCellData&                a_rho,
                                                     const Vector<EBAMRCellData*>& a_itoDensities,
                                                     const Vector<EBAMRCellData*>& a_cdrDensities) noexcept
{
  CH_TIME("ItoKMCStepper::computeSpaceChargeDensity(rho, densities)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeSpaceChargeDensity(rho, densities)" << endl;
  }

  CH_assert(a_rho.getRealm() == m_fluidRealm);
  CH_assert(a_itoDensities[0]->getRealm() == m_particleRealm);
  CH_assert(a_cdrDensities[0]->getRealm() == m_fluidRealm);

  // TLDR: a_itoDensities could be defined over the particle realm, so we use m_fluidScratch1 as a temporary storage.

  DataOps::setValue(a_rho, 0.0);

  // Alias for the plasma phase.
  EBAMRCellData rhoPhase = m_amr->alias(m_plasmaPhase, a_rho);

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    const RefCountedPtr<ItoSolver>&  solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();
    const int                        idx     = solverIt.index();
    const int                        Z       = species->getChargeNumber();

    if (Z != 0) {
      m_fluidScratch1.copy(*a_itoDensities[idx]);

      DataOps::incr(rhoPhase, m_fluidScratch1, 1.0 * Z);
    }
  }

  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    const RefCountedPtr<CdrSolver>&  solver  = solverIt();
    const RefCountedPtr<CdrSpecies>& species = solver->getSpecies();
    const int                        idx     = solverIt.index();
    const int                        Z       = species->getChargeNumber();

    if (Z != 0) {
      DataOps::incr(rhoPhase, *a_cdrDensities[idx], 1.0 * Z);
    }
  }

  DataOps::scale(a_rho, Units::Qe);

  m_amr->arithmeticAverage(a_rho, m_fluidRealm);
  m_amr->interpGhostPwl(a_rho, m_fluidRealm);

  // Interpolate to centroids.
  m_amr->interpToCentroids(rhoPhase, m_fluidRealm, m_plasmaPhase);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeConductivityCell(EBAMRCellData& a_conductivity) noexcept
{
  CH_TIME("ItoKMCStepper::computeConductivityCell(EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeConductivityCell(EBAMRCellData)" << endl;
  }

  this->computeConductivityCell(a_conductivity, m_ito->getParticles(ItoSolver::WhichContainer::Bulk));
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeConductivityCell(EBAMRCellData&                                 a_conductivity,
                                                   const Vector<ParticleContainer<ItoParticle>*>& a_particles) noexcept
{
  CH_TIME("ItoKMCStepper::computeConductivityCell(EBAMRCellData, Particles)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeConductivityCell(EBAMRCellData, Particles)" << endl;
  }

  DataOps::setValue(a_conductivity, 0.0);

  // Add contribution from particle solvers
  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();
    const int Z   = species->getChargeNumber();

    if (Z != 0 && solver->isMobile()) {
      solver->depositConductivity(m_particleScratch1, *a_particles[idx]);

      // Add to the fluid realm.
      m_fluidScratch1.copy(m_particleScratch1);
      DataOps::incr(a_conductivity, m_fluidScratch1, 1.0 * std::abs(Z));
    }
  }

  // Add contribution from CDR solvers
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    const RefCountedPtr<CdrSolver>&  solver  = solverIt();
    const RefCountedPtr<CdrSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();
    const int Z   = species->getChargeNumber();

    if (Z != 0 && solver->isMobile()) {
      const EBAMRCellData& phi      = solver->getPhi();
      const EBAMRCellData& mobility = m_cdrMobilities[idx];

      DataOps::copy(m_fluidScratch1, phi);
      DataOps::multiplyScalar(m_fluidScratch1, mobility);
      DataOps::incr(a_conductivity, m_fluidScratch1, 1.0 * std::abs(Z));
    }
  }

  DataOps::scale(a_conductivity, Units::Qe);

  m_amr->arithmeticAverage(a_conductivity, m_fluidRealm, m_plasmaPhase);
  m_amr->interpGhostPwl(a_conductivity, m_fluidRealm, m_plasmaPhase);

  // Interpolate to centroids.
  m_amr->interpToCentroids(a_conductivity, m_fluidRealm, m_plasmaPhase);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeDensityGradients() noexcept
{
  CH_TIME("ItoKMCStepper::computeDensityGradients()");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDensityGradients()" << endl;
  }

  // Compute gradients for the CDR species.
  for (auto it = m_cdr->iterator(); it.ok(); ++it) {
    const RefCountedPtr<CdrSolver>& solver = it();

    const int idx = it.index();

    // Update ghost cells and coarsenings. Then compute the gradient.
    m_fluidScratch1.copy(solver->getPhi());

    m_amr->arithmeticAverage(m_fluidScratch1, m_fluidRealm, m_plasmaPhase);
    m_amr->interpGhostPwl(m_fluidScratch1, m_fluidRealm, m_plasmaPhase);

    m_amr->computeGradient(m_fluidGradPhiCDR[idx], m_fluidScratch1, m_fluidRealm, m_plasmaPhase);
  }

  // Do the same for the Ito species.
  for (auto it = m_ito->iterator(); it.ok(); ++it) {
    const RefCountedPtr<ItoSolver>& solver = it();

    const int idx = it.index();

    // Update ghost cells and coarsenings. Then compute the gradient.
    m_fluidScratch1.copy(solver->getPhi());

    m_amr->arithmeticAverage(m_fluidScratch1, m_fluidRealm, m_plasmaPhase);
    m_amr->interpGhostPwl(m_fluidScratch1, m_fluidRealm, m_plasmaPhase);

    m_amr->computeGradient(m_fluidGradPhiIto[idx], m_fluidScratch1, m_fluidRealm, m_plasmaPhase);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeCurrentDensity(EBAMRCellData& a_J) noexcept
{
  CH_TIME("ItoKMCStepper::computeCurrentDensity(EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeCurrentDensity(EBAMRCellData)" << endl;
  }

  CH_assert(a_J[0]->nComp() == SpaceDim);

  // TLDR: a_J is defined over the fluid Realm but the computation takes place on the particle Realm.
  //       If the Realms are different we compute on a scratch storage instead

  this->computeConductivityCell(m_fluidScratch1);
  DataOps::copy(a_J, m_electricFieldFluid);

  DataOps::multiplyScalar(a_J, m_fluidScratch1);
}

template <typename I, typename C, typename R, typename F>
Real
ItoKMCStepper<I, C, R, F>::computeRelaxationTime() noexcept
{
  CH_TIME("ItoKMCStepper::computeRelaxationTime()");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeRelaxationTime()" << endl;
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

  Real min = std::numeric_limits<Real>::max();
  Real max = -std::numeric_limits<Real>::max();

  DataOps::getMaxMinNorm(max, min, relaxTime);

  return min;
}

template <typename I, typename C, typename R, typename F>
bool
ItoKMCStepper<I, C, R, F>::solvePoisson() noexcept
{
  CH_TIME("ItoKMCStepper::solvePoisson()");
  if (m_verbosity > 5) {
    pout() << m_name + "::solvePoisson()" << endl;
  }

  // Solve the Poisson equation and compute the cell-centered electric field.
  MFAMRCellData& phi   = m_fieldSolver->getPotential();
  MFAMRCellData& rho   = m_fieldSolver->getRho();
  EBAMRIVData&   sigma = m_sigmaSolver->getPhi();

  const bool converged = m_fieldSolver->solve(phi, rho, sigma, false);

  m_fieldSolver->computeElectricField();

  // Copy the electric field to appropriate data holders and perform center-to-centroid
  // interpolation.
  EBAMRCellData E;
  m_amr->allocatePointer(E, m_fluidRealm);
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

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::intersectParticles(const SpeciesSubset                     a_speciesSubset,
                                              const bool                              a_delete,
                                              const std::function<void(ItoParticle&)> a_nonDeletionModifier) noexcept
{
  CH_TIME("ItoKMCStepper::intersectParticles(SpeciesSubset, bool, std::function)");
  if (m_verbosity > 5) {
    pout() << m_name + "::intersectParticles(SpeciesSubset, bool, std::function)" << endl;
  }

  this->intersectParticles(a_speciesSubset,
                           ItoSolver::WhichContainer::Bulk,
                           ItoSolver::WhichContainer::EB,
                           ItoSolver::WhichContainer::Domain,
                           a_delete,
                           a_nonDeletionModifier);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::intersectParticles(const SpeciesSubset                     a_speciesSubset,
                                              const ItoSolver::WhichContainer         a_containerBulk,
                                              const ItoSolver::WhichContainer         a_containerEB,
                                              const ItoSolver::WhichContainer         a_containerDomain,
                                              const bool                              a_delete,
                                              const std::function<void(ItoParticle&)> a_nonDeletionModifier) noexcept
{
  CH_TIME("ItoKMCStepper::intersectParticles(SpeciesSubset, Containerx3, bool, std::function)");
  if (m_verbosity > 5) {
    pout() << m_name + "::intersectParticles(SpeciesSubset, Containerx3, bool, std::function)" << endl;
  }

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = (species->getChargeNumber() != 0);

    const EBIntersection intersectionAlgorithm = solver->getIntersectionAlgorithm();

    switch (a_speciesSubset) {
    case SpeciesSubset::All: {
      solver->intersectParticles(a_containerBulk,
                                 a_containerEB,
                                 a_containerDomain,
                                 intersectionAlgorithm,
                                 a_delete,
                                 a_nonDeletionModifier);

      break;
    }
    case SpeciesSubset::AllMobile: {
      if (mobile) {
        solver->intersectParticles(a_containerBulk,
                                   a_containerEB,
                                   a_containerDomain,
                                   intersectionAlgorithm,
                                   a_delete,
                                   a_nonDeletionModifier);
      }

      break;
    }
    case SpeciesSubset::AllDiffusive: {
      if (diffusive) {
        solver->intersectParticles(a_containerBulk,
                                   a_containerEB,
                                   a_containerDomain,
                                   intersectionAlgorithm,
                                   a_delete,
                                   a_nonDeletionModifier);
      }

      break;
    }
    case SpeciesSubset::AllMobileOrDiffusive: {
      if (mobile || diffusive) {
        solver->intersectParticles(a_containerBulk,
                                   a_containerEB,
                                   a_containerDomain,
                                   intersectionAlgorithm,
                                   a_delete,
                                   a_nonDeletionModifier);
      }

      break;
    }
    case SpeciesSubset::AllMobileAndDiffusive: {
      if (mobile && diffusive) {
        solver->intersectParticles(a_containerBulk,
                                   a_containerEB,
                                   a_containerDomain,
                                   intersectionAlgorithm,
                                   a_delete,
                                   a_nonDeletionModifier);
      }

      break;
    }
    case SpeciesSubset::Charged: {
      if (charged) {
        solver->intersectParticles(a_containerBulk,
                                   a_containerEB,
                                   a_containerDomain,
                                   intersectionAlgorithm,
                                   a_delete,
                                   a_nonDeletionModifier);
      }

      break;
    }
    case SpeciesSubset::ChargedMobile: {
      if (charged && mobile) {
        solver->intersectParticles(a_containerBulk,
                                   a_containerEB,
                                   a_containerDomain,
                                   intersectionAlgorithm,
                                   a_delete,
                                   a_nonDeletionModifier);
      }

      break;
    }
    case SpeciesSubset::ChargedDiffusive: {
      if (charged && diffusive) {
        solver->intersectParticles(a_containerBulk,
                                   a_containerEB,
                                   a_containerDomain,
                                   intersectionAlgorithm,
                                   a_delete,
                                   a_nonDeletionModifier);
      }

      break;
    }
    case SpeciesSubset::ChargedMobileOrDiffusive: {
      if (charged && (mobile || diffusive)) {
        solver->intersectParticles(a_containerBulk,
                                   a_containerEB,
                                   a_containerDomain,
                                   intersectionAlgorithm,
                                   a_delete,
                                   a_nonDeletionModifier);
      }

      break;
    }
    case SpeciesSubset::ChargedMobileAndDiffusive: {
      if (charged && (mobile && diffusive)) {
        solver->intersectParticles(a_containerBulk,
                                   a_containerEB,
                                   a_containerDomain,
                                   intersectionAlgorithm,
                                   a_delete,
                                   a_nonDeletionModifier);
      }

      break;
    }
    case SpeciesSubset::Stationary: {
      if (!mobile && !diffusive) {
        solver->intersectParticles(a_containerBulk,
                                   a_containerEB,
                                   a_containerDomain,
                                   intersectionAlgorithm,
                                   a_delete,
                                   a_nonDeletionModifier);
      }

      break;
    }
    default: {
      MayDay::Abort("ItoKMCStepper::intersectParticles - logic bust");

      break;
    }
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::removeCoveredParticles(const SpeciesSubset    a_speciesSubset,
                                                  const EBRepresentation a_representation,
                                                  const Real             a_tolerance) noexcept
{
  CH_TIME("ItoKMCStepper::removeCoveredParticles(SpeciesSubset, EBRepresentation, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::removeCoveredParticles(SpeciesSubset, EBRepresentation, Real)" << endl;
  }

  this->removeCoveredParticles(a_speciesSubset, ItoSolver::WhichContainer::Bulk, a_representation, a_tolerance);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::removeCoveredParticles(const SpeciesSubset             a_which,
                                                  const ItoSolver::WhichContainer a_container,
                                                  const EBRepresentation          a_representation,
                                                  const Real                      a_tolerance) noexcept
{
  CH_TIME("ItoKMCStepper::removeCoveredParticles(SpeciesSubset, container, EBRepresentation, tolerance)");
  if (m_verbosity > 5) {
    pout() << m_name + "::removeCoveredParticles(SpeciesSubset, container, EBRepresentation, tolerance)" << endl;
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
    case SpeciesSubset::AllMobileOrDiffusive: {
      if (mobile || diffusive) {
        solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      }

      break;
    }
    case SpeciesSubset::AllMobileAndDiffusive: {
      if (mobile && diffusive) {
        solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      }

      break;
    }
    case SpeciesSubset::Charged: {
      if (charged) {
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
    case SpeciesSubset::ChargedMobileOrDiffusive: {
      if (charged && (mobile || diffusive)) {
        solver->removeCoveredParticles(a_container, a_representation, a_tolerance);
      }

      break;
    }
    case SpeciesSubset::ChargedMobileAndDiffusive: {
      if (charged && (mobile && diffusive)) {
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
      MayDay::Abort("ItoKMCStepper::removeCoveredParticles - logic bust");

      break;
    }
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::transferCoveredParticles(const SpeciesSubset    a_speciesSubset,
                                                    const EBRepresentation a_representation,
                                                    const Real             a_tolerance) noexcept
{
  CH_TIME("ItoKMCStepper::transferCoveredParticles(SpeciesSubset, EBRepresentation, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::transferCoveredParticles(SpeciesSubset, EBRepresentation, Real)" << endl;
  }

  this->transferCoveredParticles(a_speciesSubset,
                                 ItoSolver::WhichContainer::Bulk,
                                 ItoSolver::WhichContainer::Covered,
                                 a_representation,
                                 a_tolerance);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::transferCoveredParticles(const SpeciesSubset             a_speciesSubset,
                                                    const ItoSolver::WhichContainer a_containerFrom,
                                                    const ItoSolver::WhichContainer a_containerTo,
                                                    const EBRepresentation          a_representation,
                                                    const Real                      a_tolerance) noexcept
{
  CH_TIME("ItoKMCStepper::transferCoveredParticles(SpeciesSubset, Containerx2, EBRepresentation, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::transferCoveredParticles(SpeciesSubset, Containerx2, EBRepresentation, Real)" << endl;
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
    case SpeciesSubset::AllMobileOrDiffusive: {
      if (mobile || diffusive) {
        solver->transferCoveredParticles(a_containerFrom, a_containerTo, a_representation, a_tolerance);
      }

      break;
    }
    case SpeciesSubset::AllMobileAndDiffusive: {
      if (mobile && diffusive) {
        solver->transferCoveredParticles(a_containerFrom, a_containerTo, a_representation, a_tolerance);
      }

      break;
    }
    case SpeciesSubset::Charged: {
      if (charged) {
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
    case SpeciesSubset::ChargedMobileOrDiffusive: {
      if (charged && (mobile || diffusive)) {
        solver->transferCoveredParticles(a_containerFrom, a_containerTo, a_representation, a_tolerance);
      }

      break;
    }
    case SpeciesSubset::ChargedMobileAndDiffusive: {
      if (charged && (mobile && diffusive)) {
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
      MayDay::Abort("ItoKMCStepper::transferCoveredParticles - logic bust");

      break;
    }
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::remapParticles(const SpeciesSubset a_speciesSubset) noexcept
{
  CH_TIME("ItoKMCStepper::remapParticles(SpeciesSubset)");
  if (m_verbosity > 5) {
    pout() << m_name + "::remapParticles(SpeciesSubset)" << endl;
  }

  this->remapParticles(a_speciesSubset, ItoSolver::WhichContainer::Bulk);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::remapParticles(const SpeciesSubset             a_speciesSubset,
                                          const ItoSolver::WhichContainer a_container) noexcept
{
  CH_TIME("ItoKMCStepper::remapParticles(SpeciesSubset, WhichContainer)");
  if (m_verbosity > 5) {
    pout() << m_name + "::remapParticles(SpeciesSubset, WhichContainer)" << endl;
  }

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

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
    case SpeciesSubset::AllMobileOrDiffusive: {
      if (mobile || diffusive) {
        solver->remap(a_container);
      }

      break;
    }
    case SpeciesSubset::AllMobileAndDiffusive: {
      if (mobile && diffusive) {
        solver->remap(a_container);
      }

      break;
    }
    case SpeciesSubset::Charged: {
      if (charged) {
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
    case SpeciesSubset::ChargedMobileOrDiffusive: {
      if (charged && (mobile || diffusive)) {
        solver->remap(a_container);
      }

      break;
    }
    case SpeciesSubset::ChargedMobileAndDiffusive: {
      if (charged && (mobile && diffusive)) {
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
      MayDay::Abort("ItoKMCStepper::remapParticles - logic bust");

      break;
    }
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::depositParticles(const SpeciesSubset a_speciesSubset) noexcept
{
  CH_TIME("ItoKMCStepper::depositParticles(SpeciesSubset)");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositParticles(SpeciesSubset)" << endl;
  }

  this->depositParticles(a_speciesSubset, ItoSolver::WhichContainer::Bulk);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::depositParticles(const SpeciesSubset             a_speciesSubset,
                                            const ItoSolver::WhichContainer a_container) noexcept
{
  CH_TIME("ItoKMCStepper::depositParticles(SpeciesSubset)");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositParticles(SpeciesSubset)" << endl;
  }

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

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
    case SpeciesSubset::AllMobileOrDiffusive: {
      if (mobile || diffusive) {
        solver->depositParticles(a_container);
      }

      break;
    }
    case SpeciesSubset::AllMobileAndDiffusive: {
      if (mobile && diffusive) {
        solver->depositParticles(a_container);
      }

      break;
    }
    case SpeciesSubset::Charged: {
      if (charged) {
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
    case SpeciesSubset::ChargedMobileOrDiffusive: {
      if (charged && (mobile || diffusive)) {
        solver->depositParticles(a_container);
      }

      break;
    }
    case SpeciesSubset::ChargedMobileAndDiffusive: {
      if (charged && (mobile && diffusive)) {
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
      MayDay::Abort("ItoKMCStepper::depositParticles - logic bust");

      break;
    }
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::setItoVelocityFunctions() noexcept
{
  CH_TIME("ItoKMCStepper::setItoVelocityFunctions");
  if (m_verbosity > 5) {
    pout() << m_name + "::setItoVelocityFunctions" << endl;
  }

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();
    const int                        Z       = species->getChargeNumber();

    if (solver->isMobile() && Z != 0) {
      EBAMRCellData& velocityFunction = solver->getVelocityFunction();
      velocityFunction.copy(m_electricFieldParticle);

      const int Z = species->getChargeNumber();

      if (Z < 0) {
        DataOps::scale(velocityFunction, -1.0);
      }

      // Coarsen and update ghost cells.
      m_amr->conservativeAverage(velocityFunction, m_particleRealm, m_plasmaPhase);
      m_amr->interpGhostPwl(velocityFunction, m_particleRealm, m_plasmaPhase);
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::setCdrVelocityFunctions() noexcept
{
  CH_TIME("ItoKMCStepper::setCdrVelocityFunctions");
  if (m_verbosity > 5) {
    pout() << m_name + "::setCdrVelocityFunctions" << endl;
  }

  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<CdrSolver>&        solver  = solverIt();
    const RefCountedPtr<CdrSpecies>& species = solver->getSpecies();
    const int                        Z       = species->getChargeNumber();

    if (solver->isMobile() && Z != 0) {
      EBAMRCellData& velocity = solver->getCellCenteredVelocity();
      velocity.copy(m_electricFieldFluid);

      const int Z = species->getChargeNumber();

      if (Z < 0) {
        DataOps::scale(velocity, -1.0);
      }

      // Coarsen and update ghost cells.
      m_amr->conservativeAverage(velocity, m_fluidRealm, m_plasmaPhase);
      m_amr->interpGhostPwl(velocity, m_fluidRealm, m_plasmaPhase);
    }
    else if (solver->isMobile() && Z == 0) {
      MayDay::Warning("ItoKMCStepper::setCdrVelocityFunctions -- how to handle mobile neutral species?");
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::multiplyCdrVelocitiesByMobilities() noexcept
{
  CH_TIME("ItoKMCStepper::multiplyCdrVelocitiesByMobilities()");
  if (m_verbosity > 5) {
    pout() << m_name + "::multiplyCdrVelocitiesByMobilities()" << endl;
  }

  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<CdrSolver>& solver = solverIt();
    const int                 idx    = solverIt.index();

    if (solver->isMobile()) {
      EBAMRCellData&       velocity = solver->getCellCenteredVelocity();
      const EBAMRCellData& mobility = m_cdrMobilities[idx];

      DataOps::multiplyScalar(velocity, mobility);

      // Coarsen and update ghost cells.
      m_amr->conservativeAverage(velocity, m_fluidRealm, m_plasmaPhase);
      m_amr->interpGhostPwl(velocity, m_fluidRealm, m_plasmaPhase);
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeDriftVelocities() noexcept
{
  CH_TIME("ItoKMCStepper::computeDriftVelocities()");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDriftVelocities()" << endl;
  }

  // Set velocities to be sgn(Z) * E
  this->setItoVelocityFunctions();
  this->setCdrVelocityFunctions();

  // Compute mobilities for both Ito and CDR species.
  this->computeMobilities();

  // Multiply sgn(Z) * E by the mobilities. For the CDR solvers this is just a multiplication, for the Ito solvers
  // we interpolate mu*E to the particle positions (in some form).
  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->interpolateVelocities();
  }

  this->multiplyCdrVelocitiesByMobilities();
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeMobilities() noexcept
{
  CH_TIME("ItoKMCStepper::computeMobilities()");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeMobilities()" << endl;
  }

  Vector<EBAMRCellData*> itoMobilities = m_ito->getMobilityFunctions();

  this->computeMobilities(itoMobilities, m_cdrMobilities, m_electricFieldFluid, m_time);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeMobilities(Vector<EBAMRCellData*>& a_itoMobilities,
                                             Vector<EBAMRCellData>&  a_cdrMobilities,
                                             const EBAMRCellData&    a_electricField,
                                             const Real              a_time) noexcept
{
  CH_TIME("ItoKMCStepper::computeMobilities(mobilities, E, time)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeMobilities(mobilities, E, time)" << endl;
  }

  const int numItoSpecies = m_physics->getNumItoSpecies();
  const int numCdrSpecies = m_physics->getNumCdrSpecies();

  CH_assert(a_electricField.getRealm() == m_fluidRealm);
  CH_assert(a_itoMobilities.size() == numItoSpecies);
  CH_assert(a_cdrMobilities.size() == numCdrSpecies);

  // The mesh mobilities belong on the particle realm (they are the ItoSolver mobilities) but we need to run
  // the computation on the fluid realm. So, create some transient storage for that.
  Vector<EBAMRCellData> fluidScratchMobilities(numItoSpecies);
  for (int i = 0; i < numItoSpecies; i++) {
    m_amr->allocate(fluidScratchMobilities[i], m_fluidRealm, m_plasmaPhase, 1);

    CH_assert(a_itoMobilities[i]->getRealm() == m_particleRealm);
  }

  // Now run the computation on the fluid realm, computing the mobilities into fluidScratchMobilities
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    Vector<LevelData<EBCellFAB>*> itoMobilities(numItoSpecies);
    Vector<LevelData<EBCellFAB>*> cdrMobilities(numCdrSpecies);

    for (int i = 0; i < numItoSpecies; i++) {
      itoMobilities[i] = &(*(fluidScratchMobilities[i])[lvl]);
    }

    for (int i = 0; i < numCdrSpecies; i++) {
      cdrMobilities[i] = &(*(a_cdrMobilities[i])[lvl]);
    }

    // Run the level computation, which will fill the mobilities.
    this->computeMobilities(itoMobilities, cdrMobilities, *a_electricField[lvl], lvl, a_time);
  }

  // Copy fluid realm data into particle realm and interpolate mobilities to the particle position.
  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>& solver = solverIt();

    if (solver->isMobile()) {
      const int idx = solverIt.index();

      a_itoMobilities[idx]->copy(fluidScratchMobilities[idx]);

      m_amr->conservativeAverage(*a_itoMobilities[idx], m_particleRealm, m_plasmaPhase);
      m_amr->interpGhostPwl(*a_itoMobilities[idx], m_particleRealm, m_plasmaPhase);

      solver->interpolateMobilities();
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeMobilities(Vector<LevelData<EBCellFAB>*>& a_itoMobilities,
                                             Vector<LevelData<EBCellFAB>*>& a_cdrMobilities,
                                             const LevelData<EBCellFAB>&    a_electricField,
                                             const int                      a_level,
                                             const Real                     a_time) noexcept
{
  CH_TIME("ItoKMCStepper::computeMobilities(mobilities, E, level, time)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeMobilities(mobilities, E, level, time)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_fluidRealm)[a_level];

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const EBCellFAB& E       = a_electricField[dit()];
    const Box        cellBox = dbl[dit()];

    Vector<EBCellFAB*> itoMobilities;
    Vector<EBCellFAB*> cdrMobilities;

    for (int i = 0; i < a_itoMobilities.size(); i++) {
      itoMobilities.push_back(&((*a_itoMobilities[i])[dit()]));
    }

    for (int i = 0; i < a_cdrMobilities.size(); i++) {
      cdrMobilities.push_back(&((*a_cdrMobilities[i])[dit()]));
    }

    this->computeMobilities(itoMobilities, cdrMobilities, E, a_level, dit(), cellBox, a_time);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeMobilities(Vector<EBCellFAB*>& a_itoMobilities,
                                             Vector<EBCellFAB*>& a_cdrMobilities,
                                             const EBCellFAB&    a_electricField,
                                             const int           a_level,
                                             const DataIndex     a_dit,
                                             const Box           a_box,
                                             const Real          a_time) noexcept
{
  CH_TIME("ItoKMCStepper::computeMobilities(meshMobilities, E, level, dit, box, time)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeMobilities(meshMobilities, E, level, dit, box, time)" << endl;
  }

  // TLDR: We go through each and every cell and call the physics interface. This includes cells covered by a finer grid
  //       but data is coarsened later anyways.
  const int numItoSpecies    = m_physics->getNumItoSpecies();
  const int numCdrSpecies    = m_physics->getNumCdrSpecies();
  const int numPlasmaSpecies = m_physics->getNumPlasmaSpecies();

  const Real     dx      = m_amr->getDx()[a_level];
  const RealVect probLo  = m_amr->getProbLo();
  const EBISBox& ebisbox = m_amr->getEBISLayout(m_fluidRealm, m_plasmaPhase)[a_level][a_dit];

  // Handles to regular data.
  const FArrayBox&   electricFieldReg = a_electricField.getFArrayBox();
  Vector<FArrayBox*> itoMobilitiesReg(numItoSpecies);
  Vector<FArrayBox*> cdrMobilitiesReg(numCdrSpecies);

  for (int i = 0; i < a_itoMobilities.size(); i++) {
    itoMobilitiesReg[i] = (&(a_itoMobilities[i]->getFArrayBox()));
  }

  for (int i = 0; i < a_cdrMobilities.size(); i++) {
    cdrMobilitiesReg[i] = (&(a_cdrMobilities[i]->getFArrayBox()));
  }

  // Physics interface mapping -- this maps a global index from the returned vector to an Ito or CDR solver
  const std::map<int, std::pair<SpeciesType, int>>& speciesMap = m_physics->getSpeciesMap();

  // Regular kernel
  auto regularKernel = [&](const IntVect& iv) -> void {
    const RealVect pos = m_amr->getProbLo() + dx * (RealVect(iv) + 0.5 * RealVect::Unit);
    const RealVect E   = RealVect(D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));

    // Call physics interface and compute mobilities for each species.
    const Vector<Real> mobilities = m_physics->computeMobilities(a_time, pos, E);

    CH_assert(mobilities.size() == numPlasmaSpecies);

    // Put the mobilities where they belong.
    for (const auto& s : speciesMap) {
      const int&         globalIndex = s.first;
      const SpeciesType& type        = s.second.first;
      const int&         localIndex  = s.second.second;

      if (type == SpeciesType::Ito) {
        (*itoMobilitiesReg[localIndex])(iv, 0) = mobilities[globalIndex];
      }
      else if (type == SpeciesType::CDR) {
        (*cdrMobilitiesReg[localIndex])(iv, 0) = mobilities[globalIndex];
      }
    }
  };

  // Irregular kernel.
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    const RealVect e   = RealVect(D_DECL(a_electricField(vof, 0), a_electricField(vof, 1), a_electricField(vof, 2)));
    const RealVect pos = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);

    // Call physics interface and compute mobilities for each species.
    const Vector<Real> mobilities = m_physics->computeMobilities(a_time, pos, e);

    CH_assert(mobilities.size() == numPlasmaSpecies);

    // Put the mobilities where they belong.
    for (const auto& s : speciesMap) {
      const int&         globalIndex = s.first;
      const SpeciesType& type        = s.second.first;
      const int&         localIndex  = s.second.second;

      if (type == SpeciesType::Ito) {
        (*a_itoMobilities[localIndex])(vof, 0) = mobilities[globalIndex];
      }
      else if (type == SpeciesType::CDR) {
        (*a_cdrMobilities[localIndex])(vof, 0) = mobilities[globalIndex];
      }
    }
  };

  VoFIterator& vofit = (*m_amr->getVofIterator(m_fluidRealm, m_plasmaPhase)[a_level])[a_dit];

  // Run the kernels.
  BoxLoops::loop(a_box, regularKernel);
  BoxLoops::loop(vofit, irregularKernel);

  // Covered is bogus.
  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    a_itoMobilities[solverIt.index()]->setCoveredCellVal(0.0, 0);
  }

  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    a_cdrMobilities[solverIt.index()]->setCoveredCellVal(0.0, 0);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeDiffusionCoefficients() noexcept
{
  CH_TIME("ItoKMCStepper::computeDiffusionCoefficients()");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDiffusionCoefficients()" << endl;
  }

  Vector<EBAMRCellData*> itoDiffusionCoefficients = m_ito->getDiffusionFunctions();
  Vector<EBAMRCellData*> cdrDiffusionCoefficients = m_cdr->getCellCenteredDiffusionCoefficients();

  this->computeDiffusionCoefficients(itoDiffusionCoefficients, cdrDiffusionCoefficients, m_electricFieldFluid, m_time);
  this->averageDiffusionCoefficientsCellToFace();
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeDiffusionCoefficients(Vector<EBAMRCellData*>& a_itoDiffusionCoefficients,
                                                        Vector<EBAMRCellData*>& a_cdrDiffusionCoefficients,
                                                        const EBAMRCellData&    a_electricField,
                                                        const Real              a_time) noexcept
{
  CH_TIME("ItoKMCStepper::computeDiffusionCoefficients(Vector<EBAMRCellData*>, EBAMRCellData, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDiffusionCoefficients(Vector<EBAMRCellData*>, EBAMRCellData, Real)" << endl;
  }

  const int numItoSpecies = m_physics->getNumItoSpecies();
  const int numCdrSpecies = m_physics->getNumCdrSpecies();

  CH_assert(a_electricField.getRealm() == m_fluidRealm);
  CH_assert(a_itoDiffusionCoefficients.size() == numItoSpecies);
  CH_assert(a_cdrDiffusionCoefficients.size() == numCdrSpecies);

  // Sanity check -- things need to be defined on the correct realms.
  for (int i = 0; i < numItoSpecies; i++) {
    CH_assert(a_itoDiffusionCoefficients[i]->getRealm() == m_particleRealm);
  }
  for (int i = 0; i < numCdrSpecies; i++) {
    CH_assert(a_cdrDiffusionCoefficients[i]->getRealm() == m_fluidRealm);
  }

  // The mesh diffusion coefficients belong on the particle realm (they are the ItoSolver diffusion coefficients) but we need to run
  // the computation on the fluid realm. So, create some transient storage for that.
  Vector<EBAMRCellData> fluidScratchDiffusion(numItoSpecies);
  for (int i = 0; i < numItoSpecies; i++) {
    m_amr->allocate(fluidScratchDiffusion[i], m_fluidRealm, m_plasmaPhase, 1);

    CH_assert(a_itoDiffusionCoefficients[i]->getRealm() == m_particleRealm);
  }

  // Compute mesh-based diffusion coefficients on the fluid realm.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    Vector<LevelData<EBCellFAB>*> itoDiffusionCoefficients(numItoSpecies);
    Vector<LevelData<EBCellFAB>*> cdrDiffusionCoefficients(numCdrSpecies);

    for (int i = 0; i < numItoSpecies; i++) {
      itoDiffusionCoefficients[i] = &(*(fluidScratchDiffusion[i])[lvl]);
    }
    for (int i = 0; i < numCdrSpecies; i++) {
      cdrDiffusionCoefficients[i] = &(*(*a_cdrDiffusionCoefficients[i])[lvl]);
    }

    this->computeDiffusionCoefficients(itoDiffusionCoefficients,
                                       cdrDiffusionCoefficients,
                                       *a_electricField[lvl],
                                       lvl,
                                       a_time);
  }

  // Copy the fluid realm data over to the particle realm data and then coarsen and interpolate diffusion coefficients
  // to particle positions.
  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>& solver = solverIt();

    if (solver->isDiffusive()) {
      const int idx = solverIt.index();

      a_itoDiffusionCoefficients[idx]->copy(fluidScratchDiffusion[idx]);

      m_amr->conservativeAverage(*a_itoDiffusionCoefficients[idx], m_particleRealm, m_plasmaPhase);
      m_amr->interpGhostPwl(*a_itoDiffusionCoefficients[idx], m_particleRealm, m_plasmaPhase);

      solver->interpolateDiffusion();
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeDiffusionCoefficients(Vector<LevelData<EBCellFAB>*>& a_itoDiffusionCoefficients,
                                                        Vector<LevelData<EBCellFAB>*>& a_cdrDiffusionCoefficients,
                                                        const LevelData<EBCellFAB>&    a_electricField,
                                                        const int                      a_level,
                                                        const Real                     a_time) noexcept
{
  CH_TIME("ItoKMCStepper::computeDiffusionCoefficients(Vector<LD<EBCellFAB>*>, LD<EBCellFAB>, int, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDiffusionCoefficients(Vector<LD<EBCellFAB>*>, LD<EBCellFAB>, int, Real)" << endl;
  }

  const int numItoSpecies = m_physics->getNumItoSpecies();
  const int numCdrSpecies = m_physics->getNumCdrSpecies();

  CH_assert(a_itoDiffusionCoefficients.size() == numItoSpecies);
  CH_assert(a_cdrDiffusionCoefficients.size() == numCdrSpecies);

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_fluidRealm)[a_level];

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    Vector<EBCellFAB*> itoDiffusionCoefficients(numItoSpecies);
    Vector<EBCellFAB*> cdrDiffusionCoefficients(numCdrSpecies);

    for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
      const int idx = solverIt.index();

      if (solverIt()->isDiffusive()) {
        itoDiffusionCoefficients[idx] = &(*a_itoDiffusionCoefficients[idx])[dit()];
      }
      else {
        itoDiffusionCoefficients[idx] = nullptr;
      }
    }

    for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
      const int idx = solverIt.index();

      if (solverIt()->isDiffusive()) {
        cdrDiffusionCoefficients[idx] = &(*a_cdrDiffusionCoefficients[idx])[dit()];
      }
      else {
        cdrDiffusionCoefficients[idx] = nullptr;
      }
    }

    this->computeDiffusionCoefficients(itoDiffusionCoefficients,
                                       cdrDiffusionCoefficients,
                                       a_electricField[dit()],
                                       a_level,
                                       dit(),
                                       dbl[dit()],
                                       a_time);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeDiffusionCoefficients(Vector<EBCellFAB*>& a_itoDiffusionCoefficients,
                                                        Vector<EBCellFAB*>& a_cdrDiffusionCoefficients,
                                                        const EBCellFAB&    a_electricField,
                                                        const int           a_level,
                                                        const DataIndex     a_dit,
                                                        const Box           a_box,
                                                        const Real          a_time) noexcept
{
  CH_TIME("ItoKMCStepper::computeDiffusionCoefficients(Patch)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDiffusionCoefficients(Patch)" << endl;
  }

  const int numItoSpecies    = m_physics->getNumItoSpecies();
  const int numCdrSpecies    = m_physics->getNumCdrSpecies();
  const int numPlasmaSpecies = m_physics->getNumCdrSpecies();

  CH_assert(a_electricField.nComp() == SpaceDim);
  CH_assert(a_itoDiffusionCoefficients.size() == numItoSpecies);
  CH_assert(a_cdrDiffusionCoefficients.size() == numPlasmaSpecies);

  // Geometric information that we need.
  const Real     dx      = m_amr->getDx()[a_level];
  const RealVect probLo  = m_amr->getProbLo();
  const EBISBox& ebisbox = m_amr->getEBISLayout(m_fluidRealm, m_plasmaPhase)[a_level][a_dit];

  // Handle to single-valued data.
  const FArrayBox& electricFieldReg = a_electricField.getFArrayBox();

  Vector<FArrayBox*> itoDiffCoReg(numItoSpecies, nullptr);
  Vector<FArrayBox*> cdrDiffCoReg(numCdrSpecies, nullptr);

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>& solver = solverIt();

    if (solver->isDiffusive()) {
      const int i     = solverIt.index();
      itoDiffCoReg[i] = &(a_itoDiffusionCoefficients[i]->getFArrayBox());
    }
  }

  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<CdrSolver>& solver = solverIt();

    if (solver->isDiffusive()) {
      const int i     = solverIt.index();
      cdrDiffCoReg[i] = &(a_cdrDiffusionCoefficients[i]->getFArrayBox());
    }
  }

  // Physics interface mapping -- this maps a global index from the returned vector to an Ito or CDR solver
  const std::map<int, std::pair<SpeciesType, int>>& speciesMap = m_physics->getSpeciesMap();

  // Regular kernel definition.
  auto regularKernel = [&](const IntVect& iv) -> void {
    const RealVect pos = probLo + dx * (RealVect(iv) + 0.5 * RealVect::Unit);
    const RealVect E   = RealVect(D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));

    // Compute diffusion coefficients.
    const Vector<Real> diffusionCoefficients = m_physics->computeDiffusionCoefficients(a_time, pos, E);

    CH_assert(diffusionCoefficients.size() == numPlasmaSpecies);

    // Put the diffusion coefficients in the correct solver storage.
    for (const auto& s : speciesMap) {
      const int&         globalIndex = s.first;
      const SpeciesType& type        = s.second.first;
      const int&         localIndex  = s.second.second;

      // We need an explicit check to see if
      if (type == SpeciesType::Ito) {
        if (m_ito->getSolvers()[localIndex]->isDiffusive()) {
          (*itoDiffCoReg[localIndex])(iv, 0) = diffusionCoefficients[globalIndex];
        }
      }
      else if (type == SpeciesType::CDR) {
        if (m_cdr->getSolvers()[localIndex]->isDiffusive()) {
          (*cdrDiffCoReg[localIndex])(iv, 0) = diffusionCoefficients[globalIndex];
        }
      }
    }
  };

  // Irregular kernel.
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    const RealVect E   = RealVect(D_DECL(a_electricField(vof, 0), a_electricField(vof, 1), a_electricField(vof, 2)));
    const RealVect pos = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);

    // Compute diffusion coefficients.
    const Vector<Real> diffusionCoefficients = m_physics->computeDiffusionCoefficients(a_time, pos, E);

    // Put the diffusion coefficients in the correct solver storage.
    for (const auto& s : speciesMap) {
      const int&         globalIndex = s.first;
      const SpeciesType& type        = s.second.first;
      const int&         localIndex  = s.second.second;

      // We need an explicit check to see if
      if (type == SpeciesType::Ito) {
        if (m_ito->getSolvers()[localIndex]->isDiffusive()) {
          (*a_itoDiffusionCoefficients[localIndex])(vof, 0) = diffusionCoefficients[globalIndex];
        }
      }
      else if (type == SpeciesType::CDR) {
        if (m_cdr->getSolvers()[localIndex]->isDiffusive()) {
          (*a_cdrDiffusionCoefficients[localIndex])(vof, 0) = diffusionCoefficients[globalIndex];
        }
      }
    }
  };

  // Run kernels.
  VoFIterator& vofit = (*m_amr->getVofIterator(m_fluidRealm, m_plasmaPhase)[a_level])[a_dit];

  BoxLoops::loop(a_box, regularKernel);
  BoxLoops::loop(vofit, irregularKernel);

  // Covered is bogus.
  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    if (solverIt()->isDiffusive()) {
      a_itoDiffusionCoefficients[solverIt.index()]->setCoveredCellVal(0.0, 0);
    }
  }

  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    if (solverIt()->isDiffusive()) {
      a_cdrDiffusionCoefficients[solverIt.index()]->setCoveredCellVal(0.0, 0);
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::averageDiffusionCoefficientsCellToFace() noexcept
{
  CH_TIME("ItoKMCStepper::averageDiffusionCoefficientsCellToFace");
  if (m_verbosity > 5) {
    pout() << m_name + "::averageDiffusionCoefficientsCellToFace" << endl;
  }

  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<CdrSolver>& solver = solverIt();

    if (solver->isDiffusive()) {

      EBAMRCellData& cellCenteredDiffusionCoefficient = solver->getCellCenteredDiffusionCoefficient();
      EBAMRFluxData& faceCenteredDiffusionCoefficient = solver->getFaceCenteredDiffusionCoefficient();

      CH_assert(cellCenteredDiffusionCoefficient.getRealm() == m_fluidRealm);
      CH_assert(faceCenteredDiffusionCoefficient.getRealm() == m_fluidRealm);

      DataOps::setValue(faceCenteredDiffusionCoefficient, std::numeric_limits<Real>::max());

      // Coarsen the cell-centered diffusion coefficient before averaging to faces.
      m_amr->arithmeticAverage(cellCenteredDiffusionCoefficient, m_fluidRealm, m_cdr->getPhase());
      m_amr->interpGhostPwl(cellCenteredDiffusionCoefficient, m_fluidRealm, m_cdr->getPhase());

      // Average to cell faces. Note that this call also includes one ghost face (required when there's an EBCF crossing).
      const int      tanGhost = 1;
      const Interval interv   = Interval(0, 0);
      const Average  average  = Average::Arithmetic;

      DataOps::averageCellToFace(faceCenteredDiffusionCoefficient,
                                 cellCenteredDiffusionCoefficient,
                                 m_amr->getDomains(),
                                 tanGhost,
                                 interv,
                                 interv,
                                 average);
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::getPhysicalParticlesPerCell(EBAMRCellData& a_ppc) const noexcept
{
  CH_TIME("ItoKMCStepper::getPhysicalParticlesPerCell(EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::getPhysicaParticlesPerCell(EBAMRCellData)" << endl;
  }

  CH_assert(a_ppc.getRealm() == m_particleRealm);

  for (auto it = m_ito->iterator(); it.ok(); ++it) {
    const int idx = it.index();

    EBAMRCellData ppc = m_amr->slice(a_ppc, Interval(idx, idx));

    const ParticleContainer<ItoParticle>& particles = it()->getParticles(ItoSolver::WhichContainer::Bulk);

    ParticleOps::getPhysicalParticlesPerCell<ItoParticle, &ItoParticle::weight>(ppc, particles);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeReactiveItoParticlesPerCell(EBAMRCellData& a_ppc) noexcept
{
  CH_TIME("ItoKMCStepper::computeReactiveItoParticlesPerCell(EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeReactiveItoParticlesPerCell(EBAMRCellData)" << endl;
  }

  CH_assert(a_ppc.getRealm() == m_particleRealm);

  DataOps::setValue(a_ppc, 0.0);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    this->computeReactiveItoParticlesPerCell(*a_ppc[lvl], lvl);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeReactiveItoParticlesPerCell(LevelData<EBCellFAB>& a_ppc, const int a_level) noexcept
{
  CH_TIME("ItoKMCStepper::computeReactiveItoParticlesPerCell(LD<EBCellFAB>, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeReactiveItoParticlesPerCell(LD<EBCellFAB>, int)" << endl;
  }

  const int numItoSpecies = m_physics->getNumItoSpecies();

  CH_assert(a_ppc.nComp() == numItoSpecies);

  const DisjointBoxLayout& dbl   = m_amr->getGrids(m_particleRealm)[a_level];
  const EBISLayout&        ebisl = m_amr->getEBISLayout(m_particleRealm, m_plasmaPhase)[a_level];

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box      box     = dbl[dit()];
    const EBISBox& ebisbox = ebisl[dit()];

    this->computeReactiveItoParticlesPerCell(a_ppc[dit()], a_level, dit(), box, ebisbox);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeReactiveItoParticlesPerCell(EBCellFAB&      a_ppc,
                                                              const int       a_level,
                                                              const DataIndex a_dit,
                                                              const Box       a_box,
                                                              const EBISBox&  a_ebisbox) noexcept
{
  CH_TIME("ItoKMCStepper::computeReactiveItoParticlesPerCell(EBCellFAB, int, DataIndex, Box, EBISBox)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeReactiveItoParticlesPerCell(EBCellFAB, int, DataIndex, Box, EBISBox)" << endl;
  }

  const int numItoSpecies = m_physics->getNumItoSpecies();

  CH_assert(a_ppc.nComp() == numItoSpecies);

  const Real     dx     = m_amr->getDx()[a_level];
  const RealVect probLo = m_amr->getProbLo();

  // TLDR: We go through each solver and add the number of PHYSICAL particles per cell to a_ppc.

  FArrayBox& ppcRegular = a_ppc.getFArrayBox();

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>& solver = solverIt();
    const int                 idx    = solverIt.index();

    // Get the cell-sorted particles. This will issue an error if the user has not
    // sorted the particles by cell before calling this routine.
    const ParticleContainer<ItoParticle>& particles     = solver->getParticles(ItoSolver::WhichContainer::Bulk);
    const BinFab<ItoParticle>&            cellParticles = particles.getCellParticles(a_level, a_dit);

    // Regular cells kernel.
    auto regularKernel = [&](const IntVect& iv) -> void {
      Real num = 0.0;

      if (a_ebisbox.isRegular(iv)) {
        for (ListIterator<ItoParticle> lit(cellParticles(iv, 0)); lit.ok(); ++lit) {
          num += lit().weight();
        }
      }

      ppcRegular(iv, idx) = num;
    };

    // Irregular kernel -- note that only particles that lie inside the domain get to react.
    auto irregularKernel = [&](const VolIndex& vof) -> void {
      const IntVect  iv           = vof.gridIndex();
      const RealVect normal       = a_ebisbox.normal(vof);
      const RealVect physCentroid = probLo + Location::position(Location::Cell::Boundary, vof, a_ebisbox, dx);

      Real num = 0.0;

      for (ListIterator<ItoParticle> lit(cellParticles(iv, 0)); lit.ok(); ++lit) {
        const RealVect& pos = lit().position();
        if ((pos - physCentroid).dotProduct(normal) >= 0.0) {
          num += lit().weight();
        }
      }

      ppcRegular(iv, idx) = num;
    };

    // Run the kernels.
    VoFIterator& vofit = (*m_amr->getVofIterator(m_particleRealm, m_plasmaPhase)[a_level])[a_dit];

    BoxLoops::loop(a_box, regularKernel);
    BoxLoops::loop(vofit, irregularKernel);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeReactiveCdrParticlesPerCell(EBAMRCellData& a_ppc) noexcept
{
  CH_TIME("ItoKMCStepper::computeReactiveCdrParticlesPerCell(EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeReactiveCdrParticlesPerCell(EBAMRCellData)" << endl;
  }

  CH_assert(a_ppc.getRealm() == m_fluidRealm);

  DataOps::setValue(a_ppc, 0.0);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    this->computeReactiveCdrParticlesPerCell(*a_ppc[lvl], lvl);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeReactiveCdrParticlesPerCell(LevelData<EBCellFAB>& a_ppc, const int a_level) noexcept
{
  CH_TIME("ItoKMCStepper::computeReactiveCdrParticlesPerCell(LD<EBCellFAB>, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeReactiveCdrParticlesPerCell(LD<EBCellFAB>, int)" << endl;
  }

  const int numCdrSpecies = m_physics->getNumCdrSpecies();

  CH_assert(a_ppc.nComp() == numCdrSpecies);

  if (numCdrSpecies > 0) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_fluidRealm)[a_level];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_fluidRealm, m_plasmaPhase)[a_level];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box      box     = dbl[dit()];
      const EBISBox& ebisbox = ebisl[dit()];

      this->computeReactiveCdrParticlesPerCell(a_ppc[dit()], a_level, dit(), box, ebisbox);
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeReactiveCdrParticlesPerCell(EBCellFAB&      a_ppc,
                                                              const int       a_level,
                                                              const DataIndex a_dit,
                                                              const Box       a_box,
                                                              const EBISBox&  a_ebisbox) noexcept
{
  CH_TIME("ItoKMCStepper::computeReactiveCdrParticlesPerCell(EBCellFAB, int, DataIndex, Box, EBISBox)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeReactiveCdrParticlesPerCell(EBCellFAB, int, DataIndex, Box, EBISBox)" << endl;
  }

  constexpr Real zero = 0.0;

  const int numCdrSpecies = m_physics->getNumCdrSpecies();

  CH_assert(a_ppc.nComp() == numCdrSpecies);

  const Real     dx     = m_amr->getDx()[a_level];
  const Real     vol    = std::pow(dx, SpaceDim);
  const RealVect probLo = m_amr->getProbLo();

  // TLDR: We go through each solver and add the number of PHYSICAL particles per cell to a_ppc.
  FArrayBox& ppcRegular = a_ppc.getFArrayBox();

  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<CdrSolver>& solver = solverIt();
    const int                 idx    = solverIt.index();

    const EBCellFAB& phi    = (*(solver->getPhi())[a_level])[a_dit];
    const FArrayBox& phiReg = phi.getFArrayBox();

    // Regular cells kernel.
    auto regularKernel = [&](const IntVect& iv) -> void {
      if (a_ebisbox.isRegular(iv)) {
        ppcRegular(iv, idx) = std::max(zero, std::floor(phiReg(iv, 0) * vol));
      }
    };

    // Irregular kernel -- note that only particles that lie inside the domain get to react.
    auto irregularKernel = [&](const VolIndex& vof) -> void {
      const Real kappa = a_ebisbox.volFrac(vof);

      a_ppc(vof, idx) = std::max(zero, std::floor(kappa * phi(vof, 0) * vol));
    };

    // Run the kernels.
    VoFIterator& vofit = (*m_amr->getVofIterator(m_fluidRealm, m_plasmaPhase)[a_level])[a_dit];

    BoxLoops::loop(a_box, regularKernel);
    BoxLoops::loop(vofit, irregularKernel);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeReactiveMeanEnergiesPerCell(EBAMRCellData& a_meanEnergies) noexcept
{
  CH_TIME("ItoKMCStepper::computeReactiveMaeanEnergiesPerCell(EBAMRCellData)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeReactiveMaeanEnergiesPerCell(EBAMRCellData)" << endl;
  }

  CH_assert(a_meanEnergies.getRealm() == m_particleRealm);

  DataOps::setValue(a_meanEnergies, 0.0);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    this->computeReactiveMeanEnergiesPerCell(*a_meanEnergies[lvl], lvl);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeReactiveMeanEnergiesPerCell(LevelData<EBCellFAB>& a_meanEnergies,
                                                              const int             a_level) noexcept
{
  CH_TIME("ItoKMCStepper::computeReactiveMeanEnergiesPerCell(LD<EBCellFAB>, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeReactiveMeanEnergiesPerCell(LD<EBCellFAB>, int)" << endl;
  }

  const int numPlasmaSpecies = m_physics->getNumItoSpecies();

  CH_assert(a_meanEnergies.nComp() == numPlasmaSpecies);

  const DisjointBoxLayout& dbl   = m_amr->getGrids(m_particleRealm)[a_level];
  const EBISLayout&        ebisl = m_amr->getEBISLayout(m_particleRealm, m_plasmaPhase)[a_level];

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box      box     = dbl[dit()];
    const EBISBox& ebisbox = ebisl[dit()];

    this->computeReactiveMeanEnergiesPerCell(a_meanEnergies[dit()], a_level, dit(), box, ebisbox);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeReactiveMeanEnergiesPerCell(EBCellFAB&      a_meanEnergies,
                                                              const int       a_level,
                                                              const DataIndex a_dit,
                                                              const Box       a_box,
                                                              const EBISBox&  a_ebisbox) noexcept
{
  CH_TIME("ItoKMCStepper::computeReactiveMeanEnergiesPerCell(EBCellFABint, DataIndex, Box, EBISBox)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeReactiveMeanEnergiesPerCell(EBCellFABint, DataIndex, Box, EBISBox))" << endl;
  }

  const int numPlasmaSpecies = m_physics->getNumItoSpecies();

  CH_assert(a_meanEnergies.nComp() == numPlasmaSpecies);

  const Real     dx     = m_amr->getDx()[a_level];
  const RealVect probLo = m_amr->getProbLo();

  // Get single-valued data.
  FArrayBox& meanEnergiesReg = a_meanEnergies.getFArrayBox();

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>& solver = solverIt();
    const int                 idx    = solverIt.index();

    // Get the cell-sorted particles. This will issue an error if the user has not
    // sorted the particles by cell before calling this routine.
    const ParticleContainer<ItoParticle>& particles     = solver->getParticles(ItoSolver::WhichContainer::Bulk);
    const BinFab<ItoParticle>&            cellParticles = particles.getCellParticles(a_level, a_dit);

    // Regular grid cells.
    auto regularKernel = [&](const IntVect& iv) -> void {
      if (a_ebisbox.isRegular(iv)) {
        Real totalWeight = 0.0;
        Real totalEnergy = 0.0;

        for (ListIterator<ItoParticle> lit(cellParticles(iv, 0)); lit.ok(); ++lit) {
          totalWeight += lit().weight();
          totalEnergy += lit().weight() * lit().energy();
        }

        if (totalWeight > 0.0) {
          meanEnergiesReg(iv, idx) = totalEnergy / totalWeight;
        }
        else {
          meanEnergiesReg(iv, idx) = 0.0;
        }
      }
    };

    // Irregular cells -- note that only valid particles get to play with us.
    auto irregularKernel = [&](const VolIndex& vof) -> void {
      const IntVect  iv         = vof.gridIndex();
      const RealVect normal     = a_ebisbox.normal(vof);
      const RealVect ebCentroid = probLo + Location::position(Location::Cell::Boundary, vof, a_ebisbox, dx);

      Real totalWeight = 0.0;
      Real totalEnergy = 0.0;

      for (ListIterator<ItoParticle> lit(cellParticles(iv, 0)); lit.ok(); ++lit) {
        const RealVect& pos = lit().position();

        if ((pos - ebCentroid).dotProduct(normal) >= 0.0) {
          totalWeight += lit().weight();
          totalEnergy += lit().weight() * lit().energy();
        }
      }

      if (totalWeight > 0.0) {
        meanEnergiesReg(iv, idx) = totalEnergy / totalWeight;
      }
      else {
        meanEnergiesReg(iv, idx) = 0.0;
      }
    };

    // Run the kernels.
    VoFIterator& vofit = (*m_amr->getVofIterator(m_particleRealm, m_plasmaPhase)[a_level])[a_dit];

    BoxLoops::loop(a_box, regularKernel);
    BoxLoops::loop(vofit, irregularKernel);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::advanceReactionNetwork(const Real a_dt) noexcept
{
  CH_TIME("ItoKMCStepper::advanceReactionNetwork(dt)");
  if (m_verbosity > 5) {
    pout() << m_name + "::advanceReactionNetwork(dt)" << endl;
  }

  CH_assert(a_dt > 0.0);

  this->advanceReactionNetwork(m_electricFieldFluid, a_dt);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::advanceReactionNetwork(const EBAMRCellData& a_electricField, const Real a_dt) noexcept
{
  CH_TIME("ItoKMCStepper::advanceReactionNetwork(EBAMRCellData, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::advanceReactionNetwork(EBAMRCellData, Real)" << endl;
  }

  const int numItoSpecies    = m_physics->getNumItoSpecies();
  const int numCdrSpecies    = m_physics->getNumCdrSpecies();
  const int numPhotonSpecies = m_physics->getNumPhotonSpecies();

  CH_assert(a_electricField.getRealm() == m_fluidRealm);
  CH_assert(a_dt > 0.0);

  // Compute the number of reactive particles for both Ito and CDR species. Also do a backup of the initial number
  // of particles per cell. This is required when reconciling the results lateron.
  if (numItoSpecies > 0) {
    this->computeReactiveItoParticlesPerCell(m_particleItoPPC);

    const Interval srcInterv(0, numItoSpecies - 1);
    const Interval dstInterv(0, numItoSpecies - 1);

    m_fluidPPC.copy(srcInterv, m_particleItoPPC, dstInterv);

    DataOps::copy(m_particleOldItoPPC, m_particleItoPPC);
  }
  if (numCdrSpecies > 0) {
    this->computeReactiveCdrParticlesPerCell(m_fluidCdrPPC);

    const Interval srcInterv(0, numCdrSpecies - 1);
    const Interval dstInterv(numItoSpecies, numItoSpecies + numCdrSpecies - 1);

    m_fluidPPC.copy(srcInterv, m_fluidCdrPPC, dstInterv);

    DataOps::copy(m_fluidOldCdrPPC, m_fluidCdrPPC);
  }

  DataOps::setValue(m_fluidYPC, 0.0);
  DataOps::setValue(m_particleYPC, 0.0);

  // Advance the reaction network which gives us a new number of particles per cell, as well as the number of
  // photons that need to be generated per cell.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    this->advanceReactionNetwork(*m_fluidPPC[lvl], *m_fluidYPC[lvl], *a_electricField[lvl], lvl, a_dt);
  }

  // Copy the results back to the holders that hold the number of particles per cell for Ito/Cdr solvers.
  if (numItoSpecies > 0) {
    const Interval srcInterv(0, numItoSpecies - 1);
    const Interval dstInterv(0, numItoSpecies - 1);

    m_particleItoPPC.copy(srcInterv, m_fluidPPC, dstInterv);
  }
  if (numCdrSpecies > 0) {
    const Interval srcInterv(numItoSpecies, numItoSpecies + numCdrSpecies - 1);
    const Interval dstInterv(0, numCdrSpecies - 1);

    m_fluidCdrPPC.copy(srcInterv, m_fluidPPC, dstInterv);
  }
  if (numPhotonSpecies > 0) {
    m_particleYPC.copy(m_fluidYPC);
  }

  // If we have photoionization reactions and CDR solvers we need to put the resulting particles on the mesh first,
  // and then add the result back into the CDR solvers.
  for (int i = 0; i < m_cdrPhotoiProducts.size(); i++) {
    m_cdrPhotoiProducts[i].clearParticles();
    m_cdrPhotoiProducts[i].organizeParticlesByCell();
  }

  // Reconcile the results -- for the discrete solvers we add/remove particles/photons and for the CDR solvers we update
  // the source terms.
  this->reconcileParticles(m_particleItoPPC, m_particleOldItoPPC, m_particleYPC);

  // Deposited PointParticles on the mesh and add the result to m_fluidCdrPPC. This is a bit convoluted because
  // the photoionization products are created in reconcileParticles and  deposited them on the particle realm, but
  // the resulting products need to end up on in the correct data component in m_fluidCdrPPC.
  for (int i = 0; i < m_cdrPhotoiProducts.size(); i++) {
    m_cdrPhotoiProducts[i].organizeParticlesByPatch();

    m_amr->depositParticles<PointParticle, &PointParticle::weight>(m_particleScratch1,
                                                                   m_particleRealm,
                                                                   m_plasmaPhase,
                                                                   m_cdrPhotoiProducts[i],
                                                                   DepositionType::NGP,
                                                                   CoarseFineDeposition::Halo);

    m_fluidScratch1.copy(m_particleScratch1);
    DataOps::volumeScale(m_fluidScratch1, m_amr->getDx());

    EBAMRCellData fluidCdrPPC = m_amr->slice(m_fluidCdrPPC, Interval(i, i));

    DataOps::incr(fluidCdrPPC, m_fluidScratch1, 1.0);

    m_cdrPhotoiProducts[i].clearParticles();
  }

  this->reconcileCdrDensities(m_fluidCdrPPC, m_fluidOldCdrPPC, a_dt);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::advanceReactionNetwork(LevelData<EBCellFAB>&       a_particlesPerCell,
                                                  LevelData<EBCellFAB>&       a_newPhotonsPerCell,
                                                  const LevelData<EBCellFAB>& a_electricField,
                                                  const int                   a_level,
                                                  const Real                  a_dt) noexcept
{
  CH_TIME("ItoKMCStepper::advanceReactionNetwork(LD<EBCellFAB>x3, int, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::advanceReactionNetwork(LD<EBCellFAB>x3, int, Real)" << endl;
  }

  const int numPlasmaSpecies = m_physics->getNumPlasmaSpecies();
  const int numPhotonSpecies = m_physics->getNumPhotonSpecies();

  CH_assert(a_particlesPerCell.nComp() == numPlasmaSpecies);
  CH_assert(a_newPhotonsPerCell.nComp() == numPhotonSpecies);
  CH_assert(a_electricField.nComp() == SpaceDim);

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_fluidRealm)[a_level];

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    this->advanceReactionNetwork(a_particlesPerCell[dit()],
                                 a_newPhotonsPerCell[dit()],
                                 a_electricField[dit()],
                                 a_level,
                                 dit(),
                                 dbl[dit()],
                                 m_amr->getDx()[a_level],
                                 a_dt);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::advanceReactionNetwork(EBCellFAB&       a_particlesPerCell,
                                                  EBCellFAB&       a_newPhotonsPerCell,
                                                  const EBCellFAB& a_electricField,
                                                  const int        a_level,
                                                  const DataIndex  a_dit,
                                                  const Box        a_box,
                                                  const Real       a_dx,
                                                  const Real       a_dt) noexcept
{
  CH_TIME("ItoKMCStepper::advanceReactionNetwork(EBCellFABx3, int, DataIndex, Box, Realx2)");
  if (m_verbosity > 5) {
    pout() << m_name + "::advanceReactionNetwork(EBCellFABx3, int, DataIndex, Box, Realx2)" << endl;
  }

  const int numCdrSpecies    = m_physics->getNumCdrSpecies();
  const int numItoSpecies    = m_physics->getNumItoSpecies();
  const int numPlasmaSpecies = m_physics->getNumPlasmaSpecies();
  const int numPhotonSpecies = m_physics->getNumPhotonSpecies();

  CH_assert(a_particlesPerCell.nComp() == numPlasmaSpecies);
  CH_assert(a_newPhotonsPerCell.nComp() == numPhotonSpecies);
  CH_assert(a_electricField.nComp() == SpaceDim);

  // Geometric information that we require.
  const RealVect probLo  = m_amr->getProbLo();
  const EBISBox& ebisbox = m_amr->getEBISLayout(m_fluidRealm, m_plasmaPhase)[a_level][a_dit];

  const FArrayBox& electricFieldReg = a_electricField.getFArrayBox();

  // Storage used by physics interface.
  Vector<Physics::ItoKMC::FPR> particles(numPlasmaSpecies);
  Vector<Physics::ItoKMC::FPR> newPhotons(numPhotonSpecies);
  Vector<Real>                 meanEnergies(numPlasmaSpecies);
  Vector<Real>                 energySources(numPlasmaSpecies);
  Vector<RealVect>             densityGradients(numPlasmaSpecies, RealVect::Zero);

  // Populate single-valued data.
  FArrayBox& particlesPerCellReg = a_particlesPerCell.getFArrayBox();
  FArrayBox& newPhotonsReg       = a_newPhotonsPerCell.getFArrayBox();

  // Handle to density gradients for CDR and Ito species.
  Vector<EBCellFAB*> densityGradientsCDR(numCdrSpecies);
  Vector<EBCellFAB*> densityGradientsIto(numItoSpecies);

  Vector<FArrayBox*> densityGradientsCDRReg(numCdrSpecies);
  Vector<FArrayBox*> densityGradientsItoReg(numItoSpecies);

  for (int i = 0; i < numCdrSpecies; i++) {
    densityGradientsCDR[i]    = &(*m_fluidGradPhiCDR[i][a_level])[a_dit];
    densityGradientsCDRReg[i] = &(densityGradientsCDR[i]->getFArrayBox());
  }

  for (int i = 0; i < numItoSpecies; i++) {
    densityGradientsIto[i]    = &(*m_fluidGradPhiIto[i][a_level])[a_dit];
    densityGradientsItoReg[i] = &(densityGradientsIto[i]->getFArrayBox());
  }

  // Handle to valid grid cells.
  const BaseFab<bool>& validCells = (*m_amr->getValidCells(m_fluidRealm)[a_level])[a_dit];

  // Regular cells
  auto regularKernel = [&](const IntVect& iv) -> void {
    if (ebisbox.isRegular(iv) && validCells(iv, 0)) {
      const RealVect pos = probLo + a_dx * (RealVect(iv) + 0.5 * RealVect::Unit);
      const RealVect E   = RealVect(D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));

      // Populate the data holders that the physics interface requires.
      for (int i = 0; i < numPlasmaSpecies; i++) {
        particles[i] = (long long)particlesPerCellReg(iv, i);
      }

      for (int i = 0; i < numPhotonSpecies; i++) {
        newPhotons[i] = 0LL;
      }

      // Populate gradients.
      for (int i = 0; i < numItoSpecies; i++) {
        densityGradients[i] = RealVect(D_DECL((*densityGradientsItoReg[i])(iv, 0),
                                              (*densityGradientsItoReg[i])(iv, 1),
                                              (*densityGradientsItoReg[i])(iv, 2)));
      }
      for (int i = 0; i < numCdrSpecies; i++) {
        densityGradients[numItoSpecies + i] = RealVect(D_DECL((*densityGradientsCDRReg[i])(iv, 0),
                                                              (*densityGradientsCDRReg[i])(iv, 1),
                                                              (*densityGradientsCDRReg[i])(iv, 2)));
      }

      // Do the physics advance.
      m_physics->advanceKMC(particles, newPhotons, a_dt, E, a_dx, 1.0);

      // Repopulate the input data holders with the new number of particles/photons per cell.
      for (int i = 0; i < numPlasmaSpecies; i++) {
        particlesPerCellReg(iv, i) = 1.0 * particles[i];
      }

      for (int i = 0; i < numPhotonSpecies; i++) {
        newPhotonsReg(iv, i) = 1.0 * newPhotons[i];
      }
    }
  };

  // Irregular cells
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    const IntVect iv = vof.gridIndex();

    if (validCells(iv, 0)) {
      const Real     kappa = ebisbox.volFrac(vof);
      const RealVect pos   = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, a_dx);
      const RealVect E = RealVect(D_DECL(a_electricField(vof, 0), a_electricField(vof, 1), a_electricField(vof, 2)));

      // Initialize for this cell.
      for (int i = 0; i < numPlasmaSpecies; i++) {
        particles[i] = (long long)a_particlesPerCell(vof, i);
      }

      for (int i = 0; i < numPhotonSpecies; i++) {
        newPhotons[i] = 0LL;
      }

      // Populate gradients.
      for (int i = 0; i < numItoSpecies; i++) {
        densityGradients[i] = RealVect(D_DECL((*densityGradientsIto[i])(vof, 0),
                                              (*densityGradientsIto[i])(vof, 1),
                                              (*densityGradientsIto[i])(vof, 2)));
      }
      for (int i = 0; i < numCdrSpecies; i++) {
        densityGradients[numItoSpecies + i] = RealVect(D_DECL((*densityGradientsCDR[i])(vof, 0),
                                                              (*densityGradientsCDR[i])(vof, 1),
                                                              (*densityGradientsCDR[i])(vof, 2)));
      }

      // Do the physics advance
      m_physics->advanceKMC(particles, newPhotons, a_dt, E, a_dx, kappa);

      // Repopulate the input data holders with the new number of particles/photons per cell.
      for (int i = 0; i < numPlasmaSpecies; i++) {
        a_particlesPerCell(vof, i) = 1.0 * particles[i];
      }

      for (int i = 0; i < numPhotonSpecies; i++) {
        a_newPhotonsPerCell(vof, i) = 1.0 * newPhotons[i];
      }
    }
  };

  // Run the kernels.
  VoFIterator& vofit = (*m_amr->getVofIterator(m_fluidRealm, m_plasmaPhase)[a_level])[a_dit];

  BoxLoops::loop(a_box, regularKernel);
  BoxLoops::loop(vofit, irregularKernel);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::reconcileParticles(const EBAMRCellData& a_newParticlesPerCell,
                                              const EBAMRCellData& a_oldParticlesPerCell,
                                              const EBAMRCellData& a_newPhotonsPerCell) noexcept
{
  CH_TIME("ItoKMCStepper::reconcileParticles(EBAMRCellDatax3)");
  if (m_verbosity > 5) {
    pout() << m_name + "::reconcileParticles(EBAMRCellDatax3)";
  }

  CH_assert(a_newParticlesPerCell.getRealm() == m_particleRealm);
  CH_assert(a_oldParticlesPerCell.getRealm() == m_particleRealm);
  CH_assert(a_newPhotonsPerCell.getRealm() == m_particleRealm);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    this->reconcileParticles(*a_newParticlesPerCell[lvl], *a_oldParticlesPerCell[lvl], *a_newPhotonsPerCell[lvl], lvl);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::reconcileParticles(const LevelData<EBCellFAB>& a_newParticlesPerCell,
                                              const LevelData<EBCellFAB>& a_oldParticlesPerCell,
                                              const LevelData<EBCellFAB>& a_newPhotonsPerCell,
                                              const int                   a_level) noexcept
{
  CH_TIME("ItoKMCStepper::reconcileParticles(LevelData<EBCellFAB>x3, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::reconcileParticles(LevelData<EBCellFAB>x3, int)" << endl;
  }

  const int numItoSpecies    = m_physics->getNumItoSpecies();
  const int numPhotonSpecies = m_physics->getNumPhotonSpecies();

  CH_assert(a_newParticlesPerCell.nComp() == numItoSpecies);
  CH_assert(a_oldParticlesPerCell.nComp() == numItoSpecies);
  CH_assert(a_newPhotonsPerCell.nComp() == numPhotonSpecies);

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[a_level];

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    this->reconcileParticles(a_newParticlesPerCell[dit()],
                             a_oldParticlesPerCell[dit()],
                             a_newPhotonsPerCell[dit()],
                             a_level,
                             dit(),
                             dbl[dit()],
                             m_amr->getDx()[a_level]);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::reconcileParticles(const EBCellFAB& a_newParticlesPerCell,
                                              const EBCellFAB& a_oldParticlesPerCell,
                                              const EBCellFAB& a_newPhotonsPerCell,
                                              const int        a_level,
                                              const DataIndex  a_dit,
                                              const Box        a_box,
                                              const Real       a_dx) noexcept
{
  CH_TIME("ItoKMCStepper::reconcileParticles(EBCellFABx3, int, DataIndex, Box, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::reconcileParticles(EBCellFABx3, int, DataIndex, Box, Real)" << endl;
  }

  // TLDR: This is the main routine for generating new particles/photons after the chemistry advance have finished. We have
  //       already computed the number of particles in each grid cell, and we now need to generate them. To do that we use
  //       the reconciliation routines from the physics interface, which takes the per-cell responsibility for that. The main
  //       work done in this routine is to expose the per-patch data to per-cell data that the physics interface can then use.

  const int numItoSpecies    = m_physics->getNumItoSpecies();
  const int numCdrSpecies    = m_physics->getNumCdrSpecies();
  const int numPhotonSpecies = m_physics->getNumPhotonSpecies();

  CH_assert(a_newParticlesPerCell.nComp() == numItoSpecies);
  CH_assert(a_oldParticlesPerCell.nComp() == numItoSpecies);
  CH_assert(a_newPhotonsPerCell.nComp() == numPhotonSpecies);

  // Geometric information that we need.
  const RealVect probLo  = m_amr->getProbLo();
  const EBISBox& ebisbox = m_amr->getEBISLayout(m_particleRealm, m_plasmaPhase)[a_level][a_dit];

  // List of valid grid cells
  const BaseFab<bool>& validCells = (*m_amr->getValidCells(m_particleRealm)[a_level])[a_dit];

  // Pointers to cell-centered particle storages. Note that the BinFabs hold a list of particles/photons
  // on each grid cell. We will manipulate these lists.

  Vector<BinFab<ItoParticle>*>   itoParticlesFAB(numItoSpecies);
  Vector<BinFab<PointParticle>*> cdrParticlesFAB(numCdrSpecies);
  Vector<BinFab<Photon>*>        sourcePhotonsFAB(numPhotonSpecies);
  Vector<BinFab<Photon>*>        bulkPhotonsFAB(numPhotonSpecies);

  // Get a handle to the cell-sorted particles from the Ito solvers.
  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>& solver = solverIt();

    const int idx = solverIt.index();

    ParticleContainer<ItoParticle>& solverParticles = solver->getParticles(ItoSolver::WhichContainer::Bulk);

    itoParticlesFAB[idx] = &(solverParticles.getCellParticles(a_level, a_dit));
  }

  // Get a handle to the cell-sorted particles for CDR photoionization products
  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    const int idx = solverIt.index();

    ParticleContainer<PointParticle>& pointParticles = m_cdrPhotoiProducts[idx];

    cdrParticlesFAB[idx] = &(pointParticles.getCellParticles(a_level, a_dit));
  }

  // Get a handle to the cell-sorted photons from the Monte Carlo photon solvers.
  for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<McPhoto>& solver = solverIt();

    const int idx = solverIt.index();

    ParticleContainer<Photon>& solverBulkPhotons   = solver->getBulkPhotons();
    ParticleContainer<Photon>& solverSourcePhotons = solver->getSourcePhotons();

    bulkPhotonsFAB[idx]   = &(solverBulkPhotons.getCellParticles(a_level, a_dit));
    sourcePhotonsFAB[idx] = &(solverSourcePhotons.getCellParticles(a_level, a_dit));
  }

  // The physics interface takes the physical number of particles/photons as arguments
  // to the reconciliation routines. These need to be set from the input arguments; this is the
  // storage we use in the grid cells.
  Vector<Physics::ItoKMC::FPR> numNewParticles(numItoSpecies);
  Vector<Physics::ItoKMC::FPR> numOldParticles(numItoSpecies);
  Vector<Physics::ItoKMC::FPR> numNewPhotons(numPhotonSpecies);

  // The physics interface also takes the actual particles/photons as argument to its reconciliation routines. This
  // is the storage we use for these; note that it is repopulated in every grid cell.
  Vector<List<ItoParticle>*>   itoParticles(numItoSpecies);
  Vector<List<PointParticle>*> cdrParticles(numCdrSpecies);
  Vector<List<Photon>*>        bulkPhotons(numPhotonSpecies);
  Vector<List<Photon>*>        sourcePhotons(numPhotonSpecies);

  // Regular cells
  auto regularKernel = [&](const IntVect& iv) -> void {
    if (ebisbox.isRegular(iv) && validCells(iv)) {
      const RealVect cellPos       = probLo + a_dx * (RealVect(iv) + 0.5 * RealVect::Unit);
      const RealVect centroidPos   = RealVect::Zero;
      const RealVect lo            = -0.5 * RealVect::Unit;
      const RealVect hi            = 0.5 * RealVect::Unit;
      const RealVect bndryCentroid = RealVect::Zero;
      const RealVect bndryNormal   = RealVect::Zero;
      const Real     kappa         = 1.0;

      // Populate the per-cell Ito data
      for (int i = 0; i < numItoSpecies; i++) {
        itoParticles[i]    = &((*itoParticlesFAB[i])(iv, 0));
        numNewParticles[i] = (long long)(a_newParticlesPerCell.getSingleValuedFAB()(iv, i));
        numOldParticles[i] = (long long)(a_oldParticlesPerCell.getSingleValuedFAB()(iv, i));
      }

      // Populate the per-cell CDR data
      for (int i = 0; i < numCdrSpecies; i++) {
        cdrParticles[i] = &((*cdrParticlesFAB[i])(iv, 0));
      }

      // Populate the per-cell photon data.
      for (int i = 0; i < numPhotonSpecies; i++) {
        bulkPhotons[i]   = &((*bulkPhotonsFAB[i])(iv, 0));
        sourcePhotons[i] = &((*sourcePhotonsFAB[i])(iv, 0));

        numNewPhotons[i] = (long long)(a_newPhotonsPerCell.getSingleValuedFAB()(iv, i));

        // sourcePhotons will hold the NEW number of photons to be generated -- it should already
        // have been cleared in upstream code but I'm leaving this in for safety.
        sourcePhotons[i]->clear();
      }

      // Reconcile the ItoSolver particles -- this either removes weight from the original particles (if we lost physical particles)
      // or adds new particles (if we gained physical particles)
      m_physics->reconcileParticles(itoParticles,
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

      // Reconcile the photon solver. This will generate new computational photons that are later added to the Monte Carlo photon
      // solvers.
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

      // Add the photoionization term. This will adds new particles from the photoionization reactions.
      m_physics->reconcilePhotoionization(itoParticles, cdrParticles, bulkPhotons);
    }
  };

  // Irregular cells
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    const IntVect iv = vof.gridIndex();
    if (ebisbox.isIrregular(iv) && validCells(iv, 0)) {
      const Real     kappa         = ebisbox.volFrac(vof);
      const RealVect cellPos       = probLo + Location::position(Location::Cell::Center, vof, ebisbox, a_dx);
      const RealVect centroidPos   = ebisbox.centroid(vof);
      const RealVect bndryCentroid = ebisbox.bndryCentroid(vof);
      const RealVect bndryNormal   = ebisbox.normal(vof);

      // Compute the minimum bounding box that encloses this cut-cell.
      RealVect lo = -0.5 * RealVect::Unit;
      RealVect hi = 0.5 * RealVect::Unit;
      if (kappa < 1.0) {
        DataOps::computeMinValidBox(lo, hi, bndryNormal, bndryCentroid);
      }

      // Populate the per-cell particle data.
      for (int i = 0; i < numItoSpecies; i++) {
        itoParticles[i]    = &((*itoParticlesFAB[i])(iv, 0));
        numNewParticles[i] = (long long)(a_newParticlesPerCell(vof, i));
        numOldParticles[i] = (long long)(a_oldParticlesPerCell(vof, i));
      }

      // Populate the per-cell CDR data
      for (int i = 0; i < numCdrSpecies; i++) {
        cdrParticles[i] = &((*cdrParticlesFAB[i])(iv, 0));
      }

      // Populate the per-cell photon data.
      for (int i = 0; i < numPhotonSpecies; i++) {
        bulkPhotons[i]   = &((*bulkPhotonsFAB[i])(iv, 0));
        sourcePhotons[i] = &((*sourcePhotonsFAB[i])(iv, 0));

        numNewPhotons[i] = (long long)(a_newPhotonsPerCell(vof, i));

        // sourcePhotons will hold the NEW number of photons to be generated -- it should already
        // have been cleared in upstream code but I'm leaving this in for safety.
        sourcePhotons[i]->clear();
      }

      // Reconcile the ItoSolver particles -- this either removes weight from the original particles (if we lost physical particles)
      // or adds new particles (if we gained physical particles)
      m_physics->reconcileParticles(itoParticles,
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

      // Reconcile the photon solver. This will generate new computational photons that are later added to the Monte Carlo photon
      // solvers.
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

      // Add the photoionization term. This will adds new particles from the photoionization reactions.
      m_physics->reconcilePhotoionization(itoParticles, cdrParticles, bulkPhotons);
    }
  };

  // Run the kernels.
  VoFIterator& vofit = (*m_amr->getVofIterator(m_particleRealm, m_plasmaPhase)[a_level])[a_dit];

  BoxLoops::loop(a_box, regularKernel);
  BoxLoops::loop(vofit, irregularKernel);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::reconcilePhotoionization() noexcept
{
  CH_TIME("ItoKMCStepper::reconcilePhotoionization()");
  if (m_verbosity > 5) {
    pout() << m_name + "::reconcilePhotoionization()" << endl;
  }

  const int numItoSpecies    = m_physics->getNumItoSpecies();
  const int numCdrSpecies    = m_physics->getNumCdrSpecies();
  const int numPhotonSpecies = m_physics->getNumPhotonSpecies();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {

      Vector<List<ItoParticle>*>   itoParticles(numItoSpecies);
      Vector<List<PointParticle>*> cdrParticles(numCdrSpecies);
      Vector<List<Photon>*>        photonParticles(numPhotonSpecies);

      for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
        ParticleContainer<ItoParticle>& particles = solverIt()->getParticles(ItoSolver::WhichContainer::Bulk);

        itoParticles[solverIt.index()] = &(particles[lvl][dit()].listItems());
      }

      for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
        ParticleContainer<PointParticle>& particles = m_cdrPhotoiProducts[solverIt.index()];

        cdrParticles[solverIt.index()] = &(particles[lvl][dit()].listItems());
      }

      for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
        ParticleContainer<Photon>& particles = solverIt()->getBulkPhotons();

        photonParticles[solverIt.index()] = &(particles[lvl][dit()].listItems());
      }

      m_physics->reconcilePhotoionization(itoParticles, cdrParticles, photonParticles);
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::reconcileCdrDensities(const EBAMRCellData& a_newParticlesPerCell,
                                                 const EBAMRCellData& a_oldParticlesPerCell,
                                                 const Real           a_dt) noexcept
{
  CH_TIME("ItoKMCStepper::reconcileCdrDensities(EBAMRCellDatax2, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::reconcileCdrDensities(EBAMRCellDatax2, Real)" << endl;
  }

  const int numCdrSpecies = m_physics->getNumCdrSpecies();

  CH_assert(a_newParticlesPerCell.getRealm() == m_fluidRealm);
  CH_assert(a_oldParticlesPerCell.getRealm() == m_fluidRealm);
  CH_assert(a_dt > 0.0);

  if (numCdrSpecies > 0) {

    // Increment, but don't divide by kappa.
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      this->reconcileCdrDensities(*a_newParticlesPerCell[lvl], *a_oldParticlesPerCell[lvl], lvl, a_dt);
    }

    // Redistribute if user calls for it.
    if (m_redistributeCDR) {
      for (auto it = m_cdr->iterator(); it.ok(); ++it) {
        const int idx = it.index();

        const EBAMRCellData newPPC = m_amr->slice(a_newParticlesPerCell, Interval(idx, idx));
        const EBAMRCellData oldPPC = m_amr->slice(a_oldParticlesPerCell, Interval(idx, idx));

        for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
          const DisjointBoxLayout& dbl   = m_amr->getGrids(m_fluidRealm)[lvl];
          const EBISLayout&        ebisl = m_amr->getEBISLayout(m_fluidRealm, m_plasmaPhase)[lvl];
          const Real               dx    = m_amr->getDx()[lvl];

          for (DataIterator dit(dbl); dit.ok(); ++dit) {
            const EBISBox& ebisbox = ebisl[dit()];

            BaseIVFAB<Real>& deltaMass = (*m_fluidScratchEB[lvl])[dit()];

            deltaMass.setVal(0.0);

            const EBCellFAB& newPPC = (*a_newParticlesPerCell[lvl])[dit()];
            const EBCellFAB& oldPPC = (*a_oldParticlesPerCell[lvl])[dit()];

            auto kernel = [&](const VolIndex& vof) -> void {
              const Real kappa = ebisbox.volFrac(vof);

              deltaMass(vof, 0) = (newPPC(vof, idx) - oldPPC(vof, idx)) * (1.0 - kappa) / std::pow(dx, SpaceDim);
            };

            VoFIterator& vofit = (*m_amr->getVofIterator(m_fluidRealm, m_plasmaPhase)[lvl])[dit()];

            BoxLoops::loop(vofit, kernel);
          }
        }

        const RefCountedPtr<CdrSolver>& solver = it();

        solver->redistribute(solver->getPhi(), m_fluidScratchEB);
      }
    }

    this->coarsenCDRSolvers();
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::reconcileCdrDensities(const LevelData<EBCellFAB>& a_newParticlesPerCell,
                                                 const LevelData<EBCellFAB>& a_oldParticlesPerCell,
                                                 const int                   a_level,
                                                 const Real                  a_dt) noexcept
{
  CH_TIME("ItoKMCStepper::reconcileCdrDensities(LD<EBCellFAB>x2, int, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::reconcileCdrDensities(LD<EBCellFAB>x2, int, Real)" << endl;
  }

  const int numCdrSpecies = m_physics->getNumCdrSpecies();

  CH_assert(a_newParticlesPerCell.nComp() == numCdrSpecies);
  CH_assert(a_oldParticlesPerCell.nComp() == numCdrSpecies);

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_fluidRealm)[a_level];
  const Real               dx  = m_amr->getDx()[a_level];

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    this->reconcileCdrDensities(a_newParticlesPerCell[dit()],
                                a_oldParticlesPerCell[dit()],
                                a_level,
                                dit(),
                                dbl[dit()],
                                dx,
                                a_dt);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::reconcileCdrDensities(const EBCellFAB& a_newParticlesPerCell,
                                                 const EBCellFAB& a_oldParticlesPerCell,
                                                 const int        a_level,
                                                 const DataIndex  a_dit,
                                                 const Box        a_box,
                                                 const Real       a_dx,
                                                 const Real       a_dt) noexcept
{
  CH_TIME("ItoKMCStepper::reconcileCdrDensities(EBCellFABx2, int, DataIndex, Box, Realx2)");
  if (m_verbosity > 5) {
    pout() << m_name + "::reconcileCdrDensities(EBCellFABx2, int, DataIndex, Box, Realx2)" << endl;
  }

  const int numCdrSpecies = m_physics->getNumCdrSpecies();

  CH_assert(a_newParticlesPerCell.nComp() == numCdrSpecies);
  CH_assert(a_oldParticlesPerCell.nComp() == numCdrSpecies);

  const Real volume = std::pow(a_dx, SpaceDim);

  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<CdrSolver>& solver = solverIt();
    const int                 index  = solverIt.index();

    EBCellFAB& phi = (*(solver->getPhi()[a_level]))[a_dit];
    EBCellFAB& src = (*(solver->getSource()[a_level]))[a_dit];

    // Source = (newParticles - oldParticles)/volume
    src.setVal(0.0);
    src.plus(a_newParticlesPerCell, index, 0, 1);
    src.minus(a_oldParticlesPerCell, index, 0, 1);
    src /= volume;

    // Phi += (newParticles - oldParticles)/volume
    phi += src;

    // Source = (newParticles - oldParticles)/(volume*dt)
    src /= a_dt;
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::coarsenCDRSolvers() noexcept
{
  CH_TIME("ItoKMCStepper::coarsenCDRSolvers");
  if (m_verbosity > 5) {
    pout() << m_name + "::coarsenCDRSolvers" << endl;
  }

  for (auto solverIt = this->m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    auto& solver = solverIt();

    EBAMRCellData& phi = solver->getPhi();
    EBAMRCellData& src = solver->getSource();

    this->m_amr->conservativeAverage(phi, phi.getRealm(), this->m_plasmaPhase);
    this->m_amr->conservativeAverage(src, src.getRealm(), this->m_plasmaPhase);

    this->m_amr->interpGhostPwl(phi, phi.getRealm(), this->m_plasmaPhase);
    this->m_amr->interpGhostPwl(src, src.getRealm(), this->m_plasmaPhase);

    DataOps::setCoveredValue(phi, 0.0);
    DataOps::setCoveredValue(src, 0.0);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::fillSecondaryEmissionEB(const Real a_dt) noexcept
{
  CH_TIME("ItoKMCStepper::fillSecondaryEmissionEB(Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::fillSecondaryEmissionEB(Real)" << endl;
  }

  // Particles that left the domain
  Vector<ParticleContainer<ItoParticle>*> primaryParticles;
  for (auto it = m_ito->iterator(); it.ok(); ++it) {
    ParticleContainer<ItoParticle>& intersectedParticles = it()->getParticles(ItoSolver::WhichContainer::EB);

    primaryParticles.push_back(&intersectedParticles);
  }

  // CDR solvers extrapolate their fluxes. We then copy the extrapolated fluxes to transient data holders (which are defined over the particle realm).
  EBAMRIVData tmp;
  m_amr->allocate(tmp, m_fluidRealm, m_plasmaPhase, 1);

  for (auto solverIt = m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    const int                       idx    = solverIt.index();
    const RefCountedPtr<CdrSolver>& solver = solverIt();

    EBAMRIVData& extrapFlux = m_cdrFluxesExtrap[idx];

    if (solver->isMobile()) {
      solver->extrapolateAdvectiveFluxToEB(tmp);

      extrapFlux.copy(tmp);
    }
    else {
      DataOps::setValue(extrapFlux, 0.0);
    }
  }

  // Photons that left the domain
  Vector<ParticleContainer<Photon>*> primaryPhotons;
  for (auto it = m_rte->iterator(); it.ok(); ++it) {
    ParticleContainer<Photon>& intersectedPhotons = it()->getEbPhotons();

    primaryPhotons.push_back(&intersectedPhotons);
  }

  // Call the other version.
  this->fillSecondaryEmissionEB(m_secondaryParticles,
                                m_cdrFluxes,
                                m_secondaryPhotons,
                                primaryParticles,
                                m_cdrFluxesExtrap,
                                primaryPhotons,
                                m_electricFieldParticle,
                                a_dt);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::fillSecondaryEmissionEB(Vector<ParticleContainer<ItoParticle>>&  a_secondaryParticles,
                                                   Vector<EBAMRIVData>&                     a_cdrFluxes,
                                                   Vector<ParticleContainer<Photon>>&       a_secondaryPhotons,
                                                   Vector<ParticleContainer<ItoParticle>*>& a_primaryParticles,
                                                   Vector<EBAMRIVData>&                     a_cdrFluxesExtrap,
                                                   Vector<ParticleContainer<Photon>*>&      a_primaryPhotons,
                                                   const EBAMRCellData&                     a_electricField,
                                                   const Real                               a_dt) noexcept
{
  CH_TIME("ItoKMCStepper::fillSecondaryEmissionEB(full)");
  if (m_verbosity > 5) {
    pout() << m_name + "::fillSecondaryEmissionEB(full)" << endl;
  }

  const int numItoSpecies    = m_physics->getNumItoSpecies();
  const int numCdrSpecies    = m_physics->getNumCdrSpecies();
  const int numPhotonSpecies = m_physics->getNumPhotonSpecies();

  CH_assert(a_secondaryParticles.size() == numItoSpecies);
  CH_assert(a_cdrFluxes.size() == numCdrSpecies);
  CH_assert(a_secondaryPhotons.size() == numPhotonSpecies);
  CH_assert(a_primaryParticles.size() == numItoSpecies);
  CH_assert(a_cdrFluxesExtrap.size() == numCdrSpecies);
  CH_assert(a_primaryPhotons.size() == numPhotonSpecies);
  CH_assert(a_electricField.getRealm() == m_particleRealm);
  CH_assert(a_dt >= 0.0);

  // Incoming/outgoing particle containers must be sorted by cell.
  for (int i = 0; i < numItoSpecies; i++) {
    CH_assert(a_secondaryParticles[i].getRealm() == m_particleRealm);
    CH_assert(a_primaryParticles[i]->getRealm() == m_particleRealm);

    a_secondaryParticles[i].clearParticles();
    a_secondaryParticles[i].organizeParticlesByCell();

    a_primaryParticles[i]->organizeParticlesByCell();
  }

  for (int i = 0; i < numCdrSpecies; i++) {
    CH_assert(a_cdrFluxes[i].getRealm() == m_particleRealm);
    CH_assert(a_cdrFluxesExtrap[i].getRealm() == m_particleRealm);

    DataOps::setValue(a_cdrFluxes[i], 0.0);
  }

  // Incoming/outgoing photon containers by be sorted by cell
  for (int i = 0; i < numPhotonSpecies; i++) {
    CH_assert(a_secondaryPhotons[i].getRealm() == m_particleRealm);
    CH_assert(a_primaryPhotons[i]->getRealm() == m_particleRealm);

    a_secondaryPhotons[i].clearParticles();

    a_secondaryPhotons[i].organizeParticlesByCell();
    a_primaryPhotons[i]->organizeParticlesByCell();
  }

  const RealVect probLo = m_amr->getProbLo();

  const Vector<Electrode>&  electrodes  = m_computationalGeometry->getElectrodes();
  const Vector<Dielectric>& dielectrics = m_computationalGeometry->getDielectrics();

  bool isDielectric = false;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_particleRealm)[lvl];
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_particleRealm, m_plasmaPhase)[lvl];
    const Real               dx    = m_amr->getDx()[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const EBISBox&       ebisbox       = ebisl[dit()];
      const EBCellFAB&     electricField = (*a_electricField[lvl])[dit()];
      const EBCellFAB&     oldPPC        = (*m_oldPPC[lvl])[dit()];
      const EBCellFAB&     newPPC        = (*m_newPPC[lvl])[dit()];
      const BaseFab<bool>& validCells    = (*m_amr->getValidCells(m_particleRealm)[lvl])[dit()];

      Vector<BinFab<ItoParticle>*> secondaryParticlesFAB;
      Vector<BinFab<ItoParticle>*> primaryParticlesFAB;

      Vector<BinFab<Photon>*> secondaryPhotonsFAB;
      Vector<BinFab<Photon>*> primaryPhotonsFAB;

      Vector<BaseIVFAB<Real>*> cdrFluxesFAB;
      Vector<BaseIVFAB<Real>*> cdrFluxesExtrapFAB;

      for (auto it = m_ito->iterator(); it.ok(); ++it) {
        secondaryParticlesFAB.push_back(&(a_secondaryParticles[it.index()].getCellParticles(lvl)[dit()]));
        primaryParticlesFAB.push_back(&(a_primaryParticles[it.index()]->getCellParticles(lvl)[dit()]));
      }

      for (auto it = m_cdr->iterator(); it.ok(); ++it) {
        cdrFluxesFAB.push_back(&((*(a_cdrFluxes[it.index()])[lvl])[dit()]));
        cdrFluxesExtrapFAB.push_back(&((*(a_cdrFluxes[it.index()])[lvl])[dit()]));
      }

      for (auto it = m_rte->iterator(); it.ok(); ++it) {
        secondaryPhotonsFAB.push_back(&(a_secondaryPhotons[it.index()].getCellParticles(lvl)[dit()]));
        primaryPhotonsFAB.push_back(&(a_primaryPhotons[it.index()]->getCellParticles(lvl)[dit()]));
      }

      // Kernel definition.
      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const IntVect iv = vof.gridIndex();

        if (validCells(iv)) {
          const RealVect E = RealVect(D_DECL(electricField(vof, 0), electricField(vof, 1), electricField(vof, 2)));
          const RealVect bndryNormal   = ebisbox.normal(vof);
          const RealVect bndryCentroid = ebisbox.bndryCentroid(vof);
          const RealVect cellCentroid  = ebisbox.centroid(vof);
          const RealVect cellCenter    = probLo + Location::position(Location::Cell::Center, vof, ebisbox, dx);
          const RealVect physPos       = cellCenter + bndryCentroid * dx;
          const Real     bndryArea     = ebisbox.bndryArea(vof);

          Vector<List<ItoParticle>> secondaryParticles(numItoSpecies);
          Vector<List<ItoParticle>> primaryParticles(numItoSpecies);

          Vector<Real> cdrFluxes(numCdrSpecies, 0.0);
          Vector<Real> cdrFluxesExtrap(numCdrSpecies, 0.0);

          Vector<List<Photon>> secondaryPhotons(numPhotonSpecies);
          Vector<List<Photon>> primaryPhotons(numPhotonSpecies);

          for (int i = 0; i < numItoSpecies; i++) {
            secondaryParticles[i] = (*secondaryParticlesFAB[i])(iv, 0);
            primaryParticles[i]   = (*primaryParticlesFAB[i])(iv, 0);
          }

          // Populate CDR fluxes
          for (int i = 0; i < numCdrSpecies; i++) {
            cdrFluxes[i]       = 0.0;
            cdrFluxesExtrap[i] = (*cdrFluxesExtrapFAB[i])(vof, 0);
          }

          for (int i = 0; i < numPhotonSpecies; i++) {
            secondaryPhotons[i] = (*secondaryPhotonsFAB[i])(iv, 0);
            primaryPhotons[i]   = (*primaryPhotonsFAB[i])(iv, 0);
          }

          // Figure out which material we are dealing with.
          int  matIndex = -1;
          Real minDist  = std::numeric_limits<Real>::max();

          for (int i = 0; i < electrodes.size(); i++) {
            const Real curDist = electrodes[i].getImplicitFunction()->value(physPos);

            if (std::abs(curDist) < std::abs(minDist)) {
              minDist  = curDist;
              matIndex = i;
            }
          }

          for (int i = 0; i < dielectrics.size(); i++) {
            const Real curDist = dielectrics[i].getImplicitFunction()->value(physPos);

            if (std::abs(curDist) < std::abs(minDist)) {
              minDist      = curDist;
              matIndex     = i;
              isDielectric = true;
            }
          }

          // Call the physics framework.
          m_physics->secondaryEmissionEB(secondaryParticles,
                                         cdrFluxes,
                                         secondaryPhotons,
                                         primaryParticles,
                                         cdrFluxesExtrap,
                                         primaryPhotons,
                                         E,
                                         physPos,
                                         cellCentroid,
                                         bndryCentroid,
                                         bndryNormal,
                                         bndryArea,
                                         dx,
                                         a_dt,
                                         isDielectric,
                                         matIndex);
        }
      };

      // Run the kernel.
      VoFIterator& vofit = (*m_amr->getVofIterator(m_particleRealm, m_plasmaPhase)[lvl])[dit()];
      BoxLoops::loop(vofit, irregularKernel);
    }
  }

  // Sort by patch
  for (int i = 0; i < numItoSpecies; i++) {
    a_secondaryParticles[i].organizeParticlesByPatch();
    a_primaryParticles[i]->organizeParticlesByPatch();
  }

  for (int i = 0; i < numPhotonSpecies; i++) {
    a_secondaryPhotons[i].organizeParticlesByPatch();
    a_primaryPhotons[i]->organizeParticlesByPatch();
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::resolveSecondaryEmissionEB(const Real a_dt) noexcept
{
  CH_TIME("ItoKMCStepper::resolveSecondaryEmissionEB(short)");
  if (m_verbosity > 5) {
    pout() << m_name + "::resolveSecondaryEmissionEB(short)" << endl;
  }

  Vector<ParticleContainer<ItoParticle>*> secondaryParticles;
  Vector<ParticleContainer<ItoParticle>*> primaryParticles;

  for (auto it = m_ito->iterator(); it.ok(); ++it) {
    const int idx = it.index();

    primaryParticles.push_back(&(it()->getParticles(ItoSolver::WhichContainer::EB)));
    secondaryParticles.push_back(&(m_secondaryParticles[idx]));
  }

  // Copy the CDR fluxes on the particle realm over to the fluid realm.
  Vector<EBAMRIVData*> cdrFluxes;
  for (auto it = m_cdr->iterator(); it.ok(); ++it) {
    const RefCountedPtr<CdrSolver>& solver = it();
    const int                       idx    = it.index();

    EBAMRIVData& ebFlux = solver->getEbFlux();

    ebFlux.copy(m_cdrFluxes[idx]);

    m_amr->arithmeticAverage(ebFlux, m_fluidRealm, m_plasmaPhase);

    cdrFluxes.push_back(&ebFlux);
  }

  // Handle to surface charge density.
  EBAMRIVData& surfaceChargeDensity = m_sigmaSolver->getPhi();

  this->resolveSecondaryEmissionEB(secondaryParticles, primaryParticles, cdrFluxes, surfaceChargeDensity, a_dt);

  m_sigmaSolver->resetElectrodes(0.0);
  m_amr->arithmeticAverage(surfaceChargeDensity, m_fluidRealm, m_plasmaPhase);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::resolveSecondaryEmissionEB(Vector<ParticleContainer<ItoParticle>*>& a_secondaryParticles,
                                                      Vector<ParticleContainer<ItoParticle>*>& a_primaryParticles,
                                                      Vector<EBAMRIVData*>&                    a_cdrFluxes,
                                                      EBAMRIVData&                             a_surfaceChargeDensity,
                                                      const Real                               a_dt) noexcept
{
  CH_TIME("ItoKMCStepper::resolveSecondaryEmissionEB(full)");
  if (m_verbosity > 5) {
    pout() << m_name + "::resolveSecondaryEmissionEB(full)" << endl;
  }

  const int numItoSpecies = m_physics->getNumItoSpecies();
  const int numCdrSpecies = m_physics->getNumCdrSpecies();

  CH_assert(a_secondaryParticles.size() == numItoSpecies);
  CH_assert(a_primaryParticles.size() == numItoSpecies);
  CH_assert(a_cdrFluxes.size() == numCdrSpecies);
  CH_assert(a_surfaceChargeDensity.getRealm() == m_fluidRealm);

  for (int i = 0; i < numItoSpecies; i++) {
    CH_assert(a_secondaryParticles[i]->getRealm() == m_particleRealm);
    CH_assert(a_primaryParticles[i]->getRealm() == m_particleRealm);
  }

  for (int i = 0; i < numCdrSpecies; i++) {
    CH_assert(a_secondaryParticles[i]->getRealm() == m_particleRealm);
    CH_assert(a_primaryParticles[i]->getRealm() == m_particleRealm);
  }

  // Deposit the incoming/outgoing particles on the surface and update the surface charge density.
  for (auto it = m_ito->iterator(); it.ok(); ++it) {
    const RefCountedPtr<ItoSolver>&  solver  = it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = it.index();
    const int Z   = species->getChargeNumber();

    if (Z != 0) {

      // Add charge from primary particles
      m_amr->depositParticles<ItoParticle, &ItoParticle::weight>(m_particleScratchEB,
                                                                 m_particleRealm,
                                                                 m_plasmaPhase,
                                                                 *a_primaryParticles[idx]);

      m_fluidScratchEB.copy(m_particleScratchEB);
      DataOps::incr(a_surfaceChargeDensity, m_fluidScratchEB, 1.0 * Z * Units::Qe);

      // Subtract charge from secondary particles
      m_amr->depositParticles<ItoParticle, &ItoParticle::weight>(m_particleScratchEB,
                                                                 m_particleRealm,
                                                                 m_plasmaPhase,
                                                                 *a_secondaryParticles[idx]);

      m_fluidScratchEB.copy(m_particleScratchEB);
      DataOps::incr(a_surfaceChargeDensity, m_fluidScratchEB, -1.0 * Z * Units::Qe);
    }

    // Add the secondary particles into the solvers and remove the primary particles.
    ParticleContainer<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk);

    particles.transferParticles(*a_secondaryParticles[idx]);
    a_primaryParticles[idx]->clearParticles();
  }

  // Add CDR fluxes to the CDR solvers
  for (auto it = m_cdr->iterator(); it.ok(); ++it) {
    const RefCountedPtr<CdrSolver>&  solver  = it();
    const RefCountedPtr<CdrSpecies>& species = solver->getSpecies();

    const int idx = it.index();
    const int Z   = species->getChargeNumber();

    if (Z != 0) {
      DataOps::incr(a_surfaceChargeDensity, *a_cdrFluxes[idx], Z * a_dt * Units::Qe);
    }

    // Add mass to CDR solvers -- this is an inefficient way of doing it but I don't know if it'll be a performance bottleneck as well.
    EBAMRCellData divG;
    EBAMRFluxData G;

    m_amr->allocate(divG, m_fluidRealm, m_plasmaPhase, 1);
    m_amr->allocate(G, m_fluidRealm, m_plasmaPhase, 1);

    DataOps::setValue(G, 0.0);

    solver->computeDivG(divG, G, *a_cdrFluxes[idx], false);

    EBAMRCellData& phi = solver->getPhi();
    DataOps::incr(phi, divG, -a_dt);

    m_amr->conservativeAverage(phi, m_fluidRealm, m_plasmaPhase);
    m_amr->interpGhostPwl(phi, m_fluidRealm, m_plasmaPhase);
  }

  // Conservatively coarsen the surface charge density.
  m_amr->conservativeAverage(a_surfaceChargeDensity, m_fluidRealm, m_plasmaPhase);
}

template <typename I, typename C, typename R, typename F>
Real
ItoKMCStepper<I, C, R, F>::computePhysicsDt() noexcept
{
  CH_TIME("ItoKMCStepper::computePhysicsDt()");
  if (m_verbosity > 5) {
    pout() << m_name + "::computePhysicsDt()" << endl;
  }

  const Real dt = this->computePhysicsDt(m_electricFieldFluid);

  return dt;
}

template <typename I, typename C, typename R, typename F>
Real
ItoKMCStepper<I, C, R, F>::computePhysicsDt(const EBAMRCellData& a_electricField) noexcept
{
  CH_TIME("ItoKMCStepper::computePhysicsDt(EBAMRCellFAB, Vector<EBAMRCellFAB*>)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computePhysicsDt(EBAMRCellFAB, Vector<EBAMRCellFAB*>)" << endl;
  }

  CH_assert(a_electricField.getRealm() == m_fluidRealm);

  const int numItoSpecies = m_physics->getNumItoSpecies();
  const int numCdrSpecies = m_physics->getNumCdrSpecies();

  Real minDt = std::numeric_limits<Real>::max();

  // Sort by cell so we can compute the number of particles.
  m_ito->organizeParticlesByCell(ItoSolver::WhichContainer::Bulk);

  // Compute the number of reactive particles for both Ito and CDR species.
  if (numItoSpecies > 0) {
    this->computeReactiveItoParticlesPerCell(m_particleItoPPC);

    const Interval srcInterv(0, numItoSpecies - 1);
    const Interval dstInterv(0, numItoSpecies - 1);

    m_fluidPPC.copy(srcInterv, m_particleItoPPC, dstInterv);
  }
  if (numCdrSpecies > 0) {
    this->computeReactiveCdrParticlesPerCell(m_fluidCdrPPC);

    const Interval srcInterv(0, numCdrSpecies - 1);
    const Interval dstInterv(numItoSpecies, numItoSpecies + numCdrSpecies - 1);

    m_fluidPPC.copy(srcInterv, m_fluidCdrPPC, dstInterv);
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Real levelDt = this->computePhysicsDt(*a_electricField[lvl], *m_fluidPPC[lvl], lvl);

    minDt = std::min(minDt, levelDt);
  }

  // Sort by patch.
  m_ito->organizeParticlesByPatch(ItoSolver::WhichContainer::Bulk);

  return minDt;
}

template <typename I, typename C, typename R, typename F>
Real
ItoKMCStepper<I, C, R, F>::computePhysicsDt(const LevelData<EBCellFAB>& a_electricField,
                                            const LevelData<EBCellFAB>& a_particlesPerCell,
                                            const int                   a_level) noexcept
{
  CH_TIME("ItoKMCStepper::computePhysicsDt(LD<EBCellFAB>, LD<EBCellFAB>, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computePhysicsDt(LD<EBCellFAB>, LD<EBCellFAB>, int)" << endl;
  }

  const int numPlasmaSpecies = m_physics->getNumPlasmaSpecies();

  CH_assert(a_particlesPerCell.nComp() == numPlasmaSpecies);
  CH_assert(a_electricField.nComp() == SpaceDim);

  Real minDt = std::numeric_limits<Real>::max();

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_fluidRealm)[a_level];

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Real patchDt =
      this->computePhysicsDt(a_electricField[dit()], a_particlesPerCell[dit()], a_level, dit(), dbl[dit()]);

    minDt = std::min(minDt, patchDt);
  }

  return ParallelOps::min(minDt);
}

template <typename I, typename C, typename R, typename F>
Real
ItoKMCStepper<I, C, R, F>::computePhysicsDt(const EBCellFAB& a_electricField,
                                            const EBCellFAB& a_particlesPerCell,
                                            const int        a_level,
                                            const DataIndex  a_dit,
                                            const Box        a_box) noexcept
{
  CH_TIME("ItoKMCStepper::computePhysicsDt(EBCellFAB, EBCellFAB, int, DataIndex, Box)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computePhysicsDt(EBCellFAB, EBCellFAB, int, DataIndex, Box)" << endl;
  }

  Real minDt = std::numeric_limits<Real>::max();

  const int numPlasmaSpecies = m_physics->getNumPlasmaSpecies();

  CH_assert(a_electricField.nComp() == SpaceDim);
  CH_assert(a_particlesPerCell.nComp() == numPlasmaSpecies);

  // Geometric information that we need.
  const Real     dx      = m_amr->getDx()[a_level];
  const RealVect probLo  = m_amr->getProbLo();
  const EBISBox& ebisbox = m_amr->getEBISLayout(m_fluidRealm, m_plasmaPhase)[a_level][a_dit];

  // Valid grid cells.
  const BaseFab<bool>& validCells = (*m_amr->getValidCells(m_fluidRealm)[a_level])[a_dit];

  // Handles to regular grid data.
  const FArrayBox& electricFieldReg    = a_electricField.getFArrayBox();
  const FArrayBox& particlesPerCellReg = a_particlesPerCell.getFArrayBox();

  // The kernel require storage for the per-cell number of particles. This is
  // the storage we use for that
  Vector<Physics::ItoKMC::FPR> ppc(numPlasmaSpecies);

  // Regular grid kernel.
  auto regularKernel = [&](const IntVect& iv) -> void {
    // Only regular cells not covered by a finer grid.
    if (ebisbox.isRegular(iv) && validCells(iv, 0)) {
      const RealVect E   = RealVect(D_DECL(electricFieldReg(iv, 0), electricFieldReg(iv, 1), electricFieldReg(iv, 2)));
      const RealVect pos = m_amr->getProbLo() + dx * (RealVect(iv) + 0.5 * RealVect::Unit);

      for (int i = 0; i < numPlasmaSpecies; i++) {
        ppc[i] = particlesPerCellReg(iv, i);
      }

      const Real cellDt = m_physics->computeDt(E, pos, ppc);

      minDt = std::min(minDt, cellDt);
    }
  };

  // Irregular grid kernel.
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    const IntVect& iv = vof.gridIndex();

    // Only irregular cells not covered by a finer grid.
    if (ebisbox.isIrregular(iv) && validCells(iv, 0)) {
      const RealVect E   = RealVect(D_DECL(a_electricField(vof, 0), a_electricField(vof, 1), a_electricField(vof, 2)));
      const RealVect pos = probLo + Location::position(Location::Cell::Centroid, vof, ebisbox, dx);

      for (int i = 0; i < numPlasmaSpecies; i++) {
        ppc[i] = a_particlesPerCell(vof, i);
      }

      const Real cellDt = m_physics->computeDt(E, pos, ppc);

      minDt = std::min(minDt, cellDt);
    }
  };

  // Run the kernels.
  VoFIterator& vofit = (*m_amr->getVofIterator(m_fluidRealm, m_plasmaPhase)[a_level])[a_dit];

  BoxLoops::loop(a_box, regularKernel);
  BoxLoops::loop(vofit, irregularKernel);

  return minDt;
}

template <typename I, typename C, typename R, typename F>
Real
ItoKMCStepper<I, C, R, F>::computeTotalCharge() const noexcept
{
  CH_TIME("ItoKMCStepper::computeTotalCharge()");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeTotalCharge()" << endl;
  }

  const bool kappaScale = true;

  Real totalCharge = 0.0;

  // Charge from Ito solvers.
  for (auto it = m_ito->iterator(); it.ok(); ++it) {
    const RefCountedPtr<ItoSolver>&  solver  = it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int Z = species->getChargeNumber();

    if (Z != 0) {
      const ParticleContainer<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk);

      totalCharge += Z * ParticleOps::sum<ItoParticle, &ItoParticle::weight>(particles);
    }
  }

  // Charge from CDR solvers
  for (auto it = m_cdr->iterator(); it.ok(); ++it) {
    const RefCountedPtr<CdrSolver>&  solver  = it();
    const RefCountedPtr<CdrSpecies>& species = solver->getSpecies();

    const int Z = species->getChargeNumber();

    if (Z != 0) {
      const EBAMRCellData& phi = solver->getPhi();

      totalCharge += Z * solver->computeMass(phi, kappaScale);
    }
  }

  // Surface charge
  totalCharge += m_sigmaSolver->computeMass();

  return totalCharge * Units::Qe;
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::advancePhotons(const Real a_dt) noexcept
{
  CH_TIME("ItoKMCStepper::advancePhotons(Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::advancePhotons(Real)" << endl;
  }

  // TLDR: This will add the source photons to the "bulk" photons and then advance them. If the
  //       solver is a true transient solver then the photons are moved and some of them are eventually
  //       absorbed on the mesh. If the solver is an "instanteneous" solver then all source photons
  //       are absorbed on the mesh.

  for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<McPhoto>& solver = solverIt();

    // To reiterate: photons are the photons that live in the solver and are moved around. bulkPhotons
    // are the solvers that were absorbed on the mesh, bbPhotons are the photons that collided with the EB
    // and domainPhotons are photons that moved out of the domain.
    ParticleContainer<Photon>& photons       = solver->getPhotons();
    ParticleContainer<Photon>& bulkPhotons   = solver->getBulkPhotons();
    ParticleContainer<Photon>& ebPhotons     = solver->getEbPhotons();
    ParticleContainer<Photon>& domainPhotons = solver->getDomainPhotons();
    ParticleContainer<Photon>& sourcePhotons = solver->getSourcePhotons();

    if (solver->isInstantaneous()) {
      solver->clear(photons);

      // Add source Photons
      photons.addParticles(sourcePhotons);
      solver->clear(sourcePhotons);

      // Instantaneous advance.
      solver->advancePhotonsInstantaneous(bulkPhotons, ebPhotons, domainPhotons, photons);
    }
    else {
      // Add source Photons
      photons.addParticles(sourcePhotons);
      solver->clear(sourcePhotons);

      // Stationary advance
      solver->advancePhotonsTransient(bulkPhotons, ebPhotons, domainPhotons, photons, a_dt);
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::sortPhotonsByCell(const McPhoto::WhichContainer a_which) noexcept
{
  CH_TIME("ItoKMCStepper::sortPhotonsByCell(McPhoto::WhichContainer)");
  if (m_verbosity > 5) {
    pout() << m_name + "::sortPhotonsByCell(McPhoto::WhichContainer)" << endl;
  }

  for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->sortPhotonsByCell(a_which);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::sortPhotonsByPatch(const McPhoto::WhichContainer a_which) noexcept
{
  CH_TIME("ItoKMCStepper::sortPhotonsByPatch(McPhoto::WhichContainer)");
  if (m_verbosity > 5) {
    pout() << m_name + "::sortPhotonsByPatch(McPhoto::WhichContainer)" << endl;
  }

  for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->sortPhotonsByPatch(a_which);
  }
}

template <typename I, typename C, typename R, typename F>
Vector<RefCountedPtr<ItoSolver>>
ItoKMCStepper<I, C, R, F>::getLoadBalanceSolvers() const noexcept
{
  CH_TIME("ItoKMCStepper::getLoadBalanceSolvers()");
  if (m_verbosity > 5) {
    pout() << m_name + "::getLoadBalanceSolvers()" << endl;
  }

  Vector<RefCountedPtr<ItoSolver>> lbSolvers;

  // If there's an index < 0 we load balance everything.
  bool loadBalanceAll = false;
  for (int i = 0; i < m_loadBalanceIndices.size(); i++) {
    if (m_loadBalanceIndices[i] < 0) {
      loadBalanceAll = true;
    }
  }

  if (loadBalanceAll) {
    for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
      lbSolvers.push_back(solverIt());
    }
  }
  else {
    for (int i = 0; i < m_loadBalanceIndices.size(); i++) {
      RefCountedPtr<ItoSolver>& solver = m_ito->getSolvers()[i];

      lbSolvers.push_back(solver);
    }
  }

  return lbSolvers;
}

template <typename I, typename C, typename R, typename F>
bool
ItoKMCStepper<I, C, R, F>::loadBalanceThisRealm(const std::string a_realm) const
{
  CH_TIME("TimeStepper::loadBalanceThisRealm");
  if (m_verbosity > 5) {
    pout() << "TimeStepper::loadBalanceThisRealm" << endl;
  }

  bool ret = false;

  if (a_realm == m_particleRealm && m_loadBalanceParticles) {
    ret = true;
  }
  else if (a_realm == m_fluidRealm && m_loadBalanceFluid) {
    ret = true;
  }

  return ret;
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::loadBalanceBoxes(Vector<Vector<int>>&             a_procs,
                                            Vector<Vector<Box>>&             a_boxes,
                                            const std::string                a_realm,
                                            const Vector<DisjointBoxLayout>& a_grids,
                                            const int                        a_lmin,
                                            const int                        a_finestLevel)
{
  CH_TIME("ItoKMCStepper::loadBalanceBoxes");
  if (m_verbosity > 5) {
    pout() << m_name + "::loadBalanceBoxes" << endl;
  }

  if (m_loadBalanceParticles && a_realm == m_particleRealm) {
    this->loadBalanceParticleRealm(a_procs, a_boxes, a_realm, a_grids, a_lmin, a_finestLevel);
  }
  else if (m_loadBalanceFluid && a_realm == m_fluidRealm) {
    this->loadBalanceFluidRealm(a_procs, a_boxes, a_realm, a_grids, a_lmin, a_finestLevel);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::loadBalanceParticleRealm(Vector<Vector<int>>&             a_procs,
                                                    Vector<Vector<Box>>&             a_boxes,
                                                    const std::string                a_realm,
                                                    const Vector<DisjointBoxLayout>& a_grids,
                                                    const int                        a_lmin,
                                                    const int                        a_finestLevel) noexcept
{
  CH_TIME("ItoKMCStepper::loadBalanceParticleRealm(...)");
  if (m_verbosity > 5) {
    pout() << m_name + "::loadBalanceParticleRealm(...)" << endl;
  }

  // TLDR: This is a bit involved due to the fact that the simulation lives in a state between the old grids
  //       and the new grids. We want to compute the number of computational in patch on the new grid, and
  //       use that for load balancing. We have already computed the number of computational particles per
  //       grid cell on the old grids, but the new grids are not ready (yet). We only have the proxy-grids
  //       coming in through the argument (a_grids), and our job is to take this grid and reassign the patches
  //       so that each MPI rank gets roughly the same number of computational particles. To do this we perform
  //       the following steps:
  //
  //       1. Allocate storage on the proxy grids (a_grids) so we have something to regrid into.
  //       2. Define regrid operators for going between the proxy grids and the old grids.
  //       3. Regrid the PPC on the old grids onto the proxy grids.
  //       4. Go through the patches on the proxy grids and figure out the total number of particles
  //          in each patch (there's a weird global-to-local remapping taking place through intCode()).
  //       5. Call our nifty load-balancing routines.
  //

  if (!m_loadBalanceParticles) {
    MayDay::Error("ItoKMCStepper::loadBalanceParticleRealm -- logic bust, should not have been called!");
  }

  // Get the solvers that we will use for load balancing.
  Vector<RefCountedPtr<ItoSolver>> lbSolvers = this->getLoadBalanceSolvers();

  // Decompose the DisjointBoxLayout
  a_procs.resize(1 + a_finestLevel);
  a_boxes.resize(1 + a_finestLevel);

  for (int lvl = a_lmin; lvl <= a_finestLevel; lvl++) {
    a_procs[lvl] = a_grids[lvl].procIDs();
    a_boxes[lvl] = a_grids[lvl].boxArray();
  }

  // 1. Allocate something that we can regrid the PPC for each species into, and something that holds the total
  // PPC on the new grids.
  EBAMRCellData totalPPC;
  EBAMRCellData speciesPPC;

  m_amr->allocate(totalPPC, m_particleRealm, m_plasmaPhase, 1);
  m_amr->allocate(speciesPPC, m_particleRealm, m_plasmaPhase, 1);

  DataOps::setValue(totalPPC, 0.0);
  DataOps::setValue(speciesPPC, 0.0);

  // 2. EBFineInterp is not a part of the registry for ItoSolver so we just define it here ourselves. Note that
  // it is stored on the same level that we interpolate to.
  Vector<RefCountedPtr<EBFineInterp>> interpOp(1 + a_finestLevel);
  for (int lvl = 1; lvl <= a_finestLevel; lvl++) {
    const EBLevelGrid& eblgFine = *m_amr->getEBLevelGrid(m_particleRealm, m_plasmaPhase)[lvl];
    const EBLevelGrid& eblgCoar = *m_amr->getEBLevelGrid(m_particleRealm, m_plasmaPhase)[lvl - 1];
    const int          refRat   = m_amr->getRefinementRatios()[lvl - 1];
    const int          nComp    = 1;

    interpOp[lvl] = RefCountedPtr<EBFineInterp>(
      new EBFineInterp(eblgFine, eblgCoar, refRat, nComp, eblgFine.getEBIS()));
  }

  // 3. Go through each solver and figure out the number of particles on the new grids. Add
  // these to totalPPC.
  for (int i = 0; i < lbSolvers.size(); i++) {
    const EBAMRCellData& oldData        = m_loadBalancePPC[i];
    const int            oldFinestLevel = oldData.size() - 1;

    // These levels have not changed but ownship MIGHT have changed.
    for (int lvl = 0; lvl <= std::max(0, a_lmin - 1); lvl++) {
      oldData[lvl]->copyTo(*speciesPPC[lvl]);
    }

    // These levels have changed.
    for (int lvl = std::max(1, a_lmin); lvl <= a_finestLevel; lvl++) {
      RefCountedPtr<EBFineInterp>& interpolator = interpOp[lvl];

      interpolator->regridNoSlopes(*speciesPPC[lvl], *speciesPPC[lvl - 1], Interval(0, 0));

      // There could be parts of the new grid that overlapped with the old grid (on level lvl) -- we don't want
      // to pollute the solution with interpolation there since we already have valid data.
      if (lvl <= std::min(oldFinestLevel, a_finestLevel)) {
        oldData[lvl]->copyTo(*speciesPPC[lvl]);
      }
    }

    // Add to totalPPC.
    DataOps::incr(totalPPC, speciesPPC, 1.0);
  }

  // 4. totalPPC contains the total number of computational particles per cell on the new grids,
  // we need to map this to something we can load balance.
  Vector<Vector<long int>> loads(1 + a_finestLevel, 0L);
  for (int lvl = 0; lvl <= a_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = a_grids[lvl];

    Vector<long int>& levelLoads = loads[lvl];

    levelLoads.resize(dbl.size());

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box            cellBox    = dbl[dit()];
      const EBCellFAB&     PPC        = (*totalPPC[lvl])[dit()];
      const EBISBox&       ebisbox    = PPC.getEBISBox();
      const BaseFab<bool>& validCells = (*m_amr->getValidCells(m_particleRealm)[lvl])[dit()];
      const FArrayBox&     regPPC     = PPC.getFArrayBox();

      auto regularKernel = [&](const IntVect& iv) -> void {
        if (validCells(iv, 0) && ebisbox.isRegular(iv)) {
          levelLoads[dit().intCode()] += (long int)regPPC(iv, 0);
        }
      };

      BoxLoops::loop(cellBox, regularKernel);
    }

    ParallelOps::vectorSum(levelLoads);

    // Add the "constant" load from the other PPC stuff
    for (LayoutIterator lit = dbl.layoutIterator(); lit.ok(); ++lit) {
      const Box cellBox = dbl[lit()];

      levelLoads[lit().intCode()] += (long int)m_loadPerCell * cellBox.numPts();
    }
  }

  // 5. Finally do the actual load balancing.
  LoadBalancing::sort(a_boxes, loads, m_boxSort);
  LoadBalancing::balanceLevelByLevel(a_procs, loads, a_boxes);
  //  LoadBalancing::hierarchy(a_procs, loads, a_boxes); If you want to try something crazy...
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::loadBalanceFluidRealm(Vector<Vector<int>>&             a_procs,
                                                 Vector<Vector<Box>>&             a_boxes,
                                                 const std::string                a_realm,
                                                 const Vector<DisjointBoxLayout>& a_grids,
                                                 const int                        a_lmin,
                                                 const int                        a_finestLevel) noexcept
{
  CH_TIME("ItoKMCStepper::loadBalanceFluidRealm(...)");
  if (m_verbosity > 5) {
    pout() << m_name + "::loadBalanceFluidRealm(...)" << endl;
  }

  CH_assert(m_loadBalanceFluid);
  CH_assert(a_realm == m_fluidRealm);

  // TLDR: This code tries to compute a load for each grid patch by applying a relaxation operator to each box. This means that the load
  //       should be a decent estimate that takes into account boundary conditions, coarse-fine interface arithmetic, and enlargened stencils
  //       near the embedded boundary.

  a_procs.resize(1 + a_finestLevel);
  a_boxes.resize(1 + a_finestLevel);

  // We need to make AmrMesh restore some operators that we need in order to create a multigrid object. Fortunately, FieldSolver has routines
  // for doing that but it will not know if AmrMesh has updated it's operators or not. So, we need to regrid them.
  m_amr->regridOperators(m_fluidRealm, a_lmin);

  // Field solver needs to allocate solver and set up the multigrid solver.
  m_fieldSolver->allocate();
  m_fieldSolver->setupSolver();

  // Field solver implementation gets the responsibility of computing loads on each level.
  for (int lvl = 0; lvl <= a_finestLevel; lvl++) {
    Vector<long long> loads = m_fieldSolver->computeLoads(a_grids[lvl], lvl);

    // Do the desired sorting and load balancing
    a_boxes[lvl] = a_grids[lvl].boxArray();

    LoadBalancing::sort(a_boxes[lvl], loads, m_boxSort);
    LoadBalancing::makeBalance(a_procs[lvl], loads, a_boxes[lvl]);
  }
}

template <typename I, typename C, typename R, typename F>
Vector<long int>
ItoKMCStepper<I, C, R, F>::getCheckpointLoads(const std::string a_realm, const int a_level) const
{
  CH_TIME("ItoKMCStepper::getCheckpointLoads(...)");
  if (m_verbosity > 5) {
    pout() << m_name + "::getCheckpointLoads(...)" << endl;
  }

  const DisjointBoxLayout& dbl  = m_amr->getGrids(a_realm)[a_level];
  const int                nbox = dbl.size();

  Vector<long int> loads(nbox, 0L);

  if (m_loadBalanceParticles && a_realm == m_particleRealm) {

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

template <typename I, typename C, typename R, typename F>
void
ItoKMCStepper<I, C, R, F>::computeEdotJSource(const Real a_dt) noexcept
{
  CH_TIME("ItoKMCStepper::computeEdotJSource(a_dt)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeEdotJSource(a_dt)" << endl;
  }

  CH_assert(a_dt > 0.0);

  DataOps::setValue(m_EdotJ, 0.0);

  CH_assert(m_EdotJ.getRealm() == m_fluidRealm);

  // TLDR: EdotJ is an energy term for the various species, i.e. it is the rate of energy increase as the particle moves from
  //       position A to position B, excluding friction from collision with other molecules. We compute this energy increase as
  //
  //          q * V(B) - V(A)
  //
  //       which means that the energy rate is q*(V(B) - V(A))/a_dt.
  //
  //       We simply assign this factor to the particles and then deposit them on the mesh. However, this is more complex than
  //       it sounds because the particles m_EdotJ live on different realms. The way we do this is that we copy the potential
  //       over into the particle realm and we interpolate V(B) and V(A) onto some storage in the particle container. We then
  //       assign an effective weight w * [V(B) - V(A)] to the particles which we deposit onto the mesh using the appropriate
  //       deposition scheme that the user has assgned.
  //

  // Allocate a particle data holder with three scalar storage spaces; these are the weight, V(A), and V(B). We also
  // need the positions A and B so they're in here as well.
  using CompParticle = GenericParticle<3, 1>;

  ParticleContainer<CompParticle> computationParticles;
  m_amr->allocate(computationParticles, m_particleRealm);

  // Electrostatic potential on appropriate phase. This is defined on the fluid realm
  // but we need it on the particle realm.
  const EBAMRCellData potentialPhase = m_amr->alias(m_plasmaPhase, m_fieldSolver->getPotential());
  m_particleScratch1.copy(potentialPhase);

  m_amr->conservativeAverage(m_particleScratch1, m_particleRealm, m_plasmaPhase);
  m_amr->interpGhost(m_particleScratch1, m_particleRealm, m_plasmaPhase);

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int  idx       = solverIt.index();
    const int  Z         = species->getChargeNumber();
    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    if (Z != 0 && (mobile || diffusive)) {

      const ParticleContainer<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk);

      // Copy the ItoParticles to the transient particles we use for computing these things.
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[lvl];

        for (DataIterator dit(dbl); dit.ok(); ++dit) {

          List<CompParticle>&      patchCompParticles = computationParticles[lvl][dit()].listItems();
          const List<ItoParticle>& patchItoParticles  = particles[lvl][dit()].listItems();

          for (ListIterator<ItoParticle> lit(patchItoParticles); lit.ok(); ++lit) {
            const ItoParticle& itoParticle = lit();

            CompParticle p;

            // vect<0> holds the starting position. real<0> holds the weight and
            // real<1> is used for V(A) and real<2> is used for V(B)
            p.position()         = itoParticle.position();
            p.template real<0>() = itoParticle.weight();
            p.template real<1>() = 0.0;
            p.template real<2>() = 0.0;
            p.template vect<0>() = itoParticle.oldPosition();

            patchCompParticles.add(p);
          }
        }
      }

      // Interpolate the potential to the current particle position which gives us V(B) for the particles.
      m_amr->interpolateParticles<CompParticle, &CompParticle::template real<1>>(computationParticles,
                                                                                 m_particleRealm,
                                                                                 m_plasmaPhase,
                                                                                 m_particleScratch1,
                                                                                 solver->getDeposition(),
                                                                                 false);

      // Move the particles back to their old positions and interpolate the potential there.
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[lvl];

        for (DataIterator dit(dbl); dit.ok(); ++dit) {
          List<CompParticle>& particles = computationParticles[lvl][dit()].listItems();

          for (ListIterator<CompParticle> lit(particles); lit.ok(); ++lit) {

            // Swap positions.
            const RealVect tmp = lit().position();

            // vect<0> holds the end position.
            lit().position()         = lit().template vect<0>();
            lit().template vect<0>() = tmp;
          }
        }
      }

      computationParticles.remap();

      // Interpolate the potential to the previous particle position which gives us V(A) for the particles.
      m_amr->interpolateParticles<CompParticle, &CompParticle::template real<2>>(computationParticles,
                                                                                 m_particleRealm,
                                                                                 m_plasmaPhase,
                                                                                 m_particleScratch1,
                                                                                 solver->getDeposition());

      // Move the particles back again and multiply the weight by Z * (V(A) - V(B))/a_dt
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[lvl];

        for (DataIterator dit(dbl); dit.ok(); ++dit) {
          List<CompParticle>& particles = computationParticles[lvl][dit()].listItems();

          for (ListIterator<CompParticle> lit(particles); lit.ok(); ++lit) {
            lit().position() = lit().template vect<0>();

            // Compute weight * V(B) - V(A) where V(B) is stored in real<1> and V(A) is stored in real<2>
            lit().template real<0>() *= (lit().template real<1>() - lit().template real<2>());
          }
        }
      }

      computationParticles.remap();

      // Deposit the particles
      m_amr->depositParticles<CompParticle, &CompParticle::template real<0>>(m_particleScratch1,
                                                                             m_particleRealm,
                                                                             m_plasmaPhase,
                                                                             computationParticles,
                                                                             solver->getDeposition(),
                                                                             solver->getCoarseFineDeposition());

      // Copy data back onto the fluid realm.
      m_fluidScratch1.copy(m_particleScratch1);
      DataOps::scale(m_fluidScratch1, Z * Units::Qe / a_dt);
      DataOps::plus(m_EdotJ, m_fluidScratch1, 0, idx, 1);
    }

    computationParticles.clearParticles();
  }
}

#include <CD_NamespaceFooter.H>

#endif
