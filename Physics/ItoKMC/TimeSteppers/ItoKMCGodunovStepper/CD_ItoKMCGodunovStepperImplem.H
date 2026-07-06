/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ItoKMCGodunovStepperImplem.H
  @brief  Implementation of CD_ItoKMCGodunovStepper.H
  @author Robert Marskar
*/

#ifndef CD_ITOKMCGODUNOVSTEPPERIMPLEM_H
#define CD_ITOKMCGODUNOVSTEPPERIMPLEM_H

// Std includes
#include <limits>

// Chombo includes
#include <ParmParse.H>
#include <BaseIF.H>

// Our includes
#include <CD_ItoKMCGodunovStepper.H>
#include <CD_Timer.H>
#include <CD_ParallelOps.H>
#include <CD_DataOps.H>
#include <CD_ParticleLoops.H>
#include <CD_Units.H>
#include <CD_Photon.H>
#include <CD_DischargeIO.H>
#include <CD_ParticleManagement.H>
#include <CD_NearestNeighborParticleMerge.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoKMC;

namespace {
  // Deposit an SoA particle container's weight onto a_phi using an Ito solver's configured deposition scheme.
  inline void
  depositPointParticlesLikeSolver(const RefCountedPtr<AmrMesh>&   a_amr,
                                  const std::string&              a_realm,
                                  const phase::which_phase        a_phase,
                                  const RefCountedPtr<ItoSolver>& a_solver,
                                  EBAMRCellData&                  a_phi,
                                  ParticleContainer<NoPayload>&   a_particles) noexcept
  {
    a_amr->depositWeight(a_phi,
                         a_realm,
                         a_phase,
                         a_particles,
                         a_solver->getDeposition(),
                         a_solver->getCoarseFineDeposition(),
                         false);
  }

} // namespace

template <typename I, typename C, typename R, typename F>
ItoKMCGodunovStepper<I, C, R, F>::ItoKMCGodunovStepper(RefCountedPtr<ItoKMCPhysics>& a_physics, bool a_parseOptions)
  : ItoKMCStepper<I, C, R, F>(a_physics)
{
  CH_TIME("ItoKMCGodunovStepper::ItoKMCGodunovStepper");

  this->m_name                     = "ItoKMCGodunovStepper";
  this->m_prevDt                   = 0.0;
  this->m_writeCheckpointParticles = false;
  this->m_readCheckpointParticles  = false;
  this->m_extendConductivityEB     = false;
  this->m_canRegridOnRestart       = true;
  this->m_prevDt                   = 0.0;
  this->m_maxFieldAbort            = std::numeric_limits<Real>::max();

  if (a_parseOptions) {
    this->parseOptions();
  }
}

template <typename I, typename C, typename R, typename F>
ItoKMCGodunovStepper<I, C, R, F>::~ItoKMCGodunovStepper()
{
  CH_TIME("ItoKMCGodunovStepper::~ItoKMCGodunovStepper");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCGodunovStepper::~ItoKMCGodunovStepper" << endl;
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::registerOperators() noexcept
{
  CH_TIME("ItoKMCGodunovStepper::registerOperators");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCGodunovStepper::registerOperators" << endl;
  }

  ItoKMCStepper<I, C, R, F>::registerOperators();

  // Register this because we must be able to deposit particles inside dielectrics (due to particle diffusion across the EB).
  (this->m_amr)->registerOperator(s_particle_mesh, this->m_particleRealm, phase::solid);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::registerRealms() noexcept
{
  CH_TIME("ItoKMCGodunovStepper::registerRealms");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCGodunovStepper::registerRealms" << endl;
  }

  ItoKMCStepper<I, C, R, F>::registerRealms();

  // Nearest-neighbor merge debug test (see postInitialize()) -- the ghost-mask width it needs
  // must be registered here, BEFORE Driver::setupFresh() builds the initial grids
  // (Driver::setupFresh() calls registerRealms() before m_amr->regridAmr(), but
  // registerOperators() only after -- so this is the earliest hook available). A particle ghost
  // mask is only actually built during a regrid (AmrMesh::registerParticleGhostMask() itself just
  // records that a width should be tracked from here on); registering it any later would leave
  // the mask undefined for the initial grid.
  bool debugTestNNMerge = false;

  ParmParse pp(this->m_name.c_str());
  pp.query("debug_test_nn_merge", debugTestNNMerge);

  if (debugTestNNMerge) {
    (this->m_amr)->registerParticleGhostMask(this->m_particleRealm, 1);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::allocate() noexcept
{
  CH_TIME("ItoKMCGodunovStepper::allocate");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCGodunovStepper::allocate" << endl;
  }

  ItoKMCStepper<I, C, R, F>::allocate();

  // Now allocate for the conductivity particles and rho^dagger particles. This is only done in the 'allocate' routine
  // and not in 'allocateInternals' because that would discard the particles during regrids. That has definitely never
  // happen, and there's no way I've spent countless hours tracking down such a bug.
  const int numItoSpecies = this->m_physics->getNumItoSpecies();

  m_conductivityParticles.resize(numItoSpecies);
  m_irregularParticles.resize(numItoSpecies);
  m_rhoDaggerParticles.resize(numItoSpecies);

  for (auto solverIt = (this->m_ito)->iterator(); solverIt.ok(); ++solverIt) {
    const int idx = solverIt.index();

    m_conductivityParticles[idx] = RefCountedPtr<ParticleContainer<NoPayload>>(new ParticleContainer<NoPayload>());
    m_irregularParticles[idx]    = RefCountedPtr<ParticleContainer<NoPayload>>(new ParticleContainer<NoPayload>());
    m_rhoDaggerParticles[idx]    = RefCountedPtr<ParticleContainer<NoPayload>>(new ParticleContainer<NoPayload>());

    (this->m_amr)->allocate(*m_conductivityParticles[idx], this->m_particleRealm);
    (this->m_amr)->allocate(*m_irregularParticles[idx], this->m_particleRealm);
    (this->m_amr)->allocate(*m_rhoDaggerParticles[idx], this->m_particleRealm);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::allocateInternals() noexcept
{
  CH_TIME("ItoKMCGodunovStepper::allocateInternals");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::allocateInternals" << endl;
  }

  ItoKMCStepper<I, C, R, F>::allocateInternals();

  const int numCdrSpecies = this->m_physics->getNumCdrSpecies();

  m_cdrDivD.resize(numCdrSpecies);
  for (int i = 0; i < numCdrSpecies; i++) {
    this->m_amr->allocate(m_cdrDivD[i], this->m_fluidRealm, this->m_plasmaPhase, 1);
  }

  this->m_amr->allocate(m_semiImplicitRhoCDR, this->m_fluidRealm, this->m_plasmaPhase, 1);
  this->m_amr->allocate(m_semiImplicitConductivityCDR, this->m_fluidRealm, this->m_plasmaPhase, 1);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::barrier() const noexcept
{
  CH_TIME("ItoKMCGodunovStepper::barrier");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::barrier" << endl;
  }

  if ((this->m_profile)) {
    ParallelOps::barrier();
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::parseOptions() noexcept
{
  CH_TIME("ItoKMCGodunovStepper::parseOptions");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::parseOptions" << endl;
  }

  ItoKMCStepper<I, C, R, F>::parseOptions();

  this->parseAlgorithm();
  this->parseFiltering();
  this->parseCheckpointParticles();
  this->parseSecondaryEmissionSpecification();
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::parseRuntimeOptions() noexcept
{
  CH_TIME("ItoKMCGodunovStepper::parseRuntimeOptions");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::parseRuntimeOptions" << endl;
  }

  ItoKMCStepper<I, C, R, F>::parseRuntimeOptions();

  this->parseAlgorithm();
  this->parseFiltering();
  this->parseCheckpointParticles();
  this->parseSecondaryEmissionSpecification();
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::parseAlgorithm() noexcept
{
  CH_TIME("ItoKMCGodunovStepper::parseAlgorithm");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::parseAlgorithm" << endl;
  }

  ParmParse   pp(this->m_name.c_str());
  std::string str;

  pp.get("extend_conductivity", m_extendConductivityEB);
  pp.get("algorithm", str);
  pp.get("abort_max_field", m_maxFieldAbort);

  // Get algorithm
  if (str == "euler_maruyama") {
    m_algorithm = WhichAlgorithm::EulerMaruyama;
  }
  else {
    MayDay::Abort("ItoKMCGodunovStepper::parseAlgorithm - unknown algorithm requested");
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::parseFiltering() noexcept
{
  CH_TIME("ItoKMCGodunovStepper::parseFiltering");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::parseFiltering" << endl;
  }

  ParmParse   pp(this->m_name.c_str());
  std::string str;

  m_rhoFilterNum       = -1;
  m_rhoFilterMaxStride = 1;
  m_rhoFilterAlpha     = 0.5;

  m_condFilterNum       = -1;
  m_condFilterMaxStride = 1;
  m_condFilterAlpha     = 0.5;

  pp.get("rho_filter_num", m_rhoFilterNum);
  pp.get("rho_filter_max_stride", m_rhoFilterMaxStride);
  pp.get("rho_filter_alpha", m_rhoFilterAlpha);

  pp.get("cond_filter_num", m_condFilterNum);
  pp.get("cond_filter_max_stride", m_condFilterMaxStride);
  pp.get("cond_filter_alpha", m_condFilterAlpha);

  if (m_rhoFilterAlpha <= 0.0 || m_rhoFilterAlpha >= 1.0) {
    MayDay::Abort("ItoKMCGodunovStepper::parseFiltering -- cannot have alpha <= 0 or alpha >= 1 for rho_filter");
  }
  if (m_condFilterAlpha <= 0.0 || m_condFilterAlpha >= 1.0) {
    MayDay::Abort("ItoKMCGodunovStepper::parseFiltering -- cannot have alpha <= 0 or alpha >= 1 for cond_filter");
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::parseCheckpointParticles() noexcept
{
  CH_TIME("ItoKMCGodunovStepper::parseCheckpointParticles");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::parseCheckpointParticles" << endl;
  }

  ParmParse pp(this->m_name.c_str());

  pp.query("checkpoint_particles", m_writeCheckpointParticles);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::parseSecondaryEmissionSpecification() noexcept
{
  CH_TIME("ItoKMCGodunovStepper::parseSecondaryEmissionSpecifiation");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::parseSecondaryEmissionSpecification" << endl;
  }

  ParmParse pp(this->m_name.c_str());

  std::string str;

  pp.query("secondary_emission", str);

  if (str == "before_reactions") {
    m_emitSecondaryParticlesBeforeReactions = true;
  }
  else if (str == "after_reactions") {
    m_emitSecondaryParticlesBeforeReactions = false;
  }
  else {
    std::string err;

    err = "ItoKMCGodunovStepper::parseSecondaryEmissionSpecification - expected 'before_reactions' or 'after_reactions'";
    err += "but got" + str;

    MayDay::Abort(err.c_str());
  }
}

template <typename I, typename C, typename R, typename F>
Real
ItoKMCGodunovStepper<I, C, R, F>::computeDt()
{
  CH_TIME("ItoKMCGodunovStepper::computeDt");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::computeDt" << endl;
  }

  Real dt = ItoKMCStepper<I, C, R, F>::computeDt();

  if ((this->m_maxReducedField > m_maxFieldAbort) && (m_maxFieldAbort > 0.0)) {
    pout() << this->m_name + " stopping because maximum field is too high (" << this->m_maxReducedField << ")" << endl;

    this->m_keepGoing = false;
  }

  return dt;
}

template <typename I, typename C, typename R, typename F>
Real
ItoKMCGodunovStepper<I, C, R, F>::advance(const Real a_dt)
{
  CH_TIME("ItoKMCGodunovStepper::advance");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::advance" << endl;
  }

  auto debugCharge = [this](const std::string a_message) -> void {
    const Real Qtot = this->computeTotalCharge();
    if (procID() == 0) {
      std::cout << a_message << " = " << Qtot << "\n";
    }
  };

  // Special flag for telling the class that we have the necessary things in place for doing a regrid. This is
  // an if-but-maybe situation where the user chose not to checkpoint the particles we need for regrids, yet tries
  // to restart a simulation and regrid without all the prerequisites being in place. This flag is set to true
  // because these requirements are checked during the regrid routine.
  m_canRegridOnRestart = true;

  m_timer = Timer("ItoKMCGodunovStepper::advance");

  //  debugCharge("advance");

  // Previous time step is needed when regridding.
  this->m_prevDt = a_dt;

  // Done only so we can plot the absorbed photons (advanceReactionNetwork absorbs them)
  m_timer.startEvent("Deposit photons");
  for (auto solverIt = (this->m_rte)->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<McPhoto> solver = solverIt();

    EBAMRCellData&             phi     = solver->getPhi();
    ParticleContainer<Photon>& photons = solver->getBulkPhotons();

    solver->depositPhotons(phi, photons, DepositionType::NGP);
  }
  m_timer.stopEvent("Deposit photons");

  // ====== BEGIN TRANSPORT STEP ======
  // Semi-implicitly advance the particles and the field.
  switch (m_algorithm) {
  case WhichAlgorithm::EulerMaruyama: {
    this->advanceEulerMaruyama(a_dt);

    break;
  }
  default: {
    MayDay::Abort("ItoKMCGodunovStepper::advance - logic bust");

    break;
  }
  }

  // Do intersection test and remove particles that struck the EB or domain. Transfer them to appropriate containers. Then recompute the number of particles per cell.
  this->barrier();
  m_timer.startEvent("EB/Particle intersection");
  if (m_extendConductivityEB) {

    // TLDR: This hook does a special intersection routine where instead of transferring the particles that were intersected, they
    //       are put in a separate data container. We want to do this in order to reduce artificial gradients in the particle
    //       densities near the EB. The way we do this is that we flag the particles during the intersection; particles that are flagged
    //       are transferred to a container which is added to the conductivity-related particles. The particles are later removed.
    const std::function<void(ParticleSoA<ItoParticle>&, std::size_t)> setFlag = [](ParticleSoA<ItoParticle>& leaf,
                                                                                   const std::size_t i) -> void {
      leaf.template get<&ItoParticle::scratch>(i) = 1.0;
    };
    const std::function<void(ParticleSoA<ItoParticle>&, std::size_t)> nonDeletionModifier =
      [](ParticleSoA<ItoParticle>& leaf, const std::size_t i) -> void {
      leaf.template get<&ItoParticle::scratch>(i) = -1.0;
    };

    // Set flag for identifying which particles were intersected by not removed.
    for (auto it = this->m_ito->iterator(); it.ok(); ++it) {
      ParticleOps::setData(it()->getParticles(ItoSolver::WhichContainer::Bulk), setFlag);
    }

    // Intersect particles, but don't remove them.
    const bool deleteParticles = false;
    this->intersectParticles(SpeciesSubset::AllMobileOrDiffusive, deleteParticles, nonDeletionModifier);

    // Particles that were not removed are copied to a separate data container.
    for (auto it = this->m_ito->iterator(); it.ok(); ++it) {
      const RefCountedPtr<ItoSolver>&  solver  = it();
      const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

      const int idx = it.index();
      const int Z   = species->getChargeNumber();

      ParticleContainer<NoPayload>&   irregParticles = *m_irregularParticles[idx];
      ParticleContainer<ItoParticle>& ebParticles    = solver->getParticles(ItoSolver::WhichContainer::EB);
      ParticleContainer<ItoParticle>& bulkParticles  = it()->getParticles(ItoSolver::WhichContainer::Bulk);

      irregParticles.clearParticles();

      if (Z != 0 && solver->isMobile()) {
        for (int lvl = 0; lvl <= (this->m_amr)->getFinestLevel(); lvl++) {
          const DisjointBoxLayout& dbl = (this->m_amr)->getGrids((this->m_particleRealm))[lvl];
          const DataIterator&      dit = dbl.dataIterator();

          const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
          for (int mybox = 0; mybox < nbox; mybox++) {
            const DataIndex& din = dit[mybox];

            ParticleSoA<NoPayload>&         pointParticles = irregParticles[lvl][din];
            const ParticleSoA<ItoParticle>& leaf           = bulkParticles[lvl][din];

            for (std::size_t i = 0; i < leaf.size(); i++) {
              if (leaf.template get<&ItoParticle::scratch>(i) < 0.0) {
                const RealVect pos      = leaf.position(i);
                const Real     weight   = leaf.weight(i);
                const Real     mobility = leaf.template get<&ItoParticle::mobility>(i);

                pointParticles.append(pos, weight * mobility);
              }
            }
          }
        }
      }
    }

    // Finally, delete the original particles that were flagged (scratch < 0) via swap-and-pop on the SoA leaves.
    for (auto it = this->m_ito->iterator(); it.ok(); ++it) {
      ParticleContainer<ItoParticle>& bulkParticles = it()->getParticles(ItoSolver::WhichContainer::Bulk);

      for (int lvl = 0; lvl <= (this->m_amr)->getFinestLevel(); lvl++) {
        const DisjointBoxLayout& dbl = (this->m_amr)->getGrids((this->m_particleRealm))[lvl];
        const DataIterator&      dit = dbl.dataIterator();

        const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
        for (int mybox = 0; mybox < nbox; mybox++) {
          const DataIndex& din = dit[mybox];

          ParticleSoA<ItoParticle>& leaf = bulkParticles[lvl][din];

          std::size_t i = 0;
          while (i < leaf.size()) {
            if (leaf.template get<&ItoParticle::scratch>(i) < 0.0) {
              leaf.remove(i);
            }
            else {
              i++;
            }
          }
        }
      }
    }
  }
  else {
    const bool deleteParticles = true;

    this->intersectParticles(SpeciesSubset::AllMobileOrDiffusive, deleteParticles);
  }

  // The intersection tests above may not have caught all particles -- do a cleanup sweep where particles on the wrong side of
  // the EB are put on the EB.
  for (auto it = this->m_ito->iterator(); it.ok(); ++it) {
    const RefCountedPtr<ItoSolver>& solver = it();

    ParticleContainer<ItoParticle>& ebParticles   = solver->getParticles(ItoSolver::WhichContainer::EB);
    ParticleContainer<ItoParticle>& bulkParticles = solver->getParticles(ItoSolver::WhichContainer::Bulk);

    this->m_amr->transferIrregularParticles(ebParticles, bulkParticles, this->m_plasmaPhase);
  }
  m_timer.stopEvent("EB/Particle intersection");
  // ====== END TRANSPORT STEP ======

  // Photon transport
  this->barrier();
  m_timer.startEvent("Photon transport");
  this->advancePhotons(a_dt);
  m_timer.stopEvent("Photon transport");

  // Compute the gradients of the various species densities - this is used in the KMC kernels.
  if ((this->m_physics)->needGradients()) {
    m_timer.startEvent("Gradient calculation");
    (this->m_ito)->depositParticles();
    this->computeDensityGradients();
    m_timer.stopEvent("Gradient calculation");
  }

  // Resolve secondary emission. We have filled the relevant particles in the transport step -- this can be done
  // either before or after the reactions.
  if (m_emitSecondaryParticlesBeforeReactions) {
    this->barrier();
    m_timer.startEvent("EB particle injection");
    this->fillSecondaryEmissionEB(a_dt);
    this->resolveSecondaryEmissionEB(a_dt);
    m_timer.stopEvent("EB particle injection");
  }

  // Sort the particles and photons per cell so we can call reaction algorithms
  this->barrier();
  m_timer.startEvent("Sort by cell");
  (this->m_ito)->organizeParticlesByCell(ItoSolver::WhichContainer::Bulk);
  this->sortPhotonsByCell(McPhoto::WhichContainer::Bulk);
  this->sortPhotonsByCell(McPhoto::WhichContainer::Source);
  m_timer.stopEvent("Sort by cell");

  // Run the Kinetic Monte Carlo reaction kernels.
  this->barrier();
  m_timer.startEvent("Reaction network");
  this->advanceReactionNetwork(a_dt);
  m_timer.stopEvent("Reaction network");

  // Sort particles per patch.
  this->barrier();
  m_timer.startEvent("Sort by patch");
  (this->m_ito)->organizeParticlesByPatch(ItoSolver::WhichContainer::Bulk);
  this->sortPhotonsByPatch(McPhoto::WhichContainer::Bulk);
  this->sortPhotonsByPatch(McPhoto::WhichContainer::Source);
  m_timer.stopEvent("Sort by patch");

  // Resolve secondary emission. We have filled the relevant particles in the transport step -- this can be done
  // either before or after the reactions.
  if (!m_emitSecondaryParticlesBeforeReactions) {
    this->barrier();
    m_timer.startEvent("EB particle injection");
    this->fillSecondaryEmissionEB(a_dt);
    this->resolveSecondaryEmissionEB(a_dt);
    m_timer.stopEvent("EB particle injection");
  }

  // Remove particles that are inside the EB -- this is not a part of the algorithm, just a safety measure to make sure
  // we don't do things with particles that lie inside the EB.
  this->barrier();
  m_timer.startEvent("Remove covered");
  this->removeCoveredParticles(SpeciesSubset::AllMobileOrDiffusive, EBRepresentation::Discrete, this->m_toleranceEB);
  m_timer.stopEvent("Remove covered");

  // Clear BC data holders.
  for (auto solverIt = (this->m_ito)->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->clear(ItoSolver::WhichContainer::EB);
    solverIt()->clear(ItoSolver::WhichContainer::Domain);
  }

  // Prepare for the next time step
  this->barrier();
  m_timer.startEvent("Post-compute v");
  this->computeDriftVelocities();
  m_timer.stopEvent("Post-compute v");

  this->barrier();
  m_timer.startEvent("Post-compute D");
  this->computeDiffusionCoefficients();
  m_timer.stopEvent("Post-compute D");

  this->computePhysicsDt();

  if ((this->m_profile)) {
    m_timer.eventReport(pout(), false);
  }

  m_timer.clear();

  // Compute the maximum field (in Townsend).
  this->m_maxReducedField = this->computeMaxReducedElectricField(this->m_plasmaPhase);

  return a_dt;
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::preRegrid(const int a_lmin, const int a_oldFinestLevel) noexcept
{
  CH_TIME("ItoKMCGodunovStepper::preRegrid");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCGodunovStepper::preRegrid" << endl;
  }

  const int numItoSpecies    = (this->m_physics)->getNumItoSpecies();
  const int numCdrSpecies    = (this->m_physics)->getNumCdrSpecies();
  const int numPlasmaSpecies = (this->m_physics)->getNumPlasmaSpecies();
  const int numPhotonSpecies = (this->m_physics)->getNumPhotonSpecies();

  ItoKMCStepper<I, C, R, F>::preRegrid(a_lmin, a_oldFinestLevel);

  for (auto solverIt = (this->m_ito)->iterator(); solverIt.ok(); ++solverIt) {
    const int idx = solverIt.index();

    m_conductivityParticles[idx]->preRegrid();
    m_irregularParticles[idx]->preRegrid();
    m_rhoDaggerParticles[idx]->preRegrid();
  }

  this->m_amr->allocate(m_scratchSemiImplicitRhoCDR, this->m_fluidRealm, this->m_plasmaPhase, 1);
  this->m_amr->allocate(m_scratchSemiImplicitConductivityCDR, this->m_fluidRealm, this->m_plasmaPhase, 1);

  DataOps::copy(m_scratchSemiImplicitRhoCDR, m_semiImplicitRhoCDR);
  DataOps::copy(m_scratchSemiImplicitConductivityCDR, m_semiImplicitConductivityCDR);

  // Release unnecessary storage.
  for (int i = 0; i < numCdrSpecies; i++) {
    m_cdrDivD[i].clear();
  }

  m_semiImplicitRhoCDR.clear();
  m_semiImplicitConductivityCDR.clear();
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::regrid(const int a_lmin,
                                         const int a_oldFinestLevel,
                                         const int a_newFinestLevel) noexcept
{
  CH_TIME("ItoKMCGodunovStepper::regrid");
  if (this->m_verbosity > 5) {
    pout() << "ItoKMCGodunovStepper::regrid" << endl;
  }

  m_timer = Timer("ItoKMCGodunovStepper::regrid");

  // A special flag for aborting the simulation if the user did NOT put checkpoint particles in the checkpoint
  // file but still try to regrid on restart.
  if (!m_canRegridOnRestart) {
    const std::string baseErr = "ItoKMCGodunovStepper::regrid -- can't regrid because";
    const std::string err1    = "checkpoint file does not contain particles. Set Driver.initial_regrids=0";

    pout() << baseErr + err1 << endl;

    MayDay::Error((baseErr + err1).c_str());
  }

  // Regrid solvers
  m_timer.startEvent("Regrid ItoSolver");
  (this->m_ito)->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  if (this->m_timeStep == 0) {
    // Necessary because first time step is not semi-implicit, so it will use the densities from the
    // Ito and CDr solvers. However,  ItoSolver does not re-deposit the particles during regrid, so we
    // enforce it here.
    (this->m_ito)->depositParticles();
  }
  m_timer.stopEvent("Regrid ItoSolver");

  m_timer.startEvent("Regrid CdrSolver");
  (this->m_cdr)->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_timer.stopEvent("Regrid CdrSolver");

  m_timer.startEvent("Regrid FieldSolver");
  (this->m_fieldSolver)->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_timer.stopEvent("Regrid FieldSolver");

  m_timer.startEvent("Regrid RTE");
  (this->m_rte)->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_timer.stopEvent("Regrid RTE");

  m_timer.startEvent("Regrid SurfaceODESolver");
  this->m_sigmaSolver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_timer.stopEvent("Regrid SurfaceODESolver");

  // Allocate internal memory for ItoKMCGodunovStepper now....
  m_timer.startEvent("Allocate internals");
  this->allocateInternals();
  m_timer.stopEvent("Allocate internals");

  // We need to remap/regrid the stored data required for the semi-implicit update as well.
  m_timer.startEvent("Remap algorithm-particles");
  for (auto solverIt = (this->m_ito)->iterator(); solverIt.ok(); ++solverIt) {
    const int idx = solverIt.index();
    (this->m_amr)->remapToNewGrids(*m_rhoDaggerParticles[idx], a_lmin, a_newFinestLevel);
    (this->m_amr)->remapToNewGrids(*m_conductivityParticles[idx], a_lmin, a_newFinestLevel);
    (this->m_amr)->remapToNewGrids(*m_irregularParticles[idx], a_lmin, a_newFinestLevel);
  }
  m_timer.stopEvent("Remap algorithm-particles");

  // Also regrid the space charge density contribution from the CDR equations
  this->m_amr->interpToNewGrids(m_semiImplicitRhoCDR,
                                m_scratchSemiImplicitRhoCDR,
                                this->m_plasmaPhase,
                                a_lmin,
                                a_oldFinestLevel,
                                a_newFinestLevel,
                                EBCoarseToFineInterp::Type::ConservativePWC);

  this->m_amr->interpToNewGrids(m_semiImplicitConductivityCDR,
                                m_scratchSemiImplicitConductivityCDR,
                                this->m_plasmaPhase,
                                a_lmin,
                                a_oldFinestLevel,
                                a_newFinestLevel,
                                EBCoarseToFineInterp::Type::ConservativePWC);

  // Set up the field solver with standard coefficients or with
  // modified coefficients if we are reusing data from the last time step.
  m_timer.startEvent("Setup field solver");
  (this->m_fieldSolver)->setupSolver();
  this->computeConductivities(m_conductivityParticles);
  this->setupSemiImplicitPoisson(this->m_prevDt);
  m_timer.stopEvent("Setup field solver");

  // Solve the Poisson equation.
  m_timer.startEvent("Solve Poisson");
  if (this->m_timeStep == 0) {
    this->computeSpaceChargeDensity();
  }
  else {
    this->depositPointParticles(m_rhoDaggerParticles, SpeciesSubset::All);
    this->computeSemiImplicitRho();
  }

  const bool converged = this->solvePoisson();

  if (!converged) {
    const std::string errMsg = "ItoKMCGodunovStepper::regrid - Poisson solve did not converge after regrid";

    pout() << errMsg << endl;

    if (this->m_abortOnFailure) {
      MayDay::Error(errMsg.c_str());
    }
  }
  m_timer.stopEvent("Solve Poisson");

  // Regrid superparticles.
  if (this->m_regridSuperparticles) {
    m_timer.startEvent("Make superparticles");
    (this->m_ito)->organizeParticlesByCell(ItoSolver::WhichContainer::Bulk);
    (this->m_ito)->makeSuperparticles(ItoSolver::WhichContainer::Bulk, (this->m_particlesPerCell));
    (this->m_ito)->organizeParticlesByPatch(ItoSolver::WhichContainer::Bulk);
    m_timer.stopEvent("Make superparticles");
  }

  // Now let Ihe ito solver deposit its actual particles... In the above it deposit m_rhoDaggerParticles.
  m_timer.startEvent("Deposit particles");
  (this->m_ito)->depositParticles();
  m_timer.stopEvent("Deposit particles");

  // Recompute new velocities and diffusion coefficients
  m_timer.startEvent("Prepare next step");
  this->computeDiffusionCoefficients();
  this->computeDriftVelocities();
  m_timer.stopEvent("Prepare next step");

  m_timer.eventReport(pout(), false);

  // No reason to have these lying around.
  m_scratchSemiImplicitRhoCDR.clear();
  m_scratchSemiImplicitConductivityCDR.clear();

  // Fill the neutral density on the mesh
  this->fillNeutralDensity();
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::setOldPositions() noexcept
{
  CH_TIME("ItoKMCGodunovStepper::setOldPositions");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::setOldPositions" << endl;
  }

  for (auto solverIt = (this->m_ito)->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>& solver = solverIt();

    for (int lvl = 0; lvl <= (this->m_amr)->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = (this->m_amr)->getGrids((this->m_particleRealm))[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      auto& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        ParticleSoA<ItoParticle>& leaf = particles[din];

        double* const pos[SpaceDim] = {D_DECL(leaf.positionColumn(0), leaf.positionColumn(1), leaf.positionColumn(2))};
        double* const oldPos[SpaceDim] = {D_DECL(leaf.template column<&ItoParticle::old_x>(),
                                                 leaf.template column<&ItoParticle::old_y>(),
                                                 leaf.template column<&ItoParticle::old_z>())};

        ParticleLoops::loop(leaf, [&](const std::size_t i) {
          D_DECL(oldPos[0][i] = pos[0][i], oldPos[1][i] = pos[1][i], oldPos[2][i] = pos[2][i]);
        });
      }
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::remapPointParticles(Vector<RefCountedPtr<ParticleContainer<NoPayload>>>& a_particles,
                                                      const SpeciesSubset a_subset) noexcept
{
  CH_TIME("ItoKMCGodunovStepper::remapPointParticles");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::remapPointParticles" << endl;
  }

  for (auto solverIt = (this->m_ito)->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = species->getChargeNumber() != 0;

    switch (a_subset) {
    case SpeciesSubset::All: {
      a_particles[idx]->remap();

      break;
    }
    case SpeciesSubset::AllMobile: {
      if (mobile) {
        a_particles[idx]->remap();
      }

      break;
    }
    case SpeciesSubset::AllDiffusive: {
      if (diffusive) {
        a_particles[idx]->remap();
      }

      break;
    }
    case SpeciesSubset::AllMobileOrDiffusive: {
      if (mobile || diffusive) {
        a_particles[idx]->remap();
      }

      break;
    }
    case SpeciesSubset::AllMobileAndDiffusive: {
      if (mobile && diffusive) {
        a_particles[idx]->remap();
      }

      break;
    }
    case SpeciesSubset::Charged: {
      if (charged) {
        a_particles[idx]->remap();
      }

      break;
    }
    case SpeciesSubset::ChargedMobile: {
      if (charged && mobile) {
        a_particles[idx]->remap();
      }

      break;
    }
    case SpeciesSubset::ChargedDiffusive: {
      if (charged && diffusive) {
        a_particles[idx]->remap();
      }

      break;
    }
    case SpeciesSubset::ChargedMobileOrDiffusive: {
      if (charged && (mobile || diffusive)) {
        a_particles[idx]->remap();
      }

      break;
    }
    case SpeciesSubset::ChargedMobileAndDiffusive: {
      if (charged && (mobile && diffusive)) {
        a_particles[idx]->remap();
      }

      break;
    }
    case SpeciesSubset::Stationary: {
      if (!mobile && !diffusive) {
        a_particles[idx]->remap();
      }

      break;
    }
    default: {
      MayDay::Abort("ItoKMCGodunovStepper::remapPointParticles - logic bust");

      break;
    }
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::depositPointParticles(
  const Vector<RefCountedPtr<ParticleContainer<NoPayload>>>& a_particles,
  const SpeciesSubset                                        a_subset) noexcept
{
  CH_TIME("ItoKMCGodunovStepper::depositPointParticles");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::depositPointParticles" << endl;
  }

  for (auto solverIt = (this->m_ito)->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = species->getChargeNumber() != 0;

    switch (a_subset) {
    case SpeciesSubset::All: {
      depositPointParticlesLikeSolver(this->m_amr,
                                      this->m_particleRealm,
                                      this->m_plasmaPhase,
                                      solver,
                                      solver->getPhi(),
                                      *a_particles[idx]);

      break;
    }
    case SpeciesSubset::AllMobile: {
      if (mobile) {
        depositPointParticlesLikeSolver(this->m_amr,
                                        this->m_particleRealm,
                                        this->m_plasmaPhase,
                                        solver,
                                        solver->getPhi(),
                                        *a_particles[idx]);
      }

      break;
    }
    case SpeciesSubset::AllDiffusive: {
      if (diffusive) {
        depositPointParticlesLikeSolver(this->m_amr,
                                        this->m_particleRealm,
                                        this->m_plasmaPhase,
                                        solver,
                                        solver->getPhi(),
                                        *a_particles[idx]);
      }

      break;
    }
    case SpeciesSubset::AllMobileOrDiffusive: {
      if (mobile || diffusive) {
        depositPointParticlesLikeSolver(this->m_amr,
                                        this->m_particleRealm,
                                        this->m_plasmaPhase,
                                        solver,
                                        solver->getPhi(),
                                        *a_particles[idx]);
      }
      break;
    }
    case SpeciesSubset::AllMobileAndDiffusive: {
      if (mobile && diffusive) {
        depositPointParticlesLikeSolver(this->m_amr,
                                        this->m_particleRealm,
                                        this->m_plasmaPhase,
                                        solver,
                                        solver->getPhi(),
                                        *a_particles[idx]);
      }
      break;
    }
    case SpeciesSubset::Charged: {
      if (charged) {
        depositPointParticlesLikeSolver(this->m_amr,
                                        this->m_particleRealm,
                                        this->m_plasmaPhase,
                                        solver,
                                        solver->getPhi(),
                                        *a_particles[idx]);
      }
      break;
    }
    case SpeciesSubset::ChargedMobile: {
      if (charged && mobile) {
        depositPointParticlesLikeSolver(this->m_amr,
                                        this->m_particleRealm,
                                        this->m_plasmaPhase,
                                        solver,
                                        solver->getPhi(),
                                        *a_particles[idx]);
      }
      break;
    }
    case SpeciesSubset::ChargedDiffusive: {
      if (charged && diffusive) {
        depositPointParticlesLikeSolver(this->m_amr,
                                        this->m_particleRealm,
                                        this->m_plasmaPhase,
                                        solver,
                                        solver->getPhi(),
                                        *a_particles[idx]);
      }

      break;
    }
    case SpeciesSubset::ChargedMobileOrDiffusive: {
      if (charged && (mobile || diffusive)) {
        depositPointParticlesLikeSolver(this->m_amr,
                                        this->m_particleRealm,
                                        this->m_plasmaPhase,
                                        solver,
                                        solver->getPhi(),
                                        *a_particles[idx]);
      }

      break;
    }
    case SpeciesSubset::ChargedMobileAndDiffusive: {
      if (charged && (mobile && diffusive)) {
        depositPointParticlesLikeSolver(this->m_amr,
                                        this->m_particleRealm,
                                        this->m_plasmaPhase,
                                        solver,
                                        solver->getPhi(),
                                        *a_particles[idx]);
      }

      break;
    }
    case SpeciesSubset::Stationary: {
      if (!mobile && !diffusive) {
        depositPointParticlesLikeSolver(this->m_amr,
                                        this->m_particleRealm,
                                        this->m_plasmaPhase,
                                        solver,
                                        solver->getPhi(),
                                        *a_particles[idx]);
      }

      break;
    }
    default: {
      MayDay::Abort("ItoKMCGodunovStepper::depositPointParticles - logic bust");

      break;
    }
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::clearPointParticles(
  const Vector<RefCountedPtr<ParticleContainer<NoPayload>>>& a_particles,
  const SpeciesSubset                                        a_subset) noexcept
{
  CH_TIME("ItoKMCGodunovStepper::clearPointParticles");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::clearPointParticles" << endl;
  }

  for (auto solverIt = (this->m_ito)->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = species->getChargeNumber() != 0;

    switch (a_subset) {
    case SpeciesSubset::All: {
      a_particles[idx]->clearParticles();

      break;
    }
    case SpeciesSubset::AllMobile: {
      if (mobile) {
        a_particles[idx]->clearParticles();
      }

      break;
    }
    case SpeciesSubset::AllDiffusive: {
      if (diffusive) {
        a_particles[idx]->clearParticles();
      }

      break;
    }
    case SpeciesSubset::AllMobileOrDiffusive: {
      if (mobile || diffusive) {
        a_particles[idx]->clearParticles();
      }

      break;
    }
    case SpeciesSubset::AllMobileAndDiffusive: {
      if (mobile && diffusive) {
        a_particles[idx]->clearParticles();
      }

      break;
    }
    case SpeciesSubset::Charged: {
      if (charged) {
        a_particles[idx]->clearParticles();
      }

      break;
    }
    case SpeciesSubset::ChargedMobile: {
      if (charged && mobile) {
        a_particles[idx]->clearParticles();
      }

      break;
    }
    case SpeciesSubset::ChargedDiffusive: {
      if (charged && diffusive) {
        a_particles[idx]->clearParticles();
      }

      break;
    }
    case SpeciesSubset::ChargedMobileOrDiffusive: {
      if (charged && (mobile || diffusive)) {
        a_particles[idx]->clearParticles();
      }

      break;
    }
    case SpeciesSubset::ChargedMobileAndDiffusive: {
      if (charged && (mobile && diffusive)) {
        a_particles[idx]->clearParticles();
      }

      break;
    }
    case SpeciesSubset::Stationary: {
      if (!mobile && !diffusive) {
        a_particles[idx]->clearParticles();
      }

      break;
    }
    default: {
      MayDay::Abort("ItoKMCGodunovStepper::clearPointParticles - logic bust");

      break;
    }
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::computeConductivities(
  const Vector<RefCountedPtr<ParticleContainer<NoPayload>>>& a_particles) noexcept
{
  CH_TIME("ItoKMCGodunovStepper::computeConductivities");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::computeConductivities" << endl;
  }

  this->computeCellConductivity((this->m_conductivityCell), a_particles);
  this->computeFaceConductivity();
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::computeCellConductivity(
  EBAMRCellData&                                             a_conductivityCell,
  const Vector<RefCountedPtr<ParticleContainer<NoPayload>>>& a_particles) noexcept
{
  CH_TIME("ItoKMCGodunovStepper::computeCellConductivity(EBAMRCellData, ParticleContainer)");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::computeCellConductivity(EBAMRCellData, ParticleContainer)" << endl;
  }

  DataOps::setValue(a_conductivityCell, 0.0);

  // Contribution from Ito solvers.
  for (auto solverIt = (this->m_ito)->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();
    const int Z   = species->getChargeNumber();

    if (Z != 0 && solver->isMobile()) {
      // Deposit on the particle realm.
      DataOps::setValue(this->m_particleScratch1, 0.0);
      depositPointParticlesLikeSolver(this->m_amr,
                                      this->m_particleRealm,
                                      this->m_plasmaPhase,
                                      solver,
                                      this->m_particleScratch1,
                                      *a_particles[idx]);

      // Copy to fluid realm and add to total conductivity.
      (this->m_amr)->copyData(this->m_fluidScratch1, this->m_particleScratch1);
      DataOps::incr(a_conductivityCell, this->m_fluidScratch1, 1.0 * std::abs(Z));
    }
  }

  // Contribution from CDR solvers.
  DataOps::setValue(m_semiImplicitConductivityCDR, 0.0);
  for (auto solverIt = (this->m_cdr)->iterator(); solverIt.ok(); ++solverIt) {
    const RefCountedPtr<CdrSolver>&  solver  = solverIt();
    const RefCountedPtr<CdrSpecies>& species = solver->getSpecies();

    const int index = solverIt.index();
    const int Z     = species->getChargeNumber();

    if (Z != 0 && solver->isMobile()) {
      const EBAMRCellData& phi = solver->getPhi();
      const EBAMRCellData& mu  = this->m_cdrMobilities[index];

      DataOps::copy(this->m_fluidScratch1, phi);
      DataOps::multiply(this->m_fluidScratch1, mu);

      DataOps::incr(a_conductivityCell, this->m_fluidScratch1, 1.0 * std::abs(Z));
    }
  }

  // Conductivity is mobility * weight * Q
  DataOps::scale(a_conductivityCell, Units::Qe);

  // Coarsen, update ghost cells and interpolate to centroids. User can ask for bi/tri-linear filtering.
  (this->m_amr)->arithmeticAverage(a_conductivityCell, this->m_fluidRealm, (this->m_plasmaPhase));
  (this->m_amr)->interpGhostPwl(a_conductivityCell, this->m_fluidRealm, (this->m_plasmaPhase));

  // User can choose to filter the conductivity. Need to sync after each smoothing.
  if (m_condFilterNum > 0 && m_condFilterMaxStride > 0) {
    for (int i = 0; i < m_condFilterNum; i++) {
      for (int curStride = 1; curStride <= m_condFilterMaxStride; curStride++) {
        DataOps::filterSmooth(a_conductivityCell, m_condFilterAlpha, curStride, true);

        (this->m_amr)->arithmeticAverage(a_conductivityCell, this->m_fluidRealm, (this->m_plasmaPhase));
        (this->m_amr)->interpGhostPwl(a_conductivityCell, this->m_fluidRealm, (this->m_plasmaPhase));
      }
    }
  }

  (this->m_amr)->interpToCentroids(a_conductivityCell, this->m_fluidRealm, (this->m_plasmaPhase));
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::computeFaceConductivity() noexcept
{
  CH_TIME("ItoKMCGodunovStepper::computeFaceConductivity");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::computeFaceConductivity" << endl;
  }

  DataOps::setValue((this->m_conductivityFace), 0.0);
  DataOps::setValue((this->m_conductivityEB), 0.0);

  // Average the cell-centered conductivity to faces. Note that this includes one "ghost face", which we need
  // because the multigrid solver will interpolate face-centered conductivities to face centroids.
  const Average  average  = Average::Arithmetic;
  const int      tanGhost = 1;
  const Interval interv(0, 0);

  DataOps::averageCellToFace(
    (this->m_conductivityFace),
    (this->m_conductivityCell),
    (this->m_amr)->getDomains(),
    tanGhost,
    interv,
    interv,
    average,
    (this->m_amr)->getFaceIteratorWithTangentialGhosts(this->m_fluidRealm, this->m_plasmaPhase));

  // Set the EB conductivity.
  DataOps::incr((this->m_conductivityEB),
                (this->m_conductivityCell),
                1.0,
                (this->m_amr)->getVofIterator(this->m_fluidRealm, this->m_plasmaPhase));
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::computeSemiImplicitRho() noexcept
{
  CH_TIMERS("ItoKMCGodunovStepper::computeSemiImplicitRho");
  CH_TIMER("ItoKMCGodunovStepper::computeSemiImplicitRho::plasma_phase", t1);
  CH_TIMER("ItoKMCGodunovStepper::computeSemiImplicitRho::solid_phase", t2);
  CH_TIMER("ItoKMCGodunovStepper::computeSemiImplicitRho::filter", t3);
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::computeSemiImplicitRho" << endl;
  }

  // Soft requirement?
  CH_assert(this->m_plasmaPhase == phase::gas);

  const RefCountedPtr<MultiFluidIndexSpace>& mfis        = (this->m_computationalGeometry)->getMfIndexSpace();
  const Vector<Dielectric>&                  dielectrics = (this->m_computationalGeometry)->getDielectrics();

  const bool hasDielectrics = (mfis->numPhases() > 1) && (dielectrics.size() > 0);

  MFAMRCellData& rho    = this->m_fieldSolver->getRho();
  EBAMRCellData  rhoGas = (this->m_amr)->alias(phase::gas, rho);
  EBAMRCellData  rhoSolid;

  if (hasDielectrics) {
    rhoSolid = (this->m_amr)->alias(phase::solid, rho);
  }

  DataOps::setValue(rho, 0.0);

  // Contribution from Ito solvers to the gas-side space charge density.
  CH_START(t1);
  for (auto solverIt = this->m_ito->iterator(); solverIt.ok(); ++solverIt) {
    const RefCountedPtr<ItoSolver>&  solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();
    const int                        Z       = species->getChargeNumber();

    if (Z != 0) {
      (this->m_amr)->copyData(this->m_fluidScratch1, solver->getPhi());

      DataOps::incr(rhoGas, this->m_fluidScratch1, 1.0 * Z * Units::Qe);
    }
  }

  // Contribution from CDR equations to the gas-side space charge density.
  DataOps::incr(rhoGas, m_semiImplicitRhoCDR, 1.0);
  CH_STOP(t1);

  // Contribution from gas-side particles that diffused into EBs. This might be necessary near dielectric EBs as a comparatively small
  // correction in the space charge density. There should be no contribution from the CDR solvers because the diffusion-only divergence
  // term was computed using homogeneous Neumann boundary conditions.
  if (hasDielectrics) {
    CH_START(t2);

    EBAMRCellData particleScratch;
    EBAMRCellData fluidScratch;

    (this->m_amr)->allocate(particleScratch, this->m_particleRealm, phase::solid, 1);
    (this->m_amr)->allocate(fluidScratch, this->m_fluidRealm, phase::solid, 1);

    for (auto solverIt = (this->m_ito)->iterator(); solverIt.ok(); ++solverIt) {
      const RefCountedPtr<ItoSolver>&  solver  = solverIt();
      const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();
      const int                        Z       = species->getChargeNumber();

      if (Z != 0) {
        DataOps::setValue(particleScratch, 0.0);
        DataOps::setValue(fluidScratch, 0.0);

        (this->m_amr)
          ->depositWeight(particleScratch,
                          this->m_particleRealm,
                          phase::solid,
                          *m_rhoDaggerParticles[solverIt.index()],
                          DepositionType::CIC,
                          CoarseFineDeposition::Halo,
                          true);

        (this->m_amr)->copyData(fluidScratch, particleScratch);

        DataOps::incr(rhoSolid, fluidScratch, 1.0 * Z * Units::Qe);
      }
    }

    CH_STOP(t2);
  }

  // Sync across levels
  this->m_amr->arithmeticAverage(rho, this->m_fluidRealm);
  this->m_amr->interpGhostPwl(rho, this->m_fluidRealm);

  // User can choose to filter the space charge density. Need to sync after each smoothing.
  if (m_rhoFilterNum > 0 && m_rhoFilterMaxStride > 0) {
    CH_START(t3);
    for (int i = 0; i < m_rhoFilterNum; i++) {
      for (int curStride = 1; curStride <= m_rhoFilterMaxStride; curStride++) {

        DataOps::filterSmooth(rhoGas, m_rhoFilterAlpha, curStride, true);

        this->m_amr->arithmeticAverage(rhoGas, this->m_fluidRealm, this->m_plasmaPhase);
        this->m_amr->interpGhost(rhoGas, this->m_fluidRealm, this->m_plasmaPhase);
      }
    }
    CH_STOP(t3);
  }

  // Put data on centroids.
  this->m_amr->interpToCentroids(rhoGas, this->m_fluidRealm, phase::gas);
  if (hasDielectrics) {
    this->m_amr->interpToCentroids(rhoSolid, this->m_fluidRealm, phase::solid);
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::setupSemiImplicitPoisson(const Real a_dt) noexcept
{
  CH_TIME("ItoKMCGodunovStepper::setupSemiImplicitPoisson");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::setupSemiImplicitPoisson" << endl;
  }

  // Set coefficients as usual
  (this->m_fieldSolver)->setPermittivities();

  // Get the permittivities
  MFAMRCellData& permCell = (this->m_fieldSolver)->getPermittivityCell();
  MFAMRFluxData& permFace = (this->m_fieldSolver)->getPermittivityFace();
  MFAMRIVData&   permEB   = (this->m_fieldSolver)->getPermittivityEB();

  // Get handles to the gas-phase permittivities.
  EBAMRFluxData permFaceGas = (this->m_amr)->alias((this->m_plasmaPhase), permFace);
  EBAMRIVData   permEBGas   = (this->m_amr)->alias((this->m_plasmaPhase), permEB);

  // Increment the field solver permittivities by a_factor*sigma. After this, the "permittivities" are
  // given by epsr + a_factor*sigma
  DataOps::incr(permFaceGas, (this->m_conductivityFace), a_dt / Units::eps0);
  DataOps::incr(permEBGas,
                (this->m_conductivityEB),
                a_dt / Units::eps0,
                (this->m_amr)->getVofIterator(this->m_fluidRealm, this->m_plasmaPhase));

  // Coarsen coefficients.
  (this->m_amr)->arithmeticAverage(permFaceGas, this->m_fluidRealm, (this->m_plasmaPhase));
  (this->m_amr)->arithmeticAverage(permEBGas, this->m_fluidRealm, (this->m_plasmaPhase));

  // Set up the solver with the "permittivities"
  (this->m_fieldSolver)->setSolverPermittivities(permCell, permFace, permEB);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::removeCoveredPointParticles(
  Vector<RefCountedPtr<ParticleContainer<NoPayload>>>& a_particles,
  const EBRepresentation                               a_representation,
  const Real                                           a_tolerance) const noexcept
{
  CH_TIME("ItoKMCGodunovStepper::removeCoveredPointParticles");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::removeCoveredPointParticles" << endl;
  }

  for (int i = 0; i < a_particles.size(); i++) {
    if (a_particles[i] != nullptr) {
      ParticleContainer<NoPayload>& particles = *a_particles[i];

      switch (a_representation) {
      case EBRepresentation::Discrete: {
        (this->m_amr)->removeCoveredParticlesDiscrete(particles, (this->m_plasmaPhase), a_tolerance);

        break;
      }
      case EBRepresentation::ImplicitFunction: {
        (this->m_amr)->removeCoveredParticlesIF(particles, (this->m_plasmaPhase), a_tolerance);

        break;
      }
      case EBRepresentation::Voxel: {
        (this->m_amr)->removeCoveredParticlesVoxels(particles, (this->m_plasmaPhase));

        break;
      }
      default: {
        MayDay::Error("ItoKMCGodunovStepper::removeCoveredParticles - logic bust");
      }
      }
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::copyConductivityParticles(
  Vector<RefCountedPtr<ParticleContainer<NoPayload>>>& a_conductivityParticles) noexcept
{
  CH_TIME("ItoKMCGodunovStepper::copyConductivityParticles");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::copyConductivityParticles" << endl;
  }

  // Clear particles first.
  this->clearPointParticles(a_conductivityParticles, SpeciesSubset::All);

  for (auto solverIt = (this->m_ito)->iterator(); solverIt.ok(); ++solverIt) {
    const RefCountedPtr<ItoSolver>&  solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();
    const int Z   = species->getChargeNumber();

    if (Z != 0 && solver->isMobile()) {
      const ParticleContainer<ItoParticle>& solverParticles = solver->getParticles(ItoSolver::WhichContainer::Bulk);

      for (int lvl = 0; lvl <= (this->m_amr)->getFinestLevel(); lvl++) {
        const DisjointBoxLayout& dbl = (this->m_amr)->getGrids((this->m_particleRealm))[lvl];
        const DataIterator&      dit = dbl.dataIterator();

        const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
        for (int mybox = 0; mybox < nbox; mybox++) {
          const DataIndex& din = dit[mybox];

          const ParticleSoA<ItoParticle>& leaf = solverParticles[lvl][din];

          ParticleSoA<NoPayload>& pointParticles = (*a_conductivityParticles[idx])[lvl][din];
          ParticleSoA<NoPayload>& irregParticles = (*m_irregularParticles[idx])[lvl][din];

          for (std::size_t i = 0; i < leaf.size(); i++) {
            const RealVect pos      = leaf.position(i);
            const Real     weight   = leaf.weight(i);
            const Real     mobility = leaf.template get<&ItoParticle::mobility>(i);

            pointParticles.append(pos, weight * mobility);
          }

          pointParticles.catenate(irregParticles);
        }
      }
    }
  }
}

template <typename I, typename C, typename R, typename F>
bool
ItoKMCGodunovStepper<I, C, R, F>::solvePoisson() noexcept
{
  CH_TIME("ItoKMCGodunovStepper::solvePoisson()");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::solvePoisson()" << endl;
  }

  // Solve the Poisson equation and compute the cell-centered electric field.
  MFAMRCellData& phi   = this->m_fieldSolver->getPotential();
  MFAMRCellData& rho   = this->m_fieldSolver->getRho();
  EBAMRIVData&   sigma = this->m_sigmaSolver->getPhi();

  const bool converged = (this->m_fieldSolver)->solve(phi, rho, sigma, false);

  (this->m_fieldSolver)->computeElectricField();

  // Copy the electric field to appropriate data holders and perform center-to-centroid
  // interpolation.
  EBAMRCellData E;
  (this->m_amr)->allocatePointer(E, this->m_fluidRealm);
  (this->m_amr)->alias(E, this->m_plasmaPhase, (this->m_fieldSolver)->getElectricField());

  // Fluid realm
  (this->m_amr)->copyData(this->m_electricFieldFluid, E);
  (this->m_amr)->conservativeAverage(this->m_electricFieldFluid, this->m_fluidRealm, this->m_plasmaPhase);
  (this->m_amr)->interpGhostPwl(this->m_electricFieldFluid, this->m_fluidRealm, this->m_plasmaPhase);
  (this->m_amr)->interpToCentroids(this->m_electricFieldFluid, this->m_fluidRealm, this->m_plasmaPhase);

  // Particle realm
  (this->m_amr)->copyData(this->m_electricFieldParticle, E);
  (this->m_amr)->conservativeAverage(this->m_electricFieldParticle, this->m_particleRealm, this->m_plasmaPhase);
  (this->m_amr)->interpGhostPwl(this->m_electricFieldParticle, this->m_particleRealm, this->m_plasmaPhase);
  (this->m_amr)->interpToCentroids(this->m_electricFieldParticle, this->m_particleRealm, this->m_plasmaPhase);

  return converged;
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::advanceEulerMaruyama(const Real a_dt) noexcept
{
  CH_TIME("ItoKMCGodunovStepper::advanceEulerMaruyama");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::advanceEulerMaruyama" << endl;
  }

  // Store X^k positions.
  this->setOldPositions();

  // Diffuse the particles. This copies onto m_rhoDaggerParticles and stores the hop on the full particles. We need
  // to remap the particles species that made a diffusion hop.
  this->barrier();
  m_timer.startEvent("Diffuse particles");
  this->diffuseParticlesEulerMaruyama(m_rhoDaggerParticles, a_dt);
  this->remapPointParticles(m_rhoDaggerParticles, SpeciesSubset::ChargedDiffusive);
  m_timer.stopEvent("Diffuse particles");

  // Perform the diffusive CDR advance.
  this->barrier();
  m_timer.startEvent("Diffuse CDR");
  this->computeDiffusionTermCDR(m_semiImplicitRhoCDR, a_dt);
  m_timer.stopEvent("Diffuse CDR");

  // Compute the conductivity on the mesh. This deposits q_e * Z * w * mu on the mesh.
  this->barrier();
  m_timer.startEvent("Compute conductivities");
  this->copyConductivityParticles(m_conductivityParticles);
  this->computeConductivities(m_conductivityParticles);
  m_timer.stopEvent("Compute conductivities");

  // Set up the semi-implicit Poisson solver with the computed conductivities.
  this->barrier();
  m_timer.startEvent("Setup Poisson");
  this->setupSemiImplicitPoisson(a_dt);
  m_timer.stopEvent("Setup Poisson");

  // Compute space charge density arising from the new particle positions X^k + sqrt(2*D*dt)*W. Only need to
  // do the diffusive and charged species.
  this->barrier();
  m_timer.startEvent("Deposit point particles");
  this->depositPointParticles(m_rhoDaggerParticles, SpeciesSubset::Charged);
  this->computeSemiImplicitRho();
  m_timer.stopEvent("Deposit point particles");

  // Solve the semi-implicit Poisson equation.
  this->barrier();
  m_timer.startEvent("Solve Poisson");
  const bool converged = this->solvePoisson();
  if (!converged) {
    const std::string errMsg = "ItoKMCGodunovStepper::advanceEulerMaruyama - Poisson solve did not converge";

    pout() << errMsg << endl;

    if (this->m_abortOnFailure) {
      MayDay::Error(errMsg.c_str());
    }
  }
  m_timer.stopEvent("Solve Poisson");

  // Recompute velocities with the new electric field. This interpolates the velocities to the current particle positions, i.e.
  // we compute V^(k+1)(X^k) = mu^k * E^(k+1)(X^k)
  this->barrier();
  m_timer.startEvent("Step-compute v");
#if 1 // This is what the algorithm says.
  this->setCdrVelocityFunctions();
  this->setItoVelocityFunctions();
  (this->m_ito)->interpolateVelocities();
  this->multiplyCdrVelocitiesByMobilities();
#else // Have to use this for LEA - need to debug.
  this->computeDriftVelocities();
#endif
  m_timer.stopEvent("Step-compute v");

  // Finalize the Euler-Maruyama update.
  this->barrier();
  m_timer.startEvent("Euler-Maruyama step");
  this->stepEulerMaruyamaParticles(a_dt);
  this->remapParticles(SpeciesSubset::AllMobileOrDiffusive);
  this->stepEulerMaruyamaCDR(a_dt);
  m_timer.stopEvent("Euler-Maruyama step");
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::diffuseParticlesEulerMaruyama(
  Vector<RefCountedPtr<ParticleContainer<NoPayload>>>& a_rhoDaggerParticles,
  const Real                                           a_dt) noexcept
{
  CH_TIME("ItoKMCGodunovStepper::diffuseParticlesEulerMaruyama");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::diffuseParticlesEulerMaruyama" << endl;
  }

  this->clearPointParticles(a_rhoDaggerParticles, SpeciesSubset::All);

  for (auto solverIt = (this->m_ito)->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const int  Z         = species->getChargeNumber();

    const auto& diffusionFunction = (this->m_physics)->getItoDiffusionFunctions()[idx];

    for (int lvl = 0; lvl <= (this->m_amr)->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = (this->m_amr)->getGrids((this->m_particleRealm))[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      auto& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        ParticleSoA<ItoParticle>& leaf           = particles[din];
        ParticleSoA<NoPayload>&   pointParticles = (*a_rhoDaggerParticles[idx])[lvl][din];

        for (std::size_t i = 0; i < leaf.size(); i++) {
          const Real     weight = leaf.weight(i);
          const RealVect pos    = leaf.position(i);

          // Compute a particle hop and store it on the run-time storage. The DiffusionFunction reads payload
          // columns (diffusion/velocity), so materialize the payload via gather().
          const ItoParticle p   = leaf.gather(i);
          const RealVect    hop = diffusive ? diffusionFunction(p, a_dt) : RealVect::Zero;

          D_DECL(leaf.template get<&ItoParticle::scratch_x>(i) = static_cast<ParticleReal>(hop[0]),
                 leaf.template get<&ItoParticle::scratch_y>(i) = static_cast<ParticleReal>(hop[1]),
                 leaf.template get<&ItoParticle::scratch_z>(i) = static_cast<ParticleReal>(hop[2]));

          if (Z != 0) {
            pointParticles.append(pos + hop, weight);
          }
        }
      }
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::computeDiffusionTermCDR(EBAMRCellData& a_semiImplicitRhoCDR, const Real a_dt) noexcept
{
  CH_TIME("ItoKMCGodunovStepper::diffuseCDREulerMaruyama");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::diffuseCDREulerMaruyama" << endl;
  }

  DataOps::setValue(a_semiImplicitRhoCDR, 0.0);

  for (auto solverIt = this->m_cdr->iterator(); solverIt.ok(); ++solverIt) {
    const RefCountedPtr<CdrSolver>&  solver  = solverIt();
    const RefCountedPtr<CdrSpecies>& species = solver->getSpecies();

    const int index = solverIt.index();
    const int Z     = species->getChargeNumber();

    const EBAMRCellData& phi = solver->getPhi();

    // Compute finite-volume approximation to div(D*grad(phi))
    if (solver->isDiffusive()) {
      solver->computeDivD(m_cdrDivD[index], solver->getPhi(), false, false, false);
    }

    // Update the space charge density arising from the update phi^dagger = phi^k + dt * div(D*grad(phi^k))
    if (Z != 0) {
      DataOps::incr(a_semiImplicitRhoCDR, phi, 1.0 * Z);
      if (solver->isDiffusive()) {
        DataOps::incr(a_semiImplicitRhoCDR, m_cdrDivD[index], 1.0 * Z * a_dt);
      }
    }
  }

  DataOps::scale(a_semiImplicitRhoCDR, Units::Qe);

  this->m_amr->arithmeticAverage(a_semiImplicitRhoCDR, this->m_fluidRealm, this->m_plasmaPhase);
  this->m_amr->interpGhostPwl(a_semiImplicitRhoCDR, this->m_fluidRealm, this->m_plasmaPhase);
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::stepEulerMaruyamaParticles(const Real a_dt) noexcept
{
  CH_TIME("ItoKMCGodunovStepper::stepEulerMaruyamaParticles");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::stepEulerMaruyamaParticles" << endl;
  }

  for (auto solverIt = (this->m_ito)->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>& solver = solverIt();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    const Real f = mobile ? a_dt : 0.0;
    const Real g = diffusive ? 1.0 : 0.0;

    if (mobile || diffusive) {
      for (int lvl = 0; lvl <= (this->m_amr)->getFinestLevel(); lvl++) {
        const DisjointBoxLayout& dbl = (this->m_amr)->getGrids((this->m_particleRealm))[lvl];
        const DataIterator&      dit = dbl.dataIterator();

        auto& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

        const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
        for (int mybox = 0; mybox < nbox; mybox++) {
          const DataIndex& din = dit[mybox];

          ParticleSoA<ItoParticle>& leaf = particles[din];

          double* const pos[SpaceDim] = {
            D_DECL(leaf.positionColumn(0), leaf.positionColumn(1), leaf.positionColumn(2))};
          double* const       oldPos[SpaceDim] = {D_DECL(leaf.template column<&ItoParticle::old_x>(),
                                                   leaf.template column<&ItoParticle::old_y>(),
                                                   leaf.template column<&ItoParticle::old_z>())};
          ParticleReal* const vel[SpaceDim]    = {D_DECL(leaf.template column<&ItoParticle::vx>(),
                                                      leaf.template column<&ItoParticle::vy>(),
                                                      leaf.template column<&ItoParticle::vz>())};
          ParticleReal* const hop[SpaceDim]    = {D_DECL(leaf.template column<&ItoParticle::scratch_x>(),
                                                      leaf.template column<&ItoParticle::scratch_y>(),
                                                      leaf.template column<&ItoParticle::scratch_z>())};

          ParticleLoops::loop(leaf, [&](const std::size_t i) {
            // Add in the diffusion hop and advective contribution.
            for (int dir = 0; dir < SpaceDim; dir++) {
              pos[dir][i] = oldPos[dir][i] + f * vel[dir][i] + g * hop[dir][i];
            }
          });
        }
      }
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::stepEulerMaruyamaCDR(const Real a_dt) noexcept
{
  CH_TIME("ItoKMCGodunovStepper::stepEulerMaruyamaCDR");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::stepEulerMaruyamaCDR" << endl;
  }

  for (auto solverIt = (this->m_cdr)->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<CdrSolver>& solver = solverIt();

    const int index = solverIt.index();

    EBAMRCellData& phi = solver->getPhi();

    this->m_amr->conservativeAverage(phi, this->m_fluidRealm, this->m_plasmaPhase);
    this->m_amr->interpGhostPwl(phi, this->m_fluidRealm, this->m_plasmaPhase);

    // Add in the advective term -- BC comes later.
    if (solver->isMobile()) {
      DataOps::setValue(solver->getEbFlux(), 0.0);

      // Compute the FV approximation to the advective term and do the Euler advance. If the underlying
      // CDR solver is a CTU solver, we also add in the transverse terms.
      solver->computeDivF(this->m_fluidScratch1, phi, 0.5 * a_dt, false, true, true);

      DataOps::incr(phi, this->m_fluidScratch1, -a_dt);
    }

    // Add in the diffusion term.
    if (solver->isDiffusive()) {
      DataOps::incr(phi, this->m_cdrDivD[index], a_dt);
    }

    DataOps::floor(phi, 0.0, this->m_amr->getVofIterator(this->m_fluidRealm, this->m_plasmaPhase));
  }

  this->coarsenCDRSolvers();
}

#ifdef CH_USE_HDF5
template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::writeCheckpointHeader(HDF5HeaderData& a_header) const noexcept
{
  CH_TIME("ItoKMCGodunovStepper::writeCheckpointHeader");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::writeCheckpointHeader" << endl;
  }

  a_header.m_real["prev_dt"]             = this->m_prevDt;
  a_header.m_real["physics_dt"]          = this->m_physicsDt;
  a_header.m_int["checkpoint_particles"] = m_writeCheckpointParticles ? 1 : 0;
}
#endif

#ifdef CH_USE_HDF5
template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::readCheckpointHeader(HDF5HeaderData& a_header) noexcept
{
  CH_TIME("ItoKMCGodunovStepper::readCheckpointHeader");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::readCheckpointHeader" << endl;
  }

  this->m_prevDt    = a_header.m_real["prev_dt"];
  this->m_physicsDt = a_header.m_real["physics_dt"];

  m_readCheckpointParticles = (a_header.m_int["checkpoint_particles"] != 0) ? true : false;
  m_canRegridOnRestart      = m_readCheckpointParticles;
}
#endif

#ifdef CH_USE_HDF5
template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const noexcept
{
  CH_TIME("ItoKMCGodunovStepper::writeCheckpointData");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::writeCheckpointData" << endl;
  }

  ItoKMCStepper<I, C, R, F>::writeCheckpointData(a_handle, a_lvl);

  // Write the point-particles.
  if (m_writeCheckpointParticles) {
    for (int i = 0; i < (this->m_physics)->getNumItoSpecies(); i++) {
      const std::string identifierSigma = "ItoKMCGodunovStepper::conductivityParticles_" + std::to_string(i);
      const std::string identifierRho   = "ItoKMCGodunovStepper::spaceChargeParticles_" + std::to_string(i);

      const ParticleContainer<NoPayload>& conductivityParticles = *m_conductivityParticles[i];
      const ParticleContainer<NoPayload>& rhoDaggerParticles    = *m_rhoDaggerParticles[i];

      DischargeIO::writeCheckParticlesToHDF(a_handle, conductivityParticles[a_lvl], identifierSigma);
      DischargeIO::writeCheckParticlesToHDF(a_handle, rhoDaggerParticles[a_lvl], identifierRho);
    }
  }

  // Write data required for the semi-implicit regrid.
  if (this->m_physics->getNumCdrSpecies() > 0) {
    write(a_handle, *m_semiImplicitRhoCDR[a_lvl], "ItoKMCGodunovStepper::semiImplicitRhoCDR");
    write(a_handle, *m_semiImplicitConductivityCDR[a_lvl], "ItoKMCGodunovStepper::semiImplicitConductivityCDR");
  }
}
#endif

#ifdef CH_USE_HDF5
template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::readCheckpointData(HDF5Handle& a_handle, const int a_lvl) noexcept
{
  CH_TIME("ItoKMCGodunovStepper::readCheckpointData");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::readCheckpointData" << endl;
  }

  ItoKMCStepper<I, C, R, F>::readCheckpointData(a_handle, a_lvl);

  // Write the point-particles.
  if (m_readCheckpointParticles) {
    for (int i = 0; i < (this->m_physics)->getNumItoSpecies(); i++) {
      const std::string identifierSigma = "ItoKMCGodunovStepper::conductivityParticles_" + std::to_string(i);
      const std::string identifierRho   = "ItoKMCGodunovStepper::spaceChargeParticles_" + std::to_string(i);

      ParticleContainer<NoPayload>& conductivityParticles = *m_conductivityParticles[i];
      ParticleContainer<NoPayload>& rhoDaggerParticles    = *m_rhoDaggerParticles[i];

      DischargeIO::readCheckParticlesFromHDF(a_handle, conductivityParticles[a_lvl], identifierSigma);
      DischargeIO::readCheckParticlesFromHDF(a_handle, rhoDaggerParticles[a_lvl], identifierRho);
    }
  }

  // Read in data that is required for the semi-implicit regrid.
  if (this->m_physics->getNumCdrSpecies() > 0) {
    const Interval interv(0, 0);

    read<EBCellFAB>(a_handle,
                    *m_semiImplicitRhoCDR[a_lvl],
                    "ItoKMCGodunovStepper::semiImplicitRhoCDR",
                    this->m_amr->getGrids(this->m_fluidRealm)[a_lvl],
                    interv,
                    false);

    read<EBCellFAB>(a_handle,
                    *m_semiImplicitConductivityCDR[a_lvl],
                    "ItoKMCGodunovStepper::semiImplicitConductivityCDR",
                    this->m_amr->getGrids(this->m_fluidRealm)[a_lvl],
                    interv,
                    false);
  }
}
#endif

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::postPlot() noexcept
{
  CH_TIME("ItoKMCGodunovStepper::postPlot");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::postPlot" << endl;
  }

  this->m_physicsPlotVariables.clear();

  this->plotParticles();
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::plotParticles() const noexcept
{
  CH_TIME("ItoKMCGodunovStepper::plotParticles");
  if (this->m_verbosity > 2) {
    pout() << this->m_name + "::plotParticles" << endl;
  }

  bool plotParticles = false;

  ParmParse pp(this->m_name.c_str());

  pp.query("plot_particles", plotParticles);

  if (plotParticles) {

    for (auto solverIt = (this->m_ito)->iterator(); solverIt.ok(); ++solverIt) {
      const RefCountedPtr<ItoSolver>&       solver    = solverIt();
      const ParticleContainer<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk);

      // Create the output folder
      std::string cmd     = "mkdir -p particles/" + solver->getName();
      int         success = 0;
      if (procID() == 0) {
        success = system(cmd.c_str());
      }

      if (success != 0) {
        MayDay::Error("ItoKMCGodunovStepper::plotParticles - could not create 'particles' directory");
      }

      // Set plot file name
      const std::string prefix = "./particles/" + solver->getName() + "/" + solver->getName();
      char              fileChar[1000];
      sprintf(fileChar, "%s.step%07d.%dd.h5part", prefix.c_str(), this->m_timeStep, SpaceDim);

      // Plot the particles directly from the live SoA container. The output datasets (id, position, weight
      // and the selected ItoParticle payload columns) are derived from ParticleTraits<ItoParticle>::h5PartColumns.
      DischargeIO::writeH5Part(std::string(fileChar), particles, this->m_amr->getProbLo(), this->m_time);
    }
  }
}

template <typename I, typename C, typename R, typename F>
void
ItoKMCGodunovStepper<I, C, R, F>::postInitialize() noexcept
{
  CH_TIME("ItoKMCGodunovStepper::postInitialize");
  if (this->m_verbosity > 5) {
    pout() << this->m_name + "::postInitialize" << endl;
  }

  ItoKMCStepper<I, C, R, F>::postInitialize();

  bool debugTestNNMerge = false;

  ParmParse pp(this->m_name.c_str());
  pp.query("debug_test_nn_merge", debugTestNNMerge);

  if (!debugTestNNMerge) {
    return;
  }

  int  debugNumSamplesPerRank      = 200;
  int  debugParticlesPerCellThresh = 4;
  Real debugMaxCellDistanceReal    = -1.0;
  int  debugNumMergeRounds         = 1;
  bool debugIterateToConvergence   = false;
  int  debugMaxFallbackCandidates  = 0;
  int  debugSubPartitionThresh     = std::numeric_limits<int>::max();
  pp.query("debug_nn_merge_samples_per_rank", debugNumSamplesPerRank);
  pp.query("debug_nn_merge_particles_per_cell_thresh", debugParticlesPerCellThresh);
  pp.query("debug_nn_merge_max_cell_distance", debugMaxCellDistanceReal);
  pp.query("debug_nn_merge_num_rounds", debugNumMergeRounds);
  pp.query("debug_nn_merge_iterate_to_convergence", debugIterateToConvergence);
  pp.query("debug_nn_merge_max_fallback_candidates", debugMaxFallbackCandidates);
  pp.query("debug_nn_merge_sub_partition_thresh", debugSubPartitionThresh);

  const std::optional<int> debugMaxCellDistance = (debugMaxCellDistanceReal < 0.0)
                                                    ? std::nullopt
                                                    : std::optional<int>(static_cast<int>(debugMaxCellDistanceReal));

  const RefCountedPtr<BaseIF>& implicitFunction = this->m_computationalGeometry->getImplicitFunction(
    this->m_plasmaPhase);

  // Rejects a merge whose weighted-centroid position would fall inside the embedded boundary --
  // BaseIF::value() >= 0 means inside the solid (Chombo convention: negative is inside the
  // fluid). If this geometry has no analytic implicit function (e.g. a DCEL/STL geometry), there
  // is nothing to check against, so every position is accepted rather than silently rejecting
  // everything.
  auto isPositionValid = [&implicitFunction](const RealVect& a_pos) -> bool {
    if (implicitFunction.isNull()) {
      return true;
    }
    return implicitFunction->value(a_pos) < 0.0;
  };

  // Rank-namespaced, monotonically increasing id counter -- shared between the initial sample
  // (so newly-drawn particles get real ids, since ParticleManagement::drawBoxParticles/
  // ParticleContainer::addParticlesDestructive leave id at ParticleSoA::s_invalidID) and the
  // merge's own fresh-id allocation, so the two can never collide. See
  // CD_NearestNeighborParticleMerge.H's IDAllocator docs -- there is no existing general-purpose
  // global-id allocation scheme elsewhere in the particle infrastructure to reuse instead.
  constexpr ParticleID rankStride = 100000000LL;
  ParticleID           nextID     = static_cast<ParticleID>(procID()) * rankStride;
  auto                 allocateID = [&nextID]() -> ParticleID {
    return nextID++;
  };

  for (auto solverIt = (this->m_ito)->iterator(); solverIt.ok(); ++solverIt) {
    const RefCountedPtr<ItoSolver>& solver = solverIt();

    ParticleContainer<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk);

    // Sample new ItoParticle particles uniformly across the domain.
    ParticleSoA<ItoParticle> drawn;
    ParticleManagement::drawBoxParticles(drawn,
                                         debugNumSamplesPerRank,
                                         this->m_amr->getProbLo(),
                                         this->m_amr->getProbHi());
    particles.addParticlesDestructive(drawn);

    // drawBoxParticles() samples uniformly over the domain's bounding box with no EB awareness --
    // discard whatever landed inside the solid before anything else touches these particles.
    // Reusing ItoSolver's own, already-existing removeCoveredParticles() rather than
    // reimplementing an EB check here (see also isPositionValid below, which handles the
    // SEPARATE concern of a merge RESULT's position, not the initial sample).
    solver->removeCoveredParticles(particles, EBRepresentation::ImplicitFunction, this->m_toleranceEB);

    // Assign every particle a real, globally-unique id (both pre-existing particles that may
    // still carry ParticleSoA::s_invalidID, and the newly-sampled ones above).
    for (int lvl = 0; lvl <= this->m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = this->m_amr->getGrids(solver->getRealm())[lvl];
      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        ParticleSoA<ItoParticle>& leaf = particles[lvl][dit()];
        for (std::size_t i = 0; i < leaf.size(); i++) {
          leaf.particleID(i) = allocateID();
          leaf.rankID(i)     = procID();
        }
      }
    }

    // Phase 0 standalone check for the new hybrid spatial-index design (see TODO.md's "Performance"
    // section): buildNNCellBuckets() is not consumed by anything yet, so exercise it directly here
    // against the freshly-sampled particles. blockSize/haloBlockBox are constructed inline (not via
    // MergeBlockGrid, which doesn't exist until Phase 1) purely for this standalone check.
    {
      const int          blockSize       = this->m_amr->getMinBlockSize();
      unsigned long long totalChecked    = 0;
      unsigned long long totalMismatched = 0;

      for (int lvl = 0; lvl <= this->m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout& dbl      = this->m_amr->getGrids(solver->getRealm())[lvl];
        const Real               dxScalar = this->m_amr->getDx()[lvl];
        const RealVect           dx       = RealVect(D_DECL(dxScalar, dxScalar, dxScalar));
        const RealVect           probLo   = this->m_amr->getProbLo();

        for (DataIterator dit(dbl); dit.ok(); ++dit) {
          ParticleSoA<ItoParticle>& leaf = particles[lvl][dit()];

          std::vector<ParticleManagement::NNMergeParticle<Real>> localValid;
          std::vector<ParticleManagement::NNMergeParticle<Real>> ghosts; // empty -- no ghost-fill yet
          localValid.reserve(leaf.size());
          for (std::size_t i = 0; i < leaf.size(); i++) {
            if (leaf.isGhost(i)) {
              continue;
            }
            ParticleManagement::NNMergeParticle<Real> p;
            p.position = leaf.position(i);
            p.weight   = leaf.weight(i);
            p.globalID = leaf.particleID(i);
            localValid.push_back(p);
          }

          const Box haloBlockBox = grow(coarsen(dbl[dit()], blockSize), 1);

          ParticleManagement::NNCellBuckets<Real> buckets;
          ParticleManagement::buildNNCellBuckets(localValid, ghosts, probLo, dx, blockSize, haloBlockBox, buckets);

          if (buckets.sorted.size() != localValid.size()) {
            MayDay::Error("ItoKMCGodunovStepper::postInitialize -- buildNNCellBuckets() Phase 0 "
                          "check failed: total indexed count != input count");
          }

          unsigned long long blockCountSum = 0;
          for (const int c : buckets.blockCounts) {
            blockCountSum += static_cast<unsigned long long>(c);
          }
          if (blockCountSum != localValid.size()) {
            MayDay::Error("ItoKMCGodunovStepper::postInitialize -- buildNNCellBuckets() Phase 0 "
                          "check failed: blockCounts does not sum to input count");
          }

          for (std::size_t i = 0; i < localValid.size(); i++) {
            totalChecked++;

            const IntVect cell = ParticleManagement::cellKeyOf(localValid[i].position, probLo, dx);
            const auto    it   = buckets.cellOffsets.find(cell);
            if (it == buckets.cellOffsets.end()) {
              totalMismatched++;
              continue;
            }

            bool found = false;
            for (std::size_t s = it->second.first; s < it->second.second; s++) {
              if (buckets.sorted[s].isLocal && buckets.sorted[s].idx == i) {
                found = true;
                break;
              }
            }
            if (!found) {
              totalMismatched++;
            }
          }
        }
      }

      const unsigned long long totalCheckedGlobal    = ParallelOps::sum(totalChecked);
      const unsigned long long totalMismatchedGlobal = ParallelOps::sum(totalMismatched);
      if (procID() == 0) {
        pout() << this->m_name + "::postInitialize -- buildNNCellBuckets() Phase 0 check for solver '"
               << solver->getName() << "': " << totalCheckedGlobal << " particles checked, " << totalMismatchedGlobal
               << " mismatched (should be 0)" << endl;
      }
      if (totalMismatchedGlobal != 0) {
        MayDay::Error("ItoKMCGodunovStepper::postInitialize -- buildNNCellBuckets() Phase 0 check "
                      "found mismatched particles");
      }
    }

    // Phase 1 standalone check: MergeBlockGrid has no consumer yet -- verify its own geometric
    // self-consistency directly (haloBlockBox == grow(blockBox, 1), coarsen(cellBox, blockSize) ==
    // blockBox, and that every one of this rank's own patches has a reachable entry).
    {
      const int blockSize = this->m_amr->getMinBlockSize();

      for (int lvl = 0; lvl <= this->m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout&             dbl  = this->m_amr->getGrids(solver->getRealm())[lvl];
        const RefCountedPtr<MergeBlockGrid>& grid = this->m_amr->getMergeBlockGrid(solver->getRealm())[lvl];

        if (grid->getBlockSize() != blockSize) {
          MayDay::Error("ItoKMCGodunovStepper::postInitialize -- MergeBlockGrid Phase 1 check "
                        "failed: getBlockSize() mismatch");
        }

        for (DataIterator dit(dbl); dit.ok(); ++dit) {
          const MergeBlockGrid::PatchBlocks& pb = grid->getPatchBlocks(dit());

          if (pb.cellBox != dbl[dit()]) {
            MayDay::Error("ItoKMCGodunovStepper::postInitialize -- MergeBlockGrid Phase 1 check "
                          "failed: cellBox != this patch's own box");
          }
          if (pb.blockBox != coarsen(pb.cellBox, blockSize)) {
            MayDay::Error("ItoKMCGodunovStepper::postInitialize -- MergeBlockGrid Phase 1 check "
                          "failed: blockBox != coarsen(cellBox, blockSize)");
          }
          if (pb.haloBlockBox != grow(pb.blockBox, 1)) {
            MayDay::Error("ItoKMCGodunovStepper::postInitialize -- MergeBlockGrid Phase 1 check "
                          "failed: haloBlockBox != grow(blockBox, 1)");
          }
        }
      }

      if (procID() == 0) {
        pout() << this->m_name + "::postInitialize -- MergeBlockGrid Phase 1 check for solver '" << solver->getName()
               << "': geometric self-consistency OK" << endl;
      }
    }

    // Per-level valid particle counts (each ParallelOps::sum() is a collective -- called
    // unconditionally by every rank, same iteration count everywhere, since m_amr->getFinestLevel()
    // is identical on every rank). Added specifically so multi-level behavior (e.g. "are merges
    // happening on every level, or only the coarsest") is directly visible, not just a global total.
    auto countPerLevel = [this,
                          &solver](ParticleContainer<ItoParticle>& a_particles) -> std::vector<unsigned long long> {
      std::vector<unsigned long long> counts(this->m_amr->getFinestLevel() + 1, 0);
      for (int lvl = 0; lvl <= this->m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout& dbl        = this->m_amr->getGrids(solver->getRealm())[lvl];
        unsigned long long       localCount = 0;
        for (DataIterator dit(dbl); dit.ok(); ++dit) {
          ParticleSoA<ItoParticle>& leaf = a_particles[lvl][dit()];
          for (std::size_t i = 0; i < leaf.size(); i++) {
            if (!leaf.isGhost(i)) {
              localCount++;
            }
          }
        }
        counts[lvl] = ParallelOps::sum(localCount);
      }
      return counts;
    };

    auto printPerLevel = [](const std::string& a_label, const std::vector<unsigned long long>& a_counts) {
      if (procID() == 0) {
        pout() << "  " << a_label << " = [";
        for (std::size_t lvl = 0; lvl < a_counts.size(); lvl++) {
          pout() << a_counts[lvl] << (lvl + 1 < a_counts.size() ? ", " : "");
        }
        pout() << "]" << endl;
      }
    };

    const std::vector<unsigned long long> countsBefore = countPerLevel(particles);

    const unsigned long long nBefore = particles.getNumberOfValidParticlesGlobal();
    Real                     wBefore = 0.0;
    for (int lvl = 0; lvl <= this->m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = this->m_amr->getGrids(solver->getRealm())[lvl];
      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        ParticleSoA<ItoParticle>& leaf = particles[lvl][dit()];
        for (std::size_t i = 0; i < leaf.size(); i++) {
          if (!leaf.isGhost(i)) {
            wBefore += leaf.weight(i);
          }
        }
      }
    }
    wBefore = ParallelOps::sum(wBefore);

    // Write the ORIGINAL (pre-merge) particles.
    {
      std::string cmd     = "mkdir -p particles/" + solver->getName();
      int         success = 0;
      if (procID() == 0) {
        success = system(cmd.c_str());
      }
      if (success != 0) {
        MayDay::Error("ItoKMCGodunovStepper::postInitialize -- could not create 'particles' directory");
      }

      const std::string prefix = "./particles/" + solver->getName() + "/" + solver->getName() + ".nn_merge_original";
      DischargeIO::writeH5Part(prefix + "." + std::to_string(SpaceDim) + "d.h5part",
                               particles,
                               this->m_amr->getProbLo(),
                               this->m_time);
    }

    // Ghost-only diagnostic snapshot: copy the sampled (pre-merge) particles into a SEPARATE
    // container -- never touching `particles`, which the merge below still needs completely
    // untouched -- fill ITS ghosts using the already-registered mask, strip out everything that
    // is NOT a ghost, and write what's left. Lets the actual ghost region (and how far a ghost
    // travels, e.g. across a coarse-fine boundary) be inspected directly in the visualization,
    // independent of any merge decision.
    {
      ParticleContainer<ItoParticle> ghostOnlyParticles;
      ghostOnlyParticles.define(this->m_amr->getGrids(solver->getRealm()),
                                this->m_amr->getDomains(),
                                this->m_amr->getDx(),
                                this->m_amr->getRefinementRatios(),
                                this->m_amr->getProbLo(),
                                this->m_amr->getMinBlockSize(),
                                this->m_amr->getLevelTiles(solver->getRealm()),
                                this->m_amr->getFinestLevel(),
                                solver->getRealm(),
                                nullptr);

      // Direct per-patch copy -- ghostOnlyParticles shares the EXACT SAME grid layout as
      // `particles` (defined from the same AmrMesh/realm), so this is index-for-index exact
      // (id/rank/position/weight/payload all preserved). Deliberately NOT
      // addParticlesDestructive(), which would remap particles and lose that direct
      // per-patch correspondence.
      for (int lvl = 0; lvl <= this->m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout& dbl = this->m_amr->getGrids(solver->getRealm())[lvl];
        for (DataIterator dit(dbl); dit.ok(); ++dit) {
          ghostOnlyParticles[lvl][dit()].append(particles[lvl][dit()]);
        }
      }

      ghostOnlyParticles.fillGhostParticles(this->m_amr->getParticleGhostMask(solver->getRealm(), 1),
                                            this->m_amr->getParticleGhostMaskCoarToFine(solver->getRealm(), 1),
                                            this->m_amr->getParticleGhostMaskFineToCoar(solver->getRealm(), 1));

      for (int lvl = 0; lvl <= this->m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout& dbl = this->m_amr->getGrids(solver->getRealm())[lvl];
        for (DataIterator dit(dbl); dit.ok(); ++dit) {
          ParticleSoA<ItoParticle>& leaf = ghostOnlyParticles[lvl][dit()];
          for (std::size_t i = 0; i < leaf.size();) {
            if (!leaf.isGhost(i)) {
              leaf.remove(i);
            }
            else {
              i++;
            }
          }
        }
      }

      const std::string prefix = "./particles/" + solver->getName() + "/" + solver->getName() + ".nn_merge_ghosts_only";
      DischargeIO::writeH5Part(prefix + "." + std::to_string(SpaceDim) + "d.h5part",
                               ghostOnlyParticles,
                               this->m_amr->getProbLo(),
                               this->m_time);
    }

    // The opaque reconciliation payload is just the particle's energy -- mobility/diffusion/
    // velocity are safe to drop across a merge since they are recomputed from field
    // interpolation before they are next used.
    auto gather = [](const ParticleSoA<ItoParticle>& a_leaf, std::size_t a_idx) -> Real {
      return a_leaf.template get<&ItoParticle::energy>(a_idx);
    };
    auto combine = [](const Real a_firstEnergy,
                      const Real a_firstWeight,
                      const Real a_secondEnergy,
                      const Real a_secondWeight) -> Real {
      return (a_firstWeight * a_firstEnergy + a_secondWeight * a_secondEnergy) / (a_firstWeight + a_secondWeight);
    };
    auto scatter = [](ParticleSoA<ItoParticle>& a_leaf, const ParticleManagement::NNMergeParticle<Real>& a_p) -> void {
      ItoParticle payload;
      payload.energy = a_p.payload;
      a_leaf.append(a_p.position, a_p.weight, payload);
      a_leaf.particleID(a_leaf.size() - 1) = a_p.globalID;
      a_leaf.rankID(a_leaf.size() - 1)     = a_p.ownerRank;
    };

    // Cross-level merge count -- see ParticleManagement::mergeNearestNeighborsRound()'s
    // a_crossLevelMergeCount docs. NOT reset by the library between calls, so this single counter
    // (passed by pointer to every round below) accumulates over the whole debug test. This is a
    // PER-RANK count (each rank only counts the merges IT committed, whether as a trivial-tier
    // participant or as judge) -- like nBefore/nAfter/wBefore/wAfter, it needs an explicit
    // ParallelOps::sum() (a collective, called by every rank) to get the true cross-rank total
    // before printing; printing the raw per-rank value would silently miss activity on every rank
    // but rank 0.
    unsigned long long crossLevelMergeCount = 0;

    // Per-merge diagnostic: recompute the weighted centroid from (first, second)'s OWN recorded
    // position/weight and compare against the merged particle's actual position. These must
    // match exactly (up to floating-point rounding) by construction -- any real deviation would
    // mean the merged particle's data was NOT actually derived from the two ids recorded as its
    // parents, a genuine bug distinct from "the nearest-neighbor search legitimately found a
    // distant candidate" (the centroid of two distant points still lies exactly on the segment
    // between them, so a large maxPartnerDistance with a tiny maxCenterlineDeviation means the
    // MERGE math is fine and any surprising-looking merges are a candidate-search/ghost-reach
    // question, not a centroid bug). All three are per-rank local; reduced via ParallelOps::max()
    // before printing.
    Real maxCenterlineDeviation       = 0.0;
    Real maxPartnerDistance           = 0.0;
    Real maxCrossLevelPartnerDistance = 0.0;

    auto onMergeCommitted = [&maxCenterlineDeviation, &maxPartnerDistance, &maxCrossLevelPartnerDistance](
                              const ParticleManagement::NNMergeParticle<Real>& a_first,
                              const ParticleManagement::NNMergeParticle<Real>& a_second,
                              const ParticleManagement::NNMergeParticle<Real>& a_merged) -> void {
      const Real wSum = a_first.weight + a_second.weight;

      RealVect expectedPos(D_DECL(0.0, 0.0, 0.0));
      Real     partnerDist2 = 0.0;
      for (int dir = 0; dir < SpaceDim; dir++) {
        expectedPos[dir] = (a_first.weight * a_first.position[dir] + a_second.weight * a_second.position[dir]) / wSum;
        const Real d     = a_first.position[dir] - a_second.position[dir];
        partnerDist2 += d * d;
      }

      Real deviation2 = 0.0;
      for (int dir = 0; dir < SpaceDim; dir++) {
        const Real d = expectedPos[dir] - a_merged.position[dir];
        deviation2 += d * d;
      }

      const Real deviation   = std::sqrt(deviation2);
      const Real partnerDist = std::sqrt(partnerDist2);

      maxCenterlineDeviation = std::max(maxCenterlineDeviation, deviation);
      maxPartnerDistance     = std::max(maxPartnerDistance, partnerDist);
      if (a_first.level != a_second.level) {
        maxCrossLevelPartnerDistance = std::max(maxCrossLevelPartnerDistance, partnerDist);
      }
    };

    // Repeated rounds are what RULES.md's "Single pass, caller loops for convergence" design
    // expects -- one call does a SINGLE cross-patch exchange and is not guaranteed to drain a
    // crowded cell down to the threshold; debug_nn_merge_num_rounds exists so the effect of
    // repeated calls (with no particle motion between them here, unlike a real timestep) can be
    // observed directly.
    for (int round = 0; round < debugNumMergeRounds; round++) {
      const unsigned long long crossLevelMergeCountBefore = crossLevelMergeCount;

      ParticleManagement::mergeNearestNeighborsRound<ItoParticle, Real>(*this->m_amr,
                                                                        particles,
                                                                        debugParticlesPerCellThresh,
                                                                        gather,
                                                                        combine,
                                                                        scatter,
                                                                        allocateID,
                                                                        debugIterateToConvergence,
                                                                        debugMaxFallbackCandidates,
                                                                        debugMaxCellDistance,
                                                                        isPositionValid,
                                                                        &crossLevelMergeCount,
                                                                        onMergeCommitted,
                                                                        debugSubPartitionThresh);

      // getNumberOfValidParticlesGlobal() and ParallelOps::sum() are both MPI-collectives -- must
      // be called by every rank unconditionally; only the PRINTING is rank-0-only.
      const unsigned long long nRound                          = particles.getNumberOfValidParticlesGlobal();
      const unsigned long long crossLevelMergesThisRoundGlobal = ParallelOps::sum(crossLevelMergeCount -
                                                                                  crossLevelMergeCountBefore);
      if (procID() == 0) {
        pout() << this->m_name + "::postInitialize -- nearest-neighbor merge debug test for solver '"
               << solver->getName() << "', round " << (round + 1) << "/" << debugNumMergeRounds
               << ": particles = " << nRound << ", cross-level merges this round = " << crossLevelMergesThisRoundGlobal
               << endl;
      }
    }

    // Defensive cleanup: isPositionValid above already rejects a merge whose RESULT would fall
    // inside the EB, so this should normally find nothing -- but reuse the same
    // removeCoveredParticles() here too (rather than trust that check alone) so nothing merged
    // ends up sitting inside the solid in the final, written-out state.
    const unsigned long long nBeforeEBCleanup = particles.getNumberOfValidParticlesGlobal();
    solver->removeCoveredParticles(particles, EBRepresentation::ImplicitFunction, this->m_toleranceEB);
    const unsigned long long nAfterEBCleanup = particles.getNumberOfValidParticlesGlobal();
    if (procID() == 0 && nAfterEBCleanup != nBeforeEBCleanup) {
      pout() << this->m_name + "::postInitialize -- nearest-neighbor merge debug test for solver '" << solver->getName()
             << "': removeCoveredParticles() discarded " << (nBeforeEBCleanup - nAfterEBCleanup)
             << " merged particle(s) that ended up inside the EB (isPositionValid should normally have "
                "prevented this)"
             << endl;
    }

    const std::vector<unsigned long long> countsAfter = countPerLevel(particles);

    const unsigned long long nAfter = particles.getNumberOfValidParticlesGlobal();
    Real                     wAfter = 0.0;
    for (int lvl = 0; lvl <= this->m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = this->m_amr->getGrids(solver->getRealm())[lvl];
      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        ParticleSoA<ItoParticle>& leaf = particles[lvl][dit()];
        for (std::size_t i = 0; i < leaf.size(); i++) {
          if (!leaf.isGhost(i)) {
            wAfter += leaf.weight(i);
          }
        }
      }
    }
    wAfter = ParallelOps::sum(wAfter);

    // Write the MERGED (post-merge) particles.
    {
      const std::string prefix = "./particles/" + solver->getName() + "/" + solver->getName() + ".nn_merge_merged";
      DischargeIO::writeH5Part(prefix + "." + std::to_string(SpaceDim) + "d.h5part",
                               particles,
                               this->m_amr->getProbLo(),
                               this->m_time);
    }

    // Collective -- must be called by every rank unconditionally.
    const unsigned long long crossLevelMergeCountGlobal         = ParallelOps::sum(crossLevelMergeCount);
    const Real               maxCenterlineDeviationGlobal       = ParallelOps::max(maxCenterlineDeviation);
    const Real               maxPartnerDistanceGlobal           = ParallelOps::max(maxPartnerDistance);
    const Real               maxCrossLevelPartnerDistanceGlobal = ParallelOps::max(maxCrossLevelPartnerDistance);

    if (procID() == 0) {
      const Real massError = std::abs(wAfter - wBefore);
      pout() << this->m_name + "::postInitialize -- nearest-neighbor merge debug test for solver '" << solver->getName()
             << "':" << endl;
      pout() << "  particles before merge = " << nBefore << endl;
      pout() << "  particles after merge  = " << nAfter << endl;
      pout() << "  difference             = " << (static_cast<long long>(nBefore) - static_cast<long long>(nAfter))
             << endl;
      printPerLevel("particles per level before merge", countsBefore);
      printPerLevel("particles per level after merge ", countsAfter);
      pout() << "  cross-level merges (total, all rounds) = " << crossLevelMergeCountGlobal << endl;
      pout() << "  max centerline deviation (should be ~0) = " << maxCenterlineDeviationGlobal << endl;
      pout() << "  max merge-partner distance              = " << maxPartnerDistanceGlobal << endl;
      pout() << "  max CROSS-LEVEL merge-partner distance  = " << maxCrossLevelPartnerDistanceGlobal << endl;
      pout() << "  total weight before = " << wBefore << endl;
      pout() << "  total weight after  = " << wAfter << endl;
      pout() << "  mass conservation error = " << massError << endl;

      if (massError > 1.E-6 * std::max(wBefore, 1.0)) {
        MayDay::Error(
          "ItoKMCGodunovStepper::postInitialize -- nearest-neighbor merge debug test failed mass conservation");
      }
    }
  }
}

#include <CD_NamespaceFooter.H>

#endif
