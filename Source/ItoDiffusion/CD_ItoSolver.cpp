/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ItoSolver.cpp
 * @brief  Implementation of CD_ItoSolver.H
 * @author Robert Marskar
 */

// Std includes
#include <chrono>
#include <array>
#include <unordered_set>

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>
#include <ParticleIO.H>
#include <EBCellFactory.H>

// Our includes
#include <CD_ItoSolver.H>
#include <CD_Random.H>
#include <CD_DataOps.H>
#include <CD_ParallelOps.H>
#include <CD_ParticleOps.H>
#include <CD_ParticleLoops.H>
#include <CD_DischargeIO.H>
#include <CD_ParticleManagement.H>
#include <CD_NearestNeighborParticleMerge.H>
#include <CD_KDParticleMerge.H>
#include <CD_EBParticleMesh.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

constexpr int ItoSolver::m_comp;
constexpr int ItoSolver::m_nComp;

ItoSolver::ItoSolver()
  : m_checkpointing(WhichCheckpoint::Particles),
    m_mobilityInterp(WhichMobilityInterpolation::Direct),
    m_realm(Realm::primal),
    m_phase(phase::gas),
    m_name("ItoSolver"),
    m_className("ItoSolver"),
    m_verbosity(-1),
    m_deposition(DepositionType::CIC),
    m_coarseFineDeposition(CoarseFineDeposition::Transition),
    m_plotDeposition(DepositionType::CIC)
{
  CH_TIME("ItoSolver::ItoSolver");

  // Default is to not merge particles
  m_particleCellMerger = [](ParticleSoA<ItoMergeParticle>& a_particles, const CellInfo& a_cellInfo, const int a_ppc) {

  };
}

ItoSolver::~ItoSolver()
{
  CH_TIME("ItoSolver::~ItoSolver");
}

std::string
ItoSolver::getName() const
{
  CH_TIME("ItoSolver::getName");

  return m_name;
}

std::string
ItoSolver::getRealm() const
{
  CH_TIME("ItoSolver::getRealm");

  return m_realm;
}

void
ItoSolver::setRealm(const std::string& a_realm)
{
  CH_TIME("ItoSolver::setRealm");

  m_realm = a_realm;
}

void
ItoSolver::setParticleCellMerger(
  const ParticleManagement::ParticleMerger<ItoMergeParticle>& a_particleCellMerger) noexcept
{
  CH_TIME("ItoSolver::setParticleCellMerger");

  m_particleCellMerger = a_particleCellMerger;
}

const RefCountedPtr<ItoSpecies>&
ItoSolver::getSpecies() const
{
  CH_TIME("ItoSolver::getSpecies");
  if (m_verbosity > 5) {
    pout() << m_name + "::getSpecies" << endl;
  }

  CH_assert(!m_species.isNull());

  return m_species;
}

void
ItoSolver::parseOptions()
{
  CH_TIME("ItoSolver::parseOptions");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseOptions" << endl;
  }

  this->parseVerbosity();
  this->parseRNG();
  this->parseTruncation();
  this->parsePlotVariables();
  this->parseDeposition();
  this->parseIntersectionEB();
  this->parseRedistribution();
  this->parseDivergenceComputation();
  this->parseCheckpointing();
  this->parseParticleMerger();
}

void
ItoSolver::parseRuntimeOptions()
{
  CH_TIME("ItoSolver::parseRuntimeOptions");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRuntimeOptions" << endl;
  }

  this->parseVerbosity();
  this->parsePlotVariables();
  this->parseTruncation();
  this->parseDeposition();
  this->parseIntersectionEB();
  this->parseRedistribution();
  this->parseDivergenceComputation();
  this->parseCheckpointing();
  this->parseParticleMerger();
}

void
ItoSolver::parseVerbosity()
{
  CH_TIME("ItoSolver::parseVerbosity");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseVerbosity" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("verbosity", m_verbosity);
}

void
ItoSolver::parseRNG()
{
  CH_TIME("ItoSolver::parseRNG");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRNG" << endl;
  }

  // Seed the RNG
  ParmParse pp(m_className.c_str());

  pp.get("normal_max", m_normalDistributionTruncation);
}

void
ItoSolver::parseTruncation()
{
  CH_TIME("ItoSolver::parseTruncation");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseTruncation" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("normal_max", m_normalDistributionTruncation);
}

void
ItoSolver::parsePlotVariables()
{
  CH_TIME("ItoSolver::parsePlotVariables");
  if (m_verbosity > 5) {
    pout() << m_name + "::parsePlotVariables" << endl;
  }

  m_plotPhi              = false;
  m_plotVelocity         = false;
  m_plotDiffCo           = false;
  m_plotParticles        = false;
  m_plotParticlesEB      = false;
  m_plotParticlesDomain  = false;
  m_plotParticlesSource  = false;
  m_plotParticlesCovered = false;
  m_plotEnergyDensity    = false;
  m_plotAverageEnergy    = false;

  ParmParse pp(m_className.c_str());

  const int           num = pp.countval("plt_vars");
  Vector<std::string> str(num);
  pp.getarr("plt_vars", str, 0, num);

  for (int i = 0; i < num; i++) {
    if (str[i] == "phi") {
      m_plotPhi = true;
    }
    else if (str[i] == "vel") {
      m_plotVelocity = true;
    }
    else if (str[i] == "dco") {
      m_plotDiffCo = true;
    }
    else if (str[i] == "part") {
      m_plotParticles = true;
    }
    else if (str[i] == "eb_part") {
      m_plotParticlesEB = true;
    }
    else if (str[i] == "dom_part") {
      m_plotParticlesDomain = true;
    }
    else if (str[i] == "src_part") {
      m_plotParticlesSource = true;
    }
    else if (str[i] == "covered_part") {
      m_plotParticlesCovered = true;
    }
    else if (str[i] == "energy_density") {
      m_plotEnergyDensity = true;
    }
    else if (str[i] == "average_energy") {
      m_plotAverageEnergy = true;
    }
  }
}

void
ItoSolver::parseDeposition()
{
  CH_TIME("ItoSolver::parseDeposition");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseDeposition" << endl;
  }

  ParmParse   pp(m_className.c_str());
  std::string str;

  // Deposition for particle-mesh operations
  pp.get("deposition", str);
  if (str == "ngp") {
    m_deposition = DepositionType::NGP;
  }
  else if (str == "cic") {
    m_deposition = DepositionType::CIC;
  }
  else if (str == "tsc") {
    m_deposition = DepositionType::TSC;
  }
  else {
    MayDay::Error("ItoSolver::parseDeposition - unknown deposition method requested");
  }

  // Parse coarse-fine strategy
  pp.get("deposition_cf", str);
  if (str == "interp") {
    m_coarseFineDeposition = CoarseFineDeposition::Interp;
  }
  else if (str == "halo") {
    m_coarseFineDeposition = CoarseFineDeposition::Halo;
  }
  else if (str == "halo_ngp") {
    m_coarseFineDeposition = CoarseFineDeposition::HaloNGP;
  }
  else if (str == "transition") {
    m_coarseFineDeposition = CoarseFineDeposition::Transition;
  }
  else {
    MayDay::Error("ItoSolver::parseDeposition - unknown coarse-fine deposition method requested.");
  }

  // Deposition for plotting only
  pp.get("plot_deposition", str);

  if (str == "ngp") {
    m_plotDeposition = DepositionType::NGP;
  }
  else if (str == "cic") {
    m_plotDeposition = DepositionType::CIC;
  }
  else {
    MayDay::Error("ItoSolver::parseDeposition - unknown deposition method requested");
  }

  // Mobility interpolation.
  pp.get("mobility_interp", str);
  if (str == "direct") {
    m_mobilityInterp = WhichMobilityInterpolation::Direct;
  }
  else if (str == "velocity") {
    m_mobilityInterp = WhichMobilityInterpolation::Velocity;
  }
  else {
    MayDay::Abort("ItoSolver::parseDeposition - unknown interpolation method for mobility");
  }

  pp.get("irr_ngp_deposition", m_forceIrregDepositionNGP);
  pp.get("irr_ngp_interp", m_forceIrregInterpolationNGP);
}

void
ItoSolver::parseIntersectionEB()
{
  CH_TIME("ItoSolver::parseIntersectionEB");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseIntersectionEB" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;

  pp.get("intersection_alg", str);
  pp.get("bisect_step", m_bisectionStep);

  if (str == "raycast") {
    m_intersectionAlg = EBIntersection::Raycast;
  }
  else if (str == "bisection") {
    m_intersectionAlg = EBIntersection::Bisection;
  }
  else {
    MayDay::Error("ItoSolver::parseIntersectionEB -- logic bust");
  }
}

void
ItoSolver::parseRedistribution()
{
  CH_TIME("ItoSolver::parseRedistribution");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRedistribution" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("redistribute", m_useRedistribution);
}

void
ItoSolver::parseDivergenceComputation()
{
  CH_TIME("ItoSolver::parseDivergenceComputation");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseDivergenceComputation" << endl;
  }

  ParmParse pp(m_className.c_str());

  pp.get("blend_conservation", m_blendConservation);
}

void
ItoSolver::parseCheckpointing()
{
  CH_TIME("ItoSolver::parseCheckpointing");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseCheckpointing" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("checkpointing", str);
  pp.get("ppc_restart", m_restartPPC);
  if (str == "particles") {
    m_checkpointing = WhichCheckpoint::Particles;
  }
  else if (str == "numbers") {
    m_checkpointing = WhichCheckpoint::Numbers;
  }
  else {
    MayDay::Abort("ItoSolver::parseCheckpointing - unknown checkpointing method requested");
  }
}

void
ItoSolver::parseParticleMerger()
{
  CH_TIME("ItoSolver::parseParticleMerger");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseParticleMerger" << endl;
  }

  ParmParse pp(m_className.c_str());

  std::string str;

  pp.get("merge_algorithm", str);

  // Target particle count per cell, and the regrid-time merge method. Both live here rather than on a
  // time stepper: the solver is where merging happens, so it is what owns the target to merge to and
  // whether to merge at a regrid. A stepper still owns the per-step cadence (merge_interval), which is
  // a property of the advance, not of the solver.
  m_particlesPerCell.resize(pp.countval("particles_per_cell"));
  pp.getarr("particles_per_cell", m_particlesPerCell, 0, m_particlesPerCell.size());

  // regrid_superparticles selects the merge method run inside ItoSolver::regrid(): one of the
  // merge_algorithm selectors (incl. 'none' to turn the merge off), or 'solver' to use merge_algorithm.
  std::string regridMerge;
  pp.get("regrid_superparticles", regridMerge);
  m_regridMergeMethod = (regridMerge == "solver") ? std::optional<ParticleManagement::ParticleMergeMethod>(std::nullopt)
                                                  : std::optional<ParticleManagement::ParticleMergeMethod>(
                                                      ParticleManagement::mergeMethodFromString(regridMerge));

  // The per-cell merger for each cell-based method is built locally in the corresponding
  // mergeXxx() (see mergeLite()); parseParticleMerger only records the method.
  // 'external' leaves the merger to the user (setParticleCellMerger); 'none' does nothing.
  m_mergeMethod = ParticleManagement::mergeMethodFromString(str);

  // The distributed nearest-neighbor pair merges (nn_pair_tree/nn_pair_onecell/nn_pair_hash) have
  // extra tunables. Read and validate them UNCONDITIONALLY, not only when merge_algorithm selects
  // one of them: any of them can also be selected as the regrid-time method
  // (ItoKMCStepper.regrid_superparticles), which drives makeSuperparticles() through
  // m_regridMergeMethod and bypasses m_mergeMethod entirely. Gating the read on m_mergeMethod would
  // silently drop these settings (and skip their validation) for that valid configuration, so they
  // are mandatory inputs like everything else here.
  pp.get("nn_pair_iterate", m_nnPairIterate);
  pp.get("nn_pair_fallback", m_nnPairFallback);
  pp.get("nn_pair_max_cell_dist", m_nnPairMaxCellDistance);
  pp.get("nn_pair_max_rounds", m_nnPairMaxRounds);

  // Same rationale as the nn_pair_* reads above -- kd_carve can also be selected as the
  // regrid-time method, bypassing m_mergeMethod, so this is a mandatory input like everything else
  // here rather than gated on m_mergeMethod.
  pp.get("kd_split_weight_leaf_dx", m_kdSplitWeightLeafDx);

  // A bad value here silently produces a degenerate merge, so fail loudly (in every build, not just
  // DEBUG) rather than assert.
  if (m_nnPairFallback < 0) {
    MayDay::Abort("ItoSolver::parseParticleMerger - 'nn_pair_fallback' must be >= 0");
  }
  if (m_nnPairMaxCellDistance == 0) {
    MayDay::Abort("ItoSolver::parseParticleMerger - 'nn_pair_max_cell_dist' must be >= 1 (or < 0 for unbounded)");
  }
  if (m_nnPairMaxRounds < 1) {
    MayDay::Abort("ItoSolver::parseParticleMerger - 'nn_pair_max_rounds' must be >= 1");
  }
  if (m_kdSplitWeightLeafDx < 0.0) {
    MayDay::Abort("ItoSolver::parseParticleMerger - 'kd_split_weight_leaf_dx' must be >= 0");
  }

  if (m_particlesPerCell.size() < 1) {
    MayDay::Abort("ItoSolver::parseParticleMerger - 'particles_per_cell' must have at least one entry");
  }
  for (int lvl = 0; lvl < m_particlesPerCell.size(); lvl++) {
    if (m_particlesPerCell[lvl] <= 0) {
      MayDay::Abort("ItoSolver::parseParticleMerger - 'particles_per_cell' must be > 0");
    }
  }

  // A non-uniform particles_per_cell is only meaningful to the cell-based methods; the whole-container
  // methods use one crowding threshold for the whole hierarchy and read entry 0, so a per-level vector
  // combined with one of them would be silently truncated. Fail loudly instead, for the same reason the
  // validations above do: a bad value here produces a degenerate merge rather than an error.
  // (initialData() needs no such guard -- generateParticlesFromDensity() takes the vector and resolves
  // it per level.)
  const bool uniformPPC = std::adjacent_find(m_particlesPerCell.stdVector().begin(),
                                             m_particlesPerCell.stdVector().end(),
                                             std::not_equal_to<int>()) == m_particlesPerCell.stdVector().end();

  if (!uniformPPC) {
    auto isCellBased = [](const ParticleManagement::ParticleMergeMethod a_method) -> bool {
      return a_method == ParticleManagement::ParticleMergeMethod::EqualWeightKD ||
             a_method == ParticleManagement::ParticleMergeMethod::Reinitialize ||
             a_method == ParticleManagement::ParticleMergeMethod::ReinitializeBVH ||
             a_method == ParticleManagement::ParticleMergeMethod::NnSfc ||
             a_method == ParticleManagement::ParticleMergeMethod::External;
    };

    if (!isCellBased(m_mergeMethod)) {
      MayDay::Abort("ItoSolver::parseParticleMerger - a per-level 'particles_per_cell' needs a cell-based "
                    "'merge_algorithm' (equal_weight_kd/reinitialize/reinitialize_bvh/nn_sfc/external); the "
                    "whole-container methods use one threshold for the whole hierarchy and would use entry 0 only");
    }
    if (m_regridMergeMethod.has_value() && !isCellBased(m_regridMergeMethod.value()) &&
        m_regridMergeMethod.value() != ParticleManagement::ParticleMergeMethod::None) {
      MayDay::Abort("ItoSolver::parseParticleMerger - a per-level 'particles_per_cell' needs a cell-based "
                    "'regrid_superparticles' (or 'none'); the whole-container methods would use entry 0 only");
    }
  }
}

EBIntersection
ItoSolver::getIntersectionAlgorithm() const noexcept
{
  CH_TIME("ItoSolver::getIntersectionAlgorithm");
  if (m_verbosity > 5) {
    pout() << m_name + "::getIntersectionAlgorithm" << endl;
  }

  return m_intersectionAlg;
}

const Vector<int>&
ItoSolver::getParticlesPerCell() const noexcept
{
  CH_TIME("ItoSolver::getParticlesPerCell");

  return m_particlesPerCell;
}

ParticleContainer<ItoMergeParticle>&
ItoSolver::getRegridParticles() noexcept
{
  CH_TIME("ItoSolver::getRegridParticles");

  return m_regridParticles;
}

unsigned long long
ItoSolver::getNumParticles(const WhichContainer a_whichContainer, const bool a_localOnly) const
{
  CH_TIME("ItoSolver::getNumParticles(WhichContainer, bool)");
  if (m_verbosity > 5) {
    pout() << m_name + "::getNumParticles(WhichContainer, bool)" << endl;
  }

  const ParticleContainer<ItoParticle>& particles = m_particleContainers.at(a_whichContainer);

  unsigned long long N = 0;
  if (a_localOnly) {
    N = particles.getNumberOfValidParticlesLocal();
  }
  else {
    N = particles.getNumberOfValidParticlesGlobal();
  }

  return N;
}

void
ItoSolver::setComputationalGeometry(const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry)
{
  CH_TIME("ItoSolver::setComputationalGeometry");
  if (m_verbosity > 5) {
    pout() << m_name + "::setComputationalGeometry" << endl;
  }

  CH_assert(!a_computationalGeometry.isNull());

  m_computationalGeometry = a_computationalGeometry;
}

void
ItoSolver::setAmr(const RefCountedPtr<AmrMesh>& a_amr)
{
  CH_TIME("ItoSolver::setAmr");
  if (m_verbosity > 5) {
    pout() << m_name + "::setAmr" << endl;
  }

  CH_assert(!a_amr.isNull());

  m_amr = a_amr;
}

void
ItoSolver::registerOperators() const
{
  CH_TIME("ItoSolver::registerOperators");
  if (m_verbosity > 5) {
    pout() << m_name + "::registerOperators" << endl;
  }

  if (m_amr.isNull()) {
    MayDay::Abort("CdrSolver::registerOperators - need to set AmrMesh!");
  }
  else {
    m_amr->registerOperator(s_eb_coar_ave, m_realm, m_phase);
    m_amr->registerOperator(s_eb_fill_patch, m_realm, m_phase);
    m_amr->registerOperator(s_noncons_div, m_realm, m_phase);
    m_amr->registerOperator(s_particle_mesh, m_realm, m_phase);
    m_amr->registerOperator(s_eb_multigrid, m_realm, m_phase);
    if (m_useRedistribution) {
      m_amr->registerOperator(s_eb_redist, m_realm, m_phase);
    }

    // nn_pair_tree/nn_pair_hash both read a particle ghost halo as merge candidates (see
    // mergeNnPairSearch() / ParticleManagement::mergeNearestNeighborsTree()/
    // mergeNearestNeighborsHash()). Its width is the merge distance nn_pair_max_cell_dist (>= 1), or
    // 1 when that is unbounded or neither merger is used. Register it here, alongside the other
    // operators; the same width must be filled and passed to the round -- shared by both backends
    // since they use the identical ghost-width formula.
    const int nnPairSearchGhostWidth = (m_nnPairMaxCellDistance < 0) ? 1 : m_nnPairMaxCellDistance;
    m_amr->registerParticleGhostMask(m_realm, nnPairSearchGhostWidth);

    // nn_pair_onecell's merge distance is structurally fixed at 1 (see
    // the nn_pair_onecell merge / ParticleManagement::mergeNearestNeighborsOneCell()), so it
    // always needs a width-1 mask regardless of nn_pair_max_cell_dist; register it too (a no-op if
    // nnPairSearchGhostWidth is already 1).
    m_amr->registerParticleGhostMask(m_realm, 1);

    // kd_carve and kd_skin_nn both hardcode their ghost width to 1, unlike nn_pair_tree/nn_pair_hash's
    // nn_pair_max_cell_dist-driven width. mergeKDCarve() likewise bounds a mergeable leaf's
    // own extent at a fixed one cell width, which the width-1 fill here already matches exactly:
    // buildKDQuotaLeaves() can never see a particle farther away than this fill provides, so a
    // looser leaf bound could not be honoured even if one existed. Registered unconditionally here, same as
    // the two calls above -- a no-op given the width-1 registration already present for
    // nn_pair_onecell, kept explicit so this doesn't silently depend on that other registration
    // continuing to exist.
    m_amr->registerParticleGhostMask(m_realm, 1);
  }
}

void
ItoSolver::setPhase(const phase::which_phase a_phase)
{
  CH_TIME("ItoSolver::setPhase");
  if (m_verbosity > 5) {
    pout() << m_name + "::setPhase" << endl;
  }

  m_phase = a_phase;
}

void
ItoSolver::setVerbosity(const int a_verbosity)
{
  CH_TIME("ItoSolver::setVerbosity");

  m_verbosity = a_verbosity;

  if (m_verbosity > 5) {
    pout() << m_name + "::setVerbosity" << endl;
  }
}

void
ItoSolver::setTime(const int a_step, const Real a_time, const Real a_dt)
{
  CH_TIME("ItoSolver::setTime");
  if (m_verbosity > 5) {
    pout() << m_name + "::setTime" << endl;
  }

  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;
}

void
ItoSolver::initialData()
{
  CH_TIME("ItoSolver::initialData");
  if (m_verbosity > 5) {
    pout() << m_name + "::initialData" << endl;
  }

  CH_assert(!m_species.isNull());

  // clang-format off
  // TLDR: This function will fetch the initial particles from the species and deposit them on the mesh. In most cases the various MPI ranks
  //       will have drawn a different set of initial particles (the only sane way to do it) and so those particles are put directly in
  //       the 'bulk' particle container. After that we remove the particles that fell inside the EB and deposit the particles on the mesh.
  // clang-format on

  ParticleContainer<ItoParticle>& bulkParticles = m_particleContainers.at(WhichContainer::Bulk);
  bulkParticles.clearParticles();

  // Copy the species' SoA seed particles into a buffer and ingest it. addParticlesDestructive() routes
  // each particle to its owning patch/level/rank via remap() and empties the buffer (so the species'
  // m_initialParticles is left intact for restart).
  ParticleSoA<ItoParticle> buffer;
  buffer.append(m_species->getInitialParticles());

  bulkParticles.addParticlesDestructive(buffer);

  // Generate particles from the initial density distribution.
  auto initialDensity = [&](const RealVect& x) -> Real {
    const auto& initialDensityFunc = m_species->getInitialDensity();

    return initialDensityFunc(x, m_time);
  };

  // Seed at the merge target, not at ppc_restart. This is a cold start, not a restart: ppc_restart
  // exists for redrawing particles from a 'checkpointing = numbers' checkpoint (see
  // readCheckpointLevelFluid()), and using it here meant seeding at one count and having the first
  // merge immediately take the population down to another. A cold start and the first merge now agree
  // on the target by construction.
  this->generateParticlesFromDensity(bulkParticles, initialDensity, m_particlesPerCell);

  constexpr Real tolerance = 0.0;

  // Add particles, remove the ones that are inside the EB, and then deposit
  this->removeCoveredParticles(bulkParticles, EBRepresentation::ImplicitFunction, tolerance);
  this->depositWeight(m_phi, bulkParticles, m_deposition, m_coarseFineDeposition);

  this->coarsenAndFillGhosts(m_phi);
}

void
ItoSolver::generateParticlesFromDensity(ParticleContainer<ItoParticle>&              a_particles,
                                        const std::function<Real(const RealVect x)>& a_densityFunc,
                                        const int a_maxParticlesPerCell) const noexcept
{
  CH_TIME("ItoSolver::generateParticlesFromDensity(int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::generateParticlesFromDensity(int)" << endl;
  }

  this->generateParticlesFromDensity(a_particles, a_densityFunc, Vector<int>(1, a_maxParticlesPerCell));
}

void
ItoSolver::generateParticlesFromDensity(ParticleContainer<ItoParticle>&              a_particles,
                                        const std::function<Real(const RealVect x)>& a_densityFunc,
                                        const Vector<int>& a_maxParticlesPerCell) const noexcept
{
  CH_TIME("ItoSolver::generateParticlesFromDensity(Vector<int>)");
  if (m_verbosity > 5) {
    pout() << m_name + "::generateParticlesFromDensity(Vector<int>)" << endl;
  }

  if (a_maxParticlesPerCell.size() < 1) {
    MayDay::Error("ItoSolver::generateParticlesFromDensity - empty particles-per-cell vector");
  }

  // Lambda which stochastically determines the number of particles in a cell. This is the mean number of particles,
  // plus a stochastic evaluation of whether or not to include the "fractional" particle.
  auto sampleParticles = [&](const Real a_volume, const Real a_density, const int a_ppc) -> std::vector<long long> {
    const Real meanNumParticles   = a_volume * a_density;
    const Real remainingParticles = meanNumParticles - std::floor(meanNumParticles);

    long long numParticles = llround(std::floor(meanNumParticles));
    if (Random::getUniformReal01() < remainingParticles) {
      numParticles += 1LL;
    }

    return ParticleManagement::partitionParticleWeights(numParticles, static_cast<long long>(a_ppc));
  };

  // Grid loop.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    // Same clamp applyCellMerger() uses: entry lvl when the vector is long enough, else the last one.
    const int ppc = (lvl < a_maxParticlesPerCell.size()) ? a_maxParticlesPerCell[lvl] : a_maxParticlesPerCell.back();

    const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit    = dbl.dataIterator();
    const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const Real               dx     = m_amr->getDx()[lvl];
    const Real               vol    = std::pow(dx, SpaceDim);
    const RealVect           probLo = m_amr->getProbLo();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex&          din        = dit[mybox];
      const Box&                cellbox    = dbl[din];
      const EBISBox&            ebisbox    = ebisl[din];
      const BaseFab<bool>&      validCells = (*m_amr->getValidCells(m_realm)[lvl])[din];
      ParticleSoA<ItoParticle>& particles  = a_particles[lvl][din];

      auto regularKernel = [&](const IntVect& iv) -> void {
        if (validCells(iv, 0) && ebisbox.isRegular(iv)) {
          const RealVect cellPos = probLo + (RealVect(iv) + 0.5 * RealVect::Unit) * dx;
          const Real     phi     = a_densityFunc(cellPos);

          const std::vector<long long> particleWeights = sampleParticles(vol, phi, ppc);

          const RealVect lo = probLo + (RealVect(iv)) * dx;
          const RealVect hi = probLo + (RealVect(iv) + RealVect::Unit) * dx;

          // Partition the particle weights.
          for (const auto& w : particleWeights) {
            const RealVect x = Random::randomPosition(lo, hi);

            particles.append(x, static_cast<double>(w));
          }
        }
      };

      auto irregularKernel = [&](const VolIndex& vof) -> void {
        const IntVect iv = vof.gridIndex();
        if (validCells(iv, 0) && ebisbox.isIrregular(iv)) {
          const Real     kappa         = ebisbox.volFrac(vof);
          const RealVect normal        = ebisbox.normal(vof);
          const RealVect bndryCentroid = ebisbox.bndryCentroid(vof);
          const RealVect cellPos       = probLo + (RealVect(iv) + 0.5 * RealVect::Unit) * dx;
          const Real     phi           = a_densityFunc(cellPos);

          // Compute the minimum box that encloses this cell.
          RealVect lo = -0.5 * RealVect::Unit;
          RealVect hi = +0.5 * RealVect::Unit;

          DataOps::computeMinValidBox(lo, hi, normal, bndryCentroid);

          // Partition particle weights.
          const std::vector<long long> particleWeights = sampleParticles(kappa * vol, phi, ppc);

          // Sample the particles.
          for (const auto& w : particleWeights) {
            const RealVect x = Random::randomPosition(cellPos, lo, hi, bndryCentroid, normal, dx, kappa);

            particles.append(x, static_cast<double>(w));
          }
        }
      };

      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[din];

      // Not vectorizable: per-cell variable-length particle generation (std::function density callback,
      // std::vector allocation, RNG, List append). Multi-cut N/A: the regular kernel guards isRegular and
      // cut cells go to the irregular kernel, which samples positions inside the EB-clipped cell volume.
      BoxLoops::loop<D_DECL(1, 1, 1)>(cellbox, regularKernel);
      BoxLoops::loop(vofit, irregularKernel);
    }
  }
}

void
ItoSolver::computeLoads(Vector<long int>& a_loads, const DisjointBoxLayout& a_dbl, const int a_level)
{
  CH_TIME("ItoSolver::computeLoads");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeLoads" << endl;
  }

  CH_assert(a_dbl.isClosed());
  CH_assert(a_dbl.size() > 0);

  const ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);

  a_loads.resize(a_dbl.size(), 0L);

  const DataIterator& dit = a_dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    a_loads[din.intCode()] = static_cast<long>(particles[a_level][din].size());
  }

  ParallelOps::sum(a_loads);
}

void
ItoSolver::removeCoveredParticles(const EBRepresentation a_representation, const Real a_tol)
{
  CH_TIME("ItoSolver::removeCoveredParticles(EBRepresentation, tolerance)");
  if (m_verbosity > 5) {
    pout() << m_name + "::removeCoveredParticles(EBRepresentation, tolerance)" << endl;
  }

  this->removeCoveredParticles(WhichContainer::Bulk, a_representation, a_tol);
}

void
ItoSolver::removeCoveredParticles(const WhichContainer   a_container,
                                  const EBRepresentation a_representation,
                                  const Real             a_tol)
{
  CH_TIME("ItoSolver::removeCoveredParticles(WhichContainer, EBRepresentation, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::removeCoveredParticles(WhichContainer, EBRepresentation, Real)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(a_container);

  this->removeCoveredParticles(particles, a_representation, a_tol);
}

void
ItoSolver::removeCoveredParticles(ParticleContainer<ItoParticle>& a_particles,
                                  const EBRepresentation          a_representation,
                                  const Real                      a_tol) const
{
  CH_TIME("ItoSolver::removeCoveredParticles(particles, EBRepresentation)");
  if (m_verbosity > 5) {
    pout() << m_name + "::removeCoveredParticles(particles, EBRepresentation)" << endl;
  }

  switch (a_representation) {
  case EBRepresentation::ImplicitFunction: {
    m_amr->removeCoveredParticlesIF(a_particles, m_phase, a_tol);

    break;
  }
  case EBRepresentation::Discrete: {
    m_amr->removeCoveredParticlesDiscrete(a_particles, m_phase, a_tol);

    break;
  }
  case EBRepresentation::Voxel: {
    m_amr->removeCoveredParticlesVoxels(a_particles, m_phase);

    break;
  }
  default: {
    MayDay::Error("ItoSolver::removeCoveredParticles - unsupported EB representation requested");

    break;
  }
  }
}

void
ItoSolver::transferCoveredParticles(const EBRepresentation a_representation, const Real a_tol)
{
  CH_TIME("ItoSolver::transferCoveredParticles(EBRepresentation, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::transferCoveredParticles(EBRepresentation, Real)" << endl;
  }

  this->transferCoveredParticles(WhichContainer::Bulk, WhichContainer::Covered, a_representation, a_tol);
}

void
ItoSolver::transferCoveredParticles(const WhichContainer   a_containerFrom,
                                    const WhichContainer   a_containerTo,
                                    const EBRepresentation a_representation,
                                    const Real             a_tol)
{
  CH_TIME("ItoSolver::transferCoveredParticles(WhichContainer, WhichContainer, EBRepresentation, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::transferCoveredParticles(WhichContainer, WhichContainer, EBRepresentation, Real)" << endl;
  }

  ParticleContainer<ItoParticle>& particlesFrom = this->getParticles(a_containerFrom);
  ParticleContainer<ItoParticle>& particlesTo   = this->getParticles(a_containerTo);

  this->transferCoveredParticles(particlesFrom, particlesTo, a_representation, a_tol);
}

void
ItoSolver::transferCoveredParticles(ParticleContainer<ItoParticle>& a_particlesFrom,
                                    ParticleContainer<ItoParticle>& a_particlesTo,
                                    const EBRepresentation          a_representation,
                                    const Real                      a_tol) const
{
  CH_TIME("ItoSolver::transferCoveredParticles(ParticleContainer, ParticleContainer, EBRepresentation, Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::transferCoveredParticles(ParticleContainer, ParticleContainer, EBRepresentation, Real)"
           << endl;
  }

  switch (a_representation) {
  case EBRepresentation::ImplicitFunction: {
    m_amr->transferCoveredParticlesIF(a_particlesFrom, a_particlesTo, m_phase, a_tol);

    break;
  }
  case EBRepresentation::Discrete: {
    m_amr->transferCoveredParticlesDiscrete(a_particlesFrom, a_particlesTo, m_phase, a_tol);

    break;
  }
  case EBRepresentation::Voxel: {
    m_amr->transferCoveredParticlesVoxels(a_particlesFrom, a_particlesTo, m_phase);

    break;
  }
  default: {
    MayDay::Error("ItoSolver::transferCoveredParticles -- logic bust");

    break;
  }
  }
}

void
ItoSolver::intersectParticles(const EBIntersection                                               a_ebIntersection,
                              const bool                                                         a_deleteParticles,
                              const std::function<void(ParticleSoA<ItoParticle>&, std::size_t)>& a_nonDeletionModifier)
{
  CH_TIME("ItoSolver::intersectParticles(EBIntersection, bool)");
  if (m_verbosity > 5) {
    pout() << m_name + "::intersectParticles(EBIntersection, bool)" << endl;
  }

  this->intersectParticles(WhichContainer::Bulk,
                           WhichContainer::EB,
                           WhichContainer::Domain,
                           a_ebIntersection,
                           a_deleteParticles,
                           a_nonDeletionModifier);
}

void
ItoSolver::intersectParticles(const WhichContainer                                               a_particles,
                              const WhichContainer                                               a_ebParticles,
                              const WhichContainer                                               a_domainParticles,
                              const EBIntersection                                               a_ebIntersection,
                              const bool                                                         a_deleteParticles,
                              const std::function<void(ParticleSoA<ItoParticle>&, std::size_t)>& a_nonDeletionModifier)
{
  CH_TIME("ItoSolver::intersectParticles(WhichContainerx3, EBIntersection, bool)");
  if (m_verbosity > 5) {
    pout() << m_name + "::intersectParticles(WhichContainerx3, EBIntersection, bool)" << endl;
  }

  ParticleContainer<ItoParticle>& particles       = this->getParticles(a_particles);
  ParticleContainer<ItoParticle>& ebParticles     = this->getParticles(a_ebParticles);
  ParticleContainer<ItoParticle>& domainParticles = this->getParticles(a_domainParticles);

  this->intersectParticles(particles,
                           ebParticles,
                           domainParticles,
                           a_ebIntersection,
                           a_deleteParticles,
                           a_nonDeletionModifier);
}

void
ItoSolver::intersectParticles(ParticleContainer<ItoParticle>&                                    a_particles,
                              ParticleContainer<ItoParticle>&                                    a_ebParticles,
                              ParticleContainer<ItoParticle>&                                    a_domainParticles,
                              const EBIntersection                                               a_ebIntersection,
                              const bool                                                         a_deleteParticles,
                              const std::function<void(ParticleSoA<ItoParticle>&, std::size_t)>& a_nonDeletionModifier)
{
  CH_TIME("ItoSolver::intersectParticles(ParticleContainerx3, EBIntersection, bool)");
  if (m_verbosity > 5) {
    pout() << m_name + "::intersectParticles(ParticleContainerx3, EBIntersection, bool)" << endl;
  }

  CH_assert(!a_particles.isOrganizedByCell());
  CH_assert(!a_ebParticles.isOrganizedByCell());
  CH_assert(!a_domainParticles.isOrganizedByCell());

  constexpr Real tolerance = 0.0;

  switch (a_ebIntersection) {
  case EBIntersection::Raycast: {
    m_amr->intersectParticlesRaycastIF<D_DECL(&ItoParticle::old_x, &ItoParticle::old_y, &ItoParticle::old_z)>(
      a_particles,
      a_ebParticles,
      a_domainParticles,
      m_phase,
      tolerance,
      a_deleteParticles,
      a_nonDeletionModifier);

    break;
  }
  case EBIntersection::Bisection: {
    m_amr->intersectParticlesBisectIF<D_DECL(&ItoParticle::old_x, &ItoParticle::old_y, &ItoParticle::old_z)>(
      a_particles,
      a_ebParticles,
      a_domainParticles,
      m_phase,
      m_bisectionStep,
      a_deleteParticles,
      a_nonDeletionModifier);

    break;
  }
  default: {
    MayDay::Error("ItoSolver::intersectParticles - unsupported EB intersection requested");

    break;
  }
  }
}

void
ItoSolver::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("ItoSolver::regrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::regrid" << endl;
  }

  CH_assert(a_lmin >= 0);
  CH_assert(a_oldFinestLevel >= 0);
  CH_assert(a_newFinestLevel >= 0);

  const int ncomp = 1;

  // Mesh data -- always allocate it.
  m_amr->allocate(m_phi, m_realm, m_phase, ncomp);

  // For "redistributed" particle deposition
  m_amr->allocate(m_depositionNC, m_realm, m_phase, ncomp);
  m_amr->allocate(m_massDiff, m_realm, m_phase, ncomp);

  // Per-cell scratch for the kd merges (see ParticleManagement::mergeKDCarve). One component, one
  // ghost cell -- the shape those functions document.
  m_amr->allocate(m_kdMergeCellHistogram, m_realm, 1, 1);
  m_amr->allocate(m_kdMergeLeafQuota, m_realm, 1, 1);

  // Only allocate memory for velocity if we actually have a mobile solver
  if (m_isMobile) {
    m_amr->allocate(m_mobilityFunction, m_realm, m_phase, ncomp); //
    m_amr->allocate(m_velocityFunction, m_realm, m_phase, SpaceDim);
  }
  else {
    m_amr->allocatePointer(m_mobilityFunction, m_realm);
    m_amr->allocatePointer(m_velocityFunction, m_realm);
  }

  // Only allocate memory if we actually have a diffusion solver
  if (m_isDiffusive) {
    m_amr->allocate(m_diffusionFunction, m_realm, m_phase, ncomp);
  }
  else {
    m_amr->allocatePointer(m_diffusionFunction, m_realm);
  }

  // Regrid particle containers. The bulk one is empty -- preRegrid() moved its particles into
  // m_regridParticles -- but it still has to make the trip, since ParticleContainer::regrid() is what
  // rebuilds its holders over the new layout.
  for (auto& container : m_particleContainers) {
    ParticleContainer<ItoParticle>& particles = container.second;

    m_amr->remapToNewGrids(particles, a_lmin, a_newFinestLevel);
  }

  // Redistribute the reduced particles onto the new grids, merge them there, and only then rebuild
  // ItoParticles. This ordering is the point of the exercise: the de-refinement pile-up exists between
  // the redistribution and the merge, and here it is held at 53 B per particle rather than 149 B, and
  // merged away before anything is materialised at the larger size.
  //
  // Safe to run the AMR-collective merges from here: Driver::regrid() calls AmrMesh::regridOperators()
  // before the time stepper's regrid, so the particle ghost masks exist and m_amr->getFinestLevel() is
  // already the new finest level. The kd merges' mesh scratch (m_kdMergeCellHistogram/m_kdMergeLeafQuota)
  // was allocated above.
  {
    CH_TIME("ItoSolver::regrid::mergeRegridParticles");

    m_amr->remapToNewGrids(m_regridParticles, a_lmin, a_newFinestLevel);

    // nullopt => use the solver's own merge_algorithm.
    const ParticleManagement::ParticleMergeMethod method = m_regridMergeMethod.value_or(m_mergeMethod);

    this->mergeLite(m_regridParticles, m_particlesPerCell, method);
    this->rebuildFromMergeContainer(WhichContainer::Bulk, m_regridParticles);

    // Backstop for any leaf the rebuild skipped. rebuildFromMergeContainer() releases each leaf as it
    // consumes it, but this container is a member and would otherwise hold its arenas until the next
    // regrid -- unlike the function-local merge containers, which die at scope exit.
    m_regridParticles.clearParticles();

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[lvl];
      const DataIterator&      dit  = dbl.dataIterator();
      const int                nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        m_regridParticles[lvl][dit[mybox]].shrinkToFit();
      }
    }
  }
}

void
ItoSolver::setSpecies(const RefCountedPtr<ItoSpecies>& a_species)
{
  CH_TIME("ItoSolver::setSpecies");
  if (m_verbosity > 5) {
    pout() << m_name + "::setSpecies" << endl;
  }

  m_species     = a_species;
  m_name        = a_species->getName();
  m_isDiffusive = m_species->isDiffusive();
  m_isMobile    = m_species->isMobile();
}

void
ItoSolver::allocate()
{
  CH_TIME("ItoSolver::allocate");
  if (m_verbosity > 5) {
    pout() << m_name + "::allocate" << endl;
  }

  CH_assert(!m_species.isNull());

  const int ncomp = 1;

  // Mesh data -- always allocate it.
  m_amr->allocate(m_phi, m_realm, m_phase, ncomp);

  // For "redistributed" particle deposition
  m_amr->allocate(m_depositionNC, m_realm, m_phase, ncomp);
  m_amr->allocate(m_massDiff, m_realm, m_phase, ncomp);

  // Per-cell scratch for the kd merges (see ParticleManagement::mergeKDCarve). One component, one
  // ghost cell -- the shape those functions document.
  m_amr->allocate(m_kdMergeCellHistogram, m_realm, 1, 1);
  m_amr->allocate(m_kdMergeLeafQuota, m_realm, 1, 1);

  // Only allocate memory for velocity if we actually have a mobile solver
  if (m_isMobile) {
    m_amr->allocate(m_mobilityFunction, m_realm, m_phase, ncomp); //
    m_amr->allocate(m_velocityFunction, m_realm, m_phase, SpaceDim);
  }
  else {
    m_amr->allocatePointer(m_mobilityFunction, m_realm);
    m_amr->allocatePointer(m_velocityFunction, m_realm);
  }

  // Only allocate memory if we actually have a diffusion solver
  if (m_isDiffusive) {
    m_amr->allocate(m_diffusionFunction, m_realm, m_phase, ncomp);
  }
  else {
    m_amr->allocatePointer(m_diffusionFunction, m_realm);
  }

  // ParticleContainer is move-only-incapable, so default-construct each container in place (operator[])
  // rather than emplace()ing a temporary.
  const std::array<WhichContainer, 6> containers = {WhichContainer::Bulk,
                                                    WhichContainer::EB,
                                                    WhichContainer::Domain,
                                                    WhichContainer::Source,
                                                    WhichContainer::Covered,
                                                    WhichContainer::Scratch};

  for (const auto& which : containers) {
    m_amr->allocate(m_particleContainers[which], m_realm);
  }
}

#ifdef CH_USE_HDF5
void
ItoSolver::writeCheckpointLevel(HDF5Handle& a_handle, const int a_level) const
{
  CH_TIME("ItoSolver::writeCheckpointLevel");
  if (m_verbosity > 5) {
    pout() << m_name + "::writeCheckpointLevel" << endl;
  }

  // Write state.
  write(a_handle, *m_phi[a_level], m_name);

  // Write particles.
  switch (m_checkpointing) {
  case WhichCheckpoint::Particles: {
    this->writeCheckPointLevelParticles(a_handle, a_level);

    break;
  }
  case WhichCheckpoint::Numbers: {
    this->writeCheckPointLevelFluid(a_handle, a_level);

    break;
  }
  default: {
    MayDay::Error("ItoSolver::writeCheckpointLevel -- logic bust");

    break;
  }
  }
}
#endif

#ifdef CH_USE_HDF5
void
ItoSolver::writeCheckPointLevelParticles(HDF5Handle& a_handle, const int a_level) const
{
  CH_TIME("ItoSolver::writeCheckPointLevelParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::writeCheckPointLevelParticles" << endl;
  }

  // I call this _particlesP to distinguish it from the "fluid" checkpointing method.
  const std::string str = m_name + "_particlesP";

  const ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);

  DischargeIO::writeCheckParticlesToHDF(a_handle, particles[a_level], str);
}
#endif

#ifdef CH_USE_HDF5
void
ItoSolver::writeCheckPointLevelFluid(HDF5Handle& a_handle, const int a_level) const
{
  CH_TIME("ItoSolver::writeCheckPointLevelFluid");
  if (m_verbosity > 5) {
    pout() << m_name + "::writeCheckPointLevelFluid" << endl;
  }

  // clang-format off
  // TLDR: This routine checkpoints the particle data using the number of particles in a grid cell. When the simulation is restarted, we read the
  //       number of particles per cell from the HDF5 file and re-initialize the particles. However, this function does NOT currently store the energy
  //       (or other parameters of interest) in the HDF5 file. Only the number of particles is available. I don't expect this function to be widely used
  //       by anyone.
  // clang-format on

  // I call this _particlesF to distinguish it from the "particle" checkpointing method.
  const std::string str = m_name + "_particlesF";

  // Handles to relevant grid information
  const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[a_level];
  const DataIterator&      dit    = dbl.dataIterator();
  const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[a_level];
  const RealVect           dx     = m_amr->getDx()[a_level] * RealVect::Unit;
  const RealVect           probLo = m_amr->getProbLo();

  // Handle to the particles that will be checkpointed.
  const ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);

  // Create transient storage that holds the particle numbers.
  LevelData<EBCellFAB> particleNumbers(dbl, m_nComp, IntVect::Zero, EBCellFactory(ebisl));
  DataOps::setValue(particleNumbers, 0.0);

  // Now go through the grid and add the number of particles in each cell
  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const Box cellBox = dbl[din];

    // For multi-valued cells all the particles go onto the first vof.
    BaseFab<Real>& particleNumbersFAB = particleNumbers[din].getSingleValuedFAB();

    // Accumulate the per-cell sum of particle weights directly off the SoA leaf (the leaf is const here, so
    // we compute the owning cell per particle rather than cell-sorting in place).
    const ParticleSoA<ItoParticle>& leaf = particles[a_level][din];

    const std::size_t numParticles = leaf.size();
    for (std::size_t i = 0; i < numParticles; i++) {
      const RealVect pos = leaf.position(i);
      const IntVect  iv  = ParticleOps::getParticleCellIndex(pos, probLo, dx[0]);

      if (cellBox.contains(iv)) {
        particleNumbersFAB(iv, m_comp) += leaf.weight(i);
      }
    }
  }

  // Finally, write the particle numbers to HDF5.
  write(a_handle, particleNumbers, str);
}
#endif

#ifdef CH_USE_HDF5
void
ItoSolver::readCheckpointLevel(HDF5Handle& a_handle, const int a_level)
{
  CH_TIME("ItoSolver::readCheckpointLevel");
  if (m_verbosity > 5) {
    pout() << m_name + "::readCheckpointLevel" << endl;
  }

  // Read state vector
  read<EBCellFAB>(a_handle, *m_phi[a_level], m_name, m_amr->getGrids(m_realm)[a_level], Interval(0, 0), false);

  // Instantiate the particles
  switch (m_checkpointing) {
  case WhichCheckpoint::Particles: {
    this->readCheckpointLevelParticles(a_handle, a_level);

    break;
  }
  case WhichCheckpoint::Numbers: {
    this->readCheckpointLevelFluid(a_handle, a_level);

    break;
  }
  default: {
    MayDay::Error("ItoSolver::readCheckpointLevel -- logic bust");

    break;
  }
  }
}
#endif

#ifdef CH_USE_HDF5
void
ItoSolver::readCheckpointLevelParticles(HDF5Handle& a_handle, const int a_level)
{
  CH_TIME("ItoSolver::readCheckpointLevelParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::readCheckpointLevelParticles" << endl;
  }

  // clang-format off
  // TLDR: This function is the one that reads SimpleItoParticles from the checkpoint file and instantiates full ItoParticle's from that. Recalling
  //       writeCheckpointLevelParticles we only stored the weight, position, and energy of the particles. Here we read that information back in.
  // clang-format on

  // This is the particle container that we will fill.
  ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);

  CH_assert(m_checkpointing == WhichCheckpoint::Particles);
  CH_assert(!particles.isOrganizedByCell());

  const std::string str = m_name + "_particlesP";

  DischargeIO::readCheckParticlesFromHDF(a_handle, particles[a_level], str);
}
#endif

#ifdef CH_USE_HDF5
void
ItoSolver::readCheckpointLevelFluid(HDF5Handle& a_handle, const int a_level)
{
  CH_TIME("ItoSolver::readCheckpointLevelFluid");
  if (m_verbosity > 5) {
    pout() << m_name + "::readCheckpointLevelFluid" << endl;
  }

  CH_assert(m_checkpointing == WhichCheckpoint::Numbers);

  const std::string str = m_name + "_particlesF";

  constexpr int  comp     = 0;
  constexpr int  numComp  = 1;
  constexpr int  numGhost = 0;
  const Interval interv   = Interval(comp, comp);

  // Allocate some storage that we can read into.
  EBAMRCellData particlesPerCell;
  m_amr->allocate(particlesPerCell, m_realm, m_phase, numComp, numGhost);

  read<EBCellFAB>(a_handle, *particlesPerCell[a_level], str, m_amr->getGrids(m_realm)[a_level], interv, false);

  // particlesPerCell holds the number of particles per cell -- call the other version which instantiates new particles
  // from that.
  this->drawNewParticles(*particlesPerCell[a_level], a_level, m_restartPPC);
}
#endif

void
ItoSolver::drawNewParticles(const LevelData<EBCellFAB>& a_particlesPerCell, const int a_level, const int a_newPPC)
{
  CH_TIME("ItoSolver::drawNewParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::drawNewParticles" << endl;
  }

  // Handle to grid information.
  const RealVect           probLo = m_amr->getProbLo();
  const Real               dx     = m_amr->getDx()[a_level];
  const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[a_level];
  const DataIterator&      dit    = dbl.dataIterator();
  const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[a_level];

  // Particle container that we will fill.
  ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);

  // Go through each patch and instantiate new particles.
  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const Box&           cellBox = dbl[din];
    const EBISBox&       ebisbox = ebisl[din];
    const BaseFab<Real>& ppc     = a_particlesPerCell[din].getSingleValuedFAB();

    // This should draw new particles rather than append -- so clear out any old particles.
    ParticleSoA<ItoParticle>& myParticles = particles[a_level][din];
    myParticles.clear();

    // Kernel region for cut-cells.
    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_level])[din];

    // Regular kernel
    auto regularKernel = [&](const IntVect& iv) -> void {
      // Do regular cells -- in these cells we only need to draw a random position somewhere inside the cubic cell.
      // Easy.
      if (ebisbox.isRegular(iv)) {

        // This bit of code will take the number of physical particles and divide them into a_newPPC particles with
        // approximately equal weights. It is possible that one of the particles will have a larger particle weight than
        // the others.
        const std::vector<long long> weights = ParticleManagement::partitionParticleWeights(llround(ppc(iv)),
                                                                                            (long long)a_newPPC);

        // Settings for drawing new particles in the current cell.
        const RealVect minLo = -0.5 * RealVect::Unit;
        const RealVect minHi = 0.5 * RealVect::Unit;
        const RealVect norma = RealVect::Zero;
        const RealVect centr = RealVect::Zero;
        const RealVect pos   = probLo + (RealVect(iv) + 0.5 * RealVect::Unit) * dx;
        const Real     kappa = 1.0;

        for (const auto& w : weights) {
          const RealVect x = Random::randomPosition(pos, minLo, minHi, centr, norma, dx, kappa);

          myParticles.append(x, 1.0 * static_cast<double>(w));
        }
      }
    };

    // Irregular kernel. Do the same for irregular cells. This differs from the regular-cell case only in that the
    // positions are checked against the EB.
    auto irregularKernel = [&](const VolIndex& vof) -> void {
      const IntVect  iv    = vof.gridIndex();
      const RealVect cent  = ebisbox.bndryCentroid(vof);
      const RealVect norm  = ebisbox.normal(vof);
      const RealVect pos   = probLo + dx * (RealVect(iv) + 0.5 * RealVect::Unit);
      const Real     kappa = ebisbox.volFrac(vof);

      const auto numPhysicalParticles = (unsigned long long)llround(ppc(iv));

      if (numPhysicalParticles > 0ULL) {

        // No multi-valued cells please -- I don't know how to handle them.
        CH_assert(!ebisbox.isMultiValued(iv));

        // Compute a small box that encloses the cut-cell volume
        RealVect minLo = -0.5 * RealVect::Unit;
        RealVect minHi = 0.5 * RealVect::Unit;
        if (kappa < 1.0) {
          DataOps::computeMinValidBox(minLo, minHi, norm, cent);
        }

        // This bit of code will take the number of physical particles and divide them into a_newPPC particles with
        // approximately equal weights. It is possible that one of the particles will have a larger particle weight than
        // the others.
        const std::vector<long long> weights = ParticleManagement::partitionParticleWeights(llround(ppc(iv)),
                                                                                            (long long)a_newPPC);

        for (const auto& w : weights) {
          const RealVect x = Random::randomPosition(pos, minLo, minHi, cent, norm, dx, kappa);

          myParticles.append(x, 1.0 * static_cast<double>(w));
        }
      }
    };

    // Run the kernels.
    // Not vectorizable: per-cell variable-length particle redraw (partitionParticleWeights allocation,
    // RNG, List append). Multi-cut N/A: regular kernel guards isRegular; cut cells go to the irregular
    // kernel, which clips sampling to the EB cell volume (computeMinValidBox).
    BoxLoops::loop<D_DECL(1, 1, 1)>(cellBox, regularKernel);
    BoxLoops::loop(vofit, irregularKernel);
  }
}

int
ItoSolver::getNumberOfPlotVariables() const
{
  CH_TIME("ItoSolver::getNumberOfPlotVariables");
  if (m_verbosity > 5) {
    pout() << m_name + "::getNumberOfPlotVariables" << endl;
  }

  int numPlotVars = 0;

  if (m_plotPhi) {
    numPlotVars += 1;
  }
  if (m_plotDiffCo && m_isDiffusive) {
    numPlotVars += 1;
  }
  if (m_plotVelocity && m_isMobile) {
    numPlotVars += SpaceDim;
  }
  if (m_plotParticles) {
    numPlotVars += 1;
  }
  if (m_plotParticlesEB) {
    numPlotVars += 1;
  }
  if (m_plotParticlesDomain) {
    numPlotVars += 1;
  }
  if (m_plotParticlesSource) {
    numPlotVars += 1;
  }
  if (m_plotParticlesCovered) {
    numPlotVars += 1;
  }
  if (m_plotEnergyDensity) {
    numPlotVars += 1;
  }
  if (m_plotAverageEnergy) {
    numPlotVars += 1;
  }

  return numPlotVars;
}

Vector<std::string>
ItoSolver::getPlotVariableNames() const
{
  CH_TIME("ItoSolver::getPlotVariableNames");
  if (m_verbosity > 5) {
    pout() << m_name + "::getPlotVariableNames" << endl;
  }

  Vector<std::string> names(0);

  if (m_plotPhi) {
    names.push_back(m_name + " phi");
  }
  if (m_plotDiffCo && m_isDiffusive) {
    names.push_back(m_name + " diffusion_coefficient");
  }
  if (m_plotVelocity && m_isMobile) {
    names.push_back("x-Velocity field " + m_name);
    names.push_back("y-Velocity field " + m_name);
    if (SpaceDim == 3) {
      names.push_back("z-Velocity field " + m_name);
    }
  }
  if (m_plotParticles) {
    names.push_back(m_name + " particles");
  }
  if (m_plotParticlesEB) {
    names.push_back(m_name + " eb_particles");
  }
  if (m_plotParticlesDomain) {
    names.push_back(m_name + " domain_particles");
  }
  if (m_plotParticlesSource) {
    names.push_back(m_name + " source_particles");
  }
  if (m_plotParticlesCovered) {
    names.push_back(m_name + " covered_particles");
  }
  if (m_plotEnergyDensity) {
    names.push_back(m_name + " energy * phi");
  }
  if (m_plotAverageEnergy) {
    names.push_back(m_name + " average_energy");
  }

  return names;
}

void
ItoSolver::writePlotData(LevelData<EBCellFAB>& a_output,
                         int&                  a_comp,
                         const std::string&    a_outputRealm,
                         const int             a_level) const noexcept
{
  CH_TIMERS("ItoSolver::writePlotData");
  CH_TIMER("ItoSolver::writePlotData::mesh_data", t1);
  CH_TIMER("ItoSolver::writePlotData::particle_data", t2);
  if (m_verbosity > 5) {
    pout() << m_name + "::writePlotData" << endl;
  }

  CH_assert(a_comp >= 0);
  CH_assert(a_level >= 0);
  CH_assert(a_level <= m_amr->getFinestLevel());

  LevelData<EBCellFAB> scratch;
  m_amr->allocate(scratch, m_realm, m_phase, a_level, 1);

  // Write phi
  CH_START(t1);
  if (m_plotPhi) {
    const Interval srcInterval(m_comp, m_comp);
    const Interval dstInterval(a_comp, a_comp);

    this->writeData(a_output, a_comp, m_phi, a_outputRealm, a_level, false, true);
  }

  // Plot diffusion coefficient
  if (m_plotDiffCo && m_isDiffusive) {
    const Interval srcInterval(m_comp, m_comp);
    const Interval dstInterval(a_comp, a_comp);

    this->writeData(a_output, a_comp, m_diffusionFunction, a_outputRealm, a_level, false, true);
  }

  // Write velocities
  if (m_plotVelocity && m_isMobile) {
    const Interval srcInterval(m_comp, SpaceDim - 1);
    const Interval dstInterval(a_comp, a_comp + SpaceDim - 1);

    this->writeData(a_output, a_comp, m_velocityFunction, a_outputRealm, a_level, false, true);
  }
  CH_STOP(t1);

  CH_START(t2);
  if (m_plotParticles) {
    this->depositWeightNGP(scratch, m_particleContainers.at(WhichContainer::Bulk), a_level);

    m_amr->copyData(a_output, scratch, a_level, a_outputRealm, m_realm, Interval(a_comp, a_comp), Interval(0, 0));

    a_comp++;
  }
  if (m_plotParticlesEB) {
    this->depositWeightNGP(scratch, m_particleContainers.at(WhichContainer::EB), a_level);

    m_amr->copyData(a_output, scratch, a_level, a_outputRealm, m_realm, Interval(a_comp, a_comp), Interval(0, 0));

    a_comp++;
  }
  if (m_plotParticlesDomain) {
    this->depositWeightNGP(scratch, m_particleContainers.at(WhichContainer::Domain), a_level);

    m_amr->copyData(a_output, scratch, a_level, a_outputRealm, m_realm, Interval(a_comp, a_comp), Interval(0, 0));

    a_comp++;
  }
  if (m_plotParticlesSource) {
    this->depositWeightNGP(scratch, m_particleContainers.at(WhichContainer::Source), a_level);

    m_amr->copyData(a_output, scratch, a_level, a_outputRealm, m_realm, Interval(a_comp, a_comp), Interval(0, 0));

    a_comp++;
  }
  if (m_plotParticlesCovered) {
    this->depositWeightNGP(scratch, m_particleContainers.at(WhichContainer::Covered), a_level);

    m_amr->copyData(a_output, scratch, a_level, a_outputRealm, m_realm, Interval(a_comp, a_comp), Interval(0, 0));

    a_comp++;
  }
  if (m_plotEnergyDensity) {
    this->depositGatheredNGP(scratch,
                             m_particleContainers.at(WhichContainer::Bulk),
                             a_level,
                             [](const ParticleSoA<ItoParticle>& a_leaf, const std::size_t a_i) -> Real {
                               return a_leaf.weight(a_i) * a_leaf.template get<&ItoParticle::energy>(a_i);
                             });

    m_amr->copyData(a_output, scratch, a_level, a_outputRealm, m_realm, Interval(a_comp, a_comp), Interval(0, 0));

    a_comp++;
  }
  if (m_plotAverageEnergy) {
    LevelData<EBCellFAB> weight;
    m_amr->allocate(weight, m_realm, m_phase, a_level, 1);

    this->depositWeightNGP(weight, m_particleContainers.at(WhichContainer::Bulk), a_level);

    this->depositGatheredNGP(scratch,
                             m_particleContainers.at(WhichContainer::Bulk),
                             a_level,
                             [](const ParticleSoA<ItoParticle>& a_leaf, const std::size_t a_i) -> Real {
                               return a_leaf.weight(a_i) * a_leaf.template get<&ItoParticle::energy>(a_i);
                             });

    // Set scratch = totalEnergy/totalWeight
    DataOps::divideFallback(scratch, weight, 0.0, *m_amr->getMultiCutVofIterator(m_realm, m_phase)[a_level]);

    m_amr->copyData(a_output, scratch, a_level, a_outputRealm, m_realm, Interval(a_comp, a_comp), Interval(0, 0));

    a_comp++;
  }
  CH_STOP(t2);
}

void
ItoSolver::writeData(LevelData<EBCellFAB>& a_output,
                     int&                  a_comp,
                     const EBAMRCellData&  a_data,
                     const std::string&    a_outputRealm,
                     const int             a_level,
                     const bool            a_interpToCentroids,
                     const bool            a_interpGhost) const noexcept
{
  CH_TIMERS("ItoSolver::writeData");
  CH_TIMER("ItoSolver::writeData::allocate", t1);
  CH_TIMER("ItoSolver::writeData::local_copy", t2);
  CH_TIMER("ItoSolver::writeData::interp_ghost", t3);
  CH_TIMER("ItoSolver::writeData::interp_centroid", t4);
  CH_TIMER("ItoSolver::writeData::final_copy", t5);
  if (m_verbosity > 5) {
    pout() << m_name + "::writeData" << endl;
  }

  // Number of components we are working with.
  const int numComp = a_data[a_level]->nComp();

  CH_START(t1);
  LevelData<EBCellFAB> scratch;
  m_amr->allocate(scratch, m_realm, m_phase, a_level, numComp);
  CH_STOP(t1);

  CH_START(t2);
  m_amr->copyData(scratch, *a_data[a_level], a_level, m_realm, m_realm);
  scratch.exchange();
  CH_START(t2);

  // Interpolate ghost cells
  CH_START(t3);
  if (a_level > 0 && a_interpGhost) {
    m_amr->interpGhost(scratch, *a_data[a_level - 1], a_level, m_realm, m_phase);
  }
  CH_STOP(t3);

  CH_START(t4);
  if (a_interpToCentroids) {
    m_amr->interpToCentroids(scratch, m_realm, m_phase, a_level);
  }
  CH_STOP(t4);

  DataOps::setCoveredValue(scratch, *m_amr->getCoveredCells(m_realm, m_phase)[a_level], 0.0);

  CH_START(t5);
  const Interval srcInterv(0, numComp - 1);
  const Interval dstInterv(a_comp, a_comp + numComp - 1);

  m_amr->copyData(a_output,
                  scratch,
                  a_level,
                  a_outputRealm,
                  m_realm,
                  dstInterv,
                  srcInterv,
                  CopyStrategy::ValidGhost,
                  CopyStrategy::ValidGhost);

  CH_STOP(t5);

  a_comp += numComp;
}

void
ItoSolver::depositConductivity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const
{
  CH_TIME("ItoSolver::depositConductivity(EBAMRCellData, ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositConductivity(EBAMRCellData, ParticleContainer)" << endl;
  }

  this->depositConductivity(a_phi, a_particles, m_deposition, m_coarseFineDeposition);
}

void
ItoSolver::depositConductivity(EBAMRCellData&                  a_phi,
                               ParticleContainer<ItoParticle>& a_particles,
                               const DepositionType            a_deposition,
                               const CoarseFineDeposition      a_coarseFineDeposition) const
{
  CH_TIME("ItoSolver::depositConductivity(EBAMRCellData, ParticleContainer, DepositionType)");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositConductivity(EBAMRCellData, ParticleContainer, DepositionType)" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);
  CH_assert(!a_particles.isOrganizedByCell());

  if (m_isMobile) {
    this->depositGathered(a_phi,
                          a_particles,
                          a_deposition,
                          a_coarseFineDeposition,
                          [](const ParticleSoA<ItoParticle>& a_leaf, const std::size_t a_i) -> Real {
                            return a_leaf.weight(a_i) * a_leaf.template get<&ItoParticle::mobility>(a_i);
                          });

    this->coarsenAndFillGhosts(a_phi);
  }
  else {
    DataOps::setValue(a_phi, 0.0);
  }
}

void
ItoSolver::depositDiffusivity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const
{
  CH_TIME("ItoSolver::depositDiffusivity(EBAMRCellData, ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositDiffusivity(EBAMRCellData, ParticleContainer)" << endl;
  }

  this->depositDiffusivity(a_phi, a_particles, m_deposition, m_coarseFineDeposition);
}

void
ItoSolver::depositDiffusivity(EBAMRCellData&                  a_phi,
                              ParticleContainer<ItoParticle>& a_particles,
                              const DepositionType            a_deposition,
                              const CoarseFineDeposition      a_coarseFineDeposition) const
{
  CH_TIME("ItoSolver::depositDiffusivity(EBAMRCellData, ParticleContainer, DepositionType)");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositDiffusivity(EBAMRCellData, ParticleContainer, DepositionType)" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);
  CH_assert(!a_particles.isOrganizedByCell());

  this->depositGathered(a_phi,
                        a_particles,
                        a_deposition,
                        a_coarseFineDeposition,
                        [](const ParticleSoA<ItoParticle>& a_leaf, const std::size_t a_i) -> Real {
                          return a_leaf.weight(a_i) * a_leaf.template get<&ItoParticle::diffusion>(a_i);
                        });

  this->coarsenAndFillGhosts(a_phi);
}

void
ItoSolver::depositEnergyDensity(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const
{
  CH_TIME("ItoSolver::depositEnergyDensity(EBAMRCellData, ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositEnergyDensity(EBAMRCellData, ParticleContainer)" << endl;
  }

  this->depositEnergyDensity(a_phi, a_particles, m_deposition, m_coarseFineDeposition);
}

void
ItoSolver::depositEnergyDensity(EBAMRCellData&                  a_phi,
                                ParticleContainer<ItoParticle>& a_particles,
                                const DepositionType            a_deposition,
                                const CoarseFineDeposition      a_coarseFineDeposition) const
{
  CH_TIME("ItoSolver::depositEnergyDensity(EBAMRCellData, ParticleContainer, DepositionType)");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositEnergyDensity(EBAMRCellData, ParticleContainer, DepositionType)" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);
  CH_assert(!a_particles.isOrganizedByCell());

  this->depositGathered(a_phi,
                        a_particles,
                        a_deposition,
                        a_coarseFineDeposition,
                        [](const ParticleSoA<ItoParticle>& a_leaf, const std::size_t a_i) -> Real {
                          return a_leaf.weight(a_i) * a_leaf.template get<&ItoParticle::energy>(a_i);
                        });

  this->coarsenAndFillGhosts(a_phi);
}

void
ItoSolver::computeAverageMobility(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const
{
  CH_TIME("ItoSolver::computeAverageMobility(EBAMRCellData, ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeAverageMobility(EBAMRCellData, ParticleContainer)" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);
  CH_assert(!a_particles.isOrganizedByCell());

  EBAMRCellData weight;
  m_amr->allocate(weight, m_realm, m_phase, m_nComp);

  // Deposit weight*mu and weight
  this->depositGathered(a_phi,
                        a_particles,
                        m_deposition,
                        m_coarseFineDeposition,
                        [](const ParticleSoA<ItoParticle>& a_leaf, const std::size_t a_i) -> Real {
                          return a_leaf.weight(a_i) * a_leaf.template get<&ItoParticle::mobility>(a_i);
                        });

  this->coarsenAndFillGhosts(a_phi);

  this->depositWeight(weight, a_particles, m_deposition, m_coarseFineDeposition);

  this->coarsenAndFillGhosts(weight);

  // Make averageMobility = weight*mu/weight. If there is no weight then set the value to zero.
  constexpr Real zero = 0.0;

  DataOps::divideFallback(a_phi, weight, zero, m_amr->getMultiCutVofIterator(m_realm, m_phase));
}

void
ItoSolver::computeAverageDiffusion(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const
{
  CH_TIME("ItoSolver::computeAverageDiffusion(EBAMRCellData, ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeAverageDiffusion(EBAMRCellData, ParticleContainer)" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);
  CH_assert(!a_particles.isOrganizedByCell());

  EBAMRCellData weight;
  m_amr->allocate(weight, m_realm, m_phase, m_nComp);

  // Deposit weight*D and weight
  this->depositGathered(a_phi,
                        a_particles,
                        m_deposition,
                        m_coarseFineDeposition,
                        [](const ParticleSoA<ItoParticle>& a_leaf, const std::size_t a_i) -> Real {
                          return a_leaf.weight(a_i) * a_leaf.template get<&ItoParticle::diffusion>(a_i);
                        });

  this->coarsenAndFillGhosts(a_phi);

  this->depositWeight(weight, a_particles, m_deposition, m_coarseFineDeposition);

  this->coarsenAndFillGhosts(weight);

  // Make average diffusion coefficient = weight*D/weight. If there is no weight then set the value to zero.
  constexpr Real zero = 0.0;

  DataOps::divideFallback(a_phi, weight, zero, m_amr->getMultiCutVofIterator(m_realm, m_phase));
}

void
ItoSolver::computeAverageEnergy(EBAMRCellData& a_phi, ParticleContainer<ItoParticle>& a_particles) const
{
  CH_TIME("ItoSolver::computeAverageEnergy(EBAMRCellData, ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeAverageEnergy(EBAMRCellData, ParticleContainer)" << endl;
  }

  CH_assert(a_phi[0]->nComp() == 1);
  CH_assert(!a_particles.isOrganizedByCell());

  // Need scratch storage to deposit into
  EBAMRCellData weight;
  m_amr->allocate(weight, m_realm, m_phase, m_nComp);

  // Deposit weight*energy and weight
  this->depositGathered(a_phi,
                        a_particles,
                        m_deposition,
                        m_coarseFineDeposition,
                        [](const ParticleSoA<ItoParticle>& a_leaf, const std::size_t a_i) -> Real {
                          return a_leaf.weight(a_i) * a_leaf.template get<&ItoParticle::energy>(a_i);
                        });

  this->coarsenAndFillGhosts(a_phi);

  this->depositWeight(weight, a_particles, m_deposition, m_coarseFineDeposition);

  this->coarsenAndFillGhosts(weight);

  // Make average energy = weight*energy/weight. If there is no weight then set the value to zero.
  constexpr Real zero = 0.0;

  DataOps::divideFallback(a_phi, weight, zero, m_amr->getMultiCutVofIterator(m_realm, m_phase));
}

void
ItoSolver::depositParticles()
{
  CH_TIME("ItoSolver::depositParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositParticles" << endl;
  }

  this->depositParticles(WhichContainer::Bulk);
}

void
ItoSolver::depositParticles(const WhichContainer a_container)
{
  CH_TIME("ItoSolver::depositParticles(container)");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositParticles(container)" << endl;
  }

  this->depositWeight(m_phi, m_particleContainers.at(a_container), m_deposition, m_coarseFineDeposition);

  this->coarsenAndFillGhosts(m_phi);
}

void
ItoSolver::depositWeightNGP(LevelData<EBCellFAB>&                 a_output,
                            const ParticleContainer<ItoParticle>& a_particles,
                            const int                             a_level) const noexcept
{
  CH_TIME("ItoSolver::depositWeightNGP");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositWeightNGP" << endl;
  }

  CH_assert(a_level >= 0);
  CH_assert(a_level <= m_amr->getFinestLevel());

  const ProblemDomain&     domain = m_amr->getDomains()[a_level];
  const DisjointBoxLayout& dbl    = m_amr->getGrids(a_particles.getRealm())[a_level];
  const DataIterator&      dit    = dbl.dataIterator();
  const EBISLayout&        ebisl  = m_amr->getEBISLayout(a_particles.getRealm(), m_phase)[a_level];
  const Real               dx     = m_amr->getDx()[a_level];
  const RealVect           probLo = m_amr->getProbLo();

  CH_assert(a_output.disjointBoxLayout() == dbl);

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const Box      cellBox = dbl[din];
    const EBISBox& ebisbox = ebisl[din];

    EBParticleMesh particleMesh(domain, cellBox, ebisbox, dx * RealVect::Unit, probLo);

    EBCellFAB&                      output = a_output[din];
    const ParticleSoA<ItoParticle>& leaf   = a_particles[a_level][din];

    // The per-patch deposit INCREMENTS, so start from a clean slate.
    output.setVal(0.0);

    particleMesh.depositWeight(output, leaf, DepositionType::NGP, 1.0, true);
  }
}

void
ItoSolver::coarsenAndFillGhosts(EBAMRCellData& a_phi) const
{
  CH_TIME("ItoSolver::coarsenAndFillGhosts");
  if (m_verbosity > 5) {
    pout() << m_name + "::coarsenAndFillGhosts" << endl;
  }

  m_amr->conservativeAverage(a_phi, m_realm, m_phase);
  m_amr->interpGhost(a_phi, m_realm, m_phase);
}

void
ItoSolver::redistributeAMR(EBAMRCellData& a_phi) const
{
  CH_TIME("ItoSolver::redistributeAMR");
  if (m_verbosity > 5) {
    pout() << m_name + "::redistributeAMR" << endl;
  }

  // clang-format off
  // TLDR: When we entered this routine we had a_phi = m_i/dV but we actually want to have phi = m_i/(kappa*dV) so as to have
  //       meaningful densities. Thus, we can either run with a_phi just as it is, in which case it must be interpreted as an extended
  //       state into the EB. That is perfectly fine. But we can also use O(1) accurate redistribution in order to make the scheme
  //       completely conservative, if that is important.
  //
  //       If we use redistribution then we compute a hybrid update phiH = kappa*phi = a_phi in each cell. But we are then "missing"
  //       a mass kappa*phi - kappa*phiH = a_phi(1 - kappa). This mass can be smooshed into the neighboring grid cells. The code
  //       below does even more than that -- it can compute an update phiH = kappa*phi + (1-kappa)*phiNC where phiNC is a non-conservative
  //       type of update. In this case the mass loss is just like for fluid models: dM = kappa*(1-kappa)(phiC - phiNC). But this update
  //       is not strictly non-negative.
  // clang-format on

  if (m_useRedistribution) {
    this->depositNonConservative(m_depositionNC, a_phi);    // Compute m_depositionNC = sum(kappa*Wc)/sum(kappa)
    this->depositHybrid(a_phi, m_massDiff, m_depositionNC); // Compute hybrid deposition, including mass difference

    Vector<RefCountedPtr<EBFluxRedistribution>>& redistOps = m_amr->getRedistributionOp(m_realm, m_phase);
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const Real     scale     = 1.0;
      const Interval variables = Interval(0, 0);
      const bool     hasCoar   = lvl > 0;
      const bool     hasFine   = lvl < m_amr->getFinestLevel();

      if (hasCoar) {
        redistOps[lvl]->redistributeCoar(*a_phi[lvl - 1], *m_massDiff[lvl], scale, variables);
      }

      redistOps[lvl]->redistributeLevel(*a_phi[lvl], *m_massDiff[lvl], scale, variables);

      if (hasFine) {
        redistOps[lvl]->redistributeFine(*a_phi[lvl + 1], *m_massDiff[lvl], scale, variables);
      }
    }
  }
}

void
ItoSolver::depositNonConservative(EBAMRIVData& a_depositionNC, const EBAMRCellData& a_depositionKappaC) const
{
  CH_TIME("ItoSolver::depositNonConservative");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositNonConservative" << endl;
  }

  if (m_blendConservation) {
    m_amr->nonConservativeDivergence(a_depositionNC, a_depositionKappaC, m_realm, m_phase);
  }
  else {
    DataOps::setValue(a_depositionNC, 0.0);
  }
}

void
ItoSolver::depositHybrid(EBAMRCellData&     a_depositionH,
                         EBAMRIVData&       a_massDifference,
                         const EBAMRIVData& a_depositionNC) const
{
  CH_TIME("ItoSolver::depositHybrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositHybrid" << endl;
  }

  // TLDR: Compute divH = kappa*divC + (1-kappa)*divNC on each cell. Also compute mass difference.

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl   = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit   = dbl.dataIterator();
    const EBISLayout&        ebisl = m_amr->getEBISLayout(m_realm, m_phase)[lvl];

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      // On input, divH contains kappa*depositionWeights
      EBCellFAB&             divH    = (*a_depositionH[lvl])[din];
      BaseIVFAB<Real>&       deltaM  = (*a_massDifference[lvl])[din];
      const BaseIVFAB<Real>& divNC   = (*a_depositionNC[lvl])[din];
      const EBISBox&         ebisbox = ebisl[din];

      // Iteration space.
      VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[lvl])[din];

      auto kernel = [&](const VolIndex& vof) -> void {
        const Real kappa = ebisbox.volFrac(vof);
        const Real dc    = divH(vof, m_comp);
        const Real dnc   = divNC(vof, m_comp);

        // Note that if dc - kappa*dnc can be negative, i.e. we may end up STEALING mass
        // from other cells. This is why there is a flag m_blendConservation which always
        // gives positive definite results.
        divH(vof, m_comp)   = dc + (1.0 - kappa) * dnc;         // On output, contains hybrid divergence
        deltaM(vof, m_comp) = (1 - kappa) * (dc - kappa * dnc); // Remember, dc already scaled by kappa.
      };

      BoxLoops::loop(vofit, kernel);
    }
  }
}

bool
ItoSolver::isMobile() const
{
  CH_TIME("ItoSolver::isMobile");

  return m_isMobile;
}

bool
ItoSolver::isDiffusive() const
{
  CH_TIME("ItoSolver::isDiffusive");

  return m_isDiffusive;
}

void
ItoSolver::preRegrid(const int a_lbase, const int /*a_oldFinestLevel*/)
{
  CH_TIME("ItoSolver::preRegrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::preRegrid" << endl;
  }

  CH_assert(a_lbase >= 0);

  // The bulk particles cross the regrid in reduced form. At a de-refining regrid 2^D fine cells
  // collapse into one coarse cell, so between the redistribution and the merge each rank transiently
  // holds up to 2^D times the target particles per cell -- the largest single term in the particle
  // footprint. Extract now, so that everything from here to the merge in regrid() moves 53 B particles
  // instead of 149 B ones, and so that the merge can run before the pile-up is ever rebuilt at the
  // larger size. extractIntoMergeContainer() clears the ItoParticle arenas patch by patch as it goes.
  this->extractIntoMergeContainer(WhichContainer::Bulk, m_regridParticles);

  // BOTH containers make the round trip. The bulk ItoParticle container is empty now, but
  // ParticleContainer::regrid() is the only thing that re-defines m_particles/m_grownGrids/
  // m_bufferParticles over the new layout and resets m_finestLevel. Skip its preRegrid() and the
  // rebuild in regrid() would index new-layout DataIndexes into old-layout LevelData.
  m_regridParticles.preRegrid();

  for (auto& container : m_particleContainers) {
    ParticleContainer<ItoParticle>& particles = container.second;

    particles.preRegrid();
  }

  m_phi.clear();
  m_mobilityFunction.clear();
  m_velocityFunction.clear();
  m_diffusionFunction.clear();
  m_depositionNC.clear();
  m_massDiff.clear();
  m_kdMergeCellHistogram.clear();
  m_kdMergeLeafQuota.clear();
}

ParticleContainer<ItoParticle>&
ItoSolver::getParticles(const WhichContainer a_container)
{
  CH_TIME("ItoSolver::getParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::getParticles" << endl;
  }

  return m_particleContainers.at(a_container);
}

const ParticleContainer<ItoParticle>&
ItoSolver::getParticles(const WhichContainer a_container) const
{
  CH_TIME("ItoSolver::getParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::getParticles" << endl;
  }

  return m_particleContainers.at(a_container);
}

EBAMRCellData&
ItoSolver::getPhi()
{
  CH_TIME("ItoSolver::getPhi");
  if (m_verbosity > 5) {
    pout() << m_name + "::getPhi" << endl;
  }

  return m_phi;
}

EBAMRCellData&
ItoSolver::getVelocityFunction()
{
  CH_TIME("ItoSolver::getVelocityFunction");
  if (m_verbosity > 5) {
    pout() << m_name + "::getVelocityFunction" << endl;
  }

  return m_velocityFunction;
}

const EBAMRCellData&
ItoSolver::getVelocityFunction() const
{
  CH_TIME("ItoSolver::getVelocityFunction");
  if (m_verbosity > 5) {
    pout() << m_name + "::getVelocityFunction" << endl;
  }

  return m_velocityFunction;
}

EBAMRCellData&
ItoSolver::getDiffusionFunction()
{
  CH_TIME("ItoSolver::getDiffusionFunction");
  if (m_verbosity > 5) {
    pout() << m_name + "::getDiffusionFunction" << endl;
  }

  return m_diffusionFunction;
}

const EBAMRCellData&
ItoSolver::getDiffusionFunction() const
{
  CH_TIME("ItoSolver::getDiffusionFunction");
  if (m_verbosity > 5) {
    pout() << m_name + "::getDiffusionFunction" << endl;
  }

  return m_diffusionFunction;
}

EBAMRCellData&
ItoSolver::getMobilityFunction()
{
  CH_TIME("ItoSolver::getMobilityFunction");
  if (m_verbosity > 5) {
    pout() << m_name + "::getMobilityFunction" << endl;
  }

  return m_mobilityFunction;
}

const EBAMRCellData&
ItoSolver::getMobilityFunction() const
{
  CH_TIME("ItoSolver::getMobilityFunction");
  if (m_verbosity > 5) {
    pout() << m_name + "::getMobilityFunction" << endl;
  }

  return m_mobilityFunction;
}

void
ItoSolver::setDiffusionFunction(const Real a_diffusionCoefficient)
{
  CH_TIME("ItoSolver::setDiffusionFunction");
  if (m_verbosity > 5) {
    pout() << m_name + "::setDiffusionFunction" << endl;
  }

  DataOps::setValue(m_diffusionFunction, a_diffusionCoefficient);
}

void
ItoSolver::setVelocityFunction(const RealVect& a_velocity)
{
  CH_TIME("ItoSolver::setVelocityFunction");
  if (m_verbosity > 5) {
    pout() << m_name + "::setVelocityFunction" << endl;
  }

  for (int dir = 0; dir < SpaceDim; dir++) {
    DataOps::setValue(m_velocityFunction, a_velocity[dir], dir);
  }
}

void
ItoSolver::setParticleMobility(const Real a_mobility)
{
  CH_TIME("ItoSolver::setParticleMobility");
  if (m_verbosity > 5) {
    pout() << m_name + "::setParticleMobility" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(WhichContainer::Bulk);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      ParticleSoA<ItoParticle>& leaf = particles[lvl][din];

      ParticleReal* mobility = leaf.column<&ItoParticle::mobility>();

      ParticleLoops::loop(leaf, [&](const std::size_t i) {
        mobility[i] = a_mobility;
      });
    }
  }
}

void
ItoSolver::setParticleDiffusion(const Real a_diffCo)
{
  CH_TIME("ItoSolver::setParticleDiffusion");
  if (m_verbosity > 5) {
    pout() << m_name + "::setParticleDiffusion" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(WhichContainer::Bulk);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      ParticleSoA<ItoParticle>& leaf = particles[lvl][din];

      ParticleReal* diffusion = leaf.column<&ItoParticle::diffusion>();

      ParticleLoops::loop(leaf, [&](const std::size_t i) {
        diffusion[i] = a_diffCo;
      });
    }
  }
}

void
ItoSolver::interpolateVelocities()
{
  CH_TIME("ItoSolver::interpolateVelocities");
  if (m_verbosity > 5) {
    pout() << m_name + "::interpolateVelocities" << endl;
  }

  if (m_isMobile) {
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        this->interpolateVelocities(lvl, din);
      }
    }
  }
}

void
ItoSolver::interpolateVelocities(const int a_lvl, const DataIndex& a_dit)
{
  CH_TIME("ItoSolver::interpolateVelocities");
  if (m_verbosity > 5) {
    pout() << m_name + "::interpolateVelocities" << endl;
  }

  ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);

  if (m_isMobile) {
    const ProblemDomain& domain  = m_amr->getDomains()[a_lvl];
    const Box            cellBox = m_amr->getGrids(m_realm)[a_lvl][a_dit];
    const EBISBox&       ebisbox = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl][a_dit];
    const Real           dx      = m_amr->getDx()[a_lvl];
    const RealVect       probLo  = m_amr->getProbLo();

    const EBCellFAB& velo_func = (*m_velocityFunction[a_lvl])[a_dit];

    ParticleSoA<ItoParticle>& leaf = particles[a_lvl][a_dit];

    // Interpolate the velocity function onto the particle velocity columns.
    EBParticleMesh meshInterp(domain, cellBox, ebisbox, dx * RealVect::Unit, probLo);

    meshInterp.interpolate<D_DECL(&ItoParticle::vx, &ItoParticle::vy, &ItoParticle::vz)>(leaf,
                                                                                         velo_func,
                                                                                         m_deposition,
                                                                                         m_forceIrregInterpolationNGP);

    // Set the particle velocities to velo_func * mobility.
    const ParticleReal* mobility    = leaf.column<&ItoParticle::mobility>();
    ParticleReal*       v[SpaceDim] = {
      D_DECL(leaf.column<&ItoParticle::vx>(), leaf.column<&ItoParticle::vy>(), leaf.column<&ItoParticle::vz>())};

    ParticleLoops::loop(leaf, [&](const std::size_t i) {
      for (int dir = 0; dir < SpaceDim; dir++) {
        v[dir][i] *= mobility[i];
      }
    });
  }
}

void
ItoSolver::interpolateMobilities()
{
  CH_TIME("ItoSolver::interpolateMobilities()");
  if (m_verbosity > 5) {
    pout() << m_name + "::interpolateMobilities()" << endl;
  }

  if (m_isMobile) {
    EBAMRCellData velocityMagnitude;
    m_amr->allocate(velocityMagnitude, m_realm, m_phase, 1);

    switch (m_mobilityInterp) {
    case WhichMobilityInterpolation::Velocity: {

      // Compute |v|
      DataOps::vectorLength(velocityMagnitude,
                            m_velocityFunction,
                            m_amr->getNotCoveredCells(m_realm, m_phase),
                            m_amr->getMultiCutVofIterator(m_realm, m_phase));

      m_amr->conservativeAverage(velocityMagnitude, m_realm, m_phase);
      m_amr->interpGhostPwl(velocityMagnitude, m_realm, m_phase);

      break;
    }
    default: // Do nothing
      break;
    }

    // Call the level version and interpolate the mobilities from the mesh data.
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        this->interpolateMobilities(lvl, din, (*velocityMagnitude[lvl])[din]);
      }
    }
  }
}

void
ItoSolver::interpolateMobilities(const int a_lvl, const DataIndex& a_dit, const EBCellFAB& a_velocityMagnitude) noexcept
{
  CH_TIME("ItoSolver::interpolateMobilities(lvl, patch)");
  if (m_verbosity > 5) {
    pout() << m_name + "::interpolateMobilities(lvl, patch)" << endl;
  }

  CH_assert(m_isMobile);

  switch (m_mobilityInterp) {
  case WhichMobilityInterpolation::Direct: {
    this->interpolateMobilitiesDirect(a_lvl, a_dit);

    break;
  }
  case WhichMobilityInterpolation::Velocity: {
    this->interpolateMobilitiesVelocity(a_lvl, a_dit, a_velocityMagnitude);

    break;
  }
  default: {
    MayDay::Error("ItoSolver::interpolateMobilities(int, DataIndex) - logic bust");

    break;
  }
  }
}

void
ItoSolver::interpolateMobilitiesDirect(const int a_lvl, const DataIndex& a_dit) noexcept
{
  CH_TIME("ItoSolver::interpolateMobilitiesDirect");
  if (m_verbosity > 5) {
    pout() << m_name + "::interpolateMobilitiesDirect" << endl;
  }

  CH_assert(m_isMobile);
  CH_assert(m_mobilityInterp == WhichMobilityInterpolation::Direct);

  // TLDR: This will compute the particle mobility by interpolating a scalar mobility field (stored on the mesh) to the
  // particle positions.

  ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);

  const ProblemDomain& domain  = m_amr->getDomains()[a_lvl];
  const Box            cellBox = m_amr->getGrids(m_realm)[a_lvl][a_dit];
  const EBISBox&       ebisbox = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl][a_dit];
  const Real           dx      = m_amr->getDx()[a_lvl];
  const RealVect       probLo  = m_amr->getProbLo();

  const EBCellFAB& mobilityFunction = (*m_mobilityFunction[a_lvl])[a_dit];

  ParticleSoA<ItoParticle>& leaf = particles[a_lvl][a_dit];

  // Interpolate the scalar mobility field onto the particle mobility column.
  EBParticleMesh meshInterp(domain, cellBox, ebisbox, dx * RealVect::Unit, probLo);

  meshInterp.interpolate<&ItoParticle::mobility>(leaf, mobilityFunction, m_deposition, m_forceIrregInterpolationNGP);
}

void
ItoSolver::interpolateMobilitiesVelocity(const int        a_lvl,
                                         const DataIndex& a_dit,
                                         const EBCellFAB& a_velocityMagnitude) noexcept
{
  CH_TIME("ItoSolver::interpolateMobilitiesVelocity");
  if (m_verbosity > 5) {
    pout() << m_name + "::interpolateMobilitiesVelocity" << endl;
  }

  CH_assert(m_isMobile);
  CH_assert(a_velocityMagnitude.nComp() == 1);
  CH_assert(m_mobilityInterp == WhichMobilityInterpolation::Velocity);

  // clang-format off
  // TLDR: This function computes the particle mobilities by interpolating mu*V to the particle position and then setting
  //       the mobility as mu = [mu*V(Xp)]/V(Xp).
  // clang-format on

  const ProblemDomain& domain  = m_amr->getDomains()[a_lvl];
  const Box            cellBox = m_amr->getGrids(m_realm)[a_lvl][a_dit];
  const EBISBox&       ebisbox = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl][a_dit];
  const Real           dx      = m_amr->getDx()[a_lvl];
  const RealVect       probLo  = m_amr->getProbLo();

  const EBCellFAB& mobilityFunction = (*m_mobilityFunction[a_lvl])[a_dit];

  ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);
  ParticleSoA<ItoParticle>&       leaf      = particles[a_lvl][a_dit];

  EBParticleMesh meshInterp(domain, cellBox, ebisbox, dx * RealVect::Unit, probLo);

  // Compute mu*|V| on the mesh
  EBCellFAB muV;
  muV.clone(a_velocityMagnitude);
  muV *= mobilityFunction;

  // First, interpolate |V| to the particle position, stored in the scratch column.
  meshInterp.interpolate<&ItoParticle::scratch>(leaf, a_velocityMagnitude, m_deposition, m_forceIrregInterpolationNGP);

  meshInterp.interpolate<&ItoParticle::mobility>(leaf, muV, m_deposition, m_forceIrregInterpolationNGP);

  // We now have scratch = |V(Xp)| and mobility = |mu*V|(Xp). Set mobility(Xp) = |mu*V|(Xp)/|V|(Xp).
  ParticleReal*       mobility = leaf.column<&ItoParticle::mobility>();
  const ParticleReal* scratch  = leaf.column<&ItoParticle::scratch>();

  ParticleLoops::loop(leaf, [&](const std::size_t i) {
    mobility[i] *= 1.0 / scratch[i];
  });
}

void
ItoSolver::updateMobilities()
{
  CH_TIME("ItoSolver::updateMobilities");
  if (m_verbosity > 5) {
    pout() << m_name + "::updateMobilities" << endl;
  }

  // TLDR: This routine is for computing mobilities as mu = mu(e) where e is the energy. This is done
  //       via the Ito species.

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      this->updateMobilities(lvl, din);
    }
  }
}

void
ItoSolver::updateMobilities(const int a_level, const DataIndex& a_dit)
{
  CH_TIME("ItoSolver::updateMobilities(int, DataIndex)");
  if (m_verbosity > 5) {
    pout() << m_name + "::updateMobilities(int, DataIndex)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);

  if (m_isMobile) {
    ParticleSoA<ItoParticle>& leaf = particles[a_level][a_dit];

    const std::size_t   n        = leaf.size();
    ParticleReal*       mobility = leaf.column<&ItoParticle::mobility>();
    const ParticleReal* energy   = leaf.column<&ItoParticle::energy>();

    for (std::size_t i = 0; i < n; i++) {
      mobility[i] = m_species->mobility(energy[i]);
    }
  }
}

void
ItoSolver::interpolateDiffusion()
{
  CH_TIME("ItoSolver::interpolateDiffusion");
  if (m_verbosity > 5) {
    pout() << m_name + "::interpolateDiffusion" << endl;
  }

  if (m_isDiffusive) {
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        this->interpolateDiffusion(lvl, din);
      }
    }
  }
}

void
ItoSolver::interpolateDiffusion(const int a_lvl, const DataIndex& a_dit)
{
  CH_TIME("ItoSolver::interpolateDiffusion");
  if (m_verbosity > 5) {
    pout() << m_name + "::interpolateDiffusion" << endl;
  }

  if (m_isDiffusive) {

    // These are the particles that will be interpolated.
    ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);
    ParticleSoA<ItoParticle>&       leaf      = particles[a_lvl][a_dit];

    const ProblemDomain& domain  = m_amr->getDomains()[a_lvl];
    const Box            cellBox = m_amr->getGrids(m_realm)[a_lvl][a_dit];
    const EBISBox&       ebisbox = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl][a_dit];
    const Real           dx      = m_amr->getDx()[a_lvl];
    const RealVect       probLo  = m_amr->getProbLo();

    // Create the particle interpolator and interpolate the diffusion field onto the diffusion column.
    const EBCellFAB& Dcoef = (*m_diffusionFunction[a_lvl])[a_dit];

    EBParticleMesh meshInterp(domain, cellBox, ebisbox, dx * RealVect::Unit, probLo);

    meshInterp.interpolate<&ItoParticle::diffusion>(leaf, Dcoef, m_deposition, m_forceIrregInterpolationNGP);
  }
}

void
ItoSolver::updateDiffusion()
{
  CH_TIME("ItoSolver::updateDiffusion");
  if (m_verbosity > 5) {
    pout() << m_name + "::updateDiffusion" << endl;
  }

  // TLDR: This routine is for updating the diffusion as D = D(e) where e is the energy. This is done
  //       via ItoSpecies

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      this->updateDiffusion(lvl, din);
    }
  }
}

void
ItoSolver::updateDiffusion(const int a_level, const DataIndex& a_dit)
{
  CH_TIME("ItoSolver::updateDiffusion(lvl, dit)");
  if (m_verbosity > 5) {
    pout() << m_name + "::updateDiffusion(lvl, dit)" << endl;
  }

  if (m_isDiffusive) {

    // Particles to be updated.
    ParticleContainer<ItoParticle>& particles = this->getParticles(WhichContainer::Bulk);
    ParticleSoA<ItoParticle>&       leaf      = particles[a_level][a_dit];

    const std::size_t   n         = leaf.size();
    ParticleReal*       diffusion = leaf.column<&ItoParticle::diffusion>();
    const ParticleReal* energy    = leaf.column<&ItoParticle::energy>();

    // Update diffusion coefficient.
    for (std::size_t i = 0; i < n; i++) {
      diffusion[i] = m_species->diffusion(energy[i]);
    }
  }
}

Real
ItoSolver::computeDt() const
{
  CH_TIME("ItoSolver::computeDt()");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDt()" << endl;
  }

  Real dt = std::numeric_limits<Real>::max();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Real levelDt = this->computeDt(lvl);

    dt = std::min(dt, levelDt);
  }

  return dt;
}

Real
ItoSolver::computeDt(const int a_lvl) const
{
  CH_TIME("ItoSolver::computeDt(int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDt(int)" << endl;
  }

  Real dt = std::numeric_limits<Real>::max();

  // Compute largest permitted time step on each patch.
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const DataIterator&      dit = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(min : dt)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const Real patchDt = this->computeDt(a_lvl, din);

    dt = std::min(dt, patchDt);
  }

  return ParallelOps::min(dt);
}

Real
ItoSolver::computeDt(const int a_lvl, const DataIndex& a_dit) const
{
  CH_TIME("ItoSolver::computeDt(int, DataIndex)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDt(int, DataIndex)" << endl;
  }

  // Return value.
  Real dt = std::numeric_limits<Real>::max();

  // Grid resolution on this level.
  const Real dx = m_amr->getDx()[a_lvl];

  // Particles to iterate over.
  const ParticleContainer<ItoParticle>& particles = this->getParticles(WhichContainer::Bulk);
  const ParticleSoA<ItoParticle>&       leaf      = particles[a_lvl][a_dit];

  const std::size_t   n           = leaf.size();
  const ParticleReal* v[SpaceDim] = {
    D_DECL(leaf.column<&ItoParticle::vx>(), leaf.column<&ItoParticle::vy>(), leaf.column<&ItoParticle::vz>())};
  const ParticleReal* D = leaf.column<&ItoParticle::diffusion>();

  if (m_isMobile && !m_isDiffusive) {

    // Advection but no diffusion - set the time step as dt = dx/vMax where vMax is the largest velocity component.
    for (std::size_t i = 0; i < n; i++) {
      Real vMax = 0.0;
      for (int dir = 0; dir < SpaceDim; dir++) {
        vMax = std::max(vMax, std::abs(static_cast<Real>(v[dir][i])));
      }

      const Real thisDt = (vMax > 0.0) ? dx / vMax : std::numeric_limits<Real>::max();

      dt = std::min(dt, thisDt);
    }
  }
  else if (!m_isMobile && m_isDiffusive) {
    // Diffusion but no advection -- set the time step as dt = dx*dx/(2*D)

    for (std::size_t i = 0; i < n; i++) {
      const Real Di     = static_cast<Real>(D[i]);
      const Real thisDt = (Di > 0.0) ? dx * dx / (2.0 * SpaceDim * Di) : std::numeric_limits<Real>::max();

      dt = std::min(dt, thisDt);
    }
  }
  else if (m_isMobile && m_isDiffusive) {

    // Both advection and diffusion. Compute dt = 1/(1/dtA + 1/dtD) where dtA and dtD are as in the code bits above.
    for (std::size_t i = 0; i < n; i++) {
      Real vMax = 0.0;
      for (int dir = 0; dir < SpaceDim; dir++) {
        vMax += std::abs(static_cast<Real>(v[dir][i]));
      }

      const Real Di = static_cast<Real>(D[i]);

      const Real dtAdvect  = (vMax > 0.0) ? dx / vMax : std::numeric_limits<Real>::max();
      const Real dtDiffuse = (Di > 0.0) ? dx * dx / (2.0 * SpaceDim * Di) : std::numeric_limits<Real>::max();

      const Real thisDt = 1. / (1. / dtAdvect + 1. / dtDiffuse);

      dt = std::min(dt, thisDt);
    }
  }

  return dt;
}

Real
ItoSolver::computeHopDt(const Real a_maxCellsToMove) const
{
  CH_TIME("ItoSolver::computeHopDt(Real)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeHopDt(Real)" << endl;
  }

  CH_assert(a_maxCellsToMove > 0.0);

  // clang-format off
  // TLDR: This routine computes the largest possible time step such that no particles move more than a_maxCellsToMove during
  //       a standard Ito kernel step. This is the AMR version.
  // clang-format on

  Real dt = std::numeric_limits<Real>::max();

  // Compute time steps for each grid level.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Real levelDt = this->computeHopDt(a_maxCellsToMove, lvl);

    dt = std::min(dt, levelDt);
  }

  return dt;
}

Real
ItoSolver::computeHopDt(const Real a_maxCellsToMove, const int a_lvl) const
{
  CH_TIME("ItoSolver::computeHopDt(Real, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeHopDt(Real, int)" << endl;
  }

  CH_assert(a_maxCellsToMove > 0.0);

  // clang-format off
  // TLDR: This routine computes the largest possible time step such that no particles move more than a_maxCellsToMove during
  //       a standard Ito kernel step. This is the level version.
  // clang-format on

  Real dt = std::numeric_limits<Real>::max();

  // Compute time steps for each grid patch.
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const DataIterator&      dit = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(min : dt)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const Real patchDt = this->computeHopDt(a_maxCellsToMove, a_lvl, din);

    dt = std::min(dt, patchDt);
  }

  return ParallelOps::min(dt);
}

Real
ItoSolver::computeHopDt(const Real a_maxCellsToMove, const int a_lvl, const DataIndex& a_dit) const
{
  CH_TIME("ItoSolver::computeHopDt(Real, int, DataIndex)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeHopDt(Real, int, DataIndex)" << endl;
  }

  CH_assert(a_maxCellsToMove > 0.0);

  // clang-format off
  // TLDR: This routine computes the largest possible time step such that no particles move more than a_maxCellsToMove during
  //       a standard Ito kernel step. This is the patch version.
  // clang-format on

  Real dt = std::numeric_limits<Real>::max();

  // Grid resolution on this level.
  const Real dx = m_amr->getDx()[a_lvl];

  // Some shortcuts.
  const Real dMax  = a_maxCellsToMove * dx;
  const Real dMax2 = dMax * dMax;
  const Real W0    = m_normalDistributionTruncation;
  const Real W02   = m_normalDistributionTruncation * m_normalDistributionTruncation;

  // These are the particles we will compute the time steps for.
  const ParticleContainer<ItoParticle>& particles = this->getParticles(WhichContainer::Bulk);
  const ParticleSoA<ItoParticle>&       leaf      = particles[a_lvl][a_dit];

  const std::size_t   n           = leaf.size();
  const ParticleReal* v[SpaceDim] = {
    D_DECL(leaf.column<&ItoParticle::vx>(), leaf.column<&ItoParticle::vy>(), leaf.column<&ItoParticle::vz>())};
  const ParticleReal* D = leaf.column<&ItoParticle::diffusion>();

  // Helper that returns the largest velocity component magnitude of particle i.
  auto vMaxOf = [&](const std::size_t a_i) -> Real {
    Real vMax = 0.0;
    for (int dir = 0; dir < SpaceDim; dir++) {
      vMax = std::max(vMax, std::abs(static_cast<Real>(v[dir][a_i])));
    }
    return vMax;
  };

  // Cases handled differently
  if (m_isMobile && !m_isDiffusive) {
    // Advection but no diffusion

    for (std::size_t i = 0; i < n; i++) {
      // Compute the regular time step -- recall that dMax = a_maxCellsToMove*dx
      const Real thisDt = dMax / vMaxOf(i);

      dt = std::min(dt, thisDt);
    }
  }
  else if (!m_isMobile && m_isDiffusive) {
    // Diffusion but no advection

    for (std::size_t i = 0; i < n; i++) {
      // Recall, the diffusion kernel is usually dX = sqrt(2*D*dt)*N where N is a SpaceDim vector of Gaussian numbers.
      // But we only care about not moving more than a specified number of grid cells in any one of the coordinate
      // directions so we have |dX| = sqrt(2*D*dt)*N0 where N0 is the maximum value in the Gaussian distribution (which
      // we have truncated). Solving for dt yields dt = |dX|^2/(2*D*N0^2).
      const Real Di     = static_cast<Real>(D[i]);
      const Real thisDt = dMax2 / (2.0 * Di * SpaceDim * W02);

      dt = std::min(dt, thisDt);
    }
  }
  else if (m_isMobile && m_isDiffusive) {

    // Diffusion AND advection. Much more difficult and requires us to solve a second order equation.
    for (std::size_t i = 0; i < n; i++) {
      const Real vMax = vMaxOf(i);
      const Real Di   = static_cast<Real>(D[i]);

      Real thisDt = std::numeric_limits<Real>::max();

      // This case is more complicated because we have dX = v*dt + sqrt(2*D*dt)*N0 and we need to solve for dt. This
      // yields a second order equation in the form A*dt^2 + B*dt + C = 0. We just solve for dt.
      if (vMax > 0.0) {
        const Real a = vMax;
        const Real b = W0 * sqrt(2.0 * Di);
        const Real c = dMax;

        const Real A = a * a;
        const Real B = -(b * b + 2 * a * c);
        const Real C = c * c;

        thisDt = (-B - sqrt(B * B - 4. * A * C)) / (2. * A);
      }
      else {
        if (Di > 0.0) {
          thisDt = dMax2 / (2.0 * Di * W02);
        }
      }

      dt = std::min(dt, thisDt);
    }
  }

  return dt;
}

Real
ItoSolver::computeAdvectiveDt() const
{
  CH_TIME("ItoSolver::computeAdvectiveDt");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeAdvectiveDt()" << endl;
  }

  // TLDR: We compute dt = dx/vMax for every particle.

  Real dt = std::numeric_limits<Real>::max();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Real levelDt = this->computeAdvectiveDt(lvl);

    dt = std::min(levelDt, dt);
  }

  return dt;
}

Real
ItoSolver::computeAdvectiveDt(const int a_lvl) const
{
  CH_TIME("ItoSolver::computeAdvectiveDt(int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeAdvectiveDt(int)" << endl;
  }

  CH_assert(a_lvl >= 0);

  // TLDR: We compute dt = dx/vMax on each grid patch on this level.

  Real dt = std::numeric_limits<Real>::max();

  // Iterate over patches.
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const DataIterator&      dit = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(min : dt)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const Real patchDt = this->computeAdvectiveDt(a_lvl, din);

    dt = std::min(dt, patchDt);
  }

  return ParallelOps::min(dt);
}

Real
ItoSolver::computeAdvectiveDt(const int a_lvl, const DataIndex& a_dit) const
{
  CH_TIME("ItoSolver::computeAdvectiveDt(int, DataIndex)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeAdvectiveDt(int, DataIndex, dx)" << endl;
  }

  Real dt = std::numeric_limits<Real>::max();

  if (m_isMobile) {

    const Real dx = m_amr->getDx()[a_lvl];

    // Particles that we iterate over.
    const ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);
    const ParticleSoA<ItoParticle>&       leaf      = particles[a_lvl][a_dit];

    const std::size_t   n           = leaf.size();
    const ParticleReal* v[SpaceDim] = {
      D_DECL(leaf.column<&ItoParticle::vx>(), leaf.column<&ItoParticle::vy>(), leaf.column<&ItoParticle::vz>())};

    for (std::size_t i = 0; i < n; i++) {
      // Get maximum velocity component magnitude.
      Real vMax = 0.0;
      for (int dir = 0; dir < SpaceDim; dir++) {
        vMax = std::max(vMax, std::abs(static_cast<Real>(v[dir][i])));
      }

      const Real thisDt = (vMax > 0.0) ? (dx / vMax) : std::numeric_limits<Real>::max();

      dt = std::min(dt, thisDt);
    }
  }

  return dt;
}

Real
ItoSolver::computeDiffusiveDt() const
{
  CH_TIME("ItoSolver::computeDiffusiveDt()");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDiffusiveDt()" << endl;
  }

  // TLDR: Compute dt = dx*dx/(2*D) on each grid patchon every grid level.

  Real dt = std::numeric_limits<Real>::max();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const Real levelDt = this->computeDiffusiveDt(lvl);

    dt = std::min(dt, levelDt);
  }

  return dt;
}

Real
ItoSolver::computeDiffusiveDt(const int a_lvl) const
{
  CH_TIME("ItoSolver::computeDiffusiveDt(int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDiffusiveDt(int)" << endl;
  }

  CH_assert(a_lvl >= 0);

  // TLDR: Compute dt = dx*dx/(2*D) on each grid patch.

  Real dt = std::numeric_limits<Real>::max();

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const DataIterator&      dit = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(min : dt)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const Real patchDt = this->computeDiffusiveDt(a_lvl, din);

    dt = std::min(dt, patchDt);
  }

  return ParallelOps::min(dt);
}

Real
ItoSolver::computeDiffusiveDt(const int a_lvl, const DataIndex& a_dit) const
{
  CH_TIME("ItoSolver::computeDiffusiveDt(int, DataIndex)");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeDiffusiveDt(int, DataIndex)" << endl;
  }

  CH_assert(a_lvl >= 0);

  // TLDR: Compute dt = dx*dx/(2*D) for all particles in the input grid patch.

  Real dt = std::numeric_limits<Real>::max();

  if (m_isDiffusive) {
    const Real dx  = m_amr->getDx()[a_lvl];
    const Real dx2 = dx * dx / (2.0 * SpaceDim);

    // These are the particles we iterate over.
    const ParticleContainer<ItoParticle>& particles = this->getParticles(WhichContainer::Bulk);
    const ParticleSoA<ItoParticle>&       leaf      = particles[a_lvl][a_dit];

    const std::size_t   n = leaf.size();
    const ParticleReal* D = leaf.column<&ItoParticle::diffusion>();

    for (std::size_t i = 0; i < n; i++) {
      const Real Di = static_cast<Real>(D[i]);

      const Real thisDt = (Di > 0.0) ? dx2 / Di : std::numeric_limits<Real>::max();

      dt = std::min(dt, thisDt);
    }
  }

  return dt;
}

void
ItoSolver::remapAll()
{
  CH_TIME("ItoSolver::remapAll");
  if (m_verbosity > 5) {
    pout() << m_name + "::remapAll" << endl;
  }

  for (auto& container : m_particleContainers) {
    container.second.remap();
  }
}

void
ItoSolver::remap()
{
  CH_TIME("ItoSolver::remap");
  if (m_verbosity > 5) {
    pout() << m_name + "::remap" << endl;
  }

  this->remap(WhichContainer::Bulk);
}

void
ItoSolver::remap(const WhichContainer a_container)
{
  CH_TIME("ItoSolver::remap(container)");
  if (m_verbosity > 5) {
    pout() << m_name + "::remap(container)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(a_container);

  particles.remap();
}

DepositionType
ItoSolver::getDeposition() const
{
  return m_deposition;
}

CoarseFineDeposition
ItoSolver::getCoarseFineDeposition() const
{
  return m_coarseFineDeposition;
}

phase::which_phase
ItoSolver::getPhase() const
{
  return m_phase;
}

void
ItoSolver::organizeParticlesByCell(const WhichContainer a_container)
{
  CH_TIME("ItoSolver::organizeParticlesByCell");
  if (m_verbosity > 5) {
    pout() << m_name + "::organizeParticlesByCell" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(a_container);

  particles.organizeParticlesByCell();
}

void
ItoSolver::organizeParticlesByPatch(const WhichContainer a_container)
{
  CH_TIME("ItoSolver::organizeParticlesByPatch");
  if (m_verbosity > 5) {
    pout() << m_name + "::organizeParticlesByPatch" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(a_container);

  particles.organizeParticlesByPatch();
}

void
ItoSolver::makeSuperparticles(const WhichContainer a_container)
{
  CH_TIME("ItoSolver::makeSuperparticles(WhichContainer)");
  if (m_verbosity > 5) {
    pout() << m_name + "::makeSuperparticles(WhichContainer)" << endl;
  }

  this->makeSuperparticles(a_container, m_particlesPerCell);
}

void
ItoSolver::makeSuperparticles(const WhichContainer a_container, const int a_particlesPerCell)
{
  CH_TIME("ItoSolver::makeSuperparticles(WhichContainer, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::makeSuperparticles(WhichContainer, int)" << endl;
  }

  // Uniform target across all levels -- forward to the per-level entry point.
  this->makeSuperparticles(a_container, Vector<int>(1, a_particlesPerCell));
}

void
ItoSolver::makeSuperparticles(const WhichContainer a_container, const Vector<int>& a_particlesPerCell)
{
  CH_TIME("ItoSolver::makeSuperparticles(WhichContainer, Vector<int>)");
  if (m_verbosity > 5) {
    pout() << m_name + "::makeSuperparticles(WhichContainer, Vector<int>)" << endl;
  }

  this->makeSuperparticles(a_container, a_particlesPerCell, m_mergeMethod);
}

void
ItoSolver::makeSuperparticles(const WhichContainer                          a_container,
                              const Vector<int>&                            a_particlesPerCell,
                              const ParticleManagement::ParticleMergeMethod a_method)
{
  CH_TIME("ItoSolver::makeSuperparticles(WhichContainer, Vector<int>, ParticleMergeMethod)");
  if (m_verbosity > 5) {
    pout() << m_name + "::makeSuperparticles(WhichContainer, Vector<int>, ParticleMergeMethod)" << endl;
  }

  if (a_particlesPerCell.size() < 1) {
    MayDay::Error("ItoSolver::makeSuperparticles(WhichContainer, Vector<int>) - empty particles-per-cell vector");
  }

  // 'none' short-circuits before the extract. Every other method runs on ItoMergeParticle, so the
  // pipeline is always extract -> mergeLite -> rebuild; running that round-trip for a no-op merge would
  // pay a full population copy each way and zero the ItoParticle payload for nothing. The regrid path is
  // the one place the round-trip DOES run under 'none' (ItoSolver::regrid), and there it is unavoidable
  // -- preRegrid already did the extract, so the rebuild has to happen regardless.
  if (a_method == ParticleManagement::ParticleMergeMethod::None) {
    return;
  }

  ParticleContainer<ItoMergeParticle> merge;

  this->extractIntoMergeContainer(a_container, merge);
  this->mergeLite(merge, a_particlesPerCell, a_method);
  this->rebuildFromMergeContainer(a_container, merge);
}

void
ItoSolver::mergeLite(ParticleContainer<ItoMergeParticle>&          a_particles,
                     const Vector<int>&                            a_particlesPerCell,
                     const ParticleManagement::ParticleMergeMethod a_method)
{
  CH_TIME("ItoSolver::mergeLite");
  if (m_verbosity > 5) {
    pout() << m_name + "::mergeLite" << endl;
  }

  // TLDR: this is the ONE place that knows which methods split and which do not, so no caller has to.
  //       The cell-based methods reduce a cell to at most its target and never refill an under-full one,
  //       so splitting after them would inflate the population; the whole-container methods only drain
  //       over-full cells and rely on the split to bring under-full ones back up. 'none' does neither.
  //       Getting this wrong is silent -- the population simply grows at every regrid.
  switch (a_method) {
  case ParticleManagement::ParticleMergeMethod::None: {
    // Deliberately nothing: not even a split. See the regrid path in ItoSolver::regrid().
    return;
  }

    // ---- cell-based: per-level target honoured, never splits ----
  case ParticleManagement::ParticleMergeMethod::EqualWeightKD: {
    this->mergeEqualWeightKD(a_particles, a_particlesPerCell);

    return;
  }
  case ParticleManagement::ParticleMergeMethod::Reinitialize: {
    this->mergeReinitialize(a_particles, a_particlesPerCell);

    return;
  }
  case ParticleManagement::ParticleMergeMethod::ReinitializeBVH: {
    this->mergeReinitializeBVH(a_particles, a_particlesPerCell);

    return;
  }
  case ParticleManagement::ParticleMergeMethod::NnSfc: {
    this->mergeNnSfc(a_particles, a_particlesPerCell);

    return;
  }
  case ParticleManagement::ParticleMergeMethod::External: {
    // The user supplies the per-cell merger through setParticleCellMerger().
    this->applyCellMerger(a_particles, a_particlesPerCell, m_particleCellMerger);

    return;
  }

    // ---- whole-container: one crowding threshold over the hierarchy, splits afterwards ----
  case ParticleManagement::ParticleMergeMethod::NnPairTree: {
    this->mergeNnPairSearch<NnPairSearchBackend::Tree>(a_particles, a_particlesPerCell);

    break;
  }
  case ParticleManagement::ParticleMergeMethod::NnPairHash: {
    this->mergeNnPairSearch<NnPairSearchBackend::Hash>(a_particles, a_particlesPerCell);

    break;
  }
  case ParticleManagement::ParticleMergeMethod::NnPairOneCell: {
    this->mergeNnPairOneCell(a_particles, a_particlesPerCell);

    break;
  }
  case ParticleManagement::ParticleMergeMethod::KdCarve: {
    this->mergeKDImpl(a_particles, a_particlesPerCell, true);

    break;
  }
  case ParticleManagement::ParticleMergeMethod::KdPatch: {
    this->mergeKDImpl(a_particles, a_particlesPerCell, false);

    break;
  }
  case ParticleManagement::ParticleMergeMethod::KdSkinNn: {
    this->mergeKDSkinNn(a_particles, a_particlesPerCell);

    break;
  }
  default: {
    MayDay::Abort("ItoSolver::mergeLite -- unsupported method requested");

    return;
  }
  }

  // Only the whole-container methods reach here. They mint their own rank-local ids from a counter that
  // restarts every drain round (see the allocateID() lambdas below), so what is sitting in the id column
  // now is merge scratch, not particle identity -- and it is indistinguishable from a real id downstream,
  // where it would reach H5Part output. Put the column back to "no id", which is what an append()-built
  // particle carries and what a merged particle has always effectively had.
  this->invalidateParticleIDs(a_particles);

  this->splitFromMergeContainer(a_particles, a_particlesPerCell[0]);
}

void
ItoSolver::invalidateParticleIDs(ParticleContainer<ItoMergeParticle>& a_particles) const
{
  CH_TIME("ItoSolver::invalidateParticleIDs");

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit  = dbl.dataIterator();
    const int                nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      ParticleSoA<ItoMergeParticle>& leaf = a_particles[lvl][dit[mybox]];

      for (std::size_t i = 0; i < leaf.size(); i++) {
        leaf.particleID(i) = ParticleSoA<ItoMergeParticle>::s_invalidID;
      }
    }
  }
}

void
ItoSolver::applyCellMerger(ParticleContainer<ItoMergeParticle>&                        a_particles,
                           const Vector<int>&                                          a_particlesPerCell,
                           const ParticleManagement::ParticleMerger<ItoMergeParticle>& a_merger)
{
  CH_TIME("ItoSolver::applyCellMerger");
  if (m_verbosity > 5) {
    pout() << m_name + "::applyCellMerger" << endl;
  }

  // TLDR (SoA): cell-sort so each leaf's cells own a contiguous CSR range [cellStart(c), cellStart(c+1)).
  //       For each non-empty cell we extract its particles into an SoA scratch, run a_merger on it, and
  //       collect the result; then we rebuild the leaf and return the container to a by-patch view.
  //
  // Operates on the reduced particle: everything a cell-based merge needs -- position, weight, energy --
  // is already there, and the ItoParticle columns it would otherwise discard (velocity, mobility,
  // diffusion) were never carried in the first place.
  ParticleContainer<ItoMergeParticle>& particles = a_particles;

  // A cell-based merge needs the leaves cell-sorted for the CSR ranges. Own that here rather than pushing
  // it onto every caller -- the AMR-collective nn_pair_tree/nn_pair_onecell/nn_pair_hash merges do not use
  // this path and need no sort.
  particles.organizeParticlesByCell();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const int ppc = (lvl < a_particlesPerCell.size()) ? a_particlesPerCell[lvl] : a_particlesPerCell.back();

    const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit  = dbl.dataIterator();
    const int                nbox = dit.size();

    const Real dx = m_amr->getDx()[lvl];

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      ParticleSoA<ItoMergeParticle>& leaf    = particles[lvl][din];
      const EBISBox&                 ebisbox = m_amr->getEBISLayout(m_realm, m_phase)[lvl][din];
      const Box                      cellBox = dbl[din];

      CH_assert(leaf.isSorted()); // organizeParticlesByCell() above sorted every leaf
      CH_assert(leaf.numCells() == static_cast<std::size_t>(cellBox.numPts()));

      // Accumulate the merged particles of every cell here, then swap them into the leaf. The per-cell
      // scratch is a small SoA container reused across cells (extract -> merge in place -> accumulate).
      ParticleSoA<ItoMergeParticle> merged;
      ParticleSoA<ItoMergeParticle> scratch;

      // BoxIterator visits cells in Fortran (x-fastest) order, matching sortByCell's CSR cell index.
      std::size_t cellIndex = 0;
      for (BoxIterator bit(cellBox); bit.ok(); ++bit, ++cellIndex) {
        const IntVect iv = bit();

        leaf.extractCell(cellIndex, scratch);

        if (scratch.size() == 0) {
          continue;
        }

        if (ebisbox.isRegular(iv)) {
          a_merger(scratch, CellInfo(iv, dx), ppc);
        }
        else if (ebisbox.isIrregular(iv)) {
          // Multi-valued cells are not supported -- use the (single) VoF in the cell.
          const VolIndex  vof           = ebisbox.getVoFs(iv).stdVector().front();
          const Real      kappa         = ebisbox.volFrac(vof);
          const RealVect& bndryCentroid = ebisbox.bndryCentroid(vof);
          const RealVect& bndryNormal   = ebisbox.normal(vof);

          a_merger(scratch, CellInfo(iv, dx, kappa, bndryCentroid, bndryNormal), ppc);
        }
        // Covered cells should not contain particles; if they do, they pass through unchanged.

        merged.append(scratch);
      }

      leaf.swap(merged);
    }
  }

  // Restore the by-patch view -- the merged leaves are no longer cell-sorted, and callers expect the
  // container patch-organized on return.
  particles.organizeParticlesByPatch();
}

void
ItoSolver::mergeEqualWeightKD(ParticleContainer<ItoMergeParticle>& a_particles, const Vector<int>& a_particlesPerCell)
{
  CH_TIME("ItoSolver::mergeEqualWeightKD");
  if (m_verbosity > 5) {
    pout() << m_name + "::mergeEqualWeightKD" << endl;
  }

  // Recursively partition particles into at most a_ppc equal-weight KD leaves, then reduce each leaf
  // to one particle at the weighted-centroid position. Requires particle weights >= 1 to split.
  using PType = ParticleManagement::MergeParticle<Real>; // payload = energy

  // Pack the reduced particle's columns into the per-cell AoS intermediate.
  const std::function<PType(const ParticleSoA<ItoMergeParticle>&, std::size_t)> gather =
    [](const ParticleSoA<ItoMergeParticle>& a, const std::size_t i) -> PType {
    PType p;

    p.weight   = a.weight(i);
    p.payload  = a.template get<&ItoMergeParticle::energy>(i);
    p.position = a.position(i);

    return p;
  };

  // Propagate energy to both daughters when the median particle is split across a KD boundary.
  const ParticleManagement::BinaryParticleReconcile<PType> reconcile =
    [](PType& p1, PType& p2, const PType& p0) -> void {
    p1.payload = p0.payload;
    p2.payload = p0.payload;
  };

  // Reduce each KD leaf to a single weighted-centroid particle.
  const std::function<void(ParticleSoA<ItoMergeParticle>&, const PType*, const PType*, const CellInfo&)> scatterLeaf =
    [](ParticleSoA<ItoMergeParticle>& a, const PType* first, const PType* last, const CellInfo&) -> void {
    Real     w = 0.0;
    Real     e = 0.0;
    RealVect x = RealVect::Zero;

    for (const PType* p = first; p != last; ++p) {
      w += p->weight;
      x += p->weight * p->position;
      e += p->weight * p->payload;
    }

    x *= 1.0 / w;
    e *= 1.0 / w;

    ItoMergeParticle payload;
    payload.energy = static_cast<ParticleReal>(e);
    a.append(x, w, payload);
  };

  const ParticleManagement::ParticleMerger<ItoMergeParticle>
    merger = ParticleManagement::makeEqualWeightKDMerger<PType, &PType::weight, &PType::position>(gather,
                                                                                                  reconcile,
                                                                                                  scatterLeaf);

  this->applyCellMerger(a_particles, a_particlesPerCell, merger);
}

void
ItoSolver::mergeReinitialize(ParticleContainer<ItoMergeParticle>& a_particles, const Vector<int>& a_particlesPerCell)
{
  CH_TIME("ItoSolver::mergeReinitialize");
  if (m_verbosity > 5) {
    pout() << m_name + "::mergeReinitialize" << endl;
  }

  // Sums the total number of physical particles in the cell, then redistributes them into at most a_ppc
  // computational particles with as-equal-as-possible integer weights, placed at random cell positions.
  // All output particles carry the same weight-averaged energy. Requires integer-valued weights.

  // Sum physical particle count and compute weight-averaged energy across the cell.
  const std::function<std::pair<long long, Real>(const ParticleSoA<ItoMergeParticle>&)> aggregate =
    [](const ParticleSoA<ItoMergeParticle>& a) -> std::pair<long long, Real> {
    long long numPhysical = 0LL;
    Real      E           = 0.0;

    for (std::size_t i = 0; i < a.size(); i++) {
      const Real w = a.weight(i);

      numPhysical = numPhysical + (long long)w;
      E           = E + w * a.template get<&ItoMergeParticle::energy>(i);
    }

    const Real avgE = (numPhysical > 0) ? E / static_cast<double>(numPhysical) : 0.0;

    return {numPhysical, avgE};
  };

  // Emit one new particle with the given weight and average energy at the drawn position.
  const std::function<void(ParticleSoA<ItoMergeParticle>&, const RealVect&, long long, const Real&)> emit =
    [](ParticleSoA<ItoMergeParticle>& a, const RealVect& x, const long long wt, const Real& avgE) -> void {
    ItoMergeParticle payload;

    payload.energy = static_cast<ParticleReal>(avgE);

    a.append(x, static_cast<double>(wt), payload);
  };

  const ParticleManagement::ParticleMerger<ItoMergeParticle>
    merger = ParticleManagement::makeReinitializeMerger<Real, ItoMergeParticle>(aggregate, emit, [this]() noexcept {
      return m_amr->getProbLo();
    });

  this->applyCellMerger(a_particles, a_particlesPerCell, merger);
}

void
ItoSolver::mergeReinitializeBVH(ParticleContainer<ItoMergeParticle>& a_particles, const Vector<int>& a_particlesPerCell)
{
  CH_TIME("ItoSolver::mergeReinitializeBVH");
  if (m_verbosity > 5) {
    pout() << m_name + "::mergeReinitializeBVH" << endl;
  }

  // Same KD partition as equal_weight_kd, but leaf positions are reinitialized: cut-cells use the
  // weighted centroid (to stay inside the EB), full cells draw a random point in the leaf bounding box.
  // Requires particle weights >= 1 to split.
  using PType = ParticleManagement::MergeParticle<Real>; // payload = energy

  // Pack the reduced particle's columns into the per-cell AoS intermediate.
  const std::function<PType(const ParticleSoA<ItoMergeParticle>&, std::size_t)> gather =
    [](const ParticleSoA<ItoMergeParticle>& a, const std::size_t i) -> PType {
    PType p;

    p.weight   = a.weight(i);
    p.payload  = a.template get<&ItoMergeParticle::energy>(i);
    p.position = a.position(i);

    return p;
  };

  // Propagate energy to both daughters when the median particle is split across a KD boundary.
  const ParticleManagement::BinaryParticleReconcile<PType> reconcile =
    [](PType& p1, PType& p2, const PType& p0) -> void {
    p1.payload = p0.payload;
    p2.payload = p0.payload;
  };

  // Cut-cells: weighted-centroid position to avoid placing particles outside the EB.
  // Full cells: random point in the leaf bounding box to reinitialize spatial distribution.
  const std::function<void(ParticleSoA<ItoMergeParticle>&, const PType*, const PType*, const CellInfo&)> scatterLeaf =
    [](ParticleSoA<ItoMergeParticle>& a, const PType* first, const PType* last, const CellInfo& cellInfo) -> void {
    Real w = 0.0;
    Real e = 0.0;

    if (cellInfo.getVolFrac() < 1.0) {
      RealVect x = RealVect::Zero;

      for (const PType* p = first; p != last; ++p) {
        w += p->weight;
        x += p->weight * p->position;
        e += p->weight * p->payload;
      }

      x *= 1.0 / w;
      e *= 1.0 / w;

      ItoMergeParticle payload;
      payload.energy = static_cast<ParticleReal>(e);
      a.append(x, w, payload);
    }
    else {
      RealVect xMin = +std::numeric_limits<Real>::max() * RealVect::Unit;
      RealVect xMax = -std::numeric_limits<Real>::max() * RealVect::Unit;

      for (const PType* p = first; p != last; ++p) {
        w += p->weight;
        e += p->weight * p->payload;

        const RealVect x = p->position;

        for (int dir = 0; dir < SpaceDim; dir++) {
          xMin[dir] = std::min(xMin[dir], x[dir]);
          xMax[dir] = std::max(xMax[dir], x[dir]);
        }
      }

      RealVect x;

      for (int dir = 0; dir < SpaceDim; dir++) {
        x[dir] = xMin[dir] + Random::getUniformReal01() * (xMax[dir] - xMin[dir]);
      }

      e *= 1.0 / w;

      ItoMergeParticle payload;
      payload.energy = static_cast<ParticleReal>(e);
      a.append(x, w, payload);
    }
  };

  const ParticleManagement::ParticleMerger<ItoMergeParticle>
    merger = ParticleManagement::makeEqualWeightKDMerger<PType, &PType::weight, &PType::position>(gather,
                                                                                                  reconcile,
                                                                                                  scatterLeaf);

  this->applyCellMerger(a_particles, a_particlesPerCell, merger);
}

void
ItoSolver::mergeNnSfc(ParticleContainer<ItoMergeParticle>& a_particles, const Vector<int>& a_particlesPerCell)
{
  CH_TIME("ItoSolver::mergeNnSfc");
  if (m_verbosity > 5) {
    pout() << m_name + "::mergeNnSfc" << endl;
  }

  // Sorts particles along a Hilbert curve, then merges adjacent pairs until the count reaches a_ppc.
  // Produces better spatial locality than the KD methods and does not require integer weights.
  using PType = ParticleManagement::MergeParticle<Real>; // payload = energy

  // Pack the reduced particle's columns into the per-cell AoS intermediate.
  const std::function<PType(const ParticleSoA<ItoMergeParticle>&, std::size_t)> gather =
    [](const ParticleSoA<ItoMergeParticle>& a, const std::size_t i) -> PType {
    PType p;

    p.weight   = a.weight(i);
    p.payload  = a.template get<&ItoMergeParticle::energy>(i);
    p.position = a.position(i);

    return p;
  };

  // Weighted-average position and energy when two particles are merged into one.
  const std::function<void(PType&, const PType&)> combine = [](PType& a, const PType& b) -> void {
    const Real wa  = a.weight;
    const Real wb  = b.weight;
    const Real w   = wa + wb;
    const Real inv = (w > 0.0) ? 1.0 / w : 0.0;

    a.position = (wa * a.position + wb * b.position) * inv;
    a.payload  = (wa * a.payload + wb * b.payload) * inv;
    a.weight   = w;
  };

  // Unpack the merged intermediate back into the reduced particle and append it to the SoA.
  const std::function<void(ParticleSoA<ItoMergeParticle>&, const PType&)> scatter = [](ParticleSoA<ItoMergeParticle>& a,
                                                                                       const PType& p) -> void {
    ItoMergeParticle payload;
    payload.energy = static_cast<ParticleReal>(p.payload);
    a.append(p.position, p.weight, payload);
  };

  const ParticleManagement::ParticleMerger<ItoMergeParticle>
    merger = ParticleManagement::makeSfcNearestNeighborMerger<PType, &PType::weight, &PType::position>(gather,
                                                                                                       combine,
                                                                                                       scatter);

  this->applyCellMerger(a_particles, a_particlesPerCell, merger);
}

void
ItoSolver::extractIntoMergeContainer(const WhichContainer a_container, ParticleContainer<ItoMergeParticle>& a_merge)
{
  CH_TIME("ItoSolver::extractIntoMergeContainer");
  if (m_verbosity > 5) {
    pout() << m_name + "::extractIntoMergeContainer" << endl;
  }

  // Extract into a minimal merge particle (position + weight owned by the leaf, energy the only payload).
  // ItoParticle carries far more than the merge needs; its mobility/diffusion/velocity/... are
  // recomputed from the field right after merging. Copy the few columns the merge uses into a small
  // ItoMergeParticle container on the SAME realm, then kill the ItoParticles -- so the repeated
  // cross-rank ghost exchange and container churn move only a couple of reals per particle instead of
  // the full ItoParticle.
  ParticleContainer<ItoParticle>& particles = this->getParticles(a_container);

  // Ghosts first. The loop below copies every slot with no isGhost() test and append() marks the result
  // valid, so a ghost copy would silently become an owned particle and double-count its weight. Nothing
  // ghost-fills the bulk container today, but this function is called from preRegrid() -- before the
  // ParticleContainer::preRegrid() that used to guarantee the invariant -- so enforce it here rather
  // than depend on a distant caller. preRegrid() and remap() clear defensively for the same reason.
  particles.clearGhostParticles();

  m_amr->allocate(a_merge, m_realm);

  // Extract into the SAME patch (positions unchanged => ownership unchanged => no remap needed).
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit  = dbl.dataIterator();
    const int                nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex&          din     = dit[mybox];
      ParticleSoA<ItoParticle>& itoLeaf = particles[lvl][din];

      ParticleSoA<ItoMergeParticle>& mergeLeaf = a_merge[lvl][din];

      // Exact, not geometric: growTo() doubles from 16, so a leaf built by repeated append() carries up
      // to 2x overshoot, while the final count is known here and reserve() sets the capacity exactly.
      mergeLeaf.reserve(itoLeaf.size());

      for (std::size_t i = 0; i < itoLeaf.size(); i++) {
        ItoMergeParticle p;
        p.energy = itoLeaf.template get<&ItoParticle::energy>(i);

        mergeLeaf.append(itoLeaf.position(i), itoLeaf.weight(i), p);

        // append() sets id/rank to s_invalidID/-1. Carry the id: it is preserved by every transport path
        // in ParticleContainer and is written to H5Part output, so losing it here would rename every
        // particle that merely passes through. rankID is left alone -- nothing in-tree reads it on an
        // ItoParticle, and it is already -1 after any merge.
        mergeLeaf.particleID(mergeLeaf.size() - 1) = itoLeaf.particleID(i);
      }

      // Release this patch's ItoParticle arena as soon as its reduced copy exists, rather than clearing
      // the whole container at the end. clearParticles() deliberately keeps each leaf's capacity, so the
      // 149 B arenas would otherwise stay resident alongside the 53 B copies for the whole extract --
      // and, on the regrid path, through the entire load-balancing phase. Overlap becomes one patch
      // instead of one level. ParticleContainer::gatherToPool() drains its source the same way.
      // Safe under the omp loop: each thread owns its own din, and no kernel here holds a raw column<>()
      // pointer, which shrinkToFit() would invalidate.
      itoLeaf.clear();
      itoLeaf.shrinkToFit();
    }
  }

  // The per-leaf clears above already emptied every patch; this resets the container-level bookkeeping.
  particles.clearParticles();
}

void
ItoSolver::splitFromMergeContainer(ParticleContainer<ItoMergeParticle>& a_merge, const int a_numParticlesPerCellThresh)
{
  CH_TIME("ItoSolver::splitFromMergeContainer");
  if (m_verbosity > 5) {
    pout() << m_name + "::splitFromMergeContainer" << endl;
  }

  // The merge phase only drains over-full cells. A cell BELOW the target is brought up to it by
  // splitting. The merge left the small container patch-organized; cell-sort so the per-cell CSR
  // ranges are valid.
  //
  // Greedy heaviest-split (same rule as nn_sfc): repeatedly halve the heaviest particle into two
  // co-located daughters -- floor/ceil so integer weights stay integer and both daughters keep weight
  // >= 1 -- until the cell reaches the target or no particle can be split (heaviest weight < 2).
  // Daughters sit at the parent's position and inherit its energy, so the spatial distribution and the
  // weighted centroid are preserved. No EB check is needed: a daughter is co-located with a parent that
  // already lives at a valid position.
  a_merge.organizeParticlesByCell();

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din  = dit[mybox];
      const Box        cbox = dbl[din];

      ParticleSoA<ItoMergeParticle>& leaf = a_merge[lvl][din];

      ParticleSoA<ItoMergeParticle> split;
      ParticleSoA<ItoMergeParticle> scratch;

      // BoxIterator visits cells in the same (x-fastest) order as sortByCell's CSR cell index.
      std::size_t cellIndex = 0;
      for (BoxIterator bit(cbox); bit.ok(); ++bit, ++cellIndex) {
        leaf.extractCell(cellIndex, scratch);

        if (scratch.size() == 0) {
          continue;
        }

        while (scratch.size() < static_cast<std::size_t>(a_numParticlesPerCellThresh)) {
          std::size_t hi = 0;
          for (std::size_t i = 1; i < scratch.size(); i++) {
            if (scratch.weight(i) > scratch.weight(hi)) {
              hi = i;
            }
          }

          if (scratch.weight(hi) < 2.0) {
            break;
          }

          const Real     halfWeight = std::floor(scratch.weight(hi) * 0.5);
          const RealVect x          = scratch.position(hi);

          ItoMergeParticle daughter;
          daughter.energy = scratch.template get<&ItoMergeParticle::energy>(hi);

          scratch.weight(hi) -= halfWeight;
          scratch.append(x, halfWeight, daughter);
        }

        split.append(scratch);
      }

      leaf.swap(split);
    }
  }
}

void
ItoSolver::rebuildFromMergeContainer(const WhichContainer a_container, ParticleContainer<ItoMergeParticle>& a_merge)
{
  CH_TIME("ItoSolver::rebuildFromMergeContainer");
  if (m_verbosity > 5) {
    pout() << m_name + "::rebuildFromMergeContainer" << endl;
  }

  // Rebuild ItoParticles from the reduced particles, into the same patch (no remap).
  // velocity/mobility/diffusion/... default to zero and are recomputed from the field before use.
  // particleID is carried across explicitly -- append() would set it to s_invalidID, and for the
  // methods that do not renumber (and for 'none' on the regrid path) the incoming id is the particle's
  // real identity, which reaches H5Part output. mergeLite() has already invalidated the ids of any
  // method that mints its own, so this copy is unconditional and always right.
  ParticleContainer<ItoParticle>& particles = this->getParticles(a_container);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit  = dbl.dataIterator();
    const int                nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex&               din       = dit[mybox];
      ParticleSoA<ItoMergeParticle>& mergeLeaf = a_merge[lvl][din];

      ParticleSoA<ItoParticle>& itoLeaf = particles[lvl][din];

      itoLeaf.reserve(itoLeaf.size() + mergeLeaf.size());

      for (std::size_t i = 0; i < mergeLeaf.size(); i++) {
        ItoParticle p;
        p.energy = mergeLeaf.template get<&ItoMergeParticle::energy>(i);

        itoLeaf.append(mergeLeaf.position(i), mergeLeaf.weight(i), p);

        itoLeaf.particleID(itoLeaf.size() - 1) = mergeLeaf.particleID(i);
      }

      // Hand the reduced leaf's arena back now rather than at end of scope: the ItoParticle copy of
      // this patch already exists, so from here the two would coexist for the rest of the level loop.
      mergeLeaf.clear();
      mergeLeaf.shrinkToFit();
    }
  }

  // Leave the bulk container patch-organized, matching the cell-based merge path's contract.
  particles.organizeParticlesByPatch();
}

template <ItoSolver::NnPairSearchBackend Backend>
void
ItoSolver::mergeNnPairSearch(ParticleContainer<ItoMergeParticle>& merge, const Vector<int>& a_particlesPerCell)
{
  CH_TIME("ItoSolver::mergeNnPairSearch");
  if (m_verbosity > 5) {
    pout() << m_name + "::mergeNnPairSearch" << endl;
  }

  // Reduce every over-full cell to a_particlesPerCell superparticles (and refill under-full cells) using
  // the distributed nearest-neighbor pair merge (whole-patch PointCloudBVH or PointCloudHashGrid search,
  // selected by Backend), in five steps:
  //
  //   1. extractIntoMergeContainer(): copy the ItoParticles into a minimal ItoMergeParticle container
  //      (position + weight + energy) on the same realm/patch, then clear the ItoParticles.
  //   2. Build the merge callbacks (gather/combine/scatter), the EB position-validity predicate, and the
  //      rank-namespaced fresh-id allocator.
  //   3. Merge phase: drain over-full cells by looping mergeNearestNeighborsTree()/mergeNearestNeighborsHash()
  //      until a round merges nothing (each round roughly halves a cell's surplus above the target).
  //   4-5. splitAndRebuildFromMergeContainer(): bring under-full cells up to the target by repeatedly
  //      halving the heaviest particle, then rebuild the ItoParticles.

  // The distributed nearest-neighbor merge has no per-level notion -- it uses a SINGLE crowding
  // threshold over the whole AMR hierarchy. Use the coarsest level's value; in practice this vector is
  // uniform wherever nn_pair_tree/nn_pair_hash is used.
  const int a_numParticlesPerCellThresh = a_particlesPerCell[0];

  // 2. Merge callbacks (on the small particle), plus the EB position check and id allocator below.
  auto gather = [](const ParticleSoA<ItoMergeParticle>& a_leaf, const std::size_t a_idx) -> Real {
    return a_leaf.template get<&ItoMergeParticle::energy>(a_idx);
  };

  auto combine = [](const Real a_firstEnergy,
                    const Real a_firstWeight,
                    const Real a_secondEnergy,
                    const Real a_secondWeight) -> Real {
    return (a_firstWeight * a_firstEnergy + a_secondWeight * a_secondEnergy) / (a_firstWeight + a_secondWeight);
  };

  auto scatter = [](ParticleSoA<ItoMergeParticle>& a_leaf, const ParticleManagement::MergeParticle<Real>& a_p) -> void {
    ItoMergeParticle p;
    p.energy = a_p.payload;

    a_leaf.append(a_p.position, a_p.weight, p);
    a_leaf.particleID(a_leaf.size() - 1) = a_p.globalID;
    a_leaf.rankID(a_leaf.size() - 1)     = a_p.ownerRank;
  };

  // Reject a merge whose weighted-centroid position would fall inside the embedded boundary --
  // BaseIF::value() >= 0 is inside the solid (Chombo convention: negative is inside the fluid). A
  // geometry with no analytic implicit function (e.g. DCEL/STL) has nothing to test against, so
  // every position is accepted rather than silently rejecting all merges.
  const RefCountedPtr<BaseIF>& implicitFunction = m_computationalGeometry->getImplicitFunction(m_phase);

  auto isPositionValid = [&implicitFunction](const RealVect& a_pos) -> bool {
    if (implicitFunction.isNull()) {
      return true;
    }

    return implicitFunction->value(a_pos) < 0.0;
  };

  // A negative nn_pair_max_cell_dist means "no cap" (std::nullopt).
  const std::optional<int> maxCellDistance = (m_nnPairMaxCellDistance < 0)
                                               ? std::nullopt
                                               : std::optional<int>(m_nnPairMaxCellDistance);

  // Fresh-id allocator handed to mergeNearestNeighborsTree() below (its IDAllocator argument). Why an
  // allocator at all: the merge fabricates brand-new particles (each merged pair becomes one new
  // particle) that need globally-unique ids, and it does so on every rank independently, mid-round, with
  // no cross-rank communication -- so each rank must be able to mint ids that can never collide with any
  // other rank's, without asking anyone. The scheme: give rank r the disjoint id block
  // [r * rankStride, (r+1) * rankStride) and hand out ids sequentially within it. Different ranks live in
  // different blocks, so their ids are automatically unique with zero coordination. Ids only need to be
  // unique WITHIN one round (they are a per-round labelling, not a persistent identity), so nextID is
  // reset to rankBase at the top of every drain round below. The same allocator also renumbers the
  // existing valid particles at the start of each round, so a particle and its ghosts share one id.
  //
  // rankStride caps a rank at rankStride ids per round; 1e12 is ~1e6x above any realistic per-rank count
  // (valid + merged) and stays int64-safe up to ~9.2e6 ranks (numRanks * rankStride < INT64_MAX).
  constexpr ParticleID rankStride = 1000000000000LL;
  const ParticleID     rankBase   = static_cast<ParticleID>(procID()) * rankStride;

  ParticleID nextID = rankBase;

  auto allocateID = [&nextID, rankBase]() -> ParticleID {
    // Ids must stay inside this rank's own [rankBase, rankBase + rankStride) block, or they would collide
    // with the next rank's. Fail loudly (debug builds) if that invariant is ever violated.
    CH_assert(nextID < rankBase + rankStride);

    return nextID++;
  };

  // 3. Merge phase.
  // A single mergeNearestNeighborsTree() does ONE cross-patch propose/judge/verdict exchange and
  // merges nearest-neighbor pairs; merged particles only become candidates on the NEXT round, so one
  // round roughly halves an over-full cell's surplus above the target. The loop stops early once a round
  // merges nothing (the global valid-particle count stops decreasing), i.e. once every over-full cell has
  // drained to the target. m_nnPairMaxRounds (ItoSolver.nn_pair_max_rounds) caps the number of rounds so
  // the cost stays bounded and predictable instead of running to full convergence.
  unsigned long long nAfter = merge.getNumberOfValidParticlesGlobal();

  for (int round = 0; round < m_nnPairMaxRounds; round++) {
    const unsigned long long nBefore = nAfter;

    // Fresh, globally-unique ids each round, assigned before ghosts exist so fillGhostParticles()
    // copies each ghost's id from its (now-numbered) source. rankID is this rank (all owned valid
    // particles).
    nextID = rankBase;

    if (round > 0) {
      merge.clearGhostParticles();
    }

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

      // Serial (no omp): allocateID() hands out ids from a single shared counter, so parallelising the
      // box loop would race on it and mint duplicate ids.
      for (int mybox = 0; mybox < nbox; mybox++) {
        ParticleSoA<ItoMergeParticle>& leaf = merge[lvl][dit[mybox]];

        for (std::size_t i = 0; i < leaf.size(); i++) {
          leaf.particleID(i) = allocateID();
          leaf.rankID(i)     = procID();
        }
      }
    }

    // Fill a ghost halo of the merge width (= maxCellDistance, or 1 when unbounded) so a candidate up
    // to that many cells across a patch boundary is visible to the search. This width MUST match the
    // registered mask width (registerOperators()) and the a_ghostWidth passed to the round below --
    // every rank uses the same value, which is what defines boundary-exposure and keeps the merge
    // double-merge-safe.
    const int ghostWidth = maxCellDistance.value_or(1);
    merge.fillGhostParticles(m_amr->getParticleGhostMask(m_realm, ghostWidth),
                             m_amr->getParticleGhostMaskCoarToFine(m_realm, ghostWidth),
                             m_amr->getParticleGhostMaskFineToCoar(m_realm, ghostWidth));

    // Merge particles -- this is the most expensive part of the merge algorithm. The only line in
    // this whole function that differs between the two backends.
    if constexpr (Backend == NnPairSearchBackend::Hash) {
      ParticleManagement::mergeNearestNeighborsHash<ItoMergeParticle, Real>(merge,
                                                                            *m_amr,
                                                                            a_numParticlesPerCellThresh,
                                                                            gather,
                                                                            combine,
                                                                            scatter,
                                                                            allocateID,
                                                                            m_nnPairIterate,
                                                                            m_nnPairFallback,
                                                                            maxCellDistance,
                                                                            ghostWidth,
                                                                            isPositionValid);
    }
    else {
      ParticleManagement::mergeNearestNeighborsTree<ItoMergeParticle, Real>(merge,
                                                                            *m_amr,
                                                                            a_numParticlesPerCellThresh,
                                                                            gather,
                                                                            combine,
                                                                            scatter,
                                                                            allocateID,
                                                                            m_nnPairIterate,
                                                                            m_nnPairFallback,
                                                                            maxCellDistance,
                                                                            ghostWidth,
                                                                            isPositionValid);
    }

    // getNumberOfValidParticlesGlobal() is an MPI collective -- every rank calls it and breaks on the
    // same global value, so the loop stays in lockstep across ranks.
    nAfter = merge.getNumberOfValidParticlesGlobal();

    if (nAfter >= nBefore) {
      break;
    }
  }
}

void
ItoSolver::mergeNnPairOneCell(ParticleContainer<ItoMergeParticle>& merge, const Vector<int>& a_particlesPerCell)
{
  CH_TIME("ItoSolver::mergeNnPairOneCell");
  if (m_verbosity > 5) {
    pout() << m_name + "::mergeNnPairOneCell" << endl;
  }

  // Reduce every over-full cell to a_particlesPerCell superparticles (and refill under-full cells) using
  // the distributed nearest-neighbor pair merge (per-cell PointCloudBVH search, Chebyshev cell distance
  // structurally fixed at 1); see the nn_pair_tree merge for the shared five-step outline --
  // this differs only in step 3 (no maxCellDistance, a fixed width-1 ghost halo, and
  // mergeNearestNeighborsOneCell() instead of mergeNearestNeighborsTree()).

  // The distributed nearest-neighbor merge has no per-level notion -- it uses a SINGLE crowding
  // threshold over the whole AMR hierarchy. Use the coarsest level's value; in practice this vector is
  // uniform wherever nn_pair_onecell is used.
  const int a_numParticlesPerCellThresh = a_particlesPerCell[0];

  // 2. Merge callbacks -- identical reconciliation rule to the nn_pair_tree merge.
  auto gather = [](const ParticleSoA<ItoMergeParticle>& a_leaf, const std::size_t a_idx) -> Real {
    return a_leaf.template get<&ItoMergeParticle::energy>(a_idx);
  };

  auto combine = [](const Real a_firstEnergy,
                    const Real a_firstWeight,
                    const Real a_secondEnergy,
                    const Real a_secondWeight) -> Real {
    return (a_firstWeight * a_firstEnergy + a_secondWeight * a_secondEnergy) / (a_firstWeight + a_secondWeight);
  };

  auto scatter = [](ParticleSoA<ItoMergeParticle>& a_leaf, const ParticleManagement::MergeParticle<Real>& a_p) -> void {
    ItoMergeParticle p;
    p.energy = a_p.payload;

    a_leaf.append(a_p.position, a_p.weight, p);
    a_leaf.particleID(a_leaf.size() - 1) = a_p.globalID;
    a_leaf.rankID(a_leaf.size() - 1)     = a_p.ownerRank;
  };

  const RefCountedPtr<BaseIF>& implicitFunction = m_computationalGeometry->getImplicitFunction(m_phase);

  auto isPositionValid = [&implicitFunction](const RealVect& a_pos) -> bool {
    if (implicitFunction.isNull()) {
      return true;
    }

    return implicitFunction->value(a_pos) < 0.0;
  };

  // See the nn_pair_tree merge for the id-allocator rationale.
  constexpr ParticleID rankStride = 1000000000000LL;
  const ParticleID     rankBase   = static_cast<ParticleID>(procID()) * rankStride;

  ParticleID nextID = rankBase;

  auto allocateID = [&nextID, rankBase]() -> ParticleID {
    CH_assert(nextID < rankBase + rankStride);

    return nextID++;
  };

  // 3. Merge phase -- same drain-loop structure as the nn_pair_tree merge, but the ghost width
  // is fixed at 1 (see mergeNearestNeighborsOneCell()'s own docs on why nn_pair_max_cell_dist does not
  // apply here) and mergeNearestNeighborsOneCell() replaces mergeNearestNeighborsTree().
  unsigned long long nAfter = merge.getNumberOfValidParticlesGlobal();

  for (int round = 0; round < m_nnPairMaxRounds; round++) {
    const unsigned long long nBefore = nAfter;

    nextID = rankBase;

    if (round > 0) {
      merge.clearGhostParticles();
    }

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

      // Serial (no omp): allocateID() hands out ids from a single shared counter, so parallelising the
      // box loop would race on it and mint duplicate ids.
      for (int mybox = 0; mybox < nbox; mybox++) {
        ParticleSoA<ItoMergeParticle>& leaf = merge[lvl][dit[mybox]];

        for (std::size_t i = 0; i < leaf.size(); i++) {
          leaf.particleID(i) = allocateID();
          leaf.rankID(i)     = procID();
        }
      }
    }

    merge.fillGhostParticles(m_amr->getParticleGhostMask(m_realm, 1),
                             m_amr->getParticleGhostMaskCoarToFine(m_realm, 1),
                             m_amr->getParticleGhostMaskFineToCoar(m_realm, 1));

    // Merge particles -- this is the most expensive part of the merge algorithm.
    ParticleManagement::mergeNearestNeighborsOneCell<ItoMergeParticle, Real>(merge,
                                                                             *m_amr,
                                                                             a_numParticlesPerCellThresh,
                                                                             gather,
                                                                             combine,
                                                                             scatter,
                                                                             allocateID,
                                                                             m_nnPairIterate,
                                                                             m_nnPairFallback,
                                                                             isPositionValid);

    nAfter = merge.getNumberOfValidParticlesGlobal();

    if (nAfter >= nBefore) {
      break;
    }
  }
}

void
ItoSolver::mergeKDSkinNn(ParticleContainer<ItoMergeParticle>& merge, const Vector<int>& a_particlesPerCell)
{
  CH_TIME("ItoSolver::mergeKDSkinNn");
  if (m_verbosity > 5) {
    pout() << m_name + "::mergeKDSkinNn" << endl;
  }

  // Reduce every over-full cell to a_particlesPerCell superparticles (and refill under-full cells)
  // with a whole-patch kd-tree merge whose boundary tier is the nearest-neighbor pair merge rather
  // than kd_carve's arbitration. The particles are split across TWO containers so the two tiers
  // cannot interfere: what one tier has already merged must not be merged again by the other, but
  // must still occupy room in the shared per-cell target. In six steps:
  //
  //   1. extractIntoMergeContainer(): copy the ItoParticles into a minimal ItoMergeParticle
  //      container, then clear the ItoParticles. A second, empty container of the same type takes
  //      the interior results.
  //   2. Build the merge callbacks. Two combine callbacks are needed: the kd tier collapses a whole
  //      leaf at once (N-ary) while the nn tier merges two particles at a time (pairwise).
  //   3. Interior tier: mergeKDInterior() commits every leaf holding no ghost into the interior
  //      container and leaves the contested particles -- the skin -- behind in the merge container.
  //   4. Budget: count the interior superparticles per cell so the skin tier drains each cell to
  //      what is LEFT of its target rather than to the whole of it.
  //   5. Skin tier: drain the skin with mergeNearestNeighborsOneCell(), whose Chebyshev-1 search
  //      radius matches the width-1 halo and the kd tier's one-cell leaf bound exactly.
  //   6. Fold the interior container back in and splitAndRebuildFromMergeContainer().

  // Like the other whole-container merges, this one has no per-level notion -- a SINGLE crowding
  // threshold over the whole AMR hierarchy. Use the coarsest level's value.
  const int numParticlesPerCellThresh = a_particlesPerCell[0];

  // 1. `merge` becomes the skin container once the interior tier has taken its share.

  ParticleContainer<ItoMergeParticle> interior;
  m_amr->allocate(interior, m_realm);

  // 2. Merge callbacks.
  auto gather = [](const ParticleSoA<ItoMergeParticle>& a_leaf, const std::size_t a_idx) -> Real {
    return a_leaf.template get<&ItoMergeParticle::energy>(a_idx);
  };

  // N-ary: the kd tier collapses an entire leaf in one go.
  auto kdCombine = [](const Real* a_energies, const Real* a_weights, const std::size_t a_count) -> Real {
    Real weightedSum = 0.0;
    Real totalWeight = 0.0;

    for (std::size_t i = 0; i < a_count; i++) {
      weightedSum += a_weights[i] * a_energies[i];
      totalWeight += a_weights[i];
    }

    return weightedSum / totalWeight;
  };

  // Pairwise: the nn tier merges exactly two particles. Same weighted average as above, so the two
  // tiers reconcile energy identically.
  auto nnCombine = [](const Real a_firstEnergy,
                      const Real a_firstWeight,
                      const Real a_secondEnergy,
                      const Real a_secondWeight) -> Real {
    return (a_firstWeight * a_firstEnergy + a_secondWeight * a_secondEnergy) / (a_firstWeight + a_secondWeight);
  };

  auto kdScatter = [](ParticleSoA<ItoMergeParticle>&                 a_leaf,
                      const ParticleManagement::MergeParticle<Real>& a_p) -> void {
    ItoMergeParticle p;
    p.energy = a_p.payload;

    a_leaf.append(a_p.position, a_p.weight, p);
    a_leaf.particleID(a_leaf.size() - 1) = a_p.globalID;
    a_leaf.rankID(a_leaf.size() - 1)     = a_p.ownerRank;
  };

  auto nnScatter = [](ParticleSoA<ItoMergeParticle>&                 a_leaf,
                      const ParticleManagement::MergeParticle<Real>& a_p) -> void {
    ItoMergeParticle p;
    p.energy = a_p.payload;

    a_leaf.append(a_p.position, a_p.weight, p);
    a_leaf.particleID(a_leaf.size() - 1) = a_p.globalID;
    a_leaf.rankID(a_leaf.size() - 1)     = a_p.ownerRank;
  };

  const RefCountedPtr<BaseIF>& implicitFunction = m_computationalGeometry->getImplicitFunction(m_phase);

  auto isPositionValid = [&implicitFunction](const RealVect& a_pos) -> bool {
    if (implicitFunction.isNull()) {
      return true;
    }

    return implicitFunction->value(a_pos) < 0.0;
  };

  // See the nn_pair_tree merge for the id-allocator rationale.
  constexpr ParticleID rankStride = 1000000000000LL;
  const ParticleID     rankBase   = static_cast<ParticleID>(procID()) * rankStride;

  ParticleID nextID = rankBase;

  auto allocateID = [&nextID, rankBase]() -> ParticleID {
    CH_assert(nextID < rankBase + rankStride);

    return nextID++;
  };

  // Assign fresh ids to every valid particle. Done before any ghost exists so fillGhostParticles()
  // copies each ghost's id from its (now-numbered) source and a particle and its ghosts agree.
  auto renumber = [&]() -> void {
    nextID = rankBase;

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

      // Serial (no omp): allocateID() hands out ids from a single shared counter, so parallelising
      // the box loop would race on it and mint duplicate ids.
      for (int mybox = 0; mybox < nbox; mybox++) {
        ParticleSoA<ItoMergeParticle>& leaf = merge[lvl][dit[mybox]];

        for (std::size_t i = 0; i < leaf.size(); i++) {
          leaf.particleID(i) = allocateID();
          leaf.rankID(i)     = procID();
        }
      }
    }
  };

  auto fillWidthOneGhosts = [&](ParticleContainer<ItoMergeParticle>& a_c) -> void {
    a_c.fillGhostParticles(m_amr->getParticleGhostMask(m_realm, 1),
                           m_amr->getParticleGhostMaskCoarToFine(m_realm, 1),
                           m_amr->getParticleGhostMaskFineToCoar(m_realm, 1));
  };

  // 3. Interior tier -- one local pass, no communication beyond the ghost fill it reads.
  renumber();
  fillWidthOneGhosts(merge);

  ParticleManagement::mergeKDInterior<ItoMergeParticle, Real>(merge,
                                                              interior,
                                                              m_kdMergeCellHistogram,
                                                              m_kdMergeLeafQuota,
                                                              *m_amr,
                                                              numParticlesPerCellThresh,
                                                              m_kdSplitWeightLeafDx,
                                                              gather,
                                                              kdCombine,
                                                              kdScatter,
                                                              allocateID,
                                                              isPositionValid);

  merge.clearGhostParticles();

  // 4. Per-cell budget for the skin tier. The interior superparticles are invisible to the nn tier
  // (they live in the other container) but still occupy their cells, so each one reserves a slot.
  // Ghosts are counted too, and deliberately: the nn tier consults the occupancy of cells owned by
  // a NEIGHBOURING patch when it judges a cross-patch pair, and a width-1 fill is exactly the reach
  // over which it can do so. Without this the budget would be right locally and wrong at the seam.
  ParticleManagement::NNCellBudget cellBudget(numParticlesPerCellThresh);

  fillWidthOneGhosts(interior);

  // One reservation per DISTINCT particle. A single interior super-particle is resident in its own
  // patch and simultaneously a ghost in every neighbouring patch this rank also owns, so booking it
  // once per occurrence would reserve its cell several times over and starve the skin tier. The nn
  // tier collapses the same way when it builds its live per-cell count (see gatherMergeParticles()),
  // and the two counts have to agree or the budget is measured against the wrong denominator.
  std::unordered_set<ParticleID> reservedIDs;

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

    // Serial (no omp): reservedIDs is a shared set that dedupes a particle seen from several patches,
    // and cellBudget.reserveOne() accumulates into shared per-cell counts. Both would race.
    for (int mybox = 0; mybox < nbox; mybox++) {
      const ParticleSoA<ItoMergeParticle>& leaf = interior[lvl][dit[mybox]];

      for (std::size_t i = 0; i < leaf.size(); i++) {
        if (!reservedIDs.insert(leaf.particleID(i)).second) {
          continue;
        }

        // A ghost must be booked against the level it came FROM, with that level's dx -- exactly
        // how the nn tier buckets it (see gatherMergeParticles()). Booking a coarse or fine ghost
        // under this patch's level would key it to a cell that tier never looks up, silently losing
        // the reservation at precisely the refinement boundaries where it matters most.
        int originLevel = lvl;

        if (leaf.isGhost(i)) {
          switch (leaf.ghost(i)) {
          case GhostType::Coarse: {
            originLevel = lvl - 1;

            break;
          }
          case GhostType::Fine: {
            originLevel = lvl + 1;

            break;
          }
          default: {
            originLevel = lvl;

            break;
          }
          }
        }

        const RealVect dx = m_amr->getDx()[originLevel] * RealVect::Unit;

        cellBudget.reserveOne(originLevel, leaf.position(i), m_amr->getProbLo(), dx);
      }
    }
  }

  interior.clearGhostParticles();

  // 5. Skin tier -- same drain loop as the nn_pair_onecell merge, against the budget above.
  unsigned long long nAfter = merge.getNumberOfValidParticlesGlobal();

  for (int round = 0; round < m_nnPairMaxRounds; round++) {
    const unsigned long long nBefore = nAfter;

    if (round > 0) {
      merge.clearGhostParticles();
    }

    renumber();
    fillWidthOneGhosts(merge);

    ParticleManagement::mergeNearestNeighborsOneCell<ItoMergeParticle, Real>(merge,
                                                                             *m_amr,
                                                                             cellBudget,
                                                                             gather,
                                                                             nnCombine,
                                                                             nnScatter,
                                                                             allocateID,
                                                                             m_nnPairIterate,
                                                                             m_nnPairFallback,
                                                                             isPositionValid);

    nAfter = merge.getNumberOfValidParticlesGlobal();

    if (nAfter >= nBefore) {
      break;
    }
  }

  merge.clearGhostParticles();

  // 6. Fold the interior results back in. Both containers share a layout and neither tier moved a
  // particle between patches, so this is a per-patch concatenation with no exchange. Ids may collide
  // between the two containers -- renumber() restarts this rank's block every round and the interior
  // ids were minted before that -- which is harmless: ids are a per-round labelling consumed by the
  // merge tiers, both of which are finished, and nothing downstream identifies a particle by id.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      merge[lvl][din].append(interior[lvl][din]);
    }
  }
}

void
ItoSolver::mergeKDImpl(ParticleContainer<ItoMergeParticle>& merge,
                       const Vector<int>&                   a_particlesPerCell,
                       const bool                           a_enableBoundaryCarve)
{
  CH_TIME("ItoSolver::mergeKDImpl");

  // Reduce every over-full cell to a_particlesPerCell superparticles (and refill under-full cells)
  // using the whole-patch kd-tree merge (ParticleManagement::mergeKDCarve or mergeKDPatch).
  // Unlike the nn_pair merges, there is no drain-loop round here -- both
  // are a single, fixed, non-iterative pass -- in four steps:
  //
  //   1. extractIntoMergeContainer(): copy the ItoParticles into a minimal ItoMergeParticle
  //      container (position + weight + energy) on the same realm/patch, then clear the
  //      ItoParticles.
  //   2. Build the merge callbacks (gather/N-ary combine/scatter), the EB position-validity
  //      predicate, and the rank-namespaced fresh-id allocator.
  //   3. Assign every particle a fresh id/rank, then run the selected merge once: mergeKDCarve()
  //      after filling the width-1 same-level/coarse-to-fine/fine-to-coarse particle ghost masks,
  //      or mergeKDPatch() with no ghost fill at all.
  //   4. splitAndRebuildFromMergeContainer(): bring under-full cells up to the target by repeatedly
  //      halving the heaviest particle, then rebuild the ItoParticles.

  // The whole-patch kd-tree merges have no per-level notion -- they use a SINGLE crowding threshold
  // over the whole AMR hierarchy. Use the coarsest level's value; in practice this vector is uniform
  // wherever kd_carve/kd_patch are used.
  const int a_numParticlesPerCellThresh = a_particlesPerCell[0];

  // 2. Merge callbacks -- identical reconciliation rule to the nn_pair variants, just N-ary.
  auto gather = [](const ParticleSoA<ItoMergeParticle>& a_leaf, const std::size_t a_idx) -> Real {
    return a_leaf.template get<&ItoMergeParticle::energy>(a_idx);
  };

  auto combine = [](const Real* a_energies, const Real* a_weights, const std::size_t a_count) -> Real {
    Real weightedSum = 0.0;
    Real totalWeight = 0.0;

    for (std::size_t i = 0; i < a_count; i++) {
      weightedSum += a_weights[i] * a_energies[i];
      totalWeight += a_weights[i];
    }

    return weightedSum / totalWeight;
  };

  auto scatter = [](ParticleSoA<ItoMergeParticle>& a_leaf, const ParticleManagement::MergeParticle<Real>& a_p) -> void {
    ItoMergeParticle p;
    p.energy = a_p.payload;

    a_leaf.append(a_p.position, a_p.weight, p);
    a_leaf.particleID(a_leaf.size() - 1) = a_p.globalID;
    a_leaf.rankID(a_leaf.size() - 1)     = a_p.ownerRank;
  };

  const RefCountedPtr<BaseIF>& implicitFunction = m_computationalGeometry->getImplicitFunction(m_phase);

  auto isPositionValid = [&implicitFunction](const RealVect& a_pos) -> bool {
    if (implicitFunction.isNull()) {
      return true;
    }

    return implicitFunction->value(a_pos) < 0.0;
  };

  // See the nn_pair_tree merge for the id-allocator rationale.
  constexpr ParticleID rankStride = 1000000000000LL;
  const ParticleID     rankBase   = static_cast<ParticleID>(procID()) * rankStride;

  ParticleID nextID = rankBase;

  auto allocateID = [&nextID, rankBase]() -> ParticleID {
    CH_assert(nextID < rankBase + rankStride);

    return nextID++;
  };

  // 3. Fresh ids, ghost fill, single merge pass.
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

    // Serial (no omp): allocateID() hands out ids from a single shared counter, so parallelising
    // the box loop would race on it and mint duplicate ids.
    for (int mybox = 0; mybox < nbox; mybox++) {
      ParticleSoA<ItoMergeParticle>& leaf = merge[lvl][dit[mybox]];

      for (std::size_t i = 0; i < leaf.size(); i++) {
        leaf.particleID(i) = allocateID();
        leaf.rankID(i)     = procID();
      }
    }
  }

  // Merge particles -- this is the most expensive part of the merge algorithm. Only the carve tier
  // needs a ghost halo; mergeKDPatch() requires its absence, since a ghost it merged would be
  // merged again by its true owner.
  if (a_enableBoundaryCarve) {
    merge.fillGhostParticles(m_amr->getParticleGhostMask(m_realm, 1),
                             m_amr->getParticleGhostMaskCoarToFine(m_realm, 1),
                             m_amr->getParticleGhostMaskFineToCoar(m_realm, 1));

    ParticleManagement::mergeKDCarve<ItoMergeParticle, Real>(merge,
                                                             m_kdMergeCellHistogram,
                                                             m_kdMergeLeafQuota,
                                                             *m_amr,
                                                             a_numParticlesPerCellThresh,
                                                             m_kdSplitWeightLeafDx,
                                                             gather,
                                                             combine,
                                                             scatter,
                                                             allocateID,
                                                             isPositionValid);
  }
  else {
    ParticleManagement::mergeKDPatch<ItoMergeParticle, Real>(merge,
                                                             m_kdMergeCellHistogram,
                                                             m_kdMergeLeafQuota,
                                                             *m_amr,
                                                             a_numParticlesPerCellThresh,
                                                             m_kdSplitWeightLeafDx,
                                                             gather,
                                                             combine,
                                                             scatter,
                                                             allocateID,
                                                             isPositionValid);
  }
}

void
ItoSolver::clear(const WhichContainer a_container)
{
  CH_TIME("ItoSolver::clear(string)");
  if (m_verbosity > 5) {
    pout() << m_name + "::clear(string)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(a_container);

  this->clear(particles);
}

void
ItoSolver::clear(ParticleContainer<ItoParticle>& a_particles) const
{
  CH_TIME("ItoSolver::clear(ParticleContainer)");
  if (m_verbosity > 5) {
    pout() << m_name + "::clear(ParticleContainer)" << endl;
  }

  a_particles.clearParticles();
}

#include <CD_NamespaceFooter.H>
