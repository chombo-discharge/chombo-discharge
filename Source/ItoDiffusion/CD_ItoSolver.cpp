/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoSolver.cpp
  @brief  Implementation of CD_ItoSolver.H
  @author Robert Marskar
*/

// Std includes
#include <chrono>

// Chombo includes
#include <CH_Timer.H>
#include <ParmParse.H>
#include <ParticleIO.H>

// Our includes
#include <CD_SimpleItoParticle.H>
#include <CD_NonCommParticle.H>
#include <CD_ItoSolver.H>
#include <CD_Random.H>
#include <CD_DataOps.H>
#include <CD_ParallelOps.H>
#include <CD_ParticleOps.H>
#include <CD_ParticleManagement.H>
#include <CD_BoxLoops.H>
#include <CD_Random.H>
#include <CD_NamespaceHeader.H>

constexpr int ItoSolver::m_comp;
constexpr int ItoSolver::m_nComp;

ItoSolver::ItoSolver()
{
  CH_TIME("ItoSolver::ItoSolver");

  // Default settings
  m_verbosity            = -1;
  m_name                 = "ItoSolver";
  m_className            = "ItoSolver";
  m_realm                = Realm::primal;
  m_phase                = phase::gas;
  m_haloBuffer           = 1;
  m_coarseFineDeposition = CoarseFineDeposition::Halo;
  m_deposition           = DepositionType::CIC;
  m_plotDeposition       = DepositionType::CIC;
  m_checkpointing        = WhichCheckpoint::Particles;
  m_mobilityInterp       = WhichMobilityInterpolation::Direct;

  // Default is to not merge particles
  m_particleMerger = [](List<ItoParticle>& a_particles, const CellInfo& a_cellInfo, const int a_ppc) {

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

const std::string
ItoSolver::getRealm() const
{
  CH_TIME("ItoSolver::getRealm");

  return m_realm;
}

void
ItoSolver::setRealm(const std::string a_realm)
{
  CH_TIME("ItoSolver::setRealm");

  m_realm = a_realm;
}

void
ItoSolver::setParticleMerger(const ParticleManagement::ParticleMerger<ItoParticle>& a_particleMerger) noexcept
{
  CH_TIME("ItoSolver::setParticleMerger");

  m_particleMerger = a_particleMerger;
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
  CH_TIME("McPhoto::parsePlotVariables");
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
      m_plotParticlesSource = true;
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
  if (str == "none") {
    m_particleMerger = [](List<ItoParticle>& a_particles, const CellInfo& a_cellInfo, const int a_ppc) {
    };
  }
  else if (str == "equal_weight_kd") {
    m_particleMerger = [this](List<ItoParticle>& a_particles, const CellInfo& a_cellInfo, const int a_ppc) {
      this->makeSuperparticlesEqualWeightKD(a_particles, a_cellInfo, a_ppc);
    };
  }
  else if (str == "reinitialize") {
    m_particleMerger = [this](List<ItoParticle>& a_particles, const CellInfo& a_cellInfo, const int a_ppc) {
      this->reinitializeParticles(a_particles, a_cellInfo, a_ppc);
    };
  }
  else if (str == "external") {
    // Do nothing, because the user will set the merger algorithm through setParticleMerger
  }
  else {
    MayDay::Abort("ItoSolver::parseParticleMerger - unknown particle merging algorithm requested");
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

    // Register mask for CIC deposition.
    m_amr->registerMask(s_particle_halo, m_haloBuffer, m_realm);
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

  // TLDR: This function will fetch the initial particles from the species and deposit them on the mesh. In most cases the various MPI ranks
  //       will have drawn a different set of initial particles (the only sane way to do it) and so those particles are put directly in
  //       the 'bulk' particle container. After that we remove the particles that fell inside the EB and deposit the particles on the mesh.

  ParticleContainer<ItoParticle>& bulkParticles = m_particleContainers.at(WhichContainer::Bulk);
  bulkParticles.clearParticles();
  bulkParticles.addParticles(m_species->getInitialParticles());
  bulkParticles.remap();

  constexpr Real tolerance = 0.0;

  // Add particles, remove the ones that are inside the EB, and then deposit
  this->removeCoveredParticles(bulkParticles, EBRepresentation::ImplicitFunction, tolerance);
  this->depositParticles<ItoParticle, &ItoParticle::weight>(m_phi, bulkParticles, m_deposition, m_coarseFineDeposition);
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

    a_loads[din.intCode()] = particles[a_level][din].numItems();
  }

  ParallelOps::vectorSum(a_loads);
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
ItoSolver::intersectParticles(const EBIntersection                    a_ebIntersection,
                              const bool                              a_deleteParticles,
                              const std::function<void(ItoParticle&)> a_nonDeletionModifier)
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
ItoSolver::intersectParticles(const WhichContainer                    a_particles,
                              const WhichContainer                    a_ebParticles,
                              const WhichContainer                    a_domainParticles,
                              const EBIntersection                    a_ebIntersection,
                              const bool                              a_deleteParticles,
                              const std::function<void(ItoParticle&)> a_nonDeletionModifier)
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
ItoSolver::intersectParticles(ParticleContainer<ItoParticle>&         a_particles,
                              ParticleContainer<ItoParticle>&         a_ebParticles,
                              ParticleContainer<ItoParticle>&         a_domainParticles,
                              const EBIntersection                    a_ebIntersection,
                              const bool                              a_deleteParticles,
                              const std::function<void(ItoParticle&)> a_nonDeletionModifier)
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
    m_amr->intersectParticlesRaycastIF(a_particles,
                                       a_ebParticles,
                                       a_domainParticles,
                                       m_phase,
                                       tolerance,
                                       a_deleteParticles,
                                       a_nonDeletionModifier);

    break;
  }
  case EBIntersection::Bisection: {
    m_amr->intersectParticlesBisectIF(a_particles,
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

  // Regrid particle containers.
  for (auto& container : m_particleContainers) {
    ParticleContainer<ItoParticle>& particles = container.second;

    m_amr->remapToNewGrids(particles, a_lmin, a_newFinestLevel);
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

  m_particleContainers.emplace(WhichContainer::Bulk, ParticleContainer<ItoParticle>());
  m_particleContainers.emplace(WhichContainer::EB, ParticleContainer<ItoParticle>());
  m_particleContainers.emplace(WhichContainer::Domain, ParticleContainer<ItoParticle>());
  m_particleContainers.emplace(WhichContainer::Source, ParticleContainer<ItoParticle>());
  m_particleContainers.emplace(WhichContainer::Covered, ParticleContainer<ItoParticle>());
  m_particleContainers.emplace(WhichContainer::Scratch, ParticleContainer<ItoParticle>());

  for (auto& container : m_particleContainers) {
    m_amr->allocate(container.second, m_realm);
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

  // TLDR: This routine writes particles directly to the HDF5 file -- since the ItoParticle is quite large (112 bytes, ish), we put the most
  //       important ItoParticle fields into the SimpleItoParticle and checkpoint that particle type instead. This is what happens below.

  // I call this _particlesP to distinguish it from the "fluid" checkpointing method.
  const std::string str = m_name + "_particlesP";

  // Set up a particle container with SimpleItoParticle -- this is ItoParticle's low-memory cousin.
  ParticleContainer<SimpleItoParticle> lowMemoryParticles;
  m_amr->allocate(lowMemoryParticles, m_realm);

  // Handle to the particles that will be checkpointed.
  const ParticleContainer<ItoParticle>& myParticles = this->getParticles(WhichContainer::Bulk);

  // Make ItoParticle into SimpleItoParticle. This saves a ton of disk space.
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_level];
  const DataIterator&      dit = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    // Handle to list of particles low-memory particles in the current grid patch.
    List<SimpleItoParticle>& simpleParticles = (lowMemoryParticles[a_level])[din].listItems();

    // Make the ItoParticle into a SimpleItoParticle for checkpointing purposes -- we only need weight, position, and energy.
    simpleParticles.clear();

    for (ListIterator<ItoParticle> lit(myParticles[a_level][din].listItems()); lit.ok(); ++lit) {
      const ItoParticle& p = lit();

      simpleParticles.append(SimpleItoParticle(p.weight(), p.position(), p.energy()));
    }
  }

  // Finally, write the particles to HDF5.
  writeParticlesToHDF(a_handle, lowMemoryParticles[a_level], str);
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

  // TLDR: This routine checkpoints the particle data using the number of particles in a grid cell. When the simulation is restarted, we read the
  //       number of particles per cell from the HDF5 file and re-initialize the particles. However, this function does NOT currently store the energy
  //       (or other parameters of interest) in the HDF5 file. Only the number of particles is available. I don't expect this function to be widely used
  //       by anyone.

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

    // Add the patch particles into cell-sorted particles.
    BinFab<ItoParticle> cellSortedParticles(cellBox, dx, probLo);
    cellSortedParticles.addItems(particles[a_level][din].listItems());

    // Kernel - go through the patch and get the number of particles per cell.
    auto kernel = [&](const IntVect& iv) -> void {
      const List<ItoParticle>& cellParticles = cellSortedParticles(iv, m_comp);

      // Go through the particles in the current grid cell and set the number of particles.
      particleNumbersFAB(iv, m_comp) = 0.0;
      for (ListIterator<ItoParticle> lit(cellParticles); lit.ok(); ++lit) {
        const ItoParticle& p = lit();

        particleNumbersFAB(iv, m_comp) += p.weight();
      }
    };

    // Execute kernel.
    BoxLoops::loop(cellBox, kernel);
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

  // TLDR: This function is the one that reads SimpleItoParticles from the checkpoint file and instantiates full ItoParticle's from that. Recalling
  //       writeCheckpointLevelParticles we only stored the weight, position, and energy of the particles. Here we read that information back in.

  // This is the particle container that we will fill.
  ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);

  CH_assert(m_checkpointing == WhichCheckpoint::Particles);
  CH_assert(!particles.isOrganizedByCell());

  const std::string str = m_name + "_particlesP";

  // Allocate storage for holding simple particle data and then read the particles into simpleParticles[a_level]. The other
  // grid levels are not touched.
  Vector<RefCountedPtr<ParticleData<SimpleItoParticle>>> simpleParticles;
  m_amr->allocate(simpleParticles, m_realm);

  readParticlesFromHDF(a_handle, *simpleParticles[a_level], str);

  // Go through the particles we read from the file and make them into true ItoParticles.
  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_level];
  const DataIterator&      dit = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    List<ItoParticle>&             itoParticles       = particles[a_level][din].listItems();
    const List<SimpleItoParticle>& simpleItoParticles = (*simpleParticles[a_level])[din].listItems();

    for (ListIterator<SimpleItoParticle> lit(simpleItoParticles); lit.ok(); ++lit) {
      const SimpleItoParticle& simpleParticle = lit();
      const ItoParticle        itoParticle    = ItoParticle(simpleParticle.weight(),
                                                  simpleParticle.position(),
                                                  RealVect::Zero,
                                                  0.0,
                                                  0.0,
                                                  simpleParticle.energy());

      itoParticles.append(itoParticle);
    }
  }
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

  // particlesPerCell holds the number of particles per cell -- call the other version which instantiates new particles from that.
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
    List<ItoParticle>& myParticles = particles[a_level][din].listItems();
    myParticles.clear();

    // Kernel region for cut-cells.
    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_level])[din];

    // Regular kernel
    auto regularKernel = [&](const IntVect& iv) -> void {
      // Do regular cells -- in these cells we only need to draw a random position somewhere inside the cubic cell. Easy.
      if (ebisbox.isRegular(iv)) {

        // This bit of code will take the number of physical particles and divide them into a_newPPC particles with
        // approximately equal weights. It is possible that one of the particles will have a larger particle weight than the others.
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

          myParticles.add(ItoParticle(1.0 * w, x));
        }
      }
    };

    // Irregular kernel. Do the same for irregular cells. This differs from the regular-cell case only in that the positions
    // are checked against the EB.
    auto irregularKernel = [&](const VolIndex& vof) -> void {
      const IntVect  iv    = vof.gridIndex();
      const RealVect cent  = ebisbox.bndryCentroid(vof);
      const RealVect norm  = ebisbox.normal(vof);
      const RealVect pos   = probLo + dx * (RealVect(iv) + 0.5 * RealVect::Unit);
      const Real     kappa = ebisbox.volFrac(vof);

      const unsigned long long numPhysicalParticles = (unsigned long long)llround(ppc(iv));

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
        // approximately equal weights. It is possible that one of the particles will have a larger particle weight than the others.
        const std::vector<long long> weights = ParticleManagement::partitionParticleWeights(llround(ppc(iv)),
                                                                                            (long long)a_newPPC);

        for (const auto& w : weights) {
          const RealVect x = Random::randomPosition(pos, minLo, minHi, cent, norm, dx, kappa);

          myParticles.add(ItoParticle(1.0 * w, x));
        }
      }
    };

    // Run the kernels.
    BoxLoops::loop(cellBox, regularKernel);
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
                         const std::string     a_outputRealm,
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
    this->depositParticlesNGP<ItoParticle, &ItoParticle::weight>(scratch,
                                                                 m_particleContainers.at(WhichContainer::Bulk),
                                                                 a_level);

    m_amr->copyData(a_output, scratch, a_level, a_outputRealm, m_realm, Interval(a_comp, a_comp), Interval(0, 0));

    a_comp++;
  }
  if (m_plotParticlesEB) {
    this->depositParticlesNGP<ItoParticle, &ItoParticle::weight>(scratch,
                                                                 m_particleContainers.at(WhichContainer::EB),
                                                                 a_level);

    m_amr->copyData(a_output, scratch, a_level, a_outputRealm, m_realm, Interval(a_comp, a_comp), Interval(0, 0));

    a_comp++;
  }
  if (m_plotParticlesDomain) {
    this->depositParticlesNGP<ItoParticle, &ItoParticle::weight>(scratch,
                                                                 m_particleContainers.at(WhichContainer::Domain),
                                                                 a_level);

    m_amr->copyData(a_output, scratch, a_level, a_outputRealm, m_realm, Interval(a_comp, a_comp), Interval(0, 0));

    a_comp++;
  }
  if (m_plotParticlesSource) {
    this->depositParticlesNGP<ItoParticle, &ItoParticle::weight>(scratch,
                                                                 m_particleContainers.at(WhichContainer::Source),
                                                                 a_level);

    m_amr->copyData(a_output, scratch, a_level, a_outputRealm, m_realm, Interval(a_comp, a_comp), Interval(0, 0));

    a_comp++;
  }
  if (m_plotEnergyDensity) {
    this->depositParticlesNGP<ItoParticle, &ItoParticle::totalEnergy>(scratch,
                                                                      m_particleContainers.at(WhichContainer::Bulk),
                                                                      a_level);

    m_amr->copyData(a_output, scratch, a_level, a_outputRealm, m_realm, Interval(a_comp, a_comp), Interval(0, 0));

    a_comp++;
  }
  if (m_plotAverageEnergy) {
    LevelData<EBCellFAB> weight;
    m_amr->allocate(weight, m_realm, m_phase, a_level, 1);
    this->depositParticlesNGP<ItoParticle, &ItoParticle::weight>(weight,
                                                                 m_particleContainers.at(WhichContainer::Bulk),
                                                                 a_level);

    this->depositParticlesNGP<ItoParticle, &ItoParticle::totalEnergy>(scratch,
                                                                      m_particleContainers.at(WhichContainer::Bulk),
                                                                      a_level);

    // Set scratch = totalEnergy/totalWeight
    DataOps::divideFallback(scratch, weight, 0.0);

    m_amr->copyData(a_output, scratch, a_level, a_outputRealm, m_realm, Interval(a_comp, a_comp), Interval(0, 0));

    a_comp++;
  }
  CH_STOP(t2);
}

void
ItoSolver::writeData(LevelData<EBCellFAB>& a_output,
                     int&                  a_comp,
                     const EBAMRCellData&  a_data,
                     const std::string     a_outputRealm,
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

  DataOps::setCoveredValue(scratch, 0.0);

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
    this->depositParticles<ItoParticle, &ItoParticle::conductivity>(a_phi,
                                                                    a_particles,
                                                                    a_deposition,
                                                                    a_coarseFineDeposition);
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

  this->depositParticles<ItoParticle, &ItoParticle::diffusivity>(a_phi,
                                                                 a_particles,
                                                                 a_deposition,
                                                                 a_coarseFineDeposition);
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

  this->depositParticles<ItoParticle, &ItoParticle::totalEnergy>(a_phi,
                                                                 a_particles,
                                                                 a_deposition,
                                                                 a_coarseFineDeposition);
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
  this->depositParticles<ItoParticle, &ItoParticle::conductivity>(a_phi,
                                                                  a_particles,
                                                                  m_deposition,
                                                                  m_coarseFineDeposition);
  this->depositParticles<ItoParticle, &ItoParticle::weight>(weight, a_particles, m_deposition, m_coarseFineDeposition);

  // Make averageMobility = weight*mu/weight. If there is no weight then set the value to zero.
  constexpr Real zero = 0.0;

  DataOps::divideFallback(a_phi, weight, zero);
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
  this->depositParticles<ItoParticle, &ItoParticle::diffusivity>(a_phi,
                                                                 a_particles,
                                                                 m_deposition,
                                                                 m_coarseFineDeposition);
  this->depositParticles<ItoParticle, &ItoParticle::weight>(weight, a_particles, m_deposition, m_coarseFineDeposition);

  // Make average diffusion coefficient = weight*D/weight. If there is no weight then set the value to zero.
  constexpr Real zero = 0.0;

  DataOps::divideFallback(a_phi, weight, zero);
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
  this->depositParticles<ItoParticle, &ItoParticle::totalEnergy>(a_phi,
                                                                 a_particles,
                                                                 m_deposition,
                                                                 m_coarseFineDeposition);
  this->depositParticles<ItoParticle, &ItoParticle::weight>(weight, a_particles, m_deposition, m_coarseFineDeposition);

  // Make average energy = weight*energy/weight. If there is no weight then set the value to zero.
  constexpr Real zero = 0.0;

  DataOps::divideFallback(a_phi, weight, zero);
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

  this->depositParticles<ItoParticle, &ItoParticle::weight>(m_phi,
                                                            m_particleContainers.at(a_container),
                                                            m_deposition,
                                                            m_coarseFineDeposition);
}

void
ItoSolver::redistributeAMR(EBAMRCellData& a_phi) const
{
  CH_TIME("ItoSolver::redistributeAMR");
  if (m_verbosity > 5) {
    pout() << m_name + "::redistributeAMR" << endl;
  }

  // TLDR: When we entered this routine we had a_phi = m_i/dV but we actually want to have phi = m_i/(kappa*dV) so as to have
  //       meaningful densities. Thus, we can either run with a_phi just as it is, in which case it must be interpreted as an extended
  //       state into the EB. That is perfectly fine. But we can also use O(1) accurate redistribution in order to make the scheme
  //       completely conservative, if that is important.
  //
  //       If we use redistribution then we compute a hybrid update phiH = kappa*phi = a_phi in each cell. But we are then "missing"
  //       a mass kappa*phi - kappa*phiH = a_phi(1 - kappa). This mass can be smooshed into the neighboring grid cells. The code
  //       below does even more than that -- it can compute an update phiH = kappa*phi + (1-kappa)*phiNC where phiNC is a non-conservative
  //       type of update. In this case the mass loss is just like for fluid models: dM = kappa*(1-kappa)(phiC - phiNC). But

  if (m_useRedistribution) {
    this->depositNonConservative(m_depositionNC, a_phi);    // Compute m_depositionNC = sum(kappa*Wc)/sum(kappa)
    this->depositHybrid(a_phi, m_massDiff, m_depositionNC); // Compute hybrid deposition, including mass differnce

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
    const IrregAmrStencil<NonConservativeDivergenceStencil>& stencils = m_amr->getNonConservativeDivergenceStencils(
      m_realm,
      m_phase);
    stencils.apply(a_depositionNC, a_depositionKappaC);
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
ItoSolver::preRegrid(const int a_lbase, const int a_oldFinestLevel)
{
  CH_TIME("ItoSolver::preRegrid");
  if (m_verbosity > 5) {
    pout() << m_name + "::preRegrid" << endl;
  }

  CH_assert(a_lbase >= 0);

  for (auto& container : m_particleContainers) {
    ParticleContainer<ItoParticle>& particles = container.second;

    particles.preRegrid(a_lbase);
  }

  m_phi.clear();
  m_mobilityFunction.clear();
  m_velocityFunction.clear();
  m_diffusionFunction.clear();
  m_depositionNC.clear();
  m_massDiff.clear();
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
ItoSolver::setVelocityFunction(const RealVect a_velocity)
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

  particles.setValue<&ItoParticle::mobility>(a_mobility);
}

void
ItoSolver::setParticleDiffusion(const Real a_diffCo)
{
  CH_TIME("ItoSolver::setParticleDiffusion");
  if (m_verbosity > 5) {
    pout() << m_name + "::setParticleDiffusion" << endl;
  }

  ParticleContainer<ItoParticle>& particles = this->getParticles(WhichContainer::Bulk);

  particles.setValue<&ItoParticle::diffusion>(a_diffCo);
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

  EBAMRParticleMesh& particleMesh = m_amr->getParticleMesh(m_realm, m_phase);

  if (m_isMobile) {
    const EBCellFAB& velo_func = (*m_velocityFunction[a_lvl])[a_dit];
    const RealVect   dx        = m_amr->getDx()[a_lvl] * RealVect::Unit;
    const RealVect   origin    = m_amr->getProbLo();
    const Box        box       = m_amr->getGrids(m_realm)[a_lvl][a_dit];

    List<ItoParticle>& particleList = particles[a_lvl][a_dit].listItems();

    // This interpolates the velocity function on to the particle velocities
    const EBParticleMesh& meshInterp = particleMesh.getEBParticleMesh(a_lvl, a_dit);

    meshInterp.interpolate<ItoParticle, &ItoParticle::velocity>(particleList,
                                                                velo_func,
                                                                m_deposition,
                                                                m_forceIrregInterpolationNGP);

    // Go through the particles and set their velocities to velo_func*mobility
    for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
      ItoParticle& p = lit();
      p.velocity() *= p.mobility();
    }
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
      DataOps::vectorLength(velocityMagnitude, m_velocityFunction);

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

  // TLDR: This will compute the particle mobility by interpolating a scalar mobility field (stored on the mesh) to the particle positions.

  ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);

  EBAMRParticleMesh& particleMesh = m_amr->getParticleMesh(m_realm, m_phase);

  const EBCellFAB& mobilityFunction = (*m_mobilityFunction[a_lvl])[a_dit];
  const RealVect   dx               = m_amr->getDx()[a_lvl] * RealVect::Unit;
  const RealVect   probLo           = m_amr->getProbLo();
  const Box        box              = m_amr->getGrids(m_realm)[a_lvl][a_dit];

  List<ItoParticle>& particleList = particles[a_lvl][a_dit].listItems();

  // Interpolate onto the mobility field
  const EBParticleMesh& meshInterp = particleMesh.getEBParticleMesh(a_lvl, a_dit);

  meshInterp.interpolate<ItoParticle, &ItoParticle::mobility>(particleList,
                                                              mobilityFunction,
                                                              m_deposition,
                                                              m_forceIrregInterpolationNGP);
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

  // TLDR: This function computes the particle mobilities by interpolating mu*V to the particle position and then setting
  //       the mobility as mu = [mu*V(Xp)]/V(Xp).

  const EBAMRParticleMesh& particleMesh     = m_amr->getParticleMesh(m_realm, m_phase);
  const EBCellFAB&         mobilityFunction = (*m_mobilityFunction[a_lvl])[a_dit];
  const RealVect           dx               = m_amr->getDx()[a_lvl] * RealVect::Unit;
  const RealVect           probLo           = m_amr->getProbLo();
  const Box                box              = m_amr->getGrids(m_realm)[a_lvl][a_dit];
  const EBParticleMesh&    meshInterp       = particleMesh.getEBParticleMesh(a_lvl, a_dit);

  ParticleContainer<ItoParticle>& particles    = m_particleContainers.at(WhichContainer::Bulk);
  List<ItoParticle>&              particleList = particles[a_lvl][a_dit].listItems();

  // Compute mu*|V| on the mesh
  EBCellFAB muV;
  muV.clone(a_velocityMagnitude);
  muV *= mobilityFunction;

  // First, interpolate |V| to the particle position, it will be stored on m_tmp.
  meshInterp.interpolate<ItoParticle, &ItoParticle::tmpReal>(particleList,
                                                             a_velocityMagnitude,
                                                             m_deposition,
                                                             m_forceIrregInterpolationNGP);

  meshInterp.interpolate<ItoParticle, &ItoParticle::mobility>(particleList,
                                                              muV,
                                                              m_deposition,
                                                              m_forceIrregInterpolationNGP);

  // We now have ItoParticle::tmp = |V(Xp)| and ItoParticle::mobility = |mu*V|(Xp). Now let the mobility be mu(Xp) = |mu*V|(Xp)/|V|(Xp)
  for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
    ItoParticle& p = lit();

    p.mobility() *= 1. / p.tmpReal();
  }
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
ItoSolver::updateMobilities(const int a_level, const DataIndex a_dit)
{
  CH_TIME("ItoSolver::updateMobilities(int, DataIndex)");
  if (m_verbosity > 5) {
    pout() << m_name + "::updateMobilities(int, DataIndex)" << endl;
  }

  ParticleContainer<ItoParticle>& particles = m_particleContainers.at(WhichContainer::Bulk);

  if (m_isMobile) {
    List<ItoParticle>& particleList = particles[a_level][a_dit].listItems();

    for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
      ItoParticle& p = lit();

      p.mobility() = m_species->mobility(p.energy());
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
    ParticleContainer<ItoParticle>& particles    = m_particleContainers.at(WhichContainer::Bulk);
    List<ItoParticle>&              particleList = particles[a_lvl][a_dit].listItems();

    EBAMRParticleMesh& particleMesh = m_amr->getParticleMesh(m_realm, m_phase);

    // Create the particle interpolator.
    const EBCellFAB& Dcoef = (*m_diffusionFunction[a_lvl])[a_dit];
    const RealVect   dx    = m_amr->getDx()[a_lvl] * RealVect::Unit;

    const EBParticleMesh& meshInterp = particleMesh.getEBParticleMesh(a_lvl, a_dit);

    // Interpolator
    meshInterp.interpolate<ItoParticle, &ItoParticle::diffusion>(particleList,
                                                                 Dcoef,
                                                                 m_deposition,
                                                                 m_forceIrregInterpolationNGP);
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
ItoSolver::updateDiffusion(const int a_level, const DataIndex a_dit)
{
  CH_TIME("ItoSolver::updateDiffusion(lvl, dit)");
  if (m_verbosity > 5) {
    pout() << m_name + "::updateDiffusion(lvl, dit)" << endl;
  }

  if (m_isDiffusive) {

    // Particles to be updated.
    ParticleContainer<ItoParticle>& particles    = this->getParticles(WhichContainer::Bulk);
    List<ItoParticle>&              particleList = particles[a_level][a_dit].listItems();

    // Update diffusion coefficient.
    for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
      ItoParticle& p = lit();

      p.diffusion() = m_species->diffusion(p.energy());
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
  const ParticleContainer<ItoParticle>& particles    = this->getParticles(WhichContainer::Bulk);
  const List<ItoParticle>&              particleList = particles[a_lvl][a_dit].listItems();

  // Iteratore for our particles.
  ListIterator<ItoParticle> lit(particleList);

  if (m_isMobile && !m_isDiffusive) {

    // Advection but no diffusion - set the time step as dt = dx/vMax where vMax is the largest velocity component.
    for (lit.rewind(); lit.ok(); ++lit) {
      const ItoParticle& p = lit();
      const RealVect&    v = p.velocity();

      // Get the largest velocity component
      constexpr bool doAbs  = true;
      const int      maxDir = v.maxDir(doAbs);
      const Real     vMax   = std::abs(v[maxDir]);

      // Compute this particle's time step as dt = dx/vMax
      const Real thisDt = (vMax > 0.0) ? dx / vMax : std::numeric_limits<Real>::max();

      dt = std::min(dt, thisDt);
    }
  }
  else if (!m_isMobile && m_isDiffusive) {
    // Diffusion but no advection -- set the time step as dt = dx*dx/(2*D)

    for (lit.rewind(); lit.ok(); ++lit) {
      const ItoParticle& p = lit();

      // Get the diffusion coefficient and compute dt = dx*dx/(2*D)
      const Real D      = p.diffusion();
      const Real thisDt = (D > 0.0) ? dx * dx / (2.0 * D) : std::numeric_limits<Real>::max();

      dt = std::min(dt, thisDt);
    }
  }
  else if (m_isMobile && m_isDiffusive) {

    // Both advectino and diffusion. Compute dt = 1/(1/dtA + 1/dtD) where dtA and dtD are as in the code bits above.
    for (lit.rewind(); lit.ok(); ++lit) {
      const ItoParticle& p = lit();
      const RealVect&    v = p.velocity();
      const Real&        D = p.diffusion();

      // Get the largest velocity component.
      constexpr bool doAbs  = true;
      const int      maxDir = v.maxDir(doAbs);
      const Real     vMax   = std::abs(v[maxDir]);

      const Real dtAdvect  = (vMax > 0.0) ? dx / vMax : std::numeric_limits<Real>::max();
      const Real dtDiffuse = (D > 0.0) ? dx * dx / (2.0 * D) : std::numeric_limits<Real>::max();

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

  // TLDR: This routine computes the largest possible time step such that no particles move more than a_maxCellsToMove during
  //       a standard Ito kernel step. This is the AMR version.

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

  // TLDR: This routine computes the largest possible time step such that no particles move more than a_maxCellsToMove during
  //       a standard Ito kernel step. This is the level version.

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

  // TLDR: This routine computes the largest possible time step such that no particles move more than a_maxCellsToMove during
  //       a standard Ito kernel step. This is the patch version.

  Real dt = std::numeric_limits<Real>::max();

  // Grid resolution on this level.
  const Real dx = m_amr->getDx()[a_lvl];

  // Some shortcuts.
  const Real dMax  = a_maxCellsToMove * dx;
  const Real dMax2 = dMax * dMax;
  const Real W0    = m_normalDistributionTruncation;
  const Real W02   = m_normalDistributionTruncation * m_normalDistributionTruncation;

  // These are the particles we will compute the time steps for.
  const ParticleContainer<ItoParticle>& particles    = this->getParticles(WhichContainer::Bulk);
  const List<ItoParticle>&              particleList = particles[a_lvl][a_dit].listItems();

  // Make an iterator for our particles.
  ListIterator<ItoParticle> lit(particleList);

  // Cases handled differently
  if (m_isMobile && !m_isDiffusive) {
    // Advection but no diffusion

    for (lit.rewind(); lit; ++lit) {
      const ItoParticle& p = lit();
      const RealVect&    v = p.velocity();

      // Compute the regular time step -- recall that dMax = a_maxCellsToMove*dx
      const int  maxDir = v.maxDir(true);
      const Real thisDt = dMax / std::abs(v[maxDir]);

      dt = std::min(dt, thisDt);
    }
  }
  else if (!m_isMobile && m_isDiffusive) {
    // Diffusion but no advection

    for (lit.rewind(); lit; ++lit) {
      const ItoParticle& p = lit();
      const Real&        D = p.diffusion();

      // Recall, the diffusion kernel is usually dX = sqrt(2*D*dt)*N where N is a SpaceDim vector of Gaussian numbers. But we only care about
      // not moving more than a specified number of grid cells in any one of the coordinate directions so we have
      // |dX| = sqrt(2*D*dt)*N0 where N0 is the maximum value in the Gaussian distribution (which we have truncated). Solving for dt yields dt = |dX|^2/(2*D*N0^2).
      const Real thisDt = dMax2 / (2.0 * D * SpaceDim * W02);

      dt = std::min(dt, thisDt);
    }
  }
  else if (m_isMobile && m_isDiffusive) {

    // Diffusion AND advection. Much more difficult and requires us to solve a second order equation.
    for (lit.rewind(); lit; ++lit) {
      const ItoParticle& p = lit();
      const RealVect&    v = p.velocity();
      const Real&        D = p.diffusion();

      // Get the maximum velocity component.
      const int  maxDir = v.maxDir(true);
      const Real vMax   = Abs(v[maxDir]);

      Real thisDt = std::numeric_limits<Real>::max();

      // This case is more complicated because we have dX = v*dt + sqrt(2*D*dt)*N0 and we need to solve for dt. This yields a second order equation
      // in the form A*dt^2 + B*dt + C = 0. We just solve for dt.
      if (vMax > 0.0) {
        const Real a = vMax;
        const Real b = W0 * sqrt(2.0 * D);
        const Real c = dMax;

        const Real A = a * a;
        const Real B = -(b * b + 2 * a * c);
        const Real C = c * c;

        thisDt = (-B - sqrt(B * B - 4. * A * C)) / (2. * A);
      }
      else {
        if (D > 0.0) {
          thisDt = dMax2 / (2.0 * D * W02);
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
    const ParticleContainer<ItoParticle>& particles    = m_particleContainers.at(WhichContainer::Bulk);
    const List<ItoParticle>&              particleList = particles[a_lvl][a_dit].listItems();

    for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
      const ItoParticle& p = lit();
      const RealVect&    v = p.velocity();

      // Get maximum velocity component.
      const int  maxDir = v.maxDir(true);
      const Real vMax   = std::abs(v[maxDir]);

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
    const Real dx2 = dx * dx / 2.0;

    // These are the particles we iterate over.
    const ParticleContainer<ItoParticle>& particles    = this->getParticles(WhichContainer::Bulk);
    const List<ItoParticle>&              particleList = particles[a_lvl][a_dit].listItems();

    for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
      const ItoParticle& p = lit();
      const Real&        D = p.diffusion();

      const Real thisDt = (D > 0.0) ? dx2 / D : std::numeric_limits<Real>::max();

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
ItoSolver::makeSuperparticles(const WhichContainer a_container, const int a_particlesPerCell)
{
  CH_TIME("ItoSolver::makeSuperparticles(WhichContainer, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::makeSuperparticles(WhichContainer, int)" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    this->makeSuperparticles(a_container, a_particlesPerCell, lvl);
  }
}

void
ItoSolver::makeSuperparticles(const WhichContainer a_container, const Vector<int> a_particlesPerCell)
{
  CH_TIME("ItoSolver::makeSuperparticles(WhichContainer, Vector<int>)");
  if (m_verbosity > 5) {
    pout() << m_name + "::makeSuperparticles(WhichContainer, Vector<int>)" << endl;
  }

  if (a_particlesPerCell.size() < 1) {
    MayDay::Error("ItoSolver::makeSuperParticles(Container, Vector<int>) -logic bust");
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    int ppc;

    if (lvl < a_particlesPerCell.size()) {
      ppc = a_particlesPerCell[lvl];
    }
    else {
      ppc = a_particlesPerCell.back();
    }

    this->makeSuperparticles(a_container, ppc, lvl);
  }
}

void
ItoSolver::makeSuperparticles(const WhichContainer a_container, const int a_particlesPerCell, const int a_level)
{
  CH_TIME("ItoSolver::makeSuperparticles(WhichContainer, int, int)");
  if (m_verbosity > 5) {
    pout() << m_name + "::makeSuperparticles(WhichContainer, int, int)" << endl;
  }

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_level];
  const DataIterator&      dit = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    this->makeSuperparticles(a_container, a_particlesPerCell, a_level, din);
  }
}

void
ItoSolver::makeSuperparticles(const WhichContainer a_container,
                              const int            a_particlesPerCell,
                              const int            a_level,
                              const DataIndex      a_dit)
{
  CH_TIME("ItoSolver::makeSuperparticles(WhichContainer, int, int, DataIndex)");
  if (m_verbosity > 5) {
    pout() << m_name + "::makeSuperparticles(WhichContainer, int, int, DataIndex)" << endl;
  }

  // This are the particles in the box we're currently looking at.
  ParticleContainer<ItoParticle>& particles     = this->getParticles(a_container);
  BinFab<ItoParticle>&            cellParticles = particles.getCellParticles(a_level, a_dit);

  const Real     dx      = m_amr->getDx()[a_level];
  const EBISBox& ebisbox = m_amr->getEBISLayout(m_realm, m_phase)[a_level][a_dit];

  // Kernel for particle merging in regular cells
  auto regularKernel = [&](const IntVect& iv) -> void {
    if (ebisbox.isRegular(iv)) {
      List<ItoParticle>& particles = cellParticles(iv, m_comp);

      if (particles.length() > 0) {
        m_particleMerger(particles, CellInfo(iv, dx), a_particlesPerCell);
      }
    }
  };

  // Kernel for particle merging in irregular cells
  auto irregularKernel = [&](const VolIndex& vof) -> void {
    const IntVect iv = vof.gridIndex();

    List<ItoParticle>& particles = cellParticles(iv, m_comp);

    if (particles.length() > 0) {
      const Real      kappa         = ebisbox.volFrac(vof);
      const RealVect& bndryCentroid = ebisbox.bndryCentroid(vof);
      const RealVect& bndryNormal   = ebisbox.normal(vof);

      CellInfo cellInfo(iv, dx, kappa, bndryCentroid, bndryNormal);

      m_particleMerger(particles, cellInfo, a_particlesPerCell);
    }
  };

  // Iteration space.
  const Box    cellBox = m_amr->getGrids(m_realm)[a_level][a_dit];
  VoFIterator& vofit   = (*m_amr->getVofIterator(m_realm, m_phase)[a_level])[a_dit];

  // Run kernel
  BoxLoops::loop(cellBox, regularKernel);
  BoxLoops::loop(vofit, irregularKernel);
}

void
ItoSolver::mergeParticles(List<ItoParticle>& a_particles, const CellInfo& a_cellInfo, const int a_ppc)
{
  CH_TIMERS("ItoSolver::mergeParticles");

  m_particleMerger(a_particles, a_cellInfo, a_ppc);
}

void
ItoSolver::makeSuperparticlesEqualWeightKD(List<ItoParticle>& a_particles, const CellInfo& a_cellInfo, const int a_ppc)
{
  CH_TIMERS("ItoSolver::makeSuperparticlesEqualWeightKD");
  CH_TIMER("ItoSolver::makeSuperparticlesEqualWeightKD::populate_list", t1);
  CH_TIMER("ItoSolver::makeSuperparticlesEqualWeightKD::build_kd", t2);
  CH_TIMER("ItoSolver::makeSuperparticlesEqualWeightKD::merge_particles", t3);

  using PType        = NonCommParticle<2, 1>;
  using Node         = KDNode<PType>;
  using ParticleList = KDNode<PType>::ParticleList;

  // 1. Make the input list into a vector of particles with a smaller memory footprint.
  CH_START(t1);
  Real         W = 0.0;
  ParticleList particles;
  for (ListIterator<ItoParticle> lit(a_particles); lit.ok(); ++lit) {
    PType p;

    p.template real<0>() = lit().weight();
    p.template real<1>() = lit().energy();
    p.template vect<0>() = lit().position();

    W += lit().weight();

    particles.emplace_back(p);
  }
  CH_STOP(t1);

  // Particle reconciler for manipulated two particles arising from splitting of one particle.
  auto particleReconcile = [](PType& p1, PType& p2, const PType& p0) -> void {
    p1.template real<1>() = p0.template real<1>();
    p2.template real<1>() = p0.template real<1>();
  };

  // 2. Build KD-tree.
  const std::vector<std::shared_ptr<Node>> leaves = ParticleManagement::
    recursivePartitionAndSplitEqualWeightKD<PType, &PType::template real<0>, &PType::template vect<0>>(
      particles,
      a_ppc,
      particleReconcile);

  // Merge leaves into new particles.
  CH_START(t3);
  a_particles.clear();

  for (const auto& l : leaves) {
    Real     w = 0.0;
    Real     e = 0.0;
    RealVect x = RealVect::Zero;

    for (const auto& p : l->getParticles()) {
      w += p.template real<0>();
      x += p.template real<0>() * p.template vect<0>();
      e += p.template real<0>() * p.template real<1>();
    }

    x *= 1. / w;
    e *= 1. / w;

    a_particles.add(ItoParticle(w, x, RealVect::Zero, 0.0, 0.0, e));
  }
  CH_STOP(t3);
}

void
ItoSolver::reinitializeParticles(List<ItoParticle>& a_particles, const CellInfo& a_cellInfo, const int a_ppc)
{
  CH_TIME("ItoSolver::reinitializeParticles");

  // Compute number of physical particles and partition them into a_ppc computational particles
  // with weights as equal as possible. Average energy becomes true average energy (same for all particles).
  long long numPhysicalParticles = 0LL;
  Real      averageEnergy        = 0.0;
  for (ListIterator<ItoParticle> lit(a_particles); lit.ok(); ++lit) {
    const ItoParticle& p = lit();
    const Real         w = p.weight();
    const Real         e = p.energy();

    numPhysicalParticles += (long long)w;
    averageEnergy += w * e;
  }

  if (numPhysicalParticles > 0) {
    averageEnergy *= 1.0 / numPhysicalParticles;
  }

  const std::vector<long long> weights = ParticleManagement::partitionParticleWeights(numPhysicalParticles,
                                                                                      (long long)a_ppc);

  // randomPosition wants the physical corners.
  const Real     dx            = a_cellInfo.getDx();
  const Real     kappa         = a_cellInfo.getVolFrac();
  const RealVect probLo        = m_amr->getProbLo();
  const RealVect cellPos       = probLo + dx * (a_cellInfo.getGridIndex() + 0.5 * RealVect::Unit);
  const RealVect validLo       = a_cellInfo.getValidLo();
  const RealVect validHi       = a_cellInfo.getValidHi();
  const RealVect bndryCentroid = a_cellInfo.getBndryCentroid();
  const RealVect bndryNormal   = a_cellInfo.getBndryNormal();

  a_particles.clear();

  for (int i = 0; i < weights.size(); i++) {
    const Real     w = weights[i];
    const RealVect x = Random::randomPosition(cellPos, validLo, validHi, bndryCentroid, bndryNormal, dx, kappa);

    a_particles.add(ItoParticle(w, x, RealVect::Zero, 0.0, 0.0, averageEnergy));
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

  this->clear(a_particles.getParticles());
}

void
ItoSolver::clear(AMRParticles<ItoParticle>& a_particles) const
{
  CH_TIME("ItoSolver::clear");
  if (m_verbosity > 5) {
    pout() << m_name + "::clear" << endl;
  }

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    a_particles[lvl]->clear();
  }
}

#include <CD_NamespaceFooter.H>
