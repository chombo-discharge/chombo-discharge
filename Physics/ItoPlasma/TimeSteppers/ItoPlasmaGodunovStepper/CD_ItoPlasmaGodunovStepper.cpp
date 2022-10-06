/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaGodunovStepper.cpp
  @brief  Implementation of CD_ItoPlasmaGodunovStepper.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_ItoPlasmaGodunovStepper.H>
#include <CD_Timer.H>
#include <CD_DataOps.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoPlasma;

ItoPlasmaGodunovStepper::ItoPlasmaGodunovStepper(RefCountedPtr<ItoPlasmaPhysics>& a_physics)
  : ItoPlasmaStepper(a_physics)
{
  m_name    = "ItoPlasmaGodunovStepper";
  m_physics = a_physics;

  ParmParse pp("ItoPlasmaGodunovStepper");
  pp.get("load_ppc", m_loadPerCell);
  pp.get("nwo_reactions", m_useNewReactionAlgorithm);

  m_avg_cfl = 0.0;
}

ItoPlasmaGodunovStepper::~ItoPlasmaGodunovStepper() {}

int
ItoPlasmaGodunovStepper::getNumberOfPlotVariables() const
{
  CH_TIME("ItoPlasmaGodunovStepper::getNumberOfPlotVariables");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaGodunovStepper::getNumberOfPlotVariables" << endl;
  }

  int ncomp = ItoPlasmaStepper::getNumberOfPlotVariables();

  ncomp++; // Add conductivity

  return ncomp;
}

void
ItoPlasmaGodunovStepper::writePlotData(EBAMRCellData&       a_output,
                                       Vector<std::string>& a_plotVariableNames,
                                       int&                 a_icomp) const
{
  CH_TIME("ItoPlasmaGodunovStepper::writeConductivity");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaGodunovStepper::writeConductivity" << endl;
  }

  ItoPlasmaStepper::writePlotData(a_output, a_plotVariableNames, a_icomp);

  // Do conductivity
  this->writeConductivity(a_output, a_icomp);
  a_plotVariableNames.push_back("conductivity");
}

void
ItoPlasmaGodunovStepper::writeConductivity(EBAMRCellData& a_output, int& a_icomp) const
{
  CH_TIME("ItoPlasmaStepper::writeConductivity");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaStepper::writeConductivity" << endl;
  }

  const Interval src_interv(0, 0);
  const Interval dst_interv(a_icomp, a_icomp);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
    if (m_conductivityCell.getRealm() == a_output.getRealm()) {
      m_conductivityCell[lvl]->localCopyTo(src_interv, *a_output[lvl], dst_interv);
    }
    else {
      m_conductivityCell[lvl]->copyTo(src_interv, *a_output[lvl], dst_interv);
    }
  }
  a_icomp += 1;
}

void
ItoPlasmaGodunovStepper::allocate()
{
  CH_TIME("ItoPlasmaGodunovStepper::allocate");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaGodunovStepper::allocate" << endl;
  }

  ItoPlasmaStepper::allocate();

  // Now allocate for the conductivity particles and rho^dagger particles
  const int num_ItoSpecies = m_physics->getNumItoSpecies();

  m_conductivity_particles.resize(num_ItoSpecies);
  m_rho_dagger_particles.resize(num_ItoSpecies);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    const RefCountedPtr<ItoSolver>& solver = solver_it();

    const int idx = solver_it.index();

    m_conductivity_particles[idx] = new ParticleContainer<PointParticle>();
    m_rho_dagger_particles[idx]   = new ParticleContainer<PointParticle>();

    m_amr->allocate(*m_conductivity_particles[idx], m_particleRealm);
    m_amr->allocate(*m_rho_dagger_particles[idx], m_particleRealm);
  }
}

void
ItoPlasmaGodunovStepper::parseOptions()
{
  CH_TIME("ItoPlasmaGodunovStepper::parseOptions");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseOptions" << endl;
  }

  ParmParse   pp(m_name.c_str());
  std::string str;

  pp.get("verbosity", m_verbosity);
  pp.get("ppc", m_particlesPerCell);
  pp.get("regrid_super", m_regridSuperparticles);
  pp.get("algorithm", str);
  pp.get("load_balance", m_loadBalance);
  pp.get("load_index", m_loadBalanceIndex);
  pp.get("min_dt", m_minDt);
  pp.get("max_dt", m_maxDt);
  pp.get("eb_tolerance", m_eb_tolerance);

  // Get algorithm
  if (str == "euler_maruyama") {
    m_algorithm = which_algorithm::euler_maruyama;
  }
  else if (str == "trapezoidal") {
    m_algorithm = which_algorithm::trapezoidal;
  }
  else {
    MayDay::Abort("ItoPlasmaGodunovStepper::parseOptions - unknown algorithm requested");
  }

  // Dt limitation
  pp.get("which_dt", str);
  if (str == "advection") {
    m_whichDt = which_dt::advection;
  }
  else if (str == "diffusion") {
    m_whichDt = which_dt::diffusion;
  }
  else if (str == "AdvectionDiffusion") {
    m_whichDt = which_dt::AdvectionDiffusion;
  }
  else {
    MayDay::Abort("ItoPlasmaGodunovStepper::parseOptions - unknown 'which_dt' requested");
  }

  // Setup runtime storage (requirements change with algorithm)
  this->setupRuntimeStorage();
}

void
ItoPlasmaGodunovStepper::parseRuntimeOptions()
{
  CH_TIME("ItoPlasmaGodunovStepper::parseRuntimeOptions");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRuntimeOptions" << endl;
  }

  ParmParse   pp(m_name.c_str());
  std::string str;

  pp.get("verbosity", m_verbosity);
  pp.get("ppc", m_particlesPerCell);
  pp.get("merge_interval", m_mergeInterval);
  pp.get("relax_factor", m_relaxTimeFactor);
  pp.get("regrid_super", m_regridSuperparticles);
  pp.get("algorithm", str);
  pp.get("load_balance", m_loadBalance);
  pp.get("load_index", m_loadBalanceIndex);
  pp.get("min_dt", m_minDt);
  pp.get("max_dt", m_maxDt);
  pp.get("eb_tolerance", m_eb_tolerance);

  // Get algorithm
  if (str == "euler_maruyama") {
    m_algorithm = which_algorithm::euler_maruyama;
  }
  else if (str == "trapezoidal") {
    m_algorithm = which_algorithm::trapezoidal;
  }
  else {
    MayDay::Abort("ItoPlasmaGodunovStepper::parseOptions - unknown algorithm requested");
  }

  // Dt limitation
  pp.get("which_dt", str);
  if (str == "advection") {
    m_whichDt = which_dt::advection;
  }
  else if (str == "diffusion") {
    m_whichDt = which_dt::diffusion;
  }
  else if (str == "AdvectionDiffusion") {
    m_whichDt = which_dt::AdvectionDiffusion;
  }
  else {
    MayDay::Abort("ItoPlasmaGodunovStepper::parseOptions - unknown 'which_dt' requested");
  }

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
    MayDay::Abort(
      "ItoPlasmaGodunovStepper::parseOptions - unknown box sorting method requested for argument 'BoxSorting'");
  }

  // Setup runtime storage (requirements change with algorithm)
  this->setupRuntimeStorage();

  //
  m_ito->parseRuntimeOptions();
  m_fieldSolver->parseRuntimeOptions();
  m_rte->parseRuntimeOptions();
}

void
ItoPlasmaGodunovStepper::allocateInternals()
{
  CH_TIME("ItoPlasmaGodunovStepper::allocateInternals");
  if (m_verbosity > 5) {
    pout() << m_name + "::allocateInternals" << endl;
  }

  const int num_ItoSpecies = m_physics->getNumItoSpecies();
  const int num_rtSpecies  = m_physics->getNumRtSpecies();

  m_amr->allocate(m_scratch1, m_fluidRealm, m_phase, 1);
  m_amr->allocate(m_scratch2, m_fluidRealm, m_phase, 1);

  // Allocate fluid scratch storage
  m_fscratch1.resize(num_ItoSpecies);
  m_fscratch2.resize(num_ItoSpecies);
  for (int i = 0; i < num_ItoSpecies; i++) {
    m_amr->allocate(m_fscratch1[i], m_fluidRealm, m_phase, 1);
    m_amr->allocate(m_fscratch2[i], m_fluidRealm, m_phase, 1);
  }
}

Real
ItoPlasmaGodunovStepper::advance(const Real a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::advance");
  if (m_verbosity > 5) {
    pout() << m_name + "::advance" << endl;
  }

  // Advance the particles.
  switch (m_algorithm) {
  case which_algorithm::euler_maruyama: {
    this->advanceParticlesEulerMaruyama(a_dt);

    break;
  }
  case which_algorithm::trapezoidal: {
    this->advanceParticlesTrapezoidal(a_dt);

    break;
  }
  default: {
    MayDay::Abort("ItoPlasmaGodunovStepper::advance - logic bust");

    break;
  }
  }

  // Compute current and relaxation time.
  this->computeJ(m_currentDensity, a_dt);
  const Real relaxTime = this->computeRelaxationTime(); // This is for the restricting the next step.

  // Move Photons
  this->advancePhotons(a_dt);

  // If we are using the LEA, we must compute the Ohmic heating term. This must be done
  // BEFORE sorting the particles per cell.
  if (m_physics->getCoupling() == ItoPlasmaPhysics::coupling::LEA) {
    this->computeEdotJSource(a_dt);
  }

  // Sort the particles and Photons per cell so we can call reaction algorithms
  m_ito->sortParticlesByCell(ItoSolver::WhichContainer::Bulk);
  this->sortPhotonsByCell(McPhoto::WhichContainer::Bulk);
  this->sortPhotonsByCell(McPhoto::WhichContainer::Source);

  // Chemistry kernel.
  this->advanceReactionNetwork(a_dt);

  // Make superparticles.
  if ((m_timeStep + 1) % m_mergeInterval == 0 && m_mergeInterval > 0) {
    m_ito->makeSuperparticles(ItoSolver::WhichContainer::Bulk, m_particlesPerCell);
  }

  // Sort particles per patch.
  m_ito->sortParticlesByPatch(ItoSolver::WhichContainer::Bulk);
  this->sortPhotonsByPatch(McPhoto::WhichContainer::Bulk);
  this->sortPhotonsByPatch(McPhoto::WhichContainer::Source);

  // Clear other data holders for now. BC comes later...
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    solver_it()->clear(ItoSolver::WhichContainer::EB);
    solver_it()->clear(ItoSolver::WhichContainer::Domain);
  }

  // Deposit particles
  m_ito->depositParticles();

  // Prepare next step
  this->computeItoVelocities();
  this->computeItoDiffusion();

  return a_dt;
}

Real
ItoPlasmaGodunovStepper::computeDt()
{
  CH_TIME("ItoPlasmaGodunovStepper::computeDt");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaGodunovStepper::computeDt" << endl;
  }

  Real a_dt = std::numeric_limits<Real>::max();

  if (m_whichDt == which_dt::advection) {
    a_dt = m_advectionCFL * m_ito->computeAdvectiveDt();
  }
  else if (m_whichDt == which_dt::diffusion) {
    a_dt = m_diffusionCFL * m_ito->computeDiffusiveDt();
  }
  else if (m_whichDt == which_dt::AdvectionDiffusion) {
    a_dt = m_advectionDiffusionCFL * m_ito->computeDt();
  }

  // Physics-based restriction
  const Real physicsDt = this->computePhysicsDt();
  if (physicsDt < a_dt) {
    a_dt       = physicsDt;
    m_timeCode = TimeCode::Physics;
  }

  if (a_dt < m_minDt) {
    a_dt       = m_minDt;
    m_timeCode = TimeCode::Hardcap;
  }

  if (a_dt > m_maxDt) {
    a_dt       = m_maxDt;
    m_timeCode = TimeCode::Hardcap;
  }

#if 0 // Debug code
  const Real dtCFL = m_ito->computeDt();
  m_avg_cfl += a_dt/dtCFL;
  if(procID() == 0) std::cout << "dt = " << a_dt
			      << "\t relax dt = " << relaxTime
			      << "\t factor = " << a_dt/relaxTime
			      << "\t CFL = " << a_dt/dtCFL
			      << "\t avgCFL = " << m_avg_cfl/(1+m_timeStep)
			      << std::endl;
#endif

  return a_dt;
}

void
ItoPlasmaGodunovStepper::preRegrid(const int a_lmin, const int a_oldFinestLevel)
{
  CH_TIME("ItoPlasmaGodunovStepper::preRegrid");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaGodunovStepper::preRegrid" << endl;
  }

  ItoPlasmaStepper::preRegrid(a_lmin, a_oldFinestLevel);

  // Copy conductivity to scratch storage
  const int ncomp        = 1;
  const int finest_level = m_amr->getFinestLevel();
  m_amr->allocate(m_cache, m_fluidRealm, m_phase, ncomp);
  for (int lvl = 0; lvl <= a_oldFinestLevel; lvl++) {
    m_conductivityCell[lvl]->localCopyTo(*m_cache[lvl]);
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    const int idx = solver_it.index();

    m_conductivity_particles[idx]->preRegrid(a_lmin);
    m_rho_dagger_particles[idx]->preRegrid(a_lmin);
  }
}

void
ItoPlasmaGodunovStepper::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("ItoPlasmaGodunovStepper::regrid");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaGodunovStepper::regrid" << endl;
  }

  // Regrid solvers
  m_ito->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_fieldSolver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_rte->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
  m_sigma->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);

  // Allocate internal memory for ItoPlasmaGodunovStepper now....
  this->allocateInternals();

  // We need to remap/regrid the stored particles as well.
  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    const int idx = solver_it.index();
    m_amr->remapToNewGrids(*m_rho_dagger_particles[idx], a_lmin, a_newFinestLevel);
    m_amr->remapToNewGrids(*m_conductivity_particles[idx], a_lmin, a_newFinestLevel);
  }

  // Recompute the conductivity and space charge densities.
  this->computeRegridConductivity();
  this->computeRegridRho();
  this->setupSemiImplicitPoisson(m_prevDt);

  // Solve the Poisson equation.
  const bool converged = this->solvePoisson();
  if (!converged) {
    MayDay::Abort("ItoPlasmaGodunovStepper::regrid - Poisson solve did not converge after regrid!!!");
  }

  // Regrid superparticles.
  if (m_regridSuperparticles) {
    m_ito->sortParticlesByCell(ItoSolver::WhichContainer::Bulk);
    m_ito->makeSuperparticles(ItoSolver::WhichContainer::Bulk, m_particlesPerCell);
    m_ito->sortParticlesByPatch(ItoSolver::WhichContainer::Bulk);
  }

  // Now let Ihe ito solver deposit its actual particles... In the above it deposit m_rho_dagger_particles.
  m_ito->depositParticles();

  // Recompute new velocities and diffusion coefficients
  this->computeItoVelocities();
  this->computeItoDiffusion();
}

void
ItoPlasmaGodunovStepper::setupRuntimeStorage()
{
  CH_TIME("ItoPlasmaGodunovStepper::setupRuntimeStorage");
  if (m_verbosity > 5) {
    pout() << m_name + "::setupRuntimeStorage" << endl;
  }
#if 0
  switch (m_algorithm) {
  case which_algorithm::euler_maruyama:
    ItoParticle::setNumRuntimeVectors(1);
    break;
  case which_algorithm::trapezoidal:
    ItoParticle::setNumRuntimeVectors(2); // For V^k and the diffusion hop.
    break;
  default:
    MayDay::Abort("ItoPlasmaGodunovStepper::setupRuntimeStorage - logic bust");
  }
#else
  MayDay::Error(
    "ItoPlasmaGodunovStepper::setupRuntimeStorage -- need to figure out how to add more fields to ItoParticle");
#endif
}

void
ItoPlasmaGodunovStepper::setOldPositions()
{
  CH_TIME("ItoPlasmaGodunovStepper::setOldPositions()");
  if (m_verbosity > 5) {
    pout() << m_name + "::setOldPositions()" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>& solver = solver_it();

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_particleRealm)[lvl];
      ParticleData<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

        List<ItoParticle>& particleList = particles[dit()].listItems();

        for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
          ItoParticle& p  = particleList[lit];
          p.oldPosition() = p.position();
        }
      }
    }
  }
}

void
ItoPlasmaGodunovStepper::remapGodunovParticles(Vector<ParticleContainer<PointParticle>*>& a_particles,
                                               const SpeciesSubset                       a_SpeciesSubset)
{
  CH_TIME("ItoPlasmaGodunovStepper::remapGodunovParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::remapGodunovParticles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = species->getChargeNumber() != 0;

    switch (a_SpeciesSubset) {
    case SpeciesSubset::All:
      a_particles[idx]->remap();
      break;
    case SpeciesSubset::AllMobile:
      if (mobile)
        a_particles[idx]->remap();
      break;
    case SpeciesSubset::AllDiffusive:
      if (diffusive)
        a_particles[idx]->remap();
      break;
    case SpeciesSubset::ChargedMobile:
      if (charged && mobile)
        a_particles[idx]->remap();
      break;
    case SpeciesSubset::ChargedDiffusive:
      if (charged && diffusive)
        a_particles[idx]->remap();
      break;
    case SpeciesSubset::AllMobileOrDiffusive:
      if (mobile || diffusive)
        a_particles[idx]->remap();
      break;
    case SpeciesSubset::ChargedAndMobileOrDiffusive:
      if (charged && (mobile || diffusive))
        a_particles[idx]->remap();
      break;
    case SpeciesSubset::Stationary:
      if (!mobile && !diffusive)
        a_particles[idx]->remap();
      break;
    default:
      MayDay::Abort("ItoPlasmaGodunovStepper::remapGodunovParticles - logic bust");
    }
  }
}

void
ItoPlasmaGodunovStepper::deposit_PointParticles(const Vector<ParticleContainer<PointParticle>*>& a_particles,
                                                const SpeciesSubset                             a_SpeciesSubset)
{
  CH_TIME("ItoPlasmaGodunovStepper::deposit_PointParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::deposit_PointParticles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = species->getChargeNumber() != 0;

    switch (a_SpeciesSubset) {
    case SpeciesSubset::All:
      solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      break;
    case SpeciesSubset::AllMobile:
      if (mobile)
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      break;
    case SpeciesSubset::AllDiffusive:
      if (diffusive)
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      break;
    case SpeciesSubset::ChargedMobile:
      if (charged && mobile)
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      break;
    case SpeciesSubset::ChargedDiffusive:
      if (charged && diffusive)
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      break;
    case SpeciesSubset::AllMobileOrDiffusive:
      if (mobile || diffusive)
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      break;
    case SpeciesSubset::ChargedAndMobileOrDiffusive:
      if (charged && (mobile || diffusive))
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      break;
    case SpeciesSubset::Stationary:
      if (!mobile && !diffusive)
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      break;
    default:
      MayDay::Abort("ItoPlasmaGodunovStepper::deposit_PointParticles - logic bust");
    }
  }
}

void
ItoPlasmaGodunovStepper::clearGodunovParticles(const Vector<ParticleContainer<PointParticle>*>& a_particles,
                                               const SpeciesSubset                             a_SpeciesSubset)
{
  CH_TIME("ItoPlasmaGodunovStepper::clearGodunovParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::deposit_clearParticles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = species->getChargeNumber() != 0;

    switch (a_SpeciesSubset) {
    case SpeciesSubset::All:
      a_particles[idx]->clearParticles();
      break;
    case SpeciesSubset::AllMobile:
      if (mobile)
        a_particles[idx]->clearParticles();
      break;
    case SpeciesSubset::AllDiffusive:
      if (diffusive)
        a_particles[idx]->clearParticles();
      break;
    case SpeciesSubset::ChargedMobile:
      if (charged && mobile)
        a_particles[idx]->clearParticles();
      break;
    case SpeciesSubset::ChargedDiffusive:
      if (charged && diffusive)
        a_particles[idx]->clearParticles();
      break;
    case SpeciesSubset::AllMobileOrDiffusive:
      if (mobile || diffusive)
        a_particles[idx]->clearParticles();
      break;
    case SpeciesSubset::ChargedAndMobileOrDiffusive:
      if (charged && (mobile || diffusive))
        a_particles[idx]->clearParticles();
      break;
    case SpeciesSubset::Stationary:
      if (!mobile && !diffusive)
        a_particles[idx]->clearParticles();
      break;
    default:
      MayDay::Abort("ItoPlasmaGodunovStepper::clearGodunovParticles - logic bust");
    }
  }
}

void
ItoPlasmaGodunovStepper::computeAllConductivities(const Vector<ParticleContainer<PointParticle>*>& a_particles)
{
  CH_TIME("ItoPlasmaGodunovStepper::computeAllConductivities");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeAllConductivities" << endl;
  }

  this->compute_cell_conductivity(m_conductivityCell, a_particles);

  // Now do the faces
  this->compute_face_conductivity();
}

void
ItoPlasmaGodunovStepper::compute_cell_conductivity(EBAMRCellData&                                   a_conductivity,
                                                   const Vector<ParticleContainer<PointParticle>*>& a_particles)
{
  CH_TIME("ItoPlasmaGodunovStepper::compute_cell_conductivity(conductivity, PointParticle");
  if (m_verbosity > 5) {
    pout() << m_name + "::compute_cell_conductivity(conductivity, PointParticle)" << endl;
  }

  DataOps::setValue(a_conductivity, 0.0);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();
    const int q   = species->getChargeNumber();

    if (q != 0 && solver->isMobile()) {
      DataOps::setValue(m_particleScratch1, 0.0);
#if 1 // Original code
      solver->depositParticles<PointParticle,
                               &PointParticle::weight>(m_particleScratch1,
                                                       *a_particles[idx]); // The particles should have "masses" = m*mu
#else
      const EBAMRCellData& mu  = solver->getMobilityFunction();
      const EBAMRCellData& phi = solver->getPhi();
      DataOps::copy(m_particleScratch1, mu);
      DataOps::multiply(m_particleScratch1, phi);
#endif

      // Copy to fluid Realm and add to fluid stuff
      m_fluidScratch1.copy(m_particleScratch1);
      DataOps::incr(a_conductivity, m_fluidScratch1, Abs(q));
    }
  }

  DataOps::scale(a_conductivity, Units::Qe);

  m_amr->conservativeAverage(a_conductivity, m_fluidRealm, m_phase);
  m_amr->interpGhostPwl(a_conductivity, m_fluidRealm, m_phase);

  // See if this helps....
  m_amr->interpToCentroids(a_conductivity, m_fluidRealm, m_phase);
}

void
ItoPlasmaGodunovStepper::compute_face_conductivity()
{
  CH_TIME("ItoPlasmaGodunovStepper::compute_face_conductivity");
  if (m_verbosity > 5) {
    pout() << m_name + "::compute_face_conductivity" << endl;
  }

  DataOps::setValue(m_conductivityFace, 0.0);
  DataOps::setValue(m_conductivityEB, 0.0);

  // This code does averaging from cell to face.
  DataOps::averageCellToFace(m_conductivityFace, m_conductivityCell, m_amr->getDomains());

  // This code extrapolates the conductivity to the EB. This should actually be the EB centroid but since the stencils
  // for EB extrapolation can be a bit nasty (e.g. Negative weights), we do the centroid instead and take that as an approximation.
#if 0
  const IrregAmrStencil<CentroidInterpolationStencil>& ebsten = m_amr->getCentroidInterpolationStencils(m_fluidRealm, m_phase);
  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    ebsten.apply(m_conductivityEB, m_conductivityCell, lvl);
  }
#else
  DataOps::incr(m_conductivityEB, m_conductivityCell, 1.0);
#endif
}

void
ItoPlasmaGodunovStepper::setupSemiImplicitPoisson(const Real a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::setupSemiImplicitPoisson");
  if (m_verbosity > 5) {
    pout() << m_name + "::setupSemiImplicitPoisson" << endl;
  }

  // Set coefficients as usual
  m_fieldSolver->setPermittivities();

  // Get bco and increment with mobilities
  MFAMRFluxData& bco     = m_fieldSolver->getPermittivityFace();
  MFAMRIVData&   bco_irr = m_fieldSolver->getPermittivityEB();

  EBAMRFluxData bco_gas;
  EBAMRIVData   bco_irr_gas;

  m_amr->allocatePointer(bco_gas);
  m_amr->allocatePointer(bco_irr_gas);

  m_amr->alias(bco_gas, phase::gas, bco);
  m_amr->alias(bco_irr_gas, phase::gas, bco_irr);

  DataOps::scale(m_conductivityFace, a_dt / Units::eps0);
  DataOps::scale(m_conductivityEB, a_dt / Units::eps0);

  DataOps::multiply(m_conductivityFace, bco_gas);
  DataOps::multiply(m_conductivityEB, bco_irr_gas);

  DataOps::incr(bco_gas, m_conductivityFace, 1.0);
  DataOps::incr(bco_irr_gas, m_conductivityEB, 1.0);

  m_amr->conservativeAverage(bco_gas, m_fluidRealm, phase::gas);
  m_amr->conservativeAverage(bco_irr_gas, m_fluidRealm, phase::gas);

  // Set up the solver
  m_fieldSolver->setupSolver();
}

void
ItoPlasmaGodunovStepper::setupStandardPoisson()
{
  CH_TIME("ItoPlasmaGodunovStepper::setupStandardPoisson");
  if (m_verbosity > 5) {
    pout() << m_name + "::setupStandardPoisson" << endl;
  }

  // Set coefficients as usual
  m_fieldSolver->setPermittivities();
  m_fieldSolver->setupSolver();
}

void
ItoPlasmaGodunovStepper::copyConductivityParticles(Vector<ParticleContainer<PointParticle>*>& a_conductivity_particles)
{
  CH_TIME("ItoPlasmaGodunovStepper::copyConductivityParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::copyConductivityParticles" << endl;
  }

  this->clearGodunovParticles(a_conductivity_particles, SpeciesSubset::All);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    const RefCountedPtr<ItoSolver>&  solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();
    const int q   = species->getChargeNumber();

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
        const List<ItoParticle>& ito_parts =
          solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl][dit()].listItems();
        List<PointParticle>& gdnv_parts = (*a_conductivity_particles[idx])[lvl][dit()].listItems();

        if (q != 0 && solver->isMobile()) {
          for (ListIterator<ItoParticle> lit(ito_parts); lit.ok(); ++lit) {
            const ItoParticle& p        = lit();
            const RealVect&    pos      = p.position();
            const Real&        weight   = p.weight();
            const Real&        mobility = p.mobility();

            gdnv_parts.add(PointParticle(pos, weight * mobility));
          }
        }
      }
    }
  }
}

void
ItoPlasmaGodunovStepper::copyRhoDaggerParticles(Vector<ParticleContainer<PointParticle>*>& a_rho_dagger_particles)
{
  CH_TIME("ItoPlasmaGodunovStepper::copyRhoDaggerParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::copyRhoDaggerParticles" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    const RefCountedPtr<ItoSolver>&  solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();
    const int q   = species->getChargeNumber();

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
        const List<ItoParticle>& ito_parts =
          solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl][dit()].listItems();
        List<PointParticle>& gdnv_parts = (*a_rho_dagger_particles[idx])[lvl][dit()].listItems();

        gdnv_parts.clear();

        if (q != 0) {
          for (ListIterator<ItoParticle> lit(ito_parts); lit.ok(); ++lit) {
            const ItoParticle& p      = lit();
            const RealVect&    pos    = p.position();
            const Real&        weight = p.weight();

            gdnv_parts.add(PointParticle(pos, weight));
          }
        }
      }
    }
  }
}

void
ItoPlasmaGodunovStepper::computeRegridConductivity()
{
  CH_TIME("ItoPlasmaGodunovStepper::computeRegridConductivity");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeRegridConductivity" << endl;
  }

  this->computeAllConductivities(m_conductivity_particles);
}

void
ItoPlasmaGodunovStepper::computeRegridRho()
{
  CH_TIME("ItoPlasmaGodunovStepper::computeRegridRho");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeRegridRho" << endl;
  }

  this->deposit_PointParticles(m_rho_dagger_particles, SpeciesSubset::All);
}

void
ItoPlasmaGodunovStepper::advanceParticlesEulerMaruyama(const Real a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::advanceParticlesEulerMaruyama");
  if (m_verbosity > 5) {
    pout() << m_name + "::advanceParticlesEulerMaruyama" << endl;
  }

  m_prevDt = a_dt; // Needed for regrids.

  // 1. Store X^k positions.
  this->setOldPositions();

  // 2. Diffuse the particles. This copies onto m_rho_dagger_particles and stores the hop on the full particles.
  this->diffuseParticlesEulerMaruyama(m_rho_dagger_particles, a_dt);
  this->remapGodunovParticles(m_rho_dagger_particles, SpeciesSubset::AllDiffusive);

  // 3. Solve the semi-implicit Poisson equation. Also, copy the particles used for computing the conductivity to scratch.
  this->copyConductivityParticles(m_conductivity_particles); // Sets particle "weights" = w*mu

  // Compute conductivity on mesh
  this->computeAllConductivities(m_conductivity_particles); // Deposits q_e*Z*w*mu on the mesh

  // Setup Poisson solver
  this->setupSemiImplicitPoisson(a_dt); // Multigrid setup

  // Compute space charge density
  // Diffusive should be enough because state is not changed for others.
  this->deposit_PointParticles(m_rho_dagger_particles, SpeciesSubset::AllDiffusive);

  this->solvePoisson(); // Solve the stinking equation.

  // 4. Recompute velocities with the new electric field, then do the actual semi-implicit Euler-Maruyama update.
#if 1 // This is what the algorithm says.
  this->setItoVelocityFunctions();
  m_ito->interpolateVelocities();
#else // Have to use this for LEA - need to debug.
  this->computeItoVelocities();
#endif

  this->stepEulerMaruyama(a_dt);

  this->remapParticles(SpeciesSubset::AllMobileOrDiffusive);

  // 5. Do intersection test and remove EB particles. These particles are NOT allowed to react later.
  const bool delete_eb_particles = true;
  this->intersectParticles(SpeciesSubset::AllMobileOrDiffusive,
                           EbRepresentation::ImplicitFunction,
                           delete_eb_particles);
  this->removeCoveredParticles(SpeciesSubset::AllMobileOrDiffusive,
                               EbRepresentation::ImplicitFunction,
                               m_eb_tolerance);

  // 6. Deposit particles. This shouldn't be necessary unless we want to compute (E,J)
  this->depositParticles(SpeciesSubset::AllMobileOrDiffusive);
}

void
ItoPlasmaGodunovStepper::diffuseParticlesEulerMaruyama(Vector<ParticleContainer<PointParticle>*>& a_rho_dagger,
                                                       const Real                                 a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::diffuseParticlesEulerMaruyama");
  if (m_verbosity > 5) {
    pout() << m_name + "::diffuseParticlesEulerMaruyama" << endl;
  }

  this->clearGodunovParticles(a_rho_dagger, SpeciesSubset::All);

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    const Real f = mobile ? 1.0 : 0.0;    // Multiplication factor for mobility
    const Real g = diffusive ? 1.0 : 0.0; // Multiplication factor for diffusion

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_particleRealm)[lvl];
      ParticleData<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

        List<ItoParticle>&   ItoParticles = particles[dit()].listItems();
        List<PointParticle>& gdnv_parts   = (*a_rho_dagger[idx])[lvl][dit()].listItems();

        if (diffusive) {
          for (ListIterator<ItoParticle> lit(ItoParticles); lit.ok(); ++lit) {
            ItoParticle&    p      = lit();
            const Real      factor = g * sqrt(2.0 * p.diffusion() * a_dt);
            const Real&     weight = p.weight();
            const RealVect& pos    = p.position();
#if 0
	    RealVect&       hop    = p.runtimeVector(0);
            hop                    = factor * solver->randomGaussian();


            // Add simpler particle
            gdnv_parts.add(PointParticle(pos + hop, weight));
#else
            MayDay::Error("ItoPlasmaGodunovStepper::diffuseParticlesEulerMaruayma -- runtime stuff");
#endif
          }
        }
        else { // Splitting up diffusion and non-diffusion because I dont want to generate random numbers where they're not required...
          for (ListIterator<ItoParticle> lit(ItoParticles); lit.ok(); ++lit) {
            ItoParticle&    p      = lit();
            const Real&     weight = p.weight();
            const RealVect& pos    = p.position();
#if 0
	    RealVect&       hop  = p.runtimeVector(0);
            hop                  = RealVect::Zero;

            // Add simpler particle
            gdnv_parts.add(PointParticle(pos, weight));
#else
            MayDay::Error("ItoPlasmaGodunovStepper::diffuseParticlesEulerMaruayma -- runtime stuff");
#endif
          }
        }
      }
    }
  }
}

void
ItoPlasmaGodunovStepper::stepEulerMaruyama(const Real a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::stepEulerMaruyama");
  if (m_verbosity > 5) {
    pout() << m_name + "::stepEulerMaruyama" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>& solver = solver_it();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    const Real f = mobile ? 1.0 : 0.0;
    const Real g = diffusive ? 1.0 : 0.0;

    if (mobile || diffusive) {

      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_particleRealm)[lvl];
        ParticleData<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

        for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

          List<ItoParticle>& particleList = particles[dit()].listItems();

          for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
            ItoParticle& p = lit();

            // Add diffusion hop again. The position after the diffusion hop is oldPosition() and X^k is in position()
#if 0
	    const RealVect& hop = p.runtimeVector(0);
            p.position()        = p.oldPosition() + f * p.velocity() * a_dt + g * hop;
#else
            MayDay::Error("ItoPlasmaGodunovStepper::stopEulerMaruyama -- runtime stuff");
#endif
          }
        }
      }
    }
  }
}

void
ItoPlasmaGodunovStepper::advanceParticlesTrapezoidal(const Real a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::advanceParticlesTrapezoidal");
  if (m_verbosity > 5) {
    pout() << m_name + "::advanceParticlesTrapezoidal" << endl;
  }

  m_prevDt = 0.5 * a_dt; // Needed for regrids.

  this->setOldPositions();

  // ====== PREDICTOR BEGIN ======
  this->preTrapezoidalPredictor(m_rho_dagger_particles, a_dt);
  this->remapGodunovParticles(m_rho_dagger_particles,
                              SpeciesSubset::
                                AllDiffusive); // Particles that were copied but not moved are in the right box.
  this->deposit_PointParticles(m_rho_dagger_particles, SpeciesSubset::All); // All copies need to deposit.

  this->copyConductivityParticles(m_conductivity_particles);
  this->computeAllConductivities(m_conductivity_particles);
  this->setupSemiImplicitPoisson(a_dt);
  this->solvePoisson();

  this->setItoVelocityFunctions();
  m_ito->interpolateVelocities();
  this->trapezoidalPredictor(a_dt);
  this->remapParticles(SpeciesSubset::AllMobileOrDiffusive);
  // ====== PREDICTOR END ======

  // ====== CORRECTOR BEGIN =====
  this->preTrapezoidalCorrector(m_rho_dagger_particles,
                                a_dt); // Mobile or diffusive moves to X^dagger = X^k + 0.5*dt*V^k + hop
  this->remapGodunovParticles(m_rho_dagger_particles,
                              SpeciesSubset::
                                AllMobileOrDiffusive); // Only need to remap particles that were mobile or diffusive
  this->deposit_PointParticles(m_rho_dagger_particles,
                               SpeciesSubset::All); // Everything needs to deposit...

  this->copyConductivityParticles(m_conductivity_particles);
  this->computeAllConductivities(m_conductivity_particles);
  this->setupSemiImplicitPoisson(0.5 * a_dt);
  this->solvePoisson();

  this->setItoVelocityFunctions();
  m_ito->interpolateVelocities();
  this->trapezoidalCorrector(a_dt);
  this->remapParticles(SpeciesSubset::AllMobileOrDiffusive);
  // ====== CORRECTOR END =====

  // Do particle-boundary intersection.
  this->intersectParticles(SpeciesSubset::AllMobileOrDiffusive, EbRepresentation::ImplicitFunction, true);
  this->removeCoveredParticles(SpeciesSubset::AllMobileOrDiffusive,
                               EbRepresentation::ImplicitFunction,
                               m_eb_tolerance);

  // Finally, deposit particles.
  this->depositParticles(SpeciesSubset::AllMobileOrDiffusive);
}

void
ItoPlasmaGodunovStepper::preTrapezoidalPredictor(Vector<ParticleContainer<PointParticle>*>& a_rho_dagger,
                                                 const Real                                 a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::preTrapezoidalPredictor");
  if (m_verbosity > 5) {
    pout() << m_name + "::preTrapezoidalPredictor" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    const Real f = mobile ? 1.0 : 0.0;    // Multiplication factor for mobility
    const Real g = diffusive ? 1.0 : 0.0; // Multiplication factor for diffusion

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_particleRealm)[lvl];
      ParticleData<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

        List<ItoParticle>&   ItoParticles = particles[dit()].listItems();
        List<PointParticle>& gdnv_parts   = (*a_rho_dagger[idx])[lvl][dit()].listItems();

        gdnv_parts.clear();

        // Store the diffusion hop, and add the godunov particles
        for (ListIterator<ItoParticle> lit(ItoParticles); lit.ok(); ++lit) {
          ItoParticle&    p      = lit();
          const Real      factor = sqrt(2.0 * p.diffusion() * a_dt);
          const RealVect  hop    = factor * solver->randomGaussian();
          const RealVect& Xk     = p.oldPosition();
          const Real&     weight = p.weight();

          // Store the diffusion hop and the current velocity.
#if 0
	  p.runtimeVector(0) = g * hop;
	  p.runtimeVector(1) = f * p.velocity();

          // Add simpler particle
          gdnv_parts.add(PointParticle(Xk + g * hop, weight));
#else
          MayDay::Error("ItoPlasmaGodunovStepper::preTrapezoidalPredictor -- runtime stuff");
#endif
        }
      }
    }
  }
}

void
ItoPlasmaGodunovStepper::trapezoidalPredictor(const Real a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::trapezoidalPredictor");
  if (m_verbosity > 5) {
    pout() << m_name + "::trapezoidalPredictor" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>& solver = solver_it();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    const Real f = mobile ? 1.0 : 0.0;
    const Real g = diffusive ? 1.0 : 0.0;

    if (mobile || diffusive) {

      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_particleRealm)[lvl];
        ParticleData<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

        for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

          List<ItoParticle>& particleList = particles[dit()].listItems();

          for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
            ItoParticle& p = lit();

            // Add diffusion hop again. The position after the diffusion hop is oldPosition() and X^k is in position()
#if 0
	    const RealVect& hop = p.runtimeVector(0);
	    const RealVect& Vk  = p.runtimeVector(1);

            p.position() = p.oldPosition() + f * p.velocity() * a_dt + g * hop;
#else
            MayDay::Error("ItoPlasmaGodunovStepper::trapezoidalPredictor -- runtime stuff");
#endif
          }
        }
      }
    }
  }
}

void
ItoPlasmaGodunovStepper::preTrapezoidalCorrector(Vector<ParticleContainer<PointParticle>*>& a_rho_dagger,
                                                 const Real                                 a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::preTrapezoidalCorrector");
  if (m_verbosity > 5) {
    pout() << m_name + "::preTrapezoidalCorrector" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>&        solver  = solver_it();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solver_it.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    const Real f = mobile ? 1.0 : 0.0;    // Multiplication factor for mobility
    const Real g = diffusive ? 1.0 : 0.0; // Multiplication factor for diffusion

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_particleRealm)[lvl];
      ParticleData<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

        List<ItoParticle>&   ItoParticles = particles[dit()].listItems();
        List<PointParticle>& gdnv_parts   = (*a_rho_dagger[idx])[lvl][dit()].listItems();

        gdnv_parts.clear();

        // Store the diffusion hop, and add the godunov particles
        for (ListIterator<ItoParticle> lit(ItoParticles); lit.ok(); ++lit) {
          ItoParticle& p = lit();

          const Real&     weight = p.weight();
          const RealVect& Xk     = p.oldPosition();
#if 0
	  const RealVect& hop  = p.runtimeVector(0);
	  const RealVect& Vk   = p.runtimeVector(1);

          // Move particle.
          const RealVect pos = Xk + 0.5 * a_dt * f * Vk + g * hop;
          gdnv_parts.add(PointParticle(pos, weight));
#else
          MayDay::Error("ItoPlasmaGodunovStepper::trapezoidalPredictor -- preTrapezoidalCorrector");
#endif
        }
      }
    }
  }
}

void
ItoPlasmaGodunovStepper::trapezoidalCorrector(const Real a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::trapezoidalCorrector");
  if (m_verbosity > 5) {
    pout() << m_name + "::trapezoidalCorrector" << endl;
  }

  for (auto solver_it = m_ito->iterator(); solver_it.ok(); ++solver_it) {
    RefCountedPtr<ItoSolver>& solver = solver_it();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    const Real f = mobile ? 1.0 : 0.0;
    const Real g = diffusive ? 1.0 : 0.0;

    if (mobile || diffusive) {

      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_particleRealm)[lvl];
        ParticleData<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

        for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

          List<ItoParticle>& particleList = particles[dit()].listItems();

          for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
            ItoParticle& p = lit();

            const RealVect& Xk = p.oldPosition();
#if 0
	    const RealVect& hop = p.runtimeVector(0);
	    const RealVect& Vk  = p.runtimeVector(1);
            const RealVect& Vk1 = p.velocity();

            p.position() = Xk + 0.5 * f * a_dt * (Vk + Vk1) + g * hop;
#else
            MayDay::Error("ItoPlasmaGodunovStepper::trapezoidalCorrector -- runtime stuff");
#endif
          }
        }
      }
    }
  }
}

#include <CD_NamespaceFooter.H>
