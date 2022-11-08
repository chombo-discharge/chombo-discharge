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
#include <CD_ParallelOps.H>
#include <CD_DataOps.H>
#include <CD_Units.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::ItoPlasma;

ItoPlasmaGodunovStepper::ItoPlasmaGodunovStepper(RefCountedPtr<ItoPlasmaPhysics>& a_physics)
  : ItoPlasmaStepper(a_physics)
{
  CH_TIME("ItoPlasmaGodunovStepper::ItoPlasmaGodunovStepper");

  m_name = "ItoPlasmaGodunovStepper";

  this->parseOptions();
}

ItoPlasmaGodunovStepper::~ItoPlasmaGodunovStepper()
{
  CH_TIME("ItoPlasmaGodunovStepper::~ItoPlasmaGodunovStepper");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaGodunovStepper::~ItoPlasmaGodunovStepper" << endl;
  }
}

void
ItoPlasmaGodunovStepper::allocate()
{
  CH_TIME("ItoPlasmaGodunovStepper::allocate");
  if (m_verbosity > 5) {
    pout() << "ItoPlasmaGodunovStepper::allocate" << endl;
  }

  ItoPlasmaStepper::allocate();

  // Now allocate for the conductivity particles and rho^dagger particles. This is only done in the 'allocate' routine
  // and not in 'allocateInternals' because that would discard the particles during regrids. That has definitely never
  // happen, and there's no way I've spent countless hours tracking down such a bug.
  const int numItoSpecies = m_physics->getNumPlasmaSpecies();

  m_conductivityParticles.resize(numItoSpecies);
  m_rhoDaggerParticles.resize(numItoSpecies);

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    const RefCountedPtr<ItoSolver>& solver = solverIt();

    const int idx = solverIt.index();

    m_conductivityParticles[idx] = new ParticleContainer<PointParticle>();
    m_rhoDaggerParticles[idx]    = new ParticleContainer<PointParticle>();

    m_amr->allocate(*m_conductivityParticles[idx], m_particleRealm);
    m_amr->allocate(*m_rhoDaggerParticles[idx], m_particleRealm);
  }
}

void
ItoPlasmaGodunovStepper::barrier() const noexcept
{
  CH_TIME("ItoPlasmaGodunovStepper::barrier");
  if (m_verbosity > 5) {
    pout() << m_name + "::barrier" << endl;
  }

  if (m_profile) {
    ParallelOps::barrier();
  }
}

void
ItoPlasmaGodunovStepper::parseOptions()
{
  CH_TIME("ItoPlasmaGodunovStepper::parseOptions");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseOptions" << endl;
  }

  ItoPlasmaStepper::parseOptions();

  this->parseAlgorithm();
}

void
ItoPlasmaGodunovStepper::parseRuntimeOptions()
{
  CH_TIME("ItoPlasmaGodunovStepper::parseRuntimeOptions");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseRuntimeOptions" << endl;
  }

  ItoPlasmaStepper::parseRuntimeOptions();

  this->parseAlgorithm();

  m_ito->parseRuntimeOptions();
  m_fieldSolver->parseRuntimeOptions();
  m_rte->parseRuntimeOptions();
  m_sigmaSolver->parseRuntimeOptions();
}

void
ItoPlasmaGodunovStepper::parseAlgorithm() noexcept
{
  CH_TIME("ItoPlasmaGodunovStepper::parseAlgorithm");
  if (m_verbosity > 5) {
    pout() << m_name + "::parseAlgorithm" << endl;
  }

  ParmParse   pp(m_name.c_str());
  std::string str;

  pp.get("algorithm", str);

  // Get algorithm
  if (str == "euler_maruyama") {
    m_algorithm = WhichAlgorithm::EulerMaruyama;
  }
  else {
    MayDay::Abort("ItoPlasmaGodunovStepper::parseAlgorithm - unknown algorithm requested");
  }
}

Real
ItoPlasmaGodunovStepper::advance(const Real a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::advance");
  if (m_verbosity > 5) {
    pout() << m_name + "::advance" << endl;
  }

  m_timer = Timer("ItoPlasmaGodunovStepper::advance");

  // Previous time step is needed when regridding.
  m_prevDt = a_dt;

  // Done only so we can plot the absorbed photons (advanceReactionNetwork absorbs them)
  m_timer.startEvent("Deposit photons");
  for (auto solverIt = m_rte->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->depositPhotons(solverIt()->getPhi(), solverIt()->getBulkPhotons(), DepositionType::NGP);
  }
  m_timer.stopEvent("Deposit photons");

  // ====== BEGIN TRANSPORT STEP ======
  // Semi-implicitly advance the particles and the field.
  switch (m_algorithm) {
  case WhichAlgorithm::EulerMaruyama: {
    this->advanceParticlesEulerMaruyama(a_dt);

    break;
  }
  default: {
    MayDay::Abort("ItoPlasmaGodunovStepper::advance - logic bust");

    break;
  }
  }

  // Remove the run-time configurable particle storage. It is no longer needed.
  // ====== END TRANSPORT STEP ======

  // Photon transport
  this->barrier();
  m_timer.startEvent("Photon transport");
  this->advancePhotons(a_dt);
  m_timer.stopEvent("Photon transport");

  // Sort the particles and photons per cell so we can call reaction algorithms
  this->barrier();
  m_timer.startEvent("Sort by cell");
  m_ito->organizeParticlesByCell(ItoSolver::WhichContainer::Bulk);
  this->sortPhotonsByCell(McPhoto::WhichContainer::Bulk);
  this->sortPhotonsByCell(McPhoto::WhichContainer::Source);
  m_timer.stopEvent("Sort by cell");

  // Run the Kinetic Monte Carlo reaction kernels.
  this->barrier();
  m_timer.startEvent("Reaction network");
  this->advanceReactionNetwork(a_dt);
  m_timer.stopEvent("Reaction network");

  // Build superparticles.
  if ((m_timeStep + 1) % m_mergeInterval == 0 && m_mergeInterval > 0) {
    this->barrier();
    m_timer.startEvent("Super-particle management");
    m_ito->makeSuperparticles(ItoSolver::WhichContainer::Bulk, m_particlesPerCell);
    m_timer.stopEvent("Super-particle management");
  }

  // Sort particles per patch.
  this->barrier();
  m_timer.startEvent("Sort by patch");
  m_ito->organizeParticlesByPatch(ItoSolver::WhichContainer::Bulk);
  this->sortPhotonsByPatch(McPhoto::WhichContainer::Bulk);
  this->sortPhotonsByPatch(McPhoto::WhichContainer::Source);
  m_timer.stopEvent("Sort by patch");

  // Clear other data holders for now. BC comes later...
  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->clear(ItoSolver::WhichContainer::EB);
    solverIt()->clear(ItoSolver::WhichContainer::Domain);
  }

  // Prepare for the next time step
  this->barrier();
  m_timer.startEvent("Post-compute v");
  this->computeItoVelocities();
  m_timer.stopEvent("Post-compute v");

  this->barrier();
  m_timer.startEvent("Post-compute D");
  this->computeItoDiffusion();
  m_timer.stopEvent("Post-compute D");

  this->barrier();
  m_timer.startEvent("Compute J");
  this->computeCurrentDensity(m_currentDensity);
  m_timer.stopEvent("Compute J");

  if (m_profile) {
    m_timer.eventReport(pout(), false);
  }

  m_timer.clear();

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

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    const int idx = solverIt.index();

    m_conductivityParticles[idx]->preRegrid(a_lmin);
    m_rhoDaggerParticles[idx]->preRegrid(a_lmin);
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
  m_sigmaSolver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);

  // Allocate internal memory for ItoPlasmaGodunovStepper now....
  this->allocateInternals();

  // We need to remap/regrid the stored particles as well.
  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    const int idx = solverIt.index();
    m_amr->remapToNewGrids(*m_rhoDaggerParticles[idx], a_lmin, a_newFinestLevel);
    m_amr->remapToNewGrids(*m_conductivityParticles[idx], a_lmin, a_newFinestLevel);
  }

  // Set up the field solver
  m_fieldSolver->setupSolver();

  // Recompute the conductivity and space charge densities.
  this->computeConductivities(m_conductivityParticles);
  this->depositPointParticles(m_rhoDaggerParticles, SpeciesSubset::All);
  this->setupSemiImplicitPoisson(m_prevDt);

  // Solve the Poisson equation.
  const bool converged = this->solvePoisson();
  if (!converged) {
    MayDay::Warning("ItoPlasmaGodunovStepper::regrid - Poisson solve did not converge after regrid!!!");
  }

  // Regrid superparticles.
  if (m_regridSuperparticles) {
    m_ito->organizeParticlesByCell(ItoSolver::WhichContainer::Bulk);
    m_ito->makeSuperparticles(ItoSolver::WhichContainer::Bulk, m_particlesPerCell);
    m_ito->organizeParticlesByPatch(ItoSolver::WhichContainer::Bulk);
  }

  // Now let Ihe ito solver deposit its actual particles... In the above it deposit m_rhoDaggerParticles.
  m_ito->depositParticles();

  // Recompute new velocities and diffusion coefficients
  this->computeItoVelocities();
  this->computeItoDiffusion();
}

void
ItoPlasmaGodunovStepper::setOldPositions() noexcept
{
  CH_TIME("ItoPlasmaGodunovStepper::setOldPositions");
  if (m_verbosity > 5) {
    pout() << m_name + "::setOldPositions" << endl;
  }

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>& solver = solverIt();

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_particleRealm)[lvl];
      ParticleData<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit) {

        List<ItoParticle>& particleList = particles[dit()].listItems();

        for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
          ItoParticle& p = lit();

          p.oldPosition() = p.position();
        }
      }
    }
  }
}

void
ItoPlasmaGodunovStepper::remapPointParticles(Vector<ParticleContainer<PointParticle>*>& a_particles,
                                             const SpeciesSubset                        a_subset) noexcept
{
  CH_TIME("ItoPlasmaGodunovStepper::remapPointParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::remapPointParticles" << endl;
  }

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
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
      MayDay::Abort("ItoPlasmaGodunovStepper::remapPointParticles - logic bust");

      break;
    }
    }
  }
}

void
ItoPlasmaGodunovStepper::depositPointParticles(const Vector<ParticleContainer<PointParticle>*>& a_particles,
                                               const SpeciesSubset                              a_subset) noexcept
{
  CH_TIME("ItoPlasmaGodunovStepper::depositPointParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::depositPointParticles" << endl;
  }

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();
    const bool charged   = species->getChargeNumber() != 0;

    switch (a_subset) {
    case SpeciesSubset::All: {
      solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);

      break;
    }
    case SpeciesSubset::AllMobile: {
      if (mobile) {
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      }

      break;
    }
    case SpeciesSubset::AllDiffusive: {
      if (diffusive) {
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      }

      break;
    }
    case SpeciesSubset::AllMobileOrDiffusive: {
      if (mobile || diffusive) {
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      }
      break;
    }
    case SpeciesSubset::AllMobileAndDiffusive: {
      if (mobile && diffusive) {
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      }
      break;
    }
    case SpeciesSubset::Charged: {
      if (charged) {
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      }
      break;
    }
    case SpeciesSubset::ChargedMobile: {
      if (charged && mobile) {
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      }
      break;
    }
    case SpeciesSubset::ChargedDiffusive: {
      if (charged && diffusive) {
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      }

      break;
    }
    case SpeciesSubset::ChargedMobileOrDiffusive: {
      if (charged && (mobile || diffusive)) {
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      }

      break;
    }
    case SpeciesSubset::ChargedMobileAndDiffusive: {
      if (charged && (mobile && diffusive)) {
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      }

      break;
    }
    case SpeciesSubset::Stationary: {
      if (!mobile && !diffusive) {
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      }

      break;
    }
    default: {
      MayDay::Abort("ItoPlasmaGodunovStepper::depositPointParticles - logic bust");

      break;
    }
    }
  }
}

void
ItoPlasmaGodunovStepper::clearPointParticles(const Vector<ParticleContainer<PointParticle>*>& a_particles,
                                             const SpeciesSubset                              a_subset) noexcept
{
  CH_TIME("ItoPlasmaGodunovStepper::clearPointParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::clearPointParticles" << endl;
  }

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
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
      MayDay::Abort("ItoPlasmaGodunovStepper::clearPointParticles - logic bust");

      break;
    }
    }
  }
}

void
ItoPlasmaGodunovStepper::computeConductivities(const Vector<ParticleContainer<PointParticle>*>& a_particles) noexcept
{
  CH_TIME("ItoPlasmaGodunovStepper::computeConductivities");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeConductivities" << endl;
  }

  this->computeCellConductivity(m_conductivityCell, a_particles);

  // Now do the faces
  this->computeFaceConductivity();
}

void
ItoPlasmaGodunovStepper::computeCellConductivity(EBAMRCellData&                                   a_conductivityCell,
                                                 const Vector<ParticleContainer<PointParticle>*>& a_particles) noexcept
{
  CH_TIME("ItoPlasmaGodunovStepper::computeCellConductivity(EBAMRCellData, PointParticle");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeCellConductivity(EBAMRCellData, PointParticle)" << endl;
  }

  DataOps::setValue(a_conductivityCell, 0.0);

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();
    const int Z   = species->getChargeNumber();

    if (Z != 0 && solver->isMobile()) {
      // Deposit on the particle realm.
      DataOps::setValue(m_particleScratch1, 0.0);
      solver->depositParticles<PointParticle, &PointParticle::weight>(m_particleScratch1, *a_particles[idx]);

      // Copy to fluid realm and add to total conductivity.
      m_fluidScratch1.copy(m_particleScratch1);
      DataOps::incr(a_conductivityCell, m_fluidScratch1, 1.0 * std::abs(Z));
    }
  }

  // Conductivity is mobility * weight * Q
  DataOps::scale(a_conductivityCell, Units::Qe);

  // Coarsen, update ghost cells and interpolate to centroids
  m_amr->conservativeAverage(a_conductivityCell, m_fluidRealm, m_plasmaPhase);
  m_amr->interpGhostMG(a_conductivityCell, m_fluidRealm, m_plasmaPhase);

  m_amr->interpToCentroids(a_conductivityCell, m_fluidRealm, m_plasmaPhase);
}

void
ItoPlasmaGodunovStepper::computeFaceConductivity() noexcept
{
  CH_TIME("ItoPlasmaGodunovStepper::computeFaceConductivity");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeFaceConductivity" << endl;
  }

  DataOps::setValue(m_conductivityFace, 0.0);
  DataOps::setValue(m_conductivityEB, 0.0);

  // Average the cell-centered conductivity to faces. Note that this includes one "ghost face", which we need
  // because the multigrid solver will interpolate face-centered conductivities to face centroids.
  const Average  average  = Average::Arithmetic;
  const int      tanGhost = 1;
  const Interval interv(0, 0);

  DataOps::averageCellToFace(m_conductivityFace,
                             m_conductivityCell,
                             m_amr->getDomains(),
                             tanGhost,
                             interv,
                             interv,
                             average);

  // Set the EB conductivity.
  DataOps::incr(m_conductivityEB, m_conductivityCell, 1.0);
}

void
ItoPlasmaGodunovStepper::setupSemiImplicitPoisson(const Real a_dt) noexcept
{
  CH_TIME("ItoPlasmaGodunovStepper::setupSemiImplicitPoisson");
  if (m_verbosity > 5) {
    pout() << m_name + "::setupSemiImplicitPoisson" << endl;
  }

  // Set coefficients as usual
  m_fieldSolver->setPermittivities();

  // Get the permittivities
  // Get the permittivities on the faces.
  MFAMRCellData& permCell = m_fieldSolver->getPermittivityCell();
  MFAMRFluxData& permFace = m_fieldSolver->getPermittivityFace();
  MFAMRIVData&   permEB   = m_fieldSolver->getPermittivityEB();

  // Get handles to the gas-phase permittivities.
  EBAMRFluxData permFaceGas = m_amr->alias(m_plasmaPhase, permFace);
  EBAMRIVData   permEBGas   = m_amr->alias(m_plasmaPhase, permEB);

  // Increment the field solver permittivities by a_factor*sigma. After this, the "permittivities" are
  // given by epsr + a_factor*sigma
  DataOps::incr(permFaceGas, m_conductivityFace, a_dt / Units::eps0);
  DataOps::incr(permEBGas, m_conductivityEB, a_dt / Units::eps0);

  // Coarsen coefficients.
  m_amr->arithmeticAverage(permFaceGas, m_fluidRealm, m_plasmaPhase);
  m_amr->arithmeticAverage(permEBGas, m_fluidRealm, m_plasmaPhase);

  // Set up the solver with the "permittivities"
  m_fieldSolver->setSolverPermittivities(permCell, permFace, permEB);
}

void
ItoPlasmaGodunovStepper::copyConductivityParticles(
  Vector<ParticleContainer<PointParticle>*>& a_conductivityParticles) noexcept
{
  CH_TIME("ItoPlasmaGodunovStepper::copyConductivityParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::copyConductivityParticles" << endl;
  }

  this->clearPointParticles(a_conductivityParticles, SpeciesSubset::All);

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    const RefCountedPtr<ItoSolver>&  solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();
    const int Z   = species->getChargeNumber();

    if (Z != 0 && solver->isMobile()) {
      const ParticleContainer<ItoParticle>& solverParticles = solver->getParticles(ItoSolver::WhichContainer::Bulk);

      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[lvl];

        for (DataIterator dit(dbl); dit.ok(); ++dit) {
          const List<ItoParticle>& patchParticles = solverParticles[lvl][dit()].listItems();

          List<PointParticle>& pointParticles = (*a_conductivityParticles[idx])[lvl][dit()].listItems();

          for (ListIterator<ItoParticle> lit(patchParticles); lit.ok(); ++lit) {
            const ItoParticle& p        = lit();
            const RealVect&    pos      = p.position();
            const Real&        weight   = p.weight();
            const Real&        mobility = p.mobility();

            pointParticles.add(PointParticle(pos, weight * mobility));
          }
        }
      }
    }
  }
}

void
ItoPlasmaGodunovStepper::copyRhoDaggerParticles(
  Vector<ParticleContainer<PointParticle>*>& a_rhoDaggerParticles) noexcept
{
  CH_TIME("ItoPlasmaGodunovStepper::copyRhoDaggerParticles");
  if (m_verbosity > 5) {
    pout() << m_name + "::copyRhoDaggerParticles" << endl;
  }

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    const RefCountedPtr<ItoSolver>&  solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();
    const int Z   = species->getChargeNumber();

    if (Z != 0) {
      const ParticleContainer<ItoParticle>& solverParticles = solver->getParticles(ItoSolver::WhichContainer::Bulk);

      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[lvl];

        for (DataIterator dit(dbl); dit.ok(); ++dit) {
          const List<ItoParticle>& patchParticles = solverParticles[lvl][dit()].listItems();

          List<PointParticle>& pointParticles = (*a_rhoDaggerParticles[idx])[lvl][dit()].listItems();

          pointParticles.clear();

          for (ListIterator<ItoParticle> lit(patchParticles); lit.ok(); ++lit) {
            const ItoParticle& p      = lit();
            const RealVect&    pos    = p.position();
            const Real&        weight = p.weight();

            pointParticles.add(PointParticle(pos, weight));
          }
        }
      }
    }
  }
}

void
ItoPlasmaGodunovStepper::advanceParticlesEulerMaruyama(const Real a_dt) noexcept
{
  CH_TIME("ItoPlasmaGodunovStepper::advanceParticlesEulerMaruyama");
  if (m_verbosity > 5) {
    pout() << m_name + "::advanceParticlesEulerMaruyama" << endl;
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
  m_timer.stopEvent("Deposit point particles");

  // Solve the stinking equation.
  this->barrier();
  m_timer.startEvent("Solve Poisson");
  this->solvePoisson();
  m_timer.stopEvent("Solve Poisson");

  // Recompute velocities with the new electric field. This interpolates the velocities to the current particle positions, i.e.
  // we compute V^(k+1)(X^k) = mu^k * E^(k+1)(X^k)
  this->barrier();
  m_timer.startEvent("Step-compute v");
#if 1 // This is what the algorithm says.
  this->setItoVelocityFunctions();
  m_ito->interpolateVelocities();
#else // Have to use this for LEA - need to debug.
  this->computeItoVelocities();
#endif
  m_timer.stopEvent("Step-compute v");

  // Finalize the Euler-Maruyama update.
  this->barrier();
  m_timer.startEvent("Euler-Maruyama step");
  this->stepEulerMaruyama(a_dt);
  this->remapParticles(SpeciesSubset::AllMobileOrDiffusive);
  m_timer.stopEvent("Euler-Maruyama step");

  // Do intersection test and remove EB particles. These particles are NOT allowed to react later.
  this->barrier();
  m_timer.startEvent("EB/Particle intersection");
  const bool deleteParticles = true;
  this->intersectParticles(SpeciesSubset::AllMobileOrDiffusive, EBIntersection::Bisection, deleteParticles);
  m_timer.stopEvent("EB/Particle intersection");

  this->barrier();
  m_timer.startEvent("Remove covered");
  this->removeCoveredParticles(SpeciesSubset::AllMobileOrDiffusive, EBRepresentation::ImplicitFunction, m_toleranceEB);
  m_timer.stopEvent("Remove covered");
}

void
ItoPlasmaGodunovStepper::diffuseParticlesEulerMaruyama(Vector<ParticleContainer<PointParticle>*>& a_rhoDaggerParticles,
                                                       const Real                                 a_dt) noexcept
{
  CH_TIME("ItoPlasmaGodunovStepper::diffuseParticlesEulerMaruyama");
  if (m_verbosity > 5) {
    pout() << m_name + "::diffuseParticlesEulerMaruyama" << endl;
  }

  this->clearPointParticles(a_rhoDaggerParticles, SpeciesSubset::All);

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();

    const bool diffusive = solver->isDiffusive();
    const int  Z         = species->getChargeNumber();

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_particleRealm)[lvl];
      ParticleData<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        List<ItoParticle>&   itoParticles   = particles[dit()].listItems();
        List<PointParticle>& pointParticles = (*a_rhoDaggerParticles[idx])[lvl][dit()].listItems();

        for (ListIterator<ItoParticle> lit(itoParticles); lit.ok(); ++lit) {
          ItoParticle&    p      = lit();
          const Real&     weight = p.weight();
          const RealVect& pos    = p.position();

          // Compute a particle hop and store it on the run-time storage.
          RealVect& hop = p.tmpVect();
          if (diffusive) {
            hop = sqrt(2.0 * p.diffusion() * a_dt) * solver->randomGaussian();
          }
          else {
            hop = RealVect::Zero;
          }

          if (Z != 0) {
            pointParticles.add(PointParticle(pos + hop, weight));
          }
        }
      }
    }
  }
}

void
ItoPlasmaGodunovStepper::stepEulerMaruyama(const Real a_dt) noexcept
{
  CH_TIME("ItoPlasmaGodunovStepper::stepEulerMaruyama");
  if (m_verbosity > 5) {
    pout() << m_name + "::stepEulerMaruyama" << endl;
  }

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>& solver = solverIt();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    const Real f = mobile ? a_dt : 0.0;
    const Real g = diffusive ? 1.0 : 0.0;

    if (mobile || diffusive) {
      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
        const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_particleRealm)[lvl];
        ParticleData<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

        for (DataIterator dit(dbl); dit.ok(); ++dit) {
          List<ItoParticle>& particleList = particles[dit()].listItems();

          for (ListIterator<ItoParticle> lit(particleList); lit.ok(); ++lit) {
            ItoParticle& p = lit();

            // Add in the diffusion hop and advective contribution.
            const RealVect& hop = p.tmpVect();
            p.position()        = p.oldPosition() + f * p.velocity() + g * hop;
          }
        }
      }
    }
  }
}

#include <CD_NamespaceFooter.H>
