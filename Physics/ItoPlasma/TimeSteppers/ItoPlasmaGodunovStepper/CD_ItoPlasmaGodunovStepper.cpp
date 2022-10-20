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
  CH_TIME("ItoPlasmaGodunovStepper::ItoPlasmaGodunovStepper");

  m_name = "ItoPlasmaGodunovStepper";
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

  this->allocateInternals();
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
  else if (str == "trapezoidal") {
    m_algorithm = WhichAlgorithm::Trapezoidal;
  }
  else {
    MayDay::Abort("ItoPlasmaGodunovStepper::parseAlgorithm - unknown algorithm requested");
  }
}

void
ItoPlasmaGodunovStepper::allocateInternals()
{
  CH_TIME("ItoPlasmaGodunovStepper::allocateInternals");
  if (m_verbosity > 5) {
    pout() << m_name + "::allocateInternals" << endl;
  }

  // Call parent method first
  ItoPlasmaStepper::allocateInternals();

  // Now allocate for the conductivity particles and rho^dagger particles
  const int numItoSpecies = m_physics->getNumItoSpecies();

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

Real
ItoPlasmaGodunovStepper::advance(const Real a_dt)
{
  CH_TIME("ItoPlasmaGodunovStepper::advance");
  if (m_verbosity > 5) {
    pout() << m_name + "::advance" << endl;
  }

  // Setup runtime storage (requirements change with algorithm)
  this->setRuntimeParticleStorage();

  // Semi-implicitly advance the particles and the field.
  switch (m_algorithm) {
  case WhichAlgorithm::EulerMaruyama: {
    this->advanceParticlesEulerMaruyama(a_dt);

    break;
  }
  case WhichAlgorithm::Trapezoidal: {
    this->advanceParticlesTrapezoidal(a_dt);

    break;
  }
  default: {
    MayDay::Abort("ItoPlasmaGodunovStepper::advance - logic bust");

    break;
  }
  }

  // Remove the run-time configurable particle storage.
  this->resetRuntimeParticleStorage();

  // Compute current and relaxation time.
  this->computeCurrentDensity(m_currentDensity);
  const Real relaxTime = this->computeRelaxationTime(); // This is for the restricting the next step.

  // Do the radiative transfer advance.
  this->advancePhotons(a_dt);

  // If we are using the LEA, we must compute the Ohmic heating term. This must be done
  // BEFORE sorting the particles per cell.
  if (m_physics->getCoupling() == ItoPlasmaPhysics::coupling::LEA) {
    this->computeEdotJSource(a_dt);
  }

  // Sort the particles and photons per cell so we can call reaction algorithms
  m_ito->sortParticlesByCell(ItoSolver::WhichContainer::Bulk);
  this->sortPhotonsByCell(McPhoto::WhichContainer::Bulk);
  this->sortPhotonsByCell(McPhoto::WhichContainer::Source);

  // Run the Kinetic Monte Carlo reaction kernels.
  this->advanceReactionNetwork(a_dt);

  // Build superparticles.
  if ((m_timeStep + 1) % m_mergeInterval == 0 && m_mergeInterval > 0) {
    m_ito->makeSuperparticles(ItoSolver::WhichContainer::Bulk, m_particlesPerCell);
  }

  // Sort particles per patch.
  m_ito->sortParticlesByPatch(ItoSolver::WhichContainer::Bulk);
  this->sortPhotonsByPatch(McPhoto::WhichContainer::Bulk);
  this->sortPhotonsByPatch(McPhoto::WhichContainer::Source);

  // Clear other data holders for now. BC comes later...
  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    solverIt()->clear(ItoSolver::WhichContainer::EB);
    solverIt()->clear(ItoSolver::WhichContainer::Domain);
  }

  // Deposit particles on the mesh.
  m_ito->depositParticles();

  // Prepare for the next time step
  this->computeItoVelocities();
  this->computeItoDiffusion();

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

  // Now let Ihe ito solver deposit its actual particles... In the above it deposit m_rhoDaggerParticles.
  m_ito->depositParticles();

  // Recompute new velocities and diffusion coefficients
  this->computeItoVelocities();
  this->computeItoDiffusion();
}

void
ItoPlasmaGodunovStepper::setRuntimeParticleStorage() noexcept
{
  CH_TIME("ItoPlasmaGodunovStepper::setRuntimeParticleStorage");
  if (m_verbosity > 5) {
    pout() << m_name + "::setupRuntimeParticleStorage" << endl;
  }
#if 0
  switch (m_algorithm) {
  case WhichAlgorithm::EulerMaruyama:
    ItoParticle::setNumRuntimeVectors(1);
    break;
  case WhichAlgorithm::Trapezoidal:
    ItoParticle::setNumRuntimeVectors(2); // For V^k and the diffusion hop.
    break;
  default:
    MayDay::Abort("ItoPlasmaGodunovStepper::setRuntimeParticleStorage - logic bust");
  }
#else
  MayDay::Error(
    "ItoPlasmaGodunovStepper::setRuntimeParticleStorage -- need to figure out how to add more fields to ItoParticle");
#endif
}

void
ItoPlasmaGodunovStepper::resetRuntimeParticleStorage() noexcept
{
  CH_TIME("ItoPlasmaGodunovStepper::resetRuntimeParticleStorage");
  if (m_verbosity > 5) {
    pout() << m_name + "::resetRuntimeParticleStorage" << endl;
  }
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
    case SpeciesSubset::AllMobileOrDiffusive: {
      if (mobile || diffusive) {
        a_particles[idx]->remap();
      }

      break;
    }
    case SpeciesSubset::ChargedAndMobileOrDiffusive: {
      if (charged && (mobile || diffusive)) {
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
    case SpeciesSubset::AllMobileOrDiffusive: {
      if (mobile || diffusive) {
        solver->depositParticles<PointParticle, &PointParticle::weight>(solver->getPhi(), *a_particles[idx]);
      }
      break;
    }
    case SpeciesSubset::ChargedAndMobileOrDiffusive: {
      if (charged && (mobile || diffusive)) {
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
    case SpeciesSubset::AllMobileOrDiffusive: {
      if (mobile || diffusive) {
        a_particles[idx]->clearParticles();
      }

      break;
    }
    case SpeciesSubset::ChargedAndMobileOrDiffusive: {
      if (charged && (mobile || diffusive)) {
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

  // Coarsen coefficients
  m_amr->arithmeticAverage(permFaceGas, m_fluidRealm, phase::gas);
  m_amr->conservativeAverage(permEBGas, m_fluidRealm, phase::gas);

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

    if(Z > 0 && solver->isMobile()) {
      const ParticleContainer<ItoParticle>& solverParticles = solver->getParticles(ItoSolver::WhichContainer::Bulk);

      for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
	const DisjointBoxLayout& dbl = m_amr->getGrids(m_particleRealm)[lvl];

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
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
ItoPlasmaGodunovStepper::copyRhoDaggerParticles(Vector<ParticleContainer<PointParticle>*>& a_rhoDaggerParticles) noexcept
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

	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
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
ItoPlasmaGodunovStepper::computeRegridConductivity() noexcept
{
  CH_TIME("ItoPlasmaGodunovStepper::computeRegridConductivity");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeRegridConductivity" << endl;
  }

  this->computeConductivities(m_conductivityParticles);
}

void
ItoPlasmaGodunovStepper::computeRegridRho() noexcept
{
  CH_TIME("ItoPlasmaGodunovStepper::computeRegridRho");
  if (m_verbosity > 5) {
    pout() << m_name + "::computeRegridRho" << endl;
  }

  this->depositPointParticles(m_rhoDaggerParticles, SpeciesSubset::All);
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

  // 2. Diffuse the particles. This copies onto m_rhoDaggerParticles and stores the hop on the full particles.
  this->diffuseParticlesEulerMaruyama(m_rhoDaggerParticles, a_dt);
  this->remapPointParticles(m_rhoDaggerParticles, SpeciesSubset::AllDiffusive);

  // 3. Solve the semi-implicit Poisson equation. Also, copy the particles used for computing the conductivity to scratch.
  this->copyConductivityParticles(m_conductivityParticles); // Sets particle "weights" = w*mu

  // Compute conductivity on mesh
  this->computeConductivities(m_conductivityParticles); // Deposits q_e*Z*w*mu on the mesh

  // Setup Poisson solver
  this->setupSemiImplicitPoisson(a_dt); // Multigrid setup

  // Compute space charge density
  // Diffusive should be enough because state is not changed for others.
  this->depositPointParticles(m_rhoDaggerParticles, SpeciesSubset::AllDiffusive);

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
  this->intersectParticles(SpeciesSubset::AllMobileOrDiffusive, EBIntersection::Bisection, delete_eb_particles);
  this->removeCoveredParticles(SpeciesSubset::AllMobileOrDiffusive, EBRepresentation::ImplicitFunction, m_toleranceEB);

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

  this->clearPointParticles(a_rho_dagger, SpeciesSubset::All);

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    const Real f = mobile ? 1.0 : 0.0;    // Multiplication factor for mobility
    const Real g = diffusive ? 1.0 : 0.0; // Multiplication factor for diffusion

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_particleRealm)[lvl];
      ParticleData<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

        List<ItoParticle>&   ItoParticles   = particles[dit()].listItems();
        List<PointParticle>& pointParticles = (*a_rho_dagger[idx])[lvl][dit()].listItems();

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
            pointParticles.add(PointParticle(pos + hop, weight));
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
            pointParticles.add(PointParticle(pos, weight));
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

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>& solver = solverIt();

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
  this->preTrapezoidalPredictor(m_rhoDaggerParticles, a_dt);
  this->remapPointParticles(m_rhoDaggerParticles,
                            SpeciesSubset::
                              AllDiffusive); // Particles that were copied but not moved are in the right box.
  this->depositPointParticles(m_rhoDaggerParticles, SpeciesSubset::All); // All copies need to deposit.

  this->copyConductivityParticles(m_conductivityParticles);
  this->computeConductivities(m_conductivityParticles);
  this->setupSemiImplicitPoisson(a_dt);
  this->solvePoisson();

  this->setItoVelocityFunctions();
  m_ito->interpolateVelocities();
  this->trapezoidalPredictor(a_dt);
  this->remapParticles(SpeciesSubset::AllMobileOrDiffusive);
  // ====== PREDICTOR END ======

  // ====== CORRECTOR BEGIN =====
  this->preTrapezoidalCorrector(m_rhoDaggerParticles,
                                a_dt); // Mobile or diffusive moves to X^dagger = X^k + 0.5*dt*V^k + hop
  this->remapPointParticles(m_rhoDaggerParticles,
                            SpeciesSubset::
                              AllMobileOrDiffusive); // Only need to remap particles that were mobile or diffusive
  this->depositPointParticles(m_rhoDaggerParticles,
                              SpeciesSubset::All); // Everything needs to deposit...

  this->copyConductivityParticles(m_conductivityParticles);
  this->computeConductivities(m_conductivityParticles);
  this->setupSemiImplicitPoisson(0.5 * a_dt);
  this->solvePoisson();

  this->setItoVelocityFunctions();
  m_ito->interpolateVelocities();
  this->trapezoidalCorrector(a_dt);
  this->remapParticles(SpeciesSubset::AllMobileOrDiffusive);
  // ====== CORRECTOR END =====

  // Do particle-boundary intersection.
  this->intersectParticles(SpeciesSubset::AllMobileOrDiffusive, EBIntersection::Bisection, true);
  this->removeCoveredParticles(SpeciesSubset::AllMobileOrDiffusive, EBRepresentation::ImplicitFunction, m_toleranceEB);

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

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    const Real f = mobile ? 1.0 : 0.0;    // Multiplication factor for mobility
    const Real g = diffusive ? 1.0 : 0.0; // Multiplication factor for diffusion

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_particleRealm)[lvl];
      ParticleData<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

        List<ItoParticle>&   ItoParticles   = particles[dit()].listItems();
        List<PointParticle>& pointParticles = (*a_rho_dagger[idx])[lvl][dit()].listItems();

        pointParticles.clear();

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
          pointParticles.add(PointParticle(Xk + g * hop, weight));
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

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>& solver = solverIt();

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

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>&        solver  = solverIt();
    const RefCountedPtr<ItoSpecies>& species = solver->getSpecies();

    const int idx = solverIt.index();

    const bool mobile    = solver->isMobile();
    const bool diffusive = solver->isDiffusive();

    const Real f = mobile ? 1.0 : 0.0;    // Multiplication factor for mobility
    const Real g = diffusive ? 1.0 : 0.0; // Multiplication factor for diffusion

    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout&   dbl       = m_amr->getGrids(m_particleRealm)[lvl];
      ParticleData<ItoParticle>& particles = solver->getParticles(ItoSolver::WhichContainer::Bulk)[lvl];

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

        List<ItoParticle>&   ItoParticles   = particles[dit()].listItems();
        List<PointParticle>& pointParticles = (*a_rho_dagger[idx])[lvl][dit()].listItems();

        pointParticles.clear();

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
          pointParticles.add(PointParticle(pos, weight));
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

  for (auto solverIt = m_ito->iterator(); solverIt.ok(); ++solverIt) {
    RefCountedPtr<ItoSolver>& solver = solverIt();

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
