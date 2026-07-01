/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_RadiativeTransferStepperImplem.H
  @brief  Implementation of CD_RadiativeTransferStepper.H
  @author Robert Marskar
*/

// Std includes
#include <math.h>
#include <type_traits>

// Chombo includes
#include <CH_Timer.H>
#include <DataIterator.H>
#include <BoxIterator.H>
#include <ParmParse.H>
#include <PolyGeom.H>

// Our includes
#include <CD_RadiativeTransferStepper.H>
#include <CD_RadiativeTransferSpecies.H>
#include <CD_AmrMesh.H>
#include <CD_McPhoto.H>
#include <CD_Photon.H>
#include <CD_ParticleContainer.H>
#include <CD_ParticleSoA.H>
#include <CD_Random.H>
#include <CD_NamespaceHeader.H>

using namespace Physics::RadiativeTransfer;

template <typename T>
RadiativeTransferStepper<T>::RadiativeTransferStepper()
{
  CH_TIME("RadiativeTransferStepper::RadiativeTransferStepper");

  // Get shit for setting the source term
  ParmParse    pp("RadiativeTransferStepper");
  Vector<Real> v;

  pp.get("verbosity", m_verbosity);
  pp.get("realm", m_realm);
  pp.get("dt", m_dt);

  m_forceDt = -1.0;

  // Optional ghost-particle boundary verification (McPhoto only). Off unless explicitly requested.
  m_ghostBoundaryTest          = false;
  m_ghostBoundaryTestParticles = 100000;
  m_ghostBoundaryTestWidth     = 1;
  pp.query("ghost_boundary_test", m_ghostBoundaryTest);
  pp.query("ghost_boundary_test_particles", m_ghostBoundaryTestParticles);
  pp.query("ghost_boundary_test_width", m_ghostBoundaryTestWidth);
}

template <typename T>
RadiativeTransferStepper<T>::~RadiativeTransferStepper()
{
  CH_TIME("RadiativeTransferStepper::~RadiativeTransferStepper");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::~RadiativeTransferStepper" << endl;
  }
}

template <typename T>
void
RadiativeTransferStepper<T>::setupSolvers()
{
  CH_TIME("RadiativeTransferStepper::setupSolvers");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::setupSolvers" << endl;
  }

  m_species = RefCountedPtr<RtSpecies>(new RadiativeTransferSpecies());

  // Solver setup
  m_solver = RefCountedPtr<RtSolver>(new T());
  m_solver->setVerbosity(m_verbosity);
  m_solver->setRtSpecies(m_species);
  m_solver->parseOptions();
  m_solver->setPhase(phase::gas);
  m_solver->setAmr(m_amr);
  m_solver->setComputationalGeometry(m_computationalGeometry);
  m_solver->setRealm(m_realm);
}

template <typename T>
void
RadiativeTransferStepper<T>::parseRuntimeOptions()
{
  CH_TIME("RadiativeTransferStepper::parseRuntimeOptions");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::parseRuntimeOptions" << endl;
  }

  ParmParse pp("RadiativeTransferStepper");

  pp.get("verbosity", m_verbosity);
  pp.get("dt", m_dt);

  m_solver->parseRuntimeOptions();
}

template <typename T>
void
RadiativeTransferStepper<T>::allocate()
{
  CH_TIME("RadiativeTransferStepper::allocate");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::allocate" << endl;
  }

  m_solver->allocate();
}

template <typename T>
void
RadiativeTransferStepper<T>::initialData()
{
  CH_TIME("RadiativeTransferStepper::initialData");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::initialData" << endl;
  }

  m_solver->initialData();
  this->setGaussianSource();

  if (m_solver->isStationary()) {
    m_solver->advance(0.0, false);
  }
}

template <typename T>
void
RadiativeTransferStepper<T>::postInitialize()
{
  CH_TIME("RadiativeTransferStepper::postInitialize");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::postInitialize" << endl;
  }

  if (m_ghostBoundaryTest) {
    this->ghostBoundaryTest();
  }
}

template <typename T>
void
RadiativeTransferStepper<T>::ghostBoundaryTest()
{
  CH_TIME("RadiativeTransferStepper::ghostBoundaryTest");

  // This verification needs the McPhoto particle container; it is a no-op for mesh-based solvers.
  if constexpr (std::is_same_v<T, McPhoto>) {
    ParticleContainer<Photon>& photons = static_cast<T*>(&(*m_solver))->getPhotons();

    // (1) Sample particles uniformly in the domain box (on rank 0; remap distributes them to owners).
    const RealVect probLo = m_amr->getProbLo();
    const RealVect probHi = m_amr->getProbHi();

    photons.clearParticles();

    ParticleSoA<Photon> sampled;
    if (procID() == 0) {
      for (int i = 0; i < m_ghostBoundaryTestParticles; i++) {
        sampled.append(Random::randomPosition(probLo, probHi), 1.0);
      }
    }
    photons.addParticlesDestructive(sampled); // routes every particle to its owning patch/level/rank

    // (2) Fill the ghost layer of every patch from the prebuilt Realm masks (registered width).
    const AMRParticleGhostMask& same = m_amr->getParticleGhostMask(m_realm, m_ghostBoundaryTestWidth);
    const AMRParticleGhostMask& c2f  = m_amr->getParticleGhostMaskCoarToFine(m_realm, m_ghostBoundaryTestWidth);
    const AMRParticleGhostMask& f2c  = m_amr->getParticleGhostMaskFineToCoar(m_realm, m_ghostBoundaryTestWidth);

    photons.fillGhost(same, c2f, f2c);

    // (3) Discard the valid particles, keeping only the ghosts, and turn the ghosts into ordinary
    //     (valid) particles so the subsequent remap() will route them rather than drop them.
    long long numGhost = 0;
    for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++) {
      const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[lvl];
      const DataIterator&      dit  = dbl.dataIterator();
      const int                nbox = dit.size();

      for (int mybox = 0; mybox < nbox; mybox++) {
        ParticleSoA<Photon>& leaf = photons[lvl][dit[mybox]];

        std::size_t i = 0;
        while (i < leaf.size()) {
          if (leaf.isGhost(i)) {
            leaf.ghost(i) = GhostType::Valid; // a former ghost becomes an ordinary particle
            numGhost++;
            i++;
          }
          else {
            leaf.remove(i); // discard the patch's own valid particle (swap-pop; do not advance)
          }
        }
      }
    }

    // ParallelOps::sum RETURNS the global sum (it does not modify its argument in place).
    const long long numGhostGlobal = ParallelOps::sum(numGhost);
    if (procID() == 0) {
      pout() << "RadiativeTransferStepper::ghostBoundaryTest - ghost copies filled = " << numGhostGlobal << endl;
    }

    // (4) Remap the former ghosts back to their owners. What remains are exactly the particles that lay
    //     within a ghost layer of some patch boundary, duplicated once per patch that ghosted them.
    photons.remap();

    // getNumberOfValidParticlesGlobal() is collective, so ALL ranks must call it (not just rank 0).
    const long long numBoundary = photons.getNumberOfValidParticlesGlobal();
    if (procID() == 0) {
      pout() << "RadiativeTransferStepper::ghostBoundaryTest - boundary particles after remap = " << numBoundary
             << endl;
    }
  }
  else {
    MayDay::Warning("RadiativeTransferStepper::ghostBoundaryTest - only implemented for the McPhoto solver");
  }
}

template <typename T>
void
RadiativeTransferStepper<T>::forceDt(const Real a_dt)
{
  CH_TIME("RadiativeTransferStepper::forceDt");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::forceDt" << endl;
  }

  m_forceDt = a_dt;
}

#ifdef CH_USE_HDF5
template <typename T>
void
RadiativeTransferStepper<T>::writeCheckpointData(HDF5Handle& a_handle, const int a_lvl) const
{
  CH_TIME("RadiativeTransferStepper::writeCheckpointData");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::writeCheckpointData" << endl;
  }

  m_solver->writeCheckpointLevel(a_handle, a_lvl);
}
#endif

#ifdef CH_USE_HDF5
template <typename T>
void
RadiativeTransferStepper<T>::readCheckpointData(HDF5Handle& a_handle, const int a_lvl)
{
  CH_TIME("RadiativeTransferStepper::readCheckpointData");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::readCheckpointData" << endl;
  }

  m_solver->readCheckpointLevel(a_handle, a_lvl);
}
#endif

template <typename T>
void
RadiativeTransferStepper<T>::postCheckpointSetup()
{
  CH_TIME("RadiativeTransferStepper::postCheckpointSetup");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::postCheckpointSetup" << endl;
  }

  this->setGaussianSource();
}

template <typename T>
int
RadiativeTransferStepper<T>::getNumberOfPlotVariables() const
{
  CH_TIME("RadiativeTransferStepper::getNumberOfPlotVariables");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::getNumberOfPlotVariables" << endl;
  }

  // +3 dev fields (ghost_same/ghost_fine_to_coar/ghost_coar_to_fine) visualizing the Realm ghost masks.
  return m_solver->getNumberOfPlotVariables() + (m_ghostBoundaryTest ? 3 : 0);
}

template <typename T>
Vector<std::string>
RadiativeTransferStepper<T>::getPlotVariableNames() const
{
  CH_TIME("RadiativeTransferStepper::getPlotVariableNames");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::getPlotVariableNames" << endl;
  }

  Vector<std::string> names = m_solver->getPlotVariableNames();
  if (m_ghostBoundaryTest) {
    names.push_back("ghost_same");
    names.push_back("ghost_fine_to_coar");
    names.push_back("ghost_coar_to_fine");
  }
  return names;
}

template <typename T>
void
RadiativeTransferStepper<T>::writePlotData(LevelData<EBCellFAB>& a_output,
                                           int&                  a_icomp,
                                           const std::string&    a_outputRealm,
                                           const int             a_level) const
{
  CH_TIME("RadiativeTransferStepper::writePlotData");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::writePlotData" << endl;
  }

  m_solver->writePlotData(a_output, a_icomp, a_outputRealm, a_level);

  if (m_ghostBoundaryTest) {
    this->writeGhostTargetPlotData(a_output, a_icomp, a_outputRealm, a_level);
  }
}

template <typename T>
void
RadiativeTransferStepper<T>::writeGhostTargetPlotData(LevelData<EBCellFAB>& a_output,
                                                      int&                  a_icomp,
                                                      const std::string&    a_outputRealm,
                                                      const int             a_level) const
{
  // DEV visualization: write 1 in each cell that will SCATTER its particles (has >=1 target) and 0
  // elsewhere, one field per Realm ghost mask (same / fine-to-coarse / coarse-to-fine). Assumes the
  // output realm matches the stepper realm (true for the McPhoto example).
  CH_assert(a_outputRealm == m_realm);

  const AMRParticleGhostMask& same = m_amr->getParticleGhostMask(m_realm, m_ghostBoundaryTestWidth);
  const AMRParticleGhostMask& coar = m_amr->getParticleGhostMaskFineToCoar(m_realm, m_ghostBoundaryTestWidth);
  const AMRParticleGhostMask& fine = m_amr->getParticleGhostMaskCoarToFine(m_realm, m_ghostBoundaryTestWidth);

  const DisjointBoxLayout& dbl  = m_amr->getGrids(m_realm)[a_level];
  const DataIterator&      dit  = dbl.dataIterator();
  const int                nbox = dit.size();

  const AMRParticleGhostMask* tables[3] = {&same, &coar, &fine};

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din      = dit[mybox];
    const Box        validBox = dbl[din];

    EBCellFAB&     out = a_output[din];
    BaseFab<Real>& reg = out.getSingleValuedFAB();

    for (int k = 0; k < 3; k++) {
      const int comp = a_icomp + k;

      out.setVal(comp, 0.0);

      const AMRParticleGhostMask& table = *tables[k];
      if (a_level < table.size() && !table[a_level].isNull()) {
        const ParticleGhostMask& t = (*table[a_level])[din];
        for (BoxIterator bit(validBox & t.box()); bit.ok(); ++bit) {
          if (t.numTargets(bit()) > 0) {
            reg(bit(), comp) = 1.0;
          }
        }
      }
    }
  }

  a_icomp += 3;
}

template <typename T>
Real
RadiativeTransferStepper<T>::computeDt()
{
  CH_TIME("RadiativeTransferStepper::computeDt");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::computeDt" << endl;
  }

  Real dt = m_dt;

  if (m_forceDt > 0.0) {
    dt = m_forceDt;
  }

  return dt;
}

template <typename T>
Real
RadiativeTransferStepper<T>::advance(const Real a_dt)
{
  CH_TIME("RadiativeTransferStepper::advance");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::advance" << endl;
  }

  m_solver->advance(a_dt);

  return a_dt;
}

template <typename T>
void
RadiativeTransferStepper<T>::synchronizeSolverTimes(const int a_step, const Real a_time, const Real a_dt)
{
  CH_TIME("RadiativeTransferStepper::synchronizeSolverTimes");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::synchronizeSolverTimes" << endl;
  }

  m_timeStep = a_step;
  m_time     = a_time;
  m_dt       = a_dt;

  m_solver->setTime(m_timeStep, m_time, m_dt);
}

template <typename T>
void
RadiativeTransferStepper<T>::printStepReport()
{
  CH_TIME("RadiativeTransferStepper::printStepReport");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::printStepReport" << endl;
  }
}

template <typename T>
void
RadiativeTransferStepper<T>::registerRealms()
{
  CH_TIME("RadiativeTransferStepper::registerRealms");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::registerRealms" << endl;
  }

  m_amr->registerRealm(m_realm);

  // DEV: register the ghost-mask width we visualize (see writeGhostTargetPlotData).
  if (m_ghostBoundaryTest) {
    m_amr->registerParticleGhostMask(m_realm, m_ghostBoundaryTestWidth);
  }
}

template <typename T>
void
RadiativeTransferStepper<T>::registerOperators()
{
  CH_TIME("RadiativeTransferStepper::registerOperators");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::registerOperators" << endl;
  }

  m_solver->registerOperators();
}

template <typename T>
void
RadiativeTransferStepper<T>::preRegrid(const int a_base, const int a_oldFinestLevel)
{
  CH_TIME("RadiativeTransferStepper::preRegrid");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::preRegrid" << endl;
  }

  m_solver->preRegrid(a_base, a_oldFinestLevel);
}

template <typename T>
void
RadiativeTransferStepper<T>::regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel)
{
  CH_TIME("RadiativeTransferStepper::regrid");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::regrid" << endl;
  }

  m_solver->regrid(a_lmin, a_oldFinestLevel, a_newFinestLevel);
}

template <typename T>
void
RadiativeTransferStepper<T>::postRegrid()
{
  CH_TIME("RadiativeTransferStepper::postRegrid");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::postRegrid" << endl;
  }

  this->setGaussianSource();
}

template <typename T>
const EBAMRCellData&
RadiativeTransferStepper<T>::getPhi() const
{
  CH_TIME("RadiativeTransferStepper::getPhi");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::getPhi" << endl;
  }

  return m_solver->getPhi();
}

template <typename T>
void
RadiativeTransferStepper<T>::setGaussianSource()
{
  CH_TIME("RadiativeTransferStepper::setGaussianSource");
  if (m_verbosity > 5) {
    pout() << "RadiativeTransferStepper::setGaussianSource" << endl;
  }

  ParmParse pp("RadiativeTransferStepper");

  Vector<Real> v(SpaceDim);
  RealVect     blobCenter;
  Real         blobAmp;
  Real         blobRad;

  pp.get("blob_amplitude", blobAmp);
  pp.get("blob_radius", blobRad);
  pp.getarr("blob_center", v, 0, SpaceDim);
  blobCenter = RealVect(D_DECL(v[0], v[1], v[2]));

  auto sourceFunction = [blobAmp, blobRad, blobCenter](const RealVect a_position) -> Real {
    const RealVect dist  = a_position - blobCenter;
    const Real     dist2 = dist.dotProduct(dist);
    return blobAmp * exp(-dist2 / (2 * blobRad * blobRad));
  };

  m_solver->setSource(sourceFunction);
}

#include <CD_NamespaceFooter.H>
