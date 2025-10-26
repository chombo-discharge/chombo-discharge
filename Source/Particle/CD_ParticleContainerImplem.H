/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ParticleContainerImplem.H
  @brief  Implementation of CD_ParticleContainer.H
  @author Robert Marskar
*/

#ifndef CD_ParticleContainerImplem_H
#define CD_ParticleContainerImplem_H

#warning "The particle map we send should probably be ParticleMap<std::vector<P>> rather than List<P>"

// Std includes
#include <bitset>

// Chombo includes
#include <ParmParse.H>
#include <CH_Timer.H>
#include <MPI_util.H>
#include <ParticleDataI.H>

// Our includes
#include <CD_ParticleContainer.H>
#include <CD_ParallelOps.H>
#include <CD_ParticleOps.H>
#include <CD_Timer.H>
#include <CD_BoxLoops.H>
#include <CD_OpenMP.H>
#include <CD_NamespaceHeader.H>

template <class P>
ParticleContainer<P>::ParticleContainer()
{
  m_isDefined         = false;
  m_isOrganizedByCell = false;
  m_profile           = false;
  m_debug             = false;
  m_verbose           = false;
}

template <class P>
ParticleContainer<P>::ParticleContainer(const Vector<DisjointBoxLayout>&         a_grids,
                                        const Vector<ProblemDomain>&             a_domains,
                                        const Vector<Real>&                      a_dx,
                                        const Vector<int>&                       a_refRat,
                                        const Vector<ValidMask>&                 a_validMask,
                                        const Vector<RefCountedPtr<LevelTiles>>& a_levelTiles,
                                        const RealVect&                          a_probLo,
                                        const int                                a_blockingFactor,
                                        const int                                a_finestLevel,
                                        const std::string                        a_realm)
{
  CH_TIME("ParticleContainer::ParticleContainer");

  this->define(a_grids,
               a_domains,
               a_dx,
               a_refRat,
               a_validMask,
               a_levelTiles,
               a_probLo,
               a_blockingFactor,
               a_finestLevel,
               a_realm);
}

template <class P>
ParticleContainer<P>::~ParticleContainer()
{
  CH_TIME("ParticleContainer::~ParticleContainer");
}

template <class P>
void
ParticleContainer<P>::define(const Vector<DisjointBoxLayout>&         a_grids,
                             const Vector<ProblemDomain>&             a_domains,
                             const Vector<Real>&                      a_dx,
                             const Vector<int>&                       a_refRat,
                             const Vector<ValidMask>&                 a_validMask,
                             const Vector<RefCountedPtr<LevelTiles>>& a_levelTiles,
                             const RealVect&                          a_probLo,
                             const int                                a_blockingFactor,
                             const int                                a_finestLevel,
                             const std::string                        a_realm)
{
  CH_TIME("ParticleContainer::~ParticleContainer");

  m_grids          = a_grids;
  m_domains        = a_domains;
  m_refRat         = a_refRat;
  m_validRegion    = a_validMask;
  m_probLo         = a_probLo;
  m_finestLevel    = a_finestLevel;
  m_realm          = a_realm;
  m_blockingFactor = a_blockingFactor;
  m_levelTiles     = a_levelTiles;

  m_dx.resize(1 + m_finestLevel);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    m_dx[lvl] = a_dx[lvl] * RealVect::Unit;
  }

  constexpr int base = 0;

  // Do the define stuff.
  this->setupGrownGrids(base, m_finestLevel);
  this->setupParticleData(base, m_finestLevel);

  m_isDefined         = true;
  m_isOrganizedByCell = false;
  m_profile           = false;

  ParmParse pp("ParticleContainer");
  pp.query("profile", m_profile);
  pp.query("debug", m_debug);
  pp.query("verbose", m_verbose);
}

template <class P>
void
ParticleContainer<P>::setupGrownGrids(const int a_base, const int a_finestLevel)
{
  CH_TIME("ParticleContainer::setupGrownGrids");
  if (m_verbose) {
    pout() << "ParticleContainer::setupGrownGrids" << endl;
  }

  // TLDR: This routine sets up the buffer grids which can be used when we need to fetch particles that fall off each levels' grid. This is very useful
  //       when we want to fetch particles that lie on the coarse grid of a refinement boundary.

  m_grownGrids.resize(1 + a_finestLevel);

  for (int lvl = 0; lvl <= a_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl    = m_grids[lvl];
    const ProblemDomain&     domain = m_domains[lvl];

    // Grow boxes by refinement factor on the finer levels.
    Vector<Box> boxes = dbl.boxArray();
    if (lvl > 0) {
      for (auto& box : boxes.stdVector()) {
        box.grow(m_refRat[lvl - 1]);
        box &= domain;
      }
    }

    m_grownGrids[lvl] = BoxLayout(boxes, dbl.procIDs());
  }
}

template <class P>
void
ParticleContainer<P>::setupParticleData(const int a_base, const int a_finestLevel)
{
  CH_TIME("ParticleContainer::setupParticleData");
  if (m_verbose) {
    pout() << "ParticleContainer::setupParticleData" << endl;
  }

  // TLDR: This sets up the most commonly used particle data holders for this particle AMR container. This means that
  //       we allocate:
  //
  //          1. The "normal" particle container data holder m_particles. This is defined on the DisjointBoxLayout
  //             which is the natural place for the particles to live.
  //
  //          2. A buffer particle data holder which is defined on a grown DisjointBoxLayout. This is useful when particles
  //             on the coarse level need to deposit to the fine level.
  //
  //          3. A data holder for "masked particles", providing an opportunity to copy/transfer some of the particles on a grid
  //             level to a separate data holder if they lie within a "mask". Typically used for extracting coarse-level that live
  //             just outside the fine grid (i.e. on the coarse side of the refinement boundary).
  //
  //          5. A data holder for storing cell-sorted particles. Very useful when particles need to be sorted by cell rather than patch.

  m_particles.resize(1 + a_finestLevel);
  m_bufferParticles.resize(1 + a_finestLevel);
  m_maskParticles.resize(1 + a_finestLevel);
  m_cellSortedParticles.resize(1 + a_finestLevel);

  for (int lvl = a_base; lvl <= a_finestLevel; lvl++) {
    m_particles[lvl] = RefCountedPtr<ParticleData<P>>(
      new ParticleData<P>(m_grids[lvl], m_domains[lvl], m_blockingFactor, m_dx[lvl], m_probLo));

    m_bufferParticles[lvl] = RefCountedPtr<ParticleData<P>>(
      new ParticleData<P>(m_grownGrids[lvl], m_domains[lvl], m_blockingFactor, m_dx[lvl], m_probLo));

    m_maskParticles[lvl] = RefCountedPtr<ParticleData<P>>(
      new ParticleData<P>(m_grids[lvl], m_domains[lvl], m_blockingFactor, m_dx[lvl], m_probLo));

    m_cellSortedParticles[lvl] = RefCountedPtr<LayoutData<BinFab<P>>>(new LayoutData<BinFab<P>>(m_grids[lvl]));
  }
}

template <class P>
void
ParticleContainer<P>::sortParticles() noexcept
{
  CH_TIME("ParticleContainer::sortParticles");
  if (m_verbose) {
    pout() << "ParticleContainer::sortParticles" << endl;
  }

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {

    const DisjointBoxLayout& dbl = m_grids[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      List<P>& particles = (*m_particles[lvl])[din].listItems();

      particles.sort();
    }
  }
}

template <class P>
bool
ParticleContainer<P>::isOrganizedByCell() const
{
  return m_isOrganizedByCell;
}

template <class P>
int
ParticleContainer<P>::getFinestLevel() const
{
  CH_assert(m_isDefined);

  return m_finestLevel;
}

template <class P>
const std::string
ParticleContainer<P>::getRealm() const
{
  CH_assert(m_isDefined);

  return m_realm;
}

template <class P>
const Vector<DisjointBoxLayout>&
ParticleContainer<P>::getGrids() const
{
  return m_grids;
}

template <class P>
const RealVect
ParticleContainer<P>::getProbLo() const noexcept
{
  return m_probLo;
}

template <class P>
const Vector<RealVect>
ParticleContainer<P>::getDx() const noexcept
{
  return m_dx;
}

template <class P>
AMRParticles<P>&
ParticleContainer<P>::getParticles()
{
  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer::getParticles - particles are sorted by cell!");
  }

  return m_particles;
}

template <class P>
const AMRParticles<P>&
ParticleContainer<P>::getParticles() const
{
  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {
    MayDay::Abort("ParticleContainer::getParticles - particles are sorted by cell!");
  }

  return m_particles;
}

template <class P>
AMRParticles<P>&
ParticleContainer<P>::getBufferParticles()
{
  CH_assert(m_isDefined);

  return m_bufferParticles;
}

template <class P>
const AMRParticles<P>&
ParticleContainer<P>::getBufferParticles() const
{
  CH_assert(m_isDefined);

  return m_bufferParticles;
}

template <class P>
AMRParticles<P>&
ParticleContainer<P>::getMaskParticles()
{
  CH_assert(m_isDefined);

  return m_maskParticles;
}

template <class P>
const AMRParticles<P>&
ParticleContainer<P>::getMaskParticles() const
{
  CH_assert(m_isDefined);

  return m_maskParticles;
}

template <class P>
ParticleData<P>&
ParticleContainer<P>::operator[](const int a_lvl)
{
  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer::operator[](a_lvl) - particles are sorted by cell!");
  }

  return *m_particles[a_lvl];
}

template <class P>
const ParticleData<P>&
ParticleContainer<P>::operator[](const int a_level) const
{
  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer::operator[](a_lvl) - particles are sorted by cell!");
  }

  return *m_particles[a_level];
}

template <class P>
AMRCellParticles<P>&
ParticleContainer<P>::getCellParticles()
{
  CH_assert(m_isDefined);

  if (!m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer::getCellParticles()- particles are not sorted by cell!");
  }

  return m_cellSortedParticles;
}

template <class P>
const AMRCellParticles<P>&
ParticleContainer<P>::getCellParticles() const
{
  CH_assert(m_isDefined);

  if (!m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer::getCellParticles()- particles are not sorted by cell!");
  }

  return m_cellSortedParticles;
}

template <class P>
LayoutData<BinFab<P>>&
ParticleContainer<P>::getCellParticles(const int a_level)
{
  CH_assert(m_isDefined);

  if (!m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer::getCellParticles(level)- particles are not sorted by cell!");
  }

  return *m_cellSortedParticles[a_level];
}

template <class P>
const LayoutData<BinFab<P>>&
ParticleContainer<P>::getCellParticles(const int a_level) const
{
  CH_TIME("ParticleContainer::getCellParticles(int)");

  CH_assert(m_isDefined);

  if (!m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer::getCellParticles(level)- particles are not sorted by cell!");
  }

  return *m_cellSortedParticles[a_level];
}

template <class P>
void
ParticleContainer<P>::getCellParticles(BinFab<P>& cellParticles, const int a_lvl, const DataIndex a_dit) const
{
  CH_TIME("ParticleContainer::getCellParticles(BinFab)");

  CH_assert(m_isDefined);

  if (!m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer::getCellParticles - particles are not sorted by cell!");
  }

  cellParticles.define(m_grids[a_lvl][a_dit], m_dx[a_lvl], m_probLo);
  cellParticles.addItems((*m_particles[a_lvl])[a_dit].listItems());
}

template <class P>
void
ParticleContainer<P>::getCellParticlesDestructive(BinFab<P>& cellParticles, const int a_lvl, const DataIndex a_dit)
{
  CH_TIME("ParticleContainer::getCellParticlesDestructive");

  CH_assert(m_isDefined);

  if (!m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer::getCellParticlesDestructive - particles are not sorted by cell!");
  }

  cellParticles.define(m_grids[a_lvl].get(a_dit), m_dx[a_lvl], m_probLo);
  cellParticles.addItemsDestructive((*m_particles[a_lvl])[a_dit].listItems());
}

template <class P>
BinFab<P>&
ParticleContainer<P>::getCellParticles(const int a_level, const DataIndex a_dit)
{
  CH_TIME("ParticleContainer::getCellParticles(int, DataIndex)");

  CH_assert(m_isDefined);

  if (!m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer::getCellParticles(int, dit) - particles are not sorted by cell!");
  }

  return (*m_cellSortedParticles[a_level])[a_dit];
}

template <class P>
const BinFab<P>&
ParticleContainer<P>::getCellParticles(const int a_level, const DataIndex a_dit) const
{
  CH_TIME("ParticleContainer::getCellParticles(int, DataIndex)");

  CH_assert(m_isDefined);

  if (!m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer::getCellParticles(int, dit) - particles are not sorted by cell!");
  }

  return (*m_cellSortedParticles[a_level])[a_dit];
}

template <class P>
void
ParticleContainer<P>::organizeParticlesByCell()
{
  CH_TIME("ParticleContainer::organizeParticlesByCell");
  if (m_verbose) {
    pout() << "ParticleContainer::organizeParticlesByCell" << endl;
  }

  CH_assert(m_isDefined);

  if (!m_isOrganizedByCell) {

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const DisjointBoxLayout& dbl = m_grids[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        BinFab<P>& cellParticles = (*m_cellSortedParticles[lvl])[din];

        cellParticles.define(dbl[din], m_dx[lvl], m_probLo);
        cellParticles.addItemsDestructive((*m_particles[lvl])[din].listItems());
      }
    }

    m_isOrganizedByCell = true;
  }
}

template <class P>
void
ParticleContainer<P>::organizeParticlesByPatch()
{
  CH_TIME("ParticleContainer::organizeParticlesByPatch");
  if (m_verbose) {
    pout() << "ParticleContainer::organizeParticlesByPatch" << endl;
  }

  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {

    constexpr int comp = 0;

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const DisjointBoxLayout& dbl = m_grids[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        ListBox<P>& patchParticles = (*m_particles[lvl])[din];
        BinFab<P>&  cellParticles  = (*m_cellSortedParticles[lvl])[din];

        // Kernel which adds moves particles from the cell container to the patch container.
        auto kernel = [&](const IntVect& iv) -> void {
          patchParticles.addItemsDestructive(cellParticles(iv, comp));
        };

        BoxLoops::loop(dbl[din], kernel);

        cellParticles.clear();
      }
    }

    m_isOrganizedByCell = false;
  }
}

template <class P>
void
ParticleContainer<P>::addParticles(const List<P>& a_particles)
{
  CH_TIME("ParticleContainer::addParticles");
  if (m_verbose) {
    pout() << "ParticleContainer::addParticles" << endl;
  }

  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer::addParticles(List<P>) - particles are sorted by cell!");
  }

  // TLDR: This routine is almost a precise copy of the remap routine, with the exception that the particles to be remapped
  //       are only the input particles rather than everything inside m_particles

  const uint32_t numRanks = numProc();
  const uint32_t myRank   = procID();

  using LevelAndIndex = std::pair<uint32_t, uint32_t>;

  std::vector<ParticleMap<List<P>>> particlesToSend(numRanks);

#pragma omp parallel
  {
    List<P> outcasts;

#pragma omp master
    {
      outcasts.join(a_particles);
    }

    std::vector<ParticleMap<List<P>>> threadLocalParticlesToSend(numRanks);

    // Collect particles owned by this rank/thread.
    this->transferParticlesToSingleList(outcasts, m_particles);

    // Map the particles to other grid patches
    this->mapParticlesToAMRGrid(threadLocalParticlesToSend, outcasts);

    // Catenate the thread-local particle lists onto particlesToSend
#pragma omp critical
    this->catenateParticleMaps(particlesToSend, threadLocalParticlesToSend);
  }

  // Particles that go from this rank to this rank don't need to be communicated so that we can place them directly on the correct patch.
  this->assignLocalParticles(particlesToSend[myRank], m_particles);

  // If using MPI then we have to scatter the particles across MPI ranks.
#ifdef CH_MPI
  ParticleMap<List<P>> receivedParticles;

  ParticleOps::scatterParticles(receivedParticles, particlesToSend);

  // Assign particles to the correct level and grid patch -- we iterate through receivedParticles and decode the information
  // we got from there.
  this->assignLocalParticles(receivedParticles, m_particles);
#endif

  if (m_debug) {
    this->sanityCheck();
  }
}

template <class P>
void
ParticleContainer<P>::addParticlesDestructive(List<P>& a_particles)
{
  CH_TIME("ParticleContainer::addParticlesDestructive");
  if (m_verbose) {
    pout() << "ParticleContainer::addParticlesDestructive" << endl;
  }

  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer::addParticlesDestructive(List<P>) - particles are sorted by cell!");
  }

  // TLDR: This routine is almost a precise copy of the remap routine, with the exception that the particles to be remapped
  //       are only the input particles rather than everything inside m_particles

  const uint32_t numRanks = numProc();
  const uint32_t myRank   = procID();

  using LevelAndIndex = std::pair<uint32_t, uint32_t>;

  std::vector<ParticleMap<List<P>>> particlesToSend(numRanks);

#pragma omp parallel
  {
    List<P> outcasts;

#pragma omp master
    {
      outcasts.catenate(a_particles);
    }

    std::vector<ParticleMap<List<P>>> threadLocalParticlesToSend(numRanks);

    // Collect particles owned by this rank/thread.
    this->transferParticlesToSingleList(outcasts, m_particles);

    // Map the particles to other grid patches
    this->mapParticlesToAMRGrid(threadLocalParticlesToSend, outcasts);

    // Catenate the thread-local particle lists onto particlesToSend
#pragma omp critical
    this->catenateParticleMaps(particlesToSend, threadLocalParticlesToSend);
  }

  // Particles that go from this rank to this rank don't need to be communicated so that we can place them directly on the correct patch.
  this->assignLocalParticles(particlesToSend[myRank], m_particles);

  // If using MPI then we have to scatter the particles across MPI ranks.
#ifdef CH_MPI
  ParticleMap<List<P>> receivedParticles;

  ParticleOps::scatterParticles(receivedParticles, particlesToSend);

  // Assign particles to the correct level and grid patch -- we iterate through receivedParticles and decode the information
  // we got from there.
  this->assignLocalParticles(receivedParticles, m_particles);
#endif

  if (m_debug) {
    this->sanityCheck();
  }
}

template <class P>
void
ParticleContainer<P>::addParticles(const BinFab<P>& a_particles, const int a_lvl, const DataIndex a_dit)
{
  CH_TIME("ParticleContainer::addParticles(BinFab)");
  if (m_verbose) {
    pout() << "ParticleContainer::addParticles(BinFab)" << endl;
  }

  CH_assert(m_isDefined);
  CH_assert(m_grids[a_lvl].get(a_dit) == a_particles.getRegion());

  if (!m_isOrganizedByCell) {
    MayDay::Abort("ParticleContainer::addParticles(BinFab<P>) - particles are not sorted by cell!");
  }

  constexpr int comp = 0;

  BinFab<P>& boxParticles = (*m_cellSortedParticles[a_lvl])[a_dit];

  auto kernel = [&](const IntVect& iv) -> void {
    List<P>&       myParticles    = boxParticles(iv, comp);
    const List<P>& inputParticles = a_particles(iv, comp);

    myParticles.join(inputParticles);
  };

  BoxLoops::loop(m_grids[a_lvl][a_dit], kernel);
}

template <class P>
void
ParticleContainer<P>::addParticlesDestructive(BinFab<P>& a_particles, const int a_lvl, const DataIndex a_dit)
{
  CH_TIME("ParticleContainer::addParticlesDestructive(BinFab)");
  if (m_verbose) {
    pout() << "ParticleContainer::addParticlesDestructive(BinFab)" << endl;
  }

  CH_assert(m_isDefined);
  CH_assert(m_grids[a_lvl].get(a_dit) == a_particles.getRegion());

  if (!m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer::addParticles(BinFab<P>) - particles are not sorted by cell!");
  }

  constexpr int comp = 0;

  BinFab<P>& boxParticles = (*m_cellSortedParticles[a_lvl])[a_dit];

  auto kernel = [&](const IntVect& iv) -> void {
    List<P>&       myParticles    = boxParticles(iv, comp);
    const List<P>& inputParticles = a_particles(iv, comp);

    myParticles.catenate(inputParticles);
  };

  BoxLoops::loop(m_grids[a_lvl][a_dit], kernel);
}

template <class P>
void
ParticleContainer<P>::addParticles(const ParticleContainer<P>& a_otherContainer)
{
  CH_TIME("ParticleContainer::addParticles(ParticleContainer)");
  if (m_verbose) {
    pout() << "ParticleContainer::addParticles(ParticleContainer)" << endl;
  }

  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer::addParticles(ParticleContainer<P>) - particles are sorted by cell!");
  }

  // TLDR: This routine adds particles from a different container to this container. If the containers live on the same realm
  //       we just add the particles directly. Otherwise, we have to add to the outcast list and remap.

  if (m_realm == a_otherContainer.getRealm()) {
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const DisjointBoxLayout& dbl = m_grids[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        ListBox<P>&       myParticles    = (*m_particles[lvl])[din];
        const ListBox<P>& otherParticles = a_otherContainer[lvl][din];

        myParticles.addItems(otherParticles.listItems());
      }
    }
  }
  else {
    // If we're adding particles across realms we need to collect the input particles onto lists that get remapped
    // into m_particles. A full explanation of how this works is given in the remap routin.

    const uint32_t numRanks = numProc();
    const uint32_t myRank   = procID();

    using LevelAndIndex = std::pair<uint32_t, uint32_t>;

    std::vector<ParticleMap<List<P>>> particlesToSend(numRanks);

#pragma omp parallel
    {
      List<P> outcasts;

      std::vector<ParticleMap<List<P>>> threadLocalParticlesToSend(numRanks);

      // Collect particles owned by this rank/thread, but on the other realm.
      a_otherContainer.copyParticlesToSingleList(outcasts, a_otherContainer.getParticles());

      // Map the particles to other grid patches
      this->mapParticlesToAMRGrid(threadLocalParticlesToSend, outcasts);

      // Catenate the thread-local particle lists onto particlesToSend
#pragma omp critical
      this->catenateParticleMaps(particlesToSend, threadLocalParticlesToSend);
    }

    // Particles that go from this rank to this rank don't need to be communicated so that we can place them directly on the correct patch.
    this->assignLocalParticles(particlesToSend[myRank], m_particles);

    // If using MPI then we have to scatter the particles across MPI ranks.
#ifdef CH_MPI
    ParticleMap<List<P>> receivedParticles;

    ParticleOps::scatterParticles(receivedParticles, particlesToSend);

    // Assign particles to the correct level and grid patch -- we iterate through receivedParticles and decode the information
    // we got from there.
    this->assignLocalParticles(receivedParticles, m_particles);
#endif
  }

  if (m_debug) {
    this->sanityCheck();
  }
}

template <class P>
void
ParticleContainer<P>::addParticlesDestructive(ParticleContainer<P>& a_otherContainer)
{
  CH_TIME("ParticleContainer::addParticlesDestructive(ParticleContainer)");
  if (m_verbose) {
    pout() << "ParticleContainer::addParticlesDestructive(ParticleContainer)" << endl;
  }

  // TLDR: This routine adds particles from a different container to this container. If the containers live on the same realm
  //       we just add the particles directly. Otherwise, we have to add to the outcast list and remap.

  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer::addParticles(ParticleContainer<P>) - particles are sorted by cell!");
  }

  // TLDR: This routine adds particles from a different container to this container. If the containers live on the same realm
  //       we just add the particles directly. Otherwise, we have to add to the outcast list and remap.

  if (m_realm == a_otherContainer.getRealm()) {
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const DisjointBoxLayout& dbl = m_grids[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        ListBox<P>&       myParticles    = (*m_particles[lvl])[din];
        const ListBox<P>& otherParticles = a_otherContainer[lvl][din];

        myParticles.addItems(otherParticles.listItems());
      }
    }
  }
  else {
    // If we're adding particles across realms we need to collect the input particles onto lists that get remapped
    // into m_particles. A full explanation of how this works is given in the remap routin.

    const uint32_t numRanks = numProc();
    const uint32_t myRank   = procID();

    using LevelAndIndex = std::pair<uint32_t, uint32_t>;

    std::vector<ParticleMap<List<P>>> particlesToSend(numRanks);

#pragma omp parallel
    {
      List<P> outcasts;

      std::vector<ParticleMap<List<P>>> threadLocalParticlesToSend(numRanks);

      // Collect particles owned by this rank/thread, but on the other realm.
      a_otherContainer.transferParticlesToSingleList(outcasts, a_otherContainer.getParticles());

      // Map the particles to other grid patches
      this->mapParticlesToAMRGrid(threadLocalParticlesToSend, outcasts);

      // Catenate the thread-local particle lists onto particlesToSend
#pragma omp critical
      this->catenateParticleMaps(particlesToSend, threadLocalParticlesToSend);
    }

    // Particles that go from this rank to this rank don't need to be communicated so that we can place them directly on the correct patch.
    this->assignLocalParticles(particlesToSend[myRank], m_particles);

    // If using MPI then we have to scatter the particles across MPI ranks.
#ifdef CH_MPI
    ParticleMap<List<P>> receivedParticles;

    ParticleOps::scatterParticles(receivedParticles, particlesToSend);

    // Assign particles to the correct level and grid patch -- we iterate through receivedParticles and decode the information
    // we got from there.
    this->assignLocalParticles(receivedParticles, m_particles);
#endif
  }

  if (m_debug) {
    this->sanityCheck();
  }
}

template <class P>
void
ParticleContainer<P>::transferParticles(ParticleContainer<P>& a_otherContainer)
{
  CH_TIME("ParticleContainer::transferParticles");

  CH_assert(m_isDefined);
  CH_assert(this->getRealm() == a_otherContainer.getRealm());

  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer::transferParticles(ParticleContainer<P>) - particles are sorted by cell!");
  }

  if (a_otherContainer.isOrganizedByCell()) {
    MayDay::Error("ParticleContainer::transferParticles(ParticleContainer<P>) - other particles are sorted by cell!");
  }

  if (a_otherContainer.getRealm() != m_realm) {
    MayDay::Error(
      "ParticleContainer::transferParticles(ParticleContainer<P>) - other container defined over a different realm");
  }

  this->transferParticles(a_otherContainer.getParticles());
}

template <class P>
void
ParticleContainer<P>::transferParticles(AMRParticles<P>& a_particles)
{
  CH_TIME("ParticleContainer::transferParticles(AMRParticles)");
  if (m_verbose) {
    pout() << "ParticleContainer::transferParticle(AMRParticles)" << endl;
  }

  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer::transferParticles(ParticleContainer<P>) - particles are sorted by cell!");
  }

  // Ok, transfer the particles.
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = m_grids[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      List<P>& myParticles    = (*m_particles[lvl])[din].listItems();
      List<P>& otherParticles = (*a_particles[lvl])[din].listItems();

      myParticles.catenate(otherParticles);
    }
  }
}

template <class P>
void
ParticleContainer<P>::remap()
{
  CH_TIME("ParticleContainer<P>::remap");
  if (m_verbose) {
    pout() << "ParticleContainer::remap" << endl;
  }

  // TLDR: This routine is quite long but does the full remapping on the whole hierarchy. It will discard particles that fall off the grid.
  //
  // This is done in the following steps:
  //    1) Collect all particles from this rank onto thread-local variables
  //    2) Iterate through those particles and figure out where they end up up the AMR hierarchy (level, grid index, and ownership, i.e. the MPI rank)
  //    3) Collect those particles onto a rank-only (no thread loacl) variable
  //    4) Assign particles _locally_, i.e. assign particles sent from this rank to this rank directly onto m_particles
  //    5) If using MPI, scatter the particles to the appropriate ranks.
  //    6) Assign particles locally from the particles that were scattered to this rank.

  const uint32_t numRanks = numProc();
  const uint32_t myRank   = procID();

  using LevelAndIndex = std::pair<uint32_t, uint32_t>;

  // These are the particles that get sent to each MPI rank. The pair contains the grid level and grid index.
  // Each MPI rank (and OpenMP thread) goes through it's particles and checks if they've falled off the grid. If they have,
  // then the particles are moved onto appropriate lists that will get sent to each rank.
  std::vector<ParticleMap<List<P>>> particlesToSend(numRanks);

#pragma omp parallel
  {
    List<P> outcasts;

    std::vector<ParticleMap<List<P>>> threadLocalParticlesToSend(numRanks);

    // Collect particles owned by this rank/thread.
    this->transferParticlesToSingleList(outcasts, m_particles);

    // Map the particles to other grid patches
    this->mapParticlesToAMRGrid(threadLocalParticlesToSend, outcasts);

    // Catenate the thread-local particle lists onto particlesToSend
#pragma omp critical
    this->catenateParticleMaps(particlesToSend, threadLocalParticlesToSend);
  }

  // Particles that go from this rank to this rank don't need to be communicated so that we can place them directly on the correct patch.
  this->assignLocalParticles(particlesToSend[myRank], m_particles);

  // If using MPI then we have to scatter the particles across MPI ranks.
#ifdef CH_MPI
  ParticleMap<List<P>> receivedParticles;

  ParticleOps::scatterParticles(receivedParticles, particlesToSend);

  // Assign particles to the correct level and grid patch -- we iterate through receivedParticles and decode the information
  // we got from there.
  this->assignLocalParticles(receivedParticles, m_particles);
#endif

  if (m_debug) {
    this->sanityCheck();
  }
}

template <typename P>
void
ParticleContainer<P>::sanityCheck() const noexcept
{
  CH_TIME("ParticleContainer::sanityCheck");

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = m_grids[lvl];
    const RealVect           dx  = m_dx[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const List<P>& particles = (*m_particles[lvl])[dit()].listItems();
      const Box      box       = dbl[dit()];

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        const RealVect particlePosition = lit().position();
        const IntVect  particleCell     = ParticleOps::getParticleCellIndex(particlePosition, m_probLo, dx);

        if (!(box.contains(particleCell))) {
          MayDay::Error("ParticleContainer::remap -- error2: particle is not contained in grid box!");
        }
      }
    }
  }
}

template <typename P>
inline void
ParticleContainer<P>::transferParticlesToSingleList(List<P>& a_list, AMRParticles<P>& a_particles) const noexcept
{
  CH_TIME("ParticleContainer::transferParticlesToSingleList");

  for (int lvl = 0; lvl < a_particles.size(); lvl++) {
    const BoxLayout&    dbl = a_particles[lvl]->getBoxes();
    const DataIterator& dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      List<P>& particles = (*a_particles[lvl])[din].listItems();

      a_list.catenate(particles);
    }
  }
}

template <typename P>
inline void
ParticleContainer<P>::copyParticlesToSingleList(List<P>& a_list, const AMRParticles<P>& a_particles) const noexcept
{
  CH_TIME("ParticleContainer::copyParticlesToSingleList");

  for (int lvl = 0; lvl < a_particles.size(); lvl++) {
    const BoxLayout&    dbl = a_particles[lvl]->getBoxes();
    const DataIterator& dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const List<P>& particles = (*a_particles[lvl])[din].listItems();

      a_list.join(particles);
    }
  }
}

template <typename P>
inline void
ParticleContainer<P>::mapParticlesToAMRGrid(std::vector<ParticleMap<List<P>>>& a_mappedParticles,
                                            List<P>&                           a_unmappedParticles) const noexcept
{
  CH_TIME("ParticleContainer::mapParticlesToAMRGrid");

  const uint32_t numRanks = numProc();
  const uint32_t myRank   = procID();

  std::vector<RealVect> quasiDx(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    quasiDx[lvl] = m_blockingFactor * m_dx[lvl];
  }

  // Starting on the finest level -- figure out where the particles wound up. We are going from finest to coarsest here because we want
  // the particle to end up on the finest grid level that contains them.
  for (ListIterator<P> lit(a_unmappedParticles); lit.ok();) {
    bool foundTile = false;

    for (int lvl = m_finestLevel; lvl >= 0 && !foundTile; lvl--) {
      const IntVect particleTile = locateBin(lit().position(), quasiDx[lvl], m_probLo);

      const auto& myLevelTiles    = m_levelTiles[lvl]->getMyTiles();
      const auto& otherLevelTiles = m_levelTiles[lvl]->getOtherTiles();

      if (myLevelTiles.find(particleTile) != myLevelTiles.end()) {
        // Found the particle on this level and this rank is the one who owns it.
        const uint32_t gridIndex = myLevelTiles.at(particleTile);
        const uint32_t toRank    = myRank;

        a_mappedParticles[toRank][std::pair<uint32_t, uint32_t>(lvl, gridIndex)].transfer(lit);

        foundTile = true;
      }
#ifdef CH_MPI
      else if (otherLevelTiles.find(particleTile) != otherLevelTiles.end()) {
        // Particle was found on this level but no on this rank.
        const uint32_t  gridIndex = otherLevelTiles.at(particleTile).first;
        const uint32_t& toRank    = otherLevelTiles.at(particleTile).second;

        a_mappedParticles[toRank][std::pair<uint32_t, uint32_t>(lvl, gridIndex)].transfer(lit);

        foundTile = true;
      }
#endif
    }

    // If this triggers the particle fell off the domain and just move onto the next one.
    if (!foundTile) {
      ++lit;
    }
  }

  // Rest of the particles fell of the domain and have to go.
  a_unmappedParticles.clear();
}

template <typename P>
inline void
ParticleContainer<P>::catenateParticleMaps(std::vector<ParticleMap<List<P>>>& a_globalParticles,
                                           std::vector<ParticleMap<List<P>>>& a_localParticles) const noexcept
{
  CH_TIME("ParticleContainer::catenateParticleMaps");

  const uint32_t numRanks = numProc();

  using LevelAndIndex = std::pair<uint32_t, uint32_t>;

  for (int irank = 0; irank < numRanks; irank++) {
    ParticleMap<List<P>>& globalParticles = a_globalParticles[irank];

    for (auto& m : a_localParticles[irank]) {
      globalParticles[m.first].catenate(m.second);
    }
  }
}

template <typename P>
inline void
ParticleContainer<P>::assignLocalParticles(ParticleMap<List<P>>& a_mappedParticles,
                                           AMRParticles<P>&      a_particleData) const noexcept
{
  CH_TIME("ParticleContainer::assignLocalParticles");

#pragma omp parallel
  {
#pragma omp single
    {
      for (auto mapIter = a_mappedParticles.begin(); mapIter != a_mappedParticles.end(); mapIter++) {
#pragma omp task firstprivate(mapIter)
        {
          const uint32_t  gridLevel = mapIter->first.first;
          const uint32_t  gridIndex = mapIter->first.second;
          const DataIndex din       = m_levelTiles[gridLevel]->getMyGrids().at(gridIndex);

          List<P>& particles   = mapIter->second;
          List<P>& myParticles = (*a_particleData[gridLevel])[din].listItems();

          myParticles.catenate(particles);
        }
      }
    }
  }

  a_mappedParticles.clear();
}

template <class P>
void
ParticleContainer<P>::preRegrid(const int a_lmin)
{
  CH_TIME("ParticleContainer::preRegrid");
  if (m_verbose) {
    pout() << "ParticleContainer::preRegrid" << endl;
  }

  // TLDR: This routine takes all the particles on each level and puts them in a single list (per processor). After than
  //       we can safely destruct the ParticleData on each level without losing the particles.

  CH_assert(m_isDefined);
  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer::preRegrid - particles are sorted by cell!");
  }

  // Fill cache particles on each level
  m_cacheParticles.resize(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    m_cacheParticles[lvl] = RefCountedPtr<ParticleData<P>>(
      new ParticleData<P>(m_grids[lvl], m_domains[lvl], m_blockingFactor, m_dx[lvl], m_probLo));

    const DisjointBoxLayout& dbl = m_grids[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      List<P>& cacheParticles  = (*m_cacheParticles[lvl])[din].listItems();
      List<P>& solverParticles = (*m_particles[lvl])[din].listItems();

      cacheParticles.catenate(solverParticles);
    }
  }
}

template <class P>
void
ParticleContainer<P>::regrid(const Vector<DisjointBoxLayout>&         a_grids,
                             const Vector<ProblemDomain>&             a_domains,
                             const Vector<Real>&                      a_dx,
                             const Vector<int>&                       a_refRat,
                             const Vector<ValidMask>&                 a_validMask,
                             const Vector<RefCountedPtr<LevelTiles>>& a_levelTiles,
                             const int                                a_lmin,
                             const int                                a_newFinestLevel)
{
  CH_TIME("ParticleContainer::regrid");
  if (m_verbose) {
    pout() << "ParticleContainer::regrid" << endl;
  }

  // TLDR: This is the main regrid routine for particle container. This is essentially a ::define() followed by several types of remapping steps
  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer::regrid(...) - particles are sorted by cell!");
  }

  // Update this stuff
  m_grids       = a_grids;
  m_domains     = a_domains;
  m_refRat      = a_refRat;
  m_validRegion = a_validMask;
  m_finestLevel = a_newFinestLevel;
  m_levelTiles  = a_levelTiles;

  m_dx.resize(1 + m_finestLevel);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    m_dx[lvl] = a_dx[lvl] * RealVect::Unit;
  }

  this->setupGrownGrids(a_lmin, m_finestLevel);
  this->setupParticleData(a_lmin, m_finestLevel);

  // Perform the remapping operation.
  const uint32_t numRanks = numProc();
  const uint32_t myRank   = procID();

  using LevelAndIndex = std::pair<uint32_t, uint32_t>;

  // These are the particles that get sent to each MPI rank. The pair contains the grid level and grid index.
  // Each MPI rank (and OpenMP thread) goes through it's particles and checks if they've falled off the grid. If they have,
  // then the particles are moved onto appropriate lists that will get sent to each rank.
  std::vector<ParticleMap<List<P>>> particlesToSend(numRanks);

#pragma omp parallel
  {
    List<P> outcasts;

    std::vector<ParticleMap<List<P>>> threadLocalParticlesToSend(numRanks);

    // Collect particles owned by this rank/thread.
    this->transferParticlesToSingleList(outcasts, m_cacheParticles);

    // Map the particles to other grid patches
    this->mapParticlesToAMRGrid(threadLocalParticlesToSend, outcasts);

    // Catenate the thread-local particle lists onto particlesToSend
#pragma omp critical
    this->catenateParticleMaps(particlesToSend, threadLocalParticlesToSend);
  }

  // Particles that go from this rank to this rank don't need to be communicated so that we can place them directly on the correct patch.
  this->assignLocalParticles(particlesToSend[myRank], m_particles);

  // If using MPI then we have to scatter the particles across MPI ranks.
#ifdef CH_MPI
  ParticleMap<List<P>> receivedParticles;

  ParticleOps::scatterParticles(receivedParticles, particlesToSend);

  // Assign particles to the correct level and grid patch -- we iterate through receivedParticles and decode the information
  // we got from there.
  this->assignLocalParticles(receivedParticles, m_particles);
#endif

  if (m_debug) {
    this->sanityCheck();
  }

  // Discard transient storage.
  m_cacheParticles.resize(0);
}

template <class P>
template <Real& (P::*particleScalarField)()>
void
ParticleContainer<P>::setValue(const Real a_value)
{
  CH_TIME("ParticleContainer::setValue");

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = m_grids[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      List<P>& patchParticles = (*m_particles[lvl])[din].listItems();

      for (ListIterator<P> lit(patchParticles); lit.ok(); ++lit) {
        P& p = lit();

        (p.*particleScalarField)() = a_value;
      }
    }
  }
}

template <class P>
template <RealVect& (P::*particleVectorField)()>
void
ParticleContainer<P>::setValue(const RealVect a_value)
{
  CH_TIME("ParticleContainer::setValue");

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = m_grids[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      List<P>& patchParticles = (*m_particles[lvl])[din].listItems();

      for (ListIterator<P> lit(patchParticles); lit.ok(); ++lit) {
        P& p = lit();

        (p.*particleVectorField)() = a_value;
      }
    }
  }
}

template <class P>
unsigned long long
ParticleContainer<P>::getNumberOfValidParticlesLocal() const
{
  CH_TIME("ParticleContainer::getNumberOfValidParticlesLocal");

  CH_assert(m_isDefined);

  unsigned long long n = 0;

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = m_grids[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(+ : n)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const List<P>& particles = (*m_particles[lvl])[din].listItems();

      n += (unsigned long long)particles.length();
    }
  }

  return n;
}

template <class P>
unsigned long long
ParticleContainer<P>::getNumberOfValidParticlesGlobal() const
{
  CH_TIME("ParticleContainer::getNumberOfValidParticlesGlobal");

  CH_assert(m_isDefined);

  const unsigned long long n = this->getNumberOfValidParticlesLocal();

  return ParallelOps::sum(n);
}

template <class P>
unsigned long long
ParticleContainer<P>::getNumberOfOutcastParticlesLocal() const
{
  CH_TIME("ParticleContainer::getNumberOfOutcastParticlesLocal");

  CH_assert(m_isDefined);

  unsigned long long n = 0;

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = m_grids[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(+ : n)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const List<P>& outcast = m_particles[lvl]->outcast();

      n += (unsigned long long)outcast.length();
    }
  }

  return n;
}

template <class P>
unsigned long long
ParticleContainer<P>::getNumberOfOutcastParticlesGlobal() const
{
  CH_TIME("ParticleContainer::getNumberOfOutcastParticlesGlobal");

  CH_assert(m_isDefined);

  const unsigned long long n = this->getNumberOfOutcastParticlesLocal();

  return ParallelOps::sum(n);
}

template <class P>
unsigned long long
ParticleContainer<P>::getNumberOfMaskParticlesLocal() const
{
  CH_TIME("ParticleContainer::getNumberOfMaskParticlesLocal");

  CH_assert(m_isDefined);

  unsigned long long n = 0;

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = m_grids[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(+ : n)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const List<P>& maskParticles = (*m_maskParticles[lvl])[din].listItems();

      n += (unsigned long long)maskParticles.length();
    }
  }

  return n;
}

template <class P>
unsigned long long
ParticleContainer<P>::getNumberOfMaskParticlesGlobal() const
{
  CH_TIME("ParticleContainer::getNumberOfMaskParticlesGlobal");

  CH_assert(m_isDefined);

  const unsigned long long n = this->getNumberOfMaskParticlesLocal();

  return ParallelOps::sum(n);
}

template <class P>
void
ParticleContainer<P>::copyMaskParticles(const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_mask) const
{
  CH_TIME("ParticleContainer::copyMaskParticles(amr)");

  CH_assert(m_isDefined);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    if (!a_mask[lvl].isNull()) {
      this->copyMaskParticles(lvl, *a_mask[lvl]);
    }
  }
}

template <class P>
void
ParticleContainer<P>::copyMaskParticles(const int a_level, const LevelData<BaseFab<bool>>& a_mask) const
{
  CH_TIME("ParticleContainer::copyMaskParticles(level)");
  if (m_verbose) {
    pout() << "ParticleContainer::copyMaskParticles(level)" << endl;
  }

  CH_assert(m_isDefined);

  m_maskParticles[a_level]->clear();

  const RealVect dx = m_dx[a_level];

  // Copy particles from m_particles to m_maskParticles if they lie in the input mask.
  const DisjointBoxLayout& dbl = m_grids[a_level];
  const DataIterator&      dit = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const BaseFab<bool>& mask = a_mask[din];

    if (mask.isUsable()) {
      const Box gridBox = m_grids[a_level][din];
      const Box maskBox = mask.box();

      CH_assert(gridBox == maskBox);

      List<P>&       maskParticles = (*m_maskParticles[a_level])[din].listItems();
      const List<P>& particles     = (*m_particles[a_level])[din].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        const RealVect x  = lit().position();
        const IntVect  iv = ParticleOps::getParticleCellIndex(x, m_probLo, dx);

        if (!(maskBox.contains(iv))) {
          std::cout << a_level << "\t" << iv << "\t" << x << "\t" << maskBox << std::endl;

          MayDay::Error("ParticleContainer::copyMaskParticles -- logic bust. Particle has fallen off grid");
        }
        else if (mask(iv)) {
          maskParticles.add(lit());
        }
      }
    }
  }
}

template <class P>
void
ParticleContainer<P>::transferMaskParticles(const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_mask)
{
  CH_TIME("ParticleContainer::transferMaskParticles(amr)");

  CH_assert(m_isDefined);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    if (!a_mask[lvl].isNull()) {
      this->transferMaskParticles(lvl, *a_mask[lvl]);
    }
  }
}

template <class P>
void
ParticleContainer<P>::transferMaskParticles(const int a_level, const LevelData<BaseFab<bool>>& a_mask)
{
  CH_TIME("ParticleContainer::transferMaskParticles(level)");
  if (m_verbose) {
    pout() << "ParticleContainer::transferMaskParticles(level)" << endl;
  }

  CH_assert(m_isDefined);

  const RealVect dx = m_dx[a_level];

  // Copy particles from m_particles to m_maskParticles if they lie in the input mask.
  const DisjointBoxLayout& dbl = m_grids[a_level];
  const DataIterator&      dit = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const BaseFab<bool>& mask = a_mask[din];

    if (mask.isUsable()) {
      const Box gridBox = m_grids[a_level][din];
      const Box maskBox = mask.box();

      CH_assert(gridBox == maskBox);

      List<P>&       maskParticles = (*m_maskParticles[a_level])[din].listItems();
      const List<P>& particles     = (*m_particles[a_level])[din].listItems();

      for (ListIterator<P> lit(particles); lit.ok();) {
        const RealVect x  = lit().position();
        const IntVect  iv = ParticleOps::getParticleCellIndex(x, m_probLo, dx);

        if (!(maskBox.contains(iv))) {
          //        std::cout << a_level << "\t" << iv << "\t" << x << std::endl;

          MayDay::Warning("ParticleContainer::transferMaskParticles -- logic bust. Particle has fallen off grid");

          ++lit;
        }
        else {
          if (mask(iv)) {
            maskParticles.transfer(lit);
          }
          else {
            ++lit;
          }
        }
      }
    }
  }
}

template <class P>
void
ParticleContainer<P>::clearParticles()
{
  CH_assert(m_isDefined);

  this->clear(m_particles);
}

template <class P>
void
ParticleContainer<P>::clearBufferParticles() const
{
  CH_assert(m_isDefined);

  this->clear(m_bufferParticles);
}

template <class P>
void
ParticleContainer<P>::clearMaskParticles() const
{
  CH_assert(m_isDefined);

  this->clear(m_maskParticles);
}

template <class P>
void
ParticleContainer<P>::clearOutcast() noexcept
{
  CH_assert(m_isDefined);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    ParticleData<P>& particleData = *m_particles[lvl];

    List<P>& outcast = particleData.outcast();

    outcast.clear();
  }
}

template <class P>
void
ParticleContainer<P>::clear(AMRParticles<P>& a_particles) const
{
  if (m_verbose) {
    pout() << "ParticleContainer::clear(AMRParticles)" << endl;
  }

  CH_assert(m_isDefined);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {

    ParticleData<P>& levelParticles = *a_particles[lvl];

    // Clear patch particles
    const BoxLayout&    grids = levelParticles.getBoxes();
    const DataIterator& dit   = grids.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex din = dit[mybox];

      List<P>& patchParticles = levelParticles[din].listItems();

      patchParticles.clear();
    }

    // Clear outcast
    List<P>& outcast = levelParticles.outcast();
    outcast.clear();
  }
}

#include <CD_NamespaceFooter.H>

#endif
