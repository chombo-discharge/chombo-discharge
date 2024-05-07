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
}

template <class P>
ParticleContainer<P>::ParticleContainer(const Vector<DisjointBoxLayout>& a_grids,
                                        const Vector<ProblemDomain>&     a_domains,
                                        const Vector<Real>&              a_dx,
                                        const Vector<int>&               a_refRat,
                                        const Vector<ValidMask>&         a_validMask,
                                        const RealVect&                  a_probLo,
                                        const int                        a_blockingFactor,
                                        const int                        a_finestLevel,
                                        const std::string                a_realm)
{
  CH_TIME("ParticleContainer<P>::ParticleContainer");

  this->define(a_grids, a_domains, a_dx, a_refRat, a_validMask, a_probLo, a_blockingFactor, a_finestLevel, a_realm);
}

template <class P>
ParticleContainer<P>::~ParticleContainer()
{
  CH_TIME("ParticleContainer<P>::~ParticleContainer");
}

template <class P>
void
ParticleContainer<P>::define(const Vector<DisjointBoxLayout>& a_grids,
                             const Vector<ProblemDomain>&     a_domains,
                             const Vector<Real>&              a_dx,
                             const Vector<int>&               a_refRat,
                             const Vector<ValidMask>&         a_validMask,
                             const RealVect&                  a_probLo,
                             const int                        a_blockingFactor,
                             const int                        a_finestLevel,
                             const std::string                a_realm)
{
  CH_TIME("ParticleContainer<P>::~ParticleContainer");

  m_grids          = a_grids;
  m_domains        = a_domains;
  m_refRat         = a_refRat;
  m_validRegion    = a_validMask;
  m_probLo         = a_probLo;
  m_finestLevel    = a_finestLevel;
  m_realm          = a_realm;
  m_blockingFactor = a_blockingFactor;

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
}

template <class P>
void
ParticleContainer<P>::setupGrownGrids(const int a_base, const int a_finestLevel)
{
  CH_TIME("ParticleContainer<P>::setupGrownGrids");

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
  CH_TIME("ParticleContainer<P>::setupParticleData");

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
  CH_TIME("ParticleContainer<P>::sortParticles()");

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
    MayDay::Error("ParticleContainer<P>::getParticles - particles are sorted by cell!");
  }

  return m_particles;
}

template <class P>
const AMRParticles<P>&
ParticleContainer<P>::getParticles() const
{
  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {
    MayDay::Abort("ParticleContainer<P>::getParticles - particles are sorted by cell!");
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
    MayDay::Error("ParticleContainer<P>::operator[](a_lvl) - particles are sorted by cell!");
  }

  return *m_particles[a_lvl];
}

template <class P>
const ParticleData<P>&
ParticleContainer<P>::operator[](const int a_level) const
{
  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer<P>::operator[](a_lvl) - particles are sorted by cell!");
  }

  return *m_particles[a_level];
}

template <class P>
AMRCellParticles<P>&
ParticleContainer<P>::getCellParticles()
{
  CH_assert(m_isDefined);

  if (!m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer<P>::getCellParticles()- particles are not sorted by cell!");
  }

  return m_cellSortedParticles;
}

template <class P>
const AMRCellParticles<P>&
ParticleContainer<P>::getCellParticles() const
{
  CH_assert(m_isDefined);

  if (!m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer<P>::getCellParticles()- particles are not sorted by cell!");
  }

  return m_cellSortedParticles;
}

template <class P>
LayoutData<BinFab<P>>&
ParticleContainer<P>::getCellParticles(const int a_level)
{
  CH_assert(m_isDefined);

  if (!m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer<P>::getCellParticles(level)- particles are not sorted by cell!");
  }

  return *m_cellSortedParticles[a_level];
}

template <class P>
const LayoutData<BinFab<P>>&
ParticleContainer<P>::getCellParticles(const int a_level) const
{
  CH_TIME("ParticleContainer<P>::getCellParticles(int)");

  CH_assert(m_isDefined);

  if (!m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer<P>::getCellParticles(level)- particles are not sorted by cell!");
  }

  return *m_cellSortedParticles[a_level];
}

template <class P>
void
ParticleContainer<P>::getCellParticles(BinFab<P>& cellParticles, const int a_lvl, const DataIndex a_dit) const
{
  CH_TIME("ParticleContainer<P>::getCellParticles(BinFab)");

  CH_assert(m_isDefined);

  if (!m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer<P>::getCellParticles - particles are not sorted by cell!");
  }

  cellParticles.define(m_grids[a_lvl][a_dit], m_dx[a_lvl], m_probLo);
  cellParticles.addItems((*m_particles[a_lvl])[a_dit].listItems());
}

template <class P>
void
ParticleContainer<P>::getCellParticlesDestructive(BinFab<P>& cellParticles, const int a_lvl, const DataIndex a_dit)
{
  CH_TIME("ParticleContainer<P>::getCellParticlesDestructive");

  CH_assert(m_isDefined);

  if (!m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer<P>::getCellParticlesDestructive - particles are not sorted by cell!");
  }

  cellParticles.define(m_grids[a_lvl].get(a_dit), m_dx[a_lvl], m_probLo);
  cellParticles.addItemsDestructive((*m_particles[a_lvl])[a_dit].listItems());
}

template <class P>
BinFab<P>&
ParticleContainer<P>::getCellParticles(const int a_level, const DataIndex a_dit)
{
  CH_TIME("ParticleContainer<P>::getCellParticles(int, DataIndex)");

  CH_assert(m_isDefined);

  if (!m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer<P>::getCellParticles(int, dit) - particles are not sorted by cell!");
  }

  return (*m_cellSortedParticles[a_level])[a_dit];
}

template <class P>
const BinFab<P>&
ParticleContainer<P>::getCellParticles(const int a_level, const DataIndex a_dit) const
{
  CH_TIME("ParticleContainer<P>::getCellParticles(int, DataIndex)");

  CH_assert(m_isDefined);

  if (!m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer<P>::getCellParticles(int, dit) - particles are not sorted by cell!");
  }

  return (*m_cellSortedParticles[a_level])[a_dit];
}

template <class P>
void
ParticleContainer<P>::organizeParticlesByCell()
{
  CH_TIME("ParticleContainer<P>::organizeParticlesByCell");

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
  CH_TIME("ParticleContainer<P>::organizeParticlesByPatch");

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
  CH_TIME("ParticleContainer<P>::addParticles");

  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer<P>::addParticles(List<P>) - particles are sorted by cell!");
  }

  List<P>& outcastFinest = m_particles[m_finestLevel]->outcast();
  outcastFinest.join(a_particles);

  for (int lvl = m_finestLevel; lvl >= 0; lvl--) {
    this->remapOutcast(*m_particles[lvl]);

    if (lvl > 0) {
      List<P>& outcast     = m_particles[lvl]->outcast();
      List<P>& outcastCoar = m_particles[lvl - 1]->outcast();

      outcastCoar.catenate(outcast);
    }
  }
}

template <class P>
void
ParticleContainer<P>::addParticlesDestructive(List<P>& a_particles)
{
  CH_TIME("ParticleContainer<P>::addParticlesDestructive");

  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer<P>::addParticlesDestructive(List<P>) - particles are sorted by cell!");
  }

  List<P>& outcastFinest = m_particles[m_finestLevel]->outcast();
  outcastFinest.join(a_particles);

  for (int lvl = m_finestLevel; lvl >= 0; lvl--) {
    this->remapOutcast(*m_particles[lvl]);

    if (lvl > 0) {
      List<P>& outcast     = m_particles[lvl]->outcast();
      List<P>& outcastCoar = m_particles[lvl - 1]->outcast();

      outcastCoar.catenate(outcast);
    }
  }
}

template <class P>
void
ParticleContainer<P>::addParticles(const BinFab<P>& a_particles, const int a_lvl, const DataIndex a_dit)
{
  CH_TIME("ParticleContainer<P>::addParticles(BinFab)");

  CH_assert(m_isDefined);
  CH_assert(m_grids[a_lvl].get(a_dit) == a_particles.getRegion());

  if (!m_isOrganizedByCell) {
    MayDay::Abort("ParticleContainer<P>::addParticles(BinFab<P>) - particles are not sorted by cell!");
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
  CH_TIME("ParticleContainer<P>::addParticlesDestructive(BinFab)");

  CH_assert(m_isDefined);
  CH_assert(m_grids[a_lvl].get(a_dit) == a_particles.getRegion());

  if (!m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer<P>::addParticles(BinFab<P>) - particles are not sorted by cell!");
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
  CH_TIME("ParticleContainer<P>::addParticles(ParticleContainer)");

  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");
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
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      List<P>& outcast = m_particles[lvl]->outcast();

      outcast.clear();

      // Iterate through the other container and add its particles to this container's outcast lists.
      const DisjointBoxLayout& dbl = a_otherContainer.getGrids()[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(join : outcast)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        const List<P>& otherParticles = a_otherContainer[lvl][din].listItems();

        outcast.join(otherParticles);
      }
    }

    // Remap so that particles end up in the correct location.
    this->remap();
  }
}

template <class P>
void
ParticleContainer<P>::addParticlesDestructive(ParticleContainer<P>& a_otherContainer)
{
  CH_TIME("ParticleContainer<P>::addParticlesDestructive(ParticleContainer)");

  // TLDR: This routine adds particles from a different container to this container. If the containers live on the same realm
  //       we just add the particles directly. Otherwise, we have to add to the outcast list and remap.

  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");
  }

  if (m_realm == a_otherContainer.getRealm()) {
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const DisjointBoxLayout& dbl = m_grids[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        ListBox<P>& myParticles    = (*m_particles[lvl])[din];
        ListBox<P>& otherParticles = a_otherContainer[lvl][din];

        myParticles.addItemsDestructive(otherParticles.listItems());
      }
    }
  }
  else {
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      List<P>& outcast = m_particles[lvl]->outcast();
      outcast.clear();

      // Iterate through the other container and add its particles to this container's outcast lists.
      const DisjointBoxLayout& dbl = a_otherContainer.getGrids()[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(catenate : outcast)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        List<P>& otherParticles = a_otherContainer[lvl][din].listItems();

        outcast.catenate(otherParticles);

        otherParticles.clear();
      }
    }

    // Remap so that particles end up on the correct patch.
    this->remap();
  }
}

template <class P>
void
ParticleContainer<P>::transferParticles(ParticleContainer<P>& a_otherContainer)
{
  CH_TIME("ParticleContainer<P>::transferParticles");

  CH_assert(m_isDefined);
  CH_assert(this->getRealm() == a_otherContainer.getRealm());

  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer<P>::transferParticles(ParticleContainer<P>) - particles are sorted by cell!");
  }

  if (a_otherContainer.isOrganizedByCell()) {
    MayDay::Error(
      "ParticleContainer<P>::transferParticles(ParticleContainer<P>) - other particles are sorted by cell!");
  }

  if (a_otherContainer.getRealm() != m_realm) {
    MayDay::Error(
      "ParticleContainer<P>::transferParticles(ParticleContainer<P>) - other container defined over a different realm");
  }

  this->transferParticles(a_otherContainer.getParticles());
}

template <class P>
void
ParticleContainer<P>::transferParticles(AMRParticles<P>& a_particles)
{
  CH_TIME("ParticleContainer<P>::transferParticles(AMRParticles)");

  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer<P>::transferParticles(ParticleContainer<P>) - particles are sorted by cell!");
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
ParticleContainer<P>::evictInvalidParticles(List<P>&         a_evictedParticles,
                                            ParticleData<P>& a_particles,
                                            const int        a_level)
{
  CH_TIME("ParticleContainer::evictInvalidParticles");

  const DisjointBoxLayout& dbl = m_grids[a_level];
  const RealVect           dx  = m_dx[a_level];
  const DataIterator&      dit = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const BaseFab<bool>& mask = (*m_validRegion[a_level])[din];

    for (ListIterator<P> lit(a_particles[din].listItems()); lit.ok();) {

      const RealVect x = lit().position();
      const IntVect  iv(D_DECL((int)std::floor((x[0] - m_probLo[0]) / dx[0]),
                              (int)std::floor((x[1] - m_probLo[1]) / dx[1]),
                              (int)std::floor((x[2] - m_probLo[2]) / dx[2])));

      CH_assert(mask.box().contains(iv));

      if (mask(iv)) {
        ++lit;
      }
      else {
        a_evictedParticles.transfer(lit);
      }
    }
  }
}

template <class P>
void
ParticleContainer<P>::remap()
{
  CH_TIMERS("ParticleContainer<P>::remap");
  CH_TIMER("ParticleContainer<P>::remap::catenate_to_finest", t1);
  CH_TIMER("ParticleContainer<P>::remap::catenate_to_coarser", t2);

  // This is old code.
#if 0
  CH_START(t1);
  List<P>& outcastFinest = m_particles[m_finestLevel]->outcast();
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {

    const DisjointBoxLayout& dbl = m_grids[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(catenate : outcastFinest)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      outcastFinest.catenate((*m_particles[lvl])[din].listItems());
    }
  }
  CH_STOP(t1);

  for (int lvl = m_finestLevel; lvl >= 0; lvl--) {
    this->gatherOutcast(*m_particles[lvl]);
    this->remapOutcast(*m_particles[lvl]);

    CH_START(t2);
    if (lvl > 0) {

      List<P>& outcast     = m_particles[lvl]->outcast();
      List<P>& outcastCoar = m_particles[lvl - 1]->outcast();

      outcastCoar.catenate(outcast);
    }
    CH_STOP(t2);
  }

  // Clear rest of particles -- these are no longer on the grid.
  this->clearOutcast();

#else
  // Here starts the new way.
  const unsigned int numRanks = numProc();
  const unsigned int myRank   = procID();

  // Build the tile-space. This is something that should probably be stored somewhere else.
  //
  // tileSet => Grids, viewed as a collection of tiles (one box = one tile)
  // myTiles => Tiles on each level owned by this rank. Maps a tile (an IntVect) to a grid index.
  // myGridIndices => Maps a global grid index to a local DataIndex
  // othertiles => Tiles on each level owned by another rank.
  std::vector<std::set<IntVect, CompIntVect>>               tileSet(1 + m_finestLevel);
  std::vector<std::map<IntVect, unsigned int, CompIntVect>> myTiles(1 + m_finestLevel);
  std::vector<std::map<unsigned int, DataIndex>>            myGridIndices(1 + m_finestLevel);
  std::vector<std::map<IntVect, boxids, CompIntVect>>       otherTiles(1 + m_finestLevel);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = m_grids[lvl];

    for (LayoutIterator lit = dbl.layoutIterator(); lit.ok(); ++lit) {
      const LayoutIndex  lidx   = lit();
      const IntVect      tile   = coarsen(dbl[lidx], m_blockingFactor).smallEnd();
      const unsigned int rankID = dbl.procID(lidx);
      const unsigned int tileID = dbl.index(lidx);

      // If using MPI we need to figure out who owns this tile.
#if CH_MPI
      if (myRank == rankID) {
        myTiles[lvl][tile] = tileID;
      }
      else {
        otherTiles[lvl][tile] = boxids(tileID, rankID);
      }
#else
      myTiles[lvl][tile] = tileID;
#endif

      tileSet[lvl].insert(tile);
    }

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      myGridIndices[lvl][dit().intCode()] = dit();
    }
  }

  // These routines are used for encoding a grid level and a grid index onto an unsigned long. The level is encoded onto an 8-bit representation (max 256 levels)
  // onto the left-most bits. The grid index is encoded by the remaining (hopefully 56) bits, which should be enough for anyone.
  using IntType                 = unsigned int;
  constexpr size_t numBitsLong  = sizeof(IntType) * 8;
  constexpr size_t numBitsLevel = 8;
  constexpr size_t numBitsIndex = numBitsLong - numBitsLevel;

  auto encode = [](const IntType a_level, const IntType a_index) -> IntType {
    std::bitset<numBitsLevel> levelBits(a_level);
    std::bitset<numBitsLong>  codeBits(a_index);

    for (size_t pos = 0; pos < numBitsLevel; pos++) {
      codeBits.set(numBitsIndex + pos, levelBits[pos]);
    }

    return (IntType)codeBits.to_ulong();
  };

  auto decode = [](const IntType a_code) -> std::pair<IntType, IntType> {
    std::bitset<numBitsLevel> levelBits;
    std::bitset<numBitsLong>  indexBits(a_code);

    for (size_t pos = 0; pos < numBitsLevel; pos++) {
      levelBits.set(pos, indexBits[numBitsIndex + pos]);
      indexBits.set(numBitsIndex + pos, false);
    }

    const IntType lvl = (IntType)levelBits.to_ulong();
    const IntType idx = (IntType)indexBits.to_ulong();

    return std::make_pair(lvl, idx);
  };

  // These are the particles that get sent to each MPI rank. The unsigned integer in the map holds both the grid index AND the level index,
  // which we do by using the first six bits to decode the level and the remaining 26 bits to hold the grid index.
  std::vector<std::map<IntType, List<P>>> particlesToSend(numRanks);

  // Particles received on this rank - we use the same indexing/decoding as for the sent particles
  std::map<IntType, List<P>> receivedParticles;

  // Each MPI rank (and OpenMP thread) goes through it's particles and checks if they've falled off the grid. If they have,
  // then the particles are moved onto appropriate lists that will get sent to each rank.
#pragma omp parallel
  {
    std::vector<std::map<IntType, List<P>>> threadLocalParticlesToSend(numRanks);

    // Collect the outcast list for each rank.
    List<P> outcasts;

#if 1 // Remove later
    pout() << "collecting outcasts" << endl;
#endif
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const DisjointBoxLayout& dbl = m_grids[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

#if 0 // This was the original code but I don't think it's correct -- we should just iterate through all particles
        ListBox<P>&      particles = (*m_particles[lvl])[din];

        particles.getInvalidDestructive(outcasts, dbl[din]);
#else
        List<P>& particles = (*m_particles[lvl])[din].listItems();
        outcasts.catenate(particles);
#endif
      }
    }
#if 1
    pout() << "done collecting outcasts" << endl;
#endif

    // Starting on the finest level -- figure out where the particles wound up. We are going from finest to coarsest here because we want
    // the particle to end up on the finest grid level that contains them.
#if 1 // Remove later
#pragma omp barrier
    pout() << "start outcast analysis" << endl;
#endif
    for (ListIterator<P> lit(outcasts); lit.ok();) {
      bool foundTile = false;

      for (IntType lvl = m_finestLevel; lvl >= 0; lvl--) {
        const IntVect particleTile = locateBin(lit().position(), m_dx[lvl], m_probLo);

        const auto& myLevelTiles    = myTiles[lvl];
        const auto& otherLevelTiles = otherTiles[lvl];

        if (myLevelTiles.find(particleTile) != myLevelTiles.end()) {
          // Found the particle on this level and this rank is the one who owns it.
          const IntType gridIndex = myLevelTiles.at(particleTile);
          const IntType gridCode  = encode(lvl, gridIndex);

          // threadLocalParticlesToSend[myRank][gridCode].transfer(lit);

          // foundTile = true;

          // break;
        }
        // #ifdef CH_MPI
        //         else if (otherLevelTiles.find(particleTile) != otherLevelTiles.end()) {
        //           // Particle was found on this level but no on this rank.
        //           const IntType  gridIndex = otherLevelTiles.at(particleTile).idx;
        //           const IntType  gridCode  = encode(lvl, gridIndex);
        //           const IntType& toRank    = otherLevelTiles.at(particleTile).pid;

        //           threadLocalParticlesToSend[toRank][gridCode].transfer(lit);

        //           foundTile = true;

        //           break;
        //         }
        // #endif
      }

      // If this triggers the particle fell off the domain.
      if (!foundTile) {
        ++lit;
      }
    }
#if 1 // Remove later
#pragma omp barrier
    pout() << "done collecting outcasts" << endl;
    pout() << "starting outcast collection" << endl;
#endif

#pragma omp critical
    {
      // Collect the particles to be sent to each rank.
      for (int irank = 0; irank < numRanks; irank++) {
        std::map<IntType, List<P>>& particlesToRank = particlesToSend[irank];

        for (auto& m : threadLocalParticlesToSend[irank]) {
          const IntType& gridCode  = m.first;
          List<P>&       particles = m.second;

          particlesToRank[gridCode].catenate(particles);
        }
      }
    }

    // Rest of the particles fell of the domain and have to go.
    outcasts.clear();

#if 1 // Remove later
#pragma omp barrier
    pout() << "done with omp outcast loop" << endl;
#endif
  }

  // Iterator for gonig through the particle lists.
  typename std::map<IntType, List<P>>::iterator mapIter;

  // Particles that go from this rank to this rank (only change patch) don't need to be communicated we can place them directly on the correct patch.
  std::map<IntType, List<P>>& particlesFromMeToMe = particlesToSend[myRank];

#if 1 // Remove later
  pout() << "starting to collect particles to send" << endl;
#endif
#pragma omp        parallel
#pragma omp single nowait
  {
    for (mapIter = particlesFromMeToMe.begin(); mapIter != particlesFromMeToMe.end(); mapIter++) {
#pragma omp task   firstprivate(mapIter)
      {
        const IntType gridCode  = mapIter->first;
        List<P>&      particles = mapIter->second;

        const std::pair<unsigned int, unsigned int> gridDecode = decode(gridCode);
        const unsigned int                          gridLevel  = gridDecode.first;
        const unsigned int                          gridIndex  = gridDecode.second;

        const DataIndex din = myGridIndices[gridLevel][gridIndex];

        (*m_particles[gridLevel])[din].listItems().catenate(particles);
      }
    }
  }
#if 1 // Remove later
  ParallelOps::barrier();
  pout() << "done collecting particles" << endl;
  pout() << "starting particle scatter " << endl;
#endif

#ifdef CH_MPI
  // TODO: This one should be rewritten to use unsigned long rather than unsigned int
  mpi_scatter_part(receivedParticles, particlesToSend);
#endif
#if 1 // Remove later
  ParallelOps::barrier();
  pout() << "done with particle scatter" << endl;
#endif

  // Assign particles to the correct level and grid patch -- we iterate through receivedParticles and decode the information
  // we got from there.
#pragma omp        parallel
#pragma omp single nowait
  {
    for (mapIter = receivedParticles.begin(); mapIter != receivedParticles.end(); mapIter++) {
#pragma omp task   firstprivate(mapIter)
      {
        const IntType gridCode  = mapIter->first;
        List<P>&      particles = mapIter->second;

        const std::pair<unsigned int, unsigned int> gridDecode = decode(gridCode);
        const unsigned int                          gridLevel  = gridDecode.first;
        const unsigned int                          gridIndex  = gridDecode.second;

        const DataIndex din = myGridIndices[gridLevel][gridIndex];

        (*m_particles[gridLevel])[din].listItems().catenate(particles);
      }
    }
  }
#endif
}

template <class P>
void
ParticleContainer<P>::levelRemap()
{
  CH_TIME("ParticleContainer<P>::levelRemap");

  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer<P>::addParticles(ParticleContainer<P>) - particles are sorted by cell!");
  }

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    this->levelRemap(lvl);
  }
}

template <class P>
void
ParticleContainer<P>::levelRemap(const int a_lvl)
{
  CH_TIME("ParticleContainer<P>::levelRemap(int)");

  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer<P>::levelRemap - particles are sorted by cell!");
  }

  this->gatherOutcast(*m_particles[a_lvl]);
  this->remapOutcast(*m_particles[a_lvl]);
}

template <class P>
void
ParticleContainer<P>::preRegrid(const int a_lmin)
{
  CH_TIME("ParticleContainer<P>::preRegrid");

  // TLDR: This routine takes all the particles on each level and puts them in a single list (per processor). After than
  //       we can safely destruct the ParticleData on each level without losing the particles.

  CH_assert(m_isDefined);
  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer<P>::preRegrid - particles are sorted by cell!");
  }

  // Fill cache particles on each level
  m_cacheParticles.resize(1 + m_finestLevel);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {

    List<P>& cacheParticles = m_cacheParticles[lvl];

    cacheParticles.clear();

    const DisjointBoxLayout& dbl = m_grids[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(catenate : cacheParticles)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      List<P>& p = (*m_particles[lvl])[din].listItems();

      cacheParticles.catenate(p);

      p.clear();
    }
  }
}

template <class P>
void
ParticleContainer<P>::regrid(const Vector<DisjointBoxLayout>& a_grids,
                             const Vector<ProblemDomain>&     a_domains,
                             const Vector<Real>&              a_dx,
                             const Vector<int>&               a_refRat,
                             const Vector<ValidMask>&         a_validMask,
                             const int                        a_lmin,
                             const int                        a_newFinestLevel)
{
  CH_TIME("ParticleContainer<P>::regrid");

  // TLDR: This is the main regrid routine for particle container. This is essentially a ::define() followed by several types of remapping steps
  CH_assert(m_isDefined);

  if (m_isOrganizedByCell) {
    MayDay::Error("ParticleContainer<P>::regrid(...) - particles are sorted by cell!");
  }

  // Update this stuff
  m_grids       = a_grids;
  m_domains     = a_domains;
  m_refRat      = a_refRat;
  m_validRegion = a_validMask;
  m_finestLevel = a_newFinestLevel;

  m_dx.resize(1 + m_finestLevel);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    m_dx[lvl] = a_dx[lvl] * RealVect::Unit;
  }

  this->setupGrownGrids(a_lmin, m_finestLevel);
  this->setupParticleData(a_lmin, m_finestLevel);

  const int oldFinestLevel = m_cacheParticles.size() - 1;

  // If a level was removed, the length of m_cacheParticles exceeds the size of the new data holders. We need somewhere to put the particles
  // on the grid levels that were removed -- just put them on the new finest level.
  if (m_finestLevel < oldFinestLevel) {
    List<P>& cacheFinest = m_cacheParticles[m_finestLevel];
    for (int lvl = m_finestLevel + 1; lvl <= oldFinestLevel; lvl++) {
      List<P>& p = m_cacheParticles[lvl];

      cacheFinest.catenate(p);
      p.clear();
    }
  }

  // Add the particles from m_cacheParticles to each outcast list in m_particles. After we remap, many particles will be remapped to the correct level but the
  // outcast list is probably not empty after this. For example, we could have removed a grid patch on level 'lvl' and the particles in that grid patch will
  // still be in the outcast list on level 'lvl' after remapping.
  for (int lvl = 0; lvl <= a_newFinestLevel; lvl++) {
    List<P>& outcast = m_particles[lvl]->outcast();

    // Add particles from cache
    if (lvl <= oldFinestLevel) {
      outcast.catenate(m_cacheParticles[lvl]);
      m_cacheParticles[lvl].clear();
    }

    this->remapOutcast(*m_particles[lvl]);
  }

  // If grids were removed the fine particle outcast list is definitely not empty. Transfer fine particles
  // to the outcast list below and remap the coarse level again. Some of these particles may be outside the
  // coar level PVR, but that is taken care of later. After this, only the outcast list on lvl = 0 has particles.
  for (int lvl = a_newFinestLevel; lvl > 0; lvl--) {
    List<P>& fineParticles = m_particles[lvl]->outcast();
    List<P>& coarParticles = m_particles[lvl - 1]->outcast();

    coarParticles.catenate(fineParticles);
    fineParticles.clear();

    this->remapOutcast(*m_particles[lvl - 1]);
  }

  // At this point all the particles have been moved to their correct patches, but the PVR is so far not respected. We call ::remap to ensure that
  // particles that fall off their level's PVR are put in the correct place.
  this->remap();

  m_cacheParticles.resize(0);
}

template <class P>
template <Real& (P::*particleScalarField)()>
void
ParticleContainer<P>::setValue(const Real a_value)
{
  CH_TIME("ParticleContainer<P>::setValue");

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
  CH_TIME("ParticleContainer<P>::setValue");

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
  CH_TIME("ParticleContainer<P>::getNumberOfValidParticlesLocal");

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
  CH_TIME("ParticleContainer<P>::getNumberOfValidParticlesGlobal");

  CH_assert(m_isDefined);

  const unsigned long long n = this->getNumberOfValidParticlesLocal();

  return ParallelOps::sum(n);
}

template <class P>
unsigned long long
ParticleContainer<P>::getNumberOfOutcastParticlesLocal() const
{
  CH_TIME("ParticleContainer<P>::getNumberOfOutcastParticlesLocal");

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
  CH_TIME("ParticleContainer<P>::getNumberOfOutcastParticlesGlobal");

  CH_assert(m_isDefined);

  const unsigned long long n = this->getNumberOfOutcastParticlesLocal();

  return ParallelOps::sum(n);
}

template <class P>
unsigned long long
ParticleContainer<P>::getNumberOfMaskParticlesLocal() const
{
  CH_TIME("ParticleContainer<P>::getNumberOfMaskParticlesLocal");

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
  CH_TIME("ParticleContainer<P>::getNumberOfMaskParticlesGlobal");

  CH_assert(m_isDefined);

  const unsigned long long n = this->getNumberOfMaskParticlesLocal();

  return ParallelOps::sum(n);
}

template <class P>
void
ParticleContainer<P>::copyMaskParticles(const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_mask) const
{
  CH_TIME("ParticleContainer<P>::copyMaskParticles(amr)");

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
  CH_TIME("ParticleContainer<P>::copyMaskParticles(level)");

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

    const BaseFab<bool>& mask    = a_mask[din];
    const Box            gridBox = m_grids[a_level][din];
    const Box            maskBox = mask.box();

    CH_assert(gridBox == maskBox);

    List<P>&       maskParticles = (*m_maskParticles[a_level])[din].listItems();
    const List<P>& particles     = (*m_particles[a_level])[din].listItems();

    for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
      const RealVect x  = lit().position();
      const IntVect  iv = IntVect(D_DECL(std::floor((x[0] - m_probLo[0]) / dx[0]),
                                        std::floor((x[1] - m_probLo[1]) / dx[1]),
                                        std::floor((x[2] - m_probLo[2]) / dx[2])));
      if (!(maskBox.contains(iv))) {
        std::cout << a_level << "\t" << iv << "\t" << x << "\t" << maskBox << std::endl;

        MayDay::Error("ParticleContainer<P>::copyMaskParticles -- logic bust. Particle has fallen off grid");
      }
      else if (mask(iv)) {
        maskParticles.add(lit());
      }
    }
  }
}

template <class P>
void
ParticleContainer<P>::transferMaskParticles(const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_mask)
{
  CH_TIME("ParticleContainer<P>::transferMaskParticles(amr)");

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
  CH_TIME("ParticleContainer<P>::transferMaskParticles(level)");

  CH_assert(m_isDefined);

  const RealVect dx = m_dx[a_level];

  // Copy particles from m_particles to m_maskParticles if they lie in the input mask.
  const DisjointBoxLayout& dbl = m_grids[a_level];
  const DataIterator&      dit = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const BaseFab<bool>& mask    = a_mask[din];
    const Box            gridBox = m_grids[a_level][din];
    const Box            maskBox = mask.box();

    CH_assert(gridBox == maskBox);

    List<P>&       maskParticles = (*m_maskParticles[a_level])[din].listItems();
    const List<P>& particles     = (*m_particles[a_level])[din].listItems();

    for (ListIterator<P> lit(particles); lit.ok();) {
      const RealVect x  = lit().position();
      const IntVect  iv = IntVect(D_DECL(std::floor((x[0] - m_probLo[0]) / dx[0]),
                                        std::floor((x[1] - m_probLo[1]) / dx[1]),
                                        std::floor((x[2] - m_probLo[2]) / dx[2])));
      if (!(maskBox.contains(iv))) {
        //        std::cout << a_level << "\t" << iv << "\t" << x << std::endl;

        MayDay::Warning("ParticleContainer<P>::transferMaskParticles -- logic bust. Particle has fallen off grid");

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

template <class P>
void
ParticleContainer<P>::gatherOutcast(ParticleData<P>& a_particles) noexcept
{
  CH_TIME("ParticleContainer<P>::gatherOutcast");

#if 0
  a_particles.gatherOutcast();

  return;
#endif

  List<P>& outcasts = a_particles.outcast();

  const BoxLayout&   grids = a_particles.getBoxes();
  const DataIterator dit   = grids.dataIterator();
  const RealVect&    dx    = a_particles.meshSpacing();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(catenate : outcasts)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    ListBox<P>& particles = a_particles[din];

    particles.getInvalidDestructive(outcasts, grids[din]);
  }
}

template <class P>
void
ParticleContainer<P>::remapOutcast(ParticleData<P>& a_particles) noexcept
{
  CH_TIMERS("ParticleContainer::remapOutcast");
  CH_TIMER("ParticleContainer::remapOutcast::build_map", t1);
  CH_TIMER("ParticleContainer::remapOutcast::mpi_particle_loop", t3);
  CH_TIMER("ParticleContainer::remapOutcast::omp_copy_particles", t4);
  CH_TIMER("ParticleContainer::remapOutcast::omp_particle_loop", t5);
  CH_TIMER("ParticleContainer::remapOutcast::omp_catenate", t6);
  CH_TIMER("ParticleContainer::remapOutcast::local_map", t7);
  CH_TIMER("ParticleContainer::remapOutcast::mpi_scatter", t8);
  CH_TIMER("ParticleContainer::remapOutcast::assign", t9);

  // TODO:
  // 1. Thread properly with OpenMP and create local lists
  // 2. Use unordered_map instead? Insertion should be faster.
  // 3. Think about adding another routine which does the full grid at once.

  // Uncomment this if you want to default back to Chombo's native remapping routine
#if 1
  a_particles.remapOutcast();

  return;
#endif

  const unsigned int numRanks = numProc();
  const unsigned int myRank   = procID();
  const BoxLayout&   grids    = a_particles.getBoxes();

  // Storage for particles to be delivered to each box/proc. The vector holds a map of particles to be sent to
  // the indicated rank. The map itself indicates the box id and the list of particles to go into that box.
  std::vector<std::map<unsigned int, List<P>>> particlesToSend(numRanks);

  List<P>& outcasts = a_particles.outcast();

  if (!outcasts.isEmpty()) {

    // Note, dx is a pseudo-grid spacing.
    const int      fixedBoxSize = a_particles.fixedBoxSize();
    const RealVect origin       = a_particles.origin();
    const RealVect dx           = a_particles.meshSpacing() * fixedBoxSize;

    // tileSet = List of all boxes in a "tiled" space.
    // myTiles = List of this ranks boxes
    // otherTiles = List of other ranks' boxes. boxids contain the local global box index AND the rank
    std::set<IntVect, CompIntVect>               tileSet;
    std::map<IntVect, unsigned int, CompIntVect> myTiles;
    std::map<IntVect, boxids, CompIntVect>       otherTiles;

    // Build local mapping of the grids -- this figures out which box is owned by which
    // processor.
    //
    // Note: If this becomes expensive we can probably store this structure somewhere.
    CH_START(t1);
    for (LayoutIterator lit = grids.layoutIterator(); lit.ok(); ++lit) {
      const LayoutIndex lidx = lit();
      const IntVect     tile = coarsen(grids[lidx], fixedBoxSize).smallEnd();

      // If using MPI, need to figure out who owns which tiles.
#ifdef CH_MPI
      const unsigned int rankID = grids.procID(lidx);
      const unsigned int tileID = grids.index(lidx);

      if (rankID == myRank) {
        myTiles[tile] = tileID;
      }
      else {
        otherTiles[tile] = boxids(tileID, rankID);
      }
#else
      myTiles[tile]      = grids.index(lidx);
#endif

      tileSet.insert(tile);
    }
    CH_STOP(t1);

    // Figure out which particles will go where. This switches between non-OpenMP and OpenMP versions.
#ifndef _OPENMP
    CH_START(t3);
    for (ListIterator<P> lit(outcasts); lit.ok();) {
      const IntVect tile = locateBin(lit().position(), dx, origin);

      if (myTiles.find(tile) != myTiles.end()) {
        // This particle goes from this rank to this rank.
        const auto& gridIndex = myTiles[tile];
        const auto& toRank    = myRank;

        particlesToSend[toRank][gridIndex].transfer(lit);
      }
      else {
#ifdef CH_MPI
        if (otherTiles.find(tile) != otherTiles.end()) {
          // Particle goes from this rank to another rank.
          const auto& gridIndex = otherTiles[tile].idx;
          const auto& toRank    = otherTiles[tile].pid;

          particlesToSend[toRank][gridIndex].transfer(lit);
        }
        else {
          ++lit;
        }
#else
        ++lit;
#endif
      }
    }
    CH_STOP(t3);
#else
    // Copy particles and clear outcast list.
    CH_START(t4);
    std::vector<P> allOutcasts;
    allOutcasts.reserve(outcasts.length());
    for (ListIterator<P> lit(outcasts); lit.ok(); ++lit) {
      allOutcasts.emplace_back(lit());
    }

    outcasts.clear();
    CH_STOP(t4);

    const int numOutcasts = allOutcasts.size();

#pragma omp parallel
    {
      std::vector<std::map<unsigned int, List<P>>> particlesToSendPerThread(numRanks);

      CH_START(t5);
#pragma omp for schedule(runtime) reduction(catenate : outcasts)
      for (int ipart = 0; ipart < numOutcasts; ipart++) {
        const P& particle = allOutcasts[ipart];
        const IntVect tile = locateBin(particle.position(), dx, origin);

        if (myTiles.find(tile) != myTiles.end()) {
          // This particle goes from this rank to this rank.
          const auto& gridIndex = myTiles[tile];
          const auto& toRank = myRank;

          particlesToSendPerThread[toRank][gridIndex].add(particle);
        }
        else {
#ifdef CH_MPI
          if (otherTiles.find(tile) != otherTiles.end()) {
            // Particle goes from this rank to another rank.
            const auto& gridIndex = otherTiles[tile].idx;
            const auto& toRank = otherTiles[tile].pid;

            particlesToSendPerThread[toRank][gridIndex].add(particle);
          }
          else {
            outcasts.add(particle);
          }
#else
          outcasts.add(particle);
#endif
        }
      }
      CH_STOP(t5);

#pragma omp critical
      {
        CH_START(t6);
        // Append the thread-local particle to the global particles -- this is a bit more complicated than I
        // like because we can't just append the vectors.
        for (int irank = 0; irank < numRanks; irank++) {
          std::map<unsigned int, List<P>>& particlesToRank = particlesToSend[irank];

          for (auto& m : particlesToSendPerThread[irank]) {
            const unsigned int gridIndex = m.first;
            List<P>& threadParticles = m.second;

            particlesToRank[gridIndex].catenate(threadParticles);
          }
        }

        particlesToSendPerThread.resize(0);
      }
      CH_STOP(t6);
    }
#endif
  }

  // Particles received onto this rank from all ranks.
  std::map<unsigned int, List<P>> receivedParticles;

  // Collect local particles that didn't change rank ownership so we don't have to
  // communicate them.
  CH_START(t7);
  for (auto& localParticles : particlesToSend[myRank]) {
    receivedParticles[localParticles.first].catenate(localParticles.second);
  }
  particlesToSend[myRank].clear();
  CH_STOP(t7);

  // Distribute remnant particles.
#ifdef CH_MPI
  CH_START(t8);
  mpi_scatter_part(receivedParticles, particlesToSend);
  CH_STOP(t8);
#endif

  // Assign particles to boxes
  CH_START(t9);
  const DataIterator dit = grids.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex    din       = dit[mybox];
    const unsigned int gridIndex = grids.index(din);

    if (receivedParticles.find(gridIndex) != receivedParticles.end()) {
      List<P>& particles = a_particles[din].listItems();

      particles.catenate(receivedParticles[gridIndex]);
    }

    // Debug hook which double-checks that all particles are in the correct box.
    if (m_debug) {
      const Box      box       = grids[din];
      const List<P>& particles = a_particles[din].listItems();

      for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
        const RealVect position = lit().position();

        if (!(box.contains(ParticleOps::getParticleCellIndex(lit().position(), m_probLo, a_particles.meshSpacing())))) {
          MayDay::Error("ParticleContainer::remapOutcast -- particle is not contained in grid box!");
        }
      }
    }
  }
  CH_STOP(t9);
}

#include <CD_NamespaceFooter.H>

#endif
