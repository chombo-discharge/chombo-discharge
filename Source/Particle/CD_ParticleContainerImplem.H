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

// Chombo includes
#include <ParmParse.H>
#include <CH_Timer.H>

// Our includes
#include <CD_ParticleContainer.H>
#include <CD_ParallelOps.H>
#include <CD_Timer.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

template <class P>
ParticleContainer<P>::ParticleContainer()
{
  m_isDefined         = false;
  m_isOrganizedByCell = false;
  m_profile           = false;
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
    for (DataIterator dit(m_grids[lvl]); dit.ok(); ++dit) {
      List<P>& particles = (*m_particles[lvl])[dit()].listItems();

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
  CH_TIME("ParticleContainer<P>::organizeParticlesByPatch");

  CH_assert(m_isDefined);

  if (!m_isOrganizedByCell) {

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const DisjointBoxLayout& dbl = m_grids[lvl];

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        BinFab<P>& cellParticles = (*m_cellSortedParticles[lvl])[dit()];

        cellParticles.define(dbl[dit()], m_dx[lvl], m_probLo);
        cellParticles.addItemsDestructive((*m_particles[lvl])[dit()].listItems());
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

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        ListBox<P>& patchParticles = (*m_particles[lvl])[dit()];
        BinFab<P>&  cellParticles  = (*m_cellSortedParticles[lvl])[dit()];

        // Kernel which adds moves particles from the cell container to the patch container.
        auto kernel = [&](const IntVect& iv) -> void {
          patchParticles.addItemsDestructive(cellParticles(iv, comp));
        };

        BoxLoops::loop(dbl[dit()], kernel);

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
    m_particles[lvl]->remapOutcast();

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
    m_particles[lvl]->remapOutcast();

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

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        ListBox<P>&       myParticles    = (*m_particles[lvl])[dit()];
        const ListBox<P>& otherParticles = a_otherContainer[lvl][dit()];

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
      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        const List<P>& otherParticles = a_otherContainer[lvl][dit()].listItems();

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

      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {

        ListBox<P>& myParticles    = (*m_particles[lvl])[dit()];
        ListBox<P>& otherParticles = a_otherContainer[lvl][dit()];

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
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
        List<P>& otherParticles = a_otherContainer[lvl][dit()].listItems();

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
    for (DataIterator dit(m_grids[lvl]); dit.ok(); ++dit) {
      List<P>& myParticles    = (*m_particles[lvl])[dit()].listItems();
      List<P>& otherParticles = (*a_particles[lvl])[dit()].listItems();

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

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const BaseFab<bool>& mask = (*m_validRegion[a_level])[dit()];

    for (ListIterator<P> lit(a_particles[dit()].listItems()); lit.ok();) {

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

  CH_START(t1);
  List<P>& outcastFinest = m_particles[m_finestLevel]->outcast();
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    for (DataIterator dit(m_grids[lvl]); dit.ok(); ++dit) {
      outcastFinest.catenate((*m_particles[lvl])[dit()].listItems());
    }
  }
  CH_STOP(t1);

  for (int lvl = m_finestLevel; lvl >= 0; lvl--) {
    m_particles[lvl]->gatherOutcast();
    m_particles[lvl]->remapOutcast();

    CH_START(t2);
    if (lvl > 0) {

      List<P>& outcast     = m_particles[lvl]->outcast();
      List<P>& outcastCoar = m_particles[lvl - 1]->outcast();

      outcastCoar.catenate(outcast);
    }
    CH_STOP(t2);
  }
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

  m_particles[a_lvl]->gatherOutcast();
  m_particles[a_lvl]->remapOutcast();
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
    m_cacheParticles[lvl].clear();

    for (DataIterator dit(m_grids[lvl]); dit.ok(); ++dit) {
      List<P>& p = (*m_particles[lvl])[dit()].listItems();

      m_cacheParticles[lvl].catenate(p);

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

    m_particles[lvl]->remapOutcast();
  }

  // If grids were removed the fine particle outcast list is definitely not empty. Transfer fine particles
  // to the outcast list below and remap the coarse level again. Some of these particles may be outside the
  // coar level PVR, but that is taken care of later. After this, only the outcast list on lvl = 0 has particles.
  for (int lvl = a_newFinestLevel; lvl > 0; lvl--) {
    List<P>& fineParticles = m_particles[lvl]->outcast();
    List<P>& coarParticles = m_particles[lvl - 1]->outcast();

    coarParticles.catenate(fineParticles);
    fineParticles.clear();

    m_particles[lvl - 1]->remapOutcast();
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

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<P>& patchParticles = (*m_particles[lvl])[dit()].listItems();

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

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      List<P>& patchParticles = (*m_particles[lvl])[dit()].listItems();

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
    for (DataIterator dit(m_grids[lvl]); dit.ok(); ++dit) {
      const List<P>& particles = (*m_particles[lvl])[dit()].listItems();

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
    for (DataIterator dit(m_grids[lvl]); dit.ok(); ++dit) {
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
  CH_TIME("ParticleContainer<P>::getNumberOfMaskparticlesLocal");

  CH_assert(m_isDefined);

  unsigned long long n = 0;

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    for (DataIterator dit(m_grids[lvl]); dit.ok(); ++dit) {
      const List<P>& maskParticles = (*m_maskParticles[lvl])[dit()].listItems();

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

  // 1. Copy particles from m_particles to m_maskParticles if they lie in the input mask.
  for (DataIterator dit(m_grids[a_level]); dit.ok(); ++dit) {
    const BaseFab<bool>& mask    = a_mask[dit()];
    const Box            gridBox = m_grids[a_level][dit()];
    const Box            maskBox = mask.box();

    CH_assert(gridBox == maskBox);

    List<P>&       maskParticles = (*m_maskParticles[a_level])[dit()].listItems();
    const List<P>& particles     = (*m_particles[a_level])[dit()].listItems();

    for (ListIterator<P> lit(particles); lit.ok(); ++lit) {
      const RealVect x  = lit().position();
      const IntVect  iv = IntVect(D_DECL(std::floor((x[0] - m_probLo[0]) / dx[0]),
                                        std::floor((x[1] - m_probLo[1]) / dx[1]),
                                        std::floor((x[2] - m_probLo[2]) / dx[2])));
      if (!(maskBox.contains(iv))) {
        MayDay::Warning("ParticleContainer<P>::copyMaskParticles -- logic bust. Particle has fallen off grid");
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

  // 1. Copy particles from m_particles to m_maskParticles if they lie in the input mask.
  for (DataIterator dit(m_grids[a_level]); dit.ok(); ++dit) {
    const BaseFab<bool>& mask    = a_mask[dit()];
    const Box            gridBox = m_grids[a_level][dit()];
    const Box            maskBox = mask.box();

    CH_assert(gridBox == maskBox);

    List<P>&       maskParticles = (*m_maskParticles[a_level])[dit()].listItems();
    const List<P>& particles     = (*m_particles[a_level])[dit()].listItems();

    for (ListIterator<P> lit(particles); lit.ok();) {
      const RealVect x  = lit().position();
      const IntVect  iv = IntVect(D_DECL(std::floor((x[0] - m_probLo[0]) / dx[0]),
                                        std::floor((x[1] - m_probLo[1]) / dx[1]),
                                        std::floor((x[2] - m_probLo[2]) / dx[2])));
      if (!(maskBox.contains(iv))) {
        MayDay::Warning("ParticleContainer<P>::copyMaskParticles -- logic bust. Particle has fallen off grid");
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
ParticleContainer<P>::clear(AMRParticles<P>& a_particles) const
{
  CH_assert(m_isDefined);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    a_particles[lvl]->clear();
  }
}

#include <CD_NamespaceFooter.H>

#endif
