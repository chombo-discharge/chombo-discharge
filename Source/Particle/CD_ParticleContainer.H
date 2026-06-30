/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ParticleContainer.H
  @brief  Declaration of an AMR-hierarchy container that owns per-patch ParticleSoA leaves.
  @author Robert Marskar
*/

#ifndef CD_PARTICLECONTAINER_H
#define CD_PARTICLECONTAINER_H

// Std includes
#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

// Chombo includes
#include <REAL.H>
#include <SPACE.H>
#include <RealVect.H>
#include <IntVect.H>
#include <Box.H>
#include <ProblemDomain.H>
#include <DisjointBoxLayout.H>
#include <BoxLayout.H>
#include <DataIterator.H>
#include <LayoutData.H>
#include <LevelData.H>
#include <BaseFab.H>
#include <RefCountedPtr.H>
#include <Vector.H>
#include <SPMD.H>
#include <CH_Timer.H>
#include <CH_assert.H>

// Our includes
#include <CD_LevelTiles.H>
#include <CD_ParallelOps.H>
#include <CD_ParticleSoA.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Per-level holder vector: one ParticleSoA leaf container per box, per level.
  @tparam P      Payload type.
  @tparam Traits Column descriptor.
*/
template <typename P, typename Traits>
using AMRParticlesSoA = Vector<RefCountedPtr<LayoutData<ParticleSoA<P, Traits>>>>;

/**
  @brief AMR-hierarchy container of computational particles, stored per patch in Struct-of-Arrays form.
  @details
  ParticleContainer owns particles distributed across an AMR hierarchy. For every level it holds one
  @c LayoutData<ParticleSoA<P>> indexed by @c DataIndex over that level's @c DisjointBoxLayout, so the
  particles of a single patch are reached as @c container[lvl][din] and stored contiguously in
  Struct-of-Arrays form (see ParticleSoA). Each particle carries a position and a weight (owned by the
  leaf) plus the user payload @c P, which contributes the extra per-particle columns.

  <h3>Particle holders</h3>
  Three parallel sets of per-level leaves are maintained:
  - <b>valid</b>  -- the live particles, each owned by the patch whose grid cell contains it. Reached
                     through operator[] and getParticles().
  - <b>mask</b>   -- a halo selection copied or moved out of the valid particles by a per-cell mask,
                     used for halo deposition near refinement boundaries (getMaskParticles()).
  - <b>buffer</b> -- particles living on the grown grids (each box grown by the refinement ratio), used
                     as scratch when depositing across coarse-fine boundaries (getBufferParticles()).

  <h3>Ownership and lifetime</h3>
  The per-level holders are shared @c RefCountedPtr, so the container is move-only: copy construction and
  copy assignment are deleted (a copy would alias the holders). Construct an empty container, then call
  define() with the grids and AMR metadata before any other use; isDefined() reports readiness.

  <h3>Particle ownership and routing</h3>
  A particle belongs to the patch on the finest level whose cell contains its position. After particles
  move, remap() routes each one back to its owning patch -- across MPI ranks and across levels as needed.
  On a regrid, call preRegrid() before the grids change and regrid() afterwards to carry the particles
  onto the new layout. Particles that leave the domain are dropped and reported by
  getNumberOfOutcastParticlesLocal() / getNumberOfOutcastParticlesGlobal().

  <h3>Cell sorting</h3>
  A leaf stores its particles in arbitrary (append) order by default. organizeParticlesByCell() sorts
  each valid leaf so that the particles of a cell occupy a contiguous range (queried through the leaf's
  CSR cell offsets); organizeParticlesByPatch() returns to the unsorted view. isOrganizedByCell() reports
  the current state.

  <h4>Typical usage</h4>
  @code{.cpp}
    ParticleContainer<MyPayload> particles;
    particles.define(grids, domains, dx, refRat, probLo, minBlockSize, levelTiles, finestLevel, "primal", validCells);

    // Add particles locally, then route every particle to its owning patch (and across ranks).
    particles.addParticlesLocal(lvl, pos, weight, MyPayload{...});
    particles.remap();

    // Iterate the valid particles of one patch.
    const DisjointBoxLayout& dbl  = particles.getGrids()[lvl];
    const DataIterator&      dit  = dbl.dataIterator();
    const int                nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex&        din  = dit[mybox];
      ParticleSoA<MyPayload>& leaf = particles[lvl][din];
      for (std::size_t i = 0; i < leaf.size(); i++) {
        const RealVect x = leaf.position(i);
        // ... use x, leaf.weight(i), payload columns ...
      }
    }
  @endcode

  @tparam P      User payload struct contributing the per-particle payload columns; defaults to NoPayload.
  @tparam Traits Column descriptor for the payload; defaults to ParticleTraits<P>.
*/
template <typename P = NoPayload, typename Traits = ParticleTraits<P>>
class ParticleContainer
{
public:
  /**
    @brief The per-patch leaf type owned at every (level, box).
  */
  using Leaf = ParticleSoA<P, Traits>;

  /**
    @brief The per-level holder type (one Leaf per box on the level's DisjointBoxLayout).
  */
  using LevelParticles = LayoutData<Leaf>;

  /**
    @name Construction and lifetime
  */
  ///@{

  /**
    @brief Construct an empty, undefined container. Call define() before use.
  */
  ParticleContainer()
    : m_probLo(RealVect::Zero),
      m_minBlockSize(0),
      m_finestLevel(-1),
      m_numOutcastLocal(0),
      m_isDefined(false),
      m_isOrganizedByCell(false)
  {}

  /**
    @brief Destructor.
  */
  ~ParticleContainer() = default;

  /**
    @brief Copy construction is deleted (holders are shared RefCountedPtr -- a copy would alias).
  */
  ParticleContainer(const ParticleContainer&) = delete;

  /**
    @brief Copy assignment is deleted (move-only ownership).
    @return Reference to this.
  */
  ParticleContainer&
  operator=(const ParticleContainer&) = delete;

  /**
    @brief Allocate the per-level holders and the per-level tile-ownership maps over the AMR grids.
    @details Allocates one empty Leaf per box for the valid, mask, and buffer holders on every level in
    [0, a_finestLevel], and builds a LevelTiles per level (from the grids + blocking factor) used by
    remap() to map a particle to its owning patch/rank. The mask and buffer holders start empty; they are
    populated on demand by copyMaskParticles()/transferMaskParticles() and the deposition routines. Must
    be called once before any other operation.
    @param[in] a_grids          Per-level DisjointBoxLayouts.
    @param[in] a_domains        Per-level problem domains.
    @param[in] a_dx             Per-level (isotropic) grid spacing.
    @param[in] a_refRat         Per-level refinement ratio (entry l is the ratio between level l and l+1).
    @param[in] a_probLo         Lower-left corner of the physical domain.
    @param[in] a_minBlockSize Grid blocking factor (tile size); the grids must be unions of
                                minBlockSize-sized boxes for the tile ownership to be exact.
    @param[in] a_numGhost       Width (in cells, at each destination patch's own level) of the ghost
                                layer filled by fillGhostParticles(). Zero disables ghost filling.
    @param[in] a_levelTiles     The owning Realm's per-level tile->box maps, aliased (never built locally)
                                as the patch-locating structure used by remap(). Must correspond to
                                a_grids. Not owned; must outlive this container.
    @param[in] a_finestLevel    Finest AMR level index.
    @param[in] a_realm          Realm label (diagnostics only).
    @param[in] a_validCells     Optional alias to the realm's valid-cell mask (true = cell not covered
                                by a finer level), enabling remap()'s stayer fast path. Pass nullptr to
                                always re-pool every particle. Not owned; must outlive this container.
  */
  void
  define(const Vector<DisjointBoxLayout>&                       a_grids,
         const Vector<ProblemDomain>&                           a_domains,
         const Vector<Real>&                                    a_dx,
         const Vector<int>&                                     a_refRat,
         const RealVect&                                        a_probLo,
         const int                                              a_minBlockSize,
         const int                                              a_numGhost,
         const Vector<RefCountedPtr<LevelTiles>>&               a_levelTiles,
         const int                                              a_finestLevel,
         const std::string&                                     a_realm,
         const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>* a_validCells)
  {
    CH_TIME("ParticleContainer::define");

    CH_assert(a_grids.size() > a_finestLevel);
    CH_assert(a_domains.size() > a_finestLevel);
    CH_assert(a_dx.size() > a_finestLevel);
    CH_assert(a_refRat.size() > a_finestLevel);
    CH_assert(a_levelTiles.size() > a_finestLevel);

    m_grids           = a_grids;
    m_domains         = a_domains;
    m_refRat          = a_refRat;
    m_probLo          = a_probLo;
    m_minBlockSize    = a_minBlockSize;
    m_numGhost        = a_numGhost;
    m_finestLevel     = a_finestLevel;
    m_realm           = a_realm;
    m_validCells      = a_validCells;
    m_numOutcastLocal = 0;

    m_dx.resize(a_finestLevel + 1);
    for (int lvl = 0; lvl <= a_finestLevel; lvl++) {
      m_dx[lvl] = a_dx[lvl] * RealVect::Unit;
    }

    // Alias the Realm's per-level tile->box maps (the single source of truth); never build our own.
    m_levelTiles = a_levelTiles;

    m_particles.resize(a_finestLevel + 1);
    m_maskParticles.resize(a_finestLevel + 1);
    for (int lvl = 0; lvl <= a_finestLevel; lvl++) {
      m_particles[lvl]     = RefCountedPtr<LevelParticles>(new LevelParticles(a_grids[lvl]));
      m_maskParticles[lvl] = RefCountedPtr<LevelParticles>(new LevelParticles(a_grids[lvl]));
    }

    // Buffer particles live on the grown grids (boxes grown by the refinement factor on finer levels).
    this->setupGrownGrids();
    m_bufferParticles.resize(a_finestLevel + 1);
    for (int lvl = 0; lvl <= a_finestLevel; lvl++) {
      m_bufferParticles[lvl] = RefCountedPtr<LevelParticles>(new LevelParticles(m_grownGrids[lvl]));
    }

    m_isOrganizedByCell = false;
    m_isDefined         = true;
  }

  ///@}

  /**
    @name Grid metadata
  */
  ///@{

  /**
    @brief Per-level AMR grids.
    @return The DisjointBoxLayouts, one per level.
  */
  const Vector<DisjointBoxLayout>&
  getGrids() const
  {
    return m_grids;
  }

  /**
    @brief Per-level grid spacing.
    @return The (isotropic) grid spacing per level, as RealVects.
  */
  Vector<RealVect>
  getDx() const
  {
    return m_dx;
  }

  /**
    @brief Lower-left corner of the physical domain.
    @return The probLo corner.
  */
  RealVect
  getProbLo() const
  {
    return m_probLo;
  }

  /**
    @brief Finest AMR level index.
    @return The finest level.
  */
  int
  getFinestLevel() const
  {
    return m_finestLevel;
  }

  /**
    @brief Realm label.
    @return The realm string.
  */
  std::string
  getRealm() const
  {
    return m_realm;
  }

  ///@}

  /**
    @name Particle holders
  */
  ///@{

  /**
    @brief The valid particles on all levels.
    @return Mutable per-level holder vector.
  */
  AMRParticlesSoA<P, Traits>&
  getParticles()
  {
    return m_particles;
  }

  /**
    @brief The valid particles on all levels (const).
    @return Const per-level holder vector.
  */
  const AMRParticlesSoA<P, Traits>&
  getParticles() const
  {
    return m_particles;
  }

  /**
    @brief The halo/mask particles on all levels.
    @return Mutable per-level holder vector.
  */
  AMRParticlesSoA<P, Traits>&
  getMaskParticles()
  {
    return m_maskParticles;
  }

  /**
    @brief The halo/mask particles on all levels (const).
    @return Const per-level holder vector.
  */
  const AMRParticlesSoA<P, Traits>&
  getMaskParticles() const
  {
    return m_maskParticles;
  }

  /**
    @brief The grown-grid buffer particles on all levels.
    @return Mutable per-level holder vector.
  */
  AMRParticlesSoA<P, Traits>&
  getBufferParticles()
  {
    return m_bufferParticles;
  }

  /**
    @brief The grown-grid buffer particles on all levels (const).
    @return Const per-level holder vector.
  */
  const AMRParticlesSoA<P, Traits>&
  getBufferParticles() const
  {
    return m_bufferParticles;
  }

  /**
    @brief Valid particles on one level (index by DataIndex to reach a patch's Leaf).
    @details Chain a second index by @c DataIndex to reach the patch's ParticleSoA leaf:
    @code{.cpp}
      ParticleSoA<P>& leaf = particles[lvl][din];
    @endcode
    @param[in] a_lvl Level index.
    @return Mutable per-level holder; operator[](DataIndex) yields the patch's ParticleSoA.
  */
  LevelParticles&
  operator[](const int a_lvl)
  {
    CH_assert(m_isDefined);
    CH_assert(a_lvl >= 0 && a_lvl <= m_finestLevel);

    return *m_particles[a_lvl];
  }

  /**
    @brief Valid particles on one level (const).
    @param[in] a_lvl Level index.
    @return Const per-level holder.
  */
  const LevelParticles&
  operator[](const int a_lvl) const
  {
    CH_assert(m_isDefined);
    CH_assert(a_lvl >= 0 && a_lvl <= m_finestLevel);

    return *m_particles[a_lvl];
  }

  ///@}

  /**
    @name State
  */
  ///@{

  /**
    @brief Whether the valid leaves are currently cell-sorted.
    @return True if organizeParticlesByCell() is in effect, false after organizeParticlesByPatch().
  */
  bool
  isOrganizedByCell() const
  {
    return m_isOrganizedByCell;
  }

  /**
    @brief Whether define() has been called.
    @return True iff defined.
  */
  bool
  isDefined() const
  {
    return m_isDefined;
  }

  ///@}

  /**
    @name Population
  */
  ///@{

  /**
    @brief Drop all valid particles on every level (keeps each leaf's arena capacity).
  */
  void
  clearParticles()
  {
    CH_TIME("ParticleContainer::clearParticles");

    CH_assert(m_isDefined);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      LevelParticles&          level = *m_particles[lvl];
      const DisjointBoxLayout& dbl   = m_grids[lvl];
      const DataIterator&      dit   = dbl.dataIterator();
      const int                nbox  = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        level[din].clear();
      }
    }
  }

  /**
    @brief Add a free-standing buffer of particles to the container, routing each to its owner.
    @details The buffer (e.g. filled by ParticleManagement::drawRandomParticles(ParticleSoA)) is catenated
    onto a local leaf -- O(1) when that leaf is empty -- and then remap() routes every particle (including
    across ranks and levels) to the patch/level/rank that owns its position. a_particles is left empty.
    Assumes this rank owns at least one box (true whenever level 0 covers the domain). Unlike
    addParticlesLocal(), the positions need not already belong to this rank.
    @code{.cpp}
      ParticleSoA<P> drawn;
      ParticleManagement::drawRandomParticles(drawn, ...);
      particles.addParticlesDestructive(drawn); // drawn is empty afterwards; particles are routed to owners
    @endcode
    @param[in,out] a_particles Buffer of particles to ingest (emptied on return).
  */
  void
  addParticlesDestructive(ParticleSoA<P, Traits>& a_particles)
  {
    CH_TIME("ParticleContainer::addParticlesDestructive");

    CH_assert(m_isDefined);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      DataIterator dit = m_grids[lvl].dataIterator();
      if (dit.ok()) {
        (*m_particles[lvl])[dit()].catenate(a_particles); // move the buffer onto a local leaf
        break;
      }
    }

    this->remap(); // route every ingested particle to its owning patch/level/rank
  }

  ///@}

  /**
    @name Counting
  */
  ///@{

  /**
    @brief Number of valid particles owned by this rank.
    @details Ghost particles (any non-Valid GhostType) are excluded -- only owned, valid particles
    are counted, even if a transient ghost halo is currently present in the leaves.
    @return Local valid particle count summed over all levels and patches.
  */
  unsigned long long
  getNumberOfValidParticlesLocal() const
  {
    CH_TIME("ParticleContainer::getNumberOfValidParticlesLocal");

    CH_assert(m_isDefined);

    unsigned long long n = 0;
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const LevelParticles&    level = *m_particles[lvl];
      const DisjointBoxLayout& dbl   = m_grids[lvl];
      const DataIterator&      dit   = dbl.dataIterator();
      const int                nbox  = dit.size();

#pragma omp parallel for schedule(runtime) reduction(+ : n)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din  = dit[mybox];
        const Leaf&      leaf = level[din];

        for (std::size_t i = 0; i < leaf.size(); i++) {
          if (!leaf.isGhost(i)) {
            n++;
          }
        }
      }
    }

    return n;
  }

  /**
    @brief Number of valid particles across all ranks.
    @return Global valid particle count (MPI all-reduce of the local count; == local in serial).
  */
  unsigned long long
  getNumberOfValidParticlesGlobal() const
  {
    CH_TIME("ParticleContainer::getNumberOfValidParticlesGlobal");

    CH_assert(m_isDefined);

    return ParallelOps::sum(this->getNumberOfValidParticlesLocal());
  }

  /**
    @brief Number of particles dropped off-domain by the most recent remap() on this rank.
    @return Local outcast count from the last remap.
  */
  unsigned long long
  getNumberOfOutcastParticlesLocal() const
  {
    return m_numOutcastLocal;
  }

  /**
    @brief Number of particles dropped off-domain by the most recent remap(), across all ranks.
    @return Global outcast count from the last remap (== local in serial).
  */
  unsigned long long
  getNumberOfOutcastParticlesGlobal() const
  {
    return ParallelOps::sum(m_numOutcastLocal);
  }

  ///@}

  /**
    @name Dynamics
  */
  ///@{

  /**
    @brief Redistribute every valid particle to the patch/level/rank that owns its cell.
    @details Pool model: collect all valid particles, map each to its destination via the
    finest-level tile that contains it (LevelTiles), append same-rank movers into the destination
    leaf and scatter cross-rank movers with MPI. Particles whose cell is not owned by any tile on
    any level (off-domain) are dropped and counted in getNumberOfOutcastParticles*(). Invalidates
    any cell-sort. Particle ids are preserved; rankID is set to the new owning rank.
  */
  void
  remap();

  /**
    @brief Cache the current valid particles ahead of a regrid.
    @details Snapshots the current per-level holders (and their layout) so that regrid() can
    rebuild over a new layout and redistribute the cached particles. Call exactly once before
    regrid().
  */
  void
  preRegrid();

  /**
    @brief Rebuild over a new AMR layout and redistribute the preRegrid()-cached particles onto it.
    @details Re-allocates the per-level holders and tile maps over the new grids (all levels), then
    runs the same pool->map->scatter as remap() but sourced from the cache. Positions are unchanged;
    each cached particle is routed to the patch/level/rank owning its cell on the NEW layout, with
    ids preserved. Off-domain particles are counted in getNumberOfOutcastParticles*(). Must follow a
    preRegrid().
    @param[in] a_grids          New per-level DisjointBoxLayouts.
    @param[in] a_domains        New per-level problem domains.
    @param[in] a_dx             New per-level (isotropic) grid spacing.
    @param[in] a_refRat         New per-level refinement ratios.
    @param[in] a_minBlockSize Grid blocking factor (tile size) for the new grids.
    @param[in] a_levelTiles     The owning Realm's per-level tile->box maps for the new grids, aliased
                                (never built locally). Must correspond to a_grids.
    @param[in] a_newFinestLevel New finest AMR level index.
  */
  void
  regrid(const Vector<DisjointBoxLayout>&         a_grids,
         const Vector<ProblemDomain>&             a_domains,
         const Vector<Real>&                      a_dx,
         const Vector<int>&                       a_refRat,
         const int                                a_minBlockSize,
         const Vector<RefCountedPtr<LevelTiles>>& a_levelTiles,
         const int                                a_newFinestLevel);

  ///@}

  /**
    @name Ghost particles
  */
  ///@{

  /**
    @brief Fill each patch's ghost halo with copies of particles that overlap its ghosted box.
    @details A ghost particle is a transient, non-owned copy of a particle that lives in an adjacent
    region (a neighbouring patch on the same level, the next-coarser level, or the next-finer level)
    and overlaps this patch's ghosted/grown box. The copies are appended into the SAME leaf as the
    valid particles and tagged with the appropriate GhostType (SameLevel/Coarse/Fine) so consumers
    can tell them from the owned particles via ParticleSoA::isGhost(). Ghosts are transient: this
    routine first calls clearGhostParticles() so it is safe to refill, and remap()/preRegrid()/
    regrid() drop ghosts so they are never routed or cached as if valid.
    @note The cross-patch/cross-level routing is not yet implemented -- this currently only clears any
    existing ghost halo. See the Phase-2 remap work.
  */
  void
  fillGhostParticles();

  /**
    @brief Remove every ghost particle (any non-Valid GhostType) from all levels and patches.
    @details Owned (Valid) particles are left untouched. Removal is by swap-and-pop, so particle
    order within a leaf is not preserved. Keeps leaf capacity.
  */
  void
  clearGhostParticles();

  ///@}

  /**
    @name Halo and buffer particles
  */
  ///@{

  /**
    @brief Copy the valid particles selected by a per-cell mask into the mask holder, all levels.
    @param[in] a_mask Per-level cell masks (skipped where null).
  */
  void
  copyMaskParticles(const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_mask);

  /**
    @brief Copy the valid particles whose cell the mask marks true into the mask holder, one level.
    @details Clears the level's mask leaves, then COPIES (keeps in the valid holder) every valid
    particle whose cell is true in the mask. Ids/rank are preserved.
    @param[in] a_level Level to fill.
    @param[in] a_mask  Per-cell boolean mask over that level's valid grids (box == grid box).
  */
  void
  copyMaskParticles(const int a_level, const LevelData<BaseFab<bool>>& a_mask);

  /**
    @brief Move the valid particles selected by a per-cell mask into the mask holder, all levels.
    @param[in] a_mask Per-level cell masks (skipped where null).
  */
  void
  transferMaskParticles(const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_mask);

  /**
    @brief Move the valid particles whose cell the mask marks true into the mask holder, one level.
    @details Like copyMaskParticles but MOVES (removes from the valid holder via swap-and-pop, which
    does not preserve order). Does not clear the mask holder first (accumulates). Ids/rank preserved.
    @param[in] a_level Level to fill.
    @param[in] a_mask  Per-cell boolean mask over that level's valid grids (box == grid box).
  */
  void
  transferMaskParticles(const int a_level, const LevelData<BaseFab<bool>>& a_mask);

  /**
    @brief Drop all mask/halo particles on every level (keeps capacity).
  */
  void
  clearMaskParticles()
  {
    CH_TIME("ParticleContainer::clearMaskParticles");

    CH_assert(m_isDefined);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      LevelParticles&          level = *m_maskParticles[lvl];
      const DisjointBoxLayout& dbl   = m_grids[lvl];
      const DataIterator&      dit   = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        level[din].clear();
      }
    }
  }

  /**
    @brief Drop all buffer particles on every level (keeps capacity).
  */
  void
  clearBufferParticles()
  {
    CH_TIME("ParticleContainer::clearBufferParticles");

    CH_assert(m_isDefined);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      LevelParticles&     level = *m_bufferParticles[lvl];
      const BoxLayout&    dbl   = m_grownGrids[lvl];
      const DataIterator& dit   = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        level[din].clear();
      }
    }
  }

  /**
    @brief Per-level grown grids that the buffer particles live on.
    @return The grown BoxLayouts (boxes grown by the refinement factor on finer levels).
  */
  const Vector<BoxLayout>&
  getGrownGrids() const
  {
    return m_grownGrids;
  }

  /**
    @brief Move all particles from another holder (on the same valid grids) into the valid holder.
    @details Per level/patch, catenates the source leaf into the valid leaf (a per-column bulk move
    preserving ids/rank, or an O(1) arena swap when the valid leaf is empty), leaving the source
    empty. Used to return halo/mask particles to the valid holder after a HaloNGP deposition. Assumes
    a_source is defined over this container's valid grids.
    @param[in,out] a_source Holder to drain into the valid particles.
  */
  void
  transferParticles(AMRParticlesSoA<P, Traits>& a_source)
  {
    CH_TIME("ParticleContainer::transferParticles");

    CH_assert(m_isDefined);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      LevelParticles& src = *a_source[lvl];
      LevelParticles& dst = *m_particles[lvl];

      const DisjointBoxLayout& dbl = m_grids[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        dst[din].catenate(src[din]); // bulk move (id/rank preserved); empty dst -> O(1) swap
      }
    }

    m_isOrganizedByCell = false;
  }

  ///@}

  /**
    @name Cell sorting
  */
  ///@{

  /**
    @brief Cell-sort every valid leaf (counting sort into Fortran cell order + CSR cell offsets).
    @details After this, each leaf is `isSorted()` and cell c owns the contiguous particle range
    [`cellStart(c)`, `cellStart(c+1)`) (query via the leaf accessors). Sets isOrganizedByCell().
  */
  void
  organizeParticlesByCell()
  {
    CH_TIME("ParticleContainer::organizeParticlesByCell");

    CH_assert(m_isDefined);

    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      LevelParticles& level = *m_particles[lvl];

      const DisjointBoxLayout& dbl = m_grids[lvl];
      const DataIterator&      dit = dbl.dataIterator();

      const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex& din = dit[mybox];

        level[din].sortByCell(dbl[din], m_dx[lvl], m_probLo);
      }
    }

    m_isOrganizedByCell = true;
  }

  /**
    @brief Return to a by-patch view (the leaves remain a single contiguous array per patch).
    @details Clears isOrganizedByCell(); the CSR cell offsets are no longer to be relied upon. There is no
    data movement -- the leaf already holds all particles of a patch contiguously.
  */
  void
  organizeParticlesByPatch()
  {
    CH_TIME("ParticleContainer::organizeParticlesByPatch");

    m_isOrganizedByCell = false;
  }

  ///@}

protected:
  /**
    @brief Build the grown grids (buffer-particle layout): boxes grown by the refinement factor.
  */
  void
  setupGrownGrids();

  /**
    @brief Pool key identifying a destination patch: (level, grid index within the level).
  */
  using PoolKey = std::pair<int, unsigned int>;

  /**
    @brief Per-rank pool of movers: pool[rank] maps a destination patch to a small SoA of particles.
  */
  using MoverPool = std::vector<std::map<PoolKey, Leaf>>;

  /**
    @brief Gather every particle of a source holder into a per-rank, per-destination-patch pool.
    @details Reads each particle, resolves its owning (level, grid index, rank) on the CURRENT tile
    maps, and copies it into the pool (ids preserved, rankID set to the new owner). Off-domain
    particles are counted into m_numOutcastLocal (which this resets).
    @param[in] a_source      Source holders (e.g. the valid particles, or the regrid cache).
    @param[in] a_sourceGrids The DisjointBoxLayouts the source holders are defined on.
    @return The per-rank mover pool.
  */
  MoverPool
  gatherToPool(const AMRParticlesSoA<P, Traits>& a_source, const Vector<DisjointBoxLayout>& a_sourceGrids);

  /**
    @brief Mask fast path for remap(): pool ONLY the movers, leaving stayers in place.
    @details Scans m_particles in place; a particle still inside its box AND in a cell flagged valid by
    m_validCells (not covered by a finer level) is a stayer and is left untouched. Every other particle
    is a mover: it is copied into the returned pool (ids preserved, rankID set to the new owner) and
    swap-popped out of its leaf. Off-domain movers are counted into m_numOutcastLocal (reset here).
    Requires m_validCells to be non-null and sized to m_finestLevel+1 (the caller checks this).
    @return The per-rank mover pool; stayers remain in m_particles (so distributeFromPool must append,
    not assume empty destinations, and remap() must NOT clearParticles() on this path).
  */
  MoverPool
  gatherMoversToPool();

  /**
    @brief Gather ghost copies of every valid particle into a per-rank, per-destination-patch pool.
    @details For every valid particle, finds all patches on levels [l-1, l, l+1] whose ghost layer
    (grow(box, m_numGhost), clamped to the domain) contains the particle's cell but whose own box does
    not -- i.e. the particle lies in that patch's ghost shell. A copy is binned for each such patch with
    the appropriate GhostType (SameLevel/Coarse/Fine from the destination's view) and the source's
    particleID and rankID PRESERVED (a ghost records its true owner, not the receiver). A particle may
    be a ghost of several patches, so it can be binned more than once. Non-periodic.
    @return The per-rank ghost pool (distribute it with a_setOwnerToReceiver == false).
  */
  MoverPool
  gatherGhostsToPool();

  /**
    @brief Assign a pool into the valid holders: local appends + MPI scatter.
    @param[in,out] a_pool              Per-rank pool (consumed). Same-rank bins are catenated into the
                                       destination leaf; cross-rank bins are linearized and Alltoallv'd.
    @param[in] a_setOwnerToReceiver    If true, every cross-rank particle's rankID is set to the
                                       receiving rank (remap/regrid ownership). If false, the rankID
                                       carried in the pool is preserved (ghost particles keep their owner).
  */
  void
  distributeFromPool(MoverPool& a_pool, const bool a_setOwnerToReceiver);

  /**
    @brief Per-level AMR grids.
  */
  Vector<DisjointBoxLayout> m_grids;

  /**
    @brief Per-level problem domains.
  */
  Vector<ProblemDomain> m_domains;

  /**
    @brief Per-level grid spacing (isotropic, stored as RealVect).
  */
  Vector<RealVect> m_dx;

  /**
    @brief Per-level refinement ratios (entry l is the ratio between level l and l+1).
  */
  Vector<int> m_refRat;

  /**
    @brief Per-level tile-ownership maps (cell tile -> owning grid index / rank), used by remap().
  */
  Vector<RefCountedPtr<LevelTiles>> m_levelTiles;

  /**
    @brief Lower-left corner of the physical domain.
  */
  RealVect m_probLo;

  /**
    @brief Grid blocking factor (tile size).
  */
  int m_minBlockSize;

  /**
    @brief Width (in cells, at each patch's own level) of the ghost layer filled by fillGhostParticles().
  */
  int m_numGhost;

  /**
    @brief Finest AMR level index.
  */
  int m_finestLevel;

  /**
    @brief Particles dropped off-domain by the most recent remap() (local count).
  */
  unsigned long long m_numOutcastLocal;

  /**
    @brief Realm label (diagnostics only).
  */
  std::string m_realm;

  /**
    @brief Optional alias to the realm's valid-cell mask (true = cell not covered by a finer level).
    @details Used by remap() to leave "stayer" particles in place: a particle still inside its box and
    in a cell not covered by a finer grid cannot change owner, so it need not be re-pooled. Null (or a
    stale size) makes remap() fall back to re-pooling every particle. Aliases a stable Realm member that
    is rebuilt each regrid; set at define(), never owned here.
  */
  const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>* m_validCells = nullptr;

  /**
    @brief Whether define() has been called.
  */
  bool m_isDefined;

  /**
    @brief Whether the valid leaves are currently cell-sorted.
  */
  bool m_isOrganizedByCell;

  /**
    @brief Valid particles: one ParticleSoA leaf per box, per level.
  */
  AMRParticlesSoA<P, Traits> m_particles;

  /**
    @brief Regrid cache: the pre-regrid valid holders, kept alive over the old layout by preRegrid().
  */
  AMRParticlesSoA<P, Traits> m_cacheParticles;

  /**
    @brief The old (pre-regrid) grids the cache holders are defined on (for iterating the cache).
  */
  Vector<DisjointBoxLayout> m_cacheGrids;

  /**
    @brief Halo/mask particles per level.
  */
  AMRParticlesSoA<P, Traits> m_maskParticles;

  /**
    @brief Grown-grid buffer particles per level (on m_grownGrids; filled during halo deposition).
  */
  AMRParticlesSoA<P, Traits> m_bufferParticles;

  /**
    @brief Per-level grown grids (boxes grown by the refinement factor on finer levels) for the buffer.
  */
  Vector<BoxLayout> m_grownGrids;

  /**
    @brief Cell index containing a position on a given level (floor((x - probLo)/dx)).
    @param[in] a_lvl Level index.
    @param[in] a_pos Position.
    @return The containing cell IntVect.
  */
  inline IntVect
  cellIndex(const int a_lvl, const RealVect& a_pos) const
  {
    IntVect iv;
    for (int dir = 0; dir < SpaceDim; dir++) {
      iv[dir] = static_cast<int>(std::floor((a_pos[dir] - m_probLo[dir]) / m_dx[a_lvl][dir]));
    }
    return iv;
  }

  /**
    @brief Resolved remap destination of a particle: owning (level, grid index, rank).
  */
  struct Destination
  {
    /**
      @brief Owning level (finest tile that contains the particle).
    */
    int level;

    /**
      @brief Owning grid index within that level's LevelTiles.
    */
    unsigned int gridIndex;

    /**
      @brief Owning rank.
    */
    int rank;

    /**
      @brief False if the particle is off-domain (owned by no tile).
    */
    bool found;
  };

  /**
    @brief Map a position to its owning (level, grid index, rank) via the finest containing tile.
    @param[in] a_pos Particle position.
    @return The resolved destination; found==false if no tile on any level owns the cell.
  */
  Destination
  findDestination(const RealVect& a_pos) const;
};

#include <CD_NamespaceFooter.H>

#include <CD_ParticleContainerImplem.H>

#endif
