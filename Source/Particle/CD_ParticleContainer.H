/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ParticleContainer.H
  @brief  Dev staging: AMR-hierarchy owner of per-patch ParticleSoA leaves (storage layer).
  @author Robert Marskar
  @details
  ParticleContainer is the SoA replacement for Source/Particle/CD_ParticleContainer.H. It owns
  particles across the AMR hierarchy: one @c LayoutData<ParticleSoA<P>> per level, indexed by
  @c DataIndex over that level's @c DisjointBoxLayout (mirroring Chombo's ParticleData<List<P>>).

  Built in stages. STAGE-1 (storage / ownership / accessors) + STAGE-2a (remap) are present here:
  per-level holders, accessors, valid/outcast counts, and `remap()` -- the pool-model
  redistribution of valid particles to the owning patch/level/rank using per-level LevelTiles
  (built in define from the grids + blocking factor) and an MPI Alltoallv scatter via the leaf's
  byte linearization. Still deferred to stage-2b..d: regrid (preRegrid cache + redistribute),
  halo/mask + buffer FILLING, and container-level cell-sort. See ARCHITECTURE.md.

  define() is lighter than production: it needs grids/domains/dx/refinement-ratios + probLo +
  blocking factor + finest level (no valid mask -- remap's finest-tile-wins lookup does not need
  it), and no EB geometry at all (the container never references EBLevelGrid/EBISBox).

  Named with the "SoA" suffix so it coexists with the production ParticleContainer (linked via
  libSource) until design-freeze.

  @note Move/ownership: the per-level holders are shared @c RefCountedPtr, so the container is
        non-copyable (a copy would alias the holders). Construct empty, then call define().
  @note The mask and buffer holders are allocated (so consumers can read them) but left EMPTY in
        stage-1; filling them, and giving buffer its proper grown-grid layout, is stage-2.
  @tparam P      User payload struct (extra columns only); defaults to NoPayload.
  @tparam Traits Column descriptor for the payload; defaults to ParticleTraits<P>.
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
  @brief AMR-hierarchy container of per-patch ParticleSoA leaves (stage-1 storage layer).
  @details See the file-level documentation for the data model and the stage-1/stage-2 split.
  @tparam P      User payload struct (extra columns only); defaults to NoPayload.
  @tparam Traits Column descriptor for the payload; defaults to ParticleTraits<P>.
*/
template <typename P = NoPayload, typename Traits = ParticleTraits<P>>
class ParticleContainer
{
public:
  /** @brief The per-patch leaf type owned at every (level, box). */
  using Leaf = ParticleSoA<P, Traits>;

  /** @brief The per-level holder type (one Leaf per box on the level's DisjointBoxLayout). */
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
      m_blockingFactor(0),
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
    @details Allocates one empty Leaf per box for the valid, mask, and buffer holders on every
    level in [0, a_finestLevel], and builds a LevelTiles per level (from the grids + blocking
    factor) used by remap() to map a particle to its owning patch/rank. Mask/buffer are
    placeholders (empty) until stage-2c.
    @param[in] a_grids          Per-level DisjointBoxLayouts.
    @param[in] a_domains        Per-level problem domains.
    @param[in] a_dx             Per-level (isotropic) grid spacing.
    @param[in] a_refRat         Per-level refinement ratio (entry l is the ratio between level l and l+1).
    @param[in] a_probLo         Lower-left corner of the physical domain.
    @param[in] a_blockingFactor Grid blocking factor (tile size); the grids must be unions of
                                blockingFactor-sized boxes for the tile ownership to be exact.
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
         const int                                              a_blockingFactor,
         const int                                              a_finestLevel,
         const std::string&                                     a_realm,
         const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>* a_validCells)
  {
    CH_TIME("ParticleContainer::define");

    CH_assert(a_grids.size() > a_finestLevel);
    CH_assert(a_domains.size() > a_finestLevel);
    CH_assert(a_dx.size() > a_finestLevel);
    CH_assert(a_refRat.size() > a_finestLevel);

    m_grids           = a_grids;
    m_domains         = a_domains;
    m_refRat          = a_refRat;
    m_probLo          = a_probLo;
    m_blockingFactor  = a_blockingFactor;
    m_finestLevel     = a_finestLevel;
    m_realm           = a_realm;
    m_validCells      = a_validCells;
    m_numOutcastLocal = 0;

    m_dx.resize(a_finestLevel + 1);
    for (int lvl = 0; lvl <= a_finestLevel; lvl++) {
      m_dx[lvl] = a_dx[lvl] * RealVect::Unit;
    }

    m_particles.resize(a_finestLevel + 1);
    m_maskParticles.resize(a_finestLevel + 1);
    m_levelTiles.resize(a_finestLevel + 1);
    for (int lvl = 0; lvl <= a_finestLevel; lvl++) {
      m_particles[lvl]     = RefCountedPtr<LevelParticles>(new LevelParticles(a_grids[lvl]));
      m_maskParticles[lvl] = RefCountedPtr<LevelParticles>(new LevelParticles(a_grids[lvl]));
      m_levelTiles[lvl]    = RefCountedPtr<LevelTiles>(new LevelTiles(a_grids[lvl], a_blockingFactor));
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
    @brief The halo/mask particles on all levels (empty in stage-1).
    @return Mutable per-level holder vector.
  */
  AMRParticlesSoA<P, Traits>&
  getMaskParticles()
  {
    return m_maskParticles;
  }

  /**
    @brief The halo/mask particles on all levels (const, empty in stage-1).
    @return Const per-level holder vector.
  */
  const AMRParticlesSoA<P, Traits>&
  getMaskParticles() const
  {
    return m_maskParticles;
  }

  /**
    @brief The grown-grid buffer particles on all levels (empty placeholder in stage-1).
    @return Mutable per-level holder vector.
  */
  AMRParticlesSoA<P, Traits>&
  getBufferParticles()
  {
    return m_bufferParticles;
  }

  /**
    @brief The grown-grid buffer particles on all levels (const, empty placeholder in stage-1).
    @return Const per-level holder vector.
  */
  const AMRParticlesSoA<P, Traits>&
  getBufferParticles() const
  {
    return m_bufferParticles;
  }

  /**
    @brief Valid particles on one level (index by DataIndex to reach a patch's Leaf).
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
    @brief Whether the container is cell-sorted (always false in stage-1).
    @return False (container-level cell sorting is stage-2).
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
    @name Population (stage-1: local only, no cross-rank routing)
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
      LevelParticles& level = *m_particles[lvl];
      for (DataIterator dit = m_grids[lvl].dataIterator(); dit.ok(); ++dit) {
        level[dit()].clear();
      }
    }
  }

  /**
    @brief Append one particle into the locally-owned patch that contains its cell.
    @details Finds the box on this rank that contains the particle's cell on level a_lvl and
    appends to that patch's leaf. If no local box owns the cell (it belongs to another rank or
    lies outside the domain), the call is a no-op -- there is no cross-rank routing in stage-1
    (that is remap, stage-2).
    @param[in] a_lvl     Level index.
    @param[in] a_pos     Particle position.
    @param[in] a_weight  Particle weight.
    @param[in] a_payload Payload fields.
  */
  void
  addParticlesLocal(const int a_lvl, const RealVect& a_pos, const double a_weight, const P& a_payload)
  {
    CH_assert(m_isDefined);
    CH_assert(a_lvl >= 0 && a_lvl <= m_finestLevel);

    const IntVect            iv    = this->cellIndex(a_lvl, a_pos);
    LevelParticles&          level = *m_particles[a_lvl];
    const DisjointBoxLayout& dbl   = m_grids[a_lvl];
    m_isOrganizedByCell            = false; // adding particles invalidates any cell sort
    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
      if (dbl[dit()].contains(iv)) {
        level[dit()].append(a_pos, a_weight, a_payload);
        return;
      }
    }
  }

  /**
    @brief Add a free-standing buffer of particles to the container, routing each to its owner.
    @details Destructive SoA analogue of ParticleContainer::addParticlesDestructive(List<P>): the buffer
    (e.g. filled by ParticleManagement::drawRandomParticles(ParticleSoA)) is catenated onto a local leaf
    -- O(1) when that leaf is empty -- and then remap() routes every particle (including across ranks and
    levels) to the patch/level/rank that owns its position. a_particles is left empty. Assumes this rank
    owns at least one box (true whenever level 0 covers the domain).
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

    this->remap(); // route every staged particle to its owning patch/level/rank
  }

  ///@}

  /**
    @name Counting
  */
  ///@{

  /**
    @brief Number of valid particles owned by this rank.
    @return Local valid particle count summed over all levels and patches.
  */
  unsigned long long
  getNumberOfValidParticlesLocal() const
  {
    CH_TIME("ParticleContainer::getNumberOfValidParticlesLocal");

    CH_assert(m_isDefined);

    unsigned long long n = 0;
    for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
      const LevelParticles& level = *m_particles[lvl];
      for (DataIterator dit = m_grids[lvl].dataIterator(); dit.ok(); ++dit) {
        n += static_cast<unsigned long long>(level[dit()].size());
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

    const unsigned long long local = this->getNumberOfValidParticlesLocal();
#ifdef CH_MPI
    unsigned long long global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, Chombo_MPI::comm);
    return global;
#else
    return local;
#endif
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
#ifdef CH_MPI
    unsigned long long global = 0;
    MPI_Allreduce(&m_numOutcastLocal, &global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, Chombo_MPI::comm);
    return global;
#else
    return m_numOutcastLocal;
#endif
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
    @param[in] a_blockingFactor Grid blocking factor (tile size) for the new grids.
    @param[in] a_newFinestLevel New finest AMR level index.
  */
  void
  regrid(const Vector<DisjointBoxLayout>& a_grids,
         const Vector<ProblemDomain>&     a_domains,
         const Vector<Real>&              a_dx,
         const Vector<int>&               a_refRat,
         const int                        a_blockingFactor,
         const int                        a_newFinestLevel);

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
      LevelParticles& level = *m_maskParticles[lvl];
      for (DataIterator dit = m_grids[lvl].dataIterator(); dit.ok(); ++dit) {
        level[dit()].clear();
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
      LevelParticles& level = *m_bufferParticles[lvl];
      for (DataIterator dit = m_grownGrids[lvl].dataIterator(); dit.ok(); ++dit) {
        level[dit()].clear();
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
      for (DataIterator dit = m_grids[lvl].dataIterator(); dit.ok(); ++dit) {
        dst[dit()].catenate(src[dit()]); // bulk move (id/rank preserved); empty dst -> O(1) swap
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
      for (DataIterator dit = m_grids[lvl].dataIterator(); dit.ok(); ++dit) {
        level[dit()].sortByCell(m_grids[lvl][dit()], m_dx[lvl], m_probLo);
      }
    }
    m_isOrganizedByCell = true;
  }

  /**
    @brief Return to a by-patch view (the leaves remain a single contiguous array per patch).
    @details Clears isOrganizedByCell(); the CSR cell offsets are no longer to be relied upon. There
    is no data movement -- unlike the production BinFab path, the SoA leaf already holds all particles
    of a patch contiguously.
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
    @brief Assign a mover pool into the (assumed empty) valid holders: local appends + MPI scatter.
    @param[in,out] a_pool Per-rank mover pool produced by gatherToPool (consumed).
  */
  void
  distributeFromPool(MoverPool& a_pool);

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
  int m_blockingFactor;

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
    @brief Whether the container is cell-sorted (always false in stage-1).
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
    @brief Halo/mask particles per level (allocated, empty in stage-1).
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
