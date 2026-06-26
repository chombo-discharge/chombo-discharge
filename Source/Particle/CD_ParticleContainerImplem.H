/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ParticleContainerImplem.H
  @brief  Implementation of CD_ParticleContainer.H
  @author Robert Marskar
*/

#ifndef CD_PARTICLECONTAINERIMPLEM_H
#define CD_PARTICLECONTAINERIMPLEM_H

// Std includes
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <utility>
#include <vector>

// Chombo includes
#include <CH_Timer.H>
#include <CH_assert.H>

// Our includes
#include <CD_ParticleContainer.H>
#include <CD_NamespaceHeader.H>

template <typename P, typename Traits>
inline typename ParticleContainer<P, Traits>::Destination
ParticleContainer<P, Traits>::findDestination(const RealVect& a_pos) const
{
  // Finest level whose tile contains the particle wins.
  for (int lvl = m_finestLevel; lvl >= 0; lvl--) {
    IntVect tile;

    for (int dir = 0; dir < SpaceDim; dir++) {
      tile[dir] = static_cast<int>(std::floor((a_pos[dir] - m_probLo[dir]) / (m_blockingFactor * m_dx[lvl][dir])));
    }

    const LevelTiles& tiles = *m_levelTiles[lvl];

    const auto& myTiles = tiles.getMyTiles();
    const auto  mit     = myTiles.find(tile);

    if (mit != myTiles.end()) {
      return Destination{lvl, mit->second, procID(), true};
    }

#ifdef CH_MPI
    const auto& otherTiles = tiles.getOtherTiles();
    const auto  oit        = otherTiles.find(tile);

    if (oit != otherTiles.end()) {
      return Destination{lvl, oit->second.first, static_cast<int>(oit->second.second), true};
    }
#endif
  }

  return Destination{-1, 0u, -1, false};
}

template <typename P, typename Traits>
inline typename ParticleContainer<P, Traits>::MoverPool
ParticleContainer<P, Traits>::gatherToPool(const AMRParticlesSoA<P, Traits>& a_source,
                                           const Vector<DisjointBoxLayout>&  a_sourceGrids)
{
  CH_TIME("ParticleContainer::gatherToPool");

  MoverPool pool(numProc());

  // Bin every particle of one leaf into a destination pool, returning the number of off-domain
  // (outcast) particles. Shared by the serial and OpenMP paths so the binning logic lives once.
  // findDestination() is const and reads only the (immutable-during-gather) tile maps, so it is
  // safe to call concurrently; each thread bins into its own thread-local pool.
  auto binLeaf = [this](MoverPool& a_target, const Leaf& a_leaf) -> unsigned long long {
    unsigned long long outcast = 0;

    for (std::size_t i = 0; i < a_leaf.size(); i++) {
      const RealVect    pos = a_leaf.position(i);
      const Destination dst = this->findDestination(pos);

      if (!dst.found) {
        outcast++;
        continue;
      }

      Leaf& bin = a_target[dst.rank][PoolKey(dst.level, dst.gridIndex)];

      bin.append(pos, a_leaf.weight(i), a_leaf.gather(i));

      const std::size_t k = bin.size() - 1;

      bin.particleID(k) = a_leaf.particleID(i); // preserve id
      bin.rankID(k)     = dst.rank;             // new owner
    }

    return outcast;
  };

  unsigned long long numOutcast = 0;

  for (int lvl = 0; lvl < a_source.size(); lvl++) {
    const LevelParticles& level = *a_source[lvl];
    const DataIterator&   dit   = a_sourceGrids[lvl].dataIterator();
    const int             nbox  = dit.size();

#ifdef _OPENMP
#pragma omp parallel reduction(+ : numOutcast)
    {
      MoverPool threadPool(numProc()); // thread-local; no locking while binning

#pragma omp for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        numOutcast += binLeaf(threadPool, level[dit[mybox]]);
      }

      // Merge the thread-local pool into the shared pool. First writer to a (rank,key) slot gets the
      // O(1) arena swap inside catenate(); later writers append.
#pragma omp critical(ParticleContainer_gatherMerge)
      {
        for (int r = 0; r < numProc(); r++) {
          for (auto& kv : threadPool[r]) {
            pool[r][kv.first].catenate(kv.second);
          }
        }
      }
    }
#else
    for (int mybox = 0; mybox < nbox; mybox++) {
      numOutcast += binLeaf(pool, level[dit[mybox]]);
    }
#endif
  }

  m_numOutcastLocal = numOutcast;

  return pool;
}

template <typename P, typename Traits>
inline typename ParticleContainer<P, Traits>::MoverPool
ParticleContainer<P, Traits>::gatherMoversToPool()
{
  CH_TIME("ParticleContainer::gatherMoversToPool");

  CH_assert(m_validCells != nullptr);

  MoverPool          pool(numProc());
  unsigned long long numOutcast = 0;

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    LevelParticles&                 level     = *m_particles[lvl];
    const LevelData<BaseFab<bool>>* maskLevel = (*m_validCells)[lvl].isNull() ? nullptr : &(*(*m_validCells)[lvl]);
    const DataIterator&             dit       = m_grids[lvl].dataIterator();
    const int                       nbox      = dit.size();

    // Partition one leaf in place: keep stayers (still inside their box AND in a cell flagged valid =
    // not covered by a finer level), pool the movers and swap-pop them out. a_mask == nullptr forces
    // every particle to be treated as a mover (no usable mask for this leaf). Returns outcast movers.
    // findDestination()/cellIndex() are const reads, so this is safe to call from many threads as long
    // as each leaf is touched by a single thread (the omp-for partitions boxes disjointly).
    auto partitionLeaf = [this, lvl](MoverPool&                        a_target,
                                     [[maybe_unused]] const DataIndex& a_din,
                                     Leaf&                             a_leaf,
                                     const BaseFab<bool>*              a_mask) -> unsigned long long {
      unsigned long long outcast = 0;

      const Box maskBox = (a_mask != nullptr) ? a_mask->box() : Box();

      std::size_t i = 0;

      while (i < a_leaf.size()) {
        const RealVect pos = a_leaf.position(i);
        if (a_mask != nullptr) {
          const IntVect iv = this->cellIndex(lvl, pos);

          if (maskBox.contains(iv) && (*a_mask)(iv, 0)) {
            // Stayer: in its own box and not under a finer grid, so it cannot have changed owner.
            // Debug self-check of that invariant -- the destination the gather-all path would compute
            // must be exactly this (level, grid, rank). If this ever fires, the skip is unsound.
#ifndef NDEBUG
            const Destination chk = this->findDestination(pos);
            CH_assert(chk.found && chk.level == lvl && chk.rank == procID() &&
                      m_levelTiles[lvl]->getMyGrids().at(chk.gridIndex) == a_din);
#endif
            i++;

            continue;
          }
        }

        const Destination dst = this->findDestination(pos); // mover: route it

        if (dst.found) {
          Leaf& bin = a_target[dst.rank][PoolKey(dst.level, dst.gridIndex)];

          bin.append(pos, a_leaf.weight(i), a_leaf.gather(i));

          const std::size_t k = bin.size() - 1;

          bin.particleID(k) = a_leaf.particleID(i); // preserve id
          bin.rankID(k)     = dst.rank;             // new owner
        }
        else {
          outcast++;
        }

        a_leaf.remove(i); // swap-pop; do NOT advance i (a new particle now occupies slot i)
      }

      return outcast;
    };

#ifdef _OPENMP
#pragma omp parallel reduction(+ : numOutcast)
    {
      MoverPool threadPool(numProc());

#pragma omp for schedule(runtime)
      for (int mybox = 0; mybox < nbox; mybox++) {
        const DataIndex&     din  = dit[mybox];
        const BaseFab<bool>* mask = (maskLevel != nullptr) ? &(*maskLevel)[din] : nullptr;

        if (mask != nullptr && !mask->isUsable()) {
          mask = nullptr;
        }

        numOutcast += partitionLeaf(threadPool, din, level[din], mask);
      }

#pragma omp critical(ParticleContainer_gatherMerge)
      {
        for (int r = 0; r < numProc(); r++) {
          for (auto& kv : threadPool[r]) {
            pool[r][kv.first].catenate(kv.second);
          }
        }
      }
    }
#else
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex&     din  = dit[mybox];
      const BaseFab<bool>* mask = (maskLevel != nullptr) ? &(*maskLevel)[din] : nullptr;

      if (mask != nullptr && !mask->isUsable()) {
        mask = nullptr;
      }

      numOutcast += partitionLeaf(pool, din, level[din], mask);
    }
#endif
  }

  m_numOutcastLocal = numOutcast;

  return pool;
}

template <typename P, typename Traits>
inline void
ParticleContainer<P, Traits>::distributeFromPool(MoverPool& a_pool)
{
  CH_TIME("ParticleContainer::distributeFromPool");

  const int myRank   = procID();
  const int numRanks = numProc();

  // ---- same-rank movers: catenate each pooled leaf into its destination leaf ----
  // On the gather-all path the destination leaves are empty (remap clears them; regrid allocates them
  // fresh) so catenate() takes the O(1) arena-swap fast path. On the mask fast path the destinations
  // still hold their stayers, so catenate() appends (O(movers)); a mover never returns to its own
  // origin leaf, so there is no self-append.
  for (auto& kv : a_pool[myRank]) {
    const int          lvl  = kv.first.first;
    const unsigned int gidx = kv.first.second;
    const DataIndex    din  = m_levelTiles[lvl]->getMyGrids().at(gidx);

    (*m_particles[lvl])[din].catenate(kv.second);
  }

#ifdef CH_MPI
  // ---- cross-rank movers: one Alltoallv of linearized particles ----
  constexpr std::size_t pp        = Leaf::bytesPerParticle();
  constexpr std::size_t headerLen = 2 * sizeof(std::uint32_t) + sizeof(std::uint64_t);

  std::vector<std::vector<char>> sendBuf(numRanks);
  // Each destination rank's buffer and pool map are disjoint, and linearizeParticle is a const read,
  // so packing different ranks in parallel is race-free.
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime)
#endif
  for (int r = 0; r < numRanks; r++) {
    if (r == myRank) {
      continue;
    }

    std::vector<char>& buf = sendBuf[r];

    for (auto& kv : a_pool[r]) {
      Leaf& src = kv.second;

      const std::uint32_t lvl   = static_cast<std::uint32_t>(kv.first.first);
      const std::uint32_t gidx  = static_cast<std::uint32_t>(kv.first.second);
      const std::uint64_t count = static_cast<std::uint64_t>(src.size());
      const std::size_t   off   = buf.size();

      buf.resize(off + headerLen + count * pp);

      char* p = buf.data() + off;

      std::memcpy(p, &lvl, sizeof(std::uint32_t));
      p += sizeof(std::uint32_t);

      std::memcpy(p, &gidx, sizeof(std::uint32_t));
      p += sizeof(std::uint32_t);

      std::memcpy(p, &count, sizeof(std::uint64_t));
      p += sizeof(std::uint64_t);

      for (std::size_t i = 0; i < src.size(); i++) {
        src.linearizeParticle(p, i);
        p += pp;
      }
    }
  }

  // Exchange byte counts, then the byte buffers.
  std::vector<int> sendCounts(numRanks, 0);
  std::vector<int> recvCounts(numRanks, 0);

  for (int r = 0; r < numRanks; r++) {
    CH_assert(sendBuf[r].size() <= static_cast<std::size_t>(std::numeric_limits<int>::max()));
    sendCounts[r] = static_cast<int>(sendBuf[r].size());
  }

  MPI_Alltoall(sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT, Chombo_MPI::comm);

  std::vector<int> sdispl(numRanks, 0);
  std::vector<int> rdispl(numRanks, 0);

  std::int64_t stot = 0;
  std::int64_t rtot = 0;

  for (int r = 0; r < numRanks; r++) {
    sdispl[r] = static_cast<int>(stot);
    stot += sendCounts[r];

    rdispl[r] = static_cast<int>(rtot);
    rtot += recvCounts[r];
  }

  // The accumulators are int64 so the cross-rank sum cannot overflow; MPI displacement arrays are
  // int, so assert the totals fit before allocating the flat buffers and calling MPI_Alltoallv.
  CH_assert(stot <= static_cast<std::int64_t>(std::numeric_limits<int>::max()));
  CH_assert(rtot <= static_cast<std::int64_t>(std::numeric_limits<int>::max()));

  std::vector<char> sflat(stot);
  std::vector<char> rflat(rtot);

  for (int r = 0; r < numRanks; r++) {
    if (sendCounts[r] > 0) {
      std::memcpy(sflat.data() + sdispl[r], sendBuf[r].data(), static_cast<std::size_t>(sendCounts[r]));
    }
  }

  MPI_Alltoallv(sflat.data(),
                sendCounts.data(),
                sdispl.data(),
                MPI_BYTE,
                rflat.data(),
                recvCounts.data(),
                rdispl.data(),
                MPI_BYTE,
                Chombo_MPI::comm);

  // ---- unpack: group received chunks by destination leaf, then append in parallel ----
  // Pass 1 (serial, header-only walk): record each chunk's payload pointer + count, grouped by its
  // destination (level, gridIndex). Pass 2 appends in parallel with one destination leaf per thread,
  // so no two threads ever touch the same leaf (delinearizeAndAppend mutates the leaf arena).
  struct RecvChunk
  {
    const char*   payload; // first particle of the chunk within rflat
    std::uint64_t count;
  };
  struct DestGroup
  {
    int                    level;
    unsigned int           gridIndex;
    std::uint64_t          total;
    std::vector<RecvChunk> chunks;
  };

  std::vector<DestGroup>         groups;
  std::map<PoolKey, std::size_t> keyToGroup;
  {
    const char* p   = rflat.data();
    const char* end = p + rtot;

    while (p < end) {
      std::uint32_t lvl;
      std::uint32_t gidx;
      std::uint64_t count;

      std::memcpy(&lvl, p, sizeof(std::uint32_t));
      p += sizeof(std::uint32_t);

      std::memcpy(&gidx, p, sizeof(std::uint32_t));
      p += sizeof(std::uint32_t);

      std::memcpy(&count, p, sizeof(std::uint64_t));
      p += sizeof(std::uint64_t);

      const PoolKey key(static_cast<int>(lvl), gidx);
      const auto    it = keyToGroup.find(key);

      if (it == keyToGroup.end()) {
        keyToGroup.emplace(key, groups.size());
        groups.push_back(DestGroup{static_cast<int>(lvl), gidx, count, {RecvChunk{p, count}}});
      }
      else {
        DestGroup& g = groups[it->second];
        g.total += count;
        g.chunks.push_back(RecvChunk{p, count});
      }

      p += count * pp;
    }
  }

  const int numGroups = static_cast<int>(groups.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(runtime)
#endif
  for (int d = 0; d < numGroups; d++) {
    const DestGroup& g   = groups[d];
    const DataIndex  din = m_levelTiles[g.level]->getMyGrids().at(g.gridIndex);

    Leaf& tgt = (*m_particles[g.level])[din];

    tgt.reserve(tgt.size() + g.total); // one allocation; avoids geometric regrow during the appends
    for (const RecvChunk& chunk : g.chunks) {
      const char* q = chunk.payload;

      for (std::uint64_t c = 0; c < chunk.count; c++) {
        tgt.delinearizeAndAppend(q);
        q += pp;
        tgt.rankID(tgt.size() - 1) = myRank; // owner is the receiving rank
      }
    }
  }
#endif
}

template <typename P, typename Traits>
inline void
ParticleContainer<P, Traits>::remap()
{
  CH_TIME("ParticleContainer::remap");

  CH_assert(m_isDefined);

  if (m_validCells != nullptr && static_cast<int>(m_validCells->size()) >= m_finestLevel + 1) {
    // Mask fast path: a particle still in its box and not under a finer grid cannot change owner, so we
    // leave those "stayers" in place and pool only the movers. Destinations therefore keep their stayers
    // -- do NOT clearParticles() here; distributeFromPool() appends the movers onto the existing leaves.
    MoverPool pool = this->gatherMoversToPool();

    this->distributeFromPool(pool);
  }
  else {
    // Fallback: re-pool every valid particle (used when no mask is available, e.g. addParticlesDestructive
    // staging, or before the realm's valid-cell mask has been built).
    MoverPool pool = this->gatherToPool(m_particles, m_grids);

    this->clearParticles();
    this->distributeFromPool(pool);
  }

  m_isOrganizedByCell = false;
}

template <typename P, typename Traits>
inline void
ParticleContainer<P, Traits>::preRegrid()
{
  CH_TIME("ParticleContainer::preRegrid");

  CH_assert(m_isDefined);

  // Share the current holders (RefCountedPtr) so regrid() can rebuild m_particles while the old
  // particles stay alive in the cache. Remember the old layout for iterating the cache.
  m_cacheParticles = m_particles;
  m_cacheGrids     = m_grids;
}

template <typename P, typename Traits>
inline void
ParticleContainer<P, Traits>::regrid(const Vector<DisjointBoxLayout>& a_grids,
                                     const Vector<ProblemDomain>&     a_domains,
                                     const Vector<Real>&              a_dx,
                                     const Vector<int>&               a_refRat,
                                     const int                        a_blockingFactor,
                                     const int                        a_newFinestLevel)
{
  CH_TIME("ParticleContainer::regrid");

  CH_assert(m_isDefined);

  // ---- adopt the new layout ----
  m_grids          = a_grids;
  m_domains        = a_domains;
  m_refRat         = a_refRat;
  m_blockingFactor = a_blockingFactor;
  m_finestLevel    = a_newFinestLevel;

  m_dx.resize(a_newFinestLevel + 1);
  for (int lvl = 0; lvl <= a_newFinestLevel; lvl++) {
    m_dx[lvl] = a_dx[lvl] * RealVect::Unit;
  }

  // ---- re-allocate the (empty) holders and tile maps over the new grids ----
  m_particles.resize(a_newFinestLevel + 1);
  m_maskParticles.resize(a_newFinestLevel + 1);
  m_levelTiles.resize(a_newFinestLevel + 1);

  for (int lvl = 0; lvl <= a_newFinestLevel; lvl++) {
    m_particles[lvl]     = RefCountedPtr<LevelParticles>(new LevelParticles(a_grids[lvl]));
    m_maskParticles[lvl] = RefCountedPtr<LevelParticles>(new LevelParticles(a_grids[lvl]));
    m_levelTiles[lvl]    = RefCountedPtr<LevelTiles>(new LevelTiles(a_grids[lvl], a_blockingFactor));
  }

  this->setupGrownGrids();

  m_bufferParticles.resize(a_newFinestLevel + 1);

  for (int lvl = 0; lvl <= a_newFinestLevel; lvl++) {
    m_bufferParticles[lvl] = RefCountedPtr<LevelParticles>(new LevelParticles(m_grownGrids[lvl]));
  }

  // ---- redistribute the cached particles onto the new layout ----
  MoverPool pool = this->gatherToPool(m_cacheParticles, m_cacheGrids);
  this->distributeFromPool(pool);

  // ---- release the cache ----
  m_cacheParticles.resize(0);
  m_cacheGrids.resize(0);

  m_isOrganizedByCell = false;
}

template <typename P, typename Traits>
inline void
ParticleContainer<P, Traits>::setupGrownGrids()
{
  CH_TIME("ParticleContainer::setupGrownGrids");

  m_grownGrids.resize(m_finestLevel + 1);
  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    const DisjointBoxLayout& dbl = m_grids[lvl];

    Vector<Box> boxes = dbl.boxArray();

    if (lvl > 0) {
      const Box domainBox = m_domains[lvl].domainBox();
      for (int i = 0; i < boxes.size(); i++) {
        boxes[i].grow(m_refRat[lvl - 1]);
        boxes[i] &= domainBox;
      }
    }

    m_grownGrids[lvl] = BoxLayout(boxes, dbl.procIDs());
  }
}

template <typename P, typename Traits>
inline void
ParticleContainer<P, Traits>::copyMaskParticles(const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_mask)
{
  CH_TIME("ParticleContainer::copyMaskParticles(amr)");

  CH_assert(m_isDefined);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    if (!a_mask[lvl].isNull()) {
      this->copyMaskParticles(lvl, *a_mask[lvl]);
    }
  }
}

template <typename P, typename Traits>
inline void
ParticleContainer<P, Traits>::copyMaskParticles(const int a_level, const LevelData<BaseFab<bool>>& a_mask)
{
  CH_TIME("ParticleContainer::copyMaskParticles(level)");

  CH_assert(m_isDefined);

  LevelParticles&       maskLevel = *m_maskParticles[a_level];
  const LevelParticles& srcLevel  = *m_particles[a_level];

  const DisjointBoxLayout& dbl = m_grids[a_level];
  const DataIterator&      dit = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    maskLevel[din].clear();

    const BaseFab<bool>& mask = a_mask[din];
    if (!mask.isUsable()) {
      continue;
    }
    const Box maskBox = mask.box();

    const Leaf& src = srcLevel[din];
    Leaf&       dst = maskLevel[din];

    for (std::size_t i = 0; i < src.size(); i++) {
      const IntVect iv = this->cellIndex(a_level, src.position(i));

      if (maskBox.contains(iv) && mask(iv, 0)) {
        dst.append(src.position(i), src.weight(i), src.gather(i));

        const std::size_t k = dst.size() - 1;

        dst.particleID(k) = src.particleID(i);
        dst.rankID(k)     = src.rankID(i);
      }
    }
  }
}

template <typename P, typename Traits>
inline void
ParticleContainer<P, Traits>::transferMaskParticles(const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_mask)
{
  CH_TIME("ParticleContainer::transferMaskParticles(amr)");

  CH_assert(m_isDefined);

  for (int lvl = 0; lvl <= m_finestLevel; lvl++) {
    if (!a_mask[lvl].isNull()) {
      this->transferMaskParticles(lvl, *a_mask[lvl]);
    }
  }
}

template <typename P, typename Traits>
inline void
ParticleContainer<P, Traits>::transferMaskParticles(const int a_level, const LevelData<BaseFab<bool>>& a_mask)
{
  CH_TIME("ParticleContainer::transferMaskParticles(level)");

  CH_assert(m_isDefined);

  LevelParticles& maskLevel = *m_maskParticles[a_level];
  LevelParticles& srcLevel  = *m_particles[a_level];

  const DisjointBoxLayout& dbl = m_grids[a_level];
  const DataIterator&      dit = dbl.dataIterator();

  const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
  for (int mybox = 0; mybox < nbox; mybox++) {
    const DataIndex& din = dit[mybox];

    const BaseFab<bool>& mask = a_mask[din];
    if (!mask.isUsable()) {
      continue;
    }
    const Box maskBox = mask.box();

    Leaf& src = srcLevel[din];
    Leaf& dst = maskLevel[din];

    std::size_t i = 0;
    while (i < src.size()) {
      const IntVect iv = this->cellIndex(a_level, src.position(i));

      if (maskBox.contains(iv) && mask(iv, 0)) {
        dst.append(src.position(i), src.weight(i), src.gather(i));

        const std::size_t k = dst.size() - 1;

        dst.particleID(k) = src.particleID(i);
        dst.rankID(k)     = src.rankID(i);

        src.remove(i); // swap-and-pop; the swapped-in particle now occupies index i
      }
      else {
        ++i;
      }
    }
  }
}

#include <CD_NamespaceFooter.H>

#endif
