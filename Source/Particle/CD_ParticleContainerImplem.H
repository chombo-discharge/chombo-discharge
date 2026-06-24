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

  m_numOutcastLocal = 0;

  MoverPool pool(numProc());

  for (int lvl = 0; lvl < a_source.size(); lvl++) {
    const LevelParticles& level = *a_source[lvl];
    for (DataIterator dit = a_sourceGrids[lvl].dataIterator(); dit.ok(); ++dit) {
      const Leaf& leaf = level[dit()];
      for (std::size_t i = 0; i < leaf.size(); i++) {
        const RealVect    pos = leaf.position(i);
        const Destination dst = this->findDestination(pos);
        if (!dst.found) {
          m_numOutcastLocal++;
          continue;
        }
        Leaf& bin = pool[dst.rank][PoolKey(dst.level, dst.gridIndex)];
        bin.append(pos, leaf.weight(i), leaf.gather(i));
        const std::size_t k = bin.size() - 1;
        bin.particleID(k)   = leaf.particleID(i); // preserve id
        bin.rankID(k)       = dst.rank;           // new owner
      }
    }
  }

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
  // The destination valid leaves are empty here (remap clears them; regrid allocates them fresh), so
  // catenate() takes the O(1) arena-swap fast path -- a zero-copy steal of the pooled leaf.
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
  for (int r = 0; r < numRanks; r++) {
    if (r == myRank) {
      continue;
    }
    std::vector<char>& buf = sendBuf[r];
    for (auto& kv : a_pool[r]) {
      Leaf&               src   = kv.second;
      const std::uint32_t lvl   = static_cast<std::uint32_t>(kv.first.first);
      const std::uint32_t gidx  = static_cast<std::uint32_t>(kv.first.second);
      const std::uint64_t count = static_cast<std::uint64_t>(src.size());

      const std::size_t off = buf.size();
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
    sendCounts[r] = static_cast<int>(sendBuf[r].size());
  }
  MPI_Alltoall(sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT, Chombo_MPI::comm);

  std::vector<int> sdispl(numRanks, 0);
  std::vector<int> rdispl(numRanks, 0);
  int              stot = 0;
  int              rtot = 0;
  for (int r = 0; r < numRanks; r++) {
    sdispl[r] = stot;
    stot += sendCounts[r];
    rdispl[r] = rtot;
    rtot += recvCounts[r];
  }

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

  // Unpack received movers into their destination leaves.
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

    const DataIndex din = m_levelTiles[lvl]->getMyGrids().at(gidx);
    Leaf&           tgt = (*m_particles[lvl])[din];
    for (std::uint64_t c = 0; c < count; c++) {
      tgt.delinearizeAndAppend(p);
      p += pp;
      tgt.rankID(tgt.size() - 1) = myRank; // owner is the receiving rank
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

  MoverPool pool = this->gatherToPool(m_particles, m_grids); // pool a copy of every valid particle
  this->clearParticles();                                    // empty the holders (capacity retained)
  this->distributeFromPool(pool);                            // re-assign each to its owner

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

  for (DataIterator dit = m_grids[a_level].dataIterator(); dit.ok(); ++dit) {
    maskLevel[dit()].clear();

    const BaseFab<bool>& mask = a_mask[dit()];
    if (!mask.isUsable()) {
      continue;
    }
    const Box maskBox = mask.box();

    const Leaf& src = srcLevel[dit()];
    Leaf&       dst = maskLevel[dit()];
    for (std::size_t i = 0; i < src.size(); i++) {
      const IntVect iv = this->cellIndex(a_level, src.position(i));
      if (maskBox.contains(iv) && mask(iv, 0)) {
        dst.append(src.position(i), src.weight(i), src.gather(i));
        const std::size_t k = dst.size() - 1;
        dst.particleID(k)   = src.particleID(i);
        dst.rankID(k)       = src.rankID(i);
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

  for (DataIterator dit = m_grids[a_level].dataIterator(); dit.ok(); ++dit) {
    const BaseFab<bool>& mask = a_mask[dit()];
    if (!mask.isUsable()) {
      continue;
    }
    const Box maskBox = mask.box();

    Leaf& src = srcLevel[dit()];
    Leaf& dst = maskLevel[dit()];

    std::size_t i = 0;
    while (i < src.size()) {
      const IntVect iv = this->cellIndex(a_level, src.position(i));
      if (maskBox.contains(iv) && mask(iv, 0)) {
        dst.append(src.position(i), src.weight(i), src.gather(i));
        const std::size_t k = dst.size() - 1;
        dst.particleID(k)   = src.particleID(i);
        dst.rankID(k)       = src.rankID(i);
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
