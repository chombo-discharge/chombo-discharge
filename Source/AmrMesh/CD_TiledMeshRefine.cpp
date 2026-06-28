/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_TiledMeshRefine.cpp
  @brief  Implementation of CD_TiledMeshRefine.H
  @author Robert Marskar
*/

// Std includes
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <vector>

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_ParallelOps.H>
#include <CD_TiledMeshRefine.H>
#include <CD_NamespaceHeader.H>

TiledMeshRefine::TiledMeshRefine(const ProblemDomain& a_coarsestDomain,
                                 const Vector<int>&   a_refRatios,
                                 const IntVect&       a_tileSize,
                                 const IntVect&       a_maxBlockSize) noexcept
  : m_refRatios(a_refRatios), m_tileSize(a_tileSize), m_maxBlockSize(a_maxBlockSize)
{
  CH_TIME("TiledMeshRefine::TiledMeshRefine");

  m_superVol = 1;
  for (int dir = 0; dir < SpaceDim; dir++) {
    CH_assert(a_maxBlockSize[dir] >= a_tileSize[dir]);
    CH_assert(a_maxBlockSize[dir] % a_tileSize[dir] == 0);

    m_superFactor[dir] = a_maxBlockSize[dir] / a_tileSize[dir];
    m_superVol *= m_superFactor[dir];
  }
  m_bitmaskWords = static_cast<int>((m_superVol + 63) / 64);

  m_amrDomains.resize(0);

  m_amrDomains.push_back(a_coarsestDomain);
  for (int lvl = 1; lvl < a_refRatios.size(); lvl++) {
    CH_assert(a_refRatios[lvl - 1] >= 2);
    CH_assert(a_refRatios[lvl - 1] % 2 == 0);

    const ProblemDomain domainFine = refine(m_amrDomains[lvl - 1], m_refRatios[lvl - 1]);

    m_amrDomains.push_back(domainFine);
  }
}

TiledMeshRefine::~TiledMeshRefine() noexcept
{
  CH_TIME("TiledMeshRefine::~TiledMeshRefine");
}

std::uint64_t
TiledMeshRefine::encodeSuper(const IntVect& a_super) const noexcept
{
  std::uint64_t key = 0;
  for (int dir = 0; dir < SpaceDim; dir++) {
    key |= (static_cast<std::uint64_t>(a_super[dir]) & 0x1FFFFF) << (21 * dir);
  }
  return key;
}

IntVect
TiledMeshRefine::decodeSuper(const std::uint64_t a_key) const noexcept
{
  IntVect super;
  for (int dir = 0; dir < SpaceDim; dir++) {
    super[dir] = static_cast<int>((a_key >> (21 * dir)) & 0x1FFFFF);
  }
  return super;
}

int
TiledMeshRefine::subIndex(const IntVect& a_sub) const noexcept
{
  int lin    = 0;
  int stride = 1;
  for (int dir = 0; dir < SpaceDim; dir++) {
    lin += a_sub[dir] * stride;
    stride *= m_superFactor[dir];
  }
  return lin;
}

void
TiledMeshRefine::addFineTile(SuperTiles& a_tiles, const IntVect& a_fineTile) const noexcept
{
  IntVect super;
  IntVect sub;
  for (int dir = 0; dir < SpaceDim; dir++) {
    super[dir] = a_fineTile[dir] / m_superFactor[dir];
    sub[dir]   = a_fineTile[dir] - super[dir] * m_superFactor[dir];
  }

  const std::uint64_t key = this->encodeSuper(super);
  if (a_tiles.m_full.count(key) > 0) {
    return; // already fully tagged
  }

  std::vector<std::uint64_t>& bm = a_tiles.m_partial[key];
  if (bm.empty()) {
    bm.assign(m_bitmaskWords, 0);
  }

  const int lin = this->subIndex(sub);
  bm[lin >> 6] |= (1ull << (lin & 63));
}

void
TiledMeshRefine::classify(SuperTiles& a_tiles) const noexcept
{
  for (auto it = a_tiles.m_partial.begin(); it != a_tiles.m_partial.end();) {
    long long pop = 0;
    for (const std::uint64_t w : it->second) {
      pop += __builtin_popcountll(w);
    }

    if (pop == m_superVol) {
      a_tiles.m_full.insert(it->first);
      it = a_tiles.m_partial.erase(it);
    }
    else {
      ++it;
    }
  }
}

void
TiledMeshRefine::gatherSuperTiles(SuperTiles& a_tiles) const noexcept
{
  CH_TIME("TiledMeshRefine::gatherSuperTiles");

#ifdef CH_MPI
  // Promote locally-full super-tiles first: a super-tile all of whose sub-tiles are tagged on this rank is
  // globally full (tags partition the domain by rank, so a super's tagged sub-tiles are disjoint across
  // ranks and OR-ing is exact). Those are then sent as a single key instead of an all-ones bitmask.
  this->classify(a_tiles);

  const int nProc = numProc();

  // All-gather a uint64 buffer (every rank's data, including this rank's), returning the concatenation.
  auto allGather = [&](const std::vector<unsigned long long>& a_send) {
    const int        mySendCount = static_cast<int>(a_send.size());
    int              recvCount   = mySendCount;
    std::vector<int> sendCounts(nProc);
    std::vector<int> offsets(nProc);

    MPI_Allreduce(MPI_IN_PLACE, &recvCount, 1, MPI_INT, MPI_SUM, Chombo_MPI::comm);
    MPI_Allgather(&mySendCount, 1, MPI_INT, sendCounts.data(), 1, MPI_INT, Chombo_MPI::comm);
    offsets[0] = 0;
    for (int i = 0; i < nProc - 1; i++) {
      offsets[i + 1] = offsets[i] + sendCounts[i];
    }

    std::vector<unsigned long long> recv(recvCount);
    MPI_Allgatherv(a_send.data(),
                   mySendCount,
                   MPI_UNSIGNED_LONG_LONG,
                   recv.data(),
                   sendCounts.data(),
                   offsets.data(),
                   MPI_UNSIGNED_LONG_LONG,
                   Chombo_MPI::comm);
    return recv;
  };

  // Stream 1 -- full super-tiles, one key each (the "full marker", no bitmask).
  const std::vector<unsigned long long> sendFull(a_tiles.m_full.begin(), a_tiles.m_full.end());
  const std::vector<unsigned long long> recvFull = allGather(sendFull);

  // Stream 2 -- partial super-tiles, key + sub-occupancy bitmask each.
  const int                       R = 1 + m_bitmaskWords;
  std::vector<unsigned long long> sendPartial;
  sendPartial.reserve(a_tiles.m_partial.size() * R);
  for (const auto& kv : a_tiles.m_partial) {
    sendPartial.push_back(kv.first);
    for (int w = 0; w < m_bitmaskWords; w++) {
      sendPartial.push_back(kv.second[w]);
    }
  }
  const std::vector<unsigned long long> recvPartial = allGather(sendPartial);

  // Rebuild the global representation: full markers are full; partial bitmasks are OR-ed per key, except
  // for keys already known full. (Partials that OR to all-ones across ranks are promoted by the caller's
  // subsequent classify().)
  a_tiles.m_full.clear();
  a_tiles.m_partial.clear();
  for (const unsigned long long key : recvFull) {
    a_tiles.m_full.insert(key);
  }
  for (std::size_t i = 0; i < recvPartial.size(); i += R) {
    const std::uint64_t key = recvPartial[i];
    if (a_tiles.m_full.count(key) > 0) {
      continue; // already full via a marker from some rank
    }
    std::vector<std::uint64_t>& bm = a_tiles.m_partial[key];
    if (bm.empty()) {
      bm.assign(m_bitmaskWords, 0);
    }
    for (int w = 0; w < m_bitmaskWords; w++) {
      bm[w] |= recvPartial[i + 1 + w];
    }
  }
#endif
}

void
TiledMeshRefine::nestFrom(SuperTiles&       a_tiles,
                          const SuperTiles& a_finer,
                          const int         a_refToFine,
                          const Box&        a_thisTileBox) const noexcept
{
  const Box fineTileBox = refine(a_thisTileBox, a_refToFine);

  // Add every this-level tile in coarsen(grow(a_fineBox, 1)) (clipped to the domains).
  auto addNest = [&](const Box& a_fineBox) {
    const Box grown     = grow(a_fineBox, 1) & fineTileBox;
    const Box coarsened = coarsen(grown, a_refToFine) & a_thisTileBox;
    for (BoxIterator bit(coarsened); bit.ok(); ++bit) {
      this->addFineTile(a_tiles, bit());
    }
  };

  // Full finer super-tiles: their whole block, grown and coarsened, in bulk.
  for (const std::uint64_t key : a_finer.m_full) {
    const IntVect super = this->decodeSuper(key);
    const IntVect blkLo = super * m_superFactor;
    const IntVect blkHi = blkLo + m_superFactor - IntVect::Unit;
    addNest(Box(blkLo, blkHi));
  }

  // Partial finer super-tiles: each tagged sub-tile individually (the finer boundary).
  for (const auto& kv : a_finer.m_partial) {
    const IntVect super = this->decodeSuper(kv.first);
    for (int lin = 0; lin < m_superVol; lin++) {
      if ((kv.second[lin >> 6] >> (lin & 63)) & 1ull) {
        IntVect sub;
        int     r = lin;
        for (int dir = 0; dir < SpaceDim; dir++) {
          sub[dir] = r % m_superFactor[dir];
          r /= m_superFactor[dir];
        }
        const IntVect fineTile = super * m_superFactor + sub;
        addNest(Box(fineTile, fineTile));
      }
    }
  }
}

int
TiledMeshRefine::regrid(Vector<Vector<Box>>& a_newGrids, const Vector<IntVectSet>& a_tags) const noexcept
{
  CH_TIME("TiledMeshRefine::regrid");

  // Figure out the highest level which has tags. This needs to be the same for every rank.
  int topLevel = 0;

  for (int lvl = 0; lvl < a_tags.size(); lvl++) {
    if (!a_tags[lvl].isEmpty()) {
      topLevel = lvl;
    }
  }

  topLevel = ParallelOps::max(topLevel);

  int newFinestLevel = 0;

  if (topLevel >= 0) {
    newFinestLevel = 1 + topLevel;

    // Extra (empty) finer level so the finest level can use the same makeLevelTiles call.
    std::vector<SuperTiles> amrTiles(2 + newFinestLevel);

    for (int lvl = newFinestLevel; lvl > 0; lvl--) {
      this->makeLevelTiles(amrTiles[lvl],
                           amrTiles[lvl + 1],
                           a_tags[lvl - 1],
                           m_amrDomains[lvl],
                           m_refRatios[lvl],
                           m_refRatios[lvl - 1]);
    }

    // Coarsest grid just consists of proper nesting around the finer grids. Last argument is dummy.
    this->makeLevelTiles(amrTiles[0], amrTiles[1], IntVectSet(), m_amrDomains[0], m_refRatios[0], 1);

    // Make the super-tiles into boxes.
    a_newGrids.resize(1 + newFinestLevel);
    for (int lvl = 0; lvl <= newFinestLevel; lvl++) {
      this->makeBoxesFromTiles(a_newGrids[lvl], amrTiles[lvl], m_amrDomains[lvl]);
    }
  }
  else {
    a_newGrids.resize(0);
  }

  return newFinestLevel;
}

void
TiledMeshRefine::makeLevelTiles(SuperTiles&          a_tiles,
                                const SuperTiles&    a_fineTiles,
                                const IntVectSet&    a_coarTags,
                                const ProblemDomain& a_domain,
                                const int            a_refToFine,
                                const int            a_refToCoar) const noexcept
{
  CH_TIMERS("TiledMeshRefine::makeLevelTiles");
  CH_TIMER("TiledMeshRefine::makeLevelTiles::tag_tiles", t1);
  CH_TIMER("TiledMeshRefine::makeLevelTiles::gather_tiles", t2);
  CH_TIMER("TiledMeshRefine::makeLevelTiles::add_fine_tiles", t3);

  a_tiles.m_full.clear();
  a_tiles.m_partial.clear();

  const Box tileBox = Box(IntVect::Zero, a_domain.size() / m_tileSize - IntVect::Unit);

  const ProblemDomain coarDomain = coarsen(a_domain, a_refToCoar);
  const IntVect       coarProbLo = coarDomain.domainBox().smallEnd();

  // Generate (per-rank-local) super-tiles from tags on the coarser level.
  CH_START(t1);
  for (IVSIterator ivsIt(a_coarTags); ivsIt.ok(); ++ivsIt) {
    CH_assert(a_refToCoar >= 2);
    CH_assert(a_refToCoar % 2 == 0);

    const IntVect tag = ivsIt();
    const IntVect iv  = IntVect(D_DECL((tag[0] - coarProbLo[0]) / (m_tileSize[0] / a_refToCoar),
                                      (tag[1] - coarProbLo[1]) / (m_tileSize[1] / a_refToCoar),
                                      (tag[2] - coarProbLo[2]) / (m_tileSize[2] / a_refToCoar)));

    if (tileBox.contains(iv)) {
      this->addFineTile(a_tiles, iv);
    }
  }
  CH_STOP(t1);

  // Gather the compact super-tile representation onto all ranks, then classify full vs partial globally.
  CH_START(t2);
  this->gatherSuperTiles(a_tiles);
  this->classify(a_tiles);
  CH_STOP(t2);

  // Ensure proper nesting by injecting the coarsened, grown-by-one buffer from the finer level.
  CH_START(t3);
  this->nestFrom(a_tiles, a_fineTiles, a_refToFine, tileBox);
  this->classify(a_tiles);
  CH_STOP(t3);
}

void
TiledMeshRefine::makeBoxesFromTiles(Vector<Box>&         a_boxes,
                                    const SuperTiles&    a_tiles,
                                    const ProblemDomain& a_domain) const noexcept
{
  CH_TIME("TiledMeshRefine::makeBoxesFromTiles");

  a_boxes.resize(0);

  const IntVect probLo = a_domain.domainBox().smallEnd();

  // Visit super-tiles in sorted-key order so the box list is identical on every rank.
  std::vector<std::uint64_t> fullKeys(a_tiles.m_full.begin(), a_tiles.m_full.end());
  std::sort(fullKeys.begin(), fullKeys.end());

  std::vector<std::uint64_t> partialKeys;
  partialKeys.reserve(a_tiles.m_partial.size());
  for (const auto& kv : a_tiles.m_partial) {
    partialKeys.push_back(kv.first);
  }
  std::sort(partialKeys.begin(), partialKeys.end());

  // Full super-tiles -> one (max_block_size) box each.
  for (const std::uint64_t key : fullKeys) {
    const IntVect super = this->decodeSuper(key);
    const IntVect boxLo = probLo + (super * m_superFactor) * m_tileSize;
    const IntVect boxHi = probLo + ((super + IntVect::Unit) * m_superFactor) * m_tileSize - IntVect::Unit;
    a_boxes.push_back(Box(boxLo, boxHi));
  }

  // Partial super-tiles -> pack their tagged sub-tiles into variable-sized boxes.
  std::vector<IntVect> members;
  for (const std::uint64_t key : partialKeys) {
    const std::vector<std::uint64_t>& bm    = a_tiles.m_partial.at(key);
    const IntVect                     super = this->decodeSuper(key);

    members.clear();
    for (int lin = 0; lin < m_superVol; lin++) {
      if ((bm[lin >> 6] >> (lin & 63)) & 1ull) {
        IntVect sub;
        int     r = lin;
        for (int dir = 0; dir < SpaceDim; dir++) {
          sub[dir] = r % m_superFactor[dir];
          r /= m_superFactor[dir];
        }
        members.push_back(super * m_superFactor + sub);
      }
    }

    this->packTiles(members, m_superFactor, probLo, a_boxes);
  }
}

void
TiledMeshRefine::packTiles(std::vector<IntVect>& a_tiles,
                           const IntVect&        a_maxTile,
                           const IntVect&        a_probLo,
                           Vector<Box>&          a_boxes) const noexcept
{
  CH_TIME("TiledMeshRefine::packTiles");

  if (a_tiles.empty()) {
    return;
  }

  IntVect* const T = a_tiles.data();

  // Bounding box (in tile coordinates) of a contiguous range [b, e).
  auto bbox = [&](const std::size_t b, const std::size_t e, IntVect& lo, IntVect& hi) {
    lo = T[b];
    hi = T[b];
    for (std::size_t i = b + 1; i < e; i++) {
      for (int d = 0; d < SpaceDim; d++) {
        lo[d] = std::min(lo[d], T[i][d]);
        hi[d] = std::max(hi[d], T[i][d]);
      }
    }
  };

  // Number of tagged tiles in the slab coord[dir] == lo[dir] + p over the range [b, e). The signatures
  // used by the Berger-Rigoutsos split are these slab counts as a function of p.
  std::vector<int> hist;
  auto signature = [&](const std::size_t b, const std::size_t e, const IntVect& lo, const int ext, const int dir) {
    hist.assign(ext, 0);
    for (std::size_t i = b; i < e; i++) {
      hist[T[i][dir] - lo[dir]]++;
    }
  };

  struct Work
  {
    std::size_t b, e;
    IntVect     lo, hi;
  };

  std::vector<Work> stack;
  stack.reserve(64);
  {
    IntVect lo, hi;
    bbox(0, a_tiles.size(), lo, hi);
    stack.push_back({0, a_tiles.size(), lo, hi});
  }

  while (!stack.empty()) {
    const Work w = stack.back();
    stack.pop_back();

    const std::size_t cnt = w.e - w.b;

    long long vol     = 1;
    int       overDir = -1;
    for (int d = 0; d < SpaceDim; d++) {
      const int ext = w.hi[d] - w.lo[d] + 1;
      vol *= ext;
      if (ext > a_maxTile[d]) {
        if (overDir < 0 || (w.hi[d] - w.lo[d]) > (w.hi[overDir] - w.lo[overDir])) {
          overDir = d;
        }
      }
    }

    if (overDir < 0 && vol == static_cast<long long>(cnt)) {
      // Fully tagged and within the cap -> emit one (anisotropic) box in cell coordinates.
      const IntVect boxLo = a_probLo + w.lo * m_tileSize;
      const IntVect boxHi = a_probLo + (w.hi + IntVect::Unit) * m_tileSize - IntVect::Unit;
      a_boxes.push_back(Box(boxLo, boxHi));
      continue;
    }

    // ---- choose split direction + cut plane (split into [lo, cut) and [cut, hi]) ----
    int dir = 0;
    int cut = 0;

    if (overDir >= 0) {
      // Over the cap -> split the most-over direction at a cap multiple (aligned to the tile-grid origin).
      dir = overDir;
      cut = w.lo[dir] + a_maxTile[dir];
    }
    else {
      // Berger-Rigoutsos: prefer an empty slab (zero-tag gap) in ANY direction, nearest the middle.
      int holeDir  = -1;
      int holeCut  = -1;
      int holeDist = std::numeric_limits<int>::max();
      for (int d = 0; d < SpaceDim; d++) {
        const int ext = w.hi[d] - w.lo[d] + 1;
        if (ext < 2) {
          continue;
        }
        signature(w.b, w.e, w.lo, ext, d);
        const int mid = ext / 2;
        for (int p = 1; p < ext; p++) {
          if (hist[p] == 0) {
            const int dist = std::abs(p - mid);
            if (dist < holeDist) {
              holeDist = dist;
              holeDir  = d;
              holeCut  = w.lo[d] + p;
            }
          }
        }
      }

      if (holeDir >= 0) {
        dir = holeDir;
        cut = holeCut;
      }
      else {
        // No hole: cut at the strongest signature inflection (max |second difference|) across directions.
        long long bestLap = -1;
        int       bestDir = -1;
        int       bestCut = -1;
        for (int d = 0; d < SpaceDim; d++) {
          const int ext = w.hi[d] - w.lo[d] + 1;
          if (ext < 3) {
            continue;
          }
          signature(w.b, w.e, w.lo, ext, d);
          for (int p = 1; p < ext - 1; p++) {
            const long long lap = std::llabs(static_cast<long long>(hist[p - 1]) - 2 * static_cast<long long>(hist[p]) +
                                             static_cast<long long>(hist[p + 1]));
            if (lap > bestLap) {
              bestLap = lap;
              bestDir = d;
              bestCut = w.lo[d] + p + 1; // cut just past the inflection
            }
          }
        }

        if (bestDir >= 0) {
          dir = bestDir;
          cut = bestCut;
        }
        else {
          // Degenerate (e.g. a thin region) -> longest-axis midpoint.
          for (int d = 1; d < SpaceDim; d++) {
            if ((w.hi[d] - w.lo[d]) > (w.hi[dir] - w.lo[dir])) {
              dir = d;
            }
          }
          cut = w.lo[dir] + std::max(1, (w.hi[dir] - w.lo[dir] + 1) / 2);
        }
      }
    }

    // Partition [b, e) by coord[dir] < cut.
    IntVect* const    beg  = T + w.b;
    IntVect* const    end  = T + w.e;
    IntVect* const    m    = std::partition(beg, end, [dir, cut](const IntVect& t) {
      return t[dir] < cut;
    });
    const std::size_t mIdx = static_cast<std::size_t>(m - T);

    if (mIdx > w.b) {
      IntVect lo, hi;
      bbox(w.b, mIdx, lo, hi);
      stack.push_back({w.b, mIdx, lo, hi});
    }
    if (mIdx < w.e) {
      IntVect lo, hi;
      bbox(mIdx, w.e, lo, hi);
      stack.push_back({mIdx, w.e, lo, hi});
    }
  }
}

#include <CD_NamespaceFooter.H>
