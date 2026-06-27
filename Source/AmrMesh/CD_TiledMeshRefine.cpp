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
#include <set>
#include <vector>

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_Tile.H>
#include <CD_ParallelOps.H>
#include <CD_TiledMeshRefine.H>
#include <CD_BoxLoops.H>
#include <CD_NamespaceHeader.H>

TiledMeshRefine::TiledMeshRefine(const ProblemDomain& a_coarsestDomain,
                                 const Vector<int>&   a_refRatios,
                                 const IntVect&       a_tileSize,
                                 const IntVect&       a_maxBoxSize) noexcept
  : m_refRatios(a_refRatios), m_tileSize(a_tileSize), m_maxBoxSize(a_maxBoxSize)
{
  CH_TIME("TiledMeshRefine::TiledMeshRefine");

  for (int dir = 0; dir < SpaceDim; dir++) {
    CH_assert(a_maxBoxSize[dir] >= a_tileSize[dir]);
    CH_assert(a_maxBoxSize[dir] % a_tileSize[dir] == 0);
  }

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

    // Extra level of empty tiles so we can use makeLevelTiles for all levels
    std::vector<TileSet> amrTiles(2 + newFinestLevel);

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

    // Make tiles into boxes
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
TiledMeshRefine::makeLevelTiles(TileSet&             a_tiles,
                                const TileSet&       a_fineTiles,
                                const IntVectSet&    a_coarTags,
                                const ProblemDomain& a_domain,
                                const int            a_refToFine,
                                const int            a_refToCoar) const noexcept
{
  CH_TIMERS("TiledMeshRefine::makeLevelTiles");
  CH_TIMER("TiledMeshRefine::makeLevelTiles::tag_tiles", t1);
  CH_TIMER("TiledMeshRefine::makeLevelTiles::gather_tiles", t2);
  CH_TIMER("TiledMeshRefine::makeLevelTiles::add_fine_tiles", t3);

  // Generate tiles from tags on the coarser level.
  CH_START(t1);
  const Box tileBox = Box(IntVect::Zero, a_domain.size() / m_tileSize - IntVect::Unit);

  const ProblemDomain coarDomain = coarsen(a_domain, a_refToCoar);
  const ProblemDomain fineDomain = refine(a_domain, a_refToFine);

  const IntVect coarProbLo = coarDomain.domainBox().smallEnd();
  a_tiles.clear();
  for (IVSIterator ivsIt(a_coarTags); ivsIt.ok(); ++ivsIt) {
    CH_assert(a_refToCoar >= 2);
    CH_assert(a_refToCoar % 2 == 0);

    const IntVect tag = ivsIt();
    const IntVect iv  = IntVect(D_DECL((tag[0] - coarProbLo[0]) / (m_tileSize[0] / a_refToCoar),
                                      (tag[1] - coarProbLo[1]) / (m_tileSize[1] / a_refToCoar),
                                      (tag[2] - coarProbLo[2]) / (m_tileSize[2] / a_refToCoar)));

    if (tileBox.contains(iv)) {
      a_tiles.emplace(D_DECL(iv[0], iv[1], iv[2]));
    }
  }
  CH_STOP(t1);

  // Gather tiles globally
#ifdef CH_MPI
  CH_START(t2);
  const int mySendCount  = static_cast<int>(a_tiles.size()) * SpaceDim;
  int*      mySendBuffer = new int[mySendCount];

  // Get the number of elements sent by each MPI rank and compute the offset array which is required by Allgatherv
  int  recvCount  = mySendCount;
  int* sendCounts = new int[numProc()];
  int* offsets    = new int[numProc()];

  MPI_Allreduce(MPI_IN_PLACE, &recvCount, 1, MPI_INT, MPI_SUM, Chombo_MPI::comm);
  MPI_Allgather(&mySendCount, 1, MPI_INT, sendCounts, 1, MPI_INT, Chombo_MPI::comm);
  offsets[0] = 0;
  for (int i = 0; i < numProc() - 1; i++) {
    offsets[i + 1] = offsets[i] + sendCounts[i];
  }

  // Linearize the tiles as integers onto the send buffer
  int idx = 0;
  for (const auto& t : a_tiles) {
    for (int dir = 0; dir < SpaceDim; dir++, idx++) {
      mySendBuffer[idx] = t[dir];
    }
  }
  a_tiles.clear();

  // Allocate storage and gather tiles
  int* recvBuffer = new int[recvCount];
  MPI_Allgatherv(mySendBuffer, mySendCount, MPI_INT, recvBuffer, sendCounts, offsets, MPI_INT, Chombo_MPI::comm);

  // de-linearize the received data back into tiles
  for (int i = 0; i < recvCount;) {
    a_tiles.emplace(D_DECL(recvBuffer[i], recvBuffer[i + 1], recvBuffer[i + 2]));
    i += SpaceDim;
  }

  delete[] recvBuffer;
  delete[] offsets;
  delete[] sendCounts;
  delete[] mySendBuffer;

  CH_STOP(t2);
#endif

  // Ensure proper nesting by adding coarsened tiles from the fine level. We grow by one tile (one the fine level) in order
  // to ensure that we're nesting correctly.
  CH_START(t3);

  const Box tileBoxFine = refine(tileBox, a_refToFine);

  for (const Tile& fineTile : a_fineTiles) {
    CH_assert(a_refToFine >= 2);
    CH_assert(a_refToFine % 2 == 0);

    const IntVect fineTileIV = IntVect(D_DECL(fineTile[0], fineTile[1], fineTile[2]));
    const Box     fineBox    = grow(Box(fineTileIV, fineTileIV), 1) & tileBoxFine;
    const Box     box        = coarsen(fineBox, a_refToFine);

    BoxLoops::loop<D_DECL(1, 1, 1)>(box, [&](const IntVect& iv) -> void {
      a_tiles.emplace(D_DECL(iv[0], iv[1], iv[2]));
    });
  }

  CH_STOP(t3);
}

void
TiledMeshRefine::makeBoxesFromTiles(Vector<Box>&         a_boxes,
                                    const TileSet&       a_tileSet,
                                    const ProblemDomain& a_domain) const noexcept
{
  CH_TIME("TiledMeshRefine::makeBoxesFromTiles");

  // When the cap equals the tile size there is nothing to merge -> one box per tile (legacy behaviour).
  bool merge = false;
  for (int dir = 0; dir < SpaceDim; dir++) {
    if (m_maxBoxSize[dir] > m_tileSize[dir]) {
      merge = true;
    }
  }

  if (merge) {
    this->mergeTiles(a_boxes, a_tileSet, a_domain);
    return;
  }

  a_boxes.resize(0);

  const IntVect probLo = a_domain.domainBox().smallEnd();

  for (const auto& tile : a_tileSet) {
    IntVect boxLo = probLo;
    for (int dir = 0; dir < SpaceDim; dir++) {
      boxLo[dir] += tile[dir] * m_tileSize[dir];
    }
    const IntVect boxHi = boxLo + m_tileSize - IntVect::Unit;

    a_boxes.push_back(Box(boxLo, boxHi));
  }
}

void
TiledMeshRefine::mergeTiles(Vector<Box>& a_boxes, const TileSet& a_tiles, const ProblemDomain& a_domain) const noexcept
{
  CH_TIME("TiledMeshRefine::mergeTiles");

  a_boxes.resize(0);

  if (a_tiles.empty()) {
    return;
  }

  // Cap in TILE units (the constructor guarantees m_maxBoxSize % m_tileSize == 0).
  IntVect maxTile;
  for (int dir = 0; dir < SpaceDim; dir++) {
    maxTile[dir] = m_maxBoxSize[dir] / m_tileSize[dir];
  }

  // Flatten the (sorted, globally identical) tile set into a contiguous buffer we can partition in place.
  std::vector<IntVect> tiles;
  tiles.reserve(a_tiles.size());
  for (const Tile& t : a_tiles) {
    tiles.emplace_back(D_DECL(t[0], t[1], t[2]));
  }

  IntVect* const T = tiles.data();

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

  struct Work
  {
    std::size_t b, e;
    IntVect     lo, hi;
  };

  std::vector<Work> stack;
  stack.reserve(64);
  {
    IntVect lo, hi;
    bbox(0, tiles.size(), lo, hi);
    stack.push_back({0, tiles.size(), lo, hi});
  }

  const IntVect probLo = a_domain.domainBox().smallEnd();

  std::vector<int> hist; // reused per node

  while (!stack.empty()) {
    const Work w = stack.back();
    stack.pop_back();

    const std::size_t cnt = w.e - w.b;

    long long vol  = 1;
    bool      fits = true;
    for (int d = 0; d < SpaceDim; d++) {
      const int ext = w.hi[d] - w.lo[d] + 1;
      vol *= ext;
      if (ext > maxTile[d]) {
        fits = false;
      }
    }

    if (fits && vol == static_cast<long long>(cnt)) {
      // Fully tagged and within the cap -> emit one (anisotropic) box in cell coordinates.
      const IntVect boxLo = probLo + w.lo * m_tileSize;
      const IntVect boxHi = probLo + (w.hi + IntVect::Unit) * m_tileSize - IntVect::Unit;
      a_boxes.push_back(Box(boxLo, boxHi));
      continue;
    }

    // ---- choose split direction + cut plane (split into [lo, cut) and [cut, hi]) ----
    int dir = 0;
    int cut = 0;

    // Over the cap in some direction? Split the most-over direction at a max-tile multiple.
    int overDir = -1;
    for (int d = 0; d < SpaceDim; d++) {
      if (w.hi[d] - w.lo[d] + 1 > maxTile[d]) {
        if (overDir < 0 || (w.hi[d] - w.lo[d]) > (w.hi[overDir] - w.lo[overDir])) {
          overDir = d;
        }
      }
    }

    if (overDir >= 0) {
      dir = overDir;
      cut = w.lo[dir] + maxTile[dir]; // align cap-cuts to the tile-grid origin
    }
    else {
      // Longest extent (tie -> lowest index).
      for (int d = 1; d < SpaceDim; d++) {
        if ((w.hi[d] - w.lo[d]) > (w.hi[dir] - w.lo[dir])) {
          dir = d;
        }
      }
      const int ext = w.hi[dir] - w.lo[dir] + 1;
      hist.assign(ext, 0);
      for (std::size_t i = w.b; i < w.e; i++) {
        hist[T[i][dir] - w.lo[dir]]++;
      }

      // Prefer a zero-tag gap nearest the middle (deterministic outward scan); else cut at the middle.
      const int mid  = ext / 2;
      int       best = -1;
      for (int delta = 0; delta < ext && best < 0; delta++) {
        const int kLo = mid - delta;
        const int kHi = mid + delta;
        if (kLo > 0 && kLo < ext && hist[kLo] == 0) {
          best = kLo;
        }
        else if (kHi > 0 && kHi < ext && hist[kHi] == 0) {
          best = kHi;
        }
      }
      cut = (best >= 0) ? (w.lo[dir] + best) : (w.lo[dir] + std::max(1, mid));
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
