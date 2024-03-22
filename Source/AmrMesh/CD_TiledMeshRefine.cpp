/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TiledMeshRefine.cpp
  @brief  Implementation of CD_TiledMeshRefine.H
  @author Robert Marskar
*/

// Std includes
#include <set>

// Chombo includes
#include <BoxIterator.H>
#include <CH_Timer.H>

// Our includes
#include <CD_Tile.H>
#include <CD_ParallelOps.H>
#include <CD_TiledMeshRefine.H>
#include <CD_NamespaceHeader.H>

TiledMeshRefine::TiledMeshRefine(const ProblemDomain& a_coarsestDomain,
                                 const Vector<int>&   a_refRatios,
                                 const IntVect&       a_tileSize) noexcept
{
  CH_TIME("TiledMeshRefine::TiledMeshRefine");

  m_refRatios = a_refRatios;
  m_tileSize  = a_tileSize;

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
      a_tiles.emplace(Tile(D_DECL(iv[0], iv[1], iv[2])));
    }
  }
  CH_STOP(t1);

  // Gather tiles globally
#ifdef CH_MPI
  CH_START(t2);
  const int mySendCount  = a_tiles.size() * SpaceDim;
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
  for (int i = 0; i < recvCount; i += SpaceDim) {
    a_tiles.emplace(Tile(D_DECL(recvBuffer[i], recvBuffer[i + 1], recvBuffer[i + 2])));
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

    for (BoxIterator bit(box); bit.ok(); ++bit) {
      const IntVect iv = bit();

      a_tiles.emplace(Tile(D_DECL(iv[0], iv[1], iv[2])));
    }
  }

  CH_STOP(t3);
}

void
TiledMeshRefine::makeBoxesFromTiles(Vector<Box>&         a_boxes,
                                    const TileSet&       a_tileSet,
                                    const ProblemDomain& a_domain) const noexcept
{
  CH_TIME("TiledMeshRefine::makeBoxesFromTiles");

  a_boxes.resize(0);

  const IntVect probLo = a_domain.domainBox().smallEnd();

  for (const auto& tile : a_tileSet) {
    IntVect boxLo = probLo;
    for (int dir = 0; dir < SpaceDim; dir++) {
      boxLo[dir] += tile[dir] * m_tileSize[dir];
    }
    const IntVect boxHi = boxLo + m_tileSize - IntVect::Unit;

    const Box box(boxLo, boxHi);

    a_boxes.push_back(Box(boxLo, boxHi));
  }
}

#include <CD_NamespaceFooter.H>
