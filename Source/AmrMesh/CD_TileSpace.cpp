/* chombo-discharge
 * Copyright © 2024 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TileSpace.cpp
  @brief  Implementation of CD_TileSpace.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <BoxLayout.H>
#include <LevelData.H>

// Our includes
#include <CD_TileSpace.H>
#include <CD_NamespaceHeader.H>

TileSpace::TileSpace() noexcept
{
  CH_TIME("TileSpace::TileSpace(weak)");

  m_isDefined = false;
}

TileSpace::TileSpace(const Vector<DisjointBoxLayout>& a_grids, const int a_blockingFactor) noexcept
{
  CH_TIME("TileSpace::TileSpace(full)");

  this->define(a_grids, a_blockingFactor);
}

TileSpace::~TileSpace() noexcept
{
  CH_TIME("TileSpace::~TileSpace");
}

void
TileSpace::define(const Vector<DisjointBoxLayout>& a_grids, const int a_blockingFactor) noexcept
{
  CH_TIME("TileSpace::define");

  const unsigned int numRanks = numProc();
  const unsigned int myRank   = procID();

  for (int lvl = 0; lvl < a_grids.size(); lvl++) {
    const DisjointBoxLayout& dbl = a_grids[lvl];

    CH_assert(dbl.isDefined());

    std::map<IntVect, unsigned int, TileSpace::TileComparator>      myTiles;
    std::map<IntVect, TileSpace::BoxIDs, TileSpace::TileComparator> otherTiles;
    std::map<unsigned int, DataIndex>                               myGrids;

    for (LayoutIterator lit = dbl.layoutIterator(); lit.ok(); ++lit) {
      const LayoutIndex  lidx   = lit();
      const IntVect      tile   = coarsen(dbl[lidx], a_blockingFactor).smallEnd();
      const unsigned int rankID = dbl.procID(lidx);
      const unsigned int tileID = dbl.index(lidx);

      // If using MPI we need to figure out who owns this tile.
#if CH_MPI
      if (myRank == rankID) {
        myTiles[tile] = tileID;
      }
      else {
        otherTiles[tile] = std::make_pair(tileID, rankID);
      }
#else
      myTiles[tile] = tileID;
#endif
    }

    // Figure out which global indices correspond to which local indices.
    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      myGrids[dbl.index(dit())] = dit();
    }

    m_myTiles.emplace_back(myTiles);
    m_otherTiles.emplace_back(otherTiles);
    m_myGrids.emplace_back(myGrids);
  }

  m_isDefined = true;
}

const std::vector<std::map<IntVect, unsigned int, TileSpace::TileComparator>>&
TileSpace::getMyTiles() const noexcept
{
  CH_assert(m_isDefined);

  return m_myTiles;
}

const std::vector<std::map<IntVect, TileSpace::BoxIDs, TileSpace::TileComparator>>&
TileSpace::getOtherTiles() const noexcept
{
  CH_assert(m_isDefined);

  return m_otherTiles;
}

const std::vector<std::map<unsigned int, DataIndex>>&
TileSpace::getMyGrids() const noexcept
{
  CH_assert(m_isDefined);

  return m_myGrids;
}

#include <CD_NamespaceFooter.H>
