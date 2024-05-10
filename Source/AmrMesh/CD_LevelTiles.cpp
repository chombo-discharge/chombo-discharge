/* chombo-discharge
 * Copyright © 2024 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_LevelTiles.cpp
  @brief  Implementation of CD_LevelTiles.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <BoxLayout.H>
#include <LevelData.H>

// Our includes
#include <CD_LevelTiles.H>
#include <CD_NamespaceHeader.H>

LevelTiles::LevelTiles() noexcept
{
  CH_TIME("LevelTiles::LevelTiles(weak)");

  m_isDefined = false;
}

LevelTiles::LevelTiles(const DisjointBoxLayout& a_dbl, const int a_blockingFactor) noexcept
{
  CH_TIME("LevelTiles::LevelTiles(full)");

  this->define(a_dbl, a_blockingFactor);
}

LevelTiles::~LevelTiles() noexcept
{
  CH_TIME("LevelTiles::~LevelTiles");
}

void
LevelTiles::define(const DisjointBoxLayout& a_dbl, const int a_blockingFactor) noexcept
{
  CH_TIME("LevelTiles::define");

  CH_assert(a_dbl.isDefined());

  const unsigned int numRanks = numProc();
  const unsigned int myRank   = procID();

  m_myTiles.clear();
  m_otherTiles.clear();
  m_myGrids.clear();

  for (LayoutIterator lit = a_dbl.layoutIterator(); lit.ok(); ++lit) {
    const LayoutIndex  lidx   = lit();
    const IntVect      tile   = coarsen(a_dbl[lidx], a_blockingFactor).smallEnd();
    const unsigned int rankID = a_dbl.procID(lidx);
    const unsigned int tileID = a_dbl.index(lidx);

    // If using MPI we need to figure out who owns this tile.
#if CH_MPI
    if (myRank == rankID) {
      m_myTiles[tile] = tileID;
    }
    else {
      m_otherTiles[tile] = std::make_pair(tileID, rankID);
    }
#else
    m_myTiles[tile] = tileID;
#endif
  }

  // Figure out which global indices correspond to which local indices.
  for (DataIterator dit(a_dbl); dit.ok(); ++dit) {
    m_myGrids[a_dbl.index(dit())] = dit();
  }

  m_isDefined = true;
}

const std::map<IntVect, unsigned int, LevelTiles::TileComparator>&
LevelTiles::getMyTiles() const noexcept
{
  CH_assert(m_isDefined);

  if (!m_isDefined) {
    MayDay::Abort("define snuck out");
  }

  return m_myTiles;
}

const std::map<IntVect, LevelTiles::BoxIDs, LevelTiles::TileComparator>&
LevelTiles::getOtherTiles() const noexcept
{
  CH_assert(m_isDefined);

  return m_otherTiles;
}

const std::map<unsigned int, DataIndex>&
LevelTiles::getMyGrids() const noexcept
{
  CH_assert(m_isDefined);

  return m_myGrids;
}

#include <CD_NamespaceFooter.H>
