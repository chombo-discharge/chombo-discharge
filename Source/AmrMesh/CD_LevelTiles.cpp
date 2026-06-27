/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_LevelTiles.cpp
  @brief  Implementation of CD_LevelTiles.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <BoxLayout.H>
#include <BoxIterator.H>
#include <LevelData.H>

// Our includes
#include <CD_LevelTiles.H>
#include <CD_NamespaceHeader.H>

LevelTiles::LevelTiles() noexcept : m_isDefined(false)
{
  CH_TIME("LevelTiles::LevelTiles(weak)");
}

LevelTiles::LevelTiles(const DisjointBoxLayout& a_dbl, const int a_minBlockSize) noexcept
{
  CH_TIME("LevelTiles::LevelTiles(full)");

  this->define(a_dbl, a_minBlockSize);
}

LevelTiles::~LevelTiles() noexcept
{
  CH_TIME("LevelTiles::~LevelTiles");
}

void
LevelTiles::define(const DisjointBoxLayout& a_dbl, const int a_minBlockSize) noexcept
{
  CH_TIME("LevelTiles::define");

  const unsigned int myRank = procID();

  m_myTiles.clear();
  m_otherTiles.clear();
  m_myGrids.clear();

  for (LayoutIterator lit = a_dbl.layoutIterator(); lit.ok(); ++lit) {
    const LayoutIndex& lidx   = lit();
    const unsigned int rankID = a_dbl.procID(lidx);
    const unsigned int tileID = a_dbl.index(lidx);

    // A box is a union of aligned min_block_size tiles, so coarsening by the block size gives the
    // (possibly multi-tile) range of min-tiles it covers. Register the box under every one of them.
    // When min_block_size == max_box_size this range is a single tile (the one-tile-per-box fast path).
    const Box tileRange = coarsen(a_dbl[lidx], a_minBlockSize);

    for (BoxIterator bit(tileRange); bit.ok(); ++bit) {
      const IntVect tile = bit();

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
