/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_LevelTiles.H
  @brief  Declaration of LevelTiles
  @author Robert Marskar
*/

#ifndef CD_LEVELTILES_H
#define CD_LEVELTILES_H

// Std includes
#include <cstddef>
#include <unordered_map>
#include <vector>

// Chombo includes
#include <IntVect.H>
#include <DisjointBoxLayout.H>
#include <DataIndex.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Class for storing the AMR hierarchy as a collection of tiles.
  @details This class decomposes the AMR hierarchy into uniform tiles, each tile being the cells spanned by a
  minBlockSize^SpaceDim grid patch, and maps each tile to the grid box that contains it. Boxes may be larger than
  one tile (variable-sized/anisotropic boxes that are unions of aligned minBlockSize tiles): such a box is
  registered under every tile it covers, so the tile->box map stays uniform and a point->tile->box lookup is O(1).
  When minBlockSize == maxBoxSize every box is exactly one tile (the one-tile-per-box fast path).
*/
class LevelTiles
{
public:
  /** @brief alias type alias. */
  using BoxIDs = std::pair<unsigned int, unsigned int>;

  /**
    @brief Hash functor for using IntVect tiles as unordered_map keys.
    @details Hashes all SpaceDim components (FNV-1a), with no bit-truncation, so tiles outside the
    registered range (e.g. far-outcast particles) cannot alias a registered tile.
  */
  struct TileHasher
  {
    /**
      @brief Hash operator.
      @param[in] a_tile Tile to hash.
      @return Hash value.
    */
    std::size_t
    operator()(const IntVect& a_tile) const noexcept
    {
      std::size_t hash = 1469598103934665603ULL;
      for (int dir = 0; dir < SpaceDim; dir++) {
        hash ^= static_cast<std::size_t>(static_cast<unsigned int>(a_tile[dir]));
        hash *= 1099511628211ULL;
      }
      return hash;
    }
  };

  /**
    @brief Weak constructor. Must subsequently call the define function
  */
  LevelTiles() noexcept;

  /**
    @brief Full constructor - puts object in usable state.
    @param[in] a_dbl            Grids
    @param[in] a_minBlockSize Grid size.
  */
  LevelTiles(const DisjointBoxLayout& a_dbl, const int a_minBlockSize) noexcept;

  /**
    @brief Destructor (does nothing)
  */
  virtual ~LevelTiles() noexcept;

  /**
    @brief Define function. Puts object in usable state.
    @param[in] a_dbl            Grids
    @param[in] a_minBlockSize Grid size.
  */
  virtual void
  define(const DisjointBoxLayout& a_dbl, const int a_minBlockSize) noexcept;

  /**
    @brief Get the tiles owned by this rank
    @return m_myTiles
  */
  virtual const std::unordered_map<IntVect, unsigned int, TileHasher>&
  getMyTiles() const noexcept;

  /**
    @brief Get the tiles owned by other ranks
    @return m_otherTiles
  */
  virtual const std::unordered_map<IntVect, BoxIDs, TileHasher>&
  getOtherTiles() const noexcept;

  /**
    @brief Get the mapping of grid index to DataIndex
    @return m_myGrids
  */
  virtual const std::unordered_map<unsigned int, DataIndex>&
  getMyGrids() const noexcept;

protected:
  /**
    @brief Is defined or not
  */
  bool m_isDefined;

  /**
    @brief "Tiles" owned by this rank
  */
  std::unordered_map<IntVect, unsigned int, TileHasher> m_myTiles;

  /**
    @brief Grids owned by other ranks.
  */
  std::unordered_map<IntVect, BoxIDs, TileHasher> m_otherTiles;

  /**
    @brief Mapping of grid index to DataIndex
  */
  std::unordered_map<unsigned int, DataIndex> m_myGrids;
};

#include <CD_NamespaceFooter.H>

#endif
