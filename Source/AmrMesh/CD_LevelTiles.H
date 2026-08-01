/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_LevelTiles.H
 * @brief  Declaration of LevelTiles
 * @author Robert Marskar
 */

#ifndef CD_LEVELTILES_H
#define CD_LEVELTILES_H

// Std includes
#include <cmath>
#include <cstddef>
#include <unordered_map>
#include <vector>

// Chombo includes
#include <IntVect.H>
#include <RealVect.H>
#include <Vector.H>
#include <RefCountedPtr.H>
#include <SPMD.H>
#include <DisjointBoxLayout.H>
#include <DataIndex.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
 * @brief Class for storing the AMR hierarchy as a collection of tiles.
 * @details This class decomposes the AMR hierarchy into uniform tiles, each tile being the cells spanned by a
 * minBlockSize^SpaceDim grid patch, and maps each tile to the grid box that contains it. Boxes may be larger than
 * one tile (variable-sized/anisotropic boxes that are unions of aligned minBlockSize tiles): such a box is
 * registered under every tile it covers, so the tile->box map stays uniform and a point->tile->box lookup is O(1).
 * When minBlockSize == maxBlockSize every box is exactly one tile (the one-tile-per-box fast path).
 */
class LevelTiles
{
public:
  /** @brief alias type alias. */
  using BoxIDs = std::pair<unsigned int, unsigned int>;

  /**
   * @brief Hash functor for using IntVect tiles as unordered_map keys.
   * @details Hashes all SpaceDim components (FNV-1a). Hash collisions are harmless: std::unordered_map
   * distinguishes entries by key equality (IntVect::operator==), not by hash, so a collision only
   * affects performance and never returns a wrong value.
   */
  struct TileHasher
  {
    /**
     * @brief Hash operator.
     * @param[in] a_tile Tile to hash.
     * @return Hash value.
     */
    std::size_t
    operator()(const IntVect& a_tile) const noexcept
    {
      std::size_t hash = 14695981039346656037ULL; // FNV-1a 64-bit offset basis
      for (int dir = 0; dir < SpaceDim; dir++) {
        hash ^= static_cast<std::size_t>(static_cast<unsigned int>(a_tile[dir]));
        hash *= 1099511628211ULL;
      }
      return hash;
    }
  };

  /**
   * @brief Weak constructor. Must subsequently call the define function
   */
  LevelTiles() noexcept;

  /**
   * @brief Full constructor - puts object in usable state.
   * @param[in] a_dbl            Grids
   * @param[in] a_minBlockSize Grid size.
   */
  LevelTiles(const DisjointBoxLayout& a_dbl, const int a_minBlockSize) noexcept;

  /**
   * @brief Destructor (does nothing)
   */
  virtual ~LevelTiles() noexcept;

  /**
   * @brief Define function. Puts object in usable state.
   * @param[in] a_dbl            Grids
   * @param[in] a_minBlockSize Grid size.
   */
  virtual void
  define(const DisjointBoxLayout& a_dbl, const int a_minBlockSize) noexcept;

  /**
   * @brief Get the tiles owned by this rank
   * @return m_myTiles
   */
  virtual const std::unordered_map<IntVect, unsigned int, TileHasher>&
  getMyTiles() const noexcept;

  /**
   * @brief Get the tiles owned by other ranks
   * @return m_otherTiles
   */
  virtual const std::unordered_map<IntVect, BoxIDs, TileHasher>&
  getOtherTiles() const noexcept;

  /**
   * @brief Get the mapping of grid index to DataIndex
   * @return m_myGrids
   */
  virtual const std::unordered_map<unsigned int, DataIndex>&
  getMyGrids() const noexcept;

  /**
   * @brief Result of a point->block query. See findDestination.
   */
  struct LevelAndBox
  {
    /**
     * @brief Owning AMR level (finest tile that contains the point), or -1 if not found.
     */
    int level;

    /**
     * @brief Global grid/box index of the covering box within that level's LevelTiles.
     */
    unsigned int gridIndex;

    /**
     * @brief MPI rank owning that box.
     */
    int rank;

    /**
     * @brief True if a covering tile/box was found (false when the point is off-domain).
     */
    bool valid;
  };

  /**
   * @brief Map a physical position to its owning (level, grid index, rank) via the finest containing tile.
   * @details Shared, header-inlined point->block core used by both ParticleContainer::findDestination (the
   * per-particle remap hot path) and Realm::getLevelAndBox. The position is mapped to a min-block-size tile by
   * integer division, then looked up in each level's tile->box hash, finest level first. It is O(1) per level (a
   * tile computation plus an unordered_map lookup), never a linear scan over patches, and works for
   * variable-sized/anisotropic boxes because every box is registered under all the min tiles it covers.
   * @param[in] a_pos          Physical position to locate.
   * @param[in] a_probLo       Lower-left corner of the domain.
   * @param[in] a_dx           Per-level grid spacing, indexed [level][dir] (isotropic per level).
   * @param[in] a_minBlockSize Tile size in cells (min_block_size).
   * @param[in] a_levelTiles   Per-level tile->box maps (the single source of truth).
   * @param[in] a_finestLevel  Finest AMR level to search.
   * @return {level, gridIndex, rank, valid}; valid==false (other fields unspecified) if no tile owns the cell.
   */
  static inline LevelAndBox
  findDestination(const RealVect&                          a_pos,
                  const RealVect&                          a_probLo,
                  const Vector<RealVect>&                  a_dx,
                  const int                                a_minBlockSize,
                  const Vector<RefCountedPtr<LevelTiles>>& a_levelTiles,
                  const int                                a_finestLevel) noexcept
  {
    // Finest level whose min-block tile contains the point wins.
    for (int lvl = a_finestLevel; lvl >= 0; lvl--) {
      IntVect tile;

      for (int dir = 0; dir < SpaceDim; dir++) {
        tile[dir] = static_cast<int>(std::floor((a_pos[dir] - a_probLo[dir]) / (a_minBlockSize * a_dx[lvl][dir])));
      }

      const LevelTiles& tiles = *a_levelTiles[lvl];

      const auto& myTiles = tiles.getMyTiles();
      const auto  mit     = myTiles.find(tile);

      if (mit != myTiles.end()) {
        return LevelAndBox{lvl, mit->second, procID(), true};
      }

#ifdef CH_MPI
      const auto& otherTiles = tiles.getOtherTiles();
      const auto  oit        = otherTiles.find(tile);

      if (oit != otherTiles.end()) {
        return LevelAndBox{lvl, oit->second.first, static_cast<int>(oit->second.second), true};
      }
#endif
    }

    return LevelAndBox{-1, 0u, -1, false};
  }

protected:
  /**
   * @brief Is defined or not
   */
  bool m_isDefined;

  /**
   * @brief "Tiles" owned by this rank
   */
  std::unordered_map<IntVect, unsigned int, TileHasher> m_myTiles;

  /**
   * @brief Grids owned by other ranks.
   */
  std::unordered_map<IntVect, BoxIDs, TileHasher> m_otherTiles;

  /**
   * @brief Mapping of grid index to DataIndex
   */
  std::unordered_map<unsigned int, DataIndex> m_myGrids;
};

#include <CD_NamespaceFooter.H>

#endif
