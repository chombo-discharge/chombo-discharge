/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_TiledMeshRefine.H
  @brief  Header for class TiledMeshRefine
  @author Robert Marskar
*/

#ifndef CD_TILEDMESHREFINE_H
#define CD_TILEDMESHREFINE_H

// Std includes
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Chombo includes
#include <IntVectSet.H>
#include <ProblemDomain.H>
#include <Box.H>
#include <BoxIterator.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Class for generating AMR boxes using a tiling algorithm.
  @details This class provides a scalable method for grid generation where the grids are generated in a
  pre-set tile pattern, by flagging min_block_size tiles rather than cells for refinement.

  Generation is top-down over a grid of super-tiles, where a super-tile is a max_block_size block of
  min_block_size tiles. Each level is represented compactly as the set of fully-tagged super-tiles plus,
  for partially-tagged super-tiles, a bitmask of which sub-tiles are tagged. Fine tiles are therefore
  materialized only at the boundary of the refined region (inside partial super-tiles), and only that
  compact representation is gathered across MPI ranks -- the intermediate memory scales with the refined
  surface plus the box count, not the refined volume. Full super-tiles become one (max_block_size) box each;
  partial super-tiles are packed into variable-sized boxes by packTiles. When max_block_size == min_block_size
  every super-tile is a single tile and the output is one box per tile (legacy behaviour).
*/
class TiledMeshRefine
{
public:
  /**
    @brief Default constructor -- deleted because we never need it
  */
  TiledMeshRefine() = delete;

  /**
    @brief Constructor
    @param[in] a_coarsestDomain Coarsest grid domain
    @param[in] a_refRatios      Refinement ratios
    @param[in] a_tileSize       Tile size (= minimum block size, in cells)
    @param[in] a_maxBlockSize     Maximum box size (in cells). Tiles are merged into boxes no larger than
                                this per direction. Must be a multiple of a_tileSize. When equal to
                                a_tileSize the output is one box per tile (legacy behaviour).
  */
  TiledMeshRefine(const ProblemDomain& a_coarsestDomain,
                  const Vector<int>&   a_refRatios,
                  const IntVect&       a_tileSize,
                  const IntVect&       a_maxBlockSize) noexcept;

  /**
    @brief Destructor (does nothing)
  */
  virtual ~TiledMeshRefine() noexcept;

  /**
    @brief Regrid using the tile clustering algorithm.
    @param[out] a_newGrids  The new grid boxes
    @param[in]  a_tagsLevel Grid tags on each level
    @return Returns the new finest grid level.
  */
  virtual int
  regrid(Vector<Vector<Box>>& a_newGrids, const Vector<IntVectSet>& a_tagsLevel) const noexcept;

protected:
  /**
    @brief Compact representation of one level's tiles over the super-tile grid.
    @details A super-tile is a m_superFactor block of min_block_size tiles. Fully-tagged super-tiles are
    stored by key only (m_full); partially-tagged ones additionally store a bitmask of their tagged
    sub-tiles (m_partial). Super-tile keys are packed IntVects (see encodeSuper / decodeSuper).
  */
  struct SuperTiles
  {
    /** @brief Keys of fully-tagged super-tiles. */
    std::unordered_set<std::uint64_t> m_full;

    /** @brief Partially-tagged super-tiles: key -> sub-tile occupancy bitmask. */
    std::unordered_map<std::uint64_t, std::vector<std::uint64_t>> m_partial;
  };

  /**
    @brief Computational domains on each level
    @note This are the domains for the grid.
  */
  Vector<ProblemDomain> m_amrDomains;

  /**
    @brief Refinement ratios
  */
  Vector<int> m_refRatios;

  /**
    @brief Tile size (= minimum block size, in cells)
  */
  IntVect m_tileSize;

  /**
    @brief Maximum box size (in cells). Merge cap; a multiple of m_tileSize.
  */
  IntVect m_maxBlockSize;

  /**
    @brief Number of min_block_size tiles per super-tile, per direction (= m_maxBlockSize / m_tileSize).
  */
  IntVect m_superFactor;

  /**
    @brief Number of sub-tiles in a super-tile (the product of m_superFactor).
  */
  long long m_superVol;

  /**
    @brief Number of 64-bit words in a super-tile sub-occupancy bitmask.
  */
  int m_bitmaskWords;

  /**
    @brief Pack a super-tile index into a 64-bit key (21 bits per direction).
    @param[in] a_super Super-tile index (non-negative).
    @return Packed key.
  */
  std::uint64_t
  encodeSuper(const IntVect& a_super) const noexcept;

  /**
    @brief Unpack a 64-bit super-tile key into a super-tile index.
    @param[in] a_key Packed key.
    @return Super-tile index.
  */
  IntVect
  decodeSuper(const std::uint64_t a_key) const noexcept;

  /**
    @brief Linear sub-tile index of a sub-position within a super-tile.
    @param[in] a_sub Sub-position in [0, m_superFactor) per direction.
    @return Linear index in [0, m_superVol).
  */
  int
  subIndex(const IntVect& a_sub) const noexcept;

  /**
    @brief Add a fine tile (tile coordinates) to the representation by setting its sub-tile bit.
    @details No-op if the tile's super-tile is already known to be full. Does not reclassify.
    @param[in,out] a_tiles    Representation to add to.
    @param[in]     a_fineTile Tile to add (non-negative tile coordinates).
  */
  void
  addFineTile(SuperTiles& a_tiles, const IntVect& a_fineTile) const noexcept;

  /**
    @brief Promote every partial super-tile whose bitmask is full to the full set.
    @param[in,out] a_tiles Representation to classify.
  */
  void
  classify(SuperTiles& a_tiles) const noexcept;

  /**
    @brief Gather the (per-rank-local) super-tile representation onto all ranks.
    @details Locally-full super-tiles are classified first and gathered as a single key each (a "full
    marker", no bitmask); partial super-tiles are gathered as key + bitmask and OR-ed per key. Because tags
    partition the domain by rank, a super-tile's tagged sub-tiles are disjoint across ranks, so the OR is
    exact and a locally-full super-tile is globally full. The gathered data scales with the number of
    touched super-tiles (full ones costing just a key), not the number of fine tiles. Partial super-tiles
    that OR to all-ones across ranks are promoted by the caller's subsequent classify().
    @param[in,out] a_tiles Local representation in; global representation out.
  */
  void
  gatherSuperTiles(SuperTiles& a_tiles) const noexcept;

  /**
    @brief Inject the proper-nesting buffer into a level from the (classified) finer level.
    @details Reproduces the legacy "grow each finer tile by one and coarsen" buffer exactly, but from the
    compact finer representation: full finer super-tiles contribute coarsen(grow(block,1)) in bulk, partial
    ones contribute coarsen(grow(tile,1)) per tagged sub-tile. Adds tiles to a_tiles; the caller reclassifies.
    @param[in,out] a_tiles      Level being built (this level).
    @param[in]     a_finer      Classified finer-level representation.
    @param[in]     a_refToFine  Refinement ratio from this level to the finer level.
    @param[in]     a_thisTileBox Tile-coordinate bounding box of this level (for clipping).
  */
  void
  nestFrom(SuperTiles&       a_tiles,
           const SuperTiles& a_finer,
           const int         a_refToFine,
           const Box&        a_thisTileBox) const noexcept;

  /**
    @brief Make the super-tile representation on the current level from tags and finer-level nesting.
    @param[out] a_tiles      Representation on the current level
    @param[in]  a_fineTiles  Representation on the finer level
    @param[in]  a_coarTags   Refinement tags on the current level (cells on the coarser level)
    @param[in]  a_domain     Problem domain on the current level
    @param[in]  a_refToFine  Refinement ratio to the finer level
    @param[in]  a_refToCoar  Refinement ratio to the coarser level
  */
  void
  makeLevelTiles(SuperTiles&          a_tiles,
                 const SuperTiles&    a_fineTiles,
                 const IntVectSet&    a_coarTags,
                 const ProblemDomain& a_domain,
                 const int            a_refToFine,
                 const int            a_refToCoar) const noexcept;

  /**
    @brief Turn the super-tile representation into boxes.
    @details Full super-tiles become one (max_block_size) box each; partial super-tiles are packed into
    variable-sized boxes by packTiles. When max_block_size == min_block_size every super-tile is a single
    tile and this emits one box per tile (legacy behaviour). Super-tiles are visited in sorted-key order
    so the box list is identical on every rank.
    @param[out] a_boxes  Grid boxes generated from the tiles
    @param[in]  a_tiles  Representation to convert
    @param[in]  a_domain Problem domain
  */
  void
  makeBoxesFromTiles(Vector<Box>& a_boxes, const SuperTiles& a_tiles, const ProblemDomain& a_domain) const noexcept;

  /**
    @brief Pack a set of tiles into variable-sized (anisotropic) boxes via recursive Berger-Rigoutsos
    bisection. This is the single shared box packer, called per partial super-tile.
    @details Recursive longest-extent bisection over the tile bounding box. At each node: if the node is
    fully tagged and within the per-direction cap it is emitted as one box; otherwise it is split. The
    split is tag-driven (Berger-Rigoutsos): cut at an empty slab (a zero-tag gap) in any direction if one
    exists, else at the strongest signature inflection, else at the longest-axis midpoint; a node over the
    cap is split at a cap multiple first. Deterministic (fixed axis/plane tie-breaks) and exact (a box is
    emitted only when fully tagged, so no untagged cells are covered).
    @param[in,out] a_tiles   Tiles to pack, in tile coordinates; reordered in place.
    @param[in]     a_maxTile Per-direction cap in tile units (= max_block_size / tile_size).
    @param[in]     a_probLo  Cell-space origin of the level (added when forming the cell boxes).
    @param[out]    a_boxes   Boxes (cell coordinates) appended for these tiles.
  */
  void
  packTiles(std::vector<IntVect>& a_tiles,
            const IntVect&        a_maxTile,
            const IntVect&        a_probLo,
            Vector<Box>&          a_boxes) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
