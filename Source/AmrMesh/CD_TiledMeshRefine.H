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
#include <vector>
#include <set>

// Chombo includes
#include <IntVectSet.H>
#include <ProblemDomain.H>
#include <Box.H>

// Our includes
#include <CD_Tile.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Class for generation AMR boxes using a tiling algorithm. 
  @details This class provides a scalable method for grid generation where the grids are generated in a pre-set tile pattern. This 
  work by decomposing the grid into a tiled pattern, and then flagging tiles rather than cells for refinement. 
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
    @param[in] a_maxBoxSize     Maximum box size (in cells). Tiles are merged into boxes no larger than
                                this per direction. Must be a multiple of a_tileSize. When equal to
                                a_tileSize the output is one box per tile (legacy behaviour).
  */
  TiledMeshRefine(const ProblemDomain& a_coarsestDomain,
                  const Vector<int>&   a_refRatios,
                  const IntVect&       a_tileSize,
                  const IntVect&       a_maxBoxSize) noexcept;

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
    @brief Tile representation
  */
  using Tile = TileI<int, SpaceDim>;

  /**
    @brief Tile set representation
  */
  using TileSet = std::set<Tile>;

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
  IntVect m_maxBoxSize;

  /**
    @brief Make tiles on the current level from tags and tile coarsening from finer levels
    @param[out] a_tiles      Tiles on the current level
    @param[in]  a_fineTiles  Tiles on the finer level
    @param[in]  a_coarTags   Refinement tags on the current level
    @param[in]  a_domain     Problem domain on the current level
    @param[in]  a_refToFine  Refinement ratio to the finer level
    @param[in]  a_refToCoar  Refinement ratio to the coarser level
  */
  virtual void
  makeLevelTiles(TileSet&             a_tiles,
                 const TileSet&       a_fineTiles,
                 const IntVectSet&    a_coarTags,
                 const ProblemDomain& a_domain,
                 const int            a_refToFine,
                 const int            a_refToCoar) const noexcept;

  /**
    @brief Turn tiles into boxes
    @details When m_maxBoxSize == m_tileSize this emits one box per tile (legacy behaviour). Otherwise
    it merges adjacent tiles into fewer, larger, possibly anisotropic boxes via mergeTiles.
    @param[out] a_boxes  Grid boxes generated from the tiles
    @param[in]  a_tiles  Tiles to convert
    @param[in]  a_domain Problem domain
  */
  virtual void
  makeBoxesFromTiles(Vector<Box>& a_boxes, const TileSet& a_tiles, const ProblemDomain& a_domain) const noexcept;

  /**
    @brief Merge a set of tiles into variable-sized (possibly anisotropic) boxes.
    @details KD / Berger-Rigoutsos recursive longest-axis bisection over the tile grid. Each emitted box
    is a union of aligned tiles, fully tagged, and no larger than m_maxBoxSize per direction. The result
    is a pure, deterministic function of the (globally identical) tile set and a fixed longest-axis/cut
    rule, so every rank computes the same boxes. Covers exactly the same cells as the one-box-per-tile
    output, with far fewer boxes.
    @param[out] a_boxes  Grid boxes generated from the tiles
    @param[in]  a_tiles  Tiles to merge
    @param[in]  a_domain Problem domain (supplies the cell-space origin)
  */
  virtual void
  mergeTiles(Vector<Box>& a_boxes, const TileSet& a_tiles, const ProblemDomain& a_domain) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
