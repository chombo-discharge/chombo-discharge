/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TiledMeshRefine.H
  @brief  Header for class TiledMeshRefine
  @author Robert Marskar
*/

#ifndef CD_TiledMeshRefine_H
#define CD_TiledMeshRefine_H

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

/*!
  @brief Class for generation AMR boxes using a tiling algorithm. 
  @details This class provides a scalable method for grid generation where the grids are generated in a pre-set tile pattern. This 
  work by decomposing the grid into a tiled pattern, and then flagging tiles rather than cells for refinement. 
*/
class TiledMeshRefine
{
public:
  /*!
    @brief Default constructor -- deleted because we never need it
  */
  TiledMeshRefine() = delete;

  /*!
    @brief Constructor
    @param[in] a_coarsestDomain Coarsest grid domain
    @param[in] a_refRatios      Refinement ratios
    @param[in] a_tileSize       Tile size
  */
  TiledMeshRefine(const ProblemDomain& a_coarsestDomain,
                  const Vector<int>&   a_refRatios,
                  const IntVect&       a_tileSize) noexcept;

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~TiledMeshRefine() noexcept;

  /*!
    @brief Regrid using the tile clustering algorithm.
    @param[out] a_newBoxes  The new grid boxes
    @param[in]  a_tagsLevel Grid tags on each level
    @return Returns the new finest grid level. 
  */
  virtual int
  regrid(Vector<Vector<Box>>& a_newBoxes, const Vector<IntVectSet>& a_tagsLevel) const noexcept;

protected:
  /*!
    @brief Tile representation
  */
  using Tile = TileI<int, SpaceDim>;

  /*!
    @brief Tile set representation
  */
  using TileSet = std::set<Tile>;

  /*!
    @brief Computational domains on each level
    @note This are the domains for the grid.
  */
  Vector<ProblemDomain> m_amrDomains;

  /*!
    @brief Refinement ratios
  */
  Vector<int> m_refRatios;

  /*!
    @brief Tile size
  */
  IntVect m_tileSize;

  /*!
    @brief Make tiles on the current level from tags and tile coarsening from finer levels
  */
  virtual void
  makeLevelTiles(TileSet&             a_tiles,
                 const TileSet&       a_fineTiles,
                 const IntVectSet&    a_coarTags,
                 const ProblemDomain& a_domain,
                 const int            a_refToFine,
                 const int            a_refToCoar) const noexcept;

  /*!
    @brief Turn tiles into boxes
  */
  virtual void
  makeBoxesFromTiles(Vector<Box>& a_boxes, const TileSet& a_tiles, const ProblemDomain& a_domain) const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
