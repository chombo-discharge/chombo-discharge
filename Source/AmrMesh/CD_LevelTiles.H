/* chombo-discharge
 * Copyright © 2024 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_LevelTiles.H
  @brief  Declaration of LevelTiles
  @author Robert Marskar
*/

#ifndef CD_LevelTiles_H
#define CD_LevelTiles_H

// Std includes
#include <vector>
#include <map>

// Chombo includes
#include <IntVect.H>
#include <DisjointBoxLayout.H>
#include <DataIndex.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for storing the AMR hierarchy as a collection of tiles. 
  @details This class only makes sense in terms of a grid which uses constant-size grid, which is a restriction we happily accept when using particles. This class
  decomposes the AMR hierarchy into such tiles, each tile being the cells spanned by the blockingFactor^SpaceDim grid patch. 
*/
class LevelTiles
{
public:
  using BoxIDs = std::pair<unsigned int, unsigned int>;

  /*!
    @brief For lexicographically comparing tiles
    @param[in] a First tile
    @param[in] b First tile
  */
  struct TileComparator
  {
    bool
    operator()(const IntVect& a, const IntVect& b) const
    {
      D_TERM6(if (a[0] > b[0]) return true; else if (a[0] < b[0]) return false;, if (a[1] > b[1]) return true;
              else if (a[1] < b[1]) return false;
              , if (a[2] > b[2]) return true;
              else if (a[2] < b[2]) return false;
              , if (a[3] > b[3]) return true;
              else if (a[3] < b[3]) return false;
              , if (a[4] > b[4]) return true;
              else if (a[4] < b[4]) return false;
              , if (a[5] > b[5]) return true;
              else if (a[5] < b[5]) return false;);
      return false;
    }
  };

  /*!
    @brief Weak constructor. Must subsequently call the define function
  */
  LevelTiles() noexcept;

  /*!
    @brief Full constructor - puts object in usable state. 
    @param[in] a_dbl Grids
    @param[in] a_blockingFactor Grid size. 
  */
  LevelTiles(const DisjointBoxLayout& a_grids, const int a_blockingFactor) noexcept;

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~LevelTiles() noexcept;

  /*!
    @brief Define function. Puts object in usable state. 
    @param[in] a_dbl Grids
    @param[in] a_blockingFactor Grid size. 
  */
  virtual void
  define(const DisjointBoxLayout& a_grids, const int a_blockingFactor) noexcept;

  /*!
    @brief Get the tiles own by this rank
    @return m_myTiles
  */
  virtual const std::map<IntVect, unsigned int, TileComparator>&
  getMyTiles() const noexcept;

  /*!
    @brief Get the tiles own by this rank
    @return m_otherTiles
  */
  virtual const std::map<IntVect, BoxIDs, TileComparator>&
  getOtherTiles() const noexcept;

  /*!
    @brief Get the tiles own by this rank
    @return m_myGridIndices
  */
  virtual const std::map<unsigned int, DataIndex>&
  getMyGrids() const noexcept;

protected:
  /*!
    @brief Is defined or not
  */
  bool m_isDefined;

  /*!
    @brief "Tiles" owned by this rank
  */
  std::map<IntVect, unsigned int, TileComparator> m_myTiles;

  /*!
    @brief Grids owned by other ranks. 
  */
  std::map<IntVect, BoxIDs, TileComparator> m_otherTiles;

  /*!
    @brief Mapping of grid index to DataIndex
  */
  std::map<unsigned int, DataIndex> m_myGrids;
};

#include <CD_NamespaceFooter.H>

#endif
