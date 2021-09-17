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

// Chombo includes
#include <IntVectSet.H>
#include <ProblemDomain.H>
#include <Box.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for generation AMR boxes using a tiling algorithm. 
  @details This class provides a scalable method for grid generation where the grids are generated in a pre-set tile pattern. This 
  work by decomposing the grid into a tiled pattern, and then flagging tiles rather than cells for refinement. 
*/
class TiledMeshRefine {
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
  TiledMeshRefine(const ProblemDomain& a_coarsestDomain, const Vector<int>& a_refRatios, const IntVect& a_tileSize);

  /*!
    @brief Destructor (does nothing)
  */
  ~TiledMeshRefine();

  /*!
    @brief Regrid using the tile clustering algorithm.
    @param[out] a_newBoxes  The new grid boxes
    @param[in]  a_tagsLevel Grid tags on each level
    @param[in]  a_baseLevel Coarsest level which does not change
    @param[in]  a_topLevel  Finest level which has tags
    @param[in]  a_oldBoxes  Previous grid boxes
    @return Returns the new finest grid level. 
  */
  int regrid(Vector<Vector<Box> >&       a_newBoxes,
	     const Vector<IntVectSet>&   a_tagsLevel,
	     const int                   a_baseLevel,
	     const int                   a_topLevel,
	     const Vector<Vector<Box> >& a_oldBoxes) const;

protected:

  /*!
    @brief Computational domains on each level
  */
  Vector<ProblemDomain> m_vectDomains;

  /*!
    @brief Refinement ratios
  */
  Vector<int> m_refRatios;

  /*!
    @brief Tile size
  */
  IntVect m_tileSize;

  /*!
    @brief Sanity check, performs some basic assertions
  */
  void sanityCheck() const ;

  /*!
    @brief Make tiles on the current level from a_levelTags and coarsenings of the finer level tiles
    @param[out] a_levelTiles
    @param[in]  a_fineLevelTiles Finer grid level tiles
    @param[in]  a_coarLevel      Coarser grid level tiles
    @param[in]  a_levelDomain    Computational domain for current level
    @param[in]  a_refFine        Refinement factor to fine level
    @param[in]  a_refCoar        Refinement factor to coarse level
  */
  void makeLevelTiles(IntVectSet&          a_levelTiles,
		      const IntVectSet&    a_fineLevelTiles,
		      const IntVectSet&    a_coarLevelTags,
		      const ProblemDomain& a_levelDomain,
		      const int            a_refFine,
		      const int            a_refCoar) const;

  /*!
    @brief Given tiles, generate the corresponding boxes
    @param[out] a_levelBoxes  Boxes corresponding to tiles
    @param[in]  a_levelTiles  Tiles on level
    @param[in]  a_leveLDomain Computational domain on level
  */
  void makeBoxesFromTiles(Vector<Box>& a_levelBoxes, const IntVectSet& a_levelTiles, const ProblemDomain& a_levelDomain) const;

};

#include <CD_NamespaceFooter.H>

#endif
