/*!
  @file   TiledMeshRefine.H
  @brief  Header for TiledMeshRefine
  @author Robert Marskar
*/

#ifndef _TILEDMESHREFINE_H_
#define _TILEDMESHREFINE_H_

#include <IntVectSet.H>
#include <ProblemDomain.H>
#include <Box.H>

namespace ChomboDischarge {

  /*!
    @brief Class for generation AMR boxes using tiling
  */
  class TiledMeshRefine {
  public:

    TiledMeshRefine();
    TiledMeshRefine(const ProblemDomain& a_coarsestDomain, const Vector<int> a_refRatios, const IntVect& a_tileSize);
  
    ~TiledMeshRefine();

    /*!
      @brief Regrid using the tile clustering algorithm.
    */
    int regrid(Vector<Vector<Box> >&       a_newBoxes,
	       const Vector<IntVectSet>&   a_tagsLevel,
	       const int                   a_baseLevel,
	       const int                   a_topLevel,
	       const Vector<Vector<Box> >& a_oldBoxes);

  protected:

    Vector<ProblemDomain> m_vectDomains;  // Domains
    Vector<Box>           m_tileDomains;  // Tile domains
    Vector<int>           m_refRatios;    // Refinement ratios
    IntVect               m_tileSize;     // Tile size

    void sanityCheck();

    /*!
      @brief Make tiles on the current level from a_levelTags and coarsenings of the finer level tiles
    */
    void makeLevelTiles(IntVectSet&          a_levelTiles,
			const IntVectSet&    a_fineLevelTiles,
			const IntVectSet&    a_coarLevelTags,
			const ProblemDomain& a_levelDomain,
			const int            a_refFine,
			const int            a_refCoar);
  
    void makeBoxesFromTiles(Vector<Box>& a_levelBoxes, const IntVectSet& a_levelTiles, const ProblemDomain& a_levelDomain);

  };
}

#endif
