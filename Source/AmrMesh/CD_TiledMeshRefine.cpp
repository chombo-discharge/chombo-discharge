/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TiledMeshRefine.cpp
  @brief  Implementation of CD_TiledMeshRefine.H
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_TiledMeshRefine.H>
#include <CD_NamespaceHeader.H>

TiledMeshRefine::TiledMeshRefine(const ProblemDomain& a_coarsestDomain,
                                 const Vector<int>&   a_refRatios,
                                 const IntVect&       a_tileSize)
{
  CH_TIME("TiledMeshRefine::TiledMeshRefine");

  m_refRatios = a_refRatios;
  m_tileSize  = a_tileSize;

  m_vectDomains.resize(0);

  m_vectDomains.push_back(a_coarsestDomain);
  for (int lvl = 1; lvl < a_refRatios.size(); lvl++) {
    ProblemDomain domain = m_vectDomains[lvl - 1];
    domain.refine(a_refRatios[lvl - 1]);
    m_vectDomains.push_back(domain);
  }

  this->sanityCheck();
}

TiledMeshRefine::~TiledMeshRefine() { CH_TIME("TiledMeshRefine::~TiledMeshRefine"); }

void
TiledMeshRefine::sanityCheck() const
{
  CH_TIME("TiledMeshRefine::sanityCheck");

  for (int lvl = 0; lvl < m_vectDomains.size(); lvl++) {
    for (int dir = 0; dir < SpaceDim; dir++) {
      if (m_vectDomains[lvl].domainBox().size(dir) % m_tileSize[dir] != 0) {
        MayDay::Abort("TiledMeshRefine::sanityCheck - domainBox%tileSize != 0");
      }
      for (int dir = 0; dir < SpaceDim; dir++) {
        if (m_tileSize[dir] % m_refRatios[lvl] != 0) {
          MayDay::Abort("TiledMeshRefine::sanityCheck - tileSize%refRat != 0");
        }
      }
    }
  }
}

int
TiledMeshRefine::regrid(Vector<Vector<Box>>&       a_newGrids,
                        const Vector<IntVectSet>&  a_tags,
                        const int                  a_baseLevel,
                        const int                  a_topLevel,
                        const Vector<Vector<Box>>& a_oldGrids) const
{
  CH_TIME("TiledMeshRefine::regrid");

  // CH_assert(a_topLevel >= 0 );
  // CH_assert(a_baseLevel < (a_topLevel+1) && a_baseLevel >= 0 );
  // CH_assert(a_OldGrids.size() >= a_topLevel + 1 );
  // CH_assert(m_vectDomains.size() >= a_topLevel + 1 );
  // CH_assert(m_nRefVect.size() >= a_topLevel + 1 );
  // CH_assert(a_tags.size() >= a_topLevel+1 );
  // CH_assert(a_topLevel >= 0 );
  // CH_assert(a_baseLevel < (a_topLevel+1) && a_baseLevel >= 0 );

  if (a_baseLevel > 0) {
    MayDay::Error("TiledMeshRefine::regrid - a_baseLevel>0 not yet supported");
  }

  // Set the top level to be the finest level which actually has tags. Ranks may disagree on what is the "top level".
  int myTopLevel = 0;

  for (int lvl = 0; lvl < a_tags.size(); lvl++) {
    if (!a_tags[lvl].isEmpty()) {
      myTopLevel = lvl;
    }
  }

#ifdef CH_MPI
  int topLevel = -1;
  MPI_Allreduce(&myTopLevel, &topLevel, 1, MPI_INT, MPI_MAX, Chombo_MPI::comm);
#else
  int topLevel = myTopLevel;
#endif

  int newFinestLevel = 0;
  if (topLevel >= a_baseLevel) { // We have something that can change
    newFinestLevel = 1 + topLevel;

    // Extra level of empty tiles so we can use makeLevelTiles for all levels
    Vector<IntVectSet> tiles(2 + newFinestLevel, IntVectSet());

    // Make tiles on all levels above a_baseLevel
    for (int lvl = newFinestLevel; lvl > a_baseLevel; lvl--) {
      //      pout() << "making tiles on level = " << lvl << endl;
      this->makeLevelTiles(tiles[lvl],
                           tiles[lvl + 1],
                           a_tags[lvl - 1],
                           m_vectDomains[lvl],
                           m_refRatios[lvl],
                           m_refRatios[lvl - 1]);
    }

    // Make tiles into boxes
    a_newGrids.resize(1 + newFinestLevel);
    for (int lvl = 0; lvl <= a_baseLevel; lvl++) {
      //      pout() << "copying old grids on level = " << lvl << endl;
      a_newGrids[lvl] = a_oldGrids[lvl];
    }

    for (int lvl = a_baseLevel + 1; lvl <= newFinestLevel; lvl++) {
      //      pout() << "making boxes from tiles on level = " << lvl << endl;
      this->makeBoxesFromTiles(a_newGrids[lvl], tiles[lvl], m_vectDomains[lvl]);
    }
  }
  else { // If we don't have any tags, just return the old boxes
    newFinestLevel = a_baseLevel;
    a_newGrids.resize(a_oldGrids.size());
    for (int lvl = 0; lvl <= a_baseLevel; lvl++) {
      a_newGrids[lvl] = a_oldGrids[lvl];
    }
  }

  //  MayDay::Abort("TiledMeshRefine::regrid - not implemented");
  //  pout() << "TiledMeshRefine::regrid - done" << endl;
  return newFinestLevel;
}

void
TiledMeshRefine::makeLevelTiles(IntVectSet&          a_levelTiles,
                                const IntVectSet&    a_fineLevelTiles,
                                const IntVectSet&    a_coarLevelTags,
                                const ProblemDomain& a_levelDomain,
                                const int            a_refFine,
                                const int            a_refCoar) const
{
  CH_TIME("TiledMeshRefine::makeLevelTiles");

  // Lo/Hi corners and number of tiles in each direction

  const IntVect numLevelTiles = a_levelDomain.domainBox().size() / m_tileSize;
  const Box     levelTileBox(IntVect::Zero, numLevelTiles - 1);

  // 1. Generate tiles on this level from tags on the coarser level
  Box coarBox = a_levelDomain.domainBox();
  coarBox.coarsen(a_refCoar);
  const IntVect probLo       = coarBox.smallEnd();
  IntVectSet    myLevelTiles = IntVectSet();

  for (IVSIterator ivsIt(a_coarLevelTags); ivsIt.ok(); ++ivsIt) {
    const IntVect iv = ivsIt();

    IntVect curTile;
    for (int dir = 0; dir < SpaceDim; dir++) {
      curTile[dir] = (iv[dir] - probLo[dir]) / (m_tileSize[dir] / a_refCoar);
    }
    if (levelTileBox.contains(curTile)) {
      myLevelTiles |= curTile;
    }
  }

  // 2. If the domain is periodic, add symmetry tiles
  for (int dir = 0; dir < SpaceDim; dir++) {
    if (a_levelDomain.isPeriodic(dir)) {
      IntVectSet periodicTiles;
      for (IVSIterator ivsIt(myLevelTiles); ivsIt.ok(); ++ivsIt) {
        const IntVect tile = ivsIt();

        IntVect mirrorTile = tile;
        if (tile[dir] == 0) {
          mirrorTile[dir] = numLevelTiles[dir];
          periodicTiles |= mirrorTile;
        }
        else if (tile[dir] == numLevelTiles[dir]) {
          mirrorTile[dir] = 0;
          periodicTiles |= mirrorTile;
        }
      }
      myLevelTiles |= periodicTiles;
    }
  }

  // 3. Gather tiles globally
#ifdef CH_MPI
  const int          destProc = uniqueProc(SerialTask::compute);
  Vector<IntVectSet> allTaggedTiles;
  gather(allTaggedTiles, myLevelTiles, destProc);
  if (procID() == destProc) {
    a_levelTiles.makeEmpty();
    for (int i = 0; i < allTaggedTiles.size(); i++) {
      a_levelTiles |= allTaggedTiles[i];
    }
  }
  broadcast(a_levelTiles, destProc);
#else
  a_levelTiles = myLevelTiles;
#endif

  // 4. Add finer level tiles to this level to ensure proper nesting. We do this by adding all the neighboring tiles to a tile
  //    on the finer level, coarsening all those tiles and adding them to this level
  if (a_fineLevelTiles.numPts() > 0) {
    //    MayDay::Abort("wtf");
    Box fineTileBox = levelTileBox;
    fineTileBox.refine(a_refFine);

    for (IVSIterator ivsIt(a_fineLevelTiles); ivsIt.ok(); ++ivsIt) {
      const IntVect iv = ivsIt();

      // Grow tiles and restrict to fine domain
      Box box(iv, iv);
      box.grow(1);
      IntVectSet grownFineLevelTiles(box);
      grownFineLevelTiles &= fineTileBox;

      a_levelTiles |= coarsen(grownFineLevelTiles, a_refFine);
    }
  }

  // 5. Restrict tiles to problem domains. Shouldn't be necessary
  //  a_levelTiles &= levelTileBox;
}

void
TiledMeshRefine::makeBoxesFromTiles(Vector<Box>&         a_levelBoxes,
                                    const IntVectSet&    a_levelTiles,
                                    const ProblemDomain& a_levelDomain) const
{
  CH_TIME("TiledMeshRefine::makeBoxesFromTiles");

  a_levelBoxes.resize(0);

  const IntVect probLo = a_levelDomain.domainBox().smallEnd();
  const IntVect probHi = a_levelDomain.domainBox().bigEnd();

  for (IVSIterator ivsIt(a_levelTiles); ivsIt.ok(); ++ivsIt) {
    const IntVect iv = ivsIt();

    const IntVect lo = probLo + iv * m_tileSize;
    const IntVect hi = lo + m_tileSize - IntVect::Unit;

    a_levelBoxes.push_back(Box(lo, hi));
  }
}

#include <CD_NamespaceFooter.H>
