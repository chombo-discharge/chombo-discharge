#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


/******************************************************************************/
/**
 * \file LevelGridMetrics.cpp
 *
 * \brief Non-inline definitions for multiblock in LevelMappedConsOperator.cpp
 *
 * Member functions that are only used for multiblock calculations are defined
 * here.
 *
 *//*+*************************************************************************/

#include "LevelMappedConsOperator.H"

#include "BlockRegister.H"


/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Define members required for multiblock grids
/** Such as m_facesToFill which is required for use of BlockRegister
 *//*-----------------------------------------------------------------*/

void
LevelMappedConsOperator::defineMultiblockMbrs()
{

//--Define m_facesToFill

  const MultiBlockCoordSys& coordSys = m_levelGridMetricsPtr->getCoordSys();
  const DisjointBoxLayout& grids = m_levelGridMetricsPtr->getBoxes();
  const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries =
    coordSys.boundaries();
  // Tuple< LayoutData<IntVectSet>, 2*SpaceDim> m_facesToFill;
  int faceID = 0;
  for (SideIterator sit; sit.ok(); ++sit)
    {
      Side::LoHiSide side = sit();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_facesToFill[faceID].define(grids);
          for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
            {
              const Box& baseBox = grids[dit];
              int baseBlockNum = coordSys.whichBlock(baseBox);
              const Box& blockBox = coordSys.mappingBlocks()[baseBlockNum];
              Box blockBdryBox = adjCellBox( blockBox, idir, side, 1 );
              Box clipBox = adjCellBox( baseBox, idir, side, 1 );
              clipBox &= blockBdryBox;
              if ( clipBox.isEmpty() )
                { // empty
                  m_facesToFill[faceID][dit] = IntVectSet();
                }
              else
                {
                  const BlockBoundary& bb = boundaries[baseBlockNum][faceID];
                  const IndicesTransformation& it = bb.getTransformation();
                  Box clipOtherBox = it.transform(clipBox);
                  // ivsMissing will hold all IntVects of clipBox
                  // that are NOT present in the other block.
                  IntVectSet ivsMissing(clipBox);
                  // Loop through all boxes, and remove their back transforms
                  // from ivsMissing.
                  LayoutIterator lit = grids.layoutIterator();
                  for (lit.begin(); lit.ok(); ++lit)
                    {
                      const Box& baseOtherBox = grids[lit];
                      if (baseOtherBox.intersectsNotEmpty(clipOtherBox))
                        {
                          Box baseThisBox = it.transformBack(baseOtherBox);
                          ivsMissing -= baseThisBox;
                        }
                    }
                  // Got ivsMissing.
                  // ivsPresent will hold all IntVects of clipBox
                  // that ARE present in the other block.
                  IntVectSet ivsPresent(clipBox);
                  ivsPresent -= ivsMissing;
                  m_facesToFill[faceID][dit] = IntVectSet(ivsPresent);
                }
            }
          faceID++;
        }
    }
}

/*--------------------------------------------------------------------*/
//  Sets single-valued fluxes at block boundaries
/** Using the contents of a BlockRegister, set the single-valued
 *  fluxes at block boundaries.
 *//*-----------------------------------------------------------------*/

void
LevelMappedConsOperator::setCommonFlux(
  LevelData<FluxBox>&   a_flux,
  const BlockRegister&  a_blockRegister) const
{
  CH_TIME("LevelMappedConsOperator::setCommonFlux");
  const MultiBlockCoordSys& coordSys = m_levelGridMetricsPtr->getCoordSys();
  const DisjointBoxLayout& grids = m_levelGridMetricsPtr->getBoxes();
  const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries =
    coordSys.boundaries();
  for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
  {
    FluxBox& FfaceAvg = a_flux[dit];
    const Box& baseBox = grids[dit];
    int baseBlockNum = coordSys.whichBlock(baseBox);
    int faceID = 0;
    for (SideIterator sit; sit.ok(); ++sit)
    {
      Side::LoHiSide side = sit();
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        if (a_blockRegister.hasInterface(dit(), idir, side))
        {
          // maybe better if this is done inside BlockRegister
          const BlockBoundary& bb = boundaries[baseBlockNum][faceID];
          int reorientFace = bb.reorientFace(idir);
          Box faceBox = adjCellBox(baseBox, idir, side, 1);
          // faceBox.grow(faceGrowVect);
          // if Lo, then shift +1; if Hi, then shift -1
          faceBox.shiftHalf(idir, -sign(side));
          Side::LoHiSide sideOther = flip(side);
          // Need to define these FABs.
          FArrayBox fluxThisFab(faceBox, m_numFluxes); // from this block
          FArrayBox fluxOtherFab(faceBox, m_numFluxes); // from other block
          a_blockRegister.getFlux(fluxThisFab, dit(),
              idir, side, sideOther);
          a_blockRegister.getFlux(fluxOtherFab, dit(),
              idir, side, side);
          fluxOtherFab.mult(reorientFace * 0.5);
          fluxThisFab.mult(0.5);
          fluxThisFab += fluxOtherFab;
          // FfaceAvg[idir].copy(fluxThisFab);
          FArrayBox& faceFillFab = FfaceAvg[idir];
          int ncomp = faceFillFab.nComp();
          const IntVectSet& ivs = m_facesToFill[faceID][dit];
          for (IVSIterator ivsit(ivs); ivsit.ok(); ++ivsit)
            {
              // ivClip is index of the ghost cell just outside (side, idir)
              // of the block.
              // (Store cells, not faces, because IntVectSet is for cells.)
              IntVect ivClip = ivsit();
              // ivFace is index of the cell face on the boundary,
              // or ivClip shifted by "half" in (-side, idir) direction.
              // If side is high, then ivFace is ivClip shifted down from
              // cell to face in direction idir.  So no change in coordinates.
              // If side is low, then ivFace is ivClip shifted up from
              // cell to face in direction idir.  So +1 in idir coordinate.
              IntVect ivFace = ivClip;
              if (side == Side::Lo) ivFace += BASISV(idir);
              for (int icomp = 0; icomp < ncomp; icomp++)
                faceFillFab(ivFace, icomp) = fluxThisFab(ivFace, icomp);
            }
        }
        faceID++;
      } // iterate over dimensions
    } // iterate over sides
  } // iterate over patches
}
