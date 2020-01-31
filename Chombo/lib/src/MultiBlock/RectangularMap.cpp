#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// petermc, 28 Mar 2008

#include "RectangularMap.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
RectangularMap::RectangularMap(/// physical coordinates of low corner
                               const RealVect&   a_physCornerLo,
                               /// physical coordinates of high corner
                               const RealVect&   a_physCornerHi,
                               /// CELL-centered domain in index space
                               const Box&        a_bxCells)
{
  m_physCornerLo = a_physCornerLo;
  m_physCornerHi = a_physCornerHi;
  m_bxWidth = m_physCornerHi - m_physCornerLo;
  m_bxCells = a_bxCells;
  m_cartesianOrigin = m_bxCells.smallEnd();
  m_cartesianLength = m_bxCells.size();
}


// ---------------------------------------------------------
bool RectangularMap::inThisBlock(const RealVect&   a_physPoint) const
{
  bool inRectangle = true;
  // Use "<=" and ">=" so that you don't get stuck "between" blocks.
  for (int idir = 0; idir <= 1; idir++) // Check only X and Y dimensions.
    {
      if (! ( (a_physPoint[idir] >= m_physCornerLo[idir]) &&
              (a_physPoint[idir] <= m_physCornerHi[idir]) ) )
        inRectangle = false;
    }
  return inRectangle;
}

#include "NamespaceFooter.H"
