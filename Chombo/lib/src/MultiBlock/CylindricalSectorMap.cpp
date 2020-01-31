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

#include "CylindricalSectorMap.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
CylindricalSectorMap::CylindricalSectorMap(/// physical coordinates of center of cylinder
                                           const RealVect&   a_physCenter,
                                           /// length of central cube in each dimension
                                           Real              a_bxWidth,
                                           /// outer radius of cylinder
                                           Real              a_outerRadius,
                                           /// which block:  1, 2, 3, 4 for +x, +y, -x, -y
                                           int               a_block,
                                           /// CELL-centered domain of sector in index space
                                           const Box&        a_bxCells)
{
  m_physCenter = a_physCenter;
  m_bxWidth = a_bxWidth;
  m_outerRadius = a_outerRadius;
  m_block = a_block;
  m_cartesianOrigin = a_bxCells.smallEnd();
  m_cartesianLength = a_bxCells.size();
}


// ---------------------------------------------------------
bool CylindricalSectorMap::inThisBlock(const RealVect&   a_physPoint) const
{
  Real halfWidth = m_bxWidth / 2.;
  RealVect physOffset = a_physPoint - m_physCenter;
  Real X = physOffset[0];
  Real Y = physOffset[1];
  bool inSector;
  switch (m_block)
    { // Use "<=" and ">=" so that you don't get stuck "between" blocks.
    case 1:
      { // +X
        inSector = ((X >= halfWidth) && (X >= abs(Y)));
        break;
      }
    case 2:
      { // +Y
        inSector = ((Y >= halfWidth) && (Y >= abs(X)));
        break;
      }
    case 3:
      { // -X
        inSector = ((X <= -halfWidth) && (X <= -abs(Y)));
        break;
      }
    case 4:
      { // -Y
        inSector = ((Y <= -halfWidth) && (Y <= -abs(X)));
        break;
      }
    default:
      {
        inSector = false;
      }
    }
  return inSector;
}

#include "NamespaceFooter.H"
