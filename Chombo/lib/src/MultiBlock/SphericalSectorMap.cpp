#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// petermc, 21 May 2008

#include "SphericalSectorMap.H"
#include "SphericalSectorMappingF_F.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
SphericalSectorMap::SphericalSectorMap(/// physical coordinates of center of cylinder
                                       const RealVect&   a_physCenter,
                                       /// length of central cube in each dimension
                                       Real              a_bxWidth,
                                       /// outer radius of cylinder
                                       Real              a_outerRadius,
                                       /// which block:  0, 1, 2, 3, 4, 5 for +x, +y, -x, -y, +z, -z
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
bool SphericalSectorMap::inThisBlock(const RealVect&   a_physPoint) const
{
  bool inSector = false;
#if CH_SPACEDIM == 3
  Real halfWidth = m_bxWidth / 2.;
  RealVect physOffset = a_physPoint - m_physCenter;
  Real X = physOffset[0];
  Real Y = physOffset[1];
  Real Z = physOffset[2];
  switch (m_block)
    { // Use "<=" and ">=" so that you don't get stuck "between" blocks.
    case 0:
      { // +X
        inSector = ((X >= halfWidth) && (X >= abs(Y)) && (X >= abs(Z)));
        break;
      }
    case 1:
      { // +Y
        inSector = ((Y >= halfWidth) && (Y >= abs(X)) && (Y >= abs(Z)));
        break;
      }
    case 2:
      { // -X
        inSector = ((X <= -halfWidth) && (X <= -abs(Y)) && (X <= -abs(Z)));
        break;
      }
    case 3:
      { // -Y
        inSector = ((Y <= -halfWidth) && (Y <= -abs(X)) && (Y <= -abs(Z)));
        break;
      }
    case 4:
      { // +Z
        inSector = ((Z >= halfWidth) && (Z >= abs(X)) && (Z >= abs(Y)));
        break;
      }
    case 5:
      { // -Z
        inSector = ((Z <= -halfWidth) && (Z <= -abs(X)) && (Z <= -abs(Y)));
        break;
      }
    default:
      {
        inSector = false;
      }
    }
#endif
  return inSector;
}


// ---------------------------------------------------------
void SphericalSectorMap::setPhysicalFromMap(/// physical coordinates, SpaceDim components
                                            FArrayBox&         a_physFab,
                                            /// box on which to set physical coordinates
                                            const Box&         a_bx,
                                            /// mapped coordinates, SpaceDim components
                                            const FArrayBox&   a_mapFab) const
{
  RealVect cartesianOriginRefined = m_refinement * m_cartesianOrigin;
  RealVect cartesianLengthRefined = m_refinement * m_cartesianLength;
  FORT_MAPSPHERICALSECTOR(CHF_FRA(a_physFab),
                          CHF_CONST_FRA(a_mapFab),
                          CHF_CONST_INT(m_block),
                          CHF_CONST_REALVECT(cartesianOriginRefined),
                          CHF_CONST_REALVECT(cartesianLengthRefined),
                          CHF_CONST_REALVECT(m_physCenter),
                          CHF_CONST_REAL(m_bxWidth),
                          CHF_CONST_REAL(m_outerRadius),
                          CHF_BOX(a_bx));
}


// ---------------------------------------------------------
void SphericalSectorMap::setMapFromPhysical(/// mapped coordinates, SpaceDim components
                                            FArrayBox&         a_mapFab,
                                            /// box on which to set physical coordinates
                                            const Box&         a_bx,
                                            /// physical coordinates, SpaceDim components
                                            const FArrayBox&   a_physFab) const
{
  RealVect cartesianOriginRefined = m_refinement * m_cartesianOrigin;
  RealVect cartesianLengthRefined = m_refinement * m_cartesianLength;
  FORT_MAPINVSPHERICALSECTOR(CHF_FRA(a_mapFab),
                             CHF_CONST_FRA(a_physFab),
                             CHF_CONST_INT(m_block),
                             CHF_CONST_REALVECT(cartesianOriginRefined),
                             CHF_CONST_REALVECT(cartesianLengthRefined),
                             CHF_CONST_REALVECT(m_physCenter),
                             CHF_CONST_REAL(m_bxWidth),
                             CHF_CONST_REAL(m_outerRadius),
                             CHF_BOX(a_bx));
}


// ---------------------------------------------------------
void SphericalSectorMap::setJacobian(/// jacobian on a_bx
                                     FArrayBox&         a_jacobianFab,
                                     /// box on which to set jacobian
                                     const Box&         a_bx,
                                     /// mapped coordinates, SpaceDim components
                                     const FArrayBox&   a_mapFab) const
{
  RealVect cartesianOriginRefined = m_refinement * m_cartesianOrigin;
  RealVect cartesianLengthRefined = m_refinement * m_cartesianLength;
  FORT_JACOBIANSPHERICALSECTOR(CHF_FRA1(a_jacobianFab, 0),
                               CHF_CONST_FRA(a_mapFab),
                               CHF_CONST_INT(m_block),
                               CHF_CONST_REALVECT(cartesianOriginRefined),
                               CHF_CONST_REALVECT(cartesianLengthRefined),
                               CHF_CONST_REAL(m_bxWidth),
                               CHF_CONST_REAL(m_outerRadius),
                               CHF_BOX(a_bx));
}

#include "NamespaceFooter.H"
