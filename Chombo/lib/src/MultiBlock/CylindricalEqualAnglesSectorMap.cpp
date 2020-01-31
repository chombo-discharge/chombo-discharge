#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// petermc, 31 Mar 2008

#include "CylindricalEqualAnglesSectorMap.H"
#include "CylindricalEqualAnglesSectorMappingF_F.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
CylindricalEqualAnglesSectorMap::CylindricalEqualAnglesSectorMap(/// physical coordinates of center of cylinder
                                                                 const RealVect&   a_physCenter,
                                                                 /// length of central cube in each dimension
                                                                 Real              a_bxWidth,
                                                                 /// outer radius of cylinder
                                                                 Real              a_outerRadius,
                                                                 /// which block:  1, 2, 3, 4 for +x, +y, -x, -y
                                                                 int               a_block,
                                                                 /// CELL-centered domain of sector in index space
                                                                 const Box&        a_bxCells)
  : CylindricalSectorMap(a_physCenter, a_bxWidth, a_outerRadius, a_block, a_bxCells)
{
}


// ---------------------------------------------------------
void CylindricalEqualAnglesSectorMap::setPhysicalFromMap(/// physical coordinates, SpaceDim components
                                                         FArrayBox&         a_physFab,
                                                         /// box on which to set physical coordinates
                                                         const Box&         a_bx,
                                                         /// mapped coordinates, SpaceDim components
                                                         const FArrayBox&   a_mapFab) const
{
  RealVect cartesianOriginRefined = m_refinement * m_cartesianOrigin;
  RealVect cartesianLengthRefined = m_refinement * m_cartesianLength;
  FORT_MAPCYLINDRICALEQUALANGLESSECTOR(CHF_FRA(a_physFab),
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
void CylindricalEqualAnglesSectorMap::setMapFromPhysical(/// mapped coordinates, SpaceDim components
                                                         FArrayBox&         a_mapFab,
                                                         /// box on which to set physical coordinates
                                                         const Box&         a_bx,
                                                         /// physical coordinates, SpaceDim components
                                                         const FArrayBox&   a_physFab) const
{
  RealVect cartesianOriginRefined = m_refinement * m_cartesianOrigin;
  RealVect cartesianLengthRefined = m_refinement * m_cartesianLength;
  FORT_MAPINVCYLINDRICALEQUALANGLESSECTOR(CHF_FRA(a_mapFab),
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
void CylindricalEqualAnglesSectorMap::setJacobian(/// jacobian on a_bx
                                                  FArrayBox&         a_jacobianFab,
                                                  /// box on which to set jacobian
                                                  const Box&         a_bx,
                                                  /// mapped coordinates, SpaceDim components
                                                  const FArrayBox&   a_mapFab) const
{
  RealVect cartesianOriginRefined = m_refinement * m_cartesianOrigin;
  RealVect cartesianLengthRefined = m_refinement * m_cartesianLength;
  FORT_JACOBIANCYLINDRICALEQUALANGLESSECTOR(CHF_FRA1(a_jacobianFab, 0),
                                            CHF_CONST_FRA(a_mapFab),
                                            CHF_CONST_INT(m_block),
                                            CHF_CONST_REALVECT(cartesianOriginRefined),
                                            CHF_CONST_REALVECT(cartesianLengthRefined),
                                            CHF_CONST_REAL(m_bxWidth),
                                            CHF_CONST_REAL(m_outerRadius),
                                            CHF_BOX(a_bx));
}


#include "NamespaceFooter.H"
