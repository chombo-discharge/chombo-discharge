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

#include "RectangularEqualAnglesMap.H"
#include "RectangularEqualAnglesMappingF_F.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
RectangularEqualAnglesMap::RectangularEqualAnglesMap(/// physical coordinates of low corner
                                                     const RealVect&   a_cornerLo,
                                                     /// physical coordinates of high corner
                                                     const RealVect&   a_cornerHi,
                                                     /// CELL-centered domain in index space
                                                     const Box&        a_bxCells)
  : RectangularMap(a_cornerLo, a_cornerHi, a_bxCells)
{
}


// ---------------------------------------------------------
void RectangularEqualAnglesMap::setPhysicalFromMap(/// physical coordinates, SpaceDim components
                                                   FArrayBox&         a_physFab,
                                                   /// box on which to set physical coordinates
                                                   const Box&         a_bx,
                                                   /// mapped coordinates, SpaceDim components
                                                   const FArrayBox&   a_mapFab) const
{
  RealVect cartesianOriginRefined = m_refinement * m_cartesianOrigin;
  RealVect cartesianLengthRefined = m_refinement * m_cartesianLength;
  FORT_MAPRECTANGULAREQUALANGLES(CHF_FRA(a_physFab),
                                 CHF_CONST_FRA(a_mapFab),
                                 CHF_CONST_REALVECT(cartesianOriginRefined),
                                 CHF_CONST_REALVECT(cartesianLengthRefined),
                                 CHF_CONST_REALVECT(m_physCornerLo),
                                 CHF_CONST_REALVECT(m_bxWidth),
                                 CHF_BOX(a_bx));
}


// ---------------------------------------------------------
void RectangularEqualAnglesMap::setMapFromPhysical(/// mapped coordinates, SpaceDim components
                                                   FArrayBox&         a_mapFab,
                                                   /// box on which to set physical coordinates
                                                   const Box&         a_bx,
                                                   /// physical coordinates, SpaceDim components
                                                   const FArrayBox&   a_physFab) const
{
  RealVect cartesianOriginRefined = m_refinement * m_cartesianOrigin;
  RealVect cartesianLengthRefined = m_refinement * m_cartesianLength;
  FORT_MAPINVRECTANGULAREQUALANGLES(CHF_FRA(a_mapFab),
                                    CHF_CONST_FRA(a_physFab),
                                    CHF_CONST_REALVECT(cartesianOriginRefined),
                                    CHF_CONST_REALVECT(cartesianLengthRefined),
                                    CHF_CONST_REALVECT(m_physCornerLo),
                                    CHF_CONST_REALVECT(m_bxWidth),
                                    CHF_BOX(a_bx));
}


// ---------------------------------------------------------
void RectangularEqualAnglesMap::setJacobian(/// jacobian on a_bx
                                            FArrayBox&         a_jacobianFab,
                                            /// box on which to set jacobian
                                            const Box&         a_bx,
                                            /// mapped coordinates, SpaceDim components
                                            const FArrayBox&   a_mapFab) const
{
  RealVect cartesianOriginRefined = m_refinement * m_cartesianOrigin;
  RealVect cartesianLengthRefined = m_refinement * m_cartesianLength;
  FORT_JACOBIANRECTANGULAREQUALANGLES(CHF_FRA1(a_jacobianFab, 0),
                                      CHF_CONST_FRA(a_mapFab),
                                      CHF_CONST_REALVECT(cartesianOriginRefined),
                                      CHF_CONST_REALVECT(cartesianLengthRefined),
                                      CHF_CONST_REALVECT(m_bxWidth),
                                      CHF_BOX(a_bx));
}

#include "NamespaceFooter.H"
