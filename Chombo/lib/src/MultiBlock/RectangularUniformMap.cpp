#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// petermc, 27 Mar 2008

#include "RectangularUniformMap.H"
#include "RectangularUniformMappingF_F.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
RectangularUniformMap::RectangularUniformMap(/// physical coordinates of low corner
                                             const RealVect&   a_physCornerLo,
                                             /// physical coordinates of high corner
                                             const RealVect&   a_physCornerHi,
                                             /// CELL-centered domain in index space
                                             const Box&        a_bxCells) : RectangularMap(a_physCornerLo, a_physCornerHi, a_bxCells)
{
  m_cellWidth = m_bxWidth / m_cartesianLength;
}


// ---------------------------------------------------------
void RectangularUniformMap::setPhysicalFromMap(/// physical coordinates, SpaceDim components
                                               FArrayBox&         a_physFab,
                                               /// box on which to set physical coordinates
                                               const Box&         a_bx,
                                               /// mapped coordinates, SpaceDim components
                                               const FArrayBox&   a_mapFab) const
{
  RealVect cartesianOriginRefined = m_refinement * m_cartesianOrigin;
  RealVect cellWidthRefined = m_cellWidth / (m_refinement * 1.);
  FORT_MAPRECTANGULARUNIFORM(CHF_FRA(a_physFab),
                             CHF_CONST_FRA(a_mapFab),
                             CHF_CONST_REALVECT(cartesianOriginRefined),
                             CHF_CONST_REALVECT(m_physCornerLo),
                             CHF_CONST_REALVECT(cellWidthRefined),
                             CHF_BOX(a_bx));
}


// ---------------------------------------------------------
void RectangularUniformMap::setMapFromPhysical(/// mapped coordinates, SpaceDim components
                                               FArrayBox&         a_mapFab,
                                               /// box on which to set physical coordinates
                                               const Box&         a_bx,
                                               /// physical coordinates, SpaceDim components
                                               const FArrayBox&   a_physFab) const
{
  RealVect cartesianOriginRefined = m_refinement * m_cartesianOrigin;
  RealVect cellWidthRefined = m_cellWidth / (m_refinement * 1.);
  FORT_MAPINVRECTANGULARUNIFORM(CHF_FRA(a_mapFab),
                                CHF_CONST_FRA(a_physFab),
                                CHF_CONST_REALVECT(cartesianOriginRefined),
                                CHF_CONST_REALVECT(m_physCornerLo),
                                CHF_CONST_REALVECT(cellWidthRefined),
                                CHF_BOX(a_bx));
}


// ---------------------------------------------------------
void RectangularUniformMap::setJacobian(/// jacobian on a_bx
                                        FArrayBox&         a_jacobianFab,
                                        /// box on which to set jacobian
                                        const Box&         a_bx,
                                        /// mapped coordinates, SpaceDim components
                                        const FArrayBox&   a_mapFab) const
{
  RealVect bxWidthRefined = m_bxWidth / (m_refinement * 1.);
  // I am not sure that this is correct with refinement.
  FORT_JACOBIANRECTANGULARUNIFORM(CHF_FRA1(a_jacobianFab, 0),
                                  CHF_CONST_REALVECT(bxWidthRefined),
                                  CHF_BOX(a_bx));
}

#include "NamespaceFooter.H"
