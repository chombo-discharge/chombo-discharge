#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CylinderEquiangularBlockCS.H"
#include "CylinderEquiangularMapsF_F.H"

#include "NamespaceHeader.H"

CylinderEquiangularBlockCS::CylinderEquiangularBlockCS(int a_blockNum,
                                                       const RealVect& a_dx,
                                                       const Box& a_blockBox)
{
  m_blockNum = a_blockNum;
  m_dx = a_dx;
  m_blockBox = a_blockBox;
}

CylinderEquiangularBlockCS::~CylinderEquiangularBlockCS()
{
}

RealVect
CylinderEquiangularBlockCS::realCoord(const RealVect& a_Xi) const
{
  RealVect x;
  FORT_CYLINDEREQUIANGULARREAL(CHF_REALVECT(x),
                               CHF_CONST_INT(m_blockNum),
                               CHF_CONST_REALVECT(a_Xi));
  return x;
}

RealVect
CylinderEquiangularBlockCS::mappedCoord(const RealVect& a_x) const
{
  RealVect xi;
  FORT_CYLINDEREQUIANGULARMAPPED(CHF_REALVECT(xi),
                                 CHF_CONST_INT(m_blockNum),
                                 CHF_CONST_REALVECT(a_x));
  return xi;
}

Real
CylinderEquiangularBlockCS::dXdXi(const RealVect& a_Xi,
                                  int a_dirX,
                                  int a_dirXi) const
{
  Real deriv;
  FORT_CYLINDEREQUIANGULARDERIV(CHF_REAL(deriv),
                                CHF_CONST_INT(m_blockNum),
                                CHF_CONST_REALVECT(a_Xi),
                                CHF_CONST_INT(a_dirX),
                                CHF_CONST_INT(a_dirXi));
  return deriv;
}

#include "NamespaceFooter.H"
