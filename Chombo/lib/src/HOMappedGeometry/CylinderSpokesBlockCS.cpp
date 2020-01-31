#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CylinderSpokesBlockCS.H"
#include "CylinderSpokesMapsF_F.H"

#include "NamespaceHeader.H"

CylinderSpokesBlockCS::CylinderSpokesBlockCS(int a_blockNum,
                                             const RealVect& a_dx,
                                             const Box& a_blockBox)
{
  m_blockNum = a_blockNum;
  m_dx = a_dx;
  m_blockBox = a_blockBox;
}

CylinderSpokesBlockCS::~CylinderSpokesBlockCS()
{
}

RealVect
CylinderSpokesBlockCS::realCoord(const RealVect& a_Xi) const
{
  RealVect x;
  FORT_CYLINDERSPOKESREAL(CHF_REALVECT(x),
                          CHF_CONST_INT(m_blockNum),
                          CHF_CONST_REALVECT(a_Xi));
  return x;
}

RealVect
CylinderSpokesBlockCS::mappedCoord(const RealVect& a_x) const
{
  RealVect xi;
  FORT_CYLINDERSPOKESMAPPED(CHF_REALVECT(xi),
                            CHF_CONST_INT(m_blockNum),
                            CHF_CONST_REALVECT(a_x));
  return xi;
}

Real
CylinderSpokesBlockCS::dXdXi(const RealVect& a_Xi,
                             int a_dirX,
                             int a_dirXi) const
{
  Real deriv;
  FORT_CYLINDERSPOKESDERIV(CHF_REAL(deriv),
                           CHF_CONST_INT(m_blockNum),
                           CHF_CONST_REALVECT(a_Xi),
                           CHF_CONST_INT(a_dirX),
                           CHF_CONST_INT(a_dirXi));
  return deriv;
}

#include "NamespaceFooter.H"
