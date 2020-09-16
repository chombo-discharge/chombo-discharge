#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// CylindricalHarmonicGrad.cpp
// petermc, 30 Jun 2010

// #include <math.h>
#include "CylindricalHarmonicGrad.H"
#include "FunCylindricalGradF_F.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
// default constructor
CylindricalHarmonicGrad::CylindricalHarmonicGrad()
{
  //  setDefaultValues();
}


// ---------------------------------------------------------
void CylindricalHarmonicGrad::define(int          a_funN,
                                     int          a_funK)
{
  m_funN = a_funN;
  m_funK = a_funK;
  m_isDefined = true;
}


// ---------------------------------------------------------
// destructor
CylindricalHarmonicGrad::~CylindricalHarmonicGrad()
{
}


// ---------------------------------------------------------
void CylindricalHarmonicGrad::setVectorFunctionFromPhysical(FArrayBox&         a_funFab,
                                                            const Box&         a_bx,
                                                            const FArrayBox&   a_coordsFab)
{
  CH_assert(isDefined());
  FORT_SETFUNCTIONCYLINDRICALGRAD(CHF_FRA(a_funFab),
                                  CHF_CONST_FRA(a_coordsFab),
                                  CHF_BOX(a_bx),
                                  CHF_CONST_INT(m_funN),
                                  CHF_CONST_INT(m_funK));
}


// ---------------------------------------------------------
RealVect CylindricalHarmonicGrad::setVectorFunctionMax(Real   a_bxWidth,
                                                       Real   a_outerRadius)
{
  CH_assert(isDefined());
  RealVect maxfun;
  FORT_MAXFUNCTIONCYLINDRICALGRAD(CHF_REALVECT(maxfun),
                                  CHF_CONST_REAL(a_bxWidth),
                                  CHF_CONST_REAL(a_outerRadius),
                                  CHF_CONST_INT(m_funN),
                                  CHF_CONST_INT(m_funK));
  return maxfun;
}

#include "NamespaceFooter.H"
