#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// CylindricalHarmonic.cpp
// petermc, 21 Mar 2008

// #include <math.h>
#include "CylindricalHarmonic.H"
#include "FunCylindricalF_F.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
// default constructor
CylindricalHarmonic::CylindricalHarmonic()
{
  //  setDefaultValues();
}


// ---------------------------------------------------------
void CylindricalHarmonic::define(int          a_funN,
                                 int          a_funK)
{
  m_funN = a_funN;
  m_funK = a_funK;
  m_isDefined = true;
}


// ---------------------------------------------------------
// destructor
CylindricalHarmonic::~CylindricalHarmonic()
{
}


// ---------------------------------------------------------
void CylindricalHarmonic::setFunctionFromPhysical(FArrayBox&         a_funFab,
                                                  const Box&         a_bx,
                                                  const FArrayBox&   a_coordsFab)
{
  CH_assert(isDefined());
  FORT_SETFUNCTIONCYLINDRICAL(CHF_FRA(a_funFab),
                              CHF_CONST_FRA(a_coordsFab),
                              CHF_BOX(a_bx),
                              CHF_CONST_INT(m_funN),
                              CHF_CONST_INT(m_funK));
}


// ---------------------------------------------------------
Real CylindricalHarmonic::setFunctionMax(Real   a_bxWidth,
                                         Real   a_outerRadius)
{
  CH_assert(isDefined());
  Real maxfun;
  FORT_MAXFUNCTIONCYLINDRICAL(CHF_REAL(maxfun),
                              CHF_CONST_REAL(a_bxWidth),
                              CHF_CONST_REAL(a_outerRadius),
                              CHF_CONST_INT(m_funN),
                              CHF_CONST_INT(m_funK));
  return maxfun;
}

#include "NamespaceFooter.H"
