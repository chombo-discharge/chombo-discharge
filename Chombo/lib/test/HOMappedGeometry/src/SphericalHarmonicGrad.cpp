#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// SphericalHarmonicGrad.cpp
// petermc, 30 Jul 2010

// #include <math.h>
#include "SphericalHarmonicGrad.H"
#include "FunSphericalGradF_F.H"
#include "LegendreF.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
// default constructor
SphericalHarmonicGrad::SphericalHarmonicGrad()
{
  //  setDefaultValues();
}


// ---------------------------------------------------------
void SphericalHarmonicGrad::define(int          a_funL,
                                   int          a_funM)
{
  m_funL = a_funL;
  m_funM = a_funM;
  m_legendreLen = m_funL - a_funM + 1;
  m_legendreCoeffs = new Real[m_legendreLen];
  FORT_LEGENDRECOEFFS(m_legendreCoeffs, &m_funL, &m_funM);
  m_isDefined = true;
}


// ---------------------------------------------------------
// destructor
SphericalHarmonicGrad::~SphericalHarmonicGrad()
{
  if (isDefined())
    {
      delete [] m_legendreCoeffs;
    }
}


// ---------------------------------------------------------
void SphericalHarmonicGrad::setVectorFunctionFromPhysical(FArrayBox&         a_funFab,
                                                            const Box&         a_bx,
                                                            const FArrayBox&   a_coordsFab)
{
  CH_assert(isDefined());
  FORT_SETFUNCTIONSPHERICALGRAD(CHF_FRA(a_funFab),
                                CHF_CONST_FRA(a_coordsFab),
                                CHF_BOX(a_bx),
                                CHF_CONST_R1D(m_legendreCoeffs, m_legendreLen),
                                CHF_CONST_INT(m_funL),
                                CHF_CONST_INT(m_funM));
}


// ---------------------------------------------------------
RealVect SphericalHarmonicGrad::setVectorFunctionMax(Real   a_bxWidth,
                                                       Real   a_outerRadius)
{
  CH_assert(isDefined());
  RealVect maxfun;
  FORT_MAXFUNCTIONSPHERICALGRAD(CHF_REALVECT(maxfun),
                                CHF_CONST_REAL(a_bxWidth),
                                CHF_CONST_REAL(a_outerRadius),
                                CHF_CONST_R1D(m_legendreCoeffs, m_legendreLen),
                                CHF_CONST_INT(m_funL),
                                CHF_CONST_INT(m_funM));
  return maxfun;
}

#include "NamespaceFooter.H"
