#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// SphericalHarmonic.cpp
// petermc, 21 Mar 2008

// #include <math.h>
#include "SphericalHarmonic.H"
#include "FunSphericalF_F.H"
#include "LegendreF.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
// default constructor
SphericalHarmonic::SphericalHarmonic()
{
  //  setDefaultValues();
}


// ---------------------------------------------------------
void SphericalHarmonic::define(int          a_funL,
                               int          a_funM)
{
  CH_assert(a_funL >= a_funM);
  m_funL = a_funL;
  m_funM = a_funM;
  m_legendreLen = m_funL - a_funM + 1;
  m_legendreCoeffs = new Real[m_legendreLen];
  FORT_LEGENDRECOEFFS(m_legendreCoeffs, &m_funL, &m_funM);
  m_isDefined = true;
}


// ---------------------------------------------------------
// destructor
SphericalHarmonic::~SphericalHarmonic()
{
  if (isDefined())
    {
      delete [] m_legendreCoeffs;
    }
}


// ---------------------------------------------------------
void SphericalHarmonic::setFunctionFromPhysical(FArrayBox&         a_funFab,
                                                const Box&         a_bx,
                                                const FArrayBox&   a_coordsFab)
{
  CH_assert(isDefined());
  FORT_SETFUNCTIONSPHERICAL(CHF_FRA(a_funFab),
                            CHF_CONST_FRA(a_coordsFab),
                            CHF_BOX(a_bx),
                            CHF_CONST_R1D(m_legendreCoeffs, m_legendreLen),
                            CHF_CONST_INT(m_funL),
                            CHF_CONST_INT(m_funM));
}


// ---------------------------------------------------------
Real SphericalHarmonic::setFunctionMax(Real   a_bxWidth,
                                       Real   a_outerRadius)
{
  CH_assert(isDefined());
  Real maxfun;
  FORT_MAXFUNCTIONSPHERICAL(CHF_REAL(maxfun),
                            CHF_CONST_REAL(a_bxWidth),
                            CHF_CONST_REAL(a_outerRadius),
                            CHF_CONST_R1D(m_legendreCoeffs, m_legendreLen),
                            CHF_CONST_INT(m_funL),
                            CHF_CONST_INT(m_funM));
  return maxfun;
}

#include "NamespaceFooter.H"
