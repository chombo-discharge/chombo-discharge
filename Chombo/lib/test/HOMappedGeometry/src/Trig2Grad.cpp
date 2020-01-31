#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Trig2Grad.cpp
// petermc, 30 Jun 2010

// #include <math.h>
#include "Trig2Grad.H"
#include "FunTrig2GradF_F.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
// default constructor
Trig2Grad::Trig2Grad()
{
  //  setDefaultValues();
  m_isDefined = true;
}


// ---------------------------------------------------------
// destructor
Trig2Grad::~Trig2Grad()
{
}


// ---------------------------------------------------------
void Trig2Grad::setVectorFunctionFromPhysical(FArrayBox&         a_funFab,
                                                            const Box&         a_bx,
                                                            const FArrayBox&   a_coordsFab)
{
  CH_assert(isDefined());
  FORT_SETFUNCTIONTRIG2GRAD(CHF_FRA(a_funFab),
                            CHF_CONST_FRA(a_coordsFab),
                            CHF_BOX(a_bx));
}


// ---------------------------------------------------------
RealVect Trig2Grad::setVectorFunctionMax()
{
  CH_assert(isDefined());
  RealVect maxfun = M_PI * RealVect::Unit;
  return maxfun;
}

#include "NamespaceFooter.H"
