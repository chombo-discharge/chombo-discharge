#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Trig2.cpp
// petermc, 3 Sep 2010

// #include <math.h>
#include "Trig2.H"
#include "FunTrig2F_F.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
// default constructor
Trig2::Trig2()
{
  //  setDefaultValues();
  m_isDefined = true;
}


// ---------------------------------------------------------
// destructor
Trig2::~Trig2()
{
}


// ---------------------------------------------------------
void Trig2::setFunctionFromPhysical(FArrayBox&         a_funFab,
                                    const Box&         a_bx,
                                    const FArrayBox&   a_coordsFab)
{
  CH_assert(isDefined());
  FORT_SETFUNCTIONTRIG2(CHF_FRA(a_funFab),
                        CHF_CONST_FRA(a_coordsFab),
                        CHF_BOX(a_bx));
}


// ---------------------------------------------------------
Real Trig2::setFunctionMax()
{
  CH_assert(isDefined());
  Real maxfun = Real(SpaceDim);
  return maxfun;
}

#include "NamespaceFooter.H"
