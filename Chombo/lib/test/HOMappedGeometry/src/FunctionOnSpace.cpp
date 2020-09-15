#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// FunctionOnSpace.cpp
// petermc, 7 Jun 2010

#include "FunctionOnSpace.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
Real
FunctionOnSpace::functionFromPhysical(/// physical coordinates
                                      const RealVect&   a_rv)
{
  IntVect iv = IntVect::Zero;
  Box bx(iv, iv);
  FArrayBox funFab(bx, 1);
  FArrayBox coordsFab(bx, SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      coordsFab(iv, idir) = a_rv[idir];
    }
  setFunctionFromPhysical(funFab, bx, coordsFab);
  Real funVal = funFab(iv, 0);
  return funVal;
}

#include "NamespaceFooter.H"
