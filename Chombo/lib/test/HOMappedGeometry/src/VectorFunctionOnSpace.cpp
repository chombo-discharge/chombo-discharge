#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// VectorFunctionOnSpace.cpp
// petermc, 30 Jun 2010

#include "VectorFunctionOnSpace.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
RealVect
VectorFunctionOnSpace::functionFromPhysical(/// physical coordinates
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
  setVectorFunctionFromPhysical(funFab, bx, coordsFab);
  RealVect funVal;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      funVal[idir] = funFab(iv, idir);
    }
  return funVal;
}

#include "NamespaceFooter.H"
