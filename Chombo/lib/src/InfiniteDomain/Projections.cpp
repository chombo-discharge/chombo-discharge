#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Projections.cpp
// petermc, 2 Sep 2004

#include "BoxIterator.H"
#include "Projections.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
Box coarsenOuter(const Box&   a_bx,
                 int          a_refRatio)
{
  CH_assert(a_bx.ixType() == IndexType::TheNodeType());
  Box bxCoarse = coarsen(a_bx, a_refRatio);
  if (refine(bxCoarse, a_refRatio) == a_bx)
    return bxCoarse;

  IntVect loFine = a_bx.smallEnd();
  IntVect hiFine = a_bx.bigEnd();
  IntVect loCoarse = loFine / a_refRatio;
  IntVect hiCoarse = hiFine / a_refRatio;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (loCoarse[idir] * a_refRatio > loFine[idir])
        loCoarse.setVal(idir, loCoarse[idir]-1);
      if (hiCoarse[idir] * a_refRatio < hiFine[idir])
        hiCoarse.setVal(idir, hiCoarse[idir]+1);
    }
  Box bxCoarseOuter = Box(loCoarse, hiCoarse, IndexType::TheNodeType());
  return bxCoarseOuter;
}


// ---------------------------------------------------------
Box coarsenInner(const Box&   a_bx,
                 int          a_refRatio)
{
  CH_assert(a_bx.ixType() == IndexType::TheNodeType());
  Box bxCoarse = coarsen(a_bx, a_refRatio);
  if (refine(bxCoarse, a_refRatio) == a_bx)
    return bxCoarse;

  IntVect loFine = a_bx.smallEnd();
  IntVect hiFine = a_bx.bigEnd();
  IntVect loCoarse = loFine / a_refRatio;
  IntVect hiCoarse = hiFine / a_refRatio;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (loCoarse[idir] * a_refRatio < loFine[idir])
        loCoarse.setVal(idir, loCoarse[idir]+1);
      if (hiCoarse[idir] * a_refRatio > hiFine[idir])
        hiCoarse.setVal(idir, hiCoarse[idir]-1);
    }
  Box bxCoarseInner = Box(loCoarse, hiCoarse, IndexType::TheNodeType());
  return bxCoarseInner;
}


// ---------------------------------------------------------
void projectToCoarse(FArrayBox&         a_coarse,
                     const FArrayBox&   a_fine,
                     int                a_refRatio,
                     const IntVect&     a_shift)
{
  // Only for FArrayBoxes of one component.
  CH_assert(a_coarse.nComp() == 1);
  CH_assert(a_fine.nComp() == 1);
  if (a_refRatio == 1)
    {
      a_coarse.shift(+a_shift);
      a_coarse.copy(a_fine);
      a_coarse.shift(-a_shift);
    }
  else
    {
      Box bxCoarse = a_coarse.box();
      Box bxFine = a_fine.box();
      Box bxFineCoarsenedShifted = coarsenInner(bxFine - a_shift, a_refRatio);
      Box bxCoarseIntersect = bxCoarse & bxFineCoarsenedShifted;
      for (BoxIterator bit(bxCoarseIntersect); bit.ok(); ++bit)
        {
          IntVect ivCoarse = bit();
          a_coarse(ivCoarse, 0) = a_fine(ivCoarse * a_refRatio + a_shift, 0);
        }
    }
}


// ---------------------------------------------------------
void projectToFine(FArrayBox&         a_fine,
                   const FArrayBox&   a_coarse,
                   int                a_refRatio,
                   const IntVect&     a_shift)
{
  // Only for FArrayBoxes of one component.
  CH_assert(a_coarse.nComp() == 1);
  CH_assert(a_fine.nComp() == 1);
  if (a_refRatio == 1)
    {
      a_fine.shift(-a_shift);
      a_fine.copy(a_coarse);
      a_fine.shift(+a_shift);
    }
  else
    {
      Box bxCoarse = a_coarse.box();
      Box bxFine = a_fine.box();
      Box bxFineCoarsenedShifted = coarsenInner(bxFine - a_shift, a_refRatio);
      Box bxCoarseIntersect = bxCoarse & bxFineCoarsenedShifted;
      for (BoxIterator bit(bxCoarseIntersect); bit.ok(); ++bit)
        {
          IntVect ivCoarse = bit();
          a_fine(ivCoarse * a_refRatio + a_shift, 0) = a_coarse(ivCoarse, 0);
        }
    }
}

#include "NamespaceFooter.H"
