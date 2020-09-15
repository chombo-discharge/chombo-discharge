#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "IndicesFunctions.H"

#include "NamespaceHeader.H"

IntVect inversePermutation(const IntVect& a_permutation)
{
  CH_assert(isPermutationVect(a_permutation));
  IntVect permInverse;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      permInverse[a_permutation[idir]] = idir;
    }
  return permInverse;
}


bool isPermutationVect(const IntVect& a_permutation)
{
  // First check that every component in a_permutation is in range.
  bool allInRange = true;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      allInRange &= inDimensionRange(a_permutation[idir]);
    }
  if (!allInRange)
    {
      return false;
    }
  // Now check that the components of a_permutation include the whole range.
  IntVect foundVect = IntVect::Zero;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      foundVect[a_permutation[idir]] = 1;
    }
  return (foundVect == IntVect::Unit);
}

bool isSignVect(const IntVect& a_sign)
{
  return (absolute(a_sign) == IntVect::Unit);
}

#include "NamespaceFooter.H"
