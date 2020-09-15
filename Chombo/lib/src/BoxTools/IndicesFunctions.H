#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _INDICESFUNCTIONS_H_
#define _INDICESFUNCTIONS_H_

#include "Box.H"

#include "NamespaceHeader.H"

inline size_t FortranArrayIndex(IntVect a_iv,
                                IntVect a_dims)
{
  CH_assert(a_iv >= IntVect::Zero);
  CH_assert(a_dims >= IntVect::Zero);
  // Fortran array order: last dimension changes most slowly.
  size_t ind = a_iv[CH_SPACEDIM-1];
  for (int idir = CH_SPACEDIM-2; idir >= 0; idir--)
    {
      ind *= a_dims[idir];
      ind += a_iv[idir];
    }
  return ind;
}

inline size_t FortranArrayIndex(IntVect a_iv,
                                Box a_bx)
{
  IntVect dims = a_bx.size();
  IntVect offset = a_iv - a_bx.smallEnd();
  return FortranArrayIndex(offset, dims);
}

inline IntVect FortranArrayIndex(size_t a_ind,
                                 IntVect a_dims)
{
  IntVect iv;
  // Fortran array order: first dimension changes fastest
  size_t indRemaining = a_ind;
  for (int idir = 0; idir < CH_SPACEDIM-1; idir++)
    {
      iv[idir] = indRemaining % a_dims[idir];
      indRemaining -= iv[idir];
      indRemaining /= a_dims[idir];
    }
  iv[CH_SPACEDIM-1] = indRemaining;
  return iv;
}

// Return inverse permutation of a_permutation.
IntVect inversePermutation(const IntVect& a_permutation);

// Test whether the IntVect represents a legal permutation of dimensions.
bool isPermutationVect(const IntVect& a_permutation);

// Test whether the IntVect represents a legal sign vector.
bool isSignVect(const IntVect& a_sign);

// Test whether the index is within the range of dimensions.
inline bool inDimensionRange(int a_ind)
{
  return ((a_ind >= 0) && (a_ind < CH_SPACEDIM));
}

#include "NamespaceFooter.H"

#endif
