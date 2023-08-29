/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_OpenMP.H
  @brief  Declaration of various useful OpenMP-related utilities
  @author Robert Marskar
*/

#ifndef CD_OpenMP_H
#define CD_OpenMP_H

// Chombo includes
#include <IntVectSet.H>

// Our includes
#include <CD_NamespaceHeader.H>

// Thread-safe reduction operator for taking the union of IntVectSet
#ifdef _OPENMP
void
ThreadSafeIVSUnion(IntVectSet& ivsInOut, const IntVectSet& ivsIn)
{
  ivsInOut |= ivsIn;
}
#pragma omp declare reduction(+ : IntVectSet : ThreadSafeIVSUnion(omp_out, omp_in))
#endif

#include <CD_NamespaceFooter.H>

#endif
