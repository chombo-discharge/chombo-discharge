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
#include <List.H>

// Our includes
#include <CD_NamespaceHeader.H>

#ifdef _OPENMP

// Thread-safe reduction operator for taking the union of IntVectSet
inline void
ThreadSafeIVSUnion(IntVectSet& ivsInOut, const IntVectSet& ivsIn) noexcept
{
  ivsInOut |= ivsIn;
}
#pragma omp declare reduction(+ : IntVectSet : ThreadSafeIVSUnion(omp_out, omp_in))

// Thread-safe catenation for a particle list.
template <typename P>
inline void
ThreadSafeCatenation(List<P>& particlesOut, List<P>& particlesIn) noexcept
{
  particlesOut.catenate(particlesIn);
}

// Thread-safe join for a particle list.
template <typename P>
inline void
ThreadSafeJoin(List<P>& particlesOut, List<P>& particlesIn) noexcept
{
  particlesOut.join(particlesIn);
}

// Thread-safe minimum of pair<Real, RealVect>.
inline void
ThreadSafePairMin(std::pair<Real, RealVect>& ompOut, const std::pair<Real, RealVect>& ompIn) noexcept
{
  ompOut = (ompIn.first < ompOut.first) ? ompIn : ompOut;
}
#pragma omp declare reduction(pairmin : std::pair<Real, RealVect> : ThreadSafePairMin(omp_out, omp_in)) initializer(omp_priv=omp_orig)

#endif

#include <CD_NamespaceFooter.H>

#endif
