/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_OpenMP.H
  @brief  Declaration of various useful OpenMP-related utilities
  @author Robert Marskar
*/

#ifndef CD_OPENMP_H
#define CD_OPENMP_H

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
#pragma omp declare reduction(pairmin : std::pair<Real, RealVect> : ThreadSafePairMin(omp_out, omp_in)) \
  initializer(omp_priv = omp_orig)

// Thread-safe maximum of pair<Real, RealVect>.
inline void
ThreadSafePairMax(std::pair<Real, RealVect>& ompOut, const std::pair<Real, RealVect>& ompIn) noexcept
{
  ompOut = (ompIn.first > ompOut.first) ? ompIn : ompOut;
}
#pragma omp declare reduction(pairmax : std::pair<Real, RealVect> : ThreadSafePairMax(omp_out, omp_in)) \
  initializer(omp_priv = omp_orig)

#endif

#include <CD_NamespaceFooter.H>

#endif
