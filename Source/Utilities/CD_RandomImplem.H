/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RandomImplem.H
  @brief  Implementation of CD_Random.H
  @author Robert Marskar
*/

#ifndef CD_RandomImplem_H
#define CD_RandomImplem_H

// Std includes
#include <chrono>
#include <omp.h>

// Chombo includes
#include <SPMD.H>
#include <CH_Timer.H>
#include <ParmParse.H>

// Our includes
#include <CD_Random.H>
#include <CD_NamespaceHeader.H>

inline void
Random::seed()
{
  if (!s_seeded) {
    ParmParse pp("Random");

    if (pp.contains("seed")) {
      int seed;
      pp.get("seed", seed);
      if (seed < 0) {
        Random::setRandomSeed();
      }
      else {
        Random::setSeed(seed);
      }
    }
    else {
      Random::setSeed(0);
    }

    s_seeded = true;
  }
}

inline void
Random::setSeed(const int a_seed)
{
#ifdef CH_MPI
#ifdef _OPENMP
#pragma omp parallel
  {
    const int seed = a_seed + procID() * omp_get_num_threads() + omp_get_thread_num();

    s_rng = std::mt19937_64(seed);
  }
#else
  const int seed = a_seed + procID();

  s_rng = std::mt19937_64(seed);
#endif
#else
#ifdef _OPENMP
#pragma omp parallel
  {
    const int seed = a_seed + omp_get_thread_num();

    s_rng = std::mt19937_64(seed);
  }
#else
  const int seed = a_seed;

  s_rng = std::mt19937_64(seed);
#endif
#endif

  s_seeded = true;
}

inline void
Random::setRandomSeed()
{
  int seed = std::chrono::system_clock::now().time_since_epoch().count();

  // Special hook for MPI -- master rank broadcasts the seed and everyone increments by their processor ID.
#ifdef CH_MPI
  MPI_Bcast(&seed, 1, MPI_INT, 0, Chombo_MPI::comm);
#endif

  Random::setSeed(seed);
}

template <typename T, typename>
inline T
Random::getPoisson(const Real a_mean)
{
  CH_TIME("Random::getPoisson");

  CH_assert(s_seeded);

  T ret = (T)0;

  if (a_mean < 250.0) {
    std::poisson_distribution<T> poisson(a_mean);

    ret = poisson(s_rng);
  }
  else {
    std::normal_distribution<Real> normal(a_mean, sqrt(a_mean));

    ret = (T)std::max(normal(s_rng), (Real)0.0);
  }

  return ret;
}

template <typename T, typename>
inline T
Random::getBinomial(const T a_N, const Real a_p) noexcept
{
  CH_TIME("Random::getBinomial");

  CH_assert(s_seeded);

  T ret = 0;

  const bool useNormalApprox = (a_N > (T)9.0 * ((1.0 - a_p) / a_p)) && (a_N > (T)9.0 * a_p / (1.0 - a_p));

  if (useNormalApprox) {
    const Real mean = a_N * a_p;

    std::normal_distribution<Real> normalDist(mean, mean * (1.0 - a_p));

    ret = (T)std::max(normalDist(s_rng), (Real)0.0);
  }
  else {
    std::binomial_distribution<T> binomDist(a_N, a_p);

    ret = binomDist(s_rng);
  }

  return ret;
}

inline Real
Random::getUniformReal01()
{
  CH_TIME("Random::getUniformReal01");

  CH_assert(s_seeded);

  return s_uniform01(s_rng);
}

inline Real
Random::getUniformReal11()
{
  CH_TIME("Random::getUniformReal11");

  CH_assert(s_seeded);

  return s_uniform11(s_rng);
}

inline Real
Random::getNormal01()
{
  CH_TIME("Random::getNormal01");

  CH_assert(s_seeded);

  return s_normal01(s_rng);
}

inline RealVect
Random::getDirection()
{
  CH_TIME("Random::getDirection");

  CH_assert(s_seeded);

  constexpr Real safety = 1.E-12;

#if CH_SPACEDIM == 2
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1 * x1 + x2 * x2;
  while (r >= 1.0 || r < safety) {
    x1 = Random::getUniformReal11();
    x2 = Random::getUniformReal11();
    r  = x1 * x1 + x2 * x2;
  }

  return RealVect(x1, x2) / sqrt(r);
#elif CH_SPACEDIM == 3
  Real x1 = 2.0;
  Real x2 = 2.0;
  Real r  = x1 * x1 + x2 * x2;
  while (r >= 1.0 || r < safety) {
    x1 = Random::getUniformReal11();
    x2 = Random::getUniformReal11();
    r  = x1 * x1 + x2 * x2;
  }

  const Real x = 2 * x1 * sqrt(1 - r);
  const Real y = 2 * x2 * sqrt(1 - r);
  const Real z = 1 - 2 * r;

  return RealVect(x, y, z);
#endif
}

template <typename T>
inline Real
Random::get(T& a_distribution)
{
  CH_TIME("Random::get");

  CH_assert(s_seeded);

  return a_distribution(s_rng);
}

template <typename T>
inline size_t
Random::getDiscrete(T& a_distribution)
{
  CH_TIME("Random::getDiscrete");

  CH_assert(s_seeded);

  return a_distribution(s_rng);
}

inline RealVect
Random::randomPosition(const RealVect a_cellPos,
                       const RealVect a_lo,
                       const RealVect a_hi,
                       const RealVect a_bndryCentroid,
                       const RealVect a_bndryNormal,
                       const Real     a_dx,
                       const Real     a_kappa) noexcept
{
  CH_TIME("Random::randomPosition(full)");

  RealVect pos;

  if (a_kappa < 1.0) { // Rejection sampling.
    pos = Random::randomPosition(a_lo, a_hi, a_bndryCentroid, a_bndryNormal);
  }
  else { // Regular cell. Get a position.
    pos = Random::randomPosition(a_lo, a_hi);
  }

  // Convert from unit cell coordinates to physical coordinates.
  pos = a_cellPos + pos * a_dx;

  return pos;
}

inline RealVect
Random::randomPosition(const RealVect a_lo,
                       const RealVect a_hi,
                       const RealVect a_bndryCentroid,
                       const RealVect a_bndryNormal) noexcept
{
  CH_TIME("Random::randomPosition(partial)");

  // Draw a position.
  RealVect pos   = Random::randomPosition(a_lo, a_hi);
  bool     valid = (pos - a_bndryCentroid).dotProduct(a_bndryNormal) >= 0.0;

  // Rejectino sampling.
  while (!valid) {
    pos   = Random::randomPosition(a_lo, a_hi);
    valid = (pos - a_bndryCentroid).dotProduct(a_bndryNormal) >= 0.0;
  }

  return pos;
}

inline RealVect
Random::randomPosition(const RealVect a_lo, const RealVect a_hi) noexcept
{
  CH_TIME("Random::randomPosition(box)");

  RealVect pos = RealVect::Zero;

  for (int dir = 0; dir < SpaceDim; dir++) {
    pos[dir] = a_lo[dir] + Random::getUniformReal01() * (a_hi[dir] - a_lo[dir]);
  }

  return pos;
}

#include <CD_NamespaceFooter.H>

#endif
