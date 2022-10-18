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

// Chombo includes
#include <SPMD.H>
#include <ParmParse.H>

// Our includes
#include <CD_Random.H>
#include <CD_NamespaceHeader.H>

inline void
Random::seed()
{

  // Set the RNG seed.
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

inline void
Random::setSeed(const int a_seed)
{
#ifdef CH_MPI
  const int seed = a_seed + procID();
#else
  const int seed = a_seed;
#endif

  s_rng = std::mt19937_64(seed);

  s_seeded = true;
}

inline void
Random::setRandomSeed()
{

  int seed = std::chrono::system_clock::now().time_since_epoch().count();

  // Special hook for MPI -- master rank broadcasts the seed and everyone increments by their processor ID.
#ifdef CH_MPI
  MPI_Bcast(&seed, 1, MPI_INT, 0, Chombo_MPI::comm);

  seed += procID();
#endif

  // Init the RNG
  s_rng = std::mt19937_64(seed);

  // Duh.
  s_seeded = true;
}

template <typename T>
inline T
Random::getPoisson(const Real a_mean)
{
  CH_assert(s_seeded);

  constexpr T zero = T(0);

  T ret = 0;

  if (a_mean < 250.0) {
    std::poisson_distribution<T> poisson(a_mean);

    ret = poisson(s_rng);
  }
  else {
    std::normal_distribution<Real> normal(a_mean, sqrt(a_mean));

    ret = T(normal(s_rng));

    ret = std::max(zero, ret);
  }

  return ret;
}

inline Real
Random::getUniformReal01()
{
  CH_assert(s_seeded);

  return s_uniform01(s_rng);
}

inline Real
Random::getUniformReal11()
{
  CH_assert(s_seeded);

  return s_uniform11(s_rng);
}

inline Real
Random::getNormal01()
{
  CH_assert(s_seeded);

  return s_normal01(s_rng);
}

inline RealVect
Random::getDirection()
{
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
  Real      x1   = 2.0;
  Real      x2   = 2.0;
  Real      r    = x1 * x1 + x2 * x2;
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
  RealVect pos = RealVect::Zero;

  for (int dir = 0; dir < SpaceDim; dir++) {
    pos[dir] = a_lo[dir] + 0.5 * (1.0 + Random::getUniformReal11()) * (a_hi[dir] - a_lo[dir]);
  }

  return pos;
}

#include <CD_NamespaceFooter.H>

#endif
