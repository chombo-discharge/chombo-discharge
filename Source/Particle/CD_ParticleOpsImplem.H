/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ParticleOpsImplem.H
  @brief  Implementation of CD_ParticleOps.H
  @author Robert Marskar
*/

#ifndef CD_PARTICLEOPSIMPLEM_H
#define CD_PARTICLEOPSIMPLEM_H

// Std includes
#include <cstdint>

// Chombo includes
#include <CH_Timer.H>
#include <PolyGeom.H>

// Our includes
#include <CD_ParticleOps.H>
#include <CD_ParallelOps.H>
#include <CD_PolyUtils.H>
#include <CD_Random.H>
#include <CD_NamespaceHeader.H>

inline IntVect
ParticleOps::getParticleCellIndex(const RealVect& a_particlePosition,
                                  const RealVect& a_probLo,
                                  const Real&     a_dx) noexcept
{
  return IntVect(D_DECL(std::floor((a_particlePosition[0] - a_probLo[0]) / a_dx),
                        std::floor((a_particlePosition[1] - a_probLo[1]) / a_dx),
                        std::floor((a_particlePosition[2] - a_probLo[2]) / a_dx)));
}

inline IntVect
ParticleOps::getParticleCellIndex(const RealVect& a_particlePosition,
                                  const RealVect& a_probLo,
                                  const RealVect& a_dx) noexcept
{
  return IntVect(D_DECL(std::floor((a_particlePosition[0] - a_probLo[0]) / a_dx[0]),
                        std::floor((a_particlePosition[1] - a_probLo[1]) / a_dx[1]),
                        std::floor((a_particlePosition[2] - a_probLo[2]) / a_dx[2])));
}

template <typename P, typename Traits>
inline void
ParticleOps::getPhysicalParticlesPerCell(EBAMRCellData& a_ppc, const ParticleContainerSoA<P, Traits>& a_src) noexcept
{
  CH_TIME("ParticleOps::getPhysicalParticlesPerCell(SoA)");

  const RealVect probLo = a_src.getProbLo();

  for (int lvl = 0; lvl <= a_src.getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = a_src.getGrids()[lvl];
    const DataIterator&      dit = dbl.dataIterator();
    const RealVect           dx  = a_src.getDx()[lvl];

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex&              din  = dit[mybox];
      const ParticleSoA<P, Traits>& leaf = a_src[lvl][din];

      FArrayBox& ppc = (*a_ppc[lvl])[din].getFArrayBox();
      ppc.setVal(0.0);

      for (std::size_t i = 0; i < leaf.size(); i++) {
        const IntVect iv = ParticleOps::getParticleGridCell(leaf.position(i), probLo, dx);

        ppc(iv, 0) += leaf.weight(i);
      }
    }
  }
}

template <typename P, typename Traits>
inline void
ParticleOps::getComputationalParticlesPerCell(EBAMRCellData&                         a_ppc,
                                              const ParticleContainerSoA<P, Traits>& a_src) noexcept
{
  CH_TIME("ParticleOps::getComputationalParticlesPerCell(SoA)");

  const RealVect probLo = a_src.getProbLo();

  for (int lvl = 0; lvl <= a_src.getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = a_src.getGrids()[lvl];
    const DataIterator&      dit = dbl.dataIterator();
    const RealVect           dx  = a_src.getDx()[lvl];

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex&              din  = dit[mybox];
      const ParticleSoA<P, Traits>& leaf = a_src[lvl][din];

      FArrayBox& ppc = (*a_ppc[lvl])[din].getFArrayBox();
      ppc.setVal(0.0);

      for (std::size_t i = 0; i < leaf.size(); i++) {
        const IntVect iv = ParticleOps::getParticleGridCell(leaf.position(i), probLo, dx);

        ppc(iv, 0) += 1.0;
      }
    }
  }
}

inline IntVect
ParticleOps::getParticleGridCell(const RealVect& a_particlePosition,
                                 const RealVect& a_probLo,
                                 const RealVect& a_dx) noexcept
{
  return IntVect(D_DECL(std::floor((a_particlePosition[0] - a_probLo[0]) / a_dx[0]),
                        std::floor((a_particlePosition[1] - a_probLo[1]) / a_dx[1]),
                        std::floor((a_particlePosition[2] - a_probLo[2]) / a_dx[2])));
}

inline bool
ParticleOps::domainIntersection(const RealVect& a_oldPos,
                                const RealVect& a_newPos,
                                const RealVect& a_probLo,
                                const RealVect& a_probHi,
                                Real&           a_s)
{

  // TLDR: This code does a boundary intersection test and returns where on the interval [oldPos, newPos] the intersection
  //       happened. We do this by checking if the particle moves towards a particular domain side and ends up outside of it.

  a_s = std::numeric_limits<Real>::max();

  bool crossedDomainBoundary = false;

  const RealVect path = a_newPos - a_oldPos;

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      const Side::LoHiSide side      = sit();
      const RealVect       wallPoint = (side == Side::Lo) ? a_probLo : a_probHi; // A point on the domain side
      const RealVect       n0        = sign(side) *
                          RealVect(
                            BASISV(dir)); // Normal vector pointing OUT of the domain on side sit and direction dir.
      const Real normPath = PolyGeom::dot(n0, path); // Component of the path that is normal to the domain edge/face.

      // If normPath > 0 then the particle trajectory points towards the domain edge/face and we can have an intersection.
      if (normPath > 0.0) {

        // s determines the intersection point between the particle path and the plane corresponding to the domain edge/face. Note that
        // we consider the edge/face to be an infinite plane and we just compute the intersection point between each edge/face and select the
        // closest intersection point.
        const Real s = PolyGeom::dot(wallPoint - a_oldPos, n0) / normPath;
        if (s >= 0.0 && s <= 1.0) {
          crossedDomainBoundary = true;
          if (s < a_s) {
            a_s = s;
          }
        }
      }
    }
  }

  return crossedDomainBoundary;
}

inline bool
ParticleOps::ebIntersectionBisect(const RefCountedPtr<BaseIF>& a_impFunc,
                                  const RealVect&              a_oldPos,
                                  const RealVect&              a_newPos,
                                  const Real&                  a_bisectStep,
                                  Real&                        a_s)
{

  // TLDR: We compute the intersection point using a bisection algorithm. We divide the full path into intervals and check if an interval
  //       has a root. If it does, we compute it using Brent's algorithm.

  a_s = std::numeric_limits<Real>::max();

  bool crossedEB = false;

  const Real     pathLen = (a_newPos - a_oldPos).vectorLength(); // Total path len
  const int      nsteps  = ceil(pathLen / a_bisectStep);         // Number of bisection intervals
  const RealVect dxStep  = (a_newPos - a_oldPos) / nsteps;       // Physical length of each bisection interval

  // Check each interval
  RealVect curPos = a_oldPos;
  for (int istep = 0; istep < nsteps; istep++) {
    const Real fa = a_impFunc->value(curPos); // Value of the implicit function at the start of the bisection interval
    const Real fb = a_impFunc->value(curPos +
                                     dxStep); // Value of the implicit function at the end of the bisection interval

    if (fa * fb <= 0.0) {

      // If this triggered we happen to know that f(pos+dxStep) > 0.0 and f(pos) < 0.0 and so we must have a root on the interval. We now compute the precise location
      // where the particle crossed the EB. For that we use a Brent root finder on the interval [pos, pos+dxStep]. This is a 1D problem.
      const RealVect intersectionPos = PolyUtils::brentRootFinder(a_impFunc, curPos, curPos + dxStep);
      a_s                            = (intersectionPos - a_oldPos).vectorLength() / pathLen;
      crossedEB                      = true;

      break;
    }
    else { // Move to next interval
      curPos += dxStep;
    }
  }

  return crossedEB;
}

inline bool
ParticleOps::ebIntersectionRaycast(const RefCountedPtr<BaseIF>& a_impFunc,
                                   const RealVect&              a_oldPos,
                                   const RealVect&              a_newPos,
                                   const Real&                  a_tolerance,
                                   Real&                        a_s)
{

  a_s = std::numeric_limits<Real>::max();

  bool ret = false;

  // Absolute distance to EB.
  auto dist = [&](const RealVect& x) -> Real {
    return std::abs(a_impFunc->value(x));
  };

  const Real D  = (a_newPos - a_oldPos).vectorLength(); // Total particle path length
  const Real D0 = dist(a_oldPos);                       // Distance to EB from starting position

  // If the distance to the EB from the starting position is smaller than the total path length, we need to check for intersections.
  if (D > D0) {

    const RealVect t = (a_newPos - a_oldPos) / D; // Particle trajectory.

    // Move a_oldPos along +t. If we end up too close to the boundary the particle has intersected the BC. Note that this does NOT check for whether or not
    // the particle moves tangential to the EB surface. The length of each step is the distance to the EB, so if the particle is close to the EB but moves
    // tangentially to it, this routine will be EXTREMELY slow.
    RealVect xa = a_oldPos;
    Real     r  = D;
    Real     d  = dist(xa);

    while (d < r) {

      if (d < a_tolerance) { // We collided.
        a_s = (xa - a_oldPos).vectorLength() / D;
        ret = true;
        break;
      }
      else { // We did not collide.
        xa += t * d;
        r -= d;
        d = dist(xa);
      }
    }
  }

  return ret;
}

template <typename P, typename Traits>
inline void
ParticleOps::copyDestructive(ParticleContainerSoA<P, Traits>& a_dst, ParticleContainerSoA<P, Traits>& a_src) noexcept
{
  CH_TIME("ParticleOps::copyDestructive(ParticleContainerSoA<P, Traits> x2)");

  CH_assert(a_dst.getRealm() == a_src.getRealm());

  for (int lvl = 0; lvl <= a_dst.getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = a_dst.getGrids()[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      a_dst[lvl][din].clear();                   // drop the destination leaf's particles
      a_dst[lvl][din].catenate(a_src[lvl][din]); // move the source leaf in (O(1) swap; src emptied)
    }
  }
}

template <typename P, typename Traits>
inline Real
ParticleOps::sum(const ParticleContainerSoA<P, Traits>& a_particles) noexcept
{
  CH_TIME("ParticleOps::sum(ParticleContainerSoA<P>)");

  Real particleSum = 0.0;

  for (int lvl = 0; lvl <= a_particles.getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = a_particles.getGrids()[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime) reduction(+ : particleSum)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex&              din  = dit[mybox];
      const ParticleSoA<P, Traits>& leaf = a_particles[lvl][din];

      for (std::size_t i = 0; i < leaf.size(); i++) {
        particleSum += leaf.weight(i);
      }
    }
  }

  return ParallelOps::sum(particleSum);
}

template <typename P, typename Traits>
inline void
ParticleOps::setData(ParticleContainerSoA<P, Traits>&                                 a_particles,
                     const std::function<void(ParticleSoA<P, Traits>&, std::size_t)>& a_functor) noexcept
{
  CH_TIME("ParticleOps::setData(ParticleContainerSoA<P>, std::function<void(ParticleSoA<P>&, std::size_t)>)");

  for (int lvl = 0; lvl <= a_particles.getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = a_particles.getGrids()[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      ParticleSoA<P, Traits>& leaf = a_particles[lvl][din];

      for (std::size_t i = 0; i < leaf.size(); i++) {
        a_functor(leaf, i);
      }
    }
  }
}

#include <CD_NamespaceFooter.H>

#endif
