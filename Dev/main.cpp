/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   main.cpp
  @brief  Test application exercising the basic functionality of ParticleSoA.
  @author Robert Marskar
  @details
  A chombo-discharge executable that exercises the SoA particle container against
  the intended usage patterns:

    1. defining a particle type and adding particles to the container;
    2. iterating over particles and removing some (swap-and-pop idiom);
    3. a SIMD-friendly whole-column operation (multiplying fields together) launched
       through ParticleLoops::loop;
    4. merging two particles into one while conserving the center of mass.

  This is the harness we will grow to test deposition on FArrayBoxes, EB
  intersections, and merging on real Chombo objects.
*/

// Chombo includes
#include <MayDay.H>
#include <parstream.H>

// Our includes
#include <CD_Initialize.H>
#include <CD_DemoParticle.H>
#include <CD_ParticleLoops.H>

using namespace ChomboDischarge;

using DemoSoA = ParticleSoA<DemoParticle>;

namespace {

  /** @brief Always-on check (independent of NDEBUG) -- aborts via MayDay on failure. */
  void
  require(const bool a_ok, const std::string& a_what)
  {
    if (!a_ok) {
      MayDay::Abort(("main: check failed -- " + a_what).c_str());
    }
  }

  /** @brief (1) Add a_count particles to the container. */
  void
  addParticles(DemoSoA& a_soa, const std::size_t a_count)
  {
    a_soa.reserve(a_count);
    for (std::size_t i = 0; i < a_count; i++) {
      DemoParticle p;
      p.pos      = static_cast<Real>(i) * RealVect::Unit;
      p.weight   = static_cast<Real>(1 + (i % 4)); // weights 1..4
      p.mobility = 0.5 * static_cast<Real>(i);
      p.velocity = static_cast<Real>(i % 3) * RealVect::Unit;
      a_soa.append(p);
    }
  }

  /** @brief (2) Remove every particle whose weight is below a threshold (swap-and-pop). */
  std::size_t
  removeLightParticles(DemoSoA& a_soa, const Real a_weightThreshold)
  {
    std::size_t numRemoved = 0;

    const std::vector<Real>& w = a_soa.weights();
    for (std::size_t i = 0; i < a_soa.size();) {
      if (w[i] < a_weightThreshold) {
        a_soa.remove(i); // last particle is swapped into slot i
        numRemoved++;
      }
      else {
        i++;
      }
    }
    return numRemoved;
  }

  /** @brief (3) SIMD-friendly kernel: conductivity = weight * mobility over whole columns. */
  void
  computeConductivity(DemoSoA& a_soa)
  {
    const std::vector<Real>& weight       = a_soa.column<&DemoParticle::weight>();
    const std::vector<Real>& mobility     = a_soa.column<&DemoParticle::mobility>();
    std::vector<Real>&       conductivity = a_soa.column<&DemoParticle::conductivity>();

    ParticleLoops::loop(a_soa, [&](std::size_t i) {
      conductivity[i] = weight[i] * mobility[i];
    });
  }

  /**
    @brief (4) Merge two particles into one, conserving the center of mass.
    @details Treating weight as mass: W = w_i + w_j (total weight conserved) and
    X = (w_i x_i + w_j x_j) / W (center of mass conserved). Other fields are combined
    by the same weighted average. The merged particle overwrites slot i; slot j is
    removed (swap-and-pop), so we keep i < j to avoid disturbing i.
  */
  void
  mergeConservingCenterOfMass(DemoSoA& a_soa, std::size_t a_i, std::size_t a_j)
  {
    require(a_i != a_j, "merge needs two distinct particles");
    if (a_i > a_j) {
      std::swap(a_i, a_j);
    }

    std::vector<RealVect>& pos      = a_soa.positions();
    std::vector<Real>&     weight   = a_soa.weights();
    std::vector<RealVect>& velocity = a_soa.column<&DemoParticle::velocity>();
    std::vector<Real>&     mobility = a_soa.column<&DemoParticle::mobility>();

    const Real wi = weight[a_i];
    const Real wj = weight[a_j];
    const Real W  = wi + wj;

    RealVect mergedPos;
    RealVect mergedVel;
    for (int dir = 0; dir < SpaceDim; dir++) {
      mergedPos[dir] = (wi * pos[a_i][dir] + wj * pos[a_j][dir]) / W;
      mergedVel[dir] = (wi * velocity[a_i][dir] + wj * velocity[a_j][dir]) / W;
    }

    pos[a_i]      = mergedPos;
    velocity[a_i] = mergedVel;
    mobility[a_i] = (wi * mobility[a_i] + wj * mobility[a_j]) / W;
    weight[a_i]   = W;

    a_soa.remove(a_j);
  }

  /** @brief Total weight over the container. */
  Real
  totalWeight(const DemoSoA& a_soa)
  {
    Real sum = 0.0;
    for (const Real w : a_soa.weights()) {
      sum += w;
    }
    return sum;
  }

  /** @brief Weight-weighted center of mass over the container. */
  RealVect
  centerOfMass(const DemoSoA& a_soa)
  {
    const std::vector<RealVect>& pos = a_soa.positions();
    const std::vector<Real>&     w   = a_soa.weights();

    RealVect com = RealVect::Zero;
    Real     wt  = 0.0;
    for (std::size_t i = 0; i < a_soa.size(); i++) {
      for (int dir = 0; dir < SpaceDim; dir++) {
        com[dir] += w[i] * pos[i][dir];
      }
      wt += w[i];
    }
    com *= (1.0 / wt);
    return com;
  }

} // namespace

int
main(int argc, char* argv[])
{
  initialize(argc, argv);

  {
    DemoSoA soa;

    // (1) Add particles.
    addParticles(soa, 20);
    require(soa.size() == 20, "add: expected 20 particles");
    pout() << "added                : " << soa.size() << " particles, total weight " << totalWeight(soa) << endl;

    // (2) Iterate and remove the lightest particles (weight < 2).
    const std::size_t removed = removeLightParticles(soa, 2.0);
    require(soa.size() == 20 - removed, "remove: size bookkeeping");
    for (const Real w : soa.weights()) {
      require(w >= 2.0, "remove: a light particle survived");
    }
    pout() << "removed (weight < 2) : " << removed << " -> " << soa.size() << " remain" << endl;

    // (3) SIMD-friendly field multiply (conductivity = weight * mobility).
    computeConductivity(soa);
    {
      const std::vector<Real>& w     = soa.weights();
      const std::vector<Real>& mu    = soa.column<&DemoParticle::mobility>();
      const std::vector<Real>& sigma = soa.column<&DemoParticle::conductivity>();
      for (std::size_t i = 0; i < soa.size(); i++) {
        require(sigma[i] == w[i] * mu[i], "conductivity kernel result");
      }
      pout() << "conductivity kernel  : OK (sigma[0] = " << sigma[0] << ")" << endl;
    }

    // (4) Merge two particles, conserving total weight and center of mass.
    const Real     weightBefore = totalWeight(soa);
    const RealVect comBefore    = centerOfMass(soa);

    mergeConservingCenterOfMass(soa, 0, 1);

    const Real     weightAfter = totalWeight(soa);
    const RealVect comAfter    = centerOfMass(soa);

    const Real eps = 1e-12;
    require(std::abs(weightAfter - weightBefore) < eps, "merge: total weight conserved");
    for (int dir = 0; dir < SpaceDim; dir++) {
      require(std::abs(comAfter[dir] - comBefore[dir]) < eps, "merge: center of mass conserved");
    }
    pout() << "merged particles 0,1 : " << soa.size() << " remain; weight " << weightBefore << " -> " << weightAfter
           << ", COM conserved" << endl;

    pout() << "All ParticleSoA checks passed." << endl;
  }

  return finalize();
}
