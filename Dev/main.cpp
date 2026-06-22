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
  A chombo-discharge executable that exercises the merged arena-backed SoA particle
  container against the intended usage patterns, so the class can be explored end to end:

    1. defining a payload type and adding particles (append(pos, weight, payload));
    2. iterating over particles and removing some (swap-and-pop idiom);
    3. a SIMD-friendly whole-column operation launched through ParticleLoops::loop;
    4. merging two particles into one while conserving the center of mass;
    5. cell-sorting the container (counting sort + CSR offsets) and walking the cells;
    6. an MPI-style linearize / delinearize round-trip of a single particle.

  Note how position/weight/id/rank are container-owned: the payload struct (DemoPayload)
  declares only the extra fields, and the mandatory data is reached via accessors.
*/

// Std includes
#include <cmath>

// Chombo includes
#include <Box.H>
#include <IntVect.H>
#include <MayDay.H>
#include <parstream.H>

// Our includes
#include <CD_Initialize.H>
#include <CD_DemoParticle.H>
#include <CD_ParticleLoops.H>

using namespace ChomboDischarge;

using DemoSoA = ParticleSoA<DemoPayload>;

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
      DemoPayload payload;
      payload.mobility = 0.5 * static_cast<Real>(i);
      const Real v     = static_cast<Real>(i % 3);
      D_TERM(payload.vx = v;, payload.vy = v;, payload.vz = v;);

      const RealVect pos    = static_cast<Real>(i) * RealVect::Unit;
      const Real     weight = static_cast<Real>(1 + (i % 4)); // weights 1..4
      a_soa.append(pos, weight, payload);
    }
  }

  /** @brief (2) Remove every particle whose weight is below a threshold (swap-and-pop). */
  std::size_t
  removeLightParticles(DemoSoA& a_soa, const Real a_weightThreshold)
  {
    std::size_t numRemoved = 0;
    for (std::size_t i = 0; i < a_soa.size();) {
      if (a_soa.weight(i) < a_weightThreshold) {
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
    const double*       weight       = a_soa.weightColumn();
    const ParticleReal* mobility     = a_soa.column<&DemoPayload::mobility>();
    ParticleReal*       conductivity = a_soa.column<&DemoPayload::conductivity>();

    ParticleLoops::loop(a_soa, [&](std::size_t i) {
      conductivity[i] = weight[i] * mobility[i];
    });
  }

  /**
    @brief (4) Merge two particles into one, conserving the center of mass.
    @details Treating weight as mass: W = w_i + w_j (total weight conserved) and
    X = (w_i x_i + w_j x_j) / W (center of mass conserved). Other fields are combined by
    the same weighted average. The merged particle overwrites slot i; slot j is removed
    (swap-and-pop), so we keep i < j to avoid disturbing i.
  */
  void
  mergeConservingCenterOfMass(DemoSoA& a_soa, std::size_t a_i, std::size_t a_j)
  {
    require(a_i != a_j, "merge needs two distinct particles");
    if (a_i > a_j) {
      std::swap(a_i, a_j);
    }

    D_TERM(ParticleReal* vx = a_soa.column<&DemoPayload::vx>();, ParticleReal* vy = a_soa.column<&DemoPayload::vy>();
           , ParticleReal* vz                                                     = a_soa.column<&DemoPayload::vz>(););
    ParticleReal* mobility = a_soa.column<&DemoPayload::mobility>();

    const Real wi = a_soa.weight(a_i);
    const Real wj = a_soa.weight(a_j);
    const Real W  = wi + wj;

    const RealVect xi = a_soa.position(a_i);
    const RealVect xj = a_soa.position(a_j);

    const RealVect vi = RealVect(D_DECL(vx[a_i], vy[a_i], vz[a_i]));
    const RealVect vj = RealVect(D_DECL(vx[a_j], vy[a_j], vz[a_j]));

    RealVect mergedPos;
    RealVect mergedVel;
    for (int dir = 0; dir < SpaceDim; dir++) {
      mergedPos[dir] = (wi * xi[dir] + wj * xj[dir]) / W;
      mergedVel[dir] = (wi * vi[dir] + wj * vj[dir]) / W;
    }

    a_soa.setPosition(a_i, mergedPos);
    D_TERM(vx[a_i] = mergedVel[0];, vy[a_i] = mergedVel[1];, vz[a_i] = mergedVel[2];);
    mobility[a_i]     = (wi * mobility[a_i] + wj * mobility[a_j]) / W;
    a_soa.weight(a_i) = W;

    a_soa.remove(a_j);
  }

  /** @brief Total weight over the container. */
  Real
  totalWeight(const DemoSoA& a_soa)
  {
    const double* w   = a_soa.weightColumn();
    Real          sum = 0.0;
    for (std::size_t i = 0; i < a_soa.size(); i++) {
      sum += w[i];
    }
    return sum;
  }

  /** @brief Weight-weighted center of mass over the container. */
  RealVect
  centerOfMass(const DemoSoA& a_soa)
  {
    RealVect com = RealVect::Zero;
    Real     wt  = 0.0;
    for (std::size_t i = 0; i < a_soa.size(); i++) {
      const RealVect x = a_soa.position(i);
      const Real     w = a_soa.weight(i);
      for (int dir = 0; dir < SpaceDim; dir++) {
        com[dir] += w * x[dir];
      }
      wt += w;
    }
    com *= (1.0 / wt);
    return com;
  }

  /** @brief Fortran linear cell index of a position within a box (matches ParticleSoA::sortByCell). */
  std::size_t
  cellOf(const RealVect& a_x, const Box& a_box, const RealVect& a_dx, const RealVect& a_probLo)
  {
    const IntVect loEnd   = a_box.smallEnd();
    const IntVect boxSize = a_box.size();
    std::size_t   lin     = 0;
    std::size_t   stride  = 1;
    for (int d = 0; d < SpaceDim; d++) {
      int c = static_cast<int>(std::floor((a_x[d] - a_probLo[d]) / a_dx[d])) - loEnd[d];
      if (c < 0) {
        c = 0;
      }
      if (c >= boxSize[d]) {
        c = boxSize[d] - 1;
      }
      lin += static_cast<std::size_t>(c) * stride;
      stride *= static_cast<std::size_t>(boxSize[d]);
    }
    return lin;
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
    for (std::size_t i = 0; i < soa.size(); i++) {
      require(soa.weight(i) >= 2.0, "remove: a light particle survived");
    }
    pout() << "removed (weight < 2) : " << removed << " -> " << soa.size() << " remain" << endl;

    // (3) SIMD-friendly field multiply (conductivity = weight * mobility).
    computeConductivity(soa);
    {
      const double*       w     = soa.weightColumn();
      const ParticleReal* mu    = soa.column<&DemoPayload::mobility>();
      const ParticleReal* sigma = soa.column<&DemoPayload::conductivity>();
      for (std::size_t i = 0; i < soa.size(); i++) {
        // Cast the recomputed product to ParticleReal so the comparison is exact in any payload
        // precision (the kernel stored weight*mobility into a ParticleReal column).
        require(sigma[i] == static_cast<ParticleReal>(w[i] * mu[i]), "conductivity kernel result");
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

    // (5) Cell-sort: bin a fresh container into a 4^SpaceDim box and walk the CSR ranges.
    {
      const Box      box(IntVect::Zero, 3 * IntVect::Unit); // 4 cells per direction
      const RealVect dx     = RealVect::Unit;
      const RealVect probLo = RealVect::Zero;

      DemoSoA cells;
      // Scatter a handful of particles across cells (some sharing a cell).
      for (int n = 0; n < 30; n++) {
        const RealVect pos = RealVect(D_DECL(0.5 + (n % 4), 0.5 + ((n / 4) % 4), 0.5 + ((n / 16) % 4)));
        cells.append(pos, 1.0);
      }
      require(!cells.isSorted(), "cell-sort: container is unsorted before sortByCell");

      cells.sortByCell(box, dx, probLo);
      require(cells.isSorted(), "cell-sort: container is sorted after sortByCell");
      require(cells.numCells() == box.numPts(), "cell-sort: numCells matches the box");

      // Every particle in cell c's CSR range must actually live in cell c, and the
      // per-cell counts must sum to size().
      std::size_t total    = 0;
      std::size_t occupied = 0;
      for (std::size_t c = 0; c < cells.numCells(); c++) {
        const std::size_t lo = cells.cellStart(c);
        const std::size_t hi = cells.cellStart(c + 1);
        for (std::size_t k = lo; k < hi; k++) {
          require(cellOf(cells.position(k), box, dx, probLo) == c, "cell-sort: particle in wrong cell range");
        }
        total += (hi - lo);
        if (hi > lo) {
          occupied++;
        }
      }
      require(total == cells.size(), "cell-sort: CSR ranges cover every particle");
      pout() << "cell-sorted          : " << cells.size() << " particles into " << occupied << "/" << cells.numCells()
             << " cells" << endl;
    }

    // (6) MPI-style linearize / delinearize round-trip of one particle.
    {
      std::vector<unsigned char> buffer(DemoSoA::bytesPerParticle());
      soa.particleID(2) = 4242; // give it identifiable metadata
      soa.rankID(2)     = 7;
      soa.linearizeParticle(buffer.data(), 2);

      DemoSoA received;
      received.delinearizeAndAppend(buffer.data());
      require(received.size() == 1, "linearize: one particle received");
      require(received.particleID(0) == 4242, "linearize: particleID round-trips");
      require(received.rankID(0) == 7, "linearize: rankID round-trips");
      for (int dir = 0; dir < SpaceDim; dir++) {
        require(std::abs(received.position(0)[dir] - soa.position(2)[dir]) < eps, "linearize: position round-trips");
      }
      require(received.weight(0) == soa.weight(2), "linearize: weight round-trips");
      pout() << "linearize round-trip : OK (id " << received.particleID(0) << ", rank " << received.rankID(0) << ")"
             << endl;
    }

    pout() << "All ParticleSoA checks passed." << endl;
  }

  return finalize();
}
