/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/*
  Microbenchmark for the EBParticleMeshSoA deposit-stencil tightening.

  It deposits the SAME particle weights from the SAME ParticleSoA<> through two kernels that
  differ ONLY in the iterated stencil box:
    - TIGHT: EBParticleMeshSoA::depositWeight (anchored at the cloud edges -> 2^D cells for CIC,
             3^D for TSC at the standard width);
    - LOOSE: a local replica of the previous kernel (centered 3^D box for CIC, 5^D for TSC).
  Same storage, same per-cell weight math -> the measured difference is the stencil size only.

  All-regular EBISBox so no cut-cell branch is taken; forceIrregNGP = false. 3D, default width.
*/

// Std includes
#include <chrono>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

// Chombo includes
#include <ProblemDomain.H>
#include <DisjointBoxLayout.H>
#include <BoxIterator.H>
#include <EBIndexSpace.H>
#include <EBISLevel.H>
#include <EBISLayout.H>
#include <EBISBox.H>
#include <EBCellFAB.H>
#include <AllRegularService.H>
#include <RealVect.H>
#include <IntVect.H>
#include <parstream.H>

// Our includes
#include <CD_Initialize.H>
#include <CD_DepositionType.H>
#include <CD_ParticleOps.H>
#include <CD_BoxLoops.H>
#include <CD_ParticleSoA.H>
#include <CD_EBParticleMeshSoA.H>

using namespace ChomboDischarge;

namespace {

  using BenchSoA = ParticleSoA<>;

  volatile double g_sink = 0.0; // defeats dead-code elimination of the deposited field

  /**
    @brief LOOSE deposit replica: previous centered-box kernel (CIC 3^D / TSC 5^D), same math.
    @details Identical per-cell weight formula to EBParticleMeshSoA, but over the over-provisioned
    centered box anchored at particleIV -- this is the "before" the tightening removed.
  */
  void
  looseDepositWeight(FArrayBox&           a_rho,
                     const RealVect&      a_probLo,
                     const RealVect&      a_dx,
                     const BenchSoA&      a_soa,
                     const DepositionType a_type)
  {
    Real invVol = 1.0;
    for (int dir = 0; dir < SpaceDim; dir++) {
      invVol /= a_dx[dir];
    }
    const RealVect cicWidth = 1 * RealVect::Unit;
    const RealVect tscWidth = 2 * RealVect::Unit;
    const Box      cicBox(-1 * IntVect::Unit, 1 * IntVect::Unit);
    const Box      tscBox(-2 * IntVect::Unit, 2 * IntVect::Unit);
    Real           cicFactor = invVol;
    Real           tscFactor = invVol;
    for (int dir = 0; dir < SpaceDim; dir++) {
      cicFactor *= 1.0 / cicWidth[dir];
      tscFactor *= 2.0 / tscWidth[dir];
    }

    const std::size_t n = a_soa.size();

    if (a_type == DepositionType::CIC) {
      for (std::size_t i = 0; i < n; i++) {
        const Real     w    = a_soa.weight(i);
        const RealVect pos  = a_soa.position(i);
        const IntVect  piv  = ParticleOps::getParticleCellIndex(pos, a_probLo, a_dx);
        auto           kern = [&](const IntVect& iv) -> void {
          Real weight = cicFactor;
          for (int dir = 0; dir < SpaceDim; dir++) {
            const Real a = (a_probLo[dir] - pos[dir]) / a_dx[dir] + iv[dir];
            const Real b = a + 1.0;
            const Real L = 0.5 * cicWidth[dir];
            weight *= std::max(0.0, std::min(b, L) - std::max(a, -L));
          }
          a_rho(iv, 0) += weight * w;
        };
        BoxLoops::loop<D_DECL(1, 1, 1)>(cicBox + piv, kern);
      }
    }
    else { // TSC
      for (std::size_t i = 0; i < n; i++) {
        const Real     w    = a_soa.weight(i);
        const RealVect pos  = a_soa.position(i);
        const IntVect  piv  = ParticleOps::getParticleCellIndex(pos, a_probLo, a_dx);
        auto           kern = [&](const IntVect& iv) -> void {
          Real weight = tscFactor;
          for (int dir = 0; dir < SpaceDim; dir++) {
            const Real a      = (a_probLo[dir] - pos[dir]) / a_dx[dir] + iv[dir];
            const Real b      = a + 1.0;
            const Real L      = tscWidth[dir];
            const Real alpha  = std::max(a, -0.5 * L);
            const Real beta   = std::min(b, +0.5 * L);
            const Real factor = (alpha < beta) ? 1.0 : 0.0;
            weight *= factor * ((beta - alpha) - (beta * std::abs(beta) - alpha * std::abs(alpha)) / L);
          }
          a_rho(iv, 0) += weight * w;
        };
        BoxLoops::loop<D_DECL(1, 1, 1)>(tscBox + piv, kern);
      }
    }
  }

  /** @brief Accumulate the whole FAB into the sink so the deposits are not optimized away. */
  void
  drain(const EBCellFAB& a_rho, const Box& a_box)
  {
    const FArrayBox& fb  = a_rho.getFArrayBox();
    double           sum = 0.0;
    for (BoxIterator bit(a_box); bit.ok(); ++bit) {
      sum += fb(bit(), 0);
    }
    g_sink = g_sink + sum;
  }

} // namespace

int
main(int argc, char* argv[])
{
  initialize(argc, argv);

  {
    constexpr int       nCell = 32;
    constexpr int       ppc   = 16;
    constexpr int       reps  = 100;
    const int           nPart = nCell * nCell * nCell * ppc;
    const Box           domBox(IntVect::Zero, (nCell - 1) * IntVect::Unit);
    const ProblemDomain domain(domBox);
    const Real          dx     = 1.0;
    const RealVect      dxVect = dx * RealVect::Unit;
    const RealVect      probLo = RealVect::Zero;

    // ---- all-regular EBISBox ----
    Vector<Box>       boxes(1, domBox);
    Vector<int>       procs(1, 0);
    DisjointBoxLayout dbl(boxes, procs, domain);
    AllRegularService allReg;
    EBIndexSpace*     ebis = Chombo_EBIS::instance();
    ebis->define(domain, probLo, dx, allReg, nCell, 0);
    EBISLayout ebisl;
    ebis->fillEBISLayout(ebisl, dbl, domain, 2);
    DataIterator dit = dbl.dataIterator();
    dit.reset();
    const EBISBox& ebisbox = ebisl[dit()];

    EBParticleMeshSoA mesh(domain, domBox, ebisbox, dxVect, probLo);

    // ---- particles (interior, >= 3 cells from the boundary) ----
    std::mt19937                         rng(2024);
    std::uniform_real_distribution<Real> uPos(3.0, nCell - 3.0);
    std::uniform_real_distribution<Real> uW(0.5, 1.5);
    BenchSoA                             soa;
    soa.reserve(nPart);
    for (int i = 0; i < nPart; i++) {
      RealVect pos;
      for (int dir = 0; dir < SpaceDim; dir++) {
        pos[dir] = probLo[dir] + uPos(rng) * dx;
      }
      soa.append(pos, uW(rng), NoPayload{});
    }

    auto bench = [&](const DepositionType a_type, const bool a_tight) -> double {
      EBCellFAB rho(ebisbox, domBox, 1);
      rho.setVal(0.0);
      // warmup
      for (int r = 0; r < 3; r++) {
        if (a_tight) {
          mesh.depositWeight(rho, soa, a_type, 1.0, false);
        }
        else {
          looseDepositWeight(rho.getFArrayBox(), probLo, dxVect, soa, a_type);
        }
      }
      const auto t0 = std::chrono::steady_clock::now();
      for (int r = 0; r < reps; r++) {
        if (a_tight) {
          mesh.depositWeight(rho, soa, a_type, 1.0, false);
        }
        else {
          looseDepositWeight(rho.getFArrayBox(), probLo, dxVect, soa, a_type);
        }
      }
      const auto t1 = std::chrono::steady_clock::now();
      drain(rho, domBox);
      const double ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
      return ns / (static_cast<double>(reps) * static_cast<double>(nPart));
    };

    pout() << "Deposit-stencil microbenchmark (" << nPart << " particles, " << SpaceDim << "D, " << reps << " reps)"
           << endl;
    for (const auto tc : {std::make_pair(DepositionType::CIC, "CIC"), std::make_pair(DepositionType::TSC, "TSC")}) {
      const double loose = bench(tc.first, false);
      const double tight = bench(tc.first, true);
      pout() << "  " << tc.second << " : loose " << loose << " ns/p, tight " << tight << " ns/p, speedup "
             << (loose / tight) << "x" << endl;
    }
    pout() << "  (sink " << g_sink << ")" << endl;
  }

  return finalize();
}
