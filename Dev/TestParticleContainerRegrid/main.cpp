/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/*
  Stage-2b correctness test for ParticleContainerSoA::preRegrid()/regrid().

  Populates a 2-level grid (A), then regrids onto:
    B - a 2-level grid whose fine patch covers a DIFFERENT sub-region (so many particles change
        level and box, and -- under MPI -- rank);
    C - a single-level grid (the fine level is removed).
  After each regrid it asserts: particle count is conserved, no particles were lost off-domain,
  the conserved position-sum is unchanged (positions are never modified by regrid), and every
  particle sits in the box/level owning its cell on the NEW layout.

  Run single-rank AND under mpirun -np 2/4 (round-robin box distribution forces cross-rank moves).
*/

// Std includes
#include <cmath>
#include <random>

// Chombo includes
#include <ProblemDomain.H>
#include <DisjointBoxLayout.H>
#include <DataIterator.H>
#include <LayoutIterator.H>
#include <BoxIterator.H>
#include <RealVect.H>
#include <IntVect.H>
#include <SPMD.H>
#include <MayDay.H>
#include <parstream.H>

// Our includes
#include <CD_Initialize.H>
#include <CD_ParticleSoA.H>
#include <CD_ParticleContainerSoA.H>

using namespace ChomboDischarge;

namespace {

  using PC = ParticleContainerSoA<>; // NoPayload

  void
  require(const bool a_ok, const std::string& a_what)
  {
    if (!a_ok) {
      MayDay::Abort(("TestParticleContainerRegrid: check failed -- " + a_what).c_str());
    }
  }

  /** @brief DisjointBoxLayout tiling a region into blockingFactor-sized boxes, round-robin over ranks. */
  DisjointBoxLayout
  tiledDBL(const Box& a_region, const ProblemDomain& a_domain, const int a_bf)
  {
    const IntVect lo = a_region.smallEnd();
    const IntVect hi = a_region.bigEnd();

    Vector<Box> boxes;
    const Box   tileBox(IntVect::Zero, (hi - lo) / a_bf);
    for (BoxIterator bit(tileBox); bit.ok(); ++bit) {
      const IntVect blo = lo + bit() * a_bf;
      boxes.push_back(Box(blo, blo + (a_bf - 1) * IntVect::Unit));
    }

    Vector<int> procs(boxes.size());
    for (int i = 0; i < boxes.size(); i++) {
      procs[i] = i % numProc();
    }
    return DisjointBoxLayout(boxes, procs, a_domain);
  }

  IntVect
  cellOf(const RealVect& a_pos, const RealVect& a_probLo, const Real a_dx)
  {
    IntVect iv;
    for (int dir = 0; dir < SpaceDim; dir++) {
      iv[dir] = static_cast<int>(std::floor((a_pos[dir] - a_probLo[dir]) / a_dx));
    }
    return iv;
  }

  bool
  anyBoxContains(const DisjointBoxLayout& a_dbl, const IntVect& a_iv)
  {
    for (LayoutIterator lit = a_dbl.layoutIterator(); lit.ok(); ++lit) {
      if (a_dbl[lit()].contains(a_iv)) {
        return true;
      }
    }
    return false;
  }

  /** @brief Assert every particle sits in the box/level owning its cell (in-box + finest-level). */
  void
  checkOwnership(const PC& a_pc, const Vector<Real>& a_dxScalar, const RealVect& a_probLo)
  {
    const Vector<DisjointBoxLayout>& grids = a_pc.getGrids();
    for (int lvl = 0; lvl <= a_pc.getFinestLevel(); lvl++) {
      const PC::LevelParticles& level = a_pc[lvl];
      for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
        const Box       box  = grids[lvl][dit()];
        const PC::Leaf& leaf = level[dit()];
        for (std::size_t i = 0; i < leaf.size(); i++) {
          const RealVect pos = leaf.position(i);
          require(box.contains(cellOf(pos, a_probLo, a_dxScalar[lvl])), "regrid: particle outside its box");
          for (int lf = lvl + 1; lf <= a_pc.getFinestLevel(); lf++) {
            require(!anyBoxContains(grids[lf], cellOf(pos, a_probLo, a_dxScalar[lf])),
                    "regrid: particle not on its finest owning level");
          }
        }
      }
    }
  }

  /** @brief Global sum over all particles of (sum of position components), in double -- invariant
      under regrid (positions are never modified). Sums the raw double columns, not the Real-typed
      position(), so the checksum is precise regardless of the build's Real precision. */
  double
  globalPositionSum(const PC& a_pc)
  {
    double local = 0.0;
    for (int lvl = 0; lvl <= a_pc.getFinestLevel(); lvl++) {
      const PC::LevelParticles& level = a_pc[lvl];
      for (DataIterator dit = a_pc.getGrids()[lvl].dataIterator(); dit.ok(); ++dit) {
        const PC::Leaf& leaf = level[dit()];
        for (int dir = 0; dir < SpaceDim; dir++) {
          const double* xc = leaf.positionColumn(dir);
          for (std::size_t i = 0; i < leaf.size(); i++) {
            local += xc[i];
          }
        }
      }
    }
#ifdef CH_MPI
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, Chombo_MPI::comm);
    return global;
#else
    return local;
#endif
  }

} // namespace

int
main(int argc, char* argv[])
{
  initialize(argc, argv);

  {
    constexpr int  N      = 16;
    constexpr int  bf     = 8;
    constexpr int  refRat = 2;
    const Real     dx0    = 1.0;
    const Real     dx1    = dx0 / refRat;
    const RealVect probLo = RealVect::Zero;

    const Box           domBox0(IntVect::Zero, (N - 1) * IntVect::Unit);
    const Box           domBox1(IntVect::Zero, (refRat * N - 1) * IntVect::Unit);
    const ProblemDomain domain0(domBox0);
    const ProblemDomain domain1(domBox1);

    // Two-level grid metadata shared by layouts A and B.
    Vector<ProblemDomain> domains2(2);
    domains2[0] = domain0;
    domains2[1] = domain1;
    Vector<Real> dx2(2);
    dx2[0] = dx0;
    dx2[1] = dx1;
    Vector<int> refRat2(2);
    refRat2[0] = refRat;
    refRat2[1] = refRat;

    // ---- layout A: fine patch over the lower-centre [8,23] ----
    Vector<DisjointBoxLayout> gridsA(2);
    gridsA[0] = tiledDBL(domBox0, domain0, bf);
    gridsA[1] = tiledDBL(Box((N / 2) * IntVect::Unit, ((3 * N / 2) - 1) * IntVect::Unit), domain1, bf);

    PC pc;
    pc.define(gridsA, domains2, dx2, refRat2, probLo, bf, 1, "testRealm");

    // ---- populate each local box, canonicalize ----
    std::mt19937  rng(104729 + procID());
    constexpr int perBox = 4;
    for (int lvl = 0; lvl <= 1; lvl++) {
      for (DataIterator dit = gridsA[lvl].dataIterator(); dit.ok(); ++dit) {
        const Box box = gridsA[lvl][dit()];
        for (int k = 0; k < perBox; k++) {
          RealVect x;
          for (int dir = 0; dir < SpaceDim; dir++) {
            std::uniform_real_distribution<Real> u(box.smallEnd(dir) + 0.25, box.bigEnd(dir) + 0.75);
            x[dir] = u(rng) * dx2[lvl];
          }
          pc.addParticlesLocal(lvl, x, 1.0, NoPayload{});
        }
      }
    }
    pc.remap();
    const unsigned long long baseline = pc.getNumberOfValidParticlesGlobal();
    const double             posSum   = globalPositionSum(pc);
    require(baseline > 0, "regrid: particles populated");
    checkOwnership(pc, dx2, probLo);

    // ---- regrid to layout B: fine patch over the UPPER region [16,31] ----
    Vector<DisjointBoxLayout> gridsB(2);
    gridsB[0] = tiledDBL(domBox0, domain0, bf);
    gridsB[1] = tiledDBL(Box(N * IntVect::Unit, ((refRat * N) - 1) * IntVect::Unit), domain1, bf);

    pc.preRegrid();
    pc.regrid(gridsB, domains2, dx2, refRat2, bf, 1);
    require(pc.getNumberOfValidParticlesGlobal() == baseline, "regrid B: count conserved");
    require(pc.getNumberOfOutcastParticlesGlobal() == 0, "regrid B: nothing lost off-domain");
    require(std::abs(globalPositionSum(pc) - posSum) <= 1.e-10 * (1.0 + std::abs(posSum)),
            "regrid B: positions preserved");
    checkOwnership(pc, dx2, probLo);

    // ---- regrid to layout C: single level (fine level removed) ----
    Vector<DisjointBoxLayout> gridsC(1);
    gridsC[0] = tiledDBL(domBox0, domain0, bf);
    Vector<ProblemDomain> domainsC(1);
    domainsC[0] = domain0;
    Vector<Real> dxC(1);
    dxC[0] = dx0;
    Vector<int> refRatC(1);
    refRatC[0] = refRat;

    pc.preRegrid();
    pc.regrid(gridsC, domainsC, dxC, refRatC, bf, 0);
    require(pc.getFinestLevel() == 0, "regrid C: collapsed to one level");
    require(pc.getNumberOfValidParticlesGlobal() == baseline, "regrid C: count conserved");
    require(pc.getNumberOfOutcastParticlesGlobal() == 0, "regrid C: nothing lost off-domain");
    require(std::abs(globalPositionSum(pc) - posSum) <= 1.e-10 * (1.0 + std::abs(posSum)),
            "regrid C: positions preserved");
    checkOwnership(pc, dxC, probLo);

    pout() << "All ParticleContainerSoA regrid checks passed (" << SpaceDim << "D, " << numProc()
           << " rank(s), baseline " << baseline << ")." << endl;
  }

  return finalize();
}
