/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/*
  Stage-2a correctness test for ParticleContainerSoA::remap().

  Builds a MULTI-BOX, 2-level grid (each level split into blockingFactor-sized boxes and
  load-balanced across ranks; the fine level covers a sub-region), then:
    1. populates each local box, remap()s to canonicalize, and records the baseline count;
    2. moves every particle to a new random in-domain position and remap()s, asserting count
       conservation, zero outcast, and that every particle now sits in the box/level that owns its
       cell (in-box + finest-level checks across the GLOBAL layout);
    3. moves every particle off-domain and remap()s, asserting they all become outcasts.

  Designed to run single-rank AND under mpirun -np 2/4 (load balancing distributes boxes, so the
  random moves force cross-rank transfers through the Alltoallv scatter).
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
      MayDay::Abort(("TestParticleContainerRemap: check failed -- " + a_what).c_str());
    }
  }

  /** @brief DisjointBoxLayout tiling a region into blockingFactor-sized boxes, round-robin over ranks.
      @details Built by hand (no domainSplit/LoadBalance) to keep the Dev test self-contained; the
      region extents must be multiples of the blocking factor. */
  DisjointBoxLayout
  tiledDBL(const Box& a_region, const ProblemDomain& a_domain, const int a_bf)
  {
    const IntVect lo = a_region.smallEnd();
    const IntVect hi = a_region.bigEnd();

    Vector<Box> boxes;
    const Box   tileBox(IntVect::Zero, (hi - lo) / a_bf); // one cell per blockingFactor-sized tile
    for (BoxIterator bit(tileBox); bit.ok(); ++bit) {
      const IntVect blo = lo + bit() * a_bf;
      boxes.push_back(Box(blo, blo + (a_bf - 1) * IntVect::Unit));
    }

    Vector<int> procs(boxes.size());
    for (int i = 0; i < boxes.size(); i++) {
      procs[i] = i % numProc(); // round-robin so boxes distribute across ranks under MPI
    }
    return DisjointBoxLayout(boxes, procs, a_domain);
  }

  /** @brief Cell index of a position on a level. */
  IntVect
  cellOf(const RealVect& a_pos, const RealVect& a_probLo, const Real a_dx)
  {
    IntVect iv;
    for (int dir = 0; dir < SpaceDim; dir++) {
      iv[dir] = static_cast<int>(std::floor((a_pos[dir] - a_probLo[dir]) / a_dx));
    }
    return iv;
  }

  /** @brief Whether ANY box in the (global) layout contains the cell. */
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

  /** @brief Assert every particle sits in the box/level that owns its cell (in-box + finest-level). */
  void
  checkOwnership(const PC& a_pc, const Vector<Real>& a_dxScalar, const RealVect& a_probLo)
  {
    for (int lvl = 0; lvl <= a_pc.getFinestLevel(); lvl++) {
      const PC::LevelParticles&        level = a_pc[lvl];
      const Vector<DisjointBoxLayout>& grids = a_pc.getGrids();
      for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
        const Box       box  = grids[lvl][dit()];
        const PC::Leaf& leaf = level[dit()];
        for (std::size_t i = 0; i < leaf.size(); i++) {
          const RealVect pos = leaf.position(i);

          // (a) in its own box on this level.
          require(box.contains(cellOf(pos, a_probLo, a_dxScalar[lvl])), "remap: particle outside its box");

          // (b) on the FINEST level that owns it: no finer level's box contains it.
          for (int lf = lvl + 1; lf <= a_pc.getFinestLevel(); lf++) {
            require(!anyBoxContains(grids[lf], cellOf(pos, a_probLo, a_dxScalar[lf])),
                    "remap: particle not on its finest owning level");
          }
        }
      }
    }
  }

  /** @brief Move every local particle to a fresh random position in [a_lo, a_hi)^D (physical). */
  void
  scatterPositions(PC& a_pc, std::mt19937& a_rng, const Real a_lo, const Real a_hi)
  {
    std::uniform_real_distribution<Real> u(a_lo, a_hi);
    for (int lvl = 0; lvl <= a_pc.getFinestLevel(); lvl++) {
      PC::LevelParticles& level = a_pc[lvl];
      for (DataIterator dit = a_pc.getGrids()[lvl].dataIterator(); dit.ok(); ++dit) {
        PC::Leaf& leaf = level[dit()];
        for (std::size_t i = 0; i < leaf.size(); i++) {
          RealVect x;
          for (int dir = 0; dir < SpaceDim; dir++) {
            x[dir] = u(a_rng);
          }
          leaf.setPosition(i, x);
        }
      }
    }
  }

} // namespace

int
main(int argc, char* argv[])
{
  initialize(argc, argv);

  {
    constexpr int  N      = 16;
    constexpr int  bf     = 8; // blocking factor (tile/box size)
    constexpr int  refRat = 2;
    const Real     dx0    = 1.0;
    const Real     dx1    = dx0 / refRat;
    const RealVect probLo = RealVect::Zero;

    const Box           domBox0(IntVect::Zero, (N - 1) * IntVect::Unit);
    const Box           domBox1(IntVect::Zero, (refRat * N - 1) * IntVect::Unit);
    const ProblemDomain domain0(domBox0);
    const ProblemDomain domain1(domBox1);
    const Box           fineRegion((N / 2) * IntVect::Unit, ((3 * N / 2) - 1) * IntVect::Unit); // fine cells [8,23]

    Vector<DisjointBoxLayout> grids(2);
    grids[0] = tiledDBL(domBox0, domain0, bf);
    grids[1] = tiledDBL(fineRegion, domain1, bf);

    Vector<ProblemDomain> domains(2);
    domains[0] = domain0;
    domains[1] = domain1;

    Vector<Real> dx(2);
    dx[0] = dx0;
    dx[1] = dx1;

    Vector<int> refRats(2);
    refRats[0] = refRat;
    refRats[1] = refRat;

    PC pc;
    pc.define(grids, domains, dx, refRats, probLo, bf, 1, "testRealm");

    // ---- populate: a few particles inside each LOCAL box, on each level ----
    std::mt19937  rng(7919 + procID());
    constexpr int perBox = 4;
    for (int lvl = 0; lvl <= 1; lvl++) {
      for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
        const Box box = grids[lvl][dit()];
        for (int k = 0; k < perBox; k++) {
          RealVect x;
          for (int dir = 0; dir < SpaceDim; dir++) {
            std::uniform_real_distribution<Real> u(box.smallEnd(dir) + 0.25, box.bigEnd(dir) + 0.75);
            x[dir] = u(rng) * dx[lvl]; // random position inside this box (physical)
          }
          pc.addParticlesLocal(lvl, x, 1.0, NoPayload{});
        }
      }
    }

    // ---- canonicalize, record baseline ----
    pc.remap();
    const unsigned long long baseline = pc.getNumberOfValidParticlesGlobal();
    require(baseline > 0, "remap: some particles were populated");
    require(pc.getNumberOfOutcastParticlesGlobal() == 0, "remap: no outcasts after canonicalize");
    checkOwnership(pc, dx, probLo);

    // ---- move everybody to fresh random in-domain positions, remap, re-check ----
    scatterPositions(pc, rng, 0.5, static_cast<Real>(N) - 0.5);
    pc.remap();
    require(pc.getNumberOfValidParticlesGlobal() == baseline, "remap: count conserved after move");
    require(pc.getNumberOfOutcastParticlesGlobal() == 0, "remap: no outcasts after in-domain move");
    checkOwnership(pc, dx, probLo);

    // ---- move everybody off-domain, remap, expect all outcast ----
    scatterPositions(pc, rng, -100.0, -90.0);
    pc.remap();
    require(pc.getNumberOfValidParticlesGlobal() == 0, "remap: all particles left the domain");
    require(pc.getNumberOfOutcastParticlesGlobal() == baseline, "remap: off-domain particles counted as outcast");

    pout() << "All ParticleContainerSoA::remap checks passed (" << SpaceDim << "D, " << numProc()
           << " rank(s), baseline " << baseline << ")." << endl;
  }

  return finalize();
}
