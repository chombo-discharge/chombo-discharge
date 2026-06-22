/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/*
  Stage-2c correctness test for ParticleContainerSoA halo/mask + buffer.

  Builds a multi-box 2-level grid, populates it, then exercises the mask-particle filling against a
  hand-built per-cell mask (cells with iv[0] < midpoint marked true):
    - copyMaskParticles: mask holder gets exactly the valid particles whose cell is masked, and the
      valid particles are left untouched (count + global checksum unchanged);
    - transferMaskParticles: those same particles MOVE out of the valid holder into the mask holder
      (valid loses them, total is conserved, ids preserved).
  Also checks the buffer holder lives on grown grids (boxes grown by the refinement factor on the
  fine level) and starts empty. Run single-rank AND under mpirun -np 2/4.
*/

// Std includes
#include <cmath>
#include <random>

// Chombo includes
#include <ProblemDomain.H>
#include <DisjointBoxLayout.H>
#include <BoxLayout.H>
#include <DataIterator.H>
#include <BoxIterator.H>
#include <LevelData.H>
#include <BaseFab.H>
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
      MayDay::Abort(("TestParticleContainerHalo: check failed -- " + a_what).c_str());
    }
  }

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

  unsigned long long
  globalReduce(unsigned long long a_local)
  {
#ifdef CH_MPI
    unsigned long long global = 0;
    MPI_Allreduce(&a_local, &global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, Chombo_MPI::comm);
    return global;
#else
    return a_local;
#endif
  }

  /** @brief Global number of particles in a holder (iterated over a_grids). */
  unsigned long long
  countHolder(const PC::LevelParticles& a_level, const DisjointBoxLayout& a_grids)
  {
    unsigned long long n = 0;
    for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit) {
      n += a_level[dit()].size();
    }
    return n;
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

    Vector<DisjointBoxLayout> grids(2);
    grids[0] = tiledDBL(domBox0, domain0, bf);
    grids[1] = tiledDBL(Box((N / 2) * IntVect::Unit, ((3 * N / 2) - 1) * IntVect::Unit), domain1, bf);

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

    // ---- populate + canonicalize ----
    std::mt19937  rng(15485863 + procID());
    constexpr int perBox = 4;
    for (int lvl = 0; lvl <= 1; lvl++) {
      for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
        const Box box = grids[lvl][dit()];
        for (int k = 0; k < perBox; k++) {
          RealVect x;
          for (int dir = 0; dir < SpaceDim; dir++) {
            std::uniform_real_distribution<Real> u(box.smallEnd(dir) + 0.25, box.bigEnd(dir) + 0.75);
            x[dir] = u(rng) * dx[lvl];
          }
          pc.addParticlesLocal(lvl, x, 1.0, NoPayload{});
        }
      }
    }
    pc.remap();

    // ---- build a per-cell mask: cells with iv[0] < midpoint are true ----
    const int                                       mid0          = N / 2;          // level-0 midpoint cell
    const int                                       mid1          = refRat * N / 2; // level-1 midpoint cell
    const int                                       midOfLevel[2] = {mid0, mid1};
    Vector<RefCountedPtr<LevelData<BaseFab<bool>>>> masks(2);
    for (int lvl = 0; lvl <= 1; lvl++) {
      masks[lvl] = RefCountedPtr<LevelData<BaseFab<bool>>>(new LevelData<BaseFab<bool>>(grids[lvl], 1, IntVect::Zero));
      for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
        BaseFab<bool>& fab = (*masks[lvl])[dit()];
        fab.setVal(false);
        for (BoxIterator bit(fab.box()); bit.ok(); ++bit) {
          if (bit()[0] < midOfLevel[lvl]) {
            fab(bit(), 0) = true;
          }
        }
      }
    }

    // ---- expected: count valid particles whose cell is masked (per level), and total valid ----
    unsigned long long expectMaskedLocal = 0;
    unsigned long long validBeforeLocal  = 0;
    for (int lvl = 0; lvl <= 1; lvl++) {
      const PC::LevelParticles& level = pc[lvl];
      for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
        const PC::Leaf& leaf = level[dit()];
        validBeforeLocal += leaf.size();
        for (std::size_t i = 0; i < leaf.size(); i++) {
          if (cellOf(leaf.position(i), probLo, dx[lvl])[0] < midOfLevel[lvl]) {
            expectMaskedLocal++;
          }
        }
      }
    }
    const unsigned long long expectMasked = globalReduce(expectMaskedLocal);
    const unsigned long long validBefore  = globalReduce(validBeforeLocal);
    require(expectMasked > 0 && expectMasked < validBefore, "mask selects a strict, non-empty subset");

    // ---- (A) copyMaskParticles: mask holder = masked subset; valid holder untouched ----
    pc.copyMaskParticles(masks);
    {
      unsigned long long maskLocal  = 0;
      unsigned long long validLocal = 0;
      for (int lvl = 0; lvl <= 1; lvl++) {
        maskLocal += countHolder(*pc.getMaskParticles()[lvl], grids[lvl]);
        validLocal += countHolder(pc[lvl], grids[lvl]);
        // every mask particle must be in the masked region
        const PC::LevelParticles& ml = *pc.getMaskParticles()[lvl];
        for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
          const PC::Leaf& leaf = ml[dit()];
          for (std::size_t i = 0; i < leaf.size(); i++) {
            require(cellOf(leaf.position(i), probLo, dx[lvl])[0] < midOfLevel[lvl], "copy: mask particle is masked");
          }
        }
      }
      require(globalReduce(maskLocal) == expectMasked, "copy: mask holder has the masked subset");
      require(globalReduce(validLocal) == validBefore, "copy: valid holder is untouched");
    }

    // ---- (B) transferMaskParticles: masked subset MOVES out of the valid holder ----
    pc.clearMaskParticles();
    pc.transferMaskParticles(masks);
    {
      unsigned long long maskLocal  = 0;
      unsigned long long validLocal = 0;
      for (int lvl = 0; lvl <= 1; lvl++) {
        maskLocal += countHolder(*pc.getMaskParticles()[lvl], grids[lvl]);
        validLocal += countHolder(pc[lvl], grids[lvl]);
        // no remaining valid particle may be in the masked region
        const PC::LevelParticles& vl = pc[lvl];
        for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
          const PC::Leaf& leaf = vl[dit()];
          for (std::size_t i = 0; i < leaf.size(); i++) {
            require(!(cellOf(leaf.position(i), probLo, dx[lvl])[0] < midOfLevel[lvl]),
                    "transfer: no masked particle left in the valid holder");
          }
        }
      }
      require(globalReduce(maskLocal) == expectMasked, "transfer: mask holder has the masked subset");
      require(globalReduce(validLocal) == validBefore - expectMasked, "transfer: valid holder lost the subset");
    }

    // ---- (C) grown-grid buffer: level 1 boxes grown by the refinement factor, buffer empty ----
    {
      const Vector<Box> g0 = pc.getGrownGrids()[0].boxArray();
      const Vector<Box> b0 = grids[0].boxArray();
      for (int i = 0; i < g0.size(); i++) {
        require(g0[i] == b0[i], "buffer: level 0 grids are not grown");
      }
      const Vector<Box> g1 = pc.getGrownGrids()[1].boxArray();
      const Vector<Box> b1 = grids[1].boxArray();
      for (int i = 0; i < g1.size(); i++) {
        require(g1[i].contains(b1[i]) && g1[i] != b1[i], "buffer: level 1 grids grown beyond the valid box");
      }
      unsigned long long bufLocal = 0;
      for (int lvl = 0; lvl <= 1; lvl++) {
        for (DataIterator dit = pc.getGrownGrids()[lvl].dataIterator(); dit.ok(); ++dit) {
          bufLocal += (*pc.getBufferParticles()[lvl])[dit()].size();
        }
      }
      require(globalReduce(bufLocal) == 0, "buffer: starts empty");
    }

    pout() << "All ParticleContainerSoA halo/buffer checks passed (" << SpaceDim << "D, " << numProc()
           << " rank(s), masked " << expectMasked << "/" << validBefore << ")." << endl;
  }

  return finalize();
}
