/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/*
  Stage-2d correctness test for ParticleContainerSoA cell sorting.

  Populates a multi-box 2-level grid, then organizeParticlesByCell() and verifies that every leaf is
  cell-sorted: numCells == box cells, the CSR ranges partition the leaf, and every particle in the
  range of cell c actually lives in the cell whose Fortran (column-major) index is c. Also checks the
  organize flag, count conservation, and that adding a particle invalidates the cell organization.
  Run single-rank AND under mpirun -np 2/4.
*/

// Std includes
#include <cmath>
#include <cstddef>
#include <random>

// Chombo includes
#include <ProblemDomain.H>
#include <DisjointBoxLayout.H>
#include <DataIterator.H>
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
      MayDay::Abort(("TestParticleContainerCellSort: check failed -- " + a_what).c_str());
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

  /** @brief Fortran (column-major, dir 0 fastest) linear index of a cell within a box. */
  std::size_t
  fortranIndex(const Box& a_box, const IntVect& a_iv)
  {
    const IntVect lo  = a_box.smallEnd();
    const IntVect sz  = a_box.size();
    std::size_t   idx = 0;
    std::size_t   str = 1;
    for (int dir = 0; dir < SpaceDim; dir++) {
      idx += static_cast<std::size_t>(a_iv[dir] - lo[dir]) * str;
      str *= static_cast<std::size_t>(sz[dir]);
    }
    return idx;
  }

  /** @brief Cell of sorted particle i, computed from the leaf's double position columns. */
  IntVect
  cellOfLeaf(const PC::Leaf& a_leaf, const std::size_t a_i, const RealVect& a_probLo, const Real a_dx)
  {
    IntVect iv;
    for (int dir = 0; dir < SpaceDim; dir++) {
      iv[dir] = static_cast<int>(std::floor((a_leaf.positionColumn(dir)[a_i] - a_probLo[dir]) / a_dx));
    }
    return iv;
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

    std::mt19937  rng(32452843 + procID());
    constexpr int perBox = 20; // many per box so cells hold 0..several particles
    for (int lvl = 0; lvl <= 1; lvl++) {
      for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
        const Box box = grids[lvl][dit()];
        for (int k = 0; k < perBox; k++) {
          RealVect x;
          for (int dir = 0; dir < SpaceDim; dir++) {
            std::uniform_real_distribution<Real> u(box.smallEnd(dir) + 0.05, box.bigEnd(dir) + 0.95);
            x[dir] = u(rng) * dx[lvl];
          }
          pc.addParticlesLocal(lvl, x, 1.0, NoPayload{});
        }
      }
    }
    pc.remap();
    const unsigned long long baseline = pc.getNumberOfValidParticlesGlobal();
    require(baseline > 0, "particles populated");
    require(!pc.isOrganizedByCell(), "not cell-organized before sorting");

    // ---- organize by cell, verify the CSR structure on every leaf ----
    pc.organizeParticlesByCell();
    require(pc.isOrganizedByCell(), "organized by cell after sorting");

    for (int lvl = 0; lvl <= pc.getFinestLevel(); lvl++) {
      const PC::LevelParticles& level = pc[lvl];
      for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
        const Box       box  = grids[lvl][dit()];
        const PC::Leaf& leaf = level[dit()];

        require(leaf.isSorted(), "leaf is sorted");
        require(leaf.numCells() == static_cast<std::size_t>(box.numPts()), "numCells == box cell count");

        std::size_t total = 0;
        for (std::size_t c = 0; c < leaf.numCells(); c++) {
          const std::size_t begin = leaf.cellStart(c);
          const std::size_t end   = leaf.cellStart(c + 1);
          require(end >= begin, "CSR offsets are monotone");
          require(end - begin == leaf.particlesInCell(c), "particlesInCell matches CSR range");
          total += (end - begin);

          for (std::size_t i = begin; i < end; i++) {
            require(fortranIndex(box, cellOfLeaf(leaf, i, probLo, dx[lvl])) == c, "particle is in cell c");
          }
        }
        require(total == leaf.size(), "CSR ranges partition the leaf");
      }
    }

    // ---- count is unchanged; back-to-patch clears the flag ----
    require(pc.getNumberOfValidParticlesGlobal() == baseline, "count unchanged by cell sort");
    pc.organizeParticlesByPatch();
    require(!pc.isOrganizedByCell(), "by-patch clears the cell-organized flag");
    require(pc.getNumberOfValidParticlesGlobal() == baseline, "count unchanged by organizeParticlesByPatch");

    // ---- adding a particle invalidates the cell organization ----
    pc.organizeParticlesByCell();
    require(pc.isOrganizedByCell(), "re-organized by cell");
    pc.addParticlesLocal(0, (0.5) * RealVect::Unit, 1.0, NoPayload{}); // owned by rank 0 (cell 0 of box 0)
    require(!pc.isOrganizedByCell(), "adding a particle invalidates the cell sort");

    pout() << "All ParticleContainerSoA cell-sort checks passed (" << SpaceDim << "D, " << numProc()
           << " rank(s), baseline " << baseline << ")." << endl;
  }

  return finalize();
}
