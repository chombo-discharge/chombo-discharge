/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/*
  Stage-1 correctness test for ParticleContainerSoA (AMR storage/accessor layer).

  Builds a 2-level all-regular grid hierarchy BY HAND (no AmrMesh, no EBIS -- the stage-1 container
  is geometry-light), defines the container, and exercises: grid metadata, per-level/per-patch leaf
  access via operator[], local population (addParticlesLocal + direct append), valid-particle
  counts, leaf iteration/round-trip, and clear(). Run single-rank.
*/

// Std includes
#include <cmath>

// Chombo includes
#include <ProblemDomain.H>
#include <DisjointBoxLayout.H>
#include <DataIterator.H>
#include <RealVect.H>
#include <IntVect.H>
#include <MayDay.H>
#include <parstream.H>

// Our includes
#include <CD_Initialize.H>
#include <CD_ParticleSoA.H>
#include <CD_ParticleContainerSoA.H>

using namespace ChomboDischarge;

namespace {

  using PC = ParticleContainerSoA<>; // NoPayload: position + weight + metadata

  /** @brief Always-on check. */
  void
  require(const bool a_ok, const std::string& a_what)
  {
    if (!a_ok) {
      MayDay::Abort(("TestParticleContainerSoA: check failed -- " + a_what).c_str());
    }
  }

  /** @brief Single-box DisjointBoxLayout covering a_domBox on proc 0. */
  DisjointBoxLayout
  singleBoxDBL(const Box& a_domBox, const ProblemDomain& a_domain)
  {
    Vector<Box> boxes(1, a_domBox);
    Vector<int> procs(1, 0);
    return DisjointBoxLayout(boxes, procs, a_domain);
  }

  /** @brief Sum of weights over every valid particle owned by this rank. */
  Real
  localWeightSum(const PC& a_pc)
  {
    Real sum = 0.0;
    for (int lvl = 0; lvl <= a_pc.getFinestLevel(); lvl++) {
      const PC::LevelParticles& level = a_pc[lvl];
      for (DataIterator dit = a_pc.getGrids()[lvl].dataIterator(); dit.ok(); ++dit) {
        const PC::Leaf& leaf = level[dit()];
        for (std::size_t i = 0; i < leaf.size(); i++) {
          sum += leaf.weight(i);
        }
      }
    }
    return sum;
  }

} // namespace

int
main(int argc, char* argv[])
{
  initialize(argc, argv);

  {
    constexpr int  N      = 16;
    const Real     dx0    = 1.0;
    const Real     dx1    = 0.5; // refRatio 2
    const RealVect probLo = RealVect::Zero;

    // ---- two-level all-regular hierarchy, built by hand ----
    const Box           domBox0(IntVect::Zero, (N - 1) * IntVect::Unit);
    const Box           domBox1(IntVect::Zero, (2 * N - 1) * IntVect::Unit);
    const ProblemDomain domain0(domBox0);
    const ProblemDomain domain1(domBox1);

    Vector<DisjointBoxLayout> grids(2);
    grids[0] = singleBoxDBL(domBox0, domain0);
    grids[1] = singleBoxDBL(domBox1, domain1);

    Vector<ProblemDomain> domains(2);
    domains[0] = domain0;
    domains[1] = domain1;

    Vector<Real> dx(2);
    dx[0] = dx0;
    dx[1] = dx1;

    Vector<int> refRat(2);
    refRat[0] = 2;
    refRat[1] = 2; // unused above finest, kept for shape

    PC pc;
    pc.define(grids, domains, dx, refRat, probLo, N, 1, "testRealm"); // blockingFactor = N (single box per level)

    // ---- (1) grid metadata ----
    require(pc.isDefined(), "container is defined");
    require(pc.getFinestLevel() == 1, "finest level == 1");
    require(pc.getGrids().size() == 2, "two levels of grids");
    require(pc.getDx().size() == 2, "two levels of dx");
    require(pc.getDx()[0][0] == dx0 && pc.getDx()[1][0] == dx1, "dx per level");
    for (int dir = 0; dir < SpaceDim; dir++) {
      require(pc.getProbLo()[dir] == 0.0, "probLo == 0");
    }
    require(!pc.isOrganizedByCell(), "not cell-organized in stage-1");

    // ---- (2) mask/buffer holders allocated and empty ----
    require(pc.getMaskParticles().size() == 2, "mask holder has two levels");
    require(pc.getBufferParticles().size() == 2, "buffer holder has two levels");
    {
      unsigned long long maskCount = 0;
      for (int lvl = 0; lvl <= pc.getFinestLevel(); lvl++) {
        const PC::LevelParticles& m = *pc.getMaskParticles()[lvl];
        for (DataIterator dit = pc.getGrids()[lvl].dataIterator(); dit.ok(); ++dit) {
          maskCount += m[dit()].size();
        }
      }
      require(maskCount == 0, "mask particles empty in stage-1");
    }

    // ---- (3) populate via addParticlesLocal, count ----
    const int nLvl0   = 12;
    const int nLvl1   = 8;
    Real      wExpect = 0.0;
    for (int i = 0; i < nLvl0; i++) {
      const Real     w = 1.0 + i;
      const RealVect x = (2.0 + 0.5 * i) * RealVect::Unit; // interior on level 0 (physical < N)
      pc.addParticlesLocal(0, x, w, NoPayload{});
      wExpect += w;
    }
    for (int i = 0; i < nLvl1; i++) {
      const Real     w = 100.0 + i;
      const RealVect x = (1.0 + 0.3 * i) * RealVect::Unit; // interior on level 1 (physical < N)
      pc.addParticlesLocal(1, x, w, NoPayload{});
      wExpect += w;
    }
    require(pc.getNumberOfValidParticlesGlobal() == static_cast<unsigned long long>(nLvl0 + nLvl1),
            "global count after addParticlesLocal");

    // ---- (4) direct per-patch append via operator[] ----
    {
      DataIterator dit = pc.getGrids()[0].dataIterator();
      dit.reset();
      pc[0][dit()].append(RealVect(D_DECL(3.0, 3.0, 3.0)), 7.0, NoPayload{});
      wExpect += 7.0;
    }
    require(pc.getNumberOfValidParticlesGlobal() == static_cast<unsigned long long>(nLvl0 + nLvl1 + 1),
            "global count after direct append");

    // ---- (5) round-trip: summed weight matches what we deposited ----
    require(std::abs(localWeightSum(pc) - wExpect) < 1.e-12 * wExpect, "weight round-trips through the leaves");

    // ---- (6) clear keeps capacity, zeros the count ----
    std::size_t capBefore = 0;
    {
      DataIterator dit = pc.getGrids()[0].dataIterator();
      dit.reset();
      capBefore = pc[0][dit()].capacity();
    }
    require(capBefore > 0, "level-0 leaf grew capacity");
    pc.clearParticles();
    require(pc.getNumberOfValidParticlesGlobal() == 0, "count is zero after clear");
    {
      DataIterator dit = pc.getGrids()[0].dataIterator();
      dit.reset();
      require(pc[0][dit()].capacity() == capBefore, "clear keeps arena capacity");
    }

    pout() << "All ParticleContainerSoA checks passed (" << SpaceDim << "D)." << endl;
  }

  return finalize();
}
