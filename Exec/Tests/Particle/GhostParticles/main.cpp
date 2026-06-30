/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   main.cpp
  @brief  Self-checking test of ParticleContainer::fillGhostParticles / clearGhostParticles
  @author Robert Marskar
  @details Builds a two-level AMR hierarchy by hand, scatters one particle per finest cell, fills the
  ghost layer, and verifies (independently, by brute force over the global grids) that every patch
  received exactly the ghost particles whose cell lies in its grown box but outside its own box, with
  the correct GhostType (SameLevel/Coarse/Fine), the owner's rankID preserved, and the exact position.
  Then clears the ghosts and checks the valid count is unchanged. Works in serial and under MPI.
*/

// Std includes
#include <map>
#include <vector>

// Chombo includes
#include <Box.H>
#include <BoxIterator.H>
#include <BRMeshRefine.H>
#include <DisjointBoxLayout.H>
#include <LayoutIterator.H>
#include <LoadBalance.H>
#include <ProblemDomain.H>
#include <RealVect.H>
#include <Vector.H>

// Our includes
#include <CD_Driver.H>
#include <CD_ParticleContainer.H>
#include <CD_ParticleSoA.H>
#include <CD_LevelTiles.H>

using namespace ChomboDischarge;

namespace {

  constexpr int  g_N        = 32;  // level-0 cells per direction
  constexpr int  g_minBlk   = 16;  // blocking factor (tile size)
  constexpr int  g_refRat   = 2;   // refinement ratio level 0 -> 1
  constexpr int  g_numGhost = 2;   // ghost layer width (cells)
  constexpr int  g_finest   = 1;   // finest level index
  const Real     g_L        = 1.0; // physical domain length
  const RealVect g_probLo   = RealVect::Zero;

  // Coarse-cell region (on level 0) that is refined to level 1.
  inline Box
  coarseFineRegion()
  {
    return Box(8 * IntVect::Unit, 23 * IntVect::Unit);
  }

  inline RealVect
  dxAt(const int a_lvl)
  {
    Real dx = g_L / g_N;
    for (int l = 0; l < a_lvl; l++) {
      dx /= g_refRat;
    }
    return dx * RealVect::Unit;
  }

  inline IntVect
  cellAt(const int a_lvl, const RealVect& a_pos)
  {
    const RealVect dx = dxAt(a_lvl);
    IntVect        iv;
    for (int dir = 0; dir < SpaceDim; dir++) {
      iv[dir] = static_cast<int>(std::floor((a_pos[dir] - g_probLo[dir]) / dx[dir]));
    }
    return iv;
  }

  struct BoxInfo
  {
    Box          box;
    unsigned int gridIndex;
    int          rank;
  };

  // Global box list (every rank knows all boxes -- the DBL is global metadata).
  std::vector<BoxInfo>
  globalBoxes(const DisjointBoxLayout& a_dbl)
  {
    std::vector<BoxInfo> out;
    for (LayoutIterator lit = a_dbl.layoutIterator(); lit.ok(); ++lit) {
      out.push_back(BoxInfo{a_dbl[lit()], a_dbl.index(lit()), a_dbl.procID(lit())});
    }
    return out;
  }

} // namespace

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  int nerr = 0;

  // ---------------------------------------------------------------------------
  // Build a two-level hierarchy by hand.
  // ---------------------------------------------------------------------------
  const Box           domain0(IntVect::Zero, (g_N - 1) * IntVect::Unit);
  const ProblemDomain pd0(domain0);
  const ProblemDomain pd1(refine(domain0, g_refRat));

  Vector<Box> boxes0;
  domainSplit(pd0, boxes0, g_minBlk, g_minBlk);
  Vector<int> procs0;
  LoadBalance(procs0, boxes0);
  const DisjointBoxLayout dbl0(boxes0, procs0, pd0);

  const Box   fineRegion = refine(coarseFineRegion(), g_refRat);
  Vector<Box> boxes1;
  domainSplit(fineRegion, boxes1, g_minBlk, g_minBlk);
  Vector<int> procs1;
  LoadBalance(procs1, boxes1);
  const DisjointBoxLayout dbl1(boxes1, procs1, pd1);

  Vector<DisjointBoxLayout> grids;
  grids.push_back(dbl0);
  grids.push_back(dbl1);
  Vector<ProblemDomain> domains;
  domains.push_back(pd0);
  domains.push_back(pd1);
  Vector<Real> dxs;
  dxs.push_back(g_L / g_N);
  dxs.push_back(g_L / g_N / g_refRat);
  Vector<int> refRats;
  refRats.push_back(g_refRat);
  refRats.push_back(g_refRat);

  Vector<RefCountedPtr<LevelTiles>> levelTiles;
  levelTiles.push_back(RefCountedPtr<LevelTiles>(new LevelTiles(dbl0, g_minBlk)));
  levelTiles.push_back(RefCountedPtr<LevelTiles>(new LevelTiles(dbl1, g_minBlk)));

  ParticleContainer<NoPayload> container;
  container.define(grids, domains, dxs, refRats, g_probLo, g_minBlk, g_numGhost, levelTiles, g_finest, "test", nullptr);

  // ---------------------------------------------------------------------------
  // Deterministic particle list (identical on every rank): one particle at the
  // centre of every finest cell. Level-1 cells in the refined region; level-0
  // cells everywhere else.
  // ---------------------------------------------------------------------------
  struct Part
  {
    ParticleID id;
    RealVect   pos;
  };
  std::vector<Part> allParts;
  {
    ParticleID id = 0;
    for (BoxIterator bit(fineRegion); bit.ok(); ++bit) {
      const RealVect pos = g_probLo + (RealVect(bit()) + 0.5 * RealVect::Unit) * dxAt(1);
      allParts.push_back(Part{id++, pos});
    }
    const Box cfCoarse = coarseFineRegion();
    for (BoxIterator bit(domain0); bit.ok(); ++bit) {
      if (!cfCoarse.contains(bit())) {
        const RealVect pos = g_probLo + (RealVect(bit()) + 0.5 * RealVect::Unit) * dxAt(0);
        allParts.push_back(Part{id++, pos});
      }
    }
  }

  // Stage the full list on rank 0 only (addParticlesDestructive remaps to owners).
  ParticleSoA<NoPayload> staging;
  if (procID() == 0) {
    for (const Part& p : allParts) {
      staging.append(p.pos, 1.0);
      staging.particleID(staging.size() - 1) = p.id;
    }
  }
  container.addParticlesDestructive(staging);

  const unsigned long long nValidBefore = container.getNumberOfValidParticlesGlobal();
  if (nValidBefore != allParts.size()) {
    pout() << "FAIL: valid count after add = " << nValidBefore << " expected " << allParts.size() << endl;
    nerr++;
  }

  // ---------------------------------------------------------------------------
  // Reference expectation: brute force over the global grids. owner(pos) = finest
  // box containing the cell. A particle is a ghost of every patch (within +/-1
  // level) whose grown box contains its cell but whose own box does not.
  // ---------------------------------------------------------------------------
  std::vector<std::vector<BoxInfo>> allBoxes = {globalBoxes(dbl0), globalBoxes(dbl1)};

  auto owner = [&](const RealVect& a_pos, int& a_lvl, int& a_rank) {
    for (int lvl = g_finest; lvl >= 0; lvl--) {
      const IntVect iv = cellAt(lvl, a_pos);
      for (const BoxInfo& bi : allBoxes[lvl]) {
        if (bi.box.contains(iv)) {
          a_lvl  = lvl;
          a_rank = bi.rank;
          return;
        }
      }
    }
    a_lvl  = -1;
    a_rank = -1;
  };

  struct Exp
  {
    GhostType gt;
    int       ownerRank;
    RealVect  pos;
  };
  std::map<std::pair<int, unsigned int>, std::map<ParticleID, Exp>> expected;

  for (const Part& p : allParts) {
    int ls    = -1;
    int orank = -1;
    owner(p.pos, ls, orank);
    if (ls < 0) {
      continue;
    }

    for (int lt = ls - 1; lt <= ls + 1; lt++) {
      if (lt < 0 || lt > g_finest) {
        continue;
      }

      const GhostType gt     = (ls == lt) ? GhostType::SameLevel : (ls < lt ? GhostType::Coarse : GhostType::Fine);
      const IntVect   cellLt = cellAt(lt, p.pos);
      const Box       domBox = domains[lt].domainBox();

      for (const BoxInfo& bi : allBoxes[lt]) {
        Box gb = grow(bi.box, g_numGhost);
        gb &= domBox;
        if (gb.contains(cellLt) && !bi.box.contains(cellLt)) {
          expected[{lt, bi.gridIndex}][p.id] = Exp{gt, orank, p.pos};
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Fill ghosts, then compare each LOCAL patch's ghosts to the reference.
  // ---------------------------------------------------------------------------
  container.fillGhostParticles();

  unsigned long long localGhosts = 0;
  for (int lt = 0; lt <= g_finest; lt++) {
    const DisjointBoxLayout& dbl = grids[lt];
    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const unsigned int gidx = dbl.index(dit());
      const auto&        leaf = container[lt][dit()];

      std::map<ParticleID, Exp> actual;
      for (std::size_t i = 0; i < leaf.size(); i++) {
        if (leaf.isGhost(i)) {
          localGhosts++;
          actual[leaf.particleID(i)] = Exp{leaf.ghost(i), leaf.rankID(i), leaf.position(i)};
        }
      }

      const auto&               eit = expected.find({lt, gidx});
      std::map<ParticleID, Exp> exp = (eit == expected.end()) ? std::map<ParticleID, Exp>{} : eit->second;

      if (actual.size() != exp.size()) {
        pout() << "FAIL: lvl " << lt << " grid " << gidx << " has " << actual.size() << " ghosts, expected "
               << exp.size() << endl;
        nerr++;
      }

      for (const auto& kv : exp) {
        const auto ait = actual.find(kv.first);
        if (ait == actual.end()) {
          pout() << "FAIL: lvl " << lt << " grid " << gidx << " missing ghost id " << kv.first << endl;
          nerr++;
          continue;
        }
        if (ait->second.gt != kv.second.gt) {
          pout() << "FAIL: lvl " << lt << " grid " << gidx << " id " << kv.first << " wrong GhostType" << endl;
          nerr++;
        }
        if (ait->second.ownerRank != kv.second.ownerRank) {
          pout() << "FAIL: lvl " << lt << " grid " << gidx << " id " << kv.first << " rankID " << ait->second.ownerRank
                 << " expected owner " << kv.second.ownerRank << endl;
          nerr++;
        }
        if ((ait->second.pos - kv.second.pos).vectorLength() > 1.E-12) {
          pout() << "FAIL: lvl " << lt << " grid " << gidx << " id " << kv.first << " position mismatch" << endl;
          nerr++;
        }
      }
    }
  }

  // Valid count must be unaffected by the ghost halo.
  if (container.getNumberOfValidParticlesGlobal() != allParts.size()) {
    pout() << "FAIL: valid count changed by fillGhostParticles" << endl;
    nerr++;
  }

  const unsigned long long nGhostGlobal = ParallelOps::sum(localGhosts);
  if (nGhostGlobal == 0) {
    pout() << "FAIL: no ghosts were created" << endl;
    nerr++;
  }

  // Also verify at least one ghost of each category exists (same-level, coarse, fine).
  {
    unsigned long long byType[4] = {0, 0, 0, 0};
    for (int lt = 0; lt <= g_finest; lt++) {
      for (DataIterator dit(grids[lt]); dit.ok(); ++dit) {
        const auto& leaf = container[lt][dit()];
        for (std::size_t i = 0; i < leaf.size(); i++) {
          if (leaf.isGhost(i)) {
            byType[static_cast<int>(leaf.ghost(i))]++;
          }
        }
      }
    }
    for (int t = 1; t <= 3; t++) {
      if (ParallelOps::sum(byType[t]) == 0) {
        pout() << "FAIL: no ghosts of GhostType " << t << endl;
        nerr++;
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Clear ghosts: none must remain, valid count unchanged.
  // ---------------------------------------------------------------------------
  container.clearGhostParticles();
  unsigned long long remaining = 0;
  for (int lt = 0; lt <= g_finest; lt++) {
    for (DataIterator dit(grids[lt]); dit.ok(); ++dit) {
      const auto& leaf = container[lt][dit()];
      for (std::size_t i = 0; i < leaf.size(); i++) {
        if (leaf.isGhost(i)) {
          remaining++;
        }
      }
    }
  }
  if (ParallelOps::sum(remaining) != 0) {
    pout() << "FAIL: " << ParallelOps::sum(remaining) << " ghosts remain after clearGhostParticles" << endl;
    nerr++;
  }
  if (container.getNumberOfValidParticlesGlobal() != allParts.size()) {
    pout() << "FAIL: valid count changed by clearGhostParticles" << endl;
    nerr++;
  }

  const int totalErr = static_cast<int>(ParallelOps::sum(static_cast<unsigned long long>(nerr)));
  if (procID() == 0) {
    if (totalErr == 0) {
      std::cout << "GhostParticles test: PASS (" << nGhostGlobal << " ghosts filled across the hierarchy)" << std::endl;
    }
    else {
      std::cout << "GhostParticles test: FAIL (" << totalErr << " errors)" << std::endl;
    }
  }

  ChomboDischarge::finalize();

  return (totalErr == 0) ? 0 : 1;
}
