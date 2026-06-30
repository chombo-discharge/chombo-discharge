/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   main.cpp
  @brief  Integration test of ParticleContainer::fillGhostParticles on a real AMR hierarchy
  @author Robert Marskar
  @details Plugs into the McPhoto radiative-transfer app: the Driver builds a RodDielectric geometry
  and a real AMR hierarchy (refinement ratios 4 then 2 around the dielectric). Afterwards we allocate a
  Photon ParticleContainer on that hierarchy, sample particles uniformly in the domain box, fill the
  ghost layer, and verify that (1) every ghost particle is a member of the original sampled set with an
  identical position, and (2) every ghost's recorded owner rankID, owning level and GhostType are
  consistent with the realm's LevelTiles ownership map, and that it sits in the shell of the patch that
  holds it. Works in serial and under MPI.
*/

// Std includes
#include <random>
#include <unordered_map>
#include <vector>

// Chombo includes
#include <Box.H>
#include <DataIterator.H>
#include <DisjointBoxLayout.H>
#include <ProblemDomain.H>
#include <RealVect.H>

// Our includes
#include <CD_Driver.H>
#include <CD_McPhoto.H>
#include <CD_RodDielectric.H>
#include <CD_RadiativeTransferStepper.H>
#include <CD_ParticleContainer.H>
#include <CD_ParticleSoA.H>
#include <CD_Photon.H>
#include <CD_LevelTiles.H>
#include <CD_Realm.H>
#include <CD_ParallelOps.H>

using namespace ChomboDischarge;
using namespace Physics::RadiativeTransfer;

namespace {

  // Number of particles sampled uniformly in the domain box.
  constexpr int g_numParticles = 100000;

  // One sampled particle: a globally-unique id and a position.
  struct Part
  {
    ParticleID id;
    RealVect   pos;
  };

  // Deterministically sample g_numParticles uniformly in [probLo, probHi). Identical on every rank
  // (fixed seed, fixed call order), so the full original set is known everywhere for verification.
  std::vector<Part>
  sampleUniform(const RealVect& a_probLo, const RealVect& a_probHi)
  {
    std::mt19937_64                                              rng(0x9E3779B97F4A7C15ULL);
    std::array<std::uniform_real_distribution<double>, SpaceDim> dist;
    for (int dir = 0; dir < SpaceDim; dir++) {
      dist[dir] = std::uniform_real_distribution<double>(a_probLo[dir], a_probHi[dir]);
    }

    std::vector<Part> out;
    out.reserve(g_numParticles);
    for (int i = 0; i < g_numParticles; i++) {
      RealVect pos;
      for (int dir = 0; dir < SpaceDim; dir++) {
        pos[dir] = dist[dir](rng);
      }
      out.push_back(Part{static_cast<ParticleID>(i), pos});
    }
    return out;
  }

  // Verify fillGhostParticles on the AMR hierarchy currently held by a_amr. Returns the number of
  // verification errors found on this rank.
  int
  verifyGhostParticles(AmrMesh& a_amr)
  {
    const std::string                        realm      = Realm::Primal;
    const int                                finest     = a_amr.getFinestLevel();
    const int                                minBlk     = a_amr.getBlockingFactor();
    const int                                numGhost   = a_amr.getNumberOfGhostCells();
    const RealVect                           probLo     = a_amr.getProbLo();
    const RealVect                           probHi     = a_amr.getProbHi();
    const Vector<Real>                       dx         = a_amr.getDx();
    const Vector<ProblemDomain>              domains    = a_amr.getDomains();
    const Vector<DisjointBoxLayout>&         grids      = a_amr.getGrids(realm);
    const Vector<RefCountedPtr<LevelTiles>>& levelTiles = a_amr.getLevelTiles(realm);

    int nerr = 0;

    if (finest < 1) {
      pout() << "FAIL: expected at least 2 AMR levels, got finestLevel = " << finest << endl;
      return nerr + 1;
    }

    // ---- LevelTiles oracle: finest level / owning rank of the tile containing a position ----
    auto cellAt = [&](const int a_lvl, const RealVect& a_pos) {
      IntVect iv;
      for (int dir = 0; dir < SpaceDim; dir++) {
        iv[dir] = static_cast<int>(std::floor((a_pos[dir] - probLo[dir]) / dx[a_lvl]));
      }
      return iv;
    };
    auto tileAt = [&](const int a_lvl, const RealVect& a_pos) {
      IntVect iv;
      for (int dir = 0; dir < SpaceDim; dir++) {
        iv[dir] = static_cast<int>(std::floor((a_pos[dir] - probLo[dir]) / (minBlk * dx[a_lvl])));
      }
      return iv;
    };
    // Returns {level, rank} of the finest tile owning the position, or {-1,-1} if off-domain.
    auto ownerFromTiles = [&](const RealVect& a_pos, int& a_level, int& a_rank) {
      for (int lvl = finest; lvl >= 0; lvl--) {
        const IntVect tile = tileAt(lvl, a_pos);
        const auto&   my   = levelTiles[lvl]->getMyTiles();
        if (my.find(tile) != my.end()) {
          a_level = lvl;
          a_rank  = procID();
          return;
        }
        const auto& other = levelTiles[lvl]->getOtherTiles();
        const auto  oit   = other.find(tile);
        if (oit != other.end()) {
          a_level = lvl;
          a_rank  = static_cast<int>(oit->second.second);
          return;
        }
      }
      a_level = -1;
      a_rank  = -1;
    };

    // ---- sample particles uniformly, stage on rank 0, route to owners ----
    const std::vector<Part>                  orig = sampleUniform(probLo, probHi);
    std::unordered_map<ParticleID, RealVect> origPos;
    origPos.reserve(orig.size());
    for (const Part& p : orig) {
      origPos[p.id] = p.pos;
    }

    ParticleContainer<Photon> particles;
    a_amr.allocate(particles, realm);

    ParticleSoA<Photon> staging;
    if (procID() == 0) {
      for (const Part& p : orig) {
        staging.append(p.pos, 1.0);
        staging.particleID(staging.size() - 1) = p.id;
      }
    }
    particles.addParticlesDestructive(staging);

    // Every sampled position is inside the domain, so all should be retained as valid.
    unsigned long long expectValid = 0;
    for (const Part& p : orig) {
      int l = -1;
      int r = -1;
      ownerFromTiles(p.pos, l, r);
      if (l >= 0) {
        expectValid++;
      }
    }
    const unsigned long long nValid = particles.getNumberOfValidParticlesGlobal();
    if (nValid != expectValid) {
      pout() << "FAIL: valid count = " << nValid << " expected " << expectValid << endl;
      nerr++;
    }

    // ---- fill ghosts and verify every ghost in every local patch ----
    particles.fillGhostParticles();

    unsigned long long nGhost    = 0;
    unsigned long long byType[4] = {0, 0, 0, 0};

    for (int lt = 0; lt <= finest; lt++) {
      const DisjointBoxLayout& dbl = grids[lt];
      const Box                dom = domains[lt].domainBox();

      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        const unsigned int gidx   = dbl.index(dit());
        const Box          dstBox = dbl[dit()];
        Box                grown  = grow(dstBox, numGhost);
        grown &= dom;

        const auto& leaf = particles[lt][dit()];

        for (std::size_t i = 0; i < leaf.size(); i++) {
          if (!leaf.isGhost(i)) {
            continue;
          }

          nGhost++;
          byType[static_cast<int>(leaf.ghost(i))]++;

          const ParticleID id  = leaf.particleID(i);
          const RealVect   pos = leaf.position(i);

          // (1) the ghost must be a member of the original set, with an identical position.
          const auto oit = origPos.find(id);
          if (oit == origPos.end()) {
            pout() << "FAIL: ghost id " << id << " is not in the original particle set" << endl;
            nerr++;
            continue;
          }
          if ((oit->second - pos).vectorLength() > 1.E-12) {
            pout() << "FAIL: ghost id " << id << " position differs from the original" << endl;
            nerr++;
          }

          // (2) owner rank/level and GhostType must agree with the LevelTiles ownership of the position.
          int srcLvl  = -1;
          int srcRank = -1;
          ownerFromTiles(pos, srcLvl, srcRank);
          if (srcRank != leaf.rankID(i)) {
            pout() << "FAIL: ghost id " << id << " rankID " << leaf.rankID(i) << " disagrees with LevelTiles owner "
                   << srcRank << endl;
            nerr++;
          }
          const GhostType expectType = (srcLvl == lt)  ? GhostType::SameLevel
                                       : (srcLvl < lt) ? GhostType::Coarse
                                                       : GhostType::Fine;
          if (leaf.ghost(i) != expectType) {
            pout() << "FAIL: ghost id " << id << " GhostType inconsistent with source level " << srcLvl
                   << " (dest level " << lt << ")" << endl;
            nerr++;
          }

          // The ghost must lie in the shell of the patch that holds it (in the grown box, not the box),
          // and that patch must indeed be owned by this rank.
          const IntVect cellLt = cellAt(lt, pos);
          if (!grown.contains(cellLt) || dstBox.contains(cellLt)) {
            pout() << "FAIL: ghost id " << id << " is not in the ghost shell of its host patch" << endl;
            nerr++;
          }
          if (levelTiles[lt]->getMyGrids().find(gidx) == levelTiles[lt]->getMyGrids().end()) {
            pout() << "FAIL: host patch (lvl " << lt << ", grid " << gidx << ") is not owned by this rank" << endl;
            nerr++;
          }
        }
      }
    }

    // Global tallies + sanity that the test actually exercised all three ghost categories.
    const unsigned long long nGhostGlobal = ParallelOps::sum(nGhost);
    const unsigned long long nSame        = ParallelOps::sum(byType[static_cast<int>(GhostType::SameLevel)]);
    const unsigned long long nCoarse      = ParallelOps::sum(byType[static_cast<int>(GhostType::Coarse)]);
    const unsigned long long nFine        = ParallelOps::sum(byType[static_cast<int>(GhostType::Fine)]);

    if (nGhostGlobal == 0) {
      pout() << "FAIL: no ghost particles were created" << endl;
      nerr++;
    }
    if (nCoarse == 0) {
      pout() << "FAIL: no coarse-level ghosts were created (expected at the coarse-fine boundary)" << endl;
      nerr++;
    }
    if (nFine == 0) {
      pout() << "FAIL: no fine-level ghosts were created (expected at the coarse-fine boundary)" << endl;
      nerr++;
    }

    if (procID() == 0) {
      std::cout << "GhostParticlesAmr: finestLevel=" << finest << " ghosts(global)=" << nGhostGlobal
                << " [same=" << nSame << " coarse=" << nCoarse << " fine=" << nFine << "]" << std::endl;
    }

    return nerr;
  }

} // namespace

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  auto compgeom    = RefCountedPtr<ComputationalGeometry>(new RodDielectric());
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto tagger      = RefCountedPtr<CellTagger>(nullptr);
  auto timestepper = RefCountedPtr<RadiativeTransferStepper<McPhoto>>(new RadiativeTransferStepper<McPhoto>());
  auto engine      = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

  // Build geometry + the AMR hierarchy (max_steps = 0 in the inputs => setup only, no advance).
  engine->setupAndRun();

  // Verify ghost-particle filling on the hierarchy the Driver built.
  const int nerr     = verifyGhostParticles(*amr);
  const int totalErr = static_cast<int>(ParallelOps::sum(static_cast<unsigned long long>(nerr)));

  if (procID() == 0) {
    std::cout << (totalErr == 0 ? "GhostParticlesAmr test: PASS"
                                : "GhostParticlesAmr test: FAIL (" + std::to_string(totalErr) + " errors)")
              << std::endl;
  }

  ChomboDischarge::finalize();

  return (totalErr == 0) ? 0 : 1;
}
