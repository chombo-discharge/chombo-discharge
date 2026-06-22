/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/*
  AoS-vs-SoA parity harness.

  Samples the SAME N particles, injects them into BOTH the production ParticleContainer<PointParticle>
  (AoS, List<P>) and the new ParticleContainerSoA<> (SoA), remaps both, then:

    Check A (particle identity): per box/level, the two particle SETS are BIT-FOR-BIT identical
            (positions + weights), after a canonical sort. This proves the SoA remap routes the exact
            same particles to the exact same patches/levels/ranks as production.

    Check B (deposition): deposit the weight onto the AMR mesh with both stacks
            (EBAMRParticleMesh vs EBAMRParticleMeshSoA) and compare the meshes. Because the two
            remaps yield different in-box orderings, floating-point accumulation order differs, so the
            meshes are NOT bit-identical -- they agree to ROUNDOFF. We assert max|diff| is at the
            precision/roundoff level (<= 1e-12 * total weight) and report it.

  Deposition uses NGP and CIC (TSC is excluded -- EBParticleMeshSoA carries the partition-of-unity
  fix that production TSC lacks). Real precision is the build default (double), so AoS and SoA store
  positions/weights identically. Run single-rank AND under mpirun -np 2/4.
*/

// Std includes
#include <algorithm>
#include <cmath>
#include <vector>

// Chombo includes
#include <ProblemDomain.H>
#include <DisjointBoxLayout.H>
#include <DataIterator.H>
#include <BoxIterator.H>
#include <LevelData.H>
#include <BaseFab.H>
#include <EBCellFAB.H>
#include <EBCellFactory.H>
#include <EBLevelGrid.H>
#include <EBIndexSpace.H>
#include <AllRegularService.H>
#include <List.H>
#include <RealVect.H>
#include <IntVect.H>
#include <SPMD.H>
#include <MayDay.H>
#include <parstream.H>

// Our includes
#include <CD_Initialize.H>
#include <CD_LevelTiles.H>
#include <CD_PointParticle.H>
#include <CD_ParticleContainer.H>
#include <CD_EBAMRParticleMesh.H>
#include <CD_ParticleSoA.H>
#include <CD_ParticleContainerSoA.H>
#include <CD_EBAMRParticleMeshSoA.H>

using namespace ChomboDischarge;

namespace {

  using PCSoA = ParticleContainerSoA<>; // NoPayload: position + weight

  void
  require(const bool a_ok, const std::string& a_what)
  {
    if (!a_ok) {
      MayDay::Abort(("TestRemapDepositParity: check failed -- " + a_what).c_str());
    }
  }

  DisjointBoxLayout
  tiledDBL(const Box& a_region, const ProblemDomain& a_domain, const int a_bf)
  {
    const IntVect lo = a_region.smallEnd();
    const IntVect hi = a_region.bigEnd();
    Vector<Box>   boxes;
    const Box     tileBox(IntVect::Zero, (hi - lo) / a_bf);
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

  /** @brief (position, weight) record with a canonical ordering for set comparison. */
  struct PW
  {
    RealVect x;
    Real     w;
    bool
    operator<(const PW& o) const
    {
      for (int dir = 0; dir < SpaceDim; dir++) {
        if (x[dir] != o.x[dir]) {
          return x[dir] < o.x[dir];
        }
      }
      return w < o.w;
    }
  };

} // namespace

int
main(int argc, char* argv[])
{
  initialize(argc, argv);

  {
    constexpr int  N      = 16;
    constexpr int  bf     = 8;
    constexpr int  refRat = 2;
    constexpr int  ghost  = 2;
    constexpr int  nPart  = 2000;
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

    // ---- all-regular EB hierarchy ----
    AllRegularService allReg;
    EBIndexSpace*     ebis = Chombo_EBIS::instance();
    ebis->define(domain1, probLo, dx1, allReg, refRat * N, -1);

    Vector<RefCountedPtr<EBLevelGrid>> eblgs(2);
    eblgs[0] = RefCountedPtr<EBLevelGrid>(new EBLevelGrid(grids[0], domain0, ghost, ebis));
    eblgs[1] = RefCountedPtr<EBLevelGrid>(new EBLevelGrid(grids[1], domain1, ghost, ebis));
    eblgs[0]->setMaxRefinementRatio(refRat);
    eblgs[1]->setMaxCoarseningRatio(refRat, ebis);

    // ---- production container (needs valid masks + level tiles) ----
    Vector<ParticleContainer<PointParticle>::ValidMask> validMask(2);
    Vector<RefCountedPtr<LevelTiles>>                   levelTiles(2);
    for (int lvl = 0; lvl <= 1; lvl++) {
      validMask[lvl] = ParticleContainer<PointParticle>::ValidMask(
        new LevelData<BaseFab<bool>>(grids[lvl], 1, IntVect::Zero));
      for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
        (*validMask[lvl])[dit()].setVal(true);
      }
      levelTiles[lvl] = RefCountedPtr<LevelTiles>(new LevelTiles(grids[lvl], bf));
    }

    ParticleContainer<PointParticle> pcAoS;
    pcAoS.define(grids, domains, dx, refRats, validMask, levelTiles, probLo, bf, 1, "testRealm");

    PCSoA pcSoA;
    pcSoA.define(grids, domains, dx, refRats, probLo, bf, 1, "testRealm");

    EBAMRParticleMesh amrAoS;
    amrAoS.define(eblgs, refRats, dx, probLo, ghost, 1);
    EBAMRParticleMeshSoA amrSoA;
    amrSoA.define(eblgs, refRats, dx, probLo, ghost, 1);

    // ---- the SAME N particles (deterministic on every rank), exactly representable weights ----
    std::vector<RealVect> pos(nPart);
    std::vector<Real>     wgt(nPart);
    Real                  wTotal = 0.0;
    {
      // A simple deterministic LCG so positions/weights are identical on all ranks (no <random>
      // implementation dependence). Positions span the interior so they distribute over both levels.
      unsigned long long s    = 88172645463325252ULL;
      auto               next = [&]() -> Real {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        return static_cast<Real>((s >> 11) & ((1ULL << 24) - 1)) / static_cast<Real>(1ULL << 24); // [0,1)
      };
      for (int i = 0; i < nPart; i++) {
        RealVect x;
        for (int dir = 0; dir < SpaceDim; dir++) {
          x[dir] = (0.5 + next() * (N - 1.0)); // physical [0.5, N-0.5)
        }
        pos[i] = x;
        wgt[i] = 1.0 + 0.25 * static_cast<Real>(i % 13); // exact in float and double
        wTotal += wgt[i];
      }
    }

    // ---- inject identically: rank 0 owns the global list; both stacks then route it ----
    {
      List<PointParticle> list;
      if (procID() == 0) {
        for (int i = 0; i < nPart; i++) {
          list.add(PointParticle(pos[i], wgt[i]));
        }
      }
      pcAoS.addParticles(list); // production pool->map->scatter (collective)
    }
    {
      if (procID() == 0) {
        DataIterator dit = grids[0].dataIterator();
        dit.reset();
        PCSoA::Leaf& leaf = pcSoA[0][dit()];
        for (int i = 0; i < nPart; i++) {
          leaf.append(pos[i], wgt[i], NoPayload{});
        }
      }
      pcSoA.remap(); // SoA pool->map->scatter (collective)
    }

    require(pcAoS.getNumberOfValidParticlesGlobal() == static_cast<unsigned long long>(nPart),
            "AoS kept all particles");
    require(pcSoA.getNumberOfValidParticlesGlobal() == static_cast<unsigned long long>(nPart),
            "SoA kept all particles");

    // ---- Check A: per box/level, the particle sets are bit-for-bit identical ----
    for (int lvl = 0; lvl <= 1; lvl++) {
      for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
        std::vector<PW> a;
        std::vector<PW> b;

        const List<PointParticle>& aosList = pcAoS[lvl][dit()].listItems();
        for (ListIterator<PointParticle> lit(aosList); lit.ok(); ++lit) {
          a.push_back(PW{lit().position(), lit().weight()});
        }
        const PCSoA::Leaf& soaLeaf = pcSoA[lvl][dit()];
        for (std::size_t i = 0; i < soaLeaf.size(); i++) {
          b.push_back(PW{soaLeaf.position(i), soaLeaf.weight(i)});
        }

        require(a.size() == b.size(), "same particle count per box");
        std::sort(a.begin(), a.end());
        std::sort(b.begin(), b.end());
        for (std::size_t i = 0; i < a.size(); i++) {
          bool same = (a[i].w == b[i].w);
          for (int dir = 0; dir < SpaceDim; dir++) {
            same = same && (a[i].x[dir] == b[i].x[dir]);
          }
          require(same, "particle is bit-for-bit identical in both containers");
        }
      }
    }

    // ---- Check B: deposition agrees to roundoff ----
    EBAMRCellData meshAoS;
    EBAMRCellData meshSoA;
    meshAoS.resize(2);
    meshSoA.resize(2);
    for (int lvl = 0; lvl <= 1; lvl++) {
      meshAoS[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(
        new LevelData<EBCellFAB>(grids[lvl], 1, ghost * IntVect::Unit, EBCellFactory(eblgs[lvl]->getEBISL())));
      meshSoA[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(
        new LevelData<EBCellFAB>(grids[lvl], 1, ghost * IntVect::Unit, EBCellFactory(eblgs[lvl]->getEBISL())));
    }

    auto compareMeshes = [&]() -> Real {
      Real maxAbs = 0.0;
      for (int lvl = 0; lvl <= 1; lvl++) {
        for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
          const FArrayBox& a = (*meshAoS[lvl])[dit()].getFArrayBox();
          const FArrayBox& b = (*meshSoA[lvl])[dit()].getFArrayBox();
          for (BoxIterator bit(grids[lvl][dit()]); bit.ok(); ++bit) {
            maxAbs = std::max(maxAbs, std::abs(a(bit(), 0) - b(bit(), 0)));
          }
        }
      }
#ifdef CH_MPI
      Real g = 0.0;
      MPI_Allreduce(&maxAbs, &g, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
      maxAbs = g;
#endif
      return maxAbs;
    };

    struct Case
    {
      DepositionType       dep;
      CoarseFineDeposition cf;
      const char*          name;
    };
    const Case cases[3] = {{DepositionType::NGP, CoarseFineDeposition::Interp, "NGP/Interp"},
                           {DepositionType::CIC, CoarseFineDeposition::Interp, "CIC/Interp"},
                           {DepositionType::CIC, CoarseFineDeposition::Halo, "CIC/Halo"}};

    const Real tol = 1.e-12 * wTotal;
    for (int c = 0; c < 3; c++) {
      amrAoS.deposit<PointParticle, const Real&, &PointParticle::weight>(meshAoS,
                                                                         pcAoS,
                                                                         cases[c].dep,
                                                                         cases[c].cf,
                                                                         false);
      amrSoA.depositWeight(meshSoA, pcSoA, cases[c].dep, cases[c].cf, false);

      const Real maxAbs = compareMeshes();
      require(maxAbs <= tol, std::string("deposition agrees to roundoff [") + cases[c].name + "]");
      pout() << "  [" << cases[c].name << "] max|AoS-SoA| = " << maxAbs << "  (tol " << tol << ", wTotal " << wTotal
             << ", rel " << (maxAbs / wTotal) << ")" << endl;
    }

    pout() << "All AoS-vs-SoA parity checks passed (" << SpaceDim << "D, " << numProc() << " rank(s), " << nPart
           << " particles)." << endl;
  }

  return finalize();
}
