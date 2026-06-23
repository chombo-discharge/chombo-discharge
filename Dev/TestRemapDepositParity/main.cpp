/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/*
  AoS-vs-SoA parity harness on a CUT-CELL, MULTI-LEVEL (AMR) geometry.

  Geometry: a sphere carved out of a 2-level grid via GeometryShop -> genuine cut cells on both
  levels (asserted: irregular VoF count > 0). The fine level covers a sub-region, so the coarse-fine
  machinery is exercised together with the embedded boundary.

  Samples the SAME N particles, injects them into BOTH the production ParticleContainer<ParityParticle>
  (AoS, List<P>) and ParticleContainerSoA<PhiPayload> (SoA), remaps both, then:

    Check A (particle identity): per box/level the particle SETS (position + weight) are bit-for-bit
            identical after a canonical sort -> identical remap routing.
    Check B (deposition): deposit weight with EBAMRParticleMesh vs EBAMRParticleMeshSoA and compare
            the mesh; asserted <= 1e-12*wTotal (in practice 0).
    Check C (interpolation): interpolate a known mesh field onto each particle's phi with both stacks
            and compare phi per particle; asserted <= 1e-12*max|phi| (in practice 0). Interpolation is
            per-particle (no accumulation), so equality is bit-exact through the cut-cell stencils.

  NGP and CIC only (TSC excluded -- EBParticleMeshSoA carries the partition-of-unity fix production
  TSC lacks). Real=double build, so AoS/SoA store positions/weights/phi identically. Run single-rank
  AND under mpirun -np 2/4.
*/

// Std includes
#include <algorithm>
#include <cmath>
#include <tuple>
#include <vector>

// Chombo includes
#include <ProblemDomain.H>
#include <DisjointBoxLayout.H>
#include <DataIterator.H>
#include <BoxIterator.H>
#include <LevelData.H>
#include <BaseFab.H>
#include <FArrayBox.H>
#include <EBCellFAB.H>
#include <EBCellFactory.H>
#include <EBLevelGrid.H>
#include <EBISBox.H>
#include <EBIndexSpace.H>
#include <GeometryShop.H>
#include <SphereIF.H>
#include <IntVectSet.H>
#include <List.H>
#include <RealVect.H>
#include <IntVect.H>
#include <SPMD.H>
#include <MayDay.H>
#include <parstream.H>

// Our includes
#include <CD_Initialize.H>
#include <CD_LevelTiles.H>
#include <CD_GenericParticle.H>
#include <CD_ParticleContainer.H>
#include <CD_EBAMRParticleMesh.H>
#include <CD_ParticleSoA.H>
#include <CD_ParticleContainerSoA.H>
#include <CD_EBAMRParticleMeshSoA.H>

using namespace ChomboDischarge;

namespace ChomboDischarge {

  /** @brief AoS particle: position + weight (scalar 0) + phi (scalar 1). */
  class ParityParticle : public GenericParticle<2, 0>
  {
  public:
    ParityParticle() = default;
    ParityParticle(const RealVect& a_x, const Real a_w)
    {
      this->position()         = a_x;
      this->template real<0>() = a_w;
      this->template real<1>() = 0.0;
    }
    inline Real&
    weight()
    {
      return this->template real<0>();
    }
    inline const Real&
    weight() const
    {
      return this->template real<0>();
    }
    inline Real&
    phi()
    {
      return this->template real<1>();
    }
    inline const Real&
    phi() const
    {
      return this->template real<1>();
    }
  };

  /** @brief SoA payload: a single scalar phi (weight is container-owned). */
  struct PhiPayload
  {
    Real phi = 0.0;
  };
  template <>
  struct ParticleTraits<PhiPayload>
  {
    static constexpr auto columns = std::make_tuple(&PhiPayload::phi);
  };

} // namespace ChomboDischarge

namespace {

  using PCSoA = ParticleContainerSoA<PhiPayload>;

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

  /** @brief (position, scalar) record with a canonical ordering. */
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

  Real
  globalMax(Real a_local)
  {
#ifdef CH_MPI
    Real g = 0.0;
    MPI_Allreduce(&a_local, &g, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
    return g;
#else
    return a_local;
#endif
  }

  unsigned long long
  globalSum(unsigned long long a_local)
  {
#ifdef CH_MPI
    unsigned long long g = 0;
    MPI_Allreduce(&a_local, &g, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, Chombo_MPI::comm);
    return g;
#else
    return a_local;
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

    // ---- CUT-CELL geometry: a sphere whose surface crosses both levels ----
    const RealVect sphereCenter = 0.5 * static_cast<Real>(N) * RealVect::Unit; // physical centre
    const Real     sphereRadius = 0.3 * static_cast<Real>(N);
    SphereIF       sphere(sphereRadius, sphereCenter, false);
    GeometryShop   workshop(sphere, 0, dx1 * RealVect::Unit);

    EBIndexSpace* ebis = Chombo_EBIS::instance();
    ebis->define(domain1, probLo, dx1, workshop, refRat * N, -1);

    Vector<RefCountedPtr<EBLevelGrid>> eblgs(2);
    eblgs[0] = RefCountedPtr<EBLevelGrid>(new EBLevelGrid(grids[0], domain0, ghost, ebis));
    eblgs[1] = RefCountedPtr<EBLevelGrid>(new EBLevelGrid(grids[1], domain1, ghost, ebis));
    eblgs[0]->setMaxRefinementRatio(refRat);
    eblgs[1]->setMaxCoarseningRatio(refRat, ebis);

    // ---- confirm the geometry actually has cut cells ----
    {
      unsigned long long nIrreg = 0;
      for (int lvl = 0; lvl <= 1; lvl++) {
        const EBISLayout& ebisl = eblgs[lvl]->getEBISL();
        for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
          nIrreg += static_cast<unsigned long long>(ebisl[dit()].getIrregIVS(grids[lvl][dit()]).numPts());
        }
      }
      require(globalSum(nIrreg) > 0, "geometry has cut cells");
    }

    // ---- containers ----
    Vector<ParticleContainer<ParityParticle>::ValidMask> validMask(2);
    Vector<RefCountedPtr<LevelTiles>>                    levelTiles(2);
    for (int lvl = 0; lvl <= 1; lvl++) {
      validMask[lvl] = ParticleContainer<ParityParticle>::ValidMask(
        new LevelData<BaseFab<bool>>(grids[lvl], 1, IntVect::Zero));
      for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
        (*validMask[lvl])[dit()].setVal(true);
      }
      levelTiles[lvl] = RefCountedPtr<LevelTiles>(new LevelTiles(grids[lvl], bf));
    }

    ParticleContainer<ParityParticle> pcAoS;
    pcAoS.define(grids, domains, dx, refRats, validMask, levelTiles, probLo, bf, 1, "testRealm");

    PCSoA pcSoA;
    pcSoA.define(grids, domains, dx, refRats, probLo, bf, 1, "testRealm");

    EBAMRParticleMesh amrAoS;
    amrAoS.define(eblgs, refRats, dx, probLo, ghost, 1);
    EBAMRParticleMeshSoA amrSoA;
    amrSoA.define(eblgs, refRats, dx, probLo, ghost, 1);

    // ---- the SAME N particles (deterministic on every rank) ----
    std::vector<RealVect> pos(nPart);
    std::vector<Real>     wgt(nPart);
    Real                  wTotal = 0.0;
    {
      unsigned long long s    = 88172645463325252ULL;
      auto               next = [&]() -> Real {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        return static_cast<Real>((s >> 11) & ((1ULL << 24) - 1)) / static_cast<Real>(1ULL << 24);
      };
      for (int i = 0; i < nPart; i++) {
        RealVect x;
        for (int dir = 0; dir < SpaceDim; dir++) {
          x[dir] = 0.5 + next() * (N - 1.0);
        }
        pos[i] = x;
        wgt[i] = 1.0 + 0.25 * static_cast<Real>(i % 13);
        wTotal += wgt[i];
      }
    }

    // ---- inject identically: rank 0 owns the global list; both stacks route it ----
    {
      List<ParityParticle> list;
      if (procID() == 0) {
        for (int i = 0; i < nPart; i++) {
          list.add(ParityParticle(pos[i], wgt[i]));
        }
      }
      pcAoS.addParticles(list);
    }
    {
      if (procID() == 0) {
        DataIterator dit = grids[0].dataIterator();
        dit.reset();
        PCSoA::Leaf& leaf = pcSoA[0][dit()];
        for (int i = 0; i < nPart; i++) {
          leaf.append(pos[i], wgt[i], PhiPayload{});
        }
      }
      pcSoA.remap();
    }

    require(pcAoS.getNumberOfValidParticlesGlobal() == static_cast<unsigned long long>(nPart), "AoS kept all");
    require(pcSoA.getNumberOfValidParticlesGlobal() == static_cast<unsigned long long>(nPart), "SoA kept all");

    // ---- per-particle field comparison (sorted by position) ----
    auto compareField = [&](auto a_aosGet, auto a_soaGet, const std::string& a_what) {
      for (int lvl = 0; lvl <= 1; lvl++) {
        for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
          std::vector<PW> a;
          std::vector<PW> b;

          const List<ParityParticle>& L = pcAoS[lvl][dit()].listItems();
          for (ListIterator<ParityParticle> lit(L); lit.ok(); ++lit) {
            a.push_back(PW{lit().position(), a_aosGet(lit())});
          }
          const PCSoA::Leaf& leaf = pcSoA[lvl][dit()];
          for (std::size_t i = 0; i < leaf.size(); i++) {
            b.push_back(PW{leaf.position(i), a_soaGet(leaf, i)});
          }

          require(a.size() == b.size(), a_what + ": same count per box");
          std::sort(a.begin(), a.end());
          std::sort(b.begin(), b.end());
          for (std::size_t i = 0; i < a.size(); i++) {
            bool same = (a[i].w == b[i].w);
            for (int dir = 0; dir < SpaceDim; dir++) {
              same = same && (a[i].x[dir] == b[i].x[dir]);
            }
            require(same, a_what);
          }
        }
      }
    };

    // Check A: particle identity (position + weight) after remap.
    compareField(
      [](const ParityParticle& p) {
        return p.weight();
      },
      [](const PCSoA::Leaf& leaf, std::size_t i) {
        return leaf.weight(i);
      },
      "remap: particle identity");

    // ---- Check B: deposition agrees to roundoff ----
    EBAMRCellData meshAoS;
    EBAMRCellData meshSoA;
    EBAMRCellData field;
    meshAoS.resize(2);
    meshSoA.resize(2);
    field.resize(2);
    for (int lvl = 0; lvl <= 1; lvl++) {
      const EBCellFactory fact(eblgs[lvl]->getEBISL());
      meshAoS[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(
        new LevelData<EBCellFAB>(grids[lvl], 1, ghost * IntVect::Unit, fact));
      meshSoA[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(
        new LevelData<EBCellFAB>(grids[lvl], 1, ghost * IntVect::Unit, fact));
      field[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(
        new LevelData<EBCellFAB>(grids[lvl], 1, ghost * IntVect::Unit, fact));
    }

    auto meshMaxDiff = [&]() -> Real {
      Real m = 0.0;
      for (int lvl = 0; lvl <= 1; lvl++) {
        for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
          const FArrayBox& a = (*meshAoS[lvl])[dit()].getFArrayBox();
          const FArrayBox& b = (*meshSoA[lvl])[dit()].getFArrayBox();
          for (BoxIterator bit(grids[lvl][dit()]); bit.ok(); ++bit) {
            m = std::max(m, std::abs(a(bit(), 0) - b(bit(), 0)));
          }
        }
      }
      return globalMax(m);
    };

    const DepositionType       deps[2]    = {DepositionType::NGP, DepositionType::CIC};
    const char*                depName[2] = {"NGP", "CIC"};
    const CoarseFineDeposition cfs[2]     = {CoarseFineDeposition::Interp, CoarseFineDeposition::Halo};
    const char*                cfName[2]  = {"Interp", "Halo"};

    const Real depTol = 1.e-12 * wTotal;
    for (int d = 0; d < 2; d++) {
      for (int c = 0; c < 2; c++) {
        amrAoS.deposit<ParityParticle, const Real&, &ParityParticle::weight>(meshAoS, pcAoS, deps[d], cfs[c], false);
        amrSoA.depositWeight(meshSoA, pcSoA, deps[d], cfs[c], false);
        const Real m = meshMaxDiff();
        require(m <= depTol, std::string("deposit parity [") + depName[d] + "/" + cfName[c] + "]");
        pout() << "  deposit  [" << depName[d] << "/" << cfName[c] << "] max|AoS-SoA| = " << m << "  (tol " << depTol
               << ")" << endl;
      }
    }

    // ---- Check C: interpolation agrees to roundoff ----
    // Fill an identical mesh field (valid + ghost) on both levels; both stacks read the same bits.
    for (int lvl = 0; lvl <= 1; lvl++) {
      for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
        FArrayBox& fb = (*field[lvl])[dit()].getFArrayBox();
        for (BoxIterator bit(fb.box()); bit.ok(); ++bit) {
          const IntVect iv    = bit();
          Real          v     = 1.0;
          const Real    cf[3] = {0.5, 0.25, 0.125};
          for (int dir = 0; dir < SpaceDim; dir++) {
            v += cf[dir] * static_cast<Real>(iv[dir]);
          }
          fb(iv, 0) = v;
        }
      }
    }

    Real maxPhi = 0.0;
    for (int d = 0; d < 2; d++) {
      amrAoS.interpolate<ParityParticle, Real&, &ParityParticle::phi>(pcAoS, field, deps[d], false);
      amrSoA.interpolate<&PhiPayload::phi>(pcSoA, field, deps[d], false);

      // Track scale for the tolerance.
      for (int lvl = 0; lvl <= 1; lvl++) {
        for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
          const List<ParityParticle>& L = pcAoS[lvl][dit()].listItems();
          for (ListIterator<ParityParticle> lit(L); lit.ok(); ++lit) {
            maxPhi = std::max(maxPhi, std::abs(lit().phi()));
          }
        }
      }

      // Compare interpolated phi per particle (bit-exact comparison; tolerance reported for context).
      compareField(
        [](const ParityParticle& p) {
          return p.phi();
        },
        [](const PCSoA::Leaf& leaf, std::size_t i) {
          return leaf.column<&PhiPayload::phi>()[i];
        },
        std::string("interpolate parity [") + depName[d] + "]");
      pout() << "  interp   [" << depName[d] << "] phi matched bit-for-bit (max|phi| ~ " << globalMax(maxPhi) << ")"
             << endl;
    }

    // Check D: interpolate a mesh field onto the particle WEIGHT (EBAMRParticleMeshSoA::interpolateWeight
    // vs the AoS interpolate<...,&weight>). Done last because it overwrites the weights.
    for (int d = 0; d < 2; d++) {
      amrAoS.interpolate<ParityParticle, Real&, &ParityParticle::weight>(pcAoS, field, deps[d], false);
      amrSoA.interpolateWeight(pcSoA, field, deps[d], false);

      compareField(
        [](const ParityParticle& p) {
          return p.weight();
        },
        [](const PCSoA::Leaf& leaf, std::size_t i) {
          return leaf.weight(i);
        },
        std::string("interpolateWeight parity [") + depName[d] + "]");
      pout() << "  interpW  [" << depName[d] << "] weight matched bit-for-bit" << endl;
    }

    pout() << "All AoS-vs-SoA parity checks passed on cut cells (" << SpaceDim << "D, " << numProc() << " rank(s), "
           << nPart << " particles)." << endl;
  }

  return finalize();
}
