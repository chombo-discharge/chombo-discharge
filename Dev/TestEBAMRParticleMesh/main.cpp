/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/*
  Stage-A correctness test for EBAMRParticleMeshSoA (AMR deposit/interpolate over ParticleContainerSoA).

  Builds a 2-level all-regular hierarchy where the fine level covers a SUB-region of the domain (so
  the coarse-fine mass transfer is actually exercised), then checks:
    - interpolate: a constant mesh field interpolates to that constant on every particle (both levels);
    - deposit (CoarseFineDeposition::Interp): mass is conserved -- sum over valid coarse cells (those
      not covered by the fine patch) plus all fine cells, weighted by cell volume, equals the total
      deposited weight. This exercises exchange+EBAddOp and addFineGhostsToCoarse/addInvalidCoarseToFine.

  Run single-rank.
*/

// Std includes
#include <cmath>

// Chombo includes
#include <ProblemDomain.H>
#include <DisjointBoxLayout.H>
#include <DataIterator.H>
#include <BoxIterator.H>
#include <RealVect.H>
#include <IntVect.H>
#include <LevelData.H>
#include <EBCellFAB.H>
#include <EBCellFactory.H>
#include <EBLevelGrid.H>
#include <EBIndexSpace.H>
#include <EBISLevel.H>
#include <AllRegularService.H>
#include <MayDay.H>
#include <parstream.H>

// Our includes
#include <CD_Initialize.H>
#include <CD_DataOps.H>
#include <CD_DepositionType.H>
#include <CD_CoarseFineDeposition.H>
#include <CD_ParticleSoA.H>
#include <CD_ParticleContainerSoA.H>
#include <CD_EBAMRParticleMeshSoA.H>

using namespace ChomboDischarge;

namespace ChomboDischarge {

  /** @brief Test payload with a single scalar to interpolate into. */
  struct TestPayload
  {
    Real phi = 0.0;
  };
  template <>
  struct ParticleTraits<TestPayload>
  {
    static constexpr auto columns = std::make_tuple(&TestPayload::phi);
  };

} // namespace ChomboDischarge

namespace {

  using PC = ParticleContainerSoA<TestPayload>;

  void
  require(const bool a_ok, const std::string& a_what)
  {
    if (!a_ok) {
      MayDay::Abort(("TestEBAMRParticleMesh: check failed -- " + a_what).c_str());
    }
  }

  DisjointBoxLayout
  singleBoxDBL(const Box& a_box, const ProblemDomain& a_domain)
  {
    Vector<Box> boxes(1, a_box);
    Vector<int> procs(1, 0);
    return DisjointBoxLayout(boxes, procs, a_domain);
  }

} // namespace

int
main(int argc, char* argv[])
{
  initialize(argc, argv);

  {
    constexpr int  N      = 16;
    constexpr int  refRat = 2;
    constexpr int  ghost  = 2;
    const Real     dx0    = 1.0;
    const Real     dx1    = dx0 / refRat;
    const RealVect probLo = RealVect::Zero;

    // ---- domains: level 1 is the fully-refined domain; its PATCH covers only the centre ----
    const Box           domBox0(IntVect::Zero, (N - 1) * IntVect::Unit);
    const Box           domBox1(IntVect::Zero, (refRat * N - 1) * IntVect::Unit);
    const ProblemDomain domain0(domBox0);
    const ProblemDomain domain1(domBox1);

    // Fine patch covers coarse cells [N/4, 3N/4) -> fine cells [N/2, 3N/2).
    const Box fineBox((N / 2) * IntVect::Unit, ((3 * N / 2) - 1) * IntVect::Unit);

    Vector<DisjointBoxLayout> grids(2);
    grids[0] = singleBoxDBL(domBox0, domain0);
    grids[1] = singleBoxDBL(fineBox, domain1);

    Vector<ProblemDomain> domains(2);
    domains[0] = domain0;
    domains[1] = domain1;

    Vector<Real> dx(2);
    dx[0] = dx0;
    dx[1] = dx1;

    Vector<int> refRats(2);
    refRats[0] = refRat;
    refRats[1] = refRat;

    // ---- all-regular EBIS at the finest resolution, EBLevelGrids per level ----
    AllRegularService allReg;
    EBIndexSpace*     ebis = Chombo_EBIS::instance();
    ebis->define(domain1, probLo, dx1, allReg, refRat * N, -1);

    Vector<RefCountedPtr<EBLevelGrid>> eblgs(2);
    eblgs[0] = RefCountedPtr<EBLevelGrid>(new EBLevelGrid(grids[0], domain0, ghost, ebis));
    eblgs[1] = RefCountedPtr<EBLevelGrid>(new EBLevelGrid(grids[1], domain1, ghost, ebis));
    eblgs[0]->setMaxRefinementRatio(refRat);
    eblgs[1]->setMaxCoarseningRatio(refRat, ebis);

    // ---- AMR mesh (one component) ----
    EBAMRCellData mesh;
    mesh.resize(2);
    for (int lvl = 0; lvl <= 1; lvl++) {
      mesh[lvl] = RefCountedPtr<LevelData<EBCellFAB>>(
        new LevelData<EBCellFAB>(grids[lvl], 1, ghost * IntVect::Unit, EBCellFactory(eblgs[lvl]->getEBISL())));
    }

    // ---- container + AMR particle-mesh ----
    PC pc;
    pc.define(grids, domains, dx, refRats, probLo, N, 1, "testRealm"); // blockingFactor = N (single box per level)

    EBAMRParticleMeshSoA amrPM;
    amrPM.define(eblgs, refRats, dx, probLo, ghost, 1);

    // ---- particles: coarse particles in the UNCOVERED region, fine particles in the fine patch ----
    Real wExpect = 0.0;
    int  added   = 0;
    // coarse, valid region: cells [1,3] (physical [1,4)) -- left of the covered centre
    for (int i = 0; i < 6; i++) {
      const Real     w = 1.0 + i;
      const RealVect x = (1.5 + 0.3 * i) * RealVect::Unit;
      pc.addParticlesLocal(0, x, w, TestPayload{});
      wExpect += w;
      added++;
    }
    // fine patch: physical [5, 11) -> within coarse [N/4,3N/4)=[4,12)
    for (int i = 0; i < 10; i++) {
      const Real     w = 10.0 + i;
      const RealVect x = (5.0 + 0.5 * i) * RealVect::Unit;
      pc.addParticlesLocal(1, x, w, TestPayload{});
      wExpect += w;
      added++;
    }
    require(pc.getNumberOfValidParticlesGlobal() == static_cast<unsigned long long>(added),
            "all particles landed in a local patch");

    const DepositionType types[3] = {DepositionType::NGP, DepositionType::CIC, DepositionType::TSC};
    const char*          names[3] = {"NGP", "CIC", "TSC"};

    // ---- (A) interpolate a constant field -> every particle reads the constant ----
    for (int t = 0; t < 3; t++) {
      DataOps::setValue(mesh, 7.0);
      amrPM.interpolate<&TestPayload::phi>(pc, mesh, types[t], false);

      Real maxErr = 0.0;
      for (int lvl = 0; lvl <= pc.getFinestLevel(); lvl++) {
        const PC::LevelParticles& level = pc[lvl];
        for (DataIterator dit = pc.getGrids()[lvl].dataIterator(); dit.ok(); ++dit) {
          const PC::Leaf& leaf = level[dit()];
          const Real*     phi  = leaf.column<&TestPayload::phi>();
          for (std::size_t i = 0; i < leaf.size(); i++) {
            maxErr = std::max(maxErr, std::abs(phi[i] - 7.0));
          }
        }
      }
      require(maxErr < 1.e-12, std::string("interpolate constant field [") + names[t] + "]");
    }

    // ---- (B) deposit weight (Interp) -> mass conserved across the coarse-fine boundary ----
    const Box coveredCoarse = coarsen(fineBox, refRat); // coarse cells under the fine patch

    Real volCoarse = 1.0;
    Real volFine   = 1.0;
    for (int dir = 0; dir < SpaceDim; dir++) {
      volCoarse *= dx0;
      volFine *= dx1;
    }

    // Valid mass = coarse cells NOT under the fine patch + all fine cells.
    auto depositedMass = [&]() -> Real {
      Real depMass = 0.0;
      {
        DataIterator dit = grids[0].dataIterator();
        dit.reset();
        const FArrayBox& fb = (*mesh[0])[dit()].getFArrayBox();
        for (BoxIterator bit(domBox0); bit.ok(); ++bit) {
          if (!coveredCoarse.contains(bit())) {
            depMass += fb(bit(), 0) * volCoarse;
          }
        }
      }
      {
        DataIterator dit = grids[1].dataIterator();
        dit.reset();
        const FArrayBox& fb = (*mesh[1])[dit()].getFArrayBox();
        for (BoxIterator bit(fineBox); bit.ok(); ++bit) {
          depMass += fb(bit(), 0) * volFine;
        }
      }
      return depMass;
    };

    for (int t = 0; t < 3; t++) {
      amrPM.depositWeight(mesh, pc, types[t], CoarseFineDeposition::Interp, false);
      const Real depMass = depositedMass();
      // TSC deposit is fixed (partition-of-unity) in EBParticleMeshSoA, so all schemes conserve.
      require(std::abs(depMass - wExpect) <= 1.e-9 * wExpect,
              std::string("deposit Interp mass conservation [") + names[t] + "]");
      pout() << "  [" << names[t] << "] interpolate + deposit(Interp): OK (mass " << depMass << " / " << wExpect << ")"
             << endl;
    }

    // ---- (C) Halo coarse-fine strategy: re-deposits coarse halo particles at fine resolution.
    //      Total mass is still conserved across the boundary. ----
    for (int t = 0; t < 3; t++) {
      amrPM.depositWeight(mesh, pc, types[t], CoarseFineDeposition::Halo, false);
      const Real depMass = depositedMass();
      require(std::abs(depMass - wExpect) <= 1.e-9 * wExpect,
              std::string("deposit Halo mass conservation [") + names[t] + "]");
      pout() << "  [" << names[t] << "] deposit(Halo): OK (mass " << depMass << " / " << wExpect << ")" << endl;
    }

    // ---- (D) HaloNGP coarse-fine strategy: coarse halo particles deposit with NGP (no cloud spread
    //      over the boundary). Mass is conserved, and the valid particles are restored afterwards. ----
    const unsigned long long nBefore = pc.getNumberOfValidParticlesGlobal();
    for (int t = 0; t < 3; t++) {
      amrPM.depositWeight(mesh, pc, types[t], CoarseFineDeposition::HaloNGP, false);
      const Real depMass = depositedMass();
      require(std::abs(depMass - wExpect) <= 1.e-9 * wExpect,
              std::string("deposit HaloNGP mass conservation [") + names[t] + "]");
      require(pc.getNumberOfValidParticlesGlobal() == nBefore, "HaloNGP restores the valid particles");
      pout() << "  [" << names[t] << "] deposit(HaloNGP): OK (mass " << depMass << " / " << wExpect << ")" << endl;
    }

    // ---- (E) Transition coarse-fine strategy: coarse-side transition particles deposit at fine width
    //      onto the refined-coarse grid, distributed to both levels. Mass conserved; valid restored. ----
    for (int t = 0; t < 3; t++) {
      amrPM.depositWeight(mesh, pc, types[t], CoarseFineDeposition::Transition, false);
      const Real depMass = depositedMass();
      require(std::abs(depMass - wExpect) <= 1.e-9 * wExpect,
              std::string("deposit Transition mass conservation [") + names[t] + "]");
      require(pc.getNumberOfValidParticlesGlobal() == nBefore, "Transition restores the valid particles");
      pout() << "  [" << names[t] << "] deposit(Transition): OK (mass " << depMass << " / " << wExpect << ")" << endl;
    }

    pout() << "All EBAMRParticleMeshSoA checks passed (" << SpaceDim << "D)." << endl;
  }

  return finalize();
}
