/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/*
  Dev correctness test for EBParticleMeshSoA (Dev port of EBParticleMesh).

  Strategy: build one all-regular EBISBox and feed the SAME particle set, in the SAME order,
  through both the production EBParticleMesh (List<P>) and the staged EBParticleMeshSoA
  (ParticleSoA<P>). The per-particle EB kernels are byte-identical copies and the particle
  order is identical, so results agree BITWISE -- EXCEPT TSC deposit, where EBParticleMeshSoA
  deliberately fixes a production partition-of-unity bug (production TSC over-deposits ~2.33x
  per dimension; see PORTING_EBParticleMesh.md). Accordingly:
    - deposit (weight, phi): bitwise vs production for NGP/CIC; for TSC the SoA path is validated
      by partition-of-unity (sum(rho)*cellVol == sum(strength)) since production is wrong;
    - deposit (velocity, per-component): validated by partition-of-unity for all schemes;
    - interpolate (phi, velocity): bitwise vs production for NGP/CIC/TSC (interpolation is
      unaffected by the deposit bug).

  All-regular geometry never enters the cut-cell branches; those branches are verbatim copies,
  so this isolates the actual risk surface (the SoA loop marshalling + selector).
*/

// Std includes
#include <cmath>
#include <random>
#include <vector>

// Chombo includes
#include <ProblemDomain.H>
#include <DisjointBoxLayout.H>
#include <BoxIterator.H>
#include <EBIndexSpace.H>
#include <EBISLevel.H>
#include <EBISLayout.H>
#include <EBISBox.H>
#include <EBCellFAB.H>
#include <AllRegularService.H>
#include <List.H>
#include <RealVect.H>
#include <IntVect.H>
#include <MayDay.H>
#include <parstream.H>

// Our includes
#include <CD_Initialize.H>
#include <CD_DepositionType.H>
#include <CD_GenericParticle.H>
#include <CD_EBParticleMesh.H>
#include <CD_EBParticleMeshSoA.H>
#include <CD_ParticleSoA.H>

namespace ChomboDischarge {

  /** @brief Production-style particle: weight + a scalar (phi) + a vector (velocity). */
  class TestParticle : public GenericParticle<2, 1>
  {
  public:
    Real&
    weight()
    {
      return this->real<0>();
    }
    const Real&
    weight() const
    {
      return this->real<0>();
    }
    Real&
    phi()
    {
      return this->real<1>();
    }
    const Real&
    phi() const
    {
      return this->real<1>();
    }
    RealVect&
    velocity()
    {
      return this->vect<0>();
    }
    const RealVect&
    velocity() const
    {
      return this->vect<0>();
    }
  };

  /** @brief Matching SoA payload: phi scalar + per-component velocity (position/weight are owned). */
  struct TestPayload
  {
    Real phi = 0.0;
    Real D_DECL(vx = 0.0, vy = 0.0, vz = 0.0);
  };
  template <>
  struct ParticleTraits<TestPayload>
  {
    static constexpr auto columns = std::make_tuple(&TestPayload::phi,
                                                    D_DECL(&TestPayload::vx, &TestPayload::vy, &TestPayload::vz));
  };

} // namespace ChomboDischarge

using namespace ChomboDischarge;

namespace {

  using TestSoA = ParticleSoA<TestPayload>;

  /** @brief Always-on check. */
  void
  require(const bool a_ok, const std::string& a_what)
  {
    if (!a_ok) {
      MayDay::Abort(("TestEBParticleMesh: check failed -- " + a_what).c_str());
    }
  }

  /** @brief One particle's worth of test data. */
  struct PData
  {
    RealVect pos;
    Real     weight;
    Real     phi;
    RealVect vel;
  };

  /** @brief Maximum |a-b| over all cells/components of two FABs on a_box. */
  Real
  maxFabDiff(const EBCellFAB& a_a, const EBCellFAB& a_b, const Box& a_box)
  {
    const FArrayBox& fa = a_a.getFArrayBox();
    const FArrayBox& fb = a_b.getFArrayBox();
    Real             d  = 0.0;
    for (BoxIterator bit(a_box); bit.ok(); ++bit) {
      for (int c = 0; c < a_a.nComp(); c++) {
        d = std::max(d, std::abs(fa(bit(), c) - fb(bit(), c)));
      }
    }
    return d;
  }

  /** @brief Linear mesh field f_c(x) = 1 + c + sum_dir (dir+1)*x[dir], filled per component. */
  void
  fillLinearField(EBCellFAB& a_fab, const Box& a_box, const RealVect& a_probLo, const Real a_dx)
  {
    FArrayBox& fab = a_fab.getFArrayBox();
    for (BoxIterator bit(a_box); bit.ok(); ++bit) {
      const IntVect iv = bit();
      RealVect      x;
      for (int dir = 0; dir < SpaceDim; dir++) {
        x[dir] = a_probLo[dir] + (iv[dir] + 0.5) * a_dx;
      }
      for (int c = 0; c < a_fab.nComp(); c++) {
        Real val = 1.0 + c;
        for (int dir = 0; dir < SpaceDim; dir++) {
          val += (dir + 1) * x[dir];
        }
        fab(iv, c) = val;
      }
    }
  }

  /** @brief Build a List<TestParticle> and a TestSoA holding the same data in the same order. */
  void
  buildContainers(const std::vector<PData>& a_data, List<TestParticle>& a_list, TestSoA& a_soa)
  {
    a_list.clear();
    a_soa.clear();
    a_soa.reserve(a_data.size());
    for (const PData& d : a_data) {
      TestParticle p;
      p.position() = d.pos;
      p.weight()   = d.weight;
      p.phi()      = d.phi;
      p.velocity() = d.vel;
      a_list.append(p);

      TestPayload pl;
      pl.phi = d.phi;
      D_TERM(pl.vx = d.vel[0];, pl.vy = d.vel[1];, pl.vz = d.vel[2];);
      a_soa.append(d.pos, d.weight, pl);
    }
  }

} // namespace

int
main(int argc, char* argv[])
{
  initialize(argc, argv);

  {
    constexpr int       N = 16;
    const Box           domBox(IntVect::Zero, (N - 1) * IntVect::Unit);
    const ProblemDomain domain(domBox);
    const Real          dx     = 1.0;
    const RealVect      dxVect = dx * RealVect::Unit;
    const RealVect      probLo = RealVect::Zero;
    const Real          eps    = 0.0; // expect BITWISE agreement

    // ---- all-regular EBISBox over one box ----
    Vector<Box>       boxes(1, domBox);
    Vector<int>       procs(1, 0);
    DisjointBoxLayout dbl(boxes, procs, domain);

    AllRegularService allReg;
    EBIndexSpace*     ebis = Chombo_EBIS::instance();
    ebis->define(domain, probLo, dx, allReg, N, 0);

    EBISLayout ebisl;
    ebis->fillEBISLayout(ebisl, dbl, domain, 2);

    DataIterator dit = dbl.dataIterator();
    dit.reset();
    const DataIndex din     = dit();
    const EBISBox&  ebisbox = ebisl[din];

    EBParticleMesh    meshProd(domain, domBox, ebisbox, dxVect, probLo);
    EBParticleMeshSoA meshSoA(domain, domBox, ebisbox, dxVect, probLo);

    // ---- particle data: interior positions (>= 3 cells from the boundary) ----
    std::mt19937                         rng(12345);
    std::uniform_real_distribution<Real> uPos(3.0, N - 3.0);
    std::uniform_real_distribution<Real> uVal(0.5, 2.5);

    constexpr int      nPart = 500;
    std::vector<PData> data(nPart);
    for (int i = 0; i < nPart; i++) {
      RealVect pos, vel;
      for (int dir = 0; dir < SpaceDim; dir++) {
        pos[dir] = probLo[dir] + uPos(rng) * dx;
        vel[dir] = uVal(rng);
      }
      data[i] = PData{pos, uVal(rng), uVal(rng), vel};
    }

    const DepositionType types[3] = {DepositionType::NGP, DepositionType::CIC, DepositionType::TSC};
    const char*          names[3] = {"NGP", "CIC", "TSC"};

    Real cellVol = 1.0;
    for (int dir = 0; dir < SpaceDim; dir++) {
      cellVol *= dx;
    }

    Real totalWeight = 0.0;
    for (const PData& d : data) {
      totalWeight += d.weight;
    }

    for (int t = 0; t < 3; t++) {
      List<TestParticle> list;
      TestSoA            soa;
      buildContainers(data, list, soa);

      // ---- (1) deposit weight (scalar, mandatory column) ----
      {
        EBCellFAB rhoA(ebisbox, domBox, 1);
        EBCellFAB rhoB(ebisbox, domBox, 1);
        rhoA.setVal(0.0);
        rhoB.setVal(0.0);

        meshProd.deposit<TestParticle, const Real&, &TestParticle::weight>(rhoA, list, types[t], 1.0, false);
        meshSoA.depositWeight(rhoB, soa, types[t], 1.0, false);

        // EBParticleMeshSoA FIXES the production TSC deposit partition-of-unity bug, so for TSC the
        // two paths DELIBERATELY differ -- bitwise agreement is asserted only for NGP/CIC. The SoA
        // path is validated for all schemes by the mass-conservation check below.
        if (types[t] != DepositionType::TSC) {
          require(maxFabDiff(rhoA, rhoB, domBox) <= eps, std::string("deposit weight bitwise [") + names[t] + "]");
        }

        // Mass conservation: sum(rho)*cellVol == sum(weights) (interior particles -> full stencil inside).
        const FArrayBox& fb      = rhoB.getFArrayBox();
        Real             depMass = 0.0;
        for (BoxIterator bit(domBox); bit.ok(); ++bit) {
          depMass += fb(bit(), 0);
        }
        depMass *= cellVol;
        const Real relErr = std::abs(depMass - totalWeight) / totalWeight;
        pout() << "    [" << names[t] << "] weight mass rel.err=" << relErr << endl;
        require(relErr <= 1.e-9, std::string("deposit weight mass conservation [") + names[t] + "]");
      }

      // ---- (2) deposit a payload scalar (phi) ----
      {
        EBCellFAB rhoA(ebisbox, domBox, 1);
        EBCellFAB rhoB(ebisbox, domBox, 1);
        rhoA.setVal(0.0);
        rhoB.setVal(0.0);

        meshProd.deposit<TestParticle, const Real&, &TestParticle::phi>(rhoA, list, types[t], 1.0, false);
        meshSoA.deposit<&TestPayload::phi>(rhoB, soa, types[t], 1.0, false);

        // NGP/CIC: bitwise vs production. TSC: production is buggy, so validate the SoA path by
        // partition-of-unity (sum(rho)*cellVol == sum(phi)) instead.
        if (types[t] != DepositionType::TSC) {
          require(maxFabDiff(rhoA, rhoB, domBox) <= eps, std::string("deposit phi bitwise [") + names[t] + "]");
        }

        const FArrayBox& fb     = rhoB.getFArrayBox();
        Real             depPhi = 0.0;
        for (BoxIterator bit(domBox); bit.ok(); ++bit) {
          depPhi += fb(bit(), 0);
        }
        depPhi *= cellVol;
        Real expectPhi = 0.0;
        for (const PData& d : data) {
          expectPhi += d.phi;
        }
        require(std::abs(depPhi - expectPhi) <= 1.e-9 * expectPhi,
                std::string("deposit phi partition-of-unity [") + names[t] + "]");
      }

      // ---- (3) deposit the velocity vector (SoA per-component path; self-consistency) ----
      {
        EBCellFAB rhoVec(ebisbox, domBox, SpaceDim);
        rhoVec.setVal(0.0);
        meshSoA.deposit<D_DECL(&TestPayload::vx, &TestPayload::vy, &TestPayload::vz)>(rhoVec,
                                                                                      soa,
                                                                                      types[t],
                                                                                      1.0,
                                                                                      false);

        const FArrayBox& fb = rhoVec.getFArrayBox();
        for (int c = 0; c < SpaceDim; c++) {
          Real depMass = 0.0;
          for (BoxIterator bit(domBox); bit.ok(); ++bit) {
            depMass += fb(bit(), c);
          }
          depMass *= cellVol;
          Real expect = 0.0;
          for (const PData& d : data) {
            expect += d.vel[c];
          }
          require(std::abs(depMass - expect) <= 1.e-9 * std::abs(expect),
                  std::string("deposit velocity mass conservation [") + names[t] + "]");
        }
      }

      // ---- (4) interpolate a scalar field into phi ----
      {
        EBCellFAB field(ebisbox, domBox, 1);
        fillLinearField(field, domBox, probLo, dx);

        meshProd.interpolate<TestParticle, Real&, &TestParticle::phi>(list, field, types[t], false);
        meshSoA.interpolate<&TestPayload::phi>(soa, field, types[t], false);

        int         i   = 0;
        Real        d   = 0.0;
        const Real* phi = soa.column<&TestPayload::phi>();
        for (ListIterator<TestParticle> lit(list); lit.ok(); ++lit, ++i) {
          d = std::max(d, std::abs(lit().phi() - phi[i]));
        }
        require(d <= eps, std::string("interpolate phi bitwise [") + names[t] + "]");
      }

      // ---- (5) interpolate a vector field into velocity (per-component) ----
      {
        EBCellFAB field(ebisbox, domBox, SpaceDim);
        fillLinearField(field, domBox, probLo, dx);

        meshProd.interpolate<TestParticle, RealVect&, &TestParticle::velocity>(list, field, types[t], false);
        meshSoA.interpolate<D_DECL(&TestPayload::vx, &TestPayload::vy, &TestPayload::vz)>(soa, field, types[t], false);

        D_TERM(const Real* vx = soa.column<&TestPayload::vx>();, const Real* vy = soa.column<&TestPayload::vy>();
               , const Real*                                                 vz = soa.column<&TestPayload::vz>(););
        int  i = 0;
        Real d = 0.0;
        for (ListIterator<TestParticle> lit(list); lit.ok(); ++lit, ++i) {
          const RealVect v = lit().velocity();
          D_TERM(d = std::max(d, std::abs(v[0] - vx[i]));, d = std::max(d, std::abs(v[1] - vy[i]));
                 , d                                         = std::max(d, std::abs(v[2] - vz[i])););
        }
        require(d <= eps, std::string("interpolate velocity bitwise [") + names[t] + "]");
      }

      pout() << "  [" << names[t] << "] deposit(weight,phi,vel) + interpolate(phi,vel): OK" << endl;
    }

    // ---- explicit destination-component path (core API: deposit/interpolate at a chosen comp) ----
    {
      List<TestParticle> tmp;
      TestSoA            soa;
      buildContainers(data, tmp, soa);

      // Deposit weight into component 1 of a 2-component FAB; component 0 must stay untouched.
      EBCellFAB rho2(ebisbox, domBox, 2);
      rho2.setVal(0.0);
      meshSoA.depositWeight(rho2, 1, soa, DepositionType::CIC, 1.0, false);

      const FArrayBox& fb = rho2.getFArrayBox();
      Real             m0 = 0.0;
      Real             m1 = 0.0;
      for (BoxIterator bit(domBox); bit.ok(); ++bit) {
        m0 += std::abs(fb(bit(), 0));
        m1 += fb(bit(), 1);
      }
      require(m0 == 0.0, "explicit-comp: component 0 left untouched");
      require(std::abs(m1 * cellVol - totalWeight) <= 1.e-9 * totalWeight,
              "explicit-comp: weight lands in component 1");

      // Interpolate phi from component 1 of a 2-comp field; must equal interpolating from comp 0 of
      // a 1-comp field holding the same data.
      EBCellFAB f1(ebisbox, domBox, 1);
      EBCellFAB f2(ebisbox, domBox, 2);
      fillLinearField(f1, domBox, probLo, dx);
      f2.setVal(0.0);
      for (BoxIterator bit(domBox); bit.ok(); ++bit) {
        f2.getFArrayBox()(bit(), 1) = f1.getFArrayBox()(bit(), 0);
      }

      TestSoA soaA;
      TestSoA soaB;
      buildContainers(data, tmp, soaA);
      buildContainers(data, tmp, soaB);
      meshSoA.interpolate<&TestPayload::phi>(soaA, f1, DepositionType::CIC, false);    // comp 0 (convenience)
      meshSoA.interpolate<&TestPayload::phi>(soaB, f2, 1, DepositionType::CIC, false); // comp 1 (core)

      const Real* pa = soaA.column<&TestPayload::phi>();
      const Real* pb = soaB.column<&TestPayload::phi>();
      Real        d  = 0.0;
      for (std::size_t i = 0; i < soaA.size(); i++) {
        d = std::max(d, std::abs(pa[i] - pb[i]));
      }
      require(d == 0.0, "explicit-comp: interpolate from component 1 matches component 0");

      pout() << "  [explicit-comp] deposit->comp1 + interpolate<-comp1: OK" << endl;
    }

    pout() << "All EBParticleMeshSoA checks passed (" << nPart << " particles, " << SpaceDim << "D)." << endl;
  }

  return finalize();
}
