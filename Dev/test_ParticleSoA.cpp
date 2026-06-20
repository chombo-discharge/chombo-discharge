/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   test_ParticleSoA.cpp
  @brief  Standalone driver exercising the ParticleSoA prototype against the
          common chombo-discharge particle usage patterns.
  @author Robert Marskar
  @details Build:  g++ -std=c++14 -O2 -Wall test_ParticleSoA.cpp -o test_ParticleSoA
*/

#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

#include "CD_ParticleSoA.H"

using proxy::ParticleSoA;
using proxy::ParticleTraits;
using proxy::Real;
using proxy::RealVect;
using proxy::SpaceDim;

// ---------------------------------------------------------------------------
// 1. A user-defined particle with an exotic extra field, exactly as requested.
// ---------------------------------------------------------------------------
struct MyParticle
{
  RealVect pos              = RealVect::Zero();
  Real     weight           = 0.0;
  Real     lightsaberDamage = 0.0;
};

namespace proxy {
  template <>
  struct ParticleTraits<MyParticle>
  {
    static constexpr auto columns = std::make_tuple(&MyParticle::pos,
                                                    &MyParticle::weight,
                                                    &MyParticle::lightsaberDamage);

    // Designate position/weight by MEMBER POINTER. ParticleSoA derives the indices,
    // so reordering `columns` above cannot misroute access.
    static constexpr auto positionPtr = &MyParticle::pos;
    static constexpr auto weightPtr   = &MyParticle::weight;
  };
} // namespace proxy

// ---------------------------------------------------------------------------
// 2. An ItoParticle-like type (position, weight, velocity, energy, oldPos, ...).
//    Mirrors GenericParticle<5,3> usage in the real code.
// ---------------------------------------------------------------------------
struct ItoLikeParticle
{
  RealVect pos       = RealVect::Zero();
  RealVect velocity  = RealVect::Zero();
  RealVect oldPos    = RealVect::Zero();
  Real     weight    = 0.0;
  Real     mobility  = 0.0;
  Real     diffusion = 0.0;
  Real     energy    = 0.0;
};

namespace proxy {
  template <>
  struct ParticleTraits<ItoLikeParticle>
  {
    static constexpr auto columns = std::make_tuple(&ItoLikeParticle::pos,
                                                    &ItoLikeParticle::velocity,
                                                    &ItoLikeParticle::oldPos,
                                                    &ItoLikeParticle::weight,
                                                    &ItoLikeParticle::mobility,
                                                    &ItoLikeParticle::diffusion,
                                                    &ItoLikeParticle::energy);

    static constexpr auto positionPtr = &ItoLikeParticle::pos;
    static constexpr auto weightPtr   = &ItoLikeParticle::weight;

    // HDF5 checkpoints only a SUBSET: position, weight, energy. velocity/oldPos/
    // mobility/diffusion are recomputed after restart, so they are not persisted.
    // Declared by MEMBER POINTER (no literal indices) -> follows column reorders.
    static constexpr auto h5Columns = std::make_tuple(&ItoLikeParticle::pos,
                                                      &ItoLikeParticle::weight,
                                                      &ItoLikeParticle::energy);
  };
} // namespace proxy

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------
static void
testBasic()
{
  ParticleSoA<MyParticle> soa;

  for (int i = 0; i < 10; i++) {
    MyParticle p;
    p.pos              = RealVect(static_cast<Real>(i));
    p.weight           = static_cast<Real>(2 * i);
    p.lightsaberDamage = static_cast<Real>(100 + i);
    soa.append(p);
  }

  assert(soa.size() == 10);

  // Mandatory accessors.
  assert(soa.position(3)[0] == 3.0);
  assert(soa.weight(3) == 6.0);

  // Whole-column access (the SoA win: contiguous, vectorizable).
  Real sumWeight = 0.0;
  for (const Real w : soa.weights()) {
    sumWeight += w;
  }
  assert(sumWeight == 90.0);

  // gather() round-trips the AoS view.
  const MyParticle g = soa.gather(4);
  assert(g.weight == 8.0);
  assert(g.lightsaberDamage == 104.0);

  // swap-and-pop removal.
  soa.remove(0);
  assert(soa.size() == 9);

  std::cout << "testBasic            : OK\n";
}

static void
testDepositLoop()
{
  // Mirrors EBParticleMesh::deposit: iterate particles, read position + a scalar
  // field, accumulate onto a (fake) mesh.
  ParticleSoA<ItoLikeParticle> soa;
  for (int i = 0; i < 1000; i++) {
    ItoLikeParticle p;
    p.pos    = RealVect(static_cast<Real>(i % 16));
    p.weight = 1.0;
    soa.append(p);
  }

  std::vector<Real> mesh(16, 0.0);

  const auto& pos = soa.positions();
  const auto& w   = soa.weights();

  for (std::size_t n = 0; n < soa.size(); n++) {
    const int cell = static_cast<int>(pos[n][0]);
    mesh[cell] += w[n];
  }

  Real total = 0.0;
  for (const Real m : mesh) {
    total += m;
  }
  assert(total == 1000.0);

  std::cout << "testDepositLoop      : OK\n";
}

static void
testInterpolateLoop()
{
  // Mirrors ItoSolver: interpolate a field onto particle velocity, then scale
  // velocity by mobility (an in-place per-particle update across columns).
  ParticleSoA<ItoLikeParticle> soa;
  for (int i = 0; i < 100; i++) {
    ItoLikeParticle p;
    p.pos      = RealVect(1.0);
    p.mobility = 2.0;
    soa.append(p);
  }

  auto&       vel = soa.column<&ItoLikeParticle::velocity>(); // single-token selector
  const auto& mob = soa.column<&ItoLikeParticle::mobility>();

  for (std::size_t n = 0; n < soa.size(); n++) {
    vel[n] = RealVect(3.0); // pretend this came from mesh interpolation
    vel[n] *= mob[n];
  }

  assert(vel[42][0] == 6.0);

  std::cout << "testInterpolateLoop  : OK\n";
}

static void
testLinearization()
{
  // Mirrors MPI remap: linearize a particle, ship the bytes, delinearize on the
  // far end into another container.
  ParticleSoA<ItoLikeParticle> src;
  ItoLikeParticle              p;
  p.pos      = RealVect(7.0);
  p.weight   = 3.5;
  p.energy   = 9.0;
  p.velocity = RealVect(1.0);
  src.append(p);

  std::vector<unsigned char> buffer(ParticleSoA<ItoLikeParticle>::s_bytesPerParticle);
  src.linearizeParticle(buffer.data(), 0);

  ParticleSoA<ItoLikeParticle> dst;
  dst.delinearizeAndAppend(buffer.data());

  assert(dst.size() == 1);
  const ItoLikeParticle q = dst.gather(0);
  assert(q.pos[0] == 7.0);
  assert(q.weight == 3.5);
  assert(q.energy == 9.0);
  assert(q.velocity[0] == 1.0);

  std::cout << "testLinearization    : OK  (" << ParticleSoA<ItoLikeParticle>::s_bytesPerParticle
            << " bytes/particle)\n";
}

static void
testFieldSelector()
{
  // Demonstrates the C++17 call-site ergonomics:
  //   deposit<&MyParticle::weight>(mesh, soa)
  // - field selected by a single member-pointer token (reorder-safe; index derived)
  // - field TYPE (Real here) deduced inside deposit(); caller never spells it
  ParticleSoA<MyParticle> soa;
  for (int i = 0; i < 100; i++) {
    MyParticle p;
    p.pos              = RealVect(static_cast<Real>(i % 8));
    p.weight           = 1.0;
    p.lightsaberDamage = 10.0;
    soa.append(p);
  }

  std::vector<Real> meshW(8, 0.0);
  std::vector<Real> meshD(8, 0.0);

  proxy::deposit<&MyParticle::weight>(meshW, soa);           // by member pointer, type deduced
  proxy::deposit<&MyParticle::lightsaberDamage>(meshD, soa); // different field, no other change

  Real wTot = 0.0;
  Real dTot = 0.0;
  for (int i = 0; i < 8; i++) {
    wTot += meshW[i];
    dTot += meshD[i];
  }
  assert(wTot == 100.0);
  assert(dTot == 1000.0);

  std::cout << "testFieldSelector    : OK\n";
}

static void
testUpdateNonMandatoryField()
{
  // A user iterates particles and updates a NON-mandatory field (lightsaberDamage).
  // Three supported styles, all reorder-safe (field selected by member pointer):
  using Soa = ParticleSoA<MyParticle>;

  Soa soa;
  for (int i = 0; i < 10; i++) {
    MyParticle p;
    p.lightsaberDamage = 10.0;
    soa.append(p);
  }

  // (1) Whole-column access -- RECOMMENDED. Contiguous, vectorizable, no copies.
  for (Real& d : soa.column<&MyParticle::lightsaberDamage>()) {
    d *= 2.0; // 10 -> 20
  }
  assert(soa.column<&MyParticle::lightsaberDamage>()[3] == 20.0);

  // (2) Per-particle element reference via get<&P::field>(i) -- when you need index n.
  for (std::size_t n = 0; n < soa.size(); n++) {
    soa.get<&MyParticle::lightsaberDamage>(n) += 5.0; // 20 -> 25
  }
  assert(soa.get<&MyParticle::lightsaberDamage>(3) == 25.0);

  // (3) AoS view: gather -> modify -> scatter. Most familiar, but copies the whole
  //     particle; fine for occasional per-particle work, avoid in hot loops.
  for (std::size_t n = 0; n < soa.size(); n++) {
    MyParticle p = soa.gather(n);
    p.lightsaberDamage += 100.0; // 25 -> 125
    soa.scatter(n, p);
  }
  assert(soa.gather(3).lightsaberDamage == 125.0);

  std::cout << "testUpdateField      : OK\n";
}

static void
testH5SubsetLinearization()
{
  // Mirrors GenericParticle's H5size/H5linearOut/H5linearIn: checkpoint only a
  // SUBSET of fields. Here ItoLikeParticle persists {pos, weight, energy} and
  // drops {velocity, oldPos, mobility, diffusion}.
  using Soa = ParticleSoA<ItoLikeParticle>;

  // Full MPI linearization keeps everything; H5 keeps only the subset.
  static_assert(Soa::s_bytesPerParticle == 4 * sizeof(Real) + 3 * sizeof(RealVect), "full = 4 Real + 3 RealVect");
  assert(Soa::h5BytesPerParticle() == sizeof(RealVect) + sizeof(Real) + sizeof(Real));

  Soa             src;
  ItoLikeParticle p;
  p.pos      = RealVect(7.0);
  p.weight   = 3.5;
  p.energy   = 9.0;
  p.velocity = RealVect(1.0); // NOT checkpointed
  p.mobility = 5.0;           // NOT checkpointed
  src.append(p);

  std::vector<unsigned char> buffer(Soa::h5BytesPerParticle());
  src.h5LinearizeParticle(buffer.data(), 0);

  Soa dst;
  dst.h5DelinearizeAndAppend(buffer.data());

  const ItoLikeParticle q = dst.gather(0);
  assert(q.pos[0] == 7.0);      // checkpointed -> preserved
  assert(q.weight == 3.5);      // checkpointed -> preserved
  assert(q.energy == 9.0);      // checkpointed -> preserved
  assert(q.velocity[0] == 0.0); // dropped -> default after restart
  assert(q.mobility == 0.0);    // dropped -> default after restart

  std::cout << "testH5Subset         : OK  (full " << Soa::s_bytesPerParticle << "B, H5 " << Soa::h5BytesPerParticle()
            << "B/particle)\n";
}

// Same fields as MyParticle, but the columns are registered in a DIFFERENT order
// (weight first). Nothing else changes - this proves indices are derived.
struct ReorderedParticle
{
  RealVect pos              = RealVect::Zero();
  Real     weight           = 0.0;
  Real     lightsaberDamage = 0.0;
};

namespace proxy {
  template <>
  struct ParticleTraits<ReorderedParticle>
  {
    // weight is column 0 here, position is column 1 (swapped vs MyParticle).
    static constexpr auto columns     = std::make_tuple(&ReorderedParticle::weight,
                                                    &ReorderedParticle::pos,
                                                    &ReorderedParticle::lightsaberDamage);
    static constexpr auto positionPtr = &ReorderedParticle::pos;
    static constexpr auto weightPtr   = &ReorderedParticle::weight;
  };
} // namespace proxy

static void
testReorderSafety()
{
  using A = ParticleSoA<MyParticle>;
  using B = ParticleSoA<ReorderedParticle>;

  // The physical column order differs...
  static_assert(A::weightIndex == 1, "MyParticle: weight is column 1");
  static_assert(A::positionIndex == 0, "MyParticle: position is column 0");
  static_assert(B::weightIndex == 0, "ReorderedParticle: weight is column 0");
  static_assert(B::positionIndex == 1, "ReorderedParticle: position is column 1");

  // ...yet name-based access via member pointer resolves correctly for both, with
  // no hand-edited index anywhere. A deposit-by-name still hits the weight column.
  B soa;
  for (int i = 0; i < 50; i++) {
    ReorderedParticle p;
    p.pos    = RealVect(static_cast<Real>(i % 4));
    p.weight = 2.0;
    soa.append(p);
  }
  assert(soa.weight(7) == 2.0); // mandatory accessor still correct after reorder

  std::vector<Real> mesh(4, 0.0);
  proxy::deposit<&ReorderedParticle::weight>(mesh, soa); // resolves to column 0 here

  Real tot = 0.0;
  for (const Real m : mesh) {
    tot += m;
  }
  assert(tot == 100.0);

  std::cout << "testReorderSafety    : OK\n";
}

int
main()
{
  std::cout << "SpaceDim = " << SpaceDim << "\n";
  testBasic();
  testDepositLoop();
  testInterpolateLoop();
  testLinearization();
  testFieldSelector();
  testReorderSafety();
  testUpdateNonMandatoryField();
  testH5SubsetLinearization();
  std::cout << "All ParticleSoA prototype tests passed.\n";
  return 0;
}
