/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   test_ParticleLayoutMacro.cpp
  @brief  Side-by-side comparison: explicit ParticleTraits vs CD_PARTICLE_LAYOUT.
  @author Robert Marskar
  @details Two identical 7-field Ito-like particles, one described with hand-written
           traits and one with the macro. Asserts (mostly at compile time) that the
           macro produces the SAME layout: same column count, same derived
           position/weight indices, same per-particle byte sizes, same column types.

           Build: g++ -std=c++14 -O2 -Wall test_ParticleLayoutMacro.cpp -o test_macro
*/

#include <cassert>
#include <iostream>
#include <type_traits>

#include "CD_ParticleLayoutMacro.H"

using proxy::ParticleSoA;
using proxy::Real;
using proxy::RealVect;

// ---------------------------------------------------------------------------
// Variant 1: explicit, hand-written traits (the "field list written twice" form).
// ---------------------------------------------------------------------------
struct ItoExplicit
{
  RealVect pos       = RealVect::Zero();
  Real     weight    = 0.0;
  RealVect velocity  = RealVect::Zero();
  RealVect oldPos    = RealVect::Zero();
  Real     mobility  = 0.0;
  Real     diffusion = 0.0;
  Real     energy    = 0.0;
};

namespace proxy {
  template <>
  struct ParticleTraits<ItoExplicit>
  {
    static constexpr auto columns     = std::make_tuple(&ItoExplicit::pos,
                                                    &ItoExplicit::weight,
                                                    &ItoExplicit::velocity,
                                                    &ItoExplicit::oldPos,
                                                    &ItoExplicit::mobility,
                                                    &ItoExplicit::diffusion,
                                                    &ItoExplicit::energy);
    static constexpr auto positionPtr = &ItoExplicit::pos;
    static constexpr auto weightPtr   = &ItoExplicit::weight;
  };
} // namespace proxy

// ---------------------------------------------------------------------------
// Variant 2: same particle, described by the macro (single member list).
// ---------------------------------------------------------------------------
struct ItoMacro
{
  RealVect pos       = RealVect::Zero();
  Real     weight    = 0.0;
  RealVect velocity  = RealVect::Zero();
  RealVect oldPos    = RealVect::Zero();
  Real     mobility  = 0.0;
  Real     diffusion = 0.0;
  Real     energy    = 0.0;
};

namespace proxy {
  CD_PARTICLE_LAYOUT(ItoMacro,
                     pos,    // position (required, 1st)
                     weight, // weight   (required, 2nd)
                     velocity,
                     oldPos,
                     mobility,
                     diffusion,
                     energy);
} // namespace proxy

// ---------------------------------------------------------------------------
// The macro must produce an identical layout to the explicit traits.
// ---------------------------------------------------------------------------
using SoaE = ParticleSoA<ItoExplicit>;
using SoaM = ParticleSoA<ItoMacro>;

static_assert(SoaE::NumColumns == SoaM::NumColumns, "same column count");
static_assert(SoaE::NumColumns == 7, "7 columns");

static_assert(SoaE::positionIndex == SoaM::positionIndex, "same derived position index");
static_assert(SoaE::weightIndex == SoaM::weightIndex, "same derived weight index");
static_assert(SoaM::positionIndex == 0, "position is column 0");
static_assert(SoaM::weightIndex == 1, "weight is column 1");

static_assert(SoaE::s_bytesPerParticle == SoaM::s_bytesPerParticle, "same MPI byte size");
static_assert(SoaM::s_bytesPerParticle == 4 * sizeof(Real) + 3 * sizeof(RealVect), "4 Real + 3 RealVect");

// Column types match position-by-position.
static_assert(std::is_same<SoaE::ColumnType<0>, SoaM::ColumnType<0>>::value, "col0");
static_assert(std::is_same<SoaE::ColumnType<4>, SoaM::ColumnType<4>>::value, "col4");

// A field selected by member pointer lands on the same derived column in both.
static_assert(SoaE::columnIndex(&ItoExplicit::energy) == SoaM::columnIndex(&ItoMacro::energy),
              "energy derives to the same index");

int
main()
{
  // Runtime smoke test: both behave identically.
  SoaM     soa;
  ItoMacro p;
  p.pos    = RealVect(2.0);
  p.weight = 4.0;
  p.energy = 9.0;
  soa.append(p);

  const ItoMacro q = soa.gather(0);
  assert(q.pos[0] == 2.0);
  assert(q.weight == 4.0);
  assert(q.energy == 9.0);
  assert(soa.weight(0) == 4.0);

  std::cout << "Macro vs explicit traits: identical layout (" << SoaM::NumColumns << " columns, "
            << SoaM::s_bytesPerParticle << " B/particle).\n";
  std::cout << "test_ParticleLayoutMacro: OK\n";
  return 0;
}
