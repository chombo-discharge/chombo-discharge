/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ParticleLayoutMacro.H
  @brief  Optional sugar: CD_PARTICLE_LAYOUT macro that generates a ParticleTraits
          specialization from a single member list.
  @author Robert Marskar
  @details This is the (b) half of the schema discussion - a macro that removes the
           "field list written twice" duplication of the explicit-traits approach.
           It expands to EXACTLY the same ParticleTraits a user would hand-write
           (same generated code, same performance, same reorder-safety via derived
           indices). The preprocessor plumbing is isolated here so the rest of the
           prototype stays macro-free.

           Usage (inside namespace proxy, with the particle struct already defined):

             CD_PARTICLE_LAYOUT(ItoParticle,
                                pos,      // 1st member = position (required)
                                weight,   // 2nd member = weight   (required)
                                velocity, oldPos, mobility, diffusion, energy);

           Supports up to 16 columns. Relies on the GNU ", ##__VA_ARGS__" extension
           (supported by GCC and Intel oneAPI / clang - the project's two CI
           compilers).
*/

#ifndef CD_PARTICLELAYOUTMACRO_H
#define CD_PARTICLELAYOUTMACRO_H

// Std includes
#include <tuple>

// Our includes
#include "CD_ParticleSoA.H"

// ---------------------------------------------------------------------------
// Preprocessor plumbing (count args, map a macro over a member list).
// ---------------------------------------------------------------------------
#define CD_PL_EXPAND(x) x
#define CD_PL_CAT(a, b) CD_PL_CAT_(a, b)
#define CD_PL_CAT_(a, b) a##b

// Number of variadic arguments, 1..16.
#define CD_PL_NARG(...) CD_PL_EXPAND(CD_PL_NARG_(__VA_ARGS__, CD_PL_RSEQ()))
#define CD_PL_NARG_(...) CD_PL_EXPAND(CD_PL_ARG_N(__VA_ARGS__))
#define CD_PL_ARG_N(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, N, ...) N
#define CD_PL_RSEQ() 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1

// Turn one member name into a pointer-to-member.
#define CD_PL_ADDR(TYPE, m) &TYPE::m

// Apply f(TYPE, member) across the member list, comma-separated.
#define CD_PL_MAP_1(f, T, a) f(T, a)
#define CD_PL_MAP_2(f, T, a, ...) f(T, a), CD_PL_EXPAND(CD_PL_MAP_1(f, T, __VA_ARGS__))
#define CD_PL_MAP_3(f, T, a, ...) f(T, a), CD_PL_EXPAND(CD_PL_MAP_2(f, T, __VA_ARGS__))
#define CD_PL_MAP_4(f, T, a, ...) f(T, a), CD_PL_EXPAND(CD_PL_MAP_3(f, T, __VA_ARGS__))
#define CD_PL_MAP_5(f, T, a, ...) f(T, a), CD_PL_EXPAND(CD_PL_MAP_4(f, T, __VA_ARGS__))
#define CD_PL_MAP_6(f, T, a, ...) f(T, a), CD_PL_EXPAND(CD_PL_MAP_5(f, T, __VA_ARGS__))
#define CD_PL_MAP_7(f, T, a, ...) f(T, a), CD_PL_EXPAND(CD_PL_MAP_6(f, T, __VA_ARGS__))
#define CD_PL_MAP_8(f, T, a, ...) f(T, a), CD_PL_EXPAND(CD_PL_MAP_7(f, T, __VA_ARGS__))
#define CD_PL_MAP_9(f, T, a, ...) f(T, a), CD_PL_EXPAND(CD_PL_MAP_8(f, T, __VA_ARGS__))
#define CD_PL_MAP_10(f, T, a, ...) f(T, a), CD_PL_EXPAND(CD_PL_MAP_9(f, T, __VA_ARGS__))
#define CD_PL_MAP_11(f, T, a, ...) f(T, a), CD_PL_EXPAND(CD_PL_MAP_10(f, T, __VA_ARGS__))
#define CD_PL_MAP_12(f, T, a, ...) f(T, a), CD_PL_EXPAND(CD_PL_MAP_11(f, T, __VA_ARGS__))
#define CD_PL_MAP_13(f, T, a, ...) f(T, a), CD_PL_EXPAND(CD_PL_MAP_12(f, T, __VA_ARGS__))
#define CD_PL_MAP_14(f, T, a, ...) f(T, a), CD_PL_EXPAND(CD_PL_MAP_13(f, T, __VA_ARGS__))
#define CD_PL_MAP_15(f, T, a, ...) f(T, a), CD_PL_EXPAND(CD_PL_MAP_14(f, T, __VA_ARGS__))
#define CD_PL_MAP_16(f, T, a, ...) f(T, a), CD_PL_EXPAND(CD_PL_MAP_15(f, T, __VA_ARGS__))
#define CD_PL_MAP(f, T, ...) CD_PL_EXPAND(CD_PL_CAT(CD_PL_MAP_, CD_PL_NARG(__VA_ARGS__))(f, T, __VA_ARGS__))

// ---------------------------------------------------------------------------
// The public macro.
// ---------------------------------------------------------------------------

/**
  @brief Generate a ParticleTraits specialization for TYPE.
  @details TYPE is the particle struct, POS the position member name, WEIGHT the
  weight member name; the remaining names are additional columns. Expands to the
  same explicit ParticleTraits a user would otherwise write by hand.
*/
#define CD_PARTICLE_LAYOUT(TYPE, POS, WEIGHT, ...)                                                                \
  template <>                                                                                                     \
  struct ParticleTraits<TYPE>                                                                                     \
  {                                                                                                               \
    static constexpr auto columns     = std::make_tuple(CD_PL_MAP(CD_PL_ADDR, TYPE, POS, WEIGHT, ##__VA_ARGS__)); \
    static constexpr auto positionPtr = &TYPE::POS;                                                               \
    static constexpr auto weightPtr   = &TYPE::WEIGHT;                                                            \
  }

#endif
