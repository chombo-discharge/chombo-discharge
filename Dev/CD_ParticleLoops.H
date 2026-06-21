/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ParticleLoops.H
  @brief  Declaration of a namespace for SIMD-decorated loops over SoA particles.
  @author Robert Marskar
  @details
  ParticleLoops is the particle analogue of BoxLoops: it provides @c loop overloads
  that iterate over the particles of a ParticleSoA and decorate the inner loop with
  @c CD_PRAGMA_SIMD, so elementwise per-particle kernels auto-vectorize. The pattern
  mirrors BoxLoops exactly, with the cell index replaced by the particle index:

  @code{.cpp}
    // column pointers captured once, OUTSIDE the loop
    Real*               sigma = soa.column<&P::conductivity>();
    const ParticleReal* w     = soa.weightColumn();
    const Real*         mu    = soa.column<&P::mobility>();

    ParticleLoops::loop(soa, [&](std::size_t i) { sigma[i] = w[i] * mu[i]; });
  @endcode

  <h3>Why this exists</h3>
  CD_PRAGMA_SIMD (GCC @c ivdep / clang @c vectorize(enable)) tells the compiler to
  ignore assumed loop-carried/aliasing dependencies, so the kernel vectorizes
  without the caller having to launder column pointers through @c __restrict, and
  without the compiler emitting runtime aliasing checks ("loop versioning"). ALWAYS_INLINE
  guarantees the kernel is inlined into the loop body.

  <h3>Contract (read carefully)</h3>
  - The kernel takes the particle index: @c (std::size_t i). Capture the column
    references/pointers ONCE before the loop (as above) -- do not call
    @c soa.column<>() inside the kernel.
  - Elementwise read/transform ONLY. Do NOT append(), remove(), or resize() the
    container inside the loop -- the trip count is fixed and such mutation is a data
    hazard. Structural changes belong in ordinary index loops.
  - Not for reductions. @c ivdep does not guarantee a correct/vectorized reduction;
    write sums/min/max in a plain loop (a dedicated reduce() overload can be added
    later if needed).
  - This is the INNER (per-patch, per-particle) loop. Outer OpenMP-over-patches
    parallelism is separate, exactly as BoxLoops sits inside the omp-over-boxes region.

  @note As in CD_Decorations.H, CD_PRAGMA_SIMD is disabled in debug builds (when
        NDEBUG is not defined); the loop is then an ordinary correct loop.
*/

#ifndef CD_PARTICLELOOPS_H
#define CD_PARTICLELOOPS_H

// Std includes
#include <cstddef>
#include <utility>

// Our includes
#include <CD_Decorations.H>
#include "CD_ParticleSoA.H"
#include <CD_NamespaceHeader.H>

/**
  @brief Namespace for SIMD-decorated loops over SoA particles (particle analogue of BoxLoops).
*/
namespace ParticleLoops {

  /**
    @brief Launch a kernel over [0, a_numParticles), decorating the loop with CD_PRAGMA_SIMD.
    @details The kernel is invoked as @c a_kernel(i) for each particle index i. Force-inlined
    so the kernel body fuses into the vectorizable loop.
    @param[in]     a_numParticles Number of particles to iterate over.
    @param[in,out] a_kernel       Callable taking a std::size_t particle index.
  */
  template <typename Functor>
  ALWAYS_INLINE void
  loop(const std::size_t a_numParticles, Functor&& a_kernel)
  {
    CD_PRAGMA_SIMD
    for (std::size_t i = 0; i < a_numParticles; i++) {
      a_kernel(i);
    }
  }

  /**
    @brief Convenience overload iterating over every particle in a ParticleSoA.
    @details Equivalent to loop(a_soa.size(), kernel). The container is only queried for its
    size; the kernel should capture the columns it needs from the caller's scope.
    @param[in]     a_soa    Particle container to iterate over.
    @param[in,out] a_kernel Callable taking a std::size_t particle index.
  */
  template <typename P, typename Traits, typename Functor>
  ALWAYS_INLINE void
  loop(const ParticleSoA<P, Traits>& a_soa, Functor&& a_kernel)
  {
    ParticleLoops::loop(a_soa.size(), std::forward<Functor>(a_kernel));
  }

} // namespace ParticleLoops

#include <CD_NamespaceFooter.H>

#endif
