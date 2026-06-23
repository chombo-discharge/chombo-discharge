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
  @c CD_PRAGMA_SIMD, so elementwise per-particle kernels (hopefully) auto-vectorize.
  This mirrors BoxLoops exactly, with the cell index replaced by the particle index:

  @code{.cpp}
    // column pointers captured once, OUTSIDE the loop
    ParticleReal*       sigma = soa.column<&P::conductivity>(); // payload column -> ParticleReal*
    const double*       w     = soa.weightColumn();             // weight is always double
    const ParticleReal* mu    = soa.column<&P::mobility>();     // payload column

    ParticleLoops::loop(soa, [&](std::size_t i) { sigma[i] = w[i] * mu[i]; });
  @endcode

  <h3>Why this exists</h3>
  CD_PRAGMA_SIMD (GCC @c ivdep / clang @c vectorize(enable)) tells the compiler to
  ignore assumed loop-carried/aliasing dependencies, so the kernel vectorizes
  without the caller having to launder column pointers through @c __restrict, and
  without the compiler emitting runtime aliasing checks ("loop versioning"). ALWAYS_INLINE
  guarantees the kernel is inlined into the loop body.

  <h3>Contract</h3>
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

// Our includes
#include <CD_Decorations.H>
#include <CD_ParticleSoA.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Namespace for SIMD-decorated loops over SoA particles (particle analogue of BoxLoops).
*/
namespace ParticleLoops {

  /**
    @brief Launch a kernel over every particle in a ParticleSoA, decorating the loop with CD_PRAGMA_SIMD.
    @details The kernel is invoked as @c a_kernel(i) for each particle index i in [0, a_soa.size()).
    The trip count is taken from the container, so it always matches the columns the kernel captured
    from the same container -- there is intentionally no bare-count overload to pass the wrong size to.
    Force-inlined so the kernel body fuses into the vectorizable loop.

    @note a_soa is taken by const reference because loop() itself only reads a_soa.size() -- it never
          touches particle data. Whether the kernel may WRITE columns is governed entirely by how the
          caller captured the column pointers: column()/positionColumn()/weightColumn() return mutable
          pointers from a non-const container and read-only pointers from a const one. This mirrors
          BoxLoops::loop(const Box&, ...), where the iteration domain is const and the kernel mutates
          separately-captured FArrayBox data.
    @tparam P       Payload type of the ParticleSoA.
    @tparam Traits  Column descriptor for the payload.
    @tparam Functor Callable type taking a std::size_t particle index.
    @param[in]     a_soa    Particle container to iterate over (read for its size only).
    @param[in,out] a_kernel Callable taking a std::size_t particle index.
  */
  template <typename P, typename Traits, typename Functor>
  ALWAYS_INLINE void
  loop(const ParticleSoA<P, Traits>& a_soa, Functor&& a_kernel)
  {
    const std::size_t numParticles = a_soa.size();

    CD_PRAGMA_SIMD
    for (std::size_t i = 0; i < numParticles; i++) {
      a_kernel(i);
    }
  }

} // namespace ParticleLoops

#include <CD_NamespaceFooter.H>

#endif
