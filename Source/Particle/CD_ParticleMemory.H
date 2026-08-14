/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ParticleMemory.H
 * @brief  Process-wide accounting of the memory held by particle arenas
 * @author Robert Marskar
 */

#ifndef CD_PARTICLEMEMORY_H
#define CD_PARTICLEMEMORY_H

// Std includes
#include <atomic>
#include <cstddef>

// Our includes
#include <CD_NamespaceHeader.H>

/**
 * @brief Running total of the bytes held by every ParticleSoA arena on this rank.
 * @details Chombo's memory tracking counts only allocations that pass through its Arena, i.e.
 * BaseFab-backed mesh data. ParticleSoA allocates with std::aligned_alloc, so particle storage is
 * invisible to overallMemoryUsage() -- on a particle-heavy run it is the dominant term and the
 * tracked figure can be a small fraction of the resident set. This namespace closes that gap by
 * having ParticleSoA report every arena allocation and release, so the total can be printed
 * alongside the tracked and resident figures.
 *
 * The total covers *every* live ParticleSoA, including the transient ones built during a remap (the
 * mover pool), which is the point: those transients are exactly where the peak lives, and no
 * per-solver accounting would see them.
 *
 * This is allocated bytes, not live-particle bytes. The difference is the arena slack that
 * ParticleSoA's geometric growth leaves behind, and it is the allocated figure that occupies the
 * resident set.
 */
namespace ParticleMemory {

/**
 * @brief The running total, in bytes. Use the accessors rather than touching this directly.
 * @details Atomic because arenas are grown from inside OpenMP regions -- delinearizeAndAppend()
 * mutates one leaf per thread during a remap. Relaxed ordering is sufficient: the counter is a
 * diagnostic total with no happens-before relationship to the data it describes.
 */
inline std::atomic<std::size_t> s_allocatedBytes{0};

/**
 * @brief Record an arena allocation.
 * @param[in] a_bytes Size of the newly allocated arena, in bytes.
 */
inline void
addBytes(const std::size_t a_bytes) noexcept
{
  s_allocatedBytes.fetch_add(a_bytes, std::memory_order_relaxed);
}

/**
 * @brief Record an arena release.
 * @param[in] a_bytes Size of the released arena, in bytes. Must be what addBytes() was called with.
 */
inline void
removeBytes(const std::size_t a_bytes) noexcept
{
  s_allocatedBytes.fetch_sub(a_bytes, std::memory_order_relaxed);
}

/**
 * @brief Get the bytes currently held by particle arenas on this rank.
 * @return Allocated arena bytes.
 */
inline std::size_t
getAllocatedBytes() noexcept
{
  return s_allocatedBytes.load(std::memory_order_relaxed);
}
} // namespace ParticleMemory

#include <CD_NamespaceFooter.H>

#endif
