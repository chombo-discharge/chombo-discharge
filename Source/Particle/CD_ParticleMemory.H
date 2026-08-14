/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ParticleMemory.H
 * @brief  Process-wide accounting of the memory held by particle data
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
 * @brief Running totals of the memory held by particle data on this rank.
 * @details Chombo's memory tracking counts only allocations that pass through its Arena, i.e.
 * BaseFab-backed mesh data. Particle storage never reaches it -- ParticleSoA allocates with
 * std::aligned_alloc and the exchange buffers are plain std::vector -- so on a particle-heavy run the
 * largest term in the footprint is invisible to overallMemoryUsage(). This namespace closes that gap.
 *
 * The total is split by kind, because the two answer different questions and have different remedies:
 *
 *   - Container: the arenas owned by ParticleSoA. Capacity a rank is *holding*, including the slack
 *     that geometric growth leaves behind. Released only when a container is destroyed or shrunk, so
 *     it persists between operations.
 *   - Buffer: the flat send and receive buffers of a cross-rank exchange. Payload *in flight*. Exists
 *     only for the duration of one exchange, and is sized by how much is being moved rather than by
 *     how much a rank owns.
 *
 * A peak that is mostly Container means a rank is retaining capacity, and shrinking it is the fix. A
 * peak that is mostly Buffer means too much is moved at once, and chunking the exchange is the fix.
 * The summed total cannot tell them apart, which is the reason for the split.
 *
 * Peaks are tracked per kind and for the total. The total peak is not the sum of the per-kind peaks,
 * because the two need not peak at the same moment.
 */
namespace ParticleMemory {

/**
 * @brief The kinds of particle memory tracked here. See the namespace documentation.
 */
enum class Kind
{
  Container = 0, ///< ParticleSoA arenas: capacity a rank is holding.
  Buffer    = 1, ///< Flat exchange buffers: payload in flight.
  NumKinds  = 2  ///< Number of kinds.
};

/**
 * @brief Bytes currently held, per kind. Use the accessors rather than touching this directly.
 * @details Atomic because arenas are grown from inside OpenMP regions -- delinearizeAndAppend()
 * mutates one leaf per thread during a remap. Relaxed ordering suffices: these are diagnostic totals
 * with no happens-before relationship to the data they describe.
 */
inline std::atomic<std::size_t> s_bytes[static_cast<int>(Kind::NumKinds)]{};

/**
 * @brief High-water mark per kind. See s_bytes.
 */
inline std::atomic<std::size_t> s_peak[static_cast<int>(Kind::NumKinds)]{};

/**
 * @brief Bytes currently held across all kinds.
 */
inline std::atomic<std::size_t> s_totalBytes{0};

/**
 * @brief High-water mark of s_totalBytes.
 * @details The figure that decides an out-of-memory failure, and one the step report cannot otherwise
 * see: the largest particle footprint of a step falls midway through an exchange, and by the end of
 * the step the transient is gone again.
 */
inline std::atomic<std::size_t> s_totalPeak{0};

/**
 * @brief Raise a high-water mark to a_value unless it is already at least that.
 * @details The loop reloads on contention, because a competing thread may have raised the mark in the
 * meantime; it exits as soon as the mark is at least a_value, so a thread that loses the race to a
 * larger value simply stops.
 * @param[in,out] a_mark  The mark to raise.
 * @param[in]     a_value Candidate value.
 */
inline void
raisePeak(std::atomic<std::size_t>& a_mark, const std::size_t a_value) noexcept
{
  std::size_t seen = a_mark.load(std::memory_order_relaxed);

  while (a_value > seen && !a_mark.compare_exchange_weak(seen, a_value, std::memory_order_relaxed)) {
  }
}

/**
 * @brief Record an allocation.
 * @param[in] a_kind  Which kind of particle memory this is.
 * @param[in] a_bytes Size of the allocation, in bytes.
 */
inline void
addBytes(const Kind a_kind, const std::size_t a_bytes) noexcept
{
  const int k = static_cast<int>(a_kind);

  raisePeak(s_peak[k], s_bytes[k].fetch_add(a_bytes, std::memory_order_relaxed) + a_bytes);
  raisePeak(s_totalPeak, s_totalBytes.fetch_add(a_bytes, std::memory_order_relaxed) + a_bytes);
}

/**
 * @brief Record a release.
 * @param[in] a_kind  Which kind of particle memory this is.
 * @param[in] a_bytes Size of the release, in bytes. Must be what addBytes() was called with.
 */
inline void
removeBytes(const Kind a_kind, const std::size_t a_bytes) noexcept
{
  s_bytes[static_cast<int>(a_kind)].fetch_sub(a_bytes, std::memory_order_relaxed);
  s_totalBytes.fetch_sub(a_bytes, std::memory_order_relaxed);
}

/**
 * @brief Get the bytes currently held by one kind.
 * @param[in] a_kind Which kind.
 * @return Bytes currently held.
 */
inline std::size_t
getBytes(const Kind a_kind) noexcept
{
  return s_bytes[static_cast<int>(a_kind)].load(std::memory_order_relaxed);
}

/**
 * @brief Get the high-water mark of one kind.
 * @param[in] a_kind Which kind.
 * @return Peak bytes held by that kind.
 */
inline std::size_t
getPeak(const Kind a_kind) noexcept
{
  return s_peak[static_cast<int>(a_kind)].load(std::memory_order_relaxed);
}

/**
 * @brief Get the bytes currently held by particle data of every kind.
 * @return Total bytes currently held.
 */
inline std::size_t
getAllocatedBytes() noexcept
{
  return s_totalBytes.load(std::memory_order_relaxed);
}

/**
 * @brief Get the high-water mark across all kinds.
 * @return Peak total bytes.
 */
inline std::size_t
getPeakBytes() noexcept
{
  return s_totalPeak.load(std::memory_order_relaxed);
}
} // namespace ParticleMemory

#include <CD_NamespaceFooter.H>

#endif
