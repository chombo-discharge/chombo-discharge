/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ParticleSoAArena.H
  @brief  Arena-backed (single-allocation) variant of ParticleSoA.
  @author Robert Marskar
  @details
  Same column-wise (SoA) data model and traits as ParticleSoA, but ALL columns live in
  ONE contiguous, aligned allocation (the "arena"), laid out column-major:

      [ column 0 (capacity slots) | pad | column 1 (capacity slots) | pad | ... ]

  Column base pointers are cached (recomputed only on reserve/grow), so element access
  is as cheap as the N-vector version. The point of the arena is that bulk operations
  over a *compact* container (size == capacity) are a SINGLE contiguous span:
    - build      : one allocation (not N)
    - copy / MPI pack : one memcpy of data() (glibc auto-streams; no intrinsics), or a
                        zero-copy MPI_Send of data() for a whole-container transfer.

  Cost: growth/reserve reallocates the whole arena and moves every column (vs N vectors
  growing independently), so reserve to capacity. Offsets depend on capacity, so a
  one-memcpy span is only fully tight when size == capacity (reserve-exact); with spare
  capacity the per-column used regions have gaps.
*/

#ifndef CD_PARTICLESOAARENA_H
#define CD_PARTICLESOAARENA_H

// Std includes
#include <algorithm>
#include <array>
#include <cstdlib>
#include <cstring>
#include <tuple>
#include <utility>

// Our includes
#include <CD_ParticleSoA.H> // ParticleTraits<P> + detail:: helpers
#include <CD_NamespaceHeader.H>

/**
  @brief Arena-backed SoA particle container. See file header for the layout.
  @tparam P      User particle struct.
  @tparam Traits Column descriptor; defaults to ParticleTraits<P>.
*/
template <typename P, typename Traits = ParticleTraits<P>>
class ParticleSoAArena
{
public:
  using ColumnsTuple = std::remove_cv_t<std::remove_reference_t<decltype(Traits::columns)>>;

  /** @brief Number of SoA columns. */
  static constexpr std::size_t s_numColumns = std::tuple_size<ColumnsTuple>::value;

  /** @brief Value type of column K. */
  template <std::size_t K>
  using ColumnType = typename detail::MemberPointerTraits<std::tuple_element_t<K, ColumnsTuple>>::ValueType;

  /** @brief Cache-line alignment for each column (also enables aligned SIMD). */
  static constexpr std::size_t s_align = 64;

  /** @brief Derive a column index from a pointer-to-member (reorder-safe selector). */
  template <typename M>
  static constexpr std::size_t
  columnIndex(const M a_member) noexcept
  {
    return detail::indexOf(Traits::columns, a_member);
  }

  static constexpr std::size_t s_positionIndex = detail::indexOf(Traits::columns, Traits::positionPtr);
  static constexpr std::size_t s_weightIndex   = detail::indexOf(Traits::columns, Traits::weightPtr);

  static_assert(std::is_same<ColumnType<s_positionIndex>, RealVect>::value, "positionPtr must be a RealVect member");
  static_assert(std::is_same<ColumnType<s_weightIndex>, Real>::value, "weightPtr must be a Real member");

  ParticleSoAArena() = default;

  ~ParticleSoAArena()
  {
    std::free(m_buffer);
  }

  ParticleSoAArena(const ParticleSoAArena&) = delete;
  ParticleSoAArena&
  operator=(const ParticleSoAArena&) = delete;

  /** @brief Move constructor (steals the arena). */
  ParticleSoAArena(ParticleSoAArena&& a_other) noexcept
  {
    this->steal(a_other);
  }

  /** @brief Move assignment. @return Reference to this. */
  ParticleSoAArena&
  operator=(ParticleSoAArena&& a_other) noexcept
  {
    if (this != &a_other) {
      std::free(m_buffer);
      this->steal(a_other);
    }
    return *this;
  }

  /** @brief Number of particles. */
  std::size_t
  size() const noexcept
  {
    return m_size;
  }

  /** @brief Allocated capacity in particles. */
  std::size_t
  capacity() const noexcept
  {
    return m_capacity;
  }

  /** @brief Whether the container is empty. */
  bool
  empty() const noexcept
  {
    return m_size == 0;
  }

  /** @brief Drop all particles (keeps the arena). */
  void
  clear() noexcept
  {
    m_size = 0;
  }

  /**
    @brief Ensure capacity for at least a_capacity particles (reallocates + moves columns if grown).
    @param[in] a_capacity Requested capacity.
  */
  void
  reserve(const std::size_t a_capacity)
  {
    if (a_capacity <= m_capacity) {
      return;
    }
    this->reallocate(a_capacity);
  }

  /**
    @brief Set the particle count to a_count (grows capacity if needed; new slots are uninitialized).
    @param[in] a_count New particle count.
  */
  void
  resize(const std::size_t a_count)
  {
    this->reserve(a_count);
    m_size = a_count;
  }

  /** @brief Typed base pointer of column K. */
  template <std::size_t K>
  ColumnType<K>*
  columnByIndex() noexcept
  {
    return reinterpret_cast<ColumnType<K>*>(m_colPtrs[K]);
  }

  /** @brief Const typed base pointer of column K. */
  template <std::size_t K>
  const ColumnType<K>*
  columnByIndex() const noexcept
  {
    return reinterpret_cast<const ColumnType<K>*>(m_colPtrs[K]);
  }

  /** @brief Column base pointer selected by member pointer. */
  template <auto Member>
  auto*
  column() noexcept
  {
    return this->columnByIndex<columnIndex(Member)>();
  }

  /** @brief Const column base pointer selected by member pointer. */
  template <auto Member>
  const auto*
  column() const noexcept
  {
    return this->columnByIndex<columnIndex(Member)>();
  }

  /** @brief Position of particle i. */
  RealVect&
  position(const std::size_t a_index) noexcept
  {
    return this->columnByIndex<s_positionIndex>()[a_index];
  }

  /** @brief Const position of particle i. */
  const RealVect&
  position(const std::size_t a_index) const noexcept
  {
    return this->columnByIndex<s_positionIndex>()[a_index];
  }

  /** @brief Weight of particle i. */
  Real&
  weight(const std::size_t a_index) noexcept
  {
    return this->columnByIndex<s_weightIndex>()[a_index];
  }

  /** @brief Append one particle (grows if at capacity). @param[in] a_p Particle. */
  void
  append(const P& a_p)
  {
    if (m_size == m_capacity) {
      this->reserve(m_capacity == 0 ? 16 : 2 * m_capacity);
    }
    this->appendImpl(a_p, std::make_index_sequence<s_numColumns>{});
    m_size++;
  }

  /** @brief Gather particle i into the AoS view. @param[in] a_index Index. @return Particle. */
  P
  gather(const std::size_t a_index) const
  {
    return this->gatherImpl(a_index, std::make_index_sequence<s_numColumns>{});
  }

  /**
    @brief Raw arena pointer (column-major). For a compact container this spans all data.
    @return Pointer to the arena.
  */
  const void*
  data() const noexcept
  {
    return m_buffer;
  }

  /**
    @brief Number of bytes spanned by the live data (offset of the last column + its used region).
    @details For a compact container (size == capacity) this is one contiguous span and a
    single memcpy/MPI_Send of [data(), data()+byteSpan()) transfers the whole container.
    @return Byte span.
  */
  std::size_t
  byteSpan() const noexcept
  {
    if (m_size == 0) {
      return 0;
    }
    return this->byteSpanImpl(std::make_index_sequence<s_numColumns>{});
  }

protected:
  void*                           m_buffer   = nullptr; ///< The single arena allocation.
  std::size_t                     m_capacity = 0;       ///< Capacity in particles.
  std::size_t                     m_size     = 0;       ///< Live particle count.
  std::array<void*, s_numColumns> m_colPtrs{};          ///< Cached column base pointers into m_buffer.

private:
  static std::size_t
  alignUp(const std::size_t a_x) noexcept
  {
    return (a_x + (s_align - 1)) / s_align * s_align;
  }

  void
  steal(ParticleSoAArena& a_o) noexcept
  {
    m_buffer       = a_o.m_buffer;
    m_capacity     = a_o.m_capacity;
    m_size         = a_o.m_size;
    m_colPtrs      = a_o.m_colPtrs;
    a_o.m_buffer   = nullptr;
    a_o.m_capacity = 0;
    a_o.m_size     = 0;
    a_o.m_colPtrs  = {};
  }

  template <std::size_t... I>
  std::size_t
  computeOffsets(const std::size_t                      a_capacity,
                 std::array<std::size_t, s_numColumns>& a_off,
                 std::index_sequence<I...>) const noexcept
  {
    std::size_t running = 0;
    using expander      = int[];
    (void)expander{0, (a_off[I] = running, running += alignUp(a_capacity * sizeof(ColumnType<I>)), 0)...};
    return running;
  }

  template <std::size_t... I>
  void
  moveColumns(void* a_newBuf, const std::array<std::size_t, s_numColumns>& a_off, std::index_sequence<I...>) noexcept
  {
    using expander = int[];
    (void)expander{
      0,
      (std::memcpy(static_cast<char*>(a_newBuf) + a_off[I], m_colPtrs[I], m_size * sizeof(ColumnType<I>)), 0)...};
  }

  void
  reallocate(const std::size_t a_capacity)
  {
    std::array<std::size_t, s_numColumns> off{};
    const std::size_t total  = computeOffsets(a_capacity, off, std::make_index_sequence<s_numColumns>{});
    void* const       newBuf = std::aligned_alloc(s_align, alignUp(total));
    if (m_size > 0) {
      this->moveColumns(newBuf, off, std::make_index_sequence<s_numColumns>{});
    }
    std::free(m_buffer);
    m_buffer   = newBuf;
    m_capacity = a_capacity;
    this->setColPtrs(off, std::make_index_sequence<s_numColumns>{});
  }

  template <std::size_t... I>
  void
  setColPtrs(const std::array<std::size_t, s_numColumns>& a_off, std::index_sequence<I...>) noexcept
  {
    using expander = int[];
    (void)expander{0, (m_colPtrs[I] = static_cast<char*>(m_buffer) + a_off[I], 0)...};
  }

  template <std::size_t... I>
  void
  appendImpl(const P& a_p, std::index_sequence<I...>)
  {
    constexpr auto cols = Traits::columns;
    using expander      = int[];
    (void)expander{0, ((void)(this->columnByIndex<I>()[m_size] = a_p.*std::get<I>(cols)), 0)...};
  }

  template <std::size_t... I>
  P
  gatherImpl(const std::size_t a_index, std::index_sequence<I...>) const
  {
    constexpr auto cols = Traits::columns;
    P              p;
    using expander = int[];
    (void)expander{0, ((void)(p.*std::get<I>(cols) = this->columnByIndex<I>()[a_index]), 0)...};
    return p;
  }

  template <std::size_t... I>
  std::size_t
  byteSpanImpl(std::index_sequence<I...>) const noexcept
  {
    std::size_t span = 0;
    using expander   = int[];
    // end of each column's live region; the max is the live span.
    (void)expander{0,
                   (span = std::max(span,
                                    static_cast<std::size_t>(static_cast<const char*>(m_colPtrs[I]) -
                                                             static_cast<const char*>(m_buffer)) +
                                      m_size * sizeof(ColumnType<I>)),
                    0)...};
    return span;
  }
};

#include <CD_NamespaceFooter.H>

#endif
