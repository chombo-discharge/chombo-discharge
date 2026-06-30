/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ParticleSoA.H
  @brief  Declaration of ParticleSoA, an arena-backed Struct-of-Arrays particle container.
  @author Robert Marskar
  @details ParticleSoA is the per-patch particle storage that replaces Chombo's @c List<P>. It is
  the single source of truth for one grid patch worth of particles; the AMR machinery
  (the ParticleContainer equivalent) owns a @c LayoutData<ParticleSoA<P>> per level.

  <h2>The data model </h2>

  The container stores particles column-major (Struct-of-Arrays): one contiguous array per
  field. ALL columns live in ONE aligned allocation (the "arena"), laid out as offset
  slices. Position and weight are stored in @c double precision -- they index the grid and are
  summed/conserved across the whole population, so they must not lose precision to float
  roundoff (and must stay double even in a CH_USE_FLOAT build, hence @c double, not @c Real).
  Only the user payload columns may use the narrower @c ParticleReal.

  @verbatim
      [ x | pad | y | pad | z | pad | weight | pad | id | pad | rank | pad | ghost | pad | payload... ]
  @endverbatim

  Column base pointers are cached, so element access is as cheap as a plain array. Bulk
  operations over a compact (size == capacity) container are then a SINGLE contiguous span:
  one allocation to build, one @c memcpy (or a zero-copy @c MPI_Send of data()) to transfer.

  <h2>Mandatory columns are container-owned; the user supplies only a payload</h2>

  Every particle in chombo-discharge has a position, a weight, and metadata (a global id +
  the owning rank). ParticleSoA OWNS those columns and the user does NOT declare them:

    - position : @c SpaceDim raw scalar @c double columns (NOT a RealVect column).
                 Per-component storage is what lets the field-update kernels vectorize and is
                 GPU-optimal; @c double precision is mandatory because a float global position
                 loses sub-cell resolution on fine grids (at ~10^6 cells/dim the float ULP is
                 ~0.1 cell, so small advective steps round to zero -- particles freeze).
    - weight   : one @c double column, because per-cell Sum(weight) and merge/split
                 center-of-mass conservation accumulate float roundoff over the whole
                 population, and superparticle weights span a huge dynamic range.
    - particleID / rankID : container-owned bookkeeping columns of fixed-width integer type
                 (std::int64_t / std::int32_t -- portable across the byte-blit I/O paths). They
                 travel in MPI linearization but are NOT written to HDF5 (regenerated on
                 restart).

  The user defines ONLY the extra payload fields (velocity, mobility, energy, ...) as a
  plain struct, plus a ParticleTraits specialization listing the payload columns. The
  payload may be empty, so @c ParticleSoA<> is a ready-made point/tracer particle
  (position + weight + metadata, nothing else).

  <h2>Defining a particle type (the complete recipe)</h2>

  Step 1 -- write a plain struct of the EXTRA (payload) members only. Each member becomes one
  SoA column and must be trivially copyable. Do NOT add position/weight/id/rank; the
  container already owns those. Vector fields are stored PER-COMPONENT
  as scalar columns -- there is no RealVect column type; declare x/y/z members and promote to
  RealVect in downstream code if needed. Use ParticleReal for local per-particle physics
  (velocity, mobility, energy) where float is safe; use @c double for anything position-like or
  globally reduced (e.g. an oldPosition, which subtracts from the double position).

    @code
      struct MyPayload
      {
        ParticleReal D_DECL(vx = 0.0, vy = 0.0, vz = 0.0); // velocity, per-component
        ParticleReal mobility = 0.0;
        ParticleReal energy   = 0.0;
      };
    @endcode

  Step 2 -- specialize ParticleTraits<MyPayload> with the @c columns tuple (one
  pointer-to-member per payload column; the single source of truth for the layout):

    @code
      template <>
      struct ParticleTraits<MyPayload>
      {
        static constexpr auto columns = std::make_tuple(D_DECL(&MyPayload::vx,
                                                               &MyPayload::vy,
                                                               &MyPayload::vz),
                                                        &MyPayload::mobility,
                                                        &MyPayload::energy);
      };
    @endcode

  Step 3 (optional) -- restrict which PAYLOAD columns go to HDF5 by declaring an
  @c h5Columns tuple of member pointers. Position and weight are ALWAYS checkpointed;
  id/rank NEVER are. If @c h5Columns is omitted, every payload column is checkpointed.

    @code
      static constexpr auto h5Columns = std::make_tuple(&MyPayload::energy);
    @endcode

  That is all -- @c ParticleSoA<MyPayload> now works. Payload fields are selected by member
  pointer, e.g. @c soa.column<&MyPayload::energy>(); the column index is DERIVED from the
  pointer at compile time, so reordering @c columns keeps name-based access correct.

  <h2>Precision</h2>
  Position and weight are always @c double; see the data-model note above for why. PAYLOAD
  columns use whatever type the user declares -- @c ParticleReal is a compile-time alias that
  defaults to @c double but may be set to @c float (define @c CD_PARTICLE_REAL) for
  half-memory, 8-wide-SIMD payload against a double-precision grid. position(i) returns a
  @c RealVect for cell lookup / Chombo interop (the one geometry type at the accessor boundary).

  @note ParticleSoA is move-only: copy construction and copy assignment are deleted so an O(N)
        per-column copy can never happen implicitly (matching ParticleContainer). When an independent
        copy is genuinely wanted, use the explicit deepCopy() (returns a new container) or
        deepCopyTo() (reuses an existing destination's arena).
  @tparam P      User payload struct (extra columns only); defaults to NoPayload.
  @tparam Traits Column descriptor for the payload; defaults to ParticleTraits<P>.
*/

#ifndef CD_PARTICLESOA_H
#define CD_PARTICLESOA_H

// Std includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

// Chombo includes (kept minimal -- only the grid-interop geometry types; storage uses
// fundamental C++ types so the container is not tied to Chombo's Real precision).
#include <SPACE.H>
#include <RealVect.H>
#include <IntVect.H>
#include <Box.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Compile-time precision of particle PAYLOAD data. Defaults to double; define
         CD_PARTICLE_REAL (e.g. to float) for mixed-precision payload columns.
  @note Position and weight are NOT governed by this -- they are always double; see
        ParticleSoA's data-model documentation.
*/
#ifndef CD_PARTICLE_REAL
#define CD_PARTICLE_REAL double
#endif

/**
   @brief Floating-point type a user may use for payload columns.
*/
using ParticleReal = CD_PARTICLE_REAL;

/**
   @brief Global particle identifier type (container-owned metadata column; fixed-width for I/O).
*/
using ParticleID = std::int64_t;

/**
   @brief Owning-rank identifier type (container-owned metadata column; fixed-width for I/O).
*/
using RankID = std::int32_t;

/**
   @brief Ghost designation of a particle (container-owned metadata column).
   @details A valid particle (the patch's own, owned particle) has GhostType::Valid. A ghost
   particle is a transient copy that overlaps the patch's ghosted/grown box and is owned
   elsewhere; the non-Valid enumerators record where it came from so the remap routines (and any
   per-patch consumer) know its provenance. Stored as a single byte to keep the per-particle
   footprint minimal.
*/
enum class GhostType : std::uint8_t
{
  Valid     = 0, ///< Owned, valid particle of this patch (not a ghost).
  SameLevel = 1, ///< Ghost from an adjacent patch on the same level.
  Coarse    = 2, ///< Ghost from the next-coarser level (level l-1).
  Fine      = 3  ///< Ghost from the next-finer level (level l+1).
};

/**
  @brief Traits class that a user specializes to describe the PAYLOAD columns of a
         particle type. See ParticleSoA for the full recipe.
  @note The primary template is intentionally undefined, so using ParticleSoA with a
        payload type that has no specialization is a clear hard error.
  @tparam P The user payload struct being described.
*/
template <typename P>
struct ParticleTraits;

/**
   @brief Empty payload: ParticleSoA<> is a position+weight+metadata point particle.
*/
struct NoPayload
{
};

/**
   @brief Traits for the empty payload (no payload columns).
*/
template <>
struct ParticleTraits<NoPayload>
{
  /**
     @brief No payload columns.
  */
  static constexpr auto columns = std::make_tuple();
};

namespace detail {

  /**
    @brief Primary template that extracts the class and value type from a pointer-to-member type.
    @tparam M Pointer-to-member type to decompose.
  */
  template <typename M>
  struct MemberPointerTraits;

  /**
    @brief Specialization that extracts the class and value type from a pointer-to-member type.
    @tparam C Class that owns the member.
    @tparam T Value type of the member.
  */
  template <typename C, typename T>
  struct MemberPointerTraits<T C::*>
  {
    /**
      @brief Class that owns the member.
    */
    using ClassType = C;
    /**
      @brief Value type of the member.
    */
    using ValueType = T;
  };

  /**
    @brief Type-mismatched overload of the pointer-to-member equality (always false).
    @tparam A Type of the left-hand member pointer.
    @tparam B Type of the right-hand member pointer.
    @return Always false (the two member pointers have different types).
  */
  template <typename A, typename B>
  constexpr bool
  eqMember(A, B) noexcept
  {
    return false;
  }

  /**
    @brief Same-type overload of the pointer-to-member equality.
    @tparam A Common type of both member pointers.
    @param[in] a_lhs Left-hand member pointer.
    @param[in] a_rhs Right-hand member pointer.
    @return True iff the two member pointers are equal.
  */
  template <typename A>
  constexpr bool
  eqMember(A a_lhs, A a_rhs) noexcept
  {
    return a_lhs == a_rhs;
  }

  /**
    @brief Base case of the compile-time member-pointer index search (empty sequence).
    @tparam Tup Tuple-of-member-pointer type being searched.
    @tparam M   Member-pointer type being searched for.
    @return The tuple size, signalling "not found".
  */
  template <typename Tup, typename M>
  constexpr std::size_t
  indexOfImpl(const Tup&, M, std::index_sequence<>) noexcept
  {
    return std::tuple_size<Tup>::value;
  }

  /**
    @brief Recursive step of the compile-time member-pointer index search.
    @tparam Tup  Tuple-of-member-pointer type being searched.
    @tparam M    Member-pointer type being searched for.
    @tparam I    Index currently inspected.
    @tparam Rest Remaining indices to inspect.
    @param[in] a_t      Tuple of member pointers being searched.
    @param[in] a_target Member pointer being searched for.
    @return Index of a_target in a_t, or the tuple size if absent.
  */
  template <typename Tup, typename M, std::size_t I, std::size_t... Rest>
  constexpr std::size_t
  indexOfImpl(const Tup& a_t, M a_target, std::index_sequence<I, Rest...>) noexcept
  {
    return eqMember(std::get<I>(a_t), a_target) ? I : indexOfImpl(a_t, a_target, std::index_sequence<Rest...>{});
  }

  /**
    @brief Compile-time index of a member pointer within a tuple (tuple size if absent).
    @tparam Tup Tuple-of-member-pointer type being searched.
    @tparam M   Member-pointer type being searched for.
    @param[in] a_t      Tuple of member pointers being searched.
    @param[in] a_target Member pointer being searched for.
    @return Index of a_target in a_t, or the tuple size if absent.
  */
  template <typename Tup, typename M>
  constexpr std::size_t
  indexOf(const Tup& a_t, M a_target) noexcept
  {
    return indexOfImpl(a_t, a_target, std::make_index_sequence<std::tuple_size<Tup>::value>{});
  }

  /**
    @brief std::tuple of N copies of T (used to build the SpaceDim position columns).
    @tparam T Element type to repeat.
    @tparam N Number of repetitions.
  */
  template <typename T, std::size_t N>
  struct RepeatTuple
  {
    /**
       @brief Helper alias that maps any index to T.
    */
    template <std::size_t>
    using Always = T;

    /**
      @brief Deduce the repeated-tuple type from an index sequence (declaration only).
      @tparam I Indices [0, N).
      @return A std::tuple of N copies of T.
    */
    template <std::size_t... I>
    static std::tuple<Always<I>...> make(std::index_sequence<I...>);

    /**
      @brief The resulting std::tuple type.
    */
    using Type = decltype(make(std::make_index_sequence<N>{}));
  };

  /**
    @brief Primary template mapping a tuple of pointer-to-member to the tuple of value types.
    @tparam Tuple Tuple-of-member-pointer type.
  */
  template <typename Tuple>
  struct ValueTuple;

  /**
    @brief Specialization mapping a tuple of pointer-to-member to the tuple of value types.
    @tparam Ms Member-pointer types.
  */
  template <typename... Ms>
  struct ValueTuple<std::tuple<Ms...>>
  {
    /**
      @brief Tuple of member value types.
    */
    using Type = std::tuple<typename MemberPointerTraits<Ms>::ValueType...>;
  };

  /**
    @brief Base case of the variadic constexpr size sum (no arguments).
    @return Zero.
  */
  constexpr std::size_t
  sumSizes() noexcept
  {
    return 0;
  }

  /**
    @brief Variadic constexpr sum of std::size_t values.
    @tparam Rest Types of the trailing arguments (all std::size_t).
    @param[in] a_first First value to add.
    @param[in] a_rest  Remaining values to add.
    @return The sum of all arguments.
  */
  template <typename... Rest>
  constexpr std::size_t
  sumSizes(std::size_t a_first, Rest... a_rest) noexcept
  {
    return a_first + sumSizes(a_rest...);
  }

  /**
    @brief Primary template for the compile-time sum of sizeof() over a value tuple.
    @tparam Tuple Tuple of value types.
  */
  template <typename Tuple>
  struct SumSizeof;

  /**
    @brief Specialization computing the compile-time sum of sizeof() over a value tuple.
    @tparam Ts Value types in the tuple.
  */
  template <typename... Ts>
  struct SumSizeof<std::tuple<Ts...>>
  {
    /**
      @brief Sum of sizeof(Ts).
    */
    static constexpr std::size_t s_value = sumSizes(sizeof(Ts)...);
  };

  /**
    @brief Primary template for "every type in a value tuple is trivially copyable".
    @tparam Tuple Tuple of value types.
  */
  template <typename Tuple>
  struct AllTriviallyCopyable;

  /**
    @brief Specialization testing that every type in a value tuple is trivially copyable.
    @tparam Ts Value types in the tuple.
  */
  template <typename... Ts>
  struct AllTriviallyCopyable<std::tuple<Ts...>>
  {
    /**
      @brief True iff all Ts are trivially copyable.
    */
    static constexpr bool s_value = (std::is_trivially_copyable<Ts>::value && ...);
  };

  /**
    @brief Stand-in for std::void_t (maps any type pack to void).
    @tparam ... Ignored type pack.
  */
  template <typename...>
  struct MakeVoid
  {
    /**
      @brief Always void.
    */
    using Type = void;
  };

  /**
     @brief Alias yielding void for any type pack (SFINAE helper, like std::void_t).
  */
  template <typename... Ts>
  using Void = typename MakeVoid<Ts...>::Type;

  /**
    @brief Turn a tuple of HDF5 payload member pointers into a std::index_sequence of
           DERIVED payload-column indices (each via indexOf against the payload columns).
    @tparam Traits Payload traits type providing @c columns and @c h5Columns.
    @tparam I      Indices into the @c h5Columns tuple.
    @param[in] a_seq Index sequence [0, size(h5Columns)) selecting the h5 columns.
    @return An index_sequence of the derived payload-column indices.
  */
  template <typename Traits, std::size_t... I>
  constexpr std::index_sequence<indexOf(Traits::columns, std::get<I>(Traits::h5Columns))...>
  deriveH5Seq([[maybe_unused]] std::index_sequence<I...> a_seq) noexcept
  {
    return {};
  }

  /**
    @brief Primary template resolving the HDF5 payload column subset; defaults to all columns.
    @tparam Traits   Payload traits type.
    @tparam Fallback Index sequence used when no @c h5Columns is declared.
    @tparam (anon)   SFINAE slot detecting @c Traits::h5Columns.
  */
  template <typename Traits, typename Fallback, typename = void>
  struct H5ColumnSubset
  {
    /**
      @brief The resolved subset (all payload columns).
    */
    using Type = Fallback;
  };

  /**
    @brief Specialization selecting the declared HDF5 payload subset when @c Traits::h5Columns exists.
    @tparam Traits   Payload traits type declaring @c h5Columns.
    @tparam Fallback Index sequence used when no @c h5Columns is declared (unused here).
  */
  template <typename Traits, typename Fallback>
  struct H5ColumnSubset<Traits, Fallback, Void<decltype(Traits::h5Columns)>>
  {
    /**
      @brief The declared h5Columns tuple type.
    */
    using H5Tuple = std::remove_cv_t<std::remove_reference_t<decltype(Traits::h5Columns)>>;
    /**
      @brief Derived payload-column index subset.
    */
    using Type = decltype(deriveH5Seq<Traits>(std::make_index_sequence<std::tuple_size<H5Tuple>::value>{}));
  };

  /**
    @brief memcpy one trivially-copyable value onto a byte buffer and advance the buffer pointer.
    @tparam T Trivially-copyable value type.
    @param[in,out] a_buf   Buffer pointer, advanced by sizeof(T) on return.
    @param[in]     a_value Value to copy into the buffer.
  */
  template <typename T>
  inline void
  pushBytes(unsigned char*& a_buf, const T& a_value) noexcept
  {
    static_assert(std::is_trivially_copyable<T>::value, "SoA columns must be trivially copyable");
    std::memcpy(a_buf, &a_value, sizeof(T));
    a_buf += sizeof(T);
  }

  /**
    @brief memcpy one trivially-copyable value off a byte buffer and advance the buffer pointer.
    @tparam T Trivially-copyable value type.
    @param[in,out] a_buf   Buffer pointer, advanced by sizeof(T) on return.
    @param[out]    a_value Value filled from the buffer.
  */
  template <typename T>
  inline void
  pullBytes(const unsigned char*& a_buf, T& a_value) noexcept
  {
    static_assert(std::is_trivially_copyable<T>::value, "SoA columns must be trivially copyable");
    std::memcpy(&a_value, a_buf, sizeof(T));
    a_buf += sizeof(T);
  }

} // namespace detail

/**
  @brief Arena-backed Struct-of-Arrays particle container for a single grid patch.
  @details See the file-level documentation for the data model, the container-owned
  mandatory columns, and the recipe for defining a payload type.

  <h3>Accessing data</h3>
  - Mandatory fields: position(i) (RealVect, by value), setPosition(i,x),
    positionColumn(dir) (raw double*), weight(i), weightColumn() (raw double*), particleID(i),
    rankID(i).
  - Payload fields: column<&P::field>() (raw pointer to the whole column) and
    get<&P::field>(i) (one particle's field), both reorder-safe and zero-overhead.
  - gather(i) copies particle i's payload out into a P (COPY -- keep out of hot loops).

  <h3>Mutating the population</h3>
  append(pos, weight, payload) adds a particle (grows geometrically); remove(i) deletes
  particle i by swap-and-pop (O(1), does NOT preserve order). The canonical remove loop
  does not advance i on a removal:
    @code
      for (std::size_t i = 0; i < soa.size(); ) {
        if (predicate(i)) soa.remove(i); else ++i;
      }
    @endcode
  append/remove/reserve invalidate raw pointers from column()/positionColumn(). Capacity
  never shrinks (clear() keeps it); the container is meant to be reused across steps.

  <h3>Cell sorting</h3>
  sortByCell(box,dx,probLo) counting-sorts the columns into Fortran cell order and builds
  a CSR offsets array, so cell c owns the contiguous range [cellStart(c), cellStart(c+1)).
  Append/remove invalidate the sort (isSorted()).

  <h3>Communication & I/O</h3>
  linearizeParticle()/delinearizeAndAppend() pack/unpack ALL columns (incl. id/rank) for
  MPI. h5LinearizeParticle()/h5DelinearizeAndAppend() handle the HDF5 subset (position +
  weight + the payload h5Columns; never id/rank). data()/byteSpan() expose the raw arena for
  a one-memcpy / zero-copy whole-container transfer.

  @tparam P      User payload struct (extra columns only); defaults to NoPayload.
  @tparam Traits Column descriptor for the payload; defaults to ParticleTraits<P>.
*/
template <typename P = NoPayload, typename Traits = ParticleTraits<P>>
class ParticleSoA
{
public:
  /**
    @name Compile-time layout (indices and sizes)
  */
  ///@{

  /**
     @brief The payload columns descriptor tuple type (tuple of pointer-to-member).
  */
  using PayloadColumns = std::remove_cv_t<std::remove_reference_t<decltype(Traits::columns)>>;

  /**
     @brief Number of user-declared payload columns.
  */
  static constexpr std::size_t s_numPayloadColumns = std::tuple_size<PayloadColumns>::value;

  /**
     @brief Column index of the weight (immediately after the SpaceDim position columns).
  */
  static constexpr std::size_t s_weightCol = SpaceDim;

  /**
     @brief Column index of the container-owned particleID.
  */
  static constexpr std::size_t s_idCol = SpaceDim + 1;

  /**
     @brief Column index of the container-owned rankID.
  */
  static constexpr std::size_t s_rankCol = SpaceDim + 2;

  /**
     @brief Column index of the container-owned ghost designation.
  */
  static constexpr std::size_t s_ghostCol = SpaceDim + 3;

  /**
     @brief Index of the first payload column.
  */
  static constexpr std::size_t s_payloadBegin = SpaceDim + 4;

  /**
     @brief Total number of columns (SpaceDim position + weight + id + rank + ghost + payload).
  */
  static constexpr std::size_t s_numColumns = s_payloadBegin + s_numPayloadColumns;

  /**
     @brief Cache-line alignment for each column (also enables aligned SIMD).
  */
  static constexpr std::size_t s_align = 64;

  ///@}

private:
  /**
     @brief The SpaceDim position columns, all double -- see precision rationale.
  */
  using PositionValues = typename detail::RepeatTuple<double, SpaceDim>::Type;

  /**
     @brief The mandatory column value types: positions (double), weight (double), id, rank, ghost.
  */
  using MandatoryValues = decltype(std::tuple_cat(std::declval<PositionValues>(),
                                                  std::declval<std::tuple<double, ParticleID, RankID, GhostType>>()));

  /**
     @brief The payload column value types.
  */
  using PayloadValues = typename detail::ValueTuple<PayloadColumns>::Type;

public:
  /**
    @name Compile-time layout (value types)
  */
  ///@{

  /**
     @brief The value types of every column, in storage order.
  */
  using ColumnValueTuple = decltype(std::tuple_cat(std::declval<MandatoryValues>(), std::declval<PayloadValues>()));

  /**
     @brief Value type of column K.
  */
  template <std::size_t K>
  using ColumnType = std::tuple_element_t<K, ColumnValueTuple>;

  /**
    @brief Derive the column index of a PAYLOAD member pointer (reorder-safe selector).
    @tparam M Pointer-to-member type of the payload field.
    @param[in] a_member Member pointer selecting the payload column.
    @return The storage column index of the selected payload field.
  */
  template <typename M>
  static constexpr std::size_t
  columnIndex(const M a_member) noexcept
  {
    return s_payloadBegin + detail::indexOf(Traits::columns, a_member);
  }

  /**
     @brief The HDF5 payload column subset (indices into the payload columns).
  */
  using H5PayloadColumns = typename detail::H5ColumnSubset<Traits, std::make_index_sequence<s_numPayloadColumns>>::Type;

  ///@}

  // ----- Sanity checks ------------------------------------------------------
  static_assert(std::is_default_constructible<P>::value, "payload type P must be default-constructible");
  static_assert(detail::AllTriviallyCopyable<ColumnValueTuple>::s_value,
                "every particle column must be trivially copyable for MPI/HDF5");

  /**
    @name Construction and lifetime
  */
  ///@{

  /**
     @brief Construct an empty container.
  */
  ParticleSoA() = default;

  /**
     @brief Destructor (frees the arena).
  */
  ~ParticleSoA()
  {
    std::free(m_buffer);
  }

  /**
     @brief Copy construction is deleted (move-only).
  */
  ParticleSoA(const ParticleSoA&) = delete;

  /**
    @brief Copy assignment is deleted (move-only).
    @return Reference to this.
  */
  ParticleSoA&
  operator=(const ParticleSoA&) = delete;

  /**
    @brief Move constructor (steals the arena).
    @param[in,out] a_other Source container, left empty after the move.
  */
  ParticleSoA(ParticleSoA&& a_other) noexcept
  {
    this->steal(a_other);
  }

  /**
    @brief Move assignment (steals the arena).
    @param[in,out] a_other Source container, left empty after the move.
    @return Reference to this.
  */
  ParticleSoA&
  operator=(ParticleSoA&& a_other) noexcept
  {
    if (this != &a_other) {
      std::free(m_buffer);
      this->steal(a_other);
    }
    return *this;
  }

  /**
    @brief Explicit deep copy: a new container with a private arena duplicating every particle.
    @details Copy construction and copy assignment are deleted (they would alias the single arena);
    call this when an independent copy of all the particles is genuinely wanted (rather than moving
    them). The cell-sort state (isSorted() + the CSR cell offsets) is preserved, so a cell-sorted
    container copies as cell-sorted.
    @return An independent ParticleSoA owning a copy of every particle.
  */
  [[nodiscard]] ParticleSoA
  deepCopy() const;

  /**
    @brief Deep copy into an existing destination, reusing its arena when possible.
    @details Like deepCopy() but writes into a_dst instead of returning a new container. a_dst's
    existing allocation is reused when it is already large enough (grown only if it is too small),
    which avoids a malloc when copying repeatedly into the same destination. a_dst is overwritten
    with an independent copy of every particle, and its cell-sort state is set to match this
    container. Copying onto self is a no-op.
    @param[in,out] a_dst Destination container, overwritten with a copy of every particle.
  */
  void
  deepCopyTo(ParticleSoA& a_dst) const;

  /**
    @brief Swap the entire state (arena, size, capacity, sort) with another container in O(1).
    @details Exchanges the cached column pointers and metadata only -- no particle data is moved.
    Useful for the build-then-swap idiom and for std::swap.
    @param[in,out] a_other Container to swap with.
  */
  void
  swap(ParticleSoA& a_other) noexcept
  {
    std::swap(m_buffer, a_other.m_buffer);
    std::swap(m_capacity, a_other.m_capacity);
    std::swap(m_size, a_other.m_size);
    std::swap(m_sorted, a_other.m_sorted);
    std::swap(m_colPtrs, a_other.m_colPtrs);
    m_cellStart.swap(a_other.m_cellStart);
  }

  ///@}

  /**
    @name Size and capacity
  */
  ///@{

  /**
    @brief Number of particles currently stored.
    @return The live particle count.
  */
  std::size_t
  size() const noexcept
  {
    return m_size;
  }

  /**
    @brief Allocated capacity in particles.
    @return The capacity in particles.
  */
  std::size_t
  capacity() const noexcept
  {
    return m_capacity;
  }

  /**
    @brief Whether the container holds no particles.
    @return True iff size() is zero.
  */
  bool
  empty() const noexcept
  {
    return m_size == 0;
  }

  /**
     @brief Drop all particles (keeps the arena; invalidates the cell sort).
  */
  void
  clear() noexcept
  {
    m_size   = 0;
    m_sorted = false;
  }

  /**
    @brief Ensure capacity for at least a_capacity particles (reallocates + moves on growth).
    @param[in] a_capacity Requested capacity.
  */
  void
  reserve(const std::size_t a_capacity);

  /**
    @brief Set the particle count (grows capacity if needed; new slots are uninitialized).
    @param[in] a_count New particle count.
  */
  void
  resize(const std::size_t a_count);

  /**
    @brief Reclaim unused capacity by reallocating the arena down to the current size.
    @details A no-op when already compact (capacity == size). The arena does not shrink on its own
    (clear()/remove() keep the capacity for reuse); call this to return memory after a transient
    spike. Invalidates raw column pointers. The particle order and cell-sort state are preserved.
  */
  void
  shrinkToFit();

  ///@}

  /**
    @name Population
  */
  ///@{

  /**
    @brief Append one particle with a default-constructed payload.
    @details The metadata columns (id/rank) are set to s_invalidID/-1; assign them via
    particleID(i)/rankID(i) afterwards, or let the AMR layer / delinearize set them.
    @param[in] a_position Particle position.
    @param[in] a_weight   Particle weight.
  */
  void
  append(const RealVect& a_position, const double a_weight)
  {
    this->append(a_position, a_weight, P{});
  }

  /**
    @brief Append one particle with the given payload.
    @details The metadata columns (id/rank) are set to s_invalidID/-1; assign them via
    particleID(i)/rankID(i) afterwards, or let the AMR layer / delinearize set them.
    @param[in] a_position Particle position.
    @param[in] a_weight   Particle weight.
    @param[in] a_payload  Payload columns for this particle.
  */
  void
  append(const RealVect& a_position, const double a_weight, const P& a_payload);

  /**
    @brief Append every particle of another container (bulk catenate; all columns, incl. id/rank).
    @details Copies a_other's live particles onto the end of this container with a per-column memcpy
    (one reservation, then column-wise bulk copies). Unlike the single-particle append(), this
    PRESERVES the particleID/rankID metadata as well as position/weight/payload. Invalidates the cell
    sort. Appending a container to itself is supported.
    @param[in] a_other Source container whose particles are appended.
  */
  void
  append(const ParticleSoA& a_other);

  /**
    @brief Append a single particle (all columns, incl. id/rank) copied from another container.
    @details The single-particle analogue of append(const ParticleSoA&). Pair it with remove() on the
    source to move one particle between containers. Invalidates the cell sort. Appending from self is
    supported.
    @param[in] a_src   Source container to copy from.
    @param[in] a_index Index of the particle in a_src to copy.
  */
  void
  appendParticle(const ParticleSoA& a_src, const std::size_t a_index);

  /**
    @brief Move every particle of another container into this one, leaving a_other empty (catenate).
    @details The move counterpart of append(const ParticleSoA&): appends all of a_other's particles
    (all columns, incl. id/rank) and then clears a_other. When this container is empty the two arenas
    are swapped instead, so the common "drain into an empty destination" case is O(1) with zero copy.
    Invalidates the cell sort. Catenating a container into itself is a no-op.
    @param[in,out] a_other Source container, left empty afterwards.
  */
  void
  catenate(ParticleSoA& a_other)
  {
    if (this == &a_other || a_other.m_size == 0) {
      return;
    }
    if (m_size == 0) {
      this->swap(a_other); // zero-copy steal into an empty destination
      return;
    }
    this->append(a_other); // bulk per-column copy
    a_other.clear();       // move semantics: leave the source empty
  }

  /**
    @brief Gather particle i's payload back into the AoS payload view.
    @param[in] a_index Particle index.
    @return The payload struct (mandatory fields are read via the dedicated accessors).
  */
  P
  gather(const std::size_t a_index) const
  {
    CH_assert(a_index < m_size);

    return this->gatherPayload(a_index, std::make_index_sequence<s_numPayloadColumns>{});
  }

  /**
    @brief Remove particle i using swap-and-pop (O(1), does NOT preserve order).
    @param[in] a_index Particle index.
  */
  void
  remove(const std::size_t a_index) noexcept
  {
    CH_assert(a_index < m_size);

    this->removeImpl(a_index, std::make_index_sequence<s_numColumns>{});
    m_size--;
    m_sorted = false;
  }

  ///@}

  /**
    @name Column and field access
  */
  ///@{

  /**
    @brief Typed base pointer of column K (low-level).
    @tparam K Column index.
    @return Base pointer of column K.
  */
  template <std::size_t K>
  ColumnType<K>*
  columnByIndex() noexcept
  {
    return reinterpret_cast<ColumnType<K>*>(m_colPtrs[K]);
  }

  /**
    @brief Const typed base pointer of column K (low-level).
    @tparam K Column index.
    @return Const base pointer of column K.
  */
  template <std::size_t K>
  const ColumnType<K>*
  columnByIndex() const noexcept
  {
    return reinterpret_cast<const ColumnType<K>*>(m_colPtrs[K]);
  }

  /**
    @brief Whole payload column selected by member pointer (raw pointer, vectorizable).
    @tparam Member Member pointer selecting the payload column.
    @return Base pointer of the selected payload column.
  */
  template <auto Member>
  auto*
  column() noexcept
  {
    return this->template columnByIndex<columnIndex(Member)>();
  }

  /**
    @brief Const whole payload column selected by member pointer.
    @tparam Member Member pointer selecting the payload column.
    @return Const base pointer of the selected payload column.
  */
  template <auto Member>
  const auto*
  column() const noexcept
  {
    return this->template columnByIndex<columnIndex(Member)>();
  }

  /**
    @brief One particle's payload field selected by member pointer.
    @tparam Member Member pointer selecting the payload column.
    @param[in] a_index Particle index.
    @return Reference to the selected payload field of particle a_index.
  */
  template <auto Member>
  auto&
  get(const std::size_t a_index) noexcept
  {
    CH_assert(a_index < m_size);

    return this->template columnByIndex<columnIndex(Member)>()[a_index];
  }

  /**
    @brief Const overload of get<&P::field>(i).
    @tparam Member Member pointer selecting the payload column.
    @param[in] a_index Particle index.
    @return Const reference to the selected payload field of particle a_index.
  */
  template <auto Member>
  const auto&
  get(const std::size_t a_index) const noexcept
  {
    CH_assert(a_index < m_size);

    return this->template columnByIndex<columnIndex(Member)>()[a_index];
  }

  /**
    @brief Raw position component column dir (double*, for SIMD kernels).
    @param[in] a_dir Spatial direction in [0, SpaceDim).
    @return Base pointer of the position component column for direction a_dir.
  */
  double*
  positionColumn(const int a_dir) noexcept
  {
    return reinterpret_cast<double*>(m_colPtrs[a_dir]);
  }

  /**
    @brief Const raw position component column dir.
    @param[in] a_dir Spatial direction in [0, SpaceDim).
    @return Const base pointer of the position component column for direction a_dir.
  */
  const double*
  positionColumn(const int a_dir) const noexcept
  {
    return reinterpret_cast<const double*>(m_colPtrs[a_dir]);
  }

  /**
    @brief Raw weight column (double*).
    @return Base pointer of the weight column.
  */
  double*
  weightColumn() noexcept
  {
    return this->template columnByIndex<s_weightCol>();
  }

  /**
    @brief Const raw weight column.
    @return Const base pointer of the weight column.
  */
  const double*
  weightColumn() const noexcept
  {
    return this->template columnByIndex<s_weightCol>();
  }

  ///@}

  /**
    @name Mandatory-field accessors
  */
  ///@{

  /**
    @brief Position of particle i as a RealVect (by value, assembled from the scalar columns).
    @param[in] a_index Particle index.
    @return The particle position for cell lookup / Chombo interop.
  */
  RealVect
  position(const std::size_t a_index) const noexcept
  {
    CH_assert(a_index < m_size);

    RealVect x;
    for (int dir = 0; dir < SpaceDim; dir++) {
      x[dir] = this->positionColumn(dir)[a_index];
    }
    return x;
  }

  /**
    @brief Set the position of particle i.
    @param[in] a_index    Particle index.
    @param[in] a_position New position.
  */
  void
  setPosition(const std::size_t a_index, const RealVect& a_position) noexcept
  {
    CH_assert(a_index < m_size);

    for (int dir = 0; dir < SpaceDim; dir++) {
      this->positionColumn(dir)[a_index] = a_position[dir];
    }
  }

  /**
    @brief Weight of particle i.
    @param[in] a_index Particle index.
    @return Reference to the weight.
  */
  double&
  weight(const std::size_t a_index) noexcept
  {
    CH_assert(a_index < m_size);

    return this->template columnByIndex<s_weightCol>()[a_index];
  }

  /**
    @brief Const weight of particle i.
    @param[in] a_index Particle index.
    @return Const reference to the weight.
  */
  const double&
  weight(const std::size_t a_index) const noexcept
  {
    CH_assert(a_index < m_size);

    return this->template columnByIndex<s_weightCol>()[a_index];
  }

  /**
    @brief Global id of particle i (container-owned metadata).
    @param[in] a_index Particle index.
    @return Reference to the particle id.
  */
  ParticleID&
  particleID(const std::size_t a_index) noexcept
  {
    CH_assert(a_index < m_size);

    return this->template columnByIndex<s_idCol>()[a_index];
  }

  /**
    @brief Const global id of particle i.
    @param[in] a_index Particle index.
    @return Const reference to the particle id.
  */
  const ParticleID&
  particleID(const std::size_t a_index) const noexcept
  {
    CH_assert(a_index < m_size);

    return this->template columnByIndex<s_idCol>()[a_index];
  }

  /**
    @brief Owning rank of particle i (container-owned metadata).
    @param[in] a_index Particle index.
    @return Reference to the owning rank.
  */
  RankID&
  rankID(const std::size_t a_index) noexcept
  {
    CH_assert(a_index < m_size);

    return this->template columnByIndex<s_rankCol>()[a_index];
  }

  /**
    @brief Const owning rank of particle i.
    @param[in] a_index Particle index.
    @return Const reference to the owning rank.
  */
  const RankID&
  rankID(const std::size_t a_index) const noexcept
  {
    CH_assert(a_index < m_size);

    return this->template columnByIndex<s_rankCol>()[a_index];
  }

  /**
    @brief Ghost designation of particle i (container-owned metadata).
    @param[in] a_index Particle index.
    @return Reference to the ghost designation.
  */
  GhostType&
  ghost(const std::size_t a_index) noexcept
  {
    CH_assert(a_index < m_size);

    return this->template columnByIndex<s_ghostCol>()[a_index];
  }

  /**
    @brief Const ghost designation of particle i.
    @param[in] a_index Particle index.
    @return Const reference to the ghost designation.
  */
  const GhostType&
  ghost(const std::size_t a_index) const noexcept
  {
    CH_assert(a_index < m_size);

    return this->template columnByIndex<s_ghostCol>()[a_index];
  }

  /**
    @brief Whether particle i is a ghost particle (any non-Valid designation).
    @param[in] a_index Particle index.
    @return True if the particle is a ghost, false if it is a valid/owned particle.
  */
  bool
  isGhost(const std::size_t a_index) const noexcept
  {
    return this->ghost(a_index) != GhostType::Valid;
  }

  ///@}

  /**
    @name MPI linearization
  */
  ///@{

  /**
    @brief Byte count of one fully-linearized particle (all columns, for MPI).
    @return Bytes per particle.
  */
  static constexpr std::size_t
  bytesPerParticle() noexcept
  {
    return detail::SumSizeof<ColumnValueTuple>::s_value;
  }

  /**
    @brief Linearize particle i into a byte buffer (all columns, for MPI send).
    @param[in,out] a_buffer Buffer with at least bytesPerParticle() bytes.
    @param[in]     a_index  Particle index.
  */
  void
  linearizeParticle(void* a_buffer, const std::size_t a_index) const noexcept;

  /**
    @brief Append a particle delinearized from a byte buffer (all columns, for MPI receive).
    @param[in] a_buffer Buffer with at least bytesPerParticle() bytes.
  */
  void
  delinearizeAndAppend(const void* a_buffer);

  ///@}

  /**
    @name HDF5 linearization
  */
  ///@{

  /**
    @brief Bytes for one HDF5-linearized particle (position + weight + h5 payload subset).
    @return Bytes per HDF5-checkpointed particle.
  */
  static constexpr std::size_t
  h5BytesPerParticle() noexcept
  {
    return (SpaceDim + 1) * sizeof(double) + h5PayloadSize(H5PayloadColumns{});
  }

  /**
    @brief Linearize the HDF5-checkpointed columns of particle i (no id/rank).
    @param[in,out] a_buffer Buffer with at least h5BytesPerParticle() bytes.
    @param[in]     a_index  Particle index.
  */
  void
  h5LinearizeParticle(void* a_buffer, const std::size_t a_index) const noexcept;

  /**
    @brief Append a particle whose HDF5-checkpointed columns come from a buffer.
    @details id/rank are set invalid and non-checkpointed payload columns default-constructed.
    @param[in] a_buffer Buffer with at least h5BytesPerParticle() bytes.
  */
  void
  h5DelinearizeAndAppend(const void* a_buffer);

  ///@}

  /**
    @name Raw arena and bulk transfer
  */
  ///@{

  /**
    @brief Raw arena pointer (column-major); a compact container's whole data span.
    @return Const pointer to the arena.
  */
  const void*
  data() const noexcept
  {
    return m_buffer;
  }

  /**
    @brief Mutable raw arena pointer (e.g. zero-copy MPI receive into the arena).
    @return Pointer to the arena.
  */
  void*
  data() noexcept
  {
    return m_buffer;
  }

  /**
    @brief Byte span of the live data (for a compact container, one contiguous memcpy/MPI_Send).
    @return Byte span [data(), data()+byteSpan()).
  */
  std::size_t
  byteSpan() const noexcept
  {
    if (m_size == 0) {
      return 0;
    }
    return this->byteSpanImpl(std::make_index_sequence<s_numColumns>{});
  }

  ///@}

  /**
    @name Cell sorting (CSR)
  */
  ///@{

  /**
    @brief Counting-sort the columns into Fortran cell order and build CSR cell offsets.
    @details After this, cell c owns the contiguous particle range [cellStart(c), cellStart(c+1)).
    Particles must lie within a_box; out-of-range positions are clamped to the boundary cells.
    @param[in] a_box    Patch box defining the cells.
    @param[in] a_dx     Grid spacing per direction.
    @param[in] a_probLo Physical coordinate of the lower domain corner.
  */
  void
  sortByCell(const Box& a_box, const RealVect& a_dx, const RealVect& a_probLo);

  /**
    @brief Whether the container is currently cell-sorted (false after append/remove).
    @return True iff the CSR cell offsets are valid.
  */
  bool
  isSorted() const noexcept
  {
    return m_sorted;
  }

  /**
    @brief Number of cells in the CSR sort (0 if not sorted).
    @return The number of cells.
  */
  std::size_t
  numCells() const noexcept
  {
    return m_cellStart.empty() ? 0 : m_cellStart.size() - 1;
  }

  /**
    @brief CSR start offset of cell c (valid after sortByCell; requires c <= numCells()).
    @param[in] a_cell Cell index.
    @return Index of the first particle owned by cell a_cell.
  */
  std::size_t
  cellStart(const std::size_t a_cell) const noexcept
  {
    CH_assert(m_sorted);
    CH_assert(a_cell <= this->numCells());

    return m_cellStart[a_cell];
  }

  /**
    @brief Number of particles in cell c (valid after sortByCell).
    @param[in] a_cell Cell index.
    @return Particle count in cell a_cell.
  */
  std::size_t
  particlesInCell(const std::size_t a_cell) const noexcept
  {
    CH_assert(m_sorted);
    CH_assert(a_cell < this->numCells());

    return m_cellStart[a_cell + 1] - m_cellStart[a_cell];
  }

  /**
    @brief Half-open particle index range [begin, end) owned by cell c (valid after sortByCell).
    @details Convenience over cellStart(c)/cellStart(c+1); the cell's particles are the contiguous
    indices [first, second).
    @param[in] a_cell Cell index.
    @return The pair {cellStart(c), cellStart(c+1)}.
  */
  std::pair<std::size_t, std::size_t>
  cellRange(const std::size_t a_cell) const noexcept
  {
    CH_assert(m_sorted);
    CH_assert(a_cell < this->numCells());

    return {m_cellStart[a_cell], m_cellStart[a_cell + 1]};
  }

  /**
    @brief Extract cell c's particles (the CSR range [cellStart(c), cellStart(c+1))) into a_out.
    @details a_out is cleared first, then this cell's particles are appended (all columns, including
    id/rank). Requires the container to be cell-sorted. This is the building block of the per-cell
    extract->process->rebuild idiom that replaces the AoS BinFab<List<P>> per-cell pattern: extract a
    cell into a reusable scratch container, mutate it (append/remove), accumulate it into an output
    container via append(const ParticleSoA&), and finally swap() the output into the leaf. See
    ItoSolver::makeSuperparticles for the canonical use.
    @param[in]  a_cell Cell index.
    @param[out] a_out  Destination scratch container (cleared on entry; reuse it across cells).
  */
  void
  extractCell(const std::size_t a_cell, ParticleSoA& a_out) const
  {
    CH_assert(m_sorted);

    a_out.clear();

    const std::pair<std::size_t, std::size_t> range = this->cellRange(a_cell);
    for (std::size_t i = range.first; i < range.second; i++) {
      a_out.appendParticle(*this, i);
    }
  }

  ///@}

protected:
  /**
    @brief The single arena allocation.
  */
  void* m_buffer = nullptr;

  /**
    @brief Capacity in particles.
  */
  std::size_t m_capacity = 0;

  /**
    @brief Live particle count.
  */
  std::size_t m_size = 0;

  /**
    @brief Whether m_cellStart is valid.
  */
  bool m_sorted = false;

  /**
    @brief Cached column base pointers into m_buffer.
  */
  std::array<void*, s_numColumns> m_colPtrs{};

  /**
    @brief CSR cell offsets (size numCells+1 when sorted).
  */
  std::vector<std::size_t> m_cellStart;

  /**
    @brief Sentinel for an unassigned global particle id.
  */
  static constexpr ParticleID s_invalidID = -1;

  /**
    @brief Round a_x up to the next multiple of s_align.
    @param[in] a_x Value to round up.
    @return a_x rounded up to a multiple of s_align.
  */
  static std::size_t
  alignUp(const std::size_t a_x) noexcept
  {
    return (a_x + (s_align - 1)) / s_align * s_align;
  }

  /**
    @brief Steal another container's arena (move helper).
    @param[in,out] a_o Source container, left empty afterwards.
  */
  void
  steal(ParticleSoA& a_o) noexcept
  {
    m_buffer       = a_o.m_buffer;
    m_capacity     = a_o.m_capacity;
    m_size         = a_o.m_size;
    m_sorted       = a_o.m_sorted;
    m_colPtrs      = a_o.m_colPtrs;
    m_cellStart    = std::move(a_o.m_cellStart);
    a_o.m_buffer   = nullptr;
    a_o.m_capacity = 0;
    a_o.m_size     = 0;
    a_o.m_sorted   = false;
    a_o.m_colPtrs  = {};
  }

  /** @brief Grow geometrically if the container is full (called before writing slot m_size). */
  void
  growIfFull()
  {
    if (m_size == m_capacity) {
      // Saturating geometric growth: 2 * m_capacity overflows std::size_t once m_capacity exceeds
      // SIZE_MAX/2, and reserve()'s "a_capacity <= m_capacity" guard would then silently no-op.
      const std::size_t newCap = (m_capacity == 0) ? 16
                                                   : (m_capacity <= std::numeric_limits<std::size_t>::max() / 2
                                                        ? 2 * m_capacity
                                                        : std::numeric_limits<std::size_t>::max());
      this->reserve(newCap);
    }
  }

  /**
    @brief Compute per-column byte offsets for a given capacity; returns the total bytes.
    @tparam I Column indices [0, s_numColumns).
    @param[in]  a_capacity Capacity (in particles) the columns are sized for.
    @param[out] a_off      Filled with the byte offset of each column.
    @return Total byte size of the arena for this capacity.
  */
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

  /**
    @brief Cache the column base pointers from a fresh buffer + offsets.
    @tparam I Column indices [0, s_numColumns).
    @param[in] a_off Per-column byte offsets into m_buffer.
  */
  template <std::size_t... I>
  void
  setColPtrs(const std::array<std::size_t, s_numColumns>& a_off, std::index_sequence<I...>) noexcept
  {
    using expander = int[];
    (void)expander{0, (m_colPtrs[I] = static_cast<char*>(m_buffer) + a_off[I], 0)...};
  }

  /**
    @brief Contiguously memcpy every column's live region into a freshly-allocated buffer.
    @tparam I Column indices [0, s_numColumns).
    @param[out] a_newBuf Destination arena buffer.
    @param[in]  a_off    Per-column byte offsets into a_newBuf.
  */
  template <std::size_t... I>
  void
  moveColumns(void* a_newBuf, const std::array<std::size_t, s_numColumns>& a_off, std::index_sequence<I...>) noexcept
  {
    using expander = int[];
    (void)expander{
      0,
      (std::memcpy(static_cast<char*>(a_newBuf) + a_off[I], m_colPtrs[I], m_size * sizeof(ColumnType<I>)), 0)...};
  }

  /**
    @brief Copy this container's m_size live elements of every column into another container.
    @details Per-column copy (the two arenas may have different capacities, hence different layouts,
    so a single whole-arena memcpy is not valid). a_dst must already have capacity >= m_size.
    @tparam I Column indices [0, s_numColumns).
    @param[in,out] a_dst Destination container (its column pointers are written).
  */
  template <std::size_t... I>
  void
  copyLiveColumnsTo(ParticleSoA& a_dst, std::index_sequence<I...>) const noexcept
  {
    using expander = int[];
    (void)expander{0, (std::memcpy(a_dst.m_colPtrs[I], m_colPtrs[I], m_size * sizeof(ColumnType<I>)), 0)...};
  }

  /**
    @brief Bulk-append another container's live columns at slot a_at of this container's columns.
    @details Per-column memcpy of a_other's m_size live elements into this container starting at
    element index a_at. This container must already have capacity >= a_at + a_other.size(). Disjoint
    source/destination (even when a_other == this, since the new slots lie past the existing ones).
    @tparam I Column indices [0, s_numColumns).
    @param[in] a_other Source container.
    @param[in] a_at    First destination element index (the size before appending).
  */
  template <std::size_t... I>
  void
  appendColumnsFrom(const ParticleSoA& a_other, const std::size_t a_at, std::index_sequence<I...>) noexcept
  {
    using expander = int[];
    (void)expander{0,
                   (std::memcpy(static_cast<char*>(m_colPtrs[I]) + a_at * sizeof(ColumnType<I>),
                                a_other.m_colPtrs[I],
                                a_other.m_size * sizeof(ColumnType<I>)),
                    0)...};
  }

  /**
    @brief Copy one particle's every column from a_src[a_srcIdx] into this[a_dstIdx].
    @details Per-column single-element copy. Both indices must be in range and this must have capacity
    for a_dstIdx.
    @tparam I Column indices [0, s_numColumns).
    @param[in] a_src    Source container.
    @param[in] a_srcIdx Source particle index.
    @param[in] a_dstIdx Destination particle index in this container.
  */
  template <std::size_t... I>
  void
  copyParticleColumns(const ParticleSoA& a_src,
                      const std::size_t  a_srcIdx,
                      const std::size_t  a_dstIdx,
                      std::index_sequence<I...>) noexcept
  {
    using expander = int[];
    (void)expander{0,
                   (this->template columnByIndex<I>()[a_dstIdx] = a_src.template columnByIndex<I>()[a_srcIdx], 0)...};
  }

  /**
    @brief Allocate a new arena and reallocate/move the columns into it (growth path).
    @param[in] a_capacity New capacity in particles.
  */
  void
  reallocate(const std::size_t a_capacity);

  /**
    @brief Reorder every column in place by a_dest (a_dest[i] = new slot of old particle i).
    @param[in] a_dest Destination slot for each current particle.
  */
  void
  permute(const std::vector<std::size_t>& a_dest);

  /**
    @brief Per-column scatter for permute(): newcol[dest[i]] = oldcol[i].
    @tparam K Column index to scatter.
    @param[out] a_newBuf Destination arena buffer.
    @param[in]  a_off    Byte offset of column K in a_newBuf.
    @param[in]  a_dest   Destination slot for each current particle.
  */
  template <std::size_t K>
  void
  scatterColumn(void* a_newBuf, const std::size_t a_off, const std::vector<std::size_t>& a_dest) noexcept
  {
    const ColumnType<K>* src = this->template columnByIndex<K>();
    ColumnType<K>*       dst = reinterpret_cast<ColumnType<K>*>(static_cast<char*>(a_newBuf) + a_off);
    for (std::size_t i = 0; i < m_size; i++) {
      dst[a_dest[i]] = src[i];
    }
  }

  /**
    @brief Scatter every column for permute() into a freshly-allocated buffer.
    @tparam I Column indices [0, s_numColumns).
    @param[out] a_newBuf Destination arena buffer.
    @param[in]  a_off    Per-column byte offsets into a_newBuf.
    @param[in]  a_dest   Destination slot for each current particle.
  */
  template <std::size_t... I>
  void
  scatterColumns(void*                                        a_newBuf,
                 const std::array<std::size_t, s_numColumns>& a_off,
                 const std::vector<std::size_t>&              a_dest,
                 std::index_sequence<I...>) noexcept
  {
    using expander = int[];
    (void)expander{0, ((void)this->template scatterColumn<I>(a_newBuf, a_off[I], a_dest), 0)...};
  }

  /**
    @brief Write the payload columns of slot m_size from a payload struct.
    @tparam J Payload-column indices [0, s_numPayloadColumns).
    @param[in] a_p Payload struct to read from.
  */
  template <std::size_t... J>
  void
  writePayload(const P& a_p, std::index_sequence<J...>) noexcept
  {
    constexpr auto cols = Traits::columns;
    using expander      = int[];
    (void)expander{0,
                   ((void)(this->template columnByIndex<s_payloadBegin + J>()[m_size] = a_p.*std::get<J>(cols)), 0)...};
  }

  /**
    @brief Read the payload columns of slot a_index into a payload struct.
    @tparam J Payload-column indices [0, s_numPayloadColumns).
    @param[in] a_index Particle index.
    @return The gathered payload struct.
  */
  template <std::size_t... J>
  P
  gatherPayload(const std::size_t a_index, std::index_sequence<J...>) const
  {
    constexpr auto cols = Traits::columns;
    P              p;
    using expander = int[];
    (void)expander{0,
                   ((void)(p.*std::get<J>(cols) = this->template columnByIndex<s_payloadBegin + J>()[a_index]), 0)...};
    return p;
  }

  /**
    @brief Swap-and-pop slot a_index in every column.
    @tparam I Column indices [0, s_numColumns).
    @param[in] a_index Particle index to remove.
  */
  template <std::size_t... I>
  void
  removeImpl(const std::size_t a_index, std::index_sequence<I...>) noexcept
  {
    using expander = int[];
    (void)expander{
      0,
      ((void)(this->template columnByIndex<I>()[a_index] = this->template columnByIndex<I>()[m_size - 1]), 0)...};
  }

  /**
    @brief Push every column of particle a_index onto a byte buffer (MPI).
    @tparam I Column indices [0, s_numColumns).
    @param[in,out] a_buf   Buffer pointer, advanced as columns are written.
    @param[in]     a_index Particle index.
  */
  template <std::size_t... I>
  void
  linearizeImpl(unsigned char* a_buf, const std::size_t a_index, std::index_sequence<I...>) const noexcept
  {
    using expander = int[];
    (void)expander{0, ((void)detail::pushBytes(a_buf, this->template columnByIndex<I>()[a_index]), 0)...};
  }

  /**
    @brief Pull every column of one particle off a byte buffer into slot m_size (MPI).
    @tparam I Column indices [0, s_numColumns).
    @param[in] a_buf Buffer holding one linearized particle.
  */
  template <std::size_t... I>
  void
  delinearizeImpl(const unsigned char* a_buf, std::index_sequence<I...>) noexcept
  {
    using expander = int[];
    (void)expander{0, ((void)this->template pullColumn<I>(a_buf), 0)...};
  }

  /**
    @brief Pull one column's value off a byte buffer into slot m_size and advance the pointer.
    @tparam K Column index.
    @param[in,out] a_buf Buffer pointer, advanced by the column's byte size.
  */
  template <std::size_t K>
  void
  pullColumn(const unsigned char*& a_buf) noexcept
  {
    detail::pullBytes(a_buf, this->template columnByIndex<K>()[m_size]);
  }

  /**
    @brief Compile-time byte size of the HDF5 payload subset.
    @tparam J Payload-column indices in the HDF5 subset.
    @return Total bytes of the HDF5 payload subset.
  */
  template <std::size_t... J>
  static constexpr std::size_t
  h5PayloadSize(std::index_sequence<J...>) noexcept
  {
    return detail::sumSizes(sizeof(ColumnType<s_payloadBegin + J>)...);
  }

  /**
    @brief Push the HDF5 payload subset of particle a_index onto a byte buffer.
    @tparam J Payload-column indices in the HDF5 subset.
    @param[in,out] a_buf   Buffer pointer, advanced as columns are written.
    @param[in]     a_index Particle index.
  */
  template <std::size_t... J>
  void
  h5LinearizePayload(unsigned char*& a_buf, const std::size_t a_index, std::index_sequence<J...>) const noexcept
  {
    using expander = int[];
    (void)expander{0,
                   ((void)detail::pushBytes(a_buf, this->template columnByIndex<s_payloadBegin + J>()[a_index]), 0)...};
  }

  /**
    @brief Pull the HDF5 payload subset off a byte buffer into slot m_size.
    @tparam J Payload-column indices in the HDF5 subset.
    @param[in,out] a_buf Buffer pointer, advanced as columns are read.
  */
  template <std::size_t... J>
  void
  h5DelinearizePayload(const unsigned char*& a_buf, std::index_sequence<J...>) noexcept
  {
    using expander = int[];
    (void)expander{0,
                   ((void)detail::pullBytes(a_buf, this->template columnByIndex<s_payloadBegin + J>()[m_size]), 0)...};
  }

  /**
    @brief Byte span of the live data: max over columns of (offset + used bytes).
    @tparam I Column indices [0, s_numColumns).
    @return Byte span of the live region of the arena.
  */
  template <std::size_t... I>
  std::size_t
  byteSpanImpl(std::index_sequence<I...>) const noexcept
  {
    std::size_t span = 0;
    using expander   = int[];
    (void)expander{0,
                   (span = (span > (static_cast<std::size_t>(static_cast<const char*>(m_colPtrs[I]) -
                                                             static_cast<const char*>(m_buffer)) +
                                    m_size * sizeof(ColumnType<I>)))
                             ? span
                             : (static_cast<std::size_t>(static_cast<const char*>(m_colPtrs[I]) -
                                                         static_cast<const char*>(m_buffer)) +
                                m_size * sizeof(ColumnType<I>)),
                    0)...};
    return span;
  }
};

#include <CD_NamespaceFooter.H>

#include <CD_ParticleSoAImplem.H>

#endif
