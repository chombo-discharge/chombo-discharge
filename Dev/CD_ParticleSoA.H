/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ParticleSoA.H
  @brief  Declaration of a Struct-of-Arrays (SoA) particle storage container.
  @author Robert Marskar
  @details A Struct-of-Arrays replacement for Chombo's List<P> as the per-patch
           particle storage. The central idea:

           1. The user defines a plain "Array-of-Structs" particle struct, e.g.
                struct MyParticle { RealVect pos; Real weight; Real lightsaberDamage; };
              This is the ergonomic single-particle view the rest of the code uses.

           2. The user specializes ParticleTraits<MyParticle> to describe the SoA
              columns as a compile-time tuple of pointer-to-member, and to tag which
              column is the position and which is the weight.

           3. ParticleSoA<MyParticle> then stores one contiguous std::vector per
              column (true SoA), but can still hand back / accept a MyParticle via
              gather()/append() for code that prefers the per-particle view.

           The two mandatory fields are:
             - position(i) -> RealVect& (the only mandatory vector column)
             - weight(i)   -> Real&     (the only mandatory scalar column)
           Everything else (velocity, energy, lightsaberDamage, ...) is generic
           user-defined columns.
*/

#ifndef CD_PARTICLESOA_H
#define CD_PARTICLESOA_H

// Std includes
#include <cstring>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

// Chombo includes
#include <REAL.H>
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Traits class that a user specializes to teach ParticleSoA how to
         decompose their particle struct into SoA columns.
  @details
  <h3>Defining a new particle type (the complete recipe)</h3>

  Step 1 -- write a plain struct with public data members. Each member becomes one
  SoA column. The two mandatory fields are a RealVect position and a Real weight;
  everything else is up to you. Every member must be byte-copyable (memcpy-able)
  for MPI/HDF5: trivially-copyable types qualify, and RealVect is whitelisted
  (see detail::IsBlittable); pointers and types with custom copy semantics do not.

    @code
      struct MyParticle
      {
        RealVect pos      = RealVect::Zero;
        Real     weight   = 0.0;
        RealVect velocity = RealVect::Zero;
        Real     energy   = 0.0;
      };
    @endcode

  Step 2 -- specialize ParticleTraits<MyParticle> (in the same namespace as this
  template) with three required members:

    - `columns`     : a std::make_tuple of pointer-to-member, ONE per SoA column.
                      This tuple is the single source of truth for the layout.
    - `positionPtr` : pointer to the RealVect position member.
    - `weightPtr`   : pointer to the Real weight member.

    @code
      template <>
      struct ParticleTraits<MyParticle>
      {
        static constexpr auto columns     = std::make_tuple(&MyParticle::pos,
                                                            &MyParticle::weight,
                                                            &MyParticle::velocity,
                                                            &MyParticle::energy);
        static constexpr auto positionPtr = &MyParticle::pos;
        static constexpr auto weightPtr   = &MyParticle::weight;
      };
    @endcode

  Step 3 (optional) -- restrict what gets written to HDF5 checkpoints. Declare an
  `h5Columns` tuple of member pointers naming the subset to persist; fields left
  out are recomputed after restart (and come back default-constructed). If
  `h5Columns` is omitted, ALL columns are checkpointed.

    @code
      static constexpr auto h5Columns = std::make_tuple(&MyParticle::pos,
                                                        &MyParticle::weight,
                                                        &MyParticle::energy);
    @endcode

  That is all -- `ParticleSoA<MyParticle>` now works, and call sites select fields
  by member pointer, e.g. `soa.column<&MyParticle::velocity>()`.

  <h3>Why member pointers instead of indices</h3>
  Position/weight/HDF5 columns are designated by MEMBER POINTER, never by a
  hand-written integer index. ParticleSoA derives the column indices from the
  pointers at compile time, so reordering the `columns` tuple (e.g. to tune memory
  layout) keeps every name-based access correct -- there is no second list to keep
  in sync. The one externally visible effect of reordering is that the serialized
  byte order changes, so treat a layout change as a checkpoint-format change.

  @note The primary template is intentionally left undefined, so using ParticleSoA
        with a particle type that has no ParticleTraits specialization is a hard
        compile error rather than a confusing template failure deep inside.
  @tparam P The user particle struct being described.
*/
template <typename P>
struct ParticleTraits;

namespace detail {

  /**
    @brief Extract the class type and value type from a pointer-to-member type.
  */
  template <typename M>
  struct MemberPointerTraits;

  /**
    @brief Trait: true if T may be safely byte-copied (memcpy) for linearization.
    @details std::is_trivially_copyable is too strict for Chombo's RealVect, which
    has a user-defined copy constructor (so it is not "trivially copyable") yet is
    a flat array of Real and is safe to blit. RealVect is whitelisted here.
  */
  template <typename T>
  struct IsBlittable : std::integral_constant<bool, std::is_trivially_copyable<T>::value>
  {
  };
  template <>
  struct IsBlittable<RealVect> : std::true_type
  {
  };

  template <typename C, typename T>
  struct MemberPointerTraits<T C::*>
  {
    using ClassType = C;
    using ValueType = T;
  };

  /**
    @brief Type-guarded equality of two pointer-to-member. Returns false when the
           member types differ (so heterogeneous tuples compare cleanly).
  */
  template <typename A, typename B>
  constexpr bool
  eqMember(A, B) noexcept
  {
    return false;
  }
  template <typename A>
  constexpr bool
  eqMember(A a_lhs, A a_rhs) noexcept
  {
    return a_lhs == a_rhs;
  }

  /**
    @brief Compile-time search: index of pointer-to-member a_target within a tuple
           of pointer-to-member. Returns tuple size if not found.
    @details This is what makes the layout reorder-safe: a field's column index is
    DERIVED from its member pointer rather than hand-numbered, so swapping two
    columns does not silently misroute name-based access.
  */
  template <typename Tup, typename M>
  constexpr std::size_t
  indexOfImpl(const Tup&, M, std::index_sequence<>) noexcept
  {
    return std::tuple_size<Tup>::value;
  }
  template <typename Tup, typename M, std::size_t I, std::size_t... Rest>
  constexpr std::size_t
  indexOfImpl(const Tup& a_t, M a_target, std::index_sequence<I, Rest...>) noexcept
  {
    return eqMember(std::get<I>(a_t), a_target) ? I : indexOfImpl(a_t, a_target, std::index_sequence<Rest...>{});
  }
  template <typename Tup, typename M>
  constexpr std::size_t
  indexOf(const Tup& a_t, M a_target) noexcept
  {
    return indexOfImpl(a_t, a_target, std::make_index_sequence<std::tuple_size<Tup>::value>{});
  }

  /**
    @brief Build the SoA storage type (a tuple of std::vector<columnType>) from
           the tuple of pointer-to-member types.
  */
  template <typename Tuple>
  struct SoaStorage;

  template <typename... Ms>
  struct SoaStorage<std::tuple<Ms...>>
  {
    using Type = std::tuple<std::vector<typename MemberPointerTraits<Ms>::ValueType>...>;
  };

  /** @brief Variadic constexpr sum of std::size_t values. */
  constexpr std::size_t
  sumSizes() noexcept
  {
    return 0;
  }
  template <typename... Rest>
  constexpr std::size_t
  sumSizes(std::size_t a_first, Rest... a_rest) noexcept
  {
    return a_first + sumSizes(a_rest...);
  }

  /**
    @brief Compile-time sum of sizeof() of every column value type. Used as the
           per-particle byte count when linearizing for MPI.
  */
  template <typename Tuple>
  struct LinearSize;

  template <typename... Ms>
  struct LinearSize<std::tuple<Ms...>>
  {
    static constexpr std::size_t s_value = sumSizes(sizeof(typename MemberPointerTraits<Ms>::ValueType)...);
  };

  /**
    @brief True iff every column value type is byte-copyable / blittable (required for the
           memcpy-based MPI/HDF5 linearization).
  */
  template <typename Tuple>
  struct AllBlittable;

  template <typename... Ms>
  struct AllBlittable<std::tuple<Ms...>>
  {
    static constexpr bool s_value = (IsBlittable<typename MemberPointerTraits<Ms>::ValueType>::value && ...);
  };

  /** @brief Stand-in for std::void_t. */
  template <typename...>
  struct MakeVoid
  {
    using Type = void;
  };
  template <typename... Ts>
  using Void = typename MakeVoid<Ts...>::Type;

  /**
    @brief Turn a user tuple of HDF5 member pointers into a std::index_sequence of
           DERIVED column indices (each found via indexOf against the full column
           tuple). Only the return TYPE is used (via decltype); never called.
  */
  template <typename Traits, std::size_t... I>
  constexpr std::index_sequence<indexOf(Traits::columns, std::get<I>(Traits::h5Columns))...>
  deriveH5Seq(std::index_sequence<I...>) noexcept
  {
    return {};
  }

  /**
    @brief Resolve the HDF5 column subset as an index_sequence.
    @details If the user's traits declare `static constexpr auto h5Columns =
    std::make_tuple(&P::a, &P::b, ...)` (member pointers, same single source of
    truth as `columns`), the indices are derived from those pointers - so
    reordering `columns` keeps the right fields checkpointed. If no `h5Columns` is
    declared, this falls back to "all columns". Mirrors GenericParticle's
    H5size/H5linearOut/H5linearIn (which dropped IDs and could drop scratch),
    while staying generic and reorder-safe.
  */
  template <typename Traits, typename Fallback, typename = void>
  struct H5ColumnSubset
  {
    using Type = Fallback;
  };
  template <typename Traits, typename Fallback>
  struct H5ColumnSubset<Traits, Fallback, Void<decltype(Traits::h5Columns)>>
  {
    using H5Tuple = typename std::remove_cv<typename std::remove_reference<decltype(Traits::h5Columns)>::type>::type;
    using Type    = decltype(deriveH5Seq<Traits>(std::make_index_sequence<std::tuple_size<H5Tuple>::value>{}));
  };

  /**
    @brief memcpy one trivially-copyable value onto a byte buffer and advance it.
  */
  template <typename T>
  inline void
  pushBytes(unsigned char*& a_buf, const T& a_value) noexcept
  {
    static_assert(IsBlittable<T>::value, "SoA columns must be byte-copyable (trivially copyable or RealVect)");
    std::memcpy(a_buf, &a_value, sizeof(T));
    a_buf += sizeof(T);
  }

  /**
    @brief memcpy one trivially-copyable value off a byte buffer and advance it.
  */
  template <typename T>
  inline void
  pullBytes(const unsigned char*& a_buf, T& a_value) noexcept
  {
    static_assert(IsBlittable<T>::value, "SoA columns must be byte-copyable (trivially copyable or RealVect)");
    std::memcpy(&a_value, a_buf, sizeof(T));
    a_buf += sizeof(T);
  }

} // namespace detail

/**
  @brief Struct-of-Arrays (SoA) particle container for a single grid patch.
  @details
  One ParticleSoA owns all the particles of one grid patch. It is the SoA
  replacement for the `List<P>` that today lives inside a `BinFab<P>`; the AMR
  machinery (`ParticleContainer`) owns a `LayoutData<ParticleSoA<P>>` per level.

  <h3>Storage model</h3>
  Instead of one array of particle structs (AoS), the data is stored
  column-major: one contiguous `std::vector` per particle field ("column"). For a
  particle with members {pos, weight, velocity} the storage is
  `{ vector<RealVect> pos; vector<Real> weight; vector<RealVect> velocity; }`.
  All columns have the same length, equal to size(). This layout lets per-field
  loops stream contiguous memory and auto-vectorize.

  The user never lists columns by hand here -- the particle type `P` and its
  ParticleTraits<P> specialization describe them (see ParticleTraits for the full
  recipe for defining a new particle type).

  <h3>How to access particle data</h3>
  Two complementary styles, both reorder-safe (a field is named by its member
  pointer; the column index is derived at compile time):

    1. Column access (preferred, fast). Operate on a whole field at once:
       @code
         auto&       w = soa.column<&P::weight>();     // std::vector<Real>&
         const auto& m = soa.column<&P::mobility>();
         for (std::size_t i = 0; i < soa.size(); ++i) w[i] *= m[i];
       @endcode
       column<&P::field>() returns the std::vector; .data() gives the raw pointer
       that deposition/interpolation kernels vectorize over. position(i)/weight(i)
       and positions()/weights() are convenience shortcuts for the two mandatory
       fields. get<&P::field>(i) returns one particle's field by reference.

    2. Whole-particle (AoS) view, for occasional per-particle work (splitting,
       MPI pack). gather(i) copies particle i OUT into a `P`; scatter(i, p) writes
       a `P` back IN. These COPY -- never use them in hot loops; use column access.

  <h3>Mutating the population</h3>
  append(p) adds a particle; remove(i) deletes particle i by swap-and-pop (O(1),
  does NOT preserve order -- order is not semantically meaningful). Because
  remove() swaps the last element into slot i, the canonical "iterate and remove"
  loop does not advance i on a removal:
    @code
      for (std::size_t i = 0; i < soa.size(); ) {
        if (predicate(i)) soa.remove(i);   // last element now at i
        else              ++i;
      }
    @endcode
  Note that append()/remove()/resize() invalidate references obtained from
  column()/get() (same hazard as iterators into a std::vector).

  <h3>Communication & I/O</h3>
  linearizeParticle()/delinearizeAndAppend() pack/unpack a particle to/from a byte
  buffer for MPI (all columns). h5LinearizeParticle()/h5DelinearizeAndAppend()
  do the same for the optional HDF5 column subset (see ParticleTraits::h5Columns).
  Both are generic -- no per-type linearization code is required.

  <h3>Performance contract</h3>
  - Hot loops should use column access and read the raw pointers from .data().
  - gather()/scatter() copy a whole particle and must stay out of hot loops.
  - Field selection (column<&P::f>(), get<&P::f>(), deposit<&P::f>()) is resolved
    entirely at compile time -- zero runtime overhead vs. a hand-written index.

  @note ParticleSoA is move-only (copying a patch would deep-copy every column;
        we forbid it to avoid accidental O(N) copies, matching ParticleContainer).

  @tparam P      User particle struct (the AoS single-particle view).
  @tparam Traits Column descriptor; defaults to ParticleTraits<P>.
*/
template <typename P, typename Traits = ParticleTraits<P>>
class ParticleSoA
{
public:
  /** @brief The columns descriptor tuple type (tuple of pointer-to-member). */
  using ColumnsTuple = std::remove_cv_t<std::remove_reference_t<decltype(Traits::columns)>>;

  /** @brief The SoA storage type (tuple of std::vector per column). */
  using Storage = typename detail::SoaStorage<ColumnsTuple>::Type;

  /** @brief Number of SoA columns. */
  static constexpr std::size_t s_numColumns = std::tuple_size<ColumnsTuple>::value;

  /** @brief Value type of column K (e.g. Real or RealVect). */
  template <std::size_t K>
  using ColumnType = typename detail::MemberPointerTraits<std::tuple_element_t<K, ColumnsTuple>>::ValueType;

  /**
    @brief Derive a column index from a pointer-to-member at compile time.
    @details Used to implement the member-pointer accessors. Reorder-safe; returns
             s_numColumns if the member is not a registered column.
    @param[in] a_member Pointer to a particle data member.
    @return Index of the column backing that member.
  */
  template <typename M>
  static constexpr std::size_t
  columnIndex(const M a_member) noexcept
  {
    return detail::indexOf(Traits::columns, a_member);
  }

  /** @brief Column index of the mandatory position field (DERIVED, not hand-numbered). */
  static constexpr std::size_t s_positionIndex = detail::indexOf(Traits::columns, Traits::positionPtr);

  /** @brief Column index of the mandatory weight field (DERIVED, not hand-numbered). */
  static constexpr std::size_t s_weightIndex = detail::indexOf(Traits::columns, Traits::weightPtr);

  /**
    @brief Byte count of one fully-linearized particle (all columns, for MPI).
    @return Number of bytes one particle occupies in a linearizeParticle() buffer.
  */
  static constexpr std::size_t
  bytesPerParticle() noexcept
  {
    return detail::LinearSize<ColumnsTuple>::s_value;
  }

  /**
    @brief The subset of columns written to HDF5, as a std::index_sequence.
    @details Derived from the traits' optional `static constexpr auto h5Columns`
    member-pointer tuple (reorder-safe). Defaults to "all columns" if absent.
  */
  using H5Columns = typename detail::H5ColumnSubset<Traits, std::make_index_sequence<s_numColumns>>::Type;

  // ----- Sanity checks on the user's traits ---------------------------------
  static_assert(s_numColumns >= 1, "a particle must have at least one column");
  static_assert(std::is_default_constructible<P>::value,
                "particle type P must be default-constructible (gather/HDF5 restart build a P)");
  static_assert(detail::AllBlittable<ColumnsTuple>::s_value,
                "every particle column must be byte-copyable (trivially copyable or RealVect) for MPI/HDF5");
  static_assert(s_positionIndex < s_numColumns, "ParticleTraits::positionPtr is not a registered column");
  static_assert(s_weightIndex < s_numColumns, "ParticleTraits::weightPtr is not a registered column");
  static_assert(std::is_same<ColumnType<s_positionIndex>, RealVect>::value,
                "ParticleTraits::positionPtr must designate a RealVect member");
  static_assert(std::is_same<ColumnType<s_weightIndex>, Real>::value,
                "ParticleTraits::weightPtr must designate a Real member");

  /** @brief Construct an empty container. */
  ParticleSoA() = default;

  /** @brief Copy construction is deleted (move-only: avoids accidental O(N) deep copies). */
  ParticleSoA(const ParticleSoA&) = delete;

  /** @brief Copy assignment is deleted (move-only). */
  ParticleSoA&
  operator=(const ParticleSoA&) = delete;

  /** @brief Move constructor. */
  ParticleSoA(ParticleSoA&&) noexcept = default;

  /** @brief Move assignment operator. @return Reference to this. */
  ParticleSoA&
  operator=(ParticleSoA&&) noexcept = default;

  /** @brief Destructor. */
  ~ParticleSoA() = default;

  /** @brief Number of particles currently stored. */
  std::size_t
  size() const noexcept
  {
    return std::get<0>(m_columns).size();
  }

  /** @brief Whether the container holds no particles. */
  bool
  empty() const noexcept
  {
    return this->size() == 0;
  }

  /** @brief Reserve capacity in every column. */
  void
  reserve(const std::size_t a_capacity)
  {
    this->reserveImpl(a_capacity, std::make_index_sequence<s_numColumns>{});
  }

  /** @brief Resize every column to a_count particles (new entries value-initialized). */
  void
  resize(const std::size_t a_count)
  {
    this->resizeImpl(a_count, std::make_index_sequence<s_numColumns>{});
  }

  /** @brief Drop all particles (keeps capacity). */
  void
  clear() noexcept
  {
    this->clearImpl(std::make_index_sequence<s_numColumns>{});
  }

  /**
    @brief Append a single particle (AoS view) to the back of every column.
    @param[in] a_p Particle to append.
  */
  void
  append(const P& a_p)
  {
    this->appendImpl(a_p, std::make_index_sequence<s_numColumns>{});
  }

  /**
    @brief Gather column slice i back into the AoS particle view.
    @param[in] a_index Particle index.
    @return Reconstructed particle.
  */
  P
  gather(const std::size_t a_index) const
  {
    return this->gatherImpl(a_index, std::make_index_sequence<s_numColumns>{});
  }

  /**
    @brief Remove particle i using swap-and-pop (O(1), does NOT preserve order).
    @param[in] a_index Particle index.
  */
  void
  remove(const std::size_t a_index)
  {
    this->removeImpl(a_index, std::make_index_sequence<s_numColumns>{});
  }

  /**
    @brief Access an entire column by member pointer (the public field selector).
    @details Single-token C++17 form: soa.column<&P::velocity>(). The column index
    is derived from the member pointer, so it is reorder-safe and zero-overhead.
    @return Reference to the std::vector backing that column.
  */
  template <auto Member>
  auto&
  column() noexcept
  {
    return this->columnByIndex<columnIndex(Member)>();
  }

  /** @brief Const overload of column<&P::field>(). */
  template <auto Member>
  const auto&
  column() const noexcept
  {
    return this->columnByIndex<columnIndex(Member)>();
  }

  /**
    @brief Access a field of particle i by member pointer (generalizes position(i)/weight(i)).
    @details Single-token C++17 form, e.g. soa.get<&P::lightsaberDamage>(i) *= 0.9;
    @param[in] a_index Particle index.
    @return Reference to that particle's value for the selected field.
  */
  template <auto Member>
  auto&
  get(const std::size_t a_index) noexcept
  {
    return this->columnByIndex<columnIndex(Member)>()[a_index];
  }

  /** @brief Const overload of get<&P::field>(i). */
  template <auto Member>
  const auto&
  get(const std::size_t a_index) const noexcept
  {
    return this->columnByIndex<columnIndex(Member)>()[a_index];
  }

  /**
    @brief Access an entire column by compile-time index (low-level primitive).
    @details Prefer the member-pointer column<&P::field>() in user code. This is
    used internally and by generic algorithms that already hold an index.
    @return Reference to the std::vector backing column K.
  */
  template <std::size_t K>
  std::vector<ColumnType<K>>&
  columnByIndex() noexcept
  {
    return std::get<K>(m_columns);
  }

  /** @brief Const overload of columnByIndex<K>(). */
  template <std::size_t K>
  const std::vector<ColumnType<K>>&
  columnByIndex() const noexcept
  {
    return std::get<K>(m_columns);
  }

  /**
    @brief Overwrite particle i from an AoS view (inverse of gather).
    @details Completes the AoS round-trip gather -> modify -> scatter for code that
    prefers the per-particle struct. Prefer direct column access in hot loops.
    @param[in] a_index Particle index.
    @param[in] a_p     New values for every column of particle i.
  */
  void
  scatter(const std::size_t a_index, const P& a_p) noexcept
  {
    this->scatterImpl(a_index, a_p, std::make_index_sequence<s_numColumns>{});
  }

  /** @brief The mandatory position column (whole vector). */
  std::vector<RealVect>&
  positions() noexcept
  {
    return this->columnByIndex<s_positionIndex>();
  }

  /** @brief The mandatory position column (whole vector). Const overload. */
  const std::vector<RealVect>&
  positions() const noexcept
  {
    return this->columnByIndex<s_positionIndex>();
  }

  /** @brief The mandatory weight column (whole vector). */
  std::vector<Real>&
  weights() noexcept
  {
    return this->columnByIndex<s_weightIndex>();
  }

  /** @brief The mandatory weight column (whole vector). Const overload. */
  const std::vector<Real>&
  weights() const noexcept
  {
    return this->columnByIndex<s_weightIndex>();
  }

  /** @brief Position of particle i (the one mandatory vector accessor). */
  RealVect&
  position(const std::size_t a_index) noexcept
  {
    return this->columnByIndex<s_positionIndex>()[a_index];
  }

  /** @brief Const overload of position(i). */
  const RealVect&
  position(const std::size_t a_index) const noexcept
  {
    return this->columnByIndex<s_positionIndex>()[a_index];
  }

  /** @brief Weight of particle i (the one mandatory scalar accessor). */
  Real&
  weight(const std::size_t a_index) noexcept
  {
    return this->columnByIndex<s_weightIndex>()[a_index];
  }

  /** @brief Const overload of weight(i). */
  const Real&
  weight(const std::size_t a_index) const noexcept
  {
    return this->columnByIndex<s_weightIndex>()[a_index];
  }

  /**
    @brief Linearize particle i into a byte buffer (for MPI send).
    @param[in,out] a_buffer Buffer with at least bytesPerParticle() bytes.
    @param[in]     a_index  Particle index.
  */
  void
  linearizeParticle(void* a_buffer, const std::size_t a_index) const noexcept
  {
    auto* p = static_cast<unsigned char*>(a_buffer);
    this->linearizeImpl(p, a_index, std::make_index_sequence<s_numColumns>{});
  }

  /**
    @brief Append a particle delinearized from a byte buffer (for MPI receive).
    @param[in] a_buffer Buffer with at least bytesPerParticle() bytes.
  */
  void
  delinearizeAndAppend(const void* a_buffer)
  {
    const auto* p = static_cast<const unsigned char*>(a_buffer);
    this->delinearizeImpl(p, std::make_index_sequence<s_numColumns>{});
  }

  /** @brief Bytes for one HDF5-linearized particle (only the H5Columns subset). */
  static constexpr std::size_t
  h5BytesPerParticle() noexcept
  {
    return h5SizeImpl(H5Columns{});
  }

  /**
    @brief Linearize only the HDF5-checkpointed columns of particle i.
    @details Generic over any subset; mirrors GenericParticle::H5linearOut.
    @param[in,out] a_buffer Buffer with at least h5BytesPerParticle() bytes.
    @param[in]     a_index  Particle index.
  */
  void
  h5LinearizeParticle(void* a_buffer, const std::size_t a_index) const noexcept
  {
    auto* p = static_cast<unsigned char*>(a_buffer);
    this->h5LinearizeImpl(p, a_index, H5Columns{});
  }

  /**
    @brief Append a particle whose HDF5-checkpointed columns come from a buffer.
    @details Non-checkpointed columns are left default-constructed (they are
    recomputed after restart), mirroring GenericParticle::H5linearIn.
    @param[in] a_buffer Buffer with at least h5BytesPerParticle() bytes.
  */
  void
  h5DelinearizeAndAppend(const void* a_buffer)
  {
    const auto* p = static_cast<const unsigned char*>(a_buffer);
    this->h5DelinearizeImpl(p, H5Columns{});
  }

protected:
  /** @brief The SoA storage: one std::vector per column. */
  Storage m_columns;

private:
  template <std::size_t... I>
  void
  reserveImpl(const std::size_t a_cap, std::index_sequence<I...>)
  {
    using expander = int[];
    (void)expander{0, ((void)std::get<I>(m_columns).reserve(a_cap), 0)...};
  }

  template <std::size_t... I>
  void
  resizeImpl(const std::size_t a_count, std::index_sequence<I...>)
  {
    using expander = int[];
    (void)expander{0, ((void)std::get<I>(m_columns).resize(a_count), 0)...};
  }

  template <std::size_t... I>
  void
  clearImpl(std::index_sequence<I...>) noexcept
  {
    using expander = int[];
    (void)expander{0, ((void)std::get<I>(m_columns).clear(), 0)...};
  }

  template <std::size_t... I>
  void
  appendImpl(const P& a_p, std::index_sequence<I...>)
  {
    // Local constexpr copy avoids ODR-using the static member (C++14 would
    // otherwise demand an out-of-line definition of ParticleTraits::columns).
    constexpr auto cols = Traits::columns;
    using expander      = int[];
    (void)expander{0, ((void)std::get<I>(m_columns).push_back(a_p.*std::get<I>(cols)), 0)...};
  }

  template <std::size_t... I>
  P
  gatherImpl(const std::size_t a_index, std::index_sequence<I...>) const
  {
    constexpr auto cols = Traits::columns;
    P              p;
    using expander = int[];
    (void)expander{0, ((void)(p.*std::get<I>(cols) = std::get<I>(m_columns)[a_index]), 0)...};
    return p;
  }

  template <std::size_t... I>
  void
  scatterImpl(const std::size_t a_index, const P& a_p, std::index_sequence<I...>) noexcept
  {
    constexpr auto cols = Traits::columns;
    using expander      = int[];
    (void)expander{0, ((void)(std::get<I>(m_columns)[a_index] = a_p.*std::get<I>(cols)), 0)...};
  }

  template <std::size_t... I>
  void
  removeImpl(const std::size_t a_index, std::index_sequence<I...>)
  {
    using expander = int[];
    (void)expander{0, ((void)swapPop(std::get<I>(m_columns), a_index), 0)...};
  }

  template <std::size_t... I>
  void
  linearizeImpl(unsigned char* a_buf, const std::size_t a_index, std::index_sequence<I...>) const noexcept
  {
    using expander = int[];
    (void)expander{0, ((void)detail::pushBytes(a_buf, std::get<I>(m_columns)[a_index]), 0)...};
  }

  template <std::size_t... I>
  void
  delinearizeImpl(const unsigned char* a_buf, std::index_sequence<I...>)
  {
    using expander = int[];
    (void)expander{0, ((void)appendOneFromBytes<I>(a_buf), 0)...};
  }

  template <std::size_t K>
  void
  appendOneFromBytes(const unsigned char*& a_buf)
  {
    ColumnType<K> value;
    detail::pullBytes(a_buf, value);
    std::get<K>(m_columns).push_back(value);
  }

  template <std::size_t... Is>
  static constexpr std::size_t
  h5SizeImpl(std::index_sequence<Is...>) noexcept
  {
    return detail::sumSizes(sizeof(ColumnType<Is>)...);
  }

  template <std::size_t... Is>
  void
  h5LinearizeImpl(unsigned char* a_buf, const std::size_t a_index, std::index_sequence<Is...>) const noexcept
  {
    using expander = int[];
    (void)expander{0, ((void)detail::pushBytes(a_buf, std::get<Is>(m_columns)[a_index]), 0)...};
  }

  template <std::size_t... Is>
  void
  h5DelinearizeImpl(const unsigned char* a_buf, std::index_sequence<Is...>)
  {
    constexpr auto cols = Traits::columns;
    P              p{}; // non-checkpointed columns stay default-constructed
    using expander = int[];
    (void)expander{0, ((void)readColumnFromBytes(a_buf, p, std::get<Is>(cols)), 0)...};
    this->append(p);
  }

  template <typename MemberPtr>
  static void
  readColumnFromBytes(const unsigned char*& a_buf, P& a_p, MemberPtr a_member)
  {
    typename detail::MemberPointerTraits<MemberPtr>::ValueType value;
    detail::pullBytes(a_buf, value);
    a_p.*a_member = value;
  }

  template <typename T>
  static void
  swapPop(std::vector<T>& a_vec, const std::size_t a_index) noexcept
  {
    a_vec[a_index] = a_vec.back();
    a_vec.pop_back();
  }
};

#include <CD_NamespaceFooter.H>

#endif
