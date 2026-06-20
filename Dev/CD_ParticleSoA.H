/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ParticleSoA.H
  @brief  Prototype of a Struct-of-Arrays (SoA) particle storage container.
  @author Robert Marskar
  @details This is a DESIGN PROTOTYPE for replacing Chombo's List<P> with an SoA
           layout in chombo-discharge. It is intentionally self-contained and
           C++17 (the project is bumping CXXSTD to 17 for the single-token
           member-pointer field selector, e.g. soa.column<&P::velocity>()).

           The central idea:

           1. The user defines a plain "Array-of-Structs" particle struct, e.g.
                struct MyParticle { RealVect pos; Real weight; Real lightsaberDamage; };
              This is the ergonomic single-particle view the rest of the code uses.

           2. The user specializes ParticleTraits<MyParticle> to describe the SoA
              columns as a compile-time tuple of pointer-to-member, and to tag which
              column is the position and which is the weight.

           3. ParticleSoA<MyParticle> then stores one contiguous std::vector per
              column (true SoA), but can still hand back / accept a MyParticle via
              gather()/append() for code that prefers the per-particle view.

           The two hard requirements requested are satisfied as:
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

// Prototype proxy types (swap for Chombo <REAL.H>/<RealVect.H> in the real impl).
#include "CD_RealProxy.H"

namespace proxy {

  /**
    @brief Traits class that a user specializes to teach ParticleSoA how to
           decompose their particle struct into SoA columns.
    @details A valid specialization must provide:

      static constexpr auto columns     = std::make_tuple(&P::a, &P::b, ...);
      static constexpr auto positionPtr = &P::<the RealVect position member>;
      static constexpr auto weightPtr   = &P::<the Real weight member>;

    Note that the user designates position/weight by MEMBER POINTER, not by a
    hand-written index. ParticleSoA derives the column indices from these pointers,
    so reordering the `columns` tuple cannot silently misroute access.

    Optionally, a specialization may also provide an HDF5 subset, again by member
    pointer (indices derived, so it follows column reorders):

      static constexpr auto h5Columns = std::make_tuple(&P::pos, &P::weight, ...);

    The primary template is left undefined on purpose so that using ParticleSoA
    with an undescribed particle type is a hard compile error.
  */
  template <typename P>
  struct ParticleTraits;

  namespace detail {

    /**
      @brief Extract the class type and value type from a pointer-to-member type.
    */
    template <typename M>
    struct member_pointer_traits;

    template <typename C, typename T>
    struct member_pointer_traits<T C::*>
    {
      using class_type = C;
      using value_type = T;
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
    eqMember(A a, A b) noexcept
    {
      return a == b;
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
    struct soa_storage;

    template <typename... Ms>
    struct soa_storage<std::tuple<Ms...>>
    {
      using type = std::tuple<std::vector<typename member_pointer_traits<Ms>::value_type>...>;
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
    struct linear_size;

    template <typename... Ms>
    struct linear_size<std::tuple<Ms...>>
    {
      static constexpr std::size_t value = sumSizes(sizeof(typename member_pointer_traits<Ms>::value_type)...);
    };

    /** @brief C++14 stand-in for std::void_t. */
    template <typename...>
    struct make_void
    {
      using type = void;
    };
    template <typename... Ts>
    using void_t = typename make_void<Ts...>::type;

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
    struct h5_columns
    {
      using type = Fallback;
    };
    template <typename Traits, typename Fallback>
    struct h5_columns<Traits, Fallback, void_t<decltype(Traits::h5Columns)>>
    {
      using H5Tuple = typename std::remove_cv<typename std::remove_reference<decltype(Traits::h5Columns)>::type>::type;
      using type    = decltype(deriveH5Seq<Traits>(std::make_index_sequence<std::tuple_size<H5Tuple>::value>{}));
    };

    /**
      @brief memcpy one trivially-copyable value onto a byte buffer and advance it.
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
      @brief memcpy one trivially-copyable value off a byte buffer and advance it.
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
    @brief Struct-of-Arrays particle container for a single grid patch / bin.
    @details One ParticleSoA holds all the particles for one logical region
    (analogue of today's List<P> stored inside a BinFab). The AMR/patch machinery
    (ParticleContainer) would own a collection of these, exactly as it owns
    List<P> today.

    @tparam P      User particle struct (the AoS single-particle view).
    @tparam Traits Descriptor; defaults to ParticleTraits<P>.
  */
  template <typename P, typename Traits = ParticleTraits<P>>
  class ParticleSoA
  {
  public:
    /** @brief The columns descriptor tuple type (tuple of pointer-to-member). */
    using ColumnsTuple = std::remove_cv_t<std::remove_reference_t<decltype(Traits::columns)>>;

    /** @brief The SoA storage type (tuple of std::vector per column). */
    using Storage = typename detail::soa_storage<ColumnsTuple>::type;

    /** @brief Number of SoA columns. */
    static constexpr std::size_t NumColumns = std::tuple_size<ColumnsTuple>::value;

    /** @brief Value type of column K (e.g. Real or RealVect). */
    template <std::size_t K>
    using ColumnType = typename detail::member_pointer_traits<std::tuple_element_t<K, ColumnsTuple>>::value_type;

    /**
      @brief Derive a column index from a pointer-to-member at compile time.
      @details Used to implement the member-pointer accessors. Reorder-safe; returns
               NumColumns if the member is not a registered column.
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
    static constexpr std::size_t positionIndex = detail::indexOf(Traits::columns, Traits::positionPtr);

    /** @brief Column index of the mandatory weight field (DERIVED, not hand-numbered). */
    static constexpr std::size_t weightIndex = detail::indexOf(Traits::columns, Traits::weightPtr);

    /** @brief Byte count of one fully-linearized particle (for MPI). */
    static constexpr std::size_t s_bytesPerParticle = detail::linear_size<ColumnsTuple>::value;

    /**
      @brief The subset of columns written to HDF5, as a std::index_sequence.
      @details Derived from the traits' optional `static constexpr auto h5Columns`
      member-pointer tuple (reorder-safe). Defaults to "all columns" if absent.
    */
    using H5Columns = typename detail::h5_columns<Traits, std::make_index_sequence<NumColumns>>::type;

    // ----- Sanity checks on the user's traits ---------------------------------
    static_assert(positionIndex < NumColumns, "ParticleTraits::positionPtr is not a registered column");
    static_assert(weightIndex < NumColumns, "ParticleTraits::weightPtr is not a registered column");
    static_assert(std::is_same<ColumnType<positionIndex>, RealVect>::value,
                  "ParticleTraits::positionPtr must designate a RealVect member");
    static_assert(std::is_same<ColumnType<weightIndex>, Real>::value,
                  "ParticleTraits::weightPtr must designate a Real member");

    /** @brief Construct an empty container. */
    ParticleSoA() = default;

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
      this->reserveImpl(a_capacity, std::make_index_sequence<NumColumns>{});
    }

    /** @brief Drop all particles (keeps capacity). */
    void
    clear() noexcept
    {
      this->clearImpl(std::make_index_sequence<NumColumns>{});
    }

    /**
      @brief Append a single particle (AoS view) to the back of every column.
      @param[in] a_p Particle to append.
    */
    void
    append(const P& a_p)
    {
      this->appendImpl(a_p, std::make_index_sequence<NumColumns>{});
    }

    /**
      @brief Gather column slice i back into the AoS particle view.
      @param[in] a_index Particle index.
      @return Reconstructed particle.
    */
    P
    gather(const std::size_t a_index) const
    {
      return this->gatherImpl(a_index, std::make_index_sequence<NumColumns>{});
    }

    /**
      @brief Remove particle i using swap-and-pop (O(1), does NOT preserve order).
      @param[in] a_index Particle index.
    */
    void
    remove(const std::size_t a_index)
    {
      this->removeImpl(a_index, std::make_index_sequence<NumColumns>{});
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
      this->scatterImpl(a_index, a_p, std::make_index_sequence<NumColumns>{});
    }

    /** @brief The mandatory position column (whole vector). */
    std::vector<RealVect>&
    positions() noexcept
    {
      return this->columnByIndex<positionIndex>();
    }

    /** @brief The mandatory position column (whole vector). Const overload. */
    const std::vector<RealVect>&
    positions() const noexcept
    {
      return this->columnByIndex<positionIndex>();
    }

    /** @brief The mandatory weight column (whole vector). */
    std::vector<Real>&
    weights() noexcept
    {
      return this->columnByIndex<weightIndex>();
    }

    /** @brief The mandatory weight column (whole vector). Const overload. */
    const std::vector<Real>&
    weights() const noexcept
    {
      return this->columnByIndex<weightIndex>();
    }

    /** @brief Position of particle i (the one mandatory vector accessor). */
    RealVect&
    position(const std::size_t a_index) noexcept
    {
      return this->columnByIndex<positionIndex>()[a_index];
    }

    /** @brief Const overload of position(i). */
    const RealVect&
    position(const std::size_t a_index) const noexcept
    {
      return this->columnByIndex<positionIndex>()[a_index];
    }

    /** @brief Weight of particle i (the one mandatory scalar accessor). */
    Real&
    weight(const std::size_t a_index) noexcept
    {
      return this->columnByIndex<weightIndex>()[a_index];
    }

    /** @brief Const overload of weight(i). */
    const Real&
    weight(const std::size_t a_index) const noexcept
    {
      return this->columnByIndex<weightIndex>()[a_index];
    }

    /**
      @brief Linearize particle i into a byte buffer (for MPI send).
      @param[in,out] a_buffer Buffer with at least s_bytesPerParticle bytes.
      @param[in]     a_index  Particle index.
    */
    void
    linearizeParticle(void* a_buffer, const std::size_t a_index) const noexcept
    {
      auto* p = static_cast<unsigned char*>(a_buffer);
      this->linearizeImpl(p, a_index, std::make_index_sequence<NumColumns>{});
    }

    /**
      @brief Append a particle delinearized from a byte buffer (for MPI receive).
      @param[in] a_buffer Buffer with at least s_bytesPerParticle bytes.
    */
    void
    delinearizeAndAppend(const void* a_buffer)
    {
      const auto* p = static_cast<const unsigned char*>(a_buffer);
      this->delinearizeImpl(p, std::make_index_sequence<NumColumns>{});
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
      typename detail::member_pointer_traits<MemberPtr>::value_type value;
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

  /**
    @brief Demonstrator: deposit a chosen particle field onto a 1-D mesh.
    @details Shows the C++17 call-site ergonomics: the field is selected by a single
    member-pointer token, and its TYPE (Real vs RealVect, and hence the component
    count) is DEDUCED - the caller never spells the return type. Compare to today's
    deposit<P, const Real&, &P::weight>.

    Usage:
        deposit<&MyParticle::weight>(mesh, soa);

    @tparam Member Pointer-to-member selecting the field to deposit.
    @param[in,out] a_mesh Fake 1-D mesh, indexed by floor(position[0]).
    @param[in]     a_soa  Particle container.
  */
  template <auto Member, typename P, typename Traits>
  void
  deposit(std::vector<Real>& a_mesh, const ParticleSoA<P, Traits>& a_soa)
  {
    constexpr std::size_t column = ParticleSoA<P, Traits>::columnIndex(Member);
    using FieldType              = typename ParticleSoA<P, Traits>::template ColumnType<column>;

    // Number of components is deduced, not supplied by the caller.
    constexpr int numComp = std::is_same<FieldType, Real>::value ? 1 : SpaceDim;
    static_assert(std::is_same<FieldType, Real>::value || std::is_same<FieldType, RealVect>::value,
                  "deposit field must be Real or RealVect");
    static_assert(numComp == 1, "this 1-D demonstrator only deposits scalar fields");

    const std::vector<RealVect>&  pos   = a_soa.positions();
    const std::vector<FieldType>& field = a_soa.template column<Member>();

    for (std::size_t n = 0; n < a_soa.size(); n++) {
      a_mesh[static_cast<std::size_t>(pos[n][0])] += field[n];
    }
  }

} // namespace proxy

#endif
