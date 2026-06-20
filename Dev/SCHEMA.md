# Column schema & field selector — detail for the design decision

Two linked questions:
- **Schema**: how the user declares which fields become SoA columns.
- **Selector**: how a call site names a field, e.g. `deposit<weight>(...)`.

## TL;DR

- **Runtime performance is identical across all schema options.** The choice is
  ergonomics / compile-time / tooling / error messages — not speed. (Reasoning below.)
- The **field selector can always deduce the field type** (no `Ret` template arg
  like today). Confirmed by the prototype: `deposit<MyParticle::Weight>(mesh, soa)`
  compiles and runs under `-std=c++14`, with `Real` vs `RealVect` deduced internally.
- The clean single-token selector form depends on the C++ standard:
  - **C++14**: `deposit<P::Weight>` via a named index/enum — clean. `deposit<&P::weight>`
    via a member pointer is **not** available as one arg (see below).
  - **C++17**: `deposit<&P::weight>` via member pointer (`template <auto>`) — cleanest,
    and lets the same `&P::weight` serve as both schema entry and selector.

So the real sub-decision is **stay C++14 (enum/index selector) vs. bump CXXSTD to 17
(member-pointer selector)**.

---

## Why runtime performance is identical

In all schemes the user struct `P` (e.g. `ItoParticle`) is only a **transient
Array-of-Structs view**, used by `gather()`/`append()` for single-particle work
(split/merge, MPI pack). The actual storage is always `std::tuple<std::vector<col>...>`.

Whatever names the columns — a member-pointer tuple, a macro that expands to that
tuple, or a CRTP base — the *generated storage and the field access are the same*:
`column<K>()[n]` is a direct array index, member pointers resolve to fixed offsets at
compile time, nothing is virtual. There is no runtime dispatch in any option.

The one performance rule that matters is **orthogonal to the schema**: never use the
AoS `gather()` view inside hot loops — iterate `column<K>()` directly (contiguous,
vectorizable). That holds for every schema.

---

## The three schema options, concretely (ItoParticle, 7 fields)

### (a) Traits + member-pointer tuple  — *current prototype*
```cpp
struct ItoParticle {
  RealVect pos, velocity, oldPos;
  Real weight, mobility, diffusion, energy;
  enum Field : std::size_t { Pos, Velocity, OldPos, Weight, Mobility, Diffusion, Energy };
};
template <> struct ParticleTraits<ItoParticle> {
  static constexpr auto columns = std::make_tuple(
    &ItoParticle::pos, &ItoParticle::velocity, &ItoParticle::oldPos,
    &ItoParticle::weight, &ItoParticle::mobility, &ItoParticle::diffusion, &ItoParticle::energy);
  static constexpr std::size_t positionIndex = ItoParticle::Pos;
  static constexpr std::size_t weightIndex   = ItoParticle::Weight;
};
```
- Explicit; no hidden machinery; field names are real struct members.
- Boilerplate: the member list is written twice (struct + tuple). The `enum` gives
  named selectors.
- Plays nicely with clang-format / doxygen / debuggers.

### (b) Macro-generated traits
```cpp
struct ItoParticle { RealVect pos, velocity, oldPos; Real weight, mobility, diffusion, energy; };
CD_PARTICLE_LAYOUT(ItoParticle, /*pos*/ pos, /*weight*/ weight,
                   velocity, oldPos, mobility, diffusion, energy);
```
- Expands to **exactly** option (a): identical generated code, identical performance.
- Much less typing; single source of truth for the field list.
- Downsides: macros are opaque to debuggers and to the repo's `doxygen-check` and
  `clang-format` pre-commit hooks; harder to read in error messages.

### (c) CRTP base class
```cpp
struct ItoParticle : ParticleBase<ItoParticle,
                                  Col<RealVect,Pos>, Col<Real,Weight>, Col<RealVect,Velocity>, ...> {};
```
- Can auto-generate named accessors (`p.weight()`) so the AoS view matches today's API.
- Couples every particle type to the framework base; template error messages get
  deep; storage/performance still identical.

---

## Field selector — what each standard allows

The selector reduces to "pick a column at compile time; deduce its type". The type
(and component count) is always deducible from the column, so the caller never writes
the return type. The difference is only the *spelling*:

| Form | Syntax | Standard | Notes |
|------|--------|----------|-------|
| Named index / enum | `deposit<P::Weight>` | **C++14** | Clean; needs an enum on P. Works today (tested). |
| Tag type | `deposit<Weight>` | C++14 | Clean; needs a tag struct per field. |
| Member pointer | `deposit<&P::weight>` | **C++17** | Cleanest; `template <auto>`. Unifies schema+selector. |
| Member pointer | `deposit<decltype(&P::weight), &P::weight>` or `deposit<FIELD(weight)>` | C++14 | Verbose / macro. Confirmed: single-arg member-pointer NTTP is impossible in C++14. |

Verified with the compiler:
- `deposit<MyParticle::Weight>(mesh, soa)` — **builds & runs, `-std=c++14`** (see
  `test_ParticleSoA.cpp::testFieldSelector`); the field type is deduced.
- single-arg `deposit<&P::weight>` — **C++14 reject, C++17 accept** (`template <auto>`).

If we use member pointers, we can compile-time **match** a `&P::weight` against the
traits tuple to find its column index — so the same member pointers can serve as both
the schema declaration and the selector. But the clean one-token call still needs C++17.

---

## Robustness to layout changes (e.g. swapping two columns)

What happens downstream when you reorder columns depends entirely on **how a field
is identified**, not on whether a macro was used:

1. **By name, with indices DERIVED from member pointers** *(recommended)* — reorder
   is safe, **zero downstream changes**, stays correct. The column index is computed
   at compile time by matching `&P::weight` against the column tuple, so the name
   follows the field wherever it moves. Confirmed working under C++14:

   ```cpp
   // type-guarded member-pointer equality (false if member types differ)
   template <typename A, typename B> constexpr bool eqMember(A, B)    { return false; }
   template <typename A>             constexpr bool eqMember(A a, A b){ return a == b; }
   // ... constexpr search of the tuple ...
   static constexpr std::size_t weightIndex = indexOf(columns, &P::weight); // DERIVED
   ```

2. **By name, with HAND-NUMBERED indices** (e.g. a separate `enum { Weight = 1 }`
   that must agree with the tuple by convention) — reorder is a **silent footgun**:
   if you swap the tuple entries but forget the enum, the wrong column is
   deposited/linearized with no compile error. This is the hazard to design out.

3. **By literal integer index** (`column<1>()`, or assuming `positionIndex == 0`
   somewhere) — **always breaks silently** on reorder. Forbid by convention.

### Where a macro helps, and what it cannot fix

- A `CD_PARTICLE_LAYOUT(...)` macro makes the column list a **single source of
  truth** and can emit *derived* indices (case 1) for you — so it is strictly
  **safer** than hand-written traits for reordering, because there is no duplicated
  list to drift out of sync. With the macro, swapping two columns is a one-line edit
  and all name-based call sites are unaffected.

- What no schema choice can hide: **column order is part of the serialized ABI.**
  Reordering (or adding/removing) columns changes the byte layout emitted by
  `linearizeParticle`, so **existing HDF5 checkpoints and in-flight MPI buffers
  become incompatible** unless the format is versioned or keyed by field name. Treat
  a layout change as a checkpoint-format change.

## Status in the prototype

Both halves of this discussion are now implemented and tested (`make check`):

- **Derived indices** — `ParticleTraits` no longer carries hand-written
  `positionIndex`/`weightIndex`. The user designates position/weight by member
  pointer (`positionPtr`/`weightPtr`); `ParticleSoA` derives the indices via the
  constexpr `indexOf`. `ParticleSoA::columnIndex(&P::field)` gives a reorder-safe
  selector for any field. Proven by `test_ParticleSoA.cpp::testReorderSafety`
  (`weight` is column 1 in one particle, column 0 in another; both resolve right).
- **`CD_PARTICLE_LAYOUT` macro** — `CD_ParticleLayoutMacro.H`. Generates the exact
  same traits from one member list. `test_ParticleLayoutMacro.cpp` defines the same
  7-field particle both ways and `static_assert`s identical column count, derived
  indices, byte sizes, and column types.

## Recommendation

- **Performance**: not a factor — choose on ergonomics.
- **If staying C++14**: schema (a) traits-tuple **+ an `enum Field` on the particle**,
  selector `deposit<P::Weight>`. Zero new syntax beyond the enum; type deduced.
  Optionally add macro (b) later as pure sugar over (a).
- **If bumping `CXXSTD` to 17**: use **member pointers** for both schema and selector
  → `deposit<&P::weight>`, the traits list the same `&P::weight`, one unified
  mechanism, least boilerplate. GCC and Intel oneAPI (your two CI compilers) fully
  support C++17.

Open decision to make: **C++14 + enum selector, or bump to C++17 + member-pointer
selector?**
