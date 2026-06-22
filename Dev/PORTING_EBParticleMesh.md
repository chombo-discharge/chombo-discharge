# Porting EBParticleMesh / EBAMRParticleMesh to ParticleSoA

Design + plan for moving the particle-mesh deposition/interpolation classes onto the new
`ParticleSoA` leaf. Source under review: `Source/Particle/CD_EBParticleMesh.{H,Impl}` and
`Source/Particle/CD_EBAMRParticleMesh.{H,Impl}`.

## Scope & working mode

- **All work happens in `Dev/` until design-freeze.** We stage SoA copies of these classes
  in `Dev/`, exercise them in the Dev test executable, and only replace the production
  `Source/Particle` versions at freeze.
- **Coexistence naming.** The Dev test links `libSource`, which already contains the
  production `EBParticleMesh`/`EBAMRParticleMesh` (same namespace). To avoid ODR clashes,
  the staged classes get a distinct name/file while in `Dev/` — e.g.
  `CD_EBParticleMeshSoA.H` / `class EBParticleMeshSoA` — renamed back to the canonical names
  when they replace production at freeze.

## Sequencing (decided)

1. **EBParticleMesh → Dev, tested.** Port the single-patch leaf class first; verify
   **correctness and performance** in `Dev/` against the existing kernels. It only needs the
   `ParticleSoA` leaf, which already exists.
2. **EBAMRParticleMesh → Dev.** Port the AMR wrapper next.
3. **ParticleContainer → Dev.** Finalize the container last.

> **Dependency note (sequencing tension).** EBAMRParticleMesh consumes a *container*: it pulls
> per-patch particle storage and the halo/mask/buffer holders from `ParticleContainer<P>`
> (`a_particles[lvl][din]`, `getParticles()`, `getMaskParticles()`). So step 2 cannot be fully
> exercised without at least a **minimal Dev `ParticleContainer` scaffold** (per-patch
> `LayoutData<ParticleSoA<P>>` accessor + the mask/halo holders). Practical reading of the
> plan: step 2 stands up that scaffold against the intended container interface, and step 3 is
> the *hardening/full* port (remap, regrid count→reserve→fill, MPI scatter — see
> `ARCHITECTURE.md`). Keep the container interface that EBAMRParticleMesh codes against stable
> between steps 2 and 3.

## What ports UNCHANGED (storage-agnostic)

These are the bulk of the code and need no SoA changes:

- **EBParticleMesh per-particle kernels** — `depositParticle{NGP,CIC,TSC}` and
  `interpolateParticle{NGP,CIC,TSC}`. Their boundary is already SoA-friendly:
  `const RealVect& position`, a raw `Real* strength` / `Real* particleField`, and `numComp`.
  ParticleSoA feeds exactly this. The EB/cut-cell weight math, `EBISBox` covered/irregular
  handling, and the tiny stencil `BoxLoops` stay verbatim.
- **All EBAMRParticleMesh coarse-fine machinery** — `EBCoarseFineParticleMesh`, the
  outer-halo and transition masks, `m_levelCopiers`, `EBAMRData` arithmetic, the four CF
  strategies (Interp / Halo / HaloNGP / Transition), `getTransitionMaskWidth`. None of it
  reads fields out of particles; it operates on mesh data and on *which* particles deposit.

## Incompatibilities (all in the storage / selector boundary)

1. **`List<P>` storage in every signature.** `EBParticleMesh::deposit/interpolate` take
   `List<P>&` and iterate with `ListIterator`; `EBAMRParticleMesh` pulls per-patch lists via
   `a_particles[lvl][din].listItems()` and `getParticles()/getMaskParticles()` returning
   `AMRParticles<P> = Vector<RefCountedPtr<LayoutData<ListBox<P>>>>`. → per-patch
   `ParticleSoA<P>&` leaf + the new `ParticleContainer` accessors.

2. **Field selector is a pointer-to-member-FUNCTION** — `Ret (P::*MemberFunc)() const`
   (deposit) / `Ret (P::*MemberFunc)()` (interpolate), called as `(p.*MemberFunc)()`, e.g.
   `deposit<P, const Real&, &P::weight>`. ParticleSoA payload has **no member functions** —
   fields are data columns selected by a `template <auto>` *data*-member pointer, and the
   mandatory position/weight are reached via dedicated accessors (`weightColumn()`,
   `positionColumn(d)`), which are **not** payload member pointers at all. The entire
   `<P, Ret, MemberFunc>` signature is redesigned — see "Selector design" below.

3. **Vector fields are per-component, not a contiguous `RealVect`.** The current code relies
   on contiguity: deposit does `const Ret w = (...)(); &w`; interpolate's `GetPointer` hands
   the worker `RealVect::dataPtr()` and the worker writes `a_particleField[comp]`. In
   ParticleSoA a vector payload is `vx[]`,`vy[]`,`vz[]` — **separate, non-contiguous columns**
   — so no single `Real*` spans the components and no `RealVect&` lvalue into storage exists.
   Fixed by marshalling in the loop wrapper (see "Marshalling"); the workers are unchanged.

4. **Precision (`ParticleReal` vs `Real`).** Position is now `Real` (double), so `position(i)`
   matches the workers' `const RealVect&` exactly — no change. Weight is now `Real`, so the
   charge-deposit path is clean. But payload columns may be `float` while the workers compute
   in `double` → the same marshalling step promotes on deposit and demotes on interpolate.

5. **`sanitize<Ret>` / numComp deduction** disappear — `Ret`/`numComp` were deduced from the
   accessor return type; now `numComp` is the count of selected columns and the column value
   type comes from the member pointer (`ParticleSoA::ColumnType` / `decltype`).

## Selector design — DECIDED: option (b), variadic member-pointer selector

A "vector field" is now `SpaceDim` member pointers, not one `RealVect`. We use a **variadic
pack of payload data-member pointers** (1 pointer = scalar, `SpaceDim` = vector), plus a
**dedicated path for the mandatory weight** column.

```cpp
// --- deposit ---------------------------------------------------------------
// Charge density (the dominant case): deposit the mandatory weight column.
ebpm.depositWeight(rho, soa, DepositionType::CIC, widthScale, forceIrregNGP);

// Deposit a payload scalar column (e.g. energy density strength = energy):
ebpm.deposit<&Payload::energy>(rho, soa, DepositionType::CIC, ...);

// Deposit a payload vector (pack of SpaceDim member pointers, dimension-guarded):
ebpm.deposit<D_DECL(&Payload::vx, &Payload::vy, &Payload::vz)>(rho3, soa, ...);

// --- interpolate (target is always a payload column) -----------------------
// Interpolate a scalar mesh field into a payload scalar:
ebpm.interpolate<&Payload::phi>(soa, meshScalar, DepositionType::CIC, ...);

// Interpolate a vector mesh field into a payload vector (the hot velocity path):
ebpm.interpolate<D_DECL(&Payload::vx, &Payload::vy, &Payload::vz)>(soa, meshVec, ...);
```

Rationale: preserves **single-pass** vector interpolation (the hot E→velocity path), keeps
call sites readable, and `numComp = sizeof...(Members)` is deduced. Interpolate targets are
always payload columns (you never interpolate into weight); deposit needs the extra
`depositWeight` overload because weight is container-owned, not a payload member pointer.
`static_assert(sizeof...(Members) == 1 || sizeof...(Members) == SpaceDim)` mirrors the old
"Real or RealVect" check.

## Marshalling (the loop wrapper)

The per-particle workers keep their `Real* + numComp` contract; the new loop wrapper bridges
to the (possibly `float`, possibly non-contiguous) columns via a stack scratch buffer:

```cpp
// deposit, per particle i (Members... = selected columns; here a pack):
Real strength[numComp];
// gather + promote from the per-component columns:
{ int k = 0; ( (strength[k++] = static_cast<Real>(soa.column<Members>()[i])), ... ); }
this->depositParticleCIC(rho, soa.position(i), ..., strength, numComp, forceIrregNGP);

// interpolate, per particle i:
Real field[numComp] = {0};
this->interpolateParticleCIC(field, mesh, ..., soa.position(i), numComp, forceIrregNGP);
// scatter + demote back into the per-component columns:
{ int k = 0; ( (soa.column<Members>()[i] = static_cast<ColType<Members>>(field[k++])), ... ); }
```

For `depositWeight` the scalar fast path is trivial (`Real s = soa.weightColumn()[i];` →
`&s`, `numComp = 1`) — weight is already `Real`, no promotion. Position comes from
`soa.position(i)` (already `Real`/`RealVect`, exact). Capture column base pointers ONCE
outside the loop (as in `CD_ParticleLoops.H`), not per particle.

## Testing in Dev (step 1 exit criteria)

Before moving to EBAMRParticleMesh, the Dev EBParticleMesh must demonstrate:
- **Correctness:** NGP/CIC/TSC deposit is mass-conserving and bitwise-comparable to the
  production `List<P>` path for the same particle set + order (regular and cut cells,
  `forceIrregNGP` on/off); interpolate round-trips a known mesh field within tolerance.
- **Performance:** scalar (weight) deposit/interp at parity with the benchmark numbers
  (`BENCHMARKS.md`: deposition/interpolation are grid-bound and layout-insensitive, ~tie);
  vector interpolate single-pass via the pack, no per-particle `RealVect` rebuild in the hot
  loop.

## Known bug carried by the port: production TSC deposit is not partition-of-unity

While testing the Dev port we found a pre-existing bug in production
`Source/Particle/CD_EBParticleMeshImplem.H :: depositParticleTSC`. The per-dimension cell
overlap integral is

```cpp
const Real factor = (alpha < beta) ? 1.0 : 0.0;   // out-of-support guard
weight *= factor * (beta - alpha) - (beta*|beta| - alpha*|alpha|) / L;   // BUG: factor on 1st term only
```

`factor` is meant to zero a cell that lies outside the cloud support (`alpha >= beta`), but it
multiplies only the first term — the second term `-(beta|beta| - alpha|alpha|)/L` stays active
for out-of-support cells and injects spurious mass. Result: TSC deposit over-deposits by
**~2.33x per dimension** (5.4x in 2D, 12.7x in 3D) and is not charge-conserving. The correct
form gates the whole integral:

```cpp
weight *= factor * ((beta - alpha) - (beta*|beta| - alpha*|alpha|) / L);
```

Verified: with the fix, a particle 0.3 into a cell deposits `0.245, 0.71, 0.045` (sum = 1) — the
quadratic B-spline weights the **interpolate** path already uses, so deposit and interpolate
become consistent (adjoint).

**Decision:** `EBParticleMeshSoA` ships the FIXED kernel. **Production
`CD_EBParticleMeshImplem.H` is intentionally left unpatched** — TSC deposit is unused in
production, so we avoid touching `Source/` before design-freeze; patch it (or let the SoA port
supersede it) at freeze. The Dev correctness test reflects the intentional divergence: bitwise
SoA==production for NGP/CIC deposit and for all interpolation, and partition-of-unity for TSC
deposit (where the paths now differ).

## Open items / decisions still to make

- **`depositWeight` vs a unified selector.** We use a dedicated weight overload; if a
  weight-selector *token* (so all deposits share one entry point) reads cleaner once the
  call sites are ported, revisit. Low stakes.
- **Cut-cell precision.** Cell-index + CIC/TSC weights compute in `Real` (double) from the
  double position — confirm this is what the mixed-precision plan wants at cut cells
  (`ARCHITECTURE.md` precision section); demotion happens only at the column write.
- **Container interface for step 2** — pin down the minimal per-patch accessor + mask/halo
  holder API that EBAMRParticleMesh codes against, so step 3 doesn't churn it.
