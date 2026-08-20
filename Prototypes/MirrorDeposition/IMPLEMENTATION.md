# Mirrored cut-cell deposition — implementation plan

**Design** `PLAN.md` revision 6 · **Issue** chombo-discharge#29 (parts 1 and 3) · **Branch**
`mirror_deposition` · **Depends on** PR #700

`PLAN.md` is the design and stays the authority on *why*. This file is the work order: what changes,
in what order, in which file, and what proves each step. Where it deviates from `PLAN.md` it says so
explicitly in §0.3 — those are the only deviations, and each is a correction rather than a
re-litigation.

All code citations below were re-resolved against the working tree on 2026-08-19 and hold.

---

## 0. Decisions taken before writing code

### 0.1 Staging

`PLAN.md` §6 splits the work into PR A (mechanism plus one consumer) and PR B (central wiring). PR A
as written cannot select `mirror` at all without a selector, and §5.5's enum is the thing that makes
`mirror + redistribute` unrepresentable — which §6 says PR A must guarantee. The enum therefore
leads, in its own PR:

| Stage | Content | Verified by |
|---|---|---|
| **1** | `IrregularDeposition` enum, threaded to `AmrMesh`/`EBAMRParticleMesh`; input migration; `Mirror` present but hard-errors "not implemented" | **Every regression bit-for-bit identical.** No physics in this PR, so that one property covers it |
| **2** | Per-band-cell surface data `(status, x_c, n_c, S_c)` built at regrid; the curvature fit; ghost/geometry requirements; counters | Unit test against analytic sphere/cylinder/torus curvature; no consumer yet, so regressions still bit-for-bit |
| **3** | The mirror pass: image container, scratch field, build/remap/deposit/add, guards; `ItoSolver` wired | §6.1 acceptance suite; regressions still bit-for-bit because nothing selects `mirror` |
| **4** | `PLAN.md`'s PR B: move `m_massDiff`/`m_depositionNC` onto `PhaseRealm`, move hybrid-and-redistribute into `AmrMesh`, delete the duplicated solver code, retire the ad-hoc `EBParticleMesh` constructions | Regressions bit-for-bit; the deleted code has no remaining callers |

Stages 1, 2 and 4 are each independently required to leave **every** existing regression bit-for-bit
unchanged. That is the cheap, total check `PLAN.md` §6 identifies, and splitting the enum out is what
makes it diagnostic rather than merely reassuring: a Stage-1 diff means the mechanical conversion is
wrong, a Stage-3 diff means the mirror pass leaked into a path that did not select it.

### 0.2 New types introduced by the whole feature

Per the working agreement in `CLAUDE.md`, stated up front with what already exists and why it does
not serve. **The whole feature introduces exactly one new type.**

| New | Use count | What already exists | Why it does not serve |
|---|---|---|---|
| `enum class IrregularDeposition` (`CD_IrregularDeposition.H`) | ~30 signature sites + 3 solvers | `DepositionType` (the kernel), and three independent `bool`s: `ItoSolver::m_forceIrregDepositionNGP`, `m_useRedistribution`, `m_blendConservation` | The three booleans are answers to *one* question and their cross-product is mostly meaningless. `mirror + redistribute` is representable in booleans and is a silent 1/kappa error (`PLAN.md` §6). One enum makes it unrepresentable. Follows `ParticleManagement::ParticleMergeMethod` + `mergeMethodFromString` exactly |

Everything else reuses existing containers, deliberately:

- Surface data → `EBAMRCellData` (`PLAN.md` §5.1 decided this against a `BaseIVFAB`; do not re-propose).
- Scratch field → `EBAMRCellData`, one component.
- Image container → `ParticleContainer<NoPayload>`; `NoPayload` is `CD_ParticleSoA.H:288` and is the
  container's default template argument.
- Redistribution scratch → `EBAMRIVData`, moved rather than created (Stage 4).
- **No `std::map`, no `IntVect`-keyed container, no per-cell hand-rolled store anywhere.**

Two additions that are *not* types and are called out so a reviewer can object early:

- `namespace MirrorDeposition` in `CD_MirrorDeposition.H` — component-index constants and three free
  `inline` functions (reflect, Jacobian, shape-operator lift). Same shape as `ParticleManagement`. It
  exists so §2.2/§3.1's twenty flops live in one testable place instead of being open-coded in
  `CD_AmrMeshImplem.H`. No state, no struct.
- `PhaseRealm::setRealmName` — a one-line setter (see §3.2), not a type.

### 0.3 Corrections to `PLAN.md` carried into implementation

Three, all verified against the tree.

1. **The behaviour-preserving default is `ngp`, not `native`.** §5.5 says *"The default in
   `CD_ItoSolver.options` becomes the value that reproduces today's shipped default
   (`false/false` -> `native`)"*. `CD_ItoSolver.options:21` ships `irr_ngp_deposition = true` with
   `redistribute = false`, which maps to `IrregularDeposition::NGP` by §5.5's own table. The two
   files carrying `true` are `Source/ItoDiffusion/CD_ItoSolver.options` and
   `Exec/Examples/ItoKMC/WireWire/example.inputs`. Setting the default to `native` would silently
   change every input that does not override it. **Default becomes `ngp`.** The cross-tabulation
   otherwise reproduces §5.5 exactly: 5 native / 22 redistribute / 2 ngp.

2. **The exact criterion of §1.2 is not evaluated at runtime, and must not be.** §1.2 says the
   criterion *"falls out of the deposit loop for free"*. Taken literally as a build-time test it is
   not free and not even possible: the image's cell sits up to Chebyshev 4 (CIC) / 5 (TSC) from the
   particle's cell (§1.2's `img_reach.py` table), so its cloud support reaches ~6 cells outside the
   valid box, where `EBISBox` is invalid at `eb_ghost = 2`. It is free in the other sense: an image
   whose cloud touches no `kappa > 0` cell deposits **only** into covered cells, and §4.1's estimator
   is read only at `kappa > 0`, so building it is already equivalent to discarding it.
   **Build an image for every band particle with `status >= 1` and `d > 0`; no overlap test.**

3. **The DEBUG mask-sufficiency check of §1.2 is replaced by a band-radius convergence test.** The
   outer-ring counter §1.2 asks for needs the same criterion evaluated at Chebyshev 3, and hence the
   same out-of-reach geometry. Instead the band radius becomes a runtime option
   (`AmrMesh.mirror_band_radius`, default 2, requiring `num_ghost >= radius`), and the acceptance
   suite runs the torus case at radius 2 and radius 3 and asserts the `kappa`-binned densities agree
   to round-off. That measures the thing the counter is a proxy for, needs no geometry outside the
   grown region, and costs nothing in production. See §5.4.

---

## 1. Stage 1 — `IrregularDeposition` and the input migration

**Nothing in this stage changes behaviour.** Its whole content is a signature change plus a
value-domain change in input files.

### 1.1 New header

`Source/Particle/CD_IrregularDeposition.H`, beside `CD_DepositionType.H`, with
`CD_IrregularDepositionImplem.H` included at the foot (the `CD_ParticleManagement.H:481` pattern).

```cpp
enum class IrregularDeposition
{
  Native,              ///< Deposit as-is; phi is an extended state into the EB.
  NGP,                 ///< Whole cloud into the cut cell.
  Mirror,              ///< Even extension across the EB. A cut cell then holds n, not kappa*n.
  Redistribute,        ///< Hybrid update, deltaM = (1-kappa)*dc smooshed to neighbours.
  RedistributeBlended  ///< As above, blended with the non-conservative divergence.
};

inline IrregularDeposition
irregularDepositionFromString(const std::string& a_str) noexcept;   // aborts on unknown selector
```

`Mirror`'s Doxygen carries the §7 O2 contract verbatim — *under `mirror`, a cut cell holds `n`, not
`kappa*n`* — together with the named consumers that will need attention (`DataOps::filterSmooth`,
`CD_DataOps.cpp:679`; the `arithmeticAverage` sites). That enumerator comment is the deliverable O2
asks for.

Selector strings: `native`, `ngp`, `mirror`, `redistribute`, `redistribute_blended`.

### 1.2 Signature conversion — deposition side only, ~30 sites

`bool a_forceIrregNGP` becomes `const IrregularDeposition a_irregularDeposition` in:

| File | Sites | Functions |
|---|---|---|
| `Source/AmrMesh/CD_AmrMesh.H` | 6 | `depositWeight`, `depositGathered`, `depositParticles` declarations |
| `Source/AmrMesh/CD_AmrMeshImplem.H` | 6 | the same definitions |
| `Source/Particle/CD_EBAMRParticleMesh.H` | 12 | `depositWeight`, `depositGathered`, `deposit<Members...>`, and the four `deposit*Core` helpers |
| solvers | 6 | `ItoSolver` (member, parse, 2 call sites), `McPhoto`, `TracerParticleSolver` |

**`EBParticleMesh` keeps its `bool`** (`PLAN.md` §5.5). At the point where
`EBAMRParticleMesh::depositWeight`/`depositGathered` build their `patch`/`patchNGP` lambdas
(`CD_EBAMRParticleMesh.H:261-270, 329-352`), collapse the enum:

```cpp
const bool forceIrregNGP = (a_irregularDeposition == IrregularDeposition::NGP);
```

`Mirror`, `Redistribute` and `RedistributeBlended` all map to `false` there — the mirror is a
separate pass and redistribution happens after the deposit.

**The 36 interpolation-side occurrences are untouched.** `EBParticleMesh` and `EBAMRParticleMesh`
reuse the parameter name `a_forceIrregNGP` on `interpolateWeight`/`interpolate`; `irr_ngp_interp`
keeps its boolean. Counted per file: `CD_EBParticleMesh.H` 22, `CD_EBAMRParticleMesh.H` 6,
`CD_AmrMesh.H` 4, `CD_AmrMeshImplem.H` 4. The 7 `m_forceIrregInterpolationNGP` members likewise.

### 1.3 The three hardcoded call sites

- `CD_ItoKMCStepperImplem.H:4041` passes `a_forceIrregNGP = false` -> `IrregularDeposition::Native`.
- `CD_ItoKMCStepperImplem.H:6134` — convert to whatever its literal maps to.
- `CD_ItoKMCGodunovStepperImplem.H:1338` passes `true` -> `IrregularDeposition::NGP`. This is the
  `hasDielectrics` branch of `computeSemiImplicitRho`, depositing onto `phase::solid`. **Convert the
  literal and add the comment saying why it is a literal**: `phase::solid`'s `ebisbox.normal()`
  points into the *solid*, so routing a solver's setting here would mirror across the same EB with
  the fluid side swapped and push particles deliberately placed inside the dielectric back into the
  gas. `PLAN.md` §6 PR B: *"it has to resist deciding."*

`CD_CdrSolver.cpp:1316-1319` is a `constexpr bool forceIrregNGP = true` on a direct per-patch
`EBParticleMesh::depositWeight`. It never passes through the funnel and is left alone.

### 1.4 Solver parsing

`ItoSolver::parseDeposition` (`CD_ItoSolver.cpp:251`):

- Replace `pp.get("irr_ngp_deposition", m_forceIrregDepositionNGP)` (`:320`) with
  `pp.get("irregular_deposition", str); m_irregularDeposition = irregularDepositionFromString(str);`
- Delete `pp.get("redistribute", ...)` and `pp.get("blend_conservation", ...)` and the members they
  filled; `redistributeAMR` (`:2133`) gates on
  `m_irregularDeposition == Redistribute || == RedistributeBlended`, and `depositNonConservative`
  (`:2178`) on `== RedistributeBlended`.
- **The `mirror` + NGP cross-check goes here**, thirty lines below the `deposition` parse, where it
  can name both offending keys:

```cpp
if (m_deposition == DepositionType::NGP && m_irregularDeposition == IrregularDeposition::Mirror) {
  MayDay::Error("ItoSolver::parseDeposition - '<class>.deposition = ngp' is incompatible with "
                "'<class>.irregular_deposition = mirror'; NGP puts the whole cloud in the "
                "particle's own cell, so an image in the same cell doubles it.");
}
```

  plus a cheap defensive `CH_assert` on the same condition at `AmrMesh::depositWeight` for the
  literal call sites.

- **Stage 1 only:** `Mirror` also hard-errors here with "not implemented until stage 3". Removed in
  Stage 3.

`McPhoto::parseOptions` (`CD_McPhoto.cpp:236`): `blend_conservation` -> `irregular_deposition`.
`m_blendConservation` gated `depositHybridDivergence`'s redistribution at `:1339` while the hybrid
itself ran unconditionally at `:1288/:1321` — with `blend = false` that hybrid is the identity
`divH = dc + (1-kappa)*0` and its `deltaM` is discarded, so all 14 shipped files map to `native` and
the dead work disappears. Record in the PR description that McPhoto *gains* `Redistribute`, a mode it
could not previously express, available but unselected.

### 1.5 Input migration — 29 + 14 files

| `irr_ngp_deposition` | `redistribute` | files | becomes |
|---|---|---|---|
| `false` | `true` | 22 | `irregular_deposition = redistribute` |
| `false` | `false` | 5 | `irregular_deposition = native` |
| `true` | `false` | 2 | `irregular_deposition = ngp` |

Verified by cross-tabulation over the tree; `blend_conservation` is `false` in all 29. The 22 are the
`Exec/Tests/BrownianWalker/*` pairs; the 5 are `Exec/Tests/ItoKMC/JSON/*` and three
`Exec/Examples/ItoKMC/*`; the 2 are `CD_ItoSolver.options` and
`Exec/Examples/ItoKMC/WireWire/example.inputs`.

Also: delete the now-dead `ItoSolver.redistribute` / `ItoSolver.blend_conservation` lines from every
input that sets them; rewrite `McPhoto.blend_conservation` -> `McPhoto.irregular_deposition = native`
in `CD_McPhoto.options` and the 14 inputs; update `Docs/Sphinx/source/**/Ito.rst`.

**Do not sweep in the fluid solvers.** `CdrGodunov.blend_conservation` (7 files) and
`CdrCTU.blend_conservation` (34 files) are the fluid hybrid divergence, unrelated despite the name.

Because the tree uses `pp.get` and never `pp.query`, a stale input hard-errors on the missing key.
That is the failure mode we want — **no fallback, no deprecation shim.**

### 1.6 Stage 1 verification

1. `pre-commit run --all-files`.
2. Rebuild Sphinx and check for new warnings — `.options` files are pulled in by `literalinclude`
   and the rewrites shift line counts:
   `cd Docs/Sphinx && python3 -m sphinx -W --keep-going -b html source build/html`
3. Build 2-D and 3-D, DEBUG and OPT.
4. **Run the full regression set and diff plot files bit-for-bit against the pre-branch build.** Any
   difference is a conversion bug.

---

## 2. Stage 2 — surface data at regrid

### 2.1 Component layout

`Source/Particle/CD_MirrorDeposition.H`, namespace `MirrorDeposition`:

```
comp 0                     status   0 = cannot reflect, 1 = fitted S_c, 2 = J = 1 fallback
comp 1        .. D         x_c      boundary centroid, PHYSICAL position
comp D+1      .. 2D        n_c      unit normal, INTO the fluid
comp 2D+1     .. 2D+T      S_c      symmetric DxD, upper triangle (T = 3 in 2-D, 6 in 3-D)
```

**13 components in 3-D, 8 in 2-D.**

`status` carries three values rather than the plain `ok` flag of §2.1 because §3.2/§3.3's fallback is
*"deposited with `J = 1`, not dropped"* — a crease or a rank-deficient neighbourhood must still
reflect. `status = 0` is reserved for §5.4's genuine failure: a band cell that found no *delivered*
cut cell. Counting the two separately is what makes O9's diagnostic meaningful. `S_c = 0` gives
`J = 1` exactly through §3.1's formula, so status 2 needs no branch on the hot path — it is a
diagnostic label, and the same as a plane.

Free functions in the same namespace (declared here, `inline` in `CD_MirrorDepositionImplem.H`):

```cpp
inline void liftShapeOperator(const RealVect& t1, const RealVect& t2, const Real S[2][2], Real* Sc);
inline bool reflect(const RealVect& a_pos, const Real* a_surfaceData, const Real a_maxJacobian,
                    const Real a_minDenominator, RealVect& a_image, Real& a_jacobian,
                    MirrorDeposition::Refusal& a_why);
```

`reflect` is §2.2 plus §3.1 plus §3.2's guard, and is where the acceptance evidence lives. It does
no I/O, takes no mesh, and is directly unit-testable against `mirror_frame.py`'s tables.

### 2.2 `PhaseRealm` additions

```cpp
// Regrid-lifetime, so not mutable -- the m_irregularCells convention (CD_PhaseRealm.H:600).
EBAMRCellData m_mirrorSurfaceData;   ///< 13 comps in 3-D, 8 in 2-D; see MirrorDeposition.
Vector<bool>  m_mirrorHasCutCells;   ///< Per level, ALL-REDUCED. Gates the pass; see 3.4.
```

plus a new operator string beside `s_particle_mesh` (`CD_PhaseRealm.H:60`):

```cpp
static const std::string s_mirror_deposition = "mirror_deposition";
```

registered through `PhaseRealm::registerOperator` (`CD_PhaseRealm.cpp:353`, add to the validity list
at `:362-364`), so realms that never mirror pay nothing.

**Built in `regridOperators`, not `regridBase`** (`PLAN.md` §5.1). Add
`this->defineMirrorSurfaceData(a_lmin)` beside `defineParticleMesh()` (`CD_PhaseRealm.cpp:327`),
gated on `queryOperator(s_mirror_deposition)`, with its own `timer.startEvent`/`stopEvent` and
`MemoryReport` bracket like every other operator there. `defineMasks` is unconditional and sits in
`regridBase` for a stated reason (`:182-183`) that has nothing to do with curvature.

### 2.3 Requirements checked at define time

On entry to `defineMirrorSurfaceData`, abort naming both parameters:

```
m_numEbGhostsCells >= 2                   // EBISBox valid over the 5^3 fit stencil
m_numGhostCells    >= m_mirrorBandRadius  // default 2; exchange must deliver cut-cell data
```

`CD_AmrMesh.cpp:2788/2805` validate both only as `>= 0`, and `CD_AmrMesh.options` ships 2 for each.
Every shipped input already satisfies both — the check guards a future input file, not a migration.
**Keep the two separate; do not derive one from the other.**

Move to a 7³ fit later and `eb_ghost` becomes 3 while `num_ghost` and the band stay 2.

### 2.4 The build, in order

`PhaseRealm::defineMirrorSurfaceData(int a_lmin)`, per level, allocated with `EBCellFactory(ebisl)`
and `m_numGhostCells*IntVect::Unit`, `setValue(0.0)` first so an unwritten cell is `status = 0` by
default rather than by accident.

**Pass A — seed the cut cells.** Over the **valid** box only, for each single-valued cut cell:

```cpp
n_c = ebisbox.normal(vof);
x_c = m_probLo + Location::position(Location::Cell::Boundary, vof, ebisbox, dx);   // CD_LocationImplem.H:37-51
status = 2;   // provisional: J = 1 until the fit succeeds
```

`Location::position` returns `(gridIndex + 0.5 + bndryCentroid)*dx` — the `ret *= a_dx` at `:49` is
the point of the citation. `CD_AmrMeshImplem.H:584` already calls it this way. **Do not difference
raw `bndryCentroid` values from different cells**; that drops the inter-cell offset and inflates the
fitted curvature ~10x (`PLAN.md` §2.1).

**`m_mirrorSurfaceData[lvl]->exchange()`.**

**Pass B — fit the shape operator.** Over the valid box, for each cut cell with `status >= 1`:

1. `t1, t2` orthonormal perpendicular to `n_c`, built by crossing `n_c` with the axis of smallest
   `|n_c[dir]|`. **Local to the fit and discarded here** (`PLAN.md` §2.1, §3.3).
2. Gather neighbours `j` over the 5^D window, `j != i`, taking `x_j`, `n_j` **from the field, not
   from `EBISBox`**, and accepting only `status_j >= 1`. Reading the field is what enforces §5.4's
   *delivered, not merely existing* rule: a cut cell outside the level's grids was never written and
   never exchanged, so it is `status = 0`. Skip multi-valued cells (`ebisbox.numVoFs(iv) > 1`).
3. `u_j = (t1.dx_j, t2.dx_j)`, `v_j = (t1.dn_j, t2.dn_j)`, `dx_j = x_j - x_i` (both physical),
   `dn_j = n_j - n_i`.
4. Require **>= 4** usable neighbours (this is what the harness enforces; `PLAN.md` §3.3).
5. `A = U^T U` (2x2 symmetric), `B = U^T V`. **Conditioning test on `A`** before solving: closed-form
   eigenvalues of a 2x2 symmetric, reject when `lambda_min/lambda_max < mirror_cond_tol`. Four
   neighbours can be nearly collinear in the tangent plane — a ridge, a fin, a cylinder whose cut
   cells run along the axis — and the count is necessary, not sufficient.
6. Solve `S^T = A^{-1} B`, symmetrize `S <- (S + S^T)/2`.
7. **Residual as a crease detector**: `res = dx * sqrt(sum_j |v_j - S u_j|^2 / sum_j |u_j|^2)`,
   dimensionless. Reject above `mirror_crease_tol`.
8. On success `S_c = P S P^T`, `P = [t1 t2]`, stored as the upper triangle, and `status = 1`. On any
   rejection leave `S_c = 0` and `status = 2`, and count which rejection it was.
9. Under `DEBUG`, assert `|S_c n_c| < 1e-12` — `mirror_frame.py` measures `2e-17` over 2000 draws.

> **Existing machinery considered and not used.** `LeastSquares` (`CD_LeastSquares.H`) is
> stencil-oriented: it returns `VoFStencil`s for scalar fields on cut cells, not a small-matrix fit
> of a tensor from displacement/normal-difference pairs. `LaPackUtils::computeSVD`
> (`CD_LaPackUtils.H:137`) does fit the shape, and its singular values would give the conditioning
> test directly, but it allocates `std::vector`s per call inside what is an OpenMP loop over cut
> cells. A 2x2 symmetric eigenproblem and a 2x2 inverse have exact closed forms in about fifteen
> allocation-free lines. Use the closed form; cite `computeSVD` in the comment so the next reader
> knows it was considered.

**`m_mirrorSurfaceData[lvl]->exchange()`.**

**Pass C — extend to the band.** Over the valid box, for each non-covered cell that is *not* a cut
cell and lies within `mirror_band_radius` of one:

- Search the `(2r+1)^D` window for cells with `status >= 1`; pick the one minimising
  `|cellCentre - x_c(neighbour)|^2`, where `cellCentre = probLo + (iv + 0.5)*dx`.
- **Tie-break on the lowest `IntVect` in lexicographic order**, so the result does not depend on
  iteration order.
- Copy `(x_c, n_c, S_c)` and the neighbour's `status`; on no candidate leave `status = 0` and count.

This is `argmin` over the neighbourhood, not two 1-cell propagation sweeps. §2.4's error ladder — the
1.2–2.9% that justifies the quadratic patch — was measured with `near = argmin(|x_band - x_cut|^2)`
(`mirror_source.py:48-53`). Sweeps compute a Chebyshev-nearest assignment with an arbitrary
tie-break, and the two differ exactly where more than one piece of surface is within two cells: the
O3 geometry. 125 distance comparisons per band cell once per regrid is the price of being faithful to
the evidence.

**No exchange after pass C.** Surface data is read only at a particle's own cell, and particles live
in the valid region of the patch that owns them.

**Ordering is load-bearing.** Extend first and exchange once afterwards and a band cell at a patch
edge whose nearest cut cell lies one patch over is assigned a too-far cut cell, and the exchange
cannot repair it — both patches computed from the same incomplete data and agree on the wrong answer
(`PLAN.md` §5.4). This is the one place `defineMasks`' valid-then-exchange pattern does not transfer.

Finally, set `m_mirrorHasCutCells[lvl]` from the local cut-cell count **all-reduced** with
`ParallelOps::sum` — see §3.4 for why it must be global.

### 2.5 Counters and reporting

Accumulated per thread, reduced over OpenMP, then `ParallelOps::sum` across ranks. Reported at
`m_verbose`, once per regrid:

| Counter | Meaning | Reading |
|---|---|---|
| `fitOk` | cut cells with `status = 1` | |
| `fitTooFewNeighbours` | fewer than 4 usable | under-resolved geometry |
| `fitIllConditioned` | `A` near-singular | ridge / fin / axis-aligned cylinder |
| `fitCrease` | residual over threshold | sphere union, facet edge, dielectric corner |
| `bandNoCutCell` | §5.4 / O9: valid band cell, no delivered cut cell | **diagnostic, not an abort** — says how much of the band sits at a grid edge |

A run in which the fit fallbacks are not rare has under-resolved geometry. That is a tagging problem,
not a deposition problem (`PLAN.md` §3.2).

**Report the memory** alongside the existing four masks in the memory report. At the shipped defaults
(`coarsest_domain 128^3`, `max_block_size 16`, `num_ghost 2` -> 1.91x ghost overhead) one component
costs ~33 MB per level, so 13 components cost **~430 MB per level, per phase, per realm**, against
131 MB for the four masks together. This is accepted with eyes open (`PLAN.md` §5.1) — the point is
that it appears as a labelled line rather than as unexplained growth.

### 2.6 New options

Added to `CD_AmrMesh.options` with the three-line banner intact:

```
AmrMesh.mirror_band_radius   = 2        ## Band half-width in cells. Requires num_ghost >= this.
AmrMesh.mirror_cond_tol      = 1.E-3    ## Min eigenvalue ratio of U^T U in the curvature fit.
AmrMesh.mirror_crease_tol    = 1.E-1    ## Max dimensionless fit residual before falling back to J = 1.
AmrMesh.mirror_max_jacobian  = 1.E4     ## Refuse (J = 1) above this. Must exceed ~500; see below.
AmrMesh.mirror_min_denom     = 1.E-4    ## Refuse (J = 1) below this. Must be under ~4E-3; see below.
```

**Both Jacobian guards are sized from a measurement, not a round number** (`PLAN.md` §3.2). Under
§1.2's criterion a concave cavity of `R = 3*dx` legitimately produces `J = 28.9` (CIC) and `J = 483.8`
(TSC), the latter at a denominator of `(1 - d/R)^2 = 4.1e-3`. A guard firing below `J ~ 10^3`, or on a
denominator above `~10^-3`, discards real images in a geometry this tree ships. `1e4` and `1e-4` clear
both with a decade of margin.

### 2.7 Stage 2 verification

A unit test under `Exec/Tests/ItoDiffusion/MirrorSurfaceData` (new module, see §5.1), 2-D and 3-D:

1. **Analytic curvature.** Sphere `R = 4, 6, 8 dx` convex and `R = 6, 8 dx` concave; check `2H = tr S_c`
   against `+-2/R` binned by `kappa`. Budget from `PLAN.md` §3.3: **2–3% in the `kappa <= 0.05` bin at
   tight radii, and that is systematic, not spread.**
2. **Anisotropy.** Cylinder `R = 6*dx`: principal curvatures `1/R` and `0`, so `K = 0` and
   `2H = 1/R`. This is the only case that can see a frame error — every isotropic `S` is
   frame-invariant. Cross-check against `mirror_frame.py`.
3. **`K` is the second invariant.** Assert `det S_c` is identically zero and that
   `K = (1/2)[(tr S_c)^2 - tr(S_c^2)]` recovers `1/R^2` on the sphere and `0` on the cylinder in both
   dimensions. `K = det S_c` silently degrades the exact Jacobian to the linearized one, whose error
   §3.3 measures at 48–77%.
4. **Sign convention** (`PLAN.md` §2.3). Assert `2H > 0` for a convex solid and `J < 1` there.
   Reversing the convention replaces `J` by `1/J` — 0.111 against 9.0 at `R = 4*dx`, `d = 2*dx` — with
   no crash and no NaN. One assert, one comment.
5. **Extension.** On the torus, check that band cells at Chebyshev 1 and 2 carry the nearest cut
   cell's data by the Euclidean rule and that the tie-break is deterministic across a rerun with a
   different box decomposition.
6. Regressions still bit-for-bit — nothing reads the field yet.

---

## 3. Stage 3 — the mirror pass

### 3.1 `PhaseRealm` additions

```cpp
// Per-call scratch, so mutable, following the file's own convention (eighteen such members, :605-695).
mutable EBAMRCellData                m_mirrorScratch;      ///< 1 comp; 4.2 step 3.
mutable ParticleContainer<NoPayload> m_mirrorImages;       ///< 4.2 step 1; defined at regrid.
mutable bool                         m_depositionInFlight; ///< DEBUG invariant; see below.
```

`m_mirrorScratch` is one component because `depositWeight` asserts `a_meshData[0]->nComp() == 1`
(`CD_EBAMRParticleMesh.H:639`).

**Write the reason for the scratch field where it is allocated**, because it reads like removable
weight to someone optimising later. It is not just that every AMR deposit opens with
`DataOps::setValue(a_meshData, 0.0)` (`:641, 844, 943, 1217`). The stronger reason is that the level
exchange is a *reversed* copier with an add op (`:591-592`, used at `:662`): the ghost region is the
exchange's source and the owner's valid region its destination, and **nothing writes the ghosts back
to zero**. A second deposit into the same field without a reset folds the first pass's ghost content
into the valid region again — silent double counting proportional to the ghost-region share of each
patch. The reset and the exchange are one mechanism; a `bool a_reset` parameter would be a bug.

**Also document why a ghost-width mismatch between scratch and `a_phi` is not a bug.**
`EBCellFAB::plus(src, scale)` adds over `a_src.getRegion() & getRegion()` and asserts only on
component counts, so `DataOps::incr` loses nothing in the valid region and merely skips the outer
ring — which is correct, because after step 3 only the scratch's valid region is meaningful anyway.

**`m_depositionInFlight`** is the invariant `PLAN.md` §5.1 asks for. These members are shared by every
solver on that (realm, phase), where today each solver owns its own scratch. That is safe only
because solvers deposit serially — `ItoLayout`'s iterator is a serial loop and the OpenMP parallelism
inside a deposit is over boxes. `mutable` removes the last compile-time hint that anything is being
written. Assert on entry and exit; it costs nothing under `OPT=HIGH` and fires the day someone
parallelises the solver loop.

### 3.2 Plumbing `m_mirrorImages.define`

`ParticleContainer::define` (`CD_ParticleContainer.H:193-202`) needs four things `PhaseRealm` does
not hold: `minBlockSize`, `levelTiles`, `validCells` and the realm name. All four live on `Realm`
(`CD_Realm.H:490, 550, 593`) or `AmrMesh`, all are phase-independent, and `Realm` owns the
`PhaseRealm`s. `AmrMesh::allocate(ParticleContainer<P,Traits>&, realm)` (`CD_AmrMeshImplem.H:238-261`)
shows exactly which getters to forward.

One hop, three small edits:

1. `Realm::setName(const std::string&)` + `m_name`, called from
   `AmrMesh::regridOperators(a_realm, a_lmin)` (`CD_AmrMesh.cpp:1153-1165`), which already has the
   name. `Realm` currently does not know its own label; the label is diagnostics-only
   (`ParticleContainer::m_realm` is read back only by `getRealm()` for the solvers'
   `CH_assert(getRealm() == m_realm)`).
2. `PhaseRealm::defineMirrorImages(name, minBlockSize, levelTiles, validCells)`, called from
   `Realm::regridOperators` (`CD_Realm.cpp:179-188`) after the phase loop. `Realm::regridBase` builds
   `m_validCells` and `m_levelTiles` at `:174-175`, before `regridOperators` runs, so both are current.
3. Gate on `queryOperator(s_mirror_deposition)`.

**Its lifecycle is three lines and all three are load-bearing** (`PLAN.md` §5.1):

- `define`d at regrid, because `remap()` reads `levelTiles` and `validCells`, which the regrid rebuilds.
- `clearParticles()`d at the **start** of every pass, not the end, so an aborted pass cannot leak
  images into the next deposit.
- Never `resetParticleIDs`. `ParticleSoA::append(pos, weight)` leaves the metadata columns at
  `s_invalidID`/`-1` (`CD_ParticleSoA.H:949-958`); nothing in `remap()` reads them and images need no
  identity.

### 3.3 The pass

One templated private helper in `CD_AmrMeshImplem.H`, called by both funnels:

```cpp
template <typename P, typename Traits, typename Strength>
void
AmrMesh::mirrorPass(EBAMRCellData&                      a_meshData,
                    const std::string&                  a_realm,
                    const phase::which_phase&           a_phase,
                    const ParticleContainer<P, Traits>& a_particles,
                    const DepositionType                a_depositionType,
                    const CoarseFineDeposition          a_coarseFineDeposition,
                    Strength                            a_strength);
```

`a_strength(leaf, i) -> Real` is `leaf.weight(i)` from `depositWeight` and the caller's
`a_cellGather` from `depositGathered` — the same `(const ParticleSoA<P,Traits>&, std::size_t) -> Real`
contract `EBAMRParticleMesh::depositGathered` already documents (`CD_EBAMRParticleMesh.H:305`).

Body, exactly `PLAN.md` §4.2:

```
0. if (!anyLevelHasCutCells) return;                     // 3.4 -- collective, decided at regrid
   m_mirrorImages.clearParticles();

1. for lvl, for din:                                     // VALID holder only
     for i in a_particles[lvl][din]:
       iv = ParticleOps::getParticleCellIndex(pos, probLo, dx);
       CH_assert(!ebisbox.isMultiValued(iv));            // O10
       if (status(iv) == 0) continue;
       if (!MirrorDeposition::reflect(pos, &sd(iv,0), jMax, denMin, image, J, why)) { count(why); J = 1; }
       if (d <= 0) continue;                             // 1.2 requires d_p > 0
       m_mirrorImages[lvl][din].append(image, a_strength(leaf, i) * J);

2. m_mirrorImages.remap();                               // the ORDINARY remap, no special level rule

3. particleMesh.depositWeight(m_mirrorScratch, m_mirrorImages, a_depositionType,
                              a_coarseFineDeposition, IrregularDeposition::Native);

4. DataOps::incr(a_meshData, m_mirrorScratch, 1.0);
```

**No covered-cell reset, and that is deliberate** (`PLAN.md` §4.2). Every image is inside the solid
by construction — that is what reflection means — and §4.1's estimator is defined only at
`kappa > 0`. The deposit's contract is that the field is correct where it is read, exactly as
`03cf59690` drew the line for synchronization. Four facts make deferring safe: the dominant
coarsening path is `kappa`-weighted (`ItoSolver::coarsenAndFillGhosts`, `CD_ItoSolver.cpp:2121`);
plot output resets covered cells itself (`:1813`); a covered cell has no `VolIndex`, so the elliptic
solvers never see it; and nothing zeroes covered cells after a particle deposit today either. **Do
not re-add a reset here to paper over an O2 item.**

**No per-particle normalization.** Issue #29 proposes scaling each particle by
`V = sum_{kappa>0}(w + w')`. That is wrong and must not be implemented — it reintroduces exactly the
bias the mirror removes (`PLAN.md` §4.1).

**Images are built from the valid holder only** (`PLAN.md` §4.3). `depositHaloCore:849` is a
`copyMaskParticles` — a copy — so building from a mask holder would generate a second image under
`Halo`, which is issue #29's double-counting trap in a new guise. All three coarse-fine cores restore
the valid holder before returning (`CD_EBAMRParticleMesh.H:912, 988, 1292`), so building images
*after* the real deposit is safe; the constraint is which holder, not when.

**The remap is the ordinary one.** A level-preserving remap was evaluated and refuted in both
directions (`PLAN.md` §4.2, §9): `Transition` already deposits transition particles at fine width
(`:1276`), and pinning an image to the coarse level puts it outside the 0.5–1-cell CF band where it
would be discarded. Promotion deposits the correction at the wrong bandwidth but keeps it — mass
conserved, mean correct, a variance defect.

### 3.4 The early-out must be collective

`remap()` is a collective and `DataOps::incr` is whole-domain: both are O(cells)/O(particles)
regardless of how thin the band is, and they are paid on every deposit of every species
(`PLAN.md` §7). Skipping them needs a **global** decision, or ranks deadlock in `remap()`.

`m_mirrorHasCutCells` is therefore all-reduced at regrid (§2.4), and the pass returns before step 0
when no level has cut cells anywhere. **Measure the residual fixed cost before assuming it is small
— `remap()` has been the dominant cost in this tree before.**

The reflect fractions in `PLAN.md` §7 (42–90%) are properties of the harnesses' 16³ box, not of a
simulation, and must not be read as a cost estimate. In a real run the band is two cells of a domain
hundreds of cells across and the fraction goes with surface-to-volume.

### 3.5 Wiring `ItoSolver`, and refusing the bad combination

`ItoSolver::depositWeight` is `AmrMesh::depositWeight` followed by `redistributeAMR(a_phi)`
(`CD_ItoSolverImplem.H:71, :79`); `depositGathered` does the same (`:148, :157`). `redistributeAMR`
opens with *"we had `a_phi = m_i/dV` but we actually want `phi = m_i/(kappa*dV)`"* — exactly the
assumption the mirror removes — and `depositHybrid` computes `deltaM = (1-kappa)*(dc - kappa*dnc)`
under *"Remember, `dc` already scaled by kappa"* (`CD_ItoSolver.cpp:2234`). With the mirror `dc`
is already a density, so `deltaM` is inflated by `1/kappa` — 20x at `kappa = 0.05` — and smooshed
into the neighbours.

Stage 1's enum already makes that unrepresentable: `Mirror`, `Redistribute` and
`RedistributeBlended` are one selector. **Nothing further is needed here beyond removing Stage 1's
"not implemented" error and asserting the invariant.** That is the whole reason the enum leads.

### 3.6 Documentation owed by this stage

- **`Mirror`'s enumerator Doxygen** carries the O2 contract (written in Stage 1, re-read here).
- **`depositHaloCore`'s Doxygen** records the `widthScale` mismatch (former O5): under `Halo` a
  promoted image deposits at `widthScale = 1` while its source gets `refRat`. `widthScale` is 1.0 at
  every other patch deposit in the class (`:658, 866, 966, 967, 1235, 1276`), all 60 shipped
  `deposition_cf` settings are `transition`, and the two hardcoded `Halo` sites are exactly the two
  that map to `Native`/`NGP` and never run a mirror pass. **Document it; do not engineer around it.**
- **`ItoSolver`'s plot-variable documentation** records O11: `depositWeightNGP`/`depositGatheredNGP`
  (`CD_ItoSolverImplem.H:84, :115`) build an `EBParticleMesh` per patch and never reach the funnel,
  so no mirror pass applies. They are plot paths and deliberately NGP — but the plotted particle
  density in a cut cell will then not match the `phi` the field solver sees, by exactly the factor
  the mirror introduces. **That is a diagnostic people will read as a bug. Name it.**
- **`Ito.rst`** gets the same contract statement plus the new selector's five values.

---

## 4. Stage 4 — wire it centrally (`PLAN.md`'s PR B)

**Not optional.** If the mirror pass stays the consumer's responsibility,
`irregular_deposition = mirror` becomes a setting some deposits honour and others silently ignore —
character-for-character the defect PR #700 just fixed. No assert can catch it, because a deposit
function cannot know whether its caller intends a mirror pass afterwards.

1. Move `m_massDiff` and `m_depositionNC` onto `PhaseRealm` as `mutable EBAMRIVData`. No new
   plumbing: `PhaseRealm` already owns `m_redistributionOp` and `m_nonConservativeDivergence` behind
   the `s_eb_redist` and `s_noncons_div` registrations (`CD_PhaseRealm.H:56-57`).
2. Move hybrid-and-redistribute into `AmrMesh::depositWeight`/`depositGathered`, beside the mirror
   pass, under `case Redistribute` / `case RedistributeBlended`.
3. Delete `ItoSolver::redistributeAMR`/`depositNonConservative`/`depositHybrid`
   (`CD_ItoSolver.cpp:2133, 2178, 2194`) and `McPhoto::depositHybridDivergence`/
   `depositNonConservative`/`depositHybrid` (`CD_McPhoto.cpp:1325, 1365, 1381`) — near-verbatim
   copies of each other.
4. **8 direct callers in 5 files** to re-point: `CD_McPhoto.cpp` x2 (`:1287, :1318`),
   `CD_TracerParticleSolverImplem.H:374`, `CD_ItoKMCStepperImplem.H` x2 (`:4041, :6134`),
   `CD_ItoSolverImplem.H` x2 (`:71, :148`), `CD_ItoKMCGodunovStepperImplem.H:1338`.
5. Retire the nine ad-hoc `EBParticleMesh` constructions — `McPhoto.cpp:1277, 1459`,
   `CdrSolver.cpp:1317`, `ItoSolver.cpp:2108, 2528, 2651, 2686, 2804`, `ItoSolverImplem.H:115` — by
   registering `s_particle_mesh` on those realms and using the persistent per-patch leaves. **Keep
   O11's two plot paths NGP**; retiring the construction does not change their kernel.

Regressions bit-for-bit again: `Redistribute` reproduces what `redistribute = true` did, and the
mirror is still unselected.

---

## 5. Acceptance

### 5.1 Where

New module `Exec/Tests/ItoDiffusion/` with `ItoDiffusion.ini` for `Exec/Tests/tests.py`, and:

- `Exec/Tests/ItoDiffusion/MirrorSurfaceData/` — Stage 2's curvature unit test (§2.7).
- `Exec/Tests/ItoDiffusion/MirrorDeposition/` — the three §6.1 cases, `regression2d.inputs` and
  `regression3d.inputs` each.

### 5.2 The three cases — mandatory, same PR

**Primary metric is the density binned by `kappa`.** Score cells **at least 3 cells inside** the
sampled region: a cell centre on the solid side is reconstructed almost entirely from images up to
`sqrt(3)*dx` away, and a 1-cell margin fabricates a ~1.8% deficit at small `kappa`. Margin 3 is a
defensible conservative choice, **not a measured floor** — there is no plateau at 3.

1. **A torus, and it must be the torus.** `Source/Geometry/CD_TorusSdf.H` exists. On a plane `J = 1`
   identically, so a planar-only test green-lights an implementation with no Jacobian and a 117%
   error at `R = 4*dx`. A sphere is barely better: umbilic, so `S_c` is isotropic and frame-invariant,
   and constant-curvature, so borrowing `(x_c, n_c, S_c)` from a neighbour costs nothing. A sphere
   cannot distinguish a correct anisotropic `S_c` from a wrong one — which is the gap §2.1 closes and
   the gap every harness in `Prototypes/MirrorDeposition/` shares. **Assert in the test that the
   sliver bin actually contains cells with `c_1 != c_2`**, or a coarsely-resolved torus degenerates
   to the sphere case without saying so.
2. **A concave case.** Convex-only hides the errors that no longer cancel, and the sliver bias is
   worst there.
3. **No EB at all**, asserted bit-for-bit identical to `native`.

**Report the `kappa <= 0.05` bin explicitly.** A 2–3% systematic offset there is expected and
budgeted (`PLAN.md` §3.3). A larger one, or one that grows with refinement, means the discrete
normals are worse than the noise sweep assumed.

`sum(kappa phi dV)/sum(m)` is a **secondary** check — cheap, exact, catches gross bookkeeping errors,
and reads 0.99974 on a build that is -46.5% wrong at `kappa <= 0.05`. Keep it as a second opinion, not
as the metric.

An **images-only build** (the zeroing-deposit bug) reads -97.5% to -100% in every regular cell, so
any of the three catches it; only the `kappa <= 0.05` bin is nearly blind (-1.8 to -7.7%).

### 5.3 Guard-path coverage

A tight concave cavity, `R = 3*dx`, run under both CIC and TSC, asserting the refusal counters are
zero at the §2.6 defaults. Under §1.2's criterion this geometry legitimately produces `J = 28.9`
(CIC) and `J = 483.8` (TSC) — the case that sizes the guard. A nonzero counter here means the guard
is refusing images the scheme is supposed to produce.

### 5.4 Band-radius convergence (replaces §1.2's DEBUG counter)

Run the torus case at `mirror_band_radius = 2` and `3` with `num_ghost = 3`, and assert the binned
densities agree to round-off. Chebyshev 2 is *derived* only for axis-aligned normals (O6) and
*measured* over ten tilted planar cases plus two sphere signs — good evidence, not a proof, and the
failure mode is a reflection silently not taken, which leaves no trace. This test is the trace.

See §0.3 item 3 for why this replaces the outer-ring counter rather than supplementing it.

---

## 6. Verification checklist, every stage

1. `pre-commit run --all-files` — `clang-format`, `clang-tidy`, `reuse`, `codespell`,
   `format-input-files`, `check-literalincludes`, `doxygen-check`.
2. `cd Docs/Sphinx && python3 -m sphinx -W --keep-going -b html source build/html`. The
   `sphinx-build` hook is `stages: [manual]`, so pre-commit does **not** cover this, and
   `literalinclude` line drift is partly silent. Any stage that reformats C++ or rewrites `.options`
   needs it. If a reformat shifts line counts, repair the directives by the marker-tracking procedure
   in `CLAUDE.md` §4 — do not fuzzy-match text across the reformat.
3. Build 2-D and 3-D, `DEBUG=TRUE` and `OPT=HIGH`, with and without MPI.
4. Full regression set, **diffed bit-for-bit** against the pre-branch build. Stages 1, 2 and 4 must be
   identical everywhere. Stage 3 must be identical everywhere except the new mirror tests.
5. `DEBUG` asserts exercised at least once: `|S_c n_c| < tol`, `!isMultiValued(particleIV)`,
   `!m_depositionInFlight`, and the `deposition = ngp` + `irregular_deposition = mirror` refusal.

---

## 7. Carried forward, unchanged

Open items O1, O3, O6, O7, O8, O9, O10, O11 stay as `PLAN.md` §7 dispositions them. Two are worth
restating because they land in code written here:

- **O9** — where a fine grid edge passes near the EB, an image from a fine cut cell can be demoted by
  `remap()` and a valid fine band cell can find no delivered cut cell. **Neither loses mass**:
  `addInvalidCoarseToFine` routes coarse-under-fine cloud mass back up *inside* the deposit, before
  any averaging. What is left is a bandwidth difference — the ordinary price of AMR. §2.5's
  `bandNoCutCell` counter is a **diagnostic, not an abort**.
- **O2** — `phi` changes meaning under `mirror`, and the consumers own that. What this work owes is a
  clear statement of the contract in three places (§3.6), not a fix to every reader.

And the scope limit `PLAN.md` states on all its evidence: every harness in
`Prototypes/MirrorDeposition/` is Cartesian, single-level, single-patch, analytic-geometry, with cut-cell
normals reconstructed from exact sub-sampled face area fractions. **Nothing there measures what
`GeometryShop`/`ScanShop` actually delivers.** §2.7 and §5.2 are where that number first appears; a
curvature error materially worse than the §3.3 noise sweep predicts is that gap showing, not a bug in
the fit.

**`Prototypes/MirrorDeposition/` is throwaway and must be deleted before the branch merges** — this
file included.
