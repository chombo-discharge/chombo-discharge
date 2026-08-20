# Review of PLAN.md revision 5 → hardening into revision 6 — findings

**Reviewed** `Prototypes/MirrorDeposition/PLAN.md` at revision 5, together with
`RevisionFindings.md` (review 5) · **Branch** `mirror_deposition` · **Code citations re-resolved
against** `b8721d71` · **Two harnesses changed, one added.**

This is review 6. Unlike reviews 1–5 it was asked to *harden the plan for implementation*, so its
output is a modified `PLAN.md` rather than a list of objections. This file records what changed and
why, so the next reviewer can attack the changes rather than rediscover them.

---

## 0. Verdict

**Review 5's three findings are all real. I confirmed each independently and all three are now
fixed in the plan** — two of them by a change of design rather than a change of wording.

Seven further items survive, of which four would have produced a defect in code written from
revision 5 and three are inconsistencies that would have wasted implementation time. The largest is
not a wrong number: it is that **§4.2 and §5.3 both use a mesh field that §5 never declares**, and
that the field §5 *does* declare costs ~400 MB per level per phase at the shipped defaults.

Ordered by what they cost if missed:

| # | Item | Class |
|---|---|---|
| F1 | The tangent frame is a stored quantity with no specification — **fixed by removing the frame** | blocker (review 5's Finding 2) |
| F2 | `m_mirrorScratch` is used by §4.2/§5.3 and declared nowhere | blocker |
| F3 | The surface data as a whole-domain cell field costs ~430 MB/level/phase/realm | **resolved: accepted** |
| F4 | Surface data undefined on a valid band cell at a level's grid edge | **fixed: band grown from the exchanged mask** |
| F5 | §7's cavity numbers are wrong in **both** directions | wrong number (review 5's Finding 3) |
| F6 | The `IrregularDeposition` conversion count, and where the `mirror` + NGP check goes | inconsistency |
| F7 | Step 5 cleans valid cells only; covered ghost cells stay polluted | **withdrawn — the reset came out** |
| F8 | §2.4's ladder measured a nearest-cut rule §5.4 does not propose | evidence/implementation mismatch |
| F9 | Reflect fractions quoted as a cost estimate are harness-box geometry | misleading number |
| F10 | Multi-valued band cells (**closed by invariant**), and deposits that bypass the funnel | unlisted open items |

Everything in §1 of `RevisionFindings.md` — twelve tables, fifteen counts, and the derivations — I
re-checked by sampling rather than exhaustively, and every sample reproduced. `reach_cells.py` still
gives CIC 91.596% / TSC 87.163% at Chebyshev 1 and max `d` 2.570 / 3.390. The five nits are all
correct, including the eighteen-versus-seventeen `mutable` count, which I recounted line by line.

---

## 1. The findings

### F1 — the tangent frame, removed rather than specified

Review 5's Finding 2 is right, and its own proposed fixes ("write the `t1`/`t2` rule into §2.1 as
normative", or "store `t1` explicitly, 12 comps instead of 9") both leave a frame in the format.
There is a third option that removes the question.

Store the shape operator as the symmetric **3×3 world-frame** tensor `S_c = P S P^T`, `P = [t1 t2]`.
Then

```
d    = eta + (1/2) w . (S_c w)          w = p - x_c,  eta = n_c . w
nhat = normalize( n_c + S_c w )
2H   = tr S_c        K = (1/2)[ (tr S_c)^2 - tr(S_c^2) ]
```

`S_c n_c = 0` by construction, so the normal component of `w` contributes nothing and no tangential
projection is ever taken. The tangent frame becomes local to the fit and is discarded there.

Measured (`mirror_frame.py`, new; cylinder `R = 6*dx`, principal curvatures `1/R` and `0`, 4000
draws):

| | median | p95 | max |
|---|---|---|---|
| `\|d(image)\|`, frame-free vs framed | 0.0 | 6.7e-16 | **1.8e-15** |
| `\|d(image)\|`, **mismatched** frame | **0.397 dx** | 1.36 dx | **2.24 dx** |
| `\|d(tr S)\|` under the same rotation | 2.8e-17 | 5.6e-17 | 8.3e-17 |
| `\|d(det S)\|` under the same rotation | 7.3e-20 | 1.4e-18 | 3.5e-18 |
| `\|S_c n_c\|` | 3.2e-18 | 1.0e-17 | 1.6e-17 |

Read: the reformulation is **exact**, not an approximation. And a frame mismatch moves the image by
up to 2.24 dx while `J`, and therefore §6.1's mass check, do not move at all — which is why no
existing harness could see it. Review 5 was also right that none of them *can*:
`mirror_source.py:quad_reflect` takes a scalar `kap_princ`, so its `S` is isotropic and
frame-invariant, and the anisotropic 2×2 the curvature fit produces is never fed into a reflection
anywhere in this directory.

**Two bonuses, and one new trap.** The world-frame form makes revision 5's 2-D special case —
"`S` is 1×1, set `K = 0` explicitly" — unnecessary: in 2-D `S_c = c t t^T`, so `tr S_c = c`,
`tr(S_c^2) = c^2`, and the second-invariant formula returns `K = 0` on its own. One expression, both
dimensions. The trap is that **`det S_c` is identically zero** — `S_c` is rank ≤ 2 — so carrying
`K = det S` across from the 2×2 silently degrades the exact Jacobian to the linearized one, whose
error §3.3 measures at 48–77%. That warning is now in §3.1 in bold.

Cost: 6 Reals against 3, and one fewer branch in the per-particle path.

### F2 — `m_mirrorScratch` was never declared

§4.2 step 3 deposits "into a scratch `EBAMRCellData`" and step 4 increments it into the real field;
§5.3's pseudocode names it `scratch`. §5.1's member list has four entries and none of them is it.

Now declared, `mutable`, one component — `EBAMRParticleMesh::depositWeight` asserts
`a_meshData[0]->nComp() == a_numComp` at `CD_EBAMRParticleMesh.H:639`.

While specifying it I checked the arithmetic it feeds, because a ghost mismatch between the scratch
and the caller's `a_phi` looked like a hazard. It is not one, for a reason worth writing down:
`EBCellFAB::plus(src, scale)` adds over `a_src.getRegion() & getRegion()` (Chombo
`EBCellFAB.cpp:591`) and asserts only on the component count. Both regions contain `dbl[din]`, so the
valid region is always complete; only the outer ghost ring is skipped, and after the deposit's
`exchange(..., EBAddOp())` (`CD_EBAMRParticleMesh.H:662`) the ghosts are meaningless anyway. That is
now documented in §4.2 so nobody "repairs" it.

### F3 — the surface data in a whole-domain cell field: cost quantified, container kept

> **Resolved: `EBAMRCellData` it is.** The lookup is on the per-particle hot path — one `FArrayBox`
> read at a known `IntVect`, no indirection — and that is why the container was chosen in the first
> place. The sparse alternative below is recorded so review 7 does not re-propose it; what remains
> legitimately open is the *component count*, not the container. The cost is accepted and is now
> reported in the memory report. Note also that 12 payload components is the floor for any correct
> format: revision 5's "9 Reals" was only 9 because it left the tangent frame unspecified, and the
> honest alternative (2×2 `S` plus an explicit `t1`) is 12 too.


§5.1 proposed `EBAMRCellData m_surfaceData` at 9 components, on the precedent of
`m_irregularCells`. That precedent carries **one** component.

At the shipped defaults — `AmrMesh.coarsest_domain 128 128 128`, `max_block_size 16`,
`num_ghost 2`, giving a 1.91× ghost overhead on 16³ boxes — one component is 33 MB per level. The 12
components §2.1 now needs are **~400 MB per level, per phase, per realm**; a two-phase dielectric
problem on two realms pays four times that. The four masks the realm already carries cost 131 MB per
level between them, so this triples the realm's whole mask budget to serve a set two cells thick.

The band is exactly what `BaseIVFAB` is for. `LevelData<BaseIVFAB<Real>>` over `grow(irregIVS, 2)`
with `BaseIVFactory<Real>(ebisl, bandIVS)` is the same three lines as the existing
`AmrMesh::allocate(EBAMRIVData&, ...)` at `CD_AmrMesh.cpp:361-373`, with a grown `IntVectSet` in
place of the bare `getIrregIVS`. For a sphere of `R = 32*dx` in that domain the band is ~7e4 cells,
so ~7 MB per level of payload — **well over an order of magnitude smaller**, even after
`BaseIVFAB`'s own per-patch index map — and it shrinks relative to the domain rather than tracking
it.
It is indexed by `VolIndex`, so multi-valued band cells become representable (F10), and covered
cells hold nothing at all.

The reflect test still needs to be O(1) per particle, which a `BaseIVFAB` does not give, so the band
flag stays a `LevelData<BaseFab<bool>>` — 1 byte per cell, ~4 MB per level, and the same container
`EBAMRParticleMesh` already uses for its halo and transition masks (`:511`). **No new type in either
case**, which matters under this repository's data-structure rules.

Separately: the build belongs in `regridOperators` beside `defineParticleMesh`
(`CD_PhaseRealm.cpp:327`), not in `defineMasks`. `defineMasks` sits in `regridBase` for a stated
reason — *"so the cell masks are available to load-balancing routines"* (`:182-183`) — and is
unconditional. Load balancing does not want a curvature fit, and an unconditional build makes every
realm pay the 5³ fit whether or not it deposits particles.

### F4 — surface data on a valid band cell at a level's grid edge

Review 5's Finding 1, confirmed. Its diagnosis is exactly right and I can add the cheap fix it did
not name.

The existing masks already avoid this failure by construction: `defineMasks` raises
`m_irregularCells` on `ebisbox.getIrregIVS(dbl[din])` — the **valid** box — and then exchanges
(`CD_PhaseRealm.cpp:646, 660`). A cut cell that no box on that level owns therefore never raises the
mask, even though `EBISBox` can see it. So **derive the band by growing the exchanged irregular
mask, never by asking `EBISBox::isIrregular` in the ghost region**, add the `ok` flag of §2.1, and
count the band cells that found no delivered cut cell.

The count matters more than the fallback: a nonzero count means the grids do not cover the EB band,
which is O9 seen from the other end. **O9 is now accepted rather than open** — see below — so the
`ok` flag's counter is a diagnostic, not a gate.

**Two corrections I owe O9.** First, I wrote that a demoted image is "overwritten by
`conservativeAverage`". It is not — `addInvalidCoarseToFine` interpolates coarse-grid deposition
clouds onto the fine level from inside the deposit, which is the case it exists for
(`CD_EBCoarseFineParticleMesh.H`, class doc item 3). No mass is lost; what differs is the
*bandwidth* of the cloud that delivers the correction. That is the same trade §4.2 already accepts
and documents for a promoted image, and it is the ordinary price of AMR — so O9 is accepted and
documented, not held open.

Second, on the knob: `buffer_size` is read into `m_bufferSizeBR` (`CD_AmrMesh.cpp:2704`) and reaches
exactly one place, the `BRMeshRefine` constructor in the `BergerRigoutsous` branch (`:1216-1221`).
The `Tiled` branch (`:1228`) constructs `TiledMeshRefine` with no buffer at all, and shipped inputs
split 70 `br` / 69 `tiled` — so for half the tree the parameter does nothing. Under `tiled` the fine
region is tile-aligned at `min_block_size` (16 by default), which is usually generous but can put a
tile boundary one cell from the EB with no option to widen it. The requirement therefore has to be
met by the tags, which makes §5.4's `ok` counter the only reliable detector.

### F5 — §7's cavity numbers are wrong in both directions

Review 5's Finding 3 identified the cause: `mirror_cavity.py` selected reflecting particles with
`d <= 1.5*dx*sum|n_i|`, the `3*s_max` distance band §9 retired in revision 3. It then *predicted*
the correction — `J` up to 3923 at `R = 3.5`, singular at `R = 3` under TSC — from the planar reach.

I re-ran it under the plan's real criterion instead of predicting, and **the prediction is wrong in
the more interesting direction**. The criterion is self-limiting in a cavity: the image of a deep
particle lands far inside the solid, where its cloud touches no `kappa > 0` cell, so it is never
built. Measured (`mirror_cavity.py`, rewritten):

| cavity `R/dx` | CIC reflect | CIC max `d` | CIC `Jmax` | TSC reflect | TSC max `d` | TSC `Jmax` |
|---|---|---|---|---|---|---|
| 8.0 | 42.1% | 2.289 | 3.2 | 57.9% | 3.151 | 5.3 |
| 6.0 | 47.8% | 2.328 | 5.1 | 64.9% | 3.192 | 10.7 |
| 4.0 | 65.0% | 2.133 | 10.8 | 82.8% | 2.959 | 44.7 |
| 3.5 | 75.3% | 2.234 | 20.5 | 90.6% | 2.973 | 151.0 |
| 3.0 | 73.9% | 2.059 | 28.9 | 90.5% | 2.739 | 483.8 |

So `J` reaches 28.9 at `R = 3*dx` under CIC, not 179, and 483.8 under TSC, and **it is never
singular** at any radius the harness can construct. The worst binned deviation is unchanged at
2.16%, at `R = 3.5`, so the conclusion — a NaN risk, not an accuracy risk — survives intact.

Two consequences. §3.2's guard has to be sized so that it does **not** fire at `J ~ 500`, i.e. not
on a denominator above ~1e-3; §3.2 previously said "clamp away from zero" with no number, which is
an invitation to pick 1e-2 and silently discard real images in a geometry this tree ships. And §7's
stated mechanism — *"the high-`J` images land where every stencil cell is covered and are dropped"* —
is now also the harness's mechanism, which review 5 correctly observed it was not.

### F6 — the conversion count, and where the `mirror` + NGP check goes

Two independent errors in §5.5.

**The count.** "93 C++ occurrences to convert" is a count of the *name*, not of the conversion.
`EBParticleMesh` and `EBAMRParticleMesh` use the same parameter name `a_forceIrregNGP` on the
**interpolation** path, which §5.5 explicitly excludes two paragraphs earlier ("Leave
`irr_ngp_interp` alone"). Attributing each occurrence to the function that owns it:

| File | total | interpolation | deposition |
|---|---|---|---|
| `CD_EBParticleMesh.H` | 49 | 22 | 27 |
| `CD_EBAMRParticleMesh.H` | 18 | 6 | 12 |
| `CD_AmrMesh.H` | 10 | 4 | 6 |
| `CD_AmrMeshImplem.H` | 10 | 4 | 6 |
| solvers | 6 | 0 | 6 |
| | 93 | **36** | **57** |

A mechanical conversion of all 93 crosses the scope boundary the section itself draws.

The plan should also not thread the enum to the leaf. `EBParticleMesh`'s kernels answer one
question — NGP in a cut cell, yes or no — and `Mirror`, `Redistribute` and `RedistributeBlended` all
answer *no* there, since the mirror is a separate pass and redistribution happens after the deposit.
Stopping the enum at `AmrMesh`/`EBAMRParticleMesh` and leaving the leaf a bool is **30 sites**
instead of 57, and keeps a single-purpose contract on the kernel.

**Where the check goes.** §5.5 says `Mirror` with an NGP base kernel "must still be rejected where
the kernel is chosen … because `CD_ItoKMCStepperImplem.H:4041` passes `DepositionType::NGP` as a
literal at the call site, so `parseOptions` structurally cannot see it."

That reads the call site wrong. `:4041` passes `DepositionType::NGP` **and** `a_forceIrregNGP =
false` (`:4047`), which maps to `IrregularDeposition::Native` — it can never be `Mirror`. The other
hardcoded site, `CD_ItoKMCGodunovStepperImplem.H:1338`, passes `true` (`:1344`), which maps to
`NGP`. Both are behaviour-preserving under a mechanical conversion and neither can produce the
collision. The combination that *can* arise is `deposition = ngp` with `irregular_deposition =
mirror`, and **both are options of the same solver**, parsed thirty lines apart in
`ItoSolver::parseDeposition` (`CD_ItoSolver.cpp:263` and `:320`). The check belongs there, where the
error message can name both keys.

The same misreading affects §6's `phase::solid` item: the mechanical conversion of `:1344`'s literal
`true` is `IrregularDeposition::NGP`, which preserves today's behaviour exactly. PR B does not have
to decide whether `phase::solid` deposits mirror; it has to *resist* deciding.

### F7 — the covered-cell reset, withdrawn

I raised this as a correctness gap: the reset cleans the valid box only, so covered cells stay dirty
as ghosts of the neighbouring patch, and `03cf59690` removed the ghost fill that used to erase them.

**Withdrawn — the reset itself came out instead.** The deposit's contract is that the field is
correct where it is read, and §4.1 reads it only at `kappa > 0`. Every image is inside the solid by
construction, so images landing in covered cells is the expected outcome, not something to clean up;
and ghost filling, coarsening and covered-cell presentation are the consumer's, which is the same
principle `03cf59690` established. A reset here would have been a *new* responsibility, not the
restoration of an old one — nothing zeroes covered cells after a particle deposit today either.

Four facts make the deferral safe rather than merely convenient, and they are now in §4.2:
`ItoSolver::coarsenAndFillGhosts` is `conservativeAverage` (`kappa`-weighted, so a covered fine cell
contributes exactly zero); plot output already resets covered cells at `CD_ItoSolver.cpp:1813`; a
covered cell has no `VolIndex` so the field solve never stencils it; and the pre-existing state is
that both valid and ghost covered cells hold spilled mass anyway.

What survives is the consumer-side question — which readers treat a covered cell as data — and that
is **O2**, which already names `DataOps::filterSmooth` and the `arithmeticAverage` sites. The
mechanism I documented (reversed copier at `CD_EBAMRParticleMesh.H:591-592`, so ghosts are the
exchange's source and are never written back) is worth keeping in mind when O2 is worked, but it is
not a defect in this plan.

### F8 — §2.4's ladder measured a nearest-cut rule §5.4 does not propose

§2.4's 1.2–2.9% — the number that justifies the quadratic patch over the plane, the own-foot plane
and the per-particle implicit function — is measured with `near = argmin(|x_band - x_cut|^2)` over
**all** cut cells (`mirror_source.py:48-53`). §5.4 step 3 proposes "two 1-cell propagation sweeps
with an exchange between them", which computes a *Chebyshev*-nearest assignment with an arbitrary
tie-break.

On a sphere those coincide, which is why no harness here can tell them apart. They differ exactly
where more than one piece of surface is within two cells — the O3 geometry. §5.4 now specifies
Euclidean-nearest over the radius-2 neighbourhood (125 comparisons per band cell, once per regrid)
with a fixed lexicographic tie-break, which is both faithful to the evidence and cheap.

### F9 — the reflect fractions are harness-box geometry, not a cost estimate

§7 quoted 62.7–98.0% and concluded "a vessel wall is not exotic". Every harness here puts a 16³
domain tightly around its geometry, so the band is a large share of the fluid volume *by
construction* — for a concave sphere of `R = 8*dx` with a reach of ~2.3, the shell is
`1 - (5.7/8)^3 = 64%` of the volume, and the harness duly reports 42–58%. In a real run the band is
two cells of a domain hundreds of cells across and the reflect fraction goes with
surface-to-volume.

What does **not** shrink with the band is the fixed per-deposit work: `remap()` is a collective, and
the increment and the covered-cell reset are whole-domain. Those are paid on every deposit of every
species regardless of how thin the band is, and in this tree `remap()`-class costs have dominated
before. §7 now says that, and §5.3 gains an early-out on levels with no cut cells.

### F10 — two unlisted open items

**O10, multi-valued cells on the band side — closed by invariant.** Multiply-cut cells are always
refined away in this project's workflows, so a particle's cell is single-valued and the reflect path
never has to choose between two VoFs — which it could not do anyway, since a position alone does not
determine a VoF. The invariant is a workflow invariant: I looked and **nothing in the library
enforces it** (`AmrMesh::getMultiCutVofIterator` exists precisely because the CDR solvers handle
these cells defensively, at a dozen sites). So the assertion *is* the enforcement:
`CH_assert(!ebisbox.isMultiValued(particleIV))` on the reflect path.

**O11, deposits that never reach the funnel.** `ItoSolver::depositWeightNGP` and
`depositGatheredNGP` (`CD_ItoSolverImplem.H:84`) build an `EBParticleMesh` per patch themselves
(`:115`) and never call `AmrMesh::depositWeight`, so no mirror pass can apply. They are plot paths
(`CD_ItoSolver.cpp:1705-1755`) and deliberately NGP, so this is correct — but it means the plotted
particle density in a cut cell will not match the `phi` the field solver sees, by exactly the factor
the mirror introduces. That is a diagnostic people will read as a bug.

---

## 2. Things I checked that held

Recorded so review 7 does not re-spend the time.

- **The CIC/TSC kernel facts.** `L = 0.5*a_particleWidth` top-hat at `:732`, `L = a_particleWidth`
  triangular at `:794`, `cicWidth = 1*widthScale` / `tscWidth = 2*widthScale` at `:474-475`. The TSC
  clamp is `alpha = max(a, -0.5L)`, `beta = min(b, +0.5L)`, so the support is `[-1, +1]` in cell
  units and the half-width is 1, as §1.2 says.
- **`cloudBox` (`:677-688`)** anchors the iterated box on the cloud's lower and upper edge cells, so
  it returns exactly the 2 (CIC) / 3 (TSC) cells per direction §1.2 tabulates — the plan's reach is
  the code's reach, not an over-scan.
- **`CH_assert(m_region.contains(particleIV))` at `:713`** is a second, harder reason images must be
  remapped: a patch cannot deposit a particle whose *cell* is outside its region at any ghost width.
  §1.2 now cites it.
- **`Halo` really is not the target strategy.** All **60** shipped `deposition_cf` settings are
  `transition`, and the only two hardcoded `CoarseFineDeposition::Halo` call sites are exactly the
  two that map to `Native`/`NGP` and therefore never run a mirror pass. §4.2's disposition of the
  `widthScale` mismatch is now measured rather than asserted.
- **`ParticleSoA::append(pos, weight)`** leaves id/rank at `s_invalidID`/`-1` (`:949-958`) and
  nothing in `remap()` reads them, so images need no identity — worth stating, because the obvious
  defensive move (`resetParticleIDs`) is wasted work.
- **The ghost requirement is already met everywhere.** The minimum of both shipped columns is 2
  (39×`eb_ghost=2`, 71×`num_ghost=2`), so §1.3's runtime check is a guard against a future input
  file, not a migration.
- **All five of review 5's nits**, including the eighteen `mutable` members at
  `:605, 610, 615, 622, 630, 635, 640, 645, 650, 655, 660, 665, 670, 675, 680, 685, 690, 695`.
- **`reach_cells.py`** reproduces to the printed digit: CIC 91.596% / TSC 87.163% at Chebyshev 1,
  100.000% at 2 for all ten cases both kernels, max `d` 2.570 / 3.390.

---

## 2b. Second pass — three of my own conclusions, withdrawn or downgraded

All three came from the same mistake: taking responsibility, on the deposit's behalf, for what it
hands to its consumer. The deposit's contract is that the field is correct where it is read, and
§4.1 reads it at `kappa > 0`. Everything past that — ghost cells, coarsening, covered-cell
presentation, what a downstream filter believes a cut cell means — is the consumer's, which is the
division `03cf59690` already drew for synchronization.

- **F7, the covered-cell reset — withdrawn**, and the reset came out of the plan with it. Detail
  above.
- **O9, the coarse-fine correction — accepted.** I also had a fact wrong: a demoted image is not
  overwritten by `conservativeAverage`; `addInvalidCoarseToFine` routes coarse-grid cloud mass to the
  fine level from inside the deposit. No mass is lost, the residual is bandwidth, and that is the
  trade §4.2 already accepts for a promoted image.
- **O2, the audit of `phi` readers — downgraded from blocker to contract.** What this PR owes is a
  clear statement that under `mirror` a cut cell holds `n`, not `kappa*n`, in the enumerator's
  Doxygen and in `Ito.rst`. The readers themselves belong to whoever opts in.

**What the second pass *added*, and it is worth more than what it removed.** Asked to evaluate why
two containers are needed, I found the reason the plan gave — "`depositWeight` resets the mesh
field" — is the weaker of two, and the stronger one had gone unstated through six revisions:
`DataOps::setValue` zeros ghosts too, and the level exchange is a **reversed** copier with an add op
(`CD_EBAMRParticleMesh.H:591-592`, used at `:662`), so the ghost region is the exchange's source and
is never written back. A "deposit again without resetting" would therefore fold the first pass's
ghost content into the valid region a second time — silent double counting proportional to each
patch's ghost share. So the scratch field is not ergonomics that a `bool a_reset` could remove: the
reset and the exchange are one mechanism. That is now in §4.2, next to the two alternatives I
evaluated and rejected (merging into a copy of the source container; depositing images first).

**And one measurement that made the migration safe.** The enum collapses two independent booleans
into one selector, so it is only behaviour-preserving if no shipped input uses a combination the
enum cannot express. Cross-tabulating: 22 files `false/true/false`, 5 `false/false/false`, 2
`true/false/false`, and **zero** with `irr_ngp_deposition = true` together with
`redistribute = true`. Every shipped configuration maps one-to-one, so the rename can land as a pure
signature change with every regression bit-for-bit unchanged.

## 3. What I did not do

- **I did not re-run the long harnesses.** `mirror_source.py`, `mirror_sphere_ext.py radii` and the
  curvature sweeps were re-run in full by review 5 and reproduced exactly; I sampled rather than
  repeated. If revision 7 wants a fully independent reproduction, those are the three to re-run.
- **The anisotropic composition is still not measured end to end.** `mirror_frame.py` shows the
  frame-free form is exact and shows what a mismatch costs, but no harness here yet takes a *fitted*
  anisotropic `S` from `mirror_discrete_curvature.py` and pushes it through a deposition. §6.1 now
  mandates the torus for exactly this reason, and mandates asserting that the sliver bin really
  contains cells with `c_1 != c_2` — a torus resolved too coarsely degenerates to the sphere case
  without saying so.
- **Nothing here touches what the real geometry generator delivers.** Unchanged from revision 5:
  the first number appears in phase 3.

---

## 4. Reproducing this review

From `Prototypes/MirrorDeposition/`:

```
python3 mirror_frame.py                         # F1              (new in r6, ~2 s)
python3 mirror_cavity.py                        # F5              (rewritten in r6, ~3 min)
python3 reach_cells.py                          # the r5 spot-checks
```

Counts were taken with per-file `awk` attribution to the owning function rather than raw `grep -c`,
because the deposition and interpolation paths share the parameter name — a raw count is what
produced F6. Code citations were re-resolved against `b8721d71` with `sed -n '<line>p'`.

**This file is throwaway and must be deleted before the branch merges, along with the rest of this
directory.**
