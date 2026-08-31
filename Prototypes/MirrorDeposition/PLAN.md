# Mirrored cut-cell deposition — implementation plan

**Issue** chombo-discharge#29 (parts 1 and 3) · **Branch** `mirror_deposition` · **Depends on** PR #700
· **Revision** 6, after five adversarial reviews.

> **Rebased onto `main` at `f44e5968` on 2026-08-31.** #700 is now *in* `main` (as `92448b90`), so this
> branch no longer carries its own copy of it. The design below is unchanged — no mechanism it depends
> on was touched — but a handful of line citations and shipped-input counts moved, and they have been
> re-resolved in place. `IMPLEMENTATION.md` §0.4 tabulates exactly what moved and why, and records one
> pre-existing constraint the plan had missed (`num_ghost < min_block_size`, enforced in optimized
> builds). Read §0.4 before re-checking any number here.

Particle deposition divides by `dx^D` and never by `kappa`, so a cut cell holds `kappa*n` where the
CDR solvers hold `n`. Space charge and the semi-implicit conductivity are assembled from both, so a
quasi-neutral plasma against a wall produces a full-strength spurious charge layer. This plan
removes that by depositing each particle's cloud **and** the cloud of its image reflected across the
embedded boundary, so the field is a true density up to the wall.

Reflection was chosen over hybrid-plus-redistribution (option A in the issue) because the volume
fraction `kappa` never appears in a denominator, so small `kappa` is not a singular limit. The cost
is that conservation becomes quadrature-accurate rather than exact.

---

## 0. Instructions for the next review agent

Read this section first.

**The document you are reviewing is this file.**

**What five prior reviews have already caught.** Each was a confident, plausible, wrong claim that
survived a reading of the plan and died to a measurement or to the code. Assume more remain.

| Review | What it refuted |
|---|---|
| 1 | The band formula `dx(1 + sum\|n_i\|/2)`; the "aggregate bias bounds the per-case bias" reading |
| 2 | The reflection source attribution; `sum(kappa phi dV)/sum(m)` as a primary metric; one table that did not reproduce at all |
| 3 | The band as a distance; "does not degrade at small kappa"; `bndryCentroid` units; "positivity by construction" |
| 4 | The whole post-deposit redistribution step, which the plan never mentioned; `eb_ghost` as the governing ghost count; the CIC retention column, measured with the wrong kernel; the `kappa < 0.25` threshold as kernel-independent; the level-preserving remap, which was both unnecessary and harmful; "`EBAMRParticleMesh` needs a realm" |
| 5 | Nothing quantitative — every table reproduced exactly. It refuted the *scope* of three claims: the tangent frame as an unstated accident rather than a specification; the cavity table, computed with a band §9 had already retired; and "extend to every band cell from the nearest cut cell" as a total statement |
| 6 | The tangent frame as a stored quantity at all; §7's cavity numbers, wrong in **both** directions; the 93-occurrence conversion count, 36 of which are the interpolation parameter the plan puts out of scope; "`parseOptions` structurally cannot see" `mirror` + NGP; the reflect fractions read as a cost estimate; and a mesh field the pseudocode uses that §5 never declares |
| 6, second pass | Three of revision 6's own conclusions, all from the same mistake — taking responsibility for what the deposit hands to its consumer. The covered-cell reset (withdrawn: images belong in the solid and the field is only read at `kappa > 0`); the coarse-fine correction loss (accepted: `addInvalidCoarseToFine` conserves it, the residual is bandwidth); and the audit of `phi` readers as a blocker (it is a contract to state, and the consumers' work). Revision 6 also claimed a demoted image is overwritten by `conservativeAverage`, and named `AmrMesh.buffer_size` as a knob that only the `br` grid generator reads |

**How to review this.** In descending order of yield:

1. **Re-run every table.** Revision 2 shipped a table that reproduced in one row of seven; revision 4
   shipped a table computed with a kernel the code does not use. Every number in this file is
   reproducible from this directory — see §8. Two specific traps:
   - Check the *stencil label* on curvature output: `mirror_discrete_curvature.py fit` prints a
     `stencil radius 1` (3³) block and a `stencil radius 2` (5³) block, and revision 3 quoted columns
     from both in one table.
   - Check that a harness's kernel is the code's kernel. `reach_cells.py` and
     `mirror_band_kernels.py` applied the TSC weight formula at `L = 1` and called it CIC.
     chombo-discharge's CIC is a *top-hat* cloud (`CD_EBParticleMesh.H:732`); only TSC is the
     triangular-cloud form (`:794`). Both are fixed now; the class of error is not.
2. **Check the code citations.** They are given as `file:line` throughout and were re-resolved
   against `1d86431b` in review 4. The branch moves; re-resolve them rather than trusting them.
3. **Follow the deposit past the deposit.** Review 4's largest finding was not in this plan at all —
   it was the `redistributeAMR` call one line below the `AmrMesh::depositWeight` this plan wires.
   Ask what else reads `phi` on the assumption that a cut cell holds `kappa*n`: §7 O2 lists what is
   known, and that list is not closed.
4. **Attack the derivations in §2 and §3.** The reflection point and Jacobian are where a sign or a
   frame error is silent rather than loud. Review 4 verified both against a convex sphere and a
   concave cavity and found them correct; that is one check, not a proof.
5. **Do not re-derive widths along the normal.** §1 states every width in cells. Two earlier reviews
   produced wrong ghost counts (4–5, then 3) by computing a reach in `dx` along the normal and
   re-projecting onto index axes. Work in index space.
6. **Price the containers, and check that every field the pseudocode uses is declared.** Revision 6
   found that §4.2's scratch `EBAMRCellData` appeared in no member list anywhere in §5, and that the
   surface data as a whole-domain cell field costs ~426 MB per level per phase at the shipped
   defaults — 3.2× the four masks the realm already carries. Count the fields, multiply by
   `coarsest_domain` and the ghost overhead, and cross-check §4.2 and §5.3 against §5.1 line by line.
7. **Check that a harness's *selection* is the plan's criterion, not just that its kernel is.**
   Review 5 caught `mirror_cavity.py` selecting reflecting particles with a band §9 had retired;
   revision 6 re-ran it under the real criterion and every number in §7's cavity row moved.

**Known-open items, listed so you do not spend the review rediscovering them:** §7.

**Scope limit on all the evidence.** Every harness in this directory is Cartesian, single-level,
single-patch, with analytic geometry. Cut-cell normals are reconstructed from *exact* sub-sampled
face area fractions. A real `GeometryShop`/`ScanShop` computes moments to finite order and has its
own pathologies. Nothing here measures what the real geometry generator delivers; phase 3 is where
that number first appears.

---

## 1. Required widths

All widths are stated in **cells**, in index space. This is deliberate: a width stated as a distance
along the normal has to be re-projected onto index axes to size a mask or a ghost region, and that
round trip is where three wrong answers came from.

### 1.1 What a deposition cloud reaches

Both chombo-discharge kernels are **cloud-overlap integrals**: the particle carries a cloud of side
`L*dx` centred at `x`, and a cell's weight is the overlap of that cloud with the cell volume.
`L = 1` is CIC, `L = 2` is TSC. The cells reached are exactly those overlapping
`[x - L*dx/2, x + L*dx/2]`.

The two clouds are not the same shape, and it matters for any mass-weighted statement. CIC is a
**top-hat** cloud — `weight = max(0, min(b,L) - max(a,-L))` with `L = 0.5*cicWidth`
(`CD_EBParticleMesh.H:732`). TSC is a **triangular** cloud —
`weight = (beta-alpha) - (beta|beta| - alpha|alpha|)/L` with `L = tscWidth`
(`CD_EBParticleMesh.H:794`). Both are partitions of unity after their respective factors
(`:481-482`). Only the *support* is shared machinery, and only support statements transfer between
them.

With `probLo = 0`, cell `i` spanning `[i*dx, (i+1)*dx)`, the particle's own cell `c = floor(x/dx)`
and its fractional position `f = x/dx - c` in `[0,1)`:

| Kernel | Cloud, in cell units | Cells reached, per direction | Cells in 3-D | Chebyshev reach |
|---|---|---|---|---|
| CIC, `L = 1` | `[c+f-1/2, c+f+1/2]` | `{c-1, c}` if `f < 1/2`; `{c, c+1}` if `f >= 1/2` | 2³ = 8 | **1** |
| TSC, `L = 2` | `[c+f-1, c+f+1]` | `{c-1, c, c+1}` (just `{c-1, c}` at `f = 0`) | 3³ = 27 | **1** |

(Both are also what `EBParticleMesh::cloudBox` (`:677-688`) returns: it anchors the iterated box on
the cloud's lower and upper edge cells, so the loop visits exactly the reached cells and no others.)

**Both kernels reach exactly Chebyshev distance 1 from the particle's own cell.** They differ in
*how many* of the three cells per direction they touch — two versus three — not in how far. CIC
takes its own cell plus the single neighbour the particle leans toward, so its offsets are `{-1,0}`
or `{0,+1}` and never both.

Reading the harness: `mirror_band_kernels.py:support_offsets()` returns `ceil(L/2)+1 = 2` for both
kernels. That is a safe over-scan radius for the loop, **not** the reach.

`NGP` is the third `DepositionType` (`CD_DepositionType.H:25-27`) and reaches 0. It is excluded from
`mirror` — see §5.5.

### 1.2 Band mask: Chebyshev 2, for both kernels

> Reflect only particles in cells within **Chebyshev distance 2** of a cut cell. Within the mask,
> reflect a particle if and only if its image's cloud overlaps at least one cell with `kappa > 0`,
> and `d_p > 0`.

The exact criterion falls out of the deposit loop for free and costs nothing, because the overlap
weight is already zero at the support boundary. The mask is only an early-out.

> **Check the early-out under `DEBUG`.** Chebyshev 2 is *derived* only for axis-aligned normals
> (O6) and *measured* over ten tilted planar cases plus two sphere signs. That is good evidence and
> it is not a proof, and the failure — a reflection silently not taken — leaves no trace. Under
> `DEBUG`, evaluate the exact criterion on a band grown to Chebyshev 3 and count how many particles
> in the outer ring would have reflected. The count must be zero; if it is not, the mask is short
> and the number says by how much. It costs nothing under `OPT=HIGH`.

Measured against the real overlap-integral kernels — top-hat for CIC, triangular for TSC — share of
mirrored mass retained, counting only the part of each image's cloud landing in a `kappa > 0` cell
(`reach_cells.py`, worst case over ten planar cases):

| Band mask | CIC worst | TSC worst | Verdict |
|---|---|---|---|
| Chebyshev <= 2 | **100.000%** | **100.000%** | exact, both kernels, all ten cases |
| Chebyshev <= 1 | 91.6% | 87.2% | loses 8.4% / 12.8% on `axis, s = 0.45` |
| Chebyshev <= 0 | 5.1% | 5.1% | cut cells only — hopeless |

> **Revision 4 corrected the CIC column.** Revisions 3–4 read 95.5% at Chebyshev 1 because
> `reach_cells.py` applied the TSC triangular-cloud formula at `L = 1`. With the code's top-hat CIC
> the figure is 91.6%. Every `max d` and `max cheb` entry was unaffected — those depend only on the
> support, which the two forms share — so the Chebyshev-2 result never moved, and the case for 2 got
> stronger.

Maximum Chebyshev distance from a reflected particle's cell to the nearest cut cell is **2 for both
kernels**, even though TSC's deepest reflecting particle sits at `d = 3.39*dx` against CIC's 2.57.
The mask does not depend on the kernel.

**Why 1 is not enough, and exactly when 2 is needed — this part *is* kernel-dependent.** The
Chebyshev-1 reach belongs to the *image's* cell, not the particle's, and the reflection moves the
particle across the boundary by `2d`. Take an axis-aligned EB at `x = c + t` inside cut cell `c`, so
`kappa = 1 - t`. The image must overlap cell `c`, which for a cloud of half-width `h` means
`x_img > c - h`, so

```
d_max        = t + h
deepest x_p  = c + 2t + h
```

and the particle leaves cell `c+1` when `2t + h > 2`. The half-width is the kernel's:

| Kernel | cloud half-width `h` | Chebyshev 2 needed when | i.e. |
|---|---|---|---|
| CIC | 1/2 | `t > 0.75` | **`kappa < 0.25`** |
| TSC | 1 | `t > 0.5` | **`kappa < 0.5`** |

**Chebyshev-1 suffices only if every cut cell has `kappa >= 0.25` (CIC) or `>= 0.5` (TSC)** —
precisely the regime the mirror does not care about. The sliver cells it exists to fix are the ones
that need 2. In both cases the deepest particle satisfies `x_p < c + 3`, so Chebyshev 2 is the bound
for both kernels; TSC merely reaches it over twice the `kappa` range.

The rule predicts every axis-aligned case in the harness, **once the kernel is taken into account**:

| Case | `t` | `kappa` | CIC predicted / measured | TSC predicted / measured |
|---|---|---|---|---|
| `axis-aligned x` | 0.37 | 0.63 | 1 / 1 | 1 / 1 |
| `near lo edge` | 0.40 | 0.60 | 1 / 1 | 1 / 1 |
| `axis, s=0.12` | 0.62 | 0.38 | 1 / 1 | **2 / 2** |
| `axis, s=0.45` | 0.95 | 0.05 | 2 / 2 | 2 / 2 |

Revision 4 listed `axis, s=0.12` as Chebyshev 1 for both kernels. It is 2 under TSC, and the
Chebyshev-1 retention there is 99.617%.

Worked example at `kappa = 0.05` (EB at 8.95, cut cell 8), CIC:

| `x_p` | `d` | particle cell | Cheb(p,8) | `x_img` | image cell | Cheb(img,8) | w into cell 8 |
|---|---|---|---|---|---|---|---|
| 9.20 | 0.25 | 9 | 1 | 8.70 | 8 | 0 | 0.80 |
| 9.90 | 0.95 | 9 | 1 | 8.00 | 8 | 0 | 0.50 |
| 10.00 | 1.05 | 10 | **2** | 7.90 | 7 | 1 | 0.40 |
| 10.20 | 1.25 | 10 | **2** | 7.70 | 7 | 1 | 0.20 |
| 10.40 | 1.45 | 10 | — | 7.50 | 7 | 1 | 0.00 |

**The image is never two cells away — the particle is.** The 1-D argument: `x_img >= c - h` and
`x_img < a < c+1`, so the image's cell is `c-1` or `c` for either kernel. Revision 5 measured it in
3-D rather than leaving it derived on an axis-aligned normal (`img_reach.py`):

| | max Chebyshev, image cell → nearest cut cell | max Chebyshev, image cell → source cell |
|---|---|---|
| CIC, worst of ten cases | **1** | 4 |
| TSC, worst of ten cases | **1** | 5 |

The 2 arises because the mask is indexed by the *particle's* cell — that is where the reflect
decision is taken and where `(x_c, n_c, S_c)` is read — and reflection is symmetric in *distance*
but not in *cells*. With the EB near the top face of cell `c`, walking up leaves that cell after
`1-t` and reaches `c+2` quickly, while walking down stays inside it for `t`.

The second column of that table is the reason images are remapped rather than deposited by the
source's own patch: a design in which each patch deposits its own images needs `4 + 1 = 5` ghost
cells for CIC and `5 + 1 = 6` for TSC. It is also not merely a memory argument —
`depositParticleCIC` opens with `CH_assert(m_region.contains(particleIV))`
(`CD_EBParticleMesh.H:713`), so a patch cannot deposit a particle whose *cell* lies outside its
region at all, whatever the ghost width. See §4.2.

> **Open:** the `kappa` thresholds are derived for axis-aligned normals, where `kappa` and `t` are
> the same parameter. For a tilted normal they decouple. The Chebyshev-2 *result* is measured for
> tilted normals (all ten cases, both kernels); the tilted threshold is not derived. (O6)

### 1.3 Ghost cells: `eb_ghost >= 2` **and** `num_ghost >= 2`

These are two different parameters governing two different things, and revision 4 conflated them.

- **`eb_ghost`** (`m_numEbGhostsCells`) sizes the `EBLevelGrid` — `new EBLevelGrid(m_grids[lvl],
  m_domains[lvl], m_numEbGhostsCells, &(*m_ebis))` at `CD_PhaseRealm.cpp:411`. It is how far out
  `EBISBox` geometric data (`normal`, `bndryCentroid`, `volFrac`) is valid.
- **`num_ghost`** (`m_numGhostCells`) sizes the *data*. `PhaseRealm::defineMasks(a_lmin, a_numGhost)`
  is called with `m_numGhostCells`, not `m_numEbGhostsCells` (`CD_PhaseRealm.cpp:189`), and that is
  the routine §5.1 builds the surface data in.

| Requirement | Parameter | Needs | Forced or chosen? |
|---|---|---|---|
| Curvature fit stencil (5³) reads `EBISBox` 2 cells out | `eb_ghost` | 2 | **Chosen** — 5³ was picked because 7³ does not fit in `eb_ghost = 2` |
| Band mask sees cut cells 2 out, to mark valid cells within Chebyshev 2 | `eb_ghost` | 2 | **Forced** |
| Surface data `(x_c, n_c, S_c)` exchanged into ghosts, then extended | `num_ghost` | 2 | **Forced** — a band cell's nearest cut cell is at most 2 away |

**Keep them separate in the code** — do not derive one from another. Deposition itself needs only 1,
because each image is remapped to its own owning patch (§4.2).

Move to a 7³ fit for noise headroom (3.7% against 5.0% at 10° of normal error) and `eb_ghost`
becomes 3 while `num_ghost` and the band stay 2.

The `num_ghost = 2` holds **only** under the ordering in §5.4 — exchange the fitted cut-cell values
into ghosts, then extend. Extend first and exchange once afterwards and a patch must fit `S` for cut
cells 2 cells outside its valid region, each of which needs its own radius-2 neighbourhood, reaching
4 cells out.

Neither is a safe default. `CD_AmrMesh.cpp:2772` is `pp.get("eb_ghost", ...)` validated only as
`>= 0` (`:2774`), and `:2789` is `pp.get("num_ghost", ...)` validated only as `>= 0` (`:2792`).
`CD_AmrMesh.options` ships 2 for both; shipped inputs split 40×`eb_ghost=2` / 100×`=4` and
72×`num_ghost=2` / 68×`=3`. Note also the *upper* bound already enforced in optimized builds:
`AmrMesh::sanityCheck` requires `num_ghost < min_block_size` (`:2966`), so the band radius is
bounded on both sides — see `IMPLEMENTATION.md` §0.4. **Require both `>= 2` at runtime when `mirror` is selected.**

Note what that split says: every shipped input already satisfies both requirements — the minimum in
each column is 2. The runtime check is a guard against a future input file, not a migration, and it
costs nothing today. It is also the only one of the two that a user can get wrong silently, since
`num_ghost = 1` would produce a wrong extension rather than an out-of-bounds access.

---

## 2. The reflection point

### 2.1 Stored state

Three quantities per band cell, taken from the nearest cut cell, plus a validity flag; all
computed once at regrid:

```
x_c    boundary centroid, as a PHYSICAL position    3 Reals  (2 in 2-D)
n_c    normal(vof), pointing INTO the fluid         3 Reals  (2 in 2-D)
S_c    fitted shape operator, in the WORLD frame    6 Reals  (3 in 2-D; symmetric, S_c n_c = 0)
ok     may this cell reflect at all                 1 flag   (see 5.4)
```

> **`S_c` is stored frame-free, and this is not a presentation choice.** Revision 5 stored it as the
> 3 independent entries of a 2×2 symmetric matrix "in the tangent frame", with §3.3 fitting it in
> "any orthonormal basis perpendicular to `n_i`" and §2.2 reflecting in "`t1, t2` spanning the plane
> perpendicular to `n_c`". Those have to be the *same* basis, and nothing said so. The failure mode
> is the worst available: `2H` and `K` are frame-invariant, so a mismatch leaves `J` correct, leaves
> §6.1's mass check correct, and moves only the reflection point.
>
> `mirror_frame.py`, cylinder `R = 6*dx` (principal curvatures `1/R` and `0` — every other harness
> here uses an *isotropic* `S`, which is frame-invariant and therefore cannot see this):
>
> | | median | p95 | max |
> |---|---|---|---|
> | image displacement, mismatched frame | **0.397 dx** | 1.36 dx | **2.24 dx** |
> | change in `tr S` under the same rotation | 2.8e-17 | 5.6e-17 | 8.3e-17 |
> | change in `det S` under the same rotation | 7.3e-20 | 1.4e-18 | 3.5e-18 |
>
> Storing `S_c = P S P^T` with `P = [t1 t2]` — the same operator expressed as a symmetric 3×3 in
> world coordinates — removes the question. The tangent frame becomes local to the fit (§3.3) and is
> discarded there; nothing downstream ever constructs one. Verified exact, not approximate: the
> frame-free reflection reproduces the framed one to 1.8e-15 in the image position over 4000 draws.
> It costs 6 Reals against 3, removes a branchy frame construction from the per-particle path, and
> it is what makes §3.1's 2-D special case disappear.

> **`bndryCentroid` is cell-relative — use the existing helper.** `EBISBox::bndryCentroid(vof)`
> returns coordinates in `[-0.5, 0.5]` *relative to the cell*, not a position. The tree already has
> the conversion, and it is not the arithmetic you should re-type:
>
> ```cpp
> x_c = probLo + Location::position(Location::Cell::Boundary, vof, ebisbox, dx);
> ```
>
> `Location::position` is `CD_LocationImplem.H:37-51` and returns
> `(gridIndex + 0.5 + bndryCentroid) * dx`; `CD_AmrMeshImplem.H:584` already calls it this way. The
> range must include `:49`, `ret *= a_dx` — that line is the whole point of the citation, and
> revision 5 cited `:37-41`, which stops at the `Location::Cell::Boundary` switch case.
> (`CD_EBHelmholtzEBBCImplem.H:32` open-codes the same thing; do not add a third copy.)
>
> Revision 3 wrote the curvature fit as `dx_j = bndryCentroid(j) - x_i`, differencing two
> cell-relative values from *different* cells. That drops the inter-cell offset entirely: for cut
> cells one cell apart the true tangential baseline is `~dx` and the differenced raw centroids give
> `~0.1*dx`, so `S` comes out inflated by roughly 10×, with nothing in the output looking wrong.
> Store and use `x_c` as a physical position throughout, including in `w = p - x_c`.

The boundary centroid is the only centroid the scheme uses; `Location::Cell::Centroid` (the volume
centroid) appears nowhere in it.

### 2.2 Per particle

With `w = p - x_c`:

```
eta    = n_c . w

d      = eta + (1/2) w . (S_c w)                   one-step distance to the quadratic patch
nhat   = normalize( n_c + S_c w )                  its normal at the foot point
R(p)   = p - 2 d nhat
```

`S_c n_c = 0` by construction, so the normal component of `w` contributes nothing to either
expression and no tangential projection is needed: `w . (S_c w) = xi^T S xi` and `S_c w = P S xi`
identically. **There is no tangent frame at reflect time.**

One dependent load per particle and about twenty flops. **No implicit-function evaluation
anywhere**, at regrid or per particle.

Review 4 checked the signs independently against both curvature signs, with `S` in the fit's own
convention (`S = +dn̂`, §3.3): a convex sphere of radius `R` with fluid outside gives `S = +I/R` and
patch height `h(xi) = -|xi|²/2R`, a concave cavity gives `S = -I/R` and `h(xi) = +|xi|²/2R`, and
`d = eta - h` reproduces the expression above in both. `nhat = n_c + S xi` likewise. Lifting `S` to
the world frame is an identity, not a re-derivation, so those checks carry over unchanged
(`mirror_frame.py`, rows 1–2).

`d` is the **first iterate** of the signed distance, not the signed distance: `xi` is the tangential
projection of the particle, not of the foot point. The neglected term is `O(c^2 d^2 xi)`, and at
`R = 4*dx` with `d` out to 2.6`*dx`, `c*d ~ 0.65` is not a small parameter. The end-to-end error
there measures 2.9% (§2.4), so one step is enough at the radii tested — revision 3 called it
"signed distance to the quadratic patch", which invites someone to trust it on a rod cap without
re-measuring.

### 2.3 Sign convention — assert it

Every sign depends on `ebisbox.normal()` pointing **into the fluid**: `2H > 0` for a convex solid,
`d > 0` on the fluid side, `J < 1` convex. It holds in this tree —
`CD_AmrMeshImplem.H:583-589` accepts a particle as inside the fluid when
`ebNormal.dotProduct(pos - ebCentroid) >= -tol`, and `DataOps::computeMinValidBox`
(`CD_DataOps.cpp:3687`) uses the same convention.

Reversing it replaces `J` by `1/J`: at `R = 4*dx` and `d = 2*dx` that is 0.111 against 9.0, an ~80×
error in the image weight, with no crash and no NaN. **One assert, one comment.**

### 2.4 Why the quadratic patch

Worst deviation over `kappa` bins, exact Jacobian, 3M particles, margin 3 (`mirror_source.py`):

| Sphere | nearest-cut PLANE | nearest-cut QUADRATIC | own-foot plane | per-particle IF |
|---|---|---|---|---|
| convex R = 16 | 1.3% | **1.2%** | 1.4% | 1.1% |
| convex R = 12 | 3.4% | **2.1%** | 1.5% | 1.9% |
| convex R = 8 | 7.3% | **1.8%** | 0.4% | 1.6% |
| convex R = 6 | 9.1% | **2.3%** | 1.4% | 0.9% |
| convex R = 4 | 18.2% | **2.9%** | 3.1% | 1.0% |
| concave R = 8 | 6.7% | **1.3%** | 3.1% | 1.2% |
| concave R = 6 | 9.2% | **1.9%** | 5.4% | 1.4% |

The quadratic patch is <= 2.9% everywhere, matches the per-particle IF, beats the own-foot plane on
concave geometry, and needs no implicit function. Cut-cell normals and centroids are analytic here,
so this isolates the surface *model*; §3.3 prices the discrete inputs. Together, roughly 2–4%
end to end.

---

## 3. The Jacobian

### 3.1 Form

```
J   = (1 - 2H*d + K*d^2) / (1 + 2H*d + K*d^2)
    = prod(1 - c_i d) / prod(1 + c_i d)              over principal curvatures c_i

with, from the stored world-frame S_c:

2H  = tr S_c
K   = (1/2) [ (tr S_c)^2 - tr(S_c^2) ]
```

`J` is exactly 1 on a plane, where `S_c = 0`. Sphere check: `J = (R-d)^2/(R+d)^2`, the ratio of shell
areas, as it must be.

> **`K` is NOT `det S_c`.** `S_c` is a symmetric 3×3 of rank at most 2 — it annihilates `n_c` — so
> `det S_c` is identically zero and `K = det S_c` silently degrades the exact Jacobian to the
> linearized one, whose error §3.3 measures at 48–77%. The second-invariant formula above is the
> only correct reading, and it is the same expression in both dimensions.
>
> This is what the frame-free format buys beyond §2.1's determinism. Revision 5 carried a 2-D
> special case — "`S` is 1×1, set `K = 0` explicitly, because `det` of a 1×1 `S` is `c`, not 0" —
> that is now unnecessary: in 2-D `S_c = c * t t^T`, so `tr S_c = c`, `tr(S_c^2) = c^2` and the
> formula returns `K = 0` on its own. One expression, both dimensions, no override to forget.
> Verified over 2000 random normals and curvature pairs (`mirror_frame.py`): both invariants agree
> with the 2×2 `tr`/`det` to 2e-16, and `|S_c n_c| < 2e-17`.

**Symbols.** `kappa` is the volume fraction everywhere in this document; principal curvatures are
`c_1, c_2`. Revision 3 used `kappa` for both.

Without the Jacobian the mirror carries a curvature error of +117% at `R = 4*dx`, +40% at 8, +15.3%
at 16 and +4.2% even at 40 — with *exact* reflection, so no improvement to the surface model removes
the need for it. Reflection is measure-preserving only about a plane.

Those four figures are the **worst `kappa` bin**, not the sliver bin: the worst bin is the sliver at
`R = 4, 8, 16` but the `0.15–0.30` bin at `R = 40` (+0.0424 against the sliver's +0.0219). Every
other table in this file bins by `kappa` and quotes the sliver, so the selection rule is stated here
rather than left to be rediscovered.

### 3.2 Positivity is conditional, not structural

Revision 3 claimed the scheme is "positivity-preserving by construction". It is not. The
**numerator changes sign** once `d > 1/c_max` on any non-umbilic convex surface — a negative-weight
image, i.e. `phi < 0`. (The qualifier matters: on a sphere the numerator is `(1-cd)^2` and never goes
negative, so a sphere test cannot see this.) Positivity holds for `|d| < 1/c_max` and is enforced by
the guard, not by construction.

The **denominator** vanishes at `d = -1/c_i` on concave surfaces. These are the same event at the
same threshold — a needle tip on a coarse level — and `abs()` hides both. **Guard both**: clamp
`|1 -/+ 2H*d + K*d^2|` away from zero on numerator and denominator, and signal out-of-range rather
than dividing or emitting a negative weight.

**Size the guard from a measurement, not from a round number.** The guard must not trip on images
the scheme legitimately produces. Under the criterion of §1.2 — not the retired distance band —
a concave cavity of `R = 3*dx` produces `J` up to **28.9 under CIC and 483.8 under TSC**
(`mirror_cavity.py`, revision 6), corresponding to a denominator of `(1 - d/R)^2 = 4.1e-3` at the
TSC reach. A guard that fires below `J ~ 10^3`, or on a denominator above `~10^-3`, discards real
images in a geometry this tree ships. Guard on `J` itself with a counter, so the diagnostic says how
many images were refused and at what `J`, and let the numerator sign test be separate — a negative
weight and a near-singular denominator want different messages even though they share a threshold.

**A refused image is deposited with `J = 1`, not dropped.** It is the same fallback the fit uses for
a crease and for a rank-deficient neighbourhood (§3.3), and it is the conservative choice in both
senses: the image still carries its mass, so the even extension is still present and the cut cell
still reads a density rather than `kappa*n`; and the error it carries is bounded by §3.1's
no-Jacobian column — +117% at `R = 4*dx` in the worst case — rather than being unbounded, which is
what a `J` of several hundred applied to the wrong particle would be. Dropping the image instead
would reintroduce exactly the deficit the mirror exists to remove, in precisely the cells where it is
largest. Count the fallbacks; a run in which they are not rare is a run whose geometry is
under-resolved, and that is a tagging problem, not a deposition problem.

(The separate claim that the volume fraction `kappa` never appears in a denominator is **correct**
and is the real advantage over hybrid-plus-redistribution.)

### 3.3 Curvature from the discrete normal

The discrete normal lives only on cut cells, a codimension-1 scatter, so `grad(nhat)` cannot be
taken with a Cartesian stencil. Fit the shape operator in the tangent plane, on **physical**
centroid positions:

```
for cut cell i, with normal n_i and boundary centroid position x_i:

    t1, t2 = any orthonormal basis perpendicular to n_i
    u_j    = (t1.dx_j, t2.dx_j)      dx_j = x_j - x_i        both PHYSICAL positions
    v_j    = (t1.dn_j, t2.dn_j)      dn_j = normal(j) - n_i

    S      = argmin sum_j | v_j - S u_j |^2        2x2 least squares, then symmetrized

then store  S_c = P S P^T,  P = [t1 t2]   (3x2), and DISCARD t1, t2
```

The invariants come from `S_c` by §3.1's formulas, not from this 2×2 — do not carry `2H = tr S`,
`K = det S` forward past this point, because `det` is what changes meaning between the two forms.

This makes `S` the differential of the normal, `S = +dn̂`, which is the convention §2.2 and §3.1 are
written in. No principal directions, no eigen-decomposition, no Hessian. **The tangent frame is
local to this fit and is discarded here** — `S_c` is what is stored, and §2.2 never reconstructs a
frame (§2.1). **Use a 5³ neighbourhood** and require at least **four** usable cut neighbours (this is
what the harness enforces; revision 3 said three). Skip multi-valued neighbours — a multi-valued cell
contributes two normals to the fit.

> **Four neighbours is a count, not a rank test.** Four cut neighbours can be nearly collinear in the
> tangent plane — a ridge, a thin fin, a cylinder whose cut cells run along the axis — and then
> `U^T U` is near-singular in one direction and the least squares returns an arbitrary curvature
> there. The count is necessary and not sufficient. **Test the conditioning of the 2×2 `U^T U`** and
> fall back to `J = 1` with the same counter as the crease fallback when it fails. This is separate
> from the residual test below: the residual detects a surface that is not smooth, the conditioning
> test detects a *sample* that cannot determine `S` even where the surface is smooth.

`mirror_discrete_curvature.py fit`, 5³ throughout:

| Sphere, dx = 1 | true 2H | `kappa <= 0.05` | err | `kappa > 0.75` | err | \|dJ/J\| 3³ | \|dJ/J\| 5³ | linearized |
|---|---|---|---|---|---|---|---|---|
| convex R = 8 | 0.25000 | 0.2492 | −0.32% | 0.2484 | −0.64% | 1.2% | **0.5%** | 48.4% |
| convex R = 6 | 0.33333 | 0.3310 | −0.70% | 0.3377 | +1.31% | 1.4% | **0.7%** | 62.8% |
| convex R = 4 | 0.50000 | 0.5142 | **+2.84%** | 0.5057 | +1.14% | 2.7% | **1.4%** | 77.3% |
| concave R = 6 | −0.33333 | −0.3406 | **+2.18%** | −0.3339 | +0.17% | 1.4% | **0.7%** | 26.5% |

> **What the `kappa` bins actually say.** Revision 3 read this as "it does not degrade at small
> `kappa` — the sliver bin is as good as the fat bin", tabulating the **3³** fit while shipping
> **5³**. On 5³ the sliver bin is *worse* where it matters: +2.84% against +1.14% at convex
> `R = 4*dx`, +2.18% against +0.17% concave. What *is* true is that the p05–p95 **spread** is wider
> at high `kappa` in three of the four cases (concave `R = 6` is the exception: 0.0084 in the sliver
> bin against 0.0068 at `kappa > 0.75`). Spread is not bias.
>
> So the small-`kappa` risk is resized, not retired. It is not "the fit falls apart in slivers" —
> the fit is well-conditioned there and the spread is narrow. It is "the fit carries a systematic
> 2–3% offset in slivers at tight radii". Narrow spread over 24–188 cells means it is systematic and
> will not average away. **Budget 2–3% and bin the acceptance test by `kappa`.**

Sensitivity to discrete-normal error, sphere `R = 6*dx` (`mirror_discrete_curvature.py noise`):

| normal noise | 3³ convex | 5³ convex | 7³ convex | 3³ concave | 5³ concave | 7³ concave |
|---|---|---|---|---|---|---|
| 0° | 1.4% | 0.7% | 0.5% | 1.4% | 0.7% | 0.5% |
| 1° | 2.0% | 1.0% | 0.6% | 2.1% | 1.0% | 0.6% |
| 2° | 3.3% | 1.2% | 0.7% | 3.1% | 1.3% | 0.8% |
| 5° | 6.3% | 2.5% | 1.7% | 6.9% | 2.7% | 1.7% |
| 10° | 11.7% | **5.0%** | 3.7% | 11.6% | **5.0%** | 3.7% |

Both signs are tabulated because revision 4's "concave matches convex to within 0.4 points at every
entry" is false at 5°/3³ (6.3% against 6.9%).

At 10° a 5³ fit still returns 5.0% against the linearized form's 62.8%. That margin is what the
*fit* buys: a two-point difference of normals over `2*dx` would turn 10° into a 52% curvature error.

**The fit residual is a crease detector, for free.** It is large exactly where the surface is not
locally smooth — a sphere union in `Aerosol`, a facet edge in `Tessellation`, a dielectric/electrode
corner. Threshold it, fall back to `J = 1` with a counter, and never trust a shape operator fitted
across a kink.

---

## 4. The deposition pass

### 4.1 Estimator

```
phi(x_j) = (1/dx^D) sum_p g_p [ W(x_j - x_p) + J_p * W(x_j - R(x_p)) ]     for kappa_j > 0
```

where `g_p` is whatever the caller is depositing — the weight column for `depositWeight`, the value
of the gather lambda for `depositGathered`.

The estimator is the kernel density estimate of the **even extension** of the density about the
surface. Where the true density is uniform, the even extension is uniform, so the estimate is exact
up to the wall. For a wall at `a = 0.5` in 1-D, plain CIC gives `phi_0 = 0.5n`; the mirror gives
exactly `n`. The apparent doubling of a particle on the boundary is the correction, not a double
count — in the physical measure it cancels: `kappa_0 * phi_0 * dx = 0.5 * 2m = m`.

**No per-particle normalization.** The issue proposes scaling each particle by
`V = sum_{kappa>0}(w + w')` so it deposits exactly `m`. That is wrong and must not be implemented.
`V` is 1 only for a wall on a cell **face**; for a wall at `a = 0.3` a particle at `0.4` has
`V = 1.6`, and for a wall through a cell *centre* at `a = 0.5` a particle at `0.6` has `V = 1.9`.
Those values are what makes the ensemble integral close, because a cut cell's value is a density
over the whole cell while its particles occupy only the fluid fraction. Dividing by a
position-dependent `V` reintroduces exactly the bias the mirror removes.

**Mass that lands in covered cells is ignored, not removed.** The estimator is only ever read at
`kappa_j > 0`, so the part of an image's cloud that falls in covered cells is irrelevant by
construction — and every image is inside the solid, so a large share of it falls there. The
deposition kernels do not drop it: `depositParticleCIC` and `depositParticleTSC` write
`rho(iv, ...)` for every cell in `cloudBox` with no EB test (`CD_EBParticleMesh.H:725-737,
779-798`). That is pre-existing behaviour — today's fluid-side clouds already spill into covered
cells — and the mirror raises the magnitude. **The deposit does not clean it up** (§4.2); the
question of which downstream readers treat a covered cell as data is O2.

### 4.2 The pass, step by step

```
0. clear     m_mirrorImages.clearParticles()      -- the container is persistent, not fresh

1. build     for each VALID particle inside the band, append an image at R(x_p) into
             m_mirrorImages (ParticleContainer<NoPayload>), carrying weight g_p*J_p

2. remap     m_mirrorImages.remap() -- the ordinary remap, no special level rule

3. deposit   the image container through the ORDINARY EBAMRParticleMesh::depositWeight,
             with the ordinary kernel and coarse-fine strategy, into m_mirrorScratch

4. add       DataOps::incr(real, m_mirrorScratch, 1.0)
```

**There is no covered-cell cleanup step, and that is deliberate.** Revisions 5 and 6 had one. See
below for why it came out.

**Two containers and two mesh fields — and the second mesh field is forced by more than the reset.**
The obvious reading is that `m_mirrorScratch` exists only because every AMR deposit path opens with
`DataOps::setValue(a_meshData, 0.0)` (`CD_EBAMRParticleMesh.H:641, 844, 943, 1217`), so a second
`depositWeight` into `a_phi` would erase the first. That is true, and it is the weaker of the two
reasons. The stronger one is that **the level exchange accumulates over ghost cells and never clears
them.**

`DataOps::setValue` zeros the whole `EBCellFAB`, ghosts included (`EBCellFAB::setVal`), so a deposit
starts from a clean ghost region. The level exchange is a *reversed* copier with an add op —
`m_levelCopiers[lvl].define(dbl, dbl, domain, m_ghost*IntVect::Unit, true)` then `.reverse()`
(`:591-592`), used at `:662` — so the ghost region is the exchange's **source** and the owner's
valid region its destination. The ghost cells keep what was deposited into them; nothing writes them
back to zero.

So a hypothetical "deposit again without resetting" would fold the *first* pass's ghost content into
the valid region a second time. **Double counting, silently, proportional to the ghost-region share
of each patch.** A scratch field is therefore not an ergonomic convenience that a `bool a_reset`
parameter could remove — the reset and the exchange are one mechanism, and separating them is a bug.
Say this where the scratch is allocated; it is the kind of thing that reads like dead weight to
someone optimising later.

The particle side needs its own container for a plainer reason: the source container arrives as
`const ParticleContainer<P, Traits>&` (`:251-255`) and the images are new particles. Two evaluated
alternatives, both rejected:

- **Merge images into a copy of the source container and deposit once.** Saves the second pass over
  the coarse-fine machinery and the whole-domain `incr`, but costs a copy of *every* source particle
  on *every* deposit to gain a saving proportional to the *images*, which are a thin band. It also
  forces the combined container to carry the source's payload type, so every image would carry an
  `ItoParticle`'s columns to hold one weight.
- **Deposit images into `a_phi` first, sources second.** Same double-fold problem in the other
  direction, plus it makes the mirror pass reorder the ordinary deposit.

`m_mirrorScratch` carries one component, because `depositWeight` asserts
`a_meshData[0]->nComp() == 1` (`:639`).

> **What `DataOps::incr` does when the ghosts differ.** `EBCellFAB::plus(src, scale)` adds over
> `a_src.getRegion() & getRegion()` (Chombo `EBCellFAB.cpp:591`) and asserts only that the component
> counts match. So a scratch allocated with a different ghost width than the caller's `a_phi` does
> **not** assert and does **not** lose anything in the valid region — both regions contain the same
> `dbl[din]` — it simply skips the outer ghost ring. That is correct here, because after step 3 only
> the valid region of the scratch is meaningful anyway: the deposit's `exchange(..., EBAddOp())`
> (`CD_EBAMRParticleMesh.H:662`) folds ghost mass *into* its owner and never refreshes the ghosts.
> Document this where the scratch is allocated, so nobody later "repairs" a mismatch that is not a
> bug.

> **Why there is no covered-cell reset, and why that is not a loose end.** Revisions 5 and 6 ended
> the pass with `DataOps::setCoveredValue(real, getCoveredCells(realm, phase), 0.0)`. It is gone,
> for the same reason `coarsenAndFillGhosts` is not called here either (`03cf59690`): **the deposit's
> contract is that the field is correct where it is read, and it is read at `kappa > 0`.**
> Covered-cell content, ghost-cell content and cross-level synchronization are all the consumer's,
> and the consumer is the only party that knows whether it needs them.
>
> Every image is inside the solid by construction — that is what reflection means — so images
> landing in covered cells is the expected outcome, not a defect to be cleaned up. §4.1's estimator
> is defined only at `kappa > 0`, so nothing the mirror claims depends on what a covered cell holds.
>
> Four facts make deferring it safe rather than merely convenient:
>
> 1. **The dominant coarsening path is `kappa`-weighted.** `ItoSolver::coarsenAndFillGhosts`
>    (`CD_ItoSolver.cpp:2121`) is `conservativeAverage` + `interpGhost`, and a conservative average
>    weights by `kappa`, so a covered fine cell contributes exactly zero to the coarse cut cell
>    however large its value.
> 2. **Plot output already resets covered cells itself**, at `CD_ItoSolver.cpp:1813`, immediately
>    before copying to the output realm.
> 3. **The elliptic solvers never see them.** A covered cell has no `VolIndex`, so it is not in any
>    stencil the field solve assembles.
> 4. **Nothing zeroes covered cells after a particle deposit today either.** All nineteen
>    `setCoveredValue` call sites are plot paths or the CDR solvers' own state
>    (`CD_CdrSolver.cpp:1261`, `CD_ItoKMCStepperImplem.H:4890`), never `ItoSolver`'s deposited
>    `m_phi`. A reset here would have been a new responsibility, not the restoration of an old one.
>
> What is genuinely left is the *consumer-side* question — which readers treat a covered cell as
> data — and the mirror raises the magnitude of what they would see. That is **O2**, where it
> belongs, and O2 already names `DataOps::filterSmooth` and the `arithmeticAverage` sites. Do not
> re-add a reset here to paper over an O2 item.

**Step 1 needs no kernel.** Building an image is position → cell → stored `(x_c, n_c, S_c)` → `R(x_p)`.
It is an `AmrMesh` private helper, not an `EBParticleMesh` concern.

**The image container is `ParticleContainer<NoPayload>`.** `NoPayload` (`CD_ParticleSoA.H:288`) is
the container's default template argument (`CD_ParticleContainer.H:121`) and gives position plus the
container-owned mandatory weight column and nothing else. It is an established idiom here —
`McPhoto.cpp:132`, `CdrSolver.cpp:1286`, `ItoKMCStepper.H:652`, three containers in
`ItoKMCGodunovStepper.H`. Because the mirror is linear in the deposited quantity, the gather lambda
is **evaluated once on the source particle at build time** and `J_p * g_p` stored as the image's
weight; the gather is never applied to an image. So a single scalar per image suffices, the image
container is a concrete non-template type, and step 3 is a plain call to the *existing*
`EBAMRParticleMesh::depositWeight` (`CD_EBAMRParticleMesh.H:251-255`) with a different container —
**no interface change anywhere in the deposit machinery.**

The one thing that would break this is the multi-component `AmrMesh::depositParticles<Members...>`
(`CD_AmrMeshImplem.H:410-431`). It currently has **zero callers**; if it gains one, the image
container needs `NCOMP` scalars.

**The remap is the ordinary one.** Revision 4 proposed a "level-preserving" remap that routes each
image to the finest level *at or below its source's level*, on the grounds that
`ParticleContainer::remap()` promotes an image landing under a fine grid and deposits a coarse
particle's weight through a fine-width cloud. Review 4 refuted it in both directions:

- **`Transition` already deposits at fine width.** `transferMaskParticlesTransition`
  (`CD_EBAMRParticleMesh.H:1128`) *moves* coarse particles in the transition band out of the valid
  holder, and the core deposits them on the refined-coarse grid at `widthScale = 1.0` (`:1276`,
  documented at `:1188`). Promotion is the tree's own policy for a coarse particle near the
  interface, so a promoted image simply inherits it.
- **Level preservation would silently lose the correction.** The coarse-fine masks are narrow —
  `getTransitionMaskWidth` (`:998-1014`) is `refRat/2` fine cells for CIC and `refRat` for TSC, i.e.
  0.5 and 1 coarse cells; `Halo`'s `coarseMaskWidth` is 1 coarse cell (`:847`). An image sits up to
  ~2.6 coarse cells from its source. Pin it to the coarse level while its position lies inside the
  fine region and it falls outside both bands, deposits at coarse width into coarse cells under the
  fine grid, and `ItoSolver::coarsenAndFillGhosts` (`CD_ItoSolver.cpp:2121`) then overwrites them
  via `conservativeAverage`. Plain `remap()` cannot do that, because a cell covered by a finer level
  is never a valid destination (`m_validCells`, `CD_ParticleContainer.H:928`).

Promotion deposits the correction at the wrong bandwidth but keeps it — mass conserved, mean
correct, a variance defect. Level preservation can discard it with no signal. Plain `remap()` wins,
and it needs no change to `ParticleContainer` at all.

Under `Halo` — the one strategy that deposits at `widthScale != 1`, at exactly one line
(`CD_EBAMRParticleMesh.H:897`) — a promoted image is deposited at `widthScale = 1` while its source
gets `refRat`, a bandwidth mismatch across the interface. `Halo` is not the target strategy, and that
is now measured rather than asserted: **all 62 shipped `deposition_cf` settings are `transition`**,
and the only two hardcoded `CoarseFineDeposition::Halo` call sites
(`CD_ItoKMCStepperImplem.H:4046`, `CD_ItoKMCGodunovStepperImplem.H:1407`) are exactly the two that
map to `Native` and `NGP` under §5.5 and therefore never run a mirror pass at all.
**Document this in `depositHaloCore`'s Doxygen rather than engineering around it.** The same
disposition applies to the band-mask width on that path (former O5): `widthScale` is 1.0 at every
other patch deposit in the class (`:658, :866, :966, :967, :1235, :1276`), so the reflect decision
taken at build time matches the deposit everywhere except `Halo`.

### 4.3 Build images from valid particles only

`depositHaloCore:849` is `copyMaskParticles` — a **copy**. Under `CoarseFineDeposition::Halo` a
coarse halo particle is genuinely deposited twice, once on the coarse level and once onto the
refined-coarse grid at `widthScale = refRat` (`CD_EBAMRParticleMesh.H:880`); it is not double
counting only because `conservativeAverage` later overwrites the coarse cells under the fine grid.
`HaloNGP` (`:949`) and `Transition` (`:1220`) move instead of copying.

So images must be built from the **valid** holder only. The constraint is *which holder*, not
*when*: all three coarse-fine cores restore the valid holder before returning — `depositHaloCore`
copies then clears (`CD_EBAMRParticleMesh.H:849, 912`), `depositHaloNGPCore` and
`depositTransitionCore` transfer out and transfer back (`:988`, `:1292`) — so building images after
the real deposit, as §5.3's sequence does, is safe. Revision 5 said "before any deposit", which
contradicts §5.3 and would send the next implementer to reorder correct code. Build them from a mask
holder and `Halo`'s copy generates a second image — issue #29's double-counting trap 1 in a
new guise. The image container then goes through the identical coarse-fine pass with its own mask
copy, so it must be a container **defined on the realm's grids and reused rather than reallocated
per deposit**. (Revision 4 said it must be "on the same realm with the halo masks available"; the
masks live in `EBAMRParticleMesh` and the mask *holder* lives in the container, so any container
defined on those grids qualifies.)

---

## 5. Where this slots in

The split is: **regrid-lifetime state and per-call scratch on `PhaseRealm`; orchestration in
`AmrMesh`; `EBAMRParticleMesh` unchanged.**

### 5.1 `PhaseRealm` — all the state

`m_irregularCells` is the direct precedent for the surface data:

- declared `EBAMRCellData m_irregularCells` (`CD_PhaseRealm.H:600`)
- built in `PhaseRealm::defineMasks(a_lmin, a_numGhost)` (`CD_PhaseRealm.cpp:584`) — allocate per
  level with `EBCellFactory(ebisl)`, fill on valid cells, then `exchange()` (`:657-660`)
- called from `PhaseRealm::regridBase(a_lmin)` (`:189`), after `defineEBLevelGrid` (`:163`)

Five members are added, and the file's own convention decides which are `mutable`. `PhaseRealm`
declares every operator and per-call scratch `mutable` (**eighteen** of them, `:605-695` — revision 5
said seventeen; the members are at `:605, 610, 615, 622, 630, 635, 640, 645, 650, 655, 660, 665, 670,
675, 680, 685, 690, 695`, including `m_particleMesh` at `:675`, `m_redistributionOp` at `:665` and
`m_nonConservativeDivergence` at `:695`) while regrid-built data such as `m_irregularCells` (`:600`)
is not:

```cpp
// Regrid-lifetime, so not mutable -- the m_irregularCells convention.
EBAMRCellData                        m_surfaceData;    // 13 comps: ok, x_c, n_c, S_c

// Per-call scratch, so mutable, following the file's own convention.
mutable EBAMRCellData                m_mirrorScratch;  // 1 comp, 4.2 step 3
mutable ParticleContainer<NoPayload> m_mirrorImages;
mutable EBAMRIVData                  m_massDiff;       // off the solvers, see 5.4
mutable EBAMRIVData                  m_depositionNC;
```

**`m_mirrorScratch` is the field §4.2 step 3 deposits into.** Revision 5's member list omitted it
entirely while §4.2 and §5.3 both used it. It is per-call scratch, so `mutable`, and one component.

> **`EBAMRCellData`, chosen for the lookup, with the memory cost accepted.** The reflect decision
> and the `(x_c, n_c, S_c)` read are on the per-particle hot path: one `FArrayBox` read at a known
> `IntVect`, no indirection. A `BaseIVFAB` over the band would be far smaller but is indexed through
> a per-patch offset map, and the whole point of §2.2 is that reflection is *one dependent load and
> about twenty flops*. **This was decided; do not re-propose the sparse container.** What a future
> reviewer may legitimately reopen is the *number of components*, not the container.
>
> The cost is real and is accepted with eyes open. At the shipped defaults —
> `AmrMesh.coarsest_domain 128 128 128`, `max_block_size 16`, `num_ghost 2`, so a 1.91× ghost
> overhead on 16³ boxes — one component costs 33 MB per level. Thirteen components cost **~430 MB
> per level, per phase, per realm**; a two-phase dielectric problem on two realms pays four times
> that. For comparison the four masks the realm already carries (`m_regularCells`, `m_coveredCells`,
> `m_notCoveredCells`, `m_irregularCells`) cost 131 MB per level between them. **Report it in the
> memory report** alongside the existing masks, so the number is visible rather than inferred.
>
> Note that 12 payload components is the floor for *any* correct format, not a cost of the
> frame-free one. Revision 5's "9 Reals" was only 9 because it left the tangent frame unspecified;
> the honest alternative — a 2×2 `S` plus an explicitly stored `t1` — is 3 + 3 + 3 + 3 = 12 as well.
> The frame-free format buys §2.1's determinism for nothing.

`m_surfaceData` is built **with the ordering caveat in §5.4**. Component 0 is the `ok` flag of §2.1,
zeroed at allocation, so an unwritten cell is a no-reflect cell by default rather than by accident.

> **Build it in `regridOperators`, not in `defineMasks`.** `defineMasks` sits in `regridBase`
> for a stated reason — *"so the cell masks are available to load-balancing routines, which run
> after regridBase but before regridOperators"* (`CD_PhaseRealm.cpp:182-183`) — and is
> unconditional. Nothing in load balancing wants a curvature fit, and an unconditional build makes
> every realm pay the 5³ fit and the 430 MB whether or not it deposits particles. Put the build beside
> `defineParticleMesh` (`CD_PhaseRealm.cpp:327`, from `regridOperators` at `:203`) and gate it on
> `queryOperator` like every other operator in that routine.

`m_massDiff` and `m_depositionNC` need no new plumbing: `PhaseRealm` already owns the operators that
consume them (`m_redistributionOp`, `m_nonConservativeDivergence`) behind the `s_eb_redist` and
`s_noncons_div` registrations (`:56-57`).

`m_mirrorImages` needs four things `PhaseRealm` does not currently hold — `minBlockSize`,
`levelTiles`, `validCells` and the realm name (`ParticleContainer::define`,
`CD_ParticleContainer.H:193-202`). All four live on `Realm` (`CD_Realm.H:643, 653`) or `AmrMesh`,
are phase-independent, and `Realm` owns the `PhaseRealm`s, so this is one hop of plumbing.
`AmrMesh::allocate(ParticleContainer<P,Traits>&, realm)` (`CD_AmrMeshImplem.H:238-261`) shows exactly
which getters to forward.

**Its lifecycle is three lines, and all three are load-bearing.** It is `define`d at regrid, beside
the surface data, because `remap()` reads `levelTiles` and `validCells` which the regrid rebuilds.
It is `clearParticles()`d at the *start* of every mirror pass (§4.2 step 0), not the end, so that an
aborted pass cannot leak images into the next deposit. And `ParticleSoA::append(pos, weight)` leaves
the metadata columns at `s_invalidID`/`-1` (`CD_ParticleSoA.H:949-958`), which is fine: nothing in
`remap()` reads them, and images need no identity. Do not call `resetParticleIDs` on them.

Gate all five on an operator registered through `PhaseRealm::registerOperator` (`:353`) so realms
that never deposit particles do not pay for it, in the same way `s_particle_mesh` is gated
(`CD_PhaseRealm.H:60`).

> **Invariant to assert.** These members are shared by every solver on that (realm, phase), where
> today each solver owns its own `m_massDiff`/`m_depositionNC`. That is safe only because solvers
> deposit serially — `ItoLayout`'s iterator is a serial loop and the OpenMP parallelism inside a
> deposit is over boxes. `mutable` removes the last compile-time hint that anything is being
> written, so add a `mutable bool m_depositionInFlight` asserted on entry and exit. It costs nothing
> under `OPT=HIGH` and it fires the day someone parallelises the solver loop.

### 5.2 `EBAMRParticleMesh` — unchanged

`EBAMRParticleMesh` owns the kernels, the coarse-fine strategies and the halo/transition masks, and
holds `Vector<RefCountedPtr<EBLevelGrid>> m_eblgs`, `m_refRat`, `m_dx`, `m_probLo`, `m_ghost`
(`CD_EBAMRParticleMesh.H:451-511`).

Revision 4 stated that it "**cannot allocate the image `ParticleContainer`**, which needs a realm".
That is false: `ParticleContainer::m_realm` is a pure label — assigned at `CD_ParticleContainer.H:218`,
read back only by `getRealm()` (`:300`), stored at `:1013`, used for nothing but the
`CH_assert(a_particles.getRealm() == m_realm)` checks in solvers. What the class actually lacks is
`minBlockSize`, `levelTiles` and `validCells`, and it could be given them. It should not be: once the
redistribution scratch joins the image container, route (a) means `EBAMRParticleMesh` also needs both
operators and the `BaseIVFAB` factories, at which point it carries half of `PhaseRealm`.

With the state on `PhaseRealm` and the sequencing in `AmrMesh`, **`EBAMRParticleMesh` needs no change
at all.** Its `depositWeight` (`:251-255`) already takes the container as a parameter and is not
`const`; step 3 calls it a second time with `ParticleContainer<NoPayload>`.

### 5.3 `AmrMesh` — the orchestration

`AmrMesh::depositWeight` (`CD_AmrMeshImplem.H:368`) and `AmrMesh::depositGathered` (`:388`) are the
two funnels every deposit passes through and the lowest level with both a realm and the particle
machinery. Neither is `const` (`CD_AmrMesh.H:991-997`). The full sequence lives there:

```
AmrMesh::depositWeight(...)
  ├ particleMesh.depositWeight(real, sources, ...)          // unchanged
  ├ if mode == Mirror && this (realm, phase) has cut cells:
  │   ├ m_mirrorImages.clearParticles()
  │   ├ build images into m_mirrorImages from m_surfaceData, VALID holder only
  │   ├ m_mirrorImages.remap()
  │   ├ particleMesh.depositWeight(m_mirrorScratch, m_mirrorImages, ...)
  │   └ DataOps::incr(real, m_mirrorScratch, 1.0)
  └ if mode == Redistribute | RedistributeBlended:
      └ nonConservativeDivergence + hybrid + redistribute, using m_depositionNC / m_massDiff
```

### 5.4 Ordering inside the regrid build

1. On each **valid** cut cell, read `normal(vof)` and `bndryCentroid(vof)`, **convert the centroid
   to a physical position via `Location::position`**, and fit `S` over the 5³ cut-cell
   neighbourhood. Skip multi-valued neighbours. Record the fit residual.
2. Where the fit has fewer than **four** usable neighbours, or the residual exceeds the crease
   threshold, mark the cell — those band cells get `J = 1` and a counter, never a rank-deficient `S`.
3. **Exchange the fitted cut-cell values into ghosts, then extend** `(x_c, n_c, S_c)` to every band
   cell from the nearest cut cell — or equivalently extend by two 1-cell propagation sweeps with an
   exchange between them.
4. Require `eb_ghost >= 2` **and** `num_ghost >= 2` at runtime.

> **"Nearest" means Euclidean, over the radius-2 neighbourhood, with a stated tie-break.** §2.4's
> error ladder — the 1.2–2.9% that justifies the quadratic patch over everything else — was measured
> with `near = argmin(|x_band - x_cut|^2)` over **all** cut cells (`mirror_source.py:48-53`). Two
> 1-cell propagation sweeps do not compute that: they compute a Chebyshev-nearest assignment with an
> arbitrary tie-break, and the two differ exactly where more than one piece of surface is within two
> cells — the O3 geometry. An explicit search over the 5³ neighbourhood costs 125 distance
> comparisons per band cell once per regrid and is faithful to the evidence. Fix the tie-break
> (lowest `IntVect` in lexicographic order) so the result does not depend on iteration order.

> **A band cell whose nearest cut cell is not one this level owns must not reflect.** `EBISBox` is
> valid `eb_ghost` cells out *whether or not any box on that level owns those cells* — that is what
> `new EBLevelGrid(m_grids[lvl], m_domains[lvl], m_numEbGhostsCells, &(*m_ebis))`
> (`CD_PhaseRealm.cpp:411`) buys — while `exchange()` fills ghosts only from *valid* regions of the
> same level. So on a level whose grid edge passes within two cells of the EB, a valid band cell can
> see a cut cell through `EBISBox` that no box on that level ever fits an `S` for, and no exchange
> can deliver one. Revision 5 wrote step 3 as a total statement and it is not.
>
> The cheap fix is to derive the band from data that was actually delivered rather than from
> geometry that merely exists. The existing masks already work this way: `defineMasks` raises
> `m_irregularCells` on `ebisbox.getIrregIVS(dbl[din])` — the **valid** box — and then exchanges
> (`CD_PhaseRealm.cpp:646, 660`), so a cut cell nobody owns never raises the mask. Build the band by
> growing the *exchanged* irregular mask, set the `ok` flag of §2.1 only on band cells that found a
> delivered cut cell, and **count the band cells that did not**. A nonzero count is not a warning to
> be tuned away — it means the grids do not cover the EB band, which is the same defect as O9 seen
> from the other end, and both are repaired by tagging the band rather than only the cut cells —
> see the note under O9, which corrects an earlier claim that `AmrMesh.buffer_size` was the knob.

> Revision 3 extended first and exchanged once afterwards. In that order a band cell at a patch edge
> whose nearest cut cell lies one patch over is assigned a too-far cut cell, and the exchange cannot
> repair it — both patches computed from the same incomplete data and agree on the wrong answer.
> This is the one place where `defineMasks`' valid-then-exchange pattern does **not** transfer.

### 5.5 Selection: one enum, one parser

The `forceIrregNGP` bool becomes an `IrregularDeposition` enum. It absorbs redistribution, because
`Native`, `NGP`, `Mirror` and the two redistribution schemes are five answers to the *same* question
— what do we do about a cut cell holding `kappa*n`? — and only one of them can be in force. Making
them one enum makes `mirror + redistribute` **unrepresentable** rather than something a cross-check
has to catch.

```cpp
// CD_IrregularDeposition.H, beside CD_DepositionType.H
enum class IrregularDeposition
{
  Native,              ///< Deposit as-is; phi is an extended state into the EB.
  NGP,                 ///< Whole cloud into the cut cell.
  Mirror,              ///< Even extension across the EB (this plan).
  Redistribute,        ///< Hybrid update, deltaM = (1-kappa)*dc smooshed to neighbours.
  RedistributeBlended  ///< As above, blended with the non-conservative divergence.
};

inline IrregularDeposition irregularDepositionFromString(const std::string&) noexcept;
```

This follows `ParticleManagement::ParticleMergeMethod` + `mergeMethodFromString`
(`CD_ParticleManagement.H:39-75`) exactly — documented enumerators, one `…FromString` in the same
header, abort on an unknown selector. Today every solver open-codes its own `if (str == "ngp")`
chain (`CD_ItoSolver.cpp:264` and again at `:298` for `plot_deposition`), which is how the pattern
keeps getting copied.

**How the five map onto what ships today.** `ItoSolver` gates the whole hybrid-plus-redistribute
block on `m_useRedistribution` (`CD_ItoSolver.cpp:2153`), so `blend_conservation` is only ever read
inside it; `McPhoto` has no `redistribute` option at all (`CD_McPhoto.options:8`,
`CD_McPhoto.cpp:236`) and gates its redistribution on `blend_conservation` (`:1339`). With
`blend = false` McPhoto's `depositHybrid` is an identity (`divH = dc + (1-kappa)*0`) whose `deltaM`
is discarded, so it is dead work, not a different scheme.

| mode | `divH` | `deltaM` redistributed | ItoSolver today | McPhoto today |
|---|---|---|---|---|
| `Native` | `dc` untouched | no | `redistribute=false` (7 files) | `blend=false` (14 files) |
| `Redistribute` | `dc` | `(1-kappa)*dc` | `redistribute=true, blend=false` (**23 files**) | *not expressible* |
| `RedistributeBlended` | `dc + (1-kappa)*dnc` | `(1-kappa)(dc - kappa*dnc)` | both true (0 shipped) | `blend=true` (0 shipped) |

**The migration is mechanical and behaviour-preserving, because `mirror` is not the default.** The
enum lands as a *signature* change: every shipped input is rewritten to the value that reproduces
exactly what it does today, and nothing selects `Mirror` until someone opts in deliberately. That
removes the whole class of "which of these 30 files did we mean to change" risk from this work.

The cross-product of the two flags in shipped inputs is what makes that possible, and it is not
obvious until measured — the enum collapses two independent booleans into one selector, so it can
only be behaviour-preserving if no shipped file uses a combination the enum cannot express:

| `irr_ngp_deposition` | `redistribute` | `blend_conservation` | files | becomes |
|---|---|---|---|---|
| `false` | `true` | `false` | 22 | `redistribute` |
| `false` | `false` | `false` | 5 | `native` |
| `true` | `false` | `false` | 2 | `ngp` |
| `true` | `true` | — | **0** | *not expressible — and not used* |

(29 files set `irr_ngp_deposition`; a 30th sets `redistribute = true` alone and inherits `false` for
the other from `CD_ItoSolver.options`, so it maps to `redistribute` too.) Every shipped
configuration maps one-to-one. McPhoto's 14 `blend_conservation = false` files map to `native`,
which is behaviour-preserving because with `blend = false` its `depositHybrid` is the identity
`divH = dc + (1-kappa)*0` and the resulting `deltaM` is discarded — dead work, not a different
scheme.

Two consequences to record in the PR description, not to discover later:

- **The enum makes `NGP` + `redistribute` unrepresentable.** It is representable today and no
  shipped input uses it. Dropping it is deliberate: the combination is semantically odd (put the
  whole cloud in the cut cell, then redistribute a `(1-kappa)` share of it out again).
- **McPhoto gains a mode** it cannot currently express. That is a behaviour addition, available but
  unselected.

**Do not sweep in the fluid solvers.** `CdrGodunov.blend_conservation` (7 files) and
`CdrCTU.blend_conservation` (34 files) are the *fluid* hybrid divergence and have nothing to do with
particle deposition, despite the identical option name.

**One check the enum does not subsume.** `IrregularDeposition::NGP` and `DepositionType::NGP` are
different things, and `Mirror` with an NGP base kernel is meaningless: NGP puts the whole cloud in
the particle's own cell, so an image in the same cell doubles it.

Revision 5 placed this check at `AmrMesh::depositWeight` or in `EBParticleMesh::depositCore`'s
switch, on the grounds that `CD_ItoKMCStepperImplem.H:4041` passes `DepositionType::NGP` as a literal
"so `parseOptions` structurally cannot see it". That reads the call site wrong. `:4041` passes
`DepositionType::NGP` **and** `a_forceIrregNGP = false` (`:4047`), which maps to
`IrregularDeposition::Native` — it can never be `Mirror`. The other hardcoded site,
`CD_ItoKMCGodunovStepperImplem.H:1338`, passes `true` (`:1344`), which maps to `NGP`. Both are
behaviour-preserving under a mechanical conversion and neither can produce the collision.

The combination that *can* arise is a solver configured with `deposition = ngp` and
`irregular_deposition = mirror`, and **both of those are options of the same solver**, parsed
thirty lines apart in `ItoSolver::parseDeposition` (`CD_ItoSolver.cpp:263` and `:320`). Put the check
there, where it can name both offending keys in the error message, and keep a cheap defensive
`CH_assert` at `AmrMesh::depositWeight` for the literal call sites.

**The rename.** `irr_ngp_deposition` → `irregular_deposition` across 29 `.options`/`.inputs`
occurrences in 29 files, plus `Ito.rst` and the `pp.get` at `CD_ItoSolver.cpp:320`. It is a
value-domain change, not only a key rename: `true/false` becomes one of five strings, mapped by the
table above. Because the tree uses `pp.get` and never `query`, a stale input file hard-errors on the
missing key — that is the failure mode we want, so do not add a fallback.

**No shipped input is left selecting `mirror`.** The default in `CD_ItoSolver.options` becomes the
value that reproduces today's shipped default (`false/false` → `native`), and `mirror` is opted into
by hand. Whether it later becomes a default is a separate decision, taken after the acceptance
evidence of §6.1 exists on real geometry — not in this work. Leave `irr_ngp_interp` alone; interpolation keeps its
own two-valued type, and `mirror` is meaningless when gathering.

**93 C++ occurrences across 8 files** carry the name — `CD_EBParticleMesh.H` (49),
`CD_EBAMRParticleMesh.H` (18), `CD_AmrMesh.H` (10), `CD_AmrMeshImplem.H` (10), `CD_ItoSolver.cpp`,
`CD_ItoSolver.H`, `CD_ItoSolverImplem.H` (4 between them), `CD_CdrSolver.cpp` (2) — but **93 is not
the conversion count, and a mechanical conversion of all 93 crosses the scope boundary this section
draws two paragraphs above.**

`EBParticleMesh` and `EBAMRParticleMesh` use the *same parameter name* `a_forceIrregNGP` on the
interpolation path, which "Leave `irr_ngp_interp` alone" explicitly excludes. Splitting them:

| File | total | interpolation side | deposition side |
|---|---|---|---|
| `CD_EBParticleMesh.H` | 49 | **22** | 27 |
| `CD_EBAMRParticleMesh.H` | 18 | **6** | 12 |
| `CD_AmrMesh.H` | 10 | **4** | 6 |
| `CD_AmrMeshImplem.H` | 10 | **4** | 6 |
| solvers | 6 | 0 | 6 |
| | **93** | **36** | **57** |

The 7 `m_forceIrregInterpolationNGP` members are a separate, already-excluded set.

**Stop the enum at `AmrMesh`/`EBAMRParticleMesh`; leave `EBParticleMesh` a bool.** The per-particle
kernel answers exactly one question — NGP in a cut cell, yes or no — and `Mirror`, `Redistribute`
and `RedistributeBlended` all map to *no* there, since the mirror is a separate pass (§4.2) and
redistribution happens after the deposit. Threading a five-valued enum into `depositParticleCIC`
buys nothing and forces a decision at every leaf. Converting only the deposition side of
`CD_AmrMesh.H`/`CD_AmrMeshImplem.H`/`CD_EBAMRParticleMesh.H` plus the solvers is **30 sites**, and
the leaf keeps a single-purpose contract.

---

## 6. Delivery

### PR A — mechanism, one consumer

Band mask, the per-band-cell surface data at regrid, the image container and the mirror utility.
Wire exactly one consumer — `ItoSolver::depositParticles` — as the proof.

**PR A must also deal with redistribution**, because it is one line below the deposit it wires — but
only to the extent of refusing the combination, since nothing selects `mirror` by default.
`ItoSolver::depositWeight` is `AmrMesh::depositWeight` followed by `this->redistributeAMR(a_phi)`
(`CD_ItoSolverImplem.H:71, :79`), and `depositGathered` does the same (`:148, :157`).
`redistributeAMR` (`CD_ItoSolver.cpp:2133`) opens with

> *"When we entered this routine we had `a_phi = m_i/dV` but we actually want to have
> `phi = m_i/(kappa*dV)`"*

which is exactly the assumption the mirror removes, and `depositHybrid` (`:2194`) computes
`deltaM = (1-kappa)*(dc - kappa*dnc)` under the comment *"Remember, `dc` already scaled by kappa"*
(`:2231-2234`). With the mirror, `dc` is already a density, so `deltaM` is inflated by `1/kappa`
— 20× at `kappa = 0.05` — and smooshed into the neighbouring cells. The minimum for PR A is a hard
error when `mirror` is selected with `redistribute` or `blend_conservation`; the enum of §5.5 makes
that unrepresentable instead, which is why PR B should not lag far behind.

**Done when** the acceptance test (§6.1) passes through the wired path, the `kappa <= 0.05` bin is
reported, and every other regression is bit-for-bit unchanged. The last clause is cheap to satisfy
and is the point of the whole migration strategy: nothing selects `mirror`, so *every* existing
regression must be bit-for-bit identical, and any that is not is a bug in the conversion rather than
an expected consequence of the feature.

### PR B — wire it centrally

Introduce `IrregularDeposition`, move `m_massDiff`/`m_depositionNC` onto `PhaseRealm`, move the
hybrid-and-redistribute step into `AmrMesh::depositWeight`/`depositGathered` beside the mirror pass,
and delete PR A's bespoke call and the duplicated redistribution code in the solvers
(`CD_ItoSolver.cpp:2133/2178/2193` and `CD_McPhoto.cpp:1326/1366/1382` are near-verbatim copies).

**8 direct callers in 5 files**: `CD_McPhoto.cpp` ×2 (`:1288, :1319`),
`CD_TracerParticleSolverImplem.H:374`, `CD_ItoKMCStepperImplem.H` ×2 (`:4041, :6134`),
`CD_ItoSolverImplem.H` ×2 (`:71, :148`), `CD_ItoKMCGodunovStepperImplem.H:1401`.

The eighth is the awkward one: the `hasDielectrics` branch of `computeSemiImplicitRho`
(`CD_ItoKMCGodunovStepperImplem.H:1270`), depositing onto `phase::solid` with a hardcoded
`forceIrregNGP = true` (`:1344`). The mechanical conversion of that literal is
`IrregularDeposition::NGP`, which preserves today's behaviour exactly, so **PR B does not have to
decide anything here — it has to resist deciding.** The hazard is only if someone later routes the
solver's setting into that call: `phase::solid`'s `ebisbox.normal()` points into the *solid*, so a
mirror there reflects across the same embedded boundary with the fluid side swapped, pushing
particles deliberately placed *inside* the dielectric back out into the gas. Convert the literal to
`NGP`, and say in a comment why it is a literal.

**Why B is not optional.** If the mirror pass stays the consumer's responsibility,
`irregular_deposition = mirror` becomes a setting that some deposits honour and others silently
ignore — character-for-character the defect PR #700 just fixed. No assert can catch it, because a
deposit function cannot know whether its caller intends a mirror pass afterwards.

Retire the nine ad-hoc `EBParticleMesh` constructions along the way — `McPhoto.cpp:1277,1459`,
`CdrSolver.cpp:1317`, `ItoSolver.cpp:2108,2528,2651,2686,2804`, `ItoSolverImplem.H:115` — by
registering `s_particle_mesh` on those realms and using the persistent per-patch leaves.

### 6.1 Acceptance

Three cases, mandatory, same PR. **Primary metric is the density binned by `kappa`.**

- **A curved case, and it must be the torus, not a sphere** — on a plane `J = 1` identically, so a
  planar-only test green-lights an implementation with no Jacobian and a 117% error at `R = 4*dx`.
  A sphere is barely better: it is umbilic (`c_1 = c_2`), so `S_c` is isotropic and *frame-invariant*,
  and it has constant curvature, so borrowing `(x_c, n_c, S_c)` from a neighbouring cut cell costs
  nothing. A sphere therefore cannot distinguish a correct anisotropic `S_c` from a wrong one, which
  is exactly the gap §2.1 closes and exactly the gap every harness in this directory shares. A torus
  has unequal principal curvatures **and** non-constant curvature, and it exercises the
  cut-cell-to-band-cell extension. **Assert in the test that the sliver bin actually contains cells
  with `c_1 != c_2`** — otherwise a torus resolved too coarsely degenerates to the sphere case
  without saying so.
- **A concave case** — convex-only hides the errors that no longer cancel, and the sliver bias is
  worst there.
- **A case with no EB at all**, bit-for-bit identical to `native`.

Report the `kappa <= 0.05` bin explicitly. A 2–3% systematic offset there is expected and budgeted
(§3.3); a larger one, or one that grows with refinement, is the signal that the discrete normals are
worse than the noise sweep assumed.

Score cells **at least 3 cells inside** the sampled region. A cell centre on the solid side is
reconstructed almost entirely from images, up to `sqrt(3)*dx` away; a 1-cell margin fabricates a
~1.8% deficit at small `kappa`. Margin 3 is a defensible conservative choice, **not** a measured
floor — there is no plateau at 3.

`sum(kappa phi dV)/sum(m)` is a **secondary** check. It is an exact identity, equal to
`mean_p[kappa~(x_p) + J_p*kappa~(R(x_p))]` to 2.2e-16, and because it weights each cell by `kappa`
it is least sensitive exactly where the mirror helps:

| Convex sphere, cell plane + linearized J | `sum(kappa phi dV)/sum(m)` | `kappa <= 0.05` density |
|---|---|---|
| R = 16 | 1.00027 | −0.0129 |
| R = 8 | 1.00101 | −0.0572 |
| R = 6 | 1.00019 | −0.1921 |
| R = 4 | **0.99974** | **−0.4648** |

Keep it — cheap, exact, catches gross bookkeeping errors. It is a second opinion.

An **images-only build** (the zeroing-deposit bug) reads −97.5% to −100% in every regular cell, so
the planar test does catch it; only the `kappa <= 0.05` bin is nearly blind (−1.8 to −7.7%).
Revision 2 argued the opposite.

---

## 7. Open items

| # | Item | Status |
|---|---|---|
| O1 | The mirror activates issue #29 part 4 | **Downstream, not this PR.** Today the reaction step's charge error cancels because both sides are `kappa`-weighted; with the mirror the Ito side becomes a true density and `drho = -Qe*dn*(1-kappa)` appears. The fix belongs in `ItoKMCStepper::reconcileCdrDensities` and will arrive with the coupling-algorithm rework, in consumer code. It is listed here so the connection is on record, not because this PR owes it anything: `mirror` is opt-in and nothing selects it (§5.5), so nothing regresses in the meantime. |
| O2 | `phi` changes meaning — and the consumers own that | **Consumer responsibility, documented not fixed here.** This PR builds images, scales them by `J`, and deposits them through the ordinary machinery. Where the field goes afterwards — ghost filling, coarsening, filtering, tagging, plotting — belongs to whoever asked for the deposit, which is the same division `03cf59690` drew for synchronization and §4.2 draws for covered cells. What this PR owes is a clear statement of the new contract: **under `mirror`, a cut cell holds `n`, not `kappa*n`.** Put it in the `IrregularDeposition::Mirror` enumerator's Doxygen and in `Ito.rst`. Known consumers that will need attention when they opt in, listed as a courtesy rather than as work items: `DataOps::filterSmooth` (`CD_DataOps.cpp:679`) pulls cut cells toward zero and would undo the mirror; the `arithmeticAverage` sites coarsen a density with no `kappa` weighting (`conservativeAverage`, the path `ItoSolver` actually uses, is fine); tagging and plot output change at the wall. |
| O3 | More than one surface, and creases | The reflection is through *a* surface. At a triple point, a re-entrant corner, or `Aerosol`'s sphere unions the even extension is not a single reflection. Fit residual is the detector, `J = 1` the fallback. **Decide whether that is enough.** |
| O6 | Tilted-normal band threshold | The `kappa < 0.25` (CIC) / `< 0.5` (TSC) thresholds are derived for axis-aligned normals only (§1.2). Chebyshev-2 is measured for tilted normals; the thresholds are not derived. |
| O7 | Images outside the problem domain | `remap()` drops off-domain particles with a counter. Near a domain corner where the solid extends past the boundary the image leaves and its contribution vanishes. Decide whether that wants a reflect-at-domain-boundary rule. |
| O8 | Thin solids and narrow gaps | If the solid is thinner than the reflection reach, the image lands in fluid on the far side. No criterion fixes this; start with an assert on minimum resolved solid thickness. |
| O9 | The correction degrades at an EB × coarse-fine boundary | **Accepted, and documented.** Where a fine grid edge passes near the EB, two things happen: an image from a fine cut cell can land outside the fine grids and be demoted by `remap()`, and a valid fine band cell can find no cut cell with a fitted `S` (§5.4). Neither loses mass. A demoted image deposits on the coarse level and `addInvalidCoarseToFine` routes coarse-under-fine cloud mass back up to the fine level *inside the deposit* (`CD_EBCoarseFineParticleMesh.H`, item 3 of the class doc), before any averaging — so `conservativeAverage` never sees it as something to overwrite. What is left is a **bandwidth** difference: the correction arrives through a coarse-width cloud rather than a fine-width one. That is precisely the trade §4.2 already accepts, and documents, for a promoted image, and it is the ordinary price of AMR: the answer near a refinement boundary is not the answer you would get without one. §5.4's `ok` counter stays as a **diagnostic**, not an abort — it says how much of the band was affected. |
| O10 | Multi-valued cells, on the band side | **Closed by an invariant.** Multiply-cut cells are always refined away in this project's workflows, so a particle's cell is single-valued and the reflect path never has to choose between two VoFs — which it could not do, since a position alone does not determine a VoF. The invariant is a *workflow* invariant: nothing in the library enforces it (`AmrMesh::getMultiCutVofIterator` exists precisely because the CDR solvers handle these cells defensively). **So assert it: `CH_assert(!ebisbox.isMultiValued(particleIV))` on the reflect path**, where it costs nothing under `OPT=HIGH` and fires the day a geometry violates the invariant. §5.4's "skip multi-valued neighbours" in the fit stays as belt-and-braces. |
| O11 | Deposits that never reach the funnel | `ItoSolver::depositWeightNGP` and `depositGatheredNGP` (`CD_ItoSolverImplem.H:84`) build an `EBParticleMesh` per patch themselves (`:115`) and never call `AmrMesh::depositWeight`, so no mirror pass can apply to them. They are plot paths (`CD_ItoSolver.cpp:1705-1755`) and deliberately NGP, so this is correct — but it means the plotted particle density in a cut cell will not match the `phi` the field solver sees, by exactly the factor the mirror introduces. That is a diagnostic people will read as a bug. **Name it in the plot variable's documentation.** |

> **Two corrections revision 6 owes O9.** First, it wrote that a demoted image is *"overwritten by
> `conservativeAverage`"*. It is not: `addInvalidCoarseToFine` interpolates coarse-grid deposition
> clouds onto the fine level from inside the deposit, which is exactly the case it exists for. The
> defect is bandwidth, not lost mass. Second, it named `AmrMesh.buffer_size` as the knob to raise.
> That parameter is read into `m_bufferSizeBR` (`CD_AmrMesh.cpp:2688`) and reaches only the
> `BRMeshRefine` constructor in the `BergerRigoutsous` branch (`:1217-1221`); the `Tiled` branch
> (`:1229`) takes no buffer at all, and shipped inputs split **70 `br` / 70 `tiled`**. Neither
> correction changes the disposition, which is to accept it.

**Retired since revision 4:**

- **O4 (how `J` reaches `depositGathered`)** — gone. The image is `ParticleContainer<NoPayload>` and
  the gather is evaluated on the *source* particle at build time, so no gather ever runs on an
  image and there is no linearity assumption to assert (§4.2).
- **O5 (band mask width on the `Halo` path)** — demoted to a comment in `depositHaloCore`.
  `widthScale != 1` occurs at exactly one line in `EBAMRParticleMesh` (`:897`); `Interp`, `HaloNGP`
  and `Transition` all deposit at 1.0.

**Unmeasured, honestly:**

- **What the real geometry generator delivers.** The first number appears in phase 3.
- **Extending the invariants off the cut cells.** A sphere has constant curvature, so no table here
  can see the cost of borrowing `(x_c, n_c, S_c)` from the nearest cut cell. A torus can.
- **Cost, and the reflect fraction is not the way to estimate it.** One extra `EBAMRCellData`
  (`m_mirrorScratch`), one extra `ParticleContainer`, a whole-domain increment per deposit, an extra
  `remap()`, and the per-regrid fit.

  Revision 5 quoted reflect fractions of 62.7–98.0% and concluded "a vessel wall is not exotic".
  **Those are properties of the harnesses' 16³ box, not of a simulation.** Every harness here puts a
  16³ domain tightly around its geometry, so the band is a large share of the fluid volume by
  construction — for a concave sphere of `R = 8*dx` with a reach of ~2.3, the shell is
  `1 - (5.7/8)^3 = 64%` of the volume, and the harness duly reports 42–58%. In a real run the band
  is two cells of a domain hundreds of cells across, and the reflect fraction goes with
  surface-to-volume. Under the criterion of §1.2 the corrected figures are, concave:

  | `R/dx` | CIC reflect | CIC max `d` | TSC reflect | TSC max `d` |
  |---|---|---|---|---|
  | 8 | 42.1% | 2.29 | 57.9% | 3.15 |
  | 6 | 47.8% | 2.33 | 64.9% | 3.19 |
  | 4 | 65.0% | 2.13 | 82.8% | 2.96 |
  | 3 | 73.9% | 2.06 | 90.5% | 2.74 |

  The cost that does **not** shrink with the band is the fixed per-deposit work: `remap()` is a
  collective, and the increment is whole-domain. Those are O(cells) and
  O(particles) regardless of how thin the band is, and they are paid on every deposit of every
  species. **Early-out the whole mirror pass on a level with no cut cells**, and measure the
  remaining fixed cost before assuming it is small — `remap()` has been the dominant cost in this
  tree before.

  The `EBAMRCellData` of §5.1 adds ~430 MB per level per phase per realm on top of that; report it
  in the memory report rather than letting it appear as unexplained growth.
- **Tight cavities are a NaN risk, not an accuracy risk** — and the numbers are kernel-dependent.
  Revision 5 read *"down to R = 3`*dx`, where `J` reaches 179"*, measured with the retired distance
  band of §9 and therefore on a sample the scheme does not produce. Under the criterion of §1.2
  (`mirror_cavity.py`, revision 6) the criterion itself caps the reach — the image of a deep particle
  in a small cavity lands far inside the solid where its cloud touches no `kappa > 0` cell — so `J`
  reaches **28.9 at `R = 3*dx` under CIC and 483.8 under TSC**, never singular. The worst binned
  deviation is unchanged at 2.16%, at `R = 3.5`.

  So the stated mechanism is now also the harness's mechanism, which it was not before: revision 5's
  run dropped those images by a distance test long before the covered-cell test could. The
  conclusion survives; §3.2's guard has to be sized for `J ~ 500`, not for `J ~ 179`.
- **The absorbing-wall boundary condition.** Reflection imposes `dn/dnhat = 0`, right for a uniform
  plasma and for reflecting or dielectric surfaces, wrong at an absorbing electrode. Out of scope.

---

## 8. Reproducing the numbers

Run from `Prototypes/MirrorDeposition/`. Review 4 re-ran every table in revision 4. The
reflection-source ladder (§2.4), the discrete-curvature fit (§3.3), the noise sweep, the
`sum(kappa phi dV)` table (§6.1), the zeroing table and the cavity numbers all reproduced **exactly**.
The CIC column of the band table did not, for the reason in §1.2, and is corrected here.

```
python3 reach_cells.py                          # band in cells; Chebyshev-N retention
python3 img_reach.py                            # where the IMAGE's cell sits          (new in r5)
python3 mirror_frame.py                         # the tangent frame, and the frame-free form (new in r6)
python3 band_weight.py                          # mirrored mass vs distance
python3 mirror_source.py                        # the reflection-source ladder
python3 mirror_discrete_curvature.py fit        # curvature from the discrete normal
python3 mirror_discrete_curvature.py noise      # how much normal error the fit absorbs
python3 mirror_levelset_curvature.py endtoend   # why the level-set route was dropped
python3 mirror_band_kernels.py                  # CIC and TSC bands, real kernels
python3 mirror_zeroing.py                       # the zeroing-deposit failure, on a plane
python3 mirror_cavity.py                        # tight-cavity Jacobian blow-up (rewritten in r6)
python3 mirror_sphere_ext.py radii              # the curved cross (revision 2)
python3 mirror_planar_ext.py band | criteria | margin
```

**Harness bug 6, fixed in r6.** `mirror_cavity.py` selected reflecting particles with
`d <= 1.5*dx*sum|n_i|`, the `3*s_max` distance band §9 retired two revisions earlier, and reported
the analytic amplification at `d = 1.5*sqrt(3)`. It now applies §1.2's criterion — the image's cloud
must overlap a `kappa > 0` cell — and reports the max `d` that criterion actually admits, for both
kernels. Every number in §7's cavity bullet changed. The class of error is review 5's Finding 3 and
it is one level up from harness bug 5: not *"the harness used the wrong kernel"* but *"the harness
used the right kernel with a selection rule the plan had already discarded"*. When re-running
anything here, check the **selection** as well as the kernel.

**Harness bug 5, fixed in r5.** `reach_cells.py` used the TSC triangular-cloud weight for both
kernels, so its CIC mass-retention column described a kernel the code does not have. Fixed by
branching on `L`. `mirror_band_kernels.py` had the same formula but only ever tested
`w > 1e-14`, so its support-based results were unaffected; its docstring is corrected.

**One README number does not reproduce.** The margin sweep for `near lo, tilted` is quoted there as
−0.007/−0.011/−0.019/−0.025 for margins 1–4; the current harness gives
+0.00686/+0.01108/+0.02079/+0.03559 with the sample sizes (235→24) matching exactly. The conclusion
the plan draws from it — no plateau at margin 3 — is unaffected. The other four known harness bugs
are documented in `README.md`. Do not use `mirror_test.py lattice` — its zero-standard-error output
is a sampling artifact and is the most convincing wrong number in the set.

**This whole directory is throwaway and must be deleted before the branch merges.**

---

## 9. Superseded, so nobody re-derives them

| Claim | Revision | Why it died |
|---|---|---|
| Band `= dx(1 + sum\|n_i\|/2)` | 1 | Short by up to 0.71`*dx` on tilted normals; validated only on an axis-aligned normal |
| Band `= 3*s_max` / `4*s_max`, per kernel | 2, 3 | Correct but unnecessary — stated in cells it is "grow by 2", with no dependence on the normal or the kernel |
| Reflection must be per-particle from the IF | 2 | The per-cell plane's 7–18% was the cost of *borrowing* the nearest cut cell's plane, not of storing one |
| Curvature from the level set | 1, 2 | Unspecified `h`; 231% on a `Tessellation` at `h = 1e-4`; 26.1% for a non-SDF implicit function; 4.1% of concave images through the singular guard |
| `sum(kappa phi dV)/sum(m)` as primary metric | 2 | Reads 0.99974 on a build that is −46.5% wrong at `kappa <= 0.05` |
| "Does not degrade at small `kappa`" | 3 | Tabulated 3³ while shipping 5³; conflated spread with bias |
| "Positivity-preserving by construction" | 3 | The Jacobian numerator changes sign past `d = 1/c_max` |
| `bndryCentroid` differenced directly | 3 | Cell-relative units; inflates fitted curvature ~10× |
| `eb_ghost >= 2` is the whole ghost requirement | 4 | `eb_ghost` sizes the `EBLevelGrid`; the surface data is sized by `num_ghost` (`CD_PhaseRealm.cpp:189`) |
| CIC retains 95.5% at Chebyshev 1 | 4 | Measured with the TSC weight formula at `L = 1`; the real top-hat CIC gives 91.6% |
| `kappa < 0.25` is the band threshold | 4 | CIC-only. TSC's cloud half-width is 1, giving `kappa < 0.5`, and `axis, s=0.12` measures Chebyshev 2 under TSC |
| Level-preserving remap | 4 | `Transition` already deposits transition particles at fine width (`:1276`); pinning an image to the coarse level puts it outside the 0.5–1-cell CF band, where `conservativeAverage` discards it |
| "Covered cells are never written" | 4 | The kernels write every cell in `cloudBox` with no EB test (`CD_EBParticleMesh.H:725, 779`) |
| `EBAMRParticleMesh` cannot hold the image container because it "needs a realm" | 4 | `ParticleContainer::m_realm` is a label read only by `getRealm()`; the real gaps are `minBlockSize`/`levelTiles`/`validCells`, and the container belongs on `PhaseRealm` for other reasons |
| The mirror is a self-contained change to the deposit | 1–4 | `ItoSolver::depositWeight` calls `redistributeAMR` one line later (`CD_ItoSolverImplem.H:79`), on the explicit assumption that `dc = kappa*phi` |
| `S_c` stored as a 2×2 "in the tangent frame" | 5 | The fit's frame and the reflection's frame have to be the same one and nothing said so; a mismatch moves the image by up to 2.24`*dx` while `J` and the §6.1 mass check stay exactly right (`mirror_frame.py`) |
| `K = det S` | 5 | True of the 2×2; **identically zero** for the world-frame `S_c`, which is rank ≤ 2. Use `K = ½[(tr S_c)² − tr(S_c²)]`, which also removes the 2-D special case |
| "Extend `(x_c, n_c, S_c)` to every band cell from the nearest cut cell" | 5 | Total as written. `EBISBox` sees cut cells outside the level's grids that no box on that level can fit an `S` for, and `exchange()` cannot deliver one |
| Cavity `J` reaches 179 at `R = 3*dx`; reflect fractions 62.7–98.0% | 5 | Measured with the `3*s_max` distance band this table retired in revision 3. Under §1.2's criterion: `J` = 28.9 (CIC) / 483.8 (TSC), reflect 42.1–90.5% |
| Build images "before any deposit" | 5 | All three coarse-fine cores restore the valid holder before returning (`:912, :988, :1292`); the constraint is *which holder*, not *when*, and as written it contradicts §5.3 |
| 93 occurrences to convert | 5 | 36 of the 93 are the **interpolation** parameter that §5.5 puts out of scope; the deposition side is 57, and only ~30 need the enum at all |
| `parseOptions` structurally cannot see `mirror` + NGP | 5 | The cited literal call site passes `a_forceIrregNGP = false`, i.e. `Native`, and can never collide. The collision is `deposition = ngp` with `irregular_deposition = mirror`, and both are options of the same solver |
| The scratch field of §4.2 step 3 | 5 | Used by §4.2 and §5.3, declared nowhere. It is `m_mirrorScratch` in §5.1 now |
