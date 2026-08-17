# Mirrored cut-cell deposition — implementation plan

**Issue** chombo-discharge#29 (parts 1 and 3) · **Branch** `mirror_deposition` · **Depends on** PR #700
· **Revision** 5, after four adversarial reviews.

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

**What four prior reviews have already caught.** Each was a confident, plausible, wrong claim that
survived a reading of the plan and died to a measurement or to the code. Assume more remain.

| Review | What it refuted |
|---|---|
| 1 | The band formula `dx(1 + sum\|n_i\|/2)`; the "aggregate bias bounds the per-case bias" reading |
| 2 | The reflection source attribution; `sum(kappa phi dV)/sum(m)` as a primary metric; one table that did not reproduce at all |
| 3 | The band as a distance; "does not degrade at small kappa"; `bndryCentroid` units; "positivity by construction" |
| 4 | The whole post-deposit redistribution step, which the plan never mentioned; `eb_ghost` as the governing ghost count; the CIC retention column, measured with the wrong kernel; the `kappa < 0.25` threshold as kernel-independent; the level-preserving remap, which was both unnecessary and harmful; "`EBAMRParticleMesh` needs a realm" |

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
cells for CIC and `5 + 1 = 6` for TSC. See §4.2.

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

Neither is a safe default. `CD_AmrMesh.cpp:2788` is `pp.get("eb_ghost", ...)` validated only as
`>= 0` (`:2790`), and `:2805` is `pp.get("num_ghost", ...)` validated only as `>= 0` (`:2808`).
`CD_AmrMesh.options` ships 2 for both; shipped inputs split 39×`eb_ghost=2` / 100×`=4` and
71×`num_ghost=2` / 68×`=3`. **Require both `>= 2` at runtime when `mirror` is selected.**

---

## 2. The reflection point

### 2.1 Stored state

Three quantities per band cell, taken from the nearest cut cell, computed once at regrid:

```
x_c    boundary centroid, as a PHYSICAL position    3 Reals
n_c    normal(vof)                                  3 Reals
S_c    fitted shape operator                        3 Reals (2x2 symmetric, tangent frame)
```

> **`bndryCentroid` is cell-relative — use the existing helper.** `EBISBox::bndryCentroid(vof)`
> returns coordinates in `[-0.5, 0.5]` *relative to the cell*, not a position. The tree already has
> the conversion, and it is not the arithmetic you should re-type:
>
> ```cpp
> x_c = probLo + Location::position(Location::Cell::Boundary, vof, ebisbox, dx);
> ```
>
> `Location::position` is `CD_LocationImplem.H:37-41` and returns
> `(gridIndex + 0.5 + bndryCentroid) * dx`; `CD_AmrMeshImplem.H:584` already calls it this way.
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

With `t1, t2` spanning the plane perpendicular to `n_c`, and `w = p - x_c`:

```
xi     = (t1.w, t2.w)          eta = n_c.w

d      = eta + (1/2) xi^T S xi                     one-step distance to the quadratic patch
nhat   = normalize( n_c + (S xi)_1 t1 + (S xi)_2 t2 )   its normal at the foot point
R(p)   = p - 2 d nhat
```

One dependent load per particle and about twenty flops. **No implicit-function evaluation
anywhere**, at regrid or per particle.

Review 4 checked the signs independently against both curvature signs, with `S` in the fit's own
convention (`S = +dn̂`, §3.3): a convex sphere of radius `R` with fluid outside gives `S = +I/R` and
patch height `h(xi) = -|xi|²/2R`, a concave cavity gives `S = -I/R` and `h(xi) = +|xi|²/2R`, and
`d = eta - h` reproduces the expression above in both. `nhat = n_c + S xi` likewise.

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
J   = (1 - 2H*d + K*d^2) / (1 + 2H*d + K*d^2)        2H = tr S,  K = det S
    = prod(1 - c_i d) / prod(1 + c_i d)              over principal curvatures c_i

2-D:  S is 1x1, set K = 0 explicitly, J = (1 - c d)/(1 + c d)
```

`J` is exactly 1 on a plane, where `S = 0`. Sphere check: `J = (R-d)^2/(R+d)^2`, the ratio of shell
areas, as it must be.

The 2-D line is an instruction, not a restatement: `det` of a 1×1 `S` is `c`, not 0, so
`K = det S` must be overridden rather than evaluated. The product form has one factor in 2-D and two
in 3-D.

**Symbols.** `kappa` is the volume fraction everywhere in this document; principal curvatures are
`c_1, c_2`. Revision 3 used `kappa` for both.

Without the Jacobian the mirror carries a curvature error of +117% at `R = 4*dx`, +40% at 8, +15.3%
at 16 and +4.2% even at 40 — with *exact* reflection, so no improvement to the surface model removes
the need for it. Reflection is measure-preserving only about a plane.

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
    2H     = tr S        K = det S
```

This makes `S` the differential of the normal, `S = +dn̂`, which is the convention §2.2 and §3.1 are
written in. No principal directions, no eigen-decomposition, no Hessian. **Use a 5³ neighbourhood**
and require at least **four** usable cut neighbours (this is what the harness enforces; revision 3
said three). Skip multi-valued neighbours — a multi-valued cell contributes two normals to the fit.

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

**Mass that lands in covered cells is discarded, but not by the kernel.** The estimator is only ever
read at `kappa_j > 0`, so dropping the part of an image's cloud that falls in covered cells is
correct. The deposition kernels do **not** drop it: `depositParticleCIC` and `depositParticleTSC`
write `rho(iv, ...)` for every cell in `cloudBox` with no EB test (`CD_EBParticleMesh.H:725-737,
779-798`). That is pre-existing behaviour — today's fluid-side clouds already spill into covered
cells — but the mirror puts a far larger share of mass there, and anything that reads regular data
uniformly across covered cells (plot output, `arithmeticAverage` in `EBCoarseAve`,
`DataOps::filterSmooth`) would see it. Hence step 5 below.

### 4.2 Five steps

```
1. build     for each VALID particle inside the band, append an image at R(x_p) into a
             persistent ParticleContainer<NoPayload>, carrying weight g_p*J_p

2. remap     images.remap() -- the ordinary remap, no special level rule

3. deposit   the image container through the ORDINARY EBAMRParticleMesh::depositWeight,
             with the ordinary kernel and coarse-fine strategy, into a scratch EBAMRCellData

4. add       DataOps::incr the scratch field into the real one

5. clean     DataOps::setCoveredValue(real, getCoveredCells(realm, phase), 0.0)
```

Step 4 is not optional: every AMR deposit path opens with `DataOps::setValue(a_meshData, 0.0)`
(`CD_EBAMRParticleMesh.H:641, 844, 943, 1217`), so step 3 needs a scratch field.

Step 5 must run **after** the exchange and `addFineGhostsToCoarse` inside step 3 have folded ghost
mass back to its owner, and **before** redistribution and `coarsenAndFillGhosts`. `DataOps::setCoveredValue`
(`CD_DataOps.H:1200`) with `AmrMesh::getCoveredCells(realm, phase)` is the existing pair; a covered
cell is covered in every patch that sees it, so nothing real is destroyed.

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
gets `refRat`, a bandwidth mismatch across the interface. `Halo` is not the target strategy;
**document this in `depositHaloCore`'s Doxygen rather than engineering around it.** The same
disposition applies to the band-mask width on that path (former O5): `widthScale` is 1.0 at every
other patch deposit in the class (`:658, :866, :966, :967, :1235, :1276`), so the reflect decision
taken at build time matches the deposit everywhere except `Halo`.

### 4.3 Build images from valid particles only

`depositHaloCore:849` is `copyMaskParticles` — a **copy**. Under `CoarseFineDeposition::Halo` a
coarse halo particle is genuinely deposited twice, once on the coarse level and once onto the
refined-coarse grid at `widthScale = refRat` (`CD_EBAMRParticleMesh.H:880`); it is not double
counting only because `conservativeAverage` later overwrites the coarse cells under the fine grid.
`HaloNGP` (`:949`) and `Transition` (`:1220`) move instead of copying.

So images must be built from the **valid** holder only, before any deposit. Build them from a mask
holder as well and `Halo`'s copy generates a second image — issue #29's double-counting trap 1 in a
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

Four members are added, and the file's own convention decides which are `mutable`. `PhaseRealm`
declares every operator and per-call scratch `mutable` (seventeen of them, `:605-695`, including
`m_particleMesh` at `:675`, `m_redistributionOp` at `:665` and `m_nonConservativeDivergence` at
`:695`) while regrid-built data such as `m_irregularCells` is not:

```cpp
EBAMRCellData                        m_surfaceData;    // (x_c, n_c, S_c), 9 comps or three fields
mutable ParticleContainer<NoPayload> m_mirrorImages;
mutable EBAMRIVData                  m_massDiff;       // moved off the solvers, see 5.4
mutable EBAMRIVData                  m_depositionNC;
```

`m_surfaceData` is built in `defineMasks` alongside the existing masks — **with the ordering caveat
in §5.4**, which `defineMasks` does not have to worry about because a mask is purely local per cell.
`m_massDiff` and `m_depositionNC` need no new plumbing: `PhaseRealm` already owns the operators that
consume them (`m_redistributionOp`, `m_nonConservativeDivergence`) behind the `s_eb_redist` and
`s_noncons_div` registrations (`:56-57`).

`m_mirrorImages` needs four things `PhaseRealm` does not currently hold — `minBlockSize`,
`levelTiles`, `validCells` and the realm name (`ParticleContainer::define`,
`CD_ParticleContainer.H:193-202`). All four live on `Realm` (`CD_Realm.H:643, 653`) or `AmrMesh`,
are phase-independent, and `Realm` owns the `PhaseRealm`s, so this is one hop of plumbing.
`AmrMesh::allocate(ParticleContainer<P,Traits>&, realm)` (`CD_AmrMeshImplem.H:238-261`) shows exactly
which getters to forward.

Gate all four on an operator registered through `PhaseRealm::registerOperator` (`:353`) so realms
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
  ├ if mode == Mirror:
  │   ├ build images into PhaseRealm's m_mirrorImages from m_surfaceData
  │   ├ m_mirrorImages.remap()
  │   ├ particleMesh.depositWeight(scratch, m_mirrorImages, ...)
  │   ├ DataOps::incr(real, scratch, 1.0)
  │   └ DataOps::setCoveredValue(real, getCoveredCells(realm, phase), 0.0)
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

Two consequences to decide, not to discover:

- **McPhoto gains a mode** it cannot currently express. That is a behaviour addition.
- **`mirror` collides with the dominant shipped configuration.** `CD_ItoSolver.options` ships
  `false/false`, but 23 of the 30 shipped inputs that set `redistribute` set it `true`. Every one of
  those files needs a per-case decision, not a mechanical substitution.

**Do not sweep in the fluid solvers.** `CdrGodunov.blend_conservation` (7 files) and
`CdrCTU.blend_conservation` (34 files) are the *fluid* hybrid divergence and have nothing to do with
particle deposition, despite the identical option name.

**One check the enum does not subsume.** `IrregularDeposition::NGP` and `DepositionType::NGP` are
different things. `Mirror` with an NGP base kernel must still be rejected where the kernel is
chosen — `AmrMesh::depositWeight` or `EBParticleMesh::depositCore`'s switch — because
`CD_ItoKMCStepperImplem.H:4041` passes `DepositionType::NGP` as a literal at the call site, so
`parseOptions` structurally cannot see it.

**The rename.** `irr_ngp_deposition` → `irregular_deposition` across 29 `.options`/`.inputs`
occurrences in 29 files, plus `Ito.rst` and the `pp.get` at `CD_ItoSolver.cpp:320`. It is a
value-domain change, not only a key rename: `true/false` becomes one of five strings. Because the
tree uses `pp.get` and never `query`, a stale input file hard-errors on the missing key — that is the
failure mode we want, so do not add a fallback. Leave `irr_ngp_interp` alone; interpolation keeps its
own two-valued type, and `mirror` is meaningless when gathering.

**93 C++ occurrences across 8 files** to convert: 85 `a_forceIrregNGP` + 4 `forceIrregNGP` +
4 `m_forceIrregDepositionNGP`, in `CD_EBParticleMesh.H` (49), `CD_EBAMRParticleMesh.H` (18),
`CD_AmrMesh.H` (10), `CD_AmrMeshImplem.H` (10), `CD_ItoSolver.cpp`, `CD_ItoSolver.H`,
`CD_ItoSolverImplem.H`, `CD_CdrSolver.cpp` (2). The 7 `m_forceIrregInterpolationNGP` are not in scope.

---

## 6. Delivery

### PR A — mechanism, one consumer

Band mask, the per-band-cell patch data at regrid, the image container, the mirror utility, and the
covered-cell reset. Wire exactly one consumer — `ItoSolver::depositParticles` — as the proof.

**PR A must also deal with redistribution**, because it is one line below the deposit it wires.
`ItoSolver::depositWeight` is `AmrMesh::depositWeight` followed by `this->redistributeAMR(a_phi)`
(`CD_ItoSolverImplem.H:71, :79`), and `depositGathered` does the same (`:148, :157`).
`redistributeAMR` (`CD_ItoSolver.cpp:2133`) opens with

> *"When we entered this routine we had `a_phi = m_i/dV` but we actually want to have
> `phi = m_i/(kappa*dV)`"*

which is exactly the assumption the mirror removes, and `depositHybrid` (`:2193`) computes
`deltaM = (1-kappa)*(dc - kappa*dnc)` under the comment *"Remember, `dc` already scaled by kappa"*
(`:2225-2233`). With the mirror, `dc` is already a density, so `deltaM` is inflated by `1/kappa`
— 20× at `kappa = 0.05` — and smooshed into the neighbouring cells. The minimum for PR A is a hard
error when `mirror` is selected with `redistribute` or `blend_conservation`; the enum of §5.5 makes
that unrepresentable instead, which is why PR B should not lag far behind.

**Done when** the acceptance test (§6.1) passes through the wired path, the `kappa <= 0.05` bin is
reported, and every other regression is bit-for-bit unchanged.

### PR B — wire it centrally

Introduce `IrregularDeposition`, move `m_massDiff`/`m_depositionNC` onto `PhaseRealm`, move the
hybrid-and-redistribute step into `AmrMesh::depositWeight`/`depositGathered` beside the mirror pass,
and delete PR A's bespoke call and the duplicated redistribution code in the solvers
(`CD_ItoSolver.cpp:2133/2178/2193` and `CD_McPhoto.cpp:1325/1365/1381` are near-verbatim copies).

**8 direct callers in 5 files**: `CD_McPhoto.cpp` ×2 (`:1287, :1318`),
`CD_TracerParticleSolverImplem.H:374`, `CD_ItoKMCStepperImplem.H` ×2 (`:4041, :6134`),
`CD_ItoSolverImplem.H` ×2 (`:71, :148`), `CD_ItoKMCGodunovStepperImplem.H:1338`.

The eighth is the awkward one: the `hasDielectrics` branch of `computeSemiImplicitRho`
(`CD_ItoKMCGodunovStepperImplem.H:1270`), depositing onto `phase::solid` with a hardcoded
`forceIrregNGP = true` (`:1338`). Under PR B it inherits the enum, so `mirror` would reflect across
the same embedded boundary with the fluid side swapped — pushing particles deliberately *inside* the
dielectric back out into the gas. **Decide explicitly whether `phase::solid` deposits mirror.**

**Why B is not optional.** If the mirror pass stays the consumer's responsibility,
`irregular_deposition = mirror` becomes a setting that some deposits honour and others silently
ignore — character-for-character the defect PR #700 just fixed. No assert can catch it, because a
deposit function cannot know whether its caller intends a mirror pass afterwards.

Retire the nine ad-hoc `EBParticleMesh` constructions along the way — `McPhoto.cpp:1277,1459`,
`CdrSolver.cpp:1317`, `ItoSolver.cpp:2108,2528,2651,2686,2804`, `ItoSolverImplem.H:115` — by
registering `s_particle_mesh` on those realms and using the persistent per-patch leaves.

### 6.1 Acceptance

Three cases, mandatory, same PR. **Primary metric is the density binned by `kappa`.**

- **A curved case** — on a plane `J = 1` identically, so a planar-only test green-lights an
  implementation with no Jacobian and a 117% error at `R = 4*dx`. A torus is the cheapest surface
  that also exercises unequal principal curvatures and non-constant curvature, and it exercises the
  cut-cell-to-band-cell extension that a sphere cannot.
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
| O1 | The mirror activates issue #29 part 4 | Today the reaction step's charge error cancels because both sides are `kappa`-weighted; with the mirror the Ito side becomes a true density and `drho = -Qe*dn*(1-kappa)` appears. Fix belongs in `ItoKMCStepper::reconcileCdrDensities`, not here. The acceptance test cannot see it, the obvious repair reintroduces a `kappa` denominator, and `floor(kappa*phi*vol)` already zeroes reactions in slivers. **Needs a reacting cut-cell test and a decision.** |
| O2 | `phi` changes meaning — and this list is not closed | Redistribution is handled (§5.5 makes it unrepresentable alongside `mirror`), but it was missed by three reviews, so assume peers remain. Known: `DataOps::filterSmooth` (`CD_DataOps.cpp:679`) pulls cut cells toward zero and would actively undo the mirror — reject `rho_filter`/`cond_filter` with `mirror`, or fix it. `arithmeticAverage` sites coarsen a density without `kappa` weighting. Tagging and plot output change at the wall. **Audit every reader of `phi` in a cut cell.** |
| O3 | More than one surface, and creases | The reflection is through *a* surface. At a triple point, a re-entrant corner, or `Aerosol`'s sphere unions the even extension is not a single reflection. Fit residual is the detector, `J = 1` the fallback. **Decide whether that is enough.** |
| O6 | Tilted-normal band threshold | The `kappa < 0.25` (CIC) / `< 0.5` (TSC) thresholds are derived for axis-aligned normals only (§1.2). Chebyshev-2 is measured for tilted normals; the thresholds are not derived. |
| O7 | Images outside the problem domain | `remap()` drops off-domain particles with a counter. Near a domain corner where the solid extends past the boundary the image leaves and its contribution vanishes. Decide whether that wants a reflect-at-domain-boundary rule. |
| O8 | Thin solids and narrow gaps | If the solid is thinner than the reflection reach, the image lands in fluid on the far side. No criterion fixes this; start with an assert on minimum resolved solid thickness. |
| O9 | Fine grids must cover the band | When a fine patch edge cuts through the mirror band, images from fine cut cells land outside the fine grids, are demoted by `remap()`, deposit at coarse width, and are then overwritten by `conservativeAverage`. Those fine cut cells lose their correction. This is a *grid* requirement, not a remap rule: the refined region must cover the EB plus ~3 cells. `AmrMesh.buffer_size` ships 2. **Decide whether to require more, or to accept and document.** |

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
- **Cost.** One extra `EBAMRCellData`, one extra `ParticleContainer`, a whole-domain increment and a
  covered-cell reset per deposit, an extra `remap()`, and the per-regrid fit. Reflect fraction:
  planar per-case reaches 43.4%, convex spheres 17.5–30.1%, and **concave runs 62.7% at R = 8, 71.3%
  at 6, 91.2% at 4 and 98.0% at 3**. A vessel wall is not exotic.
- **Tight cavities are a NaN risk, not an accuracy risk.** Down to R = 3`*dx`, where `J` reaches 179,
  the worst binned deviation stays under 2.2% — the high-`J` images land where every stencil cell is
  covered and are dropped.
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
python3 band_weight.py                          # mirrored mass vs distance
python3 mirror_source.py                        # the reflection-source ladder
python3 mirror_discrete_curvature.py fit        # curvature from the discrete normal
python3 mirror_discrete_curvature.py noise      # how much normal error the fit absorbs
python3 mirror_levelset_curvature.py endtoend   # why the level-set route was dropped
python3 mirror_band_kernels.py                  # CIC and TSC bands, real kernels
python3 mirror_zeroing.py                       # the zeroing-deposit failure, on a plane
python3 mirror_cavity.py                        # tight-cavity Jacobian blow-up
python3 mirror_sphere_ext.py radii              # the curved cross (revision 2)
python3 mirror_planar_ext.py band | criteria | margin
```

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
