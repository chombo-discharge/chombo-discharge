# Mirrored cut-cell deposition — implementation plan

**Issue** chombo-discharge#29 (parts 1 and 3) · **Branch** `mirror_deposition` · **Depends on** PR #700
· **Revision** 4, after three adversarial reviews.

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

**The document you are reviewing is this file.** The HTML rendering at
`scratchpad/plan_r4.html` is the same content; if they disagree, this file wins.

**What three prior reviews have already caught.** Each was a confident, plausible, wrong claim that
survived a reading of the plan and died to a measurement. Assume more remain.

| Review | What it refuted |
|---|---|
| 1 | The band formula `dx(1 + sum\|n_i\|/2)`; the "aggregate bias bounds the per-case bias" reading |
| 2 | The reflection source attribution; `sum(kappa phi dV)/sum(m)` as a primary metric; one table that did not reproduce at all |
| 3 | The band as a distance; "does not degrade at small kappa"; `bndryCentroid` units; "positivity by construction" |

**How to review this.** In descending order of yield:

1. **Re-run every table.** Revision 2 shipped a table that reproduced in one row of seven. Every
   number in this file is reproducible from this directory — see §8. Check the *stencil label* on
   curvature output: `mirror_discrete_curvature.py fit` prints a `stencil radius 1` (3³) block and a
   `stencil radius 2` (5³) block, and revision 3 quoted columns from both in one table.
2. **Check the code citations.** They are given as `file:line` throughout and were exact at
   `a00ac0d4d`. The branch moves; re-resolve them rather than trusting them.
3. **Attack the derivations in §2 and §3.** The reflection point and Jacobian are where a sign or a
   frame error is silent rather than loud.
4. **Do not re-derive widths along the normal.** §1 states every width in cells. Two earlier reviews
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

### 1.2 Band mask: Chebyshev 2, for both kernels

> Reflect only particles in cells within **Chebyshev distance 2** of a cut cell. Within the mask,
> reflect a particle if and only if its image's cloud overlaps at least one cell with `kappa > 0`,
> and `d_p > 0`.

The exact criterion falls out of the deposit loop for free and costs nothing, because the overlap
weight is already zero at the support boundary. The mask is only an early-out.

Measured against the real overlap-integral kernels, share of mirrored mass retained — counting only
the part of each image's cloud landing in a `kappa > 0` cell, since covered cells are never written
(`reach_cells.py`, worst case over ten planar cases):

| Band mask | CIC worst | TSC worst | Verdict |
|---|---|---|---|
| Chebyshev <= 2 | **100.000%** | **100.000%** | exact, both kernels, all ten cases |
| Chebyshev <= 1 | 95.5% | 87.2% | loses 4.5% / 12.8% on `axis, s = 0.45` |
| Chebyshev <= 0 | 5.1% | 5.1% | cut cells only — hopeless |

Maximum Chebyshev distance from a reflected particle's cell to the nearest cut cell is **2 for both
kernels**, even though TSC's deepest reflecting particle sits at `d = 3.39*dx` against CIC's 2.57.
The mask does not depend on the kernel.

**Why 1 is not enough, and exactly when 2 is needed.** The Chebyshev-1 reach belongs to the
*image's* cell, not the particle's, and the reflection moves the particle across the boundary by
`2d`. Take an axis-aligned EB at `x = c + t` inside cut cell `c`, so `kappa = 1 - t`. The image must
overlap cell `c`, i.e. `x_img > c - 1/2`, so

```
d_max        = t + 1/2
deepest x_p  = c + 2t + 1/2
```

which leaves cell `c+1` when `2t + 1/2 > 2`, i.e.

```
t > 0.75      <=>      kappa < 0.25        (axis-aligned)
```

**Chebyshev-1 suffices only if every cut cell has `kappa >= 0.25`** — precisely the regime the
mirror does not care about. The sliver cells it exists to fix are the ones that need 2.

Worked example at `kappa = 0.05` (EB at 8.95, cut cell 8):

| `x_p` | `d` | particle cell | Cheb(p,8) | `x_img` | image cell | Cheb(img,8) | w into cell 8 |
|---|---|---|---|---|---|---|---|
| 9.20 | 0.25 | 9 | 1 | 8.70 | 8 | 0 | 0.80 |
| 9.90 | 0.95 | 9 | 1 | 8.00 | 7 | 1 | 0.50 |
| 10.00 | 1.05 | 10 | **2** | 7.90 | 7 | 1 | 0.40 |
| 10.20 | 1.25 | 10 | **2** | 7.70 | 7 | 1 | 0.20 |
| 10.40 | 1.45 | 10 | — | 7.50 | 7 | 1 | 0.00 |

**The image is never two cells away — the particle is.** Since `x_img >= c - 1/2` and
`x_img < a < c+1`, the image's cell is always `c-1` or `c`: Chebyshev <= 1, exactly as the kernel
reach demands. Nothing reaches further than one cell. The 2 arises because the mask is indexed by
the *particle's* cell — that is where the reflect decision is taken and where `(x_c, n_c, S_c)` is
read — and reflection is symmetric in *distance* but not in *cells*. With the EB near the top face
of cell `c`, walking up leaves that cell after `1-t` and reaches `c+2` quickly, while walking down
stays inside it for `t`.

The rule predicts every axis-aligned case in the harness: `axis-aligned x` (`t = 0.37`),
`near lo edge` (0.40) and `axis, s=0.12` (0.62) measure Chebyshev 1; `axis, s=0.45` (0.95)
measures 2.

> **Open:** the `kappa < 0.25` threshold is derived for axis-aligned normals, where `kappa` and `t`
> are the same parameter. For a tilted normal they decouple. The Chebyshev-2 *result* is measured
> for tilted normals (all ten cases, both kernels); the tilted threshold is not derived.

### 1.3 Ghost cells: `eb_ghost = 2`

Three independent requirements land on 2. **Keep them separate in the code** — do not derive one
from another.

| Requirement | Needs | Forced or chosen? |
|---|---|---|
| Curvature fit stencil (5³) | 2 | **Chosen** — 5³ was picked because 7³ does not fit in `eb_ghost = 2` |
| Band mask construction | 2 | **Forced** — must see cut cells 2 out to mark valid cells within Chebyshev 2 |
| Extension of `(x_c, n_c, S_c)` | 2 | **Forced** — a band cell's nearest cut cell is at most 2 away |

Deposition itself needs only 1, because the level-preserving remap (§4.2) hands each image to its
own owning patch.

Move to a 7³ fit for noise headroom (3.7% against 5.0% at 10° of normal error) and `eb_ghost`
becomes 3 while the band stays 2.

The extension's 2 holds **only** under the ordering in §5.3 — exchange the fitted cut-cell values
into ghosts, then extend. Extend first and exchange once afterwards and a patch must fit `S` for cut
cells 2 cells outside its valid region, each of which needs its own radius-2 neighbourhood, reaching
4 cells out.

`eb_ghost` is not a default: `CD_AmrMesh.cpp:2788` is `pp.get`, validated only as `>= 0`
(`:2790`), so `eb_ghost = 1` is silently accepted today. `CD_AmrMesh.options` ships 2; shipped
inputs split 39×`=2` and 100×`=4`. **Require `>= 2` at runtime.**

---

## 2. The reflection point

### 2.1 Stored state

Three quantities per band cell, taken from the nearest cut cell, computed once at regrid:

```
x_c    boundary centroid, as a PHYSICAL position    3 Reals
n_c    normal(vof)                                  3 Reals
S_c    fitted shape operator                        3 Reals (2x2 symmetric, tangent frame)
```

> **`bndryCentroid` is cell-relative — convert it.** `EBISBox::bndryCentroid(vof)` returns
> coordinates in `[-0.5, 0.5]` *relative to the cell*, not a position. The tree's own conversion is
> `CD_EBHelmholtzEBBCImplem.H:32`:
>
> ```cpp
> x_c = probLo + (0.5*RealVect::Unit + RealVect(vof.gridIndex()) + ebisbox.bndryCentroid(vof)) * dx
> ```
>
> Revision 3 wrote the curvature fit as `dx_j = bndryCentroid(j) - x_i`, differencing two
> cell-relative values from *different* cells. That drops the inter-cell offset entirely: for cut
> cells one cell apart the true tangential baseline is `~dx` and the differenced raw centroids give
> `~0.1*dx`, so `S` comes out inflated by roughly 10×, with nothing in the output looking wrong.
> Store and use `x_c` as a physical position throughout, including in `w = p - x_c`.

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

`d` is the **first iterate** of the signed distance, not the signed distance: `xi` is the tangential
projection of the particle, not of the foot point. The neglected term is `O(c^2 d^2 xi)`, and at
`R = 4*dx` with `d` out to 2.6`*dx`, `c*d ~ 0.65` is not a small parameter. It measures at 2.9%
there, so one step is enough at the radii tested — revision 3 called it "signed distance to the
quadratic patch", which invites someone to trust it on a rod cap without re-measuring.

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

2-D:  S is 1x1, K = 0, J = (1 - c d)/(1 + c d)
```

`J` is exactly 1 on a plane, where `S = 0`. Sphere check: `J = (R-d)^2/(R+d)^2`, the ratio of shell
areas, as it must be.

**Symbols.** `kappa` is the volume fraction everywhere in this document; principal curvatures are
`c_1, c_2`. Revision 3 used `kappa` for both.

Without the Jacobian the mirror carries a curvature error of +117% at `R = 4*dx`, +40% at 8, +16%
at 16 and +4% even at 40 — with *exact* reflection, so no improvement to the surface model removes
the need for it. Reflection is measure-preserving only about a plane.

### 3.2 Positivity is conditional, not structural

Revision 3 claimed the scheme is "positivity-preserving by construction". It is not. The
**numerator changes sign** once `d > 1/c_max` on any non-umbilic convex surface — a negative-weight
image, i.e. `phi < 0`. Positivity holds for `|d| < 1/c_max` and is enforced by the guard, not by
construction.

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

No principal directions, no eigen-decomposition, no Hessian. **Use a 5³ neighbourhood** and require
at least **four** usable cut neighbours (this is what the harness enforces; revision 3 said three).
Skip multi-valued neighbours — a multi-valued cell contributes two normals to the fit.

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
> `R = 4*dx`, +2.18% against +0.17% concave. What *is* true is that the p05–p95 **spread** is widest
> at high `kappa`. Spread is not bias.
>
> So the small-`kappa` risk is resized, not retired. It is not "the fit falls apart in slivers" —
> the fit is well-conditioned there and the spread is narrow. It is "the fit carries a systematic
> 2–3% offset in slivers at tight radii". Narrow spread over 24–168 cells means it is systematic and
> will not average away. **Budget 2–3% and bin the acceptance test by `kappa`.**

Sensitivity to discrete-normal error, sphere `R = 6*dx` (`mirror_discrete_curvature.py noise`);
concave matches convex to within 0.4 points at every entry:

| normal noise | 3³ | 5³ | 7³ |
|---|---|---|---|
| 0° | 1.4% | 0.7% | 0.5% |
| 1° | 2.0% | 1.0% | 0.6% |
| 2° | 3.3% | 1.2% | 0.7% |
| 5° | 6.3% | 2.5% | 1.7% |
| 10° | 11.7% | **5.0%** | 3.7% |

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
phi(x_j) = (1/dx^D) sum_p m_p [ W(x_j - x_p) + J_p * W(x_j - R(x_p)) ]     for kappa_j > 0
```

The estimator is the kernel density estimate of the **even extension** of the density about the
surface. Where the true density is uniform, the even extension is uniform, so the estimate is exact
up to the wall. For a wall at `a = 0.5` in 1-D, plain CIC gives `phi_0 = 0.5n`; the mirror gives
exactly `n`. The apparent doubling of a particle on the boundary is the correction, not a double
count — in the physical measure it cancels: `kappa_0 * phi_0 * dx = 0.5 * 2m = m`.

**No per-particle normalization.** The issue proposes scaling each particle by
`V = sum_{kappa>0}(w + w')` so it deposits exactly `m`. That is wrong and must not be implemented.
`V` is 1 only for a face-aligned wall; for a wall at `a = 0.3` a particle at `0.4` has `V = 1.6`.
Those values are what makes the ensemble integral close, because a cut cell's value is a density
over the whole cell while its particles occupy only the fluid fraction. Dividing by a
position-dependent `V` reintroduces exactly the bias the mirror removes.

### 4.2 Four steps

```
1. build   for each VALID particle inside the band, append an image at R(x_p)
           carrying weight m_p*J_p, into a persistent image ParticleContainer

2. remap   LEVEL-PRESERVING: route each image to the patch on ITS OWN LEVEL that
           contains the image position

3. deposit the image container with the ordinary kernel and the ordinary
           coarse-fine strategy, into a scratch EBAMRCellData

4. add     DataOps::incr the scratch field into the real one
```

Step 4 is not optional: every AMR deposit path opens with `DataOps::setValue(a_meshData, 0.0)`
(`CD_EBAMRParticleMesh.H:641, 844, 943, 1217`), so step 3 needs a scratch field.

**Level-preserving remap.** `ParticleContainer::remap()` gives ownership to the patch on the
*finest* level containing the position — `LevelTiles::findDestination` loops
`lvl = finestLevel -> 0` and returns the first hit (`CD_LevelTiles.H:161-169`). For an image that is
exactly wrong: a coarse image landing under a fine grid is promoted, and a valid fine cell receives
a coarse particle's whole weight through a fine-width cloud, an `O(refRat^D)` density spike.

`findDestination` already takes the finest level to search as a parameter;
`CD_ParticleContainerImplem.H:44` simply passes `m_finestLevel`. Passing the **source** level
instead gives "the finest level at or below `lvl` whose tile owns this point" — which is also the
right rule when the image leaves `lvl`'s grids and belongs to a coarser one.

With that, containment stops being a question: an image is deposited by its own owning patch, so its
reach is 1 cell like any particle's.

### 4.3 Build images from valid particles only

`depositHaloCore:849` is `copyMaskParticles` — a **copy**. Under `CoarseFineDeposition::Halo` a
coarse halo particle is genuinely deposited twice, once on the coarse level and once onto the
refined-coarse grid at `widthScale = refRat` (`CD_EBAMRParticleMesh.H:880`); it is not double
counting only because `conservativeAverage` later overwrites the coarse cells under the fine grid.
`HaloNGP` and `Transition` move instead of copying (`:949`).

So images must be built from the **valid** holder only, before any deposit. Build them from a mask
holder as well and `Halo`'s copy generates a second image — issue #29's double-counting trap 1 in a
new guise. The image container then goes through the identical coarse-fine pass with its own mask
copy over the same masks, so it must be a **persistent container on the same realm with the halo
masks available**, not a bare scratch one.

---

## 5. Where this slots in

### 5.1 `PhaseRealm` — the surface data

`(x_c, n_c, S_c)` is per-cell mesh data with a regrid lifetime and a ghost exchange. That is exactly
`PhaseRealm`'s existing mask machinery, and `m_irregularCells` is the direct precedent:

- declared `EBAMRCellData m_irregularCells` (`CD_PhaseRealm.H:600`)
- built in `PhaseRealm::defineMasks(a_lmin, a_numGhost)` (`CD_PhaseRealm.cpp:584`) — allocate per
  level with `EBCellFactory(ebisl)`, fill on valid cells, then `exchange()` (`:657-660`)
- called from `PhaseRealm::regridBase(a_lmin)` (`:189`), after `defineEBLevelGrid` (`:163`)

Add the surface data beside it, as a 9-component `EBAMRCellData` (or three fields), built in the
same place and in the same order — **with the ordering caveat in §5.3**, which `defineMasks` does not
have to worry about because a mask is purely local per cell.

Gate it on an operator registered through `PhaseRealm::registerOperator` (`:353`) so realms that
never deposit particles do not pay for it, in the same way `s_particle_mesh` is gated
(`CD_PhaseRealm.H:60`).

### 5.2 `EBAMRParticleMesh` — the deposit path

`EBAMRParticleMesh` owns the kernels, the coarse-fine strategies and the halo masks. It is where the
image container's deposit in step 3 actually happens, and where `widthScale` lives.

**It is realm-agnostic** — it holds `Vector<RefCountedPtr<EBLevelGrid>> m_eblgs`, `m_refRat`,
`m_dx`, `m_probLo`, `m_ghost` (`CD_EBAMRParticleMesh.H:451-511`) and contains no reference to a
realm at all. So it **cannot allocate the image `ParticleContainer`**, which needs a realm, and it
cannot reach the `PhaseRealm` surface data on its own. Both have to be handed to it.

The two funnels every deposit already passes through are `AmrMesh::depositWeight` and
`AmrMesh::depositGathered` (`CD_AmrMeshImplem.H:368, 388`) — the lowest level with both a realm and
the particle-container machinery. That is where the mirror pass is invoked from; `EBAMRParticleMesh`
supplies the kernel and coarse-fine behaviour underneath it.

> **For the reviewer:** this split — data on `PhaseRealm`, kernels in `EBAMRParticleMesh`,
> orchestration in `AmrMesh` — is the one structural decision that has not been challenged by any
> review yet. The alternative (push everything into `EBAMRParticleMesh` and give it a realm) would
> make it self-contained but widens its interface considerably. Worth a look.

### 5.3 Ordering inside the regrid build

1. On each **valid** cut cell, read `normal(vof)` and `bndryCentroid(vof)`, **convert the centroid
   to a physical position**, and fit `S` over the 5³ cut-cell neighbourhood. Skip multi-valued
   neighbours. Record the fit residual.
2. Where the fit has fewer than **four** usable neighbours, or the residual exceeds the crease
   threshold, mark the cell — those band cells get `J = 1` and a counter, never a rank-deficient `S`.
3. **Exchange the fitted cut-cell values into ghosts, then extend** `(x_c, n_c, S_c)` to every band
   cell from the nearest cut cell — or equivalently extend by two 1-cell propagation sweeps with an
   exchange between them.
4. Require `eb_ghost >= 2` at runtime.

> Revision 3 extended first and exchanged once afterwards. In that order a band cell at a patch edge
> whose nearest cut cell lies one patch over is assigned a too-far cut cell, and the exchange cannot
> repair it — both patches computed from the same incomplete data and agree on the wrong answer.
> This is the one place where `defineMasks`' valid-then-exchange pattern does **not** transfer.

### 5.4 Selection and threading

`forceIrregNGP` bool becomes an `IrregularDeposition` enum: `native` / `ngp` / `mirror`. This makes
ngp-vs-mirror exclusive by construction rather than by a parse-time check someone must remember to
write. Interpolation keeps its own two-valued type — `mirror` is meaningless when gathering.

**93 occurrences across 8 files**: `CD_EBParticleMesh.H` (49), `CD_EBAMRParticleMesh.H` (18),
`CD_AmrMesh.H` (10), `CD_AmrMeshImplem.H` (10), `CD_ItoSolver.cpp`, `CD_ItoSolver.H`,
`CD_ItoSolverImplem.H` (as `m_forceIrregDepositionNGP`), and `CD_CdrSolver.cpp` (2). Rename
`irr_ngp_deposition` to `irregular_deposition` across 29 `.options`/`.inputs` occurrences in 29
files, plus `Ito.rst`. Leave `irr_ngp_interp` alone.

---

## 6. Delivery

### PR A — mechanism, one consumer

Band mask, the per-band-cell patch data at regrid, the level-preserving remap, and the mirror
utility. Wire exactly one consumer — `ItoSolver::depositParticles` — as the proof.

**Done when** the acceptance test (§6.1) passes through the wired path, the `kappa <= 0.05` bin is
reported, and every other regression is bit-for-bit unchanged.

### PR B — wire it centrally

Move the call inside `AmrMesh::depositWeight` and `depositGathered`, gated on the enum, and delete
PR A's bespoke call. **8 direct callers in 5 files**: `CD_McPhoto.cpp` ×2 (`:1287, :1318`),
`CD_TracerParticleSolverImplem.H:374`, `CD_ItoKMCStepperImplem.H` ×2 (`:4041, :6134`),
`CD_ItoSolverImplem.H` ×2 (`:71, :148`), `CD_ItoKMCGodunovStepperImplem.H:1338`.

The eighth is the awkward one: the `hasDielectrics` branch of `computeSemiImplicitRho`
(`CD_ItoKMCGodunovStepperImplem.H:1270`), depositing onto `phase::solid` with a hardcoded
`forceIrregNGP = true`. Under PR B it inherits the enum, so `mirror` would reflect across the same
embedded boundary with the fluid side swapped — pushing particles deliberately *inside* the
dielectric back out into the gas. **Decide explicitly whether `phase::solid` deposits mirror.**

Related: `CD_ItoKMCStepperImplem.H:4041` hardcodes `DepositionType::NGP` at the call site, so
"reject `mirror` with an NGP base kernel" cannot be enforced in `parseOptions`. The check has to
live where the kernel is chosen.

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
1.7% deficit at small `kappa`. Margin 3 is a defensible conservative choice, **not** a measured
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

An **images-only build** (the zeroing-deposit bug) reads −100% in every regular cell, so the planar
test does catch it; only the `kappa <= 0.05` bin is nearly blind (−1.8 to −7.7%). Revision 2 argued
the opposite.

---

## 7. Open items

| # | Item | Status |
|---|---|---|
| O1 | The mirror activates issue #29 part 4 | Today the reaction step's charge error cancels because both sides are `kappa`-weighted; with the mirror the Ito side becomes a true density and `drho = -Qe*dn*(1-kappa)` appears. Fix belongs in `ItoKMCStepper::reconcileCdrDensities`, not here. The acceptance test cannot see it, the obvious repair reintroduces a `kappa` denominator, and `floor(kappa*phi*vol)` already zeroes reactions in slivers. **Needs a reacting cut-cell test and a decision.** |
| O2 | `phi` changes meaning | `DataOps::filterSmooth` (`CD_DataOps.cpp:679`) pulls cut cells toward zero and would actively undo the mirror — reject `rho_filter`/`cond_filter` with `mirror`, or fix it. `arithmeticAverage` sites coarsen a density without `kappa` weighting. Tagging and plot output change at the wall. |
| O3 | More than one surface, and creases | The reflection is through *a* surface. At a triple point, a re-entrant corner, or `Aerosol`'s sphere unions the even extension is not a single reflection. Fit residual is the detector, `J = 1` the fallback. **Decide whether that is enough.** |
| O4 | How `J` reaches `depositGathered` | The gather lambda reads the image particle's own columns, and after the remap a side array does not survive. All five lambdas in the tree are `weight * <column>`, so scaling the image's weight is equivalent — but that is an unenforced linearity assumption on a public template. Assert it, or give the mirror path its own entry point. |
| O5 | Band mask width on the `Halo` path | The reflect decision is taken at build time at `widthScale = 1`, but coarse halo images are deposited at `widthScale = refRat`. Grow the coarse mask to `refRat + 1`, or take the decision after the strategy fixes the width. |
| O6 | Tilted-normal band threshold | `kappa < 0.25` is derived for axis-aligned normals only (§1.2). Chebyshev-2 is measured for tilted normals; the threshold is not derived. |
| O7 | Images outside the problem domain | `remap()` drops off-domain particles with a counter. Near a domain corner where the solid extends past the boundary the image leaves and its contribution vanishes. Decide whether that wants a reflect-at-domain-boundary rule. |
| O8 | Thin solids and narrow gaps | If the solid is thinner than the reflection reach, the image lands in fluid on the far side. No criterion fixes this; start with an assert on minimum resolved solid thickness. |

**Unmeasured, honestly:**

- **What the real geometry generator delivers.** The first number appears in phase 3.
- **Extending the invariants off the cut cells.** A sphere has constant curvature, so no table here
  can see the cost of borrowing `(x_c, n_c, S_c)` from the nearest cut cell. A torus can.
- **Cost.** The level-preserving remap, one extra `EBAMRCellData` and a whole-domain increment per
  deposit, and the per-regrid fit. Reflect fraction: planar per-case reaches 43.4%, convex spheres
  17.5–30.1%, and **concave runs 62.7% at R = 8, 71.3% at 6, 91.2% at 4 and 98.0% at 3**. A vessel
  wall is not exotic.
- **Tight cavities are a NaN risk, not an accuracy risk.** Down to R = 3`*dx`, where `J` reaches 179,
  the worst binned deviation stays under 2.2% — the high-`J` images land where every stencil cell is
  covered and are dropped.
- **The absorbing-wall boundary condition.** Reflection imposes `dn/dnhat = 0`, right for a uniform
  plasma and for reflecting or dielectric surfaces, wrong at an absorbing electrode. Out of scope.

---

## 8. Reproducing the numbers

Run from `Prototypes/MirrorDeposition/`. Revision 3's four tables were all re-run in review 3: the
band table, the reflection-source ladder, the `|dJ/J|` figures and the noise sweep reproduce
*exactly*; the fitted-`2H` columns did not, and are corrected in §3.3.

```
python3 reach_cells.py                          # band in cells; Chebyshev-N retention   (new in r4)
python3 band_weight.py                          # mirrored mass vs distance              (new in r4)
python3 mirror_source.py                        # the reflection-source ladder
python3 mirror_discrete_curvature.py fit        # curvature from the discrete normal
python3 mirror_discrete_curvature.py noise      # how much normal error the fit absorbs
python3 mirror_levelset_curvature.py endtoend   # why the level-set route was dropped
python3 mirror_band_kernels.py                  # CIC and TSC bands, real kernels
python3 mirror_zeroing.py                       # the zeroing-deposit failure, on a plane
python3 mirror_cavity.py                        # tight-cavity Jacobian blow-up
python3 mirror_sphere_ext.py radii              # the curved cross (revision 2)
python3 mirror_planar_ext.py band | criteria    # band and criterion tables (revision 2)
```

`reach_cells.py` and `band_weight.py` are new and **uncommitted**. Everything else, and the four
known harness bugs, is documented in `README.md`. Do not use `mirror_test.py lattice` — its
zero-standard-error output is a sampling artifact and is the most convincing wrong number in the
set.

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
