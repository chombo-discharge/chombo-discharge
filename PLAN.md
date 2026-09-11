# Cut-cell refinement — implementation plan

Refining an arbitrary single-valued cut cell, so curved geometry can carry a fine resolution
while near-planar features stay coarse, without refining the whole EB surface.

**Status.** The design is closed. This file covers **Milestone 1**: replacing Chombo's
moment machinery with ours, behind `Driver.geometry_generation = polyhedral`, with **no changes
to Chombo**. Milestones 2 (refinement) and 3 (AMR-aware generation) are in `PLAN_MILESTONES.md`. Over 48,798 Chombo-generated cells on nine geometries the
construction agrees with Chombo on every classification and every face topology and refuses
nothing; refinement to depth 2 conserves all seven moments to machine precision with zero
topology defects over 230,400 descendants. A C++ spike does the same at 2.7 µs/cell to build
and 0.9 µs per output cell to refine. Nothing has been written into Chombo yet.

---

## 1. The acceptance criteria

These are the two conditions the design is held to, and they are deliberately different from
each other:

1. **Same-level agreement.** When cells are generated on the same level, our classification
   and face topology must match Chombo's. We build from intersections and chords, Chombo does
   its own thing, and the *moments* will legitimately differ — but the topology may not.
2. **Conservation and topology through refinement.** When we refine *our* cell, every
   coarsening relation must hold and every child must be a legal cell.

Comparing our refinement against Chombo's independently generated fine level is **not** a
criterion. Chombo would be doing sub-parent-cell arithmetic there that we deliberately are
not, so a difference is expected and carries no information.

## 2. Why the obvious approaches fail

Chombo stores moments, never a surface. Refining therefore means *reconstructing* geometry
from an over-determined moment set, and every reconstruction carried its own tail of failures:

- **A plane fitted to the stored moments.** The parent's seven moment families are not the
  moments of any single plane, so children need repair fills to make each coarsening relation
  hold. Measured on curved geometry: 3.7–12 % of children came out *covered but with open
  faces* — a graph inconsistency, since a covered cell has no VoF for a face to attach to.
  Face centroids missed by up to 0.71 of a cell.
- **A quadric.** Right order (third, matching Chombo's apertures), but only one curvature sign
  can be made safe. A convex-fluid patch never produces a multi-valued child at any depth; the
  opposite sign always can, and the guard that would suppress it is exactly the guard
  refinement destroys, because refinement hunts down the point where the surface is tangent to
  a cell face.
- **A least-squares system over the child moments.** 164 unknowns against 57 exact
  constraints leaves 107 free dimensions, and realisability — "these numbers are the moments of
  an actual body" — is not a convex constraint, so the solution is free to land outside it.
  This is the fit-and-repair failure above, restated in a larger basis.
- **Re-evaluating the implicit function inside the cell.** Loses conservation immediately (the
  true child moments do not sum to the parent's stored ones) and can disagree with the coarse
  graph about connectivity, which is worse than a moment error.

## 3. The design

Store the surface explicitly and define the moments as its moments. Refinement then stops
being a reconstruction and becomes exact clipping.

### 3.1 What is stored

Per cut cell, in the cell's own `[-0.5, 0.5]^3` frame:

| | |
| --- | --- |
| **Edge crossings** | for each of the 12 edges, the parameter where the implicit function changes sign. Computed *on the edge*, which four cells share, so all four agree bit for bit. |
| **Corner values** | the implicit function at the 8 corners. Shared by the eight cells at that corner. |
| **Face bends** | up to 2 interior vertices per face, turning the face contour from a chord into a polyline. |
| **Apex** | one point on the interface, the fan centre for the interior triangulation. |

Everything else is derived. The first two are what make neighbouring cells agree by
construction rather than by tolerance.

### 3.2 What is derived

- **Face contour.** Walk each face's circuit of four corners and four edges, keeping fluid
  corners and inserting the crossing wherever an edge changes sign — the face's own marching
  squares. Bends are spliced into the chord between the two crossings they join.
- **Interface patch.** Each crossing lies on exactly two faces and each face pairs its
  crossings up, so every node has degree two and the graph is a disjoint union of cycles.
  Several cycles is legal: it just means the interface enters as several sheets. Triangulate
  each loop, fanning from the apex, falling back to the least-area triangulation where a fan
  would have to cut across a crease.
- **The body.** The six face polygons plus the patch triangles, oriented outward. Watertight
  within the cell because the patch's boundary is exactly the set of chords, and watertight
  across cells because the chords are shared.

### 3.3 The moments

| moment | computed as |
| --- | --- |
| `volFrac`, `volCentroid` | signed tetrahedral decomposition about the cell centre |
| `areaFrac`, `faceCentroid` | signed area and centroid of the face polygon |
| `bndryArea` | magnitude of the patch's area vector |
| `normal` | `-A/|A|`, Chombo's sign, pointing into the fluid |
| `bndryCentroid` | true-area-weighted centroid of the patch triangles |
| `bndryAreaTrue` | sum of the patch triangles' areas — new, see §7.1 |

`kappa = (1/3)[½ sum(alpha_f) - a_B (x_B · n)]` holds to 1e-16 by construction. Chombo's own
stored moments violate it by ~5e-3 on curved geometry.

### 3.4 Refinement

Clip the stored body by the three midplanes: Sutherland–Hodgman per polygon, then assemble
the edges left lying in the plane into the cap polygons that close each child. The eight
children exactly partition the parent, so **every coarsening relation is an identity about a
partition, not a condition to satisfy**. There are no fits, no fills and no repairs.

### 3.5 What refinement does not preserve

One chord per face is **not** an invariant — up to 20 % of child faces carry a polyline of two
to four chords, because the parent's patch is multi-triangle and its trace on a child face
bends. Two consequences for the implementation:

1. Child face polygons are often **non-convex**. All area and centroid computations must use
   *signed* areas about the polygon normal.
2. A child cannot be re-encoded as twelve edge crossings. The polygons must be carried.

The class closed under refinement is piecewise-linear, not one-chord-per-face.

## 4. Single-valued by construction

The finest level is single-valued by construction. This is permanent and it is not a
concession — it is what `GeometryShop` already does: it pushes exactly one `IrregNode` per
irregular cell with `m_cellIndex = 0` and never tests fluid connectivity.

So when a cell's fluid arrives in several disjoint pieces, it is **one VoF holding their sum**,
not a refusal and not a split. The sum is still a partition of the parent, so every coarsening
relation is untouched; and both cells adjoining a shared face see the same total aperture, so
the graph stays consistent. Multi-valued *parents* are a non-issue for the same reason: Chombo
does not generate them on the finest level, and the operator only ever runs there, since
coarser levels come from the existing coarsening chain.

### 4.1 What the collapse costs, and what refinement does to it

Collapsing several sheets — or one folded sheet — into a single `(a_B, n, x_B)` is lossy in a
specific way. `a_B = |∫dA|` under-counts whenever the patch folds or two sheets oppose; the
normal points somewhere neither sheet does; and the boundary centroid can land inside solid.
None of that breaks conservation, because the children are a literal subdivision of the parent's
triangles and the area *vectors* add whether they belong to one sheet or five. Measured: the 4
multi-piece cells among 230,400 descendants conserve at 1e-16 like everything else.

What it does mean is that the coarse-level summary is an approximation with a **bounded, and
monotonically improving, error**. Since the children partition the triangles, the triangle
inequality gives

```
a_B(parent)  <=  sum_children a_B(child) / 2^(D-1)  <=  a_B_true(parent)
```

with the same `1/2^(D-1)` factor `coarsenBoundaryAreaAndNormal` already applies. Refinement can
only move the area estimate up, toward the truth. Measured on the worst cell of each geometry:

| geometry | `a_B` parent | children, L1 | children, L2 | `a_B` true |
| --- | ---: | ---: | ---: | ---: |
| oblique plane | 1.211984 | 1.211984 | 1.211984 | 1.211984 |
| sphere | 0.025224 | 0.025312 | 0.025332 | 0.032753 |
| torus | 0.241749 | 0.241887 | 0.241939 | 0.244204 |
| swept sphere | 1.105688 | 1.324218 | 1.423017 | 1.598398 |
| rotated cube | 0.498767 | 0.609508 | 0.694278 | 0.927938 |
| rotated simplex | 0.421312 | 0.570440 | 0.699767 | 1.110899 |
| axis-aligned simplex | 0.040178 | **0.120533** | 0.120533 | 0.120533 |

Every row is squeezed between the two bounds and every row moves the right way. Three regimes
are visible:

- **Planar.** Already at the truth; refinement changes nothing.
- **Disconnected sheets.** The axis-aligned simplex reaches the exact true area in *one* level —
  a sharp edge splitting into children that each see a single flat facet. The collapse resolves
  itself.
- **A convex edge.** The rotated cube and simplex climb steadily without arriving, because the
  edge stays inside a cell at every level. These converge only as the fold within each child
  shrinks, so the collapsed normal there is a permanent compromise, not a transient one.

The ratio `a_B / a_B_true` is therefore a free per-cell detector for exactly this: 1.0 means a
single planar sheet, and anything below says how much the summary is losing. It needs no new
data — it is the same pair of quantities §7.1 already requires — and it is the sharpest
refinement criterion available (§ `PLAN_MILESTONES.md`, Milestone 3).

## 5. The rules that were expensive to find

Each of these was silent, and none is obvious from the requirement. They are the actual content
of the prototype.

### 5.1 The ambiguous four-crossing face

A face whose corners alternate fluid/solid/fluid/solid has four crossings and two ways to join
them — the textbook marching-cubes ambiguity. The **asymptotic decider** settles it from the
four corner values alone: the bilinear saddle `(f00·f11 − f10·f01)/(f00+f11−f10−f01)`. It needs
no new stored data, and both cells adjoining the face necessarily agree.

The rule is: **the saddle names the diagonal that meets through the middle; the *other*
diagonal is the pair the two chords cut off** — whether that is the fluid pair or the solid
pair. Assuming the connected pair is always the fluid one is wrong on half of all such faces.

And the face polygon must then be derived **from the pairing**, not re-read off the corner
signs. Reading it off the signs yields a geometrically correct face that silently disagrees
with the loop assembly, which does trust the pairing. Two wrongs that hide each other: the
apertures look right, the body will not close, and the residual is O(1) rather than a
tolerance issue. This was the last bug, and it accounted for a third of all refusals.

### 5.2 Classification: zeros do not make a cut

A facet lying exactly in a cell-face plane leaves every corner at or below zero with one whole
face at exactly zero. That cell is entirely fluid with that face **covered** — the cell beyond
it is solid — and it is not a cut cell at all. Reading those zeros as solid invents crossings
and reports the covered face as open: an aperture error of exactly 1.0, on 277 of 300 cells.

A cut cell needs corners **strictly** on both sides. Transcribe `insideOutsideFromNodes`
directly, including its use of `copysign`, so that negative zero counts as negative. CAD
geometry is full of axis-aligned facets; this is not exotic.

### 5.3 One fluid predicate, and it must be `copysign`

`-0.0 < 0.0` is false, and signed-distance functions do return negative zero. Corner
classification, crossing detection, and the saddle test must all use the *same* predicate, and
that predicate is `copysign(1.0, v) < 0.0`, because that is what Chombo uses. Negative zero
produced four distinct bugs in this work.

### 5.4 Small-but-not-zero decides something qualitative

The dominant failure family, stated once: a quantity that is small but not exactly zero passes
a `> 0` test and then decides something *qualitative* — a weighted-mean centroid, a connected
component, a face classification. Concretely:

- Apertures of 1e-28 dictating an entire weighted-mean centroid (three separate sites).
- Sliver triangles of ~1e-23 area left by holding crossings off edge endpoints; clipping
  detaches one and it reads as a second connected component. 125 cells were reported
  multi-valued on the strength of a speck of volume 1e-24.
- **Cull thresholds must match the dimension of the quantity.** The same sliver is an *area* of
  order 1e-24 in 3-D but a *length* of order 1e-12 in 2-D. One shared constant caught it in 3-D
  and sailed straight past it in 2-D.

Guard every magnitude comparison with a threshold, derive the threshold from the quantity's
dimension, and force covered and full faces to exactly 0 and 1.

### 5.5 And the mirror error: a threshold that is too coarse

A cell is covered when it has *no* volume, not when its volume is below a tolerance. Volume
scales as L³ and area as L², so a legitimate corner sliver of leg 1e-4 has `kappa = 1.7e-13`
with apertures of 5e-9, and a flat 1e-12 threshold applied to both declares it "covered but
with open faces". **Every** reported instance of that defect turned out to be this. Test
topology with exact zeros; use tolerances only on quantities of the same dimension.

### 5.6 The remaining robustness rules

1. **Hold crossings off edge endpoints** (1e-12). A node exactly on the interface makes the
   surface tangent to the faces meeting there and the crossing lands on a corner — a contour
   touching a face without cutting any of its edges, which marching squares cannot represent.
2. **Orient the patch by crossing identity, not position.** `np.allclose` default tolerances
   matched the wrong vertex once crossings sat a hair off a corner, silently reversing it.
3. **Snap kept vertices onto the clip plane.** A vertex retained within tolerance must be
   placed *exactly* on the plane, or its cap is no longer recognised as lying in a cell face
   and its area is charged to the interface.
4. **Assemble caps by consuming segments, not vertices.** A cross-section can have several
   components and two may meet at a shared vertex; marking that vertex used drops one loop and
   splices the other.
5. **Accumulate face area signed.** Where the interface lies in a cell face it bounds a *hole*
   in that face. Summing magnitudes gave apertures of 1.14.
6. **Interpolate edges canonically** (C++). `a + t(b−a)` and `b + t′(a−b)` differ in the last
   bits and land in different quantisation buckets, so the two cells sharing an edge disagree.
   Order the endpoints by a `lexLess` before interpolating.
7. **Do not hash vertex keys.** Collisions splice cap loops. Exact componentwise comparison took
   conservation from 5.0e-01 to 4.4e-16.

## 6. Evidence

### 6.1 Same level, against Chombo — criterion 1

Every cut cell Chombo generated on the finest level of nine geometries: five analytic (oblique
plane, two spheres, torus, swept sphere along a body diagonal) and four piecewise-planar
(rotated and axis-aligned cube, rotated and axis-aligned simplex — an STL is a piecewise-planar
surface, and a rotated polyhedron tests it far more sharply than a real tessellation because it
puts edges and corners at chosen angles rather than leaving it to chance).

| | cells | class mismatches | face-topology errors | refused |
| --- | ---: | ---: | ---: | ---: |
| all nine geometries | 48,798 | **0** | **0** (290,970 faces) | **0** |

### 6.2 Refinement — criterion 2

400 parents per geometry, refined to depth 2: 230,400 descendants. Worst error over all nine
geometries, both per level and accumulated at the root:

| moment | worst | typical |
| --- | ---: | ---: |
| `volFrac` | 9.4e-16 | 4e-16 |
| `volCentroid` | 3.2e-05 | 1e-13 |
| `areaFrac` | 3.3e-16 | 2e-16 |
| `faceCentroid` | 7.6e-16 | 2e-16 |
| `bndryArea` | 6.7e-16 | 4e-16 |
| `normal` | 1.2e-06 deg | 1e-06 |
| `bndryCentroid` | 7.2e-16 | 3e-16 |

Topology over all 230,400: zero `kappa` or `alpha` out of range, zero centroids outside their
own cell, zero isolated cells, zero covered-cells-with-faces, zero refusals. One dust cell
(`kappa` 6.2e-16, below the 1e-15 `EBISLevel` already drops) and four multi-piece cells carried
as one VoF per §4.

There is **no error growth with depth**: the accumulated error equals the single-level error,
which is what transitivity of clipping predicts. `volCentroid` outliers around 1e-5…1e-8 are
conditioning, not conservation — the centroid is a ratio and cells with `kappa` near zero
amplify rounding. It does not propagate, but solvers should not assume machine precision there.

Cross-cell consistency: fine-face apertures and face centroids agree to 1e-16 computed from
either side.

### 6.3 C++ spike

`spike/refine_spike.cpp`, storage = 12 crossings + 8 corners, body rebuilt on demand.

| | |
| --- | --- |
| refusals, all nine geometries | **0** |
| conservation | ≤ 8.3e-16 |
| build | 2.7 µs/cell (from 17.1) |
| refine | 0.9 µs/output cell (from 1.34) |

100k cells built and refined one level: 0.62 s, against a 1 s budget.

### 6.4 Two dimensions

The construction degenerates cleanly: the patch is a single chord rather than a triangle fan,
and clipping a polygon by a half-plane closes itself, so there are no cap loops. Nine
geometries to three levels, ~43,000 cells: all seven moments at machine precision, zero
refusals, zero topology defects.

**In 2-D there is no order loss at all** — the construction reproduces Chombo's own stored
moments to 1e-13, because in 2-D the interface inside a cut cell genuinely *is* a straight
chord. The third-to-second-order trade below is a 3-D-only cost.

### 6.5 Against Chombo's stored moments

For planar interfaces the construction reproduces Chombo's existing moments to 1e-15, so the
change is invisible wherever the geometry is genuinely planar. On curved geometry the apertures
drop from Chombo's third order to second — a median difference of 1.6e-2 at R/dx = 4,
converging at second order. This is the one deliberate accuracy cost of the design.

## 7. Decisions taken

### 7.1 `a_B` is stored twice — the one item that still needs Chombo

`PolyGeom::normal` and `PolyGeom::bndryArea` need `a_B = |∫dA|`, the magnitude of the area
vector, for `Σ(alpha_hi − alpha_lo) = a_B n` to hold exactly. But
`EBData::coarsenBndryCentroid` (`EBData.cpp:990`) weights child boundary centroids by that same
scalar, and a first moment only partitions when weighted by the measure it is a moment of — the
**true** area. The two are equal iff the patch is planar. The conflict is visible in the source:
the plain-sum rule sits commented out at `EBData.cpp:632`, directly below the area-vector rule
that replaced it.

Measured gap `(true − vector)/true` over all cut cells: 0 on an oblique plane; median 7.4e-4,
p99 4.8e-3 on a sphere; and up to **4.6e-1** on a rotated cube, where any cell containing an
edge has a folded patch.

**This cannot be fixed from chombo-discharge.** `EBIndexSpace` coarsens with
`EBISLevel(fineEBIS, …)` and `EBData::coarsenVoFData`, which we do not get to intercept. It
needs `bndryAreaTrue` in `BoundaryData` and a changed weight in `coarsenBndryCentroid`.

It is also independently a Chombo defect — two coarsening rules that disagree about what `a_B`
means — so it goes upstream as its own small PR, exactly as the covered-face centroid fix did.
**Milestone 1 does not block on it.** Without it, `bndryCentroid` carries a conservation error
of 8e-6…1.1e-2 on *coarsened* levels for curved geometry, which is the status quo today.

### 7.2 `bndryCentroid` on a covered face — DONE

Was a live Chombo bug, not a convention. `getFullNodeWithCoveredFace` computed the right value
and `fixRegularCellsNextToCovered` lost it for every cell an earlier covered cell had already
claimed by edge or corner. Fixed in Chombo-3.3#30, merged; superproject pin bumped. Our
`degenerate_moments` and Chombo now agree exactly.

### 7.3 What persists — SUPERSEDED: retained at generation, on our side

The earlier answer was "store the generators in `VolData`". Both halves of that were wrong.

**Recomputing is not affordable.** `BrentRootFinder` runs up to `MAXITER = 100` implicit-function
evaluations per crossing (`GeometryShop.cpp:1026`), twelve crossings per cell, and
`m_implicitFunctionGas` is a `NewIntersectionIF` that walks every electrode and dielectric on
every call. That is 10²–10³ IF evaluations per cell against a 2.7 µs body build.

**Storing it in `EBData` is not necessary.** The shop must evaluate the corners and bisect the
cut edges anyway, to build the body it takes the moments of. Keeping those values costs **zero
extra IF evaluations** — they are a by-product of work already done.

So: `CutCellSurface` is retained by the shop in a `LevelData<BaseIVFAB<CutCellSurface>>` on the
chombo-discharge side. `BaseIVFAB` is the existing sparse per-irregular-cell container, so this
stays inside the repository's working agreement on per-cell data. Cost is memory only — 23 reals
per cut cell with bends off, 46 with them on, roughly 37 MB per rank at 100k cut cells.

Milestone 1 does not need the surface retained at all; it is consumed and discarded per box.
Retention starts in Milestone 2.

### 7.4 New types — SUPERSEDED: private to `Source/Geometry/`

`CutCellSurface` is no longer a field in a Chombo struct, so the question of adding one to
`IrregNode` and `VolData` does not arise. It becomes an implementation type in
`Source/Geometry/`, alongside the shop that produces it. No new container: `BaseIVFAB` serves.

### 7.5 Why this needs no Chombo changes at all

`GeometryShop.H` is `public:` from the top of the class through line 121, then `private:`. That
line falls in exactly the useful place:

| reachable from a subclass | not reachable |
| --- | --- |
| `fillGraph` | `fillNodeValues` |
| `fixRegularCellsNextToCovered` | `insideOutsideFromNodes` |
| `getFullNodeWithCoveredFace` | `edgeData3D` / `edgeData2D` |
| `computeVoFInternals` | `BrentRootFinder` |

The two functions carrying the subtle covered-cell logic are public and can be reused verbatim.
The private ones are exactly the parts the prototype already reimplements: a node-value cache, a
classifier transcribed from `insideOutsideFromNodes`, and a bisection.

So the shop writes its own `fillGraph` and never calls `computeVoFInternals`, skipping
`Moments.cpp` and the `LSquares` machinery entirely. Same implicit-function cost as today,
cheaper quadrature. **Hypothesis to measure, not a promise:** it may be faster than `ScanShop`.

## 8. Milestone 1 — replace the moment machinery

Goal: identical topology, our moments, selectable at runtime, nothing in Chombo touched. This is
the step that lets the solvers, robustness and accuracy be tested before any refinement exists.

### 8.1 The switch

`Source/Driver/CD_Driver.options:5` currently reads:

```
Driver.geometry_generation  = chombo-discharge  ## Grid generation method, 'chombo-discharge' or 'chombo'
```

Add a third value, `polyhedral`. The existing two name *who* generates; the third names *what
the moments are*, which is the only thing that differs, and it stays accurate through Milestones
2 and 3.

`Source/Geometry/CD_ComputationalGeometry.cpp:212` already branches
`if (m_useScanShop) … else … new GeometryShop(…)`. A third branch goes there. The A/B comparison
— including a Poisson solve on identical grids — is then one input-file line.

### 8.2 New files, all in `Source/Geometry/`

| file | contents |
| --- | --- |
| `CD_CutCellSurface.H` | POD: 12 edge crossings, 8 corner values, optional per-face bends, optional apex |
| `CD_CutCellBody.H/.cpp` | body from a `CutCellSurface`; the seven moments as its integrals; the closure-and-range check |
| `CD_PolyhedralShop.H/.cpp` | `class PolyhedralShop : public ScanShop`, overriding `fillGraph` |

Deriving from `ScanShop` rather than `GeometryShop` keeps the box-level scan, so
`Driver.geometry_scan_level` still applies and `polyhedral` means "chombo-discharge's scan with
polyhedral moments".

### 8.3 What `PolyhedralShop::fillGraph` does

1. Fill a node-value cache over the ghost region — one implicit-function evaluation per node,
   shared by the eight cells meeting there.
2. Classify each cell with the transcribed corner test, including `copysign` so that negative
   zero is fluid (§5.2, §5.3).
3. Run `GeometryShop::fixRegularCellsNextToCovered` over the covered cells. It is public, it is
   now correct, and it handles the case §5.2 is about.
4. For each remaining irregular cell, bisect its cut edges through a per-direction edge cache so
   each edge is solved once and the four cells sharing it agree bit for bit (§5.6 rule 6:
   interpolate from `lexLess`-ordered endpoints).
5. Build the body and write `m_volFrac`, `m_volCentroid`, `m_bndryCentroid`, `m_areaFrac`,
   `m_faceCentroid` into the `IrregNode`.
6. Run the closure-and-range check. If a cell is declined, fall back to
   `GeometryShop::computeVoFInternals` for that one cell and count it. Currently zero cells in
   48,798, but the failure mode must be "slightly different moments here", never "broken graph".

`EBData::computeNormalsAndBoundaryAreas` then overwrites `bndryArea` and `normal` from
`PolyGeom`, which derives them from the apertures. Our construction satisfies
`Σ(alpha_hi − alpha_lo) = a_B n` exactly, so `PolyGeom` reproduces our values — a free
consistency check on every cell.

### 8.4 Cost control at generation

Two parts of the construction cost extra implicit-function evaluations and must be switchable:

- **Face bends** (§3.1). Found by pushing samples off the chord onto the zero set, so roughly
  15 samples per face plus gradients. Measured need: 0.00 bends per cell on an oblique plane and
  0.23–0.26 on rotated polyhedra, against ~7 on a sphere or torus. Default off; on for curved
  geometry.
- **Apex projection** (§3.1). A Newton projection along the numerical gradient, ~20 evaluations
  per cell. **Not required for conservation** — the depth-2 verification in §6.2 ran with the
  loop centroid as apex, not a projected one, and conserved to 9.4e-16. It improves how well
  `kappa` matches the true volume, not whether it partitions. Default off.

With both off, generation costs the same implicit-function evaluations as `ScanShop` does today.

### 8.5 Verification for this milestone

- **Topology is inherited, not asserted.** The graph, the classification and the fix-ups all come
  from Chombo, so criterion A holds by construction rather than by measurement. The 48,798-cell
  agreement result in §6.1 becomes the evidence that our body is consistent with the topology we
  adopt.
- **Per-cell identity check**, shipped, not just tested: `kappa = ⅓[½Σalpha − a_B(x_B·n)]`. It is
  one dot product per cut cell and it is what caught the Chombo covered-face bug.
- **A/B on a Poisson solve.** Same grids, same graph, only the moments differ. This is the number
  that says what the design costs, and it is the reason the switch exists.
- **The geometry test suite compares HDF5 exactly**, so `polyhedral` runs will differ wholesale
  from `chombo-discharge` runs. That is the point of the switch, not a regression — keep the
  default unchanged so the existing benchmarks stay valid.

## 9. Known limits

- **Depth 3 and beyond has not been re-verified** since the polyline representation and the
  four-crossing fix landed. It was clean at depth 3 under the earlier representation.
- **A cell with two or more covered faces.** `bndryCentroid` is a single `RealVect`, so a
  concave corner of the geometry — two covered faces meeting at one cell — cannot record both
  face centres. Phase 0 improves such a cell from "cell centre" to "one of the two faces", not
  to correct. Rare, and orthogonal to the refinement operator.
- **A cell edge cut twice.** Where a solid wedge is thinner than a cell, the implicit function
  is negative at both ends of a cell edge and positive in between, so bisection records no
  crossing. This is a resolution limit of *any* edge-crossing representation, marching cubes
  included — the input data is missing the feature — not something a better triangulation
  fixes. No instance survives in the present nine-geometry suite.
- **Second-order apertures on curved geometry**, per §6.5. Deliberate.
- **`cube_axis` exercises only the degenerate path.** An axis-aligned cube has no ordinary cut
  cells at all — every one of its 4,376 irregular cells is touch-only — so it validates §5.2
  and nothing else.

## 10. Where the milestones are

- **Milestone 1** — this file. Replace the moment machinery, behind the switch. No Chombo
  changes, no refinement, no new EBIS.
- **Milestones 2 and 3** — `PLAN_MILESTONES.md`. Refining parent cells on top of Milestone 1,
  then AMR-aware geometry generation where planar regions stay coarse and curved regions refine.

## 11. Where the code is

- `Prototypes/CutCellRefinement/main.cpp` — the exporter. Writes per-cut-cell moments, edge
  crossings, corner values, face bends and apex. Geometry by `Prototype.geometry` =
  `plane | sphere | swept_sphere | torus | polyhedron | tessellation | rough_sphere`, all
  ParmParse-configurable; `Polyhedron.shape` = `cube | simplex` with Euler angles.
- `Prototypes/CutCellRefinement/spike/refine_spike.cpp` — the C++ port, ~1200 lines,
  standalone. Reads the exporter's CSV, builds, refines to a given depth, and reports
  conservation, topology and timing.
- Python prototype, session scratchpad `scratchpad/eb/`: `patch.py` (body from crossings),
  `body.py` (clipper and moments), `samelevel.py` (criterion 1), `deepN.py` (criterion 2),
  `refusals.py` and `closure.py` (the diagnostics that found §5.1).
