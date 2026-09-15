# Lazy expansion of the EBIS level grids

Draft for revision. Everything under "Verified" was read out of the code this session; everything
under "Decisions" is open.

## Verified

1. **ScanShop's decomposition is not a tiling.** `buildCoarseLevel` splits the coarsest domain with
   `domainSplit(domain, maxGridSize)`. Every finer level is then built *from the coarser layout*:
   each coarse box is refined by two, and

   - `Covered` -> the whole refined box is kept, unsplit;
   - `Regular` -> the whole refined box is kept, unsplit;
   - `Irregular` -> `domainSplit(fineBox, maxGridSize, maxGridSize)`.

   So regular and covered boxes double in size at every level, deliberately, and only irregular
   regions are held at `maxGridSize`. `retainBox` is applied to the whole refined coarse box, so a
   box is dropped or kept entire.

2. **`TiledMeshRefine` is a different construction.** The AMR grids come from tags on a tile
   lattice set by the blocking factor. They do not align with ScanShop's boxes and never have --
   `EBISLayoutImplem::define` copies between the two through a `Copier`, which does not care.

3. **An unfilled box reads as all-regular.** `EBGraphImplem::define` sets `m_tag = AllRegular`.
   A box the level does not describe is therefore silently fluid.

4. **All-regular and all-covered boxes cost almost nothing.** They are a tag; no per-cell graph is
   allocated, and `BaseIVFAB` over an empty irregular set holds nothing. This is what makes the
   huge boxes in (1) affordable.

5. **Nothing reads the persistent level directly.** `EBISLayoutImplem::refine` reads
   `m_ebisBoxes`, built from `localGraph`/`localData` -- per-layout copies. Ratios above two read
   `m_fineLevels`, which are themselves layouts.

6. **The classification of a box that was never retained is recomputable.** `ScanShop::isRegular`
   and `isCovered` answer from the implicit function over a grown box, with no moments and no
   graph.

## The three kinds of region

Expansion is only ever needed where the default in (3) is wrong:

| region | needed? | cost | source |
| --- | --- | --- | --- |
| regular | **no** -- the default is already correct | zero | -- |
| covered | yes, or fluid appears inside the solid | a tag; big boxes fine | classification alone |
| irregular | yes | `maxGridSize` boxes, moments | cut from a coarser level's surfaces |

Only the third needs geometry. That is also the only one that needs a parent, and therefore the
only one that interacts with the surface store.

## Decisions -- settled

**D1. The level grows**, not just the layout. Growth is not needed for correctness -- (5) -- but
the layout cache is keyed by object identity (`std::map<DisjointBoxLayout, EBISLayout>`;
`BoxLayout::operator==` compares a `RefCountedPtr<int>`), so every regrid misses it even when the
grids are identical in content. Without growth the same geometry is re-cut on every regrid, for
every realm and level, for the life of the run. Treated as a cache: by (5), being wrong about what
to add costs time, not correctness.

**D2. ScanShop keeps recursing with the rule it already has.** Expansion is not "extend a
decomposition", it is "widen `retainBox`": the base predicate is unconditional (`return true`),
and the walk already knows how to span the whole domain at every level. Partial coverage exists
only because the polyhedral shop says no. Irregular boxes split as they already do. Nothing is
reconstructed or reimplemented, and the mechanism stays polyhedral-only, since the other
generators have nothing to widen.

**D3. Two grids: moments everywhere, chords only where nothing finer covers the cell.**

- a cell **covered by a finer level** carries moments alone, from coarsening. It needs no chords,
  because refining it returns the real children, which exist;
- a cell **not covered** carries chords. Refining it chords them, and each chorded child is itself
  a polyhedron, so the recursion continues without anything further being stored.

That is a total partition and it closes the refinement question: every cell refines, by one route or
the other, to any depth. Factor-4 refinement costs nothing structural -- the first halving returns
real children, the second chords them.

An earlier draft withdrew this in favour of one grid carrying `(ancestor surface, depth)`. The
withdrawal was wrong, and for a reason already admitted elsewhere in this plan: chording needs a
*polyhedron*, not a `CutCellSurface`, and a chorded child is a polyhedron. The objection -- that the
recursion would need a cell with no surface as a parent -- does not arise.

**D3a. The seam is made exact by a multichord face, not by a multichord cell.**

Where a chord-carrying cell abuts a coarsened one, its face toward that neighbour follows the
neighbour's children's chords rather than a single chord of its own -- piecewise on that face only.
Chording it then reproduces the neighbour's children across the shared face.

The whole-cell version of this does not work and was measured: unioning all 2^SpaceDim children into
one body *is* exactly invertible -- `refine(union(children)) == children` to 1e-16 on four
geometries, zero bad cells -- but it needs 45 to 58 polygons against a cap of 20, and worse, it
compounds. A coarsened cell holding such a union has 2^SpaceDim interface pieces, so its own parent
would need 2^(2*SpaceDim), and so on: the representation grows as 8^depth. Factor-4 refinement
reaches that immediately, and 20 of 145 input files in the suite use a ref_rat of 4.

Confining the multichord to seam *faces* avoids all of it: one face, one level, no interface
hierarchy. The cost is vertices on two polygons rather than polygons on a body.

**D4. Nesting is inherited, not enforced.** Nothing interpolates between EBIS levels, so there is
no proper-nesting radius to hold. The requirement that does exist -- the coarse level covering
`coarsen(fine)` plus a one-cell ring, so that coarsening can read a cell's neighbours -- is exactly
the rule `TiledMeshRefine::nestFrom` already applies when it builds the hierarchy:

    const Box grown     = grow(a_fineBox, 1) & fineTileBox;
    const Box coarsened = coarsen(grown, a_refToFine) & a_thisTileBox;

So if the retain regions are a properly nested hierarchy -- which the simulation grids are, by
construction -- then per-level retention maintains the parent invariant automatically. No upward
recursion in the expansion, and no post-condition scan: a violating state cannot be built.

The precondition therefore sits on the *input*, where it is small and cheap: assert in debug that
the region set handed in is nested, rather than scanning every box of every level afterwards. Two
things must hold, both structural: the request comes from a nested hierarchy, and the buffer is
applied per level -- uniform `n` preserves containment because
`coarsen(grow(R, n))` is inside `grow(coarsen(R), ceil(n/2))`, whereas growing only the fine level
would break it.

**D5. Box kinds are re-identified by the walk itself**, as they are on the first build. Falls out
of D2.

**D6. `ebis_buffer`, default 0.** Carries the levels wider than the simulation asks, so that
`fillEBISLayout` stays a pure copy for several regrids rather than cutting on each. Notes:

- it quantises to the box size -- `retainBox` keeps whole boxes, so a buffer only bites when it
  crosses a box boundary. Document it as cells grown before the box test;
- it applies per level, or the nesting recursion halves it away going up;
- it is measurable: count boxes cut per regrid, and zero means the copy was pure. Tune by
  experiment, not by guess.

### Where it goes wrong quietly

The failure D4 exists to prevent is silent, which is why it is worth making structural rather than
checked. A box that is never created is not misclassified -- it is absent, and absence reads as
all-regular from `EBGraphImplem::define`. The worst case is not a covered box reading as fluid but
a *cut* one: the embedded boundary vanishes over that patch, and the solver sees open fluid through
the electrode at exactly the place the physics lives.

Nothing catches it. An all-regular region is a valid graph, `checkGraph` passes, and the
prototype's divergence and face checks iterate irregular cells -- of which a hole has none. If the
hole runs consistently through every level, the coarsening check agrees with itself too.

This is the guarantee #731 had structurally, by hardwiring the AMR grids to the coverage regions
(`retainBox`: *"Nothing may ask for embedded boundary data outside what is kept"*). Phase 2 lifts
that, so the guarantee has to come from somewhere else -- and D4 puts it back where it was, in the
shape of the hierarchy rather than in a check.

## Sequence

1. D1-D5 settled.
2. Fill the remainder in `EBISLayoutImplem::define` -- the mechanism, no caching. Gate: the
   validator clean, and RodSphere 2-D `polyhedral` completing *and* validating.
3. The cache (D1), with the level growing. Gate: unchanged results, fewer cuts on a repeat regrid.
4. Nesting (D4) enforced, with a deliberately badly nested case to prove it fires.
5. Phase 2: the regrid copy path; lift the coverage-only restriction.
6. Phase 3: depth sweep.

## The rule the chain actually rests on

Stated by the author, and it is a total partition. At generation time every cell is in one of two
states:

- **no finer grid above it** -> it is described by a *polyhedron*;
- **a finer grid above it** -> its moments come through coarsening.

From which: **every cell can be refined.** Coarsened, and its children exist with `m_finerNodes`
pointing at them; polyhedral, and it can be chorded into 2^SpaceDim children. There is no third
case, and the recursion terminates because a chorded child is itself a polyhedron.

Two things recorded earlier in this plan were wrong against that rule and are withdrawn:

- **"Extended cells cannot be parents" was an artifact of the store, not a property of the design.**
  Chording needs a *polyhedron*, not a `CutCellSurface`; the surface is only the cheap thing a body
  is rebuilt from. A chorded child has populated polygons and chords again -- `refineSubtree`
  already does it. The 258-of-420 measurement showed the store was incomplete against the
  invariant, not that the invariant fails.
- **"The finer level must exist beneath every cut cell" was the wrong reading of
  `EBCoarseFineParticleMesh`.** It needs *a* refinement of each irregular coarse cell, and the rule
  above always provides one.

## The seam, measured

Chording a coarse cell rebuilds its body from the stored surface, and `reconcileSeam` never touches
that surface -- it overwrites `m_areaFrac` and re-derives the interface in `EBData` only. So a
chorded child inherits the *unreconciled* aperture while its neighbour across the seam descends
from the fine cells the reconciled value came from.

Measured at 64^3, two levels, `refine_angles 15`:

| geometry | seam faces | changed by reconciliation | worst delta |
| --- | --- | --- | --- |
| sphere | 0 | 0 | -- |
| torus | 0 | 0 | -- |
| polyhedron (cube) | 266 | 0 | 1.2e-15 |
| swept sphere | 438 | 12 | 1.3e-3 |

Small, because the seam sits where the pre-pass stopped tagging, which is where the interface is
near-planar, which is where a coarse chord and the sum of the fine chords coincide. The cube is the
clean case: its seam lies on planar faces and the agreement is exact.

Topology was never at risk -- `reconcileSeam` writes areas, not connectivity, and the graph is what
#731 certifies at `faceMismatch 0`. What is at stake is metric agreement to better than 1e-3.

## The multichord seam, tested across a real coarse-fine boundary

The first seam harness was worthless: it walked a *uniform* slab and, for every cut cell, subdivided
one of that cell's own faces. Nothing was coarsened, nothing had children, and a multichord face
always met a single-chord neighbour, so the exported surface cracked by construction. The numbers it
produced (closure, per-cell conservation) were real but answered a question nobody asked.

`validateTwoLevelSeam` builds the configuration the multichord exists for. Half the domain is carried
at `dxC`, half at `dxC/2`; the coarse column adjacent to the interface puts the abutting fine chords
on the shared face and single chords everywhere else; every other cell on both sides is
single-chorded. Both blocks together cover the whole sphere, so the union of the interface triangles
is a closed surface *if and only if* the seam agrees -- which makes edge valence the test, not
eyeballing an STL.

Sphere, radius 0.25, 16 coarse cells across `[-1,1]`, interface at `x = 0` through the centre:

| | triangles | edges | valence | open boundary |
| --- | --- | --- | --- | --- |
| multichord | 664 | 996 | all 2 | none |
| single chord | 648 | 992 | 40 at 1, rest 2 | 40 edges, total length 3.12 |

The multichord surface is watertight. All 28 edges lying in the seam plane pair one coarse triangle
with one fine triangle -- the coarse cell's chorded sub-faces and the fine cells' real faces are the
same segments, not merely close ones.

The single-chord surface tears along the full great circle: its 40 open edges have *every* endpoint
at `x = 0` exactly, and their total length 3.12 is twice the circumference 2*pi*0.25 = 1.571 -- the
coarse rim and the fine rim, both open, nowhere else damaged.

Conservation over the 12 cut cells on the interface, coarse face fraction against the sum of the four
fine sub-faces:

| | worst error |
| --- | --- |
| multichord | 2.5e-13 |
| single chord | 4.1e-2 |

So the multichord is what makes the coarse-fine boundary topologically closed and conservative; the
single chord is neither, and the failure is not small.

## The rotated cube, 1000 orientations

A sphere is the easy case: nothing is planar, no facet can land on a cell face, and the chords never
come out collinear. A cube is the hard one, and it is also the shape the polyhedral representation
exists for. `sweepTwoLevelSeam` runs the same two-level configuration over a cube swept through ten
values of each Euler angle across a quadrant -- 1000 orientations, 25527 seam cells.

Three things had to be fixed or separated before the sweep said anything.

**`mergeCoplanar` refused a seam face the body covers completely.** A cube face landing on a cell
face puts the crossings on the face's own edges, `faceWalk` returns slivers of no area, and the edge
cancellation leaves either nothing or fewer than three edges -- reported as `why 6`, merge reasons 1
and 6. An empty merged face is the answer there, not a failure: the aperture is zero on both sides
(measured `singleChord 1.0e-12`, `fineSum 0`). Refusing dropped the cell back to a single chord,
which is the crack the multichord exists to remove. It now returns success with no polygon when the
input sub-faces enclose less than a weld tolerance of area. That alone was 17 refusals in 600 seam
cells at grid-aligned placement, and 0 after.

**The STL was written at six digits.** Collinear segments that partition the same line exactly then
miss each other by ~1e-7, and a T-junction reads as a crack. At seventeen digits the same files lose
almost all of their apparent defects. Any watertightness claim measured off a six-digit export is
worthless.

**A body face lying on a cell face produces no cut cells at all**, so no interface polygon, so a hole
in any export of interface triangles. A cube of half-width 0.25 on a 0.125 mesh has all six faces on
cell planes; at uniform resolution its surface shows 48 open edges at 0 degrees and 32 at a generic
angle, total length exactly the perimeter of the missing facets. Nudging the centre to
(0.0131, 0.0172, 0.0193) takes every one of those to zero, at every angle tested. Nothing is wrong:
the covered/regular boundary carries that surface, and no cut cell is asked to.

With those separated, the sweep measures the seam and nothing else.

| | multichord | single chord |
| --- | --- | --- |
| seam cells refused, of 25527 | 0 | -- |
| worst face-fraction error, generic centre | 3.3e-16 | 4.8e-1 |
| worst face-fraction error, centre on the grid | 6.0e-11 | 3.8e-1 |
| worst closure residual | 3.9e-16 | -- |
| rotations whose surface is watertight | 854 of 1000 | 0 of 1000 |
| total unclosed area over the sweep | 7.8e-2 | 7.3e+0 |
| median unclosed area per rotation | 3.2e-4 | 6.4e-3 |

Watertight here means: every edge shared by exactly two triangles once T-junctions are resolved --
the coarse side merges collinear fine chords into one segment, which leaves the fine side's vertices
in the middle of it. That is a vertex-sharing difference, not a gap, and the partition is exact.

The 146 rotations that do not close leave slivers of one shape: a coarse segment lying along a
coarse cell edge, against two fine segments that dip past that edge by a fraction of a fine cell.
The multichord replaces the chord on the seam face only; where the contour grazes the coarse cell's
*other* faces, the coarse cell's own chord there is what the interface closes against, and the fine
cells resolve a dip it cannot see. Worst case over the sweep is 3.3e-3 of area against a coarse face
of 1.6e-2, on a surface of 1.34. The single chord fails the same test in every rotation and by two
orders of magnitude more area, so the multichord is doing its job; this is the part of the seam it
does not reach.

## What will bite

- **The prototype harness no longer measures locality.** It hands over parents gathered from local
  boxes only, so under MPI it aborts instead of counting refusals. Fix before trusting any MPI
  number.
- **`fillRefinedGraph` inherits #730's threshold/ghost inconsistency.** A cell below
  `m_volumeThreshold` is dropped only inside the valid region.
- **`classifyFromParents` still lacks the `isRecorded` check**, and cannot simply be given it:
  `defineDegenerate` populates no polygons, so `refine()` on such a body yields covered children
  silently.
- **2-D `geometry_refinement = 2` is broken**, pre-existing for the sphere; the swept sphere's
  `fm 20 divBad 20 worstDiv 0.5` has *not* been attributed and may be from the `clip` change.
- **#731 finding 1**, the `getCurvatureTags` out-of-bounds, is unverified and HIGH.
