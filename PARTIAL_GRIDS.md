# Partial grids for EB generation: what is established, and what is not

Written at the point where the grid classification needs redesigning. Everything below is either
quoted from the source or measured from a run; the last section says what is neither.

## Where the code is

Branch `cutcell-refinement` (#731). `375057688` and earlier are pushed; `8d71c2b48` and `f95a324ca`
are local.

| commit | what it does |
| --- | --- |
| `82be179de` | writes the reconstructed EB surface as STL, per phase |
| `527814b7f` | **took the AMR grids from the index space -- wrong, see below** |
| `954584947` | one STL per level beside the composite |
| `375057688` | writes the surface from the grids instead of from generation; keeps the generator alive; `checkGridsAgainstIndexSpace` |
| `8d71c2b48` | `restrictFace`, `mergeCoplanar`, `closeInterface` in `CutCellBody`, plus `detail::sameVertex`. **Nothing calls them.** Generation is unchanged. |
| `f95a324ca` | reverts `527814b7f` and fixes the check to the right invariant |

## The mistake, and what it taught

`527814b7f` made the simulation's grids equal to the index space's, on the theory that the two should
describe the same cells. That is the wrong invariant. `retainBox` reaches a ghost width past what it
was asked for precisely so that a box's ghost cells have geometry in them; grids that reach as far
have ghost cells the index space never carried.

Measured on ProfiledSurface 3D, four levels, ghost width 2:

| level | ghost cells outside what the index space carried |
| --- | --- |
| 0 | 0 |
| 1 | 77824 |
| 2 | 202848 |
| 3 | 462592 |

An unfilled graph reads as regular fluid rather than as missing, so none of that announced itself.
Reverted. `checkGridsAgainstIndexSpace` now asserts the real invariant -- `grow(grids, ebGhost)` is
contained in what the index space carries -- and errors out otherwise. It passes.

## What the EB machinery actually does with a box it does not carry

Traced through the source, not inferred.

`EBISLevel::fillEBISLayout` -> `EBISLayoutImplem::define` does two things:

```cpp
localGraph.define(a_grids, 1, ivghostgraph, graphFact);   // a graph for EVERY requested box
a_graph.copyTo(interv, localGraph, interv);               // overwrite only where the level carries
```

The requested layout is fully populated by `EBGraphFactory::create`, which is
`new EBGraph(region)` -> `EBGraphImplem::define` -> **`m_tag = AllRegular`**. The level's data is
copied over that. So:

- **a regular region may be omitted** -- the default is already the right answer;
- **a covered region may not** -- omitted, it reads as fluid, silently.

There is no hook between `localGraph.define` and `copyTo` for a caller to say "this region is
covered". `EBGraphFactory` holds only a `ProblemDomain`.

## What a carried regular or covered box costs

Nearly nothing in graph terms. `EBGraphImplem::size` returns `sizeof(int)` unless the region has
irregular cells; `linearOut` writes a single "secretCode" (0 covered, 1 regular, 2 has irregular) and
per-cell `GraphNode`s only for code 2. `linearIn` synthesises the source with `setToAllCovered()` or
`setToAllRegular()`. On processor, `EBGraphImplem::copy` returns on the tags before touching data,
and `BaseFab<GraphNode> m_graph` is only defined when the tag has to become `HasIrregular`.

So the cost of carrying a box that is not cut is its entry in the `DisjointBoxLayout` and in every
`Copier` -- box count, not storage. The inflation worth avoiding comes from *many small* boxes, and
those are the cut-cell boxes, which have to be carried anyway.

**Not checked:** whether `EBData` short-circuits on the tag the same way. If it allocates per box
regardless, a carried covered box is not as cheap as the graph analysis suggests.

## The nesting problem

`ScanShop::buildFinerLevels` recurses coarse to fine. A coarse box that was Covered or Regular is
refined whole, deliberately, so that large regular regions stay one big box and the metadata does not
inflate on large domains. A coarse box that was Irregular is `domainSplit` to `maxGridSize` and each
piece reclassified. That decomposition does not leave the fine irregular grid properly nested inside
the coarser one.

Measured on ProfiledSurface 3D, `min_block_size = max_block_size = 4`, `buffer_size = 1`, 12 ranks:

| level | cells that coarsen outside the level below | the same, grown by the buffer |
| --- | --- | --- |
| 1 | 0 | 0 |
| 2 | 0 | **972** |
| 3 | 0 | **1476** |

Contained, but not properly nested. At block size 8 both columns are zero, so it is the small blocks
that expose it.

## What the classification API offers

`ScanShop::InsideOutside(region, domain, probLo, dx, dit)` is **a cache lookup, not a classifier**:
it matches `a_domain` to a level and returns `(*m_boxMap[whichLevel])[a_dit]`, and `MayDay::Error`s
on both fallbacks. It cannot answer for a box that is not already in that level's layout.

What classifies is `ScanShop::isRegular(Box, probLo, dx)` and `isCovered(Box, probLo, dx)` --
protected, const, taking an arbitrary box at an arbitrary spacing, scanning the implicit function.
`buildFinerLevels` already uses them on freshly split boxes, so "split a box and reclassify the
pieces" is an operation that exists.

Two things that constrain any redesign using them:

- classification is done on the box **grown by `m_ebGhost`**, so a piece carved out of a large
  regular box next to a newly irregular neighbour will often come back irregular on its own, because
  the collar reaches into the feature;
- `m_boxMap` is a `LayoutData`, so a level's classification is tied to its layout: re-partitioning
  means rebuilding `m_grids[lvl]` and `m_boxMap[lvl]` together, not patching them.

And one thing that does not exist: any way to record that a box must be present *because a finer
level needs it*, as distinct from because the geometry cuts it. `GeometryService::InOut` has three
values and no room for the distinction; `retainBox` can filter boxes but cannot add them.

## The shape the constraints leave

Everything any consumer requests, grown by the ghost width, must be carried, and covered regions
inside that must be described. Regular space need not be. That does not forbid building the grids
over the cut cells with `TiledMeshRefine`; it says the carried set is cut-cell boxes plus a collar,
and the collar has to be classified honestly where it is covered.

## How big the overlap actually is

The restriction sweep computes, for each irregular box on a level, `coarsen(grow(box, 1), 2)` on its
parent. Where that lands outside the parent's own irregular boxes it overlaps boxes the parent had
classified regular or covered, and those would have to be split. Measured on ProfiledSurface in 2D,
serial so the counts are global, `min_block_size = max_block_size = 4`, `max_ebis_box = 8`,
`buffer_size = 1`, `refine_geometry = 6`:

| fine -> parent | forced | already irregular there | spilling | of which on regular/covered | parent boxes hit | cells in them |
| --- | --- | --- | --- | --- | --- | --- |
| 3 -> 4 | 180 | 160 | 20 | 20 | **6** | 384 |
| 4 -> 5 | 180 | 170 | 10 | 10 | **2** | 512 |
| 5 -> 6 | 570 | 560 | 10 | 10 | **2** | 128 |
| 6 -> 7 | 290 | 290 | 0 | 0 | 0 | 0 |
| 7 -> 8 and coarser | | all | 0 | 0 | 0 | 0 |

So the spill is 3 to 6 per cent of what the restriction forces, it stops entirely below level 6, and
it touches **two to six boxes per level**. Splitting those is cheap. The concern that the parent
would have to be repartitioned wholesale is not borne out here.

Two things to be careful of before generalising:

- **This geometry has no covered boxes at all.** Every level reports `covered 0`: a profiled plane
  cut by a half-space leaves regular and irregular boxes only. A geometry with bulk solid may spill
  onto covered boxes, which are the ones that cannot be dropped.
- **The hierarchy truncates.** Levels 0 to 2 hold no boxes: `retainBox` drops them because the
  curvature pre-pass asks for nothing at that depth. That is partial coverage working as intended,
  but it means the finest levels in these runs are not exercised, and the measurement only covers
  levels 3 and coarser.

Counting both cells and boxes mattered. The cell spill looks negligible either way, but the box
count is what costs, and it is the number that says splitting is affordable.

## How the grids reach the generator, and whose grids they are

`EBISLevel` does not receive a layout -- it asks for one:

```cpp
(const_cast<GeometryService*>(&a_geoserver))->makeGrids(a_domain, m_grids, a_nCellMax, 15);
```

and `ScanShop::makeGrids` matches `a_domain` against `m_domains` and hands back its own, by assignment:

```cpp
if (m_hasThisLevel[whichLevel]) { a_grids = m_grids[whichLevel]; }
```

So `EBISLevel::m_grids` *is* `ScanShop::m_grids[lvl]`. There is no second decomposition and nowhere for
anything to nest them on the way through; the `else` branch is a `MayDay::Warning` reading "This should
not happen!". `EBISLevel` then walks that layout calling `InsideOutside` (the `m_boxMap` cache lookup)
and `fillGraph`, which is the virtual `PolyhedralGeometryShop` overrides to build the bodies.

ScanShop decides where, PolyhedralGeometryShop decides what, EBISLevel iterates. So the polyhedra are
built on grids that are not properly nested, and nothing downstream is positioned to fix that.

The properly nested grids in the system are the curvature regions from
`ComputationalGeometry::getCurvatureTags`, which do go through `TiledMeshRefine` -- twice, with a comment
that nesting travels upward so tagging deeper can widen a level already clustered. Those never become
the EBIS layout. They reach the shop through `setCoverage` and are used only as a filter. A filter can
remove boxes from a decomposition that is not nested; it cannot make the survivors nested.

## retainBox: one call site, and what it prunes

Called in exactly one place, `ScanShop::buildFinerLevels`, and the position in the function is the point:

```cpp
const Box fineBox = refine(coarBox, 2);
const GeometryService::InOut& boxType = (*m_boxMap[coarLvl])[din];

if (!this->retainBox(fineBox, fineLvl)) {
  continue;                          // before the type switch, before domainSplit
}
```

It takes the **refined** box on the **finer** level and decides whether that refinement is created at
all. `ScanShop`'s own version returns true unconditionally; `PolyhedralGeometryShop` overrides it to test
`grow(coarsen(a_box, 2), m_ebGhost)` against `m_coverageRegions`. `buildCoarseLevel` never consults it,
so levels at and above the scan level are always built in full.

Because the test precedes the type switch, dropping an **irregular** refinement skips `domainSplit`
entirely, and the saving is double:

- the 2^D pieces the split would have produced are never created, nor anything they would have spawned:
  it prunes a subtree, not a level;
- the cut cells stay at level L's resolution rather than being re-cut at L+1, so the graph nodes and the
  moments for that region are held once at the coarse spacing instead of 2^D times at the fine one.

Same geometry, same requested depth (`max_amr_depth = 12`, ProfiledSurface 2D, serial):

| generation | cut-cell | regular | covered | total |
| --- | --- | --- | --- | --- |
| `polyhedral`, retainBox active | 83 | 190 | 0 | **273** |
| `chombo-discharge`, retainBox always true | 149121 | 199888 | 170524 | **519533** |

A factor of about 1900, entirely subtree pruning. The zero in the covered column is the same effect: the
pruning stops before reaching the depths at which the bulk solid would be decomposed at all. At ten
levels the unpruned hierarchy is 23125 cut-cell, 22022 regular and 15266 covered, 60413 boxes.

**ScanShop knows nothing about curvature.** Zero occurrences of `curvature`, `angle` or `normal` in
either of its files. Its irregular split is a uniform `domainSplit` to `maxGridSize` followed by a binary
occupancy test per piece, on the box grown by `m_ebGhost`. There are two distinct savings and only one is
curvature-driven: ScanShop keeps a Regular or Covered coarse box whole when refining it, which is driven
by emptiness; `retainBox` prunes subtrees, which is driven by the pre-pass.

The two pull against the nesting fix. The first depends on big regular boxes staying whole, which is what
the nesting bites into. The second depends on dropping refinements, which is what creates the holes that
make containment fragile. And a box the sweep forces onto a parent is one `retainBox` would have pruned,
so it re-opens a subtree beneath it unless the sweep can also say "keep pruning below this".

## What the containers will and will not do

**`NeighborIterator` is `DisjointBoxLayout`-only.** `m_neighbors` and `computeNeighbors()` are protected
members of `DisjointBoxLayout` and `NeighborIterator` is its friend. `BoxLayout`, which permits
overlapping boxes, has no neighbour machinery at all. So an overlapping layout of "regular/covered plus
grown irregular" cannot be asked for neighbours -- the one query wanted is the one the permissive
container lacks. Nor can it be promoted: `DisjointBoxLayout::define(const BoxLayout&)` checks
`isDisjoint` and throws.

`DisjointBoxLayout::closeN(neighbours)` does let a layout be closed with a neighbour list supplied from
outside, so the adjacency is not welded to the box set. It is protected.

**The neighbour reach is `grow(box, 1)`**, which is exactly what the restriction needs, since
`coarsen(grow(fineBox, 1), 2)` extends at most one coarse cell past its parent. The build is a windowed
sweep, not a quadratic scan -- boxes are x-sorted, `maxI` is the largest extent in x, and the window
advances monotonically -- so it is O(N x window). But it is built over `dataIterator()`, so `m_neighbors`
is populated only for boxes this rank owns, and it includes periodic images which must be `unshift`ed
before intersecting.

**The boxes are replicated; the classification is not.** `BoxLayout::m_boxes` is a
`RefCountedPtr<Vector<Entry>>` with `Entry = {Box box; unsigned int m_procID;}` -- every rank holds every
box, about 36 bytes each in 3D. `ScanShop::m_boxMap` is a `LayoutData<InOut>`, so a rank knows the
classification only of boxes it owns. That asymmetry is the whole source of the two-ranks-carving-the-
same-box race: rank A can see rank B's box but not that it is irregular. An `InOut` is four bytes against
the thirty-six already replicated, so closing the asymmetry costs a tenth of what the boxes cost. The
price is that every rank then redoes the whole sweep rather than its share.

**`IntVectSet` is viable as box algebra and fatal as a cell set.** A full box is a single `&full` sentinel,
so holding a huge region costs nothing and subtracting a small box grows it by O(log) nodes;
`TreeIntVectSet::createBoxes` emits one box per full node, so it grades rather than fragmenting. But
`numPts()` returns `int` -- at 500k^3 that overflows silently, and every count printed by the probes in
this session came through it -- and `IVSIterator` is per-`IntVect`. The probes used exactly the operations
that do not scale. What the sweep needs is box-minus-box, which decomposes into at most `2*SpaceDim`
boxes and enumerates no cells. No such helper exists: nothing in Chombo's `BoxTools` or in `Source`
matches `removeBoxFromBox`, `boxSubtract` or `complementOf`, and `IntVectSet::operator-=(Box)` is the
cell-set one.

## What `InOut` cannot say

Two things now, both needed by the sweep and neither expressible in a three-valued
`GeometryService::InOut`:

- that a box is carried **because a finer level needs it**, as distinct from because the geometry cuts it;
- that a box is **present but not a reason to refine below it**, so that `retainBox` keeps pruning the
  subtree it would otherwise re-open.

## How the EBIS machinery operates: five facts

**F1 -- classification is per box, and Regular/Covered boxes are never generated.** The per-box loop
in `EBISLevel`:

```cpp
inout = a_geoserver.InsideOutside(region, a_domain, a_origin, a_dx, din);
if      (inout == Regular) ebgraph.setToAllRegular();
else if (inout == Covered) ebgraph.setToAllCovered();
else                       { a_geoserver.fillGraph(...); ebgraph.buildGraph(...); }
```

A box classified Regular or Covered costs one tag and never touches the generator. `fillGraph` -- the
polyhedral path -- runs only for Irregular boxes. So the "outside" costs a tag per box, not a
computation per box.

**F2 -- but it needs a box.** `EBISLevel` classifies only the boxes `makeGrids` handed it. The
four-argument `InsideOutside` at line 288 is inside `makeBoxes`, the fallback decomposition used only
when the generator supplies no layout. Whatever is not in a box is never asked, and reads `AllRegular`.

**F3 -- a hole's parent is always generated, never coarsened.** `coarsenFrom` builds `fineCoverage`
from `coarsen(a_fineEBIS.m_grids.boxArray(), 2)` and overwrites only there, and by its own comment leaves
the outermost ring of the fine footprint alone. A hole has no fine footprint, so the coarse cell above it
keeps the body `fillGraph` built. The cuttable parent exists by construction.

**F4 -- the fill-from-parent contract is already declared, on the Chombo side.** The branch above adds
to `GeometryService`: `numSurfaceComponents`, `getSurfaces`, `refinedFillParentDx`, `fillRefinedGraph`,
all defaulting to refusal. What it never had was a guarantee of what it fills *from* (F3) or a definition
of which cells are holes.

**F5 -- the layout is assignment, not negotiation.** `makeGrids` does `a_grids = m_grids[whichLevel]`.
Whatever `Vector<DisjointBoxLayout>` is built becomes the EBIS layout verbatim. Nothing downstream
re-nests, re-splits, or fills gaps.

## The two questions, constrained by the facts

**Q1 -- classifying everything outside the nested cut-cell grids.** By F1 and F2 it splits by type:

- *regular outside* may be omitted: the default is the right answer;
- *covered outside* must be carried, **in full, at every level**. Not a collar. Three reasons, each
  sufficient: AMR level 0 is `domainSplit` of the whole domain and requests every cell; physics tags put
  fine grids inside dielectrics, which are covered to the gas phase; and in multi-fluid the solid phase's
  covered region is the gas phase's fluid, so between the two index spaces essentially everything is
  requested somewhere.

By F1 that is still cheap: interior solid is bulky and describes in few large boxes, and refining a
covered box whole keeps it one box per level regardless of depth. It is also exactly the constraint the
current `retainBox` placement violates, since it drops covered refinements on the same terms as
irregular ones. Moving the test into the Irregular branch -- `boxType == Irregular && retainBox(...)`
-- keeps every covered region at one box per level and leaves the subtree pruning where the subtree is.
A non-retained Irregular box then falls through the `else if` chain with nothing pushed: the outer
switch has no final `else`.

The one place box arithmetic happens: covered boxes must be disjoint from the nested irregular grids,
and their shared boundary is the surface. Each level's covered set is the coarser level's covered boxes
refined whole, minus the nested grids -- a few big boxes minus a rim.

**Q2 -- filling the hole.** The hole is load-bearing: it terminates the recursion. Retaining an
Irregular box whole instead does not make a hole, it makes a seed -- at the next level it is an
Irregular parent and is either split (pruning lost) or retained whole again, one full-width box per
level to the bottom, each through `fillGraph`. Absence is the only cheap terminator, and it is the
fourth state `InOut` cannot express seen from the other side.

A grid that lands in a hole must therefore not read the factory default. By F3 its parent has a body;
by F4 the contract to cut that body is declared. What #731 lacks is the *trigger*: today only
`EBGraphFactory` answers for an uncarried region, with `AllRegular`. So the actual dependency between
the PRs, stated sharply: #731 creates holes; #732 supplies the fill path; #731 shipped the holes without
it. The definition #732 needs from here: **a hole is a cell below the nested grids' resolution whose
coarse parent is Irregular.** Regular and covered parents are not holes; F1 answers for them.

This also changes what the containment check at `Driver::regrid` line 668 should assert. Not "grids
stay inside the carried region" -- unmaintainable once physics drives the tags -- but "a grid outside the
carried region sits only over holes with cuttable parents, never over nothing". #731 cannot satisfy that
alone; the assertion should fail there until #732 lands.

## Plan of attack

1. Build the properly nested `Vector<DisjointBoxLayout>` over the cut cells in `ComputationalGeometry`,
   from the tags `getCurvatureTags` already produces, through `TiledMeshRefine` as it already does. No
   EBIS contact.
2. Carry all covered regions on every level as large boxes refined whole, disjoint from the nested grids.
   Regular gets nothing.
3. Hand the result to the shop as the layout, replacing `buildFinerLevels` below the scan level. By F5
   that is the whole integration.
4. Define and record the hole -- below the nested resolution, Irregular parent -- so #732 fills against
   a definition.
5. Move the containment check to line 668 with the meaning above.

Not yet looked at: the seam at the scan level, where `buildCoarseLevel`'s full coverage stops and the
nested grids start. The two must agree there.

## The downward sweep in ScanShop, and why the clip is not quadratic

An alternative to step 3 that stays inside ScanShop: since every rank already holds every box, a
downward sweep that grows and coarsens the irregular grids onto each parent level, then clips the
parent's regular/covered boxes against the result. The worry is the clip: every regular/covered box on
a level against every cut-cell box, which at 12 levels is 370k against 149k.

It is not that. The forced set on level L is `coarsen(grow(fineIrreg, 1), 2)`, and fine irregular boxes
come *only* from splitting Irregular parents -- a Regular or Covered parent is refined whole and its
refinement is pushed with the parent's tag, never rescanned (`localRegularBoxes.push_back(fineBox)`). So
`fineIrreg` is contained in `refine(coarseIrreg, 2)` exactly, growing by one fine cell reaches at most
half a coarse cell past it, and coarsening puts the forced set inside `grow(coarseIrreg, 1)`. Therefore:

**a regular or covered box can be hit by the forced set only if it is a neighbour of an irregular box
on the same level, within one cell.**

That is precisely the relation `computeNeighbors` builds, with the `grow(box, 1)` reach, by a windowed
x-sorted sweep. The clip set is the irregular-to-non-irregular adjacency, which scales with the surface
of the irregular region in boxes, not with the product of the two counts. The measured spill -- 2 to 6
boxes per level -- is what that adjacency looks like on the profiled plane.

The argument holds for the cumulative sweep too, since the grown irregular set contains the original
and the containment is stated against whichever set is current. It does depend on `m_neighbors` being
computable for every box, which is the replicated-classification point above: the adjacency needs to
know which boxes are irregular, and today only the owning rank does.

## The compounding, made precise

Take the realistic starting point: irregular grids on consecutive levels contained in each other with
coinciding boundaries -- zero margin -- which is exactly what `buildFinerLevels` produces, since the
fine irregular boxes are splits of the coarse ones. Each level must grow by its ghost width `g` and
restrict onto the level below, which then grows by its own `g`. In each level's own cells, the excess
beyond the original irregular set is

```
g/2 + g = 1.5g,   0.75g + g = 1.75g,   0.875g + g = 1.875g,   ...  ->  2g
```

It compounds from `g` to `2g` and stops there, because the coarsening halves what is inherited every
level. It does not run away. But `2g` is the number that matters, and the consequence is this:

**the boxes a level's grown set can touch are those within `2g` cells of its original irregular region,
which is `ceil(2g / boxSize)` boxes deep.** At `g = 4, box = 8` that is one box and the original
neighbours suffice. At `g = 4, box = 4` -- the run where the nesting violations appeared -- it is two,
and the original neighbours miss the outer ring. So `m_neighbors` as `computeNeighbors` builds it, at
reach 1 from the original layout, is insufficient on both counts: wrong reach, and stale after the
first amendment. The practical constraint is boxes no smaller than `2g`, which is stronger than "larger
than the ghost width".

## Finding the overlaps: the dual BVH walk

The clip on each level -- the parent's regular/covered boxes against the freshly grown irregular boxes
-- must not be `O(N_irreg x (N_reg + N_cov))`; at 12 levels that is 149k against 370k. Two
observations fix it.

First, the grown boxes do **not** sit on the `maxGridSize` lattice. The original partition does, but
growing by `g` cells puts the grown boxes at arbitrary offsets, so any lattice-keyed scheme has to snap
them outward, at a cost of up to `maxGridSize - 1` cells per side. Without snapping, index the boxes
themselves.

Second, EBGeometry is a submodule and has what is needed. `TreeBVH<T, P, BV, K>` over arbitrary
primitives with a bounding volume, `AABBT<T>` with `intersects`, SAH builders (`SAH2WaySplit`,
`SAHKWaySplit`) that handle disparately sized primitives, `PackedBVH` for a flat copy, and node access
-- `getBoundingVolume()`, `getChildren()`, `getPrimitives()`, `isLeaf()`, `getChildOffsets()`.

**The algorithm.** Build one tree over the grown irregular boxes and one over the regular/covered
boxes, both with `P = Box` and `BV = AABBT<double>`. Then walk the two together:

```
walk(a, b):
  if (!a.bv.intersects(b.bv))       return
  if (a.isLeaf() && b.isLeaf())     test primitives pairwise with Box::intersectsNotEmpty; return
  if (a.isLeaf())                   for c in b.children: walk(a, c); return
  if (b.isLeaf())                   for c in a.children: walk(c, b); return
  descend the larger volume:        for c in children(larger): walk(c, other)
```

Cost is proportional to the number of overlapping *node pairs*, `O(output + depth)` for two spatially
coherent sets. A whole subtree of regular/covered boxes far from the surface is dismissed by a single
root-level `intersects`, not once per box -- which is where the big empty boxes are. No adjacency, no
alignment assumption, no replicated classification needed to find overlaps: the trees are built from
box lists every rank already holds. The single-tree `traverse()` EBGeometry provides is query-shaped
(prune predicate plus leaf evaluator) and would give the "index one, query with each of the other"
form instead; the dual walk is not provided and is about thirty lines on top of the node access.

**The float question, and why it is not a narrowing.** `T` must be floating point
(`static_assert(std::is_floating_point_v<T>)`). Box corners are integers, and every integer below 2^53
is exact in a `double`; at 500k cells a side the corners are ~10^6. The conversion is lossless.

**But `AABBT::intersects` is strict:** `lo < other.hi && hi > other.lo`, with the comment "touching edges
are not overlapping". That decides the construction. Build the AABB from `{smallEnd, bigEnd + 1}` --
the half-open extent -- and the strict test on integer-valued doubles *is* `Box::intersectsNotEmpty`
for cell-centred boxes, on every axis:

- `[0,7]` and `[8,15]`: AABBs `[0,8)` and `[8,16)`, `8 < 8` false, correctly not overlapping;
- `[0,8]` and `[8,15]`: AABBs `[0,9)` and `[8,16)`, overlapping, correctly detected.

Built from `bigEnd` inclusive instead, `[0,8]` against `[8,15]` gives `8 > 8` false and a real overlap
is pruned, silently. The half-open construction is required by the strict test, not a preference.

So the node-level prune is exact rather than conservative-with-slack, which matters in a dual walk
because a loose prune costs node-pair visits at every level. The `intersectsNotEmpty` at the leaves is
the authoritative statement of what is meant and insurance against the AABB ever being rebuilt from
inclusive corners; it costs integer compares on pairs that already passed an exact filter.

**Keep the primitive cell-centred.** Converting the stored boxes to node-centred gives the same AABB
numbers -- `surroundingNodes([0,7])` is `[0,8]` -- so it changes no prune decision, and it costs two
things. `intersectsNotEmpty` on node-centred boxes reports touching cell boxes as overlapping (they
share the node plane), so the leaf test would need an `enclosedCells` first to be correct. And the
decimation operates on cells, so every clip would convert back, with the half-open/closed distinction
reappearing there where it is harder to see than in one constructor. Face-centred is direction-dependent
and strictly worse. The primitive is the cell-centred `Box`, unchanged; the AABB is a derived filter
built once from `{smallEnd, bigEnd + 1}`.

**Two-dimensional builds:** `intersects` is hardcoded to three components and asserts `lo <= hi` on
each. Fill the third axis with `[0, 1)` on every box so it always overlaps and never asserts.

**What the walk does not solve:** the sliver. `R \ F` with `F`'s edge at an arbitrary offset inside `R`
can leave a remnant a cell or two thick, and no alignment fixes that once the grown boxes are off the
lattice. The rule at clip time is to absorb any remnant with a dimension below the minimum box size into
the adjacent grown box -- trading a sliver of regular area for irregular, the cheap direction -- which
keeps every box on the level at or above the minimum and is bounded by the same `2g` argument, since
absorption never grows anything by more than a sub-minimum dimension.

## Whether the decimation is needed at all (superseded -- see the next section)

The argument below was rejected. Two things it misses: (1) restricting the finer irregular boxes to the
level below does not produce a superset of that level's irregular boxes but a *different* tiling --
`2^D` boxes of size `M/2 + g` under each coarse box of size `M`, overlapping each other and the parent --
so the irregular sets of adjacent levels do not agree and the level has to be re-tiled; and (2) once the
level is re-tiled, a tile that straddles the old irregular region and a regular/covered box has one tag,
so the regular/covered box must give the tile up. Decimation is forced by the tiling, not by geometry.
The section is kept because the geometric facts in it (the band is regular/covered by construction)
still hold and still matter: they are what lets the carved tiles inherit their tag instead of being
reclassified. The design that replaces it is in the next section.


**What the buffer has to be.** The nesting requirement is that level L-1 be *generated* in a band at
least `g` wide around L's footprint, so that nothing adjacent to L is served from L-2. Generated means
present in L-1's layout with a correct classification. It does not mean tagged Irregular: a Regular tag
is a correct generated answer for a regular cell.

**What the band is made of.** By containment, every cut cell of L-1 near L's footprint is already
inside an L-1 Irregular box -- a fine irregular box comes only from splitting an Irregular parent, so
the surface never passes through a box classified Regular or Covered. The part of the `g`-band that
spills onto regular/covered boxes is therefore, by construction, geometrically regular or covered. The
spill measurement showed exactly that: 20, 10 and 10 cells landing on regular/covered boxes, none cut.

So for the spill:

- **regular** -- needs nothing; it is already correct and presence is free;
- **covered** -- must be carried, which is the constraint already established (all covered regions,
  every level, refined whole);
- **cut** -- is already in an Irregular box, **unless `retainBox` pruned that box's refinement.**

That last case is the only work: **un-prune Irregular refinements within `g` of the finer footprint.**
Those are splits of coarse Irregular boxes, already tile-aligned, already the right size, already
classified by the existing path. No regular or covered box is touched, so there is nothing to decimate,
no sliver to absorb, and no dual walk to run.

**`retainBox` already contains a `g`-buffer.** `grow(coarsen(a_box, 2), m_ebGhost)` grows the parent by
the ghost width before testing it against coverage, so retention is meant to keep exactly that band.
The nesting violations at block size 4 are then not a missing mechanism but a quantisation: retention is
decided per whole coarse box, against coverage regions that may themselves not be buffered-nested at
that tile size. A smaller problem than re-partitioning levels.

**The question that decides it:** does anything downstream require the *irregular boxes* to be nested,
as opposed to the *carried region*? Something that walks irregular boxes level to level and assumes a
parent box for each child box would need box nesting. The multichord does not -- it needs the finer
level's children of a seam cell, which is containment, not a buffer. `coarsenFrom` does not -- its ring
needs generated neighbours of any tag. If nothing does, the decimation, the slivers, the snapping and
the BVH are answering a question nobody asked. Unresolved; listed below.

**If regular/covered boxes do have to be cut, snap to the simulation's tile size.** Three reasons. The
box is the unit of graph cost -- an Irregular-tagged box allocates its `BaseFab<GraphNode>` over its
whole region however few cells are cut -- so once the minimum size is `2g` or more, free-form clipping
plus sliver absorption converges on the same box count as snapping, and snapping is deterministic.
Tile alignment with the simulation's grids makes `fillEBISLayout`'s copies box-to-box where they
coincide rather than fragmented. And `TiledMeshRefine` already made this choice for the same problem.
The overshoot is at most one tile per side, the quantisation `retainBox` already pays.

## The repartition: `TiledMeshRefine` is the tiler, and the decimation it leaves

**Why the restricted sets disagree.** One coarse Irregular box `B` of size `M` at level L-1 has `2^D`
fine Irregular boxes of size `M` under it at level L. Restricting each -- `coarsen(grow(fine_i, g), 2)`
-- gives `2^D` boxes of size `M/2 + g` at L-1, overlapping each other, overlapping `B`, and smaller than
`M`. That is not a superset of L-1's irregular set; it is a finer, non-disjoint re-tiling of nearly the
same region plus a rim. The union has to be re-tiled, and the re-tiling has to be canonical (disjoint,
sized, aligned) or a `DisjointBoxLayout` cannot be built from it. This is the repartition, per level.

**`TiledMeshRefine` already is that repartition** (`Source/AmrMesh/CD_TiledMeshRefine.{H,cpp}`). Read
from source:

- Input: `regrid(Vector<Vector<Box>>& grids, const Vector<IntVectSet>& tags)`. Tags on level `l-1`
  produce tiles on level `l` (`makeLevelTiles` maps a coarse tag to the fine tile containing it,
  `iv = (tag - coarProbLo) / (tileSize / refToCoar)`, line 338 ff). Tags are rank-local; the only
  communication is `gatherSuperTiles`, which assumes a rank's tags are disjoint from every other rank's
  -- ScanShop's rank-owned irregular boxes satisfy that.
- Representation: `SuperTiles` -- `std::unordered_set<uint64_t> m_full` of full `max_block_size`
  super-tiles keyed by `encodeSuper` (21 bits per direction) plus `m_partial` key -> sub-tile bitmask.
  This is the Morton-dedup idea from the previous section, already written: a packed key and a hash
  set, materialising fine tiles only at the boundary of the region.
- Nesting: `nestFrom` (line 224) adds `coarsen(grow(fineTile, 1), ref)` for every finer tile, in bulk
  for full super-tiles, per sub-tile for partial ones. The buffer is therefore **one fine tile** of
  `tileSize` fine cells, i.e. `tileSize / ref` coarse cells. The nesting requirement
  `T_{L-1} ⊇ coarsen(grow(T_L, g))` is met whenever `tileSize ≥ g`, which `min_block_size ≥ 2g` gives.
  The compounding of the previous sections is absorbed: each level's buffer is taken from that level's
  finished tile set, which already contains the level above's buffer.
- Output: `makeBoxesFromTiles` -- one `max_block_size` box per full super-tile, `packTiles`
  (Berger-Rigoutsos on tiles, exact, never covers an untagged tile) for partial ones. Boxes are disjoint,
  tile-aligned, between `tileSize` and `max_block_size` per direction, and emitted in sorted-key order so
  **the list is identical on every rank** without a broadcast.
- Level 0 of its output is only the nesting buffer of level 1 (`regrid` line 300); `AmrMesh` discards
  it (`newBoxes[0] = oldBoxes[0]`, `CD_AmrMesh.cpp` line 1254). For EB generation level 0 would be the
  scan level, which `buildCoarseLevel` builds whole, so the same discard applies.

Fed with `tags[l] = I_l` (the Irregular boxes ScanShop's upward build produced at level `l`, coarsened by
the refinement ratio so they land on level `l-1` where `regrid` expects them -- coarsen-then-refine
overshoots by at most `ref-1` cells, which the tile snap absorbs anyway), it returns `T_l ⊇ I_l`,
properly nested with a one-tile buffer, curvature-adapted because `I_l` is (the coverage regions have
already pruned it through `retainBox`). No new type: the rule in `CLAUDE.md` § Data structures applies
and is satisfied by using this class rather than writing a Morton set.

Two gaps, both small:

- `regrid` consumes tags with `IVSIterator`, per cell. Tagging whole irregular boxes costs the
  irregular *volume* per rank (own boxes only) -- the same volume `fillGraph` iterates anyway, so it is
  not the quadratic term, but it is `M^D` work per box for nothing. `nestFrom`'s `addNest` lambda
  already is the box-to-tiles primitive (`BoxIterator` over the tile-coordinate box); a box-tag entry
  point is that lambda exposed.
- The EBIS `maxGridSize` must be a multiple of `max_block_size`, so that no tile box straddles two
  ScanShop boxes. ScanShop's boxes are `maxGridSize`-lattice-aligned (`domainSplit` at the scan level,
  refined whole above it, irregular boxes split to `maxGridSize`); with that constraint every tile box
  lies inside exactly one ScanShop box. The tag inheritance and the subtraction below both rely on it.

**The decimation, given `T_l`.** Per level, with `S_l` the ScanShop layout (replicated `m_boxes` plus
the rank-local `m_boxMap` tags, which need a gather of one `int` per box to be usable everywhere):

1. For each tile box `t ∈ T_l`, locate the ScanShop box `B ∋ t`. Point location in a disjoint box set:
   the BVH over `S_l` from the previous section answers it in `O(log N)`; the dual walk is not needed
   because one side is now a point query. `T_l` scales with the surface, so this is `O(|T| log N)`.
2. `t` **inherits `B`'s tag.** If `B` is Irregular, `t` is Irregular. If `B` is Regular or Covered, `t`
   is Regular or Covered -- correct without any geometric query, because a sub-box of a box that
   classified Regular/Covered over its `m_ebGhost`-grown region is Regular/Covered over its own grown
   region. This is where the previous section's geometric fact does its work: the band is
   regular/covered by construction, so nothing needs `isRegular`/`isCovered` here. Mark such `B` as
   hit.
3. For each hit `B`: remainder `= B \ ∪{t ⊆ B}`, tagged as `B` was. Done in super-tile coordinates
   (`B / max_block_size` is exact under the constraint above) with `TreeIntVectSet`: `define(B)` is one
   full node (`&full` sentinel, `TreeIntVectSet.cpp` line 1744), each `-= t` descends `O(depth)` and
   splits nodes (`remove`, line 1736), `createBoxes` emits one box per remaining full node -- an
   octree-graded decomposition, every box a multiple of `max_block_size`, `O(hits_B · depth)` boxes.
   `numPts()` is never called on this path (only `operator<` calls it, line 964), so the `int` overflow
   at large boxes does not bite; the static scratch vectors (`index_local`, `parents_local`, lines
   1069-1070) make it non-reentrant, which is fine serially per rank.
4. New layout of level `l` = `T_l` with inherited tags, plus remainders, plus untouched ScanShop boxes.
   Every rank computes the same list (replicated inputs, deterministic steps), then the existing
   `LoadBalancing` assigns it. Irregular tiles that were not Irregular in `S_l` do not exist -- an
   Irregular tile is one that lies inside an `I_l` box -- so `fillGraph` runs on exactly the cut region
   plus the tile quantisation, no more.

The cost has no term in `|R_l| + |C_l|` beyond building the BVH once per level: the forbidden
`O(irregular × (regular + covered))` does not appear. What is paid is the tile overshoot
(`≤ max_block_size` per side, as `retainBox` already pays) and the octree grading of decimated
remainders, which is the same shape `TreeIntVectSet::createBoxes` produces everywhere else in Chombo.

**What this settles from the open list.** The question "must irregular *boxes* be nested, or only the
carried *region*" is decided by construction rather than by audit: the tiles are nested, the carved
tiles carry honest tags, and nothing downstream is asked to tolerate a half-and-half box. It stays in
the list below only as the audit of what downstream assumes, which is still unperformed.

**Corrections and the open item (packing the remainder).**

- Step 1 should be the dual walk after all, not a point query per tile: the dual walk is
  `O(output + depth)` because neighbouring tiles share the descent, `|T| · log N` does not. `T_l` comes
  out of `makeBoxesFromTiles` in sorted-key order, so its BVH is a build over an already coherent list.
- Steps 2-3 are the decimation. What is not settled is how tightly the remainder `B \ hits` is packed.
  The remainder is Regular/Covered, so tightness is box *count* (one `int` per box on the wire, one
  entry per box in every layout loop), not storage.
- `TreeIntVectSet::createBoxes` is not tight. It splits only at midpoints: `B` minus a one-super-tile
  slab on one face is ideally one box, the octree gives about `2^(D-1)` boxes per level of depth
  (`~4 · depth` in 3D, ~40 boxes for a `B` ten super-tiles deep). Graded, exact, but paying `O(depth)`
  boxes for a shape that needs one.
- Tight means splitting at an empty slab before splitting at the middle. That is the Berger-Rigoutsos
  rule `packTiles` already implements (`CD_TiledMeshRefine.cpp` line 420 ff: emit when the node is
  fully tagged and under the cap, else cut at a zero-signature slab nearest the middle, else at the
  strongest inflection, else at the midpoint). It runs on an explicit tile list, which the complement
  of the hits inside a large `B` cannot be -- enumerating it is `O(volume of B)`.
- The rule does not need the set explicit; it needs signatures and two predicates. Along an axis the
  complement's signature is `slab area − hit signature`; a node is fully complement iff it contains no
  hit; fully hit iff its hit count equals its volume. All three come from the hit boxes alone (each hit
  is a box of tiles, so its slab contribution is a range, `O(1)` per slab per box). BR on that implicit
  complement is `packTiles` with a different oracle: `O(hits · log)` per recursion level, never touches
  the volume of `B`, one box for the bulk when the hits form a slab, and a tightly packed layer where
  they do not.
- The cheap partial answer -- peel `B \ minBox(hits)` as at most `2D` exact slabs, then enumerate the
  complement inside `minBox(hits)` -- works when the hit layer is thin and fails when the surface
  crosses `B` corner to corner (`minBox(hits) ≈ B`). The recursion is needed for that case, which is
  the implicit-BR above.
- Hits per `B`: a big regular box next to the surface is hit along the faces that face it, so hits
  scale with `B`'s face area in tiles, and the remainder layer can leave `O(hits)` small boxes where
  the layer is bumpy. Chombo's `MergeBoxesOnLines` (`BoxTools/MergeBoxesOnLines.H`) merges a box list
  along one direction as a post-pass. Whether that layer is acceptable, or the tiles should be coarsened
  to super-tiles for the purposes of decimation to keep the layer flat, is the user's call.

**Decision: record the intersections, decimate last.** With `R`, `C` the regular and covered box
sets and `I` the irregular set, the dual walk records, for each `r ∈ R ∪ C`, the list of boxes of `I`
(the tiles) that intersect it. Nothing is cut during the walk. The decimation runs once, at the end, over
those lists -- so the packing algorithm (octree, implicit Berger-Rigoutsos, peel-and-pack, or none) is a
separate step that can be replaced or tuned without touching the walk or the tiler. Tightness is
deferred; it may not be needed.

**Where this is built.** The grid set -- properly nested irregular tiles, curvature-adapted, with `R`
and `C` decimated against them so every level is fully classified -- is built in `ComputationalGeometry`,
not in ScanShop. It will resemble ScanShop's upward build but is written from scratch with curvature in
it from the start. That re-scopes the PR stack: #731 becomes the mesh-building PR; #732, rebased later,
becomes the EBIS-generating PR that consumes those grids; a third PR on top takes the regrid path. The
work recorded in this document from #731 (surface export, generator keep-alive, the 3D restricted-face
seam machinery, this file) moves to #732, and #731 is reset to #730's head to start the mesh builder.

## The two phases

**The machinery does not couple them; the simulation grids do.** Chombo's own `MFIndexSpace::define`
would: it calls `reconcileIrreg` and `levelStitch` at every level, and both walk the two phases'
`m_grids` with two `DataIterator`s in lockstep and assert `a_otherPhase.m_grids[ditb] == region`
(`EBISLevel.cpp` lines 2032-2047 and 2093-2106). chombo-discharge does not use that path:
`MultiFluidIndexSpace::define` (`Source/Multifluid/CD_MultiFluidIndexSpace.cpp` lines 19 and 37)
defines the two `EBIndexSpace`s independently and never stitches. So nothing forces the phases onto
one layout. What couples them is that the simulation grids are shared: every box the simulation puts
on level `l`, plus ghosts, must be carried and correctly classified by *both* index spaces at `l`.

**Per-phase hierarchies break that.** The phases see different surfaces: gas IF = electrodes ∪
dielectrics, solid IF = (complement of dielectrics) ∩ electrodes (`CD_ComputationalGeometry.cpp`
lines 555-561). An electrode buried in a dielectric is a solid-only surface, an electrode in gas is
gas-only, and at a triple junction the two surfaces differ, so a curvature criterion stops at
different levels in the two phases and leaves different holes. A simulation box inside the gas's
carried region can land in a solid hole, which reads back `AllRegular` -- fluid where there should be
covered dielectric. The gas-side grid was fine; the failure is silent and on the other phase.

**Two upward builds, one tile hierarchy, two classifications.**

1. Run the upward build **per phase**, each on its own implicit function with its own curvature
   stop. This yields `R_l^p, C_l^p, I_l^p` per phase `p` and per-phase holes above each stop.
2. Tag the tiler with the **union across phases**, `tags_l = I_l^gas ∪ I_l^solid`, and run
   `TiledMeshRefine` **once**. `T_l` is common to both phases, properly nested, and is also the
   coverage the simulation regrids onto -- one set, which is what "same coverage" has to mean.
3. Classify `T_l` **per phase** against that phase's own boxes. The inherit rule gains a third case:
   - tile inside an `I^p` box: Irregular;
   - tile inside an `R^p` or `C^p` box: inherits the tag;
   - tile inside a **hole** of phase `p` (the other phase asked for this level here, `p`'s own build
     had stopped): classify geometrically on the grown tile, and fill the graph if it is cut. Its
     parent chain exists -- the level below is an `I^p` box by construction, and `T` is nested so
     every deeper tile also has a tile beneath it.
4. Decimate each phase's `R^p`/`C^p` against the common `T_l`, intersections recorded during the walk
   and cut at the end.

The outcome is the same boxes on the irregular tiles in both phases, a different decomposition of the
regular and covered remainder in each, and consistent classifications throughout. A phase carries
tiles the other asked for, but where those land on `R^p`/`C^p` they are tag-only, and where they land
in a hole they are cut cells the simulation was going to put a grid on -- the geometry that phase
would otherwise have served wrong. The union decides only *where a level exists*, not what it holds;
nothing is refined in a phase because of curvature the other phase has. At the shared dielectric
surface the two criteria agree up to box quantisation and the union takes the finer stop.

This also makes explicit what the curvature stop means: not "do not generate the surface finer than
this" -- a cut cell that curvature is happy with at level 3 is still cut at level 4 -- but "do not let
the *simulation* go finer than this here". Level 4 can be a hole there only because no grid will ask.
So the regrid PR must clip the physics tags to `T_l`, and the builder hands Driver `T_l` as the
coverage for both phases, not per-phase sets.

## The mesh builder in `ComputationalGeometry`: the consolidated plan

Built before any `EBIndexSpace` exists, on the implicit functions alone. Its product is box lists
with classifications, per phase and per level, kept as members: the shops read them through
`makeGrids`/`InsideOutside`, Driver reads the tiles as its coverage, and the irregular set will be
shown as a simulation grid for inspection later. All level indices here are AMR-indexed, coarsest
first; the record quotes EBIS code elsewhere with the opposite convention.

### Inputs

- coarsest domain, its grid spacing, and `probLo` -- the implicit functions are evaluated at physical
  points, so the builder needs the resolution, not only the index space; refinement ratio 2
  throughout;
- the start domain, from Driver (the domain ScanShop calls the scan level): the finest level built whole, level 0 here;
- the stop domain, the finest level the upward pass may reach whatever the curvature says. A domain rather
  than a depth: `max_amr_depth` counts from AMR level 0 and a depth counted from the start domain would read
  the same way, so two domains say it without a convention;
- the simulation's `min_block_size`/`max_block_size`: tile and super-tile for the tiler, and also the
  size every box the upward pass makes is split to. The EBIS `maxGridSize` does not enter: ScanShop
  used it for the start-level split and for splitting a refined irregular box, and using
  `max_block_size` for both puts every box the builder makes on the super-tile lattice, so the union
  in step 2 is exact in tile units and a tile lies in exactly one box by construction;
- `m_maxGhostEB`; `refine_angles`; the two implicit functions.

### State, per phase `p` and per level `l`

`Vector<Box>` regular, covered, irregular. Per level, common to both phases: `T_l`. After decimation,
per phase per level: the final box list with a parallel tag list. Nothing is a `DisjointBoxLayout`
until a shop load-balances it. The lists persist for the life of the object.

### Step 0 -- the start level and below

`domainSplit` of the whole domain to `max_block_size`; every box classified by cell-centre values on the
box grown by `m_maxGhostEB`: regular iff every value `< -halfDiagonal`, covered iff every value
`> halfDiagonal`, irregular otherwise, with the scan skip that assumes signed distance
(`ScanShop::isRegular`/`isCovered`, `CD_ScanShopImplem.H` lines 43-100, made callable on an implicit
function). Every level at or below the start level is whole in both phases; no hole exists there.

### Step 1 -- upward, per phase, to `the stop domain`

From level `l` to `l+1`, for each box at `l`:

- **regular or covered:** `refine(box, 2)`, pushed whole with the same tag. No scan, no split. This
  is what guarantees that nothing beneath a regular or covered box is ever a hole.
- **irregular:** the box-level curvature test, early exit. Cut cells do not exist at this point and
  are not looked for; the test runs on the implicit function. Over the cells of the box grown by one
  cell, take the cells in the band `|f(x)| ≤ √D · dx` at the cell centre (the scan skip prunes the
  rest), evaluate the normal there by central differences of `f` (no implicit function in the tree
  implements `BaseIF::derivative`; the base throws), and compare each band cell's normal with those of
  its band neighbours. The first pair whose angle exceeds `refine_angles` ends the scan: the box
  **splits** -- `refine(box, 2)`, `domainSplit` to `max_block_size`, each piece classified as in step 0
  -- and the scan moves to the next box. A box no pair in fails is a **leaf**: nothing is pushed
  beneath it, which is the hole above it.
- run to the stop domain regardless: once no irregular box remains the curvature test has nothing to do,
  but regular and covered boxes still refine whole onto every remaining level. A level above the last
  split is partial (regular/covered descendants, tiles where curvature reached, holes elsewhere), never
  missing.

Across a sharp edge the angle is the dihedral angle at every `dx` and never shrinks, so edges refine
to `the stop domain`. That is intended (it is where the multichord seam lives) and it is what
`refine_angles` does today off EBIS normals; `the stop domain` bounds the depth, `refine_angles` does
not on a geometry with edges. A difference across a kink of a min/max composite averages the two
faces, which reads as a jump against either neighbour: the same behaviour from the other side.

Work is distributed as ScanShop distributes it: each rank classifies a share of the boxes, the lists
are gathered and every rank holds all of them in a deterministic order.

### Step 2 -- the tiles, once, after both phases have finished step 1

`tags_l = I_l^gas ∪ I_l^solid` for every level above the start level. One `TiledMeshRefine` with the
start-level domain as its coarsest, ratio 2, tile `min_block_size`, super-tile `max_block_size`; each
irregular box is a super-tile already, so it enters exactly, as `coarsen(box, 2)` on the level below
or through a box-tag entry (`nestFrom`'s `addNest` lambda is that primitive). Its level 0 is the
start level, already whole, and is discarded
as `AmrMesh` discards it. The buffer is one tile per level, `≥ m_maxGhostEB` when
`min_block_size ≥ 2 · ghost`. Output `T_l`: disjoint, tile-aligned, identical on every rank. This is
the downward sweep; `regrid` descends internally and injects the buffer level by level.

No collar is added beyond the union. The simulation's ghost cells may reach past `T` into a hole;
answering for those cells is the index space's job (below), not the builder's.

### Step 3 -- the walk, per phase

BVH over the phase's `R_l ∪ C_l ∪ I_l` and BVH over `T_l` (`EBGeometry`'s `TreeBVH`, half-open
AABBs `{smallEnd, bigEnd + 1}` so the strict test equals `intersectsNotEmpty`); dual walk. Every box
the builder made is a super-tile or a whole refinement of one, so a tile meets at most one box. Every
tile at level `l`
lies in exactly one of three places, and the walk records which:

1. **inside an irregular box** of this phase at `l`: tag Irregular;
2. **inside a regular or covered box** of this phase at `l`: tag inherited, and the tile is appended
   to that box's list of intersecting tiles;
3. **inside the refinement of a leaf** of this phase at `l-1` -- a hole: no box to inherit from, so
   the tile is classified by step 0's test on the grown tile. Its parent is irregular, so any of the
   three answers can come out.

Case 3 arises within a single phase as well as across phases: the buffer from `T_{l+1}` reaches into
the refinement of a leaf at `l-1` whenever a box next to that leaf split twice. There is no fourth
place, because `T` is nested and every level at or below the start level is whole. Every tile's
answer is appended to the phase's `R`, `C` or `I` list at `l`, so each box the shop is handed has a
tag.

### Step 4 -- decimation, per phase, last

For each regular or covered box with a non-empty tile list, remainder `= box \ tiles`, tagged as the
box was. First implementation: `TreeIntVectSet` in tile coordinates (`define(box)` one node, `-=` per
tile, `createBoxes`), correct, octree-graded, not tight. The tile lists recorded in step 3 are the
interface any packer consumes, so replacing this step touches nothing else. The final list at `l`
for the phase is `T_l` with its tags, the remainders, and every box no tile touched.

### Step 5 -- the coarser levels, and irregularity pushed down

ScanShop builds every level at and below its scan level whole -- `domainSplit`, each box classified on its
own -- and no level below the scan level looks at any other (`CD_ScanShop.cpp` lines 147-156, 200-272).
The builder does the same for the levels coarser than the start domain, down to the coarsest domain that
can still be coarsened by two, with `max_block_size` as the split; no tiles, no decimation. Levels are
indexed from that coarsest domain, and the shop looks its domain up with `getLevel`.

Then the reconciliation ScanShop does not do: from the finest level down, **a box that contains a finer
irregular box is irregular**. On its own a whole level's classification agrees with the level above only
because the test is conservative for a signed-distance function; this makes it hold by construction, for
the start level as well as the coarser ones. One containment per finer irregular box: `coarsen(b, 2)` is at
most half a tile wide on the shared lattice and lies in exactly one box. The box it flips is a tile or a
super-tile, never a decimated remainder -- above the start level a finer irregular tile lies inside a tile
of this level by nesting. Runs after decimation, because the push-down reads the final lists.

### Step 6 -- what is handed over, and to whom

- To the shop (#732): per phase, per level, the final list and tags. `makeGrids` load-balances it;
  `InsideOutside` returns the recorded tag; `fillGraph` runs on Irregular boxes only.
- To Driver (the regrid PR): `T_l` as the coverage for both phases, and the physics tags clipped to
  it.
- To the index space (#732): the contract the fill-from-parent path relies on when a simulation box's
  ghost region, or a later regrid, reaches an uncarried cell:
  1. a hole exists only above an **irregular leaf of the same phase**; nothing beneath a regular or
     covered box is ever a hole, and the start level and below are whole, so a hole cell's parent
     chain is irregular down to a generated level;
  2. the leaf is generated and its surfaces stored, so the hole is filled by cutting them, at a
     resolution bounded by `the stop domain`;
  3. what the builder guarantees is the lists above; what it does not guarantee is that `T` covers
     the simulation's ghosts. To check in #732: that the fill path is wired to
     `EBISLayoutImplem::define`'s copy for a partial region (a ghost ring), not only to whole-level
     construction (`3ffd27464` cuts a level).

### Constraints on the inputs

`min_block_size ≥ 2 · m_maxGhostEB`; the start level whole in both phases; the signed-distance
assumption behind the scan skip, already made by ScanShop.

### Left open

Whether steps 3-4 run replicated on every rank (deterministic, no communication, `O(N log N)` per
rank) or distributed by box and gathered; which ranks tag which boxes for step 2's gather once the
lists are replicated (`gatherSuperTiles` is correct with duplicates, but every rank tagging every box
multiplies the gather by `nprocs`; partition the replicated list by index); the finite-difference
step relative to `dx`; what an empty level's list means to the shop (a phase whose deepest tile is at
`l` still answers `makeGrids` down to the finest domain); the packer.

## Claimed but not established

- That dropping covered regions is harmless in practice. A run on ProfiledSurface leaves 982776
  uncarried cells whose centres are inside the solid at the finest level and completes -- but that
  geometry never asks about them, so it shows nothing. A geometry with a covered interior adjacent
  to the grids, with the covered collar deliberately dropped, is the test that could fail.
- Whether anything downstream iterates `EBISLevel::m_grids` expecting the level to tile the domain.
  `EBCoarseFineParticleMesh::defineStencils` walking every irregular cell of the coarse level was one
  such surprise already; the consumers have not been audited.
- What downstream assumes about the irregular boxes level to level. The design above nests them by
  construction, so this no longer decides whether regular/covered boxes are cut; it is still the
  audit of consumers that has not been done.
- That the tile-inherits-tag step is sound when the EBIS `maxGridSize` is a multiple of
  `max_block_size`. Argued from ScanShop's lattice alignment; not checked against a run with the two
  parameters unequal.
- Why the coverage regions fail buffered nesting at block size 4 when they come out of
  `TiledMeshRefine`, whose `nestFrom` is meant to enforce it. Measured, not explained.
- The seam at the scan level between `buildCoarseLevel`'s full coverage and the nested grids.
- That the dual walk is `O(output + depth)` on these box sets. The argument is the standard one for
  spatially coherent BVHs; it has not been measured on the 12-level hierarchy.
- Whether anything downstream assumes `m_neighbors` was built by `computeNeighbors` with the
  `grow(box, 1)` reach. A list supplied through `closeN` with a different reach would change ghost
  exchange silently.
- Whether `EBData` short-circuits on the tag the way `EBGraph` does. Still unchecked, and it decides
  whether a carried covered box is as cheap as the graph analysis suggests.

## The seam, for when the grids are settled

Unrelated to the above and already measured. On ProfiledSurface with four levels: no multi-valued
cells at any level, the surface manifold within each level, 5796 open edges, every one of them on a
level-0 face plane -- the coarse-fine boundaries. Only 156 of 5952 per-level rim edges stitch. That
is one coarse chord against four fine sub-chords, which is what `8d71c2b48` was staged to fix.
