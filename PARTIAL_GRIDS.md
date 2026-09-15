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

## Claimed but not established

- That dropping covered regions is harmless in practice. A run on ProfiledSurface leaves 982776
  uncarried cells whose centres are inside the solid at the finest level and completes -- but that
  geometry never asks about them, so it shows nothing. A geometry with a covered interior adjacent
  to the grids, with the covered collar deliberately dropped, is the test that could fail.
- Whether anything downstream iterates `EBISLevel::m_grids` expecting the level to tile the domain.
  `EBCoarseFineParticleMesh::defineStencils` walking every irregular cell of the coarse level was one
  such surprise already; the consumers have not been audited.

## The seam, for when the grids are settled

Unrelated to the above and already measured. On ProfiledSurface with four levels: no multi-valued
cells at any level, the surface manifold within each level, 5796 open edges, every one of them on a
level-0 face plane -- the coarse-fine boundaries. Only 156 of 5952 per-level rim edges stitch. That
is one coarse chord against four fine sub-chords, which is what `8d71c2b48` was staged to fix.
