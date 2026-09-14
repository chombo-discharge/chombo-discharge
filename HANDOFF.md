# Cut-cell geometry — handoff

Reference for picking this up cold. `PLAN.md` is the design and Milestone 1; `PLAN_MILESTONES.md`
is Milestones 2 and 3. This file is the part that is hard to reconstruct: what is built, what is
measured, what has bitten, and what to do next.

**Where the work now lives.** These notes were written during Milestone 1 and describe it as the
whole of what exists. Milestone 3 has since been built as well. Read them as the record of how the
design was arrived at and what it cost to find, not as a description of the current tree; the pull
requests are that.

| | | |
| --- | --- | --- |
| #730 | `polyhedral-geometry-shop` | Milestone 1, the moment machinery. What §1–§6 below describe |
| #731 | `cutcell-refinement` | Milestone 3, generation over real AMR grids. Built; §7 item 4 is done |
| #732 | `polyhedral_refinement` | Milestone 2, refining a cell by cutting its body. Carries these notes |
| Chombo-3.3#31 | `generate-every-level-from-geometry` | the index-space side of #731 |

---

## 1. What exists right now

`Driver.geometry_generation = polyhedral` selects a third geometry generator alongside
`chombo-discharge` and `chombo`. It takes the same topology as `GeometryShop` — the same
classification, the same graph, the same covered-cell fix-ups — and replaces only the quadrature.
Moments become exact integrals of a closed polyhedron instead of quadratures of the implicit
function.

| file | what it is |
| --- | --- |
| `Source/Geometry/CD_CutCellSurface.{H,Implem.H}` | the values a cut cell is reconstructed from: corner values and edge crossings |
| `Source/Geometry/CD_CutCellBody.{H,cpp}` | assembles the polyhedron and takes its moments; the closure and range checks |
| `Source/Geometry/CD_PolyhedralEBUtils{.H,Implem.H}` | the `PolyhedralEB` namespace: the fluid predicate, the normals, and the `detail` helpers the bodies share |
| `Source/Geometry/CD_PolyhedralGeometryShop.{H,cpp}` | the generator; classification, edge cache, node filling, the dropped-cell sweep |
| `Source/Geometry/CD_ComputationalGeometry.{H,cpp}` | `usePolyhedralShop`, and the third branch that builds it |
| `Source/Driver/CD_Driver.{H,cpp,options}` | the switch |
| `Prototypes/CutCellRefinement/main.cpp` | the CSV exporter every measurement here came from |
| `Prototypes/CutCellRefinement/spike/refine_spike.cpp` | standalone C++ body construction **and clipper**; Milestone 2 starts here |

**Nothing in Chombo was changed for this.** That is deliberate and it is possible because
`GeometryShop.H` is `public:` from the top of the class through line 121 and `private:` after, and
that line falls in exactly the right place: `fillGraph`, `computeVoFInternals`,
`getFullNodeWithCoveredFace` and `fixRegularCellsNextToCovered` are reachable from a subclass;
`fillNodeValues`, `insideOutsideFromNodes`, `edgeData2D/3D` and `BrentRootFinder` are not. The
unreachable ones are exactly the parts the prototype already reimplements.

One Chombo change *was* made, separately and for its own reasons: Chombo-3.3#30, merged, which the
submodule pin moves to. See §5. Milestone 3 needs Chombo changes of its own and has them in
Chombo-3.3#31; that does not weaken the claim here, which is about Milestone 1.

## 2. How to run it

```bash
cd Exec/Tests/Geometry/MechanicalShaft
make -j DIM=3 MPI=TRUE OPT=HIGH
mpirun -np 8 ./main3d.*.ex regression.inputs \
    Driver.geometry_generation=polyhedral \
    Driver.geometry_only=true
```

Writes `geo/simulation.geometry.3d.hdf5` as usual. Verified in 2-D and 3-D.

The measurement harness is the prototype exporter, which writes one CSV row per cut cell with the
moments, the edge crossings and the corner values:

```bash
cd Prototypes/CutCellRefinement
make -j DIM=3 MPI=TRUE OPT=HIGH
mpirun -np 4 ./main3d.*.ex sphere3d.inputs \
    Driver.geometry_generation=polyhedral \
    Prototype.geometry=polyhedron Polyhedron.shape=cube \
    Polyhedron.angles="31.7 19.3 47.1" Polyhedron.size=0.4013
```

`Prototype.geometry` takes `plane | sphere | swept_sphere | torus | polyhedron | tessellation |
rough_sphere`; `Polyhedron.shape` is `cube | simplex` with Euler angles. Run it twice with
different `Driver.geometry_generation` and compare the CSVs cell by cell — that is how every number
below was obtained.

## 3. What is measured

Generator against generator, same runs, cell by cell. **Topology differences are zero everywhere**;
moments differ only where the design says they should.

| geometry | cells | topology | Δ volFrac | Δ areaFrac |
| --- | ---: | ---: | ---: | ---: |
| oblique plane 3-D | 8,898 | 0 | 3.8e-15 | 8.6e-15 |
| sphere 3-D | 1,220 | 0 | 4.9e-4 | 8.4e-2 |
| torus 3-D | 4,208 | 0 | 2.3e-3 | 7.8e-2 |
| rotated cube 3-D | 6,307 | 0 | 5.8e-2 | 4.9e-1 |
| rotated simplex 3-D | 9,595 | 0 | 1.5e-1 | 5.2e-1 |
| axis-aligned simplex 3-D | 9,852 | 0 | 1.2e-1 | 3.4e-1 |
| **axis-aligned cube 3-D** | 4,376 | 0 | **0.00e+00** | **0.00e+00** |
| oblique line 2-D | 105 | 0 | 3.6e-15 | 5.8e-15 |
| **grid-aligned line 2-D** | 64 | 0 | **0.00e+00** | **0.00e+00** |
| **axis-aligned box 2-D** | 108 | 0 | **0.00e+00** | **0.00e+00** |
| rotated box 2-D | 104 | 0 | 2.2e-15 | 2.2e-15 |
| circle 2-D | 64 | 0 | 2.3e-12 | 4.5e-12 |
| grid-aligned circle 2-D | 60 | 0 | 6.5e-13 | 1.2e-12 |

Read it this way. **Planar and off-grid** agrees to rounding — if a flat facet ever differs, that is
a bug. **Planar and grid-aligned** agrees exactly. **Curved** differs by the third-to-second-order
aperture trade, which is the one deliberate accuracy cost of the whole design. **Edges and corners**
differ most, because a folded patch is where a chord approximation is weakest. **All of 2-D** agrees
to rounding, because in two dimensions the interface inside a cut cell genuinely *is* the chord
between its crossings, so there is nothing being approximated; the 1e-12 residuals there are the
width crossings are held off the edge endpoints by.

Speed: sphere at 64³, **0.043 s** against **0.068 s** for `chombo-discharge`. Root finding costs the
same either way; the difference is `Moments.cpp` and the least-squares machinery, which this never
calls.

Offline, against the moments Chombo stores, over 48,798 cut cells: zero cells declined, bodies
closing to 4e-16, and `Σ(alpha_hi − alpha_lo) = a_B n` reproduced to 4.4e-16 in the area and 6.7e-10
in the normal on every geometry.

### Solvers

`Electrostatics/MechShaft` with plain multigrid, iterations for `chombo-discharge` / `polyhedral`:
**10 / 11** in two dimensions, **15 / 15** in three. `Electrostatics/ProfiledSurface` with BiCGStab:
4 / 4 and 5 / 5. Bare multigrid does not converge on that geometry for *either* generator, which is
why the shipped default is a solver chain rather than `gmg` alone.

EBIS coarsening, three levels deep: the divergence identity holds to 3.3e-16 on every level, the
volume fraction conserves to 1.8e-15 over 5,597 coarse cells, and there are no pathologies of any
kind -- no kappa or alpha out of range, no centroid outside its own cell, no isolated or ghost cells.
The numbers match the existing generator.


## 4. The rules that were expensive to find

Full list in `PLAN.md` §5. The four that will bite a reimplementation:

**One fluid predicate, and it is `copysign`.** `-0.0 < 0.0` is false and signed-distance functions
do return negative zero. Corner classification, crossing detection and the saddle test must all use
`copysign(1.0, v) < 0.0`, because that is what Chombo uses. Negative zero produced four separate
bugs in this work.

**The ambiguous four-crossing face.** The bilinear saddle names the diagonal that meets through the
middle; the *other* diagonal is the pair the two chords cut off, whether that is the fluid pair or
the solid one. And the face polygon must be derived **from the pairing**, not re-read off the corner
signs — doing the latter yields a geometrically correct face that silently disagrees with the loop
assembly, which does trust the pairing. Two wrongs that hide each other: apertures look right, the
body will not close, residual O(1).

**Small-but-not-zero deciding something qualitative.** Apertures of 1e-28 dictating a weighted-mean
centroid; slivers of 1e-23 reading as a second connected component. Guard every magnitude comparison,
and derive the threshold from the quantity's dimension — the same sliver is an *area* of 1e-24 in
3-D and a *length* of 1e-12 in 2-D. The mirror error is just as real: a cell is covered when it has
*no* volume, not when its volume is below a tolerance, and a legitimate corner sliver of leg 1e-4 has
kappa 1.7e-13 with apertures of 5e-9. Test topology with exact zeros.

**Recognise a degeneracy by where it came from, never by how big it is.** This is the single most
productive rule in the project, and it was learned five separate times.

Crossings are held off edge endpoints by `s_edgeTolerance = 1e-12` so the combinatorics stay generic.
That displacement buys genericity and pays in artifacts, and the artifacts cannot be caught by any
magnitude threshold, because **their size does not scale the way a threshold assumes**. A displaced
corner sliver is a slab, not a corner: its volume fraction goes as the displacement itself rather
than as its cube, so 1e-12 sails through a 1e-15 cut. Every attempt to tune a tolerance failed;
every fix that asked *what produced this* worked.

The five instances, in the order they appeared:

1. **A facet in a node plane** leaves one side present only as exact zeros. The sign test calls the
   cell cut; it is not. Rule: a side represented only by exact zeros encloses nothing. **This rule
   belongs in `classify`, not in the body assembly** -- the generator classifies every cell but only
   assembles the ones classification calls cut. Putting it in the wrong place made those cells be
   collected as irregular, assembled as degenerate, and then dropped for having no volume.
2. **A crossing recorded exactly at an edge endpoint** sits on a corner the interface passes through.
   Displacing it opens a face aperture of 1e-12 where the truth is zero. Rule: leave those alone.
3. **A face with no area open to flux carries no arc.** Both cells sharing a face compute the
   aperture from the same crossings, so the two sides cannot disagree about whether the face is
   there. The one exception is a regular neighbour, which is full and whose faces are open by
   definition; withholding the arc there makes EBGraph stop with "former regular vof not connected
   to anything".
4. **A segment's cell face follows from the circuit that produced it, not from where its endpoints
   sit.** In two dimensions the fluid region is one polygon whose segments lie in different faces. A
   positional test with a 1e-12 plane tolerance attributes *every* segment of a 1e-12-wide sliver to
   the same face, the interface chord included; they cancel, every aperture reads zero, the faces
   are withheld, and coarsening reports multi-valued cells. That was 114 multi-valued VoFs against
   Chombo's 6.
5. **A corner is on the interface when its crossings are pinned to it**, not when its value is small.
   A machined surface puts corners at 4e-17 rather than at exactly zero, so instance 1 misses them
   and the cell is built as a wedge of width 1e-12 with kappa around 1e-19. Ask instead whether every
   edge leaving the corner towards the other side turns over within the displacement distance. This
   subsumes instance 1, since a corner at zero pins its own crossings.

If a sixth appears, the question to ask is not "what tolerance separates these" but "what step
created this, and can I test for that step".


## 5. The Chombo bug, for the record

`GeometryShop::getFullNodeWithCoveredFace` correctly sets the boundary centroid to the covered face's
own centre. Its caller `fixRegularCellsNextToCovered` then lost it for every cell an earlier covered
cell had already claimed by edge or corner, because that cell was marked irregular first and the
later face-loop skipped it. On an axis-aligned cube, 4,050 of 4,056 cells wrong — with exactly six
correct, one per face of the cube, being the cells whose covered partner comes first in
`BoxIterator` order. Confirmed by Chombo's own divergence identity: residual exactly 0.5.

Fixed in Chombo-3.3#30, merged. The submodule pin moves to it. It was live for the lifetime of
Chombo 3.3 and survived because it only fires when a covered cell has a *regular* face neighbour —
that is, when the interface lies exactly on a grid plane, which is precisely the case nobody writes a
convergence test for.

## 6. Decisions already taken, so they are not reopened

| | |
| --- | --- |
| `a_B` | store both `\|∫dA\|` and the true area; only the latter partitions a boundary centroid. The one item still needing a Chombo change, as its own PR. `PLAN.md` §7.1 |
| single-valued | permanent, and by construction. Fluid in several pieces is one VoF holding their sum, which is what `GeometryShop` already does |
| what persists | generators, retained at generation, on our side. **Not** recomputed: `BrentRootFinder` runs up to 100 implicit-function evaluations per crossing and `NewIntersectionIF` walks every electrode per call, so recomputation is 10²–10³ evaluations per cell against a 2.7 µs body build |
| new types | `CutCellSurface` and `CutCellBody`, private to `Source/Geometry/`. No new container: `BaseIVFAB` is already the sparse per-irregular-cell store |
| declined cells | fatal. `ComputationalGeometry::s_strictGeometry`, a compile-time constant since the run-time switch was never the right thing to hand a user |

## 7. What to do next, in order

Items 1, 3 and 4 as originally written are done or superseded. What stands now:

1. **A convergence study.** Unchanged and still the gap. The solver comparison in §3 says the cost
   is near zero at the resolutions tested, but not what happens as the grid refines. The apertures
   are second order rather than third, so a study on a geometry with an analytic solution is what
   turns "no measured penalty" into a statement about the scheme.

2. **Milestone 2** — refine by cutting the body. `PLAN_MILESTONES.md`, and #732. `CutCellBody::refine`
   exists and is verified over 38,124 cut cells and all 304,992 children; what #732 still owes is
   below.

3. **Complete a partly carried index space.** #731 lets the index space be carried only where the
   geometry needs resolving, which leaves it short of a contract the rest of Chombo relies on:
   `GraphNode::refine` is documented as *"the result is only defined if this EBGraph was defined by
   coarsening"*, and in an ordinary index space that holds everywhere, because every level below the
   finest is the coarsening of the one above it over the whole domain. `PhaseRealm::defineEBLevelGrid`
   promises that capability level-wide by calling `setMaxRefinementRatio` on every level below the
   finest, and it is redeemed in `EBISLayoutImplem::setMaxRefinementRatio`, which refines the whole
   coarse layout and calls `fillEBISLayout`. Under partial coverage that asks for finer-level graph
   and data where none was carried, and `EBCoarseFineParticleMesh::defineStencils` is the first
   client to notice. `CutCellBody::refine` is what fills it, which is why this belongs to #732 and
   not to #731.

4. **Regrids.** #731 wires the AMR grids to the boxes the geometry was carried over, on both setup
   paths, through `Driver::regridAmrOntoGeometry`. Run-time regrids still cluster from tags and can
   ask for refinement outside what the index space carries. #732's job.

### The plan for refinement, as it stands

Settled, so they are not reopened:

| | |
| --- | --- |
| how a child is built | always a cut from its parent. No root finding below the level the surface was taken at, ever |
| what persists | the surface, 160 bytes a cut cell. Not the body, which is 4,544 and cheap to rebuild from the surface |
| what is transient | bodies. Materialised per parent for the length of a build or a regrid, then dropped |
| the moments | recomputed per child, inside `refine`, by integrating the child's own polygons. Not a separate step, not interpolated, and no root finding |
| a regrid | copies moments *and* surfaces where old and new grids overlap, and cuts from the parent only for cut cells that are genuinely new |

The storage decision turns on one asymmetry: reconstruction is root finding, of order a hundred
implicit-function evaluations a crossing, while building a body from a finished surface is one
arithmetic pass at 2.7 microseconds. Keeping the surface avoids the expensive half. Keeping the
body as well would save a further quarter of the per-cell time for twenty-eight times the memory,
which at six hundred thousand cut cells a level is 96 MB against 2.7 GB.

**Phase 1, the driver and the on-demand fill.** A `LevelData<BaseIVFAB<CutCellSurface>>` a level,
beside `EBData`, moving under the same `copyTo`. A driver that cuts a parent once and keeps the
whole child set, rather than cutting a path per cell. A hook on `GeometryService`, so only the
polyhedral generator implements it, letting `fillEBISLayout` serve boxes outside a level's own
grids. `m_finerNodes` on every irregular coarse cell, since `PhaseRealm::defineEBLevelGrid` promises
`setMaxRefinementRatio` for every level below the finest. And a look at whether `reconcileSeam`
becomes vacuous once the fine cells under a generated coarse cell are exact cuts of it.

**Phase 2, the regrid.** An entry in the EBIS fill path that takes the old level and the new grids
and copies what still applies; derivation only for new cut cells; and the coverage-only restriction
in `Driver::regridAmrOntoGeometry` lifted, which is what makes run-time regrids and the cell tagger
legal.

**Phase 3, the depth sweep.** Time and memory at depths one to ten. Not before Phase 1: measuring
the present path would only measure the waste described below.

### Why the level is extended, rather than the layout patched

`EBISLevel::fillEBISLayout` is a thin wrapper around `EBISLayoutImplem::define`, which builds a
local `LevelData<EBGraph>` on the requested grids and copies from the level's own. A box with no
source keeps the factory default, silently. So there is no missing-box branch to hook, and the
choice is between patching the layout after the fact and extending the level so the copy finds
something.

The level is extended, and the reason is ownership rather than taste. Cutting a cell needs its
parent's surface, and the surfaces are kept on the layouts the generator made for itself, which are
not the layouts a realm asks about: the two are load balanced separately, so the parent surface a
rank needs to fill a box will often be on another rank. Extending the level does the cutting on the
generator's own layout, where the surfaces are local, and leaves the distribution to the copy that
`EBISLayoutImplem::define` already performs.

That carries one requirement into the extension itself: a box added to a level has to be given to
the rank that owns the parent region it will be cut from, or the problem simply moves up a level.
The load balance of the extension is ours to choose, so this costs nothing but has to be chosen
deliberately.

The level only ever grows. A regrid that stops asking for a region leaves what was built for it in
place, since the same region tends to be asked for again a few steps later, and a cell that is
still there is a cell nobody has to cut twice.

### Where it will hurt

- **The cut, not the byte count, is the cost at depth.** See the note below on `buildRefinedBody`.
- **`PhaseRealm::defineEBLevelGrid` is already 11.5%** of a geometry-only run -- 4.4 s of 38.4 s on
  MechanicalShaft at 128^3 -- and it is exactly what the fill hook extends. Work added there lands
  on a path that is already the second largest single item.
- **A cut holds about 72 kB of bodies live.** `refine` keeps `buffer[2][4]` and its caller keeps
  `children[8]`, at 4,544 bytes each. Depth is iterative so it does not accumulate, but these loops
  are `#pragma omp parallel for`, so it is per thread.
- **The surface store costs an exchange.** 96 MB a level is cheap; the extra `LevelData` is added to
  a path where `copyTo` inside `coarsenFrom` is already 7.7% of the run.
- **The surface store is load-bearing, not an optimisation.** With reconstruction off the table, a
  cut cell whose surface is lost cannot be refined at all, so it has to survive a regrid and a load
  balance. It does not have to survive a checkpoint -- `setupForRestart` rebuilds the geometry --
  which makes restart pay a reconstruction, as a start-up cost.
- **The timings above are borrowed.** 2.7 microseconds to build and 0.9 an output cell to refine
  come from the spike: single threaded, and from before the polygon store was sized down. Nothing
  here rests on them except the storage decision, which has a twenty-eight-fold margin.

### Two things refinement must not be built without

Both found while sizing the body store, both cheap to design in and expensive to retrofit.

**`buildRefinedBody` throws away seven eighths of its work, at every level, for every cell.** It
walks from the ancestor down to the wanted cell, and at each step refines the parent into all eight
children, keeps one, and copies it over the parent. It then repeats the entire chain -- including a
fresh surface reconstruction, with the root finding that implies -- for the *next* fine cell in the
same parent. At refinement 2 `PLAN.md` already records the coarse surface being reconstructed nine
times per coarse cell. At ten levels it is of order eighty body refinements per fine cell, plus one
reconstruction, with the discarded children being exactly the cells that will be asked for next.

Whatever fills a refined layout should cut a parent once and keep the whole child set while the
level is being built, rather than cutting a path per cell. This is a question about the shape of
the refinement driver, not an optimisation to apply afterwards.

**A regrid must move what it already has.** When the grids change, the cells that survive keep their
geometry: their moments, and whatever representation refinement needs, should be copied rather than
rebuilt, and only genuinely new cut cells should be derived. That means the EBIS fill path needs a
notion of regridding rather than only of building -- something that can be handed the old level and
take from it what still applies. Rebuilding every cut cell on every regrid is the difference between
geometry being a start-up cost and being a per-regrid cost, and at plasma-simulation depths the
second is not affordable.

It also constrains the representation: anything moved between grids is moved between *ranks*, so a
descendant cannot hold a raw pointer into a parent's storage. Whatever a refined cell keeps has to
survive being serialised and load balanced.

**The largest unknown is no longer structural.** Whether the existing multigrid, flux registers and
redistribution stencils accept a cut cell whose neighbour is a *coarser* cut cell is still open, but
it is now a question about accuracy rather than about whether the thing runs at all: over a partly
carried index space the graphs join up, coarsening is conservative to round-off, and 3-D RodSphere
completes a Poisson solve. What has not been done is measuring whether the answer is right.

## 8. One investigation worth not repeating

Plain multigrid took 1.5 to 1.9 times as many iterations with `polyhedral` on MechShaft. The obvious
reading -- that dropping the apertures from third order to second had made the coarse-grid operators
worse -- was wrong, and the measurement that showed it was cheap: **with AMR switched off the two
generators needed the same number of iterations**, 8 and 8 in three dimensions. A genuine order
penalty appears on a single level. This one did not, it did not grow with depth, and it never
appeared on ProfiledSurface at any bottom-drop setting.

The cause was instance 5 in section 4: 224 cells built as 1e-19 wedges and then removed again by the
volume threshold, each leaving a hole where Chombo has a cell, with the holes falling along the
coarse-fine interface. Fixing the classification closed the gap to zero in three dimensions.

Keep the general lesson. Before attributing a solver symptom to the accuracy of the discretisation,
**turn off AMR**. If the symptom survives it is the operator; if it does not it is the graph, and the
graph is where the artifacts of section 4 live.

## 9. Known limits

- The four rules in §4 and the evidence in §3 are Milestone 1's. Milestone 3 added its own, and they
  are in #731 and Chombo-3.3#31 rather than here: that a face on the edge of what was coarsened has
  one area and the generated side takes its interface from its apertures; that the ring which cannot
  be coarsened has to be decided a cell at a time or the geometry carried grows by a factor tending
  to three; and that eroding a coverage set globally rather than a box at a time cost eleven times
  the rest of the build.

- `Prototypes/CutCellRefinement/spike/refine_spike.cpp` still lacks the face-bend and apex
  representation. Both default off in Milestone 1, so its feature set matches, but Milestone 2 on
  curved geometry will want them.
- A cell with two or more covered faces — a concave corner — cannot record both face centres in one
  `RealVect`. Chombo-3.3#30 improves such a cell from "cell centre" to "one of the two faces", not to
  correct.
- A cell edge cut twice, where a solid wedge is thinner than a cell, is a resolution limit of any
  edge-crossing representation, marching cubes included. No instance survives in the present suite.
- `faceCentroid`'s normal-direction component is write-only in Chombo: `FaceData::m_faceCentroid` is
  a `RealVect` but a face is (SpaceDim−1)-dimensional, every consumer skips that slot, and
  `EBISLevel::sanityCheck` only bounds-checks it against `0.5 + tolerance`, which a stale ±0.5 passes.
  Do not start reading it.
