# Bucket-tree carve superparticle merge — design notes (issue #665)

Working notes for a new particle-merge algorithm, `bucket_tree_carve`, requested in
[chombo-discharge/chombo-discharge#665](https://github.com/chombo-discharge/chombo-discharge/issues/665).
This document is the hardened design produced by an extended design review before any
implementation was written, kept in the repository so the reasoning behind the design travels with
the branch rather than living only in chat history. It is a starting point for implementation, not
a finished feature — nothing described here has been implemented yet.

## Goal and scope

Build a whole-patch spatial partition ("bucket tree") as an alternative to `nn_pair`'s k-NN graph,
with a "z-buffer carve" step resolving the boundary/skin region where neighboring patches disagree
about how to group particles. The interior of a patch should merge locally with zero communication;
only the region actually touching a patch/level interface should need cross-rank coordination.

**In scope for this PR**: `bucket_tree_carve` (enum `ParticleMergeMethod::BucketZCarve`, ParmParse
string `bucket_tree_carve`).

**Out of scope, deferred**: a sibling algorithm `bucket_tree_nn` (`bucket_tree_carve`'s planned
nearest-neighbor "skin" variant, reusing the same interior tier but resolving the boundary via
`nn_pair`'s existing candidate-search machinery instead of the carve protocol below). Both variants
were discussed together because they share the interior tier, but only the carve variant is being
built now, so they can be compared empirically later.

## Why a new tree, not an extension of `PointCloudBVH`

The build was initially expected to extend EBGeometry's `PointCloudBVH` (used by `nn_pair`) with a
custom weight-aware leaf predicate. That direction was dropped once it became clear that `ppc`
(particles-per-cell) is a **target leaf count**, not a leaf-size cap — confirmed directly from
`equal_weight_kd`'s own usage, where `a_ppc` is passed straight through as `buildEqualWeightKDLeaves`'s
`a_maxLeaves` (`Source/Particle/CD_ParticleManagementImplem.H:596`). Achieving a *global* leaf-count
budget needs breadth-first, budget-aware splitting — the shape `buildEqualWeightKDLeaves` already
has, not `PointCloudBVH::buildTree`'s depth-first, per-node local-cap recursion. So the build for
`bucket_tree_carve` is a **new sibling of `buildEqualWeightKDLeaves`**, living entirely in
`ParticleManagement`, with **no EBGeometry changes** required for this algorithm. (A per-cell build
was also considered and rejected — per-cell BVH/tree construction was tried for the unrelated
`nn_pair` prototype and lost to whole-patch batching because per-region build overhead dominates at
small N; the same lesson applies here, so this build is one whole-patch tree, not one per cell.)

## The tree build

- **Split rule**: longest (real-space) axis, data-driven, **deliberately unsnapped** to any grid —
  the whole point of this design is that particles must be able to merge across a cell face, so
  there is no grid-alignment step at all. K-ary branching via repeated binary halving, same shape as
  `PointCloudBVH::buildTree`/`buildEqualWeightKDLeaves`. No particle splitting during the build (no
  synthetic median-particle division, unlike `equal_weight_kd`). Weight and energy ride along as
  payload only and are never used to choose a split — the simplest thing that could work, to be
  measured before adding weight-awareness.
- **Stop rule**: at each node, recompute `spannedCells = product over axes of ceil(extent[d] / dx)`
  from that node's own (data-driven) bounding box. The node's local leaf budget is
  `spannedCells * ppc`, re-derived per node (not inherited as a fraction of the parent's budget) —
  this is what prevents a single, spatially-isolated hot cell in an otherwise near-empty patch from
  being allowed to fragment into far more leaves than `ppc`, which a flat `numCellsInPatch * ppc`
  ceiling alone would fail to prevent (a hot cell can trivially satisfy a ceiling sized for the whole
  patch while still needing its own, much smaller, local budget enforced). Once `spannedCells == 1`,
  behave exactly like `equal_weight_kd`'s existing single-cell logic (count-balanced split, no
  weight-splitting) until reaching `ppc` leaves or running out of particles. If a node's count
  already satisfies `ppc` but its AABB is still bigger than one cell (a sparse region), stop and mark
  it **unmergeable** — leave those particles untouched. `spannedCells`'s integer rounding (a
  barely-straddling sliver leaf counts as 2 cells) is a known, direction-safe over-estimate,
  deliberately left uncorrected until profiling says otherwise — it only ever matters when a node is
  both near a cell face *and* still holds enough particles to want splitting, and in that case a
  slightly looser budget is harmless.
- **Splitting under-full cells is not this algorithm's job.** `ItoSolver::splitAndRebuildFromMergeContainer`
  (`Source/ItoDiffusion/CD_ItoSolver.cpp`) already exists as a separate, algorithm-agnostic step,
  called unconditionally after every existing `makeSuperparticlesNnPair*` method, using the same
  `ppc` — a purely per-cell greedy-heaviest-split pass with no dependency on whatever tree/leaf shape
  the merge phase used. An "unmergeable" leaf, or a naturally-singleton leaf, is not a dead end for
  `bucket_tree_carve` to handle specially: it's simply left untouched, and this existing pass brings
  the affected cells back up to target afterward. `bucket_tree_carve`'s own
  `makeSuperparticlesBucketZCarve` should call it too, exactly like the `nn_pair` variants do. This
  is `nn_pair`'s decoupled merge/split convention, not `equal_weight_kd`'s in-tandem one, and the
  distinction matters: don't reintroduce in-tandem splitting into this algorithm's tree build.

## Leaf classification: interior vs. boundary

A leaf is **interior** (safe to merge immediately, zero communication) iff none of its *local*
(owned) members' cells are `isBoundaryExposed` — a per-cell test already used by `nn_pair`
(`isBoundaryExposed(const IntVect&, const ParticleGhostMask&)`,
`Source/Particle/CD_NearestNeighborParticleMergeImplem.H:266-269`), checked against the same-level
mask **and** the coarse-to-fine/fine-to-coarse masks together (see the existing three-mask check at
`CD_NearestNeighborParticleMergeImplem.H:1398-1406`), OR'd over every distinct cell a leaf's local
members touch (leaves are not grid-aligned and can span multiple cells). Otherwise it's a
**boundary** candidate, handed to the carve protocol below.

### An alternative rule was designed, confirmed safe, and rejected anyway — for artifact reasons, not correctness

A tempting alternative is: "a leaf is boundary iff it contains a ghost-origin particle," checked at
the leaf level rather than the cell level. It looked appealing because `isBoundaryExposed` is
conservative — since a leaf's AABB is smaller than a cell, many leaves inside an exposed cell
contain zero foreign ghost members and, on the surface, look safe to merge immediately;
`isBoundaryExposed` defers all of them to the carve pass anyway, at the cost of some extra (bounded,
single-round) communication for leaves that turn out to have no real competitor.

**This alternative rule is genuinely safe, not merely unsafe-if-you're-careless.** An initial
concern was that it allows a genuine double-merge: rank L's own leaf (ghost-free from L's own view)
gets committed immediately; rank R's independently-shaped tree build (built over a completely
different aggregate particle set) can legitimately group those same particles with R's own local
particles into a candidate box that lists them as members, since ghost-fill ships an entire cell's
contents to a neighbor regardless of how the *local* recursion happens to sub-partition it. But this
only breaks correctness if implemented naively (L's commit is unconditional and R's box has no
defined fallback on finding a listed member already gone). Requiring a box's gather step to fail
gracefully for an already-consumed member — drop it and proceed with whatever remains, via the same
chopping mechanic already used for ordinary argmin losers ("owner gets preference": a particle has
exactly one owner, so an owner's own unconditional local commit is an always-winning claim with no
ambiguity) — removes the double-processing entirely. Implementation, if this is revisited: **two
containers** — a leaf with no ghost members goes straight into a merged-particles container
(committed immediately, same as the interior tier here); a leaf with ≥1 ghost member puts only its
*local* members into a separate skin-particles container (tagged with its full local+ghost
membership, so the key and routing still work), and that container's boxes go through the same carve
protocol described below, with the one addition that a box's gather step must tolerate an
already-gone member instead of erroring.

**Why it was rejected anyway**: not safety, but the resulting merge strategy's effect on quality near
boundaries. `isBoundaryExposed` is simpler to build and verify — one build pass, uniform rule, every
exposed particle goes through the same coordinated argmin, so the grouping decision at a boundary is
always made by the argmin's global total order. The ghost-content + owner-gets-preference
alternative is faster (tighter skin, fewer leaves deferred) but makes the *actual* grouping decision
near a boundary partly an accident of which rank's independently-shaped tree happens to classify a
shared-adjacent particle as "safe" first, rather than a decision the argmin's total order would have
made uniformly. Every particle still ends up merged exactly once either way — but the alternative can
produce a different, less uniform outcome right at boundaries, and is harder to reason about under
the kind of multi-rank/multi-round stress testing this area of the code already gets. Treat this
alternative as a documented, ready-to-implement, **safe** optimization to revisit if profiling later
shows `isBoundaryExposed`'s conservatism costs something material — not a rejected-for-correctness
dead end.

## Interior tier

Commits immediately, zero communication. Guaranteed to never need cross-rank or cross-level
remapping: `isBoundaryExposed`'s same-level/C2F/F2C awareness already rules out any coarse/fine
interface proximity for an interior leaf, so its merge result is guaranteed to stay in the same
patch at the same level. The result is inserted directly into the proposer's own container, with no
destination lookup at all.

## Boundary/carve protocol

- A **box** is a boundary leaf with ≥1 local member (a leaf with zero local members isn't proposed by
  anyone — no rank has standing to claim particles it doesn't own).
- **Key**: `(AABB volume, anchor global-id)`, compared lexicographically ascending (smallest/tightest
  box wins). Anchor-id alone (lowest global particle id among members) is deterministic but
  arbitrary with respect to merge quality; volume alone is not a strict total order (every
  single-particle box has zero volume, so ties are common) — anchor-id is what makes the pair a
  strict total order and guarantees full determinism. The pair is computed once by the proposer and
  copied verbatim into every claim message, so comparing it for equality across boxes is safe
  bit-for-bit, no epsilon needed.
- **Eligibility is by membership, not geometry**: a particle may only join a box that explicitly
  listed it when built, never by raw AABB containment.
- **Dynamics — exactly two fixed communication phases, no iteration, no drain loop**, unlike
  `nn_pair`'s drain-based protocol:
  1. Every rank builds its tree, classifies leaves, commits interior merges immediately. Boundary
     leaves become boxes.
  2. Barrier.
  3. **Phase 1 — claims to owners.** Every rank ships one `(memberGlobalId, boxKey, proposerRank)`
     tuple per member per box to that member's owning rank, reusing the *same* rank-adjacency
     topology as ghost-fill itself (`ParticleGhostMask`'s target list) — the ranks that could ever
     claim one of my particles are exactly the ranks I already ship ghosts to. Same exchange shape as
     the existing `nnMergeExchangeByRank` helper.
  4. *(local, per rank as owner)* For each of its own particles, gather every claim on it — including
     its own box's claim if any, which costs no message since the rank already knows its own key —
     and take the argmin. This resolves the fate of every locally-owned member instantly, no
     round-trip. The owner can delete any particle whose winner isn't itself right away; nothing
     later changes the outcome.
  5. **Phase 2 — verdicts to proposers, foreign members only.** Every rank replies, to whoever claimed
     one of its particles, with the winner. A box already knows the fate of its own local members for
     free from step 4; the round-trip is only needed for the foreign (ghost-origin) members it
     listed.
  6. *(local, per rank as proposer)* Once all verdicts for a box's foreign members arrive, assemble
     final post-chop membership (self-owned from step 4 + foreign from step 5). Apply `PosValid`
     (reject → leave untouched, same handling as an under-threshold chop). If the surviving count is
     **≥2**, gather (all data already available locally — either truly local, or already held as a
     ghost copy from the original ghost-fill, so no new data transfer is needed here either),
     combine, insert the result via the placement rule below. If under 2 survive: leave whatever
     remains untouched.
- **Why this is bounded and non-iterative**: ghost width is hardcoded to 1 and every box is ≤1 cell
  (a direct consequence of the tree's stop rule), so only *direct* neighbors can ever hold a
  competing claim on a given particle — a rank 2+ hops away has no visibility into it at all. Fan-in/
  fan-out per particle is small and fixed, same shape as the existing ghost exchange. No third
  communication phase, no re-proposal/convergence loop.

## Template hooks (the "as general as `nn_pair`" requirement)

Five hooks, same style as `nn_pair`'s own customization surface:

- **Gather** — extracts a lightweight `{position, weight, energy, globalId}` payload from a real
  particle. Invoked once, up front, to build the whole-patch packed array that feeds the tree build;
  reused for the eventual combine, no second gather pass.
- **Combine** — **N-ary, not pairwise** (a deliberate difference from `nn_pair`'s pairwise `Combine`,
  since leaves/boxes here are groups by construction, not pairs off a k-NN graph). One call per
  leaf's/box's *final* membership range, producing one merged particle. A pairwise-fold alternative
  matching `nn_pair`'s signature exactly was considered and rejected: it would need a fixed
  deterministic fold order for floating-point reproducibility, and an N-ary
  weighted-centroid/weight-sum/energy-sum is barely more code than a pairwise one, so the
  cross-algorithm-reuse argument was judged weak.
- **Scatter** — writes the merged result back into the real container: deletes only this leaf's/box's
  *local* originals, reconstructs a full particle from the merged payload, inserts it. The same
  function serves both tiers unchanged — for an interior leaf "local originals" is the whole leaf;
  for a carve box it's only the locally-owned subset, since foreign members were already
  independently deleted by their own owning rank in step 4 of the dynamics above. `Scatter` never
  needs to know which tier it's finishing.
- **Allocator** — fresh, rank-namespaced global id for the new particle, reusing `nn_pair`'s existing
  convention as-is.
- **PosValid** — EB-aware check on the computed merge position, reusing `nn_pair`'s convention.
  Rejection handling: reject invalid merges, leave the group's members untouched, same as an
  under-threshold chop; never force an invalid merge through.

## Merge threshold and the split/merge decoupling

**≥2** surviving members to bother merging — both for a naturally-singleton interior leaf and for a
post-chop/post-`PosValid` carve box; 0 or 1 stays untouched. This is not a special carve-only rule:
it's the same check applied at two different moments (interior-commit time, immediately; carve-commit
time, once the two-phase protocol resolves a box's final membership). As covered above, a leftover
singleton or sparse "unmergeable" leaf is not a gap this algorithm needs to fill — that's exactly
what `splitAndRebuildFromMergeContainer` is for, run afterward, decoupled from this algorithm's own
tree shape.

## Result placement across the AMR hierarchy

A coarse level's own valid region is genuinely non-rectangular wherever a finer patch is nested in
it (Chombo's `DisjointBoxLayout` represents this as several boxes whose union excludes the fine
region, not one L-shaped box). So a carve box built from coarse-local + fine-ghost data (via the
C2F/F2C masks) can compute a merged centroid that lands inside the fine region rather than anywhere
in the coarse level's own boxes — a real possibility this design has to handle, unlike the interior
tier (see above, guaranteed same-patch-same-level by construction).

Handled as a cheap-first/expensive-fallback split, computing ownership only where genuine
disagreement can occur: try a trivial point-in-my-own-box test first (no lookup); only on failure
fall through to whatever general AMR-hierarchy-aware placement machinery (`ParticleContainer::findDestination`)
already handles ordinary particle drift across a rank/level boundary between timesteps — reuse, not
new code. The expensive path only ever runs for the minority of carve results that land outside the
proposer's own box.

## Naming and integration (implementation plan)

- ParmParse string `bucket_tree_carve`; enum `ParticleManagement::ParticleMergeMethod::BucketZCarve`
  (`Source/Particle/CD_ParticleManagement.H`, inserted after `NnPairHash`); string mapping in
  `mergeMethodFromString` (`Source/Particle/CD_ParticleManagementImplem.H`).
- New sibling header pair `Source/Particle/CD_BucketTreeCarveParticleMerge.H` /
  `...Implem.H`, included from the umbrella `CD_ParticleMergers.H`, structured the same way
  `CD_NearestNeighborParticleMerge.H` is (build / classify / commit-interior / propose-boxes /
  exchange-claims / resolve-argmin / exchange-verdicts / commit-carved, each its own function in
  `namespace detail`), with a top-level `ParticleManagement::mergeBucketTreeCarve(...)` entry point
  taking the five hooks above plus the three ghost masks. Single pass, **no round loop**.
- `ItoSolver::makeSuperparticlesBucketZCarve(WhichContainer, const Vector<int>&)`, modeled directly on
  the self-contained `makeSuperparticlesNnPairOneCell` (not the template-backend-split
  `makeSuperparticlesNnPairTree`/`Hash`, since there's only one build strategy here): extract →
  gather/combine(N-ary)/scatter/`isPositionValid` lambdas → rank-namespaced `allocateID` → fill
  same-level + C2F + F2C ghost masks at width 1 → call `mergeBucketTreeCarve` once → call the
  existing `splitAndRebuildFromMergeContainer`.
- Dispatch: one new `case` in `ItoSolver::makeSuperparticles(WhichContainer, Vector<int>, ParticleMergeMethod)`'s
  switch, between `NnPairHash` and `External`.
- **No new ParmParse tunables** beyond the standard `ppc` argument every merge method already takes —
  no iterate/fallback/max-cell-distance/max-rounds knobs, since there's no drain loop and no
  configurable search radius. `registerOperators()` needs one new (simpler than `nn_pair`'s)
  width-1 ghost-mask registration, conditioned on this method being active.
- New test scaffold `Exec/Tests/BrownianWalker/BucketTreeCarve/` (its own top-level name, not an
  `NNMerge*` sibling — this isn't part of the `nn_pair` family), mirroring `NNMergeOneCell/`'s
  `GNUmakefile`/`main.cpp` verbatim and its `.inputs` files with the `merge_algorithm` line changed;
  a matching block pair added to `Exec/Tests/BrownianWalker.ini`.
- Docs: a new bullet in `Docs/Sphinx/source/Solvers/Ito.rst`'s `merge_algorithm` list (and its
  AMR-level-exceptions sentence), and a new subsection in `Docs/Sphinx/source/Source/Particles.rst`
  describing this as its own distinct family (a spatial-partition/carve algorithm, not a
  nearest-neighbor pair algorithm) rather than folding it into the existing `nn_pair` bullets.
  Existing `literalinclude` ranges in `Ito.rst` will need re-verification after any line-count shift
  in `CD_ItoSolver.H`/`CD_ParticleManagement.H`, per this repo's own documented convention for that
  (`CLAUDE.md` §"Sphinx literalinclude line references").

## Verification plan

- Build 2D and 3D, `OPT=HIGH` and `DEBUG=TRUE`, with MPI.
- Run the new `BucketTreeCarve` BrownianWalker test at 1/4/8 ranks in 2D and 1/4 ranks in 3D, across
  regrids and multiple AMR levels, checking for no asserts/NaNs and that
  `BrownianWalker.verify_conservation=true` holds (total particle weight conserved across merges) —
  matching this area's existing stress-testing convention of exercising rank × round ×
  refinement-ratio combinations together, not each axis alone.
- Specifically exercise: particles clustered right at a patch face (forces boundary/carve activity,
  not just interior commits), and a coarse/fine interface with an L-shaped coarse region (exercises
  the placement fallback path).
- `pre-commit run --all-files`, plus a full Sphinx rebuild to catch literalinclude drift.
