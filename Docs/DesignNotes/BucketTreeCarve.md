# Bucket-tree carve superparticle merge — design notes (issue #665)

Working notes for a new particle-merge algorithm, `bucket_tree_carve`, requested in
[chombo-discharge/chombo-discharge#665](https://github.com/chombo-discharge/chombo-discharge/issues/665).
This document is the hardened design produced by an extended design review before any
implementation was written, kept in the repository so the reasoning behind the design travels with
the branch rather than living only in chat history. The "Goal and scope" through "Naming and
integration" sections below describe that pre-implementation design as it was locked before any code
was written; the implementation has since landed in this same PR, and the "Bugs found during initial
implementation and testing" section further down records what changed during that process.

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
- **Stop rule**: at each node, compute `fractionalSpan(node) = volume(bbox(node)) / volume(cell)`
  from that node's own (data-driven) bounding box — an exact real-valued ratio, not an integer count,
  and not a new computation: it reuses the same AABB-volume quantity the carve key already needs.
  This is deliberately *not* used as a scaled budget (`fractionalSpan * ppc`) gating when to stop —
  an earlier version of this design tried that and it's buggy: a multi-cell node containing one
  genuinely hot cell among several near-empty ones can have `count ≤ fractionalSpan * ppc` even
  though that one hot cell alone badly needs reducing, so a scaled-budget stop check can silently
  leave a crowded cell untouched (a reintroduction, one level up, of exactly the per-cell-not-per-patch
  averaging failure this design exists to avoid). The actual rule:
  - While `fractionalSpan(node) > 1` (still spanning more than one cell): keep splitting (longest
    axis, unsnapped) as long as `count(node) > ppc` — flat `ppc`, never scaled by `fractionalSpan`.
    Once `count(node) ≤ ppc` while `fractionalSpan` is still `> 1`: stop and mark the node
    **unmergeable** — leave those particles untouched (a genuinely sparse, spread-out region, not a
    hidden hot cell, since a real hot cell would still have `count > ppc` and keep splitting).
  - Once `fractionalSpan(node) ≤ 1` (the node's own footprint no longer exceeds one cell's volume):
    switch to `equal_weight_kd`'s existing single-cell logic (count-balanced split, no
    weight-splitting) targeting exactly `ppc` leaves for this node's subtree, until reaching that
    target or running out of particles (a naturally sparse sub-region just bottoms out at singleton
    leaves, which is fine — see the merge threshold below).
  Because `fractionalSpan` is an exact ratio with no rounding step, there's no "barely-straddling
  sliver" imprecision to worry about the way an integer `ceil(extent/dx)`-based cell count would
  have — a leaf at `1.02×` a cell's volume is genuinely, exactly just over threshold, not an artifact
  of rounding.
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
- **Dynamics — three fixed communication phases, no iteration, no drain loop**, unlike `nn_pair`'s
  drain-based protocol. An earlier version of this design used two phases and deleted a member the
  moment its owner learned who'd nominally won it; that is **unsafe** and was corrected (see the
  worked example below) — a box's win doesn't mean it will actually commit, since it can
  independently lose *other* members to *other*, unrelated competitors and end up under the merge
  threshold itself. Deletion must wait for that box's own outcome to be known, which needs a third,
  still-bounded round:
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
     and take the argmin. This resolves the *nominal* winner of every locally-owned member instantly,
     no round-trip — but the owner does **not** delete anything yet (see above).
  5. **Phase 2 — verdicts to proposers, foreign members only.** Every rank replies, to whoever claimed
     one of its particles, with the nominal winner. A box already knows the nominal fate of its own
     local members for free from step 4; the round-trip is only needed for the foreign (ghost-origin)
     members it listed.
  6. *(local, per rank as proposer)* Once all verdicts for a box's foreign members arrive, assemble
     its *provisional* final membership (self-owned from step 4 + foreign-won from step 5). Apply
     `PosValid` (reject → the box does not commit, same handling as an under-threshold chop below).
     If the surviving count is **≥2** and `PosValid` passes: the box commits — gather (all data
     already available locally, either truly local or already held as a ghost copy from the original
     ghost-fill, so no new data transfer is needed here), combine, insert the result via the
     placement rule below. If under 2 survive, or `PosValid` rejects: the box does not commit.
  7. **Phase 3 — commit/release, only for a box's foreign members.** Every proposer tells each foreign
     member's owner, from step 6's outcome, either "commit" (this box actually merged you — delete
     your copy now) or "release" (this box did not commit — you're untouched, exactly as if you'd
     never been claimed). Only now does an owner ever delete a particle. A locally-owned member whose
     nominal winner (step 4) was the owner's own box needs no message for this — the owner already
     knows its own box's step-6 outcome directly.
- **Worked example showing why step 7 is necessary** (three ranks at a junction; boxes' keys satisfy
  `Kx < Ky < Kz`): `Bx{p1(X-owned), p2(Y-owned)}`, `By{p2(Y-owned), p3(Z-owned)}`,
  `Bz{p3(Z-owned), p4(X-owned)}`. Argmin: `p1→Bx`, `p2→Bx` (beats `By`), `p3→By` (beats `Bz`),
  `p4→Bz`. Assembling actual memberships: `Bx` ends up with `{p1,p2}` (2 → commits); `By` loses `p2`
  to `Bx`, leaving only `{p3}` (1 → does **not** commit); `Bz` loses `p3` to `By`, leaving only `{p4}`
  (1 → does **not** commit). `p3`'s owner (Z) learns via Phase 2 that `By` nominally won `p3` — if Z
  deleted `p3` right then, `p3` would vanish with no merged particle anywhere holding its weight,
  since `By` never actually commits. Phase 3 is what tells Z "release, not commit" for `p3`, so Z
  correctly keeps it as an ordinary unmerged particle instead. This generalizes to chains of any
  length, not just three ranks — any cascade of pairwise box comparisons can produce a box that wins
  something but still ends up under threshold once its other members are independently contested
  elsewhere.
- **Options considered for closing this and why Phase 3 was chosen**: (a) collapse to winner-take-all
  per connected cluster of overlapping boxes, so a commit is unconditional once decided (like
  `nn_pair`'s pairwise accept) — rejected because resolving a cluster correctly is the same
  cascading problem again (a cluster "loser" can itself have been the effective winner over a third
  box), so it doesn't actually avoid a third round, and it merges strictly less (anything not held by
  the single cluster winner goes untouched, even members no one else contests); (b) accept the data
  loss as a documented, rare limitation — rejected outright, this is a real conservation violation in
  a codebase that explicitly tests for it (`BrownianWalker.verify_conservation=true`), not a
  cosmetic issue.
- **Why this is still bounded and non-iterative despite the third phase**: ghost width is hardcoded
  to 1 and every box is ≤1 cell (a direct consequence of the tree's stop rule), so only *direct*
  neighbors can ever hold a competing claim on a given particle — a rank 2+ hops away has no
  visibility into it at all. Fan-in/fan-out per particle is small and fixed, same shape as the
  existing ghost exchange, for all three phases. Three is still a FIXED number of rounds, not a
  drain/convergence loop — no re-proposal, nothing repeats. The added round-trip's actual cost is
  worth measuring once there's a running prototype rather than assumed.

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
  for a carve box it's only the locally-owned subset. Foreign members are never touched by `Scatter`
  at all (it can't reach another rank's container) — their deletion happens on their own owning
  rank's side, outside the `Scatter` callback entirely, triggered by that owner receiving "commit" in
  Phase 3 of the dynamics above. `Scatter` never needs to know which tier it's finishing, or anything
  about the carve protocol's phases.
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

## Bugs found during initial implementation and testing

The design above was locked before any code was written. Four genuine correctness bugs were still
found via actual multi-rank test execution (`BrownianWalker.verify_conservation=true` tripping),
not via design review — worth recording since all four trace back to the same root cause: an
insufficiently precise notion of "identity" once a rank can hold more than one competing box.

1. **Self-claims must accumulate, not overwrite.** The same particle can be a member of more than
   one of a rank's *own* boxes at once — gathered once as a local member of its home patch's build,
   and again as a ghost in a neighboring patch's independent build on that same rank. An
   `unordered_map<ParticleID, optional<BucketBoxKey>>` that assigns (rather than appends) silently
   drops the earlier claim, so only one of the two competing same-rank boxes was ever argmin'd —
   double-counting the particle's weight if both went on to commit. Fixed by making self-claims a
   `vector<pair<BucketBoxKey,int>>` per particle.
2. **"Owned by this rank" is not "physically resident in this patch."** Exposure and
   interior-vs-boundary classification must be evaluated over a leaf's *locally-resident* (non-ghost)
   members only, never over "owner == myRank" — the same particle is a ghost in every patch except
   its home one, even on the rank that owns it. Testing a ghost's cell against a *different* patch's
   own `ParticleGhostMask` is also undefined behavior (violates `numTargets()`'s precondition that
   the cell lie within that mask's own box). Fixed by adding an explicit `BucketParticle::isGhost`
   field, set at gather time, and switching every relevant check to it.
3. **A verdict needs the winning claim's own box index, not just its rank.** Two of the same rank's
   own boxes (different patches) can both claim the same foreign particle; a `{memberID, winnerRank}`
   verdict cannot tell them apart, so both claimants were told "your rank won" and both counted it.
   Fixed by adding `proposerBoxIdx` to a claim and redesigning the verdict to
   `{memberID, claimantBoxIdx, won}`, sent one per incoming claim rather than one per rank.
4. **A verdict must compare winning *identity*, not key *value*.** Comparing key values
   (`!(bestKey < claim.key)`) to decide `won` lets an exact tie between two genuinely independent
   claims — plausible here, since `BrownianWalker`'s regular initial particle distribution can give
   two independently-built boxes the exact same AABB volume — tell both tied claimants "you won".
   Fixed by comparing the winning claim's identity `(rank, boxIdx)` instead, which required Step 3 to
   record the actual winning foreign claim's own box index rather than a `-1` sentinel (needed so the
   identity comparison can ever match a foreign winner at all).

All four were caught by `verify_conservation`'s weight-drift check, confirmed root-caused via
temporary per-rank weight-diff instrumentation (removed once each fix was confirmed), and reproduced
only at multi-rank counts (bugs 1–2 also reproduce single-rank, since a single rank already owns
multiple patches; bugs 3–4 need ≥2 ranks). After all four fixes, single-rank and 4/8-rank 2D runs
pass cleanly with no asserts, NaNs, or conservation errors.

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
