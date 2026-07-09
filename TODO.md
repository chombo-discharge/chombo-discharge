# TODO — nearest-neighbor particle merge

Design and performance items that have been identified during development of
`Source/Particle/CD_NearestNeighborParticleMerge.H` / `CD_NearestNeighborParticleMergeImplem.H`
(and its ItoKMC integration), deliberately deferred rather than fixed immediately. Kept here so
neither of us has to re-derive them from scratch later.

## Performance

- ~~`findNearestNeighborCandidates()` is brute-force O(n^2) per patch.~~ **DONE.** Replaced with the
  hybrid spatial index design locked earlier in this file's history, implemented in five
  incremental, independently-verified phases (each landed as its own commit):

  1. **`MergeBlockGrid`** (new `Source/AmrMesh/CD_MergeBlockGrid.H`/`.cpp`) -- a static, purely
     geometric per-patch structure (patch box -> block-index-space bounds plus a 1-block halo),
     modeled exactly on `LevelTiles`'s lifecycle: built once per regrid in `Realm::regridBase()`
     (`Realm::defineMergeBlockGrid()`, right after `defineLevelTiles()`), aliased out via
     `Realm::getMergeBlockGrid()`/`AmrMesh::getMergeBlockGrid(realm)`.
  2. **`buildNNCellBuckets()`** (`CD_NearestNeighborParticleMerge.H`/`...Implem.H`) -- an O(n)
     counting sort (tally -> prefix-sum -> scatter) that buckets a patch's full local-valid and
     ghost sets (no cross-level carve-out) by cell into a flat `NNIndexedParticle` array
     (`NNCellBuckets`), rolling per-cell counts up into per-block counts (`blockCounts`) for free.
     Built once per patch per `mergeNearestNeighborsRound()` call (not per internal
     `a_iterateLocalTierToConvergence` pass -- those reuse it via lazy deletion against
     `a_consumedIDs` instead, since ghosts are always fresh every round and natives may have moved
     between rounds, so there was never a valid reason to cache across rounds anyway).
  3. **`findNearestNeighborCandidates()`** now walks that index via Moore-ring expansion from the
     query's own cell, pruned by point-to-AABB distance against whichever is tighter: the current
     worst-of-K distance, or a new `a_maxCellDistance` parameter (level-oblivious by construction,
     since a finer candidate only ever has a *tighter* exact radius, never looser -- the exact,
     level-aware check still runs downstream in `resolveTrivialTier()`/`judgeProposals()`,
     unchanged). Ring expansion is capped by two independent conditions: a distance-based one (pure
     speed) and a coverage-based one -- a new `NNCellBuckets::observedCellBox` field bounds every
     cell that actually has a particle, guaranteeing no candidate is ever missed even if the
     "expected" 1-block halo geometry is ever violated.
  4. **Adaptive per-cell midpoint-split sub-partition** (`NNSubNode`, `NNSpatialIndex::subPartitions`,
     built lazily via `buildNNSpatialIndex()`) -- only for cells whose occupancy exceeds a new,
     independent `a_subPartitionThresh` (a performance knob, deliberately separate from the
     physics-facing `a_numParticlesPerCellThresh`). Splits on the spatial midpoint of the longest
     axis (not a weighted median, unlike `buildEqualWeightKDLeaves()`), so a bursty cell (e.g. an
     avalanche) pays its own extra cost without inflating its neighbors'.
  5. `a_localValid`/`a_ghosts` in `findNearestNeighborCandidates()` now mean the patch's FULL set
     for the round, not a pre-filtered alive subset -- a deliberate, explicitly-documented contract
     change (the function does its own lazy-deletion filtering now, via a new `a_consumedIDs`
     parameter), which is what lets the index be built once per round instead of once per pass.

  **Verification**: acceptance bar during development was byte-for-byte identical merge edges
  (same query -> candidate/distance/fallback list) against the old brute force, not just "mass
  still conserves" -- confirmed on the smoke test at every phase, including with aggressive,
  artificially-forced sub-partitioning (threshold as low as 1) and combinations of
  `a_iterateLocalTierToConvergence`, `a_maxFallbackCandidates`, and `a_maxCellDistance`. Tested at
  1/4/8/16 MPI ranks, single- and multi-level, refRat 2 and 4 (the scenario that originally
  motivated this work). Mass conservation exact and centerline deviation ~0 throughout every
  phase. 2D/3D compiles clean, doxygen 0 warnings.

- ~~`std::map` used throughout instead of `std::unordered_map`.~~ **`liveCellCount` DONE -- and it
  turned out to be a correctness bug, not just a performance one** (see the below-threshold-merge
  entry in "Design / correctness coverage"). `particlesByID`, `exposed`, `outgoingTargetOf`, and
  `myParticlesByID` are still `std::map`, but those are keyed by `ParticleID` (an integral type with
  a genuine total order), so they only have the original, purely-performance motivation -- O(log n)
  lookup, poor cache locality, on a hot path looked up per particle/edge/fallback attempt. Should
  still become `std::unordered_map` for O(1) amortized lookups. `consumedIDs` and
  `hasOutgoingCommitment` (`std::set<ParticleID>`, same underlying issue -- ordered tree, O(log n)
  membership check) belong on this same list; not yet fixed. Confirmed via phase timing (see the
  `a_subPartitionThresh` entry below) that this doesn't matter for `a_iterateLocalTierToConvergence
  = false` (the common case so far -- `consumedIDs` stays empty for the only
  `findNearestNeighborCandidates()` call that round, so lookups are O(1) in practice regardless of
  container type), but would start costing real time once multi-pass convergence is actually used
  (`consumedIDs` grows across passes then).

- **`myParticlesByID` is a full copy of `particlesByID`, rebuilt every round.**
  (`mergeNearestNeighborsRound()`, the "5. judgeProposals()" section.) Built by iterating
  `particlesByID` and filtering out anything in `consumedIDs`, into a brand new map -- O(n log n)
  and a full memory duplication every call, purely to avoid checking `consumedIDs` at each lookup
  site inside `judgeProposals()`. Should be removed: pass `particlesByID` + `consumedIDs` through
  and check membership at the point of use instead.

- ~~`a_subPartitionThresh` defaulted to effectively "never sub-partition"
  (`std::numeric_limits<int>::max()`), and no caller had ever overridden it.~~ **DONE -- default
  changed to 32.** Triggered by a rough wall-clock comparison against `ItoSolver`'s existing
  `equal_weight_kd` merge algorithm (`makeSuperparticles()`), which came back roughly 10-18x faster
  for a single round/call on a realistic multi-level, variable-block-size, refRat-4 configuration
  (~97k particles/rank, threshold/target 16). Added temporary per-phase wall-clock timing to
  `mergeNearestNeighborsRound()` (raw `std::chrono`, deliberately NOT `CH_TIME`/`CH_TIMER` --
  Chombo's own `TraceTimer` profiler is reportedly too slow itself for a comparison like this) and
  found the local-tier pass (`findNearestNeighborCandidates()` + `resolveTrivialTier()`) accounts
  for ~90% of total time regardless. Root cause: since this function has had exactly one caller (a
  debug/test hook) so far, the adaptive midpoint-split sub-partition -- built specifically to avoid
  brute-forcing a crowded cell -- had literally never been exercised. Uniform particle sampling
  combined with unrefined coarse cells covering far more physical area than fine ones produces a
  handful of pathologically crowded cells that were each paying full brute-force cost. Repeated
  timing runs (3x each, to average out sampling noise from `drawBoxParticles()` not being seeded)
  confirmed enabling sub-partitioning at a reasonable threshold cuts total time by ~2.4x (~0.53s ->
  ~0.22s for the scenario above) with zero other code changes -- purely a matter of the config knob
  actually being used. Narrows the gap against `equal_weight_kd` from ~18x to ~7-8x for the
  single-round comparison, but does not close it: `findNearestNeighborCandidates()` runs one
  independent nearest-neighbor query per over-threshold particle (needed so the cross-rank
  propose/judge protocol always knows each particle's true nearest neighbor), whereas
  `equal_weight_kd` never computes "who is my neighbor" at all -- it recursively splits a cell's
  particles by weighted median in one shared O(k log k) pass to hit an exact target count. That is
  a structural difference between the two algorithms' approaches, not a bug, and further closing the
  remaining gap (if wanted) means rethinking that structural cost, not just tuning knobs. The
  temporary per-phase timing instrumentation (`NNPERF ...` pout() lines) and the equal_weight_kd
  comparison block (in `ItoKMCGodunovStepper::postInitialize()`'s debug hook, gated on
  `debug_test_nn_merge`) were left in place for further tuning work, at the user's request -- not
  yet removed.

- **Tried and reverted: replacing the per-cell midpoint-split sub-partition TREE
  (`NNSubNode`/`buildNNSubNodeRecursive()`/`nnMergeWalkSubNodes()`) with a flat, non-recursive
  regular sub-GRID (`gridRes^SpaceDim` sub-cells, direct O(1) arithmetic assignment instead of
  recursive `std::partition`).** Motivated by the item above's own finding that `localTier` still
  dominates even with sub-partitioning enabled, plus the observation that a midpoint split is
  purely geometric (never data-dependent), so a cell's sub-partition is mathematically identical to
  a uniform grid refinement -- no tree/recursion should be structurally necessary. Implemented,
  verified correct (exact mass conservation and `maxCenterlineDeviation ~ 0` throughout, at 1 and 8
  ranks), but consistently SLOWER than the tree, not faster -- `localTier` roughly 0.63-0.68s vs
  the tree's ~0.17s at `a_subPartitionThresh=32` on the same multi-level test config. Root cause,
  confirmed via temporary per-query instrumentation (ring/sub-cell-visit counters): the grid's
  formula splits every axis simultaneously (`gridRes = ceil((count/thresh)^(1/SpaceDim))`), so even
  a cell only barely over threshold gets ALL axes split at once (`gridRes=2` was the overwhelmingly
  common case at this test's actual population sizes, 33-256ish particles per crowded cell). The
  tree, by contrast, splits only the single longest axis per level, so a comparable "one split"
  treats "the other half" as ONE region with ONE exact AABB check; the grid instead breaks that
  same region into multiple separate sub-cells (3, for `gridRes=2`), each paying its own lookup
  cost even under a subsequently-added EXACT (not ring-count-heuristic) distance-to-boundary prune
  check -- recreating the recursion-style overhead the flat design was meant to eliminate, just
  redistributed into per-sub-cell lookups. Adding a minimum-worthwhile-`gridRes` gate (skip
  building a grid below `gridRes=4`) made it fail safe (falls back to plain brute force, ~0.49-
  0.53s) rather than actively regressing, but never actually fired at this test's population range,
  so it didn't deliver a win either. Fully reverted (tree restored byte-for-byte against the
  pre-attempt version, confirmed via diff) rather than left half-working behind a flag. The
  underlying idea (avoid data-dependent recursion since the split is geometric) may still be sound,
  but a naive "split every axis via one formula" flat grid does not reproduce the tree's actual
  advantage -- closing this gap for real would need a design closer to the tree's own adaptive,
  single-axis-at-a-time semantics (or a genuinely different approach), not just "make the tree
  flat". Left as an open idea, not a dead end -- revisit only with a concretely different design.

- **Planned: migrate the whole-patch spatial index onto `EBGeometry::BVH::PackedBVH`, once
  EBGeometry gains the constructors/leaf layout this needs.** Tracked upstream at
  [rmrsk/EBGeometry#92](https://github.com/rmrsk/EBGeometry/issues/92). `PackedBVH`'s SIMD node
  traversal (SoA child-AABB cache, dispatching to SSE/AVX/AVX-512 depending on `(T,K)`) is
  structurally the same pruned-descent walk `nnMergeWalkSubNodes()` already does by hand, just
  vectorized -- but its only build path today goes through `TreeBVH`, which is `shared_ptr`-per-node
  and `shared_ptr`-per-primitive, and its `signedDistance()` leaf evaluation pays a `std::sqrt` per
  candidate that a plain nearest-neighbor search never needs. Benchmarked standalone (not committed
  anywhere in this repo) against our own count-median build/leaf-scan to quantify exactly what's
  blocking a direct swap:
  - Building via `TreeBVH`+`pack()` was 4.5-17x slower than our own flat-array build across
    N=1,000-90,441 (worse at small N), almost entirely `shared_ptr<TreeBVH>`-per-node allocation --
    `BVCentroidPartitioner`'s count-split is the same algorithm we use, so it isn't the
    partitioning logic costing this. A prototype building straight into a `PackedBVH::Node`-shaped
    flat array (still individually `shared_ptr`-wrapping primitives) closed this to ~6-9% slower
    than our own build at realistic sizes (20K-90K particles) -- i.e. nearly the whole gap is the
    `TreeBVH` node allocation specifically, not primitive allocation.
  - Isolated at the leaf-scan level (leaf sizes 8-32, matching this tree's own typical leaf
    occupancy): `shared_ptr`-per-primitive scanning (dereference + `->signedDistance()`) was 3.4x
    slower than our current flat-array-with-indices leaf scan, and 6.4x slower than an SoA-packed
    leaf scan (`TriangleSoAT<T,W>`-shaped, position-only), at leaf size 24 -- split between scattered
    heap allocation (poor cache locality on every subsequent query, not just at build) and the
    per-candidate `sqrt` `signedDistance()`'s contract forces.
  - Net: once EBGeometry has (a) a direct-to-`PackedBVH` constructor (skip `TreeBVH` entirely), (b)
    `shared_ptr`-free primitive storage, and (c) either a squared-distance leaf contract or a
    templatized leaf-eval callback (`PackedBVH::signedDistance()`'s internal `evalLeaf` lambda is
    already factored out identically across every SIMD-width branch -- a natural seam to make it a
    template parameter instead of hardcoding `->signedDistance()`), replacing `NNSubNode`/
    `buildNNSubNodeRecursive()`/`nnMergeWalkSubNodes()` with `PackedBVH` should get SIMD node
    traversal essentially for free at current build cost, plus the leaf-scan SoA win on top. Full
    list of requested changes is on the linked issue. Not started -- blocked on that EBGeometry work
    landing (expected within days as of this writing) and updating
    `Submodules/EBGeometry` (or wherever it ends up living) to a version that has it.

  **Update, since the above was written** -- EBGeometry now has a working direct-to-`PackedBVH`
  constructor: build cost is at parity with our own build at realistic sizes (N=10K: EBGeometry
  actually faster, 795µs vs our 841µs; N=100K: 11.58ms vs our 11.43ms, ~1% slower). A real
  SAH-partitioned, SIMD-traversed (traversal + leaf eval both vectorized) prototype is measuring
  ~600-700ns/query, vs. our own production `localTier`/queries figure of ~1090ns -- roughly a 1.6-1.7x
  win, not the 3-5x initially speculated (see below for why that speculation didn't pan out).

  Two controlled experiments this branch ran to explain WHERE the win actually comes from, since
  the first two guesses were wrong:
  - **Recursion vs. iteration, isolated on our own K=2 tree**: a naive explicit-stack rewrite of
    `nnMergeWalkSubNodes()`'s exact algorithm was consistently ~18-25% SLOWER than the existing
    recursive version, not faster, at every N tested (1K-90K). Likely explanation: shallow, regular
    recursion (our tree depth is only ~12-24 levels) is well handled by the CPU's own hardware
    return-address prediction; a software array-based stack doesn't get that acceleration. So
    `PackedBVH::traverse()` being iterative is NOT where its advantage comes from -- whatever wins
    it delivers are from SIMD batching and cache-friendly layout, not iteration itself.
  - **Redundant AABB-distance recomputation** (the parent computes each child's distance to pick
    near/far order; the child then recomputes the identical distance at the top of its own visit,
    in both our code and this experiment) -- removing it gave a real, repeatable ~10-18% win,
    correctness-verified. Judged not worth applying to `nnMergeWalkSubNodes()` on its own: this
    function has already been through several bug-fix cycles this session, and if the `PackedBVH`
    migration lands in days and replaces it outright, spending verification effort on a
    soon-to-be-superseded double-digit-percent tweak isn't the right place to put it.

  **Point-cloud leaf design, agreed**: the SIMD-hot leaf primitive needs only **position + particle
  id** -- nothing else. Self-exclusion and the `consumedIDs` lazy-deletion check are both id-based
  (not index-based, unlike today's `NNIndexedParticle` scheme). Weight, payload, owner rank, AMR
  level, and even the local-vs-ghost flag are all looked up once per query AFTER a winner is
  selected, via the existing `particlesByID` map -- local-vs-ghost specifically via
  `exposed.count(id) > 0` (already built during the gather pass, populated only for locally-valid
  particles), not via anything carried in the tree at all.

  **Precision idea for the point cloud, not yet implemented**: since `CD_PARTICLE_REAL` is `double`
  throughout this codebase, storing `float` positions in the point cloud is a deliberate, scoped
  precision reduction -- but a low-risk one, since the actual merge computation (weighted centroid,
  mass) still reads the authoritative `double` record via id after the search, never the float copy.
  `float`'s ~7 decimal digits is far more than needed to correctly order candidates by distance at
  this codebase's typical domain scales. Goes further: store each particle's position NORMALIZED to
  the patch's own physical extent (patch size is already known), not the absolute domain position --
  `float` precision is relative (~7 significant digits regardless of magnitude), so normalizing to a
  patch-local, roughly-O(1) range maximizes the *absolute* precision available for
  particle-to-particle distance discrimination, regardless of where in the domain the patch sits.
  Implementation detail to get right if this is built: leaf AABBs derived from `double` positions
  truncated to `float` must round the low corner DOWN and the high corner UP (not naive
  round-to-nearest), or a genuine boundary particle could end up just outside its own leaf's
  `float`-rounded AABB and get wrongly pruned -- same safety-margin principle already used for
  `a_maxCellDistance`'s own `sqrt(SpaceDim)` bound.

  **Branching factor under consideration**: K=4 (AVX/plain, `double`) -> K=8 (AVX-512, `double`) ->
  K=16 (AVX-512, `float`) -- each step roughly matches SIMD width to element size, each halving (or
  better) tree depth for the same leaf count. Not verified: whether pruning effectiveness holds up
  as branching factor widens -- more children per level means more breadth to potentially visit,
  and while SIMD makes checking them cheap, whether `avgNodeVisitsPerQuery` actually drops
  proportionally to the depth reduction (rather than partially eating the win back via wider
  per-level breadth) hasn't been measured at K=8 or K=16 specifically.

- ~~Higher AMR refinement ratios inflate ghost candidate counts because of a ghost over-delivery
  bug in `Realm::defineParticleGhostMaskFineToCoar`.~~ **CORRECTED -- this was a misdiagnosis, not
  a bug.** Earlier investigation this branch's history found that fine-to-coarse ghosts ship ALL
  `refRat^D` fine children of a qualifying coarse halo cell, and framed this as "over-delivery" by
  comparing that reach against what a same-level, fine-level width-1 ghost would cover. That is the
  wrong baseline. The correct baseline is what the CONSUMING (coarse-side) query itself needs: its
  default search radius is 1 cell of *its own* level, i.e. 1 coarse cell -- which physically spans
  `refRat` fine cells in width. Shipping every fine particle inside that region is exactly the
  right amount, not excess; the coarse→fine direction does the mirror-image thing (delivers
  exactly 1 *fine* cell deep, matching a fine query's own 1-fine-cell radius). Both directions
  already deliver precisely "1 cell of the receiving level's own size" -- there is no missing
  acceptance-box filter to add, and adding one (as previously proposed) would be actively wrong: it
  would cut off fine particles a coarse query is legitimately entitled to see.

  What DOES still hold: at higher refinement ratios, a coarse cell's "1 cell" neighborhood
  genuinely contains more fine-level particles (`refRat^D` worth), so the candidate pool a
  coarse-side query must consider legitimately grows with `refRat` -- this is just AMR geometry,
  not a defect anywhere. This is exactly why the old O(n^2) brute-force search got so much worse at
  higher refinement ratios (confirmed directly: variable-sized-patch, 8-32 cells, `refRat=4`
  configs took much longer than an equivalent `refRat=2` config), and exactly why the spatial index
  (item above, now done) is the right fix: it makes searching a larger, legitimate candidate pool
  cheap, rather than trying to shrink the pool itself, which was never actually oversized. No
  further action needed here.

## Design / correctness coverage

- **The real C++ implementation has not been stress-tested the way the Python prototype was.** The
  Python prototype (`Scripts/ParticleMergePrototype/`) was validated across six 2D/3D domain
  topologies (corners, T-junctions, mixed same-rank/cross-rank ownership), multi-round convergence,
  and adversarially randomized processing order. The C++ implementation has so far only been
  exercised by the `Exec/Tests/ItoKMC/NNMergeSmokeTest` smoke test: uniform random sampling on a
  single-level, single-patch-per-corner-ish 32x32 grid, 1 and 4 MPI ranks. It has NOT yet been
  pushed through a deliberately adversarial scenario (3+ patch corners, mutual-match deadlock
  cases, boundary-exposed particles specifically) the way the prototype was. Worth doing before
  trusting this in a real physics run at scale.

- ~~Multi-level AMR (coarse-fine ghost masks) is untested.~~ **DONE, and it found a real bug**:
  tested by locally bumping `AmrMesh.max_amr_depth` to 1 in the smoke test and adding a per-level
  particle-count diagnostic (`ItoKMCGodunovStepper::postInitialize`'s `countPerLevel`/`printPerLevel`).
  Found that merges were only ever committing on the coarsest level. Root cause: `liveCellCount`
  was keyed by a bare `IntVect`, pooled across ALL levels, even though each cell key had been
  computed with THAT level's own dx -- a coarse-level cell key and an unrelated fine-level cell key
  could collide on the same `IntVect` by coincidence. Compounding this,
  `resolveTrivialTier()`/`judgeProposals()` hardcoded level-0's dx for every participant regardless
  of which level it actually lived on, so even without the collision, a fine-level particle's
  dynamic crowding recheck used the wrong cell key and almost always failed. Fixed by adding
  `NNMergeParticle::level`, introducing a level-qualified `NNCellKey = std::pair<int, IntVect>` for
  every pooled per-cell count, and threading a per-level `dx` array (`a_dxByLevel`) through
  `findNearestNeighborCandidates()`/`resolveTrivialTier()`/`judgeProposals()` so every cell key is
  always computed with the correct level's dx -- including for cross-level ghosts (a `GhostType::
  Coarse`/`Fine` ghost's true origin level is `lvl-1`/`lvl+1`, not the receiving patch's own level).
  `a_maxCellDistance`'s cross-level comparison (`nnMergeCrossLevelTooFar`) re-expresses both cell
  keys in the finer participant's cell-size units, reducing exactly to the original same-level
  definition when levels match. Verified: both levels now merge (e.g. level 0 5420->3407, level 1
  3413->2284 in one test), exact mass conservation preserved, on 1 and 4 MPI ranks.

- ~~Cross-rank proposals for ghost candidates were silently dropped.~~ **DONE, found via a
  cross-level merge diagnostic counter that stayed at exactly 0 for MPI runs with >1 rank while
  nonzero at 1 rank.** Root cause: `generateProposals()`'s own documented contract requires
  `a_particlesByID` to contain "every currently-alive local valid particle and every ghost", but
  `mergeNearestNeighborsRound()` only ever added NATIVE particles to that map. A ghost candidate
  can never resolve via the trivial tier (requires `candidateIsLocal`), so every cross-rank
  candidate (cross-level ones especially, since a cross-level candidate is *always* a ghost) had
  to go through `generateProposals()`, whose lookup then silently failed and dropped the proposal
  via a bare `continue` -- no error, just a lost merge. Masked itself as correct at 1 rank (every
  patch, hence every ghost's true owner too, is on rank 0, so the id is reachable via its owning
  patch's own native entry processed elsewhere in the same loop). Fixed by adding the missing
  `particlesByID[p.globalID] = p` entry to the ghost branch. Verified: the same multi-level config
  run at 1/4/8/16 ranks went from 0 cross-level merges at every rank count above 1 to a consistent
  nonzero count at every rank count.

- ~~Merged particles carried a stale, provisional `ownerRank`, causing mass creation at >1 rank
  and >1 round simultaneously.~~ **DONE.** `resolveTrivialTier()`/`judgeProposals()` set
  `NNMergeParticle::ownerRank` on a freshly merged particle to the committing parent's own rank --
  explicitly documented as provisional, since `placeMergedParticles()` determines the TRUE
  destination (level/patch/rank) from the merged position via `findDestination()` afterward. That
  provisional value was never corrected before being handed to the caller's scatter callback,
  which persists it into the new particle's own `rankID` column. Whenever the true destination
  rank differed from the provisional one (only possible with >1 rank), the newly created particle
  carried a wrong recorded owner rank -- silently wrong until that SAME particle later became a
  proposer in a SUBSEQUENT round: `judgeProposals()` routes an acceptance verdict to
  `NNMergeParticle::ownerRank`, so a stale value sends the verdict to the wrong rank, and the true
  owner never learns to remove the since-accepted source. The source survives un-consumed while
  its weight is ALSO folded into the new merged particle elsewhere -- mass is silently created.
  Reproduced by the user running the debug hook at 16 ranks with 2 rounds (mass conservation error
  of 142, vs. exact at 1 rank or 1 round). Fixed by correcting `ownerRank` to the current rank
  immediately before both scatter call sites in `mergeNearestNeighborsRound()`. Verified exact mass
  conservation at 16 ranks through 4 rounds and at ~125k particles.

- **Both bugs above were only found by testing at real physics-run scale (many ranks, many
  rounds) outside the committed smoke test, which only ever exercised 1/4 ranks and 1-2 rounds at
  small particle counts.** This is exactly the gap the "has not been stress-tested" item above
  already flagged, now with two concrete, confirmed examples of bugs that specifically required
  testing the CROSS PRODUCT of dimensions (rank count x round count), not each independently --
  1 round at 16 ranks was fine, 2 rounds at 1 rank was fine, only 2+ rounds at >1 rank failed.
  Worth keeping in mind for any future correctness pass: test combinations, not just axes.

- ~~`liveCellCount` double-counted a particle once per patch it was VISIBLE in, not once per
  distinct particle.~~ **DONE.** Found via the user reporting cells with fewer than the configured
  crowding threshold still merging. A single physical particle can be ghost-shipped to MULTIPLE
  patches this same rank owns at once (routine with tiled grid generation, which favors keeping
  neighboring patches co-located on a rank); `liveCellCount[key] += 1` was incremented once per
  occurrence in the per-patch gather loop (once per patch that saw it, native or ghost), inflating
  a crowded cell's count above its true population. `particlesByID` (keyed by `globalID`) has no
  such problem -- repeated occurrences collapse to one entry. Fixed by removing the inline
  increments and deriving `liveCellCount` in one pass, after the full gather loop, from
  `particlesByID`.

- ~~`judgeProposals()`'s dynamic crowding recheck only verified the TARGET's own cell, never the
  SOURCE's.~~ **DONE.** `resolveTrivialTier()`'s own commit condition requires BOTH participants'
  cells to dynamically recheck over threshold (`qCount > thresh && cCount > thresh`) at commit time,
  since an earlier merge within the same round can drain either side. `judgeProposals()` (the
  cross-rank propose/judge/verdict path) only ever checked the target's side -- a source particle
  whose own cell had since dropped to/below threshold could still be accepted, silently violating
  the crowding trigger on the source side. Fixed to check both symmetrically, defaulting an unknown
  source-cell count (this rank has no live ghost visibility into it) to 0 -- i.e. reject -- matching
  the existing conservative default already used for an unknown target count.

- ~~Cells could still merge below the configured crowding threshold, reproducing independent of AMR
  level, MPI rank count, and ghost/boundary exposure -- surviving both bugs immediately above.~~
  **DONE -- root cause was a fundamental key-comparison bug, not anything specific to levels, ranks,
  or ghosts.** `NNCellKey = std::pair<int, IntVect>` was used as the key of `std::map<NNCellKey,
  int> liveCellCount` (an ORDERED map, relying on `operator<`). `IntVect::operator<`
  (`Submodules/Chombo-3.3/lib/src/BoxTools/IntVect.H`) is a component-wise "strictly dominates in
  every direction" PARTIAL order -- built for box/geometry containment logic elsewhere in Chombo,
  not a valid strict-weak-ordering/total-order suitable for use as a map comparator.
  `std::pair<T1,T2>::operator<` calls `T2::operator<` directly on the second component, bypassing
  Chombo's own `std::less<IntVect>` specialization (which correctly uses `lexLT()`, a true
  lexicographic order) entirely. Result: two distinct, unrelated cells that don't dominate each
  other in every direction (the common case -- e.g. `(76,104)` vs `(80,101)`: neither `a<b` nor
  `b<a` holds) compare as neither-less-than-the-other, which `std::map` treats as an EQUAL key,
  silently pooling their live counts into one shared entry. This explains why the bug survived both
  fixes immediately above and reproduced identically at 1 rank, multi-level or single-level, deep
  in a patch's interior or at a boundary: it is a plain C++ comparator defect with nothing to do
  with any of those axes, only with which specific cells happen to collide in the map's underlying
  tree for a given insertion order -- confirmed empirically (a standalone reproduction using actual
  cell coordinates pulled from a real run's log collapsed 10 distinct cells into 2 map buckets; an
  independently-maintained ground-truth counter, checked at every commit, showed impossible
  negative populations that only make sense if multiple unrelated cells were sharing one counter).
  Fixed by adding `NNCellKeyHasher` (FNV-1a over the `IntVect`, reusing `LevelTiles::TileHasher`,
  folded with the level) and switching every `std::map<NNCellKey, int>` (`liveCellCount`, and the
  corresponding parameters in `findNearestNeighborCandidates()`/`resolveTrivialTier()`/
  `judgeProposals()`) to `std::unordered_map<NNCellKey, int, NNCellKeyHasher>`, which depends only
  on `IntVect::operator==` (correct, full component-wise equality) rather than `operator<`.
  Verified: a properly-hash-keyed ground-truth check (same idea as the standalone reproduction, but
  correct this time) found zero violations at 1/4/8 ranks, 1-3 rounds, on the exact multi-level,
  variable-block-size, refRat-4 configuration that originally reproduced the bug; exact mass
  conservation and `maxCenterlineDeviation ~ 0` throughout. `NNWalkCell`-keyed structures
  (`NNCellBuckets`, `NNSpatialIndex`) were never affected -- they already used `std::unordered_map`
  with `LevelTiles::TileHasher`, not `operator<`.
