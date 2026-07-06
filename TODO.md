# TODO — nearest-neighbor particle merge

Design and performance items that have been identified during development of
`Source/Particle/CD_NearestNeighborParticleMerge.H` / `CD_NearestNeighborParticleMergeImplem.H`
(and its ItoKMC integration), deliberately deferred rather than fixed immediately. Kept here so
neither of us has to re-derive them from scratch later.

## Performance

- **`findNearestNeighborCandidates()` is brute-force O(n^2) per patch.** It currently checks every
  local-valid particle against every other local-valid particle and every ghost, rather than using
  a spatial index (the original "kd-tree over the particles in the patch" framing this feature
  started from). Correct, but will become a real bottleneck as particles-per-patch grows. A real
  spatial index (mirroring the allocation-light, thread-local-scratch style already used by
  `ParticleManagement::buildEqualWeightKDLeaves()` in `CD_ParticleManagementImplem.H`, though that
  one builds a different kind of tree -- see `CD_NearestNeighborParticleMerge.H`'s own docs) should
  replace it. This is a pure implementation-detail swap: `findNearestNeighborCandidates()`'s
  signature and behavior are already specified independently of how the search is done, so nothing
  else needs to change.

- **`std::map` used throughout instead of `std::unordered_map`.** `particlesByID`, `liveCellCount`,
  `exposed`, `outgoingTargetOf`, and `myParticlesByID` in `mergeNearestNeighborsRound()` (and the
  corresponding parameters threaded into `resolveTrivialTier()`/`judgeProposals()`) are all
  `std::map` (red-black tree, O(log n) lookup, poor cache locality) on what is a hot path -- looked
  up per particle, per edge, per fallback attempt. Should be `std::unordered_map` throughout for
  O(1) amortized lookups. `IntVect` and `ParticleID` both need a hash function if not already
  provided elsewhere in the codebase (check before writing a new one).

- **`myParticlesByID` is a full copy of `particlesByID`, rebuilt every round.**
  (`mergeNearestNeighborsRound()`, the "5. judgeProposals()" section.) Built by iterating
  `particlesByID` and filtering out anything in `consumedIDs`, into a brand new map -- O(n log n)
  and a full memory duplication every call, purely to avoid checking `consumedIDs` at each lookup
  site inside `judgeProposals()`. Should be removed: pass `particlesByID` + `consumedIDs` through
  and check membership at the point of use instead.

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
