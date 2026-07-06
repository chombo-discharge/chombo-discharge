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

- **Multi-level AMR (coarse-fine ghost masks) is untested.** The smoke test runs with
  `AmrMesh.max_amr_depth = 0` (single level) specifically to keep the first end-to-end validation
  simple. The coarse-to-fine/fine-to-coarse particle ghost mask paths in
  `mergeNearestNeighborsRound()` (`maskC2F`/`maskF2C` in the gather loop, and the
  `isBoundaryExposed()` OR-across-three-masks logic) have never actually been exercised. Needs a
  multi-level test before relying on this near a refinement boundary.
