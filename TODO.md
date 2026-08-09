# TODO — particle merger cleanup (PR #684)

Working handoff note. **Delete this file before merging** — it is scratch, not project documentation.

Branch `merger-per-cell-cleanup`, PR
[#684](https://github.com/chombo-discharge/chombo-discharge/pull/684), based on upstream `main` at
`3d594509`. Everything below is committed and pushed.

## Done

Six commits, each verified byte-identical against its own oracle before landing:

| Commit | What |
|---|---|
| `f03108eb` | Deleted `NNParticleLocation` (a field-for-field copy of `LevelTiles::LevelAndBox`) and `cellLess` |
| `00a0b8ec` | kd per-cell data (histogram, live quota) → `EBAMRFAB`; exposure cache → `BaseFab`; added `AmrMesh::allocate(EBAMRFAB&, …)` |
| `a758aa4e` | Rebased 9 `literalinclude` ranges after the header shift |
| `d548957d` | Folded `mergeKDPatch` into `mergeKDInterior` (237 lines → 38) |
| `153398c4` | Moved 5 carve-protocol structs into `mergeKDCarve`'s body (133 lines off the header) |
| `da7b6d63` | Unified `KDParticle`/`NNMergeParticle` into `MergeParticle` (153 refs) |

Net: ten types removed, three hand-rolled per-cell stores replaced, ~200 lines of duplicated merge
body gone. KD header 7 types → 1 and 535 → 398 lines; NN header 14 types → 12.

## Next

1. **CI has not run on the final head** (`da7b6d63`). Check it. If OpenMP `Linux-GNU` jobs fail with
   exit 132 (SIGILL) across many tests, that is the known cached-binary flake — re-run the failed
   jobs before investigating. It cleared on identical code during #681.
2. **Update #682** to record that its `liveCellCount`/`NNCellBudget` conversion was attempted and
   rejected on merit, so nobody re-attempts it from the issue text.
3. Optional: `@claude review` on the PR.

## Decisions taken — do not re-litigate without new information

- **`a_liveCellCount` and `NNCellBudget` stay `std::unordered_map`.** Converting them to mesh data was
  written and abandoned. They are rank-global, cross-level state mutated inside the merge loop.
  Patch-scoping them needs a patch index on the pooled struct, a new MPI wire field so the judge can
  learn a remote cell's count, and introduces a cross-patch staleness hazard: a cell visible to
  several patches holds several independent copies of its count, and a commit decrements only one.
  Not a conservation bug (`consumedIDs` is rank-wide) but it permits over-merging a shared cell.
- **The two nearest-neighbour round drivers stay separate.** `nn_pair_tree`/`nn_pair_hash` already
  share a body via the `Cloud` template parameter. Folding in the per-cell driver needs either a new
  policy abstraction or three `if constexpr` branches in a 200-line orchestrator — it relocates
  complexity rather than removing it.
- **`EBAMRFAB`, not `EBAMRCellData`.** The latter needs a `phase` the mergers do not have, and its EB
  machinery would be discarded through `getFArrayBox()`. These are counts: interpolating or coarsening
  a particle count is not a defined operation. Switching later is one line per site.
- **`exposureCache` stays a `BaseFab`**, not a hierarchy-wide holder — it is patch-local scratch with
  no AMR meaning, and #685 proposes deleting it in favour of a regrid-lifetime mask anyway.

## How to verify anything here

- **Oracle is the integer field, never byte-identity of the default output.**
  `Driver.max_steps=5 Driver.plot_interval=1 ItoSolver.plt_vars=part`, then `h5diff` each step.
  Particles-per-cell is integer-valued so round-off cannot move it. The default `phi` differs at ~1
  ULP from particle ordering and is useless as a criterion.
- **Clean-build both sides** (`make libclean`). Incremental rebuild after a branch switch leaves stale
  objects and produced two phantom regressions during #681.
- **Determinism control first**: same binary twice. Cheap, and it invalidates everything downstream if
  it fails.
- **Build `OPT=TRUE`, not `OPT=HIGH`** — the latter defines `NDEBUG` and compiles every `CH_assert` to
  `(void)0`, so assertion-enabled runs are not what they appear.
- `verify_conservation` aborts the run, so `exit=0` already proves weight conservation at every merge.
- Free oracles: `kd_patch` and `kd_carve` never enter the NN merger, so NN-side changes must leave
  them byte-identical. A pure type change must leave *all* families byte-identical.
- Watch `literalinclude` drift. `Ito.rst`/`Particles.rst`/`MeshData.rst` cite absolute `:lines:` into
  `CD_ItoSolver.H`, `CD_ParticleManagement.H`, `CD_ParticleSoA.H` and `CD_AmrMesh.H`. Some drift is
  **silent** — Sphinx exits 0 rendering the wrong excerpt. Verify by comparing excerpt *content*
  across the change, not by trusting arithmetic. `MergeParticle` is deliberately at the END of
  `CD_ParticleManagement.H` for this reason.

## Related issues

- **#682** — the audit this PR implements the tractable part of.
- **#683** — merger per-patch loops are serial; shared accumulators block OpenMP.
- **#685** — regrid-invariant, non-templated data that should be cached rather than re-derived.
- **#679** — point-to-point MPI exchange; its prototype has no owner yet.

Open question, unresolved and probably the most valuable: none of the merger work is measured beyond
8 ranks on one machine, and whether six merge algorithms is the right surface has not been decided.
An operator-class refactor (`define()` at regrid, like `EBAMRParticleMesh`) would give #685's
invariants and #679's rank sets a home, but does not reduce the algorithm count.
