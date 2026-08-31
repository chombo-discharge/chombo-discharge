# Mirrored cut-cell deposition — where this PR stands

**Read this first.** `PLAN.md` is the design (why), `IMPLEMENTATION.md` is the work order (what, in what
order). This file is neither: it is the running state of the branch, so a session picking the work up
cold knows what is done, what is broken, and what to do next without re-deriving any of it.

**PR** chombo-discharge#703 (draft) · **Branch** `mirror_deposition`, fork `rmrsk`, based on `main` at
`f44e5968` · **Issue** chombo-discharge#29 parts 1 and 3 · **Last updated** 2026-08-31

---

## 1. One-paragraph summary

Stages 1 and 2 of `IMPLEMENTATION.md` are written and pushed. Stage 1 (the `IrregularDeposition`
selector and the input migration) is **complete and fully verified** — every shipped regression is
bit-for-bit unchanged on all four mappings. Stage 2 (the per-band-cell surface data and the
shape-operator fit at regrid) is **written, building in 2-D and 3-D, and verified in 2-D**, with the
3-D torus case also passing; but there is **one reproducible 3-D segfault at 4 MPI ranks** which is
localized and not yet root-caused. Stages 3 and 4 have not been started.

Nothing selects the mirror yet, so the branch is inert with respect to every existing simulation.

---

## 2. Commits on the branch

Oldest last. Everything before `0fce9309` is throwaway design material under `Prototypes/`.

| commit | what |
|---|---|
| `876d9aee` | **Stage 2** — surface data at regrid. WIP, see §4 |
| `137bbc3c` | What Stage 1 taught the work order (`IMPLEMENTATION.md` §1.7) |
| `38b0101a` | **Stage 1** — `IrregularDeposition` selector + input migration |
| `0fce9309` | Re-resolve the plan against `main` after the rebase (`IMPLEMENTATION.md` §0.4) |
| `ae36b7d7` … `f0eb11ba` | `PLAN.md`, `IMPLEMENTATION.md`, the review records, and the Python harnesses |

The branch was rebased onto `main` at `f44e5968`, dropping its own unsquashed copy of PR #700 (which
`main` already carries as `92448b90`). Do not resurrect those three commits.

---

## 3. What is done, and what proves it

### Stage 1 — `IrregularDeposition` (complete)

Three independent `ItoSolver` booleans (`irr_ngp_deposition`, `redistribute`, `blend_conservation`)
collapsed into one selector `Native | NGP | Mirror | Redistribute | RedistributeBlended`, which is
what makes `mirror + redistribute` — a silent `1/kappa` error — unrepresentable rather than merely
discouraged. The enum stops at `AmrMesh`/`EBAMRParticleMesh`; the leaf kernels keep their bool. The 36
interpolation-side occurrences are untouched and `irr_ngp_interp` stays a bool.

Regressions bit-for-bit on all four mappings, compared at the **dataset** level (see §6):

| test | selector exercised | result |
|---|---|---|
| `BrownianWalker/DriftDiffusion` | `redistribute` | 201/201 identical |
| `ItoKMC/JSON` | `native`, McPhoto, both hardcoded stepper sites | 2/2 identical |
| `RadiativeTransfer/McPhoto` | McPhoto `native` | 21/21 identical |
| `DriftDiffusion`, forced | `ngp` | 201/201 identical |

All four rejection paths were exercised and produce the intended message: `mirror` + `deposition =
ngp` (names both keys), `mirror` alone (not implemented), an unknown selector, and a stale input.

### Stage 2 — surface data at regrid (written; 2-D verified, 3-D partial)

- `Source/Particle/CD_MirrorDeposition.H` — component layout, the **frame-free** world-frame shape
  operator, `gaussianCurvature` as the second invariant (never `det`), `reflect()` and its guards. No
  mesh, no I/O, directly unit-testable. Verified standalone in both dimensions, worst error **1.4e-14**
  over ~2600 checks.
- `PhaseRealm::defineMirrorSurfaceData` — three passes with an exchange between them, gated on the new
  `s_mirror_deposition` operator so realms that never mirror pay nothing.
- Five `AmrMesh.mirror_*` options; three threaded to `PhaseRealm` as scalars (the shape you chose over
  a parameters struct or a foreign-prefix `ParmParse` read).
- `Exec/Tests/ItoDiffusion/MirrorSurfaceData` — curvature checked against the analytic surface, binned
  by `kappa`, plus the `|S_c n_c| = 0` and `det S_c = 0` invariants.

2-D, all four cases pass, and the error grows with curvature as `PLAN.md` §3.3 predicts — comfortably
inside its 1–3% budget:

| case | mean `\|d(2H)\|/\|2H\|` |
|---|---|
| sphere `R ~ 25.6 dx` convex | 0.016% |
| sphere `R ~ 25.6 dx` concave | 0.016% |
| sphere `R ~ 6.4 dx` convex | 0.25% |
| sphere `R ~ 4.0 dx` convex | 0.70% |

3-D: the **torus passes at 4 ranks** — 3672 cut cells, all anisotropic, 2.3% mean and 6.2% in the
`kappa <= 0.05` bin. That is the case that matters most, because a sphere is umbilic and therefore
cannot distinguish a correct anisotropic shape operator from a wrong one.

Regressions are still bit-for-bit after Stage 2 (201/201), because nothing registers the operator.

---

## 4. THE OPEN DEFECT — 3-D segfault at 4 ranks

**Reproducer.** `Exec/Tests/ItoDiffusion/MirrorSurfaceData`, `regression3d.inputs`, sphere radius
`0.4`, `mpirun -np 4`. Exit 139 inside `PhaseRealm::defineMirrorSurfaceData`, reaching
`EBISBox::bndryCentroid`, which only **pass A** calls.

**Already ruled out — do not re-test these:**

| hypothesis | result |
|---|---|
| rank count alone | 1 and 2 ranks pass; 4 crashes |
| geometry alone | the 3-D **torus passes at 4 ranks**, so it is geometry- *and* rank-dependent |
| OpenMP race | identical crash with `OMP_NUM_THREADS=1` |
| the test's checker | crashes with the checker skipped entirely |
| ScanShop / geometry generation | crashes under both `geometry_generation = chombo-discharge` and `= chombo` |
| pre-existing code | `Exec/Tests/Geometry/Aerosol` passes at 4 ranks in 3-D |
| read windows leaving the `EBISBox` region | clipping to `EBISBox::getRegion()` did not fix it |
| hand-rolled irregular-set derivation | switching to the realm's own `m_vofIter` did not fix it |

**Where to pick up.** There is still no real stack. `XTRACXXFLAGS=-g` did **not** reach the compile
line (`addr2line` returns `??:?`), and `pout()` markers are lost to buffering when the process faults.
So the next step is one of:

1. Get debug symbols in properly — find the right variable for this build system, or compile
   `CD_PhaseRealm.cpp` by hand with `-g` and relink — then `addr2line -f -C -i`.
2. Bracket pass A with **unbuffered** writes (`fprintf(stderr, ...)` + `fflush`), which survive a
   fault, and print the level, box, `IntVect` and `numVoFs` immediately before the
   `bndryCentroid` call.

The suspicion worth testing first is that some VoF is irregular in the EB *graph* but absent from the
EB *data*'s boundary IVS, so the `BaseIVFAB` lookup inside `EBData::bndryCentroid` runs out of range.
If that is it, the fix is a guard in pass A plus a counter, not a change of algorithm.

---

## 5. What is next, in order

1. **Close the §4 defect.** Stage 2 is not done until 3-D passes at 4 ranks.
2. **Finish the Stage 2 acceptance suite** (`IMPLEMENTATION.md` §2.7): the cylinder case (the anisotropy
   check that a sphere cannot provide and the torus provides only incidentally), the sign-convention
   assertion, and the extension/tie-break determinism check across a different box decomposition.
3. **Stage 3** — the mirror pass itself: image container, scratch field, build/remap/deposit/add,
   guards, `ItoSolver` wired, plus the §5.2 acceptance suite.
4. **Stage 4** — `PLAN.md`'s PR B: redistribution moved centrally, duplicated solver code deleted.

Note for stage 4, from `IMPLEMENTATION.md` §1.7: `McPhoto::depositHybridDivergence` also carries an
unconditional `conservativeAverage` and `interpGhost` (`CD_McPhoto.cpp:1369-1370`). Deleting the
function must **re-home those two calls**, not drop them.

---

## 6. Environment and method — things that cost real time to learn

**"Bit-for-bit" must mean dataset-for-dataset.** Chombo plot files are **not** byte-reproducible run
to run: a baseline compared against *itself* differs in all 201 `DriftDiffusion` files. A plain `diff`
is therefore a broken test that reports a conversion bug in every stage. Use `h5diff -c -d 0`, under
which that same self-comparison is 201/201 identical.

**The `ngp` mapping is covered by no shipped regression.** Every `BrownianWalker` test selects
`redistribute` and `ItoKMC/JSON` selects `native`, so `ngp` survives only as the `CD_ItoSolver.options`
default and the `WireWire` *example*. Force it explicitly on both sides of any change that touches the
NGP cut-cell path, and keep a control confirming the run is actually sensitive to the setting.

**New input keys are mandatory everywhere.** The tree parses with `pp.get`, never `pp.query`, so a new
`AmrMesh.*` key has to be added to every shipped input in the same commit — Stage 2's five keys meant
**139 files**. That is the intended failure mode (a stale input hard-errors rather than silently taking
a default), but budget for it.

**`literalinclude` drift is partly silent.** Stage 1's three added `#include`s shifted **30** `:lines:`
ranges. Repair them from an exact per-file diff (`difflib` opcodes mapping old line to new), never by
matching excerpt text — the procedure `CLAUDE.md` §4 prescribes — then re-render every directive and
compare at word level. `sphinx-build -W` alone does **not** catch a still-valid-looking wrong excerpt.

**Building in a git worktree.** Worktrees do not get submodules, so `Submodules/{Chombo-3.3,EBGeometry,json}`
must be symlinked to the main checkout's, and `DISCHARGE_HOME` pointed at the worktree while
`CHOMBO_HOME` points at the main checkout's Chombo. Those symlinks are deliberately **not** committed.

**Available build configurations on this machine.** Chombo is built for `2d` **and** `3d`, both at
`Linux.64.mpic++.gfortran.OPTHIGH.MPI.DOUBLE.DOUBLE`. There is **no** assertion-enabled build, so
`CH_assert` is compiled out (`OPT=HIGH` disables it even with `DEBUG=TRUE`); exercising the
`|S_c n_c|` and multi-valued-cell asserts needs a Chombo built at `OPT=TRUE`, which does not exist here
yet.

**Do not name a loop variable `wit`.** `codespell` reads it as a misspelling of "with".

---

## 7. Defects found in the plan while building it

Recorded so they are not reintroduced. `IMPLEMENTATION.md` §1.7 and §0.4 carry the fuller versions.

- **`Vector<bool>` does not compile** (`PLAN.md`/`IMPLEMENTATION.md` §2.2). Chombo's
  `Vector<T>::operator[]` returns `T&`, which the `std::vector<bool>` proxy specialization cannot
  provide, and nothing else in the tree uses it. Held as `Vector<int>`.
- **`&sd(iv,0)` is not a component array** (§3.3). A `BaseFab` stores the component index **last** —
  all cells of component 0, then all cells of component 1 — so that pointer walks into the next
  **cell**, not the next component. `reflect()` takes a gathered contiguous array and documents why.
- **The `det(S_c)` invariant must be relative.** `det` is a product of `SpaceDim` entries of order
  `1/R`, so its roundoff floor grows as the surface tightens; an absolute bound passes on a gentle
  surface and fails on a sharp one for no reason but arithmetic.
- **§1.4's "the dead work disappears" is wrong for McPhoto** — see §5 above.
- **The behaviour-preserving default is `ngp`, not `native`** (§0.3 item 1).

---

**This directory is throwaway and must be deleted before the branch merges — this file included.**
