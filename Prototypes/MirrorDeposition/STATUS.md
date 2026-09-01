# Mirrored cut-cell deposition — where this PR stands

**Read this first.** `PLAN.md` is the design (why), `IMPLEMENTATION.md` is the work order (what, in what
order). This file is neither: it is the running state of the branch, so a session picking the work up
cold knows what is done, what is broken, and what to do next without re-deriving any of it.

**PR** chombo-discharge#703 (draft) · **Branch** `mirror_deposition`, fork `rmrsk`, based on `main` at
`cd3a7d29` · **Issue** chombo-discharge#29 parts 1 and 3 · **Last updated** 2026-09-01

---

## 1. One-paragraph summary

Stages 1 and 2 of `IMPLEMENTATION.md` are written. Stage 1 (the `IrregularDeposition` selector and the
input migration) is **complete and fully verified** — every shipped regression is bit-for-bit unchanged
on all four mappings. Stage 2 (the per-band-cell surface data and the shape-operator fit at regrid) is
**written, building in 2-D and 3-D, and verified in both** at every rank count tried. The 3-D segfault
that this file previously carried as the open defect is **root-caused and fixed** — it was a
use-after-free in the test's own checker, not a defect in the library (§4). Stages 3 and 4 have not
been started.

Nothing selects the mirror yet, so the branch is inert with respect to every existing simulation.

**Rebased onto `main` at `cd3a7d29` on 2026-09-01**, which is when the 3-D defect changed shape — see
§4. The rebase pulled in main's Chombo bump (#712, #717) and #713's `diffusion_grad_drift`; all 33
conflicts and the whole 2-D verification were redone against the new base, so every number in §3 was
re-measured on this tree rather than carried over.

---

## 2. Commits on the branch

Oldest last. Everything before `f22a05a3` is throwaway design material under `Prototypes/`. **These
SHAs are from the 2026-09-01 rebase onto `cd3a7d29`**; every commit was rewritten, so SHAs quoted in
older notes (`876d9aee`, `38b0101a`, `0fce9309`, …) no longer resolve.

| commit | what |
|---|---|
| *(this file's own commit)* | §4 rewritten: the 3-D defect root-caused and closed |
| `6ecf9064` | **The §4 fix** — the test was iterating a freed `IntVectSet` |
| `c1638b84` | How to obtain a comparison baseline (§6) |
| `2833da6b` | First version of this file |
| `fa0cba28` | **Stage 2** — surface data at regrid |
| `d9be993c` | What Stage 1 taught the work order (`IMPLEMENTATION.md` §1.7) |
| `7efa02ef` | **Stage 1** — `IrregularDeposition` selector + input migration |
| `f22a05a3` | Re-resolve the plan against `main` after the first rebase (`IMPLEMENTATION.md` §0.4) |
| `67a1b260` … `d14d439d` | `PLAN.md`, `IMPLEMENTATION.md`, the review records, and the Python harnesses |

The branch has been rebased twice: onto `f44e5968`, dropping its own unsquashed copy of PR #700 (which
`main` already carries as `92448b90`) — do not resurrect those three commits — and onto `cd3a7d29`.

**What the second rebase cost, in case it has to be redone.** 33 conflicts: `CD_ItoSolver.{H,cpp,options}`
(main's #713 added `diffusion_grad_drift` exactly where Stage 1 deleted `redistribute` and
`blend_conservation`; keep the new key, drop the two dead ones, and keep main's `s_eb_gradient`
registration alongside the branch's enum condition), 29 `.inputs` files with the same collision, and
`Docs/Sphinx/source/Solvers/Ito.rst`.

The `.rst` is the one with a trap in it. Five `:lines:` directives conflicted, but **ten more merged
cleanly while carrying the branch's stale numbers** — a green merge that silently renders the wrong
excerpts. All 15 have to be rebuilt from *main's* numbers mapped through `difflib` opcodes into the
rebased header, then verified word-for-word against what main rendered. Also re-run
`pre-commit run format-input-files --all-files` afterwards: main's `diffusion_grad_drift` line uses a
narrower comment column than the branch's reflowed one, so eight input files need realigning.

And do not `git add -A` mid-rebase: it commits the `Submodules/*` symlinks from §6 as typechanges.

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

### Stage 2 — surface data at regrid (written; verified in 2-D and 3-D)

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
inside its 1–3% budget. **Re-measured on the rebased tree, 2026-09-01**, at both 1 and 2 ranks:

| case | cut cells | mean `\|d(2H)\|/\|2H\|` | `kappa <= 0.05` bin |
|---|---|---|---|
| sphere `R ~ 25.6 dx` convex | 204 | 0.016% | 0.021% |
| sphere `R ~ 25.6 dx` concave | 204 | 0.016% | 0.023% |
| sphere `R ~ 6.4 dx` convex | 52 | 0.25% | 0.24% |
| sphere `R ~ 4.0 dx` convex | 28 | 0.70% | no cells in bin |

Worst `|S_c n_c|` 1.3e-15, worst `|det S_c|/|2H|^D` 8.4e-17, no planar fallbacks, no missing data.

The 1-rank and 2-rank runs agree **to every printed digit**, which is the extension/tie-break
determinism check of `IMPLEMENTATION.md` §2.7 passing for these configurations — the Euclidean-nearest
argmin with the lexicographic tie-break really is decomposition-independent here. It is not yet checked
across *different box sizes*, which is the harder half of that item.

`0 anisotropic` in 2-D is the correct answer, not a failure: a 2-D surface is a curve with one
principal curvature.

3-D, after the §4 fix, at **1, 2, 3, 4, 5, 6 and 8 ranks** with no no-data cells at any of them:

| case | cut cells | mean `\|d(2H)\|/\|2H\|` | `kappa <= 0.05` bin |
|---|---|---|---|
| sphere `R = 0.4` | 3056 | 0.11% | 0.30% |
| torus | 3424, **all anisotropic** | 3.5% | 12.9% |

The torus is the case that matters most, because a sphere is umbilic and therefore cannot distinguish a
correct anisotropic shape operator from a wrong one. Its `kappa <= 0.05` figure sits at 12.9% against a
15% bound — passing, but with less headroom than anything else here, so treat that bound as load-bearing
rather than generous.

Existing regressions still run after Stage 2, because nothing registers the operator: `DriftDiffusion`
2-D completes under all four implemented mappings (201 plot files each), and `native` versus
`redistribute` differ in all 201 files, which is the control confirming the selector actually reaches
the deposition path. A **full** bit-for-bit comparison against `main` has *not* been redone since the
rebase — main's #713 changes Ito results by design, so it needs a fresh baseline (§6).

---

## 4. THE 3-D SEGFAULT — root-caused and FIXED (2026-09-01)

**It was a use-after-free in the test's checker, not a defect in the library.** One line:

```cpp
// WRONG -- getIrregIVS returns by value, so the temporary dies at the end of the
// for-init-statement and the iterator walks freed memory from its first ok().
for (IVSIterator ivsIt(ebisbox.getIrregIVS(dbl[it()])); ivsIt.ok(); ++ivsIt) {
```

`IVSIterator` holds `const TreeIntVectSet* m_ivs` plus a `Vector<const TreeIntVectSet::TreeNode*>` — raw
pointers **into** the set. Binding it to a temporary is a dangling read, not a copy. The fix is to name
the set:

```cpp
const IntVectSet irregIVS = ebisbox.getIrregIVS(dbl[it()]);
for (IVSIterator ivsIt(irregIVS); ivsIt.ok(); ++ivsIt) {
```

Every other one of the 21 `IVSIterator` uses in the tree already binds a named set — including
`CD_PhaseRealm.cpp:675`. This was unique to the new test.

**Everything the symptom did is explained by it.** The iterator yielded garbage `IntVect`s, so: in the
`OPTHIGH` build a garbage `VolIndex` reached `EBISBox::volFrac` and dereferenced null; in the assertion
build the same garbage tripped `BaseFab::operator()`'s bounds check first. Garbage cells that happened
to land in range read garbage `status` values, which is where the phantom "no data" cut cells came
from. And because it read freed heap, the failing cell moved run to run, the fault vanished when probe
calls perturbed the allocation, and it was sensitive to rank count and geometry — which is also why
main's Chombo bump appeared to "move" the crash from the sphere to the torus. It never moved; the heap
did.

**Verified fixed.** Both geometries, all rank counts, `OPT=HIGH` and the assertion build:

| case | ranks | result |
|---|---|---|
| sphere `R = 0.4` | 1, 2, 3, 4, 5, 6, 8 | pass, 3056 cut cells, **0 no data** at every rank count |
| torus | 1, 2, 4, 8 | pass, 3424 cut cells, **all anisotropic**, 0 no data |
| torus, `OPT=TRUE DEBUG=TRUE` | 4 | pass, no assertion |
| 2-D, all four cases | 1, 2 | unchanged, still passing |

Torus at 4 ranks: mean `|d(2H)|/|2H|` 3.5% (bound 5%), 12.9% in the `kappa <= 0.05` bin (bound 15%),
`|S_c n_c|` 1.5e-15, `|det S_c|` 4.1e-16. The anisotropic path is now verified at every rank count
tried, which it was not before.

**The test now fails on it.** `numNone > 0` used to be printed in the summary line and never checked,
so a cut cell with no surface data passed silently. It calls `expect()`-equivalent now and fails.

### What this cost, and what to reuse

The previous session spent its whole budget on this and recorded a wrong conclusion, because two
tooling beliefs were wrong. Both are worth keeping:

- **`addr2line` on the `OPTHIGH` build works.** It returns `??:?` for the *line*, which reads like
  total failure, but the binary still carries a symbol table, so
  `addr2line -f -C -i -p -e <exe> <offset>` names the *function*. Feed it the `+0x...` offsets straight
  out of the OpenMPI backtrace. No rebuild, no `XTRACXXFLAGS`. That alone would have shown
  `checkSurfaceData` rather than `defineMirrorSurfaceData` and saved the entire hunt.
- **The assertion build is what closes it.** `make DIM=3 MPI=TRUE OPT=TRUE DEBUG=TRUE USE_HDF=TRUE` —
  `OPT=TRUE`, **not** `HIGH`, which compiles `CH_assert` out even with `DEBUG=TRUE`. It turned an
  address-not-mapped into `BaseFabImplem.H:320: Assertion 'm_domain.contains(a_p)' failed` plus a
  `-g` stack naming the exact source line. It exists in the tree now; §6 used to list it as missing.

The general lesson, worth applying before instrumenting anything again: a fault that **moves when you
add an unrelated statement** is memory corruption, not a bad index. Reach for assertions and a
sanitizer at that point rather than more `fprintf`.

---

## 5. What is next, in order

1. **Finish the Stage 2 acceptance suite** (`IMPLEMENTATION.md` §2.7). What §4 left done: the
   extension/tie-break determinism check passes across *rank counts* — 1-through-8-rank runs agree to
   every printed digit in both dimensions. What remains: the same check across **different box sizes**
   (`max_block_size`), which is the harder half and the one that would catch a decomposition-dependent
   tie-break; the **cylinder** case (the anisotropy check a sphere cannot provide and a torus provides
   only incidentally); and the **sign-convention** assertion.
2. **Re-establish the bit-for-bit baseline.** It has not been redone since the rebase, and main's #713
   changes Ito results by design, so the old benchmark files are useless. Recipe in §6.
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

**How to get a baseline to compare against.** The comparison needs a *second* worktree checked out at
`origin/main`, not a stash: the baseline has to run its own pre-migration inputs (which still say
`irr_ngp_deposition` / `redistribute`), and those only exist alongside the matching code. Recipe:

```
git worktree add .claude/worktrees/mirror-baseline origin/main --detach
# symlink its Submodules the same way (see below), then build and run the same test in both trees
# and compare with:  h5diff -c -d 0 <base>/plt/<file> <new>/plt/<file>
```

Sanity-check the comparison before trusting it: run the baseline **twice** and diff it against itself.
That must come out 201/201 identical at the dataset level; if it does not, the harness is wrong, not
the code. And check the comparison is *sensitive* — an `ngp` run and a `redistribute` run of
`DriftDiffusion` must differ in all 201 files, or a passing comparison means nothing.

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

Two things destroy them, and both leave *empty directories* behind rather than an obvious error, so the
next build fails with a confusing missing-`lib` message rather than a missing-symlink one:

- `git commit` — pre-commit's stash/restore round-trip replaces them.
- any `git checkout`/`--amend` that restores the gitlinks.

Re-create them (`rmdir` then `ln -s`) after either. Check with `ls -la Submodules/`, not `ls`.

**Available build configurations on this machine.** Chombo is built for `2d` **and** `3d` at
`Linux.64.mpic++.gfortran.OPTHIGH.MPI.DOUBLE.DOUBLE`, and — since 2026-09-01 — for `3d` at
`Linux.64.mpic++.gfortran.DEBUG.OPT.MPI.DOUBLE.DOUBLE` as well. The second one is the assertion-enabled
build, and §4 is the argument for keeping it: `make DIM=3 MPI=TRUE OPT=TRUE DEBUG=TRUE USE_HDF=TRUE`.
`OPT=TRUE`, **not** `HIGH` — `HIGH` compiles `CH_assert` out even with `DEBUG=TRUE`. The 2-D assertion
build still does not exist.

**`addr2line` works on the optimised build.** It returns `??:?` for the line number, which looks like
failure, but the symbol table is there: `addr2line -f -C -i -p -e <exe> <offset>` names the function.
Feed it the `+0x...` offsets from the OpenMPI backtrace. See §4 for what assuming otherwise cost.

**A fault that moves when you add an unrelated statement is memory corruption**, not a bad index. Stop
instrumenting and reach for the assertion build or a sanitizer. See §4.

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
- **`IVSIterator` must not be bound to a temporary `IntVectSet`** (§4). It holds raw pointers into the
  set's tree nodes, so `IVSIterator it(ebisbox.getIrregIVS(box))` dangles immediately. Name the set.
  This one is not a plan defect — it is a defect the plan never thought to warn about, and it cost a
  session.

---

**This directory is throwaway and must be deleted before the branch merges — this file included.**
