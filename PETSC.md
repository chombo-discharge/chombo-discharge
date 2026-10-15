# PETSc field solver — work plan

> **This file is a planning document carried by the PETSc stack branches. It is never merged into
> `main`.** See [Keeping this file out of main](#keeping-this-file-out-of-main).
>
> Picking this up in a new session? Read [Handoff](#handoff), then
> [Building with PETSc](#building-with-petsc), then [Design](#design) and [Code map](#code-map).

## Goal

An electrostatic field solver, `FieldSolverAMG`, that assembles the composite (all AMR levels at once)
multifluid cut-cell Helmholtz operator as a PETSc matrix and solves it with a PETSc Krylov method
preconditioned by algebraic multigrid (GAMG or BoomerAMG). It must be a drop-in alternative to
`FieldSolverGMG` — same inputs, same boundary-condition options, same outputs — so the two can be
compared on the existing Electrostatics tests.

## Where things stand

What `main` already has (merged before this stack was restarted):

| Piece | Location | PR |
|-------|----------|----|
| PETSc initialisation/finalisation through `ChomboDischarge::initialize` | `Source/Driver/CD_Initialize.*` | #597, #598, #607 |
| `PetscGrid`: AMR cell ↔ PETSc row map, `Vec` create/destroy/fill, Chombo ↔ PETSc copies, HDF5 dump | `Source/AmrMesh/CD_PetscGrid.*`, `CD_PetscAMRCell.*`, `CD_PetscDOF.H` | #599 |
| `AmrMesh::getPetscGrid(realm)`, built per realm on regrid | `Source/AmrMesh/CD_AmrMesh.*`, `CD_Realm.*` | #599 |
| Explicit coarse-fine ghost stencils (`EBMultigridInterpolator::getInterpolationStencilEB/Regular`) | `Source/AmrMesh/CD_EBMultigridInterpolator.H` | #600 |
| Dual-level least squares (`LeastSquares::computeDualLevelStencils`) | `Source/Utilities/CD_LeastSquares.H` | — |
| `Exec/Examples/Electrostatics/PETSc` example | — | — |

**`main` does not build with `USE_PETSC=TRUE`.** Two independent breaks, both fixed on `field_amg`:

- `CD_PetscAMRCell.cpp` calls `detail::push/pullParticleProperty`, which #656 (SoA particles) removed.
  Nobody noticed because CI never compiles PETSc code.
- `Exec/Examples/Electrostatics/PETSc/main.cpp` includes `CD_FieldSolverAMG.H`, which only exists on
  this stack.

The two original draft PRs were restarted on 2026-10-14 by merging `main` into them:

- **#601** (`field_amg`): `PetscStencil`, a `FieldSolverAMG` skeleton whose `solve` only round-trips a
  `Vec` through `PetscGrid`, the `Linux-PETSc` CI job, and doxygen coverage of the PETSc code
  (`CH_USE_PETSC` predefined; 0 warnings). `MFHelmholtzPetsc` and the scratch
  `ProfiledSurface/dev.inputs` were removed from it.
- **#602** (`helmholtz_petsc`, stacked on #601): `PetscGrid` accessors plus per-phase iterators over
  the cut cells that are PETSc rows, and a partial `MFHelmholtzPetsc` that re-derives single-level
  least-squares EB stencils. With `CH_USE_PETSC` predefined it has 54 doxygen warnings (10 in
  `CD_PetscGrid.H`, 44 in `CD_MFHelmholtzPetsc.H`) and three `#warning`s, so its documentation CI job
  fails.

Verified: `discharge-lib` and `Exec/Examples/Electrostatics/PETSc` build with `MPI=TRUE
USE_PETSC=TRUE` against Ubuntu's PETSc 3.19 — `field_amg` in 2D and 3D, `helmholtz_petsc` in 2D. The
2D example runs to completion on 4 ranks (the solve is a stub, so it reports non-convergence).

## Design

Decided 2026-10-14.

### Our own composite classes; the GMG solvers are not touched

`EBHelmholtzOp`, `MFHelmholtzOp` and the `*EBBC`/`*JumpBC` classes are single-level and operate on data:
they fill ghost cells, apply, then reflux with `EBReflux` (or coarsen fluxes in the reflux-free path).
Near a coarse-fine interface they handle EB stencils that reach under the finer level by overwriting the
covered coarse cells with the conservative average of the fine data (`EBHelmholtzOp::coarsenCell`), and
they reduce stencil order wherever `isStencilValidCF` rejects a stencil. None of this exists as
stencils over the composite grid, and exposing AMR interfaces in those classes would be as much work as
new classes while muddying the interface of the existing solvers. So the PETSc operator is assembled by
new classes, and nothing in `Source/Elliptic` that GMG uses changes.

What is reused is utility code that is not a solver interface: `PetscGrid` (row numbering),
`LeastSquares::computeDualLevelStencils` (the fit), the ghost-cell stencils of
`EBLeastSquaresMultigridInterpolator` (#600), and the EB geometry queries.

### Which cells a stencil on level L may reference

Every PETSc row is a valid cell: on its level, not covered by a finer level, not a ghost. A stencil
built on level L resolves its entries as follows.

| Entry | Becomes |
|-------|---------|
| Valid cell on L | Its PETSc row. |
| L cell covered by L+1 | **Never used.** The neighbourhood takes the actual L+1 cells instead. |
| L+1 cell | Its PETSc row (L+1 cells used this way are valid on L+1). |
| Coarse-fine ghost cell of L | The stencil `EBLeastSquaresMultigridInterpolator` uses to fill that ghost (`getInterpolationStencilEB/Regular`), in terms of L and L-1 cells, then resolved by this table. |
| L-1 cell | Only ever reached through a ghost stencil, never directly. |

So each least-squares fit is a **two-level problem** (L and L+1), solved with
`LeastSquares::computeDualLevelStencils` (fine = L+1, coarse = L), whose separate fine/coarse output
stencils keep equal `VolIndex`es on the two levels apart. L-1 enters only through ghost cells that
another class defines explicitly. A finished matrix row can still span three levels, but only via that
ghost substitution.

### Fluxes across coarse-fine faces

The flux through a coarse face on the coarse-fine interface is the area-weighted sum of the fluxes
through the fine faces that make it up. The assembler builds the coarse row with that sum in place of
the coarse face flux, so refluxing is part of the matrix rather than a post-step. Fine-side face fluxes
that touch ghost cells resolve through the table above.

### Things to verify when building this

1. **Ghost depth.** An EB stencil on L may reach into L's ghost cells only as far as the interpolator
   fills them (its ghost width). This is what `isStencilValidCF` guards in GMG; here it is an assertion,
   or the neighbourhood is shrunk to respect it.
2. **What ghost stencils reference.** Ideally only valid L cells and L-1 cells not covered by L. If they
   reference other ghosts, the substitution recurses; if they reference L-1 cells covered by L, those
   need their own rule.
3. **Ghost stencils agree with `interpGhostMG`.** The #600 stencils must reproduce what
   `interpGhostMG` computes in cut cells, so both solvers mean the same thing by "ghost".

### Consequence for comparing against GMG

Near coarse-fine interfaces the PETSc stencils use a better neighbourhood than GMG (real L+1 cells
instead of coarsened data, no order reduction from `isStencilValidCF`). The two discretisations
therefore deliberately differ in cut cells near coarse-fine interfaces. Comparisons against
`FieldSolverGMG` are to solver tolerance away from those cells, and convergence-based near them.

## Ground rules for the stack

- Every PR must **compile with `USE_PETSC=TRUE MPI=TRUE`** in 2D and 3D. CI does not build PETSc today
  (the `Linux-PETSc` job added on `field_amg` fixes that once it lands); until it does, build locally before pushing.
- Every PR must **leave the non-PETSc build and the regression suite bit-for-bit unchanged**. All new
  code sits behind `#ifdef CH_USE_PETSC`.
- Each PR adds a test that exercises what it adds (an `Exec/Tests/Electrostatics/PETSc*` case or a
  unit-style executable), not only the final solver.
- New code does not change the GMG solver classes (see Design).
- Data-structure rules from `CLAUDE.md` §0 apply: per-cell data in mesh containers, no new type
  without naming the existing one considered, no `std::map` keyed on cells.
- `#warning` reminders in the WIP code are replaced by entries in this file as each PR lands; no PR
  is merged with a `#warning` in it.

## The stack

Each row is one PR. Each PR's branch is cut from the one above it; its GitHub base is that branch until
the one above merges, then `main`.

Agreed 2026-10-14: #601 and #602 are replaced by a single PR, cut fresh from `main`, that is the
mapping between chombo-discharge and PETSc. #601 and #602 are closed; their branches `field_amg` and
`helmholtz_petsc` stay in the fork for reference (`helmholtz_petsc` holds the partial
`MFHelmholtzPetsc` and the `PetscGrid` accessors still to be ported).

| # | Branch | PR | Scope | Status |
|---|--------|----|-------|--------|
| 1 | `petsc_mapping` | draft (replaces #601, #602) | Mapping: CI, `PetscGrid` accessors, `PetscStencil`, cell substitution rules, test | Infrastructure in; accessors, substitution and test to do |
| 2 | `petsc_regular` | — | Single-level, single-phase, no-EB operator + KSP solve; domain BCs; `FieldSolverAMG` skeleton | Not started |
| 3 | `petsc_amr` | — | Composite operator: coarse-fine fluxes as sums of fine fluxes | Not started |
| 4 | `petsc_eb` | — | Single-phase EB fluxes, dual-level least squares | Not started |
| 5 | `petsc_jump` | — | Multifluid jump condition (dielectrics, surface charge) | Not started |
| 6 | `petsc_fieldsolver` | — | Full `FieldSolverAMG`: coefficients, regrid, options, E-field, tests, docs | Not started |
| 7 | `petsc_perf` | — | Matrix reuse, stencil caching, preconditioner tuning, scaling study | Not started |

### 1 — Mapping between chombo-discharge and PETSc (`petsc_mapping`)

Done on `field_amg` and carried over:

- [x] `Linux-PETSc` CI job: 2D and 3D, `MPI=TRUE USE_PETSC=TRUE`, Ubuntu `petsc-dev`; builds
      `discharge-lib` and runs the `PetscMapping` test; `CI-passed` depends on it. Skips draft PRs, like the other
      build jobs.
- [x] `CH_USE_PETSC` in doxygen `PREDEFINED`; all warnings this exposed in code on `main` fixed.
- [x] Build fix for `PetscAMRCell` after #656.
- [x] `PetscStencil` documented, including why it is not a `VoFStencil` (entries are global PETSc rows,
      which is what `MatSetValues` takes; a `VolIndex` is only meaningful with a level and box) and why
      `add` is a linear scan (stencils have tens of entries, ~100 at most).

Added since:

- [x] `PetscGrid` accessors from #602: `getNumPhases`, `getMFLevelGrids{,CoFi,FiCo}`,
      `getAMRToPetsc{,CoFi,FiCo}`, `getVoFIterator(phase)`. Kept only the `mutable LayoutData<VoFIterator>`
      (the codebase idiom, e.g. `EBHelmholtzOp`): multi-valued cells abort, so the cut-cell `IntVectSet` was
      redundant. Dropped `getGrids()` (no callers, duplicates `getMFLevelGrids`) and the `#warning`.
- [x] Global-row fix: rows are offset by `m_localRowBegin` and `m_amrToPetsc` is exchanged afterwards, so
      ghost cells and the CoFi/FiCo copies hold correct global rows on any rank count.
- [x] Substitution: `PetscGrid::getPetscStencil(VoFStencil, level, din, phase, interpolator)`. Valid cell
      (own or neighbouring box) → row; CF ghost → the interpolator's stencil, fine part as level-L rows and
      coarse part as level-(L-1) rows via the CoFi map; anything else aborts with the reason. Two-level
      overload takes `(fine, coarse)` from `computeDualLevelStencils`, fine part resolved through the FiCo map.
      Recursion turned out unnecessary: ghost stencils only reference valid cells (verification item 2).
- [x] `EBMultigridInterpolator::getGhostStencil(vof, din)`: returns the stencil `coarseFineInterp` actually
      uses for a ghost cell (EB stencil if the cell is in the EB ghost set, otherwise the regular one).
- [x] Verification items 2 and 3, as the test below. Results:
  - A 3D bug on `main`: `getInterpolationStencilRegular` summed the mixed-derivative term over the
    tangential displacements, while `regularCoarseFineInterp` (what GMG runs) multiplies them. Fixed in the
    stencil; with the old formula the test fails in 3D (`interpGhostMG` 0.18400 vs stencil 0.18343).
  - Item 2 found a bug on `main` in `CoarseInterpQuadCF::defineStencils`/`defineMixedDerivStencils`: the
    valid-cell set is in coarse index space but the fine boxes were subtracted without coarsening, so the
    regular CF interpolation could difference through coarse cells covered by the fine level. Fixed by
    coarsening the boxes. **This changes GMG results at coarse-fine interfaces** (same convergence rate on
    `RodSphere`; residual history differs), so AMR benchmark files change. Kept as its own commit so it
    can go to `main` as a separate PR.
  - With both fixes, every ghost stencil references only rows, and reproduces `interpGhostMG` to round-off.
- [x] Test `Exec/Tests/Electrostatics/PetscMapping` (RodSphere setup, two phases, three levels in 2D): round
      trip exact; ghost cells overlapping neighbour boxes hold the neighbour's global row; two-level stencils
      find the level-(L+1) rows under covered cells (average equals `conservativeAverage`); every CF ghost
      cell's PETSc stencil reproduces `interpGhostMG` (max error 7e-16). Passes on 1 and 4 ranks in 2D and
      3D. The `Linux-PETSc` CI job builds and runs it on 2 ranks instead of building the example. It is not in
      `Electrostatics.ini`, since `tests.py` also runs non-PETSc builds.
- [x] `Exec/Examples/Electrostatics/PETSc`: switched to `FieldSolverGMG`, so `main` builds with PETSc
      without `FieldSolverAMG`. Builds in 2D and 3D; the 2D example runs on 4 ranks.
- [x] `FieldSolverAMG` skeleton is not in this PR; it moves to PR 2 (source on `field_amg`).
      `MFHelmholtzPetsc` from #602 is not carried over (its neighbourhood gathering is single-level).

### 2 — Regular operator and first solve (`petsc_regular`)

- [ ] Assemble `alpha*a*phi + beta*div(b grad phi)` on a single level, single phase, no EB, into an
      `MPIAIJ` matrix via `PetscGrid` rows; preallocate from stencil widths.
- [ ] Domain BCs (Dirichlet / Neumann / Robin, constant and function-valued) moved into the RHS, in the
      new classes.
- [ ] `KSPSolve` with `-ksp_type cg|gmres -pc_type gamg|hypre` taken from options; residual reported in
      the same norm as `FieldSolverGMG`.
- [ ] `FieldSolverAMG` skeleton from #601 lands here, solving this case.
- [ ] Test: manufactured solution on a box, second-order convergence over three resolutions;
      `FieldSolverGMG` and PETSc agree to solver tolerance.

### 3 — Composite AMR operator (`petsc_amr`)

- [ ] Coarse cells adjacent to a coarse-fine interface: the coarse face flux is replaced by the
      area-weighted sum of the fine face fluxes (refluxing folded into the matrix).
- [ ] Fine face fluxes that touch ghost cells go through the PR 1 substitution.
- [ ] Covered coarse cells are not rows; after the solve their values are set by `conservativeAverage`.
- [ ] Test: 2- and 3-level static refinement without EB matches `FieldSolverGMG` to solver tolerance.

### 4 — Embedded boundaries, single phase (`petsc_eb`)

- [ ] Composite neighbourhood selection for cut cells on L: valid L cells, L ghosts (within the
      interpolator's ghost width, verification item 1), and real L+1 cells instead of covered L cells.
      `EBGradient` near coarse-fine interfaces is the precedent.
- [ ] EB flux stencils for Dirichlet / Neumann / Robin from `computeDualLevelStencils`, mapped with the
      PR 1 rules. Fallbacks (order reduction, neighbourhood growth) are decided here, in the new classes.
- [ ] Test: `Exec/Tests/Electrostatics/RodSphere` (single phase) agrees with GMG away from coarse-fine
      cut cells; convergence study near them.

### 5 — Jump condition (`petsc_jump`)

- [ ] Two-phase rows: each multifluid cut cell has one row per phase; the jump condition
      `eps_g dphi_g/dn - eps_s dphi_s/dn = sigma` couples them. Gradient stencils from the same
      dual-level machinery as PR 4.
- [ ] Both `natural` and `saturation_charge` jump types.
- [ ] Test: `ProfileCircle`/`ProfileSquare`/`MechShaft` agree with GMG (same caveat as PR 4).

### 6 — `FieldSolverAMG` (`petsc_fieldsolver`)

- [ ] `setupSolver`, `setPermittivities`/`setSolverPermittivities`, `regrid`, `parseOptions` with the
      same option names as `FieldSolverGMG` (BCs, `jump_bc`, `kappa_source`, filtering, plot vars) plus
      `FieldSolverAMG.ksp_*` / `pc_*`.
- [ ] Rebuild the matrix on regrid; on time-dependent permittivity, update coefficients only.
- [ ] Restore `Exec/Examples/Electrostatics/PETSc` on `FieldSolverAMG`; `ProfiledSurface` can use it;
      `.ini` regression entries gated on a PETSc build.
- [ ] Sphinx: a `FieldSolverAMG` section next to `FieldSolverGMG`, with PETSc build instructions.

### 7 — Performance (`petsc_perf`)

- [ ] Cache EB/jump stencils so a coefficient change does not redo least squares.
- [ ] `MatSetValues` batching, preallocation accuracy, `MAT_NEW_NONZERO_ALLOCATION_ERR` on.
- [ ] Weak/strong scaling against `FieldSolverGMG` on `MechShaft3d`.

## Open decisions

1. **Fallbacks for EB stencils near coarse-fine interfaces** (PR 4): what to do when the composite
   neighbourhood is too small for the requested order, now that GMG's choices are not inherited.
2. **Preconditioner default**: GAMG (always in PETSc) or BoomerAMG (needs hypre).

Settled 2026-10-14: the discretisation is assembled by new classes rather than by translating GMG
stencils or extending the GMG/BC classes (see Design); #601/#602 collapsed into `petsc_mapping`.

## Handoff

State on 2026-10-14:

- `petsc_mapping` (#741, draft, on top of `main` at `0800b4d45`) is the live branch. Its PR description
  is filled in from the template; the checklist is unticked.
- #601 and #602 are closed. Their branches `field_amg` (`dddc6ffe`) and `helmholtz_petsc`
  (`25536f6b`) stay in the fork for reference; both are synced with `main` at `0800b4d45`.
- The `PetscGrid` accessors and the global-row fix from `helmholtz_petsc` have been ported into
  `petsc_mapping` by hand. Nothing else on those branches is needed for PR 1.

## Building with PETSc

Recipe that worked on Ubuntu 24.04 (the cloud container):

```bash
apt-get install -y petsc-dev libopenmpi-dev openmpi-bin gfortran csh   # csh: Chombo's build scripts
git submodule update --init          # if Submodules/json fails, fetch the pinned commit inside it by hand
export DISCHARGE_HOME=$PWD CHOMBO_HOME=$PWD/Submodules/Chombo-3.3/lib PETSC_DIR=/usr/lib/petsc
sed 's/USE_PETSC      = FALSE/USE_PETSC      = TRUE/' Lib/Local/Make.defs.MPI.GNU > $CHOMBO_HOME/mk/Make.defs.local
make -j4 discharge-lib DIM=2         # the first build (Chombo included) is slow
cd Exec/Examples/Electrostatics/PETSc && make -j4 DIM=2
mpirun --oversubscribe -np 4 ./main2d.*.ex example2d.inputs Driver.max_steps=0
```

- `Make.defs.local` goes in `$CHOMBO_HOME/mk/`, not `Lib/Local/` (that is what CI does too).
- `USE_PETSC=TRUE` requires `MPI=TRUE`; `CD_Initialize.H` errors out otherwise.
- PETSc's `petscvariables` brings in hypre, so `-pc_type hypre` is available with this package.
- A clone made with `--depth` breaks `git merge-base`; run `git fetch --unshallow` first.

## Code map

Facts about existing code that the to-do items build on. Where `main` and
`helmholtz_petsc` differ, both are given.

### Rows: `PetscGrid`

- One PETSc row per (valid cell, phase). Valid = on its level, not covered by a finer level, not
  covered by the EB. Multi-valued cells abort (`definePetscRows`).
- `m_amrToPetsc[lvl]` is a `LevelData<BaseFab<PetscAMRCell>>` on level `lvl`'s grids with ghosts.
  `PetscAMRCell::getPetscRow(phase)` gives the global row (`-1` if not a row), plus the flags
  `isCoveredByFinerGrid`, `isCoarCF`, `isFineCF`, `isGhostCF`, `isDomainBoundaryCell`. Ghost cells
  overlapping a neighbouring box hold that box's rows; ghost cells outside the level keep `isGhostCF`.
  (On `main` the row is rank-local and ghosts hold none; PR 1 fixes that.)
- `m_amrToPetscCoFi[lvl]` (`lvl < finest`): the level-`lvl` map copied onto *coarsen(level `lvl+1`)*
  grids. These share `DataIndex` with level `lvl+1`, so a fine box can look up the coarse rows under
  it. Note the offset: the CoFi data for fine level `l` lives at index `l-1`.
- `m_amrToPetscFiCo[lvl]` (`lvl > 0`): the level-`lvl` map copied onto *refine(level `lvl-1`)* grids,
  so a coarse box can look up the fine rows on top of it (the CF flux sum in PR 3).
- Accessors (PR 1): `getAMRToPetsc*`, `getMFLevelGrids*`, `getNumPhases`, and per-phase
  `getVoFIterator` over cut cells that are rows.

### Coarse-fine ghost stencils (#600)

- `EBMultigridInterpolator::getInterpolationStencilEB(ghostVoF, din)` returns a `VoFStencil` pair and
  `getInterpolationStencilRegular(ghostIV, din, dir, side)` a `DerivStencil` (IntVect-keyed) pair, each
  `(fine, coarse)`. The coarse part is in the **CoFi index space** (coarsened fine grids, same `din`),
  so resolve it through `m_amrToPetscCoFi[l-1]`.
- The only implementation is `EBLeastSquaresMultigridInterpolator`, which aborts if asked for a
  non-ghost cell. `getGhostStencil(vof, din)` (PR 1) picks whichever of the two `coarseFineInterp` uses.
- Regular ghost cells are one layer across each face (no corners); EB ghost cells reach `mg_interp_radius`
  deep near cut cells.

### GMG stencils (reference only)

The Design section rules out reusing these, but they are where to look when checking what GMG does
in the same cell, e.g. for the comparison tests:

- Interior regular cells: `EBHelmholtzOp` uses Fortran kernels, not stencils.
- Irregular cells: `EBHelmholtzOp::m_relaxStencils` and `m_centroidFluxStencil[dir]` (protected, no
  accessors).
- EB fluxes: `EBHelmholtzEBBC::getGradPhiStencils()`, `MFHelmholtzEBBC::getLeastSquaresBoundaryGradStencil`.
- Jump cells: `MFHelmholtzJumpBC::getGradPhiStencils()`.
- Configuration from inputs: `FieldSolverGMG::setupHelmholtzFactory`.

### Not carried over from #601/#602

- `MFHelmholtzPetsc` (#602): `computeEBDirichletStencil` gathers a single-level neighbourhood and
  returns an empty stencil, and every other `compute*Stencil` is empty. It also has a bug in the
  `hasFine` branch: `eblgFine = mflg.getEBLevelGrid(...)` uses the coarse grid, not `mflgFine`.
  Superseded by the dual-level design.
- `FieldSolverAMG` (#601, lands in PR 2): `solve()` round-trips a `Vec` filled with the rank number
  through `PetscGrid` and returns `false`. `setupSolver()` builds `MFHelmholtzPetsc` with hard-coded
  orders and weights (drop that), the permittivity setters are empty, the E-field computation is
  copied from `FieldSolverGMG`, and it guards with `#if CH_USE_PETSC` instead of `#ifdef`.

## Keeping this file out of main

`main` is updated by squash merges, so a file reaches `main` only if it is in the PR's final tree.

- `PETSC.md` lives on every stack branch.
- When a PR in the stack is ready to merge, its **last commit removes `PETSC.md`** (`git rm PETSC.md`);
  then squash-merge. `main` never sees the file.
- The next branch up does **not** merge that removal from its parent. After the parent lands,
  retarget the next PR to `main` and merge `main` into it; since `main` never touched `PETSC.md`, the
  file stays and there is nothing to resolve.
- Update the stack table above in the same branch that changes the status, and merge the change
  upward (`git merge <parent>` on each child) so every branch carries the current plan.
