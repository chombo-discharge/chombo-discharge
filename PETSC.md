# PETSc field solver — work plan

> **This file is a planning document carried by the PETSc stack branches. It is never merged into
> `main`.** See [Keeping this file out of main](#keeping-this-file-out-of-main).

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
| Explicit coarse-fine interpolation stencils (`EBMultigridInterpolator::getInterpolationStencilEB/Regular`) | `Source/AmrMesh/CD_EBMultigridInterpolator.H` | #600 |
| `Exec/Examples/Electrostatics/PETSc` example | — | — |

**`main` does not build with `USE_PETSC=TRUE`.** Two independent breaks, both fixed on `field_amg`:

- `CD_PetscAMRCell.cpp` calls `detail::push/pullParticleProperty`, which #656 (SoA particles) removed.
  Nobody noticed because CI never compiles PETSc code.
- `Exec/Examples/Electrostatics/PETSc/main.cpp` includes `CD_FieldSolverAMG.H`, which only exists on
  this stack.

The two original draft PRs, restarted on 2026-10-14 by merging `main` (80 commits ahead) into them:

- **#601** (`field_amg`): `PetscStencil` (row/weight list), a `FieldSolverAMG` skeleton whose `solve`
  only round-trips a `Vec` through `PetscGrid`, and the `MFHelmholtzPetsc` declaration.
- **#602** (`helmholtz_petsc`, stacked on #601): `PetscGrid` accessors plus per-phase iterators over
  the cut cells that are PETSc rows; `MFHelmholtzPetsc` with constructor/setters, the per-level loop
  for EB Dirichlet stencils, and the start of the neighbourhood gathering for one cut cell. Every other
  stencil function returns an empty stencil.

The sync merges were conflict-light (the only conflicts were the include-guard / Doxygen-style
reformat from #673 and the `BoxLoops::loop<...>` signature change from #643). Follow-up commits:

- moved the PETSc files to main's conventions (SPDX header, leading-star Doxygen, uppercase guards,
  clang-format);
- fixed the PETSc build (the `PetscAMRCell` break above; #601 also called `PetscGrid::isDefined()`,
  which only #602 added, so that accessor moved down into #601; missing explicit includes).

Verified after the sync: `discharge-lib` and `Exec/Examples/Electrostatics/PETSc` build in 2D with
`MPI=TRUE USE_PETSC=TRUE` against Ubuntu's PETSc 3.19, for both branches, and the example runs to
completion on 4 ranks (the solve is still a stub, so it reports non-convergence). 3D is not yet
checked.

## Ground rules for the stack

- Every PR must **compile with `USE_PETSC=TRUE MPI=TRUE`** in 2D and 3D. CI does not build PETSc today
  (PR 1 below fixes that); until it does, build locally before pushing.
- Every PR must **leave the non-PETSc build and the regression suite bit-for-bit unchanged**. All new
  code sits behind `#ifdef CH_USE_PETSC`.
- Each PR adds a test that exercises what it adds (an `Exec/Tests/Electrostatics/PETSc*` case or a
  unit-style executable), not only the final solver.
- Data-structure rules from `CLAUDE.md` §0 apply: per-cell data in mesh containers, no new type
  without naming the existing one considered, no `std::map` keyed on cells.
- `#warning` reminders in the WIP code are replaced by entries in this file as each PR lands; no PR
  is merged with a `#warning` in it.

## The stack

Each row is one PR. Each PR's branch is cut from the one above it; its GitHub base is that branch until
the one above merges, then `main`.

| # | Branch | PR | Scope | Status |
|---|--------|----|-------|--------|
| 1 | `field_amg` | #601 | Infrastructure: PETSc CI build, `PetscStencil`, `FieldSolverAMG` skeleton | Checklist done; needs retitle and PR description |
| 2 | `helmholtz_petsc` | #602 | `PetscGrid` accessors; VoF-stencil → PETSc-row translation | Synced; needs re-scope (see below) |
| 3 | `petsc_regular` | — | Single-level, single-phase, no-EB operator assembly + KSP solve; domain BCs | Not started |
| 4 | `petsc_amr` | — | Composite operator across coarse-fine interfaces | Not started |
| 5 | `petsc_eb` | — | Single-phase EB boundary conditions (Dirichlet / Neumann / Robin) | Not started |
| 6 | `petsc_jump` | — | Multifluid jump condition (dielectrics, surface charge) | Not started |
| 7 | `petsc_fieldsolver` | — | Full `FieldSolverAMG`: coefficients, regrid, options, E-field, tests, docs | Not started |
| 8 | `petsc_perf` | — | Matrix reuse, stencil caching, preconditioner tuning, scaling study | Not started |

### 1 — Infrastructure (`field_amg`, #601)

Retitle: *"PETSc infrastructure for an algebraic multigrid field solver"*.

- [x] Add a CI job that builds with `USE_PETSC=TRUE MPI=TRUE` (Ubuntu `petsc-dev`, `libopenmpi-dev`,
      `csh`), 2D and 3D, `discharge-lib` + `Exec/Examples/Electrostatics/PETSc`. Without this every
      PETSc file is compiled by nobody, which is how the two draft PRs drifted. → `Linux-PETSc` in
      `CI.yml`, a `CI-passed` dependency. Like the other build jobs it skips draft PRs.
- [x] Add `CH_USE_PETSC` to `PREDEFINED` in `Docs/doxygen.conf` so the PETSc classes are documented
      and the doxygen warnings check actually sees them. This exposed 37 warnings, most in PETSc code
      already on `main` (`PetscGrid`, `PetscAMRCell`, `getPetscGrid`, `initialize`); all fixed, the
      check is at 0. `PetscGrid::create/destroy` parameters renamed `x` → `a_x` on the way.
- [ ] `PetscStencil`: justify against `VoFStencil` in the PR description (it is keyed on a global
      PETSc row instead of a `VolIndex` + level/box, which is what `MatSetValues` needs). *PR
      description still to do; the same justification is now in the class Doxygen.*
- [x] `PetscStencil::add`: kept as a linear scan; the expected width (tens of entries, ~100 at most
      for a 3D least-squares neighbourhood plus CF interpolation) is documented on `add`. Doxygen for
      all members.
- [x] `FieldSolverAMG`: keep as a skeleton whose `solve` round-trips through `PetscGrid` (useful as a
      smoke test). Remove the dead `setupSolver` body that constructs `MFHelmholtzPetsc` — that moves
      to PR 7. Fix `#if CH_USE_PETSC` → `#ifdef`. *`setupSolver` was already empty after the sync;
      the `#if 1` around the smoke test is gone and the class Doxygen says it is a skeleton.*
- [x] Move the `MFHelmholtzPetsc` declaration out of this PR (it belongs with its implementation in
      PR 2/3), or keep it but with no `#warning`s. *Removed from `field_amg`; `helmholtz_petsc` keeps
      its own copy (the merge upward resolves the modify/delete in its favour).*
- [x] `Exec/Examples/Electrostatics/PETSc` compiles and runs (`mpirun -np 4`), 2D. 3D unchecked.
- [x] Drop `Exec/Examples/Electrostatics/ProfiledSurface/dev.inputs` (a scratch file; #602 already
      deletes it).

### 2 — Row translation (`helmholtz_petsc`, #602)

Retitle: *"Translate cut-cell stencils to PETSc rows"*.

The central design question for the whole stack is settled here. **Proposal: do not re-derive the
discretisation for PETSc.** `MFHelmholtzOp` and its BC classes already build, per level, every stencil
the operator needs as `VoFStencil`s over level cells *including ghost cells*:
`MFHelmholtzEBBC::getLeastSquaresBoundaryGradStencil`, `MFHelmholtzJumpBC::getGradPhiStencils`, the
interior face fluxes, and the domain BCs. What PETSc needs that GMG does not is the *composite* form:
every column must be a PETSc row, i.e. a valid (uncovered, non-ghost) cell on some level. So the job
of this PR is a translator:

- [ ] `PetscGrid` accessors from #602 (`getMFLevelGrids*`, `getAMRToPetsc*`, valid cut-cell iterators,
      `getIrregIVS`). Check the `LayoutData<IntVectSet>` + `mutable LayoutData<VoFIterator>` pair
      against `CLAUDE.md` §0 and drop whichever is redundant.
- [ ] A function `VoFStencil` (level `l`, box `din`) → `PetscStencil` that resolves each entry:
  - valid cell on level `l` → its row;
  - coarse-fine ghost cell on level `l` → substitute the interpolation stencil from
    `EBMultigridInterpolator::getInterpolationStencil*` (fine part on `l`, coarse part on `l-1`),
    recursively resolved;
  - cell covered by level `l+1` → only legal inside a CF flux, handled in PR 4; abort otherwise.
- [ ] Unit test: for random `phi` on a 2-level grid, `apply(PetscStencil, phi_petsc)` equals
      `apply(VoFStencil, phi_with_ghosts_filled_by_interpGhostMG)` to round-off, for every regular
      and cut cell.
- [ ] Remove the hand-rolled least-squares neighbourhood gathering in
      `MFHelmholtzPetsc::computeEBDirichletStencil` if the translator approach is accepted; otherwise
      finish it and record why reuse was rejected here.

### 3 — Regular operator and first solve (`petsc_regular`)

- [ ] Assemble `alpha*a*phi + beta*div(b grad phi)` on a single level, single phase, no EB, into an
      `MPIAIJ` matrix via `PetscGrid` rows; preallocate from stencil widths.
- [ ] Domain BCs (Dirichlet / Neumann / Robin, constant and function-valued) moved into the RHS,
      reusing `MFHelmholtzDomainBCFactory` stencils.
- [ ] `KSPSolve` with `-ksp_type cg|gmres -pc_type gamg|hypre` taken from options; residual
      reported in the same norm as `FieldSolverGMG`.
- [ ] Test: manufactured solution on a box, second-order convergence over three resolutions;
      `FieldSolverGMG` and PETSc agree to solver tolerance.

### 4 — Composite AMR operator (`petsc_amr`)

- [ ] Coarse cells adjacent to a CF interface: replace the coarse face flux by the area-weighted sum
      of fine face fluxes (refluxing folded into the matrix).
- [ ] Fine fluxes that touch CF ghosts go through the PR 2 translator.
- [ ] Covered coarse cells are not rows; their values are set from the solution by
      `conservativeAverage` after the solve.
- [ ] Test: 2- and 3-level static refinement, matches `FieldSolverGMG` to solver tolerance;
      PETSc residual of the GMG solution is at GMG's tolerance.

### 5 — Embedded boundaries, single phase (`petsc_eb`)

- [ ] Dirichlet / Neumann / Robin EB fluxes from the existing `MFHelmholtz*EBBC` stencils, translated.
- [ ] Handle the fallbacks GMG uses (order reduction, neighbourhood priority list) without
      duplicating them — expose them from the BC classes if needed.
- [ ] Cut cells next to CF interfaces (the case GMG handles with `isStencilValidCF`).
- [ ] Test: `Exec/Tests/Electrostatics/RodSphere` (single phase) agrees with GMG.

### 6 — Jump condition (`petsc_jump`)

- [ ] Two-phase rows: each multifluid cut cell has one row per phase; the jump condition
      `eps_g dphi_g/dn - eps_s dphi_s/dn = sigma` couples them. Reuse `MFHelmholtzJumpBC` stencils.
- [ ] Both `natural` and `saturation_charge` jump types.
- [ ] Test: `ProfileCircle`/`ProfileSquare`/`MechShaft` agree with GMG.

### 7 — `FieldSolverAMG` (`petsc_fieldsolver`)

- [ ] `setupSolver`, `setPermittivities`/`setSolverPermittivities`, `regrid`, `parseOptions` with the
      same option names as `FieldSolverGMG` (BCs, `jump_bc`, `kappa_source`, filtering, plot vars)
      plus `FieldSolverAMG.ksp_*` / `pc_*`.
- [ ] Rebuild the matrix on regrid; on time-dependent permittivity, update coefficients only.
- [ ] `Exec/Examples/Electrostatics/PETSc` and `ProfiledSurface` use it; add `.ini` regression
      entries gated on a PETSc build.
- [ ] Sphinx: a `FieldSolverAMG` section next to `FieldSolverGMG`, with build instructions for PETSc.

### 8 — Performance (`petsc_perf`)

- [ ] Cache EB/jump stencils (the `#warning` in the original `MFHelmholtzPetsc.H`) so a coefficient
      change does not redo least squares.
- [ ] `MatSetValues` batching, preallocation accuracy, `MAT_NEW_NONZERO_ALLOCATION_ERR` on.
- [ ] Weak/strong scaling against `FieldSolverGMG` on `MechShaft3d`.

## Open decisions

1. **Reuse vs. re-derive the discretisation** (PR 2). The plan assumes reuse of the GMG stencils via a
   translator. #602 had started re-deriving least squares inside `MFHelmholtzPetsc`.
2. **`PetscStencil` as a new type** vs. a `VoFStencil` plus a separate row lookup (PR 1).
3. **Does `MFHelmholtzPetsc` survive?** With the translator it may reduce to an assembler living in
   `FieldSolverAMG`, or a thin `MFHelmholtzPetscAssembler`.
4. **Preconditioner default**: GAMG (always in PETSc) or BoomerAMG (needs hypre).

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
