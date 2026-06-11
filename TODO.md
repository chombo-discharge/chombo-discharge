# SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

# Ongoing documentation and compliance PR — task list

This file tracks the remaining work for the documentation/REUSE/linting overhaul.

## Completed

- [x] Full REUSE compliance — all 1327 files covered (inline SPDX headers or `REUSE.toml` bulk annotations)
- [x] `codespell` clean on `Source/`, `Physics/`, `Exec/`, `Docs/`
- [x] `Tesselation` → `Tessellation` class rename everywhere (files, content, `.options`)
- [x] All `@file` directive mismatches fixed (18 files)
- [x] All `@param[inout]` → `@param[in,out]` fixed (84 files, batch sed)
- [x] `@paramo` → `@param` typos fixed (6 files)
- [x] `@detail` → `@details` typos fixed (3 files)
- [x] `\int` in Doxygen comments wrapped in `@f$…@f$` math delimiters
- [x] Local `using` type aliases in `Geometries/` `.cpp` files wrapped in `@cond DOXYGEN_SKIP`
- [x] Missing `@return` and wrong `@param` names fixed in:
  - `Source/CellTagger/CD_CellTagger.H`
  - `Source/Driver/CD_Driver.H`, `CD_TimeStepper.H`, `CD_Initialize.H`
  - `Source/TracerParticles/CD_TracerParticle.H`, `CD_TracerParticleSolver.H`
  - `Source/MeshODESolver/CD_MeshODESolver.H`
  - `Source/SurfaceODESolver/CD_SurfaceODESolver.H`
  - `Source/ConvectionDiffusionReaction/CD_CdrSolver.H`, `CD_CdrGodunov.H`, `CD_CdrCTU.H`,
    `CD_CdrMultigrid.H`, `CD_CdrIterator.H`, `CD_CdrLayout.H`, `CD_CdrDomainBC.H`
  - `Source/KineticMonteCarlo/CD_KMCDualState.H`, `CD_KMCDualStateReaction.H`, `CD_KMCSingleState.H`, `CD_KMCSingleStateReaction.H`, `CD_KMCSolver.H` (full reformat: SPDX, `/**`, uppercase guards, all warnings)
  - `Source/Geometry/CD_ComputationalGeometry.H`
  - `Source/RadiativeTransfer/CD_RtSolver.H` (partial)
  - `Source/Multifluid/CD_MultiFluidIndexSpace.H`
  - Various `Geometries/` headers (`Rod`, `RodDielectric`, `RodPlaneProfile`, `MechanicalShaft`,
    `WireWire`, `RandomInterface`)
- [x] Doxygen generates HTML (exits non-zero only due to undocumented `Source/` entities)
- [x] Pre-commit hooks: `clang-format`, `clang-tidy` (CI skip), `reuse`, `codespell`,
  `format-input-files`, `check-literalincludes`, `doxygen-check`, `cppcheck`
- [x] CI: added `REUSE`, `Codespell`, `Clang-tidy`, and `Cppcheck` jobs; build jobs depend on them
- [x] `doxygen.conf`: `WARN_AS_ERROR = FAIL_ON_WARNINGS`, `WARN_NO_PARAMDOC = YES`
- [x] Sphinx: `-W --keep-going` in `SPHINXOPTS`
- [x] `Docs/check_literalincludes.py` validator added and wired into CI
- [x] `clang-tidy-cache` integrated into `run_clang_tidy.sh` and the CI `Clang-tidy` job
- [x] `cppcheck` added to pre-commit (staged `.cpp` files, no compilation database required) and to CI (with Bear compilation database); two cppcheck findings fixed: dead self-assignment in `CD_CdrPlasmaImExSdcStepper.cpp`, `unknownMacro` suppression for Chombo `D_TERM6` in `CD_LoadBalancing.cpp`
- [x] All doxygen `@param` name mismatches fixed across 46 header files — `doxygen-check` pre-commit hook passes with zero warnings
- [x] `bugprone-narrowing-conversions`: all 157 warnings fixed with `static_cast` (Python script + manual fixes for compound expressions)

## Still to do

### ~~Doxygen warnings~~ — DONE

`pre-commit run doxygen-check --all-files` passes with zero warnings.
All `@param` name mismatches, missing `@return` tags, and undocumented defaulted special
members have been fixed across the codebase.  The per-file checklist below is retained for
reference when adding new code:

1. Convert every `/*!` doc comment to `/**`.
2. Add `@return`/`@param[in]`/`@param[out]`/`@param[in,out]` where missing.
3. Document all non-trivial protected/private member variables with at least `@brief`.

### SPDX headers — remaining files

Files under `Source/` and `Geometries/` that still carry the old `/* chombo-discharge / Copyright © */`
block. Replace with the standard SPDX block. The `REUSE.toml` catch-all already makes `reuse lint`
pass, but inline headers are the long-term goal so the catch-all can eventually be removed.

### clang-tidy warnings — files to fix

All 46 per-file boxes below are ticked; the bulk-fix plan in the section that follows covers
the remaining work (3a–3k).

**Per-file checklist** — when addressing clang-tidy warnings in any file, also do ALL of the following:
1. Ensure correct Doxygen formatting per the checklist in the "Doxygen warnings" section above.
2. Check whether the file is referenced by a `.. literalinclude::` directive in any RST file
   (run `grep -r "literalinclude" Docs/Sphinx/source/ | grep <filename>`).
   If a reference exists, verify that the RST file still points to the same code block after your changes.

**Verification steps** — after all files are done, run:
1. `pre-commit run clang-tidy --all-files` — must produce no warnings
2. `pre-commit run doxygen-check --all-files` — must produce no warnings
3. For each changed file, check for RST literalincludes:
   `grep -r "literalinclude" Docs/Sphinx/source/ | grep <filename>`
   and verify that each matching RST block still references valid code after edits.

- [x] `Exec/Convergence/AdvectionDiffusion/C2/main.cpp`
- [x] `Exec/Convergence/KineticMonteCarlo/C1/main.cpp`
- [x] `Exec/Convergence/RadiativeTransfer/C2/main.cpp`
- [x] `Exec/Tests/Utilities/LookupTable/main.cpp`
- [x] `Geometries/NoisePlane/CD_NoisePlane.H`
- [x] `Geometries/SpherePlane/CD_SpherePlane.H`
- [x] `Physics/AdvectionDiffusion/CD_AdvectionDiffusionTagger.cpp`
- [x] `Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.H`
- [x] `Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStorage.H`
- [x] `Physics/ItoKMC/CD_ItoKMCSurfaceReactionsImplem.H`
- [x] `Source/AmrMesh/CD_EBAMRData.H`
- [x] `Source/AmrMesh/CD_EBAMRDataImplem.H`
- [x] `Source/AmrMesh/CD_EBCentroidInterpolation.cpp`
- [x] `Source/AmrMesh/CD_LinearStencil.cpp`
- [x] `Source/Electrostatics/CD_MFHelmholtzElectrostaticEBBC.H`
- [x] `Source/Electrostatics/CD_MFHelmholtzElectrostaticEBBCImplem.H`
- [x] `Source/Elliptic/CD_EBHelmholtzOpFactory.H`
- [x] `Source/Elliptic/CD_MFHelmholtzNeumannEBBC.cpp`
- [x] `Source/ImplicitFunctions/CD_HyperboloidIF.cpp`
- [x] `Source/ImplicitFunctions/CD_SphereArray.H`
- [x] `Source/KineticMonteCarlo/CD_KMCSingleStateReaction.H`
- [x] `Source/KineticMonteCarlo/CD_KMCSingleStateReactionImplem.H`
- [x] `Source/KineticMonteCarlo/CD_KMCSolver.H`
- [x] `Source/KineticMonteCarlo/CD_KMCSolverImplem.H`
- [x] `Source/MeshODESolver/CD_MeshODESolver.H`
- [x] `Source/MeshODESolver/CD_MeshODESolverImplem.H`
- [x] `Source/Particle/CD_EBAMRParticleMesh.H`
- [x] `Source/Particle/CD_EBAMRParticleMeshImplem.H`
- [x] `Source/Particle/CD_GenericParticle.H`
- [x] `Source/Particle/CD_GenericParticleImplem.H`
- [x] `Source/Particle/CD_ParticleContainer.H`
- [x] `Source/Particle/CD_ParticleContainerImplem.H`
- [x] `Source/Particle/CD_PointParticle.H`
- [x] `Source/Particle/CD_PointParticleImplem.H`
- [x] `Source/Utilities/CD_DischargeIO.H`
- [x] `Source/Utilities/CD_DischargeIOImplem.H`
- [x] `Source/Utilities/CD_LeastSquares.H`
- [x] `Source/Utilities/CD_LeastSquaresImplem.H`
- [x] `Source/Utilities/CD_Location.H`
- [x] `Source/Utilities/CD_LocationImplem.H`
- [x] `Source/Utilities/CD_LookupTable.H`
- [x] `Source/Utilities/CD_LookupTable1D.H`
- [x] `Source/Utilities/CD_ParallelOps.H`

### Cleanup
After completing the above checklists, `TODO.md` should be removed before merging (`git rm TODO.md`).

---

## clang-tidy bulk-fix plan

Generated from `clang-tidy.log` (debug-mode build, ~1500 warnings in Source/Geometries/Physics/Exec).
Three-step plan: disable noisy checks, auto-fix mechanically, then fix manually.

### Step 1 — Disable 5 checks in `.clang-tidy` (~175 warnings)

Add to the disabled list:

| Check | Count | Reason |
|-------|-------|--------|
| `misc-confusable-identifiers` | 128 | `I` (current) and `l` (length) are standard physics symbols |
| `clang-analyzer-optin.cplusplus.VirtualCall` | 30 | Chombo `define()` from ctor is an embedded design pattern |
| `misc-use-anonymous-namespace` | 10 | Static functions in headers cannot go in anon namespaces (ODR) |
| `clang-diagnostic-c++17-attribute-extensions` | 4 | `[[nodiscard]]` is intentional; works fine in C++14 mode |
| `bugprone-unhandled-exception-at-new` | 3 | Chombo/scientific code never exception-handles `new` |

- [x] Step 1 done

### Step 2 — Auto-fix via `run-clang-tidy --fix` (~1046 warnings, plus additional performance warnings exposed by fixed HeaderFilterRegex)

Run from `$DISCHARGE_HOME` (where `compile_commands.json` lives). One check at a time to avoid fix conflicts:

```bash
for CHECK in \
    "readability-avoid-const-params-in-decls" \
    "readability-braces-around-statements" \
    "cppcoreguidelines-explicit-virtual-functions,modernize-use-override" \
    "cppcoreguidelines-prefer-member-initializer" \
    "modernize-use-emplace" \
    "modernize-use-equals-default" \
    "modernize-use-auto" \
    "readability-container-size-empty" \
    "readability-qualified-auto" \
    "readability-redundant-string-init" \
    "readability-redundant-string-cstr" \
    "readability-redundant-casting" \
    "readability-redundant-member-init" \
    "readability-redundant-control-flow" \
    "readability-redundant-inline-specifier" \
    "modernize-deprecated-headers" \
    "modernize-use-nullptr" \
    "modernize-make-shared" \
    "modernize-loop-convert"; do
  run-clang-tidy -p . -j$(nproc) -fix \
      -checks="-*,$CHECK" \
      'Source/|Exec/|Physics/|Geometries/'
done
```

After the loop, reformat all modified files:
```bash
find Source Physics Geometries Exec \( -name "*.H" -o -name "*.cpp" \) \
    -exec clang-format -i {} +
```

- [x] Step 2 done (also fixed: performance-for-range-copy, performance-faster-string-find, performance-trivially-destructible, performance-inefficient-string-concatenation, performance-inefficient-vector-operation, readability-string-compare, readability-static-accessed-through-instance, readability-duplicate-include, readability-const-return-type)
- Note: `performance-avoid-endl` was auto-applied then reverted; disabled in `.clang-tidy` (bare `endl` is correct for this codebase)

### Step 3 — Manual fixes — DONE

`./run_clang_tidy.sh` produces **zero warnings** in `Source/`, `Physics/`, `Geometries/`,
and `Exec/`.  The 6 remaining diagnostics are all in `Submodules/Chombo-3.3` (outside
project scope: 2× `clang-analyzer-core.StackAddressEscape` in `LevelData.H`, 4×
`clang-analyzer-cplusplus.NewDeleteLeaks` in `SPMDI.H`).

All sub-items (3a–3k) are resolved:
- **3a** `misc-unused-parameters` — unused parameter names commented out
- **3b** `cppcoreguidelines-special-member-functions` — `= delete`/`= default` added
- **3c** `bugprone-narrowing-conversions` — 157 `static_cast` fixes
- **3d** `clang-diagnostic-unused-variable` — all dead variables removed
- **3e** `readability-inconsistent-declaration-parameter-name` — all param names synced
- **3f** `clang-analyzer-deadcode.DeadStores` — dead stores removed
- **3g** `performance-unnecessary-value-param` — pass-by-value parameters converted to `const&`
- **3h** `readability-convert-member-functions-to-static` — static functions marked `static`
- **3i** `clang-diagnostic-overloaded-virtual` — hidden base overloads brought into scope
- **3j** Small-count fixes — `noexcept`, wrong-base copy-ctor, `const` member, unused variable
- **3k** Real clang-analyzer bugs — all fixed (including the missing `break` in `CD_CdrMultigrid.cpp`)

---

## PR finalization checklist

Complete every item below **before** merging. Run them in order — later steps depend on
earlier ones passing.

### 1. Remove non-PR files

`compile_commands.json` is already gitignored.

- [ ] `git rm TODO.md` and commit the removal (`run_clang_tidy.sh` is now a first-class CI script — keep it)
- [ ] Delete `clang-tidy.log` from the working tree (untracked; add to `.gitignore` if it
  keeps reappearing)

### 2. CI dry-run — mirrors each job in `.github/workflows/CI.yml`

Run locally in the same order the CI jobs execute.

**Formatting** (`clang-format` job)
- [ ] `pre-commit run clang-format --all-files` — must produce an empty diff

**REUSE** (`REUSE` job)
- [ ] `reuse lint` — must exit 0

**Codespell** (`Codespell` job)
- [ ] `codespell Source/ Physics/ Exec/ Docs/ Geometries/` — must exit 0

**Build** (`Linux-GNU` job — 2-D, no MPI, GCC)
- [ ] `make -j$(nproc)` — must compile without errors

### 3. Documentation (`Build-documentation` job)

- [x] **Doxygen**: `doxygen Docs/doxygen.conf` — exits 0 with zero warnings
- [ ] **Sphinx**: `cd Docs/Sphinx && make html SPHINXOPTS="-W --keep-going"` — must exit 0
- [ ] **literalinclude validator**: `python3 Docs/check_literalincludes.py` — all paths
  resolve

### 4. Test executables

- [ ] Build at least one representative test:
  ```bash
  cd Exec/Tests/Electrostatics/RodSphere && make -j$(nproc) DIM=2
  ```
- [ ] Run with the regression inputs and confirm the simulation completes without assertion
  failures or NaNs:
  ```bash
  ./main.*.ex regression2d.inputs
  ```

### 5. Benchmark / regression inputs

- [ ] `git diff HEAD -- '*.inputs' '*.options'` — confirm no accidental modifications to
  any `.inputs` or `.options` file; every change must be intentional and reviewed

### 6. Full pre-commit sweep

- [ ] `pre-commit run --all-files` — every hook must pass (clang-format, reuse, codespell,
  cppcheck, format-input-files, doxygen-check all run automatically)
  - Note: the `clang-tidy` hook requires `compile_commands.json`; run
    `./run_clang_tidy.sh 2>&1 | grep -c "warning:"` manually to verify zero warnings in
    `Source/`, `Geometries/`, `Physics/`, `Exec/`
