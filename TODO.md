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
  `format-input-files`, `check-literalincludes`, `doxygen-check`
- [x] CI: added `REUSE` and `Codespell` jobs; build jobs depend on them
- [x] `doxygen.conf`: `WARN_AS_ERROR = FAIL_ON_WARNINGS`, `WARN_NO_PARAMDOC = YES`
- [x] Sphinx: `-W --keep-going` in `SPHINXOPTS`
- [x] `Docs/check_literalincludes.py` validator added and wired into CI

## Still to do

### Doxygen warnings — Source/ directories (ordered by warning count)

**Per-file checklist** — when touching any `.H` or `.cpp` file, do ALL of the following:
1. Replace the legacy `/* chombo-discharge … */` block with the SPDX header:
   ```
   /*
    * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
    *
    * SPDX-License-Identifier: GPL-3.0-or-later
    */
   ```
2. Convert every `/*!` doc comment to `/**`.
3. Uppercase the header guard (e.g. `CD_FooBar_H` → `CD_FOOBAR_H`).
4. Add `@return`/`@param[in]`/`@param[out]`/`@param[in,out]` where missing.
5. Document all non-trivial protected/private member variables with at least `@brief`.
6. Wrap any local `using` / `typedef` aliases in `.cpp` files with `/// @cond DOXYGEN_SKIP` / `/// @endcond`.

Approximate warning counts (will decrease as files are processed):

- [ ] `Geometries/` — unknown warning count (not yet run through doxygen)
- [ ] `Source/AmrMesh/` — ~293 warnings
- [ ] `Source/Elliptic/` — ~244 warnings
- [ ] `Source/Utilities/` — ~202 warnings
- [ ] `Source/Particle/` — ~127 warnings
- [ ] `Source/ConvectionDiffusionReaction/` — ~112 warnings (partially done)
- [ ] `Source/RadiativeTransfer/` — ~100 warnings (partially done)
- [ ] `Source/ImplicitFunctions/` — ~100 warnings
- [ ] `Source/Geometry/` — ~82 warnings (partially done)
- [ ] `Source/ItoDiffusion/` — ~76 warnings
- [ ] `Source/Electrostatics/` — ~64 warnings
- [ ] `Source/Multifluid/` — ~53 warnings (partially done)
- [ ] `Source/KineticMonteCarlo/` — ~45 warnings (partially done)
- [ ] `Source/Driver/` — ~24 warnings (partially done)
- [ ] `Source/CellTagger/` — ~24 warnings (partially done)
- [ ] `Source/TracerParticles/` — ~18 warnings (partially done)
- [ ] `Source/SurfaceODESolver/` — ~17 warnings (partially done)
- [ ] `Source/MeshODESolver/` — ~9 warnings (partially done)

### SPDX headers — remaining files

Files under `Source/` and `Geometries/` that still carry the old `/* chombo-discharge / Copyright © */`
block. Replace with the standard SPDX block. The `REUSE.toml` catch-all already makes `reuse lint`
pass, but inline headers are the long-term goal so the catch-all can eventually be removed.

### clang-tidy warnings — files to fix

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

### clang-tidy CI

The `cmake-generate-compile-commands` hook is listed in `ci: skip` but not yet defined.
The project uses GNU makefiles, not CMake, so a `compile_commands.json` must be generated
via `bear` or an equivalent wrapper. This is a separate task.

### Cleanup
After completing the above checklists, warn the user about various stubs that are still present in this branch.
*  TODO.md should not be a part of the PR, and clang-tidy must be integrated into the CI pipeline.

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

### Step 3 — Manual fixes (~300 warnings, pending user review)

#### 3a. `misc-unused-parameters` (103 warnings)
Comment out unused parameter names in definitions (Chombo convention: `const int /*a_level*/`).
Heavy files: `Physics/CdrPlasma/` (~48), `Source/Electrostatics/` (~23), `Source/Elliptic/` (8).

- [ ] `Physics/CdrPlasma/` (~48)
- [ ] `Source/Electrostatics/` (~23)
- [ ] `Source/Elliptic/` (8)
- [ ] `Physics/Electrostatics/` (7), `Physics/ItoKMC/` (6), `Exec/` (4), `Source/AmrMesh/` (3), others (4)

#### 3b. `cppcoreguidelines-special-member-functions` (59 warnings)
For each class that declares a custom dtor but not all five special members, add `= delete` (polymorphic)
or `= default` (value). ~30 header files including:
- `Geometries/CoaxialCable/`, `DoubleRod/`, `Rod/`, `Vessel/`, `WireWire/`
- `Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.H`
- `Physics/CdrPlasma/CD_CdrPlasmaFieldTagger.H`, `…StreamerTagger.H`
- `Physics/ItoKMC/CD_ItoKMCFieldTagger.H`, `…GodunovStepper.H`
- `Source/AmrMesh/CD_EBAMRData.H`, `CD_EBLeastSquaresMultigridInterpolator.H`, and more

- [ ] `Geometries/` classes (5 headers)
- [ ] `Physics/` classes (~8 headers)
- [ ] `Source/` classes (~17 headers)

#### 3c. `bugprone-narrowing-conversions` (56 warnings)
Add `static_cast<TargetType>(expr)` for `size_t`→`int` and `long long`→`double`.
- [ ] `Physics/ItoKMC/CD_ItoKMCStepperImplem.H` (9 warnings)
- [ ] `Source/Driver/CD_Driver.cpp` (10 warnings)
- [ ] Other files (1–3 each)

#### 3d. `clang-diagnostic-unused-variable` (30+ warnings)
Remove or `(void)var;` unused variables.
- [x] `Source/AmrMesh/CD_EBCoarseToFineInterp.cpp` — `domainFine`, `ebisBoxFine` (×3 functions), `ebisBoxCoar`, `coarBox`, `volFactor`
- [x] `Source/AmrMesh/CD_EBGhostCellInterpolator.cpp` — `domainCoar` (in `define`), `fineDomainBox` (in `interpolateIrregular`); note: `domainFine` and `coarDomainBox` are used and were kept
- [x] `Source/AmrMesh/CD_EBFluxRedistribution.cpp` — `cellBox`, `domainCoar`, `domain`
- [x] `Source/AmrMesh/CD_LevelTiles.cpp` — `numRanks`
- [x] `Source/AmrMesh/CD_EBReflux.cpp` — `domainCoFi`, `boxFine`, `graphFine`, `dxFine`, etc.
- [x] `Source/AmrMesh/CD_CellCentroidInterpolation.cpp` — `domainBox` (×3)
- [x] `Source/AmrMesh/CD_PhaseRealm.cpp` — `comps` (×2)
- [x] `Source/AmrMesh/CD_EBLeastSquaresMultigridInterpolator.cpp` — `nboxCoar`, `domainCoar`, `fineBox`
- [x] `Source/AmrMesh/CD_EBGradient.cpp` — `cellBox`, `irregularKernel`, `domainFine`, `levelStencils`
- [x] `Source/AmrMesh/CD_EBAMRParticleMesh.cpp` — `domainFine`, `domainCoar`, `numBoxesCoar`
- [x] `Source/Driver/CD_Driver.cpp` — `probLo`
- [x] `Source/Utilities/CD_PolyUtils.cpp` — `x`, `fx`
- [x] `Source/Utilities/CD_DischargeIO.cpp` — `numCompTotal`
- [x] `Source/Utilities/CD_DataOps.cpp` — `irregularKernel`, `dbl` (×6)
- [x] `Source/Elliptic/CD_MFHelmholtzOp.cpp` — `isIrregular`, `homogeneousCFBC` (×2)
- [x] `Source/Elliptic/CD_MFHelmholtzOpFactory.cpp` — `found`
- [x] `Source/Elliptic/CD_EBHelmholtzOp.cpp` — `tol`
- [x] `Source/Elliptic/CD_MFHelmholtzRobinEBBC.cpp` — `box`
- [x] `Source/Elliptic/CD_MFHelmholtzJumpBC.cpp` — `box` (×3), `ebisBoxPhase0`, `ebisBoxPhase1`, `avgStencilsPhase0`, `avgStencilsPhase1`, `denomFactorPhase1`
- [x] `Source/RadiativeTransfer/CD_McPhoto.cpp` — `factor` (dead store), `box`
- [x] `Source/RadiativeTransfer/CD_EddingtonSP1.cpp` — `helmAcoBox`, `helmBcoBox`, `finestLevel`
- [x] `Source/ConvectionDiffusionReaction/CD_CdrGodunov.cpp` — dangling `:` in constructor
- [x] `Source/ConvectionDiffusionReaction/CD_CdrSolver.cpp` — `cellBox` (×3)
- [x] `Source/ItoDiffusion/CD_ItoSolver.cpp` — `dx`, `origin`, `box`, `probLo` (×3 functions)
- [x] `Source/Particle/CD_EBCoarseFineParticleMesh.cpp` — `nboxCoFi`
- [x] `Source/ImplicitFunctions/CD_NeedleIF.cpp` — `tipLength`
- [x] `Geometries/GECReferenceCell/CD_GECReferenceCell.cpp` — `h`
- [x] `Geometries/RodNeedleDisk/CD_RodNeedleDisk.cpp` — `tipLength`
- [x] `Physics/AdvectionDiffusion/CD_AdvectionDiffusionStepper.cpp` — `fluxFunc` (×2)
- [x] `Physics/BrownianWalker/CD_BrownianWalkerStepper.cpp` — `nComp`
- [x] `Physics/CdrPlasma/CD_CdrPlasmaStepper.cpp` — `normal`, `domain`, `electricFieldCell`, `ebgraph`, `cellBox`
- [x] `Physics/CdrPlasma/Timesteppers/CdrPlasmaImExSdcStepper/CD_CdrPlasmaImExSdcStepper.cpp` — `sigma_p`, `safety`, `max_dt_cfl`, `t`
- [x] `Physics/ItoKMC/PlasmaModels/ItoKMCJSON/CD_ItoKMCJSON.cpp` — `plasmaProducts`, `N`
- [ ] `Physics/ItoKMC/CD_ItoKMCStepperImplem.H` (8 warnings — header-only template file)
- [ ] `Physics/ItoKMC/TimeSteppers/ItoKMCGodunovStepper/CD_ItoKMCGodunovStepperImplem.H` (1)

#### 3e. `readability-inconsistent-declaration-parameter-name` (25 warnings)
Sync parameter names in `.H` declarations to match the `.cpp`/`Implem.H` definitions.
- [ ] `Source/Particle/CD_ParticleContainer.H`
- [ ] `Source/Electrostatics/CD_ElectrostaticDomainBc.H`, `CD_FieldSolverGMG.H`
- [ ] `Source/Geometry/CD_ComputationalGeometry.H`
- [ ] `Source/Utilities/CD_Location.H`, `CD_DischargeIO.H`
- [ ] `Source/AmrMesh/CD_LevelTiles.H`
- [ ] `Source/Multifluid/CD_MFLevelGrid.H`
- [ ] `Source/RadiativeTransfer/CD_RtLayout.H`
- [ ] `Source/Elliptic/CD_MFHelmholtzOp.H`
- [ ] `Source/ImplicitFunctions/CD_ProfilePlaneIF.H`
- [ ] `Physics/ItoKMC/…/CD_ItoKMCGodunovStepper.H`

#### 3f. `clang-analyzer-deadcode.DeadStores` (10+ warnings)
Remove or use dead-stored variables.
- [x] `Source/Driver/CD_Driver.cpp` — `numThreads` (wrapped in `#ifdef _OPENMP`), `lastOutputTime`, `probLo`
- [x] `Source/Utilities/CD_PolyUtils.cpp` — `c` (dead store before unreachable), `x`, `fx`
- [x] `Source/RadiativeTransfer/CD_McPhoto.cpp` — `factor = 0.0` before `MayDay::Error`
- [ ] `Physics/ItoKMC/CD_ItoKMCStepperImplem.H` — various (header-only template file)

#### 3g. `performance-unnecessary-value-param` (417 warnings)
Pass parameters by `const&` instead of by value where the parameter is not mutated and copy is not
intentional. Requires changing both the `.H` declaration and the `.cpp`/`Implem.H` definition.
Heavy files: `Physics/CdrPlasma/` (~130), `Physics/ItoKMC/` (~90), `Source/Elliptic/` (~50).
Not auto-fixed in Step 2 because it changes function signatures.
- [ ] `Physics/CdrPlasma/` (~130 warnings)
- [ ] `Physics/ItoKMC/` (~90 warnings)
- [ ] `Source/Elliptic/` (~50 warnings)
- [ ] Remaining files (~147 warnings spread across Source/ and other Physics/)

#### 3h. `readability-convert-member-functions-to-static` (8 warnings)
Verify function doesn't use `this`, then add `static`.
- [ ] `Source/Geometry/CD_ScanShopImplem.H:61` — `getSortedBoxesAndTypes`
- [ ] Other files

#### 3i. `clang-diagnostic-overloaded-virtual` (2 warnings)
Add `using Base::method;` to bring hidden base overloads into scope.
- [ ] `Physics/CdrPlasma/Timesteppers/CdrPlasmaGodunovStepper/CD_CdrPlasmaGodunovStepper.H`
  - `computeCdrDomainFluxes` (line 389)
  - `computeCdrDriftVelocities` (line 459)

#### 3j. Small-count fixes
- [ ] `Source/Particle/CD_ParticleContainer.H:96` — add `noexcept` to move assignment
- [ ] `Source/ItoDiffusion/CD_ItoSolver.H:89` — add `noexcept` to move assignment
- [ ] `Source/ItoDiffusion/CD_ItoParticleImplem.H:41` — copy ctor calls wrong base
- [ ] `Physics/ItoKMC/CD_ItoKMCSurfaceReactionSet.H:73` — `const int max` member → non-const
- [ ] `Physics/ItoKMC/CD_ItoKMCStepperImplem.H:403` — remove unused `std::string str`

#### 3k. Real bugs from clang-analyzer (priority — investigate before committing)

These emerged from the full run after `HeaderFilterRegex` was fixed to cover all four directories.
Run `./run_clang_tidy.sh 2>&1 | grep -E "clang-analyzer|clang-diagnostic-(sometimes|uninitialized|vla)"` to get per-file details.

- [x] `Source/ConvectionDiffusionReaction/CD_CdrMultigrid.cpp` — missing `break` in GMRES case of
  bottom-solver `switch`; without it, GMRES fell through to `default: MayDay::Error(...)` and
  aborted the program whenever GMRES was the selected solver. Fixed by adding `break;`.

- [ ] `clang-analyzer-core.NullDereference` (~49 warnings): potential null pointer dereferences.
  Likely concentrated in `Physics/ItoKMC/` and `Physics/CdrPlasma/` — needs per-file investigation.

- [ ] `clang-analyzer-core.StackAddressEscape` (~3 warnings): function returns pointer/reference to
  local variable — undefined behaviour on any caller that uses the returned value.

- [ ] `clang-diagnostic-sometimes-uninitialized` (~3 warnings): variable used uninitialised on some
  code paths. Known instance: `Source/Geometry/CD_ScanShop.cpp:158` — `whichLevel` uninitialized
  when for-loop exits with condition false, then used at line 165.

- [ ] `clang-analyzer-core.uninitialized.Assign` (~2 warnings): value assigned to field/variable
  before it is initialised.

- [ ] `clang-analyzer-cplusplus.NewDelete` (~3 warnings): memory management errors (leak or
  use-after-free pattern detected by the analyser).

- [ ] `clang-analyzer-unix.MismatchedDeallocator` (1 warning): `malloc`/`new` paired with wrong
  deallocation (`free`/`delete` mismatch).

- [ ] `clang-diagnostic-vla-cxx-extension` (~3 warnings): variable-length arrays — a GCC extension,
  not standard C++14. Replace with `std::vector` or fixed-size array.

- [ ] `clang-diagnostic-uninitialized-const-reference` (~3 warnings): uninitialized object passed
  as `const T&` — compiler can generate a temporary copy but the value is garbage.

- [ ] `readability-suspicious-call-argument` (2 warnings): `Source/Elliptic/CD_MFHelmholtzOp.cpp`
  lines 1099 and 1183 — `a_phiFine`/`a_phi`/`a_phiCoar` arguments may be swapped. Verify order
  against the called function's signature before concluding it is a real swap.

- [ ] `bugprone-suspicious-include` (1 warning): an `#include` directive that looks like it may
  include a `.cpp` file by mistake.
