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
- [ ] `Source/KineticMonteCarlo/CD_KMCSolver.H`
- [ ] `Source/KineticMonteCarlo/CD_KMCSolverImplem.H`
- [ ] `Source/MeshODESolver/CD_MeshODESolver.H`
- [x] `Source/MeshODESolver/CD_MeshODESolverImplem.H`
- [x] `Source/Particle/CD_EBAMRParticleMesh.H`
- [x] `Source/Particle/CD_EBAMRParticleMeshImplem.H`
- [x] `Source/Particle/CD_GenericParticle.H`
- [x] `Source/Particle/CD_GenericParticleImplem.H`
- [x] `Source/Particle/CD_ParticleContainer.H`
- [ ] `Source/Particle/CD_ParticleContainerImplem.H`
- [x] `Source/Particle/CD_PointParticle.H`
- [ ] `Source/Particle/CD_PointParticleImplem.H`
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
* Check if ANY /*! block appears in the code, and fix the files that contain such blocks.
* We must reset the H5 files to not use the cycle -- this causes a bug int he VisIt reader