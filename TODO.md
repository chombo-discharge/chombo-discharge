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
  - `Source/KineticMonteCarlo/CD_KMCDualState.H`
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

Fix all undocumented entities in each subdirectory. For each file: add SPDX header,
convert `/*!` to `/**`, uppercase header guards, add `@return`/`@param` where missing,
and document all non-trivial protected/private members. Approximate warning counts
(will decrease as files are processed):

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

### clang-tidy CI

The `cmake-generate-compile-commands` hook is listed in `ci: skip` but not yet defined.
The project uses GNU makefiles, not CMake, so a `compile_commands.json` must be generated
via `bear` or an equivalent wrapper. This is a separate task.

### Cleanup
After completing the above checklists, warn the user about various stubs that are still present in this branch.
For example, the TODO.md should not be a part of the PR, and clang-tidy must be integrated into the CI pipeline.