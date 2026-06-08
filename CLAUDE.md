# chombo-discharge — Claude Code guidelines

## 1. Building the code

chombo-discharge uses the Chombo GNU makefile build system. Two environment variables must be
set before any build:

```bash
export DISCHARGE_HOME=/path/to/chombo-discharge
export CHOMBO_HOME=/path/to/Chombo/lib
```

A per-machine `Make.defs.local` file in `Lib/Local/` controls compiler, MPI, HDF5, and
optimisation flags. Copy the template and edit it for the local machine:

```bash
cp Lib/Local/Make.defs.local.template Lib/Local/Make.defs.<hostname>
```

### Building the core library

```bash
# Build Chombo + chombo-discharge Source + Geometries + all Physics modules
make -j$(nproc)

# Build only the core discharge library (Source + Geometries, no Physics)
make discharge-lib -j$(nproc)

# Build a single Physics module, e.g. Electrostatics
make electrostatics -j$(nproc)
```

Common make flags (append to any of the above):

| Flag | Values | Notes |
|------|--------|-------|
| `DIM` | `2`, `3` | Spatial dimension (default 2) |
| `DEBUG` | `TRUE`, `FALSE` | Enables assertions and debug symbols |
| `OPT` | `HIGH`, `FALSE` | Optimisation level |
| `MPI` | `TRUE`, `FALSE` | Build with MPI |
| `USE_HDF` | `TRUE`, `FALSE` | Build with HDF5 I/O |

Example: 3-D MPI optimised build:

```bash
make -j$(nproc) DIM=3 MPI=TRUE OPT=HIGH
```

### Building a test executable

Each test lives under `Exec/Tests/<Module>/<TestName>/` and has its own `GNUmakefile`.
Building a test also triggers the required library builds via `make dependencies`:

```bash
cd Exec/Tests/Electrostatics/RodSphere
make -j$(nproc) DIM=2
```

The compiled binary is named `main.<config>.ex` where `<config>` encodes the build flags
(e.g. `main.Linux.64.g++.gfortran.DEBUG.OPT.MPI.ex`).

### Template-only code

Several modules (`MeshODESolver`, `SurfaceODESolver`, `TracerParticles`, …) are
header-only (`.H` + `*Implem.H`). They are compiled implicitly when a test that uses them
is built. There is no standalone library build target for these modules.

### Cleaning

```bash
make libclean    # Remove discharge library objects/archives
make allclean    # Remove all discharge objects/archives
make pristine    # Full clean including Chombo
```

---

## 2. Testing the code

### Running a test

After building, run the executable from its directory (it reads a `.inputs` file):

```bash
cd Exec/Tests/Electrostatics/RodSphere
./main.Linux.64.g++.gfortran.DEBUG.OPT.MPI.ex regression2d.inputs
```

With MPI:

```bash
mpirun -np 4 ./main.Linux.64.g++.gfortran.DEBUG.OPT.MPI.ex regression2d.inputs
```

Output goes to `pout.<rank>` files and HDF5 plot/checkpoint files in `plt/`, `chk/`, etc.

### Regression tests

Each test directory contains `regression2d.inputs` and/or `regression3d.inputs`. These are
the canonical reference inputs. Running with these inputs and checking that the simulation
completes without assertion failures or NaNs is the baseline regression check.

There is no automated regression-comparison framework in the repository yet; correctness is
verified by visual inspection of the output or by comparison against reference solutions in
the Sphinx documentation.

### CI tests

The GitHub Actions CI (`CI.yml`) builds with GNU and oneAPI compilers in 2-D without MPI.
It does not run the executables; it only checks that compilation succeeds.

---

## 3. Code style

### Naming conventions

- Classes: `PascalCase` (`AmrMesh`, `CdrSolver`).
- Member variables: `m_camelCase` (`m_phi`, `m_amrMesh`).
- Local variables and parameters: `a_camelCase` for function parameters, `camelCase` for locals.
- Constants and enumerators: as in the surrounding Chombo/C++ style.

### File naming

- Header files use `.H` (capital H), not `.h` or `.hpp`.
- Implementation files use `.cpp`.
- Template implementations that are `#include`-d from a header use `*Implem.H`.
- All files are prefixed with `CD_` (e.g. `CD_AmrMesh.H`, `CD_AmrMesh.cpp`).

### Header guards

Header guards must be fully uppercased. A file `CD_FooBar.H` uses:

```cpp
#ifndef CD_FOOBAR_H
#define CD_FOOBAR_H
// ...
#endif
```

### Include ordering

Group and order includes as follows; precede each non-empty group with a comment label:

```cpp
// Std includes
#include <iostream>

// Chombo includes
#include <ParmParse.H>
#include <EBISBox.H>

// Our includes
#include <CD_Location.H>
#include <CD_NamespaceHeader.H>
```

Omit a group label if that group has no entries for the file in question.

### Namespace

All chombo-discharge code lives inside the `CH_SPACEDIM`-aware Chombo namespace, opened
and closed with:

```cpp
#include <CD_NamespaceHeader.H>
// ... declarations ...
#include <CD_NamespaceFooter.H>
```

---

## 4. Documentation style

### File-level Doxygen block

Every `.H` and `.cpp` file must begin with:

1. An SPDX copyright block (see § Copyright headers below).
2. A `/**` Doxygen block containing exactly `@file`, `@brief`, and `@author`.

For implementation files `@brief` must read `"Implementation of <header>.H"`:

```cpp
/**
  @file   CD_Foo.cpp
  @brief  Implementation of CD_Foo.H
  @author Robert Marskar
*/
```

### Doxygen comment style

- Use `/**` for all Doxygen comment blocks, **never** `/*!`.
- Every function — public, protected, or private — must have at minimum a `@brief`.
- Expand into `@brief` + `@details` where the behaviour is non-obvious.
- Document **all** parameters with `@param[in]`, `@param[out]`, or `@param[in,out]`.
  The direction tag is mandatory; `@param[inout]` is **invalid** — use `@param[in,out]`.
- Document return values with `@return`.
- Document non-trivial protected and private member variables with at least a `@brief`.

### Copyright headers

Every `.H` and `.cpp` file must carry a REUSE-compliant SPDX block at the very top:

```cpp
/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
```

Use `reuse annotate` to apply it, then remove any leftover legacy `/* chombo-discharge … */`
block manually.

**Do not** add SPDX comment headers to `.options` or `.inputs` files — `REUSE.toml` covers
them via glob patterns.

### Options and inputs file headers

Every `.options` and `.inputs` file must begin with a three-line banner:

```
# ====================================================================================================
# ClassName class options
# ====================================================================================================
```

`ClassName` must match the C++ class name exactly. Nothing may appear before this block.

---

## 5. Pre-commit and CI

### Running pre-commit locally

```bash
pip install pre-commit
pre-commit install          # install hooks into .git/hooks
pre-commit run --all-files  # run all hooks on every file
```

### Active hooks (`.pre-commit-config.yaml`)

| Hook | What it checks |
|------|----------------|
| `clang-format` | C++ formatting (`.clang-format` config) |
| `clang-tidy` | Static analysis (requires `compile_commands.json`; skipped in CI) |
| `reuse` | REUSE/SPDX licence compliance |
| `codespell` | Spelling in `Source/`, `Docs/`, `Exec/`, `Physics/`, `Geometries/` |
| `format-input-files` | Banner comment format in `.options`/`.inputs` |
| `check-literalincludes` | Validates all `.. literalinclude::` paths in RST files |
| `doxygen-check` | Runs `doxygen Docs/doxygen.conf`; fails if warnings remain |

### CI jobs (`.github/workflows/CI.yml`)

| Job | Runs on | Purpose |
|-----|---------|---------|
| `Formatting` | Ubuntu | `clang-format` diff check |
| `REUSE` | Ubuntu | `reuse lint` |
| `Codespell` | Ubuntu | `codespell` spelling check |
| `Linux-GNU` | Ubuntu | Full 2-D build with GCC (needs Formatting + REUSE + Codespell) |
| `Linux-oneAPI` | Ubuntu | Full 2-D build with Intel oneAPI (same dependencies) |
| `Build-documentation` | Ubuntu | Doxygen HTML + Sphinx HTML; literalinclude validation |
| `CI-passed` | — | Final gate job that all others must satisfy |

### Codespell ignore list (`.codespellignore`)

Entries must be **lowercase** (codespell lowercases tokens before lookup). Current minimum
set needed for `Source/`, `Physics/`, `Exec/`, `Docs/`, `Geometries/`:

```
ans       # LevelData<EBCellFAB> ans variable in CD_MFHelmholtzOp.cpp
ba        # const int BA = p[B] + Z in CD_PerlinSdf.cpp
fpr       # FPR type alias in ItoKMC
hashi     # bool hasHi = ... (lowercased: hashi) in several AmrMesh/Elliptic files
inout     # GeometryService::InOut Chombo type in CD_ScanShop
lod       # const Real loD = ... in CD_CdrSolver.cpp
visiter   # EBGeometry::BVH::Visiter<> external API type
```
