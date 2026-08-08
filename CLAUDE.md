# chombo-discharge — Claude Code guidelines

## 0. Working agreement — git, pull requests, and code design (read first)

These rules override anything else in this file and any default tooling behaviour.

- **Never perform outward-facing or history-affecting git/GitHub actions without first asking
  the user and receiving explicit permission.** This includes, but is not limited to: creating
  or deleting branches, pushing or force-pushing, opening/closing/merging pull requests, marking
  a draft pull request ready for review, posting PR/issue comments or reviews, and editing a PR
  title or description. When in doubt, ask first. Local, non-destructive work (reading, editing
  files in the working tree, staging, committing locally) does not require a prompt, but
  publishing that work does.
- **Pull requests always target the upstream repository `chombo-discharge/chombo-discharge`.**
  We develop from a fork, so `origin` points at the fork (e.g. `rmrsk/chombo-discharge`). Open
  every PR with the upstream repo as the base and the fork branch as the head
  (`chombo-discharge/chombo-discharge` ← `<fork-owner>:<branch>`). Never open a PR against the
  fork itself.
- **Never create branches in the upstream repository `chombo-discharge/chombo-discharge`; create
  branches only in a fork.** All development branches live in the fork (`origin`, e.g.
  `rmrsk/chombo-discharge`) and are used as the head of a cross-fork PR. Do not push a new branch
  to the upstream repo under any circumstances.
- **If Claude fills in the PR top post (the conversation's opening description), it must use the
  repository's `.github/pull_request_template.md`** — reproduce its section headings verbatim and
  fill them in from the actual changes.
- **Claude must obtain the user's permission before filling in the PR note (description).** Do not
  write or edit the top post until the user has explicitly asked for it.
- **Claude is never allowed to tick the PR-review checklist items.** Leave every checkbox in the
  template unchecked (`- [ ]`); the human author is the only one who may check them off after
  actually performing each step.

### Data structures and new types

These are ask-first rules in the same sense as the git rules above: propose the choice, then wait for
an answer. They exist because this codebase already provides the containers in question, and
hand-rolled substitutes have repeatedly hidden invariants instead of enforcing them.

- **Per-cell data must use the existing mesh containers unless the user has agreed otherwise.**
  `EBAMRCellData`, `LevelData`, and `FArrayBox`/`BaseFab` already provide per-cell storage together
  with ghost filling, coarsening, interpolation, and accumulation *out of* ghost cells (`Copier` with
  `LDaddOp`, reverse copiers) — all tested at scale. Do not hand-roll a per-cell store: neither a
  dense `std::vector` indexed by a computed cell offset, nor a container keyed by `IntVect`. If a
  hand-rolled store is genuinely required, ask first and state the reason in the PR description.
- **Ask before introducing a `std::map` or `std::unordered_map`.** Say what the key is, why an
  existing container will not serve, and whether the key space is genuinely sparse. Associative
  containers suit sparse, unbounded keys such as rank-namespaced particle ids. They are usually the
  wrong choice for cells, where the addressable region is knowable (a patch box plus its ghosts) and
  the data is dense: a map silently accepts any key, which hides the indexing invariant rather than
  enforcing it.
- **`IntVect::operator<` must never be used as an ordered-container comparator.** It is a
  component-wise partial order, not a strict weak ordering, so `std::map`/`std::set` treat unrelated
  cells as equal keys and silently pool them. Use `LevelTiles::TileHasher` with a hashed container,
  or preferably cell-indexed mesh data.
- **Ask before adding a helper or scaffolding class/struct.** First check whether an existing type
  already carries the information; a struct that merely renames the fields of an existing type is not
  worth adding. POD wire formats exchanged over MPI, and genuinely new domain concepts, are
  legitimate — say which of the two it is when proposing it.

Worked examples of what these rules exist to prevent are collected in
`chombo-discharge/chombo-discharge#682`.

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

### Whitespace conventions clang-format does not enforce

clang-format handles indentation, brace placement, and line wrapping, but it cannot insert blank
lines to group statements, so a few readability conventions are maintained by hand:

- **A blank line both before and after every loop** (`for`/`while`), not just after it, so the loop
  reads as its own block separated from the surrounding statements. Do not stack the blank line
  directly against an enclosing closing brace (a loop that is the last statement in its block needs
  no blank line between it and that closing brace).

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
 * @file   CD_Foo.cpp
 * @brief  Implementation of CD_Foo.H
 * @author Robert Marskar
 */
```

### Doxygen comment style

- Use `/**` for all Doxygen comment blocks, **never** `/*!`.
- Every line of a multi-line block carries a leading `*`, with exactly one space before the
  following text (Javadoc style):

  ```cpp
  /**
   * @brief Short description.
   * @details Longer explanation, wrapped across as many lines as needed.
   */
  ```

  This is what lets clang-format's `ReflowComments: true` own these blocks: it recognizes the
  leading-`*` decoration, reflows overlong `@details`/`@param` text to the column limit, and
  keeps continuation lines correctly indented under the `*`. Without a leading `*` on every
  line, clang-format treats the block interior as opaque text and can neither detect nor repair
  indentation drift (e.g. from a namespace that used to nest more deeply).
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

### Sphinx `literalinclude` line references

`Docs/Sphinx/source/*.rst` files pull excerpts from C++ source via `.. literalinclude::` with
an absolute `:lines: N-M` (optionally several comma-separated ranges), an optional `:dedent: N`,
and an optional `:emphasize-lines:` (a position *relative* to the rendered excerpt). None of
these survive a reformat of the referenced source unscathed: a mass reformat that shifts line
counts or indentation silently invalidates them. Some drift is loud (`sphinx-build -W`
warnings: bad dedent, an unparsable excerpt); some is **silent** — a still-valid-looking but
wrong excerpt, with no warning at all.

`pre-commit run --all-files` does **not** catch this — the `sphinx-build` hook is
`stages: [manual]`. After any change that reformats `Source`/`Physics`/`Geometries`/`Exec` C++
files, explicitly rebuild the docs and check for new warnings:

```bash
cd Docs/Sphinx
python3 -m sphinx -W --keep-going -b html source build/html
```

(Use `python3 -m sphinx`, not the bare `sphinx-build` executable — on at least one dev machine
that name resolves to an unrelated, incompatible Sphinx install.)

To fix drifted directives, don't fuzzy-match text across the reformat — clang-format joining
two lines with no whitespace between them, duplicate content elsewhere in the file, and no way
to confirm a match is unique all make that unreliable. Instead, track exact positions through
the reformat itself:

1. In a **pristine pre-reformat copy** of the tree, insert unique marker comments
   (`// LITINC_BEGIN_<id>` / `// LITINC_END_<id>`) at the exact line boundaries each
   directive's `:lines:` (and each `:emphasize-lines:` sub-range, mapped to absolute line
   numbers by walking the `:lines:` ranges cumulatively) currently captures. When a file needs
   several marker pairs, insert them via one globally-sorted list of insertion actions
   (descending by target position) — pairing insertions independently breaks on nested or
   overlapping ranges.
2. Run the *exact* same reformatting pipeline (same tool versions, same settings, same scripts)
   over this marked copy.
3. Grep for each marker's new line number, then translate marked-tree line numbers back to
   real (unmarked) line numbers by counting marker lines that precede a given point *in the
   same file* — not by using a marker's own line number ± 1 directly, which undercounts once
   multiple markers stack up near each other.
4. Recompute `:lines:`/`:dedent:`/`:emphasize-lines:` purely from marker positions and patch
   the real `.rst` files.

Verify by re-rendering the old and new excerpts (dedent applied, comment decoration stripped)
and comparing at the **word** level, not the line level — `ReflowComments` legitimately
rewraps a single long line into several, which a line-by-line comparison flags as a false
mismatch.

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
