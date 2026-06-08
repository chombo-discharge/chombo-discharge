# chombo-discharge — Claude Code guidelines

## File-level Doxygen block

Every `.H` and `.cpp` file must have a Doxygen block near the top containing exactly `@file`, `@brief`,
and `@author`. For implementation files the `@brief` must say "Implementation of `<header>.H`" where
`<header>.H` is the corresponding header file that actually exists in the same directory. Example:

```cpp
/**
  @file   CD_Foo.cpp
  @brief  Implementation of CD_Foo.H
  @author Robert Marskar
*/
```

## Doxygen documentation

When revising or adding Doxygen documentation to files in `Source/` or `Physics/`:

- Use `/**` for all Doxygen comment blocks, never `/*!`.
- Every function — public, protected, or private — must have at minimum a `@brief` description.
- Expand brief descriptions into meaningful `@brief` + `@details` blocks where appropriate.
- Document **all** parameters with `@param[in]`, `@param[out]`, or `@param[in,out]`, and return values with `@return`.
- Document non-trivial protected and private member variables.

## Copyright headers

Every file must carry a REUSE-compliant SPDX header instead of the old hand-written block:

```cpp
/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
```

Use `reuse annotate` to apply the header, then remove any leftover legacy copyright lines manually.

**Do not** add SPDX comment headers to `.options` or `.inputs` files. Setup scripts merge these files
and inline headers cause duplication and clutter. Instead, `REUSE.toml` at the repository root
declares their copyright and license in bulk via glob patterns — no per-file action needed.

## Options and inputs file headers

Every `.options` and `.inputs` file must begin with a banner comment that names the class whose
options are contained in that file. The format is exactly:

```
# ====================================================================================================
# ClassName class options
# ====================================================================================================
```

where `ClassName` matches the C++ class name exactly (e.g. `FieldSolverGMG`, `CdrPlasmaGodunovStepper`).
No other text may appear before this three-line block.

## Header guards

Header guards must be fully uppercased. For example, a file named `CD_FooBar.H` uses:

```cpp
#ifndef CD_FOOBAR_H
#define CD_FOOBAR_H
// ...
#endif
```

## Include ordering

Includes must be grouped and ordered as follows: standard library / third-party first, then Chombo, then our own headers. Each group is preceded by a comment label. Example:

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
