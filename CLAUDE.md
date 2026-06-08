# chombo-discharge — Claude Code guidelines

## Doxygen documentation

When revising or adding Doxygen documentation to files in `Source/` or `Physics/`:

- Use `/**` for all Doxygen comment blocks, never `/*!`.
- Expand brief descriptions into meaningful `@brief` + `@details` blocks where appropriate.
- Document all parameters with `@param[in]`, `@param[out]`, or `@param[in,out]`, and return values with `@return`.
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
