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

## Header guards

Header guards must be fully uppercased. For example, a file named `CD_FooBar.H` uses:

```cpp
#ifndef CD_FOOBAR_H
#define CD_FOOBAR_H
// ...
#endif
```
