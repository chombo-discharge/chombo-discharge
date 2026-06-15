/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_Decorations.H
  @brief  Declaration of various decoration macros. 
  @author Robert Marskar
*/

#ifndef CD_DECORATIONS_H
#define CD_DECORATIONS_H

#include <CD_NamespaceHeader.H>

/// @cond DOXYGEN_SKIP
// SIMD instructions
#ifndef NDEBUG
#define CD_PRAGMA_SIMD

#elif defined(__INTEL_COMPILER)
#define CD_PRAGMA_SIMD _Pragma("ivdep")

#elif defined(__INTEL_LLVM_COMPILER)
#define CD_PRAGMA_SIMD _Pragma("SIMD")

#elif defined(__clang__)
#define CD_PRAGMA_SIMD _Pragma("clang loop vectorize(enable)")

#elif defined(__GNUC__)
#define CD_PRAGMA_SIMD _Pragma("GCC ivdep")

#else
#define CD_PRAGMA_SIMD

#endif

// Force inline.
#if defined(__GNUC__)
#define ALWAYS_INLINE __attribute__((always_inline)) __inline__
#else
#define ALWAYS_INLINE inline
#endif
/// @endcond

#include <CD_NamespaceFooter.H>
#endif
