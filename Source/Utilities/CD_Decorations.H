/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_SIMD.H
  @brief  Declaration of various decoration macros. 
  @author Robert Marskar
*/

#ifndef CD_Decoration_H_
#define CD_Decoration_H_

#include <CD_NamespaceHeader.H>

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

#include <CD_NamespaceFooter.H>
#endif
