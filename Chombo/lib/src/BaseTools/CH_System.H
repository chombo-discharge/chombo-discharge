#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CH_SYSTEM_H_
#define _CH_SYSTEM_H_

/******************************************************************************/
/**
 * \file
 *
 * \brief System dependent functions
 *
 *//*+*************************************************************************/

#include <cstring>

#include "CH_config.H"
#include "BaseNamespaceHeader.H"

namespace CHSystem
{
  /// Check if a file exists
  int fileExists(const char *const a_filename);

  /// Allocate aligned memory
  int memalign(void **a_memptr, size_t a_alignment, size_t a_size);

/*--------------------------------------------------------------------*/
//  Find the number of bits set to 1 in a 32 bit int
/**
 *  \param[in]  a_i     Integer
 *  \return             Number of bits set to 1
 *//*-----------------------------------------------------------------*/

  inline int popcnt(unsigned i)
  {
#if CHDEF_SYSTEM_POPCNT_METHOD == 1
    int ret;
    __asm__ ("popcnt %1, %0" : "=r" (ret) : "r" (i));
    return ret;
#else
    // Known as a parallel SWAR algorithm
    i = i - ((i >> 1) & 0x55555555);
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
    return (((i + (i >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24;
#endif
  }
}

// BOOST_PP_STRINGIZE(expr) - Return expr as a string literal.
#define BOOST_PP_STRINGIZE_(expr) #expr
#define BOOST_PP_STRINGIZE(expr) BOOST_PP_STRINGIZE_(expr)

// BOOST_PRAGMA(args) - Emits a pragma.
#define BOOST_PRAGMA(args) _Pragma(BOOST_PP_STRINGIZE(args))
#if   defined(__INTEL_COMPILER)
    #define BOOST_BUILTIN_ASSUME(expr)                                         \
        __assume(expr)                                                         \
        /**/
    #define BOOST_BUILTIN_ASSUME_ALIGNED(ptr, align)                           \
        __assume_aligned(ptr, align)                                           \
        /**/
#elif defined(__clang__)
    #define BOOST_BUILTIN_ASSUME(expr)                                         \
        __builtin_assume(expr)                                                 \
        /**/
    #define BOOST_BUILTIN_ASSUME_ALIGNED(ptr, align)                           \
        {                                                                      \
            ptr = reinterpret_cast<decltype(ptr)>(                             \
                __builtin_assume_aligned(ptr, align)                           \
            );                                                                 \
        }                                                                      \
        /**/
#else // GCC
    #define BOOST_BUILTIN_ASSUME(expr)                                         \
        do { if (!(expr)) __builtin_unreachable(); } while (0)                 \
        /**/
    #define BOOST_BUILTIN_ASSUME_ALIGNED(ptr, align)                           \
        {                                                                      \
            ptr = reinterpret_cast<decltype(ptr)>(                             \
                __builtin_assume_aligned(ptr, align)                           \
            );                                                                 \
        }                                                                      \
        /**/
#endif
#if   defined(__INTEL_COMPILER)
    #define BOOST_DEMAND_VECTORIZATION                                         \
        BOOST_PRAGMA(simd)                                                     \
        /**/
#elif defined(__clang__)
    #define BOOST_DEMAND_VECTORIZATION                                         \
        BOOST_PRAGMA(clang loop vectorize(enable) interleave(enable))          \
        /**/
#else
    #define BOOST_DEMAND_VECTORIZATION                                         \
        BOOST_PRAGMA(GCC ivdep)                                                \
        /**/
#endif

#include "BaseNamespaceFooter.H"
#endif
