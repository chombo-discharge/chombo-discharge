#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CHARRAY_FWD_H_
#define _CHARRAY_FWD_H_

/******************************************************************************/
/**
 * \file
 *
 * \brief Forward declarations for class CHArray
 *
 * These class declarations have to be in a separate file to ensure that the
 * default template parameters are only specified once.  Include this file like
 * normal if all you want is a foward declaration of class CHArray.
 *
 *//***************************************************************************/

#ifdef CH_SPACEDIM
// In Chombo land
  #include "NamespaceHeader.H"
#endif

/// Array configuration parameters
enum
{
  ArZeroRow,                          ///< Zero lower bound and row-storage
  ArZeroCol,                          ///< Zero lower bound and column-storage
  ArRangeRow,                         ///< Range of subscripts and row-storage
  ArRangeCol                          ///< Range of subscripts and column-
                                      ///< storage
};

/// Array limits
enum
{
  CHArrayMaximumRank = 7              ///< Maximum supported rank
};

namespace ArSp  // CHArray support
{

// An index refers to a single dimension and a size to the overall
// size of the array.  Except for 1D arrays, indexes are smaller than
// sizes by division by at least 1 stride.

  typedef int      IIx_t;             ///< Type of int for an index.  Must be
                                      ///< signed to allow for lower bounds !=
                                      ///< 0
  typedef unsigned USz_t;             ///< Type of unsigned int for a size

  template <typename T> class DefaultArrayAlloc;
}
template <typename T,
          unsigned Rank,
          int ArConf = ArZeroCol,
          typename Alloc = ArSp::DefaultArrayAlloc<T> >
class CHArray;

#ifdef CH_SPACEDIM
  #include "NamespaceFooter.H"
#endif
#endif
