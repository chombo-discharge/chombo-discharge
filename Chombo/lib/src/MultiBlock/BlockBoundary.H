#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// bvs  06/06/06

#ifndef _BLOCKBOUNDARY_H_
#define _BLOCKBOUNDARY_H_


#include <iostream>
#include "IntVect.H"
#include "Box.H"
#include "Misc.H"
#include "SPACE.H"
#include <map>
#include <string>


#include "NamespaceHeader.H"

// ---------------------------------------------------------
/// Class to describe the multi-block topology between 2 blocks, or placehold a BC
class BlockBoundary
{
public:

  /// null constructor leaves object in type=UNKNOWN state. It can be overridden by MappedDomain
  BlockBoundary();

  /// the boundary is an external boundary condition.
  /**
     Boundary condition label.  Application is responsible for interpreting
     the pointer.
  */
  // BlockBoundary(void*  a_boundaryClass);
  BlockBoundary(int a_boundaryClass);

  ///
  /*
    Same as constructor.
   */
  // void define(void*  a_boundaryClass);
  void define(int a_boundaryClass);

  /// the boundary is a conformal connection between two blocks.
  /**
     The BlockMap between these two MappedBlocks are conformal, hence no spatial
     interpolation needs to be performed. However, the index spaces may require a
     transform.

     these topological connections are redundant between two blocks. the user
     can choose to specify all connections, and have MappedDomain verify the
     topology, or the user can provide just the sufficient non-redundant set
     and have MappedDomain construct the mirrors.

     IntVect pold:  indices of a cell in current block
     IntVect pnew:  indices of the same cell in the neighboring block

     Then for each direction idir:
     pnew[idir] = m_sign[idir]*pold[m_permutation[idir]] + m_translation[idir]
     and hence:
     pold[m_permutation[idir]] = m_sign[idir]*pnew[idir]
                               - m_sign[idir]*m_translation[idir]
  */
  BlockBoundary(const IntVect& a_permutation,
                const IntVect& a_sign,
                const IntVect& a_translation,
                int a_neighbor,
                int a_face);

  ///
  /*
    Same as constructor.
   */
  void define(const IntVect& a_permutation,
              const IntVect& a_sign,
              const IntVect& a_translation,
              int a_neighbor,
              int a_face);

  ///
  /**
     conformal mapping with unchanged coordinate system
   */
  void defineConformal(int a_neighbor,
                       int a_face);

  // default copy constructor and assign are OK.

  /// non-conformal block mating constructor
  /**
  */
  BlockBoundary(int a_neighbor,
                int a_face);

  ///
  /*
    Same as constructor.
   */
  void define(int a_neighbor,
              int a_face);

  friend BlockBoundary refine(const BlockBoundary&   a_bb,
                              int                    a_refinement_ratio);

  enum btype
  {
    BOUNDARY,
    CONFORMAL,
    MAPPED,
    UNKNOWN
  };

  btype type() const
  {
    return m_type;
  }

  int neighbor() const
  {
    return m_neighbor;
  }

  IntVect convertOldToNew(const IntVect&   a_ivOld) const;

  IntVect convertNewToOld(const IntVect&   a_ivNew) const;

  ///
  /** convert CELL-centered box
   */
  Box convertOldToNew(const Box&   a_bxOld) const;

  Box convertNewToOld(const Box&   a_bxNew) const;

  //============= serialization routines =============

  IntVect getPermutation() const
  {
    return m_permutation;
  }
  IntVect getSign() const
  {
    return m_sign;
  }
  IntVect getTranslation() const
  {
    return m_translation;
  }
protected:

  friend class MappedDomain;

  IntVect m_permutation;
  IntVect m_sign;
  IntVect m_translation;
  // void*   m_bc;
  int     m_bc;
  int     m_neighbor;
  int     m_face;
  btype   m_type;

};

#include "NamespaceFooter.H"

#endif // include guard
