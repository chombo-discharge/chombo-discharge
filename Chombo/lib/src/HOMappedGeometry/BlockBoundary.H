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
#include "IndicesTransformation.H"

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
  */
  BlockBoundary(const IndicesTransformation& a_transformation,
                int a_neighbor);

  ///
  /*
    Same as constructor.
   */
  void define(const IndicesTransformation& a_transformation,
              int a_neighbor);

  ///
  /**
     conformal mapping with unchanged coordinate system
   */
  void defineConformal(int a_neighbor);

  ///
  /**
     Returns true if boundary is an interface between two blocks or false
     otherwise (e.g., a domain boundary)
   */
  bool isInterface() const
  {
    return ( (type() == BlockBoundary::CONFORMAL) ||
             (type() == BlockBoundary::MAPPED) );
  }

  friend BlockBoundary coarsen(const BlockBoundary&   a_bb,
                               int                    a_refinement_ratio);

  friend BlockBoundary refine(const BlockBoundary&   a_bb,
                              int                    a_refinement_ratio);


  int neighbor() const
  {
    return m_neighbor;
  }

  IndicesTransformation getTransformation() const
  {
    return m_transformation;
  }

  bool isDomainBoundary() const
  {
    return (m_type == BOUNDARY);
  }

  bool isConformal() const
  {
    return (m_type == CONFORMAL);
  }

  bool isMapped() const
  {
    return (m_type == MAPPED);
  }

  ///
  /**
     Returns 1 if flux in direction a_dir points the same way on the
     other side of the block boundary, or -1 if it points the opposite way.
   */
  int reorientFace(int a_dir) const;

  ///
  /**
     Returns the direction in the block on the other side of the
     block boundary corresponding to a_dir in this block.
   */
  int dirOther(int a_dir) const;

  enum btype
  {
    BOUNDARY, CONFORMAL, MAPPED, UNKNOWN
  };

  btype type() const
  {
    return m_type;
  }

protected:

  IndicesTransformation m_transformation;

  // void*   m_bc;
  int     m_bc;
  int     m_neighbor;
  btype   m_type;

};

#include "NamespaceFooter.H"

#endif // include guard
