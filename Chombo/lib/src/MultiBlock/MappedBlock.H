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

#ifndef _MAPPEDBLOCK_H_
#define _MAPPEDBLOCK_H_

#include <iostream>
#include "Tuple.H"
#include "Box.H"
#include "Misc.H"
#include "SPACE.H"
#include "BlockBoundary.H"
#include "BlockMap.H"
#include <map>
#include <string>

#include "NamespaceHeader.H"

// class BlockMap;

// class BlockBoundary;

// ---------------------------------------------------------
///
class MappedBlock
{
public:

  ///
  MappedBlock();


  ///
  /**
     Requires 2*CH_SPACEDIM BlockBoundary objects in a_boundaries.

     a_map cannot be NULL.  if you wish to use the identity mapping, you
     can use &BlockMap::Identity for this argument.
  */
  MappedBlock(const Box& a_domain,
              const BlockMap* a_map,
              const Tuple<BlockBoundary, 2*SpaceDim>& a_boundaries);

  void define(const Box& a_domain,
              const BlockMap* a_map,
              const Tuple<BlockBoundary, 2*SpaceDim>& a_boundaries);

  ///
  bool hasMap() const
  { return m_map != NULL;}

  const BlockMap& map() const
  {
    CH_assert(m_map != NULL);
    return *m_map;
  }


  ///
  Box box() const
  {return m_domain;}

  ///
  const BlockBoundary& boundary(int face) const
  { return m_boundaries[face]; }

  ///
  friend MappedBlock refine(const MappedBlock& a_mb,
                            int   a_refinement_ratio);

protected:

  friend class MappedDomain;

  Box m_domain;

  const BlockMap* m_map;

  Tuple<BlockBoundary, 2*SpaceDim> m_boundaries;
};

#include "NamespaceFooter.H"

#endif // include guard
