#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// petermc, 22 May 2008

#ifndef _SETMAPPEDDOMAIN_H_
#define _SETMAPPEDDOMAIN_H_

#include <iostream>
#include "Vector.H"
#include "Tuple.H"
#include "MappedDomain.H"
#include "BlockMap.H"
#include "BlockBoundary.H"
#include "Misc.H"
#include "SPACE.H"
#include <map>
#include <string>

class MappedBlock;

#include "NamespaceHeader.H"

// ---------------------------------------------------------
class SetMappedDomain
{
public:

  // Constructors

  /// The default constructor.  The constructed as an empty domain.
  /**
   */
  SetMappedDomain();

  ///
  /** Destructor.
   */
  virtual ~SetMappedDomain();

  //  Accessors

  /// Returns true if this SetMappedDomain is empty or undefined.
  /**
  */
  bool isDefined() const;

  ///
  /** Returns the MappedDomain.
   */
  const MappedDomain& mappedDomain() const;

protected:

  virtual void defineMappedDomain();

  /** number of blocks:  was const static int in derived class
   */
  int m_nblocks;

  /** defined?
   */
  bool m_isDefined;

  /** CELL-centered block boxes, length m_nblocks
   */
  Vector<Box> m_boxes;

  /** CELL-centered full computational domain
   */
  Box m_domainBox;

  /** does the actual work
   */
  MappedDomain m_mappedDomain;

  /** length m_nblocks, used in m_mappedDomain
   */
  Vector<BlockMap*> m_blockMaps;

  /** length m_nblocks, used in m_mappedDomain
   */
  Vector<MappedBlock> m_mappedBlocks;

  /** length m_nblocks
   */
  Vector< Tuple<BlockBoundary, 2*SpaceDim> > m_blockBoundaries;
};

#include "NamespaceFooter.H"

#endif // include guard
