#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// TripleCartesianCS.H
// petermc, 28 Jan 2011

#ifndef _TRIPLECARTESIANCS_H_
#define _TRIPLECARTESIANCS_H_

#include "MultiBlockCoordSys.H"
#include "NamespaceHeader.H"

/// MultiBlockCoordSys consisting of three blocks in each dimension.
class TripleCartesianCS : public MultiBlockCoordSys
{
public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor.
  */
  TripleCartesianCS();

  ///
  /** Destructor.
   */
  virtual ~TripleCartesianCS();

  virtual void define(const ProblemDomain& a_levelDomain,
                      const RealVect& a_dx);

  /*@}*/

  /// block mapping conversion function
  /** given xi_k in block k, returns which block (k_valid)
      contains that location in its valid region and its
      location in that block's computational coordinate space.
  */
  virtual void blockRemapping(RealVect& a_xi_valid,
                              int& a_n_valid,
                              const RealVect& a_xiSrc,
                              int a_nSrc) const;

  virtual Vector<RealVect> displacements(const Vector<RealVect>&   a_dstCoords,
                                         const Vector<int>&        a_dstBlocks,
                                         const RealVect&           a_srcCoords,
                                         int                       a_srcBlock) const;

protected:

  const static int NUMBLOCKS = D_TERM6(3, *3, *3, *3, *3, *3);

  Tuple<IntVect, NUMBLOCKS> m_origin;

  Box m_blockIndicesBox;

  BaseFab<int> m_blockIndices;

  // define m_boundaries
  void defineBoundaries();

private:

  void operator = (const TripleCartesianCS&)
  {
  }

  TripleCartesianCS(const TripleCartesianCS&)
  {
  }
};


class TripleCartesianCSFactory : public MultiBlockCoordSysFactory
{

public:

  TripleCartesianCSFactory()
  {
  }

  virtual ~TripleCartesianCSFactory()
  {
  }

  virtual MultiBlockCoordSys* getCoordSys(const ProblemDomain& a_levelDomain,
                                          const RealVect& a_dx) const;
};

#include "NamespaceFooter.H"

#endif
