#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// MultiCartesianCS.H
// petermc, 25 Jul 2016

#ifndef _MULTICARTESIANCS_H_
#define _MULTICARTESIANCS_H_

#include "MultiBlockCoordSys.H"
#include "NamespaceHeader.H"

/// MultiBlockCoordSys consisting of rectangular array of boxes, with Cartesian mapping in each
class MultiCartesianCS : public MultiBlockCoordSys
{
public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor.
  */
  MultiCartesianCS();

  ///
  /** Constructor.
  */
  MultiCartesianCS(IntVect a_count,
                   IntVect a_wraparound);

  ///
  /** Destructor.
   */
  virtual ~MultiCartesianCS();

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

  IntVect m_count;

  IntVect m_wraparound;

  Box m_blockIndicesBox;

  // number of cells in m_blockIndicesBox
  int m_numBlocks;

  // lives on m_blockIndicesBox; entry is the block index
  BaseFab<int> m_blockIndices;

  Vector<IntVect> m_origin;

  // define m_boundaries
  void defineBoundaries();

private:

  void operator = (const MultiCartesianCS&)
  {
  }

  MultiCartesianCS(const MultiCartesianCS&)
  {
  }
};


class MultiCartesianCSFactory : public MultiBlockCoordSysFactory
{

public:

  MultiCartesianCSFactory()
  {
    m_count = 2*IntVect::Unit;
    m_wraparound = IntVect::Zero;
  }

  virtual ~MultiCartesianCSFactory()
  {
  }

  virtual MultiBlockCoordSys* getCoordSys(const ProblemDomain& a_levelDomain,
                                          const RealVect& a_dx) const;

  void setCount(const IntVect& a_count)
  {
    m_count = a_count;
  }

  void setWraparound(const IntVect& a_wraparound)
  {
    m_wraparound = a_wraparound;
  }

protected:

  IntVect m_count;

  IntVect m_wraparound;
};

#include "NamespaceFooter.H"

#endif
