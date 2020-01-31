#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DoubleCartesianRotateCS.H
// petermc, 7 Jul 2014

#ifndef _DOUBLECARTESIANROTATECS_H_
#define _DOUBLECARTESIANROTATECS_H_

#include "MultiBlockCoordSys.H"

#include "NamespaceHeader.H"

/// MultiBlockCoordSys consisting of two blocks in each dimension, which are rotated from mapped to real space
class DoubleCartesianRotateCS : public MultiBlockCoordSys
/** MultiBlockCoordSys with a central cubic box and four sectors.
 */
{
public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor.
  */
  DoubleCartesianRotateCS();

  ///
  /** Destructor.
   */
  virtual ~DoubleCartesianRotateCS();

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

  /// Have separate volume flux on each block.  Normally does nothing, but double-cartesian coordinate system needs it.
  virtual void separateVolFlux(LevelData<FluxBox>& a_flux) const;
 
protected:

  const static int NUMBLOCKS = D_TERM6(2, *2, *2, *2, *2, *2);

  Tuple<IntVect, NUMBLOCKS> m_origin;

  Box m_blockIndicesBox;

  BaseFab<int> m_blockIndices;

  BaseFab<IntVect> m_perm;

  BaseFab<IntVect> m_sign;

  BaseFab<IndicesTransformation> m_indTransform;

  BaseFab<VectorTransformation> m_vecTransform;

  // define m_boundaries
  void defineBoundaries();

  int findDim(const RealVect& a_rv) const;

  int findSign(const RealVect& a_rv) const;

private:

  void operator = (const DoubleCartesianRotateCS&)
  {
  }

  DoubleCartesianRotateCS(const DoubleCartesianRotateCS&)
  {
  }
};


class DoubleCartesianRotateCSFactory : public MultiBlockCoordSysFactory
{

public:

  DoubleCartesianRotateCSFactory()
  {
  }

  virtual ~DoubleCartesianRotateCSFactory()
  {
  }

  virtual MultiBlockCoordSys* getCoordSys(const ProblemDomain& a_levelDomain,
                                          const RealVect& a_dx) const;
};

#include "NamespaceFooter.H"

#endif
