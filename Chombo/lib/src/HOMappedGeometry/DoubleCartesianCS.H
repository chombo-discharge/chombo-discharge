#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DoubleCartesianCS.H
// petermc, 3 Sep 2010

#ifndef _DOUBLECARTESIANCS_H_
#define _DOUBLECARTESIANCS_H_

#include "MultiBlockCoordSys.H"
#include "NamespaceHeader.H"

/// MultiBlockCoordSys consisting of two blocks in each dimension.
class DoubleCartesianCS : public MultiBlockCoordSys
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
  DoubleCartesianCS();

  ///
  /** Destructor.
   */
  virtual ~DoubleCartesianCS();

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

  // define m_boundaries
  void defineBoundaries();

private:

  void operator = (const DoubleCartesianCS&)
  {
  }

  DoubleCartesianCS(const DoubleCartesianCS&)
  {
  }
};


class DoubleCartesianCSFactory : public MultiBlockCoordSysFactory
{

public:

  DoubleCartesianCSFactory()
  {
  }

  virtual ~DoubleCartesianCSFactory()
  {
  }

  virtual MultiBlockCoordSys* getCoordSys(const ProblemDomain& a_levelDomain,
                                          const RealVect& a_dx) const;
};

#include "NamespaceFooter.H"

#endif
