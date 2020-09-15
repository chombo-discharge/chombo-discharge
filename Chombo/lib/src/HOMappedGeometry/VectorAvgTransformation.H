#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _VECTORAVGTRANSFORMATION_H_
#define _VECTORAVGTRANSFORMATION_H_

#include "VectorTransformation.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
/// Class to describe transformation of SpaceDim-vectors, typically from one block to another.
class VectorAvgTransformation
{
public:

  /// null constructor leaves object in undefined state
  VectorAvgTransformation()
  { m_interpDims = 0; }

  /// null constructor leaves object in undefined state
  VectorAvgTransformation(const VectorTransformation& a_tfmAvg,
                          const Vector<VectorTransformation>& a_tfmGrad);

  /// apply this tranformation to vector
  RealVect transform(const RealVect& a_vec,
                     const Vector<RealVect>& a_gradVec) const;

  /// increment a vector-transformed cell-averaged FAB
  void addTransformFab(FArrayBox& a_dstFab,
                       const FArrayBox& a_srcFab,
                       const FArrayBox& a_gradFab,
                       const Box& a_bx,
                       Real a_weight = 1.) const;

protected:

  int m_interpDims;

  Vector<Real> m_coeffsAll;
};

#include "NamespaceFooter.H"

#endif // include guard
