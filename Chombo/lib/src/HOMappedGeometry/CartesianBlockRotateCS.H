#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CARTESIANBLOCKROTATECS_H_
#define _CARTESIANBLOCKROTATECS_H_

// #include <utility>
// using std::pair;

#include "NewFourthOrderCoordSys.H"
#include "IndicesTransformation.H"
#include "VectorTransformation.H"

#include "NamespaceHeader.H"

/// Single Cartesian block
/**
   The CartesianBlockRotateCS class implements the
   NewFourthOrderCoordSys class for a cartesian block.
*/
class CartesianBlockRotateCS : public NewFourthOrderCoordSys
{
public:
  /// constructor
  /**
     a_blockNum:  block number
     a_dx:  cell spacing
     a_baseCell:  index of bottom-left corner of this block in mapped space.
  */
  CartesianBlockRotateCS(int a_blockNum,
                         const IntVect& a_blockIndex,
                         const RealVect& a_dx,
                         const Box& a_blockBox,
                         IntVect a_perm = IntVect(D_DECL6(0,1,2,3,4,5)),
                         IntVect a_sign = IntVect::Unit);

  /**
     Destructor.
   */
  virtual ~CartesianBlockRotateCS();

  /// given coordinate in mapped space, return its location in real space
  RealVect realCoord(const RealVect& a_Xi) const;

  /// given coordinate in real space, return its location in the mapped space
  RealVect mappedCoord(const RealVect& a_x) const;

  // Calculate the derivative of each coordinate vector
  virtual Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;

  virtual void volFlux(FluxBox& a_volFlux,
                       const FluxBox& a_Nt,
                       const Box& a_box) const;

  // From mapped space to rotated mapped space.
  RealVect rotateFwd(const RealVect& a_xi) const;

  // From rotated mapped space back to mapped space.
  RealVect rotateBack(const RealVect& a_xi) const;

  // accessor
  VectorTransformation vecTransform() const
  { return m_vecTransform; }

protected:

  // Block id (this object has different behaviour depending on the number)
  int m_blockNum;

  // Block index in [0:1]^D
  IntVect m_blockIndex;

  // Box containing the block
  Box m_blockBox;

  // permutation of this block
  IntVect m_perm;

  // signs of this block
  IntVect m_sign;

  // shift of center of block
  RealVect m_blockCenter;

  // transformation from mapped space to rotated mapped space
  VectorTransformation m_vecTransform;
};

#include "NamespaceFooter.H"

#endif // include guard
