#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CARTESIANBLOCKCS_H_
#define _CARTESIANBLOCKCS_H_

// #include <utility>
// using std::pair;

#include "NewFourthOrderCoordSys.H"

#include "NamespaceHeader.H"

/// Single Cartesian block
/**
   The CartesianBlockCS class implements the
   NewFourthOrderCoordSys class for a cartesian block.
*/
class CartesianBlockCS : public NewFourthOrderCoordSys
{
public:
  /// constructor
  /**
     a_blockNum:  block number
     a_dx:  cell spacing
     a_baseCell:  index of bottom-left corner of this block in mapped space.
  */
  CartesianBlockCS(int a_blockNum,
                   const IntVect& a_blockIndex,
                   const RealVect& a_dx,
                   const Box& a_blockBox);

  /**
     Destructor.
   */
  virtual ~CartesianBlockCS();

  /// given coordinate in mapped space, return its location in real space
  RealVect realCoord(const RealVect& a_Xi) const;

  /// given coordinate in real space, return its location in the mapped space
  RealVect mappedCoord(const RealVect& a_x) const;

  // Calculate the derivative of each coordinate vector
  virtual Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;

  virtual void volFlux(FluxBox& a_volFlux,
                       const FluxBox& a_Nt,
                       const Box& a_box) const;

protected:

  // Block id (this object has different behaviour depending on the number)
  int m_blockNum;

  // Block index in [0:1]^D for DoubleCartesian, 0:count-1 generally
  IntVect m_blockIndex;

  // Box containing the block in mapped space (where blocks are separated)
  Box m_blockBox;
};

#include "NamespaceFooter.H"

#endif // include guard
