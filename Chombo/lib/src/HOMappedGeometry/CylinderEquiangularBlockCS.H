#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CYLINDEREQUIANGULARBLOCKCS_H_
#define _CYLINDEREQUIANGULARBLOCKCS_H_

// #include <utility>
// using std::pair;

#include "NewFourthOrderCoordSys.H"

#include "NamespaceHeader.H"

/// Basic single-panel cubed sphere Coordinates
/**
   The CylinderEquiangularBlockCS class implements the
   NewFourthOrderCoordSys class for cylindrical panels in equiangular
   coordinates.
*/
class CylinderEquiangularBlockCS : public NewFourthOrderCoordSys
{
public:
  /// constructor
  /**
     a_blockNum:  block number
     a_dx:  cell spacing
     a_baseCell:  index of bottom-left corner of this block in mapped space.
  */
  CylinderEquiangularBlockCS(int a_blockNum,
                             const RealVect& a_dx,
                             const Box& a_blockBox);

  /**
     Destructor.
   */
  virtual ~CylinderEquiangularBlockCS();

  /// given coordinate in mapped space, return its location in real space
  RealVect realCoord(const RealVect& a_Xi) const;

  /// given coordinate in real space, return its location in the mapped space
  RealVect mappedCoord(const RealVect& a_x) const;

  // Calculate the derivative of each coordinate vector
  virtual Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;

protected:

  // Block id (this object has different behaviour depending on the number)
  int m_blockNum;

  // Box containing the block
  Box m_blockBox;
};

#include "NamespaceFooter.H"

#endif // include guard
