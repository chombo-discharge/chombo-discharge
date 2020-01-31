#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CYLINDERSPOKESBLOCKCS_H_
#define _CYLINDERSPOKESBLOCKCS_H_

// #include <utility>
// using std::pair;

#include "NewFourthOrderCoordSys.H"

#include "NamespaceHeader.H"

/// Basic single-panel cubed sphere Coordinates
/**
   The CylinderSpokesBlockCS class implements the
   NewFourthOrderCoordSys class for cylindrical panels in equiangular
   coordinates.
*/
class CylinderSpokesBlockCS : public NewFourthOrderCoordSys
{
public:
  /// constructor
  /**
     a_blockNum:  block number
     a_dx:  cell spacing
     a_baseCell:  index of bottom-left corner of this block in mapped space.
  */
  CylinderSpokesBlockCS(int a_blockNum,
                        const RealVect& a_dx,
                        const Box& a_blockBox);

  /**
     Destructor.
   */
  virtual ~CylinderSpokesBlockCS();

  /// given coordinate in mapped space, return its location in real space
  RealVect realCoord(const RealVect& a_Xi) const;

  /// given coordinate in real space, return its location in the mapped space
  RealVect mappedCoord(const RealVect& a_x) const;

  // Calculate the derivative of each coordinate vector
  virtual Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;

protected:

  const static int CENTRAL = 0;
  const static int XPOS = 1;
  const static int YPOS = 2;
  const static int XNEG = 3;
  const static int YNEG = 4;

  // Block id (this object has different behaviour depending on the number)
  int m_blockNum;

  // Box containing the block
  Box m_blockBox;

  // coordinates of low corner in mapped space
  RealVect m_mappedLo;

  // index (0 or 1) of real variable in radial direction
  int m_rCoord;

  // index (0 or 1) of real variable multiplied by a
  int m_aCoord;

  // signs of real variables
  RealVect m_signCoord;
};

#include "NamespaceFooter.H"

#endif // include guard
