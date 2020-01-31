#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CUBEDSPHERESOLIDBLOCKCS_H_
#define _CUBEDSPHERESOLIDBLOCKCS_H_

// #include <utility>
// using std::pair;

#include "REAL.H"
#include "RealVect.H"
#include "NewFourthOrderCoordSys.H"

#include "NamespaceHeader.H"

/// Basic single-panel cubed sphere Coordinates
/// Multiblock coordinate system for single block of solid cubed sphere.
/**
   The CubedSphereSolidBlockCS class implements the NewFourthOrderCoordSys class
   for solid cubed sphere block.
*/
class CubedSphereSolidBlockCS : public NewFourthOrderCoordSys
{
public:

  /// constructor
  CubedSphereSolidBlockCS(int a_blockID,
                          const RealVect& a_dxVect,
                          const Box& a_bx);

  void setR0(Real a_r0)
  {
    CH_assert(a_r0 > 0.);
    m_r0 = a_r0;
  }

  Real getR0() const
  {
    CH_assert(m_r0 > 0.);
    return m_r0;
  }

  void setR1(Real a_r1)
  {
    CH_assert(a_r1 > 0.);
    m_r1 = a_r1;
  }

  Real getR1() const
  {
    CH_assert(m_r1 > 0.);
    return m_r1;
  }

  /**
     Destructor.
   */
  virtual ~CubedSphereSolidBlockCS();

  /// given coordinate in mapped space, return its location in real space
  RealVect realCoord(const RealVect& a_Xi) const;

  //  void realCoord(FArrayBox& a_x,
  //                 const FArrayBox& a_Xi,
  //                 const Box& a_box) const;

  /// given coordinate in real space, return its location in the mapped space
  RealVect mappedCoord(const RealVect& a_x) const;

  // Calculate the derivative of each coordinate vector
  virtual Real dXdXi(const RealVect& a_Xi,
                     int a_dirX,
                     int a_dirXi) const;

  /// access function to block id
  int blockID() const
  {
    return m_blockID;
  }

protected:

  const static int XPOS = 0;
  const static int YPOS = 1;
  const static int XNEG = 2;
  const static int YNEG = 3;
  const static int ZPOS = 4;
  const static int ZNEG = 5;
  const static int CENTRAL = 6;

  // block id (this object has different behavior depending on the panel)
  int m_blockID;

  // half-width of inner cube in physical space
  Real m_r0;

  // radius of sphere in physical space
  Real m_r1;

  // coordinates of low corner in mapped space
  RealVect m_mappedLo;

  // index (0 or 1 or 2) of real variable in radial direction
  int m_rCoord;

  // index (0 or 1 or 2) of real variable multiplied by a
  int m_aCoord;

  // index (0 or 1 or 2) of real variable multiplied by b
  int m_bCoord;

  // signs of real variables
  RealVect m_signCoord;
};

#include "NamespaceFooter.H"

#endif
