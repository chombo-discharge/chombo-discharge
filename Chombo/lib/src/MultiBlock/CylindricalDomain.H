#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// petermc, 1 Apr 2008

#ifndef _CYLINDRICALDOMAIN_H_
#define _CYLINDRICALDOMAIN_H_

#include <iostream>
#include "SetMappedDomain.H"
#include "RealVect.H"
#include "Misc.H"
#include "SPACE.H"
#include <map>
#include <string>

#include "NamespaceHeader.H"

/** Derived classes implement
    CYLINDERSPOKES, CYLINDEREQUALANGLES, and CYLINDERTRANSITION.
 */

// ---------------------------------------------------------
class CylindricalDomain : public SetMappedDomain
{
public:

  // Constructors

  /// The default constructor.  The constructed as an empty domain.
  /**
   */
  CylindricalDomain();

  /// Construct CylindricalDomain with a single \em a_domBox as domain of central cube
  /**
     This constructor defaults to non-periodic domain.

             +-------+
             |   2   |
     +-------+-------+-------+
     |   3   |   0   |   1   |
     +-------+-------+-------+
             |   4   |
             +-------+
  */
  CylindricalDomain(const Box&        a_centralBox,
                    const RealVect&   a_center,
                    Real              a_bxWidth,
                    Real              a_outerRadius);

  ///
  /** Destructor.
   */
  virtual ~CylindricalDomain();

  /// Construct CylindricalDomain with \em a_domBox as computational domain
  /**
     This constructor defaults to non-periodic domain
  */
  virtual void define(const Box&        a_centralBox,
                      const RealVect&   a_center,
                      Real              a_bxWidth,
                      Real              a_outerRadius);

protected:

  virtual void setCubeAndOuterBoundaries(int a_neighbor,
                                         int a_face);

  virtual void setInnerBoundary(int a_srcBlock,
                                int a_dstBlock,
                                int a_srcFace,
                                int a_dstFace,
                                Side::LoHiSide a_corner0side,
                                Side::LoHiSide a_corner1side);

  /** index of central cubic block
   */
  const static int CUBE = 0;

  /** width of central cubic box
   */
  Real m_bxWidth;

  /** outer radius of cylinder
   */
  Real m_outerRadius;

  /** center of cylinder and of cube
   */
  RealVect m_center;

  /** physical corners of central cubic block
   */
  RealVect m_centralCornerLo, m_centralCornerHi;
};

#include "NamespaceFooter.H"

#endif // include guard
