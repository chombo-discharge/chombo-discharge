#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// petermc, 20 May 2008

#ifndef _SPHERICALDOMAIN_H_
#define _SPHERICALDOMAIN_H_

#include <iostream>
#include "SetMappedDomain.H"
#include "RealVect.H"
#include "Misc.H"
#include "SPACE.H"
#include <map>
#include <string>

#include "NamespaceHeader.H"

/** This class implements the cubed-sphere map.
 */

// ---------------------------------------------------------
class SphericalDomain : public SetMappedDomain
{
public:

  // Constructors

  /// The default constructor.  The constructed as an empty domain.
  /**
   */
  SphericalDomain();

  /// Construct SphericalDomain with a single \em a_domBox as domain of central cube
  /**
     This constructor defaults to non-periodic domain.

     block 0:  +x, 0Hi from CUBE
     block 1:  +y, 1Hi from CUBE
     block 2:  -x, 0Lo from CUBE
     block 3:  -y, 1Lo from CUBE
     block 4:  +z, 2Hi from CUBE
     block 5:  -z, 2Lo from CUBE
     block 6:  CUBE
  */
  SphericalDomain(const Box&        a_centralBox,
                  const RealVect&   a_center,
                  Real              a_bxWidth,
                  Real              a_outerRadius);

  ///
  /** Destructor.
   */
  virtual ~SphericalDomain();

  /// Construct SphericalDomain with \em a_domBox as computational domain
  /**
     This constructor defaults to non-periodic domain
  */
  virtual void define(const Box&        a_centralBox,
                      const RealVect&   a_center,
                      Real              a_bxWidth,
                      Real              a_outerRadius);

protected:

  /** index of central cubic block
   */
  const static int CUBE = 6;

  /** axis (0 or 1 or 2) of each sectoral block
   */
  Tuple<int, 2*SpaceDim> m_axis;

  /** direction (Hi or Lo) of axis of each sectoral block
   */
  Tuple<Side::LoHiSide, 2*SpaceDim> m_direction;

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
