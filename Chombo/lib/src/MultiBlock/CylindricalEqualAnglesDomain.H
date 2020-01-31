#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// petermc, 19 May 2008

#ifndef _CYLINDRICALEQUALANGLESDOMAIN_H_
#define _CYLINDRICALEQUALANGLESDOMAIN_H_

#include "CylindricalDomain.H"

#include "NamespaceHeader.H"

/** This implements the CYLINDEREQUALANGLES map.
 */

// ---------------------------------------------------------
class CylindricalEqualAnglesDomain : public CylindricalDomain
{
public:

  // Constructors

  /// The default constructor.  The constructed as an empty domain.
  /**
   */
  CylindricalEqualAnglesDomain();

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
  CylindricalEqualAnglesDomain(const Box&        a_centralBox,
                               const RealVect&   a_center,
                               Real              a_bxWidth,
                               Real              a_outerRadius);

  ///
  /** Destructor.
   */
  virtual ~CylindricalEqualAnglesDomain();

  /// Construct CylindricalEqualAnglesDomain with \em a_domBox as computational domain
  /**
     This constructor defaults to non-periodic domain
  */
  virtual void define(const Box&        a_centralBox,
                      const RealVect&   a_center,
                      Real              a_bxWidth,
                      Real              a_outerRadius);
};

#include "NamespaceFooter.H"

#endif // include guard
