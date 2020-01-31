#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// CylindricalHarmonic.H
// petermc, 21 Mar 2008

#ifndef _CYLINDRICALHARMONIC_H_
#define _CYLINDRICALHARMONIC_H_

#include "FunctionOnSpace.H"

#include "NamespaceHeader.H"

/// Class for cylindrical harmonic function
class CylindricalHarmonic : public FunctionOnSpace
/** In 2D, computes r^n * cos(n*theta);
    in 3D, computes J_n(k*rho) * cos(n*theta) * exp(-k*z).
 */
{

public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor leaves CylindricalHarmonic undefined.
  */
  CylindricalHarmonic();

  ///
  /**
     Full constructor.
  */
  CylindricalHarmonic(/// n for angle in function
                      int          a_funN,
                      /// k for exponent in cylindrical
                      int          a_funK)
  {
    define(a_funN, a_funK);
  }

  ///
  /** Destructor.
   */
  virtual ~CylindricalHarmonic();

  virtual void define(/// n for angle in function
                      int          a_funN,
                      /// k for exponent in cylindrical
                      int          a_funK);

  ///
  /** Set a_funFab on a_bx to function values at points with
      physical coordinates a_coordsFab.
   */
  virtual void setFunctionFromPhysical(/// function values, 1 component
                                       FArrayBox&         a_funFab,
                                       /// indices on which to compute function
                                       const Box&         a_bx,
                                       /// physical coordinates, SpaceDim components
                                       const FArrayBox&   a_coordsFab);

  ///
  /** Return maximum value of function on the multiblock domain.
   */
  virtual Real setFunctionMax(/// physical side length of central cube
                              Real   a_bxWidth,
                              /// outer radius of circle
                              Real   a_outerRadius);

protected:

  int m_funN, m_funK;

private:

  void operator = (const CylindricalHarmonic&)
  {
  }

  CylindricalHarmonic(const CylindricalHarmonic&)
  {
  }
};

#include "NamespaceFooter.H"

#endif
