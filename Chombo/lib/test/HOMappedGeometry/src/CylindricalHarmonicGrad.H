#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// CylindricalHarmonicGrad.H
// petermc, 30 Jun 2010

#ifndef _CYLINDRICALHARMONICGRAD_H_
#define _CYLINDRICALHARMONICGRAD_H_

#include "VectorFunctionOnSpace.H"

#include "NamespaceHeader.H"

/// Class for cylindrical harmonic function
class CylindricalHarmonicGrad : public VectorFunctionOnSpace
/** In 2D, computes gradient of r^n * cos(n*theta);
    in 3D, computes gradient of J_n(k*rho) * cos(n*theta) * exp(-k*z).
 */
{

public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor leaves CylindricalHarmonicGrad undefined.
  */
  CylindricalHarmonicGrad();

  ///
  /**
     Full constructor.
  */
  CylindricalHarmonicGrad(/// n for angle in function
                          int          a_funN,
                          /// k for exponent in cylindrical
                          int          a_funK)
  {
    define(a_funN, a_funK);
  }

  ///
  /** Destructor.
   */
  virtual ~CylindricalHarmonicGrad();

  virtual void define(/// n for angle in function
                      int          a_funN,
                      /// k for exponent in cylindrical
                      int          a_funK);

  ///
  /** Set a_funFab on a_bx to function values at points with
      physical coordinates a_coordsFab.
   */
  virtual void setVectorFunctionFromPhysical(/// function values, Spacedim components
                                             FArrayBox&         a_funFab,
                                             /// indices on which to compute function
                                             const Box&         a_bx,
                                             /// physical coordinates, SpaceDim components
                                             const FArrayBox&   a_coordsFab);

  ///
  /** Return maximum value of function on the multiblock domain.
   */
  virtual RealVect setVectorFunctionMax(/// physical side length of central cube
                                        Real   a_bxWidth,
                                        /// outer radius of circle
                                        Real   a_outerRadius);

protected:

  int m_funN, m_funK;

private:

  void operator = (const CylindricalHarmonicGrad&)
  {
  }

  CylindricalHarmonicGrad(const CylindricalHarmonicGrad&)
  {
  }
};

#include "NamespaceFooter.H"

#endif
