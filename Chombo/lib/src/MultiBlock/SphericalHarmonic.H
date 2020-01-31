#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// SphericalHarmonic.H
// petermc, 21 Mar 2008

#ifndef _SPHERICALHARMONIC_H_
#define _SPHERICALHARMONIC_H_

#include "FunctionOnSpace.H"

#include "NamespaceHeader.H"

/// Class for cylindrical harmonic function
class SphericalHarmonic : public FunctionOnSpace
/** In 3D, computes r^l * cos(m*theta) * P^m_l (cos phi).
 */
{

public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor leaves SphericalHarmonic undefined.
  */
  SphericalHarmonic();

  ///
  /**
     Full constructor.
  */
  SphericalHarmonic(/// l for Legendre and power of r
                    int          a_funL,
                    /// m for Legendre and multiple of theta
                    int          a_funM)
  {
    define(a_funL, a_funM);
  }

  ///
  /** Destructor.
   */
  virtual ~SphericalHarmonic();

  virtual void define(/// l for Legendre and power of r
                      int          a_funL,
                      /// m for Legendre and multiple of theta
                      int          a_funM);

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

  int m_funL, m_funM, m_legendreLen;

  Real* m_legendreCoeffs;

private:

  void operator = (const SphericalHarmonic&)
  {
  }

  SphericalHarmonic(const SphericalHarmonic&)
  {
  }
};

#include "NamespaceFooter.H"

#endif
