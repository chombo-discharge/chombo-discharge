#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// SphericalHarmonicGrad.H
// petermc, 30 Jul 2010

#ifndef _SPHERICALHARMONICGRAD_H_
#define _SPHERICALHARMONICGRAD_H_

#include "VectorFunctionOnSpace.H"

#include "NamespaceHeader.H"

/// Class for cylindrical harmonic function
class SphericalHarmonicGrad : public VectorFunctionOnSpace
/** In 3D, computes gradient of r^l * cos(m*theta) * P^m_l (cos phi).
 */
{

public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor leaves SphericalHarmonicGrad undefined.
  */
  SphericalHarmonicGrad();

  ///
  /**
     Full constructor.
  */
  SphericalHarmonicGrad(/// l for Legendre and power of r
                        int          a_funL,
                        /// m for Legendre and multiple of theta
                        int          a_funM)
  {
    define(a_funL, a_funM);
  }

  ///
  /** Destructor.
   */
  virtual ~SphericalHarmonicGrad();

  virtual void define(/// l for Legendre and power of r
                      int          a_funL,
                      /// m for Legendre and multiple of theta
                      int          a_funM);

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

  int m_funL, m_funM, m_legendreLen;

  Real* m_legendreCoeffs;

private:

  void operator = (const SphericalHarmonicGrad&)
  {
  }

  SphericalHarmonicGrad(const SphericalHarmonicGrad&)
  {
  }
};

#include "NamespaceFooter.H"

#endif
