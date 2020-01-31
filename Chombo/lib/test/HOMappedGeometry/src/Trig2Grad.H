#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Trig2Grad.H
// petermc, 3 Sep 2010

#ifndef _TRIG2GRAD_H_
#define _TRIG2GRAD_H_

#include "VectorFunctionOnSpace.H"

#include "NamespaceHeader.H"

/// Class for cylindrical harmonic function
class Trig2Grad : public VectorFunctionOnSpace
/** computes gradient of cos(pi*x) + cos(pi*y) + cos(pi*z).
 */
{

public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor leaves Trig2Grad undefined.
  */
  Trig2Grad();

  ///
  /** Destructor.
   */
  virtual ~Trig2Grad();

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
  virtual RealVect setVectorFunctionMax();

private:

  void operator = (const Trig2Grad&)
  {
  }

  Trig2Grad(const Trig2Grad&)
  {
  }
};

#include "NamespaceFooter.H"

#endif
