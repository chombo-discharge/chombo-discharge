#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Trig2.H
// petermc, 3 Sep 2010

#ifndef _TRIG2_H_
#define _TRIG2_H_

#include "FunctionOnSpace.H"

#include "NamespaceHeader.H"

/// Class for sum of trigonometric functions
class Trig2 : public FunctionOnSpace
/** computes cos(x*pi) + cos(y*pi) + cos(z*pi).
 */
{

public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor.
  */
  Trig2();

  ///
  /** Destructor.
   */
  virtual ~Trig2();

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
  virtual Real setFunctionMax();

private:

  void operator = (const Trig2&)
  {
  }

  Trig2(const Trig2&)
  {
  }
};

#include "NamespaceFooter.H"

#endif
