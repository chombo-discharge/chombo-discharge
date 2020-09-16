#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// VectorFunctionOnSpace.H
// petermc, 30 Jun 2010

#ifndef _VECTORFUNCTIONONSPACE_H_
#define _VECTORFUNCTIONONSPACE_H_

#include "FArrayBox.H"
#include "RealVect.H"

#include "NamespaceHeader.H"

/// Class for functions to be used with mapped grids.
class VectorFunctionOnSpace
/** Class for functions to be used with mapped grids.
 */
{

public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor leaves VectorFunctionOnSpace undefined.
  */
  VectorFunctionOnSpace()
  {
    m_isDefined = false;
  }

  ///
  /** Destructor.
   */
  virtual ~VectorFunctionOnSpace()
  {
  }

  ///
  /** Set a_funFab on a_bx to function values at points with
      physical coordinates a_coordsFab.
   */
  virtual void setVectorFunctionFromPhysical(/// function values, SpaceDim components
                                             FArrayBox&         a_funFab,
                                             /// indices on which to compute function
                                             const Box&         a_bx,
                                             /// physical coordinates, SpaceDim components
                                             const FArrayBox&   a_coordsFab) = 0;

  RealVect functionFromPhysical(/// physical coordinates
                                const RealVect&   a_rv);

  ///
  /** Return maximum value of function on the multiblock domain.
   */
  //  virtual RealVect setVectorFunctionMax(/// physical side length of central cube
  //                                        Real   a_bxWidth,
  //                                        /// outer radius of circle
  //                                        Real   a_outerRadius) = 0;

  ///
  /** Whether the function has been defined yet.
   */
  virtual bool isDefined()
  {
    return m_isDefined;
  }

protected:

  bool m_isDefined;

private:

  void operator = (const VectorFunctionOnSpace&)
  {
  }

  VectorFunctionOnSpace(const VectorFunctionOnSpace&)
  {
  }
};

#include "NamespaceFooter.H"

#endif
