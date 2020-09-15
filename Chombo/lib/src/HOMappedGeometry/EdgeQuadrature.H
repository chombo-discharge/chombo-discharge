#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EDGEQUADRATURE_H_
#define _EDGEQUADRATURE_H_

#include <utility>
using std::pair;

#include "REAL.H"
#include "RealVect.H"
#include "Vector.H"

#include "NamespaceHeader.H"

struct QuadratureElement
{
  RealVect location;
  Real weight;
};

/// Virtual base class for defining edge-based quadratures
/**
   Note that an "edge" is something which is codimension 2 in the
   SpaceDim space (in 2D, it's a point, in 3D, it's an edge, and in 4D,
   it's a face.  This API is designed to seamlessly work in any
   dimensional space, as long as the coefficients have been defined
   correctly in the derived class.
*/
class EdgeQuadrature
{
public:

  /**
     Destructor.
   */
  virtual ~EdgeQuadrature()
  {;}

  /// returns quadrature coefficients -- intervals are from [-1,1]
  /** dim1 and dim2 define the directions in which the edge coordinate
      do not vary. In other words, in 3D, if dir1=0 and dir2=1, then
      the edge will be constant in x and y, and will vary in z.
      The elements of the "location" vector in dir1 and dir2 are zeros.
  */
  virtual const Vector<QuadratureElement > coefficients(int dir1,
                                                        int dir2) const =0;

  /// scale to apply to weight
  /* this is essentially the length (area, volume, etc) of an edge
     given interval length a_length and dir1, dir2 the same as in the
     coefficients function
  */
  virtual Real weightMult(const RealVect& a_length,
                          int dir1, int dir2) const = 0;

  /// make a copy of myself
  virtual EdgeQuadrature* new_quadrature() const = 0;

};

#include "NamespaceFooter.H"

#endif
