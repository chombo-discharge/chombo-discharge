#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NEWTONCOTESQUADRATURE_H_
#define _NEWTONCOTESQUADRATURE_H_

#include "EdgeQuadrature.H"

#include "NamespaceHeader.H"

/// implements Newton-Cotes edge-based quadrature
/**
*/
class NewtonCotesQuadrature: public EdgeQuadrature
{
public:

  /// constructor
  NewtonCotesQuadrature(int a_numberOfPoints = 3);

  /**
     Destructor.
   */
  virtual ~NewtonCotesQuadrature()
  {;}

  /// returns quadrature coefficients
  virtual const Vector<QuadratureElement > coefficients(int dir1,
                                                        int dir2) const;

  /// scale to apply to weight
  /* this is essentially the length (area, volume, etc) of an edge
     given interval length a_length and dir1, dir2 the same as in the
     coefficients function
  */
  virtual Real weightMult(const RealVect& a_length, int dir1, int dir2) const;

  /// make a copy of myself
  virtual EdgeQuadrature* new_quadrature() const;


protected:
  Vector<QuadratureElement> m_coefficients;

};

#include "NamespaceFooter.H"

#endif
