#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Normals.H
// petermc, 19 Aug 2004

#ifndef _NORMALS_H_
#define _NORMALS_H_

#include <cstdlib>
#include <iostream>
#include "SPACE.H"
#include <cmath>
#include "FArrayBox.H"
#include "RealVect.H"
#include "Vector.H"
#include "Tuple.H"
#include "REAL.H"

#include "NamespaceHeader.H"

/// Computes normal derivative on faces of a rectangular domain.
class Normals
/** Computes normal derivative on faces of a rectangular domain.
*/
{

public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor leaves Normals undefined.
  */
  Normals();

  ///
  /**
     Full constructor.
     Calls full define() function with identical arguments.
  */
  Normals(/// NODE-centered domain
          const Box&        a_bx,
          /// grid spacing
          const RealVect&   a_dx,
          /// degree of normal derivative approximation
          int               a_degree);

  ///
  /** Destructor.
   */
  virtual ~Normals();

  ///
  /**
     Full define function.
  */
  virtual void define(/// NODE-centered domain
                      const Box&        a_bx,
                      /// grid spacing
                      const RealVect&   a_dx,
                      /// degree of normal derivative approximation
                      int               a_degree);

  /*@}*/

  /**
     \name Access functions
  */
  /*@{*/

  ///
  /** Returns <tt>true</tt> if this object was created with the defining
      constructor or if define() has been called.
  */
  bool isDefined() const;


  /*@}*/

  /**
     \name Parameter-setting functions
  */
  /*@{*/

  ///
  /** Set verbosity.  Default is 0.
   */
  void setVerbose(int a_verbose);

  /*@}*/

  /**
     \name Data modification functions
  */
  /*@{*/

  ///
  /** Computes the discrete normal derivative on all faces of the domain.
  */
  void evalNormal(/// normal derivatives.  order xlo, ylo, lo, xhi, yhi, zhi
                  Vector<FArrayBox*>&   a_deriv,
                  /// data from which to take normal derivatives
                  const FArrayBox&      a_phi);

  ///
  /** Computes the discrete normal derivative on one face of the domain.
  */
  void evalNormal(/// normal derivatives on one face
                  FArrayBox&              a_deriv,
                  /// data from which to take normal derivative
                  const FArrayBox&        a_phi,
                  /// face
                  int                     a_idir,
                  /// face
                  const Side::LoHiSide&   a_side);

  ///
  /** Computes the discrete normal derivative on all faces of the domain.
  */
  void evalNormalExtra(/// normal derivatives.  order xlo, ylo, lo, xhi, yhi, zhi
                       Vector<FArrayBox*>&   a_deriv,
                       /// data from which to take normal derivatives
                       const FArrayBox&      a_phi,
                       /// extra layer of data beyond all faces.  order xlo, ylo, lo, xhi, yhi, zhi
                       const Vector<FArrayBox*>&   a_phiExtra);

  ///
  /** Computes the discrete normal derivative on one face of the domain.
  */
  void evalNormalExtra(/// normal derivatives on one face
                       FArrayBox&              a_deriv,
                       /// data from which to take normal derivative
                       const FArrayBox&        a_phi,
                       /// extra layer of data beyond face
                       const FArrayBox&        a_phiExtra,
                       /// face
                       int                     a_idir,
                       /// face
                       const Side::LoHiSide&   a_side);

  /*@}*/

protected:

  /** has this Normals been defined yet?
   */
  bool m_isDefined;

  /** NODE-centered grid
   */
  Box m_bx;

  /** mesh spacings
   */
  RealVect m_dx;

  /** verbosity; default 0
   */
  int m_verbose;

  /** degree of normal derivative approximation
   */
  int m_degree;

  /** coefficients for one-sided discrete normal derivatives
   */
  Tuple<Real*, SpaceDim> m_coeffs;

  /** coefficients for discrete normal derivatives with an extra point
   */
  Tuple<Real*, SpaceDim> m_coeffsExtra;

private:

  void operator = (const Normals&)
  {
  }

  Normals(const Normals&)
  {
  }
};

#include "NamespaceFooter.H"

#endif
