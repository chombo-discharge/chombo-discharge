#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// PoissonDirichlet.H
// petermc, 19 Aug 2004

#ifndef _POISSONDIRICHLET_H_
#define _POISSONDIRICHLET_H_

#include <cstdlib>
#include <iostream>
#include "SPACE.H"
#include <cmath>
#include "FArrayBox.H"
#include "RealVect.H"
#include "REAL.H"

#include "NamespaceHeader.H"

/// Solves Poisson's equation on a rectangular domain with Dirichlet boundary conditions.
class PoissonDirichlet
/** PoissonDirichlet solves Poisson's equation on a rectangular domain with Dirichlet boundary conditions.
*/
{

public:

  ///
  /**
     SPECTRAL is the spectral operator.<br>
     POINT19 is the 19-point operator.<br>
     POINT7 is the 7-point operator.<br>
     POINT27 is the 27-point operator.
   */
  enum OperatorType
  {
    SPECTRAL = 0,
    POINT19  = 1,
    POINT7   = 2,
    POINT27  = 3
  };

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor leaves PoissonDirichlet undefined.
  */
  PoissonDirichlet();

  ///
  /**
     Full constructor.
     Calls full define() function with identical arguments.
  */
  PoissonDirichlet(/// NODE-centered domain
                   const Box&                         a_bx,
                   /// grid spacing
                   const RealVect&                    a_dx,
                   /// which operator in OperatorType
                   OperatorType                       a_op);

  ///
  /** Destructor.
   */
  virtual ~PoissonDirichlet();

  ///
  /**
     Full define function.
  */
  virtual void define(/// NODE-centered domain
                      const Box&                         a_bx,
                      /// grid spacing
                      const RealVect&                    a_dx,
                      /// which operator in OperatorType
                      OperatorType                       a_op);


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
  /** Solve with homogeneous Dirichlet boundary conditions.
  */
  void solveHomogeneous(FArrayBox&         a_phi,
                        const FArrayBox&   a_rhs);

  ///
  /** Solve in place with homogeneous Dirichlet boundary conditions.
  */
  void solveHomogeneousInPlace(FArrayBox&   a_phirhs);

  ///
  /** Solve with homogeneous Dirichlet boundary conditions.
  */
  void solveInhomogeneous(FArrayBox&         a_phi,
                          const FArrayBox&   a_rhs);

  ///
  /** Solve in place with homogeneous Dirichlet boundary conditions.
  */
  void solveInhomogeneousInPlace(FArrayBox&   a_phirhs);

  ///
  /** Evaluate operator.
  */
  void evalOperator(FArrayBox&         a_eval,
                    const FArrayBox&   a_phi);

  /*@}*/

protected:

  /** has this PoissonDirichlet been defined yet?
   */
  bool m_isDefined;

  /** NODE-centered grid
   */
  Box m_bx;

  /** mesh spacing at this level
   */
  RealVect m_dx;

  /** verbosity; default 0
   */
  int m_verbose;

  /** type of discrete Laplacian operator
   */
  OperatorType m_op;

private:

  void operator = (const PoissonDirichlet&)
  {
  }

  PoissonDirichlet(const PoissonDirichlet&)
  {
  }
};

#include "NamespaceFooter.H"

#endif
