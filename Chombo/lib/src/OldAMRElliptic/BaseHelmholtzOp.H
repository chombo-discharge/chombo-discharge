#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DMartin, Weds, July 21, 1999

#ifndef _BASEHELMHOLTZOP_H_
#define _BASEHELMHOLTZOP_H_

#include "LevelOp.H"
#include "GhostBC.H"
#include "NamespaceHeader.H"

/// Pure-virtual base Helmholtz operator class (alpha*I + div*beta*grad)
/** Virtual base class to solve Helmholtz-type problems.  The raison
    d'Etre for this class is to allow for scaling of Alpha and Beta
    for any Helmholtz-type operator
 */
class BaseHelmholtzOp: public LevelOp
{
public:
  /// Scale the helmholtz coeficient alpha in (alpha*I + beta*laplacian)phi
  /**
     (alpha := a_scale*alpha)
  */
  virtual void scaleAlpha(Real a_scale) = 0;

  /// Scale the helmholtz coefficient beta in (alpha*I + beta*laplacian)phi
  /**
     (beta := a_scale*beta)
  */
  virtual void scaleBeta(Real a_scale) =0;

  /// set (or reset) the physical boundary conditions
  /**
   */
  virtual void setDomainGhostBC(const DomainGhostBC& a_dombc) = 0;

  /// Null Constructor
  BaseHelmholtzOp()
  {
  }

  /// Destructor
  virtual ~BaseHelmholtzOp()
  {
  }
};

#include "NamespaceFooter.H"
#endif
