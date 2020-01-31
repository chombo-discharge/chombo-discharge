#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DFMartin, Sun, May 5, 2002

#ifndef _BICGSTABSMOOTHER_H_
#define _BICGSTABSMOOTHER_H_

#include <iostream>
#include <cmath>
#include "SPACE.H"
#include <cstdlib>

#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "ProblemDomain.H"

#include "BaseBottomSmoother.H"
#include "NamespaceHeader.H"

/// Implements BiCGStab as a bottom smoother for LevelOp-derived classes
/**
    The BiCGStabSmoother class implements the BiCGStab algorithm for
    LevelOp-derived classes
*/
class BiCGStabSmoother: public BaseBottomSmoother
{
public:
  /// Null constructor
  BiCGStabSmoother();

  /// Destructor
  virtual ~BiCGStabSmoother();

  /// "Virtual constructor" workaround.
  /**
      This gets around the "no virtual constructor" rule
   */
  virtual BaseBottomSmoother* new_bottomSmoother() const;

  /// Set maximum number of BiCGStab iterations.
  void setMaxIter(int a_maxIter);

  /// Set solver tolerance
  /**
      Set solver tolerance
   */
  void setSolverTol(Real a_solverTol);

  /// Set whether to enable restart capabillity
  /**
      If this is true, if BiCGStab reaches a hang point, the
      solver will recompute the residual and attempt to restart the
      computation
  */
  void setEnableRestart(bool a_enableRestart);

  /// Set verbosity level
  /**
      Set whether i/o is done when solving (default is false)
  */
  void setVerbose(bool a_verbose);

  /// Do smoothing
  /**
      Performs BiCGStab smoothing of L(a_phi) = a_rhs; L is the operator
      defined in a_levelopPtr, and the preconditioner is
      a_levelopPtr->smooth
  */
  virtual void doBottomSmooth(LevelData<FArrayBox>&       a_phi,
                              const LevelData<FArrayBox>& a_rhs,
                              LevelOp*                    a_levelopPtr);

  /// set metric for judging convergence (if not set, defaults to 1.0)
  virtual void setConvergenceMetric(Real a_metric, int a_comp);

protected:
  /**
     Maximum number of BiCGStab iterations
  */
  int m_maxIter;

  /**
      Solver convergence metrics
  */
  Vector<Real> m_convergenceMetrics;

  /**
     Solver tolerance (we declare the problem solved if the residual is
     smaller than solverTol*convergenceMetric
  */
  Real m_solverTol;

  /**
     Small parameter to check for division by zero
  */
  Real m_small;

  /**
     Small parameter to check for size of correction (to prevent hanging)
  */
  Real m_convergeSmall;

  /**
     If true, enable restarting capability
  */
  bool m_enableRestart;

  /**
     Default is false
  */
  bool m_verbose;
};

#include "NamespaceFooter.H"
#endif
