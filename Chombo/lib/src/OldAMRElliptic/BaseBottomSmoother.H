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

#ifndef _BASEBOTTOMSMOOTHER_H_
#define _BASEBOTTOMSMOOTHER_H_

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
#include "LevelOp.H"
#include "NamespaceHeader.H"

/// Pure base class to encapsulate bottom smooth for LevelOp-derived objects.
/** The BaseBottomSmoother class is a virtual base class for implementing
    different bottomSmoother algorithms for LevelOp-derived classes. It
    takes a LevelOp pointer as an argument during its definition, in order
    to access the problem-dependent (equation-dependent) functionality of
    the specific LevelOp derivative class.
*/
class BaseBottomSmoother
{
public:
  /// Null constructor
  BaseBottomSmoother()
  {
  }

  /// Destructor
  virtual ~BaseBottomSmoother()
  {
  }

  /// This gets around the "no virtual constructor" rule
  virtual BaseBottomSmoother* new_bottomSmoother() const = 0;

  /// Do smoothing
  virtual void doBottomSmooth(LevelData<FArrayBox>&       a_phi,
                              const LevelData<FArrayBox>& a_rhs,
                              LevelOp*                    a_levelopPtr) = 0;

  /// sets metric to judge convergence (normally norm(rhs) of outer problem
  /** To make the iterative solver convergence more consistent, the solver
      will consider itself "solved" if the residual is reduced to
      epsilon*metric, where epsilon is the solver tolerance (set by
      setSolverTol).  A good starting point for the metric is the norm of
      the RHS of the original problem (i.e. the rhs sent in to LevelSolver
      or AMRSolver).  Default value is 1.0.
   */
  virtual void setConvergenceMetric(Real a_metric, int a_comp=0) = 0;

  /// set solver tolerance
  virtual void setSolverTol(Real a_tol) = 0;

};

#include "NamespaceFooter.H"
#endif
