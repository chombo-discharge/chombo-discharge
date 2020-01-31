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

#ifndef _NOOPSMOOTHER_H_
#define _NOOPSMOOTHER_H_

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

/// Does nothing as a bottom smoother.
/**
    The NoOpSmoother class is a bottom smoother for LevelOp-derived
    classes which does nothing (mostly as a debugging test option)
*/
class NoOpSmoother: public BaseBottomSmoother
{
public:
  /// Null constructor
  NoOpSmoother();

  /// Destructor
  virtual ~NoOpSmoother();

  /// This gets around the "no virtual constructor" rule
  /**
   */
  virtual BaseBottomSmoother* new_bottomSmoother() const;

  /// Do (no) smoothing
  /**
      True to its name, this class does nothing in this functino
  */
  virtual void doBottomSmooth(LevelData<FArrayBox>&       a_phi,
                              const LevelData<FArrayBox>& a_rhs,
                              LevelOp*                    a_levelopPtr);

  // this doesn't really do anything, but is required by BaseBottomSmoother
  virtual void setConvergenceMetric(Real a_metric, int a_comp=0)
  {
  }

  /// this doesn't really do anything, but is required by BaseBottomSmoother
  virtual void setSolverTol(Real a_tol)
  {
  }

};

#include "NamespaceFooter.H"
#endif
