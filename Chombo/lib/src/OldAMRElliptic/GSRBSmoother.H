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

#ifndef _GSRBSMOOTHER_H_
#define _GSRBSMOOTHER_H_

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

/// GSRB smoothing as a bottom smoother.
/** The GSRBSmoother class is a bottom smoother for LevelOp-derived
    classes which does nothing but a whole lotta GSRB smoothing
*/
class GSRBSmoother: public BaseBottomSmoother
{
public:
  /// Null constructor
  GSRBSmoother();

  /// Destructor
  virtual ~GSRBSmoother();

  ///
  /**
      This gets around the "no virtual constructor" rule
   */
  virtual BaseBottomSmoother* new_bottomSmoother() const;

  /// Do GSRB smoothing
  /**
      True to its name, this class does GSRB smoothing in this function
  */
  virtual void doBottomSmooth(LevelData<FArrayBox>&       a_phi,
                              const LevelData<FArrayBox>& a_rhs,
                              LevelOp*                    a_levelopPtr);

  /// Set the number of smoother passes to perform
  virtual void setNumSmooth(int a_numSmooth);

  /// Returns the number of smoother passes which are performed
  virtual int numSmooth() const;

  // this doesn't really do anything, but is required by BaseBottomSmoother
  virtual void setConvergenceMetric(Real a_metric, int a_comp=0)
  {
  }

  /// this doesn't really do anything, but is required by BaseBottomSmoother
  virtual void setSolverTol(Real a_tol)
  {
  }

protected:
  /// Number of GSRB smoother passes to do.
  int m_numSmooth;
};

#include "NamespaceFooter.H"
#endif
