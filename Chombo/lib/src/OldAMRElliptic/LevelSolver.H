#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Fri, July 23, 1999

#ifndef _LEVELSOLVER_H_
#define _LEVELSOLVER_H_

#include <cstdlib>
#include <iostream>
#include "SPACE.H"
#include <cmath>

#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "LevelMG.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

/// Solves elliptic equations on a level
/**
   LevelSolver Solves elliptic equations on a level
   using multigrid.
*/
class LevelSolver
{
public:
  /// Null constructor, leaves LevelSolver undefined.
  /**
  */
  LevelSolver();

  /// Destructor
  virtual ~LevelSolver();

  /// Full constructor.  Calls full define function.
  /**
      Calls corresponding full define function.

      \b Inputs:
      - \b a_grids: Grid on which solution will be computed.
      - \b a_baseGrids: pointer to next coarser grids.
           Should be == NULL if there is no coarser level.
      - \b a_domain: domain of solution level
      - \b a_dxLevel: grid spacing at solution level
      - \b a_nRefCrse: refinement ratio between two levels.
           Not used if there is no finer level
      - \b a_opin: operator L used in solution of Lphi = rho
      - \b a_limitCoarsening: if true, only do multigrid coarsening down
           to next coarser AMR level (only coarsen by a_nRefCrse).  If false,
           coarsen as far as possible.
  */
  LevelSolver(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout* a_baseGrids,
              const ProblemDomain&     a_domain,
              Real                     a_dxLevel,
              int                      a_nRefCrse,
              const LevelOp* const     a_opin,
              int                      a_ncomp=1,
              bool                     a_limitCoarsening = false);

  /// Full constructor (Deprecated -- use ProblemDomain instead)
  /**
      Calls corresponding full define function.

      \b Inputs:
      - \b a_grids: Grid on which solution will be computed.
      - \b a_baseGrids: pointer to next coarser grids.
           Should be == NULL if there is no coarser level.
      - \b a_domain: domain of solution level
      - \b a_dxLevel: grid spacing at solution level
      - \b a_nRefCrse: refinement ratio between two levels.
           Not used if there is no finer level
      - \b a_opin: operator L used in solution of Lphi = rho
      - \b a_limitCoarsening: if true, only do multigrid coarsening down
           to next coarser AMR level (only coarsen by a_nRefCrse).  If false,
           coarsen as far as possible.
  */
  LevelSolver(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout* a_baseGrids,
              const Box&               a_domain,
              Real                     a_dxLevel,
              int                      a_nRefCrse,
              const LevelOp* const     a_opin,
              int                      a_ncomp = 1,
              bool                     a_limitCoarsening = false);

  /// Full define function.
  /**
      \b Inputs:
      - \b a_grids: Grid on which solution will be computed.
      - \b a_baseGrids: pointer to next coarser grids.
           Should be == NULL if there is no coarser level.
      - \b a_domain: domain of solution level
      - \b a_dxLevel: grid spacing at solution level
      - \b a_nRefCrse: refinement ratio between two levels.
           Not used if there is no finer level
      - \b a_opin: operator L used in solution of Lphi = rho
      - \b a_limitCoarsening: if true, only do multigrid coarsening down
           to next coarser AMR level (only coarsen by a_nRefCrse).  If false,
           coarsen as far as possible.
  */
  virtual void define(const DisjointBoxLayout& a_grids,
                      const DisjointBoxLayout* a_baseGrids,
                      const ProblemDomain&     a_domain,
                      Real                     a_dxLevel,
                      int                      a_nRefCrse,
                      const LevelOp* const     a_opin,
                      int                      a_ncomp = 1,
                      bool                     a_limitCoarsening = false);

  /// Full define function. (Deprecated -- use ProblemDomain version instead)
  /**
      \b Inputs:
      - \b a_grids: Grid on which solution will be computed.
      - \b a_baseGrids: pointer to next coarser grids.
           Should be == NULL if there is no coarser level.
      - \b a_domain: domain of solution level
      - \b a_dxLevel: grid spacing at solution level
      - \b a_nRefCrse: refinement ratio between two levels.
           Not used if there is no finer level
      - \b a_opin: operator L used in solution of Lphi = rho
      - \b a_limitCoarsening: if true, only do multigrid coarsening down
           to next coarser AMR level (only coarsen by a_nRefCrse).  If false,
           coarsen as far as possible.
  */
  virtual void define(const DisjointBoxLayout& a_grids,
                      const DisjointBoxLayout* a_baseGrids,
                      const Box&               a_domain,
                      Real                     a_dxLevel,
                      int                      a_nRefCrse,
                      const LevelOp* const     a_opin,
                      int                      a_ncomp=1,
                      bool                     a_limitCoarsening = false);

  /// Has this object been defined already?
  bool isDefined() const;

  /// Clears storage
  void clearMemory();

  /// Internal use function to set default values for solver parameters
  void setDefaultValues();

  /// Do a level solve with homogeneous BCs on this level using multigrid
  /**
      Uses homogeneous bcs at coarse/fine interfaces
  */
  void levelSolveH(/// initial guess at solution, returns solution
                   LevelData<FArrayBox>&       a_phi,
                   /// right-hand side
                   const LevelData<FArrayBox>& a_rhs,
                   /// if true, set phi to 0 at start (o/w use initial guess)
                   bool                        a_initializePhiToZero=true);

  /// Do a level solve with inhomogeneous BCs on this level using multigrid
  /**
      Uses inhomogeneous bcs at interfaces.  phic should be set to NULL
      if there is no coarser level.
  */
  void levelSolve(/// initial guess at solution, returns solution
                  LevelData<FArrayBox>&       a_phi,
                  /// solution at next coarser level
                  const LevelData<FArrayBox>* a_phic,
                  /// right-hand side
                  const LevelData<FArrayBox>& a_rhs,
                  /// if true, set phi to 0 at start (o/w use initial guess)
                  bool                        a_initializePhiToZero=true);

  /// Set number of multigrid smoothings on way up v-cycle. Default is 4
  void setNumSmoothUp(int a_numSmoothUp);

  /// Set number of multigrid smoothings on way down v-cycle; Default is 4.
  void setNumSmoothDown(int a_numSmoothDown);

  /// Set solution tolerance
  void setTolerance(Real a_tolerance);

  /// Set "operator tolerance" of iterative solution.
  /**
      Iteration will stop if
      (new_residual/old_residual) > 1-operatorTolerance (and
      at least minIter iterations have been performed).
      Default is 1.0e-5.
  **/
  void setOperatorTolerance(Real a_operatorTolerance);

  /// Set how much info is printed out during solve
  void setVerbose(bool a_verbose);

  /// Set max number multigrid vcycles to use in solution
  void setMaxIter(int a_maxIter);

  /// Set min number of iterations.
  /**
      Only relevant when residual is
      not decreasing fast enough to satisfy the "operator tolerance".
      Default is 4.
  */
  // added by Dan M, (9/20/02)
  void setMinIter(int a_minIter);

  /// set the convergence metrics (if appropriate) for this problem
  /**
      sets metric against which we define convergence (for example,
      the norm of the RHS).  Also passes a_metric through to the
      LevelMG (and eventually to a bottom solver.)
  */
  virtual void setConvergenceMetric(Real a_metric, int a_comp);

  /// rest the convergence metrics to an undefined state (for solver re-use)
  virtual void resetConvergenceMetrics();

  /// set type of norm to use to judge convergence
  void setNormType(int a_normType);

protected:
  /**
     LevelMG level solver object to relax on this level
  */
  LevelMG m_levelMG;

  /**
     LevelMG level solver object to relax on this level
  */
  LevelOp* m_levelOpPtr;

  /**
     Has this LevelSolver been defined yet?
  */
  bool m_isDefined;

  /**
     Grids at this level
  */
  DisjointBoxLayout  m_grids;

  /**
     Problem domain at this level
  */
  ProblemDomain m_domain;

  /**
     Residual
  */
  LevelData<FArrayBox> m_resid;

  /**
     Useful scratch space.
  */
  LevelData<FArrayBox> m_scratch;

  /**
     Correction
  */
  LevelData<FArrayBox> m_corr;

  /**
     Refinement ratio to coarser level
   */
  int m_nRefCrse;

  /**
     Cell spacing on this level
  */
  Real m_dxLevel;

  /**
     Maximum number of solver iterations
  */
  int m_maxIter;

  /**
     Minimum number of iterations
  */
  int m_minIter;

  /**
     Number of smoother iterations on up cycle
  */
  int m_numSmoothUp;

  /**
     Number of smoother iterations on down cycle
  */
  int m_numSmoothDown;

  /**
     Solver tolerance
  */
  Real m_tolerance;

  /**
     type of norm to use for evaluating convergence
  */
  int m_normType;

  /**
     metrics used to judge convergence
  */
  Vector<Real> m_convergenceMetrics;

  /**
     Tolerance used to stop solver iteration if convergence hangs.
  */
  Real m_operatorTolerance;

  /**
     How verbose shold the solve be?
  */
  bool m_verbose;

private:
  void operator = (const LevelSolver&)
  {
  }

  LevelSolver(const LevelSolver&)
  {
  }
};

#include "NamespaceFooter.H"
#endif
