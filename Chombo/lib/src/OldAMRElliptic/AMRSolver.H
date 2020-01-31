#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Tues, July 6, 1999

#ifndef _AMRSOLVER_H_
#define _AMRSOLVER_H_

#include <cstdlib>
#include <cmath>
#include <iostream>
#include "SPACE.H"

#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "AMRLevelMG.H"
#include "LevelSolver.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

/// Class which manages grid hierarchy and composite elliptic solution
/**
   AMRSolver manages the AMR/multigrid solution to the elliptic equation on
   a multiple level grid that satisfies certain proper nesting conditions.
   It can be used either as a solver, or to apply the AMR / multigrid
   V-cycle as a preconditioner for some other iterative method, such as
   a Krylov method.
 */
class AMRSolver
{
public:
  friend class AMRLevelMG;

  /// Creates a AMRSolver whose internal state is  undefined.
  /**
     Need to call define(...) function to use any of
     the functionality of the class
  */
  AMRSolver();

  /// Creates a fully-defined AMRSolver.
  /**
     Calls define function with identical arguments
  */
  AMRSolver(const Vector<DisjointBoxLayout>& a_gridsLevel,
            const Vector<Box>&               a_domainLevel,
            const Vector<Real>&              a_dxLevel,
            const Vector<int>&               a_refRatio,
            int                              a_numLevels,
            int                              a_lBase,
            const LevelOp* const             a_opin,
            int                              a_ncomp = 1,
            bool                             a_limitCoarsening = false);

  /// Creates a fully-defined AMRSolver.
  /**
     Calls define function with identical arguments
  */
  AMRSolver(const Vector<DisjointBoxLayout>& a_gridsLevel,
            const Vector<ProblemDomain>&     a_domainLevel,
            const Vector<Real>&              a_dxLevel,
            const Vector<int>&               a_refRatio,
            int                              a_numLevels,
            int                              a_lBase,
            const LevelOp* const             a_opin,
            int                              a_ncomp = 1,
            bool                             a_limitCoarsening = false);

  ///      Defines AMRSolver object's internal state
  /**
     \b Inputs:
      - \b gridsLevel -- The grids at all levels.
        Each element in the vector is a level's worth of grids.
        gridsLevel[0] are grids for level 0 and so forth.
        Vector index corresponds to level number.
      - \b domainLevel -- The domains at all levels.
        Each element in the vector is a level's domain.
        domainLevel[0] is the domain for level 0 and so forth.
        Vector index corresponds to level number.
      - \b dxLevel -- The grid spacing at all levels.
        Each element in the vector is a level's $\dx$.
        dxLevel[0] -- is $\dx$ for level 0 and so forth.
        Vector index corresponds to level number.
      - \b refRatio -- The refinement ratio between all levels.
        refRatio[0] is the refinement ratio between
        level 0 and level 1;
        Vector index corresponds to level number.
      - \b numlevels -- The number of AMR levels in the calculation.
        The length of the Vector arguments has to be at least numlevels.
      - \b lBase -- coarsest level on which solution is to be computed.
        This needs to be set at the time of definition, in order to
        build the bottom LevelSolver.
      - \b a_opin -- The levelop  to use in the solution.
      - \b a_limitCoarsening: if true, only do multigrid coarsening down
           to next coarser AMR level (only coarsen by a_nRefCrse).  If false,
           coarsen as far as possible. Only relevant when lBase > 0.
  */
  virtual void define(const Vector<DisjointBoxLayout>& a_gridsLevel,
                      const Vector<Box>&                    a_domainLevel,
                      const Vector<Real>&                   a_dxLevel,
                      const Vector<int>&                    a_refRatio,
                      int                                   a_numLevels,
                      int                                   a_lBase,
                      const LevelOp* const                  a_opin,
                      int                                   a_ncomp=1,
                      bool                                  a_limitCoarsening = false);

  /// Defines AMRSolver object's internal state
  /**
     \b Inputs:
      - \b gridsLevel -- The grids at all levels.
        Each element in the vector is a level's worth of grids.
        gridsLevel[0] are grids for level 0 and so forth.
        Vector index corresponds to level number.
      - \b domainLevel -- The domains at all levels.
        Each element in the vector is a level's domain.
        domainLevel[0] is the domain for level 0 and so forth.
        Vector index corresponds to level number.
      - \b dxLevel -- The grid spacing at all levels.
        Each element in the vector is a level's $\dx$.
        dxLevel[0] -- is $\dx$ for level 0 and so forth.
        Vector index corresponds to level number.
      - \b refRatio -- The refinement ratio between all levels.
        refRatio[0] is the refinement ratio between
        level 0 and level 1;
        Vector index corresponds to level number.
      - \b numlevels -- The number of AMR levels in the calculation.
        The length of the Vector arguments has to be at least numlevels.
      - \b lBase -- coarsest level on which solution is to be computed.
        This needs to be set at the time of definition, in order to
        build the bottom LevelSolver.
      - \b a_opin -- The levelop  to use in the solution.
      - \b a_limitCoarsening: if true, only do multigrid coarsening down
           to next coarser AMR level (only coarsen by a_nRefCrse).  If false,
           coarsen as far as possible.  Only relevant when lBase > 0.
  */
  virtual void
  define(const Vector<DisjointBoxLayout>& a_gridsLevel,
         const Vector<ProblemDomain>&     a_domainLevel,
         const Vector<Real>&              a_dxLevel,
         const Vector<int>&               a_refRatio,
         int                              a_numLevels,
         int                              a_lBase,
         const LevelOp* const             a_opin,
         int                              a_ncomp=1,
         bool                             a_limitCoarsening = false);

  /// Destructor
  virtual ~AMRSolver();

  /// Returns true if full define function has been called.
  bool isDefined() const;

  /// Set number of multigrid smoothings on way up v-cycle. Default is 4
  void setNumSmoothUp(int a_numSmoothUp);

  /// Set number of multigrid smoothings on way down v-cycle; Default is 4.
  void setNumSmoothDown(int a_numSmoothDown);

  /// Set tolerance of iterative solution.  Default is 1.0e-10.
  void setTolerance(Real a_tolerance);

  /// returns tolerance for iterative solution
  inline Real Tolerance() const
  {
    return m_tolerance;
  }

  /// Set "operator tolerance" of iterative solution.
  /**
     Iteration will stop if (new_residual/old_residual) > 1-operatorTolerance
     (and at least minIter iterations have been performed). Default is 1.0e-5.
  **/
  void setOperatorTolerance(Real a_operatorTolerance);

  /// Set max number of iterations.  Default is 42.
  void setMaxIter(int a_maxIter);

  /// Set min number of iterations.
  /**
     Only relevant when residual is
     not decreasing fast enough to satisfy the "operator tolerance".
     Default is 4.
  */
  // added by Dan M, 3/24/2000
  void setMinIter(int a_minIter);

  /// Set the number of V-cycles performed at the bottom of the AMR v-cycle
  /**
     Set the number of v-cycles performed from the base level to the maximum
     coarsening of the base level during one overall v-cycle on the AMR
     hierarchy.  Default is 1.
  */
  void setNumVCyclesBottom(int a_numVCycleBottom);

  /// set the convergence metrics (if appropriate) for this problem
  /**
      sets metric against which we define convergence (for example,
      the norm of the RHS).  Also passes a_metric through to the
      LevelMG (and eventually to a bottom solver.)
  */
  void setConvergenceMetric(Real a_metric, int a_comp);

  /// set norm type to use for evaluating convergence
  void setNormType(int a_normType);

  /// Solves on hierarchy to tolerance m_tolerance
  /**
     Solves the elliptic equation over the
     hierarchy of levels lBase ... lMax where $lMax = numlevels-1$.
     If  lBase > 0, then the data at level lBase - 1 is used
     to interpolate boundary conditions at boundary cells that are not
     adjacent to the domain boundary.

     \b Inputs:

     - \b phiLevel - pointers to current guess at the
       solution values for levels (lMin = max(lBase-1,0)) ... lMax.
       Vector index corresponds to level number.
     - \b rhsLevel - pointers to right-hand side for levels lmin ... lMax.
        Vector index corresponds to level number.

     \b Outputs:
     - phiLevel -  LevelData<FArrayBox>s pointed to for levels
       lMin,...,lMax are updated in place.
  */
  void solveAMR(Vector<LevelData<FArrayBox> *>&       a_phiLevel,
                const Vector<LevelData<FArrayBox> *>& a_rhsLevel,
                bool a_initializePhiToZero = true);

  /// Does one relaxation V-cycle using a MG solver.
  /**
     Problem
     is assumed to already be in residual-correction form.

     \b Inputs:
     - \b corrLevel - pointers to current guess at the
       correction values for levels lMin = max(lBase-1,0) ... lMax.
       Vector index corresponds to level number.
     - \b residLevel - pointers to AMR residual for levels lMin ... lMax.
        Vector index corresponds to level number.

     \b Outputs:
     - \b corrLevel - LevelData<FArrayBox>s pointed to for levels
       lMin,...,lMax are updated in place.
  */
  virtual void AMRVCycleMG(Vector<LevelData<FArrayBox> *>&       a_corrLevel,
                           const Vector<LevelData<FArrayBox> *>& a_residLevel);

  /// Calculate norm of multilevel residual on levels lBase to lmax.
  /**
     Does not include data covered by finer levels.
  */
  Vector<Real> computeResidualNorm(Vector<LevelData<FArrayBox> *>&       a_phiLevel,
                                   const Vector<LevelData<FArrayBox> *>& a_rhsLevel,
                                   int                                   a_normType);

  /// Calculate norm of multilevel residual on levels lBase to lmax.
  /**
     Does not include data covered by finer levels. Also returns
     residual in argument a_residLevel.
  */
  Vector<Real> computeResidualNorm(Vector<LevelData<FArrayBox> *>&       a_residLevel,
                                   Vector<LevelData<FArrayBox> *>&       a_phiLevel,
                                   const Vector<LevelData<FArrayBox> *>& a_rhsLevel,
                                   int                                   a_normType);

  /// Calculate multilevel residual on level ilev.
  void computeAMRResidual(Vector<LevelData<FArrayBox> *>&       a_phiLevel,
                          const Vector<LevelData<FArrayBox> *>& a_rhsLevel,
                          LevelData<FArrayBox>&                 a_res,
                          int                                   a_ilev);

  /// Calculate multilevel L(phi).  includes refluxing and all that
  /**
     Includes refluxing and all that. This is the three-level operator.

  */
  void applyAMROperator(Vector<LevelData<FArrayBox> *>& a_phiLevel,
                        LevelData<FArrayBox>&           a_LofPhi,
                        int                             a_ilev);

  /// Calculate multilevel L(phi) with homogeneous physical BCs
  /**
     Calculate multilevel L(phi) with homogeneous physical BCs
     (but with inhomogeneous C/F BCs).  includes refluxing and all that
     This is the three-level operator.
  */
  void applyAMROperatorHphys(Vector<LevelData<FArrayBox> *>& a_phiLevel,
                             LevelData<FArrayBox>&           a_LofPhi,
                             int                             a_ilev);

  /// set whether the solver does i/o while solving. default is true
  void setVerbose(bool a_verbose);

protected:
  /**
      Set data to default values.  not for external use
   */
  virtual void setDefaultValues();

  /**
     Clear internal storage
   */
  virtual void clear();

  /**
      compute multilevel norm on levels lBase to lMax for all
      components in a_phiVect.  Multilevel norm does not include
      covered regions.
  */
  Vector<Real> computeNorm(const Vector<LevelData<FArrayBox>* >& a_phiVect,
                           int a_normType);

  // accessor for  the residual in AMRMGLevel
  LevelData<FArrayBox>& getResid(int a_level);

  // accessor for correction in AmrMGLevel
  LevelData<FArrayBox>& getCorr(int a_level);

  /**
     Error tolerance for error estimation
  */
  Real m_errorTolerance;

  /**
     Solver tolerance
  */
  Real m_tolerance;

  /**
     metrics used to judge convergence
  */
  Vector<Real> m_convergenceMetrics;

  /**
     Tolerance on maxOperator
  */
  Real m_operatorTolerance;

  /**
     Number of levels allowed
  */
  int m_numLevels;

  /**
     Maximum number of working level is m_numLevels - 1
  */
  int m_finestLevel;

  /**
     Refinement ratio.
  */
  Vector<int> m_refRatio;

  /**
     Grid spacing.
  */
  Vector<Real> m_dxLevel;

  /**
     Grids.
  */
  Vector<DisjointBoxLayout> m_gridsLevel;

  /**
     Domains.
  */
  Vector<ProblemDomain> m_domainLevel;

  /**
     Maximum number of solver iterations
  */
  int m_maxIter;

  /**
     Minimum number of solver iterations before "hung converence" criteria
     is valid
  */
  int m_minIter;

  /**
     Number of smoother iterations on up cycle
  */
  int m_numSmoothUp;

  /**
     Number of VCycles for max iter of levelsolver
  */
  int m_numVCyclesBottom;

  /**
     Number of smoother iterations on down cycle
  */
  int m_numSmoothDown;

  /**
      Norm type to use for evaluatiing convergence
  */
  int m_normType;

  /**
     Array of AMRLevelMGs containing data at each refinement level
  */
  Vector<AMRLevelMG*> m_amrmgLevel;

  /**
     Has this solver been defined?
   */
  bool m_isDefined;

  /**
   */
  LevelSolver m_levelSolver;

  /**
   */
  int m_lBase;

  /**
   */
  int m_ncomp;

  /**
     default is true
  */
  bool m_verbose;

private:
  AMRSolver(const AMRSolver&);
  AMRSolver& operator=(const AMRSolver&);
};

#include "NamespaceFooter.H"
#endif
