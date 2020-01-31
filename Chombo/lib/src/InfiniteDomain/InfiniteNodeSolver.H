#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// InfiniteNodeSolver.H
// petermc, 26 Jan 2004

#ifndef _INFINITENODESOLVER_H_
#define _INFINITENODESOLVER_H_

#include <cstdlib>
#include <iostream>
#include "SPACE.H"
#include <cmath>
#include "REAL.H"
#include "FArrayBox.H"
#include "RealVect.H"
#include "PoissonDirichlet.H"
#include "Normals.H"
#include "Multipoles.H"

#include "NamespaceHeader.H"

/// Solves elliptic equations on a level.
class InfiniteNodeSolver
/** InfiniteNodeSolver solves elliptic equations on a level using multigrid.
*/
{

public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor leaves InfiniteNodeSolver undefined.
  */
  InfiniteNodeSolver();

  ///
  /**
     Full constructor.
     Calls full define() function with identical arguments.
  */
  InfiniteNodeSolver(/// NODE-centered domain
                     const Box&        a_domain,
                     /// size of buffer from rhs domain to D1
                     int               a_s1,
                     /// size of buffer from D1 to D2, the solution domain
                     int               a_s2,
                     /// patch size on faces of D1
                     int               a_patchSize,
                     /// coarsening ratio on faces of D2, normally same as a_patchSize
                     int               a_dstFaceCoarsening,
                     /// mesh spacing
                     const RealVect&   a_dx,
                     /// which operator in OperatorType
                     PoissonDirichlet::OperatorType      a_op,
                     /// degree of multipole expansion
                     int               a_multipoleOrder,
                     /// degree of normal derivative expansion
                     int               a_degreeNormalDerivative,
                     /// number of extra points to use in polynomial interpolation (even)
                     int               a_interpBorder,
                     /// run in parallel?
                     bool              a_parallel = false,
                     /// solve for coarser-level solution also?
                     bool              a_getOuterCoarse = false,
                     /// refinement ratio of coarse and fine solutions
                     int               a_refToCoarse = 1,
                     /// radius of additional points in coarse solution beyond fine solution
                     int               a_coarseAddRadius = 0);


  ///
  /**
     Abbreviated constructor.
     Calls abbreviated define() function with identical arguments.
  */
  InfiniteNodeSolver(/// NODE-centered domain
                     const Box&        a_domain,
                     /// mesh spacing
                     const RealVect&   a_dx,
                     /// order of accuracy, 2 for O(h^2) or 4 for O(h^4)
                     int               a_order,
                     /// run in parallel?
                     bool              a_parallel = false,
                     /// solve for coarser-level solution also?
                     bool              a_getOuterCoarse = false,
                     /// refinement ratio of coarse and fine solutions
                     int               a_refToCoarse = 1,
                     /// radius of additional points in coarse solution beyond fine solution
                     int               a_coarseAddRadius = 0);

  ///
  ///
  /** Destructor.
   */
  virtual ~InfiniteNodeSolver();

  ///
  /**
     Full define function.
  */
  virtual void define(/// NODE-centered domain
                     const Box&        a_domain,
                     /// size of buffer from rhs domain to D1
                     int               a_s1,
                     /// size of buffer from D1 to D2, the solution domain
                     int               a_s2,
                     /// patch size and also coarsening ratio on D2 faces
                     int               a_patchSize,
                     /// coarsening ratio on faces of D2, normally same as a_patchSize
                     int               a_dstFaceCoarsening,
                     /// mesh spacing
                     const RealVect&   a_dx,
                     /// which operator in OperatorType
                     PoissonDirichlet::OperatorType      a_op,
                     /// degree of multipole expansion
                     int               a_multipoleOrder,
                     /// degree of normal derivative expansion
                     int               a_degreeNormalDerivative,
                     /// number of extra points to use in polynomial interpolation (even)
                     int               a_interpBorder,
                     /// run in parallel?
                     bool              a_parallel = false,
                     /// solve for coarser-level solution also?
                     bool              a_getOuterCoarse = false,
                     /// refinement ratio of coarse and fine solutions
                     int               a_refToCoarse = 1,
                     /// radius of additional points in coarse solution beyond fine solution
                     int               a_coarseAddRadius = 0);


  ///
  /**
     Abbreviated define function.
  */
  virtual void define(/// NODE-centered domain
                     const Box&        a_domain,
                     /// mesh spacing
                     const RealVect&   a_dx,
                     /// order of accuracy, 2 for O(h^2) or 4 for O(h^4)
                     int               a_order,
                     /// run in parallel?
                     bool              a_parallel = false,
                     /// solve for coarser-level solution also?
                     bool              a_getOuterCoarse = false,
                     /// refinement ratio of coarse and fine solutions
                     int               a_refToCoarse = 1,
                     /// radius of additional points in coarse solution beyond fine solution
                     int               a_coarseAddRadius = 0);

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

  /** Delete internal memory and set pointers to NULL.
   */
  void clearMemory();

  /** Set default values.  Does not deal with memory.
   */
  void setDefaultValues();

  ///
  /** Solve for a_phi on m_dstBox.
  */
  void solve(FArrayBox&         a_phi,
             const FArrayBox&   a_rhs);

  ///
  /** Solve for a_phi on m_dstBox,
      and also for a_phiCoarse.
  */
  void solveWithCoarse(FArrayBox&         a_phi,
                       FArrayBox&         a_phiCoarse,
                       const FArrayBox&   a_rhs);

  ///
  /** Solve for a_phi on m_dstBox,
      and if m_getOuterCoarse then also for a_phiCoarse.
  */
  void solveUnified(FArrayBox&         a_phi,
                    FArrayBox*         a_phiCoarsePtr,
                    const FArrayBox&   a_rhs);

  ///
  /** Begin timers for solve.
  */
  void beginTimers();

  ///
  /** End timers for solve.
  */
  void endTimers();

  /*@}*/

protected:

  /** has this InfiniteNodeSolver been defined yet?
   */
  bool m_isDefined;

  /** solver for initial Poisson-Dirichlet solve
   */
  PoissonDirichlet m_solverPoissonSrc;

  /** solver for final Poisson-Dirichlet solve
   */
  PoissonDirichlet m_solverPoissonDst;

  /** for finding normal derivatives
   */
  Normals m_normals;

  /** for computing multipole coefficients and evaluating multipole expansions
   */
  Multipoles m_multipoles;

  /** run in parallel?
   */
  bool m_parallel;

  /** normal derivatives on each face of D1 (m_srcBox):  length 2*SpaceDim
   */
  Vector<FArrayBox*> m_derivs;

  /** NODE-centered grid on which rhs lives
   */
  Box m_domain;

  /** mesh spacing at this level
   */
  RealVect m_dx;

  /** verbosity; default 0
   */
  int m_verbose;

  /** type of discrete Laplacian operator
   */
  PoissonDirichlet::OperatorType m_op;

  /** NODE-centered box D1 on which charge sits
   */
  Box m_srcBox;

  /** NODE-centered box D2 on which potential is evaluated
   */
  Box m_dstBox;

  /** fab on m_srcBox (D1)
   */
  FArrayBox m_phirhsSrcFab;

  /** fab on m_dstBox (D2)
   */
  FArrayBox m_phirhsDstFab;

  /** size of buffer from m_domain to m_srcBox (D1)
   */
  int m_s1;

  /** size of buffer from m_srcBox (D1) to m_dstBox (D2)
   */
  int m_s2;

  /** length of each patch on faces of m_srcBox (D1)
   */
  int m_patchSize;

  /**  coarsening ratio on faces of m_dstBox (D2)
   */
  int m_dstFaceCoarsening;

  /** degree of multipole expansions in m_multipoles
   */
  int m_multipoleOrder;

  /** degree of normal derivative approximation in m_normals
   */
  int m_degreeNormalDerivative;

  /** width of layer for polynomial interpolation on faces of D2 in m_multipoles
   */
  int m_interpBorder;

  /** solve for coarser-level solution also?
   */
  bool m_getOuterCoarse;

  /** refinement ratio to next coarser level
   */
  int m_refToCoarse;

  /** radius of additional points in coarse solution beyond fine solution
   */
  int m_coarseAddRadius;

  /** number of uses of this solver
   */
  int m_uses;

private:

  void operator = (const InfiniteNodeSolver&)
  {
  }

  InfiniteNodeSolver(const InfiniteNodeSolver&)
  {
  }
};

#include "NamespaceFooter.H"

#endif
