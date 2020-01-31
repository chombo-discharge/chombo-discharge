#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// InfiniteMLCSolver.H
// petermc, 11 Jun 2004

#ifndef _INFINITEMLCSOLVER_H_
#define _INFINITEMLCSOLVER_H_

#include <cstdlib>
#include <iostream>
#include "SPACE.H"
#include <cmath>
#include "REAL.H"
#include "NodeFArrayBox.H"
#include "InfiniteNodeSolver.H"
#include "LevelData.H"

#include "NamespaceHeader.H"

/// Solves elliptic equations on a decomposed domain.
class InfiniteMLCSolver
/** InfiniteMLCSolver solves Poisson's equation L(phi) = rhs
    with infinite-domain boundary conditions, on a domain split into patches.
*/
{

public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor leaves InfiniteMLCSolver undefined.
  */
  InfiniteMLCSolver();

  ///
  /**
     Full constructor.
     Calls full define() function with identical arguments.
  */
  InfiniteMLCSolver(/// CELL-centered base box in all grids
                    const Box&                 a_fineBaseBox,
                    /// CELL-centered physical domain
                    const Box&                 a_fineDomain,
                    /// refinement ratio to level for coarse global solve
                    int                        a_refToCoarse,
                    /// mesh spacing
                    const RealVect&            a_dx,
                    /// size of buffer from local domain to local D1
                    int                        a_s1Local,
                    /// size of buffer from local D1 to local D2
                    int                        a_s2Local,
                    /// patch size on faces of local D1
                    int                        a_patchSizeLocal,
                    /// coarsening ratio on faces of local D2, normally same as a_patchSizeLocal
                    int                        a_dstFaceCoarseningLocal,
                    /// operator for initial local solve
                    PoissonDirichlet::OperatorType a_opLocalInitial,
                    /// degree of multipole expansion for initial local solve
                    int                        a_multipoleOrderLocal,
                    /// degree of normal derivative expansion for initial local solve
                    int                        a_degreeNormalDerivativeLocal,
                    /// width of layer to use in polynomial interpolation in initial local solve
                    int                        a_interpBorderLocal,
                    /// size of buffer from global domain to global D1
                    int                        a_s1Global,
                    /// size of buffer from global D1 to global D2
                    int                        a_s2Global,
                    /// patch size on faces of global D1
                    int                        a_patchSizeGlobal,
                    /// coarsening ratio on faces of global D2, normally same as a_patchSizeGlobal
                    int                        a_dstFaceCoarseningGlobal,
                    /// operator for global coarse solve
                    PoissonDirichlet::OperatorType a_opGlobal,
                    /// degree of multipole expansion for global coarse solve
                    int                        a_multipoleOrderGlobal,
                    /// degree of normal derivative expansion for global coarse solve
                    int                        a_degreeNormalDerivativeGlobal,
                    /// width of layer to use in polynomial interpolation in coarse global solve
                    int                        a_interpBorderGlobal,
                    /// width of extra layer to use in coarse/fine polynomial interpolation
                    int                        a_interpBorderCF,
                    /// minimum width, in coarse cells, of buffer between local grids
                    int                        a_minBufferCoarse,
                    /// operator for final local solve
                    PoissonDirichlet::OperatorType a_opLocalFinal);

  ///
  /**
     Abbreviated constructor.
     Calls abbreviated define() function with identical arguments.
  */
  InfiniteMLCSolver(/// CELL-centered base box in all grids
                    const Box&                 a_fineBaseBox,
                    /// CELL-centered physical domain
                    const Box&                 a_fineDomain,
                    /// refinement ratio to level for coarse global solve
                    int                        a_refToCoarse,
                    /// mesh spacing
                    const RealVect&            a_dx);

  ///
  /** Destructor.
   */
  virtual ~InfiniteMLCSolver();

  ///
  /**
     Full define function.
  */
  virtual void define(/// CELL-centered base box in all grids
                      const Box&                 a_fineBaseBox,
                      /// CELL-centered physical domain
                      const Box&                 a_fineDomain,
                      /// refinement ratio to level for coarse global solve
                      int                        a_refToCoarse,
                      /// mesh spacing
                      const RealVect&            a_dx,
                      /// size of buffer from local domain to local D1
                      int                        a_s1Local,
                      /// size of buffer from local D1 to local D2
                      int                        a_s2Local,
                      /// patch size on faces of local D1
                      int                        a_patchSizeLocal,
                      /// coarsening ratio on faces of local D2, normally same as a_patchSizeLocal
                      int                        a_dstFaceCoarseningLocal,
                      /// operator for initial local solve
                      PoissonDirichlet::OperatorType a_opLocalInitial,
                      /// degree of multipole expansion for initial local solve
                      int                        a_multipoleOrderLocal,
                      /// degree of normal derivative expansion for initial local solve
                      int                        a_degreeNormalDerivativeLocal,
                      /// width of layer to use in polynomial interpolation in initial local solve
                      int                        a_interpBorderLocal,
                      /// size of buffer from global domain to global D1
                      int                        a_s1Global,
                      /// size of buffer from global D1 to global D2
                      int                        a_s2Global,
                      /// patch size on faces of global D1
                      int                        a_patchSizeGlobal,
                      /// coarsening ratio on faces of global D2, normally same as a_patchSizeLocal
                      int                        a_dstFaceCoarseningGlobal,
                      /// operator for global coarse solve
                      PoissonDirichlet::OperatorType a_opGlobal,
                      /// degree of multipole expansion for global coarse solve
                      int                        a_multipoleOrderGlobal,
                      /// degree of normal derivative expansion for global coarse solve
                      int                        a_degreeNormalDerivativeGlobal,
                      /// width of layer to use in polynomial interpolation in coarse global solve
                      int                        a_interpBorderGlobal,
                      /// width of layer to use in coarse/fine polynomial interpolation
                      int                        a_interpBorderCF,
                      /// minimum width, in coarse cells, of buffer between local grids
                      int                        a_minBufferCoarse,
                      /// operator for final local solve
                      PoissonDirichlet::OperatorType a_opLocalFinal);

  ///
  /**
     Abbreviated define function.
  */
  virtual void define(/// CELL-centered base box in all grids
                      const Box&                 a_fineBaseBox,
                      /// CELL-centered physical domain
                      const Box&                 a_fineDomain,
                      /// refinement ratio to level for coarse global solve
                      int                        a_refToCoarse,
                      /// mesh spacing
                      const RealVect&            a_dx);

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
  /** Solve.
  */
  void solve(/// solution, lives on fineGrids, no ghosts
             LevelData<NodeFArrayBox>&    a_phi,
             /// right-hand side, lives on fineGrids, no ghosts.  is not changed, but isn't const because it must be shifted.
             LevelData<NodeFArrayBox>&    a_rhs);

  /*@}*/

protected:

  /** has this InfiniteNodeSolver been defined yet?
   */
  bool m_isDefined;

  /** verbosity; default 0
   */
  int m_verbose;

  /** size of buffer from a box in m_fineGrids to corresponding local D1
   */
  int m_s1Local;

  /** size of buffer from local D1 to local D2
   */
  int m_s2Local;

  /** length of each patch on faces of local D1
   */
  int m_patchSizeLocal;

  /** coarsening ratio on faces of local D2
   */
  int m_dstFaceCoarseningLocal;

  /** size of buffer from m_coarseGlobalDomain to global D1
   */
  int m_s1Global;

  /** size of buffer from global D1 to global D2
   */
  int m_s2Global;

  /** length of each patch on faces of global D1
   */
  int m_patchSizeGlobal;

  /** coarsening ratio on faces of global D2
   */
  int m_dstFaceCoarseningGlobal;

  /** CELL-centered base box of which all fine grids are translations;
      its lower corner is at zero.
   */
  Box m_fineBaseBox;

  /** CELL-centered box covering all of the domains of the fine grids
   */
  Box m_fineDomain;

  /** refinement ratio between this level and the coarser level
   */
  int m_refToCoarse;

  /** mesh spacing at this level
   */
  RealVect m_dx;

  /** mesh spacing at global coarse level
   */
  RealVect m_dxCoarse;

  /** Poisson-Dirichlet solver for final local solve
   */
  PoissonDirichlet m_poisson;

  /** operator for initial local solve:  should be Mehrstellen
   */
  PoissonDirichlet::OperatorType m_opLocalInitial;

  /** degree of multipole expansions for initial local solve
   */
  int m_multipoleOrderLocal;

  /** degree of normal derivative expansion for initial local solve
   */
  int m_degreeNormalDerivativeLocal;

  /** width of layer to use in polynomial interpolation in initial local solve
   */
  int m_interpBorderLocal;

  /** operator for coarse global solve:  should be Mehrstellen
   */
  PoissonDirichlet::OperatorType m_opGlobal;

  /** degree of multipole expansions for coarse global solve
   */
  int m_multipoleOrderGlobal;

  /** degree of normal derivative expansion for coarse global solve
   */
  int m_degreeNormalDerivativeGlobal;

  /** width of layer to use in polynomial interpolation in coarse global solve
   */
  int m_interpBorderGlobal;

  /** width of extra layer to use in coarse/fine polynomial interpolation
   */
  int m_interpLayer;

  /** minimum width, in coarse cells, of buffer between local grids
   */
  int m_minBufferCoarse;

  /** operator for final local Poisson-Dirichlet solve
   */
  PoissonDirichlet::OperatorType m_opLocalFinal;

  /** class for initial local infinite-domain solve
   */
  InfiniteNodeSolver m_solverLocalInitial;

  /** class for global infinite-domain solve
   */
  InfiniteNodeSolver m_solverGlobal;

  /** width of layer of additional points in coarse local solution beyond
      fine local solution;
      equal to m_interpLayer + m_bufferLocal / m_refToCoarse
   */
  int m_coarseAddRadiusLocal;

  /** size of buffer for m_phiLocal;
      equal to m_minBufferCoarse * m_refToCoarse.
   */
  int m_bufferLocal;

  /** size of buffer beyond m_fineGrids on which solution is computed;
      equal to m_s1Local + m_s2Local.
  */
  int m_bufferLocalSoln;

  /** CELL-centered domain for coarse global solve;
      global residual lives on its surrounding nodes
   */
  Box m_coarseGlobalDomain;

  /** CELL-centered domain on surrounding nodes of which
      coarse global solution lives;
      it is m_coarseGlobalDomain grown by m_s1Global + m_s2Global
      in all directions.
   */
  Box m_coarseGlobalDomainPhi;

  /** coefficients for interpolating from coarse global solution
      to local solution
   */
  Real* m_interpCoeffs;

  /** CELL-centered fine local grids m_fineGrids grown by m_bufferLocal;
      m_phiLocal lives on these.
   */
  // BoxLayout m_fineGridsBuffer;

  /** CELL-centered fine local grids m_fineGrids coarsened by m_refToCoarse
      and then grown by m_bufferLocal / m_refToCoarse - 1;
      m_residLocalCoarse lives on these.
   */
  // BoxLayout m_coarseLocalGridsBufferm1;

  /** CELL-centered fine local grids m_fineGrids coarsened by m_refToCoarse
      and then grown by m_coarseAddRadiusLocal;
      m_phiLocalCoarse lives on these.
   */
  // BoxLayout m_coarseLocalGridsBufferInterp;

  /** CELL-centered fine local grids m_fineGrids coarsened by m_refToCoarse
      and then grown by m_interpLayer;
      intCoarseData lives on these.
   */
  // BoxLayout m_coarseLocalGridsInterp;

  /** initial local solution,;
      lives on nodes surrounding boxes of m_fineGridsBuffer
   */
  // BoxLayoutData<NodeFArrayBox> m_phiLocal;

  /** sampling of initial local solution;
      lives on nodes surrounding boxes of m_coarseLocalGridsBufferInterp
   */
  // BoxLayoutData<NodeFArrayBox> m_phiLocalCoarse;

  /** residual of sampled initial local solution;
      lives on nodes surrounding boxes of m_coarseLocalGridsBufferm1
   */
  // BoxLayoutData<NodeFArrayBox> m_residLocalCoarse;

private:

  void operator = (const InfiniteMLCSolver&)
  {
  }

  InfiniteMLCSolver(const InfiniteMLCSolver&)
  {
  }

  ///
  /** Return sum of a_data in a_sum.
   */
  void plusReduce(NodeFArrayBox&                                a_sum,
                  const Vector< RefCountedPtr<NodeFArrayBox> >& a_data);

  ///
  /** Along each edge of the box of a_fab that is not an edge of a_bxFace,
      halve the values of a_fab.
   */
  void halveInternalEdges(FArrayBox&  a_fab,
                          const Box&  a_bxFace);
};

#include "NamespaceFooter.H"

#endif
