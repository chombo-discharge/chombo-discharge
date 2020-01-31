#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Mon, July 26, 1999

#ifndef _LEVELOP_H_
#define _LEVELOP_H_

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
#include "LoHiSide.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

/// Pure base class to encapsulate level operations API for elliptic solvers.
/**
    LevelOp is a pure base class to
    encapsulate level operations API for elliptic solvers.
    All elliptic operators used by AMRSolver inherit the LevelOp interface.
*/
class LevelOp
{
public:
  /// Null constructor
  LevelOp()
  {
  }

  /// Destructor
  virtual ~LevelOp()
  {
  }

  /// Full define function.
  /**
      Makes all coarse-fine
      information and sets internal variables
  */
  virtual void define(/// Grids at this level
                      const DisjointBoxLayout& a_grids,
                      /// Coarser-level grids
                      const DisjointBoxLayout* a_baseBAPtr,
                      /// Cell spacing at this level
                      Real                     a_dxLevel,
                      /// Refinement ratio btwn this level and coarser level
                      int                      a_refRatio,
                      /// Physical problem domain
                      const ProblemDomain&     a_domain,
                      /// Will we only need homogeneous C/F BC's?
                      bool                     a_homogeneousOnly,
                      /// Number of components
                      int                      a_ncomp = 1) = 0;

  /// Full define function. (Deprecated -- use the ProblemDomain version)
  /**
      Makes all coarse-fine
      information and sets internal variables
  */
  virtual void define(const DisjointBoxLayout& a_grids,
                      const DisjointBoxLayout* a_baseBAPtr,
                      Real                     a_dxLevel,
                      int                      a_refRatio,
                      const Box&               a_domain,
                      bool                     a_homogeneousOnly,
                      int                      a_ncomp = 1) = 0;

  /// Full define function.
  /**
      Define from finer levelop.
      reftofine is the refinement ratio
      to get from this LevelOp to \em opfine
  */
  virtual void define(/// Finer levelOp from which this will be defined
                      const LevelOp* a_opfine,
                      /// Amount this LevelOp will be coarser than \em opfine
                      int            a_reftoFine) = 0;

  /// Coarse / Fine interpolation operator.
  /**
      Fill ghost cells of phi by interpolating using coarse-level and
      fine-level data.
  */
  virtual void CFInterp(/// Data at this level
                        LevelData<FArrayBox>&       a_phif,
                        /// Coarser-level data
                        const LevelData<FArrayBox>& a_phic) = 0;

  /// Homogeneous coarse/fine interpolation operator
  /**
  */
  virtual void homogeneousCFInterp(LevelData<FArrayBox>& a_phif) = 0;

  /// "Virtual constructor" workaround"
  /**
  */
  virtual LevelOp* new_levelop() const = 0;

  /// Has this object been fully defined?
  virtual bool isDefined() const = 0;

  /// Smoother.
  /**
     Assumes that problem has already been put in
     residual correction form,
     so that C/F boundary conditions are homogeneous.
  */
  virtual void smooth(/// Current guess at solution
                      LevelData<FArrayBox>&       a_phi,
                      /// Right-hand side of problem we're solving
                      const LevelData<FArrayBox>& a_rhs) = 0;

  /// Apply preconditioning to the solution on this level.
  /**
      Given rhshat, returns phihat from M(phihat) = rhshat
  */
  virtual void levelPreconditioner(LevelData<FArrayBox>&       a_phihat,
                                   const LevelData<FArrayBox>& a_rhshat) = 0;

  /// Evaluate operator, inhomogeneous boundary conditions
  /**
       Uses inhomogeneous form of both C/F and physical boundary
       conditions.
  **/
  virtual void applyOpI(LevelData<FArrayBox>&       a_phi,
                        const LevelData<FArrayBox>* a_phicPtr,
                        LevelData<FArrayBox>&       a_lofPhi) = 0;

  /// Evaluate operator
  /**
       Uses inhomogeneous form of C/F boundary conditions.
       homogeneous form of phys boundary conditions.
  **/
  virtual void applyOpIcfHphys(LevelData<FArrayBox>&       a_phi,
                               const LevelData<FArrayBox>* a_phicPtr,
                               LevelData<FArrayBox>&       a_lofPhi)  = 0;

  /// Evaluate operator
  /**
       Uses homogeneous form of C/F boundary conditions
       and inhomogeneous form of physical boundary conditions.
  */
  virtual void applyOpHcfIphys(LevelData<FArrayBox>& a_phi,
                               LevelData<FArrayBox>& a_lofPhi) = 0;

  /// Evaluate Operator with homogeneous boundary conditions
  /**
       Uses homogeneous form of both  C/F and physical boundary
       conditions.
  */
  virtual void applyOpH(LevelData<FArrayBox>& a_phi,
                        LevelData<FArrayBox>& a_lofPhi) = 0;

  /// Smoother at bottom level
  /**
  */
  virtual void bottomSmoother(LevelData<FArrayBox>&       a_phi,
                              const LevelData<FArrayBox>& a_rhs) = 0;

  /// Get flux at this level in the direction dir over the given Box
  /**
      The fluxes live on the cell faces with direction dir.
      Fluxes are computed for all interior edges of \em data.
      The flux fab is resized inside the routine.
  */
  virtual void getFlux(/// Face-centered Fluxes in direction dir (resized in function)
                       FArrayBox&       a_fineFlux,
                       /// Cell-centered box over which fluxes will be computed
                       const FArrayBox& a_data,
                       /// Data index of grid where flux is located
                       const DataIndex& a_datInd,
                       /// Direction of fluxes to be computed
                       int              a_dir) = 0;

  /// set the convergence metrics (if appropriate) for this problem
  /** If this LevelOp is being used in the context of a solve, we may
      wish to set a metric against which to define convergence (for example,
      the norm of the RHS).  In many cases, this will simply be a pass-through
      to a bottom solver.
  */
  virtual void setConvergenceMetric(Real a_metric, int a_comp) = 0;
};

#include "NamespaceFooter.H"
#endif
