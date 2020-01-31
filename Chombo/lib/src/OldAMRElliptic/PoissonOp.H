#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// DTGraves, Weds, July 21, 1999

#ifndef _POISSONOP_H_
#define _POISSONOP_H_

#include <iostream>
#include <cmath>
#include "SPACE.H"
#include <cstdlib>

#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BaseFab.H"
#include "FArrayBox.H"
#include "QuadCFInterp.H"
#include "Copier.H"

#include "GhostBC.H"
#include "LevelOp.H"
#include "CFIVS.H"
#include "BaseBottomSmoother.H"
#include "NamespaceHeader.H"

///  one-sided differencing version of poisson operator
/** this class contains stuff to solve  div (grad phi) = rhs
 */
class PoissonOp: public LevelOp
{
public:
  /// Default constructor
  PoissonOp();

  /// Destructor
  ~PoissonOp();

  /// Virtual constructor workaround
  /**
  */
  virtual LevelOp* new_levelop() const;

  ///
  /**
  */
  void levelGSMC(LevelData<FArrayBox>&       a_phi,
                 const LevelData<FArrayBox>& a_rhs);

  /// Has this object been defined?
  virtual bool isDefined() const;

  /// Full define function.
  /**
      Makes all coarse-fine information and sets internal variables
  */
  void define(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout* a_baseBAPtr,
              Real                     a_dxLevel,
              int                      a_refRatio,
              const ProblemDomain&     a_domain,
              bool                     a_homogeneousOnly=false,
              int                      a_ncomp = 1);

  /// Full define function. Define from finer levelop.
  /**
  */
  void define(const LevelOp* a_opfine,
              int            a_reftoFine);

  /// Full define function (Deprecated -- use ProblemDomain version)
  /**
      Full define function.  makes all coarse-fine information and sets
      internal variables
  */
  void define(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout* a_baseBAPtr,
              Real                     a_dxLevel,
              int                      a_refRatio,
              const Box&               a_domain,
              bool                     a_homogeneousOnly=false,
              int                      a_ncomp = 1);

  /// Smoother.
  /**
      Assumes that problem has already been put in residual correction form,
      so that C/F boundary conditions are homogeneous.
  */
  void smooth(LevelData<FArrayBox>&       a_phi,
              const LevelData<FArrayBox>& a_rhs);

  /// Apply preconditioning to the solution on this level.
  /**
      Given rhshat, returns phihat from M(phihat) = rhshat
  */
  virtual void levelPreconditioner(LevelData<FArrayBox>&       a_phihat,
                                   const LevelData<FArrayBox>& a_rhshat);

  /// Evaluate Operator, inhomogeneous boundary conditions.
  /**
      Uses inhomogeneous C/F boundary conditions.
      and inhomogeneous phys boundary conditions
  */
  void applyOpI(LevelData<FArrayBox>&       a_phi,
                const LevelData<FArrayBox>* a_phicPtr,
                LevelData<FArrayBox>&       a_lofPhi);

  /// Evaluate Operator, mixed boundary conditions.
  /**
      Uses inhomogeneous form of C/F boundary conditions.
      and homogeneous physical boundary conditions
  */
  void applyOpIcfHphys(LevelData<FArrayBox>&       a_phi,
                       const LevelData<FArrayBox>* a_phicPtr,
                       LevelData<FArrayBox>&       a_lofPhi);

  /// Evaluate Operator using homogeneous boundary conditions.
  /**
      Uses homogeneous form of both coarse-fine and physical  boundary
      conditions.
  */
  void applyOpH(LevelData<FArrayBox>& a_phi,
                LevelData<FArrayBox>& a_lofPhi);

  /// Coarse / Fine interpolation operator.
  /**
  */
  void CFInterp(LevelData<FArrayBox>&       a_phif,
                const LevelData<FArrayBox>& a_phic);

  /// Evaluate Operator
  /**
      Uses homogeneous form of C/F boundary conditions.
      with inhomogeneous physical boundary conditions
  */
  void applyOpHcfIphys(LevelData<FArrayBox>& a_phi,
                       LevelData<FArrayBox>& a_lOfPhi);

  /// Homogeneous coarse/fine interpolation operator
  /**
  */
  void homogeneousCFInterp(LevelData<FArrayBox>& a_phif);

  /// Set physical boundary condition types
  void setDomainGhostBC(const DomainGhostBC& a_dombcin);

  /// Set the bottom smoother to use in solve
  void setBottomSmoother(const BaseBottomSmoother& a_bottomSmoother);

  /// Does homogeneous coarse/fine interpolation
  /**
  */
  void
  homogeneousCFInterp(LevelData<FArrayBox>& a_phif,
                      const DataIndex&      a_datInd,
                      int                   a_idir,
                      Side::LoHiSide        a_hiorlo);

  /// This performs the smoothing at the bottom of a multigrid V-cycle.
  /**
  */
  void bottomSmoother(LevelData<FArrayBox>&       a_phi,
                      const LevelData<FArrayBox>& a_rhs);

  /// Get flux on this level in the direction dir on the given Box
  /**
      The fluxes live on the cell faces with direction a_dir.
      Fluxes are computed for all interior edges of a_data.
      The flux fab is resized inside the routine.
  */
  void getFlux(FArrayBox&       a_fineFlux,
               const FArrayBox& a_data,
               const DataIndex& a_datInd,
               int              a_dir);

  /// pass-through to bottom smoother
  virtual void setConvergenceMetric(Real a_metric, int a_comp);

  /// number of GSRB's to do in the LevelPreconditioner call
  void numPrecondSmooth(int a_numSmooth);

protected:
  /**
     Boundary conditons
  */
  DomainGhostBC m_domghostbc;

  /**
     Domain of fine grid
  */
  ProblemDomain m_domain;

  /**
     The grids at the current level
  */
  DisjointBoxLayout m_grids;

  /**
     Copier object that knows how to perform exchange operations on
     LevelData objects defined on m_grids
  */
  Copier m_exchangeCopier;

  /**
     DisjointBoxLayout for next coarser level
     this is == NULL if there is no coarser level
  */
  DisjointBoxLayout m_baseBA;

  /**
     Number of components
  */
  int m_ncomp;

  /**
     Refinement ratio between this and the next coarser level
  */
  int m_refRatio;

  /**
     Mesh spacing at this level
  */
  Real m_dxLevel;

  /**
     Mesh spacing at next coarser level
  */
  Real m_dxCrse;

  /*
    Object which performs Coarse-fine interpolation between this level
    and the next coarser level
  */
  QuadCFInterp m_quadCFI;

  /**
     Bottom smoother object
  */
  BaseBottomSmoother* m_bottomSmootherPtr;

  /**
     Has full define function been called?
  */
  bool m_isDefined;

  /**
     Has bc been defined?
  */
  bool m_isBCDefined;

  /**
     Is inhomogeneous CF interpolation possible?
  */
  bool m_ihcfiEnabled;

  /// number of GSRB to do in LevelPreconditioner
  int m_numPrecondSmooth;

  /**
     Coarse-fine stencils for homogeneous CFInterp
  */
  LayoutData<CFIVS> m_loCFIVS[SpaceDim];
  /**
     Coarse-fine stencils for homogeneous CFInterp
  */
  LayoutData<CFIVS> m_hiCFIVS[SpaceDim];

private:
  // Internally useful functions
  void levelGSRB(LevelData<FArrayBox>&       a_phi,
                 const LevelData<FArrayBox>& a_rhs);

//   void interpOnIVS(LevelData<FArrayBox>& a_phif,
//                    const FArrayBox&      a_phistar,
//                    const DataIndex&      a_datInd,
//                    const int             a_idir,
//                    const Side::LoHiSide  a_hiorlo,
//                    const IntVectSet&     a_interpIVS);

  void interpOnIVSHomo(LevelData<FArrayBox>& a_phif,
                       const DataIndex&      a_datInd,
                       const int             a_idir,
                       const Side::LoHiSide  a_hiorlo,
                       const IntVectSet&     a_interpIVS);

  void clearMemory();

  void setDefaultValues();
};

#include "NamespaceFooter.H"
#endif
