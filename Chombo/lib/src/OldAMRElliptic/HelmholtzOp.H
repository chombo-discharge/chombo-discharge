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

#ifndef _HELMHOLTZOP_H_
#define _HELMHOLTZOP_H_

#include <iostream>
#include <cmath>
#include "SPACE.H"
#include <cstdlib>

#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "QuadCFInterp.H"
#include "Copier.H"

#include "GhostBC.H"
#include "LevelOp.H"
#include "CFIVS.H"
#include "BaseBottomSmoother.H"
#include "BaseHelmholtzOp.H"
#include "NamespaceHeader.H"

/// Helmholtz operator class
/**
    Operator to solve Helmholtz equation (ie (alpha*I+beta*laplacian)phi = rho)
 */
class HelmholtzOp: public BaseHelmholtzOp
{
public:
  /// Set helmholtz coefficient beta in (alpha*I +beta*laplacian)phi = rho.
  void setHelmCoeff(Real a_helmcoeff);

  /// Scale the helmholtz coefficient beta in (alpha*I + beta*laplacian)phi
  /**
     beta := a_scale*beta.  This function exists for backward compatibility
     with the old helmholtzOp.  in reality, this just calls scaleBeta
  */
  void scaleHelmCoeff(Real a_scale);

  /// Set helmholtz coefficient alpha in (alpha*I +beta*laplacian)phi = rho.
  void setAlpha(Real a_alpha);

  /// Function to scale the coefficient alpha
  /**
      Scale the helmholtz coefficient alpha in (alpha*I + beta*laplacian)phi.
  */
  void scaleAlpha(Real a_scale);

  /// Scale the helmholtz coefficient beta in (alpha*I + beta*laplacian)phi
  void scaleBeta(Real a_scale);

  /// Virtual constructor workaround
  /**
      Coefficient and bcs must be defined before called.
  */
  LevelOp* new_levelop() const;

  /// Full define function (Deprecated -- use ProblemDomain)
  /**
      Makes all coarse-fine information and sets internal variables
  */
  void define(const DisjointBoxLayout& a_grids,
              const DisjointBoxLayout* a_baseBAPtr,
              Real                     a_dxLevel,
              int                      a_refRatio,
              const Box&               a_domain,
              bool                     a_homogeneousOnly= false,
              int                      a_ncomp = 1);

  /// Full define function.
  /**
      Makes all coarse-fine information and sets internal variables
  */
  void define(/// Grids at this level
              const DisjointBoxLayout& a_grids,
              /// Coarser-level grids
              const DisjointBoxLayout* a_baseBAPtr,
              /// Cell spacing at this level
              Real                     a_dxLevel,
              /// Refinement ratio between this level and coarser
              int                      a_refRatio,
              /// Problem domain at this level
              const ProblemDomain&     a_domain,
              /// Only need homogeneous BCs?
              bool                     a_homogeneousOnly= false,
              /// number of components
              int                      a_ncomp = 1);

  /// Define from finer levelop.
  /**
  */
  void define(const LevelOp* a_opfine,
              int            a_refToFine);

  /// Default constructor
  HelmholtzOp();

  /// Destructor
  ~HelmholtzOp();

  /// Smoother.
  /**
     Assumes that problem has already been put in
     residual correction form,
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

  /// Evaluate Operator with inhomogeneous BCs
  /**
     Uses inhomogeneous form of both C/F boundary conditions
     and phys boundary conditions.
  */
  void applyOpI(
                LevelData<FArrayBox>& phi,
                const LevelData<FArrayBox>* phiCoarse,
                LevelData<FArrayBox>& LOfPhi
                );

  /// Evaluate operator
  /**
      Uses inhomogeneous C/F boundary conditions.
      homogeneous domain boundary conditions.
  */
  void applyOpIcfHphys(
                       LevelData<FArrayBox>& phi,
                       const LevelData<FArrayBox>* phiCoarse,
                       LevelData<FArrayBox>& LOfPhi
                       );

  /// Evaluate operator with homogeneous boundary conditions
  /**
     Uses homogeneous form of both C/F boundary conditions
     and phys boundary conditions.
  */
  void applyOpH(
                LevelData<FArrayBox>& phi,
                LevelData<FArrayBox>& LOfPhi
                );

  /// Evaluate Operator.
  /**
     Uses homogeneous C/F boundary conditions.
     inhomogeneous domain boundary conditions
  */
  void applyOpHcfIphys(
                       LevelData<FArrayBox>& phi,
                       LevelData<FArrayBox>& LOfPhi
                       );

  /// Coarse / Fine interpolation operator. fills ghost cells of phi.
  /**
  */
  void CFInterp(
                LevelData<FArrayBox>& phi,
                const LevelData<FArrayBox>& phiCoarse
                );

  /// homogeneous coarse/fine interpolation operator
  /**
   */
  void homogeneousCFInterp(LevelData<FArrayBox>& phi);

  /// Sets the physical boundary conditions for the operator.
  /**
   */
  void setDomainGhostBC(const DomainGhostBC& a_dombcin);

  /// Sets the bottom smoother which will be used.
  /**
     Default is BiCGStabSmoother
  */
  void setBottomSmoother(const BaseBottomSmoother& a_bottomSmoother);

  /// does homogeneous coarse/fine interpolation
  /**
  */
  void
  homogeneousCFInterp(LevelData<FArrayBox>& phif,
                      const DataIndex& a_datIndex,
                      int idir,
                      Side::LoHiSide hiorlo);

  /// This sends the problem to the bottom smoother
  /**
  */
  void bottomSmoother(LevelData<FArrayBox>& phi,
                      const LevelData<FArrayBox>& rhs);

  /// Has full define function been called?
  bool isDefined() const;

  /// Fills flux with the flux at this level
  /**
      The fluxes live on the cell faces with direction dir.
      Fluxes are computed for all interior edges of data.
      The flux fab is resized inside the routine.
  */
  void getFlux(
               FArrayBox& flux,
               const FArrayBox& data,
               const DataIndex& a_datInd,
               int dir);

  /// pass-through to bottom smoother
  virtual void setConvergenceMetric(Real a_metric, int a_comp);

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
     Refinement ratio between this and the next coarser level
  */
  int m_refRatio;

  /**
     Mesh spacing at this level
  */
  Real m_dxLevel;

  /**
     Mesh spacing at the next coarser level
  */
  Real m_dxCrse;

  /**
     Number of components
  */
  int m_ncomp;

  /**
    Coarse-fine information between this level
    and the next coarser level
  */
  QuadCFInterp m_quadCFI;

  ///
  /**
      Bottom smoother object
  */
  BaseBottomSmoother* m_bottomSmootherPtr;

  /**
     Has full define function been called?
  */
  bool m_isDefined;

  /**
     Has the boundary condition been set?
  */
  bool m_isBCDefined;

  /**
     Is inhomogeneous CF interpolation possible?
  */
  bool m_ihcfiEnabled;

  /**
     Is the helmholtz coefficient defined?
  */
  bool m_hcoeffDefined;

  /**
     Intvectset for coarse-fine interpolation
  */
  LayoutData<CFIVS> m_loCFIVS[SpaceDim];
  /**
     Intvectset for coarse-fine interpolation
  */
  LayoutData<CFIVS> m_hiCFIVS[SpaceDim];

  /**
     The beta in (alpha*I+beta*laplacian)phi = rho
  */
  Real m_helmCoeff;

  /**
      The alpha in (alpha*I +beta*laplacian)phi = rho
  */
  Real m_alphaCoeff;

private:
  /**
     Internally useful -- performs a GSRB smoothing
  */
  void levelGSRB(LevelData<FArrayBox>&       a_phi,
                 const LevelData<FArrayBox>& a_rhs);

//  void interpOnIVS(LevelData<FArrayBox>& a_phif,
//                   const FArrayBox&      a_phistar,
//                   const DataIndex&      a_datInd,
//                   const int             a_idir,
//                   const Side::LoHiSide  a_hiorlo,
//                   const IntVectSet&     a_interpIVS);

  /**
     Does the homoegeneous CF interpolation
  */
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
