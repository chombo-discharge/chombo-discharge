#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TENSOROP_H_
#define _TENSOROP_H_

#include <iostream>
#include <cmath>
#include "SPACE.H"
#include <cstdlib>
#include <REAL.H>
#include <IntVect.H>
#include <Box.H>
#include <DisjointBoxLayout.H>
#include <LevelData.H>
#include <BaseFab.H>
#include "FArrayBox.H"
#include "FluxBox.H"
#include "GhostBC.H"
#include "LevelOp.H"
#include "TensorCFInterp.H"
#include "CFIVS.H"
#include "Copier.H"
#include "BaseBottomSmoother.H"
#include "NamespaceHeader.H"

///  LevelOp for Tensor solves
/** this class contains stuff to solve  div (grad U) + grad(div U) = rhs
 */
class TensorOp: public LevelOp
{
public:

  ///
  /**
     virtual constructor workaround
  */
  virtual LevelOp* new_levelop() const;

  ///
  virtual bool isDefined() const;

  ///
  /**
     Full define function.  makes all coarse-fine
     information and sets internal variables
  */
  void define(
              const DisjointBoxLayout& Ba,
              const DisjointBoxLayout* base_ba,
              Real  DxLevel,
              int refratio,
              const Box& domf,
              bool a_homogeneousOnly=false,
              int ncomp = 1);

  ///
  /**
     Full define function.  makes all coarse-fine
     information and sets internal variables
  */
  void define(
              const DisjointBoxLayout& Ba,
              const DisjointBoxLayout* base_ba,
              Real  DxLevel,
              int refratio,
              const ProblemDomain& domf,
              bool a_homogeneousOnly=false,
              int ncomp = 1);

  ///
  /**
     Full define function.  define
     from finer levelop.
  */
  void define(
              const LevelOp* opfine,
              int refratio);

  ///
  TensorOp();

  ///
  ~TensorOp();

  ///
  /**
     Smoother.
     Assumes that problem has already been put in
     residual correction form,
     so that C/F boundary conditions are homogeneous.
  */
  void smooth(
              LevelData<FArrayBox>& a_phi,
              const LevelData<FArrayBox>& a_rhs);



  ///
  /** apply preconditioning to the solution on this level.
      Given rhshat, returns phihat from M(phihat) = rhshat
  */
  virtual void levelPreconditioner(
                                   LevelData<FArrayBox>& a_phihat,
                                   const LevelData<FArrayBox>& a_rhshat
                                   );

  ///
  /**
     Evaluate Operator, inhomogeneous C/F boundary conditions.
     and inhomogeneous phys boundary conditions
  */
  void applyOpI(
                LevelData<FArrayBox>& a_phi,
                const LevelData<FArrayBox>* a_phiCoarse,
                LevelData<FArrayBox>& a_LOfPhi);

  ///
  /**
     Evaluate Operator, inhomogeneous C/F boundary conditions.
     with homogeneous physical boundary conditions
  */
  void applyOpIcfHphys(
                       LevelData<FArrayBox>& a_phi,
                       const LevelData<FArrayBox>* a_phiCoarse,
                       LevelData<FArrayBox>& a_LOfPhi);

  ///
  /**
     Evaluate Operator, homogeneous C/F boundary conditions.
     and homogeneous physical boundary conditions
  */
  void applyOpH(
                LevelData<FArrayBox>& a_phi,
                LevelData<FArrayBox>& a_LOfPhi);

  ///
  /**
     Coarse / Fine interpolation operator.
  */
  void CFInterp(LevelData<FArrayBox>& a_phi,
                const LevelData<FArrayBox>& a_phiCoarse
                );

  ///
  /**
      Evaluate Operator, homogeneous C/F boundary conditions.
      with inhomogeneous physical boundary conditions
  */
  void applyOpHcfIphys(
                       LevelData<FArrayBox>& a_phi,
                       LevelData<FArrayBox>& a_lOfPhi);


  ///
  /**
     homogeneous coarse/fine interpolation operator
  */
  void homogeneousCFInterp(LevelData<FArrayBox>& a_phi);

  ///
  void setDomainGhostBC(const DomainGhostBC& a_dombcin);

  /// set boundary conditions for tangential gradients (for divergence)
  void setTanGradBC(const DomainGhostBC& a_dombcIn);

  ///
  void setBottomSmoother(const BaseBottomSmoother& a_bottomSmoother);

  ///
  /**
     does homogeneous coarse/fine interpolation for phi
  */
  void
  homogeneousCFInterpPhi(LevelData<FArrayBox>& a_phif,
                         const DataIndex& a_datInd,
                         int a_idir,
                         Side::LoHiSide a_hiorlo);

  ///
  /** does homogeneous coarse/fine interpolation for tangential gradient
      (needs phi ghost cells to be filled in first, so should call
      homogeneousCFInterpPhi first)
  */
  void homogeneousCFInterpTanGrad(LevelData<FArrayBox>& a_tanGrad,
                                  const LevelData<FArrayBox>& a_phi,
                                  const DataIndex& a_datInd,
                                  int a_idir,
                                  Side::LoHiSide a_hiorlo);

  ///
  /**
     This does a GSRB Pre/Conditioned BiCGStab on a level
     for the bottom solver.
  */
  void bottomSmoother(
                      LevelData<FArrayBox>& a_phi,
                      const LevelData<FArrayBox>& a_rhs);

  ///
  /**
      get flux( == flux at THIS level)
      The fluxes live on the cell faces with direction dir.
      Fluxes are computed for all interior edges of data.
      The flux fab is resized inside the routine.
  */
  void getFlux(
               FArrayBox& flux,
               const FArrayBox& a_data,
               const DataIndex& a_datInd,
               int a_dir);

  /// pass-through to bottom smoother
  virtual void setConvergenceMetric(Real a_metric, int a_comp);


  /// access function
  const LevelData<FArrayBox>& getTanGrad() const;

protected:

  // boundary conditons
  DomainGhostBC m_domghostbc;

  // boundary conditions for tangential component of gradient
  DomainGhostBC m_tangradbc;

  //domain of fine grid
  ProblemDomain m_domain;

  //the grids at the current level
  DisjointBoxLayout m_grids;

  // Copier object that knows how to perform exchange operations on
  // LevelData objects defined on m_grids
  Copier m_exchangeCopier;

  /* DisjointBoxLayout for next coarser level
     this is == NULL if there is no coarser level */
  DisjointBoxLayout m_baseBA;

  // number of components
  int m_ncomp;

  // refinement ratio between this and the next coarser level
  int m_refRatio;

  //mesh spacing at this level
  Real m_dxLevel;

  //mesh spacing at next coarser level
  Real m_dxCrse;

  /*
    coarse-fine information between this level
    and the next coarser level
  */
  TensorCFInterp m_quadCFI;


  ///
  /* bottom smoother object
   */
  BaseBottomSmoother* m_bottom_smoother_ptr;

  //contains tangential gradient info
  LevelData<FArrayBox> m_tanGrad;

  //has full define function been called?
  bool m_isDefined;

  //has bc been defined?
  bool m_isBCDefined;

  // has bc for tangential gradients been defined?
  bool m_isGradBCDefined;

  //is inhomogeneous CF interpolation possible?
  bool m_ihcfiEnabled;

  LayoutData<CFIVS> m_loCFIVS[SpaceDim];
  LayoutData<CFIVS> m_hiCFIVS[SpaceDim];

  // will need these to do tangential gradient computations
  LayoutData<TensorFineStencilSet> m_hiTanStencilSets[SpaceDim];
  LayoutData<TensorFineStencilSet> m_loTanStencilSets[SpaceDim];

private:
  //internally useful functions
  void
  levelGSRB(
            LevelData<FArrayBox>& a_phi,
            const LevelData<FArrayBox>& a_rhs);
  void
  interpPhiOnIVS(
                 LevelData<FArrayBox>& a_phif,
                 const FArrayBox& a_phistar,
                 const DataIndex& dFine,
                 const int a_idir,
                 const Side::LoHiSide a_hiorlo,
                 const IntVectSet& a_interpIVS);

  void
  computeFaceDiv(LevelData<FluxBox>& a_div,
                 const LevelData<FArrayBox>& a_vel,
                 const LevelData<FArrayBox>& a_tanGrad);

  void
  computeFaceDiv(FArrayBox& a_div,
                 const FArrayBox& a_vel,
                 const FArrayBox& a_tanGrad,
                 const DataIndex& a_dataInd,
                 int faceDir);

  // computes tangential gradients of a_vel and stores in m_tanGrad
  void
  computeTanGradInterior(const LevelData<FArrayBox>& a_vel);

  // computes tangential gradients of a_vel and stores in m_tanGrad
  void
  computeTanGradInterior(const FArrayBox& a_vel,
                         const DataIndex& a_dataInd);


  // computes gradient of a_vel in a_gradDir and stores in m_tanGrad
  void
  computeTanGradInterior(const FArrayBox& a_vel,
                         const DataIndex& a_dataInd,
                         int a_gradDir);

  void clearMemory();

  void setDefaultValues();

};

#include "NamespaceFooter.H"
#endif
