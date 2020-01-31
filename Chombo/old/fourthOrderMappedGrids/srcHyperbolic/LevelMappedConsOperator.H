#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELMAPPEDCONSOPERATOR_H_
#define _LEVELMAPPEDCONSOPERATOR_H_

// Original design:                     Phil Colella
// Original implementation:
// Implementation date:

#include "LevelConsOperator.H"
#include "MOLPhysics.H"
#include "MultiBlockCoordSys.H"
#include "MultiBlockLevelExchangeAverage.H"
#include "BlockRegister.H"
#include "LevelSourceTerm.H"

#include "NamespaceHeader.H"

/// Level Cons Operator - a class for evaluating the Laplacian.

/**
 */
class LevelMappedConsOperator : public LevelConsOperator
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
   */
  LevelMappedConsOperator();

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  ~LevelMappedConsOperator();

  /// set coordinate system object
  void setCoordSys(MultiBlockCoordSys* a_coordSysPtr);

  /// set multiblock level exchange object
  void setLevelExchange(MultiBlockLevelExchangeAverage* a_mblexPtr);

  /// sets whether to use source term
  void useSourceTerm(bool a_useSourceTerm);

  /// set source term object
  void setSourceTerm(LevelSourceTerm* a_sourceTermPtr);

  void cellUJToCellU(LevelData<FArrayBox>& a_Uavg,
                     const LevelData<FArrayBox>& a_UJavg) const;

  void cellUToCellUJ(LevelData<FArrayBox>& a_JUavg,
                     const LevelData<FArrayBox>& a_Uavg) const;

  virtual void evalRHS(
                       LevelData<FArrayBox>&       a_LofU,
                       LevelData<FArrayBox>&       a_U,
                       LevelFluxRegister&          a_finerFluxRegister,
                       LevelFluxRegister&          a_coarserFluxRegister,
                       const LevelData<FArrayBox>& a_UcoarseOld,
                       const Real&                 a_timeCoarseOld,
                       const LevelData<FArrayBox>& a_UcoarseNew,
                       const Real&                 a_timeCoarseNew,
                       Real                        a_time,
                       Real                        a_weight);

  virtual void evalRHSpatches(
                              LevelData<FArrayBox>&       a_LofU,
                              const LevelData<FArrayBox>& a_U,
                              LevelFluxRegister&          a_finerFluxRegister,
                              LevelFluxRegister&          a_coarserFluxRegister,
                              Real                        a_weight);

  virtual void exchangeGhosts(LevelData<FArrayBox>&   a_U);

protected:

  // use source term?  Default false.  If true, in m_sourceTermPtr.
  bool m_useSourceTerm;

  // put this in PatchMappedConsOperator
  // as well as computeFaceAreas computeFaceGradU computeNJInvGradUOnBox
  // computeNTNJInvGradUOnBox computeNormalDiffusiveFluxVector
  // secondOrderMappedGridDivergence divideOutCellVolume secondOrderLaplacian

  virtual void computeCompFaceFluxes(LevelData<FluxBox>& uTimesV,
                                     const LevelData<FluxBox>& u,
                                     const LevelData<FluxBox>& v) const;

  //  virtual void computePhysFaceFluxes(LevelData<FluxBox>&,
  //                                     const LevelData<FluxBox>&);

  virtual void getPhysicalCellVolumes(LevelData<FArrayBox>&) const;

  virtual void computeFaceAverages(LevelData<FluxBox>& face_data,
                                   const LevelData<FArrayBox>& cell_data ) const;

  virtual void computeFaceCenters(LevelData<FluxBox>& face_data,
                                  const LevelData<FArrayBox>& cell_data ) const;

  virtual void setCommonFlux(LevelData<FluxBox>&   a_flux,
                             const BlockRegister&  a_blockRegister) const;

  //  virtual void secondOrderLaplacian( LevelData<FArrayBox>&       a_lap_u,
  //                                     const LevelData<FArrayBox>& a_u,
  //                                     const bool a_div_cell_volume = true );

  //  virtual void secondOrderMappedGridDivergence( LevelData<FArrayBox>&      a_divF,
  //                                                const LevelData<FluxBox>&  a_F);

  //  virtual void computeNormalDiffusiveFluxVector( LevelData<FluxBox>& a_flux,
  //                                                 const LevelData<FArrayBox>& a_u );

  //  virtual void computeNJInvGradUOnBox( FluxBox&                   a_NJInvGradU,
  //                                       const FluxBox&             a_NJInv,
  //                                       const FluxBox&             a_gradU,
  //                                       const int                  a_comp );

  //  virtual void computeNTNJInvGradUOnBox( FluxBox&                   a_fluxNormal,
  //                                         const FluxBox&             a_N,
  //                                         const FluxBox&             a_NJInvGradU,
  //                                         const int                  a_comp );

  //  virtual void computeFaceGradU( LevelData<FluxBox>&         a_grad_u,
  //                                 const LevelData<FArrayBox>& a_u );

  //  NewFourthOrderCoordSys* m_coordSysPtr;
  MultiBlockCoordSys* m_coordSysPtr;

  MultiBlockLevelExchangeAverage* m_mblexPtr;

  // Source term doohickey.
  LevelSourceTerm* m_sourceTermPtr;

  // BlockRegister m_blockRegister;
};

#include "NamespaceFooter.H"

#endif
