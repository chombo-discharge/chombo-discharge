/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_EBHelmholtzOp.H
  @brief  Declaration of Helmholtz multigrid operators. 
  @author Robert Marskar
*/

#ifndef CD_EBHELMHOLTZOP_H
#define CD_EBHELMHOLTZOP_H

// Std includes
#include <map>

// Chombo includes
#include <BaseEBBC.H>
#include <LevelTGA.H>
#include <BaseIVFAB.H>
#include <VCAggStencil.H>

// Our includes
#include <CD_Timer.H>
#include <CD_Location.H>
#include <CD_EBMultigridInterpolator.H>
#include <CD_EBCoarAve.H>
#include <CD_EBReflux.H>
#include <CD_EBMGRestrict.H>
#include <CD_EBMGProlong.H>
#include <CD_EBHelmholtzEBBC.H>
#include <CD_EBHelmholtzDomainBC.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Helmholtz operator for equations like alpha*a(x)*phi(x) + beta*div(b(x)*grad(phi(x))) = rho. 
  @details This can be used with TGA time stepping. 
*/
class EBHelmholtzOp : public LevelTGAHelmOp<LevelData<EBCellFAB>, EBFluxFAB>
{
public:
  /**
    @brief Relaxation method for the operators
  */
  enum class Smoother
  {
    NoRelax,
    PointJacobi,
    GauSaiRedBlack,
    GauSaiMultiColor,
    Chebyshev,
    RestrictedAdditiveSchwarz,
  };

  /**
    @brief Disallowed default constructor
  */
  EBHelmholtzOp() = delete;

  /**
    @brief Disallowed copy constructor
    @param[in] a_other Other operator
  */
  EBHelmholtzOp(const EBHelmholtzOp& a_other) = delete;

  /**
    @brief Disallowed move constructor
    @param[in] a_other Other operator
  */
  EBHelmholtzOp(const EBHelmholtzOp&& a_other) = delete;

  /**
    @brief Full constructor
    @param[in] a_dataLocation  Data location, either cell center or cell centroid
    @param[in] a_eblgFine      Fine grids
    @param[in] a_eblg          Grids on this level
    @param[in] a_eblgCoFi      Coarsening of fine level grids
    @param[in] a_eblgCoar      Coarse grids
    @param[in] a_eblgCoarMG    Multigrid-grids
    @param[in] a_validCells    Valid grid cells
    @param[in] a_interpolator  Interpolator
    @param[in] a_fluxReg       Flux register
    @param[in] a_coarAve       Coarsener
    @param[in] a_domainBC      Domain BC
    @param[in] a_ebBC          Boundary conditions on EBs
    @param[in] a_probLo        Lower-left corner of computational domain
    @param[in] a_dx            Grid resolution
    @param[in] a_refToFine     Refinement ratio to fine level
    @param[in] a_refToCoar     Refinement ratio to coarse level
    @param[in] a_hasFine       Has fine level or not
    @param[in] a_hasCoar       Has coarse level or not
    @param[in] a_hasMGObjects  Has multigrid-objects (special objects between AMR levels, or below the AMR levels)
    @param[in] a_alpha         Operator alpha
    @param[in] a_beta          Operator beta
    @param[in] a_Acoef         Operator A-coefficient
    @param[in] a_Bcoef         Operator B-coefficient
    @param[in] a_ghostPhi      Ghost cells in solution vector
    @param[in] a_ghostRhs      Ghost cells in right-hand side
    @param[in] a_smoother      Which smoother to use
    @param[in] a_relaxFactor   Relaxation factor
    @param[in] a_chebyOrder    Chebyshev polynomial degree (only used when a_smoother is Smoother::Chebyshev)
    @param[in] a_chebyEigRatio Chebyshev eigenvalue ratio lambda_max/lambda_min (only used when a_smoother is Smoother::Chebyshev)
    @param[in] a_rasInnerSweeps Number of inner frozen-ghost sweeps per restricted-additive-Schwarz iteration (only used when a_smoother is Smoother::RestrictedAdditiveSchwarz)
    @param[in] a_refluxFree    If true, use the reflux-free AMR operator (skips coarse-fine flux refluxing)
    @param[in] a_BcoefIrreg    B-coefficient on EB faces
  */
  EBHelmholtzOp(const Location::Cell                             a_dataLocation,
                const EBLevelGrid&                               a_eblgFine,
                const EBLevelGrid&                               a_eblg,
                const EBLevelGrid&                               a_eblgCoFi,
                const EBLevelGrid&                               a_eblgCoar,
                const EBLevelGrid&                               a_eblgCoarMG,
                const RefCountedPtr<LevelData<BaseFab<bool>>>&   a_validCells,
                const RefCountedPtr<EBMultigridInterpolator>&    a_interpolator,
                const RefCountedPtr<EBReflux>&                   a_fluxReg,
                const RefCountedPtr<EBCoarAve>&                  a_coarAve,
                const RefCountedPtr<EBHelmholtzDomainBC>&        a_domainBC,
                const RefCountedPtr<EBHelmholtzEBBC>&            a_ebBC,
                const RealVect&                                  a_probLo,
                const Real&                                      a_dx,
                const int&                                       a_refToFine,
                const int&                                       a_refToCoar,
                const bool&                                      a_hasFine,
                const bool&                                      a_hasCoar,
                const bool&                                      a_hasMGObjects,
                const Real&                                      a_alpha,
                const Real&                                      a_beta,
                const RefCountedPtr<LevelData<EBCellFAB>>&       a_Acoef,
                const RefCountedPtr<LevelData<EBFluxFAB>>&       a_Bcoef,
                const RefCountedPtr<LevelData<BaseIVFAB<Real>>>& a_BcoefIrreg,
                const IntVect&                                   a_ghostPhi,
                const IntVect&                                   a_ghostRhs,
                const Smoother&                                  a_smoother,
                const Real&                                      a_relaxFactor,
                const int&                                       a_chebyOrder,
                const Real&                                      a_chebyEigRatio,
                const int&                                       a_rasInnerSweeps,
                const bool                                       a_refluxFree = false);

  /**
    @brief Dtor
  */
  virtual ~EBHelmholtzOp();

  /**
    @brief No copy assignment allowed
    @param[in] a_oper Oper
  */
  EBHelmholtzOp&
  operator=(const EBHelmholtzOp& a_oper) = delete;

  /**
    @brief No move assignment allowed
    @param[in] a_oper Oper
  */
  EBHelmholtzOp&
  operator=(const EBHelmholtzOp&& a_oper) = delete;

  /**
    @brief Turn off BCs
  */
  void
  turnOffCFInterp();

  /**
    @brief Turn on BCs
  */
  void
  turnOnCFInterp();

  /**
    @brief Turn off exchange operation
    @note Done when MFHelmholtzOp does the exchange and we don't want EBHelmholtzOp to redo it.
  */
  void
  turnOffExchange();

  /**
    @brief Turn on exchange operation
    @note Done when MFHelmholtzOp does the exchange and we don't want EBHelmholtzOp to redo it.
  */
  void
  turnOnExchange();

  /**
    @brief Turn off coarsening operation
    @note Done when MFHelmholtzOp does the coarsening and we don't want EBHelmholtzOp to redo it.
  */
  void
  turnOffCoarsening();

  /**
    @brief Turn on coarsening operation
    @note Done when MFHelmholtzOp does the coarsening and we don't want EBHelmholtzOp to redo it.
  */
  void
  turnOnCoarsening();

  /**
    @brief Update with new A and B coefficients.
    @note This will call defineStencils but will not recompute the stencils from the BC objects. 
    @param[in] a_Acoef         Operator A-coefficient
    @param[in] a_Bcoef         Operator B-coefficient
    @param[in] a_BcoefIrreg    Operator B-coefficient (on EB faces)
  */
  void
  setAcoAndBco(const RefCountedPtr<LevelData<EBCellFAB>>&       a_Acoef,
               const RefCountedPtr<LevelData<EBFluxFAB>>&       a_Bcoef,
               const RefCountedPtr<LevelData<BaseIVFAB<Real>>>& a_BcoefIrreg);

  /**
    @brief Get the Helmholtz A-coefficient on cell centers
    @return m_Acoef
  */
  const RefCountedPtr<LevelData<EBCellFAB>>&
  getAcoef();

  /**
    @brief Get the Helmholtz B-coefficient on faces
    @return m_Bcoef
  */
  const RefCountedPtr<LevelData<EBFluxFAB>>&
  getBcoef();

  /**
    @brief Get the Helmholtz B-coefficient on the EB
    @return m_Bcoef
  */
  const RefCountedPtr<LevelData<BaseIVFAB<Real>>>&
  getBcoefIrreg();

  /**
    @brief Coarsen data from fine to coar level
    @param[in] a_phi     Data on this level
    @param[in] a_phiFine Data on finer level
    @note This routine is used in AMROperator to ensure that EB fluxes which are described by stencils get consistent data under the fine level. 
  */
  void
  coarsenCell(LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_phiFine);

  /**
    @brief Coarsen fluxes on the fine level onto this level. 
    @details User must update fluxes before calling this routine. 
    @param[in] a_flux     Flux on this level. 
    @param[in] a_fineFlux Flux on finer level
    @note This routine is used in the reflux-free AMROperator
  */
  void
  coarsenFlux(LevelData<EBFluxFAB>& a_flux, const LevelData<EBFluxFAB>& a_fineFlux);

  /**
    @brief Point Jacobi kernel. 
    @param[in,out] a_Lcorr      Storage for computing L(a_corr)
    @param[in,out] a_corr       Correction
    @param[in]    a_resid      Residual
    @param[in]    a_Acoef      A-coefficient
    @param[in]    a_Bcoef      B-coefficient
    @param[in]    a_BcoefIrreg B-coefficient on EB faces
    @param[in]    a_cellBox    Grid box
    @param[in]    a_dit        Data index
  */
  void
  pointJacobiKernel(EBCellFAB&             a_Lcorr,
                    EBCellFAB&             a_corr,
                    const EBCellFAB&       a_resid,
                    const EBCellFAB&       a_Acoef,
                    const EBFluxFAB&       a_Bcoef,
                    const BaseIVFAB<Real>& a_BcoefIrreg,
                    const Box&             a_cellBox,
                    const DataIndex&       a_dit) const noexcept;

  /**
    @brief Chebyshev smoother kernel.
    @details Same as pointJacobiKernel but scales the update by a_omega (the Chebyshev step size).
    @param[in,out] a_Lcorr      Storage for computing L(a_corr)
    @param[in,out] a_corr       Correction
    @param[in]    a_resid      Residual
    @param[in]    a_Acoef      A-coefficient
    @param[in]    a_Bcoef      B-coefficient
    @param[in]    a_BcoefIrreg B-coefficient on EB faces
    @param[in]    a_cellBox    Grid box
    @param[in]    a_dit        Data index
    @param[in]    a_omega      Chebyshev step-size scaling factor
  */
  void
  chebyshevKernel(EBCellFAB&             a_Lcorr,
                  EBCellFAB&             a_corr,
                  const EBCellFAB&       a_resid,
                  const EBCellFAB&       a_Acoef,
                  const EBFluxFAB&       a_Bcoef,
                  const BaseIVFAB<Real>& a_BcoefIrreg,
                  const Box&             a_cellBox,
                  const DataIndex&       a_dit,
                  const Real             a_omega) const noexcept;

  /**
    @brief Red-black Gauss-Seidel kernel
    @param[in,out] a_Lcorr      Storage for computing L(a_corr)
    @param[in,out] a_corr       Correction
    @param[in]    a_resid      Residual
    @param[in]    a_Acoef      A-coefficient
    @param[in]    a_Bcoef      B-coefficient
    @param[in]    a_BcoefIrreg B-coefficient on EB faces
    @param[in]    a_redBlack   Red or black
    @param[in]    a_cellBox    Grid box
    @param[in]    a_dit        Data index
  */
  void
  gauSaiRedBlackKernel(EBCellFAB&             a_Lcorr,
                       EBCellFAB&             a_corr,
                       const EBCellFAB&       a_resid,
                       const EBCellFAB&       a_Acoef,
                       const EBFluxFAB&       a_Bcoef,
                       const BaseIVFAB<Real>& a_BcoefIrreg,
                       const Box&             a_cellBox,
                       const DataIndex&       a_dit,
                       const int&             a_redBlack) const noexcept;

  /**
    @brief Multi-color Gauss-Seidel kernel
    @param[in,out] a_Lcorr      Storage for computing L(a_corr)
    @param[in,out] a_corr       Correction
    @param[in]    a_resid      Residual
    @param[in]    a_Acoef      A-coefficient
    @param[in]    a_Bcoef      B-coefficient
    @param[in]    a_BcoefIrreg B-coefficient on EB faces
    @param[in]    a_cellBox    Grid box
    @param[in]    a_dit        Data index
    @param[in]    a_color      Color
  */
  void
  gauSaiMultiColorKernel(EBCellFAB&             a_Lcorr,
                         EBCellFAB&             a_corr,
                         const EBCellFAB&       a_resid,
                         const EBCellFAB&       a_Acoef,
                         const EBFluxFAB&       a_Bcoef,
                         const BaseIVFAB<Real>& a_BcoefIrreg,
                         const Box&             a_cellBox,
                         const DataIndex&       a_dit,
                         const IntVect&         a_color) const noexcept;

  /**
    @brief Compute residual on this level. 
    @param[out] a_residual          Residual rhs - L(phi)
    @param[in]  a_phi               phi
    @param[in]  a_rhs               Right-hand side of system
    @param[in] a_homogeneousPhysBc Homogeneous phys bc
  */
  void
  residual(LevelData<EBCellFAB>&       a_residual,
           const LevelData<EBCellFAB>& a_phi,
           const LevelData<EBCellFAB>& a_rhs,
           const bool                  a_homogeneousPhysBc) override;

  /**
    @brief Precondition system before bottom solve
    @param[in] a_corr     Correction
    @param[in] a_residual Residual
    @details This just runs a few relaxations. 
  */
  void
  preCond(LevelData<EBCellFAB>& a_corr, const LevelData<EBCellFAB>& a_residual) override final;

  /**
    @brief Apply coarse-fine boundary conditions
    @param[in,out] a_phiFine         Fine data
    @param[in,out] a_phiCoar         Coarse data
    @param[in]    a_homogeneousCFBC Homogeneous CF or not (i.e. coar data is zero);
  */
  void
  interpolateCF(LevelData<EBCellFAB>& a_phiFine, const LevelData<EBCellFAB>* a_phiCoar, const bool a_homogeneousCFBC);

  /**
    @brief Do homogeneous coarse-fine interpolation
    @param[in] a_phi Data
  */
  void
  homogeneousCFInterp(LevelData<EBCellFAB>& a_phi);

  /**
    @brief Inhomogeneous coarse-fine interpolation
    @param[in,out] a_phiCoar Coarse data. 
    @param[in] a_phi Phi
  */
  void
  inhomogeneousCFInterp(LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_phiCoar);

  /**
    @brief Multi-colored Gauss-Seidel kernel. Public because MFHelmholtzOp may want to use use. 
    @param[in] a_phi Correction
    @param[in] a_rhs Residual
    @param[in] a_color "Color":
    @param[in] a_Lphi Description
  */
  void
  gauSaiMultiColor(LevelData<EBCellFAB>&       a_phi,
                   const LevelData<EBCellFAB>& a_Lphi,
                   const LevelData<EBCellFAB>& a_rhs,
                   const IntVect&              a_color) const;

  /**
    @brief Apply operator on this level. This is a more general version which can turn on/off homogeneous and CF bcs.
    @param[out] a_Lphi              L(phi)
    @param[out] a_phi               Phi on this level
    @param[out] a_phiCoar           Coarse-level phi. If you have a coar this 
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
    @param[in]  a_homogeneousCFBC   Use homogeneous coarse-fine bcs or not
    @note a_phi is not really const and we do a nasty cast. Leaving it as const because we only modify the ghost cells!
  */
  void
  applyOp(LevelData<EBCellFAB>&             a_Lphi,
          const LevelData<EBCellFAB>&       a_phi,
          const LevelData<EBCellFAB>* const a_phiCoar,
          const bool                        a_homogeneousPhysBC,
          const bool                        a_homogeneousCFBC);

  /**
    @brief Apply operator 
    @param[out] a_Lphi              L(phi)
    @param[in]  a_phi               Phi
    @param[in]  a_homogeneousPhysBc Homogeneous physical BCs or not
    @details This computes a_Lphi = L(a_phi) using homogeneous physical BCs or not
  */
  void
  applyOp(LevelData<EBCellFAB>& a_Lphi, const LevelData<EBCellFAB>& a_phi, bool a_homogeneousPhysBc) override final;

  /**
    @brief Apply operator in a grid box. 
    @param[out] a_Lphi              L(phi)
    @param[in]  a_phi               Phi
    @param[in]  a_cellBox           Grid box
    @param[in]  a_dit               Data indxe
    @param[in]  a_homogeneousPhysBC Homogeneous physical BCs or not
    @param[in] a_Acoef Description
    @param[in] a_Bcoef Description
    @param[in] a_BcoefIrreg Description
  */
  void
  applyOp(EBCellFAB&             a_Lphi,
          EBCellFAB&             a_phi,
          const EBCellFAB&       a_Acoef,
          const EBFluxFAB&       a_Bcoef,
          const BaseIVFAB<Real>& a_BcoefIrreg,
          const Box&             a_cellBox,
          const DataIndex&       a_dit,
          const bool             a_homogeneousPhysBC) const noexcept;

  /**
    @brief Apply operator in regular cells.
    @param[out] a_Lphi              L(phi)
    @param[in]  a_phi               Phi
    @param[in]  a_cellBox           Grid box
    @param[in]  a_dit               Data indxe
    @param[in]  a_homogeneousPhysBC Homogeneous physical BCs or not
    @param[in] a_Acoef Description
    @param[in] a_Bcoef Description
    @param[in] a_BcoefIrreg Description
  */
  void
  applyOpRegular(EBCellFAB&             a_Lphi,
                 EBCellFAB&             a_phi,
                 const EBCellFAB&       a_Acoef,
                 const EBFluxFAB&       a_Bcoef,
                 const BaseIVFAB<Real>& a_BcoefIrreg,
                 const Box&             a_cellBox,
                 const DataIndex&       a_dit,
                 const bool             a_homogeneousPhysBC) const noexcept;

  /**
    @brief Apply domain flux. 
    @param[in,out] a_phi               Cell data
    @param[in]    a_Bcoef             Helmholtz B-coefficient
    @param[in]    a_cellBox           Computation box
    @param[in]    a_dit               Data index
    @note This sets the ghost cells of phi to phi(i-1) = phi(i) + F*dx/bco such that when we take the centered difference we get
    -bco*(phi(i) - phi(i-1))/dx = F. 
    @param[in] a_homogeneousPhysBc Homogeneous phys bc
  */
  void
  applyDomainFlux(EBCellFAB&       a_phi,
                  const EBFluxFAB& a_Bcoef,
                  const Box&       a_cellBox,
                  const DataIndex& a_dit,
                  const bool       a_homogeneousPhysBc) const noexcept;

  /**
    @brief Fill domain flux. This fills the flux on the domain face using centered differencing ala applyDomainFlux
    @details a_flux is replaced by the user-specified flux on the domain faces.
    @param[in,out] a_flux                  Flux data holder.
    @param[in,out] a_phi                   Cell-centered data
    @param[in]    a_cellBox               Computation box
    @param[in]    a_dit                   Data index
    @param[in]    a_homogeneousPhysBC     If true, use homogeneous physical BCs
  */
  void
  fillDomainFlux(EBFluxFAB&       a_flux,
                 const EBCellFAB& a_phi,
                 const Box&       a_cellBox,
                 const DataIndex& a_dit,
                 const bool       a_homogeneousPhysBC);

  /**
    @brief Apply operator in irregular cells
    @param[out] a_Lphi              L(phi)
    @param[in]  a_phi               Cell-centered data
    @param[in]  a_cellBox           Cell box
    @param[in]  a_dit               Grid index
    @param[in]  a_homogeneousPhysBC Homogeneous physical BCs or not
    @param[in] a_Acoef Description
    @param[in] a_Bcoef Description
    @param[in] a_BcoefIrreg Description
    @param[in] a_alphaDiagWeight Description
  */
  void
  applyOpIrregular(EBCellFAB&             a_Lphi,
                   const EBCellFAB&       a_phi,
                   const EBCellFAB&       a_Acoef,
                   const EBFluxFAB&       a_Bcoef,
                   const BaseIVFAB<Real>& a_BcoefIrreg,
                   const BaseIVFAB<Real>& a_alphaDiagWeight,
                   const Box&             a_cellBox,
                   const DataIndex&       a_dit,
                   const bool             a_homogeneousPhysBC) const noexcept;

  /**
    @brief Create data which clones the layout of the other
    @param[out] a_lhs Data clone (returned data is not initialized)
    @param[out] a_rhs Data layout to be cloned
  */
  void
  create(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) override final;

  /**
    @brief Assign data.
    @details This does a copy from rhs to lhs
    @param[in,out] a_lhs Description
    @param[in,out] a_rhs Description
  */
  void
  assign(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) override final;

  /**
    @brief Assign lhs
    @details This is the version that is called by AMRMultiGrid::VCycle. Note that the other version
    might be called by other operators (e.g., EBBackwardEuler)
    @param[out] a_lhs    Outgoing data
    @param[in]  a_rhs    Incoming data
    @param[in]  a_copier Copier
  */
  void
  assignCopier(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Copier& a_copier) override final;

  /**
    @brief Local assignment function
    @param[out] a_lhs Equal to a_rhs on output
    @param[in] a_rhs Data
  */
  void
  assignLocal(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) override final;

  /**
    @brief Build copier
    @param[out] a_copier Copier for copying between a_lhs and a_rhs
    @param[in]  a_lhs    Copying from
    @param[in]  a_rhs    Copying to
  */
  void
  buildCopier(Copier& a_copier, const LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) override;

  /**
    @brief Compute the dot product??
    @return Return value
    @param[in,out] a_lhs Description
    @param[in,out] a_rhs Description
  */
  Real
  dotProduct(const LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) override final;

  /**
    @brief Increment operator
    @param[in,out] a_lhs   Data to be incremented
    @param[in]    a_rhs   Incrementation data
    @param[in]    a_scale Scaling factor
  */
  void
  incr(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Real a_scale) override final;

  /**
    @brief Set a_lhs = a*x + b*y
    @param[out] a_lhs Result data
    @param[in]  a_x   x-data
    @param[in]  a_y   y-data
    @param[in]  a_a   Scaling factor
    @param[in]  a_b   Scaling factor
  */
  void
  axby(LevelData<EBCellFAB>&       a_lhs,
       const LevelData<EBCellFAB>& a_x,
       const LevelData<EBCellFAB>& a_y,
       const Real                  a_a,
       const Real                  a_b) override final;

  /**
    @brief Scale data. Returns a_lhs = a_lhs*a_scale
    @param[in,out] a_lhs   Data to be scaled
    @param[in]     a_scale Scaling factor
  */
  void
  scale(LevelData<EBCellFAB>& a_lhs, const Real& a_scale) override final;

  /**
    @brief Local scaling function. Multiplies the left-hand side by the right-hand side.
    @param[out] a_lhs Equal to a_lhs*a_rhs on output
    @param[in] a_rhs Scaling data
  */
  static void
  scaleLocal(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) noexcept;

  /**
    @brief Compute norm of data
    @param[in] a_rhs   Data to compute norm of
    @param[in] a_order Norm order (currently ignored)
    @details a_order is currently ignored as we call a static routine in EBAMRPoissonOp for now.
    @return Norm value
  */
  Real
  norm(const LevelData<EBCellFAB>& a_rhs, const int a_order) override final;

  /**
    @brief Set data to zero
    @param[in,out] a_lhs Data to be set to zero. 
  */
  void
  setToZero(LevelData<EBCellFAB>& a_lhs) override final;

  /**
    @brief Create coarsened data
    @param[out] a_coarse  Coarse data
    @param[in]  a_fine    Fine data
    @param[in]  a_ghosted Include ghost cells or nto
  */
  void
  createCoarser(LevelData<EBCellFAB>& a_coarse, const LevelData<EBCellFAB>& a_fine, bool a_ghosted) override final;

  /**
    @brief Relaxation method. This does smoothing for the system L(correction) = residual
    @param[in,out] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of iterations
  */
  void
  relax(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, int a_iterations) override final;

  /**
    @brief Restrict residual onto coarse level
    @param[in,out] a_resCoar      Coarse residual
    @param[in,out] a_phi Phi on this level
    @param[in]    a_rhs Rhs on this level
  */
  void
  restrictResidual(LevelData<EBCellFAB>&       a_resCoar,
                   LevelData<EBCellFAB>&       a_phi,
                   const LevelData<EBCellFAB>& a_rhs) override final;

  /**
    @brief Prolongation method. 
    @param[out] a_phi  Correction on this level
    @param[out] a_correctCoarse Correction on coarse level
  */
  void
  prolongIncrement(LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_correctCoarse) override final;

  /**
    @brief Return coarsening factor to coarser level (1 if there is no coarser level);
    @return Return value
  */
  int
  refToCoarser() override final;

  /**
    @brief Apply the AMR operator, i.e. compute L(phi) in an AMR context. 
    @param[in]  a_phiFine           Phi on fine level
    @param[in]  a_phi               Phi on this level
    @param[in]  a_phiCoar           Phi on coar level
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
    @param[in]  a_finerOp           Finer operator
    @details This involves ghost cell interpolation if there's a coarse level, and refluxing if there's a fine level.
    @param[in] a_Lphi Lphi
  */
  void
  AMROperator(LevelData<EBCellFAB>&             a_Lphi,
              const LevelData<EBCellFAB>&       a_phiFine,
              const LevelData<EBCellFAB>&       a_phi,
              const LevelData<EBCellFAB>&       a_phiCoar,
              const bool                        a_homogeneousPhysBC,
              AMRLevelOp<LevelData<EBCellFAB>>* a_finerOp) override final;

  /**
    @brief Apply the AMR operator, i.e. compute L(phi) in an AMR context. 
    @details This is the reflux-free version of AMROperator.
    @param[in]  a_phiFine           Phi on fine level
    @param[in]  a_phi               Phi on this level
    @param[in]  a_phiCoar           Phi on coar level
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
    @param[in]  a_finerOp           Finer operator
    @details This involves ghost cell interpolation if there's a coarse level, and refluxing if there's a fine level.
    @param[in] a_Lphi Lphi
  */
  void
  refluxFreeAMROperator(LevelData<EBCellFAB>&             a_Lphi,
                        const LevelData<EBCellFAB>&       a_phiFine,
                        const LevelData<EBCellFAB>&       a_phi,
                        const LevelData<EBCellFAB>&       a_phiCoar,
                        const bool                        a_homogeneousPhysBC,
                        AMRLevelOp<LevelData<EBCellFAB>>* a_finerOp);

  /**
    @brief Compute residual on this level. AMR version. 
    @param[out] a_residual          Residual on this level
    @param[in]  a_phiFine           Phi on fine level
    @param[in]  a_phi               Phi on this level
    @param[in]  a_phiCoar           Phi on coar level
    @param[in]  a_rhs               Right-hand side on this level
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
    @param[in]  a_finerOp           Finer operator
  */
  void
  AMRResidual(LevelData<EBCellFAB>&             a_residual,
              const LevelData<EBCellFAB>&       a_phiFine,
              const LevelData<EBCellFAB>&       a_phi,
              const LevelData<EBCellFAB>&       a_phiCoar,
              const LevelData<EBCellFAB>&       a_rhs,
              bool                              a_homogeneousPhysBC,
              AMRLevelOp<LevelData<EBCellFAB>>* a_finerOp) override final;

  /**
    @brief Compute AMR residual on finest AMR level. 
    @param[out] a_residual          Residual on this level
    @param[in]  a_phi               Phi on this level
    @param[in]  a_phiCoar           Phi on coar level
    @param[in]  a_rhs               Right-hand side on this level
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
  */
  void
  AMRResidualNF(LevelData<EBCellFAB>&       a_residual,
                const LevelData<EBCellFAB>& a_phi,
                const LevelData<EBCellFAB>& a_phiCoar,
                const LevelData<EBCellFAB>& a_rhs,
                bool                        a_homogeneousPhysBC) override final;

  /**
    @brief Compute AMR residual on coarsest
    @param[out] a_residual Residual on this level
    @param[in]  a_phiFine           Phi on fine level
    @param[in]  a_phi               Phi on this level
    @param[in]  a_rhs               Right-hand side on this level
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
    @param[in]  a_finerOp           Finer operator
  */
  void
  AMRResidualNC(LevelData<EBCellFAB>&             a_residual,
                const LevelData<EBCellFAB>&       a_phiFine,
                const LevelData<EBCellFAB>&       a_phi,
                const LevelData<EBCellFAB>&       a_rhs,
                bool                              a_homogeneousPhysBC,
                AMRLevelOp<LevelData<EBCellFAB>>* a_finerOp) override final;

  /**
    @brief Apply the AMR operator, i.e. compute L(phi) in an AMR context, assuming no finer levels. 
    @param[out] a_Lphi              L(phi)
    @param[in]  a_phi               Phi on this level
    @param[in]  a_phiCoar           Phi on coar level
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
    @details This involves ghost cell interpolation if there's a coarse level, and refluxing if there's a fine level.
  */
  void
  AMROperatorNF(LevelData<EBCellFAB>&       a_Lphi,
                const LevelData<EBCellFAB>& a_phi,
                const LevelData<EBCellFAB>& a_phiCoar,
                bool                        a_homogeneousPhysBC) override final;

  /**
    @brief Apply the AMR operator, i.e. compute L(phi) in an AMR context, assuming no coarser AMR levels.
    @param[out] a_Lphi              L(phi)
    @param[in]  a_phiFine           Phi on fine level
    @param[in]  a_phi               Phi on this level
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
    @param[in]  a_finerOp           Finer operator
    @details This involves ghost cell interpolation if there's a coarse level, and refluxing if there's a fine level.
  */
  void
  AMROperatorNC(LevelData<EBCellFAB>&             a_Lphi,
                const LevelData<EBCellFAB>&       a_phiFine,
                const LevelData<EBCellFAB>&       a_phi,
                bool                              a_homogeneousPhysBC,
                AMRLevelOp<LevelData<EBCellFAB>>* a_finerOp) override final;

  /**
    @brief Restrict residual
    @param[out] a_residualCoarse Coarse residual
    @param[out] a_residual         Residual
    @param[out] a_correction       Correction on this level 
    @param[out] a_coarseCorrection Coarse level correction
    @param[in]  a_skip_res         I have no idea what this one is supposed to do. 
  */
  void
  AMRRestrict(LevelData<EBCellFAB>&       a_residualCoarse,
              const LevelData<EBCellFAB>& a_residual,
              const LevelData<EBCellFAB>& a_correction,
              const LevelData<EBCellFAB>& a_coarseCorrection,
              bool                        a_skip_res) override final;

  /**
    @brief Prolongation onto AMR level
    @param[out] a_correction       Interpolated correction
    @param[in]  a_coarseCorrection Correction on coarse level
  */
  void
  AMRProlong(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_coarseCorrection) override final;

  /**
    @brief Update AMR residual
    @param[in] a_residual         Residual
    @param[in] a_correction       Correction
    @param[in] a_coarseCorrection Coarse-level correction
  */
  void
  AMRUpdateResidual(LevelData<EBCellFAB>&       a_residual,
                    const LevelData<EBCellFAB>& a_correction,
                    const LevelData<EBCellFAB>& a_coarseCorrection) override final;

  /**
    @brief Set alpha coefficient and beta coefficient (can change as diffusion solvers progress)
    @param[in] a_alpha Alpha-coefficient
    @param[in] a_beta Beta-coefficient
  */
  void
  setAlphaAndBeta(const Real& a_alpha, const Real& a_beta) override final;

  /**
    @brief Create coarsening of data holder
    @param[out] a_lhs    Coarsened data
    @param[in]  a_rhs    Fine data
    @param[in]  a_refRat Coarsening factor
  */
  void
  createCoarsened(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const int& a_refRat) override final;

  /**
    @brief Do diagonal scaling
    @param[in,out] a_rhs           Data to be scaled
    @param[in]    a_kappaWeighted Use kappa-weighted scaling or not
  */
  void
  diagonalScale(LevelData<EBCellFAB>& a_rhs, bool a_kappaWeighted) override final;

  /**
    @brief Divide by the a-coefficient
    @param[in,out] a_rhs Divided data
  */
  void
  divideByIdentityCoef(LevelData<EBCellFAB>& a_rhs) override final;

  /**
    @brief Apply operator but turn off all BCs.
    @param[in,out] a_Lphi L(phi) on output
    @param[in]     a_phi  Phi
  */
  void
  applyOpNoBoundary(LevelData<EBCellFAB>& a_Lphi, const LevelData<EBCellFAB>& a_phi) override final;

  /**
    @brief Not called, I think
    @param[in] a_phi Phi
  */
  void
  fillGrad(const LevelData<EBCellFAB>& a_phi) override final;

  /**
    @brief Fill flux
    @param[out] a_flux  Flux
    @param[in]  a_data  Data for which we will compute the flux
    @param[in]  a_grid  Grid
    @param[in]  a_dit   Corresponding data index. 
    @param[in]  a_scale Scaling factor
  */
  void
  getFlux(EBFluxFAB&                  a_flux,
          const LevelData<EBCellFAB>& a_data,
          const Box&                  a_grid,
          const DataIndex&            a_dit,
          Real                        a_scale) override final;

  /**
    @brief Allocate m_flux
  */
  void
  allocateFlux() const noexcept;

  /**
    @brief Deallocate m_flux
  */
  void
  deallocateFlux() const noexcept;

  /**
    @brief Returns m_flux. This is used in refluxing routines. 
    @return Flux
  */
  LevelData<EBFluxFAB>&
  getFlux() const;

  /**
    @brief Returns m_relCoef. Used by MFHelmholtzOp
    @return Relaxation coeff
  */
  const LevelData<EBCellFAB>&
  getRelaxationCoeff() const noexcept;

  /**
    @brief Returns m_chebyOrder. Used by MFHelmholtzOp to replicate Chebyshev parameters.
    @return Chebyshev polynomial degree
  */
  int
  getChebyOrder() const noexcept;

  /**
    @brief Returns m_chebyEigRatio. Used by MFHelmholtzOp to replicate Chebyshev parameters.
    @return Chebyshev eigenvalue ratio (lambda_max / lambda_min)
  */
  Real
  getChebyEigRatio() const noexcept;

  /**
    @brief Returns m_spectralRadius. Used by MFHelmholtzOp to replicate Chebyshev parameters.
    @return Gershgorin bound on rho(D^{-1}A)
  */
  Real
  getSpectralRadius() const noexcept;

  /**
    @brief Returns m_rasInnerSweeps. Used by MFHelmholtzOp to replicate the Schwarz block solve.
    @return Number of inner (frozen-ghost) red-black sweeps per outer Schwarz iteration
  */
  int
  getRasInnerSweeps() const noexcept;

  /**
    @brief Return all-cut VoFIterators (one per patch). Used by MFHelmholtzOp.
    @return All-cut VoF iterators over irregular cells on this level.
  */
  LayoutData<VoFIterator>&
  getVofIterIrreg() const noexcept;

protected:
  /**
    @brief Number of components that we solve for (always one..)
  */
  static constexpr int m_nComp = 1;

  /**
    @brief Component that we solve for
  */
  static constexpr int m_comp = 0;

  /**
    @brief Timer so user can profile.
  */
  Timer m_timer;

  /**
    @brief Interval
  */
  Interval m_interval;

  /**
    @brief Relaxation method
  */
  Smoother m_smoother;

  /**
    @brief Data centering
  */
  Location::Cell m_dataLocation;

  /**
    @brief Use reflux-free formulation or not
  */
  bool m_refluxFree;

  /**
    @brief Profile the operator
  */
  bool m_profile;

  /**
    @brief True if there is a multigrid level below this operator
  */
  bool m_hasMGObjects;

  /**
    @brief True if there's a finer level
  */
  bool m_hasFine;

  /**
    @brief True if there's a coarser level
  */
  bool m_hasCoar;

  /**
    @brief Refinement factor to coarse level
  */
  int m_refToCoar;

  /**
    @brief Refinement factor to fine level
  */
  int m_refToFine;

  /**
    @brief Number of smoothings in the preconditioner
  */
  int m_numSmoothPreCond;

  /**
    @brief Do coarse-fine interpolation or not. 
    @note Used by MFHelmholtzOp for optimizing how calls in EBHelmholtzOp are made
  */
  bool m_doInterpCF;

  /**
    @brief Turn on/off exchange operation.
    @note Used by MFHelmholtzOp for optimizing how calls in EBHelmholtzOp are made
  */
  bool m_doExchange;

  /**
    @brief Turn on/off exchange operation.
    @note Used by MFHelmholtzOp for optimizing how calls in EBHelmholtzOp are made
  */
  bool m_doCoarsen;

  /**
    @brief Ghost cells for phi
  */
  IntVect m_ghostPhi;

  /**
    @brief Ghost cells for rhs (note, the operator rhs)
  */
  IntVect m_ghostRhs;

  /**
    @brief Alpha-coefficient
  */
  Real m_alpha;

  /**
    @brief Beta-coefficient
  */
  Real m_beta;

  /**
    @brief Grid resolution;
  */
  Real m_dx;

  /**
    @brief Successive over-relaxation factor (Must be >= 1).
  */
  Real m_relaxFactor;

  /**
    @brief Number of Chebyshev polynomial steps per smoother invocation.
  */
  int m_chebyOrder;

  /**
    @brief Ratio lambda_max/lambda_min for the Chebyshev eigenvalue window.
    @details Targets the upper (1 - 1/m_chebyEigRatio) fraction of the spectrum.
  */
  Real m_chebyEigRatio;

  /**
    @brief Spectral radius estimate for the Jacobi-preconditioned operator D^{-1}A.
    @details Computed via the Gershgorin bound after defineStencils(). Used as lambda_max in the Chebyshev smoother.
  */
  Real m_spectralRadius;

  /**
    @brief Number of inner (frozen-ghost) red-black sweeps per outer restricted-additive-Schwarz iteration.
  */
  int m_rasInnerSweeps;

  /**
    @brief Vector resolution
  */
  RealVect m_vecDx;

  /**
    @brief Lower-left corner of domain
  */
  RealVect m_probLo;

  /**
    @brief Pre-built exchange copier
  */
  Copier m_exchangeCopier;

  /**
    @brief Pre-built exchange copier
  */
  Copier m_exchangeCopierFine;

  /**
    @brief Fine level grid (if the operator has a fine level)
  */
  EBLevelGrid m_eblgFine;

  /**
    @brief Grid
  */
  EBLevelGrid m_eblg;

  /**
    @brief Coarsened of m_eblg
  */
  EBLevelGrid m_eblgCoFi;

  /**
    @brief Coarse level grid (if the operator has a coarse level)
  */
  EBLevelGrid m_eblgCoar;

  /**
    @brief Coarser grids (multigrid level)
  */
  EBLevelGrid m_eblgCoarMG;

  /**
    @brief Restriction operator for AMR levels.
  */
  EBMGRestrict m_restrictOp;

  /**
    @brief Restriction operator if this is an MG level
  */
  EBMGRestrict m_restrictOpMG;

  /**
    @brief Prolongation operator for AMR levels
  */
  EBMGProlong m_prolongOp;

  /**
    @brief Prolongation operator for MG levels
  */
  EBMGProlong m_prolongOpMG;

  /**
    @brief Valid grid cells (might be nullptr on MG levels)
  */
  RefCountedPtr<LevelData<BaseFab<bool>>> m_validCells;

  /**
    @brief Domain bc object
  */
  RefCountedPtr<EBHelmholtzDomainBC> m_domainBc;

  /**
    @brief Domain bc object
  */
  RefCountedPtr<EBHelmholtzEBBC> m_ebBc;

  /**
    @brief Interpolator object
  */
  RefCountedPtr<EBMultigridInterpolator> m_interpolator;

  /**
    @brief Flux register
  */
  RefCountedPtr<EBReflux> m_fluxReg;

  /**
    @brief Conservative coarsener
  */
  RefCountedPtr<EBCoarAve> m_coarAve;

  /**
    @brief A-coefficient in Helmholtz equation
  */
  RefCountedPtr<LevelData<EBCellFAB>> m_Acoef;

  /**
    @brief B-coefficient in Helmholtz equation
  */
  RefCountedPtr<LevelData<EBFluxFAB>> m_Bcoef;

  /**
    @brief B-coefficient in Helmholtz equation, but on EB faces
  */
  RefCountedPtr<LevelData<BaseIVFAB<Real>>> m_BcoefIrreg;

  /**
    @brief Relaxation coefficient
  */
  LevelData<EBCellFAB> m_relCoef;

  /**
    @brief For holding fluxes
  */
  mutable LevelData<EBFluxFAB>* m_flux;

  /**
    @brief Interpolant for when we want centroid fluxes
  */
  LayoutData<BaseIFFAB<Real>> m_interpolant[SpaceDim];

  /**
    @brief Face centroid interpolation stencil
  */
  LayoutData<BaseIFFAB<FaceStencil>> m_interpStencil[SpaceDim];

  /**
    @brief Face centroid flux stencil. Defined on all faces connecting one or more irregular vofs
  */
  LayoutData<BaseIFFAB<VoFStencil>> m_centroidFluxStencil[SpaceDim];

  /**
    @brief Operator stencils in irregular cells (and ones that border irregular cells if using a centroid discretization).
    @details This stencil is => sum(fluxes)/dx, not including boundary faces or EB faces. I.e. this is the same as
    kappa*div(F) with the exclusion of faces where we have boundary conditions. 
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_relaxStencils;

  /**
    @brief For making irregular stencil applications go faster.
    @details This wraps m_relaxStencils in VCAggStencil (which computes explicit stencil offsets)
  */
  LayoutData<RefCountedPtr<VCAggStencil>> m_aggRelaxStencil;

  /**
    @brief Weights of diagonal alpha terms
  */
  LayoutData<BaseIVFAB<Real>> m_alphaDiagWeight;

  /**
    @brief Weights of diagonal beta terms
  */
  LayoutData<BaseIVFAB<Real>> m_betaDiagWeight;

  /**
    @brief VoFIterator for irregular cells
  */
  mutable LayoutData<VoFIterator> m_vofIterIrreg;

  /**
    @brief VoFIterator for "multi-cells". 
  */
  mutable LayoutData<VoFIterator> m_vofIterMulti;

  /**
    @brief VoFIterator which iterates over all cells that are 1) a cut-cell or 2) borders a cut-cell. 
  */
  mutable LayoutData<VoFIterator> m_vofIterStenc;

  /**
    @brief VoF iterators for lo domain side. 
  */
  mutable LayoutData<VoFIterator> m_vofIterDomLo[SpaceDim];

  /**
    @brief VoF iterators for hi domain side
  */
  mutable LayoutData<VoFIterator> m_vofIterDomHi[SpaceDim];

  /**
    @brief Domain boxes on each side
  */
  std::map<std::pair<int, Side::LoHiSide>, Box> m_sideBox;

  /**
    @brief Cached domain-boundary boxes for a single grid patch.
    @details Precomputed in EBHelmholtzOp::defineStencils so that applyDomainFlux and fillDomainFlux
             can avoid recomputing EBArith::loHi (and its associated timer/box-arithmetic overhead) on
             every V-cycle. The boxes depend only on the problem domain and the grid layout, both of
             which are fixed once the operator has been defined.
  */
  struct DomainBndryBoxes
  {
    /**
      @brief True if the patch abuts the domain boundary on any side, i.e. if any boundary-flux work is needed at all.
    */
    bool touchesDomain = false;

    /**
      @brief Per-direction flag: 1 if the patch abuts the lo domain face in that direction, 0 otherwise.
    */
    int hasLo[SpaceDim] = {D_DECL(0, 0, 0)};

    /**
      @brief Per-direction flag: 1 if the patch abuts the hi domain face in that direction, 0 otherwise.
    */
    int hasHi[SpaceDim] = {D_DECL(0, 0, 0)};

    /**
      @brief Lo-side cell-centered box abutting the domain. Valid only when hasLo[dir] == 1.
    */
    Box loBox[SpaceDim];

    /**
      @brief Hi-side cell-centered box abutting the domain. Valid only when hasHi[dir] == 1.
    */
    Box hiBox[SpaceDim];
  };

  /**
    @brief Precomputed domain-boundary boxes for each grid patch.
    @details Filled in EBHelmholtzOp::defineStencils and consumed by applyDomainFlux/fillDomainFlux.
  */
  LayoutData<DomainBndryBoxes> m_domainBndryBoxes;

  /**
    @brief "Colors" for the relaxation methods
  */
  Vector<IntVect> m_colors;

  /**
    @brief Jacobi relaxation
    @param[in,out] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of iterations
  */
  void
  relaxPointJacobi(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, const int a_iterations);

  /**
    @brief Jacobi relaxation
    @param[in,out] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of iterations
  */
  void
  relaxGSRedBlack(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, const int a_iterations);

  /**
    @brief Multi-colored gauss-seidel relaxation
    @param[in,out] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of iterations
  */
  void
  relaxGSMultiColor(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, const int a_iterations);

  /**
    @brief Chebyshev polynomial smoother.
    @details Applies m_chebyOrder Chebyshev-Richardson steps of the form
             correction += omega_i * relCoef * (residual - L(correction)).
             The step sizes omega_i are derived from Chebyshev nodes on [lambdaMax/m_chebyEigRatio, lambdaMax]
             where lambdaMax = 2 (safe bound for the Jacobi-preconditioned Laplacian-type operator).
    @param[in,out] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of outer smoother invocations (each applies m_chebyOrder steps)
  */
  void
  relaxChebyshev(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, const int a_iterations);

  /**
    @brief Restricted additive Schwarz smoother (block smoother).
    @details Each outer iteration exchanges ghost cells once and then performs m_rasInnerSweeps
             red-black Gauss-Seidel sweeps per patch with the ghost cells held frozen (Dirichlet).
             The blocks are the disjoint patches, so the update is automatically restricted to the
             valid region and no additive damping is required (overlap = 0). Because the ghost
             exchange is amortised over the inner sweeps, this smoother performs fewer halo exchanges
             per multigrid cycle than point relaxation.
    @param[in,out] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of outer Schwarz iterations
  */
  void
  relaxRestrictedAdditiveSchwarz(LevelData<EBCellFAB>&       a_correction,
                                 const LevelData<EBCellFAB>& a_residual,
                                 const int                   a_iterations);

  /**
    @brief Calculate the weight of the diagonal term
  */
  void
  computeDiagWeight();

  /**
    @brief Calculate relaxation coefficient
  */
  void
  computeRelaxationCoefficient();

  /**
    @brief Estimate the spectral radius of D^{-1}A (Jacobi-preconditioned operator) via the Gershgorin bound.
    @details rho_i = 2 - alpha * A_i * relCoef_i. Equals 2 for pure Laplacian and decreases with alpha > 0.
             Must be called after computeRelaxationCoefficient() and computeDiagWeight().
  */
  void
  computeSpectralRadius();

  /**
    @brief Compute aggregated stencils
  */
  void
  makeAggStencil();

  /**
    @brief Define stencils
  */
  void
  defineStencils();

  /**
    @brief Get the face-centered flux stencil
    @param[in]  a_face Face
    @param[in]  a_dit  Data index. Need because we multiply by B-coefficient
    @return Returns a stencil for the face-centered flux for the given face
    @note If the face is a boundary face the returned stencil is empty. 
  */
  VoFStencil
  getFaceCenterFluxStencil(const FaceIndex& a_face, const DataIndex& a_dit) const;

  /**
    @brief Get the face-centroid flux stencil
    @param[in]  a_face Face
    @param[in]  a_dit  Data index. Need because we multiply by B-coefficient
    @return Returns a stencil for the face-centroid flux for the given face
    @note The current version just interpolates the centered fluxes. 
  */
  VoFStencil
  getFaceCentroidFluxStencil(const FaceIndex& a_face, const DataIndex& a_dit) const;

  /**
    @brief Compute centroid fluxes in a direction.
    @param[out] a_fluxCentroid Flux on face centroid
    @param[in]  a_phi          Cell-centered data
    @param[in]  a_cellBox      Cell-centered compute box.
    @param[in]  a_dit          Data index
    @param[in]  a_dir          Coordinate direction
    @note This computes the flux for all dir-oriented faces in a_cellBox except for boundary faces.
  */
  void
  computeFlux(EBFaceFAB&       a_fluxCentroid,
              const EBCellFAB& a_phi,
              const Box&       a_cellBox,
              const DataIndex& a_dit,
              const int        a_dir);

  /**
    @brief Compute face-centered fluxes
    @param[out] a_fluxCenter Flux on center
    @param[in]  a_phi        Cell-centered data
    @param[in]  a_cellBox    Cell-centered compute box.
    @param[in]  a_dit        Data index
    @param[in]  a_dir        Coordinate direction
    @note This computes the flux for all dir-oriented faces in a_cellBox except for boundary faces.
  */
  void
  computeFaceCenteredFlux(EBFaceFAB&       a_fluxCenter,
                          const EBCellFAB& a_phi,
                          const Box&       a_cellBox,
                          const DataIndex& a_dit,
                          const int        a_dir);

  /**
    @brief Compute face-centroid fluxes
    @param[out] a_flux    Flux on face centroid
    @param[in]  a_phi     Cell-centered data
    @param[in]  a_cellBox Cell-centered compute box.
    @param[in]  a_dit     Data index
    @param[in]  a_dir     Coordinate direction
    @note This computes the flux for all dir-oriented faces in a_cellBox except for boundary faces.
  */
  void
  computeFaceCentroidFlux(EBFaceFAB&       a_flux,
                          const EBCellFAB& a_phi,
                          const Box&       a_cellBox,
                          const DataIndex& a_dit,
                          const int        a_dir);

  /**
    @brief Fill face centroid-centered fluxes on this level
    @param[in] a_phi Cell-centered data. Must have updated ghost cells for this to make any sense. 
    @note This computes the flux onto m_flux everywhere except on boundary faces. 
  */
  void
  computeFlux(const LevelData<EBCellFAB>& a_phi);

  /**
    @brief Reflux algorithm.
    @param[in,out] a_Lphi L(phiFine, phiCoar). This is modified along the CF interface to consistently account for fine-level fluxes. 
    @param[in]    a_phiFine Data on finer level
    @param[in]    a_phi     Data on this operator's level
    @param[in]    a_finerOp Operator on fine level
  */
  void
  reflux(LevelData<EBCellFAB>&             a_Lphi,
         const LevelData<EBCellFAB>&       a_phiFine,
         const LevelData<EBCellFAB>&       a_phi,
         AMRLevelOp<LevelData<EBCellFAB>>& a_finerOp);
};

#include <CD_NamespaceFooter.H>

#endif
