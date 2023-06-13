/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBHelmholtzOp.H
  @brief  Declaration of Helmholtz multigrid operators. 
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzOp_H
#define CD_EBHelmholtzOp_H

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

/*!
  @brief Helmholtz operator for equations like alpha*a(x)*phi(x) + beta*div(b(x)*grad(phi(x))) = rho. 
  @details This can be used with TGA time stepping. 
*/
class EBHelmholtzOp : public LevelTGAHelmOp<LevelData<EBCellFAB>, EBFluxFAB>
{
public:
  /*!
    @brief Relaxation method for the operators
  */
  enum class Smoother
  {
    NoRelax,
    PointJacobi,
    GauSaiRedBlack,
    GauSaiMultiColor,
  };

  /*!
    @brief Disallowed default constructor
  */
  EBHelmholtzOp() = delete;

  /*!
    @brief Disallowed copy constructor
    @param[in] a_other Other operator
  */
  EBHelmholtzOp(const EBHelmholtzOp& a_other) = delete;

  /*!
    @brief Disallowed move constructor
    @param[in] a_other Other operator
  */
  EBHelmholtzOp(const EBHelmholtzOp&& a_other) = delete;

  /*!
    @brief Full constructor
    @param[in] a_dataLocation  Data location, either cell center or cell centroid
    @param[in] a_eblgFine      Fine grids
    @param[in] a_eblg          Grids on this level
    @param[in] a_eblgCoFi      Coarsening of fine level grids
    @param[in] a_eblgCoar      Coarse grids
    @param[in] a_eblgCoarMG    Multigrid-grids
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
    @param[in] a_BcoefIrreg    Operator B-coefficient (on EB faces)
    @param[in] a_ghostCellsPhi Ghost cells in data holders
    @param[in] a_ghostCellsPhi Ghost cells in data holders
    @param[in] a_smoother      Which smoother to use
  */
  EBHelmholtzOp(const Location::Cell                             a_dataLocation,
                const EBLevelGrid&                               a_eblgFine,
                const EBLevelGrid&                               a_eblg,
                const EBLevelGrid&                               a_eblgCoFi,
                const EBLevelGrid&                               a_eblgCoar,
                const EBLevelGrid&                               a_eblgCoarMG,
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
                const RefCountedPtr<LevelData<BaseIVFAB<Real>>>& a_BcoIrreg,
                const IntVect&                                   a_ghostCellsPhi,
                const IntVect&                                   a_ghostCellsRHS,
                const Smoother&                                  a_smoother);

  /*!
    @brief Dtor
  */
  virtual ~EBHelmholtzOp();

  /*!
    @brief No copy assigment allowed
    @param[in] a_other Other operator
  */
  EBHelmholtzOp&
  operator=(const EBHelmholtzOp& a_oper) = delete;

  /*!
    @brief No move assigment allowed
    @param[in] a_other Other operator
  */
  EBHelmholtzOp&
  operator=(const EBHelmholtzOp&& a_oper) = delete;

  /*!
    @brief Turn off BCs
  */
  void
  turnOffCFInterp();

  /*!
    @brief Turn on BCs
  */
  void
  turnOnCFInterp();

  /*!
    @brief Turn off exchange operation
    @note Done when MFHelmholtzOp does the exchange and we don't want EBHelmholtzOp to redo it.
  */
  void
  turnOffExchange();

  /*!
    @brief Turn on exchange operation
    @note Done when MFHelmholtzOp does the exchange and we don't want EBHelmholtzOp to redo it.
  */
  void
  turnOnExchange();

  /*!
    @brief Turn off coarsening operation
    @note Done when MFHelmholtzOp does the coarsening and we don't want EBHelmholtzOp to redo it.
  */
  void
  turnOffCoarsening();

  /*!
    @brief Turn on coarsening operation
    @note Done when MFHelmholtzOp does the coarsening and we don't want EBHelmholtzOp to redo it.
  */
  void
  turnOnCoarsening();

  /*!
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

  /*!
    @brief Get the Helmholtz A-coefficient on cell centers
    @return m_Acoef
  */
  const RefCountedPtr<LevelData<EBCellFAB>>&
  getAcoef();

  /*!
    @brief Get the Helmholtz B-coefficient on faces
    @return m_Bcoef
  */
  const RefCountedPtr<LevelData<EBFluxFAB>>&
  getBcoef();

  /*!
    @brief Get the Helmholtz B-coefficient on the EB
    @return m_Bcoef
  */
  const RefCountedPtr<LevelData<BaseIVFAB<Real>>>&
  getBcoefIrreg();

  /*!
    @brief Coarsen data from fine to coar level
    @param[in] a_phi     Data on this level
    @param[in] a_phiFine Data on finer level
    @note This routine is used in AMROperator to ensure that EB fluxes which are described by stencils get consistent data under the fine level. 
  */
  void
  coarsenCell(LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_phiFine);

  /*!
    @brief Coarsen fluxes on the fine level onto this level. 
    @details User must update fluxes before calling this routine. 
    @param[in] a_flux     Flux on this level. 
    @param[in] a_fineFlux Flux on finer level
    @note This routine is used in the reflux-free AMROperator
  */
  void
  coarsenFlux(LevelData<EBFluxFAB>& a_flux, const LevelData<EBFluxFAB>& a_fineFlux);

  /*!
    @brief Point Jacobi kernel. 
    @param[inout] a_Lcorr Storage for computing L(a_corr)
    @param[inout] a_corr    Correction
    @param[in]    a_resid   Residual
    @param[in]    a_cellBox Grid box
    @param[in]    a_dit     Data index
  */
  void
  pointJacobiKernel(EBCellFAB&       a_Lcorr,
                    EBCellFAB&       a_corr,
                    const EBCellFAB& a_resid,
                    const Box&       a_cellBox,
                    const DataIndex& a_dit);

  /*!
    @brief Red-black Gauss-Seidel kernel
    @param[inout] a_Lcorr    Storage for computing L(a_corr)
    @param[inout] a_corr     Correction
    @param[in]    a_resid    Residual
    @param[in]    a_redBlack Red or black
    @param[in]    a_cellBox  Grid box
    @param[in]    a_dit      Data index
  */
  void
  gauSaiRedBlackKernel(EBCellFAB&       a_Lcorr,
                       EBCellFAB&       a_corr,
                       const EBCellFAB& a_resid,
                       const Box&       a_cellBox,
                       const DataIndex& a_dit,
                       const int&       a_redBlack);

  /*!
    @brief Multi-color Gauss-Seidel kernel
    @param[inout] a_Lcorr    Storage for computing L(a_corr)
    @param[inout] a_corr     Correction
    @param[in]    a_resid    Residual
    @param[in]    a_color    Color
  */
  void
  gauSaiMultiColorKernel(EBCellFAB&       a_Lcorr,
                         EBCellFAB&       a_corr,
                         const EBCellFAB& a_resid,
                         const Box&       a_cellBox,
                         const DataIndex& a_dit,
                         const IntVect&   a_color);

  /*!
    @brief Compute residual on this level. 
    @param[out] a_residual          Residual rhs - L(phi)
    @param[in]  a_phi               phi
    @param[in]  a_rhs               Right-hand side of system
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
  */
  void
  residual(LevelData<EBCellFAB>&       a_residual,
           const LevelData<EBCellFAB>& a_phi,
           const LevelData<EBCellFAB>& a_rhs,
           const bool                  a_homogeneousPhysBc);

  /*!
    @brief Precondition system before bottom solve
    @param[in] a_corr     Correction
    @param[in] a_residual Residual
    @details This just runs a few relaxations. 
  */
  void
  preCond(LevelData<EBCellFAB>& a_corr, const LevelData<EBCellFAB>& a_residual) override final;

  /*!
    @brief Apply coarse-fine boundary conditions
    @param[inout] a_phiFine         Fine data
    @param[inout] a_phiCoar         Coarse data
    @param[in]    a_homogeneousCFBC Homogeneous CF or not (i.e. coar data is zero);
  */
  void
  interpolateCF(LevelData<EBCellFAB>& a_phiFine, const LevelData<EBCellFAB>* a_phiCoar, const bool a_homogeneousCFBC);

  /*!
    @brief Do homogeneous coarse-fine interpolation
    @param[in] a_phi Data
  */
  void
  homogeneousCFInterp(LevelData<EBCellFAB>& a_phi);

  /*!
    @brief Inhomogeneous coarse-fine interpolation
    @param[inout] a_phiFine Fine data. Ghost cells will be filled.
    @param[inout] a_phiCoar Coarse data. 
  */
  void
  inhomogeneousCFInterp(LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_phiCoar);

  /*!
    @brief Multi-colored Gauss-Seidel kernel. Public because MFHelmholtzOp may want to use use. 
    @param[in] a_phi Correction
    @param[in] a_phi Lphi L(corrcetion)
    @param[in] a_rhs Residual
    @param[in] a_color "Color":
  */
  void
  gauSaiMultiColor(LevelData<EBCellFAB>&       a_phi,
                   const LevelData<EBCellFAB>& a_Lphi,
                   const LevelData<EBCellFAB>& a_rhs,
                   const IntVect&              a_color) const;

  /*!
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

  /*!
    @brief Apply operator 
    @param[out] a_Lphi              L(phi)
    @param[in]  a_phi               Phi
    @param[in]  a_homogeneousPhysBc Homogeneous physical BCs or not
    @details This computes a_Lphi = L(a_phi) using homogeneous physical BCs or not
  */
  void
  applyOp(LevelData<EBCellFAB>& a_Lphi, const LevelData<EBCellFAB>& a_phi, bool a_homogeneousPhysBc) override final;

  /*!
    @brief Apply operator in a grid box. 
    @param[out] a_Lphi              L(phi)
    @param[in]  a_phi               Phi
    @param[in]  a_cellBox           Grid box
    @param[in]  a_dit               Data indxe
    @param[in]  a_homogeneousPhysBC Homogeneous physical BCs or not
  */
  void
  applyOp(EBCellFAB&       a_Lphi,
          EBCellFAB&       a_phi,
          const Box&       a_cellBox,
          const DataIndex& a_dit,
          const bool       a_homogeneousPhysBC);

  /*!
    @brief Apply operator in regular cells.
    @param[out] a_Lphi              L(phi)
    @param[in]  a_phi               Phi
    @param[in]  a_cellBox           Grid box
    @param[in]  a_dit               Data indxe
    @param[in]  a_homogeneousPhysBC Homogeneous physical BCs or not
  */
  void
  applyOpRegular(EBCellFAB&       a_Lphi,
                 EBCellFAB&       a_phi,
                 const Box&       a_cellBox,
                 const DataIndex& a_dit,
                 const bool       a_homogeneousPhysBC);

  /*!
    @brief Apply domain flux. 
    @param[inout] a_phi               Cell data
    @param[in]    a_cellBox           Computation box
    @param[in]    a_dit               Data index
    @param[in]    a_homogeneousPhysBC Homogeneous BC or not
    @note This sets the ghost cells of phi to phi(i-1) = phi(i) + F*dx/bco such that when we take the centered difference we get
    -bco*(phi(i) - phi(i-1))/dx = F. 
  */
  void
  applyDomainFlux(EBCellFAB& a_phi, const Box& a_cellBox, const DataIndex& a_dit, const bool a_homogeneousPhysBc);

  /*!
    @brief Fill domain flux. This fills the flux on the domain face using centered differencing ala applyDomainFlux
    @details a_flux is replaced by the user-specified flux on the domain faces. 
    @param[inout] a_flux              Flux data holder.
    @param[inout] a_phi               Cell-centered data
    @param[in]    a_cellBox           Computation box
    @param[in]    a_dit               Data index
  */
  void
  fillDomainFlux(EBFluxFAB& a_flux, const EBCellFAB& a_phi, const Box& a_cellBox, const DataIndex& a_dit);

  /*!
    @brief Apply operator in irregular cells
    @param[out] a_Lphi              L(phi)
    @param[in]  a_phi               Cell-centered data
    @param[in]  a_cellBox           Cell box
    @param[in]  a_dit               Grid index
    @param[in]  a_homogeneousPhysBC Homogeneous physical BCs or not
  */
  void
  applyOpIrregular(EBCellFAB&       a_Lphi,
                   const EBCellFAB& a_phi,
                   const Box&       a_cellBox,
                   const DataIndex& a_dit,
                   const bool       a_homogeneousPhysBC);

  /*!
    @brief Create data which clones the layout of the other
    @param[out] a_lhs Data clone (returned data is not initialized)
    @param[out] a_rhs Data layout to be cloned
  */
  void
  create(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) override final;

  /*!
    @brief Assign data.
    @details This does a copy from rhs to lhs
  */
  void
  assign(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) override final;

  /*!
    @brief Assign lhs
    @details This is the version that is called by AMRMultiGrid::VCycle. Note that the other version
    might be called by other operators (e.g., EBBackwardEuler)
    @param[out] a_lhs    Outgoing data
    @param[in]  a_rhs    Incoming data
    @param[in]  a_copier Copier
  */
  void
  assignCopier(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Copier& a_copier) override final;

  /*!
    @brief Local assignment function
    @param[out] a_lhs. Equal to a_rhs on output
    @param[in]  a_rhs. Data
  */
  void
  assignLocal(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) override final;

  /*!
    @brief Build copier
    @param[out] a_copier Copier for copying between a_lhs and a_rhs
    @param[in]  a_lhs    Copying from
    @param[in]  a_rhs    Copying to
  */
  void
  buildCopier(Copier& a_copier, const LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) override;

  /*!
    @brief Compute the dot product??
  */
  Real
  dotProduct(const LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) override final;

  /*!
    @brief Increment operator
    @param[inout] a_lhs   Data to be incremented
    @param[in]    a_rhs   Incrementation data
    @param[in]    a_scale Scaling factor
  */
  void
  incr(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Real a_scale) override final;

  /*!
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

  /*!
    @brief Scale data. Returns a_lhs = a_lhs*a_scale
    @param[inout] a_lhs   Data to be scaled
    @param[i]     a_scale Scaling factor
  */
  void
  scale(LevelData<EBCellFAB>& a_lhs, const Real& a_scale) override final;

  /*!
    @brief Compute norm of data
    @param[in]
    @details a_order is currently ignored as we call a static routine in EBAMRPoissonOp for now. 
  */
  Real
  norm(const LevelData<EBCellFAB>& a_rhs, const int a_order) override final;

  /*!
    @brief Set data to zero
    @param[inout] a_lhs Data to be set to zero. 
  */
  void
  setToZero(LevelData<EBCellFAB>& a_lhs) override final;

  /*!
    @brief Create coarsened data
    @param[out] a_coarse  Coarse data
    @param[in]  a_fine    Fine data
    @param[in]  a_ghosted Include ghost cells or nto
  */
  void
  createCoarser(LevelData<EBCellFAB>& a_coarse, const LevelData<EBCellFAB>& a_fine, bool a_ghosted) override final;

  /*!
    @brief Relaxation method. This does smoothing for the system L(correction) = residual
    @param[inout] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of iterations
  */
  void
  relax(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, int a_iterations) override final;

  /*!
    @brief Restrict residual onto coarse level
    @param[inout] a_resCoar      Coarse residual
    @param[inout] a_phi Phi on this level
    @param[in]    a_rhs Rhs on this level
  */
  void
  restrictResidual(LevelData<EBCellFAB>&       a_resCoar,
                   LevelData<EBCellFAB>&       a_phi,
                   const LevelData<EBCellFAB>& a_rhs) override final;

  /*!
    @brief Prolongation method. 
    @param[out] a_phi  Correction on this level
    @param[out] a_correctCoarse Correction on coarse level
  */
  void
  prolongIncrement(LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_correctCoarse) override final;

  /*!
    @brief Return coarsening factor to coarser level (1 if there is no coarser level);
  */
  int
  refToCoarser() override final;

  /*!
    @brief Apply the AMR operator, i.e. compute L(phi) in an AMR context. 
    @param[out] a_residual          Residual on this level
    @param[in]  a_phiFine           Phi on fine level
    @param[in]  a_phi               Phi on this level
    @param[in]  a_phiCoar           Phi on coar level
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
    @param[in]  a_finerOp           Finer operatator
    @details This involves ghost cell interpolation if there's a coarse level, and refluxing if there's a fine level.
  */
  void
  AMROperator(LevelData<EBCellFAB>&             a_Lphi,
              const LevelData<EBCellFAB>&       a_phiFine,
              const LevelData<EBCellFAB>&       a_phi,
              const LevelData<EBCellFAB>&       a_phiCoar,
              const bool                        a_homogeneousPhysBC,
              AMRLevelOp<LevelData<EBCellFAB>>* a_finerOp);

  /*!
    @brief Apply the AMR operator, i.e. compute L(phi) in an AMR context. 
    @details This is the reflux-free version of AMROperator.
    @param[out] a_residual          Residual on this level
    @param[in]  a_phiFine           Phi on fine level
    @param[in]  a_phi               Phi on this level
    @param[in]  a_phiCoar           Phi on coar level
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
    @param[in]  a_finerOp           Finer operatator
    @details This involves ghost cell interpolation if there's a coarse level, and refluxing if there's a fine level.
  */
  void
  refluxFreeAMROperator(LevelData<EBCellFAB>&             a_Lphi,
                        const LevelData<EBCellFAB>&       a_phiFine,
                        const LevelData<EBCellFAB>&       a_phi,
                        const LevelData<EBCellFAB>&       a_phiCoar,
                        const bool                        a_homogeneousPhysBC,
                        AMRLevelOp<LevelData<EBCellFAB>>* a_finerOp);

  /*!
    @brief Compute residual on this level. AMR version. 
    @param[out] a_residual          Residual on this level
    @param[in]  a_phiFine           Phi on fine level
    @param[in]  a_phi               Phi on this level
    @param[in]  a_phiCoar           Phi on coar level
    @param[in]  a_rhs               Right-hand side on this level
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
    @param[in]  a_finerOp           Finer operatator
  */
  void
  AMRResidual(LevelData<EBCellFAB>&             a_residual,
              const LevelData<EBCellFAB>&       a_phiFine,
              const LevelData<EBCellFAB>&       a_phi,
              const LevelData<EBCellFAB>&       a_phiCoar,
              const LevelData<EBCellFAB>&       a_rhs,
              bool                              a_homogeneousPhysBC,
              AMRLevelOp<LevelData<EBCellFAB>>* a_finerOp) override final;

  /*!
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

  /*!
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

  /*!
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

  /*!
    @brief Apply the AMR operator, i.e. compute L(phi) in an AMR context, assuming no coarser AMR levels. 
    @param[out] a_Lphi              L(phi)
    @param[in]  a_phiFine           Phi on finer level
    @param[in]  a_phi               Phi on this level
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
    @details This involves ghost cell interpolation if there's a coarse level, and refluxing if there's a fine level.
  */
  void
  AMROperatorNC(LevelData<EBCellFAB>&             a_Lphi,
                const LevelData<EBCellFAB>&       a_phi,
                const LevelData<EBCellFAB>&       a_phiCoar,
                bool                              a_homogeneousPhysBC,
                AMRLevelOp<LevelData<EBCellFAB>>* a_finerOp) override final;

  /*!
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

  /*!
    @brief Prolongation onto AMR level
    @param[out] a_correction       Interpolated correction
    @param[in]  a_coarseCorrection Correction on coarse level
  */
  void
  AMRProlong(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_coarseCorrection) override final;

  /*!
    @brief Update AMR residual
    @param[in] a_residual         Residual
    @param[in] a_correction       Correction
    @param[in] a_coarseCorrection Coarse-level correction
  */
  void
  AMRUpdateResidual(LevelData<EBCellFAB>&       a_residual,
                    const LevelData<EBCellFAB>& a_correction,
                    const LevelData<EBCellFAB>& a_coarseCorrection) override final;

  /*!
    @brief Set alpha coefficient and beta coefficient (can change as diffusion solvers progress)
    @param[in] a_alpha Alpha-coefficient
    @param[in] a_beta Beta-coefficient
  */
  void
  setAlphaAndBeta(const Real& a_alpha, const Real& a_beta) override final;

  /*!
    @brief Create coarsening of data holder
    @param[out] a_lhs    Coarsened data
    @param[in]  a_rhs    Fine data
    @param[in]  a_refRat Coarsening factor
  */
  void
  createCoarsened(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const int& a_refRat) override final;

  /*!
    @brief Do diagonal scaling
    @param[inout] a_rhs           Data to be scaled
    @param[in]    a_kappaWeighted Use kappa-weighted scaling or not
  */
  void
  diagonalScale(LevelData<EBCellFAB>& a_rhs, bool a_kappaWeighted);

  /*!
    @brief Divide by the a-coefficient
    @param[inout] a_rhs Divided data
  */
  void
  divideByIdentityCoef(LevelData<EBCellFAB>& a_rhs) override final;

  /*!
    @brief Apply operator but turn off all BCs. 
  */
  void
  applyOpNoBoundary(LevelData<EBCellFAB>& a_ans, const LevelData<EBCellFAB>& a_phi) override final;

  /*!
    @brief Not called, I think
    @param[in] a_phi Phi
  */
  void
  fillGrad(const LevelData<EBCellFAB>& a_phi) override final;

  /*!
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

  /*!
    @brief Allocate m_flux
  */
  void
  allocateFlux() const noexcept;

  /*!
    @brief Deallocate m_flux
  */
  void
  deallocateFlux() const noexcept;

  /*!
    @brief Returns m_flux. This is used in refluxing routines. 
  */
  LevelData<EBFluxFAB>&
  getFlux() const;

protected:
  /*!
    @brief Number of components that we solve for (always one..)
  */
  static constexpr int m_nComp = 1;

  /*!
    @brief Component that we solve for
  */
  static constexpr int m_comp = 0;

  /*!
    @brief Timer so user can profile.
  */
  Timer m_timer;

  /*!
    @brief Interval
  */
  Interval m_interval;

  /*!
    @brief Relaxation method
  */
  Smoother m_smoother;

  /*!
    @brief Data centering
  */
  Location::Cell m_dataLocation;

  /*!
    @brief Use reflux-free formulation or not
  */
  bool m_refluxFree;

  /*!
    @brief Profile the operator
  */
  bool m_profile;

  /*!
    @brief True if there is a multigrid level below this operator
  */
  bool m_hasMGObjects;

  /*!
    @brief True if there's a finer level
  */
  bool m_hasFine;

  /*!
    @brief True if there's a coarser level
  */
  bool m_hasCoar;

  /*!
    @brief Refinement factor to coarse level
  */
  int m_refToCoar;

  /*!
    @brief Refinement factor to fine level
  */
  int m_refToFine;

  /*!
    @brief Do coarse-fine interpolation or not. 
    @note Used by MFHelmholtzOp for optimizing how calls in EBHelmholtzOp are made
  */
  bool m_doInterpCF;

  /*!
    @brief Turn on/off exchange operation.
    @note Used by MFHelmholtzOp for optimizing how calls in EBHelmholtzOp are made
  */
  bool m_doExchange;

  /*!
    @brief Turn on/off exchange operation.
    @note Used by MFHelmholtzOp for optimizing how calls in EBHelmholtzOp are made
  */
  bool m_doCoarsen;

  /*!
    @brief Ghost cells for phi
  */
  IntVect m_ghostPhi;

  /*!
    @brief Ghost cells for rhs (note, the operator rhs)
  */
  IntVect m_ghostRhs;

  /*!
    @brief Alpha-coefficient
  */
  Real m_alpha;

  /*!
    @brief Beta-coefficient
  */
  Real m_beta;

  /*!
    @brief Grid resolution;
  */
  Real m_dx;

  /*!
    @brief Vector resolution
  */
  RealVect m_vecDx;

  /*!
    @brief Lower-left corner of domain
  */
  RealVect m_probLo;

  /*!
    @brief Fine level grid (if the operator has a fine level)
  */
  EBLevelGrid m_eblgFine;

  /*!
    @brief Grid
  */
  EBLevelGrid m_eblg;

  /*!
    @brief Coarsened of m_eblg
  */
  EBLevelGrid m_eblgCoFi;

  /*!
    @brief Coarse level grid (if the operator has a coarse level)
  */
  EBLevelGrid m_eblgCoar;

  /*!
    @brief Coarser grids (multigrid level)
  */
  EBLevelGrid m_eblgCoarMG;

  /*!
    @brief Restriction operator for AMR levels.
  */
  EBMGRestrict m_restrictOp;

  /*!
    @brief Restriction operator if this is an MG level
  */
  EBMGRestrict m_restrictOpMG;

  /*!
    @brief Prolongation operator for AMR levels
  */
  EBMGProlong m_prolongOp;

  /*!
    @brief Prolongation operator for MG levels
  */
  EBMGProlong m_prolongOpMG;

  /*!
    @brief Domain bc object
  */
  RefCountedPtr<EBHelmholtzDomainBC> m_domainBc;

  /*!
    @brief Domain bc object
  */
  RefCountedPtr<EBHelmholtzEBBC> m_ebBc;

  /*!
    @brief Interpolator object
  */
  RefCountedPtr<EBMultigridInterpolator> m_interpolator;

  /*!
    @brief Flux register
  */
  RefCountedPtr<EBReflux> m_fluxReg;

  /*!
    @brief Conservative coarsener
  */
  RefCountedPtr<EBCoarAve> m_coarAve;

  /*!
    @brief A-coefficient in Helmholtz equation
  */
  RefCountedPtr<LevelData<EBCellFAB>> m_Acoef;

  /*!
    @brief B-coefficient in Helmholtz equation
  */
  RefCountedPtr<LevelData<EBFluxFAB>> m_Bcoef;

  /*!
    @brief B-coefficient in Helmholtz equation, but on EB faces
  */
  RefCountedPtr<LevelData<BaseIVFAB<Real>>> m_BcoefIrreg;

  /*!
    @brief Relaxation coefficient
  */
  LevelData<EBCellFAB> m_relCoef;

  /*!
    @brief For holding fluxes
  */
  mutable LevelData<EBFluxFAB>* m_flux;

  /*!
    @brief Interpolant for when we want centroid fluxes
  */
  LayoutData<BaseIFFAB<Real>> m_interpolant[SpaceDim];

  /*!
    @brief Face centroid interpolation stencil
  */
  LayoutData<BaseIFFAB<FaceStencil>> m_interpStencil[SpaceDim];

  /*!
    @brief Face centroid flux stencil. Defined on all faces connecting one or more irregular vofs
  */
  LayoutData<BaseIFFAB<VoFStencil>> m_centroidFluxStencil[SpaceDim];

  /*!
    @brief Operator stencils in irregular cells (and ones that border irregular cells if using a centroid discretization).
    @details This stencil is => sum(fluxes)/dx, not including boundary faces or EB faces. I.e. this is the same as
    kappa*div(F) with the exclusion of faces where we have boundary conditions. 
  */
  LayoutData<BaseIVFAB<VoFStencil>> m_relaxStencils;

  /*!
    @brief For making irregular stencil applications go faster.
    @details This wraps m_relaxStencils in VCAggStencil (which computes explicit stencil offsets)
  */
  LayoutData<RefCountedPtr<VCAggStencil>> m_aggRelaxStencil;

  /*!
    @brief Weights of diagonal alpha terms
  */
  LayoutData<BaseIVFAB<Real>> m_alphaDiagWeight;

  /*!
    @brief Weights of diagonal beta terms
  */
  LayoutData<BaseIVFAB<Real>> m_betaDiagWeight;

  /*!
    @brief VoFIterator for irregular cells
  */
  mutable LayoutData<VoFIterator> m_vofIterIrreg;

  /*!
    @brief VoFIterator for "multi-cells". 
  */
  mutable LayoutData<VoFIterator> m_vofIterMulti;

  /*!
    @brief VoFIterator which iterates over all cells that are 1) a cut-cell or 2) borders a cut-cell. 
  */
  mutable LayoutData<VoFIterator> m_vofIterStenc;

  /*!
    @brief VoF iterators for lo domain side. 
  */
  mutable LayoutData<VoFIterator> m_vofIterDomLo[SpaceDim];

  /*!
    @brief VoF iterators for hi domain side
  */
  mutable LayoutData<VoFIterator> m_vofIterDomHi[SpaceDim];

  /*!
    @brief Domain boxes on each side
  */
  std::map<std::pair<int, Side::LoHiSide>, Box> m_sideBox;

  /*!
    @brief "Colors" for the relaxation methods
  */
  Vector<IntVect> m_colors;

  /*!
    @brief Jacobi relaxation
    @param[inout] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of iterations
  */
  void
  relaxPointJacobi(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, const int a_iterations);

  /*!
    @brief Jacobi relaxation
    @param[inout] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of iterations
  */
  void
  relaxGSRedBlack(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, const int a_iterations);

  /*!
    @brief Multi-colored gauss-seidel relaxation
    @param[inout] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of iterations
  */
  void
  relaxGSMultiColor(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, const int a_iterations);

  /*!
    @brief Calculate the weight of the diagonal term
  */
  void
  computeDiagWeight();

  /*!
    @brief Calculate relaxation coefficient
  */
  void
  computeRelaxationCoefficient();

  /*!
    @brief Compute aggregated stencils
  */
  void
  makeAggStencil();

  /*!
    @brief Define stencils
  */
  void
  defineStencils();

  /*!
    @brief Get the face-centered flux stencil
    @param[in]  a_face Face
    @param[in]  a_dit  Data index. Need because we multiply by B-coefficient
    @return Returns a stencil for the face-centered flux for the given face
    @note If the face is a boundary face the returned stencil is empty. 
  */
  VoFStencil
  getFaceCenterFluxStencil(const FaceIndex& a_face, const DataIndex& a_dit) const;

  /*!
    @brief Get the face-centroid flux stencil
    @param[in]  a_face Face
    @param[in]  a_dit  Data index. Need because we multiply by B-coefficient
    @return Returns a stencil for the face-centroid flux for the given face
    @note The current version just interpolates the centered fluxes. 
  */
  VoFStencil
  getFaceCentroidFluxStencil(const FaceIndex& a_face, const DataIndex& a_dit) const;

  /*!
    @brief Compute centroid fluxes in a direction.
    @param[out] a_fluxCentroid   Flux on face centroid
    @param[in   a_phi          Cell-centered data
    @param[in]  a_cellBox      Cell-centered compute box.
    @param[in]  a_dit          Data index
    @note This computes the flux for all dir-oriented faces in a_cellBox except for boundary faces. 
  */
  void
  computeFlux(EBFaceFAB&       a_fluxCentroid,
              const EBCellFAB& a_phi,
              const Box&       a_cellBox,
              const DataIndex& a_dit,
              const int        a_dir);

  /*!
    @brief Compute face-centered fluxes
    @param[out] a_fluxCenter   Flux on center
    @param[in   a_phi          Cell-centered data
    @param[in]  a_cellBox      Cell-centered compute box.
    @param[in]  a_dit          Data index
    @note This computes the flux for all dir-oriented faces in a_cellBox except for boundary faces. 
  */
  void
  computeFaceCenteredFlux(EBFaceFAB&       a_fluxCenter,
                          const EBCellFAB& a_phi,
                          const Box&       a_cellBox,
                          const DataIndex& a_dit,
                          const int        a_dir);

  /*!
    @brief Compute face-centered fluxes
    @param[out] a_fluxCenter   Flux on center
    @param[in   a_phi          Cell-centered data
    @param[in]  a_cellBox      Cell-centered compute box.
    @param[in]  a_dit          Data index
    @note This computes the flux for all dir-oriented faces in a_cellBox except for boundary faces. 
  */
  void
  computeFaceCentroidFlux(EBFaceFAB&       a_flux,
                          const EBCellFAB& a_phi,
                          const Box&       a_cellBox,
                          const DataIndex& a_dit,
                          const int        a_dir);

  /*!
    @brief Fill face centroid-centered fluxes on this level
    @param[in] a_phi Cell-centered data. Must have updated ghost cells for this to make any sense. 
    @note This computes the flux onto m_flux everywhere except on boundary faces. 
  */
  void
  computeFlux(const LevelData<EBCellFAB>& a_phi);

  /*!
    @brief Reflux algorithm.
    @param[inout] a_Lphi L(phiFine, phiCoar). This is modified along the CF interface to consistently account for fine-level fluxes. 
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
