/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzOp.H
  @brief  Declaration of a class for solving multiphase Helmholtz equations. 
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzOp_H
#define CD_MFHelmholtzOp_H

// Std includes
#include <map>
#include <chrono>

// Chombo includes
#include <MFCellFAB.H>
#include <MFFluxFAB.H>

// Our includes
#include <CD_Location.H>
#include <CD_MFMultigridInterpolator.H>
#include <CD_MFCoarAve.H>
#include <CD_MFReflux.H>
#include <CD_MFLevelGrid.H>
#include <CD_MFBaseIVFAB.H>
#include <CD_MFHelmholtzJumpBC.H>
#include <CD_EBHelmholtzOp.H>
#include <CD_MFHelmholtzDomainBCFactory.H>
#include <CD_MFHelmholtzEBBCFactory.H>
#include <CD_MFHelmholtzJumpBCFactory.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Operator for solving multifluid Helmholtz on a grid level. 
  @details This operator should typically not be constructed directly by the user but used together with Chombo's AMRMultiGrid. In this
  case the user will want to constructor these operators through the factor.
*/
class MFHelmholtzOp : public LevelTGAHelmOp<LevelData<MFCellFAB>, MFFluxFAB>
{
public:
  /*!
    @brief Relaxation methods for this operator
  */
  enum class Smoother
  {
    PointJacobi,
    GauSaiRedBlack,
    GauSaiMultiColor,
  };

  /*!
    @brief Constructor. Must subsequently call define(...)
  */
  MFHelmholtzOp() = delete;

  /*!
    @brief No copy construction allowed
  */
  MFHelmholtzOp(const MFHelmholtzOp& a_op) = delete;

  /*!
    @brief No move construction allowed.
  */
  MFHelmholtzOp(const MFHelmholtzOp&& a_op) = delete;

  /*!
    @brief Full constructor
    @param[in] a_dataLocation  Data location, either cell center or cell centroid
    @param[in] a_mflgFine      Fine grids
    @param[in] a_mflg          Grids on this level
    @param[in] a_mflgCoFi      Coarsening of fine level grids
    @param[in] a_mflgCoar      Coarse grids
    @param[in] a_mflgCoarMG    Multigrid-grids
    @param[in] a_validCells    Valid grid cells (will be null on MG levels)
    @param[in] a_interpolator  Interpolator
    @param[in] a_fluxReg       Flux register
    @param[in] a_coarAve       Coarsener
    @param[in] a_domainBC      Domain BC
    @param[in] a_ebBC          Boundary conditions on EBs
    @param[in] a_jumpBcFactory Factory class for making jump BCs
    @param[in] a_probLo        Lower-left corner of computational domain
    @param[in] a_dx            Grid resolution
    @param[in] a_refToFine     Refinement ratio to fine level
    @param[in] a_refToCoar     Refinement ratio to coarse level
    @param[in] a_hasFine       Has fine level or not
    @param[in] a_hasCoar       Has coarse level or not
    @param[in] a_hasMGObjects  Has multigrid-objects (special objects between AMR levels, or below the AMR levels)
    @param[in] a_isMGOperator  Is MG operator or not
    @param[in] a_alpha         Operator alpha
    @param[in] a_beta          Operator beta
    @param[in] a_Acoef         Operator A-coefficient
    @param[in] a_Bcoef         Operator B-coefficient
    @param[in] a_BcoefIrreg    Operator B-coefficient (on EB faces)
    @param[in] a_ghostCellsPhi Ghost cells in data holders
    @param[in] a_ghostCellsPhi Ghost cells in data holders
    @param[in] a_jumpOrder     Stencil order to use on jump cells
    @param[in] a_jumpWeights   Weights for least squares stencils on jump cells
    @param[in] a_smoothPreCond Number of smoothings in the preconditioner
    @param[in] a_smoother      Which smoother to use
    @param[in] a_relaxFactor   Relaxation factor
  */
  MFHelmholtzOp(const Location::Cell                             a_dataLocation,
                const MFLevelGrid&                               a_mflgFine,
                const MFLevelGrid&                               a_mflg,
                const MFLevelGrid&                               a_mflgCoFi,
                const MFLevelGrid&                               a_mflgCoar,
                const MFLevelGrid&                               a_mflgCoarMG,
                const MFMultigridInterpolator&                   a_interpolator,
                const MFReflux&                                  a_fluxReg,
                const MFCoarAve&                                 a_coarAve,
                const RefCountedPtr<LevelData<BaseFab<bool>>>&   a_validCells,
                const RefCountedPtr<MFHelmholtzDomainBCFactory>& a_domainBC,
                const RefCountedPtr<MFHelmholtzEBBCFactory>&     a_ebBC,
                const RefCountedPtr<MFHelmholtzJumpBCFactory>&   a_jumpBcFactory,
                const RealVect&                                  a_probLo,
                const Real&                                      a_dx,
                const int&                                       a_refToFine,
                const int&                                       a_refToCoar,
                const bool&                                      a_hasFine,
                const bool&                                      a_hasCoar,
                const bool&                                      a_hasMGObjects,
                const bool&                                      a_isMGOperator,
                const Real&                                      a_alpha,
                const Real&                                      a_beta,
                const RefCountedPtr<LevelData<MFCellFAB>>&       a_Acoef,
                const RefCountedPtr<LevelData<MFFluxFAB>>&       a_Bcoef,
                const RefCountedPtr<LevelData<MFBaseIVFAB>>&     a_BcoefIrreg,
                const IntVect&                                   a_ghostPhi,
                const IntVect&                                   a_ghostRhs,
                const int&                                       a_jumpOrder,
                const int&                                       a_jumpWeight,
                const int&                                       a_preCondSmooth,
                const Smoother&                                  a_smoother,
                const Real&                                      a_relaxFactor);

  /*!
    @brief Destructor. 
  */
  virtual ~MFHelmholtzOp();

  /*!
    @brief No copy assigment allowed
  */
  MFHelmholtzOp&
  operator=(const MFHelmholtzOp& a_oper) = delete;

  /*!
    @brief No move assigment allowed
  */
  MFHelmholtzOp&
  operator=(const MFHelmholtzOp&& a_oper) = delete;

  /*!
    @brief Update operators with new coefficients.
    @param[in] a_Acoef         Operator A-coefficient
    @param[in] a_Bcoef         Operator B-coefficient
    @param[in] a_BcoefIrreg    Operator B-coefficient (on EB faces)
  */
  void
  setAcoAndBco(const RefCountedPtr<LevelData<MFCellFAB>>&   a_Acoef,
               const RefCountedPtr<LevelData<MFFluxFAB>>&   a_Bcoef,
               const RefCountedPtr<LevelData<MFBaseIVFAB>>& a_BcoefIrreg);

  /*!
    @brief Get the Helmholtz A-coefficient on cell centers
    @return m_Acoef
  */
  const RefCountedPtr<LevelData<MFCellFAB>>&
  getAcoef();

  /*!
    @brief Get the Helmholtz B-coefficient on faces
    @return m_Bcoef
  */
  const RefCountedPtr<LevelData<MFFluxFAB>>&
  getBcoef();

  /*!
    @brief Get the Helmholtz B-coefficient on the EB
    @return m_Bcoef
  */
  const RefCountedPtr<LevelData<MFBaseIVFAB>>&
  getBcoefIrreg();

  /*!
    @brief Set the jump boundary condition
  */
  void
  setJump(RefCountedPtr<LevelData<BaseIVFAB<Real>>>& a_jump);

  /*!
    @brief Get Helmholtz operator
  */
  RefCountedPtr<EBHelmholtzOp>&
  getHelmOp(const int a_phase);

  /*!
    @brief Return coarsening factor to coarser level (1 if there is no coarser level);
  */
  int
  refToCoarser() override final;

  /*!
    @brief Set alpha and beta
    @param[in] a_alpha  Alpha
    @param[in] a_beta   Beta
  */
  void
  setAlphaAndBeta(const Real& a_alpha, const Real& a_beta) override final;

  /*!
    @brief Divide by the a-coefficient
    @param[inout] a_rhs Divided data
  */
  void
  divideByIdentityCoef(LevelData<MFCellFAB>& a_rhs) override final;

  /*!
    @brief Apply operator but turn off all BCs. 
  */
  void
  applyOpNoBoundary(LevelData<MFCellFAB>& a_ans, const LevelData<MFCellFAB>& a_phi) override final;

  /*!
    @brief Not called, I think
    @param[in] a_phi Phi
  */
  void
  fillGrad(const LevelData<MFCellFAB>& a_phi) override final;

  /*!
    @brief Fill flux
    @param[out] a_flux  Flux
    @param[in]  a_data  Data for which we will compute the flux
    @param[in]  a_grid  Grid
    @param[in]  a_dit   Corresponding data index. 
    @param[in]  a_scale Scaling factor
  */
  void
  getFlux(MFFluxFAB&                  a_flux,
          const LevelData<MFCellFAB>& a_data,
          const Box&                  a_grid,
          const DataIndex&            a_dit,
          Real                        a_scale) override final;

  /*!
    @brief Compute residual on this level. 
    @param[out] a_residual          Residual rhs - L(phi)
    @param[in]  a_phi               phi
    @param[in]  a_rhs               Right-hand side of system
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
  */
  void
  residual(LevelData<MFCellFAB>&       a_residual,
           const LevelData<MFCellFAB>& a_phi,
           const LevelData<MFCellFAB>& a_rhs,
           const bool                  a_homogeneousPhysBc) override final;

  /*!
    @brief Precondition system before bottom solve
    @param[in] a_corr     Correction
    @param[in] a_residual Residual
    @details This just runs a few relaxations. 
  */
  void
  preCond(LevelData<MFCellFAB>& a_corr, const LevelData<MFCellFAB>& a_residual) override final;

  /*!
    @brief Apply operator 
    @param[out] a_Lphi              L(phi)
    @param[in]  a_phi               Phi
    @param[in]  a_homogeneousPhysBc Homogeneous physical BCs or not
    @details This computes a_Lphi = L(a_phi) using homogeneous physical BCs or not
  */
  void
  applyOp(LevelData<MFCellFAB>& a_Lphi, const LevelData<MFCellFAB>& a_phi, bool a_homogeneousPhysBc) override final;

  /*!
    @brief Time the applyOp routine. Template parameter is std::chrono duration. E.g. std::chrono::microseconds
    @param[in] a_phi Dummy cell-centered data. must have the correct number of ghost cells
    @param[in] a_numApply Number of times we apply the operator. 
    @return List of execution times in microseconds, ordered along the input DBLs boxArray()
  */
  template <typename Duration = std::chrono::microseconds>
  Vector<long long>
  computeOperatorLoads(LevelData<MFCellFAB>& a_phi, const int a_numApply) noexcept;

  /*!
    @brief Relaxation method. This does smoothing for the system L(correction) = residual
    @param[inout] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of iterations
  */
  void
  relax(LevelData<MFCellFAB>& a_correction, const LevelData<MFCellFAB>& a_residual, int a_iterations) override final;

  /*!
    @brief Jacobi relaxation
    @param[inout] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of iterations
  */
  void
  relaxPointJacobi(LevelData<MFCellFAB>& a_correction, const LevelData<MFCellFAB>& a_residual, const int a_iterations);

  /*!
    @brief Jacobi relaxation
    @param[inout] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of iterations
  */
  void
  relaxGSRedBlack(LevelData<MFCellFAB>& a_correction, const LevelData<MFCellFAB>& a_residual, const int a_iterations);

  /*!
    @brief Multi-colored gauss-seidel relaxation
    @param[inout] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of iterations
  */
  void
  relaxGSMultiColor(LevelData<MFCellFAB>& a_correction, const LevelData<MFCellFAB>& a_residual, const int a_iterations);

  /*!
    @brief Create method
    @param[out] a_lhs Clone
    @param[out] a_rhs Original data
  */
  void
  create(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs) override final;

  /*!
    @brief Create coarsened data
    @param[out] a_coarse  Coarse data
    @param[in]  a_fine    Fine data
    @param[in]  a_ghosted Include ghost cells or nto
  */
  void
  createCoarser(LevelData<MFCellFAB>& a_coarse, const LevelData<MFCellFAB>& a_fine, bool a_ghosted) override final;

  /*!
    @brief Create coarsening of data holder
    @param[out] a_lhs    Coarsened data
    @param[in]  a_rhs    Fine data
    @param[in]  a_refRat Coarsening factor
  */
  void
  createCoarsened(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs, const int& a_refRat) override final;

  /*!
    @brief Increment function
    @param[inout] a_lhs   On output, a_lhs = a_lhs + a_rhs*a_scale
    @param[in]    a_rhs   Data
    @param[in]    a_scale Scaling factor
  */
  void
  incr(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs, Real a_scale) override final;

  /*!
    @brief Dot product.
  */
  Real
  dotProduct(const LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_2) override final;

  /*!
    @brief Scale function
    @param[inout] a_lhs   On output, a_lhs = a_lhs*a_scale
    @param[in]    a_scale Scaling factor
  */
  void
  scale(LevelData<MFCellFAB>& a_lhs, const Real& a_scale) override final;

  /*!
    @brief Set to zero
    @param[out] a_lhs Data
  */
  void
  setToZero(LevelData<MFCellFAB>& a_lhs) override final;

  /*!
    @brief Assignment fucntion
    @param[out] a_lhs. Equal to a_rhs on output
    @param[in]  a_rhs. Data
  */
  void
  assign(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs) override final;

  /*!
    @brief Assign lhs
    @details This is the version that is called by AMRMultiGrid::VCycle. Note that the other version
    might be called by other operators (e.g., EBBackwardEuler)
    @param[out] a_lhs    Outgoing data
    @param[in]  a_rhs    Incoming data
    @param[in]  a_copier Copier
  */
  void
  assignCopier(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs, const Copier& a_copier) override final;

  /*!
    @brief Local assignment function
    @param[out] a_lhs. Equal to a_rhs on output
    @param[in]  a_rhs. Data
  */
  void
  assignLocal(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs) override final;

  /*!
    @brief Build copier
    @param[out] a_copier Copier for copying between a_lhs and a_rhs
    @param[in]  a_lhs    Copying from
    @param[in]  a_rhs    Copying to
  */
  void
  buildCopier(Copier& a_copier, const LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs) override;

  /*!
    @brief Compute solution norm. 
    @param[in] a_lhs   Data
    @param[in] a_order Norm order. Not used. 
  */
  Real
  norm(const LevelData<MFCellFAB>& a_lhs, int a_order) override final;

  /*!
    @brief Set a_lhs = a*x + b*y
    @param[out] a_lhs Result data
    @param[in]  a_x   x-data
    @param[in]  a_y   y-data
    @param[in]  a_a   Scaling factor
    @param[in]  a_b   Scaling factor
  */
  void
  axby(LevelData<MFCellFAB>&       a_lhs,
       const LevelData<MFCellFAB>& a_x,
       const LevelData<MFCellFAB>& a_y,
       const Real                  a_a,
       const Real                  a_b) override final;

  /*!
    @brief Restrict residual onto coarse level
    @param[inout] a_resCoar      Coarse residual
    @param[inout] a_phi Phi on this level
    @param[in]    a_rhs Rhs on this level
  */
  void
  restrictResidual(LevelData<MFCellFAB>&       a_resCoar,
                   LevelData<MFCellFAB>&       a_phi,
                   const LevelData<MFCellFAB>& a_rhs) override final;

  /*!
    @brief Prolongation method. 
    @param[out] a_phi  Correction on this level
    @param[out] a_correctCoarse Correction on coarse level
  */
  void
  prolongIncrement(LevelData<MFCellFAB>& a_phi, const LevelData<MFCellFAB>& a_correctCoarse) override final;

  /*!
    @brief Update AMR residual
    @param[in] a_residual         Residual
    @param[in] a_correction       Correction
    @param[in] a_coarseCorrection Coarse-level correction
  */
  void
  AMRUpdateResidual(LevelData<MFCellFAB>&       a_residual,
                    const LevelData<MFCellFAB>& a_correction,
                    const LevelData<MFCellFAB>& a_coarseCorrection) override final;

  /*!
    @brief Restrict residual
    @param[out] a_residualCoarse Coarse residual
    @param[out] a_residual         Residual
    @param[out] a_correction       Correction on this level 
    @param[out] a_coarseCorrection Coarse level correction
    @param[in]  a_skip_res         I have no idea what this one is supposed to do. 
  */
  void
  AMRRestrict(LevelData<MFCellFAB>&       a_residualCoarse,
              const LevelData<MFCellFAB>& a_residual,
              const LevelData<MFCellFAB>& a_correction,
              const LevelData<MFCellFAB>& a_coarseCorrection,
              bool                        a_skip_res) override final;

  /*!
    @brief Prolongation onto AMR level
    @param[out] a_correction       Interpolated correction
    @param[in]  a_coarseCorrection Correction on coarse level
  */
  void
  AMRProlong(LevelData<MFCellFAB>& a_correction, const LevelData<MFCellFAB>& a_coarseCorrection) override final;

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
  AMRResidual(LevelData<MFCellFAB>&             a_residual,
              const LevelData<MFCellFAB>&       a_phiFine,
              const LevelData<MFCellFAB>&       a_phi,
              const LevelData<MFCellFAB>&       a_phiCoar,
              const LevelData<MFCellFAB>&       a_rhs,
              bool                              a_homogeneousPhysBC,
              AMRLevelOp<LevelData<MFCellFAB>>* a_finerOp) override final;

  /*!
    @brief Compute AMR residual on finest AMR level. 
    @param[out] a_residual          Residual on this level
    @param[in]  a_phi               Phi on this level
    @param[in]  a_phiCoar           Phi on coar level
    @param[in]  a_rhs               Right-hand side on this level
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
  */
  void
  AMRResidualNF(LevelData<MFCellFAB>&       a_residual,
                const LevelData<MFCellFAB>& a_phi,
                const LevelData<MFCellFAB>& a_phiCoar,
                const LevelData<MFCellFAB>& a_rhs,
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
  AMRResidualNC(LevelData<MFCellFAB>&             a_residual,
                const LevelData<MFCellFAB>&       a_phiFine,
                const LevelData<MFCellFAB>&       a_phi,
                const LevelData<MFCellFAB>&       a_rhs,
                bool                              a_homogeneousPhysBC,
                AMRLevelOp<LevelData<MFCellFAB>>* a_finerOp) override final;

  /*!
    @brief Apply the AMR operator, i.e. compute L(phi) in an AMR context, assuming no finer levels. 
    @param[out] a_Lphi              L(phi)
    @param[in]  a_phi               Phi on this level
    @param[in]  a_phiCoar           Phi on coar level
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
    @details This involves ghost cell interpolation if there's a coarse level, and refluxing if there's a fine level.
  */
  void
  AMROperatorNF(LevelData<MFCellFAB>&       a_Lphi,
                const LevelData<MFCellFAB>& a_phi,
                const LevelData<MFCellFAB>& a_phiCoar,
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
  AMROperatorNC(LevelData<MFCellFAB>&             a_Lphi,
                const LevelData<MFCellFAB>&       a_phi,
                const LevelData<MFCellFAB>&       a_phiCoar,
                bool                              a_homogeneousPhysBC,
                AMRLevelOp<LevelData<MFCellFAB>>* a_finerOp) override final;

protected:
  /*!
    @brief Component that we solve for
  */
  static constexpr int m_comp = 0;

  /*!
    @brief Number of components that we solve for.
  */
  static constexpr int m_nComp = 1;

  /*!
    @brief Interpretation of data. Either on cell center or on cell centroid
  */
  Location::Cell m_dataLocation;

  /*!
    @brief Relaxation method
  */
  Smoother m_smoother;

  /*!
    @brief "Colors" for the multi-coloered relaxation method
  */
  Vector<IntVect> m_colors;

  /*!
    @brief Helmholtz operators on each phase. Note that I'm using int rather than Phase as identifier because that is the standard terminology for MFCellFAB. 
  */
  std::map<int, RefCountedPtr<EBHelmholtzOp>> m_helmOps;

  /*!
    @brief BC jump object. This is the one that has the stencils and can compute derivatives. 
  */
  RefCountedPtr<MFHelmholtzJumpBC> m_jumpBC;

  /*!
    @brief Actual BC jump in data-based format. This is the right-hand side of dphi/dn1 + dphi/dn2 = jump
  */
  RefCountedPtr<LevelData<BaseIVFAB<Real>>> m_jump;

  /*!
    @brief Valid grid cells (will be nullpointer on MG levels)
  */
  RefCountedPtr<LevelData<BaseFab<bool>>> m_validCells;

  /*!
    @brief Dirichlet BC values for each phase
  */
  std::map<int, RefCountedPtr<LevelData<BaseIVFAB<Real>>>> m_dirichletBcValues;

  /*!
    @brief Level grid
  */
  MFLevelGrid m_mflg;

  /*!
    @brief Fine level grid
  */
  MFLevelGrid m_mflgFine;

  /*!
    @brief Coarse grid
  */
  MFLevelGrid m_mflgCoar;

  /*!
    @brief Coarsened version of this grid
  */
  MFLevelGrid m_mflgCoFi;

  /*!
    @brief Coarse multigrid-grid
  */
  MFLevelGrid m_mflgCoarMG;

  /*!
    @brief Multi grid interpolator
  */
  MFMultigridInterpolator m_interpolator;

  /*!
    @brief Coarsener
  */
  MFCoarAve m_coarAve;

  /*!
    @brief Helmholtz A-coefficeint
  */
  RefCountedPtr<LevelData<MFCellFAB>> m_Acoef;

  /*!
    @brief Helmholtz B-coefficeint
  */
  RefCountedPtr<LevelData<MFFluxFAB>> m_Bcoef;

  /*!
    @brief Helmholtz B-coefficeint on EB
  */
  RefCountedPtr<LevelData<MFBaseIVFAB>> m_BcoefIrreg;

  /*!
    @brief Copier for exchange operation
  */
  Copier m_exchangeCopier;

  /*!
    @brief Multifluid operator or not. 
  */
  bool m_multifluid;

  /*!
    @brief Has MG objects or not.
  */
  bool m_hasMGObjects;

  /*!
    @brief Number of smoothings in the preconditioner
  */
  int m_numSmoothPreCond;

  /*!
    @brief Number of phases
  */
  int m_numPhases;

  /*!
    @brief Refinement factor to coarser AMR level
  */
  int m_refToCoar;

  /*!
    @brief Number of ghost cells
  */
  IntVect m_ghostPhi;

  /*!
    @brief Number of ghost cells
  */
  IntVect m_ghostRhs;

  /*!
    @brief True if there is a coarser AMR level
  */
  bool m_hasCoar;

  /*!
    @brief True if there is a finer AMR level
  */
  bool m_hasFine;

  /*!
    @brief True if there are multigrid levels
  */
  bool m_hasMGObjcts;

  /*!
    @brief Update the jump condition. 
    @param[in] a_phi Data
    @param[in] a_homogeneousPhysBC homogeneous physical BC or not
  */
  void
  updateJumpBC(const LevelData<MFCellFAB>& a_phi, const bool a_homogeneousPhysBC);

  /*!
    @brief Perform an exchange operation, event if the data is const.
    @param[in] a_phi Data 
  */
  void
  exchangeGhost(const LevelData<MFCellFAB>& a_phi) const;

  /*!
    @brief Do coarse-fine interpolation
    @param[inout] a_phi Fine-level data
    @param[in]    a_phiCoar Coarse-level data
    @param[in]    a_homogeneousCF Homogeneous interpolation or not. 
    @note Only ghost cells in a_phi are touched. 
  */
  void
  interpolateCF(const LevelData<MFCellFAB>& a_phi, const LevelData<MFCellFAB>* a_phiCoar, const bool a_homogeneousCF);

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
  applyOp(LevelData<MFCellFAB>&             a_Lphi,
          const LevelData<MFCellFAB>&       a_phi,
          const LevelData<MFCellFAB>* const a_phiCoar,
          const bool                        a_homogeneousPhysBC,
          const bool                        a_homogeneousCFBC);

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
  AMROperator(LevelData<MFCellFAB>&             a_Lphi,
              const LevelData<MFCellFAB>&       a_phiFine,
              const LevelData<MFCellFAB>&       a_phi,
              const LevelData<MFCellFAB>&       a_phiCoar,
              const bool                        a_homogeneousPhysBC,
              AMRLevelOp<LevelData<MFCellFAB>>* a_finerOp) override;
};

#include <CD_NamespaceFooter.H>

#include <CD_MFHelmholtzOpImplem.H>

#endif
