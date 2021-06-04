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

// Chombo includes
#include <ConductivityBaseDomainBC.H>
#include <BaseEBBC.H>
#include <EBMGAverage.H>
#include <EBMGInterp.H>
#include <LevelTGA.H>

// Our includes
#include <CD_EBMultigridInterpolator.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Factory class for making variable-coefficient Helmholtz operators. This is meant to be used with LevelTGA so the operator can be used for TGA time stepping. 
*/
class EBHelmholtzOp : public LevelTGAHelmOp<LevelData<EBCellFAB>, EBFluxFAB> {
public:

  using HelmholtzDomainBc = ConductivityBaseDomainBC;
  using HelmholtzEbBc     = ConductivityBaseEBBC;

  /*!
    @brief Relaxation method for the operators
  */
  enum class RelaxationMethod {
    Jacobi,
    GSRB,
    GSRBFast,
  };

  /*!
    @brief Disallowed default constructor
  */
#if 0 // Default, just trying to get this shit to compile
  EBHelmholtzOp() = delete;
#else
  EBHelmholtzOp(){

  }
#endif

  /*!
    @brief Disallowed copy constructor
  */
  EBHelmholtzOp(const EBHelmholtzOp& a_other) = delete;

  /*!
    @brief Full constructor
  */
  EBHelmholtzOp(const EBLevelGrid &                                  a_eblgFine,
		const EBLevelGrid &                                  a_eblg,
		const EBLevelGrid &                                  a_eblgCoar,
		const EBLevelGrid &                                  a_eblgCoarMG,
		const RefCountedPtr<EBMultigridInterpolator>&        a_quadCFI,
		const RefCountedPtr<HelmholtzDomainBc>&              a_domainBC,
		const RefCountedPtr<HelmholtzEbBc>&                  a_ebBC,
		const Real    &                                      a_dx,
		const Real    &                                      a_dxCoar,
		const int&                                           a_refToFine,
		const int&                                           a_refToCoar,
		const bool&                                          a_hasFine,
		const bool&                                          a_hasCoar,
		const bool&                                          a_hasMGObjects,
		const bool&                                          a_layoutChanged,
		const Real&                                          a_alpha,
		const Real&                                          a_beta,
		const RefCountedPtr<LevelData<EBCellFAB> >&          a_acoef,
		const RefCountedPtr<LevelData<EBFluxFAB> >&          a_bcoef,
		const RefCountedPtr<LevelData<BaseIVFAB<Real> > >&   a_bcoIrreg,
		const IntVect&                                       a_ghostCellsPhi,
		const IntVect&                                       a_ghostCellsRHS,
		const RelaxationMethod&                              a_relaxationMethod);

  /*!
    @brief Dtor
  */
  ~EBHelmholtzOp();

  /*!
    @brief Compute residual on this level. 
    @param[out] a_residual          Residual rhs - L(phi)
    @param[in]  a_phi               phi
    @param[in]  a_rhs               Right-hand side of system
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
  */
  void residual(LevelData<EBCellFAB>& a_residual, const LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_rhs, const bool a_homogeneousPhysBc);

  /*!
    @brief Precondition system before bottom solve
    @param[in] a_corr     Correction
    @param[in] a_residual Residual
    @details This just runs a few relaxations. 
  */
  void preCond(LevelData<EBCellFAB>& a_corr, const LevelData<EBCellFAB>& a_residual) override final;

  /*!
    @brief Apply operator 
    @param[out] a_Lphi              L(phi)
    @param[in]  a_phi               Phi
    @param[in]  a_homogeneousPhysBc Homogeneous physical BCs or not
    @details This computes a_Lphi = L(a_phi) using homogeneous physical BCs or not
  */
  void applyOp(LevelData<EBCellFAB>& a_Lphi, const LevelData<EBCellFAB>& a_phi, bool a_homogeneousPhysBc) override final;

  /*!
    @brief Create data which clones the layout of the other
    @param[out] a_lhs Data clone (returned data is not initialized)
    @param[out] a_rhs Data layout to be cloned
  */
  void create(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) override final;

  /*!
    @brief Assign data. 
    @details This does a local copy from rhs to lhs
  */
  void assign(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) override final;

  /*!
    @brief Compute the dot product??
  */
  Real dotProduct(const LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) override final;

  /*!
    @brief Increment operator
    @param[inout] a_lhs   Data to be incremented
    @param[in]    a_rhs   Incrementation data
    @param[in]    a_scale Scaling factor
  */
  void incr(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Real a_scale) override final;

  /*!
    @brief Set a_lhs = a*x + b*y
    @param[out] a_lhs Result data
    @param[in]  a_x   x-data
    @param[in]  a_y   y-data
    @param[in]  a_a   Scaling factor
    @param[in]  a_b   Scaling factor
  */
  void axby(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_x, const LevelData<EBCellFAB>& a_y, const Real a_a, const Real a_b) override final;

  /*!
    @brief Scale data. Returns a_lhs = a_lhs*a_scale
    @param[inout] a_lhs   Data to be scaled
    @param[i]     a_scale Scaling factor
  */
  void scale(LevelData<EBCellFAB>& a_lhs, const Real& a_scale) override final;

  /*!
    @brief Compute norm of data
    @param[in]
    @details a_order is currently ignored as we call a static routine in EBAMRPoissonOp for now. 
  */
  Real norm(const LevelData<EBCellFAB>& a_rhs, const int a_order) override final;

  /*!
    @brief Set data to zero
    @param[inout] a_lhs Data to be set to zero. 
  */
  void setToZero(LevelData<EBCellFAB>& a_lhs) override final;

  /*!
    @brief Create coarsened data
    @param[out] a_coarse  Coarse data
    @param[in]  a_fine    Fine data
    @param[in]  a_ghosted Include ghost cells or nto
  */
  void createCoarser(LevelData<EBCellFAB>& a_coarse, const LevelData<EBCellFAB>& a_fine, bool a_ghosted) override final;

  /*!
    @brief Relaxation method. This does smoothing for the system L(correction) = residual
    @param[inout] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of iterations
  */
  void relax(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, int a_iterations) override final;

  /*!
    @brief Restrict residual onto coarse level
    @param[inout] a_resCoar      Coarse residual
    @param[inout] a_phi Phi on this level
    @param[in]    a_rhs Rhs on this level
  */
  void restrictResidual(LevelData<EBCellFAB>& a_resCoar, LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_rhs) override final;

  /*!
    @brief Prolongation method. 
    @param[out] a_phi  Correction on this level
    @param[out] a_correctCoarse Correction on coarse level
  */
  void prolongIncrement(LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_correctCoarse) override final;

  /*!
    @brief Return coarsening factor to coarser level (1 if there is no coarser level);
  */
  int refToCoarser() override final;

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
  void AMRResidual(LevelData<EBCellFAB>&              a_residual,
		   const LevelData<EBCellFAB>&        a_phiFine,
		   const LevelData<EBCellFAB>&        a_phi,
		   const LevelData<EBCellFAB>&        a_phiCoar,
		   const LevelData<EBCellFAB>&        a_rhs,
		   bool                               a_homogeneousPhysBC,
		   AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp) override final;

  /*!
    @brief Compute AMR residual on finest AMR level. 
    @param[out] a_residual          Residual on this level
    @param[in]  a_phi               Phi on this level
    @param[in]  a_phiCoar           Phi on coar level
    @param[in]  a_rhs               Right-hand side on this level
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
  */
  void AMRResidualNF(LevelData<EBCellFAB>&              a_residual,
		     const LevelData<EBCellFAB>&        a_phi,
		     const LevelData<EBCellFAB>&        a_phiCoar,
		     const LevelData<EBCellFAB>&        a_rhs,
		     bool                               a_homogeneousPhysBC) override final;

  /*!
    @brief Compute AMR residual on coarsest
    @param[out] a_residual Residual on this level
    @param[in]  a_phiFine           Phi on fine level
    @param[in]  a_phi               Phi on this level
    @param[in]  a_rhs               Right-hand side on this level
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
    @param[in]  a_finerOp           Finer operator
  */
  void AMRResidualNC(LevelData<EBCellFAB>&              a_residual,
		     const LevelData<EBCellFAB>&        a_phiFine,
		     const LevelData<EBCellFAB>&        a_phi,
		     const LevelData<EBCellFAB>&        a_rhs,
		     bool                               a_homogeneousPhysBC,
		     AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp) override final;

  /*!
    @brief Apply the AMR operator, i.e. compute L(phi) in an AMR context, assuming no finer levels. 
    @param[out] a_Lphi              L(phi)
    @param[in]  a_phi               Phi on this level
    @param[in]  a_phiCoar           Phi on coar level
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
    @details This involves ghost cell interpolation if there's a coarse level, and refluxing if there's a fine level.
  */
  void AMROperatorNF(LevelData<EBCellFAB>&       a_Lphi,
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
  void AMROperatorNC(LevelData<EBCellFAB>&              a_Lphi,
		     const LevelData<EBCellFAB>&        a_phi,
		     const LevelData<EBCellFAB>&        a_phiCoar,
		     bool                               a_homogeneousPhysBC,
		     AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp) override final;

  /*!
    @brief Restrict residual
    @param[out] a_residualCoarse Coarse residual
    @param[out] a_residual         Residual
    @param[out] a_correction       Correction on this level 
    @param[out] a_coarseCorrection Coarse level correction
    @param[in]  a_skip_res         I have no idea what this one is supposed to do. 
  */
  void AMRRestrict(LevelData<EBCellFAB>&       a_residualCoarse,
		   const LevelData<EBCellFAB>& a_residual,
		   const LevelData<EBCellFAB>& a_correction,
		   const LevelData<EBCellFAB>& a_coarseCorrection,
		   bool                        a_skip_res) override final;

protected:

  /*!
    @brief Relaxation method
  */
  RelaxationMethod m_relaxationMethod;

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
    @brief Refinement factor to fine level
  */
  int m_refToFine;

  /*!
    @brief Refinement factor to fine level
  */
  int m_refToCoar;

  /*!
    @brief Fine level grid (if the operator has a fine level)
  */
  EBLevelGrid m_eblgFine;

  /*!
    @brief Grid
  */
  EBLevelGrid m_eblg;

  /*!
    @brief Coarse level grid (if the operator has a coarse level)
  */
  EBLevelGrid m_eblgCoar;

  /*!
    @brief Coarser grids (multigrid level)
  */
  EBLevelGrid m_eblgCoarMG;

  /*!
    @brief Averaging operator to MG level
  */
  EBMGAverage m_ebAverageMG;

  /*!
    @brief Averaging operator to AMR level
  */
  EBMGAverage m_ebAverage;

  /*!
    @brief Prolongation operator from MG level
  */
  EBMGInterp m_ebInterpMG;

  /*!
    @brief Prolongation operator from AMR level
  */
  EBMGInterp m_ebInterp;


  /*!
    @brief Apply operator on this level. This is a more general version which can turn on/off homogeneous and CF bcs.
    @param[out] a_Lphi              L(phi)
    @param[out] a_phi               Phi on this level
    @param[out] a_phiCoar           Coarse-level phi. If you have a coar this 
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
    @param[in]  a_homogeneousCFBC   Use homogeneous coarse-fine bcs or not
    @note This routine does a nasty cast of a_phi in order to work the ghost cells. 
  */
  void applyOp(LevelData<EBCellFAB>&             a_Lphi,
	       const LevelData<EBCellFAB>&       a_phi,
	       const LevelData<EBCellFAB>* const a_phiCoar,
	       const bool                        a_homogeneousPhysBC,
	       const bool                        a_homogeneousCFBC);

  /*!
    @brief Get factor for relaxation coefficient
    @details Returns 0.5 if using Jacobi and 1 otherwise. 
  */
  Real getSafety() const;

  /*!
    @brief Jacobi relaxation
    @param[inout] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of iterations
  */
  void relaxJacobi(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, const int a_iterations);

  /*!
    @brief Gauss-seidel relaxation
    @param[inout] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of iterations
  */
  void relaxGauSai(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, const int a_iterations);

  /*!
    @brief Fast Gauss-seidel relaxation
    @param[inout] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of iterations
  */
  void relaxGSRBFast(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, const int a_iterations);

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
  void AMROperator(LevelData<EBCellFAB>&              a_Lphi,
		   const LevelData<EBCellFAB>&        a_phiFine,
		   const LevelData<EBCellFAB>&        a_phi,
		   const LevelData<EBCellFAB>&        a_phiCoar,
		   const bool                         a_homogeneousPhysBC,
		   AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp);

  
};

#include <CD_NamespaceFooter.H>

#endif
