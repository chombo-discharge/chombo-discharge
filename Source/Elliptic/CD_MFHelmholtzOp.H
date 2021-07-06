/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzOp.H
  @brief  Declaration of a factory class for making multifluid Helmholtz operators for multigrid. 
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzOp_H
#define CD_MFHelmholtzOp_H

// Std includes
#include <map>

// Chombo includes
#include <MFCellFAB.H>
#include <MFFluxFAB.H>

// Our includes
#include <CD_MFMultigridInterpolator.H>
#include <CD_MFCoarAve.H>
#include <CD_MFFluxReg.H>
#include <CD_MFLevelGrid.H>
#include <CD_MFBaseIVFAB.H>
#include <CD_JumpBC.H>
#include <CD_EBHelmholtzOp.H>
#include <CD_MFHelmholtzDomainBCFactory.H>
#include <CD_MFHelmholtzEBBCFactory.H>
#include <CD_NamespaceHeader.H>

class MFHelmholtzOp : public AMRLevelOp<LevelData<MFCellFAB> > {
public:

  /*!
    @brief Relaxation methods for this operator
  */
  enum class RelaxationMethod {
    PointJacobi,
    GauSaiRedBlack,
    GauSaiMultiColor,
  };

  /*!
    @brief Constructor. Must subsequently call define(...)
  */
  MFHelmholtzOp();

  /*!
    @brief Full constructor
  */
  MFHelmholtzOp(const MFLevelGrid&                               a_mflgFine,
		const MFLevelGrid&                               a_mflg,
		const MFLevelGrid&                               a_mflgCoFi,		
		const MFLevelGrid&                               a_mflgCoar,
		const MFLevelGrid&                               a_mflgCoarMG,
		const MFMultigridInterpolator&                   a_interpolator,
		const MFFluxReg&                                 a_fluxReg,
		const MFCoarAve&                                 a_coarAve,
		const RefCountedPtr<MFHelmholtzDomainBCFactory>& a_domainBcFactory,
		const RefCountedPtr<MFHelmholtzEBBCFactory>&     a_ebBcFactory,
		const RealVect&                                  a_probLo,
		const Real&                                      a_dx,
		const int&                                       a_refToCoar,
		const bool&                                      a_hasFine,
		const bool&                                      a_hasCoar,		
		const bool&                                      a_hasMGObjects,
		const bool&                                      a_isMGOperator,
		const Real&                                      a_alpha,
		const Real&                                      a_beta,
		const RefCountedPtr<LevelData<MFCellFAB> >&      a_Acoef,
		const RefCountedPtr<LevelData<MFFluxFAB> >&      a_Bcoef,
		const RefCountedPtr<LevelData<MFBaseIVFAB> >&    a_BcoefIrreg,
		const IntVect&                                   a_ghostPhi,
		const IntVect&                                   a_ghostRhs,
		const int&                                       a_jumpOrder,
		const int&                                       a_jumpWeight,
		const RelaxationMethod&                          a_relaxType);

  /*!
    @brief Destructor. 
  */
  ~MFHelmholtzOp();

  /*!
    @brief Set the jump boundary condition
  */
  void setJump(const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_jump);

  /*!
    @brief Get Helmholtz operator
  */
  RefCountedPtr<EBHelmholtzOp>& getHelmOp(const int a_phase);

  /*!
    @brief Return coarsening factor to coarser level (1 if there is no coarser level);
  */
  int refToCoarser() override final;

  /*!
    @brief Compute residual on this level. 
    @param[out] a_residual          Residual rhs - L(phi)
    @param[in]  a_phi               phi
    @param[in]  a_rhs               Right-hand side of system
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
  */
  void residual(LevelData<MFCellFAB>& a_residual, const LevelData<MFCellFAB>& a_phi, const LevelData<MFCellFAB>& a_rhs, const bool a_homogeneousPhysBc);

  /*!
    @brief Precondition system before bottom solve
    @param[in] a_corr     Correction
    @param[in] a_residual Residual
    @details This just runs a few relaxations. 
  */
  void preCond(LevelData<MFCellFAB>& a_corr, const LevelData<MFCellFAB>& a_residual) override final;

  /*!
    @brief Apply operator 
    @param[out] a_Lphi              L(phi)
    @param[in]  a_phi               Phi
    @param[in]  a_homogeneousPhysBc Homogeneous physical BCs or not
    @details This computes a_Lphi = L(a_phi) using homogeneous physical BCs or not
  */
  void applyOp(LevelData<MFCellFAB>& a_Lphi, const LevelData<MFCellFAB>& a_phi, bool a_homogeneousPhysBc) override final;

  /*!
    @brief Relaxation method. This does smoothing for the system L(correction) = residual
    @param[inout] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of iterations
  */
  void relax(LevelData<MFCellFAB>& a_correction, const LevelData<MFCellFAB>& a_residual, int a_iterations) override final;

  /*!
    @brief Jacobi relaxation
    @param[inout] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of iterations
  */
  void relaxPointJacobi(LevelData<MFCellFAB>& a_correction, const LevelData<MFCellFAB>& a_residual, const int a_iterations);

  /*!
    @brief Jacobi relaxation
    @param[inout] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of iterations
  */
  void relaxGSRedBlack(LevelData<MFCellFAB>& a_correction, const LevelData<MFCellFAB>& a_residual, const int a_iterations);

  /*!
    @brief Multi-colored gauss-seidel relaxation
    @param[inout] a_correction Correction
    @param[in]    a_residual   Residual
    @param[in]    a_iterations Number of iterations
  */
  void relaxGSMultiColor(LevelData<MFCellFAB>& a_correction, const LevelData<MFCellFAB>& a_residual, const int a_iterations);

  /*!
    @brief Create method
    @param[out] a_lhs Clone
    @param[out] a_rhs Original data
  */
  void create(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs) override final;

  /*!
    @brief Create coarsened data
    @param[out] a_coarse  Coarse data
    @param[in]  a_fine    Fine data
    @param[in]  a_ghosted Include ghost cells or nto
  */
  void createCoarser(LevelData<MFCellFAB>& a_coarse, const LevelData<MFCellFAB>& a_fine, bool a_ghosted) override final;

  /*!
    @brief Create coarsening of data holder
    @param[out] a_lhs    Coarsened data
    @param[in]  a_rhs    Fine data
    @param[in]  a_refRat Coarsening factor
  */
  void createCoarsened(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs, const int& a_refRat) override final;
  
  /*!
    @brief Increment function
    @param[inout] a_lhs   On output, a_lhs = a_lhs + a_rhs*a_scale
    @param[in]    a_rhs   Data
    @param[in]    a_scale Scaling factor
  */
  void incr(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs, Real a_scale) override final;

  /*!
    @brief Dot product.
  */
  Real dotProduct(const LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_2) override final;

  /*!
    @brief Scale function
    @param[inout] a_lhs   On output, a_lhs = a_lhs*a_scale
    @param[in]    a_scale Scaling factor
  */
  void scale(LevelData<MFCellFAB>& a_lhs, const Real& a_scale) override final;

  /*!
    @brief Set to zero
    @param[out] a_lhs Data
  */
  void setToZero(LevelData<MFCellFAB>& a_lhs) override final;

  /*!
    @brief Assignment fucntion
    @param[out] a_lhs. Equal to a_rhs on output
    @param[in]  a_rhs. Data
  */
  void assign(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs) override final;

  /*!
    @brief Compute solution norm. 
    @param[in] a_lhs   Data
    @param[in] a_order Norm order. Not used. 
  */
  Real norm(const LevelData<MFCellFAB>& a_lhs, int a_order) override final;

  /*!
    @brief Set a_lhs = a*x + b*y
    @param[out] a_lhs Result data
    @param[in]  a_x   x-data
    @param[in]  a_y   y-data
    @param[in]  a_a   Scaling factor
    @param[in]  a_b   Scaling factor
  */
  void axby(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_x, const LevelData<MFCellFAB>& a_y, const Real a_a, const Real a_b) override final;

  /*!
    @brief Restrict residual onto coarse level
    @param[inout] a_resCoar      Coarse residual
    @param[inout] a_phi Phi on this level
    @param[in]    a_rhs Rhs on this level
  */
  void restrictResidual(LevelData<MFCellFAB>& a_resCoar, LevelData<MFCellFAB>& a_phi, const LevelData<MFCellFAB>& a_rhs) override final;

  /*!
    @brief Prolongation method. 
    @param[out] a_phi  Correction on this level
    @param[out] a_correctCoarse Correction on coarse level
  */
  void prolongIncrement(LevelData<MFCellFAB>& a_phi, const LevelData<MFCellFAB>& a_correctCoarse) override final;

  /*!
    @brief Update AMR residual
    @param[in] a_residual         Residual
    @param[in] a_correction       Correction
    @param[in] a_coarseCorrection Coarse-level correction
  */
  void AMRUpdateResidual(LevelData<MFCellFAB>&       a_residual,
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
  void AMRRestrict(LevelData<MFCellFAB>&       a_residualCoarse,
		   const LevelData<MFCellFAB>& a_residual,
		   const LevelData<MFCellFAB>& a_correction,
		   const LevelData<MFCellFAB>& a_coarseCorrection,
		   bool                        a_skip_res) override final;

  /*!
    @brief Prolongation onto AMR level
    @param[out] a_correction       Interpolated correction
    @param[in]  a_coarseCorrection Correction on coarse level
  */
  void AMRProlong(LevelData<MFCellFAB>& a_correction, const LevelData<MFCellFAB>& a_coarseCorrection) override final;

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
  void AMRResidual(LevelData<MFCellFAB>&              a_residual,
		   const LevelData<MFCellFAB>&        a_phiFine,
		   const LevelData<MFCellFAB>&        a_phi,
		   const LevelData<MFCellFAB>&        a_phiCoar,
		   const LevelData<MFCellFAB>&        a_rhs,
		   bool                               a_homogeneousPhysBC,
		   AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp) override final;

  /*!
    @brief Compute AMR residual on finest AMR level. 
    @param[out] a_residual          Residual on this level
    @param[in]  a_phi               Phi on this level
    @param[in]  a_phiCoar           Phi on coar level
    @param[in]  a_rhs               Right-hand side on this level
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
  */
  void AMRResidualNF(LevelData<MFCellFAB>&              a_residual,
		     const LevelData<MFCellFAB>&        a_phi,
		     const LevelData<MFCellFAB>&        a_phiCoar,
		     const LevelData<MFCellFAB>&        a_rhs,
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
  void AMRResidualNC(LevelData<MFCellFAB>&              a_residual,
		     const LevelData<MFCellFAB>&        a_phiFine,
		     const LevelData<MFCellFAB>&        a_phi,
		     const LevelData<MFCellFAB>&        a_rhs,
		     bool                               a_homogeneousPhysBC,
		     AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp) override final;

  /*!
    @brief Apply the AMR operator, i.e. compute L(phi) in an AMR context, assuming no finer levels. 
    @param[out] a_Lphi              L(phi)
    @param[in]  a_phi               Phi on this level
    @param[in]  a_phiCoar           Phi on coar level
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
    @details This involves ghost cell interpolation if there's a coarse level, and refluxing if there's a fine level.
  */
  void AMROperatorNF(LevelData<MFCellFAB>&       a_Lphi,
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
  void AMROperatorNC(LevelData<MFCellFAB>&              a_Lphi,
		     const LevelData<MFCellFAB>&        a_phi,
		     const LevelData<MFCellFAB>&        a_phiCoar,
		     bool                               a_homogeneousPhysBC,
		     AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp) override final;
  
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
    @brief Relaxation method
  */
  RelaxationMethod m_relaxType;

  /*!
    @brief "Colors" for the multi-coloered relaxation method
  */
  Vector<IntVect> m_colors;

  /*!
    @brief Helmholtz operators on each phase. Note that I'm using int rather than Phase as identifier because that is the standard terminology for MFCellFAB. 
  */
  std::map<int, RefCountedPtr<EBHelmholtzOp> > m_helmOps;

  /*!
    @brief BC jump object. This is the one that has the stencils and can compute derivatives. 
  */
  RefCountedPtr<JumpBC> m_jumpBC;

  /*!
    @brief Actual BC jump in data-based format. This is the right-hand side of dphi/dn1 + dphi/dn2 = jump
  */
  RefCountedPtr<LevelData<BaseIVFAB<Real> > > m_jump;

  /*!
    @brief Dirichlet BC values for each phase
  */
  std::map<int, RefCountedPtr<LevelData<BaseIVFAB<Real> > > > m_dirichletBcValues;

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
    @brief Multifluid operator or not. 
  */
  bool m_multifluid;

  /*!
    @brief Has MG objects or not.
  */
  bool m_hasMGObjects;

  /*!
    @brief Number of phases
  */
  int m_numPhases;

  /*!
    @brief Refinement factor to coarser AMR level
  */
  int m_refToCoar;

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
    @note Note sure how to do this function just yet. 
  */
  void updateJumpBC(const LevelData<MFCellFAB>& a_phi, const bool a_homogeneousPhysBC);

  /*!
    @brief Do coarse-fine interpolation
    @param[inout] a_phi Fine-level data
    @param[in]    a_phiCoar Coarse-level data
    @param[in]    a_homogeneousCF Homogeneous interpolation or not. 
    @note Only ghost cells in a_phi are touched. 
  */
  void interpolateCF(const LevelData<MFCellFAB>& a_phi, const LevelData<MFCellFAB>* a_phiCoar, const bool a_homogeneousCF);

  /*!
    @brief Apply operator on this level. This is a more general version which can turn on/off homogeneous and CF bcs.
    @param[out] a_Lphi              L(phi)
    @param[out] a_phi               Phi on this level
    @param[out] a_phiCoar           Coarse-level phi. If you have a coar this 
    @param[in]  a_homogeneousPhysBC Use homogeneous physical BCs or not
    @param[in]  a_homogeneousCFBC   Use homogeneous coarse-fine bcs or not
    @note a_phi is not really const and we do a nasty cast. Leaving it as const because we only modify the ghost cells!
  */
  void applyOp(LevelData<MFCellFAB>&             a_Lphi,
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
  void AMROperator(LevelData<MFCellFAB>&              a_Lphi,
		   const LevelData<MFCellFAB>&        a_phiFine,
		   const LevelData<MFCellFAB>&        a_phi,
		   const LevelData<MFCellFAB>&        a_phiCoar,
		   const bool                         a_homogeneousPhysBC,
		   AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp);

};

#include <CD_NamespaceFooter.H>

#endif
