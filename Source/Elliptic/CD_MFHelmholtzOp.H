/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzOp.H
  @brief  Declaration of a factory class for making multifluid Helmholtz operators for multigrid. 
  @author Robert Marskar
n*/

#ifndef CD_MFHelmholtzOp_H
#define CD_MFHelmholtzOp_H

// Std includes
#include <map>

// Chombo includes
#include <MFCellFAB.H>
#include <MFFluxFAB.H>

// Our includes
#include <CD_Phases.H>
#include <CD_MFBaseIVFAB.H>
#include <CD_EBHelmholtzOp.H>
#include <CD_EBHelmholtzDomainBCFactory.H>
#include <CD_EBHelmholtzEBBCFactory.H>
#include <CD_NamespaceHeader.H>

class MFHelmholtzOp : public AMRLevelOp<LevelData<MFCellFAB> > {
public:

  /*!
    @brief Constructor. Must subsequently call define(...)
  */
  MFHelmholtzOp();

  /*!
    @brief Destructor. 
  */
  ~MFHelmholtzOp();

  /*!
    @brief Define function. Not implemented (yet). 
  */
  void define();

  /*!
    @brief Return coarsening factor to coarser level (1 if there is no coarser level);
  */
  int refToCoarser() override final;

  /*!
    @brief Weird but true, coarsening between VCycles give better convergence
  */
  unsigned int orderOfAccuracy(void) const override final;

  /*!
    @brief Enforce CF consistency
    @note This does conservative averaging on each phase. 
  */
  void enforceCFConsistency(LevelData<MFCellFAB>& a_coarCorr, const LevelData<MFCellFAB>& a_fineCorr) override final;

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
    @brief Create method
    @param[out] a_lhs Clone
    @param[out] a_rhs Original data
  */
  void create(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs) override final;

  /*!
    @brief Increment function
    @param[inout] a_lhs   On output, a_lhs = a_lhs + a_rhs*a_scale
    @param[in]    a_rhs   Data
    @param[in]    a_scale Scaling factor
  */
  void incr(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs, Real a_scale) override final;

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
  
protected:

  /*!
    @brief Number of aliases needed. 
  */
  static constexpr int m_numAlias = 6;

  /*!
    @brief Helmholtz operators on each phase. Note that I'm using int rather than Phase as identifier because that is the standard terminology for MFCellFAB. 
  */
  std::map<int, RefCountedPtr<EBHelmholtzOp> > m_helmOps;

  /*!
    @brief Alias data
  */
  std::vector<LevelData<EBCellFAB>* > m_alias;

  /*!
    @brief Arihmetic multifluid operators
  */
  LevelDataOps<MFCellFAB> m_ops;

  /*!
    @brief Refinement factor to coarser AMR level
  */
  int m_refToCoar;
};

#include <CD_NamespaceFooter.H>

#endif
