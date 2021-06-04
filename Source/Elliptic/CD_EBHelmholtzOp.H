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
    @param[out] a_residual Residual rhs - L(phi)
    @param[in]  a_phi      phi
    @param[in]  a_rhs      Right-hand side of system
    @param[in]  a_homogeneous
  */
  void residual(LevelData<EBCellFAB>& a_residual, const LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_rhs, bool a_homogeneous) override final;

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

protected:

  /*!
    @brief Relaxation method
  */
  RelaxationMethod m_relaxationMethod;

  /*!
    @brief Grid
  */
  EBLevelGrid m_eblg;

  /*!
    @brief Get factor for relaxation coefficient
    @details Returns 0.5 if using Jacobi and 1 otherwise. 
  */
  Real getSafety() const;

};

#include <CD_NamespaceFooter.H>

#endif
