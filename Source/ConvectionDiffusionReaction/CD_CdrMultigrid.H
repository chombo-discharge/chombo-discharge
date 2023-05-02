/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrMultigrid.H
  @brief  Extensions of CdrSolver which use EBHelmholtzOp for diffusion solves
  @author Robert Marskar
*/

#ifndef CD_CdrMultigrid_H
#define CD_CdrMultigrid_H

// Std includes
#include <random>
#include <time.h>
#include <chrono>

// Chombo includes
#include <AMRMultiGrid.H>
#include <BiCGStabSolver.H>
#include <EBSimpleSolver.H>
#include <GMRESSolver.H>

// Our includes
#include <CD_EBHelmholtzOpFactory.H>
#include <CD_CdrSolver.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Extension class of CdrSolver that uses multigrid for diffusion part. Can also solve for stochastic diffusion.
*/
class CdrMultigrid : public CdrSolver
{
protected:
  /*!
    @brief Enum class for supported bottom solvers in multigrid.
  */
  enum class BottomSolverType
  {
    Simple,
    BiCGStab,
    GMRES
  };

  /*!
    @brief Enum for multigrid cycle types. 
  */
  enum class MultigridType
  {
    VCycle,
    WCycle,
  };

public:
  /*!
    @brief Constructor
  */
  CdrMultigrid();

  /*!
    @brief Destructor
  */
  virtual ~CdrMultigrid();

  /*!
    @brief Parse class options
  */
  virtual void
  parseOptions() = 0;

  /*!
    @brief Parse class options
  */
  virtual void
  parseRuntimeOptions() = 0;

  /*!
    @brief Register operator
  */
  virtual void
  registerOperators() override;

  /*!
    @brief Allocate internal storage
  */
  virtual void
  allocate() override;

  /*!
    @brief Compute div(J) explicitly, where J = nV - D*grad(n)
    @param[out] a_divJ             Divergence term, i.e. finite volume approximation to 
    @param[in]  a_phi              Cell-centered state
    @param[in]  a_conservativeOnly If true, we compute div(J) = 1/dx*sum(fluxes), which does not involve redistribution.
    @param[in]  a_extrapDt         Extrapolation in time, i.e. shifting of div(J) towards e.g. half time step. Only affects the advective term. 
    @param[in]  a_ebFlux           If true, the embedded boundary flux will injected and included in div(J)
    @note a_phi is non-const because ghost cells will be re-filled
  */
  virtual void
  computeDivJ(EBAMRCellData& a_divJ,
              EBAMRCellData& a_phi,
              const Real     a_extrapDt,
              const bool     a_conservativeOnly,
              const bool     a_ebFlux,
              const bool     a_domainFlux) override final;

  /*!
    @brief Compute div(v*phi) explicitly
    @param[out] a_divF             Divergence term, i.e. finite volume approximation to Div(v*phi), including redistribution magic. 
    @param[in]  a_phi              Cell-centered state
    @param[in]  a_extrapDt         Extrapolation in time, i.e. shifting of div(F) towards e.g. half time step. Only affects the advective term. 
    @param[in]  a_conservativeOnly If true, we compute div(F) = 1/dx*sum(fluxes), which does not involve redistribution.
    @param[in]  a_domainBc         How to set domain fluxes
    @param[in]  a_ebFlux           If true, the embedded boundary flux will injected be included in div(F)
    @note a_phi is non-const because ghost cells will be interpolated in this routine. Valid data in a_phi is not touched. 
  */
  virtual void
  computeDivF(EBAMRCellData& a_divF,
              EBAMRCellData& a_phi,
              const Real     a_extrapDt,
              const bool     a_conservativeOnly,
              const bool     a_ebFlux,
              const bool     a_domainFlux) override final;

  /*!
    @brief Compute div(D*grad(phi)) explicitly
    @param[out] a_divF             Divergence term, including redistribution magic. 
    @param[in]  a_phi              Cell-centered state
    @param[in]  a_conservativeOnly If true, we compute div(D) = 1/dx*sum(fluxes), which does not involve redistribution.
    @param[in]  a_domainBc         How to set domain fluxes
    @param[in]  a_useEbFlux        If true, the embedded boundary flux will be included in div(D)
    @note a_phi is non-const because ghost cells will be interpolated in this routine. Valid data in a_phi is not touched. 
  */
  virtual void
  computeDivD(EBAMRCellData& a_divD,
              EBAMRCellData& a_phi,
              const bool     a_conservativeOnly,
              const bool     a_ebFlux,
              const bool     a_domainFlux) override final;

  /*!
    @brief Implicit diffusion Euler advance with source term. 
    @details This implementation assumes that the incoming previous solution has not been weighted by the volume fraction, and
    that the incoming source term a_source HAS been weighted by the volume fraction. 
    @param[inout] a_newPhi Solution at time t + dt
    @param[in]    a_oldPhi Solution at time t. Should be unweighted by the volume fraction. 
    @param[in]    a_source Source term. Should be weighted by the volume fraction. 
    @param[in]    a_dt     Time step.
  */
  virtual void
  advanceEuler(EBAMRCellData&       a_newPhi,
               const EBAMRCellData& a_oldPhi,
               const EBAMRCellData& a_source,
               const Real           a_dt) override final;

  /*!
    @brief Implicit diffusion Crank-Nicholson advance with source term.
    @details This implementation assumes that the incoming previous solution has not been weighted by the volume fraction, and
    that the incoming source term a_source HAS been weighted by the volume fraction. 
    @param[inout] a_newPhi Solution at time t + dt
    @param[in]    a_oldPhi Solution at time t
    @param[in]    a_source Source term. 
    @param[in]    a_dt     Time step
  */
  virtual void
  advanceCrankNicholson(EBAMRCellData&       a_newPhi,
                        const EBAMRCellData& a_oldPhi,
                        const EBAMRCellData& a_source,
                        const Real           a_dt) override final;

protected:
  /*!
    @brief Relaxation type for gmg
  */
  EBHelmholtzOp::Smoother m_smoother;

  /*!
    @brief GMG multigrid type
  */
  MultigridType m_multigridType;

  /*!
    @brief Geometric multigrid solver
  */
  RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB>>> m_multigridSolver;

  /*!
    @brief Operator factory
  */
  RefCountedPtr<EBHelmholtzOpFactory> m_helmholtzOpFactory;

  /*!
    @brief Conjugate gradient solver bottom MG level
  */
  BiCGStabSolver<LevelData<EBCellFAB>> m_bicgstab;

  /*!
    @brief Simple solver
  */
  EBSimpleSolver m_simpleSolver;

  /*!
    @brief GMRES solver
  */
  GMRESSolver<LevelData<EBCellFAB>> m_gmres;

  /*!
    @brief Data which is always zero. 
  */
  EBAMRCellData m_zero;

  /*!
    @brief Storage for Helmholtz a-coefficient. Always 1. 
  */
  EBAMRCellData m_helmAcoef;

  /*!
    @brief Storage for multigrid residual
  */
  EBAMRCellData m_residual;

  /*!
    @brief Verbosity for geometric multigrid
  */
  int m_multigridVerbosity;

  /*!
    @brief Number of smoothings before averaging
  */
  int m_multigridPreSmooth;

  /*!
    @brief Number of smoothings before averaging
  */
  int m_multigridPostSmooth;

  /*!
    @brief Number of smoothing before bottom solver
  */
  int m_multigridBottomSmooth;

  /*!
    @brief Maximum number of iterations
  */
  int m_multigridMaxIterations;

  /*!
    @brief Minimum number of iterations
  */
  int m_multigridMinIterations;

  /*!
    @brief Bottom solver type
  */
  BottomSolverType m_bottomSolverType;

  /*!
    @brief Number of smoothing for bottom solver
  */
  int m_numSmoothingsForSimpleSolver;

  /*!
    @brief Set bottom drop depth
  */
  int m_minCellsBottom;

  /*!
    @brief Multigrid tolerance
  */
  Real m_multigridExitTolerance;

  /*!
    @brief Multigrid hand
  */
  Real m_multigridExitHang;

  /*!
    @brief Advection-only extrapolation to faces
  */
  virtual void
  advectToFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_phi, const Real a_extrapDt) = 0;

  /*!
    @brief Set up diffusion solver
  */
  virtual void
  setupDiffusionSolver();

  /*!
    @brief Setup the operator factory
  */
  virtual void
  setupHelmholtzFactory();

  /*!
    @brief Setup multigrid
  */
  virtual void
  setupMultigrid();

  /*!
    @brief Reset alpha and beta-coefficients in the multigrid solvers
    @param[in] a_alpha Alpha-coefficient for EBHelmholtzOp
    @param[in] a_beta  Beta-coefficient for EBHelmholtzOp
  */
  virtual void
  resetAlphaAndBeta(const Real a_alpha, const Real a_beta);

  /*!
    @brief Compute kappa * L(phi) using the multigrid operator. This is different from computeDivD.
    @param[in] a_kappaLphi kappa*L(phi)
    @param[in] a_phi        Cell-centered state
    @note This will reset the alpha and beta coefficients. 
  */
  virtual void
  computeKappaLphi(EBAMRCellData& a_kappaLphi, const EBAMRCellData& a_phi);

  /*!
    @brief Parse solver settings for geometric multigrid
  */
  virtual void
  parseMultigridSettings();
};

#include <CD_NamespaceFooter.H>

#endif
