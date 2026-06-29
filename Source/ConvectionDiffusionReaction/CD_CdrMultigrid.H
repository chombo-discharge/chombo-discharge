/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_CdrMultigrid.H
  @brief  Extensions of CdrSolver which use EBHelmholtzOp for diffusion solves
  @author Robert Marskar
*/

#ifndef CD_CDRMULTIGRID_H
#define CD_CDRMULTIGRID_H

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
#include <CD_KrylovMultigrid.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Extension class of CdrSolver that uses multigrid for diffusion part. Can also solve for stochastic diffusion.
*/
class CdrMultigrid : public CdrSolver
{
protected:
  /**
    @brief Enum class for supported bottom solvers in multigrid.
  */
  enum class BottomSolverType
  {
    Simple,
    BiCGStab,
    GMRES
  };

  /**
    @brief Enum for multigrid cycle types. 
  */
  enum class MultigridType
  {
    VCycle,
    WCycle,
  };

public:
  /**
    @brief Constructor
  */
  CdrMultigrid();

  /**
    @brief Destructor
  */
  virtual ~CdrMultigrid();

  /**
    @brief Parse class options
  */
  virtual void
  parseOptions() override = 0;

  /**
    @brief Parse class options
  */
  virtual void
  parseRuntimeOptions() override = 0;

  /**
    @brief Register operator
  */
  virtual void
  registerOperators() override;

  /**
    @brief Allocate internal storage
  */
  virtual void
  allocate() override;

  /**
    @brief Perform pre-regrid operations. 
    @param[in] a_lbase          Coarsest level that changed during regrid. 
    @param[in] a_oldFinestLevel Finest grid level before the regrid operation. 
    @note This stores m_phi and m_source. 
  */
  virtual void
  preRegrid(const int a_lbase, const int a_oldFinestLevel) override;

  /**
    @brief Compute div(J) explicitly, where J = nV - D*grad(n).
    @param[out] a_divJ             Divergence term, i.e. finite volume approximation to div(J).
    @param[in]  a_phi              Cell-centered state.
    @param[in]  a_extrapDt         Extrapolation in time, i.e. shifting of div(J) towards e.g. half time step. Only affects the advective term.
    @param[in]  a_conservativeOnly If true, we compute div(J) = 1/dx*sum(fluxes), which does not involve redistribution.
    @param[in]  a_ebFlux           If true, the embedded boundary flux will injected and included in div(J).
    @param[in]  a_domainFlux       If true, domain boundary fluxes are included.
    @note a_phi is non-const because ghost cells will be re-filled.
  */
  virtual void
  computeDivJ(EBAMRCellData& a_divJ,
              EBAMRCellData& a_phi,
              const Real     a_extrapDt,
              const bool     a_conservativeOnly,
              const bool     a_ebFlux,
              const bool     a_domainFlux) override final;

  /**
    @brief Compute div(v*phi) explicitly.
    @param[out] a_divF             Divergence term, i.e. finite volume approximation to Div(v*phi), including redistribution.
    @param[in]  a_phi              Cell-centered state.
    @param[in]  a_extrapDt         Extrapolation in time, i.e. shifting of div(F) towards e.g. half time step. Only affects the advective term.
    @param[in]  a_conservativeOnly If true, we compute div(F) = 1/dx*sum(fluxes), which does not involve redistribution.
    @param[in]  a_ebFlux           If true, the embedded boundary flux will be included in div(F).
    @param[in]  a_domainFlux       If true, domain boundary fluxes are included.
    @note a_phi is non-const because ghost cells will be interpolated in this routine. Valid data in a_phi is not touched.
  */
  virtual void
  computeDivF(EBAMRCellData& a_divF,
              EBAMRCellData& a_phi,
              const Real     a_extrapDt,
              const bool     a_conservativeOnly,
              const bool     a_ebFlux,
              const bool     a_domainFlux) override final;

  /**
    @brief Compute div(D*grad(phi)) explicitly.
    @param[out] a_divD             Divergence term, including redistribution.
    @param[in]  a_phi              Cell-centered state.
    @param[in]  a_conservativeOnly If true, we compute div(D) = 1/dx*sum(fluxes), which does not involve redistribution.
    @param[in]  a_ebFlux           If true, the embedded boundary flux will be included in div(D).
    @param[in]  a_domainFlux       If true, domain boundary fluxes are included.
    @note a_phi is non-const because ghost cells will be interpolated in this routine. Valid data in a_phi is not touched.
  */
  virtual void
  computeDivD(EBAMRCellData& a_divD,
              EBAMRCellData& a_phi,
              const bool     a_conservativeOnly,
              const bool     a_ebFlux,
              const bool     a_domainFlux) override final;

  /**
    @brief Implicit diffusion Euler advance with source term. 
    @details This implementation assumes that the incoming previous solution has not been weighted by the volume fraction, and
    that the incoming source term a_source HAS been weighted by the volume fraction. 
    @param[in,out] a_newPhi Solution at time t + dt
    @param[in]    a_oldPhi Solution at time t. Should be unweighted by the volume fraction. 
    @param[in]    a_source Source term. Should be weighted by the volume fraction. 
    @param[in]    a_dt     Time step.
  */
  virtual void
  advanceEuler(EBAMRCellData&       a_newPhi,
               const EBAMRCellData& a_oldPhi,
               const EBAMRCellData& a_source,
               const Real           a_dt) override final;

  /**
    @brief Implicit diffusion Crank-Nicholson advance with source term.
    @details This implementation assumes that the incoming previous solution has not been weighted by the volume fraction, and
    that the incoming source term a_source HAS been weighted by the volume fraction. 
    @param[in,out] a_newPhi Solution at time t + dt
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
  /**
    @brief Relaxation type for gmg
  */
  EBHelmholtzOp::Smoother m_smoother;

  /**
    @brief Chebyshev polynomial degree (only used when m_smoother is Smoother::Chebyshev)
  */
  int m_chebyOrder;

  /**
    @brief Chebyshev eigenvalue ratio lambda_max/lambda_min (only used when m_smoother is Smoother::Chebyshev)
  */
  Real m_chebyEigRatio;

  /**
    @brief Number of inner frozen-ghost sweeps per restricted-additive-Schwarz iteration (only used when m_smoother is Smoother::RestrictedAdditiveSchwarz)
  */
  int m_rasInnerSweeps;

  /**
    @brief GMG multigrid type
  */
  MultigridType m_multigridType;

  /**
    @brief Geometric multigrid solver
  */
  RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB>>> m_multigridSolver;

  /**
    @brief Operator factory
  */
  RefCountedPtr<EBHelmholtzOpFactory> m_helmholtzOpFactory;

  /**
    @brief Conjugate gradient solver bottom MG level
  */
  BiCGStabSolver<LevelData<EBCellFAB>> m_bicgstab;

  /**
    @brief Simple solver
  */
  EBSimpleSolver m_simpleSolver;

  /**
    @brief GMRES solver
  */
  GMRESSolver<LevelData<EBCellFAB>> m_gmres;

  /**
    @brief Outer-Krylov settings (selects gmg/gmres/bicgstab and the krylov_* parameters)
  */
  KrylovMultigrid::Settings m_krylovSettings;

  /**
    @brief Fallback-chain latch policy (persists across regrids; re-probes GMG every krylov_retry_interval solves)
  */
  KrylovMultigrid::FallbackPolicy m_fallbackPolicy;

  /**
    @brief Outer-Krylov driver (V-cycle-preconditioned BiCGStab/GMRES with persistent work vectors)
  */
  KrylovMultigrid::Solver<EBCellFAB> m_krylov;

  /**
    @brief Is multigrid solver defined or not?
  */
  bool m_hasMultigridSolver;

  /**
    @brief Use reflux-free AMR operator in GMG
  */
  bool m_multigridRefluxFree;

  /**
    @brief Storage for Helmholtz a-coefficient. Always 1. 
  */
  EBAMRCellData m_helmAcoef;

  /**
    @brief Storage for multigrid residual
  */
  EBAMRCellData m_residual;

  /**
    @brief Verbosity for geometric multigrid
  */
  int m_multigridVerbosity;

  /**
    @brief Number of smoothings before averaging
  */
  int m_multigridPreSmooth;

  /**
    @brief Number of smoothings before averaging
  */
  int m_multigridPostSmooth;

  /**
    @brief Number of smoothing before bottom solver
  */
  int m_multigridBottomSmooth;

  /**
    @brief Maximum number of iterations
  */
  int m_multigridMaxIterations;

  /**
    @brief Minimum number of iterations
  */
  int m_multigridMinIterations;

  /**
    @brief Bottom solver type
  */
  BottomSolverType m_bottomSolverType;

  /**
    @brief Number of smoothing for bottom solver
  */
  int m_numSmoothingsForSimpleSolver;

  /**
    @brief Set bottom drop depth
  */
  int m_minCellsBottom;

  /**
    @brief Multigrid tolerance
  */
  Real m_multigridExitTolerance;

  /**
    @brief Multigrid hand
  */
  Real m_multigridExitHang;

  /**
    @brief Advection-only extrapolation to faces
    @param[out] a_facePhi Face-centered states
    @param[in] a_phi Cell-centered states
    @param[in] a_extrapDt Extrapolating time step.
  */
  virtual void
  advectToFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_phi, const Real a_extrapDt) override = 0;

  /**
    @brief Set up diffusion solver
  */
  virtual void
  setupDiffusionSolver();

  /**
    @brief Setup the operator factory
  */
  virtual void
  setupHelmholtzFactory();

  /**
    @brief Setup multigrid
  */
  virtual void
  setupMultigrid();

  /**
    @brief Set the multigrid solver coefficients. 
    @details Useful when coefficients change underneath us, but we don't want to set up multigrid again.
  */
  virtual void
  setMultigridSolverCoefficients();

  /**
    @brief Reset alpha and beta-coefficients in the multigrid solvers
    @param[in] a_alpha Alpha-coefficient for EBHelmholtzOp
    @param[in] a_beta  Beta-coefficient for EBHelmholtzOp
  */
  virtual void
  resetAlphaAndBeta(const Real a_alpha, const Real a_beta);

  /**
    @brief Compute kappa * L(phi) using the multigrid operator. This is different from computeDivD.
    @param[in] a_kappaLphi kappa*L(phi)
    @param[in] a_phi        Cell-centered state
    @note This will reset the alpha and beta coefficients. 
  */
  virtual void
  computeKappaLphi(EBAMRCellData& a_kappaLphi, const EBAMRCellData& a_phi);

  /**
    @brief Parse solver settings for geometric multigrid
  */
  virtual void
  parseMultigridSettings();
};

#include <CD_NamespaceFooter.H>

#endif
