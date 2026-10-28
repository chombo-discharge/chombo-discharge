/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_FieldSolverAMG.H
 * @brief  Declaration of the FieldSolverAMG class
 * @author Robert Marskar
 */

#ifdef CH_USE_PETSC

#ifndef CD_FIELDSOLVERAMG_H
#define CD_FIELDSOLVERAMG_H

// Std includes
#include <string>

// Petsc includes
#include <petsc.h>

// Our includes
#include <CD_FieldSolver.H>
#include <CD_PetscHelmholtzOp.H>
#include <CD_PetscHelmholtzEBBC.H>
#include <CD_PetscHelmholtzJumpBC.H>
#include <CD_NamespaceHeader.H>

/**
 * @brief Electrostatic field solver that assembles the operator as a PETSc matrix and solves it with a PETSc Krylov
 * method, preconditioned by algebraic multigrid.
 * @details This is the PETSc counterpart of FieldSolverGMG, with the same inputs and outputs. The operator is
 * -div(eps * grad(phi)) = rho/eps0, assembled by PetscHelmholtzOp. The operator spans all AMR levels and both phases,
 * with electrodes on the EB (PetscHelmholtzElectrostaticEBBC) and the natural or saturation-charge condition on the
 * interface between the gas and the dielectrics. Periodic domains and multi-valued cells are not supported.
 */
class FieldSolverAMG : public FieldSolver
{
public:
  /**
   * @brief Default constructor.
   */
  FieldSolverAMG();

  /**
   * @brief Disallowed copy constructor.
   * @param[in] a_other Other solver
   */
  FieldSolverAMG(const FieldSolverAMG& a_other) = delete;

  /**
   * @brief Disallowed move constructor.
   * @param[in] a_other Other solver
   */
  FieldSolverAMG(const FieldSolverAMG&& a_other) = delete;

  /**
   * @brief Disallowed copy assignment operator.
   * @param[in] a_other Other solver
   * @return This solver.
   */
  FieldSolverAMG&
  operator=(const FieldSolverAMG& a_other) = delete;

  /**
   * @brief Disallowed move assignment operator.
   * @param[in] a_other Other solver
   * @return This solver.
   */
  FieldSolverAMG&
  operator=(const FieldSolverAMG&& a_other) = delete;

  /**
   * @brief Destructor.
   */
  virtual ~FieldSolverAMG();

  /**
   * @brief Solve the Poisson equation.
   * @param[in,out] a_potential Potential
   * @param[in]     a_rho       Space charge density
   * @param[in]     a_sigma     Surface charge density on the interface between the gas and the dielectrics
   * @param[in]     a_zeroPhi   Use zero as the initial guess or not
   * @return True if the solver converged.
   */
  virtual bool
  solve(MFAMRCellData&       a_potential,
        const MFAMRCellData& a_rho,
        const EBAMRIVData&   a_sigma,
        const bool           a_zeroPhi = false) override;

  /**
   * @brief Parse class options.
   */
  virtual void
  parseOptions() override;

  /**
   * @brief Parse class options that can change at run time.
   */
  virtual void
  parseRuntimeOptions() override;

  /**
   * @brief Compute the cell-centered electric field.
   * @param[out] a_E         Electric field
   * @param[in]  a_potential Potential
   */
  virtual void
  computeElectricField(MFAMRCellData& a_E, const MFAMRCellData& a_potential) const override;

  /**
   * @brief Compute the face-centered electric field.
   * @param[out] a_E         Electric field
   * @param[in]  a_potential Potential
   */
  virtual void
  computeElectricField(MFAMRFluxData& a_E, const MFAMRCellData& a_potential) const override;

  /**
   * @brief Compute the cell-centered electric field on a phase.
   * @param[out] a_E         Electric field
   * @param[in]  a_phase     Phase
   * @param[in]  a_potential Potential
   */
  virtual void
  computeElectricField(EBAMRCellData&           a_E,
                       const phase::which_phase a_phase,
                       const MFAMRCellData&     a_potential) const override;

  /**
   * @brief Compute the face-centered electric field on a phase.
   * @param[out] a_E         Electric field
   * @param[in]  a_phase     Phase
   * @param[in]  a_potential Potential
   */
  virtual void
  computeElectricField(EBAMRFluxData&           a_E,
                       const phase::which_phase a_phase,
                       const MFAMRCellData&     a_potential) const override;

  /**
   * @brief Regrid. The operator is rebuilt on the next solve.
   * @param[in] a_lmin           Coarsest level that changed
   * @param[in] a_oldFinestLevel Finest level before the regrid
   * @param[in] a_newFinestLevel Finest level after the regrid
   */
  virtual void
  regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

  /**
   * @brief Register the AMR operators this class needs.
   */
  virtual void
  registerOperators() override;

  /**
   * @brief Assemble the operator and set up the PETSc solver.
   */
  virtual void
  setupSolver() override;

  /**
   * @brief Change the permittivities used by the solver, e.g. for a semi-implicit scheme.
   * @details The solver must be set up first. The operator is refilled with the new coefficients, and the
   * preconditioner is rebuilt at the next solve. The coefficients stay in use until the solver is set up again, e.g.
   * after a regrid.
   * @param[in] a_permittivityCell Cell-centered permittivity
   * @param[in] a_permittivityFace Face-centered permittivity
   * @param[in] a_permittivityEB   EB-centered permittivity
   */
  virtual void
  setSolverPermittivities(const MFAMRCellData& a_permittivityCell,
                          const MFAMRFluxData& a_permittivityFace,
                          const MFAMRIVData&   a_permittivityEB) override;

protected:
  /**
   * @brief Operator alpha. There is no diagonal term.
   */
  static constexpr Real m_alpha = 0.0;

  /**
   * @brief Operator beta. The operator is -div(eps * grad(phi)).
   */
  static constexpr Real m_beta = -1.0;

  /**
   * @brief Is the solver set up or not.
   */
  bool m_isSolverSetup;

  /**
   * @brief PETSc Krylov solver type, e.g. 'gmres', 'fgmres', 'bcgs' (also given as 'bicgstab') or 'cg'.
   */
  std::string m_kspType;

  /**
   * @brief Restart length of GMRES and FGMRES.
   */
  int m_kspGmresRestart;

  /**
   * @brief KSP type of the smoother on each level of a multigrid preconditioner, e.g. 'chebyshev' or 'richardson'.
   */
  std::string m_amgSmootherKsp;

  /**
   * @brief PC type of the smoother on each level of a multigrid preconditioner, e.g. 'sor' or 'jacobi'.
   */
  std::string m_amgSmootherPc;

  /**
   * @brief Number of smoothing steps on each level of a multigrid preconditioner.
   */
  int m_amgSmootherIts;

  /**
   * @brief GAMG graph threshold: connections weaker than this are dropped when aggregating. Negative keeps them all.
   */
  Real m_gamgThreshold;

  /**
   * @brief Number of levels that GAMG coarsens aggressively (square graph).
   */
  int m_gamgAggressiveCoarsening;

  /**
   * @brief PETSc preconditioner type, e.g. 'gamg' or 'hypre'.
   */
  std::string m_pcType;

  /**
   * @brief Relative residual tolerance.
   */
  Real m_kspRtol;

  /**
   * @brief Absolute residual tolerance.
   */
  Real m_kspAtol;

  /**
   * @brief Maximum number of Krylov iterations.
   */
  int m_kspMaxIter;

  /**
   * @brief Conditions on the interface between the gas and the dielectrics.
   */
  enum class JumpBCType
  {
    Natural,
    SaturationCharge
  };

  /**
   * @brief Interface condition.
   */
  JumpBCType m_jumpBcType;

  /**
   * @brief Impose the interface conditions at the points on the dielectric surface closest to the EB face centroids,
   * rather than at the EB face centroids.
   */
  bool m_interfacePoints;

  /**
   * @brief Highest order of the EB face stencils, for the electrode and the interface faces.
   */
  int m_ebOrder;

  /**
   * @brief Least-squares weighting of the EB face stencils.
   */
  int m_ebWeight;

  /**
   * @brief Weight the space charge density by the volume fraction in solve, or not, as in FieldSolverGMG.
   */
  bool m_kappaSource;

  /**
   * @brief Boundary conditions on the electrode faces.
   */
  RefCountedPtr<PetscHelmholtzEBBC> m_petscEbBc;

  /**
   * @brief Conditions on the interface between the gas and the dielectrics.
   */
  RefCountedPtr<PetscHelmholtzJumpBC> m_petscJumpBc;

  /**
   * @brief The operator.
   */
  RefCountedPtr<PetscHelmholtzOp> m_petscOp;

  /**
   * @brief The PETSc solver.
   */
  KSP m_ksp;

  /**
   * @brief Verbosity of the PETSc solve.
   * @details At 1 and above, the residual is printed at every Krylov iteration. At 2 and above, the true residual is
   * printed too. At 3 and above, the solver and the preconditioner are viewed after each solve.
   */
  int m_petscVerbosity;

  /**
   * @brief The 2-norm of the zero-residual b - A*0 = b of the current solve, which the true residual is printed
   * relative to.
   */
  Real m_monitorZeroResidNorm;

  /**
   * @brief Residual norm of the KSP at its first iteration of the current solve.
   */
  Real m_monitorFirstNorm;

  /**
   * @brief Residual norm of the KSP at the previous iteration of the current solve.
   */
  Real m_monitorPreviousNorm;

  /**
   * @brief Parse the PETSc solver settings.
   */
  virtual void
  parseSolverSettings();

  /**
   * @brief KSP monitor that prints the residual at every iteration, as the Krylov solvers of FieldSolverGMG do.
   * @details The residual is in the norm the KSP tests convergence on (the preconditioned residual for
   * left-preconditioned KSPs), relative to that norm at the first iteration. The rate is the reduction factor since
   * the previous iteration. With petsc_verbosity 2 and above the true residual b - A*x is also printed, relative to the
   * zero-residual b.
   * @param[in] a_ksp          The KSP
   * @param[in] a_iteration    Iteration number
   * @param[in] a_residualNorm Residual norm of the KSP
   * @param[in] a_context      The FieldSolverAMG instance
   * @return PETSc error code
   */
  static PetscErrorCode
  monitorResidual(KSP a_ksp, PetscInt a_iteration, PetscReal a_residualNorm, void* a_context);

  /**
   * @brief Destroy the PETSc solver and the operator.
   */
  virtual void
  destroySolver();
};

#include <CD_NamespaceFooter.H>

#endif

#endif
