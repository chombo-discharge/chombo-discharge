/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_FieldSolverMG.H
 * @brief  Declaration of FieldSolverMG
 * @author Robert Marskar
 */

#ifndef CD_FIELDSOLVERMG_H
#define CD_FIELDSOLVERMG_H

// Std includes
#include <string>

// Our includes
#include <CD_FieldSolver.H>
#include <CD_FieldSolverGMG.H>
#ifdef CH_USE_PETSC
#include <CD_FieldSolverAMG.H>
#endif
#include <CD_NamespaceHeader.H>

/**
 * @brief Field solver that solves with geometric multigrid, and falls back to algebraic multigrid when that fails.
 * @details This class composes a FieldSolverGMG and, in builds with PETSc, a FieldSolverAMG. Each solve goes through
 * the geometric multigrid solver chain first (e.g. 'gmg bicgstab'). If that does not converge, the solve is repeated
 * with PETSc's algebraic multigrid, starting from the same initial guess. The algebraic multigrid solver is only set up
 * the first time it is needed after a regrid. Once it has rescued a solve, the following solves go straight to it, and
 * geometric multigrid is re-tried every amg_retry_interval solves.
 *
 * The two solvers share this solver's data: they read and write its potential, space and surface charge, electric
 * field, residue and permittivities. They also read their settings from this class's input options, so that the
 * boundary conditions, the plot variables, the interface condition and the EB stencils are given once, and agree. The
 * EB stencil order and weight (bc_order, bc_weight) apply to electrodes and dielectric interfaces in both solvers.
 *
 * Without PETSc, the algebraic multigrid fallback is not compiled, and amg_fallback is ignored with a warning. The
 * fallback is also skipped on periodic domains, and when the EB stencils reach further than the coarse-fine
 * interpolation radius, since algebraic multigrid does not support these.
 *
 * Geometric and algebraic multigrid assemble the operator differently, so a solve that falls back gives a slightly
 * different discrete solution.
 */
class FieldSolverMG : public FieldSolver
{
public:
  /**
   * @brief Constructor. Creates the geometric multigrid solver.
   */
  FieldSolverMG();

  /**
   * @brief Disallowed copy constructor.
   * @param[in] a_other Other solver
   */
  FieldSolverMG(const FieldSolverMG& a_other) = delete;

  /**
   * @brief Disallowed move constructor.
   * @param[in] a_other Other solver
   */
  FieldSolverMG(const FieldSolverMG&& a_other) = delete;

  /**
   * @brief Disallowed copy assignment operator.
   * @param[in] a_other Other solver
   * @return This solver.
   */
  FieldSolverMG&
  operator=(const FieldSolverMG& a_other) = delete;

  /**
   * @brief Disallowed move assignment operator.
   * @param[in] a_other Other solver
   * @return This solver.
   */
  FieldSolverMG&
  operator=(const FieldSolverMG&& a_other) = delete;

  /**
   * @brief Destructor.
   */
  virtual ~FieldSolverMG();

  /**
   * @brief Solve the Poisson equation, with geometric multigrid first and algebraic multigrid as fallback.
   * @param[in,out] a_phi     Potential
   * @param[in]     a_rho     Space charge density
   * @param[in]     a_sigma   Surface charge density on the interface between the gas and the dielectrics
   * @param[in]     a_zeroPhi Use zero as the initial guess or not
   * @return True if one of the solvers converged.
   */
  virtual bool
  solve(MFAMRCellData& a_phi, const MFAMRCellData& a_rho, const EBAMRIVData& a_sigma, const bool a_zeroPhi) override;

  /**
   * @brief Parse class options, including those of the geometric and algebraic multigrid solvers.
   */
  virtual void
  parseOptions() override;

  /**
   * @brief Parse class options that can change at run time.
   */
  virtual void
  parseRuntimeOptions() override;

  /**
   * @brief Register the AMR operators that the geometric and algebraic multigrid solvers need.
   */
  virtual void
  registerOperators() override;

  /**
   * @brief Allocate the data, and share it with the solvers.
   */
  virtual void
  allocate() override;

  /**
   * @brief Release the data, also in the solvers.
   */
  virtual void
  deallocate() override;

  /**
   * @brief Set up geometric multigrid. Algebraic multigrid is set up when it is first needed.
   */
  virtual void
  setupSolver() override;

  /**
   * @brief Discard the discretizations and the solver hierarchies of both solvers.
   */
  virtual void
  invalidateSolver() override;

  /**
   * @brief Change the permittivities used by the solvers, e.g. for a semi-implicit scheme.
   * @details The solver must be set up first. The coefficients stay in use until the solver is set up again, e.g. after
   * a regrid. They are also given to algebraic multigrid if it is set up later.
   * @param[in] a_permittivityCell Cell-centered permittivity
   * @param[in] a_permittivityFace Face-centered permittivity
   * @param[in] a_permittivityEB   EB-centered permittivity
   */
  virtual void
  setSolverPermittivities(const MFAMRCellData& a_permittivityCell,
                          const MFAMRFluxData& a_permittivityFace,
                          const MFAMRIVData&   a_permittivityEB) override;

  /**
   * @brief Cache the potential before a regrid, and discard both solvers.
   * @param[in] a_lbase          Coarsest level which changes during regrid.
   * @param[in] a_oldFinestLevel Finest AMR level before regrid.
   */
  virtual void
  preRegrid(const int a_lbase, const int a_oldFinestLevel) override;

  /**
   * @brief Regrid. The solvers are set up again on the next solve.
   * @param[in] a_lmin           Coarsest level that changed
   * @param[in] a_oldFinestLevel Finest level before the regrid
   * @param[in] a_newFinestLevel Finest level after the regrid
   */
  virtual void
  regrid(const int a_lmin, const int a_oldFinestLevel, const int a_newFinestLevel) override;

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
   * @brief Compute grid loads, from the geometric multigrid relaxation.
   * @param[in] a_dbl   Grids on input level
   * @param[in] a_level Grid level
   * @return Estimated loads, indexed as DisjointBoxLayout::boxArray()
   */
  virtual Vector<long long>
  computeLoads(const DisjointBoxLayout& a_dbl, const int a_level) override;

protected:
  /**
   * @brief The geometric multigrid solver.
   */
  RefCountedPtr<FieldSolverGMG> m_gmg;

  /**
   * @brief Is the solver set up or not.
   */
  bool m_isSolverSetup;

  /**
   * @brief Order of the EB stencils on electrodes and dielectric interfaces, in both solvers.
   */
  int m_bcOrder;

  /**
   * @brief Least-squares weight of the EB stencils on electrodes and dielectric interfaces, in both solvers.
   */
  int m_bcWeight;

  /**
   * @brief True if setSolverPermittivities was called after the last setup.
   */
  bool m_hasSolverPermittivities;

  /**
   * @brief Cell-centered permittivity from the last setSolverPermittivities. Shares the caller's storage.
   */
  MFAMRCellData m_solverPermCell;

  /**
   * @brief Face-centered permittivity from the last setSolverPermittivities. Shares the caller's storage.
   */
  MFAMRFluxData m_solverPermFace;

  /**
   * @brief EB-centered permittivity from the last setSolverPermittivities. Shares the caller's storage.
   */
  MFAMRIVData m_solverPermEB;

#ifdef CH_USE_PETSC
  /**
   * @brief The algebraic multigrid solver. Null unless amg_fallback is true.
   */
  RefCountedPtr<FieldSolverAMG> m_amg;

  /**
   * @brief Initial guess of the current solve. Algebraic multigrid starts from it when geometric multigrid fails.
   */
  MFAMRCellData m_phiStart;

  /**
   * @brief Fall back to algebraic multigrid or not.
   */
  bool m_amgFallback;

  /**
   * @brief True after algebraic multigrid rescued a solve, until geometric multigrid converges again.
   */
  bool m_amgLatched;

  /**
   * @brief While latched, geometric multigrid is re-tried every this many solves.
   */
  int m_amgRetryInterval;

  /**
   * @brief Number of solves since geometric multigrid was last tried.
   */
  int m_solvesSinceProbe;
#endif

  /**
   * @brief Parse the settings of both solvers, the EB stencils and the fallback.
   * @details The geometric multigrid settings have the prefix 'gmg_'. The algebraic multigrid settings have the
   * names that FieldSolverAMG uses.
   */
  virtual void
  parseSolverSettings();

  /**
   * @brief Share this solver's data and configuration with the solvers that exist.
   */
  virtual void
  syncSolvers();

#ifdef CH_USE_PETSC
  /**
   * @brief Check if algebraic multigrid supports the current domain and grids.
   * @param[out] a_reason Why it is not supported, if it is not
   * @return True if algebraic multigrid can be used
   */
  virtual bool
  isAMGSupported(std::string& a_reason) const;

  /**
   * @brief Solve with algebraic multigrid, setting it up first if needed.
   * @param[in,out] a_phi   Potential, which holds the initial guess on input
   * @param[in]     a_rho   Space charge density
   * @param[in]     a_sigma Surface charge density
   * @return True if algebraic multigrid converged.
   */
  virtual bool
  solveAMG(MFAMRCellData& a_phi, const MFAMRCellData& a_rho, const EBAMRIVData& a_sigma);
#endif
};

#include <CD_NamespaceFooter.H>

#endif
