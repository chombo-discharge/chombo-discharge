/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_KrylovMultigrid.H
  @brief  Outer Krylov solve (GMRES/BiCGStab) preconditioned by an AMR-multigrid V-cycle
  @author Robert Marskar
*/

#ifndef CD_KrylovMultigrid_H
#define CD_KrylovMultigrid_H

// Chombo includes
#include <AMRMultiGrid.H>
#include <LevelData.H>
#include <Vector.H>
#include <ParmParse.H>

// Our includes
#include <CD_AMRMultigridKrylovOp.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Support code for using an AMR-multigrid V-cycle as a preconditioner for an outer Krylov solver.
  @details The V-cycle itself (smoother, cycle type, bottom solver, BC stencils, ...) is configured by the usual
           gmg_* options of the host solver. The settings here control the outer Krylov wrapper.
*/
namespace KrylovMultigrid {

  /*!
    @brief Top-level solve method.
  */
  enum class SolverType
  {
    GMG,     /*!< Stand-alone multigrid V-cycling (no outer Krylov). */
    GMRES,   /*!< GMRES preconditioned by the V-cycle. */
    BiCGStab /*!< BiCGStab preconditioned by the V-cycle. */
  };

  /*!
    @brief Outer-Krylov configuration. Defaults reproduce stand-alone multigrid (solvers == {GMG}).
  */
  struct Settings
  {
    /*!
      @brief Solvers attempted in order; the first that converges wins. E.g. {GMG, GMRES} means "try stand-alone
             multigrid, fall back to GMRES if it does not converge".
    */
    Vector<SolverType> solvers = Vector<SolverType>(1, SolverType::GMG);

    Real eps        = 1.E-10; /*!< Relative residual tolerance. */
    int  maxIter    = 50;     /*!< Maximum outer iterations. */
    int  restart    = 30;     /*!< GMRES restart length (unused for BiCGStab). */
    int  numVCycles = 1;      /*!< V-cycles applied per preconditioner invocation. */
    int  verbosity  = -1;     /*!< Outer Krylov verbosity. */

    /*!
      @brief Whether any solver in the chain is a Krylov method (and hence needs the preconditioner adapter).
      @return True if the chain contains GMRES or BiCGStab.
    */
    bool
    usesKrylov() const
    {
      for (int i = 0; i < solvers.size(); i++) {
        if (solvers[i] != SolverType::GMG) {
          return true;
        }
      }
      return false;
    }
  };

  /*!
    @brief Parse the solver path and Krylov settings from an input group.
    @details Reads the 'solver' key (gmg|gmres|bicgstab) and the 'krylov_*' keys from the given input group. All are
             optional; absent keys keep the defaults so existing input scripts default to stand-alone multigrid.
    @param[in]  a_pp       ParmParse scoped to the host solver's class name
    @param[out] a_settings Settings to fill
  */
  inline void
  parseSettings(ParmParse& a_pp, Settings& a_settings)
  {
    const int num = a_pp.contains("solver") ? a_pp.countval("solver") : 0;
    if (num > 0) {
      a_settings.solvers.resize(num);
      for (int i = 0; i < num; i++) {
        std::string str;
        a_pp.get("solver", str, i);
        if (str == "gmg") {
          a_settings.solvers[i] = SolverType::GMG;
        }
        else if (str == "gmres") {
          a_settings.solvers[i] = SolverType::GMRES;
        }
        else if (str == "bicgstab") {
          a_settings.solvers[i] = SolverType::BiCGStab;
        }
        else {
          MayDay::Error("KrylovMultigrid::parseSettings - 'solver' entries must be 'gmg', 'gmres', or 'bicgstab'");
        }
      }
    }

    a_pp.query("krylov_eps", a_settings.eps);
    a_pp.query("krylov_max_iter", a_settings.maxIter);
    a_pp.query("krylov_restart", a_settings.restart);
    a_pp.query("krylov_vcycles", a_settings.numVCycles);
    a_pp.query("krylov_verbosity", a_settings.verbosity);
  }

  /*!
    @brief Solve L(phi) = rhs with an outer Krylov method preconditioned by the multigrid V-cycle.
    @details Runs in residual-correction form so the inhomogeneous boundary contribution (and any jump/surface-charge
             contribution) enters only through the initial residual r0 = rhs - L(phi); the matrix-vector product and
             preconditioner are applied with homogeneous boundary conditions. Warm starts are honoured: when
             a_zeroPhi is false the incoming a_phi is used as the initial guess.
    @param[in]     a_mg       Fully set-up multigrid solver (reused for matvec and preconditioning)
    @param[in,out] a_op       Adapter, already define()d over a_mg for this grid hierarchy
    @param[in,out] a_phi      Solution (in: initial guess unless a_zeroPhi; out: solution)
    @param[in]     a_rhs      Right-hand side
    @param[in]     a_lbase    Coarsest level
    @param[in]     a_lmax     Finest level
    @param[in]     a_zeroPhi  If true, start from phi = 0
    @param[in]     a_type     Which Krylov method to run (GMRES or BiCGStab)
    @param[in]     a_settings Outer-Krylov settings (tolerance, iterations, restart, ...)
    @return True if the Krylov solver reported convergence.
  */
  template <class T>
  bool
  solve(AMRMultiGrid<LevelData<T>>*  a_mg,
        AMRMultigridKrylovOp<T>&     a_op,
        Vector<LevelData<T>*>&       a_phi,
        const Vector<LevelData<T>*>& a_rhs,
        const int                    a_lbase,
        const int                    a_lmax,
        const bool                   a_zeroPhi,
        const SolverType             a_type,
        const Settings&              a_settings);

} // namespace KrylovMultigrid

#include <CD_NamespaceFooter.H>

#include <CD_KrylovMultigridImplem.H>

#endif
