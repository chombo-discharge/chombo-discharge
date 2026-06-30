/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_KrylovMultigrid.H
  @brief  Outer Krylov solve (GMRES/BiCGStab) preconditioned by an AMR-multigrid cycle
  @author Robert Marskar
*/

#ifndef CD_KrylovMultigrid_H
#define CD_KrylovMultigrid_H

// Std includes
#include <algorithm>
#include <string>

// Chombo includes
#include <AMRMultiGrid.H>
#include <LevelData.H>
#include <Vector.H>
#include <parstream.H>

// Our includes
#include <CD_AMRMultigridKrylovOp.H>
#include <CD_KrylovBiCGStab.H>
#include <CD_KrylovGMRES.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Support code for using an AMR-multigrid cycle as a preconditioner for an outer Krylov solver.
  @details The multigrid cycle itself (smoother, cycle type, bottom solver, BC stencils, ...) is configured by the usual
           gmg_* options of the host solver. The settings here control the outer Krylov wrapper.
*/
namespace KrylovMultigrid {

  /**
    @brief Top-level solve method.
  */
  enum class SolverType
  {
    /**
      @brief Stand-alone multigrid V-cycling (no outer Krylov).
    */
    GMG,

    /**
      @brief GMRES preconditioned by the multigrid cycle.
    */
    GMRES,

    /**
      @brief BiCGStab preconditioned by the multigrid cycle.
    */
    BiCGStab
  };

  /**
    @brief Outer-Krylov configuration. Defaults reproduce stand-alone multigrid (solvers == {GMG}).
  */
  struct Settings
  {
    /**
      @brief Solvers attempted in order; the first that converges wins. E.g. {GMG, GMRES} means "try stand-alone
             multigrid, fall back to GMRES if it does not converge".
    */
    Vector<SolverType> solvers = Vector<SolverType>(1, SolverType::GMG);

    /**
      @brief Relative residual tolerance.
    */
    Real eps = 1.E-10;

    /**
      @brief Maximum outer iterations.
    */
    int maxIter = 50;

    /**
      @brief GMRES restart length (unused for BiCGStab).
    */
    int restart = 30;

    /**
      @brief Number of multigrid cycles applied per preconditioner invocation. The cycle type (V/W) is deferred
             to the multigrid configuration via gmg_cycle; it should normally be a V-cycle for efficiency.
    */
    int numVCycles = 1;

    /**
      @brief Outer Krylov verbosity.
    */
    int verbosity = -1;

    /**
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

  /**
    @brief Stateful policy for the ``solver`` fallback chain.
    @details When the leading stand-alone multigrid (GMG) entry of a chain such as ``gmg gmres`` fails, the chain
             falls back to the Krylov solver. Without memory the next solve pays a full failed GMG attempt again, which
             is wasteful when the failure is persistent (e.g. a fixed, hard field operator). This policy latches onto
             the Krylov fallback after a GMG failure and only re-probes GMG every ``retryInterval`` solves, reclaiming
             multigrid's speed whenever the problem becomes tractable again. ``retryInterval == 1`` reproduces the
             memoryless "always retry GMG first" behaviour; a large value approaches a permanent latch. Optionally,
             ``reprobeIters`` arms an *adaptive* early re-probe: if the Krylov fallback converges in very few iterations
             the preconditioner is clearly strong, so GMG alone is likely to succeed and is re-probed on the next solve
             without waiting for ``retryInterval``.
    @note The state is deliberately decoupled from regridding: hold an instance as a host-solver member and do *not*
          reset it when the operators are rebuilt (regrids may be far more frequent than the desired re-probe cadence).
          The counter advances per solve, so the re-probe cadence is measured in solves, not necessarily time steps.
  */
  struct FallbackPolicy
  {
    /**
      @brief Re-probe GMG every this many solves while latched. Clamped to >= 1; 1 reproduces always-retry behaviour.
    */
    int retryInterval = 1;

    /**
      @brief Adaptive early re-probe: re-probe GMG on the next solve whenever the Krylov fallback converged in this many
             iterations or fewer (a strong preconditioner is a sign that GMG alone may now succeed). 0 disables it, so
             only the periodic retryInterval applies.
    */
    int reprobeIters = 0;

    /**
      @brief Start index into a_solvers for this solve, advancing the re-probe counter. Call exactly once per solve.
      @details Returns 0 (include the leading GMG) when unlatched, when the periodic re-probe is due, or when an
               adaptive early re-probe was armed by a cheap Krylov solve; otherwise returns the index of the first
               non-GMG solver so the persistently-failing GMG attempt is skipped.
      @param[in] a_solvers The configured fallback chain.
      @return Start index into a_solvers.
    */
    int
    startIndex(const Vector<SolverType>& a_solvers)
    {
      if (!m_latched) {
        return 0;
      }

      m_sinceProbe++;
      if (m_reprobeNext || m_sinceProbe >= std::max(1, retryInterval)) {
        m_sinceProbe  = 0;
        m_reprobeNext = false;

        return 0; // re-probe GMG (periodic or adaptive)
      }

      return this->firstNonGmg(a_solvers);
    }

    /**
      @brief Record the outcome of this solve and update the latch.
      @details Latching (a failed GMG attempt) is already announced by the host's fallback-chain message; this method
               additionally announces the *release* (a successful re-probe) at verbosity >= 1, since that transition is
               otherwise silent. While latched and skipping GMG, a Krylov solve that converged in <= reprobeIters
               iterations arms an adaptive early re-probe of GMG on the next solve.
      @param[in] a_gmgAttempted Whether a stand-alone GMG attempt was made this solve.
      @param[in] a_gmgConverged Whether that GMG attempt converged.
      @param[in] a_krylovIters  Iterations taken by the Krylov fallback this solve, or < 0 if none ran.
      @param[in] a_verbosity    Host verbosity; a release message is printed when >= 1.
      @param[in] a_name         Host solver name, used to prefix the release message.
    */
    void
    recordOutcome(const bool         a_gmgAttempted,
                  const bool         a_gmgConverged,
                  const int          a_krylovIters = -1,
                  const int          a_verbosity   = -1,
                  const std::string& a_name        = "")
    {
      if (a_gmgAttempted) {
        const bool wasLatched = m_latched;

        m_latched = !a_gmgConverged;
        if (m_latched) {
          m_sinceProbe = 0;
        }

        if (a_verbosity >= 1 && wasLatched && !m_latched) {
          pout() << a_name << " - multigrid re-probe succeeded; releasing the Krylov fallback latch" << endl;
        }
      }
      else if (m_latched && reprobeIters > 0 && a_krylovIters >= 0 && a_krylovIters <= reprobeIters) {
        // The Krylov fallback found the problem easy, so GMG alone may now succeed -- re-probe it next solve.
        m_reprobeNext = true;
      }
    }

    /**
      @brief Whether the policy is currently latched onto the Krylov fallback.
      @return True if latched.
    */
    bool
    isLatched() const
    {
      return m_latched;
    }

  protected:
    /**
      @brief Whether a failed GMG attempt has latched us onto the Krylov fallback.
    */
    bool m_latched = false;

    /**
      @brief Solves since the last GMG attempt while latched.
    */
    int m_sinceProbe = 0;

    /**
      @brief Whether a cheap Krylov solve has armed an adaptive early re-probe of GMG on the next solve.
    */
    bool m_reprobeNext = false;

    /**
      @brief Index of the first non-GMG solver in the chain, or 0 if there is none.
      @param[in] a_solvers The fallback chain.
      @return Start index that skips a leading GMG.
    */
    int
    firstNonGmg(const Vector<SolverType>& a_solvers) const
    {
      for (int i = 0; i < a_solvers.size(); i++) {
        if (a_solvers[i] != SolverType::GMG) {
          return i;
        }
      }

      return 0;
    }
  };

  /**
    @brief Map a solver-type string to the corresponding enum.
    @details Pure string-to-enum conversion; the host solver is responsible for reading the input database (the
             Krylov module deliberately does not call ParmParse). Aborts on an unrecognised string.
    @param[in] a_str One of 'gmg', 'gmres', or 'bicgstab'.
    @return The corresponding SolverType.
  */
  inline SolverType
  toSolverType(const std::string& a_str)
  {
    if (a_str == "gmg") {
      return SolverType::GMG;
    }
    else if (a_str == "gmres") {
      return SolverType::GMRES;
    }
    else if (a_str == "bicgstab") {
      return SolverType::BiCGStab;
    }
    else {
      MayDay::Error("KrylovMultigrid::toSolverType - solver must be 'gmg', 'gmres', or 'bicgstab'");
    }

    return SolverType::GMG;
  }

  /**
    @brief Human-readable name for a solver type, for diagnostic output.
    @param[in] a_type Solver type.
    @return "gmg", "gmres", or "bicgstab".
  */
  inline std::string
  solverTypeName(const SolverType a_type)
  {
    switch (a_type) {
    case SolverType::GMG:
      return "gmg";
    case SolverType::GMRES:
      return "gmres";
    case SolverType::BiCGStab:
      return "bicgstab";
    }

    return "unknown";
  }

  /**
    @brief Stateful owner of the outer-Krylov solver path.
    @details Holds the multigrid-cycle adapter (AMRMultigridKrylovOp), the chombo-discharge BiCGStab/GMRES solvers, and the
             persistent residual-correction work vectors. define() is called once per regrid; solve() reuses the
             adapter, the (lazily defined) custom solvers and the work vectors across all subsequent solves, so there
             is no per-solve work-vector allocation. solve() runs in residual-correction form so the inhomogeneous
             boundary/jump contribution enters only through the initial residual r0 = rhs - L(phi); the matvec and
             preconditioner are applied with homogeneous boundary conditions, and warm starts are honoured.
    @note T is the FAB type (EBCellFAB or MFCellFAB).
  */
  template <class T>
  class Solver
  {
  public:
    /**
      @brief Multilevel data type.
    */
    using Vec = Vector<LevelData<T>*>;

    /**
      @brief Default constructor. Call define() before use.
    */
    Solver() noexcept = default;

    /**
      @brief Destructor. Frees the persistent work vectors and custom solvers.
    */
    ~Solver() noexcept;

    /**
      @brief Define the multigrid-preconditioner adapter and reset the persistent work. Call once per regrid.
      @param[in] a_mg         Fully set-up multigrid solver (reused for matvec and preconditioning)
      @param[in] a_validCells Per-level valid-cell mask (true where not covered by a finer level), e.g. from
                              AmrMesh::getValidCells
      @param[in] a_dx         Grid spacing on each level (index 0..a_lmax)
      @param[in] a_lbase      Coarsest level in the solve
      @param[in] a_lmax       Finest level in the solve
      @param[in] a_numVCycles Number of multigrid cycles per preconditioner invocation (cycle type from gmg_cycle; normally a V-cycle)
    */
    void
    define(AMRMultiGrid<LevelData<T>>*                            a_mg,
           const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_validCells,
           const Vector<Real>&                                    a_dx,
           const int                                              a_lbase,
           const int                                              a_lmax,
           const int                                              a_numVCycles);

    /**
      @brief Solve L(phi) = rhs with the chosen outer Krylov method preconditioned by the multigrid cycle.
      @param[in,out] a_phi      Solution (in: initial guess unless a_zeroPhi; out: solution)
      @param[in]     a_rhs      Right-hand side
      @param[in]     a_zeroPhi  If true, start from phi = 0
      @param[in]     a_type     Which Krylov method to run (GMRES or BiCGStab)
      @param[in]     a_settings Outer-Krylov settings (tolerance, iterations, restart, verbosity)
      @return True if the Krylov solver reported convergence.
    */
    bool
    solve(Vec& a_phi, const Vec& a_rhs, const bool a_zeroPhi, const SolverType a_type, const Settings& a_settings);

    /**
      @brief Iteration count of the most recent Krylov solve (GMRES/BiCGStab).
      @return Number of iterations, or -1 if no Krylov solve has run.
    */
    int
    lastKrylovIterations() const
    {
      return m_lastKrylovIters;
    }

    /**
      @brief Snapshot the chain's starting guess and its residual for the fallback warm-start rule.
      @details Call once, before the first attempt in a fallback chain. Stores a copy of the starting solution (the
               "initial-residual solution") so keepOrRestore() can fall back to it.
      @param[in,out] a_phi The starting guess (after any zeroing the host has applied); its ghost cells may be filled.
      @param[in]     a_rhs Right-hand side, used to evaluate the initial residual.
    */
    void
    snapshotGuess(Vec& a_phi, const Vec& a_rhs);

    /**
      @brief Warm-start rule applied between fallback attempts.
      @details If the just-failed attempt left a residual smaller than the snapshot's initial residual, a_phi is kept
               (the failed solve is a better warm start); otherwise a_phi is restored to the snapshot (the
               initial-residual solution). snapshotGuess() must have been called first.
      @param[in,out] a_phi Current solution; restored to the snapshot if it did not improve on the initial residual.
      @param[in]     a_rhs Right-hand side, used to evaluate the current residual.
      @return True if the current solution was kept, false if it was restored.
    */
    bool
    keepOrRestore(Vec& a_phi, const Vec& a_rhs);

  protected:
    /**
      @brief Lazily allocate the persistent residual-correction vectors r0/dphi from a template hierarchy.
      @param[in] a_template A hierarchy whose layout the work vectors mirror (e.g. the rhs).
    */
    void
    allocateWork(const Vec& a_template);

    /**
      @brief Free the persistent work vectors and undefine the custom solvers.
    */
    void
    freeWork();

    /**
      @brief Multigrid solver (not owned).
    */
    AMRMultiGrid<LevelData<T>>* m_mg = nullptr;

    /**
      @brief Multigrid-cycle-as-preconditioner adapter.
    */
    AMRMultigridKrylovOp<T> m_op;

    /**
      @brief Custom BiCGStab solver (lazily defined on first BiCGStab solve).
    */
    KrylovBiCGStab<T> m_bicg;

    /**
      @brief Custom GMRES solver (lazily defined on first GMRES solve).
    */
    KrylovGMRES<T> m_gmres;

    /**
      @brief Iteration count of the most recent Krylov solve, or -1 if none has run.
    */
    int m_lastKrylovIters = -1;

    /**
      @brief Coarsest level.
    */
    int m_lbase = 0;

    /**
      @brief Finest level.
    */
    int m_lmax = 0;

    /**
      @brief Whether define() has been called.
    */
    bool m_isDefined = false;

    /**
      @brief Whether r0/dphi are allocated for the current grid.
    */
    bool m_workAllocated = false;

    /**
      @brief Persistent initial residual.
    */
    Vec m_r0;

    /**
      @brief Persistent correction.
    */
    Vec m_dphi;

    /**
      @brief Snapshot of the chain's starting guess (the initial-residual solution) for the warm-start rule.
    */
    Vec m_phiInit;

    /**
      @brief Residual of the snapshot (the initial residual of the fallback chain).
    */
    Real m_initResidual = 0.0;
  };

} // namespace KrylovMultigrid

#include <CD_NamespaceFooter.H>

#include <CD_KrylovMultigridImplem.H>

#endif
