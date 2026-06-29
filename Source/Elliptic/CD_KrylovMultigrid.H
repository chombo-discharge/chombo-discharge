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

// Std includes
#include <string>

// Chombo includes
#include <AMRMultiGrid.H>
#include <LevelData.H>
#include <Vector.H>

// Our includes
#include <CD_AMRMultigridKrylovOp.H>
#include <CD_KrylovBiCGStab.H>
#include <CD_KrylovGMRES.H>
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

  /*!
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

  /*!
    @brief Stateful owner of the outer-Krylov solver path.
    @details Holds the V-cycle adapter (AMRMultigridKrylovOp), the chombo-discharge BiCGStab/GMRES solvers, and the
             persistent residual-correction work vectors. define() is called once per regrid; solve() reuses the
             adapter, the (lazily defined) custom solvers and the work vectors across all subsequent solves, so there
             is no per-solve work-vector allocation. solve() runs in residual-correction form so the inhomogeneous
             boundary/jump contribution enters only through the initial residual r0 = rhs - L(phi); the matvec and
             preconditioner are applied with homogeneous boundary conditions, and warm starts are honoured.
    @note T is the FAB type (EBCellFAB or MFCellFAB).
  */
  template <class T>
  class Driver
  {
  public:
    using Vec = Vector<LevelData<T>*>; /*!< Multilevel data type. */

    /*! @brief Default constructor. Call define() before use. */
    Driver() noexcept = default;

    /*! @brief Destructor. Frees the persistent work vectors and custom solvers. */
    ~Driver() noexcept;

    /*!
      @brief Define the V-cycle adapter and reset the persistent work. Call once per regrid.
      @param[in] a_mg         Fully set-up multigrid solver (reused for matvec and preconditioning)
      @param[in] a_grids      Grids on each level (index 0..a_lmax)
      @param[in] a_refRat     Refinement ratio to the next finer level (index 0..a_lmax)
      @param[in] a_dx         Grid spacing on each level (index 0..a_lmax)
      @param[in] a_lbase      Coarsest level in the solve
      @param[in] a_lmax       Finest level in the solve
      @param[in] a_numVCycles V-cycles per preconditioner invocation
    */
    void
    define(AMRMultiGrid<LevelData<T>>*      a_mg,
           const Vector<DisjointBoxLayout>& a_grids,
           const Vector<int>&               a_refRat,
           const Vector<Real>&              a_dx,
           const int                        a_lbase,
           const int                        a_lmax,
           const int                        a_numVCycles);

    /*!
      @brief Solve L(phi) = rhs with the chosen outer Krylov method preconditioned by the V-cycle.
      @param[in,out] a_phi      Solution (in: initial guess unless a_zeroPhi; out: solution)
      @param[in]     a_rhs      Right-hand side
      @param[in]     a_zeroPhi  If true, start from phi = 0
      @param[in]     a_type     Which Krylov method to run (GMRES or BiCGStab)
      @param[in]     a_settings Outer-Krylov settings (tolerance, iterations, restart, verbosity)
      @return True if the Krylov solver reported convergence.
    */
    bool
    solve(Vec& a_phi, const Vec& a_rhs, const bool a_zeroPhi, const SolverType a_type, const Settings& a_settings);

  protected:
    /*!
      @brief Lazily allocate the persistent residual-correction vectors r0/dphi from a template hierarchy.
      @param[in] a_template A hierarchy whose layout the work vectors mirror (e.g. the rhs).
    */
    void
    allocateWork(const Vec& a_template);

    /*! @brief Free the persistent work vectors and undefine the custom solvers. */
    void
    freeWork();

    /*! @brief Multigrid solver (not owned). */
    AMRMultiGrid<LevelData<T>>* m_mg = nullptr;

    /*! @brief V-cycle-as-preconditioner adapter. */
    AMRMultigridKrylovOp<T> m_op;

    /*! @brief Custom BiCGStab solver (lazily defined on first BiCGStab solve). */
    KrylovBiCGStab<T> m_bicg;

    /*! @brief Custom GMRES solver (lazily defined on first GMRES solve). */
    KrylovGMRES<T> m_gmres;

    int  m_lbase         = 0;     /*!< Coarsest level. */
    int  m_lmax          = 0;     /*!< Finest level. */
    bool m_isDefined     = false; /*!< Whether define() has been called. */
    bool m_workAllocated = false; /*!< Whether r0/dphi are allocated for the current grid. */

    Vec m_r0;   /*!< Persistent initial residual. */
    Vec m_dphi; /*!< Persistent correction. */
  };

} // namespace KrylovMultigrid

#include <CD_NamespaceFooter.H>

#include <CD_KrylovMultigridImplem.H>

#endif
