/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_KrylovMultigridImplem.H
  @brief  Implementation of CD_KrylovMultigrid.H
  @author Robert Marskar
*/

#ifndef CD_KrylovMultigridImplem_H
#define CD_KrylovMultigridImplem_H

// Std includes
#include <algorithm>

// Chombo includes
#include <CH_Timer.H>
#include <parstream.H>

// Our includes
#include <CD_KrylovMultigrid.H>
#include <CD_NamespaceHeader.H>

namespace KrylovMultigrid {

  template <class T>
  Solver<T>::~Solver() noexcept
  {
    this->freeWork();
  }

  template <class T>
  void
  Solver<T>::freeWork()
  {
    m_bicg.undefine();
    m_gmres.undefine();

    if (m_workAllocated) {
      m_op.clear(m_r0);
      m_op.clear(m_dphi);
      m_op.clear(m_phiInit);
      m_workAllocated = false;
    }
  }

  template <class T>
  void
  Solver<T>::allocateWork(const Vec& a_template)
  {
    m_op.create(m_r0, a_template);
    m_op.create(m_dphi, a_template);
    m_op.create(m_phiInit, a_template);
    m_workAllocated = true;
  }

  template <class T>
  void
  Solver<T>::snapshotGuess(Vec& a_phi, const Vec& a_rhs)
  {
    CH_TIME("KrylovMultigrid::Solver::snapshotGuess");

    if (!m_workAllocated) {
      this->allocateWork(a_phi);
    }

    m_op.assign(m_phiInit, a_phi);
    m_initResidual = m_mg->computeAMRResidual(a_phi, a_rhs, m_lmax, m_lbase);
  }

  template <class T>
  bool
  Solver<T>::keepOrRestore(Vec& a_phi, const Vec& a_rhs)
  {
    CH_TIME("KrylovMultigrid::Solver::keepOrRestore");

    const Real r = m_mg->computeAMRResidual(a_phi, a_rhs, m_lmax, m_lbase);

    // Keep the just-failed solve only if it improved on the chain's initial residual; otherwise revert to the
    // initial-residual solution.
    if (r >= m_initResidual) {
      m_op.assign(a_phi, m_phiInit);

      return false;
    }

    return true;
  }

  template <class T>
  void
  Solver<T>::define(AMRMultiGrid<LevelData<T>>*                            a_mg,
                    const Vector<RefCountedPtr<LevelData<BaseFab<bool>>>>& a_validCells,
                    const Vector<Real>&                                    a_dx,
                    const int                                              a_lbase,
                    const int                                              a_lmax,
                    const int                                              a_numVCycles)
  {
    CH_TIME("KrylovMultigrid::Solver::define");

    // Release any storage tied to the previous grid before redefining the adapter for the new one. The work vectors
    // are re-allocated lazily on the first solve (when a template hierarchy is available).
    this->freeWork();

    m_mg    = a_mg;
    m_lbase = a_lbase;
    m_lmax  = a_lmax;

    m_op.define(a_mg, a_validCells, a_dx, a_lbase, a_lmax, a_numVCycles);

    m_isDefined = true;
  }

  template <class T>
  bool
  Solver<T>::solve(Vec&             a_phi,
                   const Vec&       a_rhs,
                   const bool       a_zeroPhi,
                   const SolverType a_type,
                   const Settings&  a_settings)
  {
    CH_TIME("KrylovMultigrid::Solver::solve");

    CH_assert(m_isDefined);

    // Allocate the persistent residual-correction vectors once per grid (first solve after a define()/regrid).
    if (!m_workAllocated) {
      this->allocateWork(a_rhs);
    }

    if (a_zeroPhi) {
      m_op.setToZero(a_phi);
    }

    // The preconditioner drives AMRMultiGrid::AMRVCycle directly, which needs the solver's *intermediate* multigrid
    // depths. AMRMultiGrid sets these up in init() and tears them down in revert(); the host solvers revert() after
    // every solve (and the stand-alone GMG path re-init()s every solve), so we re-initialise here -- once per field
    // solve, exactly as the GMG path does.
    m_mg->init(a_phi, a_rhs, m_lmax, m_lbase);
    m_mg->setBottomSolver(m_lmax, m_lbase);

    // Zero-residual ||rhs - L(0)|| in the adapter norm. The printed residuals are normalized against this so they are
    // on a common, solver-independent scale (a relative residual), comparable to the convergence tolerance. m_dphi is
    // the (zero) correction accumulator, reused here as the zero solution.
    m_op.setToZero(m_dphi);
    m_mg->computeAMRResidual(m_r0, m_dphi, a_rhs, m_lmax, m_lbase, false, false);
    const Real zeroResid     = m_op.norm(m_r0, 2);
    const Real residualScale = (zeroResid > 0.0) ? zeroResid : 1.0;

    // Initial residual r0 = rhs - L_inhomogeneous(phi). This is the only inhomogeneous evaluation; the Krylov matvec
    // and preconditioner run with homogeneous boundary conditions. computeAMRResidual returns rhs - L(phi) in
    // Chombo's convention, so it must NOT be negated.
    m_mg->computeAMRResidual(m_r0, a_phi, a_rhs, m_lmax, m_lbase, false, false);

    const Real rNorm0 = m_op.norm(m_r0, 2);

    bool converged    = false;
    int  exitStatus   = -1;
    Real residualNorm = rNorm0;

    switch (a_type) {
    case SolverType::GMRES: {
      const int restart = std::max(1, a_settings.restart);
      if (!m_gmres.isDefined() || m_gmres.m_restart != restart) {
        m_gmres.m_restart = restart;
        m_gmres.define(&m_op, a_rhs);
      }
      m_gmres.m_eps           = a_settings.eps;
      m_gmres.m_maxIter       = a_settings.maxIter;
      m_gmres.m_verbosity     = a_settings.verbosity;
      m_gmres.m_residualScale = residualScale;

      converged         = m_gmres.solve(m_dphi, m_r0);
      exitStatus        = m_gmres.m_exitStatus;
      residualNorm      = m_gmres.m_residualNorm;
      m_lastKrylovIters = m_gmres.m_iterations;

      break;
    }
    case SolverType::BiCGStab: {
      if (!m_bicg.isDefined()) {
        m_bicg.define(&m_op, a_rhs);
      }
      m_bicg.m_eps           = a_settings.eps;
      m_bicg.m_maxIter       = a_settings.maxIter;
      m_bicg.m_verbosity     = a_settings.verbosity;
      m_bicg.m_residualScale = residualScale;

      converged         = m_bicg.solve(m_dphi, m_r0);
      exitStatus        = m_bicg.m_exitStatus;
      residualNorm      = m_bicg.m_residualNorm;
      m_lastKrylovIters = m_bicg.m_iterations;

      break;
    }
    default: {
      MayDay::Error("KrylovMultigrid::Solver::solve - SolverType must be GMRES or BiCGStab");

      break;
    }
    }

    // One-line convergence summary (the per-iteration history, if any, is printed by the solver itself). The residual
    // is reported relative to the zero-residual (rel.res = ||r|| / ||r_zero||), the same solver-independent scale used
    // in the per-iteration output, with the absolute norms kept for reference. The verbosity is the host's
    // gmg_verbosity, reused for the Krylov path.
    if (a_settings.verbosity >= 1) {
      const Real relRes = residualNorm / residualScale;

      pout() << "KrylovMultigrid - " << solverTypeName(a_type) << (converged ? " converged" : " did NOT converge")
             << " (rel.res = " << relRes << ", |r| = " << residualNorm << ", |r_zero| = " << zeroResid
             << ", iters = " << m_lastKrylovIters << ", status = " << exitStatus << ")" << endl;
    }

    // phi += dphi.
    m_op.incr(a_phi, m_dphi, 1.0);

    return converged;
  }

} // namespace KrylovMultigrid

#include <CD_NamespaceFooter.H>

#endif
