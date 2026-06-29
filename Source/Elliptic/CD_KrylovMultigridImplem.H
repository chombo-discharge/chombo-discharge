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
#include <GMRESSolver.H>
#include <BiCGStabSolver.H>

// Our includes
#include <CD_KrylovMultigrid.H>
#include <CD_NamespaceHeader.H>

template <class T>
bool
KrylovMultigrid::solve(AMRMultiGrid<LevelData<T>>*  a_mg,
                       AMRMultigridKrylovOp<T>&     a_op,
                       Vector<LevelData<T>*>&       a_phi,
                       const Vector<LevelData<T>*>& a_rhs,
                       const int                    a_lbase,
                       const int                    a_lmax,
                       const bool                   a_zeroPhi,
                       const SolverType             a_type,
                       const Settings&              a_settings)
{
  CH_TIME("KrylovMultigrid::solve");

  using Vec = Vector<LevelData<T>*>;

  if (a_zeroPhi) {
    a_op.setToZero(a_phi);
  }

  // The preconditioner drives AMRMultiGrid::AMRVCycle directly, which needs the solver's multigrid holders and its
  // *intermediate* multigrid depths. AMRMultiGrid sets these up in init() and tears them down in revert(); the host
  // solvers call revert() after every solve (and the stand-alone GMG path re-init()s every solve via
  // solveNoInitResid). We must therefore (re-)initialise here before every Krylov solve -- otherwise the second and
  // later solves on a given grid run the V-cycle with reverted (full-coarsening) depths and return a wrong
  // correction. This mirrors what AMRMultiGrid::solve()/solveNoInitResid() do internally for the stand-alone path.
  a_mg->init(a_phi, a_rhs, a_lmax, a_lbase);
  a_mg->setBottomSolver(a_lmax, a_lbase);

  // Initial residual r0 = rhs - L_inhomogeneous(phi). This is the only place the inhomogeneous boundary/jump
  // contribution enters; the Krylov matvec and preconditioner run homogeneously. computeAMRResidual already returns
  // the residual in Chombo's convention (rhs - L(phi)), so it must NOT be negated -- negating it flips the sign of
  // the correction and yields phi = 2*phi0 - phi* for a warm start (a wrong field; only the cold phi0 = 0 case
  // survives, as -phi* happens to have the right field magnitude).
  Vec r0;
  a_op.create(r0, a_rhs);
  a_mg->computeAMRResidual(r0, a_phi, a_rhs, a_lmax, a_lbase, false, false); // r0 = rhs - L(phi)

  // Initial residual norm (same norm the Krylov solver uses internally) -- kept for the convergence-metric printout.
  const Real rNorm0 = a_op.norm(r0, 2);

  // Correction dphi, solved homogeneously: L_homogeneous(dphi) = r0.
  Vec dphi;
  a_op.create(dphi, a_phi);
  a_op.setToZero(dphi);

  bool converged  = false;
  int  exitStatus = -1;

  // Keep the Chombo solver output terse. At verbosity >= 4 both GMRES and BiCGStab print a one-line-per-iteration
  // convergence history (akin to AMRMultiGrid's residual/rate report) plus restart messages; at >= 5 BiCGStab also
  // dumps its internal scalars (alpha/beta/omega/rho) every iteration. Cap at 4 so we only ever get the terse history.
  const int solverVerbosity = std::min(a_settings.verbosity, 4);

  switch (a_type) {
  case SolverType::GMRES: {
    GMRESSolver<Vec> gmres;
    gmres.define(&a_op, true);
    gmres.m_verbosity = solverVerbosity;
    gmres.m_imax      = a_settings.maxIter;
    gmres.m_reps      = a_settings.eps;
    gmres.setRestartLen(a_settings.restart);

    gmres.solve(dphi, r0);

    exitStatus = gmres.m_exitStatus;
    converged  = (exitStatus == 1);

    break;
  }
  case SolverType::BiCGStab: {
    BiCGStabSolver<Vec> bicg;
    bicg.define(&a_op, true);
    bicg.m_verbosity = solverVerbosity;
    bicg.m_imax      = a_settings.maxIter;
    bicg.m_eps       = a_settings.eps;

    bicg.solve(dphi, r0);

    exitStatus = bicg.m_exitStatus;
    converged  = (exitStatus == 1);

    break;
  }
  default: {
    MayDay::Error("KrylovMultigrid::solve - SolverType must be GMRES or BiCGStab");

    break;
  }
  }

  // Convergence metrics. The verbosity is the host's gmg_verbosity (reused for the Krylov solvers). At this level we
  // print a one-line summary: the residual reduction in the Krylov norm and the solver exit status. The full
  // per-iteration history is emitted by the Chombo solver itself at higher verbosity.
  if (a_settings.verbosity >= 1) {
    Vec res;
    a_op.create(res, r0);
    a_op.applyOp(res, dphi, true); // res = L_homogeneous(dphi)
    a_op.scale(res, -1.0);
    a_op.incr(res, r0, 1.0); // res = r0 - L_homogeneous(dphi) = final Krylov residual
    const Real rNorm = a_op.norm(res, 2);
    a_op.clear(res);

    const Real reduction = (rNorm0 > 0.0) ? rNorm / rNorm0 : 0.0;

    pout() << "KrylovMultigrid::solve - " << solverTypeName(a_type) << (converged ? " converged" : " did NOT converge")
           << " (|r0| = " << rNorm0 << ", |r| = " << rNorm << ", reduction = " << reduction
           << ", status = " << exitStatus << ")" << endl;
  }

  // phi += dphi.
  a_op.incr(a_phi, dphi, 1.0);

  a_op.clear(r0);
  a_op.clear(dphi);

  return converged;
}

#include <CD_NamespaceFooter.H>

#endif
