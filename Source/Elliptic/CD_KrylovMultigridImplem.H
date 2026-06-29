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

// Chombo includes
#include <CH_Timer.H>
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
                       const Settings&              a_settings)
{
  CH_TIME("KrylovMultigrid::solve");

  using Vec = Vector<LevelData<T>*>;

  if (a_zeroPhi) {
    a_op.setToZero(a_phi);
  }

  // Initial residual r0 = rhs - L_inhomogeneous(phi). This is the only place the inhomogeneous boundary/jump
  // contribution enters; the Krylov matvec and preconditioner run homogeneously.
  Vec r0;
  a_op.create(r0, a_rhs);
  a_mg->computeAMRResidual(r0, a_phi, a_rhs, a_lmax, a_lbase, false, false); // r0 = L(phi) - rhs
  a_op.scale(r0, -1.0);                                                      // r0 = rhs - L(phi)

  // Correction dphi, solved homogeneously: L_homogeneous(dphi) = r0.
  Vec dphi;
  a_op.create(dphi, a_phi);
  a_op.setToZero(dphi);

  bool converged = false;

  switch (a_settings.type) {
  case SolverType::GMRES: {
    GMRESSolver<Vec> gmres;
    gmres.define(&a_op, true);
    gmres.m_verbosity = a_settings.verbosity;
    gmres.m_imax      = a_settings.maxIter;
    gmres.m_reps      = a_settings.eps;
    gmres.setRestartLen(a_settings.restart);

    gmres.solve(dphi, r0);

    converged = (gmres.m_exitStatus == 1);

    break;
  }
  case SolverType::BiCGStab: {
    BiCGStabSolver<Vec> bicg;
    bicg.define(&a_op, true);
    bicg.m_verbosity = a_settings.verbosity;
    bicg.m_imax      = a_settings.maxIter;
    bicg.m_eps       = a_settings.eps;

    bicg.solve(dphi, r0);

    converged = (bicg.m_exitStatus == 1);

    break;
  }
  default: {
    MayDay::Error("KrylovMultigrid::solve - SolverType must be GMRES or BiCGStab");

    break;
  }
  }

  // phi += dphi.
  a_op.incr(a_phi, dphi, 1.0);

  a_op.clear(r0);
  a_op.clear(dphi);

  return converged;
}

#include <CD_NamespaceFooter.H>

#endif
