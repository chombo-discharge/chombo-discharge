/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_KrylovBiCGStabImplem.H
  @brief  Implementation of CD_KrylovBiCGStab.H
  @author Robert Marskar
*/

#ifndef CD_KrylovBiCGStabImplem_H
#define CD_KrylovBiCGStabImplem_H

// Chombo includes
#include <CH_Timer.H>
#include <parstream.H>

// Our includes
#include <CD_KrylovBiCGStab.H>
#include <CD_NamespaceHeader.H>

template <class T>
KrylovBiCGStab<T>::~KrylovBiCGStab() noexcept
{
  this->undefine();
}

template <class T>
void
KrylovBiCGStab<T>::define(Op* a_op, const Vec& a_template) noexcept
{
  CH_TIME("KrylovBiCGStab::define");

  this->undefine();

  m_op = a_op;

  // Allocate the work vectors once; they are reused across all solves until the next define() (regrid).
  m_op->create(m_r, a_template);
  m_op->create(m_rTilde, a_template);
  m_op->create(m_e, a_template);
  m_op->create(m_p, a_template);
  m_op->create(m_pTilde, a_template);
  m_op->create(m_sTilde, a_template);
  m_op->create(m_t, a_template);
  m_op->create(m_v, a_template);

  m_isDefined = true;
}

template <class T>
void
KrylovBiCGStab<T>::undefine() noexcept
{
  if (m_isDefined) {
    m_op->clear(m_r);
    m_op->clear(m_rTilde);
    m_op->clear(m_e);
    m_op->clear(m_p);
    m_op->clear(m_pTilde);
    m_op->clear(m_sTilde);
    m_op->clear(m_t);
    m_op->clear(m_v);

    m_isDefined = false;
  }
}

template <class T>
bool
KrylovBiCGStab<T>::isDefined() const noexcept
{
  return m_isDefined;
}

template <class T>
bool
KrylovBiCGStab<T>::solve(Vec& a_phi, const Vec& a_rhs) noexcept
{
  CH_TIME("KrylovBiCGStab::solve");

  CH_assert(m_isDefined);

  Op& op = *m_op;

  m_exitStatus = -1;

  // Initialize the (homogeneous) residual r = rhs - A*phi and the BiCGStab state. The correction is accumulated in
  // m_e and added to a_phi at the very end, so warm starts (a_phi != 0) are honoured.
  auto computeResidual = [&]() -> void {
    op.applyOp(m_r, a_phi, true);        // r = A*phi
    op.axby(m_r, a_rhs, m_r, 1.0, -1.0); // r = rhs - A*phi
  };

  computeResidual();
  op.assign(m_rTilde, m_r);
  op.setToZero(m_e);
  op.setToZero(m_p);
  op.setToZero(m_v);

  const Real rnorm0 = op.norm(m_r, 2);

  if (rnorm0 == 0.0) {
    m_exitStatus   = 1;
    m_residualNorm = 0.0;
    m_iterations   = 0;
    return true;
  }

  Real rnorm    = rnorm0;
  Real rho      = 1.0;
  Real alpha    = 1.0;
  Real omega    = 1.0;
  int  iter     = 0;
  int  restarts = 0;

  // Fold the work-so-far into a_phi and restart the iteration from the current solution. Used on breakdown.
  auto restart = [&]() -> void {
    op.incr(a_phi, m_e, 1.0);
    op.setToZero(m_e);
    computeResidual();
    op.assign(m_rTilde, m_r);
    op.setToZero(m_p);
    op.setToZero(m_v);
    rho   = 1.0;
    alpha = 1.0;
    omega = 1.0;
    rnorm = op.norm(m_r, 2);
    restarts++;
  };

  // Treat a breakdown: restart if we still can, otherwise stop with status 2.
  auto breakdown = [&]() -> bool {
    if (restarts < m_numRestarts) {
      restart();
      return false; // continue
    }
    m_exitStatus = 2;
    return true; // stop
  };

  while (iter < m_maxIter && rnorm > m_eps * rnorm0) {
    const Real rhoNew = op.dotProduct(m_rTilde, m_r);
    if (rhoNew == 0.0) {
      if (breakdown()) {
        break;
      }
      continue;
    }

    // p = r + beta*(p - omega*v). With the initial/post-restart state (p = v = 0, rho = alpha = omega = 1) this
    // correctly reduces to p = r, so no special first-iteration branch is needed.
    const Real beta = (rhoNew / rho) * (alpha / omega);
    op.incr(m_p, m_v, -omega);
    op.scale(m_p, beta);
    op.incr(m_p, m_r, 1.0);

    op.preCond(m_pTilde, m_p);       // p_tilde = K^{-1} p
    op.applyOp(m_v, m_pTilde, true); // v = A p_tilde

    const Real rtv = op.dotProduct(m_rTilde, m_v);
    if (rtv == 0.0) {
      if (breakdown()) {
        break;
      }
      continue;
    }
    alpha = rhoNew / rtv;

    op.incr(m_e, m_pTilde, alpha); // e += alpha*p_tilde
    op.incr(m_r, m_v, -alpha);     // r = s = r - alpha*v

    rnorm = op.norm(m_r, 2);
    if (rnorm <= m_eps * rnorm0) { // converged on the half-step s
      m_exitStatus = 1;
      break;
    }

    op.preCond(m_sTilde, m_r);       // s_tilde = K^{-1} s
    op.applyOp(m_t, m_sTilde, true); // t = A s_tilde

    // omega = (t,s)/(t,t) -- both inner products in a single fused MPI reduction.
    Real      dots[2];
    const Vec rhsPair[2] = {m_r, m_t};
    op.mDotProduct(m_t, 2, rhsPair, dots);
    const Real tDotT = dots[1];
    if (tDotT == 0.0) {
      if (breakdown()) {
        break;
      }
      continue;
    }
    omega = dots[0] / tDotT;

    op.incr(m_e, m_sTilde, omega); // e += omega*s_tilde
    op.incr(m_r, m_t, -omega);     // r = s - omega*t

    const Real rnormPrev = rnorm;
    rnorm                = op.norm(m_r, 2);

    rho = rhoNew;

    if (m_verbosity >= 4) {
      pout() << "      KrylovBiCGStab:: iteration = " << iter << ", residual = " << rnorm
             << ", rate = " << ((rnormPrev > 0.0) ? rnorm / rnormPrev : 0.0) << endl;
    }

    if (rnorm <= m_eps * rnorm0) {
      m_exitStatus = 1;
      break;
    }

    if (omega == 0.0) {
      if (breakdown()) {
        break;
      }
      continue;
    }

    iter++;
  }

  if (m_exitStatus == -1) {
    m_exitStatus = (rnorm <= m_eps * rnorm0) ? 1 : 3; // converged on the loop guard, else max iterations
  }

  m_residualNorm = rnorm;
  m_iterations   = iter;

  op.incr(a_phi, m_e, 1.0); // commit the accumulated correction

  return (m_exitStatus == 1);
}

#include <CD_NamespaceFooter.H>

#endif
