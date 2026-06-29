/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_KrylovGMRESImplem.H
  @brief  Implementation of CD_KrylovGMRES.H
  @author Robert Marskar
*/

#ifndef CD_KrylovGMRESImplem_H
#define CD_KrylovGMRESImplem_H

// Std includes
#include <algorithm>
#include <cmath>

// Chombo includes
#include <CH_Timer.H>
#include <parstream.H>

// Our includes
#include <CD_KrylovGMRES.H>
#include <CD_NamespaceHeader.H>

template <class T>
KrylovGMRES<T>::~KrylovGMRES() noexcept
{
  this->undefine();
}

template <class T>
void
KrylovGMRES<T>::define(Op* a_op, const Vec& a_template) noexcept
{
  CH_TIME("KrylovGMRES::define");

  this->undefine();

  m_op           = a_op;
  m_allocRestart = std::max(1, m_restart);

  m_basis.resize(m_allocRestart + 1);
  for (int i = 0; i < m_basis.size(); i++) {
    m_op->create(m_basis[i], a_template);
  }
  m_op->create(m_w, a_template);
  m_op->create(m_z, a_template);
  m_op->create(m_c, a_template);
  m_op->create(m_r, a_template);

  m_H    = Vector<Real>((m_allocRestart + 1) * m_allocRestart, 0.0);
  m_cs   = Vector<Real>(m_allocRestart, 0.0);
  m_sn   = Vector<Real>(m_allocRestart, 0.0);
  m_g    = Vector<Real>(m_allocRestart + 1, 0.0);
  m_y    = Vector<Real>(m_allocRestart, 0.0);
  m_hcol = Vector<Real>(m_allocRestart, 0.0);
  m_corr = Vector<Real>(m_allocRestart, 0.0);

  m_isDefined = true;
}

template <class T>
void
KrylovGMRES<T>::undefine() noexcept
{
  if (m_isDefined) {
    for (int i = 0; i < m_basis.size(); i++) {
      m_op->clear(m_basis[i]);
    }
    m_basis.clear();

    m_op->clear(m_w);
    m_op->clear(m_z);
    m_op->clear(m_c);
    m_op->clear(m_r);

    m_isDefined = false;
  }
}

template <class T>
bool
KrylovGMRES<T>::isDefined() const noexcept
{
  return m_isDefined;
}

template <class T>
bool
KrylovGMRES<T>::solve(Vec& a_phi, const Vec& a_rhs) noexcept
{
  CH_TIME("KrylovGMRES::solve");

  CH_assert(m_isDefined);

  Op& op = *m_op;

  m_exitStatus = -1;

  const int  m    = m_allocRestart; // restart length
  const int  ldH  = m + 1;          // Hessenberg leading dimension
  const auto Hidx = [ldH](const int r, const int c) -> int {
    return r + c * ldH;
  };

  // True (right-preconditioning preserves it) initial residual r = rhs - A*phi.
  op.applyOp(m_r, a_phi, true);
  op.axby(m_r, a_rhs, m_r, 1.0, -1.0);

  const Real beta0 = op.norm(m_r, 2);
  if (beta0 == 0.0) {
    m_exitStatus   = 1;
    m_residualNorm = 0.0;
    m_iterations   = 0;
    return true;
  }

  Real rnorm = beta0;
  int  iter  = 0;

  while (iter < m_maxIter && rnorm > m_eps * beta0) {

    // Start a restart cycle: v_0 = r / ||r||.
    const Real beta = op.norm(m_r, 2);
    op.assign(m_basis[0], m_r);
    op.scale(m_basis[0], 1.0 / beta);

    for (int i = 0; i <= m; i++) {
      m_g[i] = 0.0;
    }
    m_g[0] = beta;

    int jUsed = 0;

    for (int j = 0; j < m && iter < m_maxIter; j++) {
      iter++;

      // Arnoldi step on A*K^{-1}: w = A (K^{-1} v_j). The preconditioned vector is not stored; the solution update
      // is preconditioned once per restart cycle below.
      op.preCond(m_z, m_basis[j]);
      op.applyOp(m_w, m_z, true);

      // Classical Gram-Schmidt against v_0..v_j -- one fused reduction (mDotProduct) instead of j+1 dot products.
      // m_hcol is persistent scratch; only entries 0..j are written and read this step.
      op.mDotProduct(m_w, j + 1, &m_basis[0], &m_hcol[0]);
      for (int i = 0; i <= j; i++) {
        op.incr(m_w, m_basis[i], -m_hcol[i]);
        m_H[Hidx(i, j)] = m_hcol[i];
      }

      if (m_reorthogonalize) {
        op.mDotProduct(m_w, j + 1, &m_basis[0], &m_corr[0]);
        for (int i = 0; i <= j; i++) {
          op.incr(m_w, m_basis[i], -m_corr[i]);
          m_H[Hidx(i, j)] += m_corr[i];
        }
      }

      const Real hjp1     = op.norm(m_w, 2);
      m_H[Hidx(j + 1, j)] = hjp1;
      if (hjp1 != 0.0) {
        op.assign(m_basis[j + 1], m_w);
        op.scale(m_basis[j + 1], 1.0 / hjp1);
      }

      // Apply existing Givens rotations to the new Hessenberg column.
      for (int k = 0; k < j; k++) {
        const Real t        = m_cs[k] * m_H[Hidx(k, j)] + m_sn[k] * m_H[Hidx(k + 1, j)];
        m_H[Hidx(k + 1, j)] = -m_sn[k] * m_H[Hidx(k, j)] + m_cs[k] * m_H[Hidx(k + 1, j)];
        m_H[Hidx(k, j)]     = t;
      }

      // New Givens rotation eliminating H[j+1][j].
      const Real hjj   = m_H[Hidx(j, j)];
      const Real hj1j  = m_H[Hidx(j + 1, j)];
      const Real denom = std::sqrt(hjj * hjj + hj1j * hj1j);
      if (denom == 0.0) {
        break;
      }
      m_cs[j]             = hjj / denom;
      m_sn[j]             = hj1j / denom;
      m_H[Hidx(j, j)]     = denom;
      m_H[Hidx(j + 1, j)] = 0.0;

      const Real gtemp = m_cs[j] * m_g[j];
      m_g[j + 1]       = -m_sn[j] * m_g[j];
      m_g[j]           = gtemp;

      rnorm = std::abs(m_g[j + 1]);
      jUsed = j + 1;

      if (m_verbosity >= 4) {
        pout() << "      KrylovGMRES:: iteration = " << iter << ", residual = " << rnorm << endl;
      }

      if (rnorm <= m_eps * beta0 || hjp1 == 0.0) {
        break;
      }
    }

    // Back-substitute the (jUsed x jUsed) upper-triangular system H y = g.
    for (int i = jUsed - 1; i >= 0; i--) {
      Real s = m_g[i];
      for (int k = i + 1; k < jUsed; k++) {
        s -= m_H[Hidx(i, k)] * m_y[k];
      }
      m_y[i] = (m_H[Hidx(i, i)] != 0.0) ? s / m_H[Hidx(i, i)] : 0.0;
    }

    // Form the unpreconditioned update c = sum_k y_k v_k, precondition once, and apply: phi += K^{-1} c.
    op.setToZero(m_c);
    for (int k = 0; k < jUsed; k++) {
      op.incr(m_c, m_basis[k], m_y[k]);
    }
    op.preCond(m_z, m_c);
    op.incr(a_phi, m_z, 1.0);

    // Recompute the true residual for the restart/convergence decision.
    op.applyOp(m_r, a_phi, true);
    op.axby(m_r, a_rhs, m_r, 1.0, -1.0);
    rnorm = op.norm(m_r, 2);

    if (jUsed == 0) { // no progress (degenerate) -- avoid spinning
      break;
    }
  }

  m_exitStatus   = (rnorm <= m_eps * beta0) ? 1 : 3;
  m_residualNorm = rnorm;
  m_iterations   = iter;

  return (m_exitStatus == 1);
}

#include <CD_NamespaceFooter.H>

#endif
