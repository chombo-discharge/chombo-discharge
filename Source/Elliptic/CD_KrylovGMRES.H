/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_KrylovGMRES.H
  @brief  Restarted, right-preconditioned GMRES over an AMR hierarchy, with persistent basis and fused reductions
  @author Robert Marskar
*/

#ifndef CD_KrylovGMRES_H
#define CD_KrylovGMRES_H

// Std includes
#include <vector>

// Chombo includes
#include <LevelData.H>
#include <Vector.H>

// Our includes
#include <CD_AMRMultigridKrylovOp.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Restarted GMRES(m) acting on an AMR hierarchy through AMRMultigridKrylovOp.
  @details A chombo-discharge replacement for Chombo's GMRESSolver on the outer Krylov path, written for scale-out.
           The Arnoldi basis (m_restart + 1 hierarchies) and the few scratch vectors are allocated once in define()
           and reused across solves (no per-solve allocation). Orthogonalization uses classical Gram-Schmidt issued
           through the adapter's batched mDotProduct, so each Arnoldi step costs a single MPI reduction (plus one for
           the norm) regardless of the AMR level count -- as opposed to the modified Gram-Schmidt / per-vector dot
           products that scale poorly. Right preconditioning keeps the tracked residual the true residual; the
           preconditioned solution update is formed once per restart cycle (one extra preCond), so the basis vectors
           need not be stored in preconditioned form.
  @note T is the FAB type (EBCellFAB or MFCellFAB).
*/
template <class T>
class KrylovGMRES
{
public:
  using Op  = AMRMultigridKrylovOp<T>; /*!< Operator/adapter type. */
  using Vec = Vector<LevelData<T>*>;   /*!< Multilevel data type the solver acts on. */

  /*! @brief Relative residual tolerance (against the initial residual norm). */
  Real m_eps = 1.E-10;

  /*! @brief Maximum number of iterations (total Arnoldi steps across restart cycles). */
  int m_maxIter = 50;

  /*! @brief Restart length m in GMRES(m). Must be set before define() (sizes the persistent basis). */
  int m_restart = 30;

  /*! @brief Whether to do one classical Gram-Schmidt reorthogonalization pass (robustness; on by default). */
  bool m_reorthogonalize = true;

  /*! @brief Verbosity; at >= 4 a terse one-line-per-iteration residual history is printed. */
  int m_verbosity = -1;

  /*! @brief Exit status after solve(): 1 = converged, 3 = max iterations, -1 = not run. */
  int m_exitStatus = -1;

  /*! @brief Final residual norm after solve() (in the adapter inner-product norm). */
  Real m_residualNorm = -1.0;

  /*! @brief Default constructor. Call define() before use. */
  KrylovGMRES() noexcept = default;

  /*! @brief Destructor. Frees the persistent basis and scratch. */
  ~KrylovGMRES() noexcept;

  /*!
    @brief Allocate the persistent basis (m_restart + 1 vectors) and scratch. Call once per regrid.
    @param[in] a_op       The (already define()d) operator/adapter to solve against
    @param[in] a_template A hierarchy whose layout the work vectors mirror (e.g. the rhs)
  */
  void
  define(Op* a_op, const Vec& a_template) noexcept;

  /*! @brief Free the persistent basis and scratch. */
  void
  undefine() noexcept;

  /*!
    @brief Whether define() has been called.
    @return True if defined.
  */
  bool
  isDefined() const noexcept;

  /*!
    @brief Solve a_op(a_phi) = a_rhs (homogeneous operator) with restarted GMRES(m_restart).
    @param[in,out] a_phi Solution (in: initial guess; out: solution). Warm starts are honoured.
    @param[in]     a_rhs Right-hand side.
    @return True if converged (m_exitStatus == 1).
  */
  bool
  solve(Vec& a_phi, const Vec& a_rhs) noexcept;

protected:
  /*! @brief Operator/adapter (not owned). */
  Op* m_op = nullptr;

  /*! @brief Whether the basis/scratch are allocated. */
  bool m_isDefined = false;

  /*! @brief Restart length the persistent storage was sized for. */
  int m_allocRestart = 0;

  /*! @brief Arnoldi basis, size m_allocRestart + 1. */
  std::vector<Vec> m_basis;

  Vec m_w; /*!< A*K^{-1}*v_j scratch. */
  Vec m_z; /*!< K^{-1}*v_j and final preconditioned-update scratch. */
  Vec m_c; /*!< Accumulated unpreconditioned solution update (sum y_k v_k). */
  Vec m_r; /*!< Residual. */

  std::vector<Real> m_H;  /*!< Hessenberg matrix, column-major, leading dimension m_allocRestart + 1. */
  std::vector<Real> m_cs; /*!< Givens cosines. */
  std::vector<Real> m_sn; /*!< Givens sines. */
  std::vector<Real> m_g;  /*!< Rotated residual right-hand side. */
  std::vector<Real> m_y;  /*!< Least-squares solution. */
};

#include <CD_NamespaceFooter.H>

#include <CD_KrylovGMRESImplem.H>

#endif
