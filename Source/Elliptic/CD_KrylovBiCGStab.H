/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_KrylovBiCGStab.H
  @brief  Preconditioned BiCGStab over an AMR hierarchy, with persistent work vectors and fused reductions
  @author Robert Marskar
*/

#ifndef CD_KrylovBiCGStab_H
#define CD_KrylovBiCGStab_H

// Chombo includes
#include <LevelData.H>
#include <Vector.H>

// Our includes
#include <CD_AMRMultigridKrylovOp.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Preconditioned BiCGStab solver acting on an AMR hierarchy through AMRMultigridKrylovOp.
  @details A chombo-discharge replacement for Chombo's BiCGStabSolver on the outer Krylov path, written for
           scale-out: the eight work vectors are allocated once in define() and reused across solves (no per-solve
           allocation), and the omega numerator/denominator inner products are fused into a single mDotProduct so the
           adapter issues one MPI reduction instead of two. The matrix-vector product and preconditioner are the
           homogeneous AMR operator and V-cycle supplied by the adapter; solve() updates the solution in place and
           honours a non-zero initial guess (warm start).
  @note T is the FAB type (EBCellFAB or MFCellFAB).
*/
template <class T>
class KrylovBiCGStab
{
public:
  using Op  = AMRMultigridKrylovOp<T>; /*!< Operator/adapter type. */
  using Vec = Vector<LevelData<T>*>;   /*!< Multilevel data type the solver acts on. */

  /*! @brief Relative residual tolerance (against the initial residual norm). */
  Real m_eps = 1.E-10;

  /*! @brief Maximum number of iterations. */
  int m_maxIter = 50;

  /*! @brief Maximum number of restarts on breakdown. */
  int m_numRestarts = 5;

  /*! @brief Verbosity; at >= 4 a terse one-line-per-iteration residual/rate history is printed. */
  int m_verbosity = -1;

  /*! @brief Exit status after solve(): 1 = converged, 2 = breakdown, 3 = max iterations, -1 = not run. */
  int m_exitStatus = -1;

  /*! @brief Final residual norm after solve() (in the adapter inner-product norm). */
  Real m_residualNorm = -1.0;

  /*! @brief Default constructor. Call define() before use. */
  KrylovBiCGStab() noexcept = default;

  /*! @brief Destructor. Frees the persistent work vectors. */
  ~KrylovBiCGStab() noexcept;

  /*!
    @brief Allocate the persistent work vectors. Call once per regrid.
    @param[in] a_op       The (already define()d) operator/adapter to solve against
    @param[in] a_template A hierarchy whose layout the work vectors mirror (e.g. the rhs)
  */
  void
  define(Op* a_op, const Vec& a_template) noexcept;

  /*! @brief Free the persistent work vectors. */
  void
  undefine() noexcept;

  /*!
    @brief Whether define() has been called and the work vectors are allocated.
    @return True if defined.
  */
  bool
  isDefined() const noexcept;

  /*!
    @brief Solve a_op(a_phi) = a_rhs (homogeneous operator) with preconditioned BiCGStab.
    @param[in,out] a_phi Solution (in: initial guess; out: solution). Warm starts are honoured.
    @param[in]     a_rhs Right-hand side.
    @return True if converged (m_exitStatus == 1).
  */
  bool
  solve(Vec& a_phi, const Vec& a_rhs) noexcept;

protected:
  /*! @brief Operator/adapter (not owned). */
  Op* m_op = nullptr;

  /*! @brief Whether the work vectors are allocated. */
  bool m_isDefined = false;

  Vec m_r;      /*!< Residual / s. */
  Vec m_rTilde; /*!< Shadow residual. */
  Vec m_e;      /*!< Accumulated correction (added to a_phi at the end). */
  Vec m_p;      /*!< Search direction. */
  Vec m_pTilde; /*!< Preconditioned search direction. */
  Vec m_sTilde; /*!< Preconditioned s. */
  Vec m_t;      /*!< A*s_tilde. */
  Vec m_v;      /*!< A*p_tilde. */
};

#include <CD_NamespaceFooter.H>

#include <CD_KrylovBiCGStabImplem.H>

#endif
