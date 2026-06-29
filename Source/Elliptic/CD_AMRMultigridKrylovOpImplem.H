/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_AMRMultigridKrylovOpImplem.H
  @brief  Implementation of CD_AMRMultigridKrylovOp.H
  @author Robert Marskar
*/

#ifndef CD_AMRMultigridKrylovOpImplem_H
#define CD_AMRMultigridKrylovOpImplem_H

// Chombo includes
#include <CH_Timer.H>
#include <LayoutIterator.H>

// Our includes
#include <CD_AMRMultigridKrylovOp.H>
#include <CD_ParallelOps.H>
#include <CD_NamespaceHeader.H>

template <class T>
AMRMultigridKrylovOp<T>::~AMRMultigridKrylovOp()
{
  if (m_precondInit) {
    this->clear(m_precondTmp);
    this->clear(m_precondDcor);
  }
}

template <class T>
void
AMRMultigridKrylovOp<T>::define(AMRMultiGrid<LevelData<T>>*      a_mg,
                                const Vector<DisjointBoxLayout>& a_grids,
                                const Vector<int>&               a_refRat,
                                const Vector<Real>&              a_dx,
                                const int                        a_lbase,
                                const int                        a_lmax,
                                const int                        a_numVCycles)
{
  CH_TIME("AMRMultigridKrylovOp::define");

  m_mg           = a_mg;
  m_amrOps       = a_mg->getAMROperators();
  m_lbase        = a_lbase;
  m_lmax         = a_lmax;
  m_numVCycles   = std::max(1, a_numVCycles);
  m_precondInit  = false;
  m_volumeCached = false; // Geometry changed; the cached per-level volumes are stale.

  // AMR inner-product weight dx^SpaceDim per level.
  m_dxScale.resize(a_lmax + 1, 1.0);
  for (int lvl = 0; lvl <= a_lmax; lvl++) {
    m_dxScale[lvl] = D_TERM(a_dx[lvl], *a_dx[lvl], *a_dx[lvl]);
  }

  // Precompute a per-level "valid" mask: true everywhere, then false on regions covered by a finer level. This is
  // the allocation-free analogue of MultilevelLinearOp, which instead zeroes covered data before each dot product.
  m_validMask.resize(a_lmax + 1);
  for (int lvl = a_lbase; lvl <= a_lmax; lvl++) {
    m_validMask[lvl] = RefCountedPtr<LevelData<BaseFab<bool>>>(
      new LevelData<BaseFab<bool>>(a_grids[lvl], 1, IntVect::Zero));

    const DataIterator& dit = a_grids[lvl].dataIterator();
    for (int mybox = 0; mybox < dit.size(); mybox++) {
      (*m_validMask[lvl])[dit[mybox]].setVal(true);
    }

    if (lvl < a_lmax) {
      const DisjointBoxLayout& finerGrids = a_grids[lvl + 1];
      const int                nRefFine   = a_refRat[lvl];

      for (int mybox = 0; mybox < dit.size(); mybox++) {
        const DataIndex& din     = dit[mybox];
        const Box        thisBox = a_grids[lvl][din];
        BaseFab<bool>&   mask    = (*m_validMask[lvl])[din];

        for (LayoutIterator lit = finerGrids.layoutIterator(); lit.ok(); ++lit) {
          Box coarFine = finerGrids[lit];
          coarFine.coarsen(nRefFine);
          coarFine &= thisBox;
          if (!coarFine.isEmpty()) {
            mask.setVal(false, coarFine, 0, 1);
          }
        }
      }
    }
  }
}

template <class T>
void
AMRMultigridKrylovOp<T>::applyOp(Vec& a_lhs, const Vec& a_phi, bool a_homogeneous)
{
  m_mg->computeAMROperator(a_lhs, const_cast<Vec&>(a_phi), m_lmax, m_lbase, a_homogeneous);
}

template <class T>
void
AMRMultigridKrylovOp<T>::residual(Vec& a_lhs, const Vec& a_phi, const Vec& a_rhs, bool a_homogeneous)
{
  m_mg->computeAMRResidual(a_lhs, const_cast<Vec&>(a_phi), a_rhs, m_lmax, m_lbase, a_homogeneous, false);
}

template <class T>
void
AMRMultigridKrylovOp<T>::preCond(Vec& a_cor, const Vec& a_residual)
{
  CH_TIME("AMRMultigridKrylovOp::preCond");

  if (!m_precondInit) {
    this->create(m_precondTmp, a_residual);
    this->create(m_precondDcor, a_residual);
    m_precondInit = true;
  }

  this->setToZero(a_cor);

  for (int iter = 0; iter < m_numVCycles; iter++) {
    // The V-cycle needs (residual - L(cor)) as its input and mangles that argument, so we feed it a scratch copy.
    // On the first cycle cor == 0, so the input is just a_residual -- skip the (expensive) residual evaluation. Only
    // multi-cycle preconditioning (numVCycles > 1) needs the recompute, which keeps the common single-cycle case as
    // cheap as one V-cycle (one extra matvec per cycle would otherwise roughly double the work per Krylov iteration).
    if (iter == 0) {
      this->assign(m_precondTmp, a_residual);
    }
    else {
      this->residual(m_precondTmp, a_cor, a_residual, true);
      this->scale(m_precondTmp, -1.0);
    }

    this->setToZero(m_precondDcor);
    m_mg->AMRVCycle(m_precondDcor, m_precondTmp, m_lmax, m_lmax, m_lbase);

    this->incr(a_cor, m_precondDcor, 1.0);
  }
}

template <class T>
void
AMRMultigridKrylovOp<T>::create(Vec& a_lhs, const Vec& a_rhs)
{
  a_lhs.resize(a_rhs.size(), nullptr);
  for (int lvl = m_lbase; lvl < a_rhs.size(); lvl++) {
    a_lhs[lvl] = new LevelData<T>();
    m_amrOps[lvl]->create(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

template <class T>
void
AMRMultigridKrylovOp<T>::clear(Vec& a_lhs)
{
  for (int lvl = 0; lvl < a_lhs.size(); lvl++) {
    if (a_lhs[lvl] != nullptr) {
      delete a_lhs[lvl];
      a_lhs[lvl] = nullptr;
    }
  }
}

template <class T>
void
AMRMultigridKrylovOp<T>::assign(Vec& a_lhs, const Vec& a_rhs)
{
  for (int lvl = m_lbase; lvl < a_rhs.size(); lvl++) {
    m_amrOps[lvl]->assign(*a_lhs[lvl], *a_rhs[lvl]);
  }
}

template <class T>
Real
AMRMultigridKrylovOp<T>::dotProduct(const Vec& a_1, const Vec& a_2)
{
  CH_TIME("AMRMultigridKrylovOp::dotProduct");

  using OpType = typename AMRMultigridKrylovOpTraits<T>::OpType;

  const int  nLevels = m_lmax - m_lbase + 1;
  const bool needVol = !m_volumeCached;

  // Rank-local partials packed into one buffer: per-level numerators, then (first call only) per-level volumes.
  // A single MPI_Allreduce of this buffer replaces the 2*nLevels separate reductions of the old per-level path.
  Vector<Real> partials(needVol ? 2 * nLevels : nLevels, 0.0);

  for (int lvl = m_lbase; lvl <= m_lmax; lvl++) {
    auto* op  = static_cast<OpType*>(m_amrOps[lvl]);
    Real  xy  = 0.0;
    Real  vol = 0.0;
    op->dotProductMaskedLocal(xy, vol, *a_1[lvl], *a_2[lvl], *m_validMask[lvl], needVol);

    partials[lvl - m_lbase] = xy;
    if (needVol) {
      partials[nLevels + (lvl - m_lbase)] = vol;
    }
  }

  ParallelOps::sum(partials);

  if (needVol) {
    m_volume.resize(m_lmax + 1, 0.0);
    for (int lvl = m_lbase; lvl <= m_lmax; lvl++) {
      m_volume[lvl] = partials[nLevels + (lvl - m_lbase)];
    }
    m_volumeCached = true;
  }

  Real prod = 0.0;
  for (int lvl = m_lbase; lvl <= m_lmax; lvl++) {
    if (m_volume[lvl] > 0.0) {
      prod += (partials[lvl - m_lbase] / m_volume[lvl]) * m_dxScale[lvl];
    }
  }

  return prod;
}

template <class T>
void
AMRMultigridKrylovOp<T>::mDotProduct(const Vec& a_1, const int a_sz, const Vec a_2[], Real a_mdots[])
{
  CH_TIME("AMRMultigridKrylovOp::mDotProduct");

  using OpType = typename AMRMultigridKrylovOpTraits<T>::OpType;

  const int  nLevels = m_lmax - m_lbase + 1;
  const int  nNum    = a_sz * nLevels;
  const bool needVol = !m_volumeCached;

  // All a_sz*nLevels rank-local numerators (plus per-level volumes on the first call) in one buffer -> one reduction.
  Vector<Real> partials(needVol ? nNum + nLevels : nNum, 0.0);

  for (int lvl = m_lbase; lvl <= m_lmax; lvl++) {
    auto* op  = static_cast<OpType*>(m_amrOps[lvl]);
    Real  vol = 0.0;

    for (int k = 0; k < a_sz; k++) {
      Real       xy          = 0.0;
      Real       v           = 0.0;
      const bool needVolHere = needVol && (k == 0); // Volume is geometry-only; compute it once per level.

      op->dotProductMaskedLocal(xy, v, *a_1[lvl], *(a_2[k][lvl]), *m_validMask[lvl], needVolHere);

      partials[k * nLevels + (lvl - m_lbase)] = xy;
      if (needVolHere) {
        vol = v;
      }
    }

    if (needVol) {
      partials[nNum + (lvl - m_lbase)] = vol;
    }
  }

  ParallelOps::sum(partials);

  if (needVol) {
    m_volume.resize(m_lmax + 1, 0.0);
    for (int lvl = m_lbase; lvl <= m_lmax; lvl++) {
      m_volume[lvl] = partials[nNum + (lvl - m_lbase)];
    }
    m_volumeCached = true;
  }

  for (int k = 0; k < a_sz; k++) {
    Real prod = 0.0;
    for (int lvl = m_lbase; lvl <= m_lmax; lvl++) {
      if (m_volume[lvl] > 0.0) {
        prod += (partials[k * nLevels + (lvl - m_lbase)] / m_volume[lvl]) * m_dxScale[lvl];
      }
    }
    a_mdots[k] = prod;
  }
}

template <class T>
void
AMRMultigridKrylovOp<T>::incr(Vec& a_lhs, const Vec& a_x, Real a_scale)
{
  for (int lvl = m_lbase; lvl < a_lhs.size(); lvl++) {
    m_amrOps[lvl]->incr(*a_lhs[lvl], *a_x[lvl], a_scale);
  }
}

template <class T>
void
AMRMultigridKrylovOp<T>::axby(Vec& a_lhs, const Vec& a_x, const Vec& a_y, Real a_a, Real a_b)
{
  for (int lvl = m_lbase; lvl < a_lhs.size(); lvl++) {
    m_amrOps[lvl]->axby(*a_lhs[lvl], *a_x[lvl], *a_y[lvl], a_a, a_b);
  }
}

template <class T>
void
AMRMultigridKrylovOp<T>::scale(Vec& a_lhs, const Real& a_scale)
{
  for (int lvl = m_lbase; lvl < a_lhs.size(); lvl++) {
    m_amrOps[lvl]->scale(*a_lhs[lvl], a_scale);
  }
}

template <class T>
Real
AMRMultigridKrylovOp<T>::norm(const Vec& a_rhs, int a_ord)
{
  if (a_ord == 2) {
    return std::sqrt(this->dotProduct(a_rhs, a_rhs));
  }

  Real maxNorm = 0.0;
  for (int lvl = m_lbase; lvl <= m_lmax; lvl++) {
    maxNorm = std::max(maxNorm, m_amrOps[lvl]->norm(*a_rhs[lvl], a_ord));
  }

  return maxNorm;
}

template <class T>
void
AMRMultigridKrylovOp<T>::setToZero(Vec& a_lhs)
{
  for (int lvl = m_lbase; lvl < a_lhs.size(); lvl++) {
    m_amrOps[lvl]->setToZero(*a_lhs[lvl]);
  }
}

#include <CD_NamespaceFooter.H>

#endif
