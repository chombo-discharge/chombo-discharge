/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_MirrorDepositionImplem.H
 * @brief  Implementation of CD_MirrorDeposition.H
 * @author Robert Marskar
 */

#ifndef CD_MIRRORDEPOSITIONIMPLEM_H
#define CD_MIRRORDEPOSITIONIMPLEM_H

// Std includes
#include <cmath>

// Our includes
#include <CD_MirrorDeposition.H>
#include <CD_NamespaceHeader.H>

namespace MirrorDeposition {

inline int
shapeIndex(const int a_i, const int a_j) noexcept
{
  // Upper triangle in row-major order: (0,0),(0,1),...,(0,D-1),(1,1),... The row offset is the number of entries
  // in the preceding rows, i.e. sum_{k<i} (D - k).
  const int lo = (a_i < a_j) ? a_i : a_j;
  const int hi = (a_i < a_j) ? a_j : a_i;

  return lo * SpaceDim - (lo * (lo - 1)) / 2 + (hi - lo);
}

inline void
liftShapeOperator(const RealVect (&a_tangents)[numTangents],
                  const Real (&a_S)[numTangents][numTangents],
                  Real (&a_shape)[numShapeComp]) noexcept
{
  // S_c = P S P^T with P = [t_1 ... t_(D-1)]. Only the upper triangle is stored; the tensor is symmetric because
  // the fit symmetrizes S before it gets here.
  for (int i = 0; i < SpaceDim; i++) {
    for (int j = i; j < SpaceDim; j++) {
      Real sum = 0.0;

      for (int a = 0; a < numTangents; a++) {
        for (int b = 0; b < numTangents; b++) {
          sum += a_tangents[a][i] * a_S[a][b] * a_tangents[b][j];
        }
      }

      a_shape[shapeIndex(i, j)] = sum;
    }
  }
}

inline Real
shapeEntry(const Real* a_shape, const int a_i, const int a_j) noexcept
{
  return a_shape[shapeIndex(a_i, a_j)];
}

inline RealVect
applyShapeOperator(const Real* a_shape, const RealVect& a_vec) noexcept
{
  RealVect ret = RealVect::Zero;

  for (int i = 0; i < SpaceDim; i++) {
    Real sum = 0.0;

    for (int j = 0; j < SpaceDim; j++) {
      sum += shapeEntry(a_shape, i, j) * a_vec[j];
    }

    ret[i] = sum;
  }

  return ret;
}

inline Real
meanCurvatureTimesTwo(const Real* a_shape) noexcept
{
  Real tr = 0.0;

  for (int i = 0; i < SpaceDim; i++) {
    tr += shapeEntry(a_shape, i, i);
  }

  return tr;
}

inline Real
gaussianCurvature(const Real* a_shape) noexcept
{
  // K = (1/2)[(tr S)^2 - tr(S^2)]. See the header for why this is NOT det(S).
  const Real tr = meanCurvatureTimesTwo(a_shape);

  Real trSq = 0.0;

  for (int i = 0; i < SpaceDim; i++) {
    for (int j = 0; j < SpaceDim; j++) {
      const Real sij = shapeEntry(a_shape, i, j);

      trSq += sij * sij;
    }
  }

  return 0.5 * (tr * tr - trSq);
}

inline Real
jacobian(const Real a_twoH, const Real a_K, const Real a_d, Real& a_denominator) noexcept
{
  const Real dSq = a_d * a_d;

  a_denominator = 1.0 + a_twoH * a_d + a_K * dSq;

  return (1.0 - a_twoH * a_d + a_K * dSq) / a_denominator;
}

inline bool
reflect(const RealVect& a_pos,
        const Real*     a_surfaceData,
        const Real      a_maxJacobian,
        const Real      a_minDenominator,
        RealVect&       a_image,
        Real&           a_jacobian,
        Refusal&        a_refusal) noexcept
{
  RealVect centroid;
  RealVect normal;

  for (int dir = 0; dir < SpaceDim; dir++) {
    centroid[dir] = a_surfaceData[compCentroid + dir];
    normal[dir]   = a_surfaceData[compNormal + dir];
  }

  const Real* shape = a_surfaceData + compShape;

  // The local quadratic patch. Note that d is the first iterate of the signed distance, not the signed distance.
  const RealVect w   = a_pos - centroid;
  const RealVect Sw  = applyShapeOperator(shape, w);
  const Real     eta = normal.dotProduct(w);
  const Real     d   = eta + 0.5 * w.dotProduct(Sw);

  // The patch normal at the foot point. S_c annihilates n_c, so the normal component of w contributes nothing here
  // and no tangential projection of w is needed.
  RealVect   nhat = normal + Sw;
  const Real len  = nhat.vectorLength();

  if (len > 0.0) {
    nhat /= len;
  }
  else {
    nhat = normal;
  }

  a_image = a_pos - 2.0 * d * nhat;

  Real       denominator = 0.0;
  const Real J           = jacobian(meanCurvatureTimesTwo(shape), gaussianCurvature(shape), d, denominator);

  a_refusal = Refusal::None;

  if (std::abs(denominator) < a_minDenominator) {
    a_refusal = Refusal::SmallDenominator;
  }
  else if (!(J > 0.0)) {
    a_refusal = Refusal::NonPositiveJacobian;
  }
  else if (J > a_maxJacobian) {
    a_refusal = Refusal::LargeJacobian;
  }

  // A refused image still deposits, with unit weight -- dropping it would remove the correction precisely where it
  // is largest. See the header.
  a_jacobian = (a_refusal == Refusal::None) ? J : 1.0;

  return a_refusal == Refusal::None;
}

} // namespace MirrorDeposition

#include <CD_NamespaceFooter.H>

#endif
