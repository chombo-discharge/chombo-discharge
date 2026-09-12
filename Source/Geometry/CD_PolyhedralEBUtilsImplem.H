/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PolyhedralEBUtilsImplem.H
 * @brief  Implementation of CD_PolyhedralEBUtils.H
 * @author Robert Marskar
 */

#ifndef CD_POLYHEDRALEBUTILSIMPLEM_H
#define CD_POLYHEDRALEBUTILSIMPLEM_H

// Std includes
#include <cmath>

// Our includes
#include <CD_PolyhedralEBUtils.H>
#include <CD_NamespaceHeader.H>

namespace PolyhedralEB {

inline bool
isFluid(const Real a_value) noexcept
{
  return std::copysign(1.0, a_value) < 0.0;
}

inline RealVect
gradientNormal(const CutCellSurface& a_surface) noexcept
{
  RealVect gradient = RealVect::Zero;

  for (int c = 0; c < CutCellSurface::s_numCorners; c++) {
    for (int d = 0; d < SpaceDim; d++) {
      if (((c >> d) & 1) == 0) {
        gradient[d] += a_surface.m_corner[c | (1 << d)] - a_surface.m_corner[c];
      }
    }
  }

  const Real length = gradient.vectorLength();

  if (length <= 0.0) {
    return RealVect::Zero;
  }

  // The fluid is where the implicit function is negative, so the normal runs against its gradient.
  return -gradient / length;
}

inline RealVect
crossingNormal(const CutCellSurface& a_surface) noexcept
{
  RealVect point[CutCellSurface::s_numEdges];

  int numPoints = 0;

  for (int e = 0; e < CutCellSurface::s_numEdges; e++) {
    if (a_surface.hasCrossing(e)) {
      point[numPoints++] = detail::crossingPosition(a_surface, e, 0.0);
    }
  }

  if (numPoints < SpaceDim) {
    return RealVect::Zero;
  }

  RealVect centre = RealVect::Zero;

  for (int i = 0; i < numPoints; i++) {
    centre += point[i];
  }

  centre /= static_cast<Real>(numPoints);

  Real covariance[SpaceDim][SpaceDim] = {{0.0}};

  for (int i = 0; i < numPoints; i++) {
    const RealVect offset = point[i] - centre;

    for (int d = 0; d < SpaceDim; d++) {
      for (int e = 0; e < SpaceDim; e++) {
        covariance[d][e] += offset[d] * offset[e];
      }
    }
  }

  // The adjugate of a singular matrix has its columns in the null space, so for points that lie
  // in a plane it is built from the plane's normal. The longest column is the best conditioned.
  Real adjugate[SpaceDim][SpaceDim];

#if CH_SPACEDIM == 3
  const Real a = covariance[0][0];
  const Real b = covariance[0][1];
  const Real c = covariance[0][2];
  const Real d = covariance[1][1];
  const Real e = covariance[1][2];
  const Real f = covariance[2][2];

  adjugate[0][0] = d * f - e * e;
  adjugate[0][1] = c * e - b * f;
  adjugate[0][2] = b * e - c * d;
  adjugate[1][0] = adjugate[0][1];
  adjugate[1][1] = a * f - c * c;
  adjugate[1][2] = b * c - a * e;
  adjugate[2][0] = adjugate[0][2];
  adjugate[2][1] = adjugate[1][2];
  adjugate[2][2] = a * d - b * b;
#else
  adjugate[0][0] = covariance[1][1];
  adjugate[0][1] = -covariance[0][1];
  adjugate[1][0] = -covariance[0][1];
  adjugate[1][1] = covariance[0][0];
#endif

  RealVect normal      = RealVect::Zero;
  Real     bestLength2 = 0.0;

  for (int col = 0; col < SpaceDim; col++) {
    RealVect candidate;

    for (int row = 0; row < SpaceDim; row++) {
      candidate[row] = adjugate[row][col];
    }

    const Real length2 = candidate.dotProduct(candidate);

    if (length2 > bestLength2) {
      bestLength2 = length2;
      normal      = candidate;
    }
  }

  if (bestLength2 <= 0.0) {
    return RealVect::Zero;
  }

  normal /= std::sqrt(bestLength2);

  // The fit fixes the plane but not which way round it faces; the corner values do that.
  if (normal.dotProduct(gradientNormal(a_surface)) < 0.0) {
    normal = -normal;
  }

  return normal;
}

inline Real
maxPrincipalCurvature(const Real a_stencil[s_curvatureStencilSize], const Real a_spacing) noexcept
{
  // Offsets are packed base three, so stepping one lattice point in direction d changes the
  // index by 3^d.
  int stride[SpaceDim];
  int centre = 0;

  for (int d = 0; d < SpaceDim; d++) {
    stride[d] = (d == 0) ? 1 : stride[d - 1] * 3;
    centre += stride[d];
  }

  const Real f0 = a_stencil[centre];

  Real gradient[SpaceDim];
  Real hessian[SpaceDim][SpaceDim];

  for (int d = 0; d < SpaceDim; d++) {
    const Real fp = a_stencil[centre + stride[d]];
    const Real fm = a_stencil[centre - stride[d]];

    gradient[d]   = 0.5 * (fp - fm) / a_spacing;
    hessian[d][d] = (fp - 2.0 * f0 + fm) / (a_spacing * a_spacing);
  }

  for (int d = 0; d < SpaceDim; d++) {
    for (int e = d + 1; e < SpaceDim; e++) {
      const Real fpp = a_stencil[centre + stride[d] + stride[e]];
      const Real fpm = a_stencil[centre + stride[d] - stride[e]];
      const Real fmp = a_stencil[centre - stride[d] + stride[e]];
      const Real fmm = a_stencil[centre - stride[d] - stride[e]];

      hessian[d][e] = 0.25 * (fpp - fpm - fmp + fmm) / (a_spacing * a_spacing);
      hessian[e][d] = hessian[d][e];
    }
  }

  Real gradientSquared = 0.0;

  for (int d = 0; d < SpaceDim; d++) {
    gradientSquared += gradient[d] * gradient[d];
  }

  // A vanishing gradient is a point the level set does not pass through cleanly, and the
  // curvature there is not defined. Reporting zero leaves such a cell unrefined by this
  // criterion rather than refined on a divided-by-nothing value.
  if (gradientSquared <= 0.0) {
    return 0.0;
  }

  const Real gradientLength = std::sqrt(gradientSquared);

  Real trace = 0.0;

  for (int d = 0; d < SpaceDim; d++) {
    trace += hessian[d][d];
  }

  Real quadratic = 0.0;

  for (int d = 0; d < SpaceDim; d++) {
    for (int e = 0; e < SpaceDim; e++) {
      quadratic += gradient[d] * hessian[d][e] * gradient[e];
    }
  }

  // Sum of the principal curvatures: the divergence of the unit normal.
  const Real curvatureSum = (gradientSquared * trace - quadratic) / (gradientSquared * gradientLength);

#if CH_SPACEDIM == 3
  // The Gaussian curvature is the same quadratic form taken against the adjugate of the Hessian,
  // which for a symmetric matrix is its cofactor matrix. With the sum and the product of the two
  // principal curvatures in hand they follow from the quadratic they are the roots of.
  const Real a = hessian[0][0];
  const Real b = hessian[0][1];
  const Real c = hessian[0][2];
  const Real d = hessian[1][1];
  const Real e = hessian[1][2];
  const Real f = hessian[2][2];

  Real adjugate[3][3];

  adjugate[0][0] = d * f - e * e;
  adjugate[0][1] = c * e - b * f;
  adjugate[0][2] = b * e - c * d;
  adjugate[1][0] = adjugate[0][1];
  adjugate[1][1] = a * f - c * c;
  adjugate[1][2] = b * c - a * e;
  adjugate[2][0] = adjugate[0][2];
  adjugate[2][1] = adjugate[1][2];
  adjugate[2][2] = a * d - b * b;

  Real adjugateForm = 0.0;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      adjugateForm += gradient[i] * adjugate[i][j] * gradient[j];
    }
  }

  const Real gaussian = adjugateForm / (gradientSquared * gradientSquared);
  const Real mean     = 0.5 * curvatureSum;

  // Equal principal curvatures make the discriminant zero, and rounding can carry it just below.
  const Real discriminant = std::max(mean * mean - gaussian, 0.0);

  return std::abs(mean) + std::sqrt(discriminant);
#else
  return std::abs(curvatureSum);
#endif
}

namespace detail {

inline int
edgeDirection(const int a_edge) noexcept
{
  return a_edge / (CutCellSurface::s_numEdges / SpaceDim);
}

inline void
edgeOrigin(const int a_edge, int a_offset[SpaceDim]) noexcept
{
  const int dir   = edgeDirection(a_edge);
  const int local = a_edge % (CutCellSurface::s_numEdges / SpaceDim);

  for (int d = 0; d < SpaceDim; d++) {
    a_offset[d] = 0;
  }

#if CH_SPACEDIM == 3
  a_offset[s_transverse[dir][0]] = local & 1;
  a_offset[s_transverse[dir][1]] = (local >> 1) & 1;
#else
  a_offset[1 - dir] = local & 1;
#endif

  a_offset[dir] = 0;
}

inline void
edgeCorners(const int a_edge, int& a_lo, int& a_hi) noexcept
{
  const int dir = edgeDirection(a_edge);

  int offset[SpaceDim];
  edgeOrigin(a_edge, offset);

  a_lo = 0;

  for (int d = 0; d < SpaceDim; d++) {
    a_lo |= offset[d] << d;
  }

  a_hi = a_lo | (1 << dir);
}

inline RealVect
cornerPosition(const int a_corner) noexcept
{
  RealVect x = RealVect::Zero;

  for (int d = 0; d < SpaceDim; d++) {
    x[d] = -0.5 + static_cast<Real>((a_corner >> d) & 1);
  }

  return x;
}

inline void
faceCorners(const int a_dir, const int a_side, int a_corner[1 << (SpaceDim - 1)]) noexcept
{
#if CH_SPACEDIM == 3
  const int t0 = s_transverse[a_dir][0];
  const int t1 = s_transverse[a_dir][1];

  const int ring[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};

  for (int i = 0; i < 4; i++) {
    a_corner[i] = (a_side << a_dir) | (ring[i][0] << t0) | (ring[i][1] << t1);
  }
#else
  const int t = 1 - a_dir;

  for (int i = 0; i < 2; i++) {
    a_corner[i] = (a_side << a_dir) | (i << t);
  }
#endif
}

inline int
edgeIndex(const int a_dir, const int a_offset[SpaceDim]) noexcept
{
#if CH_SPACEDIM == 3
  return 4 * a_dir + ((a_offset[s_transverse[a_dir][0]] & 1) | ((a_offset[s_transverse[a_dir][1]] & 1) << 1));
#else
  return 2 * a_dir + (a_offset[1 - a_dir] & 1);
#endif
}

inline void
faceEdges(const int a_dir, const int a_side, int a_edge[2 * (SpaceDim - 1)]) noexcept
{
#if CH_SPACEDIM == 3
  const int t0 = s_transverse[a_dir][0];
  const int t1 = s_transverse[a_dir][1];

  // each entry is the t0 and t1 offset of the edge's low corner, and the direction it runs
  const int spec[4][3] = {{0, 0, t0}, {1, 0, t1}, {1, 1, t0}, {0, 1, t1}};

  for (int i = 0; i < 4; i++) {
    int offset[SpaceDim];

    for (int d = 0; d < SpaceDim; d++) {
      offset[d] = 0;
    }

    offset[a_dir] = a_side;
    offset[t0]    = spec[i][0];
    offset[t1]    = spec[i][1];

    const int run = spec[i][2];

    offset[run] = 0;

    a_edge[i] = edgeIndex(run, offset);
  }
#else
  (void)a_dir;
  (void)a_side;
  (void)a_edge;
#endif
}

inline RealVect
crossingPosition(const CutCellSurface& a_surface, const int a_edge, const Real a_tolerance) noexcept
{
  const int dir = edgeDirection(a_edge);

  int offset[SpaceDim];
  edgeOrigin(a_edge, offset);

  Real t = a_surface.m_crossing[a_edge];

  // A crossing recorded exactly at an endpoint sits on a corner the interface passes through,
  // so it is already where it belongs and displacing it would open a sliver of the
  // displacement's own width. Every other crossing is held off the endpoints, which is what
  // keeps the combinatorics generic.
  if (t != 0.0 && t != 1.0) {
    t = std::max(t, a_tolerance);
    t = std::min(t, 1.0 - a_tolerance);
  }

  RealVect x = RealVect::Zero;

  for (int d = 0; d < SpaceDim; d++) {
    x[d] = -0.5 + static_cast<Real>(offset[d]);
  }

  x[dir] = -0.5 + t;

  return x;
}

inline void
polygonMoments(const RealVect* a_vertex,
               const int       a_num,
               Real&           a_area,
               RealVect&       a_vector,
               RealVect&       a_centroid) noexcept
{
  a_area     = 0.0;
  a_vector   = RealVect::Zero;
  a_centroid = RealVect::Zero;

  if (a_num < 3) {
    return;
  }

#if CH_SPACEDIM == 3
  for (int i = 1; i < a_num - 1; i++) {
    const RealVect u = a_vertex[i] - a_vertex[0];
    const RealVect v = a_vertex[i + 1] - a_vertex[0];

    a_vector += 0.5 * RealVect(D_DECL(u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]));
  }

  a_area = a_vector.vectorLength();

  if (a_area <= 0.0) {
    return;
  }

  const RealVect unit = a_vector / a_area;

  for (int i = 1; i < a_num - 1; i++) {
    const RealVect u = a_vertex[i] - a_vertex[0];
    const RealVect v = a_vertex[i + 1] - a_vertex[0];

    const RealVect n = RealVect(
      D_DECL(u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]));

    const Real w = 0.5 * n.dotProduct(unit);

    a_centroid += w * (a_vertex[0] + a_vertex[i] + a_vertex[i + 1]) / 3.0;
  }

  a_centroid /= a_area;
#else
  (void)a_vertex;
#endif
}

#if CH_SPACEDIM == 3
inline int
facePairs(const int a_dir, const int a_side, const CutCellSurface& a_surface, int a_pair[2][2]) noexcept
{
  int faceEdge[4];
  int faceCorner[4];

  faceEdges(a_dir, a_side, faceEdge);
  faceCorners(a_dir, a_side, faceCorner);

  int hit[4];
  int numHit = 0;

  for (int i = 0; i < 4; i++) {
    if (a_surface.hasCrossing(faceEdge[i])) {
      hit[numHit++] = i;
    }
  }

  if (numHit == 0) {
    return 0;
  }

  if (numHit == 2) {
    a_pair[0][0] = faceEdge[hit[0]];
    a_pair[0][1] = faceEdge[hit[1]];

    return 1;
  }

  if (numHit != 4) {
    return -1;
  }

  const Real f00 = a_surface.m_corner[faceCorner[0]];
  const Real f10 = a_surface.m_corner[faceCorner[1]];
  const Real f11 = a_surface.m_corner[faceCorner[2]];
  const Real f01 = a_surface.m_corner[faceCorner[3]];

  const Real den    = f00 + f11 - f10 - f01;
  const Real saddle = (std::abs(den) > 0.0) ? (f00 * f11 - f10 * f01) / den : (f00 + f11);

  if (isFluid(saddle) == isFluid(f00)) {
    a_pair[0][0] = faceEdge[0]; // chords cut off face corners 1 and 3
    a_pair[0][1] = faceEdge[1];
    a_pair[1][0] = faceEdge[2];
    a_pair[1][1] = faceEdge[3];
  }
  else {
    a_pair[0][0] = faceEdge[3]; // chords cut off face corners 0 and 2
    a_pair[0][1] = faceEdge[0];
    a_pair[1][0] = faceEdge[1];
    a_pair[1][1] = faceEdge[2];
  }

  return 2;
}

inline int
crossingLoops(const CutCellSurface& a_surface,
              int                   a_loop[CutCellSurface::s_numEdges],
              int                   a_start[CutCellSurface::s_numEdges + 1]) noexcept
{
  constexpr int numEdges = CutCellSurface::s_numEdges;

  int adjacent[numEdges][2];
  int degree[numEdges] = {0};

  for (int d = 0; d < SpaceDim; d++) {
    for (int side = 0; side < 2; side++) {
      int       pair[2][2];
      const int numPairs = facePairs(d, side, a_surface, pair);

      if (numPairs < 0) {
        return -1;
      }

      for (int k = 0; k < numPairs; k++) {
        const int x = pair[k][0];
        const int y = pair[k][1];

        if (degree[x] > 1 || degree[y] > 1) {
          return -1;
        }

        adjacent[x][degree[x]++] = y;
        adjacent[y][degree[y]++] = x;
      }
    }
  }

  int numCrossings = 0;

  for (int e = 0; e < numEdges; e++) {
    if (a_surface.hasCrossing(e)) {
      numCrossings++;

      if (degree[e] != 2) {
        return -1;
      }
    }
  }

  if (numCrossings < 3) {
    return -1;
  }

  bool used[numEdges] = {false};
  int  numLoops       = 0;
  int  put            = 0;

  a_start[0] = 0;

  for (int e = 0; e < numEdges; e++) {
    if (!a_surface.hasCrossing(e) || used[e]) {
      continue;
    }

    const int begin = put;

    int previous = -1;
    int current  = e;

    while (true) {
      used[current] = true;
      a_loop[put++] = current;

      const int next = (adjacent[current][0] != previous) ? adjacent[current][0] : adjacent[current][1];

      if (next == e) {
        break;
      }

      if (used[next] || put > numEdges) {
        return -1;
      }

      previous = current;
      current  = next;
    }

    if (put - begin < 3) {
      return -1;
    }

    a_start[++numLoops] = put;
  }

  return numLoops;
}
#endif

inline VertexKey
vertexKey(const RealVect& a_x) noexcept
{
  VertexKey key = {};

  for (int d = 0; d < SpaceDim; d++) {
    key.m_q[d] = std::lround(a_x[d] * 1.0E11);
  }

  return key;
}

inline bool
lexLess(const RealVect& a_a, const RealVect& a_b) noexcept
{
  for (int d = 0; d < SpaceDim - 1; d++) {
    if (a_a[d] != a_b[d]) {
      return a_a[d] < a_b[d];
    }
  }

  return a_a[SpaceDim - 1] < a_b[SpaceDim - 1];
}

inline Real
polygonArea(const RealVect* a_vertex, const int a_num) noexcept
{
  if (a_num < 3) {
    return 0.0;
  }

#if CH_SPACEDIM == 3
  Real     unused = 0.0;
  RealVect vector;
  RealVect centroid;

  polygonMoments(a_vertex, a_num, unused, vector, centroid);

  return unused;
#else
  Real twice = 0.0;

  for (int i = 0; i < a_num; i++) {
    const RealVect& a = a_vertex[i];
    const RealVect& b = a_vertex[(i + 1) % a_num];

    twice += a[0] * b[1] - b[0] * a[1];
  }

  return 0.5 * std::abs(twice);
#endif
}

} // namespace detail
} // namespace PolyhedralEB

#include <CD_NamespaceFooter.H>

#endif
