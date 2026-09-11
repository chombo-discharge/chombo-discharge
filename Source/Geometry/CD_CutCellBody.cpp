/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_CutCellBody.cpp
 * @brief  Implementation of CD_CutCellBody.H
 * @author Robert Marskar
 */

// Std includes
#include <cmath>

// Our includes
#include <CD_CutCellBody.H>
#include <CD_NamespaceHeader.H>

namespace {

/**
 * @brief The two directions transverse to each coordinate direction.
 */
#if CH_SPACEDIM == 3
constexpr int s_transverse[3][2] = {{1, 2}, {0, 2}, {0, 1}};
#endif

/**
 * @brief Direction a cell edge runs along.
 */
inline int
edgeDirection(const int a_edge) noexcept
{
  return a_edge / (CutCellSurface::s_numEdges / SpaceDim);
}

/**
 * @brief Corner offsets of a cell edge's low end.
 */
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

/**
 * @brief The corners a cell edge joins, low end first.
 */
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

/**
 * @brief Position of a cell corner in the cell's own frame.
 */
inline RealVect
cornerPosition(const int a_corner) noexcept
{
  RealVect x = RealVect::Zero;

  for (int d = 0; d < SpaceDim; d++) {
    x[d] = -0.5 + static_cast<Real>((a_corner >> d) & 1);
  }

  return x;
}

/**
 * @brief The corners of a cell face, in circuit order around the face.
 */
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

/**
 * @brief Index of the edge running along a_dir whose low corner has the given offsets.
 */
inline int
edgeIndex(const int a_dir, const int a_offset[SpaceDim]) noexcept
{
#if CH_SPACEDIM == 3
  return 4 * a_dir + ((a_offset[s_transverse[a_dir][0]] & 1) | ((a_offset[s_transverse[a_dir][1]] & 1) << 1));
#else
  return 2 * a_dir + (a_offset[1 - a_dir] & 1);
#endif
}

/**
 * @brief The edges of a cell face, edge i joining face corners i and i+1.
 */
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

/**
 * @brief Position of an edge crossing in the cell's own frame.
 */
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

/**
 * @brief Area, area vector and centroid of a simple planar polygon in circuit order.
 * @details Areas are taken signed about the polygon's own normal. A face contour that is a
 * polyline rather than a single chord leaves a non-convex polygon, and summing unsigned
 * triangle areas over a fan over-counts those.
 */
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
/**
 * @brief The chords on a cell face, each an ordered pair of edge indices.
 * @details A face with two crossings carries one chord. A face whose corners alternate fluid
 * and solid carries four crossings and there are two ways to join them; the bilinear
 * interpolant of the four corner values decides it. The saddle shares its side with exactly
 * one of the two diagonals: that diagonal meets through the middle, and the other is the pair
 * the chords cut off. Both cells adjoining the face hold the same four values, so they cannot
 * disagree.
 * @return Number of chords, or -1 if the face carries an impossible number of crossings.
 */
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

  if (CutCellBody::isFluid(saddle) == CutCellBody::isFluid(f00)) {
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

/**
 * @brief Order the crossings into closed loops, or fail if they do not form clean cycles.
 * @details Every crossing lies on exactly two faces and every face pairs its crossings up, so
 * each node has degree two and the graph is a disjoint union of cycles. Several cycles is
 * legal: the interface simply enters the cell as several sheets.
 * @return Number of loops, or -1 if the crossings do not form clean cycles.
 */
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
} // namespace

CutCellBody::CutCellBody() noexcept
{
  m_numPolygons      = 0;
  m_volumeFraction   = 0.0;
  m_volumeCentroid   = RealVect::Zero;
  m_boundaryArea     = 0.0;
  m_trueBoundaryArea = 0.0;
  m_normal           = RealVect::Zero;
  m_boundaryCentroid = RealVect::Zero;
  m_closure          = RealVect::Zero;

  for (int f = 0; f < s_numFaces; f++) {
    m_areaFraction[f] = 0.0;
    m_faceCentroid[f] = RealVect::Zero;
  }
}

bool
CutCellBody::isFluid(const Real a_value) noexcept
{
  return std::copysign(1.0, a_value) < 0.0;
}

CutCellBody::Kind
CutCellBody::classify(const CutCellSurface& a_surface) noexcept
{
  const Real first     = a_surface.m_corner[0];
  const Real firstSign = std::copysign(1.0, first);

  bool cut = false;

  for (int c = 0; c < CutCellSurface::s_numCorners; c++) {
    const Real value = a_surface.m_corner[c];
    const Real sign  = std::copysign(1.0, value);

    if ((value == 0.0 || first == 0.0) && sign * firstSign < 0.0) {
      cut = true;
    }

    if (value * first < 0.0) {
      cut = true;
    }
  }

  // A facet lying in a node plane leaves every corner on one side at exactly zero. The signs
  // then disagree and the sign test alone reads the cell as cut, but a side represented only by
  // exact zeros encloses nothing: the cell is entirely on the other side, with the face in that
  // plane covered. Building it instead yields a sliver as wide as the crossings are held off
  // the edge endpoints, which is dust of a size no volume threshold is scaled to catch.
  if (cut) {
    if (touchesOnly(a_surface, true)) {
      return Kind::Covered;
    }

    if (touchesOnly(a_surface, false)) {
      return Kind::Regular;
    }

    return Kind::Cut;
  }

  return (firstSign < 0.0) ? Kind::Regular : Kind::Covered;
}

bool
CutCellBody::touchesOnly(const CutCellSurface& a_surface, const bool a_fluidSide) noexcept
{
  bool any = false;

  for (int c = 0; c < CutCellSurface::s_numCorners; c++) {
    const Real value = a_surface.m_corner[c];

    if (isFluid(value) == a_fluidSide) {
      any = true;

      if (value != 0.0) {
        return false;
      }
    }
  }

  return any;
}

void
CutCellBody::orientOutward(Polygon& a_polygon, const int a_dir, const int a_side) const noexcept
{
  RealVect outward = RealVect::Zero;
  outward[a_dir]   = (a_side == 0) ? -1.0 : 1.0;

  Real     area = 0.0;
  RealVect vector;
  RealVect centroid;

  polygonMoments(a_polygon.m_vertex, a_polygon.m_numVertices, area, vector, centroid);

  if (vector.dotProduct(outward) < 0.0) {
    for (int i = 0; i < a_polygon.m_numVertices / 2; i++) {
      const int j = a_polygon.m_numVertices - 1 - i;

      std::swap(a_polygon.m_vertex[i], a_polygon.m_vertex[j]);
      std::swap(a_polygon.m_vertexEdge[i], a_polygon.m_vertexEdge[j]);
    }
  }
}

int
CutCellBody::faceWalk(const int a_dir, const int a_side, const CutCellSurface& a_surface, Polygon* a_out) const noexcept
{
#if CH_SPACEDIM == 3
  int faceEdge[4];
  int faceCorner[4];

  faceEdges(a_dir, a_side, faceEdge);
  faceCorners(a_dir, a_side, faceCorner);

  int numHit = 0;

  for (int i = 0; i < 4; i++) {
    numHit += a_surface.hasCrossing(faceEdge[i]) ? 1 : 0;
  }

  int       pair[2][2] = {{-1, -1}, {-1, -1}};
  const int numPairs   = facePairs(a_dir, a_side, a_surface, pair);

  if (numPairs < 0) {
    return -1;
  }

  if (numHit != 4) {
    Polygon polygon;
    polygon.m_numVertices = 0;
    polygon.m_face        = 2 * a_dir + a_side;

    for (int i = 0; i < 4; i++) {
      if (isFluid(a_surface.m_corner[faceCorner[i]])) {
        polygon.m_vertexEdge[polygon.m_numVertices] = -1;
        polygon.m_vertex[polygon.m_numVertices++]   = cornerPosition(faceCorner[i]);
      }

      if (a_surface.hasCrossing(faceEdge[i])) {
        polygon.m_vertexEdge[polygon.m_numVertices] = faceEdge[i];
        polygon.m_vertex[polygon.m_numVertices++]   = crossingPosition(a_surface, faceEdge[i], s_edgeTolerance);
      }
    }

    if (polygon.m_numVertices < 3) {
      return 0;
    }

    this->orientOutward(polygon, a_dir, a_side);

    Real     area = 0.0;
    RealVect vector;
    RealVect centroid;

    polygonMoments(polygon.m_vertex, polygon.m_numVertices, area, vector, centroid);

    if (area <= 0.0) {
      return 0;
    }

    a_out[0] = polygon;

    return 1;
  }

  // Four crossings: the chords cut off one diagonal pair of corners, and which pair that is is
  // the decider's whole output. Reading the fluid region off the corner signs instead lets the
  // face disagree with the loop assembly, which trusts the pairing.
  bool isolated[4] = {false, false, false, false};

  for (int k = 0; k < 2; k++) {
    for (int i = 0; i < 4; i++) {
      const int previous = faceEdge[(i + 3) % 4];

      if ((pair[k][0] == previous && pair[k][1] == faceEdge[i]) ||
          (pair[k][0] == faceEdge[i] && pair[k][1] == previous)) {
        isolated[i] = true;
      }
    }
  }

  int first       = -1;
  int numIsolated = 0;

  for (int i = 0; i < 4; i++) {
    if (isolated[i]) {
      numIsolated++;

      if (first < 0) {
        first = i;
      }
    }
  }

  if (numIsolated != 2) {
    return -1;
  }

  if (!isFluid(a_surface.m_corner[faceCorner[first]])) {
    // the cut-off corners are solid, so the fluid is one polygon wrapping around both
    Polygon polygon;
    polygon.m_numVertices = 0;
    polygon.m_face        = 2 * a_dir + a_side;

    for (int i = 0; i < 4; i++) {
      if (isFluid(a_surface.m_corner[faceCorner[i]])) {
        polygon.m_vertexEdge[polygon.m_numVertices] = -1;
        polygon.m_vertex[polygon.m_numVertices++]   = cornerPosition(faceCorner[i]);
      }

      polygon.m_vertexEdge[polygon.m_numVertices] = faceEdge[i];
      polygon.m_vertex[polygon.m_numVertices++]   = crossingPosition(a_surface, faceEdge[i], s_edgeTolerance);
    }

    this->orientOutward(polygon, a_dir, a_side);

    Real     area = 0.0;
    RealVect vector;
    RealVect centroid;

    polygonMoments(polygon.m_vertex, polygon.m_numVertices, area, vector, centroid);

    if (area <= 0.0) {
      return 0;
    }

    a_out[0] = polygon;

    return 1;
  }

  // the cut-off corners are the fluid ones, so the fluid is two corner triangles
  int numOut = 0;

  for (int i = 0; i < 4; i++) {
    if (!isolated[i]) {
      continue;
    }

    const int previous = faceEdge[(i + 3) % 4];

    Polygon polygon;
    polygon.m_numVertices = 0;
    polygon.m_face        = 2 * a_dir + a_side;

    polygon.m_vertexEdge[polygon.m_numVertices] = previous;
    polygon.m_vertex[polygon.m_numVertices++]   = crossingPosition(a_surface, previous, s_edgeTolerance);
    polygon.m_vertexEdge[polygon.m_numVertices] = -1;
    polygon.m_vertex[polygon.m_numVertices++]   = cornerPosition(faceCorner[i]);
    polygon.m_vertexEdge[polygon.m_numVertices] = faceEdge[i];
    polygon.m_vertex[polygon.m_numVertices++]   = crossingPosition(a_surface, faceEdge[i], s_edgeTolerance);

    this->orientOutward(polygon, a_dir, a_side);

    Real     area = 0.0;
    RealVect vector;
    RealVect centroid;

    polygonMoments(polygon.m_vertex, polygon.m_numVertices, area, vector, centroid);

    if (area > 0.0) {
      a_out[numOut++] = polygon;
    }
  }

  return numOut;
#else
  (void)a_dir;
  (void)a_side;
  (void)a_surface;
  (void)a_out;

  return -1;
#endif
}

void
CutCellBody::defineDegenerate(const CutCellSurface& a_surface, const Kind a_kind) noexcept
{
  m_volumeFraction = (a_kind == Kind::Regular) ? 1.0 : 0.0;

  for (int d = 0; d < SpaceDim; d++) {
    for (int side = 0; side < 2; side++) {
      const int face = 2 * d + side;

      if (a_kind == Kind::Covered) {
        m_areaFraction[face] = 0.0;

        continue;
      }

      int faceCorner[1 << (SpaceDim - 1)];
      faceCorners(d, side, faceCorner);

      bool allZero = true;

      for (int k = 0; k < (1 << (SpaceDim - 1)); k++) {
        allZero = allZero && (a_surface.m_corner[faceCorner[k]] == 0.0);
      }

      // the face lies in the interface exactly when all of its corners do
      m_areaFraction[face] = allZero ? 0.0 : 1.0;
    }
  }

  RealVect areaVector = RealVect::Zero;

  for (int d = 0; d < SpaceDim; d++) {
    areaVector[d] = m_areaFraction[2 * d + 1] - m_areaFraction[2 * d];
  }

  const Real length = areaVector.vectorLength();

  if (length > 0.0) {
    m_boundaryArea     = length;
    m_trueBoundaryArea = length;
    m_normal           = areaVector / length;

    // the boundary here is the covered face, so its centroid is that face's own centre
    if (a_kind == Kind::Regular) {
      for (int d = 0; d < SpaceDim; d++) {
        for (int side = 0; side < 2; side++) {
          if (m_areaFraction[2 * d + side] == 0.0) {
            m_boundaryCentroid    = RealVect::Zero;
            m_boundaryCentroid[d] = -0.5 + static_cast<Real>(side);
          }
        }
      }
    }
  }
}

void
CutCellBody::accumulateMoments() noexcept
{
  Real     faceArea[s_numFaces] = {0.0};
  RealVect faceMoment[s_numFaces];

  for (int f = 0; f < s_numFaces; f++) {
    faceMoment[f] = RealVect::Zero;
  }

  RealVect boundaryVector = RealVect::Zero;
  RealVect boundaryMoment = RealVect::Zero;
  Real     boundaryArea   = 0.0;
  RealVect volumeMoment   = RealVect::Zero;
  Real     volume         = 0.0;

  m_closure = RealVect::Zero;

#if CH_SPACEDIM == 2
  // the fluid region is a single polygon, and each of its segments is either an aperture or the
  // chord the interface follows
  if (m_numPolygons == 1) {
    const Polygon& polygon = m_polygon[0];

    Real twiceArea = 0.0;

    for (int i = 0; i < polygon.m_numVertices; i++) {
      const RealVect& a = polygon.m_vertex[i];
      const RealVect& b = polygon.m_vertex[(i + 1) % polygon.m_numVertices];

      const Real wedge = a[0] * b[1] - b[0] * a[1];

      twiceArea += wedge;
      volumeMoment += wedge * (a + b);

      const RealVect outward = RealVect(D_DECL(b[1] - a[1], a[0] - b[0], 0.0));
      const Real     length  = (b - a).vectorLength();

      m_closure += outward;

      if (length <= 0.0) {
        continue;
      }

      const RealVect midpoint = 0.5 * (a + b);

      const int face = polygon.m_segmentFace[i];

      if (face >= 0) {
        const Real signedLength = outward[face / 2] * ((face % 2 == 0) ? -1.0 : 1.0);

        faceArea[face] += signedLength;
        faceMoment[face] += signedLength * midpoint;
      }
      else {
        boundaryArea += length;
        boundaryVector += outward;
        boundaryMoment += length * midpoint;
      }
    }

    volume       = 0.5 * twiceArea;
    volumeMoment = volumeMoment / 6.0;
  }
#else
  for (int i = 0; i < m_numPolygons; i++) {
    const Polygon& polygon = m_polygon[i];

    Real     area = 0.0;
    RealVect vector;
    RealVect centroid;

    polygonMoments(polygon.m_vertex, polygon.m_numVertices, area, vector, centroid);

    m_closure += vector;

    if (area > 0.0) {
      if (polygon.m_face >= 0) {
        // signed: where the interface lies in a cell face it bounds a hole in that face, and
        // the aperture is the net area open to flux
        const int  d          = polygon.m_face / 2;
        const Real signedArea = vector[d] * ((polygon.m_face % 2 == 0) ? -1.0 : 1.0);

        faceArea[polygon.m_face] += signedArea;
        faceMoment[polygon.m_face] += signedArea * centroid;
      }
      else {
        boundaryArea += area;
        boundaryVector += vector;
        boundaryMoment += area * centroid;
      }
    }

    for (int k = 1; k < polygon.m_numVertices - 1; k++) {
      const RealVect& x = polygon.m_vertex[0];
      const RealVect& y = polygon.m_vertex[k];
      const RealVect& z = polygon.m_vertex[k + 1];

      const RealVect c = RealVect(
        D_DECL(y[1] * z[2] - y[2] * z[1], y[2] * z[0] - y[0] * z[2], y[0] * z[1] - y[1] * z[0]));

      const Real tet = x.dotProduct(c) / 6.0;

      volume += tet;
      volumeMoment += (0.25 * tet) * (x + y + z);
    }
  }
#endif

  m_volumeFraction = volume;

  if (std::abs(volume) > 0.0) {
    m_volumeCentroid = volumeMoment / volume;
  }

  for (int f = 0; f < s_numFaces; f++) {
    if (std::abs(faceArea[f]) <= s_nullArea) {
      m_areaFraction[f] = 0.0;
      m_faceCentroid[f] = RealVect::Zero;
    }
    else {
      m_areaFraction[f]        = faceArea[f];
      m_faceCentroid[f]        = faceMoment[f] / faceArea[f];
      m_faceCentroid[f][f / 2] = 0.0;
    }
  }

  m_trueBoundaryArea = boundaryArea;

  // EBData stores the magnitude of the area vector, not the area of the patch: that is what
  // makes sum(alpha_hi - alpha_lo) equal a_B*n exactly for any closed body
  m_boundaryArea = boundaryVector.vectorLength();

  if (boundaryArea > 0.0) {
    m_boundaryCentroid = boundaryMoment / boundaryArea;
  }

  if (m_boundaryArea > 0.0) {
    m_normal = -boundaryVector / m_boundaryArea;
  }
}

#if CH_SPACEDIM == 2
bool
CutCellBody::defineCut(const CutCellSurface& a_surface) noexcept
{
  // In two dimensions a cell's faces and its edges are the same four segments, so the fluid
  // region is one polygon and the circuit is walked once. Each segment of that polygon either
  // lies in a cell face, where it is an aperture, or it is the chord the interface follows.
  constexpr int ringCorner[4] = {0, 1, 3, 2};
  constexpr int ringEdge[4]   = {0, 3, 1, 2};

  Polygon polygon;
  polygon.m_numVertices = 0;
  polygon.m_face        = -1;

  // a vertex carries the circuit position it came from, negative for a crossing and
  // non-negative for a corner, so that each segment can be attributed to the cell face whose
  // stretch of the circuit produced it
  int ring[s_maxVertices];

  for (int i = 0; i < 4; i++) {
    if (isFluid(a_surface.m_corner[ringCorner[i]])) {
      ring[polygon.m_numVertices]                 = i;
      polygon.m_vertexEdge[polygon.m_numVertices] = -1;
      polygon.m_vertex[polygon.m_numVertices++]   = cornerPosition(ringCorner[i]);
    }

    if (a_surface.hasCrossing(ringEdge[i])) {
      ring[polygon.m_numVertices]                 = -(i + 1);
      polygon.m_vertexEdge[polygon.m_numVertices] = ringEdge[i];
      polygon.m_vertex[polygon.m_numVertices++]   = crossingPosition(a_surface, ringEdge[i], s_edgeTolerance);
    }
  }

  // a segment lies in a cell face unless both its ends are crossings, which is the interface.
  // The face is the one holding the stretch of circuit the segment came from.
  for (int i = 0; i < polygon.m_numVertices; i++) {
    const int here = ring[i];
    const int next = ring[(i + 1) % polygon.m_numVertices];

    if (here < 0 && next < 0) {
      polygon.m_segmentFace[i] = -1;
    }
    else {
      const int position = (here >= 0) ? here : (-here - 1);
      const int edge     = ringEdge[position];
      const int dir      = edgeDirection(edge);

      int offset[SpaceDim];
      edgeOrigin(edge, offset);

      polygon.m_segmentFace[i] = 2 * (1 - dir) + offset[1 - dir];
    }
  }

  if (polygon.m_numVertices < 3) {
    return false;
  }

  // orient counter-clockwise, so that the outward normal of the segment from a to b is
  // (b[1] - a[1], a[0] - b[0])
  Real twiceArea = 0.0;

  for (int i = 0; i < polygon.m_numVertices; i++) {
    const RealVect& a = polygon.m_vertex[i];
    const RealVect& b = polygon.m_vertex[(i + 1) % polygon.m_numVertices];

    twiceArea += a[0] * b[1] - b[0] * a[1];
  }

  if (twiceArea < 0.0) {
    const int n = polygon.m_numVertices;

    int reversed[s_maxVertices];

    // reversing the circuit turns the segment leaving vertex i into the one arriving at it
    for (int i = 0; i < n; i++) {
      reversed[i] = polygon.m_segmentFace[(n - 2 - i + n) % n];
    }

    for (int i = 0; i < n / 2; i++) {
      const int j = n - 1 - i;

      std::swap(polygon.m_vertex[i], polygon.m_vertex[j]);
      std::swap(polygon.m_vertexEdge[i], polygon.m_vertexEdge[j]);
    }

    for (int i = 0; i < n; i++) {
      polygon.m_segmentFace[i] = reversed[i];
    }
  }

  m_polygon[0]  = polygon;
  m_numPolygons = 1;

  return true;
}
#else
bool
CutCellBody::defineCut(const CutCellSurface& a_surface) noexcept
{
  m_numPolygons = 0;

  for (int d = 0; d < SpaceDim; d++) {
    for (int side = 0; side < 2; side++) {
      if (m_numPolygons + 2 > s_maxPolygons) {
        return false;
      }

      const int numWalked = this->faceWalk(d, side, a_surface, &m_polygon[m_numPolygons]);

      if (numWalked < 0) {
        return false;
      }

      m_numPolygons += numWalked;
    }
  }

  // the face polygons occupy the front of the body, so a loop is oriented against them without
  // walking the faces a second time
  const int numFacePolygons = m_numPolygons;

  int loop[CutCellSurface::s_numEdges];
  int start[CutCellSurface::s_numEdges + 1];

  const int numLoops = crossingLoops(a_surface, loop, start);

  if (numLoops <= 0) {
    return false;
  }

  for (int l = 0; l < numLoops; l++) {
    const int begin = start[l];
    const int end   = start[l + 1];

    // a chord is shared by the patch and one face polygon, and a closed surface traverses a
    // shared edge once each way, so the face fixes the loop's direction
    bool oriented = false;

    for (int f = 0; f < numFacePolygons && !oriented; f++) {
      const Polygon& face = m_polygon[f];

      for (int i = 0; i < face.m_numVertices && !oriented; i++) {
        const int x = face.m_vertexEdge[i];
        const int y = face.m_vertexEdge[(i + 1) % face.m_numVertices];

        if (x < 0 || y < 0) {
          continue;
        }

        int indexX = -1;
        int indexY = -1;

        for (int j = begin; j < end; j++) {
          if (loop[j] == x) {
            indexX = j;
          }

          if (loop[j] == y) {
            indexY = j;
          }
        }

        if (indexX < 0 || indexY < 0) {
          continue;
        }

        const int next = begin + ((indexY - begin + 1) % (end - begin));

        if (loop[next] != x) {
          for (int q = 0; q < (end - begin) / 2; q++) {
            std::swap(loop[begin + q], loop[end - 1 - q]);
          }
        }

        oriented = true;
      }
    }

    RealVect apex = RealVect::Zero;

    for (int i = begin; i < end; i++) {
      apex += crossingPosition(a_surface, loop[i], s_edgeTolerance);
    }

    apex /= static_cast<Real>(end - begin);

    for (int i = begin; i < end; i++) {
      const int nextInLoop = begin + ((i - begin + 1) % (end - begin));

      Polygon triangle;
      triangle.m_numVertices = 3;
      triangle.m_face        = -1;
      triangle.m_vertex[0]   = apex;
      triangle.m_vertex[1]   = crossingPosition(a_surface, loop[i], s_edgeTolerance);
      triangle.m_vertex[2]   = crossingPosition(a_surface, loop[nextInLoop], s_edgeTolerance);

      for (int k = 0; k < 3; k++) {
        triangle.m_vertexEdge[k] = -1;
      }

      Real     area = 0.0;
      RealVect vector;
      RealVect centroid;

      polygonMoments(triangle.m_vertex, triangle.m_numVertices, area, vector, centroid);

      if (area >= s_nullArea) {
        if (m_numPolygons >= s_maxPolygons) {
          return false;
        }

        m_polygon[m_numPolygons++] = triangle;
      }
    }
  }

  return true;
}
#endif

bool
CutCellBody::define(const CutCellSurface& a_surface) noexcept
{
  *this = CutCellBody();

  const Kind kind = CutCellBody::classify(a_surface);

  if (kind != Kind::Cut) {
    this->defineDegenerate(a_surface, kind);

    return true;
  }

  if (!this->defineCut(a_surface)) {
    return false;
  }

  this->accumulateMoments();

  // verify rather than assume: a folded patch leaves the body open or the volume outside its
  // range, and every individual moment can still look plausible when it does
  const bool closed  = this->closureResidual() <= 1.0E-9;
  const bool inRange = m_volumeFraction >= -1.0E-12 && m_volumeFraction <= 1.0 + 1.0E-12;

  return closed && inRange;
}

Real
CutCellBody::volumeFraction() const noexcept
{
  return m_volumeFraction;
}

const RealVect&
CutCellBody::volumeCentroid() const noexcept
{
  return m_volumeCentroid;
}

Real
CutCellBody::areaFraction(const int a_dir, const Side::LoHiSide a_side) const noexcept
{
  return m_areaFraction[2 * a_dir + ((a_side == Side::Lo) ? 0 : 1)];
}

const RealVect&
CutCellBody::faceCentroid(const int a_dir, const Side::LoHiSide a_side) const noexcept
{
  return m_faceCentroid[2 * a_dir + ((a_side == Side::Lo) ? 0 : 1)];
}

Real
CutCellBody::boundaryArea() const noexcept
{
  return m_boundaryArea;
}

Real
CutCellBody::trueBoundaryArea() const noexcept
{
  return m_trueBoundaryArea;
}

const RealVect&
CutCellBody::normal() const noexcept
{
  return m_normal;
}

const RealVect&
CutCellBody::boundaryCentroid() const noexcept
{
  return m_boundaryCentroid;
}

Real
CutCellBody::closureResidual() const noexcept
{
  return m_closure.vectorLength();
}

Real
CutCellBody::divergenceResidual() const noexcept
{
  RealVect apertureVector = RealVect::Zero;

  for (int d = 0; d < SpaceDim; d++) {
    apertureVector[d] = m_areaFraction[2 * d + 1] - m_areaFraction[2 * d];
  }

  RealVect residual = apertureVector - m_boundaryArea * m_normal;

  return residual.vectorLength();
}

#include <CD_NamespaceFooter.H>
