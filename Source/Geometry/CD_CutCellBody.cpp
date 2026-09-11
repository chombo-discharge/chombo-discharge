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
#include <CD_PolyhedralEBUtils.H>
#include <CD_NamespaceHeader.H>

namespace PolyhedralEB {

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
    if (isFluid(a_surface.m_corner[c]) != a_fluidSide) {
      continue;
    }

    any = true;

    // the corner is on the interface when every edge leaving it towards the other side turns
    // over within the distance crossings are held off the endpoints by. Asking where the
    // crossing is, rather than how small the corner value is, keeps this a question about the
    // surface rather than about a magnitude
    for (int d = 0; d < SpaceDim; d++) {
      const int other = c ^ (1 << d);

      if (isFluid(a_surface.m_corner[other]) == a_fluidSide) {
        continue;
      }

      const int low = c & ~(1 << d);

      int offset[SpaceDim];

      for (int k = 0; k < SpaceDim; k++) {
        offset[k] = (low >> k) & 1;
      }

      const int  edge = detail::edgeIndex(d, offset);
      const Real t    = a_surface.m_crossing[edge];

      if (!a_surface.hasCrossing(edge)) {
        return false;
      }

      const Real distance = (((c >> d) & 1) == 1) ? (1.0 - t) : t;

      if (distance > s_edgeTolerance) {
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

  detail::polygonMoments(a_polygon.m_vertex, a_polygon.m_numVertices, area, vector, centroid);

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

  detail::faceEdges(a_dir, a_side, faceEdge);
  detail::faceCorners(a_dir, a_side, faceCorner);

  int numHit = 0;

  for (int i = 0; i < 4; i++) {
    numHit += a_surface.hasCrossing(faceEdge[i]) ? 1 : 0;
  }

  int       pair[2][2] = {{-1, -1}, {-1, -1}};
  const int numPairs   = detail::facePairs(a_dir, a_side, a_surface, pair);

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
        polygon.m_vertex[polygon.m_numVertices++]   = detail::cornerPosition(faceCorner[i]);
      }

      if (a_surface.hasCrossing(faceEdge[i])) {
        polygon.m_vertexEdge[polygon.m_numVertices] = faceEdge[i];
        polygon.m_vertex[polygon.m_numVertices++]   = detail::crossingPosition(a_surface, faceEdge[i], s_edgeTolerance);
      }
    }

    if (polygon.m_numVertices < 3) {
      return 0;
    }

    this->orientOutward(polygon, a_dir, a_side);

    Real     area = 0.0;
    RealVect vector;
    RealVect centroid;

    detail::polygonMoments(polygon.m_vertex, polygon.m_numVertices, area, vector, centroid);

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
        polygon.m_vertex[polygon.m_numVertices++]   = detail::cornerPosition(faceCorner[i]);
      }

      polygon.m_vertexEdge[polygon.m_numVertices] = faceEdge[i];
      polygon.m_vertex[polygon.m_numVertices++]   = detail::crossingPosition(a_surface, faceEdge[i], s_edgeTolerance);
    }

    this->orientOutward(polygon, a_dir, a_side);

    Real     area = 0.0;
    RealVect vector;
    RealVect centroid;

    detail::polygonMoments(polygon.m_vertex, polygon.m_numVertices, area, vector, centroid);

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
    polygon.m_vertex[polygon.m_numVertices++]   = detail::crossingPosition(a_surface, previous, s_edgeTolerance);
    polygon.m_vertexEdge[polygon.m_numVertices] = -1;
    polygon.m_vertex[polygon.m_numVertices++]   = detail::cornerPosition(faceCorner[i]);
    polygon.m_vertexEdge[polygon.m_numVertices] = faceEdge[i];
    polygon.m_vertex[polygon.m_numVertices++]   = detail::crossingPosition(a_surface, faceEdge[i], s_edgeTolerance);

    this->orientOutward(polygon, a_dir, a_side);

    Real     area = 0.0;
    RealVect vector;
    RealVect centroid;

    detail::polygonMoments(polygon.m_vertex, polygon.m_numVertices, area, vector, centroid);

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
      detail::faceCorners(d, side, faceCorner);

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

    detail::polygonMoments(polygon.m_vertex, polygon.m_numVertices, area, vector, centroid);

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
      polygon.m_vertex[polygon.m_numVertices++]   = detail::cornerPosition(ringCorner[i]);
    }

    if (a_surface.hasCrossing(ringEdge[i])) {
      ring[polygon.m_numVertices]                 = -(i + 1);
      polygon.m_vertexEdge[polygon.m_numVertices] = ringEdge[i];
      polygon.m_vertex[polygon.m_numVertices++]   = detail::crossingPosition(a_surface, ringEdge[i], s_edgeTolerance);
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
      const int dir      = detail::edgeDirection(edge);

      int offset[SpaceDim];
      detail::edgeOrigin(edge, offset);

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

  const int numLoops = detail::crossingLoops(a_surface, loop, start);

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
      apex += detail::crossingPosition(a_surface, loop[i], s_edgeTolerance);
    }

    apex /= static_cast<Real>(end - begin);

    for (int i = begin; i < end; i++) {
      const int nextInLoop = begin + ((i - begin + 1) % (end - begin));

      Polygon triangle;
      triangle.m_numVertices = 3;
      triangle.m_face        = -1;
      triangle.m_vertex[0]   = apex;
      triangle.m_vertex[1]   = detail::crossingPosition(a_surface, loop[i], s_edgeTolerance);
      triangle.m_vertex[2]   = detail::crossingPosition(a_surface, loop[nextInLoop], s_edgeTolerance);

      for (int k = 0; k < 3; k++) {
        triangle.m_vertexEdge[k] = -1;
      }

      Real     area = 0.0;
      RealVect vector;
      RealVect centroid;

      detail::polygonMoments(triangle.m_vertex, triangle.m_numVertices, area, vector, centroid);

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

} // namespace PolyhedralEB

#include <CD_NamespaceFooter.H>
