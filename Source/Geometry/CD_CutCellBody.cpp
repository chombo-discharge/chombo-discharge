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
#include <algorithm>
#include <cmath>
#include <cstddef>

// Chombo includes
#include <CH_assert.H>
#include <PolyGeom.H>

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
  // one predicate decides every corner, and an exact zero is on the solid side of it. A cell with
  // a face in the interface therefore reads as cut from the fluid side -- full, with that face as
  // its interface, the crossings on the edges leaving the plane placed exactly on the corners --
  // and as covered from the solid side
  const bool firstFluid = isFluid(a_surface.m_corner[0]);

  bool cut = false;

  for (int c = 1; c < CutCellSurface::s_numCorners && !cut; c++) {
    cut = isFluid(a_surface.m_corner[c]) != firstFluid;
  }

  if (!cut) {
    return firstFluid ? Kind::Regular : Kind::Covered;
  }

  return Kind::Cut;
}

int
CutCellBody::numSheets(const CutCellSurface& a_surface) noexcept
{
#if CH_SPACEDIM == 2
  // Two dimensions: a cell's faces and its edges are the same four segments, and the crossings on them pair
  // into chords, one around each run of solid corners met going round the cell.
  int numCrossings = 0;

  for (int e = 0; e < CutCellSurface::s_numEdges; e++) {
    if (a_surface.hasCrossing(e)) {
      numCrossings++;
    }
  }

  if (numCrossings % 2 != 0) {
    return -1;
  }

  if (numCrossings == 0) {
    return 0;
  }

  constexpr int circuit[4] = {0, 1, 3, 2};

  // start the walk at a fluid corner, so that every run of solid corners is met whole
  int first = 0;

  while (!isFluid(a_surface.m_corner[circuit[first]])) {
    first++;
  }

  int  numSheets = 0;
  bool inRun     = false;
  bool allZero   = true;

  for (int i = 1; i <= 4; i++) {
    const Real value = a_surface.m_corner[circuit[(first + i) % 4]];

    if (!isFluid(value)) {
      allZero = (inRun ? allZero : true) && (value == 0.0);
      inRun   = true;
    }
    else if (inRun) {
      numSheets += allZero ? 0 : 1;
      inRun = false;
    }
  }

  return numSheets;
#else
  int loop[CutCellSurface::s_numEdges];
  int start[CutCellSurface::s_numEdges + 1];

  const int numLoops = detail::crossingLoops(a_surface, loop, start);

  if (numLoops <= 0) {
    return numLoops;
  }

  // A loop whose every edge has, for its solid end, a corner at exactly zero bounds nothing.
  int numSheets = 0;

  for (int l = 0; l < numLoops; l++) {
    bool degenerate = true;

    for (int k = start[l]; k < start[l + 1] && degenerate; k++) {
      int low  = -1;
      int high = -1;

      detail::edgeCorners(loop[k], low, high);

      const Real solidEnd = isFluid(a_surface.m_corner[low]) ? a_surface.m_corner[high] : a_surface.m_corner[low];

      degenerate = (solidEnd == 0.0);
    }

    numSheets += degenerate ? 0 : 1;
  }

  return numSheets;
#endif
}

#if CH_SPACEDIM == 3
void
CutCellBody::orientOutward(Polygon& a_polygon, const int a_dir, const int a_side) const noexcept
{
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);
  CH_assert(a_side == 0 || a_side == 1);
  CH_assert(a_polygon.m_numVertices >= 3);

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
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);
  CH_assert(a_side == 0 || a_side == 1);
  CH_assert(a_out != nullptr);

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

    CH_assert(polygon.m_numVertices <= s_maxVertices);

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

  CH_assert(numIsolated == 2);

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

    CH_assert(polygon.m_numVertices <= s_maxVertices);

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
}
#endif

void
CutCellBody::defineDegenerate(const CutCellSurface& a_surface, const Kind a_kind) noexcept
{
  CH_assert(a_kind != Kind::Cut);

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

    // the boundary here is the covered face, so its centroid is that face's own centre. A
    // regular cell has at most one such face: two would leave a corner with no edge towards
    // the fluid side, which classify reads as cut.
    if (a_kind == Kind::Regular) {
      int numCovered = 0;

      for (int d = 0; d < SpaceDim; d++) {
        for (int side = 0; side < 2; side++) {
          if (m_areaFraction[2 * d + side] == 0.0) {
            numCovered++;

            m_boundaryCentroid    = RealVect::Zero;
            m_boundaryCentroid[d] = -0.5 + static_cast<Real>(side);
          }
        }
      }

      CH_assert(numCovered == 1);
    }
  }
}

#if CH_SPACEDIM == 3
bool
CutCellBody::mergeCoplanar(const Polygon* a_in, const int a_num, Polygon* a_out, const int a_maxOut, int& a_numOut)
  const noexcept
{
  CH_assert(a_in != nullptr);
  CH_assert(a_out != nullptr);
  CH_assert(a_num >= 0);
  CH_assert(a_maxOut > 0);

  a_numOut = 0;

  // A seam face that the body covers completely carries no polygon, and that is an answer, not a
  // failure. It happens wherever a face of the geometry lands on a cell face: the crossings sit on
  // the face's own edges, faceWalk returns slivers of no area, and everything below cancels them
  // against each other. Refusing here would drop the cell back to a single chord, which is the
  // crack the multichord exists to remove.
  Real inputArea = 0.0;

  for (int ip = 0; ip < a_num; ip++) {
    RealVect twice = RealVect::Zero;

    for (int i = 0; i < a_in[ip].m_numVertices; i++) {
      const RealVect& a = a_in[ip].m_vertex[i];
      const RealVect& b = a_in[ip].m_vertex[(i + 1) % a_in[ip].m_numVertices];

      twice += PolyGeom::cross(a, b);
    }

    inputArea += 0.5 * twice.vectorLength();
  }

  // the face is one unit square in these coordinates, so this is a sliver a weld tolerance wide
  if (inputArea <= PolyhedralEB::detail::s_weldTolerance) {
    return true;
  }

  RealVect from[4 * s_maxVertices];
  RealVect to[4 * s_maxVertices];

  int numEdges = 0;

  for (int ip = 0; ip < a_num; ip++) {
    for (int i = 0; i < a_in[ip].m_numVertices; i++) {
      const RealVect& a = a_in[ip].m_vertex[i];
      const RealVect& b = a_in[ip].m_vertex[(i + 1) % a_in[ip].m_numVertices];

      if (detail::sameVertex(a, b)) {
        continue;
      }

      bool interior = false;

      for (int jp = 0; jp < a_num && !interior; jp++) {
        if (jp == ip) {
          continue;
        }

        for (int j = 0; j < a_in[jp].m_numVertices && !interior; j++) {
          const RealVect& c = a_in[jp].m_vertex[j];
          const RealVect& d = a_in[jp].m_vertex[(j + 1) % a_in[jp].m_numVertices];

          interior = detail::sameVertex(a, d) && detail::sameVertex(b, c);
        }
      }

      if (!interior) {
        if (numEdges >= 4 * s_maxVertices) {
          return false;
        }

        from[numEdges] = a;
        to[numEdges]   = b;
        numEdges++;
      }
    }
  }

  if (numEdges < 3) {
    return false;
  }

  bool used[4 * s_maxVertices] = {false};

  // The union need not be one connected region. A cube's edge grazing the corner of a quadrant
  // leaves a sliver detached from the main patch, and a face is allowed to carry a polygon for
  // each: accumulateMoments sums over polygons, so nothing has to be joined that geometry has
  // separated. Refusing here instead would throw away exactly the cells the multichord is for.
  for (int seed = 0; seed < numEdges; seed++) {
    if (used[seed]) {
      continue;
    }

    RealVect walk[4 * s_maxVertices];

    int numWalk = 0;

    used[seed]      = true;
    walk[numWalk++] = from[seed];

    RealVect       current = to[seed];
    const RealVect end     = from[seed];

    bool closed = false;

    for (int guard = 0; guard <= numEdges && !closed; guard++) {
      if (detail::sameVertex(current, end)) {
        closed = true;

        break;
      }

      int next = -1;

      for (int j = 0; j < numEdges && next < 0; j++) {
        if (!used[j] && detail::sameVertex(from[j], current)) {
          next = j;
        }
      }

      if (next < 0 || numWalk >= 4 * s_maxVertices) {
        return false;
      }

      used[next]      = true;
      walk[numWalk++] = current;
      current         = to[next];
    }

    if (!closed) {
      return false;
    }

    if (a_numOut >= a_maxOut) {
      return false;
    }

    // drop vertices that sit on the straight line between their neighbours
    Polygon& out = a_out[a_numOut];

    out.m_numVertices = 0;
    out.m_face        = a_in[0].m_face;

    for (int i = 0; i < numWalk; i++) {
      const RealVect& prev = walk[(i + numWalk - 1) % numWalk];
      const RealVect& here = walk[i];
      const RealVect& next = walk[(i + 1) % numWalk];

      const RealVect back  = here - prev;
      const RealVect ahead = next - here;

      const Real backLength  = back.vectorLength();
      const Real aheadLength = ahead.vectorLength();

      bool straight = false;

      if (backLength > 0.0 && aheadLength > 0.0) {
        const RealVect unit = back / backLength;

        const Real along = ahead.dotProduct(unit);

        straight = (along > 0.0) && ((ahead - along * unit).vectorLength() <= 1.0E-11 * aheadLength);
      }

      // A chord vertex on the boundary between two children is a vertex of the cells on the other side of the
      // seam, and dropping it would leave their two segments meeting the middle of one of ours: watertight, but
      // not a shared edge. The children sit at plus and minus a quarter, so their boundaries in the face are at
      // exactly zero. Such a vertex is kept even where it falls on a straight run, and even where it sits on the
      // face's own boundary and so splits an edge that the neighbouring face carries whole -- weldTJunctions
      // puts it into that face too before the interface is closed.
      const int faceDir = a_in[0].m_face / 2;

      bool onChildBoundary = false;

      for (int d = 0; d < SpaceDim; d++) {
        if (d == faceDir) {
          continue;
        }

        onChildBoundary = onChildBoundary || (std::abs(here[d]) <= detail::s_weldTolerance);
      }

      if (straight && !onChildBoundary) {
        continue;
      }

      if (out.m_numVertices >= s_maxVertices) {
        return false;
      }

      out.m_vertexEdge[out.m_numVertices]  = -1;
      out.m_segmentFace[out.m_numVertices] = a_in[0].m_face;
      out.m_vertex[out.m_numVertices++]    = here;
    }

    // a loop that collapses under the collinear pass enclosed nothing
    if (out.m_numVertices >= 3) {
      a_numOut++;
    }
  }

  if (a_numOut == 0) {
    return false;
  }

  return true;
}

bool
CutCellBody::restrictFace(const CutCellSurface* a_children,
                          const bool*           a_closedQuadrant,
                          const int             a_dir,
                          const int             a_side) noexcept
{
  CH_assert(a_children != nullptr);
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);
  CH_assert(a_side == 0 || a_side == 1);
  CH_assert(m_numPolygons >= 0 && m_numPolygons <= s_maxPolygons);

  const int face = 2 * a_dir + a_side;

  // this face's chord goes, and so does the interface, which was built to meet it
  int kept = 0;

  for (int ip = 0; ip < m_numPolygons; ip++) {
    if (m_polygon[ip].m_face != face && m_polygon[ip].m_face >= 0) {
      m_polygon[kept++] = m_polygon[ip];
    }
  }

  m_numPolygons = kept;

  Polygon sub[4 * (1 << (SpaceDim - 1))];

  int numSub = 0;

  for (int q = 0; q < (1 << (SpaceDim - 1)); q++) {
    // the child of this cell that covers quadrant q of the face
    int which = 0;
    int bit   = 0;

    for (int d = 0; d < SpaceDim; d++) {
      if (d == a_dir) {
        which |= a_side << d;
      }
      else {
        which |= ((q >> bit) & 1) << d;
        bit++;
      }
    }

    RealVect origin;

    for (int d = 0; d < SpaceDim; d++) {
      origin[d] = -0.25 + 0.5 * static_cast<Real>((which >> d) & 1);
    }

    Polygon walked[2];

    const int numWalked = this->faceWalk(a_dir, a_side, a_children[q], walked);

    if (numWalked < 0) {
      return false;
    }

    if (a_closedQuadrant != nullptr && a_closedQuadrant[q]) {
      continue;
    }

    for (int n = 0; n < numWalked; n++) {
      if (numSub >= 4 * (1 << (SpaceDim - 1))) {
        return false;
      }

      Polygon& to = sub[numSub];

      to = walked[n];

      for (int iv = 0; iv < to.m_numVertices; iv++) {
        to.m_vertex[iv] = origin + 0.5 * walked[n].m_vertex[iv];
      }

      numSub++;
    }
  }

  // one polygon for the face, not one per child
  if (numSub > 0) {
    Polygon merged[1 << (SpaceDim - 1)];

    int numMerged = 0;

    if (!this->mergeCoplanar(sub, numSub, merged, 1 << (SpaceDim - 1), numMerged)) {
      return false;
    }

    for (int n = 0; n < numMerged; n++) {
      if (m_numPolygons >= s_maxPolygons) {
        return false;
      }

      m_polygon[m_numPolygons++] = merged[n];
    }
  }

  return true;
}

void
CutCellBody::defineWhole() noexcept
{
  *this = CutCellBody();

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (int side = 0; side < 2; side++) {
      int faceCorner[1 << (SpaceDim - 1)];

      detail::faceCorners(dir, side, faceCorner);

      Polygon& polygon = m_polygon[m_numPolygons];

      polygon               = Polygon();
      polygon.m_numVertices = 0;
      polygon.m_face        = 2 * dir + side;

      for (int i = 0; i < (1 << (SpaceDim - 1)); i++) {
        polygon.m_vertexEdge[polygon.m_numVertices] = -1;
        polygon.m_vertex[polygon.m_numVertices++]   = detail::cornerPosition(faceCorner[i]);
      }

      this->orientOutward(polygon, dir, side);

      m_numPolygons++;
    }
  }

  this->accumulateMoments();
}

void
CutCellBody::recordFace(const int a_face, const int a_reason, CutCellFaceOverrides& a_overrides) const
{
  CH_assert(a_face >= 0 && a_face < s_numFaces);

  a_overrides.beginFace(a_face, a_reason);

  if (a_reason == CutCellFaceOverrides::s_closed) {
    return;
  }

  for (int ip = 0; ip < m_numPolygons; ip++) {
    const Polygon& polygon = m_polygon[ip];

    if (polygon.m_face == a_face) {
      a_overrides.addPolygon(polygon.m_vertex, polygon.m_vertexEdge, polygon.m_numVertices);
    }
  }
}

bool
CutCellBody::replaceFace(const CutCellFaceOverrides& a_overrides, const int a_entry) noexcept
{
  CH_assert(a_overrides.reason(a_entry) == CutCellFaceOverrides::s_finer);
  CH_assert(m_numPolygons >= 0 && m_numPolygons <= s_maxPolygons);

  const int face = a_overrides.face(a_entry);

  // this face's chord goes, and so does the interface, which was built to meet it
  int kept = 0;

  for (int ip = 0; ip < m_numPolygons; ip++) {
    if (m_polygon[ip].m_face != face && m_polygon[ip].m_face >= 0) {
      m_polygon[kept++] = m_polygon[ip];
    }
  }

  m_numPolygons = kept;

  int polyBegin = 0;
  int polyEnd   = 0;

  a_overrides.polygons(a_entry, polyBegin, polyEnd);

  for (int p = polyBegin; p < polyEnd; p++) {
    int vertBegin = 0;
    int vertEnd   = 0;

    a_overrides.vertices(p, vertBegin, vertEnd);

    if (m_numPolygons >= s_maxPolygons || vertEnd - vertBegin > s_maxVertices) {
      return false;
    }

    Polygon& polygon = m_polygon[m_numPolygons++];

    polygon.m_face        = face;
    polygon.m_numVertices = vertEnd - vertBegin;

    for (int v = vertBegin; v < vertEnd; v++) {
      polygon.m_vertex[v - vertBegin]      = a_overrides.vertex(v);
      polygon.m_vertexEdge[v - vertBegin]  = a_overrides.vertexEdge(v);
      polygon.m_segmentFace[v - vertBegin] = -1;
    }
  }

  return true;
}

bool
CutCellBody::snapFace(const int a_dir, const int a_side, const bool a_neighbourIsFluid) noexcept
{
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);
  CH_assert(a_side == 0 || a_side == 1);

  const int face = 2 * a_dir + a_side;

  // this face's polygon goes, and so does the interface, which was built to meet its chord
  int kept = 0;

  for (int ip = 0; ip < m_numPolygons; ip++) {
    if (m_polygon[ip].m_face != face && m_polygon[ip].m_face >= 0) {
      m_polygon[kept++] = m_polygon[ip];
    }
  }

  m_numPolygons = kept;

  // A neighbour that holds no solid says the whole face is open; one that holds no fluid leaves it closed, and
  // then the face has no polygon at all.
  if (a_neighbourIsFluid) {
    if (m_numPolygons >= s_maxPolygons) {
      return false;
    }

    int faceCorner[1 << (SpaceDim - 1)];

    detail::faceCorners(a_dir, a_side, faceCorner);

    Polygon& polygon = m_polygon[m_numPolygons];

    polygon               = Polygon();
    polygon.m_numVertices = 0;
    polygon.m_face        = face;

    for (int i = 0; i < (1 << (SpaceDim - 1)); i++) {
      polygon.m_vertexEdge[polygon.m_numVertices] = -1;
      polygon.m_vertex[polygon.m_numVertices++]   = detail::cornerPosition(faceCorner[i]);
    }

    this->orientOutward(polygon, a_dir, a_side);

    m_numPolygons++;
  }

  if (!this->closeInterface()) {
    return false;
  }

  const bool closed  = this->closureResidual() <= 1.0E-9;
  const bool inRange = m_volumeFraction >= -1.0E-12 && m_volumeFraction <= 1.0 + 1.0E-12;

  return closed && inRange;
}

bool
CutCellBody::weldTJunctions() noexcept
{
  CH_assert(m_numPolygons >= 0 && m_numPolygons <= s_maxPolygons);

  for (int ip = 0; ip < m_numPolygons; ip++) {
    Polygon& p = m_polygon[ip];

    for (int i = 0; i < p.m_numVertices; i++) {
      const RealVect a = p.m_vertex[i];
      const RealVect b = p.m_vertex[(i + 1) % p.m_numVertices];

      const RealVect along  = b - a;
      const Real     length = along.vectorLength();

      if (length <= detail::s_weldTolerance) {
        continue;
      }

      // Every vertex of another polygon that lies on this edge, strictly between its ends.
      RealVect inside[s_maxVertices];
      Real     where[s_maxVertices];

      int numInside = 0;

      for (int jp = 0; jp < m_numPolygons; jp++) {
        if (jp == ip) {
          continue;
        }

        const Polygon& q = m_polygon[jp];

        for (int j = 0; j < q.m_numVertices; j++) {
          const RealVect& v = q.m_vertex[j];

          if (detail::sameVertex(v, a) || detail::sameVertex(v, b)) {
            continue;
          }

          const Real t = (v - a).dotProduct(along) / (length * length);

          if (t <= 0.0 || t >= 1.0) {
            continue;
          }

          if (((v - a) - t * along).vectorLength() > detail::s_weldTolerance) {
            continue;
          }

          bool have = false;

          for (int k = 0; k < numInside && !have; k++) {
            have = detail::sameVertex(inside[k], v);
          }

          if (have) {
            continue;
          }

          if (numInside >= s_maxVertices) {
            return false;
          }

          inside[numInside] = v;
          where[numInside]  = t;

          numInside++;
        }
      }

      if (numInside == 0) {
        continue;
      }

      // in order along the edge, so that the split reads as one walk from a to b
      for (int m = 1; m < numInside; m++) {
        const RealVect v = inside[m];
        const Real     t = where[m];

        int n = m - 1;

        while (n >= 0 && where[n] > t) {
          inside[n + 1] = inside[n];
          where[n + 1]  = where[n];

          n--;
        }

        inside[n + 1] = v;
        where[n + 1]  = t;
      }

      if (p.m_numVertices + numInside > s_maxVertices) {
        return false;
      }

      for (int n = p.m_numVertices - 1; n > i; n--) {
        p.m_vertex[n + numInside]      = p.m_vertex[n];
        p.m_vertexEdge[n + numInside]  = p.m_vertexEdge[n];
        p.m_segmentFace[n + numInside] = p.m_segmentFace[n];
      }

      for (int n = 0; n < numInside; n++) {
        p.m_vertex[i + 1 + n]      = inside[n];
        p.m_vertexEdge[i + 1 + n]  = -1;
        p.m_segmentFace[i + 1 + n] = p.m_segmentFace[i];
      }

      p.m_numVertices += numInside;

      i += numInside;
    }
  }

  return true;
}

bool
CutCellBody::closeInterface() noexcept
{
  if (!this->closeBoundary(-1)) {
    return false;
  }

  this->accumulateMoments();

  return true;
}

bool
CutCellBody::closeBoundary(const int a_face) noexcept
{
  CH_assert(a_face >= -1 && a_face < s_numFaces);
  CH_assert(m_numPolygons >= 0 && m_numPolygons <= s_maxPolygons);

  if (!this->weldTJunctions()) {
    return false;
  }

  RealVect from[s_maxPolygons * s_maxVertices];
  RealVect to[s_maxPolygons * s_maxVertices];

  int numOpen = 0;

  for (int ip = 0; ip < m_numPolygons; ip++) {
    const Polygon& p = m_polygon[ip];

    for (int i = 0; i < p.m_numVertices; i++) {
      const RealVect& a = p.m_vertex[i];
      const RealVect& b = p.m_vertex[(i + 1) % p.m_numVertices];

      if (detail::sameVertex(a, b)) {
        continue;
      }

      bool shared = false;

      for (int jp = 0; jp < m_numPolygons && !shared; jp++) {
        if (jp == ip) {
          continue;
        }

        const Polygon& q = m_polygon[jp];

        for (int j = 0; j < q.m_numVertices && !shared; j++) {
          const RealVect& c = q.m_vertex[j];
          const RealVect& d = q.m_vertex[(j + 1) % q.m_numVertices];

          shared = detail::sameVertex(a, d) && detail::sameVertex(b, c);
        }
      }

      if (!shared) {
        CH_assert(numOpen < s_maxPolygons * s_maxVertices);

        from[numOpen] = b;
        to[numOpen]   = a;
        numOpen++;
      }
    }
  }

  if (numOpen == 0) {
    return true;
  }

  bool used[s_maxPolygons * s_maxVertices] = {false};

  for (int s0 = 0; s0 < numOpen; s0++) {
    if (used[s0]) {
      continue;
    }

    used[s0] = true;

    RealVect loop[s_maxVertices];

    int numLoop = 0;

    loop[numLoop++] = from[s0];

    RealVect       current = to[s0];
    const RealVect end     = from[s0];

    bool closed = false;

    for (int guard = 0; guard <= numOpen && !closed; guard++) {
      if (detail::sameVertex(current, end)) {
        closed = true;

        break;
      }

      int next = -1;

      for (int j = 0; j < numOpen && next < 0; j++) {
        if (!used[j] && detail::sameVertex(from[j], current)) {
          next = j;
        }
      }

      if (next < 0 || numLoop >= s_maxVertices) {
        break;
      }

      used[next] = true;

      loop[numLoop++] = current;
      current         = to[next];
    }

    if (!closed || numLoop < 3) {
      return false;
    }

    // A patch closing an opening cut by a plane lies in that plane, so the loop is a polygon already and is
    // kept as one: fanning it would be exact too, but it costs a polygon per vertex rather than one in total,
    // and three cuts of a body that already holds a fanned interface do not fit. The interface's own loop is
    // not planar in general, which is what the fan is for.
    if (a_face >= 0) {
      if (m_numPolygons >= s_maxPolygons || numLoop > s_maxVertices) {
        return false;
      }

      Polygon& patch = m_polygon[m_numPolygons];

      patch.m_numVertices = numLoop;
      patch.m_face        = a_face;

      for (int i = 0; i < numLoop; i++) {
        patch.m_vertex[i]      = loop[i];
        patch.m_vertexEdge[i]  = -1;
        patch.m_segmentFace[i] = -1;
      }

      this->orientOutward(patch, a_face / 2, a_face % 2);

      m_numPolygons++;

      continue;
    }

    RealVect apex = RealVect::Zero;

    for (int i = 0; i < numLoop; i++) {
      apex += loop[i];
    }

    apex /= static_cast<Real>(numLoop);

    for (int i = 0; i < numLoop; i++) {
      if (m_numPolygons >= s_maxPolygons) {
        return false;
      }

      Polygon& t = m_polygon[m_numPolygons];

      t.m_numVertices = 3;
      t.m_face        = -1;

      t.m_vertex[0] = apex;
      t.m_vertex[1] = loop[i];
      t.m_vertex[2] = loop[(i + 1) % numLoop];

      for (int k = 0; k < 3; k++) {
        t.m_vertexEdge[k]  = -1;
        t.m_segmentFace[k] = -1;
      }

      m_numPolygons++;
    }
  }

  return true;
}

void
CutCellBody::printPolygons(std::ostream& a_out) const noexcept
{
  a_out << "polygons " << m_numPolygons << std::endl;

  for (int ip = 0; ip < m_numPolygons; ip++) {
    const Polygon& p = m_polygon[ip];

    a_out << "  polygon " << ip << " face " << p.m_face << " vertices " << p.m_numVertices << ":";

    for (int i = 0; i < p.m_numVertices; i++) {
      a_out << " (" << p.m_vertex[i][0] << "," << p.m_vertex[i][1] << "," << p.m_vertex[i][2] << ")";
    }

    a_out << std::endl;
  }
}

void
CutCellBody::appendInterfaceFacets(Vector<Real>&   a_facets,
                                   const IntVect&  a_cell,
                                   const RealVect& a_probLo,
                                   const Real      a_dx) const noexcept
{
  CH_assert(a_dx > 0.0);

  // a body that is not cut holds no interface polygon, so the loop appends nothing for it
  for (int ip = 0; ip < m_numPolygons; ip++) {
    const Polygon& p = m_polygon[ip];

    if (p.m_face >= 0 || p.m_numVertices < 3) {
      continue;
    }

    // the interface is already fanned, but a polygon carrying more than three vertices is fanned
    // again here rather than left for the reader to triangulate. A triangle without area is not
    // written: a surface tangent to a cell edge leaves the cells beside it an interface that is
    // that edge and nothing more, and its fan is a set of degenerate triangles carrying no moment
    for (int v = 1; v + 1 < p.m_numVertices; v++) {
      const RealVect* corner[3] = {&p.m_vertex[0], &p.m_vertex[v], &p.m_vertex[v + 1]};

      if (PolyGeom::cross(*corner[1] - *corner[0], *corner[2] - *corner[0]).vectorLength() <= s_nullArea) {
        continue;
      }

      for (int k = 0; k < 3; k++) {
        for (int d = 0; d < SpaceDim; d++) {
          // a face vertex sits at local 0.5 exactly, so this is an integer times the spacing from
          // either side of the face; a shared edge crossing has the same local offset in both cells
          a_facets.push_back(a_probLo[d] + a_dx * (static_cast<Real>(a_cell[d]) + ((*corner[k])[d] + 0.5)));
        }
      }
    }
  }
}

#endif

#if CH_SPACEDIM == 2
void
CutCellBody::defineWhole() noexcept
{
  *this = CutCellBody();

  // the four corners counter-clockwise, and the face each segment leaving one of them lies in
  constexpr int ringCorner[4] = {0, 1, 3, 2};
  constexpr int ringFace[4]   = {2, 1, 3, 0};

  Polygon& polygon = m_polygon[0];

  polygon.m_numVertices = 4;
  polygon.m_face        = -1;

  for (int i = 0; i < 4; i++) {
    polygon.m_vertex[i]      = detail::cornerPosition(ringCorner[i]);
    polygon.m_vertexEdge[i]  = -1;
    polygon.m_segmentFace[i] = ringFace[i];
  }

  m_numPolygons = 1;

  this->accumulateMoments();
}

void
CutCellBody::recordFace(const int a_face, const int a_reason, CutCellFaceOverrides& a_overrides) const
{
  CH_assert(a_face >= 0 && a_face < s_numFaces);

  // a face in two dimensions is changed in place, and what changed it is all there is to record
  a_overrides.beginFace(a_face, a_reason);
}

bool
CutCellBody::snapFace(const int a_dir, const int a_side, const bool a_neighbourIsFluid) noexcept
{
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);
  CH_assert(a_side == 0 || a_side == 1);
  CH_assert(m_numPolygons == 1);

  // Only closing: a face onto a cell that holds no fluid keeps its stretch of the polygon, which becomes interface
  // lying in the face.
  if (a_neighbourIsFluid) {
    return false;
  }

  const int face = 2 * a_dir + a_side;

  Polygon& polygon = m_polygon[0];

  for (int i = 0; i < polygon.m_numVertices; i++) {
    if (polygon.m_segmentFace[i] == face) {
      polygon.m_segmentFace[i] = -1;
    }
  }

  this->accumulateMoments();

  const bool closed  = this->closureResidual() <= 1.0E-9;
  const bool inRange = m_volumeFraction >= -1.0E-12 && m_volumeFraction <= 1.0 + 1.0E-12;

  return closed && inRange;
}

bool
CutCellBody::closeHalfFace(const int a_dir, const int a_side, const int a_half) noexcept
{
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);
  CH_assert(a_side == 0 || a_side == 1);
  CH_assert(a_half == 0 || a_half == 1);
  CH_assert(m_numPolygons == 1);

  const int face    = 2 * a_dir + a_side;
  const int tangent = 1 - a_dir;

  Polygon& polygon = m_polygon[0];

  // a stretch of the face that runs past the midpoint is split there, so that each half is its own segments
  for (int i = 0; i < polygon.m_numVertices; i++) {
    if (polygon.m_segmentFace[i] != face) {
      continue;
    }

    const RealVect a = polygon.m_vertex[i];
    const RealVect b = polygon.m_vertex[(i + 1) % polygon.m_numVertices];

    if (!((a[tangent] < 0.0 && b[tangent] > 0.0) || (a[tangent] > 0.0 && b[tangent] < 0.0))) {
      continue;
    }

    if (polygon.m_numVertices >= s_maxVertices) {
      return false;
    }

    for (int k = polygon.m_numVertices; k > i + 1; k--) {
      polygon.m_vertex[k]      = polygon.m_vertex[k - 1];
      polygon.m_vertexEdge[k]  = polygon.m_vertexEdge[k - 1];
      polygon.m_segmentFace[k] = polygon.m_segmentFace[k - 1];
    }

    RealVect midpoint = a;

    midpoint[tangent] = 0.0;

    polygon.m_vertex[i + 1]      = midpoint;
    polygon.m_vertexEdge[i + 1]  = -1;
    polygon.m_segmentFace[i + 1] = face;

    polygon.m_numVertices++;

    i++;
  }

  for (int i = 0; i < polygon.m_numVertices; i++) {
    if (polygon.m_segmentFace[i] != face) {
      continue;
    }

    const RealVect& a = polygon.m_vertex[i];
    const RealVect& b = polygon.m_vertex[(i + 1) % polygon.m_numVertices];

    const Real middle = 0.5 * (a[tangent] + b[tangent]);

    if ((a_half == 0 && middle < 0.0) || (a_half == 1 && middle > 0.0)) {
      polygon.m_segmentFace[i] = -1;
    }
  }

  this->accumulateMoments();

  const bool closed  = this->closureResidual() <= 1.0E-9;
  const bool inRange = m_volumeFraction >= -1.0E-12 && m_volumeFraction <= 1.0 + 1.0E-12;

  return closed && inRange;
}
#endif

int
CutCellBody::numPolygons() const noexcept
{
  return m_numPolygons;
}

int
CutCellBody::widestPolygon() const noexcept
{
  int widest = 0;

  for (int ip = 0; ip < m_numPolygons; ip++) {
    widest = std::max(widest, m_polygon[ip].m_numVertices);
  }

  return widest;
}

void
CutCellBody::accumulateMoments() noexcept
{
  CH_assert(m_numPolygons >= 0 && m_numPolygons <= s_maxPolygons);

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
  CH_assert(m_numPolygons == 1);

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

      CH_assert(face >= -1 && face < s_numFaces);

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

    CH_assert(polygon.m_face >= -1 && polygon.m_face < s_numFaces);

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

    CH_assert(std::abs(m_normal.vectorLength() - 1.0) <= 1.0E-10);
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

  CH_assert(polygon.m_numVertices <= s_maxVertices);

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

bool
CutCellBody::define(const CutCellSurface&       a_surface,
                    const CutCellFaceOverrides& a_overrides,
                    const IntVect&              a_cell) noexcept
{
  if (CutCellBody::classify(a_surface) == Kind::Regular) {
    this->defineWhole();
  }
  else if (!this->define(a_surface)) {
    return false;
  }

  const int entry = a_overrides.find(a_cell);

  if (entry < 0) {
    return true;
  }

  int faceBegin = 0;
  int faceEnd   = 0;

  a_overrides.faces(entry, faceBegin, faceEnd);

#if CH_SPACEDIM == 3
  bool restricted = false;

  for (int f = faceBegin; f < faceEnd; f++) {
    if (a_overrides.reason(f) != CutCellFaceOverrides::s_finer) {
      continue;
    }

    if (!this->replaceFace(a_overrides, f)) {
      return false;
    }

    restricted = true;
  }

  if (restricted && !this->closeInterface()) {
    return false;
  }
#else
  for (int f = faceBegin; f < faceEnd; f++) {
    const int reason = a_overrides.reason(f);

    if (reason != CutCellFaceOverrides::s_closedLowHalf && reason != CutCellFaceOverrides::s_closedHighHalf) {
      continue;
    }

    const int face = a_overrides.face(f);

    if (!this->closeHalfFace(face / 2, face % 2, (reason == CutCellFaceOverrides::s_closedLowHalf) ? 0 : 1)) {
      return false;
    }
  }
#endif

  for (int f = faceBegin; f < faceEnd; f++) {
    if (a_overrides.reason(f) != CutCellFaceOverrides::s_closed) {
      continue;
    }

    const int face = a_overrides.face(f);

    if (!this->snapFace(face / 2, face % 2, false)) {
      return false;
    }
  }

  return true;
}

#if CH_SPACEDIM == 3
bool
CutCellBody::subdivide(CutCellBody* a_children) const noexcept
{
  CH_assert(a_children != nullptr);
  CH_assert(m_numPolygons >= 0 && m_numPolygons <= s_maxPolygons);

  constexpr int numChildren = 1 << SpaceDim;

  for (int c = 0; c < numChildren; c++) {
    a_children[c] = CutCellBody();
  }

  if (m_numPolygons == 0) {
    return true;
  }

  for (int c = 0; c < numChildren; c++) {
    CutCellBody& child = a_children[c];

    child.m_numPolygons = m_numPolygons;

    for (int ip = 0; ip < m_numPolygons; ip++) {
      child.m_polygon[ip] = m_polygon[ip];
    }

    // This body, cut down to the child's quadrant one direction at a time. Cutting leaves it open along the
    // plane, and the patch closing that opening lies in a face of the child -- the one between it and the
    // sibling across the cut, which is the child's high face where the child is low.
    for (int dir = 0; dir < SpaceDim; dir++) {
      const int  side     = (c >> dir) & 1;
      const bool keepHigh = (side == 1);
      const int  cutFace  = 2 * dir + (1 - side);

      const auto inside = [&](const RealVect& a_x) -> bool {
        return keepHigh ? (a_x[dir] >= 0.0) : (a_x[dir] <= 0.0);
      };

      int kept = 0;

      for (int ip = 0; ip < child.m_numPolygons; ip++) {
        const Polygon& in = child.m_polygon[ip];

        Polygon out;
        out.m_numVertices = 0;
        out.m_face        = in.m_face;

        for (int i = 0; i < in.m_numVertices; i++) {
          const RealVect& a = in.m_vertex[i];
          const RealVect& b = in.m_vertex[(i + 1) % in.m_numVertices];

          const bool aIn = inside(a);
          const bool bIn = inside(b);

          RealVect add[2];

          int numAdd = 0;

          if (aIn) {
            add[numAdd++] = a;
          }

          if (aIn != bIn) {
            add[numAdd++] = a + (b - a) * (a[dir] / (a[dir] - b[dir]));
          }

          for (int k = 0; k < numAdd; k++) {
            // A vertex on the cut is inside for both children, so a crossing there repeats the vertex it came
            // from; the repeat carries no length and is dropped rather than leaving a zero-length edge, which
            // the closure walk would read as an opening.
            if (out.m_numVertices > 0 && detail::sameVertex(out.m_vertex[out.m_numVertices - 1], add[k])) {
              continue;
            }

            if (out.m_numVertices >= s_maxVertices) {
              return false;
            }

            out.m_vertexEdge[out.m_numVertices]  = -1;
            out.m_segmentFace[out.m_numVertices] = -1;
            out.m_vertex[out.m_numVertices++]    = add[k];
          }
        }

        // the first and last can meet the same way round the circuit
        while (out.m_numVertices > 1 && detail::sameVertex(out.m_vertex[0], out.m_vertex[out.m_numVertices - 1])) {
          out.m_numVertices--;
        }

        if (out.m_numVertices >= 3) {
          child.m_polygon[kept++] = out;
        }
      }

      child.m_numPolygons = kept;

      if (child.m_numPolygons == 0) {
        break;
      }

      if (!child.closeBoundary(cutFace)) {
        return false;
      }
    }

    if (child.m_numPolygons == 0) {
      child = CutCellBody();

      continue;
    }

    // Into the child's own frame: its centre sits a quarter of a cell from this one's, and a length here is
    // half of one there.
    RealVect centre;

    for (int d = 0; d < SpaceDim; d++) {
      centre[d] = (((c >> d) & 1) == 0) ? -0.25 : 0.25;
    }

    for (int ip = 0; ip < child.m_numPolygons; ip++) {
      Polygon& polygon = child.m_polygon[ip];

      for (int i = 0; i < polygon.m_numVertices; i++) {
        polygon.m_vertex[i] = 2.0 * (polygon.m_vertex[i] - centre);
      }
    }

    child.accumulateMoments();
  }

  return this->partitions(a_children);
}

#endif

#if CH_SPACEDIM == 2
bool
CutCellBody::subdivide(CutCellBody* a_children) const noexcept
{
  CH_assert(a_children != nullptr);
  CH_assert(m_numPolygons >= 0 && m_numPolygons <= s_maxPolygons);

  constexpr int numChildren = 1 << SpaceDim;

  for (int c = 0; c < numChildren; c++) {
    a_children[c] = CutCellBody();
  }

  // Nothing to cut: every child is as this cell is, and a cell holding no fluid leaves them empty.
  if (m_numPolygons == 0) {
    return true;
  }

  for (int c = 0; c < numChildren; c++) {
    Polygon polygon = m_polygon[0];

    // The fluid of this cell, cut down to the child's quadrant one direction at a time. Clipping a closed
    // polygon against a half-plane closes it again along the cut, so the segment the cut leaves is already
    // part of the answer and only has to be told which face it lies in: the face between this child and the
    // sibling on the other side of the cut, which is the child's high face where the child is low.
    bool alive = true;

    for (int dir = 0; dir < SpaceDim && alive; dir++) {
      const int  side     = (c >> dir) & 1;
      const bool keepHigh = (side == 1);
      const int  cutFace  = 2 * dir + (1 - side);

      Polygon out;
      out.m_numVertices = 0;

      const auto inside = [&](const RealVect& a_x) -> bool {
        return keepHigh ? (a_x[dir] >= 0.0) : (a_x[dir] <= 0.0);
      };

      const auto emit = [&](const RealVect& a_x, const int a_face) -> bool {
        // A vertex on the cut is inside for both children, so an entry or exit there repeats the vertex it
        // came from; the repeat carries no length and is dropped rather than closing a zero-area segment.
        if (out.m_numVertices > 0 && detail::sameVertex(out.m_vertex[out.m_numVertices - 1], a_x)) {
          out.m_segmentFace[out.m_numVertices - 1] = a_face;

          return true;
        }

        if (out.m_numVertices >= s_maxVertices) {
          return false;
        }

        out.m_vertexEdge[out.m_numVertices]  = -1;
        out.m_segmentFace[out.m_numVertices] = a_face;
        out.m_vertex[out.m_numVertices++]    = a_x;

        return true;
      };

      for (int i = 0; i < polygon.m_numVertices && alive; i++) {
        const RealVect& a     = polygon.m_vertex[i];
        const RealVect& b     = polygon.m_vertex[(i + 1) % polygon.m_numVertices];
        const int       tagAB = polygon.m_segmentFace[i];

        const bool aIn = inside(a);
        const bool bIn = inside(b);

        if (aIn && bIn) {
          alive = emit(a, tagAB);
        }
        else if (aIn) {
          // leaving: the segment from here runs along the cut until the fluid comes back
          alive = emit(a, tagAB) && emit(a + (b - a) * (a[dir] / (a[dir] - b[dir])), cutFace);
        }
        else if (bIn) {
          // returning: the segment from the crossing to b is the stretch of a -> b that survived
          alive = emit(a + (b - a) * (a[dir] / (a[dir] - b[dir])), tagAB);
        }
      }

      if (!alive) {
        return false;
      }

      // The cut may only enter and leave once. More than one stretch of it means the fluid of this cell meets
      // the child's quadrant in pieces that do not touch, which one polygon cannot describe.
      int runs = 0;

      for (int i = 0; i < out.m_numVertices; i++) {
        const int previous = out.m_segmentFace[(i + out.m_numVertices - 1) % out.m_numVertices];

        if (out.m_segmentFace[i] == cutFace && previous != cutFace) {
          runs++;
        }
      }

      if (runs > 1) {
        return false;
      }

      if (out.m_numVertices < 3) {
        out.m_numVertices = 0;
      }

      out.m_face = -1;
      polygon    = out;
    }

    if (polygon.m_numVertices < 3) {
      continue;
    }

    // Into the child's own frame: its centre sits a quarter of a cell from this one's, and a length here is
    // half of one there.
    RealVect centre;

    for (int d = 0; d < SpaceDim; d++) {
      centre[d] = (((c >> d) & 1) == 0) ? -0.25 : 0.25;
    }

    for (int i = 0; i < polygon.m_numVertices; i++) {
      polygon.m_vertex[i] = 2.0 * (polygon.m_vertex[i] - centre);
    }

    a_children[c].m_polygon[0]  = polygon;
    a_children[c].m_numPolygons = 1;

    a_children[c].accumulateMoments();
  }

  return this->partitions(a_children);
}
#endif

bool
CutCellBody::isConnected() const noexcept
{
  CH_assert(m_numPolygons >= 0 && m_numPolygons <= s_maxPolygons);

  if (m_numPolygons <= 1) {
    return true;
  }

  // Whether two polygons carry the same edge, which on a closed surface they do in opposite directions. The
  // welding closeBoundary does first is what makes the two sides agree segment for segment, so a shared edge
  // is one pair of vertices and not a stretch of one against several of the other.
  const auto adjacent = [&](const int a_ip, const int a_jp) -> bool {
    const Polygon& p = m_polygon[a_ip];
    const Polygon& q = m_polygon[a_jp];

    for (int i = 0; i < p.m_numVertices; i++) {
      const RealVect& a = p.m_vertex[i];
      const RealVect& b = p.m_vertex[(i + 1) % p.m_numVertices];

      for (int j = 0; j < q.m_numVertices; j++) {
        const RealVect& c = q.m_vertex[j];
        const RealVect& d = q.m_vertex[(j + 1) % q.m_numVertices];

        if (detail::sameVertex(a, d) && detail::sameVertex(b, c)) {
          return true;
        }
      }
    }

    return false;
  };

  bool reached[s_maxPolygons] = {false};
  int  pending[s_maxPolygons];

  int top   = 0;
  int found = 1;

  reached[0]     = true;
  pending[top++] = 0;

  while (top > 0) {
    const int ip = pending[--top];

    for (int jp = 0; jp < m_numPolygons; jp++) {
      if (reached[jp] || !adjacent(ip, jp)) {
        continue;
      }

      reached[jp]    = true;
      pending[top++] = jp;

      found++;
    }
  }

  return found == m_numPolygons;
}

bool
CutCellBody::interfaceIsPlanar() const noexcept
{
  CH_assert(m_numPolygons >= 0 && m_numPolygons <= s_maxPolygons);

#if CH_SPACEDIM == 2
  // In two dimensions the fluid region is one polygon whose segments lie in different cell faces, and the
  // interface is those segments belonging to none of them. They share a line exactly when every endpoint of
  // every one of them lies on the line of the first, which is the two-dimensional reading of one plane.
  if (m_numPolygons == 0) {
    return true;
  }

  const Polygon& polygon = m_polygon[0];

  int first = -1;

  for (int i = 0; i < polygon.m_numVertices && first < 0; i++) {
    if (polygon.m_segmentFace[i] < 0) {
      first = i;
    }
  }

  // A body with no interface is one the surface does not enter, whose fluid is the whole cell.
  if (first < 0) {
    return true;
  }

  const RealVect& base      = polygon.m_vertex[first];
  const RealVect  direction = polygon.m_vertex[(first + 1) % polygon.m_numVertices] - base;

  const Real length = direction.vectorLength();

  // A first segment of no length gives no line to measure the others against, so nothing is claimed.
  if (length <= s_edgeTolerance) {
    return false;
  }

  const RealVect unit = direction / length;

  for (int i = 0; i < polygon.m_numVertices; i++) {
    if (polygon.m_segmentFace[i] >= 0) {
      continue;
    }

    for (int k = 0; k < 2; k++) {
      const RealVect& x = polygon.m_vertex[(i + k) % polygon.m_numVertices];

      const RealVect offset = x - base;

      if (std::abs(offset[0] * unit[1] - offset[1] * unit[0]) > s_edgeTolerance) {
        return false;
      }
    }
  }

  return true;
#else
  int first = -1;

  for (int ip = 0; ip < m_numPolygons && first < 0; ip++) {
    if (m_polygon[ip].m_face < 0) {
      first = ip;
    }
  }

  // A body with no interface is one the surface does not enter, whose fluid is the whole cell.
  if (first < 0) {
    return true;
  }

  Real     area = 0.0;
  RealVect vector;
  RealVect centroid;

  detail::polygonMoments(m_polygon[first].m_vertex, m_polygon[first].m_numVertices, area, vector, centroid);

  // A first patch of no area gives no plane to measure the others against, so nothing is claimed.
  if (area <= s_nullArea) {
    return false;
  }

  const RealVect normal = vector / area;

  for (int ip = first + 1; ip < m_numPolygons; ip++) {
    const Polygon& p = m_polygon[ip];

    if (p.m_face >= 0) {
      continue;
    }

    for (int iv = 0; iv < p.m_numVertices; iv++) {
      if (std::abs((p.m_vertex[iv] - centroid).dotProduct(normal)) > s_edgeTolerance) {
        return false;
      }
    }
  }

  return true;
#endif
}

bool
CutCellBody::hasMultiValuedChildren(const int a_refRat) const noexcept
{
  CH_assert(a_refRat >= 2);
  CH_assert((a_refRat & (a_refRat - 1)) == 0);
  CH_assert(m_numPolygons >= 0 && m_numPolygons <= s_maxPolygons);

  constexpr int numChildren = 1 << SpaceDim;

  if (this->interfaceIsPlanar()) {
    return false;
  }

  CutCellBody children[numChildren];

  // The children are built before the moments are checked, so they can be asked even when the cut is refused
  // for a reason of its own. A child the cut never reached is left empty, which reads as connected.
  const bool cut = this->subdivide(children);

  for (int c = 0; c < numChildren; c++) {
    if (!children[c].isConnected()) {
      return true;
    }
  }

  if (!cut || a_refRat == 2) {
    return false;
  }

  for (int c = 0; c < numChildren; c++) {
    if (children[c].volumeFraction() <= 0.0) {
      continue;
    }

    if (children[c].hasMultiValuedChildren(a_refRat / 2)) {
      return true;
    }
  }

  return false;
}

bool
CutCellBody::partitions(const CutCellBody* a_children) const noexcept
{
  CH_assert(a_children != nullptr);

  constexpr int  numChildren = 1 << SpaceDim;
  constexpr int  numShared   = 1 << (SpaceDim - 1);
  constexpr Real volumeScale = 1.0 / static_cast<Real>(numChildren);
  constexpr Real areaScale   = 1.0 / static_cast<Real>(numShared);

  constexpr Real tolerance = 1.0E-12;

  Real volume = 0.0;

  for (int c = 0; c < numChildren; c++) {
    // Refining a cell whose fluid is in one piece must not leave a child whose fluid is in two. If it does, the
    // child is a cell this generator cannot describe, and it was produced rather than encountered, so it is a
    // fault here rather than a geometry to be refused.
    // A singly cut cell can genuinely refine into a multi-cut one: the fluid region is not convex once the
    // interface has a crease, and a non-convex region can meet an octant in two pieces. It is reported rather
    // than refused here, because hasMultiValuedChildren asks exactly this question and has to get an answer.
    if (!a_children[c].isConnected()) {
      return false;
    }

    if (a_children[c].divergenceResidual() > tolerance) {
      return false;
    }

    volume += a_children[c].volumeFraction();
  }

  if (std::abs(volumeScale * volume - m_volumeFraction) > tolerance) {
    return false;
  }

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (int side = 0; side < 2; side++) {
      Real aperture = 0.0;

      for (int c = 0; c < numChildren; c++) {
        if (((c >> dir) & 1) == side) {
          aperture += a_children[c].areaFraction(dir, (side == 0) ? Side::Lo : Side::Hi);
        }
      }

      if (std::abs(areaScale * aperture - m_areaFraction[2 * dir + side]) > tolerance) {
        return false;
      }
    }

    // The face between two children is the high one of the low child and the low one of the high child.
    for (int c = 0; c < numChildren; c++) {
      if (((c >> dir) & 1) != 0) {
        continue;
      }

      const Real mine   = a_children[c].areaFraction(dir, Side::Hi);
      const Real theirs = a_children[c | (1 << dir)].areaFraction(dir, Side::Lo);

      if (std::abs(mine - theirs) > tolerance) {
        return false;
      }
    }
  }

  RealVect area = RealVect::Zero;

  for (int c = 0; c < numChildren; c++) {
    area += a_children[c].boundaryArea() * a_children[c].normal();
  }

  if ((areaScale * area - m_boundaryArea * m_normal).vectorLength() > tolerance) {
    return false;
  }

  return true;
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
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  return m_areaFraction[2 * a_dir + ((a_side == Side::Lo) ? 0 : 1)];
}

const RealVect&
CutCellBody::faceCentroid(const int a_dir, const Side::LoHiSide a_side) const noexcept
{
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

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

bool
CutCellBody::identical(const CutCellBody& a_other) const noexcept
{
  const auto sameVector = [](const RealVect& a_first, const RealVect& a_second) -> bool {
    for (int d = 0; d < SpaceDim; d++) {
      if (a_first[d] != a_second[d]) {
        return false;
      }
    }

    return true;
  };

  if (m_numPolygons != a_other.m_numPolygons) {
    return false;
  }

  for (int ip = 0; ip < m_numPolygons; ip++) {
    const Polygon& mine   = m_polygon[ip];
    const Polygon& theirs = a_other.m_polygon[ip];

    if (mine.m_face != theirs.m_face || mine.m_numVertices != theirs.m_numVertices) {
      return false;
    }

    for (int iv = 0; iv < mine.m_numVertices; iv++) {
      if (!sameVector(mine.m_vertex[iv], theirs.m_vertex[iv]) || mine.m_vertexEdge[iv] != theirs.m_vertexEdge[iv]) {
        return false;
      }
    }
  }

  for (int face = 0; face < s_numFaces; face++) {
    if (m_areaFraction[face] != a_other.m_areaFraction[face] ||
        !sameVector(m_faceCentroid[face], a_other.m_faceCentroid[face])) {
      return false;
    }
  }

  return m_volumeFraction == a_other.m_volumeFraction && m_boundaryArea == a_other.m_boundaryArea &&
         m_trueBoundaryArea == a_other.m_trueBoundaryArea && sameVector(m_volumeCentroid, a_other.m_volumeCentroid) &&
         sameVector(m_normal, a_other.m_normal) && sameVector(m_boundaryCentroid, a_other.m_boundaryCentroid) &&
         sameVector(m_closure, a_other.m_closure);
}

std::uint64_t
CutCellBody::fingerprint() const noexcept
{
  // FNV-1a over the bytes of every field identical compares, in the same order
  std::uint64_t hash = 14695981039346656037ULL;

  const auto mix = [&hash](const void* a_data, const std::size_t a_bytes) {
    const unsigned char* bytes = static_cast<const unsigned char*>(a_data);

    for (std::size_t i = 0; i < a_bytes; i++) {
      hash ^= bytes[i];
      hash *= 1099511628211ULL;
    }
  };

  mix(&m_numPolygons, sizeof(int));

  for (int ip = 0; ip < m_numPolygons; ip++) {
    const Polygon& polygon = m_polygon[ip];

    mix(&polygon.m_face, sizeof(int));
    mix(&polygon.m_numVertices, sizeof(int));

    for (int iv = 0; iv < polygon.m_numVertices; iv++) {
      mix(&polygon.m_vertex[iv], sizeof(RealVect));
      mix(&polygon.m_vertexEdge[iv], sizeof(int));
    }
  }

  mix(m_areaFraction, sizeof(m_areaFraction));
  mix(m_faceCentroid, sizeof(m_faceCentroid));
  mix(&m_volumeFraction, sizeof(Real));
  mix(&m_boundaryArea, sizeof(Real));
  mix(&m_trueBoundaryArea, sizeof(Real));
  mix(&m_volumeCentroid, sizeof(RealVect));
  mix(&m_normal, sizeof(RealVect));
  mix(&m_boundaryCentroid, sizeof(RealVect));
  mix(&m_closure, sizeof(RealVect));

  return hash;
}

} // namespace PolyhedralEB

#include <CD_NamespaceFooter.H>
