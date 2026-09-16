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

  for (int c = 1; c < CutCellSurface::s_numCorners; c++) {
    if (isFluid(a_surface.m_corner[c]) != firstFluid) {
      return Kind::Cut;
    }
  }

  return firstFluid ? Kind::Regular : Kind::Covered;
}

#if CH_SPACEDIM == 3
void
CutCellBody::orientOutward(Polygon& a_polygon, const int a_dir, const int a_side) const noexcept
{
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);
  CH_assert(a_side == 0 || a_side == 1);

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

      // a chord vertex on the boundary between two children is a vertex of the cells on the other
      // side of the seam, and dropping it would leave their two segments meeting the middle of one
      // of ours: watertight, but not a shared edge. The children sit at plus and minus a quarter, so
      // their boundaries in the face are at exactly zero. Only chord vertices are kept: a vertex on
      // the face's own boundary has to go, or the edge it splits no longer matches the neighbouring
      // face's whole one and closeInterface reads the face boundary as open.
      bool onChildBoundary = false;
      bool onFaceBoundary  = false;

      for (int d = 0; d < SpaceDim; d++) {
        if (d != a_in[0].m_face / 2) {
          onChildBoundary = onChildBoundary || (std::abs(here[d]) <= detail::s_weldTolerance);
          onFaceBoundary  = onFaceBoundary || (std::abs(std::abs(here[d]) - 0.5) <= detail::s_weldTolerance);
        }
      }

      if (straight && !(onChildBoundary && !onFaceBoundary)) {
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
CutCellBody::restrictFace(const CutCellSurface* a_children, const int a_dir, const int a_side) noexcept
{
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

bool
CutCellBody::closeInterface() noexcept
{

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
  // a body that is not cut holds no interface polygon, so the loop appends nothing for it
  for (int ip = 0; ip < m_numPolygons; ip++) {
    const Polygon& p = m_polygon[ip];

    if (p.m_face >= 0 || p.m_numVertices < 3) {
      continue;
    }

    // the interface is already fanned, but a polygon carrying more than three vertices is fanned
    // again here rather than left for the reader to triangulate
    for (int v = 1; v + 1 < p.m_numVertices; v++) {
      const RealVect* corner[3] = {&p.m_vertex[0], &p.m_vertex[v], &p.m_vertex[v + 1]};

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

} // namespace PolyhedralEB

#include <CD_NamespaceFooter.H>
