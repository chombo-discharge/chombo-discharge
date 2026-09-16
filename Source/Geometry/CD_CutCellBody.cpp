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
  m_kind             = Kind::Covered;
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
CutCellBody::coarsenKind(const Kind a_kinds[1 << SpaceDim]) noexcept
{
  constexpr int numChildren = 1 << SpaceDim;

  bool anyCut     = false;
  bool anyRegular = false;
  bool anyCovered = false;

  for (int c = 0; c < numChildren; c++) {
    switch (a_kinds[c]) {
    case Kind::Cut: {
      anyCut = true;

      break;
    }
    case Kind::Regular: {
      anyRegular = true;

      break;
    }
    default: {
      anyCovered = true;

      break;
    }
    }
  }

  if (anyCut || (anyRegular && anyCovered)) {
    return Kind::Cut;
  }

  return anyRegular ? Kind::Regular : Kind::Covered;
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
CutCellBody::interfaceLiesInFace(const CutCellSurface& a_surface) noexcept
{
  for (int d = 0; d < SpaceDim; d++) {
    for (int side = 0; side < 2; side++) {
      int faceCorner[1 << (SpaceDim - 1)];

      detail::faceCorners(d, side, faceCorner);

      // Every corner of the face is on the solid side. In a cell classify has called regular
      // that can only be a face the interface lies in, since the only solid corners such a cell
      // has are ones touchesOnly found on the interface. A cell whose corners are all fluid has
      // no such face, so the test costs those cells the loop and nothing else.
      bool allSolid = true;

      for (int k = 0; k < (1 << (SpaceDim - 1)); k++) {
        allSolid = allSolid && !isFluid(a_surface.m_corner[faceCorner[k]]);
      }

      if (allSolid) {
        return true;
      }
    }
  }

  return false;
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
    bool reachesOtherSide = false;

    for (int d = 0; d < SpaceDim; d++) {
      const int other = c ^ (1 << d);

      if (isFluid(a_surface.m_corner[other]) == a_fluidSide) {
        continue;
      }

      reachesOtherSide = true;

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

    // A corner every one of whose neighbours is on its own side has no edge to read a crossing
    // from, and nothing above has tested it. It is enclosed by its own side rather than sitting
    // on the interface, so the side it belongs to encloses something and this is a cut cell. The
    // shape the rule is written for, a facet lying in a node plane, leaves every corner on that
    // side facing one across the cell, so it is unaffected.
    if (!reachesOtherSide) {
      return false;
    }
  }

  return any;
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

void
CutCellBody::appendInterfaceFacets(Vector<Real>& a_facets, const RealVect& a_cellCentre, const Real a_dx) const noexcept
{
  if (m_kind != Kind::Cut) {
    return;
  }

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
          a_facets.push_back(a_cellCentre[d] + a_dx * (*corner[k])[d]);
        }
      }
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

      if (straight) {
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

    m_kind = kind;

    return true;
  }

  if (!this->defineCut(a_surface)) {
    return false;
  }

  this->accumulateMoments();

  m_kind = this->kindFromPolygons();

  // verify rather than assume: a folded patch leaves the body open or the volume outside its
  // range, and every individual moment can still look plausible when it does
  const bool closed  = this->closureResidual() <= 1.0E-9;
  const bool inRange = m_volumeFraction >= -1.0E-12 && m_volumeFraction <= 1.0 + 1.0E-12;

  return closed && inRange;
}

bool
CutCellBody::clip(const int a_dir, const Real a_coordinate, const bool a_keepLow, CutCellBody& a_out) const noexcept
{
  a_out = CutCellBody();

#if CH_SPACEDIM == 3
  // Where the cut crosses each polygon it leaves a segment in the cutting plane, and those
  // segments bound the cap the cut creates. Three dimensions only: in two the cut leaves an edge
  // of the clipped polygon itself, tagged with the face it lies in, and there is no separate cap
  // to assemble -- so a walk over these segments would never close, and must not be asked to.
  RealVect segmentFrom[s_maxPolygons];
  RealVect segmentTo[s_maxPolygons];
  int      segmentFace[s_maxPolygons];

  int numSegments = 0;
#endif

  const Real sign = a_keepLow ? 1.0 : -1.0;

  for (int ip = 0; ip < m_numPolygons; ip++) {
    const Polygon& polygon = m_polygon[ip];

    if (a_out.m_numPolygons >= s_maxPolygons) {
      return false;
    }

    Polygon& cut = a_out.m_polygon[a_out.m_numPolygons];

    cut.m_numVertices = 0;
    cut.m_face        = polygon.m_face;

    for (int i = 0; i < polygon.m_numVertices; i++) {
      const RealVect& a = polygon.m_vertex[i];
      const RealVect& b = polygon.m_vertex[(i + 1) % polygon.m_numVertices];

      const Real fa = sign * (a[a_dir] - a_coordinate);
      const Real fb = sign * (b[a_dir] - a_coordinate);

      if (cut.m_numVertices + 2 > s_maxVertices) {
        return false;
      }

      if (fa <= detail::s_clipTolerance) {
        RealVect keep = a;

        if (fa >= -detail::s_clipTolerance) {
          keep[a_dir] = a_coordinate;
        }

        cut.m_segmentFace[cut.m_numVertices] = polygon.m_segmentFace[i];
        cut.m_vertexEdge[cut.m_numVertices]  = polygon.m_vertexEdge[i];
        cut.m_vertex[cut.m_numVertices++]    = keep;
      }

      if ((fa < -detail::s_clipTolerance && fb > detail::s_clipTolerance) ||
          (fb < -detail::s_clipTolerance && fa > detail::s_clipTolerance)) {
        // interpolate from the lexicographically lower end, whichever way this polygon walks
        // the edge, so that the two polygons sharing it land on the same point bit for bit
        const bool      ordered = detail::lexLess(a, b);
        const RealVect& first   = ordered ? a : b;
        const RealVect& second  = ordered ? b : a;
        const Real      f0      = ordered ? fa : fb;
        const Real      f1      = ordered ? fb : fa;

        RealVect x = first + (f0 / (f0 - f1)) * (second - first);

        x[a_dir] = a_coordinate;

        // Leaving the half-space, the segment starting here runs along the cut and so belongs
        // to the face the cut creates. Entering it, the segment continues along the edge it
        // came from and keeps that edge's face. Only the two-dimensional moments read this,
        // where the cut leaves an edge of the clipped polygon rather than a separate cap.
        const bool leaving = fa <= detail::s_clipTolerance;

        cut.m_segmentFace[cut.m_numVertices] = leaving ? (2 * a_dir + (a_keepLow ? 1 : 0)) : polygon.m_segmentFace[i];
        cut.m_vertexEdge[cut.m_numVertices]  = -1;
        cut.m_vertex[cut.m_numVertices++]    = x;
      }
    }

    if (cut.m_numVertices < 3) {
      continue;
    }

    if (detail::polygonArea(cut.m_vertex, cut.m_numVertices) < s_nullArea) {
      continue;
    }

    a_out.m_numPolygons++;

#if CH_SPACEDIM == 3
    bool whollyInPlane = true;

    for (int i = 0; i < cut.m_numVertices; i++) {
      whollyInPlane = whollyInPlane && (std::abs(cut.m_vertex[i][a_dir] - a_coordinate) < detail::s_clipTolerance);
    }

    if (whollyInPlane) {
      continue;
    }

    for (int i = 0; i < cut.m_numVertices; i++) {
      const RealVect& a = cut.m_vertex[i];
      const RealVect& b = cut.m_vertex[(i + 1) % cut.m_numVertices];

      if (numSegments >= s_maxPolygons) {
        return false;
      }

      const bool bothInPlane = std::abs(a[a_dir] - a_coordinate) < detail::s_clipTolerance &&
                               std::abs(b[a_dir] - a_coordinate) < detail::s_clipTolerance;

      if (bothInPlane && !detail::sameVertex(a, b)) {
        segmentFrom[numSegments] = b;
        segmentTo[numSegments]   = a;
        segmentFace[numSegments] = 2 * a_dir + (a_keepLow ? 1 : 0);
        numSegments++;
      }
    }
#endif
  }

#if CH_SPACEDIM == 3
  bool used[s_maxPolygons] = {false};

  for (int s0 = 0; s0 < numSegments; s0++) {
    if (used[s0]) {
      continue;
    }

    used[s0] = true;

    if (a_out.m_numPolygons >= s_maxPolygons) {
      return false;
    }

    Polygon& loop = a_out.m_polygon[a_out.m_numPolygons];

    loop.m_numVertices = 0;
    loop.m_face        = segmentFace[s0];

    loop.m_segmentFace[loop.m_numVertices] = segmentFace[s0];
    loop.m_vertexEdge[loop.m_numVertices]  = -1;
    loop.m_vertex[loop.m_numVertices++]    = segmentFrom[s0];

    RealVect        current = segmentTo[s0];
    const RealVect& end     = segmentFrom[s0];

    bool closed = false;

    for (int guard = 0; guard <= numSegments + 1; guard++) {
      if (detail::sameVertex(current, end)) {
        closed = true;

        break;
      }

      // the nearest segment starting where this one ended, rather than the first within reach,
      // so that a vertex two loops pass close to does not splice them together
      int  next = -1;
      Real best = detail::s_weldTolerance;

      for (int j = 0; j < numSegments; j++) {
        if (used[j]) {
          continue;
        }

        const Real distance = (segmentFrom[j] - current).vectorLength();

        if (distance <= best) {
          best = distance;
          next = j;
        }
      }

      if (next < 0 || loop.m_numVertices >= s_maxVertices) {
        break;
      }

      used[next] = true;

      loop.m_segmentFace[loop.m_numVertices] = segmentFace[s0];
      loop.m_vertexEdge[loop.m_numVertices]  = -1;
      loop.m_vertex[loop.m_numVertices++]    = current;

      current = segmentTo[next];
    }

    // A cap that does not close is a face the body will be left without, and a body missing a
    // face still has moments and still looks like an answer. Refusing here is what turns that
    // into a cell the caller can decline rather than one it cannot tell from a good one.
    if (!closed) {
      return false;
    }

    if (loop.m_numVertices < 3) {
      continue;
    }

    if (detail::polygonArea(loop.m_vertex, loop.m_numVertices) >= s_nullArea) {
      a_out.m_numPolygons++;
    }
  }
#endif

  return true;
}

bool
CutCellBody::refine(CutCellBody a_children[1 << SpaceDim]) const noexcept
{
  constexpr int numChildren = 1 << SpaceDim;

  // Cut as a tree. Splitting on the first direction gives two bodies rather than 2^SpaceDim,
  // and splitting those on the second gives four, so each plane is applied once per body it
  // actually divides instead of once per child. The last direction writes the children
  // directly, so the working buffers only ever hold half of them.
  CutCellBody buffer[2][numChildren / 2];

  buffer[0][0] = *this;

  int source = 0;
  int count  = 1;

  for (int d = 0; d < SpaceDim; d++) {
    const bool last = (d == SpaceDim - 1);

    for (int i = 0; i < count; i++) {
      for (int side = 0; side < 2; side++) {
        // the child's index carries the side it took in bit d, which is the same convention the
        // frame shift below reads
        const int    destination = i + side * count;
        CutCellBody& target      = last ? a_children[destination] : buffer[1 - source][destination];

        if (!buffer[source][i].clip(d, 0.0, side == 0, target)) {
          return false;
        }
      }
    }

    source = 1 - source;
    count *= 2;
  }

  for (int c = 0; c < numChildren; c++) {
    CutCellBody& child = a_children[c];

    for (int ip = 0; ip < child.m_numPolygons; ip++) {
      for (int iv = 0; iv < child.m_polygon[ip].m_numVertices; iv++) {
        for (int d = 0; d < SpaceDim; d++) {
          const Real centre = -0.25 + 0.5 * static_cast<Real>((c >> d) & 1);

          child.m_polygon[ip].m_vertex[iv][d] = 2.0 * (child.m_polygon[ip].m_vertex[iv][d] - centre);
        }
      }
    }

    child.accumulateMoments();

    child.m_kind = child.kindFromPolygons();
  }

  return true;
}

bool
CutCellBody::coarsen(const CutCellBody a_children[1 << SpaceDim]) noexcept
{
  constexpr int  numChildren = 1 << SpaceDim;
  constexpr Real volumeScale = 1.0 / static_cast<Real>(numChildren);
  constexpr Real areaScale   = 2.0 * volumeScale;

  *this = CutCellBody();

  // The cell is cut if the surface passes through any of the cells partitioning it, whatever its
  // own corners say. This is the whole point of building it from them: a feature finer than the
  // cell leaves every corner on one side, and reading the corners would call the cell whole.
  Kind kinds[numChildren];

  for (int c = 0; c < numChildren; c++) {
    kinds[c] = a_children[c].kind();
  }

  m_kind = CutCellBody::coarsenKind(kinds);

  RealVect volumeMoment   = RealVect::Zero;
  RealVect boundaryVector = RealVect::Zero;
  RealVect boundaryMoment = RealVect::Zero;

  Real     faceArea[s_numFaces] = {0.0};
  RealVect faceMoment[s_numFaces];

  for (int f = 0; f < s_numFaces; f++) {
    faceMoment[f] = RealVect::Zero;
  }

  for (int c = 0; c < numChildren; c++) {
    const CutCellBody& child = a_children[c];

    // Where the child sits in this cell's frame, and how a point in the child's own frame maps
    // into it.
    RealVect origin;

    for (int d = 0; d < SpaceDim; d++) {
      origin[d] = -0.25 + 0.5 * static_cast<Real>((c >> d) & 1);
    }

    const Real volume = volumeScale * child.m_volumeFraction;

    m_volumeFraction += volume;
    volumeMoment += volume * (origin + 0.5 * child.m_volumeCentroid);

    const RealVect vector = -child.m_boundaryArea * child.m_normal;
    const Real     patch  = areaScale * child.m_trueBoundaryArea;

    boundaryVector += areaScale * vector;
    m_trueBoundaryArea += patch;
    boundaryMoment += patch * (origin + 0.5 * child.m_boundaryCentroid);

    m_closure += areaScale * child.m_closure;

    // Only the children lying against a face of this cell carry any of its aperture.
    for (int d = 0; d < SpaceDim; d++) {
      for (int side = 0; side < 2; side++) {
        if (((c >> d) & 1) != side) {
          continue;
        }

        const int  f    = 2 * d + side;
        const Real area = areaScale * child.m_areaFraction[f];

        faceArea[f] += area;
        faceMoment[f] += area * (origin + 0.5 * child.m_faceCentroid[f]);
      }
    }
  }

  if (m_volumeFraction > 0.0) {
    m_volumeCentroid = volumeMoment / m_volumeFraction;
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

  // The interface is taken from the apertures rather than from the children's own interfaces,
  // which is what the divergence identity says it is: what leaves through the faces has to arrive
  // through the boundary. Summing the children's interfaces instead only agrees when the faces
  // they share cancel, and a child the corner values call covered sitting against one that calls
  // the face between them half open does not cancel. The identity is exact whether or not the
  // children agree, so it is the thing to build on.
  RealVect apertureVector = RealVect::Zero;

  for (int d = 0; d < SpaceDim; d++) {
    apertureVector[d] = m_areaFraction[2 * d + 1] - m_areaFraction[2 * d];
  }

  m_boundaryArea = apertureVector.vectorLength();

  if (m_trueBoundaryArea > 0.0) {
    m_boundaryCentroid = boundaryMoment / m_trueBoundaryArea;
  }

  if (m_boundaryArea > 0.0) {
    m_normal = apertureVector / m_boundaryArea;
  }
  else if (boundaryVector.vectorLength() > 0.0) {
    // Nothing leaves through the faces, so the identity fixes neither the area nor the direction.
    // The children's interfaces still say which way it faces, and interfaceIsOneSided will refuse
    // the cell on the finest level, where such a thing is not allowed to stand.
    m_normal = -boundaryVector / boundaryVector.vectorLength();
  }

  return this->divergenceResidual() <= s_nullArea + s_edgeTolerance;
}

bool
CutCellBody::interfaceIsOneSided() const noexcept
{
  if (m_trueBoundaryArea <= s_nullArea) {
    return true;
  }

  return m_trueBoundaryArea <= s_maxAreaSpread * m_boundaryArea;
}

CutCellBody::Kind
CutCellBody::kind() const noexcept
{
  return m_kind;
}

CutCellBody::Kind
CutCellBody::kindFromPolygons() const noexcept
{
  if (m_numPolygons == 0) {
    return Kind::Covered;
  }

  for (int i = 0; i < m_numPolygons; i++) {
    if (m_polygon[i].m_face < 0) {
      return Kind::Cut;
    }
  }

  return Kind::Regular;
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
