/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// Exports the per-cut-cell EB moments, together with the true edge crossings of the implicit
// function, so that the cut-cell refinement census can be run offline.
//
// The census asks whether the edge crossings that bound the interface patch inside a cut cell can
// be recovered from the stored moments alone -- the face apertures and face centroids -- without
// calling the implicit function. If they can, a refinement operator built on a triangulation
// through those crossings is a pure function of the stored data.

#include <fstream>
#include <iomanip>

#include <CD_Driver.H>
#include <CD_RoughSphere.H>
#include <CD_Tessellation.H>
#include <PlaneIF.H>
#include <IntersectionIF.H>
#include <CD_TorusSdf.H>
#include <CD_SphereSdf.H>
#include <CD_CylinderSdf.H>
#include <CD_Electrode.H>
#include <CD_GeometryStepper.H>
#include <CD_PolyhedralEBUtils.H>
#include <CD_CutCellBody.H>
#include <BRMeshRefine.H>
#include <LoadBalance.H>
#include <CD_PolyhedralGeometryShop.H>

using namespace ChomboDischarge;
using namespace Physics::Geometry;

// A convex body as an intersection of half-spaces: the zero set is exactly the polyhedron's
// surface, so its edges and corners are sharp at whatever orientation the normals are given
// in. Positive inside the solid and negative in the fluid, matching the other geometries.
class ConvexBody : public BaseIF
{
public:
  ConvexBody(const Vector<RealVect>& a_normals, const Vector<Real>& a_offsets, const RealVect& a_center)
    : m_normals(a_normals), m_offsets(a_offsets), m_center(a_center)
  {}

  ConvexBody(const ConvexBody& a_other)
    : m_normals(a_other.m_normals), m_offsets(a_other.m_offsets), m_center(a_other.m_center)
  {}

  virtual Real
  value(const RealVect& a_point) const override
  {
    Real f = -std::numeric_limits<Real>::max();

    for (int i = 0; i < m_normals.size(); i++) {
      Real d = 0.0;
      for (int k = 0; k < SpaceDim; k++) {
        d += m_normals[i][k] * (a_point[k] - m_center[k]);
      }

      f = std::max(f, d - m_offsets[i]);
    }

    return -f;
  }

  virtual BaseIF*
  newImplicitFunction() const override
  {
    return static_cast<BaseIF*>(new ConvexBody(*this));
  }

protected:
  Vector<RealVect> m_normals;
  Vector<Real>     m_offsets;
  RealVect         m_center;
};

// Rotation by the given Euler angles, in degrees. In 2D only the first angle is used.
inline RealVect
rotate(const RealVect& a_v, const RealVect& a_deg)
{
  RealVect v = a_v;

#if CH_SPACEDIM == 2
  const Real c = cos(a_deg[0] * M_PI / 180.0);
  const Real s = sin(a_deg[0] * M_PI / 180.0);
  const Real x = v[0];
  v[0]         = c * x - s * v[1];
  v[1]         = s * x + c * v[1];
#else
  for (int axis = 0; axis < 3; axis++) {
    const Real c = cos(a_deg[axis] * M_PI / 180.0);
    const Real s = sin(a_deg[axis] * M_PI / 180.0);
    const int  i = (axis + 1) % 3;
    const int  j = (axis + 2) % 3;
    const Real a = v[i];
    v[i]         = c * a - s * v[j];
    v[j]         = s * a + c * v[j];
  }
#endif

  return v;
}

// Cell edges: one per direction, per corner of the transverse directions. Four in 2D, twelve in 3D.
constexpr int numCellEdges = SpaceDim * (1 << (SpaceDim - 1));

// Transverse directions of a_dir, in increasing order.
inline void
transverseDirs(const int a_dir, int a_trans[SpaceDim - 1])
{
  int j = 0;
  for (int d = 0; d < SpaceDim; d++) {
    if (d != a_dir) {
      a_trans[j++] = d;
    }
  }
}

// Lower node of cell edge a_edge, expressed as an offset in {0,1}^SpaceDim from the cell's lower
// corner, together with the direction the edge runs along.
inline void
edgeGeometry(const int a_edge, int& a_dir, IntVect& a_loOffset)
{
  a_dir = a_edge / (1 << (SpaceDim - 1));

  const int corner = a_edge % (1 << (SpaceDim - 1));

  int trans[SpaceDim - 1];
  transverseDirs(a_dir, trans);

  a_loOffset = IntVect::Zero;
  for (int j = 0; j < SpaceDim - 1; j++) {
    a_loOffset[trans[j]] = (corner >> j) & 1;
  }
}

// Bisect the implicit function along an edge. Returns the crossing as a fraction of the edge
// length, or a negative value when the endpoints do not bracket a crossing.
inline Real
edgeCrossing(const BaseIF& a_implicitFunction, const RealVect& a_loPt, const RealVect& a_hiPt)
{
  const Real fLo = a_implicitFunction.value(a_loPt);
  const Real fHi = a_implicitFunction.value(a_hiPt);

  if (fLo * fHi > 0.0) {
    return -1.0;
  }

  Real tLo = 0.0;
  Real tHi = 1.0;

  for (int iter = 0; iter < 60; iter++) {
    const Real     tMid = 0.5 * (tLo + tHi);
    const RealVect xMid = a_loPt + tMid * (a_hiPt - a_loPt);
    const Real     fMid = a_implicitFunction.value(xMid);

    if (fMid * fLo <= 0.0) {
      tHi = tMid;
    }
    else {
      tLo = tMid;
    }
  }

  return 0.5 * (tLo + tHi);
}

// Pull a point onto the zero level set along the implicit function's gradient. A cell holding a
// sharp edge has its crease inside the cell, where the edge crossings cannot see it; anchoring
// the interface patch on the surface itself puts the fan apex on or near that crease instead of
// floating off it.
inline RealVect
projectToSurface(const BaseIF& a_implicitFunction, const RealVect& a_point, const Real a_dx)
{
  RealVect   x = a_point;
  const Real h = 1.0E-6 * a_dx;

  for (int iter = 0; iter < 4; iter++) {
    const Real f = a_implicitFunction.value(x);

    RealVect grad;
    for (int d = 0; d < SpaceDim; d++) {
      RealVect xp = x;
      RealVect xm = x;
      xp[d] += h;
      xm[d] -= h;

      grad[d] = (a_implicitFunction.value(xp) - a_implicitFunction.value(xm)) / (2.0 * h);
    }

    const Real g2 = grad.dotProduct(grad);
    if (g2 <= 0.0) {
      break;
    }

    x -= (f / g2) * grad;
  }

  return x;
}

// Trace the interface across one cell face and report where it bends.
//
// The two endpoints are the face's edge crossings, which are shared. Between them the contour
// is sampled by pushing points off the straight chord onto the zero set along the in-face
// perpendicular. Samples that stand far enough off the chord are kept as contour vertices, in
// order from the first endpoint to the second. Those vertices lie in the INTERIOR of the face,
// so no edge crossing can see them, and a chord-only description throws them away: one bend is
// an edge crossing the face, two is a corner where three surfaces meet.
inline int
faceBends(const BaseIF&   a_implicitFunction,
          const RealVect& a_p0,
          const RealVect& a_p1,
          const int       a_faceDir,
          const Real      a_dx,
          RealVect        a_bend[2])
{
  RealVect   chord = a_p1 - a_p0;
  const Real len   = chord.vectorLength();
  if (len <= 0.0) {
    return 0;
  }
  chord /= len;

  RealVect faceNormal   = RealVect::Zero;
  faceNormal[a_faceDir] = 1.0;

  RealVect   perp = PolyGeom::cross(faceNormal, chord);
  const Real pl   = perp.vectorLength();
  if (pl <= 0.0) {
    return 0;
  }
  perp /= pl;

  constexpr int numSamples = 15;

  Real     offset[numSamples];
  RealVect point[numSamples];

  for (int k = 0; k < numSamples; k++) {
    const Real     u = Real(k + 1) / Real(numSamples + 1);
    const RealVect q = a_p0 + (u * len) * chord;

    offset[k] = 0.0;
    point[k]  = q;

    Real lo  = -0.5 * a_dx;
    Real hi  = 0.5 * a_dx;
    Real flo = a_implicitFunction.value(q + lo * perp);
    Real fhi = a_implicitFunction.value(q + hi * perp);

    if (flo * fhi > 0.0) {
      continue;
    }

    for (int iter = 0; iter < 50; iter++) {
      const Real mid = 0.5 * (lo + hi);
      const Real fm  = a_implicitFunction.value(q + mid * perp);

      if (fm * flo <= 0.0) {
        hi  = mid;
        fhi = fm;
      }
      else {
        lo  = mid;
        flo = fm;
      }
    }

    offset[k] = 0.5 * (lo + hi);
    point[k]  = q + offset[k] * perp;
  }

  // a smooth arc bows by O(dx^2/R) and is well served by the chord; a bend is not
  const Real tol = 1.0E-3 * a_dx;

  int first = -1;
  for (int k = 0; k < numSamples; k++) {
    if (std::abs(offset[k]) > tol && (first < 0 || std::abs(offset[k]) > std::abs(offset[first]))) {
      first = k;
    }
  }

  if (first < 0) {
    return 0;
  }

  // the largest remaining deviation on either side of the first is the second vertex
  int second = -1;
  for (int k = 0; k < numSamples; k++) {
    if (std::abs(k - first) < 3) {
      continue;
    }
    if (std::abs(offset[k]) > tol && (second < 0 || std::abs(offset[k]) > std::abs(offset[second]))) {
      second = k;
    }
  }

  int n = 0;
  if (second >= 0 && second < first) {
    a_bend[n++] = point[second];
    a_bend[n++] = point[first];
  }
  else {
    a_bend[n++] = point[first];
    if (second >= 0) {
      a_bend[n++] = point[second];
    }
  }

  return n;
}

// Reconstruct a cut cell's surface from the implicit function alone: corner values, and the
// position of the sign change along each edge that has one. No moments and no interior, which is
// all a refinement criterion needs.
inline PolyhedralEB::CutCellSurface
sampleSurface(const BaseIF&   a_implicitFunction,
              const IntVect&  a_cell,
              const RealVect& a_probLo,
              const Real      a_dx,
              const bool      a_root)
{
  PolyhedralEB::CutCellSurface surface;

  RealVect corner[PolyhedralEB::CutCellSurface::s_numCorners];

  for (int c = 0; c < PolyhedralEB::CutCellSurface::s_numCorners; c++) {
    RealVect x = a_probLo;

    for (int d = 0; d < SpaceDim; d++) {
      x[d] += a_dx * static_cast<Real>(a_cell[d] + ((c >> d) & 1));
    }

    corner[c]           = x;
    surface.m_corner[c] = a_implicitFunction.value(x);
  }

  for (int e = 0; e < PolyhedralEB::CutCellSurface::s_numEdges; e++) {
    int lo = 0;
    int hi = 0;

    PolyhedralEB::detail::edgeCorners(e, lo, hi);

    if (PolyhedralEB::isFluid(surface.m_corner[lo]) != PolyhedralEB::isFluid(surface.m_corner[hi])) {
      surface.m_crossing[e] = a_root ? edgeCrossing(a_implicitFunction, corner[lo], corner[hi])
                                     : surface.m_corner[lo] / (surface.m_corner[lo] - surface.m_corner[hi]);
    }
  }

  return surface;
}

// Sample the implicit function on the lattice the curvature estimator expects, centred on a_point.
inline void
curvatureStencil(const BaseIF&   a_implicitFunction,
                 const RealVect& a_point,
                 const Real      a_spacing,
                 Real            a_stencil[PolyhedralEB::s_curvatureStencilSize])
{
  for (int i = 0; i < PolyhedralEB::s_curvatureStencilSize; i++) {
    RealVect x = a_point;

    int packed = i;

    for (int d = 0; d < SpaceDim; d++) {
      x[d] += a_spacing * static_cast<Real>((packed % 3) - 1);
      packed /= 3;
    }

    a_stencil[i] = a_implicitFunction.value(x);
  }
}

// Run the curvature pre-pass and write the cells it tags, so that the tags it produces from the
// implicit function alone can be set against the ones the embedded boundary produces.
void
exportCurvatureTags(const RefCountedPtr<ComputationalGeometry>& a_compgeom,
                    const RefCountedPtr<AmrMesh>&               a_amr,
                    const Real                                  a_angle,
                    const int                                   a_growth,
                    const int                                   a_maxDepth,
                    const std::string&                          a_fileName)
{
  // The pre-pass goes as deep as it is asked to, which is not tied to how deep this run refines.
  Vector<int> refRatios = a_amr->getRefinementRatios();

  while (static_cast<int>(refRatios.size()) < a_maxDepth) {
    refRatios.push_back(2);
  }

  Vector<Vector<Box>> regions;

  const Vector<IntVectSet> tags = a_compgeom->getCurvatureTags(regions,
                                                               a_amr->getDomains()[0],
                                                               refRatios,
                                                               a_amr->getBlockingFactor() * IntVect::Unit,
                                                               a_amr->getMaxBoxSize() * IntVect::Unit,
                                                               a_amr->getProbLo(),
                                                               a_amr->getDx()[0],
                                                               a_angle,
                                                               a_growth,
                                                               a_maxDepth);

  std::ofstream out(a_fileName);

  out << "level";
  for (int d = 0; d < SpaceDim; d++) {
    out << ",iv" << d;
  }
  out << "\n";

  for (int lvl = 0; lvl < tags.size(); lvl++) {
    for (IVSIterator ivsIt(tags[lvl]); ivsIt.ok(); ++ivsIt) {
      const IntVect iv = ivsIt();

      out << lvl;
      for (int d = 0; d < SpaceDim; d++) {
        out << "," << iv[d];
      }
      out << "\n";
    }
  }
}

// Measure the seam that a partially refined index space would carry.
//
// Where a fine level exists, the coarse cells under it come through coarsening and agree with
// their children by construction. Where it does not, the coarse cells are reconstructed at their
// own resolution. Along the boundary between the two, a coarse cell built one way sits next to
// one built the other, and a face they share carries a single stored aperture. This reports how
// far apart the two constructions are: for every cut cell, the moments taken from the coarse
// reconstruction against the same moments summed from eight independently reconstructed children.
// Build a cell's body from the cells that partition it, reconstructed a level finer, and those
// from theirs, down to a_depth. At depth zero this is just the reconstruction on the cell itself.
// This is the operation a coarse level performs where a finer one covers it: its geometry comes
// up through coarsening rather than being reconstructed a second time at its own resolution.
bool
buildAggregated(const BaseIF&              a_implicitFunction,
                PolyhedralEB::CutCellBody& a_body,
                const IntVect&             a_cell,
                const RealVect&            a_probLo,
                const Real                 a_dx,
                const int                  a_depth)
{
  if (a_depth <= 0) {
    return a_body.define(sampleSurface(a_implicitFunction, a_cell, a_probLo, a_dx, true));
  }

  constexpr int numChildren = 1 << SpaceDim;

  PolyhedralEB::CutCellBody child[numChildren];

  for (int c = 0; c < numChildren; c++) {
    IntVect fine = 2 * a_cell;

    for (int d = 0; d < SpaceDim; d++) {
      fine[d] += (c >> d) & 1;
    }

    if (!buildAggregated(a_implicitFunction, child[c], fine, a_probLo, 0.5 * a_dx, a_depth - 1)) {
      return false;
    }
  }

  return a_body.coarsen(child);
}

// Total up the geometry the whole surface carries, at a range of aggregation depths, so that the
// figures can be set against the ones the geometry is known to have.
// Look for cells holding more than one sheet of interface.
//
// A cell built from the cells partitioning it can end up with two patches of interface facing
// opposite ways -- a feature thinner than the cell passing through it leaves one on each side.
// EBData stores one plane per cell, so such a cell cannot be described: the area vectors cancel,
// leaving the magnitude of their sum near zero while the interface is really there. Worse, a body
// whose interface is wholly inside it touches none of its faces, so the divergence identity reads
// zero equals zero and the body is accepted. This reports what the moments actually come to.
// Find where building a cell from the cells that partition it leaves the divergence identity
// unsatisfied, and show the parent against its children so the face that fails to cancel is
// visible. Two siblings share an internal face and it should contribute nothing to the parent,
// which only holds if they agree about it.
// Check the invariant a partially covered index space rests on: every grid a run uses has to lie
// inside what the index space actually holds at that level. Nothing enforces it -- the grids come
// from tags and the coverage from the geometry -- so it is worth asking rather than assuming.
// Cuts a level the generator never built, and checks what comes out against the level it was cut
// from.
//
// The shop is driven directly rather than through the index space. The surfaces it keeps live on
// the layouts it makes for itself, and until EBISLevel learns to ask for a level it does not
// hold, going in the front door would only reach the levels that were generated.
//
// Three things are asked of the result. The volume of the fine cells under a coarse one has to
// come to the coarse one's, which is the property the whole approach rests on: a child is a piece
// of its parent, not a fresh reconstruction that happens to sit inside it. Two fine cells sharing
// a face have to agree on its area. And the fill has to answer at all -- a refusal is counted
// rather than ignored, because a box the generator cannot cut is a box the index space would be
// left without.
void
validateRefinedFill(const RefCountedPtr<ComputationalGeometry>& a_compgeom,
                    const RefCountedPtr<AmrMesh>&               a_amr,
                    const int                                   a_depth,
                    const int                                   a_coarseBoxSize)
{
  CH_TIME("validateRefinedFill");

  if (a_depth <= 0) {
    return;
  }

  const RefCountedPtr<BaseIF>& implicitFunction = a_compgeom->getGasImplicitFunction();

  if (implicitFunction.isNull()) {
    return;
  }

  const RealVect probLo  = a_amr->getProbLo();
  const int      ebGhost = 1;

  // The level that is generated, and the finer one cut from it. The shop coarsens down from the
  // finest domain it is handed, so the finest of its levels is the one being cut, not the one
  // being generated.
  const ProblemDomain coarseDomain = a_amr->getDomains()[0];
  const Real          coarseDx     = a_amr->getDx()[0];

  const int refinement = 1 << a_depth;

  ProblemDomain fineDomain = coarseDomain;
  fineDomain.refine(refinement);

  const Real fineDx = coarseDx / static_cast<Real>(refinement);

  PolyhedralGeometryShop shop(*implicitFunction, 0, fineDx, probLo, fineDomain, fineDomain, ebGhost, 0.0, true, 1);

  Vector<Box> boxes;
  domainSplit(coarseDomain, boxes, a_coarseBoxSize, 1);

  Vector<int> ranks;
  LoadBalance(ranks, boxes);

  DisjointBoxLayout dbl(boxes, ranks, coarseDomain);

  static_cast<GeometryService&>(shop).postMakeBoxLayout(dbl, coarseDx * RealVect::Unit);

  // The generated-at-its-own-resolution comparison below fills the graph at the fine spacing, and
  // the generator records a surface for every cut cell it makes, so the fine level's store has to
  // exist before it is asked to.
  DisjointBoxLayout fineDbl;
  refine(fineDbl, dbl, refinement);

  static_cast<GeometryService&>(shop).postMakeBoxLayout(fineDbl, fineDx * RealVect::Unit);

  LayoutData<Vector<IrregNode>> coarseNodes(dbl);
  LayoutData<BaseFab<int>>      coarseKinds(dbl);

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
    const Box valid = dbl[dit()];

    Box ghost = grow(valid, ebGhost);
    ghost &= coarseDomain.domainBox();

    shop.fillGraph(coarseKinds[dit()], coarseNodes[dit()], valid, ghost, coarseDomain, probLo, coarseDx, dit());
  }

  // The parents are handed to the fill rather than looked up by it, since in the index space they
  // live on the level and not on the generator. They are assembled here the same way EBISLevel
  // assembles them, so that this exercises the interface the index space will use.
  const int numComponents = shop.numSurfaceComponents();

  LayoutData<RefCountedPtr<EBGraph>>         coarseGraphs(dbl);
  LayoutData<RefCountedPtr<BaseIVFAB<Real>>> coarseSurfaces(dbl);

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
    const Box valid = dbl[dit()];

    Box ghost = grow(valid, ebGhost);
    ghost &= coarseDomain.domainBox();

    // A BaseIVFAB asks its graph only how many volumes each of its cells holds, and the generator
    // mandates single-valued cut cells, so an all-regular graph over the region gives the same one
    // volume per cell that the index space's own graph does. Building the real graph here would
    // cover the valid region alone -- buildGraph redefines itself to the region it is given --
    // and the parents wanted reach past it.
    coarseGraphs[dit()] = RefCountedPtr<EBGraph>(new EBGraph(ghost));
    coarseGraphs[dit()]->setToAllRegular();

    // Reaching past the box, since cutting a fine box needs the parents of its ghost cells too.
    // The index space does this with a copyTo; here the boxes this rank holds are gathered by
    // hand, which is why a box whose ghost parents sit on another rank is refused below.
    Vector<IntVect> cells;
    Vector<Real>    values;

    for (DataIterator source = dbl.dataIterator(); source.ok(); ++source) {
      Vector<IntVect> theirCells;
      Vector<Real>    theirValues;

      shop.getSurfaces(theirCells, theirValues, ghost & dbl[source()], coarseDx);

      for (int n = 0; n < theirCells.size(); n++) {
        cells.push_back(theirCells[n]);

        for (int comp = 0; comp < numComponents; comp++) {
          values.push_back(theirValues[n * numComponents + comp]);
        }
      }
    }

    IntVectSet ivs;

    for (int n = 0; n < cells.size(); n++) {
      ivs |= cells[n];
    }

    coarseSurfaces[dit()] = RefCountedPtr<BaseIVFAB<Real>>(
      new BaseIVFAB<Real>(ivs, *coarseGraphs[dit()], numComponents));

    for (int n = 0; n < cells.size(); n++) {
      const VolIndex vof(cells[n], 0);

      for (int comp = 0; comp < numComponents; comp++) {
        (*coarseSurfaces[dit()])(vof, comp) = values[n * numComponents + comp];
      }
    }
  }

  long long refused         = 0;
  long long parents         = 0;
  long long compared        = 0;
  long long eitherCut       = 0;
  long long volumeDiffers   = 0;
  long long topologyDiffers = 0;
  Real      worstCutVsGen   = 0.0;
  long long coarseCut       = 0;
  long long volumeBad       = 0;
  long long faceMismatch    = 0;

  Real worstVolume = 0.0;
  Real worstFace   = 0.0;

  const Real tolerance = 1.0E-12;

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit) {
    const Box coarseBox = dbl[dit()];
    const Box fineBox   = refine(coarseBox, refinement);

    // What the coarse level says each of its cells holds, which is what the fine cells under it
    // have to come to.
    BaseFab<Real> coarseVolume(coarseBox, 1);

    for (BoxIterator bit(coarseBox); bit.ok(); ++bit) {
      coarseVolume(bit(), 0) = (coarseKinds[dit()](bit(), 0) > 0) ? 1.0 : 0.0;
    }

    const Vector<IrregNode>& nodes = coarseNodes[dit()];

    for (int n = 0; n < nodes.size(); n++) {
      if (coarseBox.contains(nodes[n].m_cell)) {
        coarseVolume(nodes[n].m_cell, 0) = nodes[n].m_volFrac;
      }
    }

    // The fine cells, as the fill hands them back. Apertures are kept for the whole refinement of
    // the coarse box rather than one call at a time, so that faces two separate calls produced
    // are checked against each other as well.
    BaseFab<Real> fineVolume(fineBox, 1);
    BaseFab<Real> fineAperture(fineBox, 2 * SpaceDim);

    fineVolume.setVal(-1.0);
    fineAperture.setVal(-1.0);

    // one call per parent, which is what the driver is meant to do: cut a parent once and keep
    // the whole child set rather than cutting a path per output cell
    for (BoxIterator parentIt(coarseBox); parentIt.ok(); ++parentIt) {
      Box chunk(parentIt(), parentIt());
      chunk.refine(refinement);

      Box chunkGhost = grow(chunk, ebGhost);
      chunkGhost &= fineDomain.domainBox();

      BaseFab<int>      kinds;
      Vector<IrregNode> fineNodes;

      parents++;

      if (!shop.fillRefinedGraph(kinds,
                                 fineNodes,
                                 chunk,
                                 chunkGhost,
                                 fineDomain,
                                 probLo,
                                 fineDx,
                                 *coarseSurfaces[dit()],
                                 coarseDx)) {
        refused++;

        continue;
      }

      for (BoxIterator bit(chunk); bit.ok(); ++bit) {
        const Real full = (kinds(bit(), 0) > 0) ? 1.0 : 0.0;

        fineVolume(bit(), 0) = full;

        for (int c = 0; c < 2 * SpaceDim; c++) {
          fineAperture(bit(), c) = full;
        }
      }

      for (int n = 0; n < fineNodes.size(); n++) {
        const IntVect& iv = fineNodes[n].m_cell;

        if (!chunk.contains(iv)) {
          continue;
        }

        fineVolume(iv, 0) = fineNodes[n].m_volFrac;

        for (int dir = 0; dir < SpaceDim; dir++) {
          for (SideIterator sit; sit.ok(); ++sit) {
            const int index = fineNodes[n].index(dir, sit());

            Real aperture = 0.0;

            for (int f = 0; f < fineNodes[n].m_areaFrac[index].size(); f++) {
              aperture += fineNodes[n].m_areaFrac[index][f];
            }

            fineAperture(iv, 2 * dir + ((sit() == Side::Lo) ? 0 : 1)) = aperture;
          }
        }
      }
    }

    // The same cells, generated at their own resolution instead of cut from their parent. If the
    // two agree then it does not matter which way a refined layout is filled; if they do not, the
    // difference is the information coarsening threw away.
    {
      const Box fineValid = fineBox;

      Box fineGhost = grow(fineValid, 1);
      fineGhost &= fineDomain.domainBox();

      BaseFab<int>      genKinds;
      Vector<IrregNode> genNodes;

      shop.fillGraph(genKinds, genNodes, fineValid, fineGhost, fineDomain, probLo, fineDx, dit());

      BaseFab<Real> genVolume(fineValid, 1);

      for (BoxIterator bit(fineValid); bit.ok(); ++bit) {
        genVolume(bit(), 0) = (genKinds(bit(), 0) > 0) ? 1.0 : 0.0;
      }

      for (int n = 0; n < genNodes.size(); n++) {
        if (fineValid.contains(genNodes[n].m_cell)) {
          genVolume(genNodes[n].m_cell, 0) = genNodes[n].m_volFrac;
        }
      }

      for (BoxIterator bit(fineValid); bit.ok(); ++bit) {
        if (fineVolume(bit(), 0) < 0.0) {
          continue;
        }

        const Real cut = fineVolume(bit(), 0);
        const Real gen = genVolume(bit(), 0);

        const bool cutIsCut = (cut > 0.0 && cut < 1.0);
        const bool genIsCut = (gen > 0.0 && gen < 1.0);

        compared++;

        if (!cutIsCut && !genIsCut) {
          continue;
        }

        eitherCut++;

        if (cutIsCut != genIsCut) {
          topologyDiffers++;
        }

        const Real d = std::abs(cut - gen);

        worstCutVsGen = std::max(worstCutVsGen, d);

        if (d > 1.0E-12) {
          volumeDiffers++;
        }
      }
    }

    // the volume of the cells under a coarse one has to come to the coarse one's
    for (BoxIterator bit(coarseBox); bit.ok(); ++bit) {
      Box under(bit(), bit());
      under.refine(refinement);

      Real sum = 0.0;
      bool got = true;

      for (BoxIterator sub(under); sub.ok(); ++sub) {
        if (fineVolume(sub(), 0) < 0.0) {
          got = false;

          break;
        }

        sum += fineVolume(sub(), 0);
      }

      if (!got) {
        continue;
      }

      if (coarseKinds[dit()](bit(), 0) == 0) {
        coarseCut++;
      }

      const Real mismatch = std::abs(sum / static_cast<Real>(under.numPts()) - coarseVolume(bit(), 0));

      worstVolume = std::max(worstVolume, mismatch);

      if (mismatch > tolerance) {
        volumeBad++;
      }
    }

    // two fine cells sharing a face have to agree on its area
    for (int dir = 0; dir < SpaceDim; dir++) {
      Box interior = fineBox;
      interior.growHi(dir, -1);

      for (BoxIterator bit(interior); bit.ok(); ++bit) {
        const IntVect lo = bit();
        const IntVect hi = lo + BASISV(dir);

        if (fineAperture(lo, 2 * dir + 1) < 0.0 || fineAperture(hi, 2 * dir) < 0.0) {
          continue;
        }

        const Real difference = std::abs(fineAperture(lo, 2 * dir + 1) - fineAperture(hi, 2 * dir));

        worstFace = std::max(worstFace, difference);

        if (difference > tolerance) {
          faceMismatch++;
        }
      }
    }
  }

  pout() << "CUTvsGEN depth " << a_depth << ": cutCells " << eitherCut << " volumeDiffers " << volumeDiffers
         << " topologyDiffers " << topologyDiffers << " worst " << worstCutVsGen << endl;

  pout() << "REFINEFILL depth " << a_depth << " dx " << fineDx << ": parents " << parents << " refused " << refused
         << " coarseCutCells " << coarseCut << " volumeBad " << volumeBad << " worstVolume " << worstVolume
         << " faceMismatch " << faceMismatch << " worstFace " << worstFace << endl;
}

// How many vertices a multichord face polygon needs.
//
// A chord-carrying cell whose neighbour is coarsened takes that face's chords from the neighbour's
// children instead of drawing one of its own. The face is a square; under one halving it is a 3x3
// node grid with twelve sub-edges, and the fluid region on it is bounded by the nodes that lie in
// the fluid together with the crossings on those sub-edges. That count is what the polygon has to
// hold, and it is what decides whether the representation can carry a multichord seam.
//
// Counted here rather than built, because the count is all that is in question: the cap is on
// vertices, and a walk of the boundary would visit exactly these.
void
validateMultichordFace(const RefCountedPtr<ComputationalGeometry>& a_compgeom,
                       const RefCountedPtr<AmrMesh>&               a_amr,
                       const int                                   a_numCells)
{
  CH_TIME("validateMultichordFace");

  const RefCountedPtr<BaseIF>& implicitFunction = a_compgeom->getGasImplicitFunction();

  if (implicitFunction.isNull()) {
    return;
  }

  const RealVect probLo = a_amr->getProbLo();
  const Real     dx     = a_amr->getDx()[0];

  long long faces      = 0;
  int       worstOne   = 0;
  int       worstMulti = 0;
  long long over14     = 0;
  long long over20     = 0;

  Box slab = a_amr->getDomains()[0].domainBox();

  for (int d = 0; d < SpaceDim; d++) {
    if (slab.size(d) > a_numCells) {
      slab.setBig(d, slab.smallEnd(d) + a_numCells - 1);
    }
  }

  for (BoxIterator bit(slab); bit.ok(); ++bit) {
    const IntVect iv = bit();

    RealVect lo = probLo;

    for (int d = 0; d < SpaceDim; d++) {
      lo[d] += dx * iv[d];
    }

    for (int dir = 0; dir < SpaceDim; dir++) {
      for (int side = 0; side < 2; side++) {
        int trans[SpaceDim - 1];
        transverseDirs(dir, trans);

        // the nine nodes of the halved face, and the values there
        Real value[3][3];

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            RealVect x = lo;
            x[dir] += side * dx;
            x[trans[0]] += 0.5 * dx * i;
#if CH_SPACEDIM == 3
            x[trans[1]] += 0.5 * dx * j;
#endif
            value[i][j] = implicitFunction->value(x);
          }
        }

        // is this face cut at all, at the coarse spacing?
        const bool coarseCut = !((value[0][0] < 0.0) == (value[2][0] < 0.0) &&
                                 (value[0][0] < 0.0) == (value[0][2] < 0.0) &&
                                 (value[0][0] < 0.0) == (value[2][2] < 0.0));

        bool anyFine = false;

        for (int i = 0; i < 3 && !anyFine; i++) {
          for (int j = 0; j < 3 && !anyFine; j++) {
            anyFine = anyFine || ((value[i][j] < 0.0) != (value[0][0] < 0.0));
          }
        }

        if (!coarseCut && !anyFine) {
          continue;
        }

        faces++;

        // one chord: the four corners in the fluid, plus the crossings on the four coarse edges
        int one = 0;

        for (int i = 0; i < 3; i += 2) {
          for (int j = 0; j < 3; j += 2) {
            if (value[i][j] < 0.0) {
              one++;
            }
          }
        }

        if ((value[0][0] < 0.0) != (value[2][0] < 0.0))
          one++;
        if ((value[0][2] < 0.0) != (value[2][2] < 0.0))
          one++;
        if ((value[0][0] < 0.0) != (value[0][2] < 0.0))
          one++;
        if ((value[2][0] < 0.0) != (value[2][2] < 0.0))
          one++;

        // multichord: every node in the fluid, plus every crossing on the twelve sub-edges
        int multi = 0;

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            if (value[i][j] < 0.0) {
              multi++;
            }
          }
        }

        for (int i = 0; i < 2; i++) {
          for (int j = 0; j < 3; j++) {
            if ((value[i][j] < 0.0) != (value[i + 1][j] < 0.0)) {
              multi++;
            }
          }
        }

        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 2; j++) {
            if ((value[i][j] < 0.0) != (value[i][j + 1] < 0.0)) {
              multi++;
            }
          }
        }

        worstOne   = std::max(worstOne, one);
        worstMulti = std::max(worstMulti, multi);

        if (multi > 14) {
          over14++;
        }
        if (multi > 20) {
          over20++;
        }
      }
    }
  }

  pout() << "MULTIFACE faces " << faces << " worstSingleChord " << worstOne << " worstMultichord " << worstMulti
         << " over14 " << over14 << " over20 " << over20 << endl;
}

// Build one cell's surface straight from the implicit function, in its own [-1/2,1/2] frame.
static void
makeSurface(const BaseIF&                 a_implicitFunction,
            const RealVect&               a_cellCentre,
            const Real                    a_dx,
            PolyhedralEB::CutCellSurface& a_surface)
{
  for (int c = 0; c < PolyhedralEB::CutCellSurface::s_numCorners; c++) {
    RealVect x = a_cellCentre;

    for (int d = 0; d < SpaceDim; d++) {
      x[d] += a_dx * (((c >> d) & 1) - 0.5);
    }

    a_surface.m_corner[c] = a_implicitFunction.value(x);
  }

  for (int e = 0; e < PolyhedralEB::CutCellSurface::s_numEdges; e++) {
    int     dir;
    IntVect loOffset;
    edgeGeometry(e, dir, loOffset);

    RealVect lo = a_cellCentre;

    for (int d = 0; d < SpaceDim; d++) {
      lo[d] += a_dx * (loOffset[d] - 0.5);
    }

    RealVect hi = lo;
    hi[dir] += a_dx;

    // Bracket and bisect on the same predicate the body classifies corners with. A test on the
    // product of the endpoint values disagrees with it wherever a value is a negative zero, and
    // the surface is then internally inconsistent: a corner the body calls fluid sits on an edge
    // this says never crosses.
    const Real fLo = a_implicitFunction.value(lo);
    const Real fHi = a_implicitFunction.value(hi);

    if (PolyhedralEB::isFluid(fLo) == PolyhedralEB::isFluid(fHi)) {
      a_surface.m_crossing[e] = PolyhedralEB::CutCellSurface::s_noCrossing;

      continue;
    }

    Real tLo = 0.0;
    Real tHi = 1.0;

    for (int iter = 0; iter < 60 && (tHi - tLo) > 1.0E-15; iter++) {
      const Real     tMid = 0.5 * (tLo + tHi);
      const RealVect xMid = lo + tMid * (hi - lo);

      if (PolyhedralEB::isFluid(a_implicitFunction.value(xMid)) == PolyhedralEB::isFluid(fLo)) {
        tLo = tMid;
      }
      else {
        tHi = tMid;
      }
    }

    a_surface.m_crossing[e] = 0.5 * (tLo + tHi);
  }
}

// A cell one of whose faces carries the chords of the cells on the other side of it, instead of a
// single chord of its own.
//
// The face polygons for that face are not assembled by a wider face walk -- faceWalk is built
// around a face having four edges and four corners, and generalising it means redoing the
// four-crossing saddle logic, which is the most error-prone thing in this whole representation.
// They are the face polygons of the cells abutting the face, obtained from the same faceWalk at
// their own spacing and mapped into this cell's frame. The rest of the body is untouched.
class SeamBody : public PolyhedralEB::CutCellBody
{
public:
  // Returns false if the assembly does not fit. a_gap is the closure residual, which is what the
  // interface loop would have to absorb: the body's other faces and its interface still describe
  // the single chord this face no longer has.
  // Attach the interface to whatever the face polygons leave open.
  //
  // In a closed body every edge is shared by exactly two polygons, traversed once each way. The
  // face polygons share their cell-edge segments with each other; every edge they do not share is
  // where the interface has to attach. So the patch can be read off the faces rather than built
  // from the cell's crossing loop, which is what lets a face carry any number of chords: the
  // construction never asks how many.
  //
  // Orientation comes free. A face traverses its open edge one way, so the interface traverses it
  // the other, which is exactly the rule a closed surface obeys.
  // Merge coplanar polygons that share edges into the single polygon bounding their union.
  //
  // Their shared edges are traversed once each way, so dropping every edge that has a reverse
  // among the others leaves exactly the union's boundary, and walking what remains gives it in
  // order. Orientation is inherited: the surviving edges keep the direction they had.
  //
  // Collinear vertices are then dropped, which is what makes the merged boundary match the
  // neighbouring faces. A sub-face's edge along a cell edge is half the length of the coarse face's
  // edge beside it, so without this the two never compare equal, and every cell edge is mistaken
  // for a chord.
  bool
  mergeCoplanar(const Polygon* a_in, const int a_num, Polygon* a_out, const int a_maxOut, int& a_numOut, int& a_why)
    const
  {
    a_why    = 0;
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

        if (PolyhedralEB::detail::sameVertex(a, b)) {
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

            interior = PolyhedralEB::detail::sameVertex(a, d) && PolyhedralEB::detail::sameVertex(b, c);
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
      a_why = 1;

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
        if (PolyhedralEB::detail::sameVertex(current, end)) {
          closed = true;

          break;
        }

        int next = -1;

        for (int j = 0; j < numEdges && next < 0; j++) {
          if (!used[j] && PolyhedralEB::detail::sameVertex(from[j], current)) {
            next = j;
          }
        }

        if (next < 0 || numWalk >= 4 * s_maxVertices) {
          a_why = 2;

          return false;
        }

        used[next]      = true;
        walk[numWalk++] = current;
        current         = to[next];
      }

      if (!closed) {
        a_why = 3;

        return false;
      }

      if (a_numOut >= a_maxOut) {
        a_why = 7;

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
          a_why = 5;

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
      a_why = 6;

      return false;
    }

    return true;
  }

  // Replace one face's chord with the chords of the four children that cover it. This is what
  // buildSeam does to the seam face, applied to an arbitrary face, and it also drops the interface
  // because the interface was built to meet the chord being replaced.
  bool
  restrictFace(const BaseIF&   a_implicitFunction,
               const RealVect& a_centre,
               const Real      a_dx,
               const int       a_dir,
               const int       a_side,
               int&            a_why,
               int&            a_mergeWhy)
  {
    const int face = 2 * a_dir + a_side;

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
      RealVect childCentre = a_centre;

      for (int d = 0; d < SpaceDim; d++) {
        origin[d] = -0.25 + 0.5 * static_cast<Real>((which >> d) & 1);
        childCentre[d] += 0.5 * a_dx * (((which >> d) & 1) - 0.5);
      }

      PolyhedralEB::CutCellSurface fine;
      makeSurface(a_implicitFunction, childCentre, 0.5 * a_dx, fine);

      Polygon walked[2];

      const int numWalked = this->faceWalk(a_dir, a_side, fine, walked);

      if (numWalked < 0) {
        a_why = 3;

        return false;
      }

      for (int n = 0; n < numWalked; n++) {
        if (numSub >= 4 * (1 << (SpaceDim - 1))) {
          a_why = 4;

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

    if (numSub > 0) {
      Polygon merged[1 << (SpaceDim - 1)];

      int numMerged = 0;

      if (!this->mergeCoplanar(sub, numSub, merged, 1 << (SpaceDim - 1), numMerged, a_mergeWhy)) {
        a_why = 6;

        return false;
      }

      for (int n = 0; n < numMerged; n++) {
        if (m_numPolygons >= s_maxPolygons) {
          a_why = 4;

          return false;
        }

        m_polygon[m_numPolygons++] = merged[n];
      }
    }

    return true;
  }

  // Pass B. Restrict the seam face, and every face that meets it across an edge the coarse cell
  // gets wrong, so that all faces sharing such an edge describe it the same way. Restricting only
  // the seam face is what leaves the two descriptions in contradiction.
  bool
  buildRestricted(const BaseIF&                       a_implicitFunction,
                  const PolyhedralEB::CutCellSurface& a_coarse,
                  const RealVect&                     a_centre,
                  const Real                          a_dx,
                  const int                           a_seamDir,
                  const int                           a_seamSide,
                  const bool*                         a_marked,
                  int&                                a_facesRestricted,
                  int&                                a_why,
                  int&                                a_mergeWhy)
  {
    a_facesRestricted = 0;
    a_why             = 0;
    a_mergeWhy        = 0;

    if (!this->define(a_coarse)) {
      a_why = 1;

      return false;
    }

    if (m_kind != Kind::Cut) {
      a_why = 2;

      return false;
    }

    if (!this->restrictFace(a_implicitFunction, a_centre, a_dx, a_seamDir, a_seamSide, a_why, a_mergeWhy)) {
      return false;
    }

    a_facesRestricted++;

    for (int e = 0; e < PolyhedralEB::CutCellSurface::s_numEdges; e++) {
      if (!a_marked[e]) {
        continue;
      }

      int     dir;
      IntVect loOffset;
      edgeGeometry(e, dir, loOffset);

      // only edges of the seam face have a finer description to restrict from
      if (dir == a_seamDir || loOffset[a_seamDir] != a_seamSide) {
        continue;
      }

      for (int d = 0; d < SpaceDim; d++) {
        if (d == dir || d == a_seamDir) {
          continue;
        }

        if (!this->restrictFace(a_implicitFunction, a_centre, a_dx, d, loOffset[d], a_why, a_mergeWhy)) {
          return false;
        }

        a_facesRestricted++;
      }
    }

    int loops = 0, needFan = 0, needFlat = 0, flatLoops = 0;

    Real bend = 0.0;

    if (!this->closeInterface(loops, needFan, needFlat, flatLoops, bend)) {
      a_why = 5;

      return false;
    }

    this->accumulateMoments();

    return true;
  }

  // Interface segments lying along a cell edge, which the embedded boundary never does except
  // where two of this body's faces disagree about that edge.
  int
  cellEdgeInterfaceSegments() const
  {
    int found = 0;

    for (int ip = 0; ip < m_numPolygons; ip++) {
      if (m_polygon[ip].m_face >= 0) {
        continue;
      }

      for (int i = 0; i < m_polygon[ip].m_numVertices; i++) {
        const RealVect& a = m_polygon[ip].m_vertex[i];
        const RealVect& b = m_polygon[ip].m_vertex[(i + 1) % m_polygon[ip].m_numVertices];

        int fixed = 0;

        for (int d = 0; d < SpaceDim; d++) {
          if (std::abs(a[d] - b[d]) <= 1.0E-12 && std::abs(std::abs(a[d]) - 0.5) <= 1.0E-12) {
            fixed++;
          }
        }

        if (fixed >= 2) {
          found++;
        }
      }
    }

    return found;
  }

  // The face polygons this body holds, in cell-relative coordinates, where a cell edge is where two
  // coordinates are both +/- 0.5.
  void
  dumpPolygons(const char* a_tag) const
  {
    pout() << "DUMP " << a_tag << " polygons " << m_numPolygons << endl;

    for (int ip = 0; ip < m_numPolygons; ip++) {
      pout() << "  poly " << ip << " face " << m_polygon[ip].m_face << " verts " << m_polygon[ip].m_numVertices;

      for (int i = 0; i < m_polygon[ip].m_numVertices; i++) {
        pout() << " (" << m_polygon[ip].m_vertex[i][0] << "," << m_polygon[ip].m_vertex[i][1] << ","
               << m_polygon[ip].m_vertex[i][2] << ")";
      }

      pout() << endl;
    }
  }

  // Write this body's interface patch as STL facets, in physical coordinates.
  //
  // Only the interface is written: the cell-face polygons are not part of the embedded boundary
  // and would bury it. Polygons are already triangles where the fan produced them, but anything
  // wider is fanned about its own vertex mean so the file is triangles throughout.
  void
  appendInterfaceSTL(std::ofstream& a_file, const RealVect& a_centre, const Real a_dx, long long& a_count) const
  {
    for (int ip = 0; ip < m_numPolygons; ip++) {
      const Polygon& p = m_polygon[ip];

      if (p.m_face >= 0 || p.m_numVertices < 3) {
        continue;
      }

      RealVect mean = RealVect::Zero;

      for (int i = 0; i < p.m_numVertices; i++) {
        mean += p.m_vertex[i];
      }

      mean /= static_cast<Real>(p.m_numVertices);

      const int numTriangles = (p.m_numVertices == 3) ? 1 : p.m_numVertices;

      for (int t = 0; t < numTriangles; t++) {
        RealVect v[3];

        if (p.m_numVertices == 3) {
          v[0] = p.m_vertex[0];
          v[1] = p.m_vertex[1];
          v[2] = p.m_vertex[2];
        }
        else {
          v[0] = mean;
          v[1] = p.m_vertex[t];
          v[2] = p.m_vertex[(t + 1) % p.m_numVertices];
        }

        for (int k = 0; k < 3; k++) {
          v[k] = a_centre + a_dx * v[k];
        }

        const RealVect e1 = v[1] - v[0];
        const RealVect e2 = v[2] - v[0];

        RealVect n = RealVect(
          D_DECL(e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]));

        const Real length = n.vectorLength();

        if (length > 0.0) {
          n /= length;
        }

        a_file << "  facet normal " << n[0] << " " << n[1] << " " << n[2] << "\n";
        a_file << "    outer loop\n";

        for (int k = 0; k < 3; k++) {
          a_file << "      vertex " << v[k][0] << " " << v[k][1] << " " << v[k][2] << "\n";
        }

        a_file << "    endloop\n  endfacet\n";

        a_count++;
      }
    }
  }

  // The body as the generator makes it today, one chord per face, for comparison.
  bool
  buildPlain(const PolyhedralEB::CutCellSurface& a_coarse)
  {
    return this->define(a_coarse) && m_kind == Kind::Cut;
  }

  // Counts first: how many loops the faces leave open, how planar they are, and what the polygon
  // total would be if every loop were fanned against if only the non-planar ones were. A planar
  // loop is already a polygon and needs no fan at all.
  bool
  closeInterface(int& a_loops, int& a_neededFan, int& a_neededPlanar, int& a_planarLoops, Real& a_worstBend)
  {
    a_loops        = 0;
    a_neededFan    = m_numPolygons;
    a_neededPlanar = m_numPolygons;
    a_planarLoops  = 0;
    a_worstBend    = 0.0;

    RealVect from[s_maxPolygons * s_maxVertices];
    RealVect to[s_maxPolygons * s_maxVertices];

    int numOpen = 0;

    for (int ip = 0; ip < m_numPolygons; ip++) {
      const Polygon& p = m_polygon[ip];

      for (int i = 0; i < p.m_numVertices; i++) {
        const RealVect& a = p.m_vertex[i];
        const RealVect& b = p.m_vertex[(i + 1) % p.m_numVertices];

        if (PolyhedralEB::detail::sameVertex(a, b)) {
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

            shared = PolyhedralEB::detail::sameVertex(a, d) && PolyhedralEB::detail::sameVertex(b, c);
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
        if (PolyhedralEB::detail::sameVertex(current, end)) {
          closed = true;

          break;
        }

        int next = -1;

        for (int j = 0; j < numOpen && next < 0; j++) {
          if (!used[j] && PolyhedralEB::detail::sameVertex(from[j], current)) {
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

      a_loops++;

      // how far the loop departs from a plane, measured against its own size so that a big loop
      // and a small one are judged the same way
      RealVect areaVector = RealVect::Zero;

      for (int i = 0; i < numLoop; i++) {
        const RealVect& u = loop[i];
        const RealVect& v = loop[(i + 1) % numLoop];

        areaVector += 0.5 *
                      RealVect(D_DECL(u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]));
      }

      const Real areaLength = areaVector.vectorLength();

      Real bend = 0.0;

      if (areaLength > 1.0E-30) {
        const RealVect normal = areaVector / areaLength;

        for (int i = 0; i < numLoop; i++) {
          bend = std::max(bend, std::abs((loop[i] - loop[0]).dotProduct(normal)));
        }

        bend /= std::sqrt(areaLength);
      }

      a_worstBend = std::max(a_worstBend, bend);

      const bool planar = (bend < 1.0E-10);

      a_planarLoops += planar ? 1 : 0;

      a_neededFan += numLoop;
      a_neededPlanar += planar ? 1 : numLoop;

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

  bool
  buildSeam(const BaseIF&                       a_implicitFunction,
            const PolyhedralEB::CutCellSurface& a_coarse,
            const RealVect&                     a_centre,
            const Real                          a_dx,
            const int                           a_seamDir,
            const int                           a_seamSide,
            int&                                a_polygons,
            int&                                a_widest,
            Real&                               a_gap,
            int&                                a_why,
            int&                                a_mergeWhy,
            int&                                a_leftover,
            int&                                a_subFaces,
            int&                                a_loops,
            int&                                a_neededFan,
            int&                                a_neededPlanar,
            int&                                a_planarLoops,
            Real&                               a_worstBend,
            Real&                               a_singleChord,
            Real&                               a_multiChord,
            Real&                               a_fineSum,
            int&                                a_components)
  {
    a_why        = 0;
    a_mergeWhy   = 0;
    a_leftover   = 0;
    a_subFaces   = 0;
    a_components = 0;

    if (!this->define(a_coarse)) {
      a_why = 1;

      return false;
    }

    if (m_kind != Kind::Cut) {
      a_why = 2;

      return false;
    }

    // what the single chord said this face's aperture was
    a_singleChord = this->areaFraction(a_seamDir, (a_seamSide == 0) ? Side::Lo : Side::Hi);
    a_fineSum     = 0.0;

    const Real before = this->closureResidual();

    // drop what the single chord put on this face
    const int seamFace = 2 * a_seamDir + a_seamSide;

    // the seam face's chord goes, and so does the interface, which was built to meet it
    int kept = 0;

    for (int ip = 0; ip < m_numPolygons; ip++) {
      if (m_polygon[ip].m_face != seamFace && m_polygon[ip].m_face >= 0) {
        m_polygon[kept++] = m_polygon[ip];
      }
    }

    m_numPolygons = kept;

    // and put the abutting cells' chords there instead
    Polygon sub[4 * (1 << (SpaceDim - 1))];

    int numSub = 0;

    for (int q = 0; q < (1 << (SpaceDim - 1)); q++) {
      // the quadrant of this cell touching the seam face
      int which = 0;
      int bit   = 0;

      for (int d = 0; d < SpaceDim; d++) {
        if (d == a_seamDir) {
          which |= a_seamSide << d;
        }
        else {
          which |= ((q >> bit) & 1) << d;
          bit++;
        }
      }

      RealVect origin;
      RealVect childCentre = a_centre;

      for (int d = 0; d < SpaceDim; d++) {
        origin[d] = -0.25 + 0.5 * static_cast<Real>((which >> d) & 1);
        childCentre[d] += 0.5 * a_dx * (((which >> d) & 1) - 0.5);
      }

      PolyhedralEB::CutCellSurface fine;
      makeSurface(a_implicitFunction, childCentre, 0.5 * a_dx, fine);

      Polygon walked[2];

      const int numWalked = this->faceWalk(a_seamDir, a_seamSide, fine, walked);

      if (numWalked < 0) {
        a_why = 3;

        return false;
      }

      // what the cell on the other side of this quadrant holds for the shared face
      PolyhedralEB::CutCellBody fineBody;

      if (fineBody.define(fine)) {
        a_fineSum += fineBody.areaFraction(a_seamDir, (a_seamSide == 0) ? Side::Lo : Side::Hi);
      }

      for (int n = 0; n < numWalked; n++) {
        if (numSub >= 4 * (1 << (SpaceDim - 1))) {
          a_why = 4;

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

    // one polygon for the face, not one per quadrant
    if (numSub > 0) {
      Polygon merged[1 << (SpaceDim - 1)];

      int numMerged = 0;
      int mergeWhy  = 0;

      if (!this->mergeCoplanar(sub, numSub, merged, 1 << (SpaceDim - 1), numMerged, mergeWhy)) {
        a_why      = 6;
        a_mergeWhy = mergeWhy;
        a_subFaces = numSub;

        return false;
      }

      // The four fine faces carry four apertures and the coarse face carries one. Where their union
      // is disconnected no single aperture describes it, and a coarse cell that is required to be
      // single-valued cannot hold the result. That is a refusal, not something to patch.
      a_components = numMerged;

      for (int n = 0; n < numMerged; n++) {
        if (m_numPolygons >= s_maxPolygons) {
          a_why = 4;

          return false;
        }

        m_polygon[m_numPolygons++] = merged[n];
      }
    }

    if (!this->closeInterface(a_loops, a_neededFan, a_neededPlanar, a_planarLoops, a_worstBend)) {
      a_why = 5;

      return false;
    }

    this->accumulateMoments();

    a_polygons   = m_numPolygons;
    a_widest     = this->widestPolygon();
    a_gap        = this->closureResidual();
    a_multiChord = this->areaFraction(a_seamDir, (a_seamSide == 0) ? Side::Lo : Side::Hi);

    (void)before;

    return true;
  }
};

// Measure what a multichord seam face costs and what it leaves for the interface to absorb.
void
validateSeamFace(const RefCountedPtr<ComputationalGeometry>& a_compgeom,
                 const RefCountedPtr<AmrMesh>&               a_amr,
                 const int                                   a_numCells)
{
  CH_TIME("validateSeamFace");

  const RefCountedPtr<BaseIF>& implicitFunction = a_compgeom->getGasImplicitFunction();

  if (implicitFunction.isNull()) {
    return;
  }

  const RealVect probLo = a_amr->getProbLo();
  const Real     dx     = a_amr->getDx()[0];

  long long cells            = 0;
  long long refused          = 0;
  long long why1             = 0;
  long long why2             = 0;
  long long why3             = 0;
  long long why4             = 0;
  long long why5             = 0;
  long long why6             = 0;
  long long mergeWhyCount[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  int       reportedRefusals = 0;
  int       maxNeedFan       = 0;
  int       maxNeedFlat      = 0;
  long long totalLoops       = 0;
  long long planarLoops      = 0;
  Real      worstBend        = 0.0;
  long long conserveBad      = 0;
  Real      worstConserve    = 0.0;
  Real      worstChordShift  = 0.0;
  int       maxPolys         = 0;
  int       maxVerts         = 0;
  Real      worstGap         = 0.0;

  Box slab = a_amr->getDomains()[0].domainBox();

  for (int d = 0; d < SpaceDim; d++) {
    if (slab.size(d) > a_numCells) {
      slab.setBig(d, slab.smallEnd(d) + a_numCells - 1);
    }
  }

  int stlFace = -1;
  {
    ParmParse pp("Prototype");
    pp.query("seamface_stl_face", stlFace);
  }

  std::ofstream multiFile;
  std::ofstream plainFile;

  long long multiFacets = 0;
  long long plainFacets = 0;

  if (stlFace >= 0 && procID() == 0) {
    multiFile.open("seam_multichord.stl");
    plainFile.open("seam_singlechord.stl");

    multiFile << std::scientific;
    plainFile << std::scientific;

    multiFile << "solid multichord\n";
    plainFile << "solid singlechord\n";
  }

  for (BoxIterator bit(slab); bit.ok(); ++bit) {
    RealVect centre = probLo;

    for (int d = 0; d < SpaceDim; d++) {
      centre[d] += dx * (bit()[d] + 0.5);
    }

    PolyhedralEB::CutCellSurface coarse;
    makeSurface(*implicitFunction, centre, dx, coarse);

    if (PolyhedralEB::CutCellBody::classify(coarse) != PolyhedralEB::CutCellBody::Kind::Cut) {
      continue;
    }

    if (stlFace >= 0 && plainFile.is_open()) {
      SeamBody plain;

      if (plain.buildPlain(coarse)) {
        plain.appendInterfaceSTL(plainFile, centre, dx, plainFacets);
      }
    }

    for (int dir = 0; dir < SpaceDim; dir++) {
      for (int side = 0; side < 2; side++) {
        SeamBody body;

        int  polys     = 0;
        int  verts     = 0;
        Real gap       = 0.0;
        int  why       = 0;
        int  loops     = 0;
        int  needFan   = 0;
        int  needFlat  = 0;
        int  flatLoops = 0;
        Real bend      = 0.0;
        Real single    = 0.0;
        Real multi     = 0.0;
        Real fine      = 0.0;

        int mergeWhy   = 0;
        int leftover   = 0;
        int subFaces   = 0;
        int components = 0;

        const bool ok = body.buildSeam(*implicitFunction,
                                       coarse,
                                       centre,
                                       dx,
                                       dir,
                                       side,
                                       polys,
                                       verts,
                                       gap,
                                       why,
                                       mergeWhy,
                                       leftover,
                                       subFaces,
                                       loops,
                                       needFan,
                                       needFlat,
                                       flatLoops,
                                       bend,
                                       single,
                                       multi,
                                       fine,
                                       components);

        maxNeedFan  = std::max(maxNeedFan, needFan);
        maxNeedFlat = std::max(maxNeedFlat, needFlat);
        totalLoops += loops;
        planarLoops += flatLoops;
        worstBend = std::max(worstBend, bend);

        if (!ok) {
          refused++;
          why1 += (why == 1);
          why2 += (why == 2);
          why3 += (why == 3);
          why4 += (why == 4);
          why5 += (why == 5);
          why6 += (why == 6);

          if (why == 6) {
            mergeWhyCount[mergeWhy]++;

            if (mergeWhy == 4 && reportedRefusals < 3) {
              reportedRefusals++;

              pout() << "  DUMP cell " << bit() << " face " << (2 * dir + side) << " subFacePolygons " << subFaces
                     << " leftoverEdges " << leftover << endl;

              // what the coarse cell looks like
              PolyhedralEB::CutCellBody coarseBody;
              coarseBody.define(coarse);

              pout() << "    coarse corners";
              for (int c = 0; c < PolyhedralEB::CutCellSurface::s_numCorners; c++) {
                pout() << " " << (PolyhedralEB::isFluid(coarse.m_corner[c]) ? "F" : "S");
              }
              pout() << "  kappa " << coarseBody.volumeFraction() << " oneSided " << coarseBody.interfaceIsOneSided()
                     << " bndryArea " << coarseBody.boundaryArea() << " trueArea " << coarseBody.trueBoundaryArea()
                     << endl;

              // and each quadrant of the face
              for (int q = 0; q < (1 << (SpaceDim - 1)); q++) {
                int which = 0;
                int bitq  = 0;
                for (int d = 0; d < SpaceDim; d++) {
                  if (d == dir) {
                    which |= side << d;
                  }
                  else {
                    which |= ((q >> bitq) & 1) << d;
                    bitq++;
                  }
                }

                RealVect childCentre = centre;
                for (int d = 0; d < SpaceDim; d++) {
                  childCentre[d] += 0.5 * dx * (((which >> d) & 1) - 0.5);
                }

                PolyhedralEB::CutCellSurface fineSurf;
                makeSurface(*implicitFunction, childCentre, 0.5 * dx, fineSurf);

                PolyhedralEB::CutCellBody fineBody;
                const bool                fineOk = fineBody.define(fineSurf);

                const char* kindName = "?";
                if (fineOk) {
                  switch (fineBody.kind()) {
                  case PolyhedralEB::CutCellBody::Kind::Regular:
                    kindName = "Regular";
                    break;
                  case PolyhedralEB::CutCellBody::Kind::Covered:
                    kindName = "Covered";
                    break;
                  default:
                    kindName = "Cut";
                    break;
                  }
                }

                pout() << "    quadrant " << q << " kind " << kindName << " kappa "
                       << (fineOk ? fineBody.volumeFraction() : -1.0) << " faceAperture "
                       << (fineOk ? fineBody.areaFraction(dir, (side == 0) ? Side::Lo : Side::Hi) : -1.0) << endl;
              }
            }
          }

          continue;
        }

        cells++;

        // Conservation of the seam face: the multichord aperture has to be what the cells on the
        // other side hold, averaged over the quadrants, or the two sides of the seam disagree
        // about how much is open.
        const Real conserve = std::abs(multi - 0.25 * fine);

        worstConserve = std::max(worstConserve, conserve);

        if (conserve > 1.0E-12) {
          conserveBad++;
        }

        worstChordShift = std::max(worstChordShift, std::abs(multi - single));

        if (stlFace >= 0 && multiFile.is_open() && (2 * dir + side) == stlFace) {
          body.appendInterfaceSTL(multiFile, centre, dx, multiFacets);
        }

        maxPolys = std::max(maxPolys, polys);
        maxVerts = std::max(maxVerts, verts);
        worstGap = std::max(worstGap, gap);
      }
    }
  }

  pout() << "SEAMFACE refusals: defineFailed " << why1 << " notCut " << why2 << " faceWalkFailed " << why3
         << " overPolygonCap " << why4 << " interfaceFailed " << why5 << " mergeFailed " << why6 << endl;
  pout() << "SEAMFACE mergeReasons: tooFewEdges " << mergeWhyCount[1] << " walkStuck " << mergeWhyCount[2]
         << " didNotClose " << mergeWhyCount[3] << " moreThanOneLoop " << mergeWhyCount[4] << " vertexOverflow "
         << mergeWhyCount[5] << " degenerate " << mergeWhyCount[6] << endl;
  // The cells the surface lives in, as a rectilinear mesh VisIt can overlay on the STL. Cell data
  // says what the generator made of each one, so a facet can be read against the cell that
  // produced it.
  if (stlFace >= 0 && procID() == 0) {
    std::ofstream grid("seam_grid.vtk");

    grid << std::scientific;

    const IntVect lo = slab.smallEnd();
    const IntVect hi = slab.bigEnd();

    grid << "# vtk DataFile Version 3.0\n";
    grid << "Cartesian cells carrying the cut-cell surface\n";
    grid << "ASCII\n";
    grid << "DATASET RECTILINEAR_GRID\n";
    grid << "DIMENSIONS " << (hi[0] - lo[0] + 2) << " " << (hi[1] - lo[1] + 2) << " "
         << (SpaceDim == 3 ? (hi[SpaceDim - 1] - lo[SpaceDim - 1] + 2) : 1) << "\n";

    const char* axis[3] = {"X_COORDINATES", "Y_COORDINATES", "Z_COORDINATES"};

    for (int d = 0; d < 3; d++) {
      if (d < SpaceDim) {
        grid << axis[d] << " " << (hi[d] - lo[d] + 2) << " double\n";

        for (int i = lo[d]; i <= hi[d] + 1; i++) {
          grid << (probLo[d] + dx * i) << " ";
        }

        grid << "\n";
      }
      else {
        grid << axis[d] << " 1 double\n0\n";
      }
    }

    const long long numCells = slab.numPts();

    grid << "CELL_DATA " << numCells << "\n";
    grid << "SCALARS kind int 1\nLOOKUP_TABLE default\n";

    for (BoxIterator bit(slab); bit.ok(); ++bit) {
      RealVect centre = probLo;

      for (int d = 0; d < SpaceDim; d++) {
        centre[d] += dx * (bit()[d] + 0.5);
      }

      PolyhedralEB::CutCellSurface surface;
      makeSurface(*implicitFunction, centre, dx, surface);

      const PolyhedralEB::CutCellBody::Kind kind = PolyhedralEB::CutCellBody::classify(surface);

      grid << ((kind == PolyhedralEB::CutCellBody::Kind::Regular)
                 ? 1
                 : ((kind == PolyhedralEB::CutCellBody::Kind::Covered) ? -1 : 0))
           << "\n";
    }

    grid << "SCALARS kappa double 1\nLOOKUP_TABLE default\n";

    for (BoxIterator bit(slab); bit.ok(); ++bit) {
      RealVect centre = probLo;

      for (int d = 0; d < SpaceDim; d++) {
        centre[d] += dx * (bit()[d] + 0.5);
      }

      PolyhedralEB::CutCellSurface surface;
      makeSurface(*implicitFunction, centre, dx, surface);

      PolyhedralEB::CutCellBody body;

      grid << (body.define(surface) ? body.volumeFraction() : -1.0) << "\n";
    }

    grid.close();

    pout() << "SEAMFACE wrote seam_grid.vtk over " << numCells << " cells" << endl;
  }

  if (multiFile.is_open()) {
    multiFile << "endsolid multichord\n";
    plainFile << "endsolid singlechord\n";

    multiFile.close();
    plainFile.close();

    pout() << "SEAMFACE stl: multichord facets " << multiFacets << " singlechord facets " << plainFacets << endl;
  }

  pout() << "SEAMFACE loops " << totalLoops << " planar " << planarLoops << " worstBend " << worstBend
         << " polygonsNeeded fan " << maxNeedFan << " planarKept " << maxNeedFlat << endl;
  pout() << "SEAMFACE conservation: bad " << conserveBad << " worst " << worstConserve << " chordShift "
         << worstChordShift << endl;

  pout() << "SEAMFACE built " << cells << " refused " << refused << " maxPolygons " << maxPolys << " (cap 20)"
         << " maxVertices " << maxVerts << " (cap 20) worstClosureGap " << worstGap << endl;
}

// Pass A of the seam restriction: which of a coarse cell's edges does endpoint bracketing get wrong?
//
// An edge with fluid at both ends and solid across the middle carries two crossings, and the sign
// test between the ends reports none. Split in half, each piece holds one crossing and the level
// below finds both. Where the two counts differ, the coarse cell's 1-D description is not the
// restriction of the fine one, and every face meeting that edge inherits the error: the faces
// disagree about whether the edge borders fluid, and closeInterface renders the disagreement as
// interface lying along the edge.
//
// This only detects. It marks the edges a patch would have to rebuild.
int
markCoarseEdges(const BaseIF&   a_implicitFunction,
                const RealVect& a_centre,
                const Real      a_dx,
                const int       a_seamDir,
                const int       a_seamSide,
                bool*           a_marked,
                int&            a_onSeamFace,
                int&            a_elsewhere)
{
  a_onSeamFace = 0;
  a_elsewhere  = 0;

  int total = 0;

  for (int e = 0; e < PolyhedralEB::CutCellSurface::s_numEdges; e++) {
    int     dir;
    IntVect loOffset;
    edgeGeometry(e, dir, loOffset);

    RealVect lo = a_centre;

    for (int d = 0; d < SpaceDim; d++) {
      lo[d] += a_dx * (loOffset[d] - 0.5);
    }

    RealVect hi = lo;
    hi[dir] += a_dx;

    RealVect mid = lo;
    mid[dir] += 0.5 * a_dx;

    const bool fLo  = PolyhedralEB::isFluid(a_implicitFunction.value(lo));
    const bool fHi  = PolyhedralEB::isFluid(a_implicitFunction.value(hi));
    const bool fMid = PolyhedralEB::isFluid(a_implicitFunction.value(mid));

    const int coarse   = (fLo != fHi) ? 1 : 0;
    const int restrict = ((fLo != fMid) ? 1 : 0) + ((fMid != fHi) ? 1 : 0);

    a_marked[e] = (coarse != restrict);

    if (a_marked[e]) {
      total++;

      // the seam face is the one whose fine neighbours supply the restriction; an edge lies in it
      // when it runs transverse to the seam and sits on the seam side
      if (dir != a_seamDir && loOffset[a_seamDir] == a_seamSide) {
        a_onSeamFace++;
      }
      else {
        a_elsewhere++;
      }
    }
  }

  return total;
}

// A real coarse-fine boundary, and what the seam looks like across it.
//
// Half the domain is carried at one spacing and half at twice the resolution, which is the
// arrangement the multichord exists for and the one the earlier harness never had: there, every
// cell subdivided a face of its own, so a multichord face always met a single-chord neighbour and
// the surface cracked by construction.
//
// Both blocks together cover the whole body, so the union of the interface triangles is a closed
// surface exactly when the seam agrees, and edge valence in the exported STL decides it. Two
// surfaces are written where a path is given: they differ only in what the last column of coarse
// cells puts on the face it shares with the fine block -- one chord of its own, or the chords of
// the cells on the other side.
void
twoLevelSeam(const BaseIF&      a_implicitFunction,
             const RealVect&    a_probLo,
             const Real         a_dxCoarse,
             const int          a_numCoarse,
             const int          a_split,
             const std::string& a_stlPrefix,
             const IntVect&     a_dumpCell,
             const bool         a_passB,
             int&               a_coarseCut,
             int&               a_coarsePlainRefused,
             int&               a_fineCut,
             int&               a_finePlainRefused,
             int&               a_seamCells,
             int&               a_seamRefused,
             int*               a_seamWhy,
             int*               a_mergeWhy,
             const int          a_numWhy,
             Real&              a_worstMulti,
             Real&              a_worstSingle,
             Real&              a_worstClosure,
             int&               a_markedEdges,
             int&               a_markedOnSeamFace,
             int&               a_cellsMarked,
             int&               a_cellsMarkedOnSeamFace,
             int&               a_cellsTorn,
             int&               a_predicted,
             int&               a_missed,
             int&               a_facesRestricted,
             int&               a_restrictRefused,
             int&               a_blindNeighbours,
             int&               a_multiValued,
             int&               a_multiValuedTorn)
{
  CH_TIME("twoLevelSeam");

  const Real dxC = a_dxCoarse;
  const Real dxF = 0.5 * dxC;

  a_coarseCut             = 0;
  a_coarsePlainRefused    = 0;
  a_fineCut               = 0;
  a_finePlainRefused      = 0;
  a_seamCells             = 0;
  a_seamRefused           = 0;
  a_worstMulti            = 0.0;
  a_worstSingle           = 0.0;
  a_worstClosure          = 0.0;
  a_markedEdges           = 0;
  a_markedOnSeamFace      = 0;
  a_cellsMarked           = 0;
  a_cellsMarkedOnSeamFace = 0;
  a_cellsTorn             = 0;
  a_predicted             = 0;
  a_missed                = 0;

  for (int i = 0; i < a_numWhy; i++) {
    a_seamWhy[i]  = 0;
    a_mergeWhy[i] = 0;
  }

  const bool writing = !a_stlPrefix.empty();

  std::ofstream multi;
  std::ofstream single;

  long long multiFacets  = 0;
  long long singleFacets = 0;

  if (writing) {
    multi.open(a_stlPrefix + "_multichord.stl");
    single.open(a_stlPrefix + "_singlechord.stl");

    multi << std::scientific << std::setprecision(17) << "solid multichord\n";
    single << std::scientific << std::setprecision(17) << "solid singlechord\n";
  }

  const Box coarseBox(IntVect::Zero, (a_numCoarse - 1) * IntVect::Unit);

  for (BoxIterator bit(coarseBox); bit.ok(); ++bit) {
    if (bit()[0] >= a_split) {
      continue;
    }

    RealVect centre = a_probLo;

    for (int d = 0; d < SpaceDim; d++) {
      centre[d] += dxC * (bit()[d] + 0.5);
    }

    PolyhedralEB::CutCellSurface coarse;
    makeSurface(a_implicitFunction, centre, dxC, coarse);

    const bool isSeamColumn = (bit()[0] == a_split - 1);

    if (PolyhedralEB::CutCellBody::classify(coarse) != PolyhedralEB::CutCellBody::Kind::Cut) {
      // A cell the corner test calls uniform can still have a face the level below cuts: its edges
      // carry an even number of crossings, so every corner reports the same side. It is skipped
      // here, which is why restricting its cut neighbour's shared face only moves the hole.
      if (isSeamColumn) {
        bool marked[PolyhedralEB::CutCellSurface::s_numEdges];

        int onSeamFace = 0;
        int elsewhere  = 0;

        markCoarseEdges(a_implicitFunction, centre, dxC, 0, 1, marked, onSeamFace, elsewhere);

        if (onSeamFace > 0) {
          a_blindNeighbours++;
        }
      }

      continue;
    }

    a_coarseCut++;

    SeamBody plain;

    if (!plain.buildPlain(coarse)) {
      a_coarsePlainRefused++;

      continue;
    }

    if (writing) {
      plain.appendInterfaceSTL(single, centre, dxC, singleFacets);
    }

    bool wroteMulti = false;

    if (isSeamColumn) {
      SeamBody seam;

      int  polys = 0, verts = 0, why = 0, mergeWhy = 0, leftover = 0, subFaces = 0;
      int  loops = 0, needFan = 0, needFlat = 0, flatLoops = 0;
      Real gap = 0.0, bend = 0.0, singleChord = 0.0, multiChord = 0.0, fineSum = 0.0;

      int components = 0;

      const bool dumping = writing && (bit() == a_dumpCell);

      if (dumping) {
        SeamBody before;

        if (before.buildPlain(coarse)) {
          before.dumpPolygons("plain");
        }
      }

      bool marked[PolyhedralEB::CutCellSurface::s_numEdges];

      int onSeamFace = 0;
      int elsewhere  = 0;

      const int numMarked = markCoarseEdges(a_implicitFunction, centre, dxC, 0, 1, marked, onSeamFace, elsewhere);

      a_markedEdges += numMarked;
      a_markedOnSeamFace += onSeamFace;

      if (numMarked > 0) {
        a_cellsMarked++;
      }

      if (onSeamFace > 0) {
        a_cellsMarkedOnSeamFace++;
      }

      const bool ok = seam.buildSeam(a_implicitFunction,
                                     coarse,
                                     centre,
                                     dxC,
                                     0,
                                     1,
                                     polys,
                                     verts,
                                     gap,
                                     why,
                                     mergeWhy,
                                     leftover,
                                     subFaces,
                                     loops,
                                     needFan,
                                     needFlat,
                                     flatLoops,
                                     bend,
                                     singleChord,
                                     multiChord,
                                     fineSum,
                                     components);

      a_seamCells++;

      if (dumping) {
        seam.dumpPolygons("seam");
      }

      if (ok && a_passB) {
        // Pass B: rebuild the seam face and every face meeting it across a marked edge, so the
        // faces sharing that edge cannot contradict one another
        SeamBody restricted;

        int faces = 0, whyB = 0, mergeWhyB = 0;

        if (restricted.buildRestricted(a_implicitFunction, coarse, centre, dxC, 0, 1, marked, faces, whyB, mergeWhyB)) {
          seam = restricted;

          a_facesRestricted += faces;
        }
        else {
          a_restrictRefused++;

          if (writing) {
            pout() << "PASSB refused " << bit() << " why " << whyB << " mergeWhy " << mergeWhyB << " faces " << faces
                   << endl;
          }
        }
      }

      if (ok) {
        const int torn = seam.cellEdgeInterfaceSegments();

        if (components > 1) {
          a_multiValued++;

          if (torn > 0) {
            a_multiValuedTorn++;
          }
        }

        if (torn > 0) {
          a_cellsTorn++;

          if (onSeamFace > 0) {
            a_predicted++;
          }
          else {
            a_missed++;

            if (writing) {
              pout() << "PASSA missed " << bit() << " tornSegments " << torn << " markedEdges " << numMarked
                     << " onSeamFace " << onSeamFace << " elsewhere " << elsewhere << endl;
            }
          }
        }

        if (writing) {
          seam.appendInterfaceSTL(multi, centre, dxC, multiFacets);
        }

        wroteMulti = true;

        a_worstMulti   = std::max(a_worstMulti, std::abs(multiChord - 0.25 * fineSum));
        a_worstSingle  = std::max(a_worstSingle, std::abs(singleChord - 0.25 * fineSum));
        a_worstClosure = std::max(a_worstClosure, gap);
      }
      else {
        a_seamRefused++;

        if (writing) {
          pout() << "TWOLEVEL refused " << bit() << " why " << why << " mergeWhy " << mergeWhy << " subFaces "
                 << subFaces << " leftover " << leftover << " singleChord " << singleChord << " fineSum " << fineSum
                 << endl;
        }

        if (why >= 0 && why < a_numWhy) {
          a_seamWhy[why]++;
        }

        if (mergeWhy > 0 && mergeWhy < a_numWhy) {
          a_mergeWhy[mergeWhy]++;
        }
      }
    }

    if (writing && !wroteMulti) {
      plain.appendInterfaceSTL(multi, centre, dxC, multiFacets);
    }
  }

  // the fine block, identical in both files
  const Box fineBox(IntVect(D_DECL(2 * a_split, 0, 0)),
                    IntVect(D_DECL(2 * a_numCoarse - 1, 2 * a_numCoarse - 1, 2 * a_numCoarse - 1)));

  for (BoxIterator bit(fineBox); bit.ok(); ++bit) {
    RealVect centre = a_probLo;

    for (int d = 0; d < SpaceDim; d++) {
      centre[d] += dxF * (bit()[d] + 0.5);
    }

    PolyhedralEB::CutCellSurface fine;
    makeSurface(a_implicitFunction, centre, dxF, fine);

    if (PolyhedralEB::CutCellBody::classify(fine) != PolyhedralEB::CutCellBody::Kind::Cut) {
      continue;
    }

    a_fineCut++;

    SeamBody body;

    if (!body.buildPlain(fine)) {
      a_finePlainRefused++;

      continue;
    }

    if (writing) {
      body.appendInterfaceSTL(multi, centre, dxF, multiFacets);
      body.appendInterfaceSTL(single, centre, dxF, singleFacets);
    }
  }

  if (writing) {
    multi << "endsolid multichord\n";
    single << "endsolid singlechord\n";

    multi.close();
    single.close();
  }
}

// The two blocks as meshes VisIt can overlay.
void
writeTwoLevelGrids(const std::string& a_prefix,
                   const RealVect&    a_probLo,
                   const Real         a_dxCoarse,
                   const int          a_numCoarse,
                   const int          a_split)
{
  for (int which = 0; which < 2; which++) {
    const bool isFine = (which == 1);
    const Real dx     = isFine ? 0.5 * a_dxCoarse : a_dxCoarse;
    const int  iLo    = isFine ? 2 * a_split : 0;
    const int  iHi    = isFine ? 2 * a_numCoarse : a_split;
    const int  jHi    = isFine ? 2 * a_numCoarse : a_numCoarse;

    std::ofstream grid(a_prefix + (isFine ? "_fine.vtk" : "_coarse.vtk"));

    grid << std::scientific;
    grid << "# vtk DataFile Version 3.0\n" << (isFine ? "fine block\n" : "coarse block\n") << "ASCII\n";
    grid << "DATASET RECTILINEAR_GRID\n";
    grid << "DIMENSIONS " << (iHi - iLo + 1) << " " << (jHi + 1) << " " << (jHi + 1) << "\n";

    grid << "X_COORDINATES " << (iHi - iLo + 1) << " double\n";
    for (int i = iLo; i <= iHi; i++) {
      grid << (a_probLo[0] + dx * i) << " ";
    }
    grid << "\n";

    const char* rest[2] = {"Y_COORDINATES", "Z_COORDINATES"};

    for (int d = 0; d < 2; d++) {
      grid << rest[d] << " " << (jHi + 1) << " double\n";
      for (int j = 0; j <= jHi; j++) {
        grid << (a_probLo[d + 1] + dx * j) << " ";
      }
      grid << "\n";
    }

    grid.close();
  }
}

// The seam on the geometry supplied by the inputs file, exported for inspection.
void
validateTwoLevelSeam(const RefCountedPtr<ComputationalGeometry>& a_compgeom,
                     const RefCountedPtr<AmrMesh>&               a_amr,
                     const int                                   a_numCoarse,
                     const int                                   a_split)
{
  CH_TIME("validateTwoLevelSeam");

  const RefCountedPtr<BaseIF>& implicitFunction = a_compgeom->getGasImplicitFunction();

  if (implicitFunction.isNull() || procID() != 0) {
    return;
  }

  const RealVect probLo = a_amr->getProbLo();
  const Real     dxC    = a_amr->getDx()[0] * (a_amr->getDomains()[0].domainBox().size(0) / a_numCoarse);

  constexpr int numWhy = 8;

  IntVect dumpCell = IntVect::Unit * (-1);
  {
    Vector<int> v;
    ParmParse   pp("Prototype");

    if (pp.contains("twolevel_dump_cell")) {
      pp.getarr("twolevel_dump_cell", v, 0, SpaceDim);

      for (int d = 0; d < SpaceDim; d++) {
        dumpCell[d] = v[d];
      }
    }
  }

  int why[numWhy];
  int mergeWhy[numWhy];
  int coarseCut = 0, coarsePlain = 0, fineCut = 0, finePlain = 0, seamCells = 0, seamRefused = 0;
  int markedEdges = 0, markedOnSeamFace = 0, cellsMarked = 0, cellsMarkedOnSeamFace = 0;
  int cellsTorn = 0, predicted = 0, missed = 0, facesRestricted = 0, restrictRefused = 0;
  int blindNeighbours = 0, multiValued = 0, multiValuedTorn = 0;

  bool passB = false;
  {
    ParmParse pp("Prototype");
    int       v = 0;
    pp.query("twolevel_passb", v);
    passB = (v > 0);
  }

  Real worstMulti = 0.0, worstSingle = 0.0, worstClosure = 0.0;

  twoLevelSeam(*implicitFunction,
               probLo,
               dxC,
               a_numCoarse,
               a_split,
               "seam2",
               dumpCell,
               passB,
               coarseCut,
               coarsePlain,
               fineCut,
               finePlain,
               seamCells,
               seamRefused,
               why,
               mergeWhy,
               numWhy,
               worstMulti,
               worstSingle,
               worstClosure,
               markedEdges,
               markedOnSeamFace,
               cellsMarked,
               cellsMarkedOnSeamFace,
               cellsTorn,
               predicted,
               missed,
               facesRestricted,
               restrictRefused,
               blindNeighbours,
               multiValued,
               multiValuedTorn);

  writeTwoLevelGrids("seam2", probLo, dxC, a_numCoarse, a_split);

  pout() << "TWOLEVEL seamCells " << seamCells << " refused " << seamRefused << " worstMultichord " << worstMulti
         << " worstSingleChord " << worstSingle << " worstClosure " << worstClosure << endl;
  pout() << "TWOLEVEL coarseCut " << coarseCut << " coarseRefused " << coarsePlain << " fineCut " << fineCut
         << " fineRefused " << finePlain << endl;
  pout() << "PASSA markedEdges " << markedEdges << " onSeamFace " << markedOnSeamFace << " cellsMarked " << cellsMarked
         << " cellsMarkedOnSeamFace " << cellsMarkedOnSeamFace << " cellsTorn " << cellsTorn << " predicted "
         << predicted << " missed " << missed << endl;
  pout() << "PASSB facesRestricted " << facesRestricted << " refused " << restrictRefused << " blindNeighbours "
         << blindNeighbours << endl;
}

// The same seam, swept over the Euler angles of a rotated cube. A cube is the hard case: its edges
// and corners land on cell centres and cell faces, which is where the chords degenerate, and
// rotating it walks those degeneracies through every relative orientation.
void
sweepTwoLevelSeam(const RefCountedPtr<AmrMesh>& a_amr,
                  const int                     a_numCoarse,
                  const int                     a_split,
                  const int                     a_samples,
                  const Real                    a_span,
                  const Real                    a_size,
                  const RealVect&               a_center,
                  const bool                    a_write,
                  const bool                    a_passB)
{
  CH_TIME("sweepTwoLevelSeam");

  if (procID() != 0) {
    return;
  }

  const RealVect probLo = a_amr->getProbLo();
  const Real     dxC    = a_amr->getDx()[0] * (a_amr->getDomains()[0].domainBox().size(0) / a_numCoarse);

  constexpr int numWhy = 8;

  int  totalRotations = 0, totalSeam = 0, totalRefused = 0, badRotations = 0;
  int  sumMarkedEdges = 0, sumOnSeamFace = 0, sumCellsMarked = 0, sumCellsOnSeamFace = 0;
  int  sumCellsTorn = 0, sumPredicted = 0, sumMissed = 0, sumFaces = 0, sumRestrictRefused = 0, sumBlind = 0;
  int  sumMultiValued = 0, sumMultiValuedTorn = 0;
  int  totalWhy[numWhy];
  int  totalMergeWhy[numWhy];
  Real worstMulti = 0.0, worstSingle = 0.0, worstClosure = 0.0;

  RealVect worstAngles = RealVect::Zero;

  for (int i = 0; i < numWhy; i++) {
    totalWhy[i]      = 0;
    totalMergeWhy[i] = 0;
  }

  for (int ia = 0; ia < a_samples; ia++) {
    for (int ib = 0; ib < a_samples; ib++) {
      for (int ic = 0; ic < a_samples; ic++) {
        const RealVect angles(D_DECL(a_span * ia / a_samples, a_span * ib / a_samples, a_span * ic / a_samples));

        Vector<RealVect> normals;
        Vector<Real>     offsets;

        for (int d = 0; d < SpaceDim; d++) {
          for (int s = 0; s < 2; s++) {
            RealVect n = RealVect::Zero;
            n[d]       = (s == 0) ? -1.0 : 1.0;

            normals.push_back(rotate(n, angles));
            offsets.push_back(a_size);
          }
        }

        const ConvexBody cube(normals, offsets, a_center);

        int  why[numWhy];
        int  mergeWhy[numWhy];
        int  coarseCut = 0, coarsePlain = 0, fineCut = 0, finePlain = 0, seamCells = 0, seamRefused = 0;
        int  markedEdges = 0, markedOnSeamFace = 0, cellsMarked = 0, cellsMarkedOnSeamFace = 0;
        int  cellsTorn = 0, predicted = 0, missed = 0, facesRestricted = 0, restrictRefused = 0;
        int  blindNeighbours = 0, multiValued = 0, multiValuedTorn = 0;
        Real multi = 0.0, singleChord = 0.0, closure = 0.0;

        std::string prefix;

        if (a_write) {
          char name[64];
          snprintf(name, sizeof(name), "sweep_%02d_%02d_%02d", ia, ib, ic);
          prefix = name;
        }

        twoLevelSeam(cube,
                     probLo,
                     dxC,
                     a_numCoarse,
                     a_split,
                     prefix,
                     IntVect::Unit * (-1),
                     a_passB,
                     coarseCut,
                     coarsePlain,
                     fineCut,
                     finePlain,
                     seamCells,
                     seamRefused,
                     why,
                     mergeWhy,
                     numWhy,
                     multi,
                     singleChord,
                     closure,
                     markedEdges,
                     markedOnSeamFace,
                     cellsMarked,
                     cellsMarkedOnSeamFace,
                     cellsTorn,
                     predicted,
                     missed,
                     facesRestricted,
                     restrictRefused,
                     blindNeighbours,
                     multiValued,
                     multiValuedTorn);

        totalRotations++;
        totalSeam += seamCells;
        totalRefused += seamRefused;
        sumMarkedEdges += markedEdges;
        sumOnSeamFace += markedOnSeamFace;
        sumCellsMarked += cellsMarked;
        sumCellsOnSeamFace += cellsMarkedOnSeamFace;
        sumCellsTorn += cellsTorn;
        sumPredicted += predicted;
        sumMissed += missed;
        sumFaces += facesRestricted;
        sumRestrictRefused += restrictRefused;
        sumBlind += blindNeighbours;
        sumMultiValued += multiValued;
        sumMultiValuedTorn += multiValuedTorn;

        for (int i = 0; i < numWhy; i++) {
          totalWhy[i] += why[i];
          totalMergeWhy[i] += mergeWhy[i];
        }

        const bool bad = (seamRefused > 0) || (coarsePlain > 0) || (finePlain > 0) || (multi > 1.0E-10);

        if (bad) {
          badRotations++;

          pout() << "SWEEP bad angles " << angles << " seamCells " << seamCells << " refused " << seamRefused
                 << " coarseRefused " << coarsePlain << " fineRefused " << finePlain << " worstMultichord " << multi
                 << " worstClosure " << closure << endl;
        }

        if (multi > worstMulti) {
          worstMulti  = multi;
          worstAngles = angles;
        }

        worstSingle  = std::max(worstSingle, singleChord);
        worstClosure = std::max(worstClosure, closure);
      }
    }
  }

  pout() << "SWEEP rotations " << totalRotations << " bad " << badRotations << " seamCells " << totalSeam << " refused "
         << totalRefused << endl;
  pout() << "SWEEP why";
  for (int i = 0; i < numWhy; i++) {
    pout() << " " << i << ":" << totalWhy[i];
  }
  pout() << " mergeWhy";
  for (int i = 0; i < numWhy; i++) {
    pout() << " " << i << ":" << totalMergeWhy[i];
  }
  pout() << endl;
  pout() << "SWEEP worstMultichord " << worstMulti << " at " << worstAngles << " worstSingleChord " << worstSingle
         << " worstClosure " << worstClosure << endl;
  pout() << "PASSA markedEdges " << sumMarkedEdges << " onSeamFace " << sumOnSeamFace << " cellsMarked "
         << sumCellsMarked << " cellsMarkedOnSeamFace " << sumCellsOnSeamFace << " cellsTorn " << sumCellsTorn
         << " predicted " << sumPredicted << " missed " << sumMissed << endl;
  pout() << "PASSB facesRestricted " << sumFaces << " refused " << sumRestrictRefused << " blindNeighbours " << sumBlind
         << endl;
  pout() << "MULTIVALUED seamFaces " << sumMultiValued << " ofWhichTorn " << sumMultiValuedTorn << endl;
}

// Check that a partially carried index space is sound before anything is asked to solve on it.
//
// Four things, reported per level rather than asserted, so that one run says everything that is
// wrong rather than the first thing. Faces first, since a face the two cells sharing it disagree
// about is what every failure so far has come down to.
void
validateIndexSpace(const RefCountedPtr<AmrMesh>& a_amr, const RefCountedPtr<ComputationalGeometry>& a_compgeom)
{
  const RefCountedPtr<EBIndexSpace>& ebis   = a_compgeom->getMfIndexSpace()->getEBIndexSpace(phase::gas);
  const Vector<DisjointBoxLayout>&   grids  = a_amr->getGrids(Realm::primal);
  const Vector<EBISLayout>&          ebisl  = a_amr->getEBISLayout(Realm::primal, phase::gas);
  const Vector<Real>&                dx     = a_amr->getDx();
  const Vector<int>&                 refRat = a_amr->getRefinementRatios();

  for (int lvl = 0; lvl <= a_amr->getFinestLevel(); lvl++) {
    long int faceMismatch    = 0;
    long int openIntoCovered = 0;
    long int divergenceBad   = 0;
    long int badFullCells    = 0;
    long int acrossCovered   = 0;
    long int acrossRegular   = 0;
    long int acrossIrregular = 0;
    long int acrossOutside   = 0;
    Real     worstDivergence = 0.0;
    int      reported        = 0;

    // Where the next level down sits, so a fault can be placed relative to the coarse-fine
    // boundary rather than only counted.
    IntVectSet finerCoverage;
    if (lvl < a_amr->getFinestLevel()) {
      const Vector<Box> finerBoxes = grids[lvl + 1].boxArray();
      for (int i = 0; i < finerBoxes.size(); i++) {
        finerCoverage |= coarsen(finerBoxes[i], refRat[lvl]);
      }
    }

    IntVectSet levelCoverage;
    {
      const Vector<Box> levelBoxes = grids[lvl].boxArray();
      for (int i = 0; i < levelBoxes.size(); i++) {
        levelCoverage |= levelBoxes[i];
      }
    }

    for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
      const EBISBox&   ebisBox = ebisl[lvl][dit()];
      const Box&       box     = grids[lvl][dit()];
      const IntVectSet irreg   = ebisBox.getIrregIVS(box);

      for (IVSIterator ivsIt(irreg); ivsIt.ok(); ++ivsIt) {
        const Vector<VolIndex> vofs = ebisBox.getVoFs(ivsIt());

        for (int iv = 0; iv < vofs.size(); iv++) {
          const VolIndex& vof = vofs[iv];

          RealVect apertureVector = RealVect::Zero;

          for (int dir = 0; dir < SpaceDim; dir++) {
            for (SideIterator sit; sit.ok(); ++sit) {
              const Vector<FaceIndex> faces = ebisBox.getFaces(vof, dir, sit());

              Real area = 0.0;
              for (int f = 0; f < faces.size(); f++) {
                area += ebisBox.areaFrac(faces[f]);
              }

              apertureVector[dir] += (sit() == Side::Hi) ? area : -area;

              // The cell across each face has to agree that the face is there.
              const IntVect other = vof.gridIndex() + sign(sit()) * BASISV(dir);

              if (!box.contains(other)) {
                continue;
              }

              if (ebisBox.isCovered(other) && area > 0.0) {
                openIntoCovered++;
              }

              if (!ebisBox.isCovered(other)) {
                const Vector<VolIndex> otherVoFs = ebisBox.getVoFs(other);
                int                    backFaces = 0;
                for (int j = 0; j < otherVoFs.size(); j++) {
                  backFaces += ebisBox.getFaces(otherVoFs[j], dir, flip(sit())).size();
                }
                if (backFaces < faces.size()) {
                  faceMismatch++;

                  if (reported < 0) {
                    reported++;

                    pout() << "   ASYM lvl " << lvl << " " << vof.gridIndex() << " dir " << dir << " side "
                           << (sit() == Side::Lo ? "lo" : "hi") << " faces " << faces.size() << " back " << backFaces
                           << " area " << area << " | this kappa " << ebisBox.volFrac(vof) << " bndry "
                           << ebisBox.bndryArea(vof) << " n " << ebisBox.normal(vof) << " | other " << other << " nvof "
                           << otherVoFs.size() << " kappa " << ebisBox.volFrac(otherVoFs[0]) << " bndry "
                           << ebisBox.bndryArea(otherVoFs[0]) << " n " << ebisBox.normal(otherVoFs[0]) << endl;
                  }
                }
              }
            }
          }

          const RealVect residual = apertureVector - ebisBox.bndryArea(vof) * ebisBox.normal(vof);
          const Real     r        = residual.vectorLength();

          worstDivergence = std::max(worstDivergence, r);
          if (r > 1.0E-9) {
            divergenceBad++;

            // What sits across each face the cell does not have. A full cell bordering a covered
            // one has to carry that face as its boundary; bordering a regular one it does not.
            for (int dir = 0; dir < SpaceDim; dir++) {
              for (SideIterator sit; sit.ok(); ++sit) {
                if (ebisBox.getFaces(vof, dir, sit()).size() > 0) {
                  continue;
                }

                const IntVect other = vof.gridIndex() + sign(sit()) * BASISV(dir);

                if (!ebisBox.getRegion().contains(other)) {
                  acrossOutside++;
                }
                else if (ebisBox.isCovered(other)) {
                  acrossCovered++;
                }
                else if (ebisBox.isRegular(other)) {
                  acrossRegular++;
                }
                else {
                  acrossIrregular++;
                }
              }
            }

            if (ebisBox.isRegular(vof.gridIndex())) {
              badFullCells++;
            }

            if (reported < 300) {
              reported++;

              // Is the cell at the edge of what this level carries?
              bool onLevelEdge = false;
              for (int dir = 0; dir < SpaceDim; dir++) {
                for (SideIterator sit; sit.ok(); ++sit) {
                  if (!levelCoverage.contains(vof.gridIndex() + sign(sit()) * BASISV(dir))) {
                    onLevelEdge = true;
                  }
                }
              }

              pout() << "   BAD lvl " << lvl << " " << vof.gridIndex() << " res " << r << " nvof " << vofs.size()
                     << " kappa " << ebisBox.volFrac(vof) << " bndryArea " << ebisBox.bndryArea(vof) << " normal "
                     << ebisBox.normal(vof) << " aperture " << apertureVector << " underFiner "
                     << finerCoverage.contains(vof.gridIndex()) << " onLevelEdge " << onLevelEdge << endl;
            }
          }
        }
      }
    }

    pout() << "VALIDATE lvl " << lvl << " dx " << dx[lvl] << ": faceMismatch " << faceMismatch << " openIntoCovered "
           << openIntoCovered << " divergenceBad " << divergenceBad << " worstDivergence " << worstDivergence << endl;

    if (divergenceBad > 0) {
      pout() << "VALIDATE lvl " << lvl << " of the " << divergenceBad << " bad cells, " << badFullCells
             << " are full cells; the faces they lack sit across " << acrossCovered << " covered, " << acrossRegular
             << " regular, " << acrossIrregular << " irregular and " << acrossOutside << " outside" << endl;
    }
  }

  // Coarsening has to conserve: a coarse cell holds what its fine cells hold.
  for (int lvl = 0; lvl < a_amr->getFinestLevel(); lvl++) {
    const int ratio = refRat[lvl];

    Real     worstVolume = 0.0;
    long int checked     = 0;

    for (DataIterator dit = grids[lvl + 1].dataIterator(); dit.ok(); ++dit) {
      const EBISBox&   fineBox = ebisl[lvl + 1][dit()];
      const IntVectSet irreg   = fineBox.getIrregIVS(grids[lvl + 1][dit()]);

      IntVectSet coarseSeen;
      for (IVSIterator ivsIt(irreg); ivsIt.ok(); ++ivsIt) {
        IntVect civ = ivsIt();
        civ.coarsen(ratio);
        coarseSeen |= civ;
      }

      for (IVSIterator ivsIt(coarseSeen); ivsIt.ok(); ++ivsIt) {
        const Box fineCells = refine(Box(ivsIt(), ivsIt()), ratio);

        Real fineVolume = 0.0;
        bool complete   = true;

        for (BoxIterator bit(fineCells); bit.ok(); ++bit) {
          if (!grids[lvl + 1][dit()].contains(bit())) {
            complete = false;
            break;
          }
          const Vector<VolIndex> fv = fineBox.getVoFs(bit());
          for (int j = 0; j < fv.size(); j++) {
            fineVolume += fineBox.volFrac(fv[j]);
          }
        }

        if (!complete) {
          continue;
        }

        fineVolume /= std::pow(Real(ratio), SpaceDim);

        // The coarse cell may live on another rank's box; only check where we hold it.
        for (DataIterator cdit = grids[lvl].dataIterator(); cdit.ok(); ++cdit) {
          if (!grids[lvl][cdit()].contains(ivsIt())) {
            continue;
          }
          const EBISBox&         coarBox    = ebisl[lvl][cdit()];
          Real                   coarVolume = 0.0;
          const Vector<VolIndex> cv         = coarBox.getVoFs(ivsIt());
          for (int j = 0; j < cv.size(); j++) {
            coarVolume += coarBox.volFrac(cv[j]);
          }
          worstVolume = std::max(worstVolume, std::abs(coarVolume - fineVolume));
          checked++;
        }
      }
    }

    pout() << "VALIDATE coarsening " << lvl + 1 << "->" << lvl << ": cells " << checked << " worstVolumeMismatch "
           << worstVolume << endl;
  }
}

void
reportCoverage(const RefCountedPtr<AmrMesh>& a_amr, const RefCountedPtr<ComputationalGeometry>& a_compgeom)
{
  const RefCountedPtr<EBIndexSpace>& ebis = a_compgeom->getMfIndexSpace()->getEBIndexSpace(phase::gas);

  const Vector<DisjointBoxLayout>& grids   = a_amr->getGrids(Realm::primal);
  const Vector<ProblemDomain>&     domains = a_amr->getDomains();

  for (int lvl = 0; lvl <= a_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& ebisGrids = ebis->getGrids(domains[lvl]);

    IntVectSet covered;
    for (int i = 0; i < ebisGrids.boxArray().size(); i++) {
      covered |= ebisGrids.boxArray()[i];
    }

    long int inside  = 0;
    long int outside = 0;

    for (int i = 0; i < grids[lvl].boxArray().size(); i++) {
      if (covered.contains(grids[lvl].boxArray()[i])) {
        inside++;
      }
      else {
        outside++;
      }
    }

    pout() << "COVERAGE level " << lvl << ": amr boxes " << (inside + outside) << ", inside the index space " << inside
           << ", outside " << outside << ", index space boxes " << ebisGrids.boxArray().size() << endl;
  }
}

void
exportAggregationFailures(const RefCountedPtr<AmrMesh>& a_amr, const int a_depth, const std::string& a_fileName)
{
  const RefCountedPtr<BaseIF>& implicitFunction = a_amr->getBaseImplicitFunction(phase::gas);
  const Vector<Real>&          dx               = a_amr->getDx();
  const RealVect               probLo           = a_amr->getProbLo();
  const ProblemDomain&         domain           = a_amr->getDomains()[0];

  std::ofstream out(a_fileName);
  out << std::setprecision(10);

  constexpr int numChildren = 1 << SpaceDim;

  long int failures = 0;

  for (BoxIterator bit(domain.domainBox()); bit.ok(); ++bit) {
    const IntVect iv = bit();

    RealVect centre = probLo;
    for (int d = 0; d < SpaceDim; d++) {
      centre[d] += dx[0] * (static_cast<Real>(iv[d]) + 0.5);
    }

    if (std::abs(implicitFunction->value(centre)) > dx[0] * std::sqrt(1.0 * SpaceDim)) {
      continue;
    }

    PolyhedralEB::CutCellBody body;

    if (buildAggregated(*implicitFunction, body, iv, probLo, dx[0], a_depth)) {
      continue;
    }

    failures++;

    if (failures > 3) {
      continue;
    }

    out << "cell " << iv << " dx " << dx[0] << " divergence " << body.divergenceResidual() << " volFrac "
        << body.volumeFraction() << " aB " << body.boundaryArea() << " trueA " << body.trueBoundaryArea() << " normal "
        << body.normal() << "\n";

    for (int d = 0; d < SpaceDim; d++) {
      out << "   parent aperture dir " << d << " lo " << body.areaFraction(d, Side::Lo) << " hi "
          << body.areaFraction(d, Side::Hi) << "\n";
    }

    for (int c = 0; c < numChildren; c++) {
      IntVect fine = 2 * iv;
      for (int d = 0; d < SpaceDim; d++) {
        fine[d] += (c >> d) & 1;
      }

      PolyhedralEB::CutCellBody child;
      const bool                ok = buildAggregated(*implicitFunction, child, fine, probLo, 0.5 * dx[0], a_depth - 1);

      out << "   child " << c << " built " << ok << " kind " << static_cast<int>(child.kind()) << " volFrac "
          << child.volumeFraction() << " aB " << child.boundaryArea() << " divergence " << child.divergenceResidual();
      for (int d = 0; d < SpaceDim; d++) {
        out << " | dir" << d << " lo " << child.areaFraction(d, Side::Lo) << " hi " << child.areaFraction(d, Side::Hi);
      }
      out << "\n";
    }
  }

  out << "total failures " << failures << "\n";
}

void
exportSheetCancellation(const RefCountedPtr<AmrMesh>& a_amr, const int a_maxDepth, const std::string& a_fileName)
{
  const RefCountedPtr<BaseIF>& implicitFunction = a_amr->getBaseImplicitFunction(phase::gas);
  const Vector<Real>&          dx               = a_amr->getDx();
  const RealVect               probLo           = a_amr->getProbLo();
  const ProblemDomain&         domain           = a_amr->getDomains()[0];

  std::ofstream out(a_fileName);
  out << std::setprecision(17);
  out << "depth,cells,cancelling,worstRatio,worstVolFrac,worstBndryArea,worstTrueArea,worstResidual,multiSheet\n";

  for (int depth = 0; depth <= a_maxDepth; depth++) {
    long int cells      = 0;
    long int cancelling = 0;
    long int multiSheet = 0;

    Real worstRatio    = 0.0;
    Real worstVol      = 0.0;
    Real worstBndry    = 0.0;
    Real worstTrue     = 0.0;
    Real worstResidual = 0.0;

    for (BoxIterator bit(domain.domainBox()); bit.ok(); ++bit) {
      const IntVect iv = bit();

      // Only cells the interface could reach are worth building.
      RealVect centre = probLo;
      for (int d = 0; d < SpaceDim; d++) {
        centre[d] += dx[0] * (static_cast<Real>(iv[d]) + 0.5);
      }

      if (std::abs(implicitFunction->value(centre)) > dx[0] * std::sqrt(1.0 * SpaceDim)) {
        continue;
      }

      PolyhedralEB::CutCellBody body;

      if (!buildAggregated(*implicitFunction, body, iv, probLo, dx[0], depth)) {
        continue;
      }

      if (body.kind() != PolyhedralEB::CutCellBody::Kind::Cut) {
        continue;
      }

      cells++;

      // Two sheets facing opposite ways cancel in the sum while both count in the total.
      const Real ratio = (body.boundaryArea() > 0.0) ? (body.trueBoundaryArea() / body.boundaryArea())
                                                     : ((body.trueBoundaryArea() > 0.0) ? 1.0E30 : 1.0);

      if (ratio > 1.5) {
        cancelling++;
      }

      // How many separate sheets of interface the cell holds, at the resolution the body is
      // built on. More than one is a cell the single-valued construction cannot honour.
#if CH_SPACEDIM == 3
      {
        const PolyhedralEB::CutCellSurface surface = sampleSurface(*implicitFunction, iv, probLo, dx[0], true);

        int loop[PolyhedralEB::CutCellSurface::s_numEdges];
        int start[PolyhedralEB::CutCellSurface::s_numEdges + 1];

        if (PolyhedralEB::detail::crossingLoops(surface, loop, start) > 1) {
          multiSheet++;
        }
      }
#endif

      if (ratio > worstRatio) {
        worstRatio    = ratio;
        worstVol      = body.volumeFraction();
        worstBndry    = body.boundaryArea();
        worstTrue     = body.trueBoundaryArea();
        worstResidual = body.divergenceResidual();
      }
    }

    out << depth << "," << cells << "," << cancelling << "," << worstRatio << "," << worstVol << "," << worstBndry
        << "," << worstTrue << "," << worstResidual << "," << multiSheet << "\n";
  }
}

void
exportAggregationTotals(const RefCountedPtr<AmrMesh>& a_amr, const int a_maxDepth, const std::string& a_fileName)
{
  const RefCountedPtr<BaseIF>&     implicitFunction = a_amr->getBaseImplicitFunction(phase::gas);
  const Vector<DisjointBoxLayout>& grids            = a_amr->getGrids(Realm::primal);
  const Vector<EBISLayout>&        ebisl            = a_amr->getEBISLayout(Realm::primal, phase::gas);
  const Vector<Real>&              dx               = a_amr->getDx();
  const RealVect                   probLo           = a_amr->getProbLo();

  std::ofstream out(a_fileName);
  out << std::setprecision(17);
  out << "depth,dx,cells,cutVolume,boundaryArea,trueBoundaryArea\n";

  const int lvl = 0;

  for (int depth = 0; depth <= a_maxDepth; depth++) {
    Real     volume   = 0.0;
    Real     boundary = 0.0;
    Real     patch    = 0.0;
    long int cells    = 0;

    for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
      const EBISBox&   ebisBox  = ebisl[lvl][dit()];
      const IntVectSet irregIVS = ebisBox.getIrregIVS(grids[lvl][dit()]);

      for (IVSIterator ivsIt(irregIVS); ivsIt.ok(); ++ivsIt) {
        PolyhedralEB::CutCellBody body;

        if (!buildAggregated(*implicitFunction, body, ivsIt(), probLo, dx[lvl], depth)) {
          continue;
        }

        cells++;
        volume += body.volumeFraction();
        boundary += body.boundaryArea();
        patch += body.trueBoundaryArea();
      }
    }

    const Real cellVolume = std::pow(dx[lvl], SpaceDim);
    const Real cellArea   = std::pow(dx[lvl], SpaceDim - 1);

    out << depth << "," << dx[lvl] << "," << cells << "," << volume * cellVolume << "," << boundary * cellArea << ","
        << patch * cellArea << "\n";
  }
}

void
exportCoarseningSeam(const RefCountedPtr<AmrMesh>& a_amr, const std::string& a_fileName)
{
  const RefCountedPtr<BaseIF>&     implicitFunction = a_amr->getBaseImplicitFunction(phase::gas);
  const Vector<DisjointBoxLayout>& grids            = a_amr->getGrids(Realm::primal);
  const Vector<EBISLayout>&        ebisl            = a_amr->getEBISLayout(Realm::primal, phase::gas);
  const Vector<Real>&              dx               = a_amr->getDx();
  const RealVect                   probLo           = a_amr->getProbLo();

  std::ofstream out(a_fileName);
  out << std::setprecision(17);
  out << "level,dx,volumeCoarse,volumeFine,boundaryCoarse,boundaryFine,worstAperture,roundTrip\n";

  constexpr int numChildren = 1 << SpaceDim;
  constexpr int numPerFace  = 1 << (SpaceDim - 1);

  for (int lvl = 0; lvl <= a_amr->getFinestLevel(); lvl++) {
    for (DataIterator dit = grids[lvl].dataIterator(); dit.ok(); ++dit) {
      const EBISBox&   ebisBox  = ebisl[lvl][dit()];
      const IntVectSet irregIVS = ebisBox.getIrregIVS(grids[lvl][dit()]);

      for (IVSIterator ivsIt(irregIVS); ivsIt.ok(); ++ivsIt) {
        const IntVect& iv = ivsIt();

        PolyhedralEB::CutCellBody coarse;

        if (!coarse.define(sampleSurface(*implicitFunction, iv, probLo, dx[lvl], true))) {
          continue;
        }

        PolyhedralEB::CutCellBody child[numChildren];

        bool built = true;

        for (int c = 0; c < numChildren && built; c++) {
          IntVect fine = 2 * iv;

          for (int d = 0; d < SpaceDim; d++) {
            fine[d] += (c >> d) & 1;
          }

          built = child[c].define(sampleSurface(*implicitFunction, fine, probLo, 0.5 * dx[lvl], true));
        }

        if (!built) {
          continue;
        }

        Real volumeFine   = 0.0;
        Real boundaryFine = 0.0;

        for (int c = 0; c < numChildren; c++) {
          volumeFine += child[c].volumeFraction();
          boundaryFine += child[c].boundaryArea();
        }

        volumeFine /= static_cast<Real>(numChildren);
        boundaryFine /= static_cast<Real>(numPerFace);

        Real worst = 0.0;

        for (int d = 0; d < SpaceDim; d++) {
          for (SideIterator sit; sit.ok(); ++sit) {
            const int bit = (sit() == Side::Lo) ? 0 : 1;

            Real fine = 0.0;

            for (int c = 0; c < numChildren; c++) {
              if (((c >> d) & 1) == bit) {
                fine += child[c].areaFraction(d, sit());
              }
            }

            fine /= static_cast<Real>(numPerFace);

            worst = std::max(worst, std::abs(fine - coarse.areaFraction(d, sit())));
          }
        }

        // refine then coarsen has to give this body back: the children partition it, so the
        // moments they are summed from are the ones they were cut from.
        Real roundTrip = -1.0;

        PolyhedralEB::CutCellBody split[numChildren];

        if (coarse.refine(split)) {
          PolyhedralEB::CutCellBody rebuilt;

          rebuilt.coarsen(split);

          roundTrip = std::abs(rebuilt.volumeFraction() - coarse.volumeFraction());
          roundTrip = std::max(roundTrip, std::abs(rebuilt.boundaryArea() - coarse.boundaryArea()));

          for (int d = 0; d < SpaceDim; d++) {
            for (SideIterator sit; sit.ok(); ++sit) {
              roundTrip = std::max(roundTrip, std::abs(rebuilt.areaFraction(d, sit()) - coarse.areaFraction(d, sit())));
            }

            roundTrip = std::max(roundTrip, std::abs(rebuilt.volumeCentroid()[d] - coarse.volumeCentroid()[d]));
            roundTrip = std::max(roundTrip, std::abs(rebuilt.normal()[d] - coarse.normal()[d]));
            roundTrip = std::max(roundTrip, std::abs(rebuilt.boundaryCentroid()[d] - coarse.boundaryCentroid()[d]));
          }
        }

        out << lvl << "," << dx[lvl] << "," << coarse.volumeFraction() << "," << volumeFine << ","
            << coarse.boundaryArea() << "," << boundaryFine << "," << worst << "," << roundTrip << "\n";
      }
    }
  }
}

void
exportCutCells(const RefCountedPtr<AmrMesh>& a_amr, const std::string& a_fileName)
{
  const RefCountedPtr<BaseIF>&     implicitFunction = a_amr->getBaseImplicitFunction(phase::gas);
  const Vector<DisjointBoxLayout>& grids            = a_amr->getGrids(Realm::primal);
  const Vector<EBISLayout>&        ebisl            = a_amr->getEBISLayout(Realm::primal, phase::gas);
  const Vector<Real>&              dx               = a_amr->getDx();
  const RealVect                   probLo           = a_amr->getProbLo();

  std::ofstream out(a_fileName);
  out << std::setprecision(17);

  // Header. Face and edge blocks are emitted in a fixed order so the reader can index them.
  out << "level,dx";
  for (int d = 0; d < SpaceDim; d++) {
    out << ",iv" << d;
  }
  out << ",numVoFs,volFrac";
  for (int d = 0; d < SpaceDim; d++) {
    out << ",volCentroid" << d;
  }
  out << ",bndryArea";
  for (int d = 0; d < SpaceDim; d++) {
    out << ",normal" << d;
  }
  for (int d = 0; d < SpaceDim; d++) {
    out << ",bndryCentroid" << d;
  }
  out << ",kappa,kappaDx";
  for (int d = 0; d < SpaceDim; d++) {
    out << ",xNormal" << d;
  }
  for (int d = 0; d < SpaceDim; d++) {
    out << ",gNormal" << d;
  }
  for (int d = 0; d < SpaceDim; d++) {
    out << ",lNormal" << d;
  }
  for (int face = 0; face < 2 * SpaceDim; face++) {
    out << ",numFaces" << face << ",areaFrac" << face;
    for (int d = 0; d < SpaceDim; d++) {
      out << ",faceCentroid" << face << "_" << d;
    }
  }
  for (int edge = 0; edge < numCellEdges; edge++) {
    out << ",edgeCut" << edge;
  }
  for (int corner = 0; corner < (1 << SpaceDim); corner++) {
    out << ",corner" << corner;
  }
  for (int d = 0; d < SpaceDim; d++) {
    out << ",apex" << d;
  }
  for (int face = 0; face < 2 * SpaceDim; face++) {
    out << ",kink" << face;
    for (int b = 0; b < 2; b++) {
      for (int d = 0; d < SpaceDim; d++) {
        out << ",kink" << face << "_" << b << "_" << d;
      }
    }
  }
  out << "\n";

  for (int lvl = 0; lvl <= a_amr->getFinestLevel(); lvl++) {
    const DisjointBoxLayout& dbl = grids[lvl];
    const DataIterator&      dit = dbl.dataIterator();

    for (int mybox = 0; mybox < dit.size(); mybox++) {
      const DataIndex& din     = dit[mybox];
      const Box&       box     = dbl[din];
      const EBISBox&   ebisBox = ebisl[lvl][din];

      const IntVectSet irregIVS = ebisBox.getIrregIVS(box);

      for (IVSIterator ivsIt(irregIVS); ivsIt.ok(); ++ivsIt) {
        const IntVect& iv = ivsIt();

        const Vector<VolIndex> vofs = ebisBox.getVoFs(iv);

        // Multi-valued cells are outside the scope of the refinement operator, but they are
        // counted rather than dropped so the census can report how often they occur.
        const VolIndex& vof = vofs[0];

        out << lvl << "," << dx[lvl];
        for (int d = 0; d < SpaceDim; d++) {
          out << "," << iv[d];
        }
        out << "," << vofs.size();
        out << "," << ebisBox.volFrac(vof);

        const RealVect volCentroid = ebisBox.centroid(vof);
        for (int d = 0; d < SpaceDim; d++) {
          out << "," << volCentroid[d];
        }

        out << "," << ebisBox.bndryArea(vof);

        const RealVect normal = ebisBox.normal(vof);
        for (int d = 0; d < SpaceDim; d++) {
          out << "," << normal[d];
        }

        const RealVect bndryCentroid = ebisBox.bndryCentroid(vof);
        for (int d = 0; d < SpaceDim; d++) {
          out << "," << bndryCentroid[d];
        }

        // Measured on the surface itself, at a lattice spacing of one cell width: the bend of
        // the surface across this cell rather than its pointwise differential curvature.
        {
          RealVect x = probLo;

          for (int d = 0; d < SpaceDim; d++) {
            x[d] += dx[lvl] * (static_cast<Real>(iv[d]) + 0.5 + bndryCentroid[d]);
          }

          x = projectToSurface(*implicitFunction, x, dx[lvl]);

          Real stencil[PolyhedralEB::s_curvatureStencilSize];
          curvatureStencil(*implicitFunction, x, dx[lvl], stencil);

          const Real kappa = PolyhedralEB::maxPrincipalCurvature(stencil, dx[lvl]);

          out << "," << kappa << "," << kappa * dx[lvl];
        }

        {
          const PolyhedralEB::CutCellSurface surface = sampleSurface(*implicitFunction, iv, probLo, dx[lvl], true);

          const RealVect crossingN = PolyhedralEB::crossingNormal(surface);
          const RealVect gradientN = PolyhedralEB::gradientNormal(surface);

          for (int d = 0; d < SpaceDim; d++) {
            out << "," << crossingN[d];
          }
          for (int d = 0; d < SpaceDim; d++) {
            out << "," << gradientN[d];
          }

          const PolyhedralEB::CutCellSurface linear  = sampleSurface(*implicitFunction, iv, probLo, dx[lvl], false);
          const RealVect                     linearN = PolyhedralEB::crossingNormal(linear);

          for (int d = 0; d < SpaceDim; d++) {
            out << "," << linearN[d];
          }
        }

        for (int dir = 0; dir < SpaceDim; dir++) {
          for (SideIterator sit; sit.ok(); ++sit) {
            const Vector<FaceIndex> faces = ebisBox.getFaces(vof, dir, sit());

            const int faceIndex = 2 * dir + ((sit() == Side::Lo) ? 0 : 1);
            (void)faceIndex;

            out << "," << faces.size();

            if (faces.size() == 1) {
              out << "," << ebisBox.areaFrac(faces[0]);

              const RealVect faceCentroid = ebisBox.centroid(faces[0]);
              for (int d = 0; d < SpaceDim; d++) {
                out << "," << faceCentroid[d];
              }
            }
            else {
              out << "," << 0.0;
              for (int d = 0; d < SpaceDim; d++) {
                out << "," << 0.0;
              }
            }
          }
        }

        // The true crossings, from the implicit function itself. These are the anchors a
        // triangulation would use, and the reference the recovery is judged against.
        RealVect crossingSum  = RealVect::Zero;
        int      numCrossings = 0;

        for (int edge = 0; edge < numCellEdges; edge++) {
          int     edgeDir;
          IntVect loOffset;
          edgeGeometry(edge, edgeDir, loOffset);

          IntVect hiOffset = loOffset;
          hiOffset[edgeDir] += 1;

          RealVect loPt;
          RealVect hiPt;
          for (int d = 0; d < SpaceDim; d++) {
            loPt[d] = probLo[d] + dx[lvl] * (iv[d] + loOffset[d]);
            hiPt[d] = probLo[d] + dx[lvl] * (iv[d] + hiOffset[d]);
          }

          const Real t = edgeCrossing(*implicitFunction, loPt, hiPt);

          out << "," << t;

          if (t >= 0.0) {
            crossingSum += loPt + t * (hiPt - loPt);
            numCrossings++;
          }
        }

        // Corner values of the implicit function. These fix which side of the interface is fluid,
        // and are shared by every cell meeting at the corner.
        for (int corner = 0; corner < (1 << SpaceDim); corner++) {
          RealVect pt;
          for (int d = 0; d < SpaceDim; d++) {
            pt[d] = probLo[d] + dx[lvl] * (iv[d] + ((corner >> d) & 1));
          }

          out << "," << implicitFunction->value(pt);
        }

        // The interface patch is anchored here rather than at the centroid of the crossings,
        // which for a cell holding a crease floats off the surface entirely.
        RealVect apex = RealVect::Zero;
        if (numCrossings > 0) {
          const RealVect centre = crossingSum / Real(numCrossings);
          const RealVect onSurf = projectToSurface(*implicitFunction, centre, dx[lvl]);

          for (int d = 0; d < SpaceDim; d++) {
            const Real rel = (onSurf[d] - probLo[d]) / dx[lvl] - Real(iv[d]) - 0.5;

            apex[d] = std::max(-0.5, std::min(0.5, rel));
          }
        }

        for (int d = 0; d < SpaceDim; d++) {
          out << "," << apex[d];
        }

        // Where the interface bends inside each face. Both cells sharing a face compute this
        // from the same face and the same implicit function, so they cannot disagree.
        for (int faceDir = 0; faceDir < SpaceDim; faceDir++) {
          for (int faceSide = 0; faceSide < 2; faceSide++) {
            std::vector<RealVect> ends;

            for (int edge = 0; edge < numCellEdges; edge++) {
              int     edgeDir;
              IntVect loOffset;
              edgeGeometry(edge, edgeDir, loOffset);

              if (edgeDir == faceDir || loOffset[faceDir] != faceSide) {
                continue;
              }

              RealVect loPt;
              RealVect hiPt;
              for (int d = 0; d < SpaceDim; d++) {
                loPt[d] = probLo[d] + dx[lvl] * (iv[d] + loOffset[d]);
                hiPt[d] = loPt[d];
              }
              hiPt[edgeDir] += dx[lvl];

              const Real t = edgeCrossing(*implicitFunction, loPt, hiPt);
              if (t >= 0.0) {
                ends.push_back(loPt + t * (hiPt - loPt));
              }
            }

            RealVect bend[2];
            int      numBends = 0;
            if (ends.size() == 2) {
              // order the endpoints by edge index so both cells sharing the face agree on
              // which way along the contour the vertices are listed
              numBends = faceBends(*implicitFunction, ends[0], ends[1], faceDir, dx[lvl], bend);

              // a vertex outside the face is not on this contour: the perpendicular search can
              // find a root beyond the face, and inserting it makes the face polygon a bowtie
              for (int b = 0; b < numBends; b++) {
                for (int d = 0; d < SpaceDim; d++) {
                  const Real rel = (bend[b][d] - probLo[d]) / dx[lvl] - Real(iv[d]) - 0.5;

                  if (rel < -0.5 - 1.0E-9 || rel > 0.5 + 1.0E-9) {
                    numBends = 0;
                  }
                }
              }
            }

            out << "," << numBends;
            for (int b = 0; b < 2; b++) {
              for (int d = 0; d < SpaceDim; d++) {
                const Real rel = (b < numBends) ? (bend[b][d] - probLo[d]) / dx[lvl] - Real(iv[d]) - 0.5 : 0.0;

                out << "," << rel;
              }
            }
          }
        }

        out << "\n";
      }
    }
  }

  out.close();
}

int
main(int argc, char* argv[])
{
  ChomboDischarge::initialize(argc, argv);

  // the census is run on a sphere (smooth, known curvature) and on a tessellation (piecewise
  // planar, so the interface inside a triangle is exactly a plane)
  std::string whichGeom = "rough_sphere";
  {
    ParmParse pp("Prototype");
    pp.query("geometry", whichGeom);
  }

  // Reads a RealVect from the input file, leaving the supplied value in place when the key is absent.
  auto queryVect = [](const char* a_prefix, const char* a_key, RealVect& a_value) {
    ParmParse         pp(a_prefix);
    std::vector<Real> v(SpaceDim);
    for (int d = 0; d < SpaceDim; d++) {
      v[d] = a_value[d];
    }

    pp.queryarr(a_key, v, 0, SpaceDim);

    for (int d = 0; d < SpaceDim; d++) {
      a_value[d] = v[d];
    }
  };

  // A plane at a configurable angle: the interface inside every cut cell is exactly planar, so any
  // residual is the reconstruction's own error rather than a curvature effect. Sweeping the normal
  // towards a coordinate axis drives the cut towards tangency with the cell faces.
  class TiltedPlane : public ComputationalGeometry
  {
  public:
    TiltedPlane(const RealVect& a_normal, const RealVect& a_point)
    {
      RefCountedPtr<BaseIF> plane = RefCountedPtr<BaseIF>(new PlaneIF(a_normal, a_point, true));

      m_electrodes.push_back(Electrode(plane, true));
    }
  };

  // A tilted plate of a set thickness. Below one cell thick every corner of every cell it passes
  // through lands on the same side of it, so the plate is invisible to the corner signs and to
  // anything built from them -- which is the case a bend criterion cannot see.
  class Slab : public ComputationalGeometry
  {
  public:
    Slab(const RealVect& a_normal, const RealVect& a_point, const Real a_thickness)
    {
      RealVect unit = a_normal;
      unit /= unit.vectorLength();

      Vector<BaseIF*> planes;

      // An electrode's implicit function is negative inside the solid, and IntersectionIF takes
      // the larger of the two, so each plane is oriented negative on the slab's side of itself.
      planes.push_back(new PlaneIF(unit, a_point - 0.5 * a_thickness * unit, true));
      planes.push_back(new PlaneIF(unit, a_point + 0.5 * a_thickness * unit, false));

      RefCountedPtr<BaseIF> slab = RefCountedPtr<BaseIF>(new IntersectionIF(planes));

      m_electrodes.push_back(Electrode(slab, true));

      for (int i = 0; i < planes.size(); i++) {
        delete planes[i];
      }
    }
  };

  // A sphere at a configurable centre. Shifting the centre by a fraction of a cell sweeps the phase
  // of the surface relative to the grid, which is what exposes tangency and near-degenerate cuts.
  class Sphere : public ComputationalGeometry
  {
  public:
    Sphere(const RealVect& a_center, const Real a_radius)
    {
      RefCountedPtr<BaseIF> sphere = RefCountedPtr<BaseIF>(new SphereSdf(a_center, a_radius, false));

      m_electrodes.push_back(Electrode(sphere, true));
    }
  };

  // A sphere swept along a segment. The barrel has one vanishing principal curvature and one equal
  // to the inverse radius, so unlike the sphere it never has two equal curvatures, and the segment
  // can be laid along a coordinate axis or along a diagonal.
  class SweptSphere : public ComputationalGeometry
  {
  public:
    SweptSphere(const RealVect& a_center1, const RealVect& a_center2, const Real a_radius)
    {
      RefCountedPtr<BaseIF> cyl = RefCountedPtr<BaseIF>(new CylinderSdf(a_center1, a_center2, a_radius, false));

      m_electrodes.push_back(Electrode(cyl, true));
    }
  };

  // A torus: the two principal curvatures differ everywhere, and the inner equator is a saddle
  // (Gaussian curvature negative), which an isotropic paraboloid cannot represent at all.
  class Torus : public ComputationalGeometry
  {
  public:
    Torus()
    {
      Real major = 0.5;
      Real minor = 0.15;
      {
        ParmParse pp("Torus");
        pp.query("major_radius", major);
        pp.query("minor_radius", minor);
      }

      RefCountedPtr<BaseIF> torus = RefCountedPtr<BaseIF>(new TorusSdf(RealVect::Zero, major, minor, false));

      m_electrodes.push_back(Electrode(torus, true));
    }
  };

  // A rotated cube or simplex: piecewise planar with sharp edges and corners, which is what a
  // tessellated surface looks like inside a cell once the facets are larger than the mesh.
  class Polyhedron : public ComputationalGeometry
  {
  public:
    Polyhedron(const std::string& a_shape, const RealVect& a_angles, const Real a_size, const RealVect& a_center)
    {
      Vector<RealVect> normals;
      Vector<Real>     offsets;

      if (a_shape == "simplex") {
#if CH_SPACEDIM == 2
        for (int i = 0; i < 3; i++) {
          const Real t = 2.0 * M_PI * i / 3.0;
          normals.push_back(RealVect(D_DECL(cos(t), sin(t), 0.0)));
        }
#else
        const Real r = 1.0 / sqrt(3.0);
        normals.push_back(RealVect(D_DECL(r, r, r)));
        normals.push_back(RealVect(D_DECL(r, -r, -r)));
        normals.push_back(RealVect(D_DECL(-r, r, -r)));
        normals.push_back(RealVect(D_DECL(-r, -r, r)));
#endif
      }
      else {
        for (int d = 0; d < SpaceDim; d++) {
          for (int s = 0; s < 2; s++) {
            RealVect n = RealVect::Zero;
            n[d]       = (s == 0) ? -1.0 : 1.0;
            normals.push_back(n);
          }
        }
      }

      for (int i = 0; i < normals.size(); i++) {
        normals[i] = rotate(normals[i], a_angles);
        offsets.push_back(a_size);
      }

      RefCountedPtr<BaseIF> body = RefCountedPtr<BaseIF>(new ConvexBody(normals, offsets, a_center));

      m_electrodes.push_back(Electrode(body, true));
    }
  };

  RefCountedPtr<ComputationalGeometry> compgeom;
  if (whichGeom == "polyhedron") {
    std::string shape  = "cube";
    RealVect    angles = RealVect(D_DECL(31.7, 19.3, 47.1));
    RealVect    center = RealVect::Zero;
    Real        size   = 0.4;
    {
      ParmParse pp("Polyhedron");
      pp.query("shape", shape);
      pp.query("size", size);
    }
    queryVect("Polyhedron", "angles", angles);
    queryVect("Polyhedron", "center", center);
    compgeom = RefCountedPtr<ComputationalGeometry>(new Polyhedron(shape, angles, size, center));
  }
  else if (whichGeom == "sphere") {
    RealVect center = RealVect::Zero;
    Real     radius = 0.25;
    queryVect("Sphere", "center", center);
    {
      ParmParse pp("Sphere");
      pp.query("radius", radius);
    }
    compgeom = RefCountedPtr<ComputationalGeometry>(new Sphere(center, radius));
  }
  else if (whichGeom == "swept_sphere") {
    RealVect center1 = -0.5 * RealVect::Unit;
    RealVect center2 = 0.5 * RealVect::Unit;
    Real     radius  = 0.25;
    queryVect("SweptSphere", "center1", center1);
    queryVect("SweptSphere", "center2", center2);
    {
      ParmParse pp("SweptSphere");
      pp.query("radius", radius);
    }
    compgeom = RefCountedPtr<ComputationalGeometry>(new SweptSphere(center1, center2, radius));
  }
  else if (whichGeom == "torus") {
    compgeom = RefCountedPtr<ComputationalGeometry>(new Torus());
  }
  else if (whichGeom == "plane") {
    RealVect normal = RealVect(D_DECL(0.3721, -0.5839, 0.7214));
    RealVect point  = RealVect(D_DECL(0.013, -0.021, 0.037));
    queryVect("Plane", "normal", normal);
    queryVect("Plane", "point", point);
    compgeom = RefCountedPtr<ComputationalGeometry>(new TiltedPlane(normal, point));
  }
  else if (whichGeom == "slab") {
    RealVect normal    = RealVect(D_DECL(0.3721, -0.5839, 0.7214));
    RealVect point     = RealVect(D_DECL(0.013, -0.021, 0.037));
    Real     thickness = 0.01;
    queryVect("Slab", "normal", normal);
    queryVect("Slab", "point", point);
    {
      ParmParse pp("Slab");
      pp.query("thickness", thickness);
    }
    compgeom = RefCountedPtr<ComputationalGeometry>(new Slab(normal, point, thickness));
  }
  else if (whichGeom == "tessellation") {
    compgeom = RefCountedPtr<ComputationalGeometry>(new Tessellation());
  }
  else {
    compgeom = RefCountedPtr<ComputationalGeometry>(new RoughSphere());
  }
  auto amr         = RefCountedPtr<AmrMesh>(new AmrMesh());
  auto tagger      = RefCountedPtr<CellTagger>(nullptr);
  auto timestepper = RefCountedPtr<GeometryStepper>(new GeometryStepper());
  auto engine      = RefCountedPtr<Driver>(new Driver(compgeom, timestepper, amr, tagger));

  engine->setupAndRun();

  reportCoverage(amr, compgeom);

  validateIndexSpace(amr, compgeom);

  {
    int      coarse = 0;
    int      split  = 0;
    int      sweep  = 0;
    int      write  = 0;
    int      passB  = 0;
    Real     span   = 90.0;
    Real     size   = 0.25;
    RealVect center = RealVect::Zero;
    {
      ParmParse pp("Prototype");
      pp.query("twolevel_coarse", coarse);
      pp.query("twolevel_split", split);
      pp.query("twolevel_sweep", sweep);
      pp.query("twolevel_span", span);
      pp.query("twolevel_size", size);
      pp.query("twolevel_write", write);
      pp.query("twolevel_passb", passB);
    }

    queryVect("Prototype", "twolevel_center", center);

    // split 0 leaves only the fine block and split == coarse only the coarse one, which is how the
    // seam is told apart from what the triangulation does at a single resolution
    if (coarse > 0 && split >= 0) {
      if (sweep > 0) {
        sweepTwoLevelSeam(amr, coarse, split, sweep, span, size, center, write > 0, passB > 0);
      }
      else {
        validateTwoLevelSeam(compgeom, amr, coarse, split);
      }
    }
  }

  {
    int cells = 32;
    {
      ParmParse pp("Prototype");
      pp.query("seamface_cells", cells);
    }

    validateSeamFace(compgeom, amr, cells);
  }

  {
    int cells = 32;
    {
      ParmParse pp("Prototype");
      pp.query("multichord_cells", cells);
    }

    validateMultichordFace(compgeom, amr, cells);
  }

  {
    int depth         = 2;
    int coarseBoxSize = 8;
    {
      ParmParse pp("Prototype");
      pp.query("refined_fill_depth", depth);
      pp.query("refined_fill_box", coarseBoxSize);
    }

    validateRefinedFill(compgeom, amr, depth, coarseBoxSize);
  }

  char fileName[256];
  snprintf(fileName, sizeof(fileName), "cutcells.%dd.%d.csv", SpaceDim, procID());
  {
    Real angle    = 15.0;
    int  growth   = 0;
    int  maxDepth = 3;
    {
      ParmParse pp("Prototype");
      pp.query("curvature_angle", angle);
      pp.query("curvature_growth", growth);
      pp.query("curvature_max_depth", maxDepth);
    }

    char tagFile[256];
    snprintf(tagFile, sizeof(tagFile), "curvaturetags.%dd.%d.csv", SpaceDim, procID());

    exportCurvatureTags(compgeom, amr, angle, growth, maxDepth, std::string(tagFile));
  }

  {
    char seamFile[256];
    snprintf(seamFile, sizeof(seamFile), "seam.%dd.%d.csv", SpaceDim, procID());

    exportCoarseningSeam(amr, std::string(seamFile));
  }

  {
    int aggregationDepth = 0;
    {
      ParmParse pp("Prototype");
      pp.query("aggregation_depth", aggregationDepth);
    }

    if (aggregationDepth > 0) {
      char aggFile[256];
      snprintf(aggFile, sizeof(aggFile), "aggregation.%dd.%d.csv", SpaceDim, procID());

      exportAggregationTotals(amr, aggregationDepth, std::string(aggFile));

      char sheetFile[256];
      snprintf(sheetFile, sizeof(sheetFile), "sheets.%dd.%d.csv", SpaceDim, procID());

      exportSheetCancellation(amr, aggregationDepth, std::string(sheetFile));

      char failFile[256];
      snprintf(failFile, sizeof(failFile), "aggfail.%dd.%d.csv", SpaceDim, procID());

      exportAggregationFailures(amr, aggregationDepth, std::string(failFile));
    }
  }

  exportCutCells(amr, std::string(fileName));

  ChomboDischarge::finalize();
}
