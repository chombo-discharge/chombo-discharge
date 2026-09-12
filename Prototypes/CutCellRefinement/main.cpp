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

  const Vector<IntVectSet> tags = a_compgeom->getCurvatureTags(a_amr->getDomains()[0],
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
    }
  }

  exportCutCells(amr, std::string(fileName));

  ChomboDischarge::finalize();
}
