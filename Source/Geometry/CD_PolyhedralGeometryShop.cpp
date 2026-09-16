/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PolyhedralGeometryShop.cpp
 * @brief  Implementation of CD_PolyhedralGeometryShop.H
 * @author Robert Marskar
 */

// Std includes
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>

// Chombo includes
#include <BoxIterator.H>
#include <CH_assert.H>
#include <IntVectSet.H>
#include <LoHiSide.H>
#include <MayDay.H>

// Our includes
#include <CD_PolyhedralEBUtils.H>
#include <CD_PolyhedralGeometryShop.H>
#include <CD_NamespaceHeader.H>

PolyhedralGeometryShop::PolyhedralGeometryShop(const BaseIF&        a_localGeom,
                                               const int            a_verbosity,
                                               const Real           a_dx,
                                               const RealVect&      a_probLo,
                                               const ProblemDomain& a_finestDomain,
                                               const ProblemDomain& a_scanLevel,
                                               const int            a_ebGhost,
                                               const Real           a_thrshdVoF,
                                               const bool           a_strict,
                                               const int            a_refinement)
  : ScanShop(a_localGeom, a_verbosity, a_dx, a_probLo, a_finestDomain, a_scanLevel, a_ebGhost, a_thrshdVoF)
{
  CH_assert(a_dx > 0.0);
  CH_assert(a_thrshdVoF >= 0.0 && a_thrshdVoF <= 1.0);

  m_strict          = a_strict;
  m_volumeThreshold = a_thrshdVoF;
  m_refinement      = std::max(1, a_refinement);
  m_coverageLevel   = -1;
  m_coverageBuffer  = 0;
}

PolyhedralGeometryShop::~PolyhedralGeometryShop()
{}

Real
PolyhedralGeometryShop::edgeRoot(const IntVect&  a_edgeIV,
                                 const int       a_dir,
                                 const Real      a_loValue,
                                 const RealVect& a_probLo,
                                 const Real&     a_dx) const noexcept
{
  RealVect lowPoint = a_probLo;

  for (int d = 0; d < SpaceDim; d++) {
    lowPoint[d] += a_dx * static_cast<Real>(a_edgeIV[d]);
  }

  Real lo      = 0.0;
  Real hi      = 1.0;
  Real loValue = a_loValue;

  for (int iter = 0; iter < 100; iter++) {
    const Real mid = 0.5 * (lo + hi);

    RealVect x = lowPoint;
    x[a_dir] += a_dx * mid;

    const Real value = m_baseIF->value(x);

    if (PolyhedralEB::isFluid(value) == PolyhedralEB::isFluid(loValue)) {
      lo      = mid;
      loValue = value;
    }
    else {
      hi = mid;
    }

    if (hi - lo < 1.0E-15) {
      break;
    }
  }

  return 0.5 * (lo + hi);
}

void
PolyhedralGeometryShop::interfaceFacets(Vector<Real>&   a_facets,
                                        const IntVect&  a_cell,
                                        const RealVect& a_probLo,
                                        const Real      a_dx) const noexcept
{
  PolyhedralEB::CutCellSurface surface;

  for (int c = 0; c < PolyhedralEB::CutCellSurface::s_numCorners; c++) {
    RealVect x = a_probLo;

    for (int d = 0; d < SpaceDim; d++) {
      x[d] += a_dx * static_cast<Real>(a_cell[d] + ((c >> d) & 1));
    }

    surface.m_corner[c] = m_baseIF->value(x);
  }

  for (int e = 0; e < PolyhedralEB::CutCellSurface::s_numEdges; e++) {
    int low  = -1;
    int high = -1;

    PolyhedralEB::detail::edgeCorners(e, low, high);

    const Real loValue = surface.m_corner[low];
    const Real hiValue = surface.m_corner[high];

    if (PolyhedralEB::isFluid(loValue) != PolyhedralEB::isFluid(hiValue)) {
      if (loValue == 0.0) {
        surface.m_crossing[e] = 0.0;
      }
      else if (hiValue == 0.0) {
        surface.m_crossing[e] = 1.0;
      }
      else {
        const int dir = PolyhedralEB::detail::edgeDirection(e);

        int offset[SpaceDim];
        PolyhedralEB::detail::edgeOrigin(e, offset);

        IntVect edgeIV = a_cell;

        for (int d = 0; d < SpaceDim; d++) {
          edgeIV[d] += offset[d];
        }

        surface.m_crossing[e] = this->edgeRoot(edgeIV, dir, loValue, a_probLo, a_dx);
      }
    }
  }

  PolyhedralEB::CutCellBody body;

  if (body.define(surface)) {
    RealVect centre = a_probLo;

    for (int d = 0; d < SpaceDim; d++) {
      centre[d] += a_dx * (static_cast<Real>(a_cell[d]) + 0.5);
    }

    body.appendInterfaceFacets(a_facets, centre, a_dx);
  }
}

int
PolyhedralGeometryShop::levelFromDx(const Real a_dx) const noexcept
{
  for (int lvl = 0; lvl < static_cast<int>(m_dx.size()); lvl++) {
    if (std::abs(m_dx[lvl] - a_dx) <= 1.0E-12 * m_dx[lvl]) {
      return lvl;
    }
  }

  return -1;
}

void
PolyhedralGeometryShop::postMakeBoxLayout(const DisjointBoxLayout& a_dbl, const RealVect& a_dx)
{
  CH_TIME("PolyhedralGeometryShop::postMakeBoxLayout");

  ScanShop::postMakeBoxLayout(a_dbl, a_dx);

  const int level = this->levelFromDx(a_dx[0]);

  if (level < 0) {
    return;
  }

  if (static_cast<int>(m_surfaces.size()) < static_cast<int>(m_dx.size())) {
    m_surfaceCells.resize(m_dx.size());
    m_surfaces.resize(m_dx.size());
  }

  // The surfaces are recorded as the graph is filled, which happens once per box of this layout
  // and after this call, so the store is only sized here
  m_surfaceCells[level] = RefCountedPtr<LayoutData<Vector<IntVect>>>(new LayoutData<Vector<IntVect>>(a_dbl));
  m_surfaces[level]     = RefCountedPtr<LayoutData<Vector<PolyhedralEB::CutCellSurface>>>(
    new LayoutData<Vector<PolyhedralEB::CutCellSurface>>(a_dbl));
}

void
PolyhedralGeometryShop::fillNodeValues(BaseFab<Real>&  a_nodeValues,
                                       const Box&      a_region,
                                       const RealVect& a_probLo,
                                       const Real&     a_dx) const
{
  CH_assert(a_dx > 0.0);

  Box nodeBox = a_region;
  nodeBox.surroundingNodes();

  a_nodeValues.define(nodeBox, 1);

  for (BoxIterator bit(nodeBox); bit.ok(); ++bit) {
    const IntVect iv = bit();

    RealVect x = a_probLo;

    for (int d = 0; d < SpaceDim; d++) {
      x[d] += a_dx * static_cast<Real>(iv[d]);
    }

    a_nodeValues(iv, 0) = m_baseIF->value(x);
  }
}

Real
PolyhedralGeometryShop::edgeCrossing(BaseFab<Real>   a_intercept[SpaceDim],
                                     const IntVect&  a_cell,
                                     const int       a_edge,
                                     const Real      a_lo,
                                     const RealVect& a_probLo,
                                     const Real&     a_dx) const
{
  CH_assert(a_dx > 0.0);

  const int dir = PolyhedralEB::detail::edgeDirection(a_edge);

  int offset[SpaceDim];
  PolyhedralEB::detail::edgeOrigin(a_edge, offset);

  // the edge is addressed by the node at its low end
  IntVect edgeIV = a_cell;

  for (int d = 0; d < SpaceDim; d++) {
    edgeIV[d] += offset[d];
  }

  if (a_intercept[dir].box().contains(edgeIV) &&
      a_intercept[dir](edgeIV, 0) != PolyhedralEB::CutCellSurface::s_noCrossing) {
    return a_intercept[dir](edgeIV, 0);
  }

  // bisection between the two endpoints. The endpoints are taken from the shared node values
  // and the edge is addressed by its own index, so every cell reaching this edge hands the
  // solver the same interval and gets the same root back
  const Real root = this->edgeRoot(edgeIV, dir, a_lo, a_probLo, a_dx);

  if (a_intercept[dir].box().contains(edgeIV)) {
    a_intercept[dir](edgeIV, 0) = root;
  }

  return root;
}

void
PolyhedralGeometryShop::fillCorners(PolyhedralEB::CutCellSurface& a_surface,
                                    const BaseFab<Real>&          a_nodeValues,
                                    const IntVect&                a_cell) const
{
  for (int c = 0; c < PolyhedralEB::CutCellSurface::s_numCorners; c++) {
    IntVect node = a_cell;

    for (int d = 0; d < SpaceDim; d++) {
      node[d] += (c >> d) & 1;
    }

    a_surface.m_corner[c] = a_nodeValues(node, 0);
  }
}

void
PolyhedralGeometryShop::buildSurface(BaseFab<Real>                 a_intercept[SpaceDim],
                                     PolyhedralEB::CutCellSurface& a_surface,
                                     const BaseFab<Real>&          a_nodeValues,
                                     const IntVect&                a_cell,
                                     const RealVect&               a_probLo,
                                     const Real&                   a_dx) const
{
  a_surface = PolyhedralEB::CutCellSurface();

  this->fillCorners(a_surface, a_nodeValues, a_cell);

  for (int e = 0; e < PolyhedralEB::CutCellSurface::s_numEdges; e++) {
    int low  = -1;
    int high = -1;

    PolyhedralEB::detail::edgeCorners(e, low, high);

    CH_assert(low >= 0 && low < PolyhedralEB::CutCellSurface::s_numCorners);
    CH_assert(high >= 0 && high < PolyhedralEB::CutCellSurface::s_numCorners);

    const Real loValue = a_surface.m_corner[low];
    const Real hiValue = a_surface.m_corner[high];

    // an edge carries a crossing exactly when its two ends disagree under the one predicate the
    // corners are classified by, so the number of crossings on a face counts sign changes
    if (PolyhedralEB::isFluid(loValue) != PolyhedralEB::isFluid(hiValue)) {
      // a corner at exactly zero is on the interface, so the crossing is that corner rather
      // than a root to be searched for
      if (loValue == 0.0) {
        a_surface.m_crossing[e] = 0.0;
      }
      else if (hiValue == 0.0) {
        a_surface.m_crossing[e] = 1.0;
      }
      else {
        a_surface.m_crossing[e] = this->edgeCrossing(a_intercept, a_cell, e, loValue, a_probLo, a_dx);
      }
    }
  }
}

void
PolyhedralGeometryShop::classifyFromParents(BaseFab<int>&   a_regIrregCovered,
                                            IntVectSet&     a_irregularCells,
                                            const Box&      a_validRegion,
                                            const Box&      a_ghostRegion,
                                            const RealVect& a_probLo,
                                            const Real&     a_dx) const
{
  const Real coarseDx = a_dx * static_cast<Real>(m_refinement);

  Box coarseRegion = a_ghostRegion;
  coarseRegion.coarsen(m_refinement);

  BaseFab<Real> nodeValues;
  this->fillNodeValues(nodeValues, coarseRegion, a_probLo, coarseDx);

  for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit) {
    const IntVect fine = bit();

    IntVect coarse = fine;
    coarse.coarsen(m_refinement);

    PolyhedralEB::CutCellSurface surface;

    for (int c = 0; c < PolyhedralEB::CutCellSurface::s_numCorners; c++) {
      IntVect node = coarse;

      for (int d = 0; d < SpaceDim; d++) {
        node[d] += (c >> d) & 1;
      }

      surface.m_corner[c] = nodeValues(node, 0);
    }

    // a coarse cell the interface never enters gives the same answer to every cell below it,
    // and no body has to be built to find that out
    const PolyhedralEB::CutCellBody::Kind coarseKind = PolyhedralEB::CutCellBody::classify(surface);

    if (coarseKind != PolyhedralEB::CutCellBody::Kind::Cut) {
      a_regIrregCovered(fine, 0) = (coarseKind == PolyhedralEB::CutCellBody::Kind::Regular) ? 1 : -1;

      continue;
    }

    PolyhedralEB::CutCellBody body;

    if (!this->buildRefinedBody(body, fine, a_probLo, a_dx)) {
      if (m_strict) {
        std::ostringstream message;

        message << "PolyhedralGeometryShop::classifyFromParents - could not cut the body down to cell " << fine;

        MayDay::Error(message.str().c_str());
      }

      a_regIrregCovered(fine, 0) = -1;

      continue;
    }

    switch (body.kind()) {
    case PolyhedralEB::CutCellBody::Kind::Covered: {
      a_regIrregCovered(fine, 0) = -1;

      break;
    }
    case PolyhedralEB::CutCellBody::Kind::Regular: {
      a_regIrregCovered(fine, 0) = 1;

      break;
    }
    default: {
      a_regIrregCovered(fine, 0) = 0;

      if (a_validRegion.contains(fine)) {
        a_irregularCells |= fine;
      }

      break;
    }
    }
  }
}

bool
PolyhedralGeometryShop::buildRefinedBody(PolyhedralEB::CutCellBody& a_body,
                                         const IntVect&             a_cell,
                                         const RealVect&            a_probLo,
                                         const Real&                a_dx) const
{
  // The cell's ancestor at the resolution the surface is reconstructed on, and where the cell
  // sits inside it.
  IntVect ancestor = a_cell;
  IntVect offset   = IntVect::Zero;

  int stride = 1;

  for (int r = 1; r < m_refinement; r *= 2) {
    for (int d = 0; d < SpaceDim; d++) {
      offset[d] += (((ancestor[d] % 2) + 2) % 2) * stride;
      ancestor[d] = (ancestor[d] >= 0) ? (ancestor[d] / 2) : -((-ancestor[d] + 1) / 2);
    }

    stride *= 2;
  }

  const Real coarseDx = a_dx * static_cast<Real>(m_refinement);

  // the ancestor's surface is reconstructed at its own resolution, with its own edge cache, so
  // that nothing here depends on the fine level's caches
  BaseFab<Real> intercept[SpaceDim];

  for (int d = 0; d < SpaceDim; d++) {
    Box edgeBox(ancestor, ancestor);

    edgeBox.surroundingNodes();
    edgeBox.enclosedCells(d);

    intercept[d].define(edgeBox, 1);
    intercept[d].setVal(PolyhedralEB::CutCellSurface::s_noCrossing);
  }

  BaseFab<Real> nodeValues;
  this->fillNodeValues(nodeValues, Box(ancestor, ancestor), a_probLo, coarseDx);

  PolyhedralEB::CutCellSurface surface;
  this->buildSurface(intercept, surface, nodeValues, ancestor, a_probLo, coarseDx);

  PolyhedralEB::CutCellBody coarse;

  if (!coarse.define(surface)) {
    return false;
  }

  // cut down to the cell, one level at a time
  for (int r = m_refinement; r > 1; r /= 2) {
    PolyhedralEB::CutCellBody children[1 << SpaceDim];

    if (!coarse.refine(children)) {
      return false;
    }

    int which = 0;

    for (int d = 0; d < SpaceDim; d++) {
      which |= (((offset[d] / (r / 2)) & 1) << d);
    }

    for (int d = 0; d < SpaceDim; d++) {
      offset[d] %= (r / 2);
    }

    coarse = children[which];
  }

  a_body = coarse;

  return true;
}

void
PolyhedralGeometryShop::fillNode(IrregNode&                       a_node,
                                 const PolyhedralEB::CutCellBody& a_body,
                                 const BaseFab<int>&              a_regIrregCovered,
                                 const IntVect&                   a_cell,
                                 const ProblemDomain&             a_domain) const
{
  CH_assert(a_regIrregCovered.box().contains(a_cell));
  CH_assert(a_regIrregCovered(a_cell, 0) == 0);

  a_node.m_cell          = a_cell;
  a_node.m_cellIndex     = 0;
  a_node.m_volFrac       = a_body.volumeFraction();
  a_node.m_volCentroid   = a_body.volumeCentroid();
  a_node.m_bndryCentroid = a_body.boundaryCentroid();

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (SideIterator sit; sit.ok(); ++sit) {
      const IntVect shifted   = a_cell + sign(sit()) * BASISV(dir);
      const int     nodeIndex = a_node.index(dir, sit());

      Vector<int>      arc;
      Vector<Real>     areaFrac;
      Vector<RealVect> faceCentroid;

      // a neighbour inside the domain is inside the ghost region the flags were built on
      CH_assert(!a_domain.contains(shifted) || a_regIrregCovered.box().contains(shifted));

      // the face polygons are oriented outward before their areas are summed, so a net
      // aperture is never negative
      CH_assert(a_body.areaFraction(dir, sit()) >= 0.0);

      // the arcs are topology: they follow the covered set and the domain, not the moments, so
      // the graph is the one GeometryShop would have built
      if (!a_domain.contains(shifted)) {
        arc.resize(1, -1);
      }
      else if (a_regIrregCovered(shifted, 0) < 0) {
        arc.resize(0);
      }
      else {
        arc.resize(1, 0);
      }

      // A face with no area open to flux carries no arc either. The aperture is forced to an
      // exact zero rather than a small number, and two cut cells sharing a face compute it from
      // the same crossings, so they cannot disagree about whether the face is there. The one
      // exception is a regular neighbour, which is full and whose faces are open by definition:
      // withholding the arc on that side leaves the two cells disagreeing about the face and
      // EBGraph stops.
      const bool neighbourIsRegular = a_domain.contains(shifted) && a_regIrregCovered(shifted, 0) > 0;
      const bool faceIsOpen         = neighbourIsRegular || a_body.areaFraction(dir, sit()) > 0.0;

      if (arc.size() > 0 && faceIsOpen) {
        areaFrac.resize(1, a_body.areaFraction(dir, sit()));
        faceCentroid.resize(1, a_body.faceCentroid(dir, sit()));
      }
      else {
        arc.resize(0);
      }

      a_node.m_arc[nodeIndex]          = arc;
      a_node.m_areaFrac[nodeIndex]     = areaFrac;
      a_node.m_faceCentroid[nodeIndex] = faceCentroid;
    }
  }
}

bool
PolyhedralGeometryShop::retainBox(const Box& a_box, const int a_level) const noexcept
{
  if (m_coverageLevel < 0) {
    return true;
  }

  const int which = m_coverageLevel - a_level;

  if (which < 0 || which >= static_cast<int>(m_coverageRegions.size())) {
    return true;
  }

  // Reach past what the region asks for. Coarsening a cell reads what lies under the cells around
  // it, so the outermost cells carried are the ones the cells inside them read and do not come up
  // through coarsening themselves. One cell would do; the ghost region the generator already
  // wants is wider than that, so nothing extra is carried for it.
  const Box parent = grow(coarsen(a_box, 2), m_ebGhost + m_coverageBuffer);

  for (int i = 0; i < m_coverageRegions[which].size(); i++) {
    if (parent.intersectsNotEmpty(coarsen(m_coverageRegions[which][i], 2))) {
      return true;
    }
  }

  return false;
}

void
PolyhedralGeometryShop::setCoverage(const Vector<Vector<Box>>& a_regions,
                                    const ProblemDomain&       a_coarsestDomain,
                                    const int                  a_buffer) noexcept
{
  CH_TIME("PolyhedralGeometryShop::setCoverage");

  CH_assert(a_buffer >= 0);

  m_coverageRegions = a_regions;
  m_coverageBuffer  = a_buffer;
  m_coverageLevel   = -1;

  for (int lvl = 0; lvl < static_cast<int>(m_domains.size()); lvl++) {
    if (m_domains[lvl].domainBox() == a_coarsestDomain.domainBox()) {
      m_coverageLevel = lvl;

      break;
    }
  }

  if (m_coverageLevel < 0 && a_regions.size() > 0) {
    MayDay::Error("PolyhedralGeometryShop::setCoverage - the boxes do not sit on a level of this hierarchy");
  }
}

void
PolyhedralGeometryShop::fillGraph(BaseFab<int>&        a_regIrregCovered,
                                  Vector<IrregNode>&   a_nodes,
                                  const Box&           a_validRegion,
                                  const Box&           a_ghostRegion,
                                  const ProblemDomain& a_domain,
                                  const RealVect&      a_probLo,
                                  const Real&          a_dx,
                                  const DataIndex&     a_di) const
{
  CH_TIME("PolyhedralGeometryShop::fillGraph");

  const int level = (static_cast<int>(m_surfaces.size()) > 0) ? this->levelFromDx(a_dx) : -1;

  CH_assert(a_domain.contains(a_ghostRegion));
  CH_assert(a_ghostRegion.contains(a_validRegion));
  CH_assert(a_dx > 0.0);

  a_regIrregCovered.resize(a_ghostRegion, 1);
  a_nodes.resize(0);

  // Every cell is written below, from the parent that holds it. Starting from a value no cell can
  // end on is what turns a parent this rank turned out not to hold into a refusal rather than
  // into a hole in the graph that nothing downstream would question.
  a_regIrregCovered.setVal(s_unclassified);

  BaseFab<Real> nodeValues;
  this->fillNodeValues(nodeValues, a_ghostRegion, a_probLo, a_dx);

  IntVectSet irregularCells;

  if (m_refinement > 1) {
    this->classifyFromParents(a_regIrregCovered, irregularCells, a_validRegion, a_ghostRegion, a_probLo, a_dx);
  }

  for (BoxIterator bit(a_ghostRegion); bit.ok() && m_refinement == 1; ++bit) {
    const IntVect iv = bit();

    PolyhedralEB::CutCellSurface surface;
    this->fillCorners(surface, nodeValues, iv);

    if (PolyhedralGeometryShop::isRecorded(surface)) {
      a_regIrregCovered(iv, 0) = 0;

      if (a_validRegion.contains(iv)) {
        irregularCells |= iv;
      }
    }
    else {
      a_regIrregCovered(iv,
                        0) = (PolyhedralEB::CutCellBody::classify(surface) == PolyhedralEB::CutCellBody::Kind::Regular)
                               ? 1
                               : -1;
    }
  }

  // a regular cell bordering a covered one is a full cell with a covered face, and the node for
  // it is the one GeometryShop builds
  for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit) {
    if (a_regIrregCovered(bit(), 0) == -1) {
      GeometryShop::fixRegularCellsNextToCovered(a_nodes, a_regIrregCovered, a_validRegion, a_domain, bit(), a_dx);
    }
  }

  BaseFab<Real> intercept[SpaceDim];

  for (int dir = 0; dir < SpaceDim; dir++) {
    Box edgeBox = a_validRegion;
    edgeBox.grow(1);
    edgeBox &= a_ghostRegion;
    edgeBox.surroundingNodes();
    edgeBox.enclosedCells(dir);

    intercept[dir].define(edgeBox, 1);
    intercept[dir].setVal(PolyhedralEB::CutCellSurface::s_noCrossing);
  }

  IntVectSet droppedCells;

  for (IVSIterator ivsIt(irregularCells); ivsIt.ok(); ++ivsIt) {
    const IntVect iv = ivsIt();

    PolyhedralEB::CutCellSurface surface;
    PolyhedralEB::CutCellBody    body;

    bool built = false;

    if (m_refinement > 1) {
      built = this->buildRefinedBody(body, iv, a_probLo, a_dx);
    }
    else {
      this->buildSurface(intercept, surface, nodeValues, iv, a_probLo, a_dx);

      built = body.define(surface);
    }

    // A cell holding a feature finer than itself has no single plane to stand for its interface.
    // Only worth refusing where the cells are claimed to be single valued, which is the finest
    // level: every coarser one is going to fail to resolve something, and is not asked to.
    const bool finestLevel = std::abs(a_dx - m_dx[0]) <= 1.0E-12 * a_dx;

    if (built && finestLevel && !body.interfaceIsOneSided()) {
      built = false;
    }

    if (!built) {
      if (m_strict) {
        std::ostringstream message;

        // The two ways a cell is refused read very differently, and saying which is which is the
        // difference between looking for a bug and reaching for more resolution.
        if (!body.interfaceIsOneSided()) {
          message << "PolyhedralGeometryShop::fillGraph - the interface folds back on itself in cell " << iv
                  << ", so no single plane stands for it (interface area " << body.trueBoundaryArea() << " against "
                  << body.boundaryArea()
                  << " carrying flux). The cell holds a feature finer than itself and wants resolving.";
        }
        else {
          // A body built here is asked to close; one carried up from below is asked to satisfy the
          // divergence identity instead, so both are reported rather than guessing which applies.
          message << "PolyhedralGeometryShop::fillGraph - could not build the body in cell " << iv
                  << " (closure residual " << body.closureResidual() << ", divergence residual "
                  << body.divergenceResidual() << ", volume fraction " << body.volumeFraction() << ")";
        }

        MayDay::Error(message.str().c_str());
      }

      continue;
    }

    if (m_volumeThreshold > 0.0 && body.volumeFraction() < m_volumeThreshold) {
      droppedCells |= iv;

      a_regIrregCovered(iv, 0) = -1;

      continue;
    }

    if (m_strict && body.divergenceResidual() > s_divergenceTolerance) {
      std::ostringstream message;

      message << "PolyhedralGeometryShop::fillGraph - moments in cell " << iv
              << " do not satisfy sum(alpha_hi - alpha_lo) = a_B*n (residual " << body.divergenceResidual() << ")";

      MayDay::Error(message.str().c_str());
    }

    IrregNode node;
    this->fillNode(node, body, a_regIrregCovered, iv, a_domain);

    a_nodes.push_back(node);

    // Keep what the cell was reconstructed from, so that a finer cell can be had by cutting this
    // one rather than by finding its roots again. Only on the path where this cell's own surface
    // is what it was built from: where the body came from cutting a coarser one, the surface that
    // matters is that ancestor's and is already kept against it.
    if (m_refinement == 1 && level >= 0) {
      (*m_surfaceCells[level])[a_di].push_back(iv);
      (*m_surfaces[level])[a_di].push_back(surface);
    }
  }

  this->dropCells(a_nodes, a_regIrregCovered, droppedCells, a_validRegion, a_domain, a_dx);
}

void
PolyhedralGeometryShop::dropCells(Vector<IrregNode>&   a_nodes,
                                  BaseFab<int>&        a_regIrregCovered,
                                  const IntVectSet&    a_droppedCells,
                                  const Box&           a_validRegion,
                                  const ProblemDomain& a_domain,
                                  const Real&          a_dx) const
{
  CH_TIME("PolyhedralGeometryShop::dropCells");

  if (a_droppedCells.isEmpty()) {
    return;
  }

  // Where each cell's node sits in a_nodes. Rebuilt from the nodes rather than carried along by
  // the caller, and appended to below as the regular cells this converts get nodes of their own,
  // which a later dropped cell may need to find.
  BaseFab<int> nodeIndex(a_validRegion, 1);
  nodeIndex.setVal(-1);

  auto indexNodes = [&](const int a_from) -> void {
    for (int n = a_from; n < static_cast<int>(a_nodes.size()); n++) {
      const IntVect& cell = a_nodes[n].m_cell;

      CH_assert(a_validRegion.contains(cell));
      CH_assert(nodeIndex(cell, 0) < 0);

      nodeIndex(cell, 0) = n;
    }
  };

  indexNodes(0);

  for (IVSIterator ivsIt(a_droppedCells); ivsIt.ok(); ++ivsIt) {
    const IntVect iv = ivsIt();

    for (int dir = 0; dir < SpaceDim; dir++) {
      for (SideIterator sit; sit.ok(); ++sit) {
        const IntVect other = iv + sign(sit()) * BASISV(dir);

        if (!a_validRegion.contains(other) || a_regIrregCovered(other, 0) != 0) {
          continue;
        }

        const int n = nodeIndex(other, 0);

        if (n < 0) {
          MayDay::Error("PolyhedralGeometryShop::dropCells - an irregular neighbour has no node");
        }

        const int arcIndex = a_nodes[n].index(dir, flip(sit()));

        a_nodes[n].m_arc[arcIndex].resize(0);
        a_nodes[n].m_areaFrac[arcIndex].resize(0);
        a_nodes[n].m_faceCentroid[arcIndex].resize(0);
      }
    }

    const int numNodes = static_cast<int>(a_nodes.size());

    GeometryShop::fixRegularCellsNextToCovered(a_nodes, a_regIrregCovered, a_validRegion, a_domain, iv, a_dx);

    indexNodes(numNodes);
  }
}

bool
PolyhedralGeometryShop::isRecorded(const PolyhedralEB::CutCellSurface& a_surface) noexcept
{
  const PolyhedralEB::CutCellBody::Kind kind = PolyhedralEB::CutCellBody::classify(a_surface);

  if (kind == PolyhedralEB::CutCellBody::Kind::Covered) {
    return false;
  }

  if (kind == PolyhedralEB::CutCellBody::Kind::Regular) {
    return PolyhedralEB::CutCellBody::interfaceLiesInFace(a_surface);
  }

  return true;
}

int
PolyhedralGeometryShop::cuttableParent(const Box& a_ghostRegion, const int a_level) const noexcept
{
  for (int lvl = a_level + 1; lvl < static_cast<int>(m_surfaces.size()); lvl++) {
    if (m_surfaces[lvl].isNull()) {
      continue;
    }

    const Box coarseRegion = coarsen(a_ghostRegion, 1 << (lvl - a_level));

    // The boxes of a layout are disjoint, so what the level describes is the sum of the overlaps
    // and it describes the whole region exactly when that comes to all of it. Every box of the
    // level is walked, not only this rank's, so that every rank reaches the same answer.
    const BoxLayout& dbl = m_surfaces[lvl]->boxLayout();

    long long covered = 0;

    for (LayoutIterator lit = dbl.layoutIterator(); lit.ok(); ++lit) {
      const Box overlap = dbl[lit()] & coarseRegion;

      covered += overlap.numPts();
    }

    if (covered == coarseRegion.numPts()) {
      return lvl;
    }
  }

  return -1;
}

int
PolyhedralGeometryShop::numSurfaceComponents() const
{
  return (m_refinement == 1) ? PolyhedralEB::CutCellSurface::s_numValues : 0;
}

void
PolyhedralGeometryShop::getSurfaces(Vector<IntVect>& a_cells,
                                    Vector<Real>&    a_values,
                                    const Box&       a_region,
                                    const Real&      a_dx) const
{
  CH_TIME("PolyhedralGeometryShop::getSurfaces");

  a_cells.resize(0);
  a_values.resize(0);

  const int level = this->levelFromDx(a_dx);

  if (level < 0 || level >= static_cast<int>(m_surfaces.size()) || m_surfaces[level].isNull()) {
    return;
  }

  // Every box this object holds of the level is looked through, rather than the one the caller
  // names: the level that asks may have grown since it was made, and its boxes are then not the
  // ones the surfaces were recorded against.
  const BoxLayout& dbl = m_surfaces[level]->boxLayout();

  for (DataIterator dit = m_surfaces[level]->dataIterator(); dit.ok(); ++dit) {
    if (!dbl[dit()].intersectsNotEmpty(a_region)) {
      continue;
    }

    const Vector<IntVect>&                      cells    = (*m_surfaceCells[level])[dit()];
    const Vector<PolyhedralEB::CutCellSurface>& surfaces = (*m_surfaces[level])[dit()];

    for (int n = 0; n < cells.size(); n++) {
      if (!a_region.contains(cells[n])) {
        continue;
      }

      a_cells.push_back(cells[n]);

      const int offset = a_values.size();
      a_values.resize(offset + PolyhedralEB::CutCellSurface::s_numValues);

      surfaces[n].store(&a_values[offset], 1);
    }
  }
}

Real
PolyhedralGeometryShop::refinedFillParentDx(const Box& a_ghostRegion, const Real& a_dx) const
{
  if (m_refinement != 1) {
    return -1.0;
  }

  const int level = this->levelFromDx(a_dx);

  if (level < 0) {
    return -1.0;
  }

  const int parent = this->cuttableParent(a_ghostRegion, level);

  return (parent < 0) ? -1.0 : m_dx[parent];
}

bool
PolyhedralGeometryShop::fillRefinedGraph(BaseFab<int>&          a_regIrregCovered,
                                         Vector<IrregNode>&     a_nodes,
                                         const Box&             a_validRegion,
                                         const Box&             a_ghostRegion,
                                         const ProblemDomain&   a_domain,
                                         const RealVect&        a_probLo,
                                         const Real&            a_dx,
                                         const BaseIVFAB<Real>& a_parents,
                                         const Real&            a_parentDx) const
{
  CH_TIME("PolyhedralGeometryShop::fillRefinedGraph");

  CH_assert(a_domain.contains(a_ghostRegion));
  CH_assert(a_ghostRegion.contains(a_validRegion));
  CH_assert(a_dx > 0.0);

  // Only the path that reconstructs each cell's surface on the cell itself records one, so it is
  // the only path a finer cell can be cut from.
  if (m_refinement != 1) {
    return false;
  }

  const int level = this->levelFromDx(a_dx);

  if (level < 0) {
    return false;
  }

  // Which level the parents came from is the caller's word, since the caller is what fetched
  // them. It has to be one of this hierarchy's and coarser than the level being filled.
  const int parent = this->levelFromDx(a_parentDx);

  if (parent <= level) {
    return false;
  }

  const int  depth      = parent - level;
  const int  refinement = 1 << depth;
  const Real coarseDx   = a_dx * static_cast<Real>(refinement);

  const Box coarseRegion = coarsen(a_ghostRegion, refinement);

  // A parent the interface never entered has no surface stored against it, and its corner values
  // are enough to say which side of the interface everything under it lies on. Wanted by the
  // locality test as well as by the fill, so taken once.
  BaseFab<Real> nodeValues;
  this->fillNodeValues(nodeValues, coarseRegion, a_probLo, coarseDx);

  a_regIrregCovered.resize(a_ghostRegion, 1);
  a_nodes.resize(0);

  // The bodies of the cut cells this box owns, kept until the classification of the whole ghost
  // region is known, since a node reads the classification of its neighbours. Keeping them costs
  // the surface of the box; cutting the parents a second time to avoid it would cost its volume.
  Vector<PolyhedralEB::CutCellBody> bodies;

  BaseFab<int> bodyIndex(a_validRegion, 1);
  bodyIndex.setVal(-1);

  // The parent whose refinement is being walked, so that a child that comes out wrong says which
  // cut produced it rather than only where it landed.
  IntVect parentCell = IntVect::Zero;

  auto visit = [&](const IntVect& a_cell, const PolyhedralEB::CutCellBody& a_body, const int a_depthRemaining) -> void {
    Box sub(a_cell, a_cell);
    sub.refine(1 << a_depthRemaining);
    sub &= a_ghostRegion;

    if (sub.isEmpty()) {
      return;
    }

    switch (a_body.kind()) {
    case PolyhedralEB::CutCellBody::Kind::Covered: {
      a_regIrregCovered.setVal(-1, sub, 0, 1);

      break;
    }
    case PolyhedralEB::CutCellBody::Kind::Regular: {
      a_regIrregCovered.setVal(1, sub, 0, 1);

      break;
    }
    default: {
      // only a cell followed to the bottom is still cut, so this is one cell
      CH_assert(a_depthRemaining == 0);

      a_regIrregCovered(a_cell, 0) = 0;

      if (a_validRegion.contains(a_cell)) {
        // A body cut from its parent is asked to satisfy the divergence identity rather than to
        // close: it is a piece of a closed body, and the identity is what the coarsening of it
        // will be read through. Asked here rather than where the node is written, so that the
        // parent it came from is still known.
        if (m_strict && a_body.divergenceResidual() > s_divergenceTolerance) {
          std::ostringstream message;

          message << "PolyhedralGeometryShop::fillRefinedGraph - cell " << a_cell << ", cut " << depth
                  << " levels out of cell " << parentCell << " on level " << parent
                  << ", does not satisfy sum(alpha_hi - alpha_lo) = a_B*n (residual " << a_body.divergenceResidual()
                  << ", volume fraction " << a_body.volumeFraction() << ")";

          MayDay::Error(message.str().c_str());
        }

        bodyIndex(a_cell, 0) = static_cast<int>(bodies.size());

        bodies.push_back(a_body);
      }

      break;
    }
    }
  };

  const IntVectSet& held = a_parents.getIVS();

  for (BoxIterator bit(coarseRegion); bit.ok(); ++bit) {
    const IntVect iv = bit();

    if (!held.contains(iv)) {
      // Nothing was kept for this parent, so it was whole or empty and everything under it is the
      // same. A parent whose corners say otherwise was declined when the level was generated, or
      // the caller handed over an incomplete set; either way there is nothing to cut it from.
      PolyhedralEB::CutCellSurface surface;
      this->fillCorners(surface, nodeValues, iv);

      if (PolyhedralGeometryShop::isRecorded(surface)) {
        if (m_strict) {
          std::ostringstream message;

          message << "PolyhedralGeometryShop::fillRefinedGraph - cell " << iv << " on level " << parent
                  << " is cut but no surface was handed over for it";

          MayDay::Error(message.str().c_str());
        }

        return false;
      }

      Box sub(iv, iv);
      sub.refine(refinement);
      sub &= a_ghostRegion;

      const PolyhedralEB::CutCellBody::Kind kind = PolyhedralEB::CutCellBody::classify(surface);

      a_regIrregCovered.setVal((kind == PolyhedralEB::CutCellBody::Kind::Regular) ? 1 : -1, sub, 0, 1);

      continue;
    }

    // The generator mandates single-valued cut cells, so a parent is one volume of fluid and its
    // surface sits at the first index.
    const VolIndex vof(iv, 0);

    // Read out a component at a time rather than off one pointer: a BaseIVFAB holds a component
    // at a time, so consecutive components of one cell are a whole level's cut cells apart.
    Real values[PolyhedralEB::CutCellSurface::s_numValues];

    for (int c = 0; c < PolyhedralEB::CutCellSurface::s_numValues; c++) {
      values[c] = a_parents(vof, c);
    }

    PolyhedralEB::CutCellSurface surface;
    surface.load(values, 1);

    PolyhedralEB::CutCellBody body;

    if (!body.define(surface)) {
      if (m_strict) {
        std::ostringstream message;

        message << "PolyhedralGeometryShop::fillRefinedGraph - the surface handed over for cell " << iv << " on level "
                << parent << " does not close (residual " << body.closureResidual() << ")";

        MayDay::Error(message.str().c_str());
      }

      return false;
    }

    parentCell = iv;

    if (!PolyhedralEB::refineSubtree(body, iv, depth, visit)) {
      if (m_strict) {
        std::ostringstream message;

        message << "PolyhedralGeometryShop::fillRefinedGraph - could not cut cell " << iv << " on level " << parent
                << " down " << depth << " levels";

        MayDay::Error(message.str().c_str());
      }

      return false;
    }
  }

  for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit) {
    if (a_regIrregCovered(bit(), 0) == s_unclassified) {
      if (m_strict) {
        std::ostringstream message;

        message << "PolyhedralGeometryShop::fillRefinedGraph - no parent on level " << parent << " answered for cell "
                << bit();

        MayDay::Error(message.str().c_str());
      }

      return false;
    }
  }

  // a regular cell bordering a covered one is a full cell with a covered face, and the node for
  // it is the one GeometryShop builds
  for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit) {
    if (a_regIrregCovered(bit(), 0) == -1) {
      GeometryShop::fixRegularCellsNextToCovered(a_nodes, a_regIrregCovered, a_validRegion, a_domain, bit(), a_dx);
    }
  }

  IntVectSet droppedCells;

  for (BoxIterator bit(a_validRegion); bit.ok(); ++bit) {
    const IntVect iv = bit();

    const int n = bodyIndex(iv, 0);

    if (n < 0) {
      continue;
    }

    const PolyhedralEB::CutCellBody& body = bodies[n];

    if (m_volumeThreshold > 0.0 && body.volumeFraction() < m_volumeThreshold) {
      droppedCells |= iv;

      a_regIrregCovered(iv, 0) = -1;

      continue;
    }

    IrregNode node;
    this->fillNode(node, body, a_regIrregCovered, iv, a_domain);

    a_nodes.push_back(node);
  }

  this->dropCells(a_nodes, a_regIrregCovered, droppedCells, a_validRegion, a_domain, a_dx);

  return true;
}

#include <CD_NamespaceFooter.H>
