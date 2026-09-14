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
                                               const bool           a_strict)
  : ScanShop(a_localGeom, a_verbosity, a_dx, a_probLo, a_finestDomain, a_scanLevel, a_ebGhost, a_thrshdVoF)
{
  CH_assert(a_dx > 0.0);
  CH_assert(a_thrshdVoF >= 0.0 && a_thrshdVoF <= 1.0);

  m_strict          = a_strict;
  m_volumeThreshold = a_thrshdVoF;
}

PolyhedralGeometryShop::~PolyhedralGeometryShop()
{}

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
  RealVect lowPoint = a_probLo;

  for (int d = 0; d < SpaceDim; d++) {
    lowPoint[d] += a_dx * static_cast<Real>(edgeIV[d]);
  }

  Real lo      = 0.0;
  Real hi      = 1.0;
  Real loValue = a_lo;

  for (int iter = 0; iter < 100; iter++) {
    const Real mid = 0.5 * (lo + hi);

    RealVect x = lowPoint;
    x[dir] += a_dx * mid;

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

  const Real root = 0.5 * (lo + hi);

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

  CH_assert(a_domain.contains(a_ghostRegion));
  CH_assert(a_ghostRegion.contains(a_validRegion));
  CH_assert(a_dx > 0.0);

  a_regIrregCovered.resize(a_ghostRegion, 1);
  a_nodes.resize(0);

  BaseFab<Real> nodeValues;
  this->fillNodeValues(nodeValues, a_ghostRegion, a_probLo, a_dx);

  IntVectSet irregularCells;

  for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit) {
    const IntVect iv = bit();

    PolyhedralEB::CutCellSurface surface;
    this->fillCorners(surface, nodeValues, iv);

    switch (PolyhedralEB::CutCellBody::classify(surface)) {
    case PolyhedralEB::CutCellBody::Kind::Covered: {
      a_regIrregCovered(iv, 0) = -1;

      break;
    }
    case PolyhedralEB::CutCellBody::Kind::Regular: {
      // A cell the interface only grazes is full, which is why it is classified regular, but the
      // face the interface lies in is closed. A regular cell has every face open, so this one is
      // given a body instead: leaving it regular has the cut cell across that face read the face
      // from the interface, find nothing of it open, and drop its side of a face this cell keeps.
      if (PolyhedralEB::CutCellBody::interfaceLiesInFace(surface)) {
        a_regIrregCovered(iv, 0) = 0;

        if (a_validRegion.contains(iv)) {
          irregularCells |= iv;
        }
      }
      else {
        a_regIrregCovered(iv, 0) = 1;
      }

      break;
    }
    default: {
      a_regIrregCovered(iv, 0) = 0;

      if (a_validRegion.contains(iv)) {
        irregularCells |= iv;
      }

      break;
    }
    }
  }

  // A node is looked up again by its cell when a neighbour of it is dropped, so the position of
  // each cell's node in a_nodes is kept per cell. Nodes are appended in three places, and each
  // registers what it appended.
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

  // a regular cell bordering a covered one is a full cell with a covered face, and the node for
  // it is the one GeometryShop builds
  for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit) {
    if (a_regIrregCovered(bit(), 0) == -1) {
      GeometryShop::fixRegularCellsNextToCovered(a_nodes, a_regIrregCovered, a_validRegion, a_domain, bit(), a_dx);
    }
  }

  indexNodes(0);

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
    this->buildSurface(intercept, surface, nodeValues, iv, a_probLo, a_dx);

    PolyhedralEB::CutCellBody body;

    if (!body.define(surface)) {
      if (m_strict) {
        std::ostringstream message;

        message << "PolyhedralGeometryShop::fillGraph - could not close the body in cell " << iv
                << " (closure residual " << body.closureResidual() << ", volume fraction " << body.volumeFraction()
                << ")";

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

    indexNodes(static_cast<int>(a_nodes.size()) - 1);
  }

  // A cell dropped for carrying no fluid worth keeping has just become covered, so the faces its
  // irregular neighbours point back with no longer lead anywhere, and the regular ones among
  // them are now full cells with a covered face.
  for (IVSIterator ivsIt(droppedCells); ivsIt.ok(); ++ivsIt) {
    const IntVect iv = ivsIt();

    for (int dir = 0; dir < SpaceDim; dir++) {
      for (SideIterator sit; sit.ok(); ++sit) {
        const IntVect other = iv + sign(sit()) * BASISV(dir);

        if (!a_validRegion.contains(other) || a_regIrregCovered(other, 0) != 0) {
          continue;
        }

        const int n = nodeIndex(other, 0);

        if (n < 0) {
          MayDay::Error("PolyhedralGeometryShop::fillGraph - an irregular neighbour has no node");
        }

        const int arcIndex = a_nodes[n].index(dir, flip(sit()));

        a_nodes[n].m_arc[arcIndex].resize(0);
        a_nodes[n].m_areaFrac[arcIndex].resize(0);
        a_nodes[n].m_faceCentroid[arcIndex].resize(0);
      }
    }

    // the regular neighbours this converts get nodes of their own, which a later dropped cell
    // may need to find
    const int numNodes = static_cast<int>(a_nodes.size());

    GeometryShop::fixRegularCellsNextToCovered(a_nodes, a_regIrregCovered, a_validRegion, a_domain, iv, a_dx);

    indexNodes(numNodes);
  }

  (void)a_di;
}

#include <CD_NamespaceFooter.H>
