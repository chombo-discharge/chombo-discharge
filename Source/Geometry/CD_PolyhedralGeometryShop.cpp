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
#include <IntVectSet.H>
#include <LoHiSide.H>
#include <MayDay.H>

// Our includes
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
                                     const Real      a_hi,
                                     const RealVect& a_probLo,
                                     const Real&     a_dx) const
{
  const int dir = a_edge / (CutCellSurface::s_numEdges / SpaceDim);

  IntVect edgeIV = a_cell;

  {
    const int local = a_edge % (CutCellSurface::s_numEdges / SpaceDim);

#if CH_SPACEDIM == 3
    constexpr int transverse[3][2] = {{1, 2}, {0, 2}, {0, 1}};

    edgeIV[transverse[dir][0]] += local & 1;
    edgeIV[transverse[dir][1]] += (local >> 1) & 1;
#else
    edgeIV[1 - dir] += local & 1;
#endif
  }

  if (a_intercept[dir].box().contains(edgeIV) && a_intercept[dir](edgeIV, 0) != CutCellSurface::s_noCrossing) {
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

    if (CutCellBody::isFluid(value) == CutCellBody::isFluid(loValue)) {
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

  (void)a_hi;

  return root;
}

void
PolyhedralGeometryShop::buildSurface(BaseFab<Real>        a_intercept[SpaceDim],
                                     CutCellSurface&      a_surface,
                                     const BaseFab<Real>& a_nodeValues,
                                     const IntVect&       a_cell,
                                     const RealVect&      a_probLo,
                                     const Real&          a_dx) const
{
  a_surface = CutCellSurface();

  for (int c = 0; c < CutCellSurface::s_numCorners; c++) {
    IntVect node = a_cell;

    for (int d = 0; d < SpaceDim; d++) {
      node[d] += (c >> d) & 1;
    }

    a_surface.m_corner[c] = a_nodeValues(node, 0);
  }

  for (int e = 0; e < CutCellSurface::s_numEdges; e++) {
    const int dir   = e / (CutCellSurface::s_numEdges / SpaceDim);
    const int local = e % (CutCellSurface::s_numEdges / SpaceDim);

    int low = 0;

#if CH_SPACEDIM == 3
    constexpr int transverse[3][2] = {{1, 2}, {0, 2}, {0, 1}};

    low |= (local & 1) << transverse[dir][0];
    low |= ((local >> 1) & 1) << transverse[dir][1];
#else
    low |= (local & 1) << (1 - dir);
#endif

    const int  high    = low | (1 << dir);
    const Real loValue = a_surface.m_corner[low];
    const Real hiValue = a_surface.m_corner[high];

    // an edge carries a crossing exactly when its two ends disagree under the one predicate the
    // corners are classified by, so the number of crossings on a face counts sign changes
    if (CutCellBody::isFluid(loValue) != CutCellBody::isFluid(hiValue)) {
      // a corner at exactly zero is on the interface, so the crossing is that corner rather
      // than a root to be searched for
      if (loValue == 0.0) {
        a_surface.m_crossing[e] = 0.0;
      }
      else if (hiValue == 0.0) {
        a_surface.m_crossing[e] = 1.0;
      }
      else {
        a_surface.m_crossing[e] = this->edgeCrossing(a_intercept, a_cell, e, loValue, hiValue, a_probLo, a_dx);
      }
    }
  }
}

void
PolyhedralGeometryShop::fillNode(IrregNode&           a_node,
                                 const CutCellBody&   a_body,
                                 const BaseFab<int>&  a_regIrregCovered,
                                 const IntVect&       a_cell,
                                 const ProblemDomain& a_domain) const
{
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

      // A face with no area open to flux carries no arc either, but only where the cell across
      // it is itself cut. The aperture is forced to an exact zero rather than a small number,
      // and two cut cells sharing a face compute it from the same crossings, so they cannot
      // disagree about whether the face is there. A regular neighbour has no crossings to agree
      // with and its faces are open by definition, so the arc stays whatever the aperture says.
      const bool neighbourIsCut = a_domain.contains(shifted) && a_regIrregCovered(shifted, 0) == 0;
      const bool faceIsOpen     = !neighbourIsCut || a_body.areaFraction(dir, sit()) > 0.0;

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

  a_regIrregCovered.resize(a_ghostRegion, 1);
  a_nodes.resize(0);

  BaseFab<Real> nodeValues;
  this->fillNodeValues(nodeValues, a_ghostRegion, a_probLo, a_dx);

  IntVectSet irregularCells;

  for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit) {
    const IntVect iv = bit();

    CutCellSurface surface;

    for (int c = 0; c < CutCellSurface::s_numCorners; c++) {
      IntVect node = iv;

      for (int d = 0; d < SpaceDim; d++) {
        node[d] += (c >> d) & 1;
      }

      surface.m_corner[c] = nodeValues(node, 0);
    }

    switch (CutCellBody::classify(surface)) {
    case CutCellBody::Kind::Covered: {
      a_regIrregCovered(iv, 0) = -1;

      break;
    }
    case CutCellBody::Kind::Regular: {
      a_regIrregCovered(iv, 0) = 1;

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
    intercept[dir].setVal(CutCellSurface::s_noCrossing);
  }

  IntVectSet droppedCells;

  for (IVSIterator ivsIt(irregularCells); ivsIt.ok(); ++ivsIt) {
    const IntVect iv = ivsIt();

    CutCellSurface surface;
    this->buildSurface(intercept, surface, nodeValues, iv, a_probLo, a_dx);

    CutCellBody body;

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

        bool found = false;

        for (int n = 0; n < a_nodes.size() && !found; n++) {
          if (a_nodes[n].m_cell == other) {
            const int arcIndex = a_nodes[n].index(dir, flip(sit()));

            a_nodes[n].m_arc[arcIndex].resize(0);
            a_nodes[n].m_areaFrac[arcIndex].resize(0);
            a_nodes[n].m_faceCentroid[arcIndex].resize(0);

            found = true;
          }
        }

        if (!found) {
          MayDay::Error("PolyhedralGeometryShop::fillGraph - an irregular neighbour has no node");
        }
      }
    }

    GeometryShop::fixRegularCellsNextToCovered(a_nodes, a_regIrregCovered, a_validRegion, a_domain, iv, a_dx);
  }

  (void)a_di;
}

#include <CD_NamespaceFooter.H>
