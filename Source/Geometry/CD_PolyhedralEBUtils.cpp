/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PolyhedralEBUtils.cpp
 * @brief  Implementation of CD_PolyhedralEBUtils.H
 * @author Robert Marskar
 */

// Std includes
#include <cmath>

// Chombo includes
#include <BoxIterator.H>
#include <CH_Timer.H>
#include <CH_assert.H>

// Our includes
#include <CD_PolyhedralEBUtils.H>
#include <CD_CutCellBody.H>
#include <CD_PolyUtils.H>
#include <CD_NamespaceHeader.H>

namespace PolyhedralEB {

void
fillNodeValues(const BaseIF&   a_function,
               BaseFab<Real>&  a_nodeValues,
               const Box&      a_region,
               const RealVect& a_probLo,
               const Real&     a_dx)
{
  CH_TIME("PolyhedralEB::fillNodeValues");

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

    a_nodeValues(iv, 0) = snappedValue(a_function, x, a_dx);
  }
}

void
defineIntercepts(BaseFab<Real> a_intercept[SpaceDim], const Box& a_region)
{
  CH_TIME("PolyhedralEB::defineIntercepts");

  for (int dir = 0; dir < SpaceDim; dir++) {
    Box edgeBox = a_region;
    edgeBox.surroundingNodes();
    edgeBox.enclosedCells(dir);

    a_intercept[dir].define(edgeBox, 1);
    a_intercept[dir].setVal(PolyhedralEB::CutCellSurface::s_noCrossing);
  }
}

Real
edgeCrossing(const BaseIF&   a_function,
             BaseFab<Real>   a_intercept[SpaceDim],
             const IntVect&  a_cell,
             const int       a_edge,
             const Real      a_lo,
             const Real      a_hi,
             const RealVect& a_probLo,
             const Real&     a_dx)
{
  CH_assert(a_dx > 0.0);
  CH_assert(a_edge >= 0 && a_edge < PolyhedralEB::CutCellSurface::s_numEdges);
  CH_assert(PolyhedralEB::isFluid(a_lo) != PolyhedralEB::isFluid(a_hi));

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
  const Real root = edgeRoot(a_function, edgeIV, dir, a_lo, a_hi, a_probLo, a_dx);

  if (a_intercept[dir].box().contains(edgeIV)) {
    a_intercept[dir](edgeIV, 0) = root;
  }

  return root;
}

bool
isDust(const PolyhedralEB::CutCellBody& a_body) noexcept
{
  // Exactly, not within a tolerance: a degenerate body has no polygon of positive area at all, so both of these
  // are zero by construction rather than by luck.
  return (1.0 - a_body.volumeFraction() <= 0.0) && (a_body.boundaryArea() <= 0.0);
}

Real
snappedValue(const BaseIF& a_function, const RealVect& a_point, const Real a_dx) noexcept
{
  CH_assert(a_dx > 0.0);

  const Real value = a_function.value(a_point);

  if (std::abs(value) > s_snapTolerance * a_dx) {
    return value;
  }

  return resolveTangency(a_function, a_point, a_dx);
}

Real
resolveTangency(const BaseIF& a_function, const RealVect& a_point, const Real a_dx) noexcept
{
  CH_assert(a_dx > 0.0);

  const Real delta    = s_probeSpacing * a_dx;
  const Real zeroBand = s_snapTolerance * a_dx;

  Real probe[2 * SpaceDim];

  for (int dir = 0; dir < SpaceDim; dir++) {
    const RealVect e = BASISREALV(dir);

    probe[2 * dir]     = a_function.value(a_point - delta * e);
    probe[2 * dir + 1] = a_function.value(a_point + delta * e);
  }

  Real reading = 0.0;

  for (int dir = 0; dir < SpaceDim; dir++) {
    const Real lo = probe[2 * dir];
    const Real hi = probe[2 * dir + 1];

    if ((lo < 0.0 && hi > 0.0) || (lo > 0.0 && hi < 0.0)) {
      return 0.0;
    }

    if (std::abs(lo) > zeroBand && std::abs(hi) > zeroBand) {
      if (std::abs(lo) > std::abs(reading)) {
        reading = lo;
      }
      if (std::abs(hi) > std::abs(reading)) {
        reading = hi;
      }
    }
  }

  return reading;
}

Real
edgeRoot(const BaseIF&   a_function,
         const IntVect&  a_edgeIV,
         const int       a_dir,
         const Real      a_loValue,
         const Real      a_hiValue,
         const RealVect& a_probLo,
         const Real      a_dx) noexcept
{
  CH_assert(a_dx > 0.0);
  CH_assert(a_dir >= 0 && a_dir < SpaceDim);

  RealVect lowPoint = a_probLo;

  for (int d = 0; d < SpaceDim; d++) {
    lowPoint[d] += a_dx * static_cast<Real>(a_edgeIV[d]);
  }

  // An edge whose ends are both away from zero and of opposite sign carries one plain root, and Brent's method
  // brackets it in a handful of evaluations where bisection needs tens. An end that reads as exactly zero is a
  // different problem: the function may be zero along a stretch of the edge, and the crossing is then where it
  // leaves that stretch rather than where it reaches it. That is a step in the fluid predicate, which has no
  // continuous root to interpolate, so those edges are walked by bisecting the predicate below.
  if (a_loValue != 0.0 && a_hiValue != 0.0 && ((a_loValue < 0.0) != (a_hiValue < 0.0))) {
    const auto alongEdge = [&](const Real a_t) -> Real {
      RealVect x = lowPoint;

      x[a_dir] += a_dx * a_t;

      return snappedValue(a_function, x, a_dx);
    };

    const Real brent = PolyUtils::brentSolve(0.0, 1.0, alongEdge);

    if (brent >= 0.0 && brent <= 1.0) {
      return brent;
    }
  }

  Real lo      = 0.0;
  Real hi      = 1.0;
  Real loValue = a_loValue;

  for (int iter = 0; iter < 100; iter++) {
    const Real mid = 0.5 * (lo + hi);

    RealVect x = lowPoint;
    x[a_dir] += a_dx * mid;

    const Real value = snappedValue(a_function, x, a_dx);

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

  Real root = 0.5 * (lo + hi);

  // an endpoint exactly on the interface owns a crossing that lands close to it
  if (a_loValue == 0.0 && root < s_rootSnap) {
    root = 0.0;
  }
  else if (a_hiValue == 0.0 && root > 1.0 - s_rootSnap) {
    root = 1.0;
  }

  return root;
}

void
fillCorners(PolyhedralEB::CutCellSurface& a_surface, const BaseFab<Real>& a_nodeValues, const IntVect& a_cell)
{
  CH_assert(a_nodeValues.box().contains(surroundingNodes(Box(a_cell, a_cell))));

  for (int c = 0; c < PolyhedralEB::CutCellSurface::s_numCorners; c++) {
    IntVect node = a_cell;

    for (int d = 0; d < SpaceDim; d++) {
      node[d] += (c >> d) & 1;
    }

    a_surface.m_corner[c] = a_nodeValues(node, 0);
  }
}

void
buildSurface(const BaseIF&                 a_function,
             BaseFab<Real>                 a_intercept[SpaceDim],
             PolyhedralEB::CutCellSurface& a_surface,
             const BaseFab<Real>&          a_nodeValues,
             const IntVect&                a_cell,
             const RealVect&               a_probLo,
             const Real&                   a_dx)
{
  CH_assert(a_dx > 0.0);

  a_surface = PolyhedralEB::CutCellSurface();

  fillCorners(a_surface, a_nodeValues, a_cell);

  for (int e = 0; e < PolyhedralEB::CutCellSurface::s_numEdges; e++) {
    int low  = -1;
    int high = -1;

    PolyhedralEB::detail::edgeCorners(e, low, high);

    CH_assert(low >= 0 && low < PolyhedralEB::CutCellSurface::s_numCorners);
    CH_assert(high >= 0 && high < PolyhedralEB::CutCellSurface::s_numCorners);

    const Real loValue = a_surface.m_corner[low];
    const Real hiValue = a_surface.m_corner[high];

    // an edge carries a crossing exactly when its two ends disagree under the one predicate the
    // corners are classified by, so the number of crossings on a face counts sign changes. A corner
    // at exactly zero is not by itself the crossing: the function may be zero along part of the
    // edge, and edgeRoot finds where it stops being so
    if (PolyhedralEB::isFluid(loValue) != PolyhedralEB::isFluid(hiValue)) {
      a_surface.m_crossing[e] = edgeCrossing(a_function, a_intercept, a_cell, e, loValue, hiValue, a_probLo, a_dx);
    }
  }

  decideFaces(a_function, a_surface, a_cell, a_probLo, a_dx);
}

void
decideFaces(const BaseIF&                 a_function,
            PolyhedralEB::CutCellSurface& a_surface,
            const IntVect&                a_cell,
            const RealVect&               a_probLo,
            const Real                    a_dx) noexcept
{
  CH_assert(a_dx > 0.0);

  using PolyhedralEB::isFluid;

  a_surface.m_fluidJoined = 0;

  // which way one face goes, from its four corners in circuit order and its centre
  auto fluidJoined = [&](const Real      a_f00,
                         const Real      a_f10,
                         const Real      a_f11,
                         const Real      a_f01,
                         const RealVect& a_centre,
                         bool&           a_alternating) -> bool {
    a_alternating = isFluid(a_f00) == isFluid(a_f11) && isFluid(a_f10) == isFluid(a_f01) &&
                    isFluid(a_f00) != isFluid(a_f10);

    if (!a_alternating) {
      return false;
    }

    const Real value = snappedValue(a_function, a_centre, 0.5 * a_dx);

    if (value != 0.0) {
      return isFluid(value);
    }

    const Real den = a_f00 + a_f11 - a_f10 - a_f01;

    return isFluid((std::abs(den) > 0.0) ? (a_f00 * a_f11 - a_f10 * a_f01) / den : (a_f00 + a_f11));
  };

#if CH_SPACEDIM == 3
  for (int dir = 0; dir < SpaceDim; dir++) {
    for (int side = 0; side < 2; side++) {
      int faceCorner[4];

      PolyhedralEB::detail::faceCorners(dir, side, faceCorner);

      RealVect centre = a_probLo;

      for (int d = 0; d < SpaceDim; d++) {
        const Real offset = (d == dir) ? static_cast<Real>(side) : 0.5;

        centre[d] += a_dx * (static_cast<Real>(a_cell[d]) + offset);
      }

      bool alternating = false;

      const bool joined = fluidJoined(a_surface.m_corner[faceCorner[0]],
                                      a_surface.m_corner[faceCorner[1]],
                                      a_surface.m_corner[faceCorner[2]],
                                      a_surface.m_corner[faceCorner[3]],
                                      centre,
                                      alternating);

      if (alternating && joined) {
        a_surface.m_fluidJoined |= 1 << (2 * dir + side);
      }
    }
  }
#else
  RealVect centre = a_probLo;

  for (int d = 0; d < SpaceDim; d++) {
    centre[d] += a_dx * (static_cast<Real>(a_cell[d]) + 0.5);
  }

  bool alternating = false;

  const bool joined = fluidJoined(a_surface.m_corner[0],
                                  a_surface.m_corner[1],
                                  a_surface.m_corner[3],
                                  a_surface.m_corner[2],
                                  centre,
                                  alternating);

  if (alternating && joined) {
    a_surface.m_fluidJoined = 1;
  }
#endif
}

} // namespace PolyhedralEB

#include <CD_NamespaceFooter.H>
