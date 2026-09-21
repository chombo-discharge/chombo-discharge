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
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

// Chombo includes
#include <BoxIterator.H>
#include <CH_assert.H>
#include <IntVectSet.H>
#include <LoHiSide.H>
#include <MayDay.H>
#include <ParmParse.H>
#include <PolyGeom.H>
#include <SPMD.H>
#include <TreeIntVectSet.H>

// Our includes
#include <CD_PolyhedralEBUtils.H>
#include <CD_PolyhedralGeometryShop.H>
#include <CD_ComputationalGeometry.H>
#include <CD_Timer.H>
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
  m_compGeom        = nullptr;
  m_phase           = phase::gas;
  m_writeSTL        = false;
  m_sanityCheck     = true;
  m_profile         = false;

  // Hidden options, as ScanShop keeps its own.
  ParmParse pp("PolyhedralGeometryShop");

  pp.query("write_stl", m_writeSTL);
  pp.query("sanity_check", m_sanityCheck);
  pp.query("profile", m_profile);
}

void
PolyhedralGeometryShop::setGrids(const ComputationalGeometry& a_compGeom, const phase::which_phase a_phase) noexcept
{
  m_compGeom = &a_compGeom;
  m_phase    = a_phase;
}

void
PolyhedralGeometryShop::verifySurface() const
{
  CH_TIME("PolyhedralGeometryShop::verifySurface");

  if (m_compGeom == nullptr) {
    MayDay::Error("PolyhedralGeometryShop::verifySurface - setGrids has not been called");
  }

  // The surface pass is three-dimensional: the seam restriction and the facets are written for polygons, and
  // in two dimensions the interface of a cell is one chord, which is checked by the graph as it is built.
  if (SpaceDim != 3) {
    return;
  }

  if (!m_sanityCheck && !m_writeSTL) {
    return;
  }

  const std::string phaseName = (m_phase == phase::gas) ? "gas" : "solid";

  Timer timer("PolyhedralGeometryShop::verifySurface (" + phaseName + ")");

  const int numLevels = m_compGeom->getNumGridLevels();

  Vector<Real> composite;

  for (int lvl = 0; lvl < numLevels; lvl++) {
    Vector<Real> levelFacets;

    timer.startEvent("Build polyhedra, level " + std::to_string(lvl));
    this->collectFacets(levelFacets, lvl);
    timer.stopEvent("Build polyhedra, level " + std::to_string(lvl));

    if (m_writeSTL) {
      this->writeSTL("surface_mesh_" + phaseName + ".level" + std::to_string(lvl) + ".stl",
                     phaseName + "_level" + std::to_string(lvl),
                     levelFacets);
    }

    composite.append(levelFacets);
  }

  if (m_writeSTL) {
    timer.startEvent("Write STL");
    this->writeSTL("surface_mesh_" + phaseName + ".stl", phaseName, composite);
    timer.stopEvent("Write STL");

    // The boxes themselves, one file per level, for reading the surface against the grids it was built on.
    if (procID() == 0) {
      for (int lvl = 0; lvl < numLevels; lvl++) {
        std::ofstream out("surface_mesh_boxes.level" + std::to_string(lvl) + ".txt");

        const Vector<Box>&                    boxes      = m_compGeom->getBoxes(lvl);
        const Vector<GeometryService::InOut>& gasTypes   = m_compGeom->getTypes(phase::gas, lvl);
        const Vector<GeometryService::InOut>& solidTypes = m_compGeom->getTypes(phase::solid, lvl);

        for (int i = 0; i < boxes.size(); i++) {
          out << boxes[i].smallEnd() << " " << boxes[i].bigEnd() << " gas " << gasTypes[i] << " solid " << solidTypes[i]
              << "\n";
        }

        // and the boxes that split on this level, with why
        std::ofstream splits("surface_mesh_splits.level" + std::to_string(lvl) + ".txt");

        Vector<int> reasons;

        const Vector<Box>& splitBoxes = m_compGeom->getSplitBoxes(lvl, reasons);

        for (int i = 0; i < splitBoxes.size(); i++) {
          splits << splitBoxes[i].smallEnd() << " " << splitBoxes[i].bigEnd() << " reason " << reasons[i] << "\n";
        }
      }
    }
  }

  if (m_sanityCheck) {
    timer.startEvent("Sanity check");
    this->sanityCheck(composite);
    timer.stopEvent("Sanity check");
  }

  if (m_profile) {
    timer.eventReport(pout(), false);
  }
}

void
PolyhedralGeometryShop::collectFacets(Vector<Real>& a_facets, const int a_level) const
{
  CH_TIME("PolyhedralGeometryShop::collectFacets");

  using PolyhedralEB::CutCellBody;
  using PolyhedralEB::CutCellSurface;

  const BaseIF& f = *m_baseIF;

  const ComputationalGeometry& compGeom = *m_compGeom;

  {
    const int  lvl    = a_level;
    const int  finest = compGeom.getNumGridLevels() - 1;
    const Real dx     = compGeom.getDx(lvl);
    const Box& domain = compGeom.getDomain(lvl).domainBox();

    const Vector<Box>&                    boxes = compGeom.getBoxes(lvl);
    const Vector<GeometryService::InOut>& types = compGeom.getTypes(m_phase, lvl);

    // The cells the next finer level carries, on this level's index space. A cell among them is written from
    // the finer level; a cell next to one of them is on the level boundary and has that face restricted.
    TreeIntVectSet covered;

    if (lvl < finest) {
      const Vector<Box>& finerBoxes = compGeom.getBoxes(lvl + 1);

      for (int j = 0; j < finerBoxes.size(); j++) {
        covered |= coarsen(finerBoxes[j], 2);
      }
    }

    for (int i = procID(); i < boxes.size(); i += numProc()) {
      if (types[i] != GeometryService::Irregular) {
        continue;
      }

      // node values once per node and each crossed edge bisected once, shared by the cells of the box
      BaseFab<Real> nodeValues;
      BaseFab<Real> intercept[SpaceDim];

      PolyhedralGeometryShop::fillNodeValues(f, nodeValues, boxes[i], m_probLo, dx);
      PolyhedralGeometryShop::defineIntercepts(intercept, boxes[i]);

      for (BoxIterator bit(boxes[i]); bit.ok(); ++bit) {
        const IntVect iv = bit();

        if (covered.contains(iv)) {
          continue;
        }

        CutCellSurface surface;

        PolyhedralGeometryShop::buildSurface(f, intercept, surface, nodeValues, iv, m_probLo, dx);

        const CutCellBody::Kind kind = CutCellBody::classify(surface);

        if (kind != CutCellBody::Kind::Cut) {
          continue;
        }

        CutCellBody body;

        if (!body.define(surface)) {
          pout() << "PolyhedralGeometryShop::collectFacets - cell " << iv << " on level " << lvl << " did not close"
                 << endl;

          MayDay::Error("PolyhedralGeometryShop::collectFacets - a cut cell's body did not close");
        }

        // the same thresholds the graph applies: too little fluid is a covered cell, too little solid a
        // regular one, and neither has a surface to write
        if (m_volumeThreshold > 0.0 && body.volumeFraction() < m_volumeThreshold) {
          continue;
        }

        if (PolyhedralGeometryShop::isDust(body, m_volumeThreshold)) {
          continue;
        }

#if CH_SPACEDIM == 3
        // A face shared with a cell the finer level carries is described at the finer level's resolution. The
        // children are this cell's own refinement, reconstructed from the implicit function at the finer
        // spacing; their faces on the shared plane coincide, edge for edge and root for root, with those of the
        // finer cells across it.
        bool restricted = false;

        CutCellSurface children[CutCellSurface::s_numCorners];
        bool           haveChildren = false;

        for (int dir = 0; dir < SpaceDim; dir++) {
          for (int side = 0; side < 2; side++) {
            const IntVect neighbour = iv + (2 * side - 1) * BASISV(dir);

            if (!domain.contains(neighbour) || !covered.contains(neighbour)) {
              continue;
            }

            if (!haveChildren) {
              // the children share their nodes and edges among themselves at the finer spacing
              const Box fineBox = refine(Box(iv, iv), 2);

              BaseFab<Real> fineNodeValues;
              BaseFab<Real> fineIntercept[SpaceDim];

              PolyhedralGeometryShop::fillNodeValues(f, fineNodeValues, fineBox, m_probLo, 0.5 * dx);
              PolyhedralGeometryShop::defineIntercepts(fineIntercept, fineBox);

              for (int c = 0; c < CutCellSurface::s_numCorners; c++) {
                IntVect child = 2 * iv;

                for (int d = 0; d < SpaceDim; d++) {
                  child[d] += (c >> d) & 1;
                }

                PolyhedralGeometryShop::buildSurface(f,
                                                     fineIntercept,
                                                     children[c],
                                                     fineNodeValues,
                                                     child,
                                                     m_probLo,
                                                     0.5 * dx);
              }

              haveChildren = true;
            }

            // A coarse edge of this face whose ends agree, but whose two finer halves each carry a crossing,
            // is a chord this cell cannot represent: two crossings on one edge is a multi-valued coarse cell.
            for (int e = 0; e < CutCellSurface::s_numEdges; e++) {
              int low  = -1;
              int high = -1;

              PolyhedralEB::detail::edgeCorners(e, low, high);

              const bool onFace = (((low >> dir) & 1) == side) && (((high >> dir) & 1) == side);

              if (!onFace) {
                continue;
              }

              const int edgeDir = PolyhedralEB::detail::edgeDirection(e);

              int offset[SpaceDim];
              PolyhedralEB::detail::edgeOrigin(e, offset);

              // The two halves of this edge, as the children see them: which halves carry a crossing.
              bool halfCrossed[2] = {false, false};

              for (int half = 0; half < 2; half++) {
                int which = 0;

                for (int d = 0; d < SpaceDim; d++) {
                  which |= ((d == edgeDir) ? half : offset[d]) << d;
                }

                halfCrossed[half] = (children[which].m_crossing[e] != CutCellSurface::s_noCrossing);
              }

              // A coarse edge whose ends agree must have no crossing in either half; one whose ends disagree
              // has one crossing, and the half the children put it in must be the half this cell's crossing
              // lies in. Anything else is a surface crossing the edge more than once at the finer spacing --
              // a feature thinner than this cell, which no single chord can represent -- and it stops the run.
              bool consistent = true;

              if (surface.m_crossing[e] == CutCellSurface::s_noCrossing) {
                consistent = !halfCrossed[0] && !halfCrossed[1];
              }
              else {
                const int half = (surface.m_crossing[e] < 0.5) ? 0 : 1;

                consistent = halfCrossed[half] && !halfCrossed[1 - half];
              }

              if (!consistent) {
                pout() << std::setprecision(17) << "PolyhedralGeometryShop::collectFacets - cell " << iv << " on level "
                       << lvl << ": edge " << e << " of face " << dir << "/" << side << " has crossing "
                       << surface.m_crossing[e] << " but its halves at the finer level have " << halfCrossed[0] << "/"
                       << halfCrossed[1] << endl;

                MayDay::Error("PolyhedralGeometryShop::collectFacets - a coarse edge on a level boundary is crossed "
                              "more than once by the finer level");
              }
            }

            // restrictFace takes the children on this face in quadrant order: quadrant q's bits fill the
            // directions other than dir, and dir takes the side.
            CutCellSurface faceChildren[1 << (SpaceDim - 1)];

            for (int q = 0; q < (1 << (SpaceDim - 1)); q++) {
              int which = 0;
              int bit   = 0;

              for (int d = 0; d < SpaceDim; d++) {
                if (d == dir) {
                  which |= side << d;
                }
                else {
                  which |= ((q >> bit) & 1) << d;
                  bit++;
                }
              }

              faceChildren[q] = children[which];
            }

            if (!body.restrictFace(faceChildren, dir, side)) {
              pout() << "PolyhedralGeometryShop::collectFacets - cell " << iv << " on level " << lvl
                     << " could not take face " << dir << "/" << side << " from the finer level" << endl;

              MayDay::Error("PolyhedralGeometryShop::collectFacets - a face could not be restricted");
            }

            restricted = true;
          }
        }

        if (restricted && !body.closeInterface()) {
          pout() << "PolyhedralGeometryShop::collectFacets - cell " << iv << " on level " << lvl
                 << " did not close after its faces were restricted" << endl;

          pout() << std::setprecision(17) << "  corners:";

          for (int c = 0; c < CutCellSurface::s_numCorners; c++) {
            pout() << " " << surface.m_corner[c];
          }

          pout() << endl << "  crossings:";

          for (int e = 0; e < CutCellSurface::s_numEdges; e++) {
            pout() << " " << surface.m_crossing[e];
          }

          pout() << endl;

          for (int c = 0; c < CutCellSurface::s_numCorners; c++) {
            pout() << "  child " << c << " corners:";

            for (int k = 0; k < CutCellSurface::s_numCorners; k++) {
              pout() << " " << children[c].m_corner[k];
            }

            pout() << " crossings:";

            for (int e = 0; e < CutCellSurface::s_numEdges; e++) {
              pout() << " " << children[c].m_crossing[e];
            }

            pout() << endl;
          }

          body.printPolygons(pout());

          MayDay::Error("PolyhedralGeometryShop::collectFacets - a restricted cell's interface did not close");
        }

        body.appendInterfaceFacets(a_facets, iv, m_probLo, dx);
#endif
      }
    }
  }
}

void
PolyhedralGeometryShop::writeSTL(const std::string&  a_fileName,
                                 const std::string&  a_name,
                                 const Vector<Real>& a_facets) const
{
  CH_TIME("PolyhedralGeometryShop::writeSTL");

  Vector<Vector<Real>> everyone;

  gather(everyone, a_facets, 0);

  if (procID() != 0) {
    return;
  }

  std::ofstream out(a_fileName);

  if (!out.good()) {
    MayDay::Error("PolyhedralGeometryShop::writeSTL - could not open the file");
  }

  out << std::scientific << std::setprecision(17) << "solid " << a_name << "\n";

  for (int rank = 0; rank < everyone.size(); rank++) {
    const Vector<Real>& rankFacets = everyone[rank];

    for (int i = 0; i + 9 <= rankFacets.size(); i += 9) {
      const RealVect a(D_DECL(rankFacets[i + 0], rankFacets[i + 1], rankFacets[i + 2]));
      const RealVect b(D_DECL(rankFacets[i + 3], rankFacets[i + 4], rankFacets[i + 5]));
      const RealVect c(D_DECL(rankFacets[i + 6], rankFacets[i + 7], rankFacets[i + 8]));

      RealVect n = PolyGeom::cross(b - a, c - a);

      if (n.vectorLength() > 0.0) {
        n /= n.vectorLength();
      }

      out << "  facet normal";

      for (int d = 0; d < 3; d++) {
        out << " " << ((d < SpaceDim) ? n[d] : 0.0);
      }

      out << "\n    outer loop\n";

      for (int v = 0; v < 3; v++) {
        out << "      vertex";

        for (int d = 0; d < 3; d++) {
          out << " " << ((d < SpaceDim) ? rankFacets[i + 3 * v + d] : 0.0);
        }

        out << "\n";
      }

      out << "    endloop\n  endfacet\n";
    }
  }

  out << "endsolid " << a_name << "\n";
}

void
PolyhedralGeometryShop::sanityCheck(const Vector<Real>& a_facets) const
{
  CH_TIME("PolyhedralGeometryShop::sanityCheck");

  Vector<Vector<Real>> everyone;

  gather(everyone, a_facets, 0);

  if (procID() != 0) {
    return;
  }

  // Vertices are welded onto a lattice a small fraction of the finest spacing wide, and an edge is the
  // ordered pair of its two welded vertices. Sorting the edges then puts the uses of one edge together.
  const int  finest  = m_compGeom->getNumGridLevels() - 1;
  const Real spacing = s_weldSpacing * m_compGeom->getDx(finest);

  const Box&     finestBox = m_compGeom->getDomain(finest).domainBox();
  const RealVect probLo    = m_probLo;

  RealVect probHi = m_probLo;

  for (int d = 0; d < SpaceDim; d++) {
    probHi[d] += m_compGeom->getDx(finest) * static_cast<Real>(finestBox.size(d));
  }

  auto weld = [&](const Real* a_x) -> std::array<long long, 3> {
    std::array<long long, 3> key = {0, 0, 0};

    for (int d = 0; d < SpaceDim; d++) {
      key[d] = std::llround((a_x[d] - probLo[d]) / spacing);
    }

    return key;
  };

  // The planes of the domain boundary a welded vertex lies in, one bit per plane: 2 * d for the low face in
  // direction d and 2 * d + 1 for the high one. An edge is on the boundary when its two ends share a plane.
  auto boundaryPlanes = [&](const std::array<long long, 3>& a_key) -> int {
    int planes = 0;

    for (int d = 0; d < SpaceDim; d++) {
      const long long hi = std::llround((probHi[d] - probLo[d]) / spacing);

      if (a_key[d] == 0) {
        planes |= 1 << (2 * d);
      }

      if (a_key[d] == hi) {
        planes |= 1 << (2 * d + 1);
      }
    }

    return planes;
  };

  // One entry per triangle edge: the two welded vertices in a fixed order.
  std::vector<std::array<long long, 6>> edges;

  for (int rank = 0; rank < everyone.size(); rank++) {
    const Vector<Real>& facets = everyone[rank];

    for (int i = 0; i + 9 <= facets.size(); i += 9) {
      std::array<long long, 3> v[3];

      for (int k = 0; k < 3; k++) {
        v[k] = weld(&facets[i + 3 * k]);
      }

      // a triangle two of whose vertices weld together has collapsed to a line and bounds nothing
      if (v[0] == v[1] || v[1] == v[2] || v[2] == v[0]) {
        continue;
      }

      for (int k = 0; k < 3; k++) {
        std::array<long long, 3> a = v[k];
        std::array<long long, 3> b = v[(k + 1) % 3];

        if (b < a) {
          std::swap(a, b);
        }

        edges.push_back({a[0], a[1], a[2], b[0], b[1], b[2]});
      }
    }
  }

  std::sort(edges.begin(), edges.end());

  long long numOpen     = 0;
  long long numOverused = 0;
  int       numReported = 0;

  for (std::size_t i = 0; i < edges.size();) {
    std::size_t j = i;

    while (j < edges.size() && edges[j] == edges[i]) {
      j++;
    }

    const std::size_t uses = j - i;

    if (uses != 2) {
      const std::array<long long, 3> a = {edges[i][0], edges[i][1], edges[i][2]};
      const std::array<long long, 3> b = {edges[i][3], edges[i][4], edges[i][5]};

      const bool boundary = (boundaryPlanes(a) & boundaryPlanes(b)) != 0;

      if (!boundary) {
        if (uses == 1) {
          numOpen++;
        }
        else {
          numOverused++;
        }

        if (numReported < 10) {
          pout() << std::setprecision(17) << "PolyhedralGeometryShop::sanityCheck - edge used " << uses << " times:";

          for (int d = 0; d < SpaceDim; d++) {
            pout() << " " << probLo[d] + spacing * static_cast<Real>(a[d]);
          }

          pout() << " ->";

          for (int d = 0; d < SpaceDim; d++) {
            pout() << " " << probLo[d] + spacing * static_cast<Real>(b[d]);
          }

          pout() << endl;

          numReported++;
        }
      }
    }

    i = j;
  }

  pout() << "PolyhedralGeometryShop::sanityCheck - " << edges.size() / 3 << " triangles, " << numOpen
         << " interior edges open, " << numOverused << " interior edges used more than twice" << endl;

  if (numOpen > 0 || numOverused > 0) {
    MayDay::Error("PolyhedralGeometryShop::sanityCheck - the interface is not closed away from the domain boundary");
  }
}

PolyhedralGeometryShop::~PolyhedralGeometryShop()
{}
void
PolyhedralGeometryShop::fillNodeValues(const BaseIF&   a_function,
                                       BaseFab<Real>&  a_nodeValues,
                                       const Box&      a_region,
                                       const RealVect& a_probLo,
                                       const Real&     a_dx)
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

    a_nodeValues(iv, 0) = PolyhedralGeometryShop::snappedValue(a_function, x, a_dx);
  }
}

void
PolyhedralGeometryShop::defineIntercepts(BaseFab<Real> a_intercept[SpaceDim], const Box& a_region)
{
  for (int dir = 0; dir < SpaceDim; dir++) {
    Box edgeBox = a_region;
    edgeBox.surroundingNodes();
    edgeBox.enclosedCells(dir);

    a_intercept[dir].define(edgeBox, 1);
    a_intercept[dir].setVal(PolyhedralEB::CutCellSurface::s_noCrossing);
  }
}

Real
PolyhedralGeometryShop::edgeCrossing(const BaseIF&   a_function,
                                     BaseFab<Real>   a_intercept[SpaceDim],
                                     const IntVect&  a_cell,
                                     const int       a_edge,
                                     const Real      a_lo,
                                     const Real      a_hi,
                                     const RealVect& a_probLo,
                                     const Real&     a_dx)
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
  const Real root = PolyhedralGeometryShop::edgeRoot(a_function, edgeIV, dir, a_lo, a_hi, a_probLo, a_dx);

  if (a_intercept[dir].box().contains(edgeIV)) {
    a_intercept[dir](edgeIV, 0) = root;
  }

  return root;
}

bool
PolyhedralGeometryShop::isDust(const PolyhedralEB::CutCellBody& a_body, const Real a_threshold) noexcept
{
  return (a_threshold > 0.0) && (1.0 - a_body.volumeFraction() < a_threshold) && (a_body.boundaryArea() < a_threshold);
}

Real
PolyhedralGeometryShop::snappedValue(const BaseIF& a_function, const RealVect& a_point, const Real a_dx) noexcept
{
  const Real value = a_function.value(a_point);

  return (std::abs(value) <= s_snapTolerance * a_dx) ? 0.0 : value;
}

Real
PolyhedralGeometryShop::edgeRoot(const BaseIF&   a_function,
                                 const IntVect&  a_edgeIV,
                                 const int       a_dir,
                                 const Real      a_loValue,
                                 const Real      a_hiValue,
                                 const RealVect& a_probLo,
                                 const Real      a_dx) noexcept
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

    const Real value = PolyhedralGeometryShop::snappedValue(a_function, x, a_dx);

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
PolyhedralGeometryShop::fillCorners(PolyhedralEB::CutCellSurface& a_surface,
                                    const BaseFab<Real>&          a_nodeValues,
                                    const IntVect&                a_cell)
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
PolyhedralGeometryShop::buildSurface(const BaseIF&                 a_function,
                                     BaseFab<Real>                 a_intercept[SpaceDim],
                                     PolyhedralEB::CutCellSurface& a_surface,
                                     const BaseFab<Real>&          a_nodeValues,
                                     const IntVect&                a_cell,
                                     const RealVect&               a_probLo,
                                     const Real&                   a_dx)
{
  a_surface = PolyhedralEB::CutCellSurface();

  PolyhedralGeometryShop::fillCorners(a_surface, a_nodeValues, a_cell);

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
      a_surface.m_crossing[e] = PolyhedralGeometryShop::edgeCrossing(a_function,
                                                                     a_intercept,
                                                                     a_cell,
                                                                     e,
                                                                     loValue,
                                                                     hiValue,
                                                                     a_probLo,
                                                                     a_dx);
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
  PolyhedralGeometryShop::fillNodeValues(*m_baseIF, nodeValues, a_ghostRegion, a_probLo, a_dx);

  IntVectSet irregularCells;

  for (BoxIterator bit(a_ghostRegion); bit.ok(); ++bit) {
    const IntVect iv = bit();

    PolyhedralEB::CutCellSurface surface;
    PolyhedralGeometryShop::fillCorners(surface, nodeValues, iv);

    switch (PolyhedralEB::CutCellBody::classify(surface)) {
    case PolyhedralEB::CutCellBody::Kind::Covered: {
      a_regIrregCovered(iv, 0) = -1;

      break;
    }
    case PolyhedralEB::CutCellBody::Kind::Regular: {
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

  PolyhedralGeometryShop::defineIntercepts(intercept, grow(a_validRegion, 1) & a_ghostRegion);

  IntVectSet droppedCells;

  for (IVSIterator ivsIt(irregularCells); ivsIt.ok(); ++ivsIt) {
    const IntVect iv = ivsIt();

    PolyhedralEB::CutCellSurface surface;
    PolyhedralGeometryShop::buildSurface(*m_baseIF, intercept, surface, nodeValues, iv, a_probLo, a_dx);

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

    // the mirror image: a body with next to no solid in it is a regular cell
    if (PolyhedralGeometryShop::isDust(body, m_volumeThreshold)) {
      a_regIrregCovered(iv, 0) = 1;

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
