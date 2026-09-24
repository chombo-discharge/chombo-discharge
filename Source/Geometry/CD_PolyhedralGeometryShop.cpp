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
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <unordered_map>
#include <vector>

#ifdef CH_USE_HDF5
#include <hdf5.h>
#endif

// Chombo includes
#include <BRMeshRefine.H>
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
#include <CD_PolyUtils.H>
#include <CD_PolyhedralGeometryShop.H>
#include <CD_ComputationalGeometry.H>
#include <CD_ParallelOps.H>
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
  m_writeSurface    = false;
#ifndef NDEBUG
  m_sanityCheck = true;
#else
  m_sanityCheck = false;
#endif
  m_profile  = false;
  m_testCopy = false;
  m_verbose  = false;

  // Hidden options, as ScanShop keeps its own.
  ParmParse pp("PolyhedralGeometryShop");

  pp.query("write_surface", m_writeSurface);
  pp.query("sanity_check", m_sanityCheck);
  pp.query("profile", m_profile);
  pp.query("test_copy", m_testCopy);
  pp.query("verbose", m_verbose);

  if (m_verbose) {
    pout() << "PolyhedralGeometryShop::PolyhedralGeometryShop()" << endl;
  }
}

void
PolyhedralGeometryShop::setGrids(const ComputationalGeometry& a_compGeom, const phase::which_phase a_phase) noexcept
{
  if (m_verbose) {
    pout() << "PolyhedralGeometryShop::setGrids" << endl;
  }

  m_compGeom = &a_compGeom;
  m_phase    = a_phase;
}

void
PolyhedralGeometryShop::buildGraphs()
{
  CH_TIME("PolyhedralGeometryShop::buildGraphs");

  if (m_verbose) {
    pout() << "PolyhedralGeometryShop::buildGraphs" << endl;
  }

  if (m_compGeom == nullptr) {
    MayDay::Error("PolyhedralGeometryShop::buildGraphs - setGrids has not been called");
  }

  const int numLevels  = m_compGeom->getNumGridLevels();
  const int startLevel = m_compGeom->getStartLevel();

  m_graphs.resize(numLevels);

  for (int lvl = 0; lvl < numLevels; lvl++) {
    m_graphs[lvl] = RefCountedPtr<PolyhedralEBGraph>(new PolyhedralEBGraph());
  }

  Timer timer("PolyhedralGeometryShop::buildGraphs (" + std::string((m_phase == phase::gas) ? "gas" : "solid") + ")");

  // Each level with cut tiles on its own, from this phase's implicit function.
  for (int lvl = startLevel; lvl < numLevels; lvl++) {
    const Vector<Box>& tiles = m_compGeom->getCutTiles(lvl);

    if (tiles.size() == 0) {
      continue;
    }

    // One ghost cell: the graph's own needs end at the ring around each box, through which it reads its
    // neighbours' surfaces; a consumer wanting a wider ring of cell states fills it from the geometry's
    // classification, which is the answer outside the tiles in any case.
    timer.startEvent("Define level " + std::to_string(lvl));
    m_graphs[lvl]->define(*m_baseIF, tiles, m_compGeom->getDomain(lvl), m_probLo, m_compGeom->getDx(lvl), 1);
    timer.stopEvent("Define level " + std::to_string(lvl));
  }

  // Then every level to the one above it, so the coarse side of each level boundary knows the fine side.
  for (int lvl = startLevel; lvl + 1 < numLevels; lvl++) {
    if (!m_graphs[lvl]->isDefined() || !m_graphs[lvl + 1]->isDefined()) {
      continue;
    }

    timer.startEvent("Link levels " + std::to_string(lvl) + "/" + std::to_string(lvl + 1));
    PolyhedralEBGraph::link(*m_graphs[lvl], *m_graphs[lvl + 1]);
    timer.stopEvent("Link levels " + std::to_string(lvl) + "/" + std::to_string(lvl + 1));
  }

  if (m_profile) {
    timer.eventReport(pout(), false);
  }
}

const PolyhedralEBGraph&
PolyhedralGeometryShop::getGraph(const int a_level) const noexcept
{
  if (m_verbose) {
    pout() << "PolyhedralGeometryShop::getGraph" << endl;
  }

  if (a_level < 0 || a_level >= m_graphs.size()) {
    MayDay::Error("PolyhedralGeometryShop::getGraph - no such level, or buildGraphs has not run");
  }

  return *m_graphs[a_level];
}

void
PolyhedralGeometryShop::verifySurface() const
{
  CH_TIME("PolyhedralGeometryShop::verifySurface");

  if (m_verbose) {
    pout() << "PolyhedralGeometryShop::verifySurface" << endl;
  }

  if (m_compGeom == nullptr) {
    MayDay::Error("PolyhedralGeometryShop::verifySurface - setGrids has not been called");
  }

  if (!m_sanityCheck && !m_writeSurface && !m_testCopy) {
    return;
  }

  const std::string phaseName = (m_phase == phase::gas) ? "gas" : "solid";

  Timer timer("PolyhedralGeometryShop::verifySurface (" + phaseName + ")");

  const int numLevels = m_compGeom->getNumGridLevels();

  if (m_sanityCheck) {
    timer.startEvent("Sanity check");
    this->sanityCheck();
    timer.stopEvent("Sanity check");
  }

  if (m_testCopy) {
    timer.startEvent("Copy test");
    this->testGraphCopy();
    timer.stopEvent("Copy test");
  }

  // The surface itself is written in three dimensions only: the facets are polygons, and in two dimensions the
  // interface of a cell is one chord, which the graph checks as it builds it. What sanityCheck says about cells
  // that share a face holds in both dimensions, and so does the copy test.
#if CH_SPACEDIM == 3
  if (m_writeSurface) {
    Vector<Vector<Real>> facets(numLevels);

    for (int lvl = 0; lvl < numLevels; lvl++) {
      timer.startEvent("Build polyhedra, level " + std::to_string(lvl));
      this->collectFacets(facets[lvl], lvl);
      timer.stopEvent("Build polyhedra, level " + std::to_string(lvl));
    }

    timer.startEvent("Write surface");
    this->writeSurface("surface_mesh_" + phaseName, facets);
    timer.stopEvent("Write surface");
  }

#endif

  if (m_profile) {
    timer.eventReport(pout(), false);
  }
}

void
PolyhedralGeometryShop::collectFacets(Vector<Real>& a_facets, const int a_level) const
{
  CH_TIME("PolyhedralGeometryShop::collectFacets");

  if (m_verbose) {
    pout() << "PolyhedralGeometryShop::collectFacets" << endl;
  }

  using PolyhedralEB::CutCellBody;
  using PolyhedralEB::CutCellSurface;

  // The levels below the start level have no graph and no cut cell of their own to write: the start level covers
  // them whole.
  if (a_level >= m_graphs.size() || !m_graphs[a_level]->isDefined()) {
    return;
  }

#if CH_SPACEDIM == 3
  const PolyhedralEBGraph& graph = *m_graphs[a_level];

  const Real dx = graph.getDx();

  const DisjointBoxLayout&                 grids    = graph.getGrids();
  const LayoutData<IntVectSet>&            cutCells = graph.getCutCells();
  const LevelData<IVSFAB<CutCellSurface>>& surfaces = graph.getSurfaces();
  const LevelData<BaseFab<signed char>>&   refined  = graph.getRefinedMask();

  // Every cut cell of this rank's tiles that the finer level does not carry.
  for (DataIterator dit(grids); dit.ok(); ++dit) {
    const Box                     box        = grids[dit()];
    const IntVectSet&             cut        = cutCells[dit()];
    const IVSFAB<CutCellSurface>& stored     = surfaces[dit()];
    const BaseFab<signed char>&   refinedFab = refined[dit()];

    // in the order a BoxIterator meets the cells, so that the surface's order does not depend on how the set
    // happens to be stored
    for (BoxIterator bit(box); bit.ok(); ++bit) {
      const IntVect iv = bit();

      if (!cut.contains(iv) || refinedFab(iv, 0) != 0) {
        continue;
      }

      CutCellBody body;

      this->defineBody(body, graph, stored, refinedFab, iv);

      body.appendInterfaceFacets(a_facets, iv, m_probLo, dx);
    }
  }
#endif
}

#if CH_SPACEDIM == 3
void
PolyhedralGeometryShop::defineBody(PolyhedralEB::CutCellBody&                  a_body,
                                   const PolyhedralEBGraph&                    a_graph,
                                   const IVSFAB<PolyhedralEB::CutCellSurface>& a_surfaces,
                                   const BaseFab<signed char>&                 a_refined,
                                   const IntVect&                              a_cell) const
{
  CH_TIME("PolyhedralGeometryShop::defineBody");

  using PolyhedralEB::CutCellBody;
  using PolyhedralEB::CutCellSurface;

  const BaseIF& f = *m_baseIF;

  const Real dx     = a_graph.getDx();
  const Box& domain = a_graph.getDomain().domainBox();

  CutCellSurface surface = a_surfaces(a_cell, 0);

  // Which of this cell's edges the finer level describes. An edge is shared by the cells that meet along it, and
  // the crossing on it has to be one point for all of them: a cell that finds it by halving its own edge and a
  // cell that finds it by halving the finer edge land on the same root only to the accuracy of the bisection,
  // which is the square root of the machine epsilon where the surface grazes the edge rather than crossing it. So
  // the question is asked of the edge and not of the cell -- if any cell meeting the edge is refined, every cell
  // meeting it reads the crossing from the finer spacing -- and the four of them then agree by construction. This
  // takes in the cell that meets the finer level along an edge alone, whose faces all lie at this level and which
  // would otherwise never look at the finer spacing at all.
  bool edgeRefined[CutCellSurface::s_numEdges];

  bool anyEdgeRefined = false;

  for (int e = 0; e < CutCellSurface::s_numEdges; e++) {
    const int edgeDir = PolyhedralEB::detail::edgeDirection(e);

    int offset[SpaceDim];
    PolyhedralEB::detail::edgeOrigin(e, offset);

    edgeRefined[e] = false;

    for (int share = 0; share < (1 << (SpaceDim - 1)); share++) {
      IntVect jv = a_cell;

      int bit = 0;

      for (int d = 0; d < SpaceDim; d++) {
        if (d == edgeDir) {
          continue;
        }

        if ((share >> bit) & 1) {
          jv[d] += (offset[d] == 0) ? -1 : 1;
        }

        bit++;
      }

      if (domain.contains(jv) && a_refined.box().contains(jv) && a_refined(jv, 0) != 0) {
        edgeRefined[e] = true;
      }
    }

    anyEdgeRefined = anyEdgeRefined || edgeRefined[e];
  }

  // The children are this cell's own refinement, reconstructed from the implicit function at the finer spacing;
  // their faces on a shared plane coincide, edge for edge and root for root, with those of the finer cells across
  // it, and their edges on this cell's edges are the segments the finer level bisects.
  CutCellSurface children[CutCellSurface::s_numCorners];

  const bool haveChildren = anyEdgeRefined;

  if (haveChildren) {
    // the children share their nodes and edges among themselves at the finer spacing
    const Box fineBox = refine(Box(a_cell, a_cell), 2);

    BaseFab<Real> fineNodeValues;
    BaseFab<Real> fineIntercept[SpaceDim];

    PolyhedralGeometryShop::fillNodeValues(f, fineNodeValues, fineBox, m_probLo, 0.5 * dx);
    PolyhedralGeometryShop::defineIntercepts(fineIntercept, fineBox);

    for (int c = 0; c < CutCellSurface::s_numCorners; c++) {
      IntVect child = 2 * a_cell;

      for (int d = 0; d < SpaceDim; d++) {
        child[d] += (c >> d) & 1;
      }

      PolyhedralGeometryShop::buildSurface(f, fineIntercept, children[c], fineNodeValues, child, m_probLo, 0.5 * dx);
    }

    // Take the crossing on every edge the finer level describes from the half the finer level put it in, in this
    // cell's own parameterisation. An edge the finer level crosses in both halves is a feature no single chord
    // can carry; it keeps this level's crossing here and is reported where the face it belongs to is restricted.
    for (int e = 0; e < CutCellSurface::s_numEdges; e++) {
      if (!edgeRefined[e]) {
        continue;
      }

      const int edgeDir = PolyhedralEB::detail::edgeDirection(e);

      int offset[SpaceDim];
      PolyhedralEB::detail::edgeOrigin(e, offset);

      // A surface carries a crossing on an edge exactly when its two ends disagree, and the body is built on that.
      // The finer level can see a crossing on an edge this level's ends agree about -- the function leaves and
      // re-enters on the way between them -- and that is a pair of crossings, which no single chord carries and
      // which the tiler refines away. Taking the finer level's word for it here would write a crossing onto an
      // edge that cannot hold one and hand the body a surface it cannot close.
      int low  = -1;
      int high = -1;

      PolyhedralEB::detail::edgeCorners(e, low, high);

      if (PolyhedralEB::isFluid(surface.m_corner[low]) == PolyhedralEB::isFluid(surface.m_corner[high])) {
        surface.m_crossing[e] = CutCellSurface::s_noCrossing;

        continue;
      }

      Real half[2] = {CutCellSurface::s_noCrossing, CutCellSurface::s_noCrossing};

      for (int h = 0; h < 2; h++) {
        int which = 0;

        for (int d = 0; d < SpaceDim; d++) {
          which |= ((d == edgeDir) ? h : offset[d]) << d;
        }

        half[h] = children[which].m_crossing[e];
      }

      const bool lowCrossed  = half[0] != CutCellSurface::s_noCrossing;
      const bool highCrossed = half[1] != CutCellSurface::s_noCrossing;

      if (lowCrossed && !highCrossed) {
        surface.m_crossing[e] = 0.5 * half[0];
      }
      else if (highCrossed && !lowCrossed) {
        surface.m_crossing[e] = 0.5 * (1.0 + half[1]);
      }

      // Ends that disagree must carry a crossing; if the finer level found none, this level keeps its own.
    }
  }

  if (!a_body.define(surface)) {
    pout() << "PolyhedralGeometryShop::defineBody - cell " << a_cell << " did not close" << endl;

    MayDay::Error("PolyhedralGeometryShop::defineBody - a cut cell's body did not close");
  }

  // A face shared with a cell the finer level carries is described at the finer level's resolution.
  bool restricted = false;

  for (int dir = 0; dir < SpaceDim; dir++) {
    for (int side = 0; side < 2; side++) {
      const IntVect neighbour = a_cell + (2 * side - 1) * BASISV(dir);

      if (!domain.contains(neighbour) || a_refined(neighbour, 0) == 0) {
        continue;
      }

      // A coarse edge of this face whose ends agree, but whose two finer halves each carry a crossing, is a
      // chord this cell cannot represent: two crossings on one edge is a multi-valued coarse cell.
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

        // A coarse edge whose ends agree must have no crossing in either half; one whose ends disagree has one
        // crossing, and the half the children put it in must be the half this cell's crossing lies in. Anything
        // else is a surface crossing the edge more than once at the finer spacing -- a feature thinner than this
        // cell, which no single chord can represent -- and it stops the run.
        bool consistent = true;

        if (surface.m_crossing[e] == CutCellSurface::s_noCrossing) {
          consistent = !halfCrossed[0] && !halfCrossed[1];
        }
        else {
          const int half = (surface.m_crossing[e] < 0.5) ? 0 : 1;

          consistent = halfCrossed[half] && !halfCrossed[1 - half];
        }

        if (!consistent) {
          pout() << std::setprecision(17) << "PolyhedralGeometryShop::defineBody - cell " << a_cell << ": edge " << e
                 << " of face " << dir << "/" << side << " has crossing " << surface.m_crossing[e]
                 << " but its halves at the finer level have " << halfCrossed[0] << "/" << halfCrossed[1] << endl;

          MayDay::Error("PolyhedralGeometryShop::defineBody - a coarse edge on a level boundary is crossed more "
                        "than once by the finer level");
        }
      }

      // restrictFace takes the children on this face in quadrant order: quadrant q's bits fill the directions
      // other than dir, and dir takes the side.
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

      if (!a_body.restrictFace(faceChildren, dir, side)) {
        pout() << "PolyhedralGeometryShop::defineBody - cell " << a_cell << " could not take face " << dir << "/"
               << side << " from the finer level" << endl;

        MayDay::Error("PolyhedralGeometryShop::defineBody - a face could not be restricted");
      }

      restricted = true;
    }
  }

  if (restricted && !a_body.closeInterface()) {
    pout() << "PolyhedralGeometryShop::defineBody - cell " << a_cell << " did not close after its faces were restricted"
           << endl;

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

    a_body.printPolygons(pout());

    MayDay::Error("PolyhedralGeometryShop::defineBody - a restricted cell's interface did not close");
  }
}
#endif

void
PolyhedralGeometryShop::indexFacets(const Vector<Real>& a_facets,
                                    const Real          a_tolerance,
                                    std::vector<Real>&  a_vertices,
                                    std::vector<int>&   a_connectivity)
{
  CH_TIME("PolyhedralGeometryShop::indexFacets");

  const int numPoints = static_cast<int>(a_facets.size() / 3);

  a_vertices.clear();
  a_connectivity.assign(numPoints, -1);

  // Two positions closer than the tolerance are one vertex. Candidates are found by the cell of a lattice of that
  // spacing that a position falls in: anything within the tolerance of it lies in that cell or in one of the
  // twenty-six around it, so those are the only cells to look in. The lattice narrows the search and nothing
  // more -- every candidate it offers is still measured against, which is what keeps two positions on either
  // side of a lattice wall from being told apart when they are a hair from one another.
  //
  // The cells are held under a mix of their three indices rather than under the indices themselves. A mix
  // collides, and collisions cost nothing here: a cell that answers for two places offers candidates from both,
  // and the distance test throws out the ones that do not belong.
  std::unordered_map<long long, std::vector<int>> buckets;

  buckets.reserve(2 * numPoints);

  const auto cellOf = [&](const Real a_x) -> long long {
    return static_cast<long long>(std::floor(a_x / a_tolerance));
  };

  const auto keyOf = [](const long long a_i, const long long a_j, const long long a_k) -> long long {
    // three odd multipliers, so that neighbouring cells land far apart in the table
    return a_i * 0x9E3779B97F4A7C15LL ^ a_j * 0xC2B2AE3D27D4EB4FLL ^ a_k * 0x165667B19E3779F9LL;
  };

  for (int i = 0; i < numPoints; i++) {
    const Real x = a_facets[3 * i];
    const Real y = a_facets[3 * i + 1];
    const Real z = a_facets[3 * i + 2];

    const long long ci = cellOf(x);
    const long long cj = cellOf(y);
    const long long ck = cellOf(z);

    int found = -1;

    for (int di = -1; di <= 1 && found < 0; di++) {
      for (int dj = -1; dj <= 1 && found < 0; dj++) {
        for (int dk = -1; dk <= 1 && found < 0; dk++) {
          const auto bucket = buckets.find(keyOf(ci + di, cj + dj, ck + dk));

          if (bucket == buckets.end()) {
            continue;
          }

          for (const int candidate : bucket->second) {
            if (std::abs(x - a_vertices[3 * candidate]) <= a_tolerance &&
                std::abs(y - a_vertices[3 * candidate + 1]) <= a_tolerance &&
                std::abs(z - a_vertices[3 * candidate + 2]) <= a_tolerance) {
              found = candidate;

              break;
            }
          }
        }
      }
    }

    if (found < 0) {
      found = static_cast<int>(a_vertices.size() / 3);

      a_vertices.push_back(x);
      a_vertices.push_back(y);
      a_vertices.push_back(z);

      buckets[keyOf(ci, cj, ck)].push_back(found);
    }

    a_connectivity[i] = found;
  }
}

void
PolyhedralGeometryShop::writeSurface(const std::string& a_fileName, const Vector<Vector<Real>>& a_facets) const
{
  CH_TIME("PolyhedralGeometryShop::writeSurface");

  if (m_verbose) {
    pout() << "PolyhedralGeometryShop::writeSurface - writing " << a_fileName << endl;
  }

#ifdef CH_USE_HDF5
  const int numLevels = a_facets.size();

  std::string directory = ".";
  {
    ParmParse pp("Driver");

    pp.query("output_directory", directory);
  }

  const std::string stem = directory + "/geo/" + a_fileName;

  // Every rank writes its own share of every dataset, so the file is opened for parallel access and the writes
  // are collective. Vertices are identified within a rank and not across them: a vertex on the boundary between
  // two ranks is written by both, which costs a little size and saves the communication a global identification
  // would need. Nothing reading an indexed mesh requires otherwise.
  hid_t access = H5Pcreate(H5P_FILE_ACCESS);

#ifdef CH_MPI
  H5Pset_fapl_mpio(access, Chombo_MPI::comm, MPI_INFO_NULL);
#endif

  const hid_t file = H5Fcreate((stem + ".h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, access);

  H5Pclose(access);

  if (file < 0) {
    MayDay::Error("PolyhedralGeometryShop::writeSurface - could not open the file");
  }

  hid_t transfer = H5Pcreate(H5P_DATASET_XFER);

#ifdef CH_MPI
  // Creating the datasets is collective, writing into them is not. Every rank owns a contiguous stretch of each
  // one and nothing else touches it, so there is nothing for a collective transfer to coordinate -- and asking
  // for one is harmful here: a rank whose stretch is empty makes the library decide collective access is not
  // possible for that write, which it then does independently while the others do not, and the next metadata
  // call finds them out of step and blocks.
  H5Pset_dxpl_mpio(transfer, H5FD_MPIO_INDEPENDENT);
#endif

  const auto writeAttribute = [&](const hid_t        a_where,
                                  const std::string& a_name,
                                  const hid_t        a_type,
                                  const hsize_t      a_count,
                                  const void*        a_data) -> void {
    const hid_t space     = (a_count == 1) ? H5Screate(H5S_SCALAR) : H5Screate_simple(1, &a_count, nullptr);
    const hid_t attribute = H5Acreate2(a_where, a_name.c_str(), a_type, space, H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(attribute, a_type, a_data);
    H5Aclose(attribute);
    H5Sclose(space);
  };

  // One dataset, sized by what every rank holds together, with this rank's rows written into its own stretch of
  // it. A rank with nothing to write selects nothing and takes part in the call all the same, which is what
  // collective access asks of it.
  const auto writeSlab = [&](const hid_t        a_where,
                             const std::string& a_name,
                             const hid_t        a_type,
                             const hsize_t      a_rows,
                             const hsize_t      a_offset,
                             const hsize_t      a_total,
                             const hsize_t      a_columns,
                             const void*        a_data) -> void {
    if (a_total == 0) {
      return;
    }

    const hsize_t dims[2] = {a_total, a_columns};
    const hsize_t mine[2] = {a_rows, a_columns};
    const hsize_t at[2]   = {a_offset, 0};

    const int rank = (a_columns > 1) ? 2 : 1;

    const hid_t space   = H5Screate_simple(rank, dims, nullptr);
    const hid_t dataset = H5Dcreate2(a_where, a_name.c_str(), a_type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // A rank with no rows describes that with the null dataspace and selects nothing in the file. It cannot
    // describe it with a simple dataspace of zero extent, which is not one: the call fails, and a rank that
    // fails on its way into a collective write leaves every other rank waiting in it.
    const hid_t memory = (a_rows > 0) ? H5Screate_simple(rank, mine, nullptr) : H5Screate(H5S_NULL);

    if (a_rows > 0) {
      H5Sselect_hyperslab(space, H5S_SELECT_SET, at, nullptr, mine, nullptr);
    }
    else {
      H5Sselect_none(space);
    }

    H5Dwrite(dataset, a_type, memory, space, transfer, a_data);

    H5Sclose(memory);
    H5Dclose(dataset);
    H5Sclose(space);
  };

  const std::string phaseName = (m_phase == phase::gas) ? "gas" : "solid";

  double probLo[3] = {0.0, 0.0, 0.0};

  for (int d = 0; d < SpaceDim; d++) {
    probLo[d] = m_probLo[d];
  }

  writeAttribute(file, "probLo", H5T_NATIVE_DOUBLE, 3, probLo);
  writeAttribute(file, "numLevels", H5T_NATIVE_INT, 1, &numLevels);

  const hid_t nameType = H5Tcopy(H5T_C_S1);
  H5Tset_size(nameType, phaseName.size());
  writeAttribute(file, "phase", nameType, 1, phaseName.c_str());
  H5Tclose(nameType);

  std::vector<long long> totalTriangles(numLevels, 0);
  std::vector<long long> totalVertices(numLevels, 0);

  for (int lvl = 0; lvl < numLevels; lvl++) {
    const std::string name = "level" + std::to_string(lvl);

    const hid_t group = H5Gcreate2(file, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    const double dx        = m_compGeom->getDx(lvl);
    const double tolerance = s_weldSpacing * dx;

    writeAttribute(group, "dx", H5T_NATIVE_DOUBLE, 1, &dx);
    writeAttribute(group, "weldTolerance", H5T_NATIVE_DOUBLE, 1, &tolerance);

    std::vector<Real> vertices;
    std::vector<int>  connectivity;

    PolyhedralGeometryShop::indexFacets(a_facets[lvl], tolerance, vertices, connectivity);

    // A triangle two of whose vertices are the same vertex has collapsed to a line and bounds nothing.
    int collapsed = 0;

    {
      std::vector<int> kept;

      for (size_t i = 0; i + 3 <= connectivity.size(); i += 3) {
        const int a = connectivity[i];
        const int b = connectivity[i + 1];
        const int c = connectivity[i + 2];

        if (a == b || b == c || c == a) {
          collapsed++;

          continue;
        }

        kept.push_back(a);
        kept.push_back(b);
        kept.push_back(c);
      }

      connectivity.swap(kept);
    }

    const long long mineVertices  = vertices.size() / 3;
    const long long mineTriangles = connectivity.size() / 3;

    long long vertexOffset   = 0;
    long long triangleOffset = 0;

    totalVertices[lvl]  = mineVertices;
    totalTriangles[lvl] = mineTriangles;

#ifdef CH_MPI
    MPI_Exscan(&mineVertices, &vertexOffset, 1, MPI_LONG_LONG, MPI_SUM, Chombo_MPI::comm);
    MPI_Exscan(&mineTriangles, &triangleOffset, 1, MPI_LONG_LONG, MPI_SUM, Chombo_MPI::comm);

    if (procID() == 0) {
      vertexOffset   = 0;
      triangleOffset = 0;
    }

    MPI_Allreduce(&mineVertices, &totalVertices[lvl], 1, MPI_LONG_LONG, MPI_SUM, Chombo_MPI::comm);
    MPI_Allreduce(&mineTriangles, &totalTriangles[lvl], 1, MPI_LONG_LONG, MPI_SUM, Chombo_MPI::comm);
#endif

    // This rank's vertices sit at vertexOffset in the file, so its triangles point there too.
    for (size_t i = 0; i < connectivity.size(); i++) {
      connectivity[i] += static_cast<int>(vertexOffset);
    }

    const int totalCollapsed = ParallelOps::sum(collapsed);

    writeAttribute(group, "numCollapsed", H5T_NATIVE_INT, 1, &totalCollapsed);

    writeSlab(group, "vertices", H5T_NATIVE_DOUBLE, mineVertices, vertexOffset, totalVertices[lvl], 3, vertices.data());
    writeSlab(group,
              "connectivity",
              H5T_NATIVE_INT,
              mineTriangles,
              triangleOffset,
              totalTriangles[lvl],
              3,
              connectivity.data());

    // The boxes of the level and the boxes that split on it are the same on every rank, so the master writes them
    // and the others take part with nothing selected.
    const Vector<Box>&                    boxes      = m_compGeom->getBoxes(lvl);
    const Vector<GeometryService::InOut>& gasTypes   = m_compGeom->getTypes(phase::gas, lvl);
    const Vector<GeometryService::InOut>& solidTypes = m_compGeom->getTypes(phase::solid, lvl);

    Vector<int> reasons;

    const Vector<Box>& splitBoxes = m_compGeom->getSplitBoxes(lvl, reasons);

    std::vector<int> corners;
    std::vector<int> gasType;
    std::vector<int> solidType;
    std::vector<int> splitCorners;
    std::vector<int> splitReason;

    const auto appendBox = [](std::vector<int>& a_into, const Box& a_box) -> void {
      for (int d = 0; d < SpaceDim; d++) {
        a_into.push_back(a_box.smallEnd()[d]);
      }
      for (int d = 0; d < SpaceDim; d++) {
        a_into.push_back(a_box.bigEnd()[d]);
      }
    };

    if (procID() == 0) {
      for (int i = 0; i < boxes.size(); i++) {
        appendBox(corners, boxes[i]);

        gasType.push_back(static_cast<int>(gasTypes[i]));
        solidType.push_back(static_cast<int>(solidTypes[i]));
      }

      for (int i = 0; i < splitBoxes.size(); i++) {
        appendBox(splitCorners, splitBoxes[i]);

        splitReason.push_back(reasons[i]);
      }
    }

    const hsize_t mineBoxes  = (procID() == 0) ? boxes.size() : 0;
    const hsize_t mineSplits = (procID() == 0) ? splitBoxes.size() : 0;

    writeSlab(group, "boxes", H5T_NATIVE_INT, mineBoxes, 0, boxes.size(), 2 * SpaceDim, corners.data());
    writeSlab(group, "boxTypeGas", H5T_NATIVE_INT, mineBoxes, 0, boxes.size(), 1, gasType.data());
    writeSlab(group, "boxTypeSolid", H5T_NATIVE_INT, mineBoxes, 0, boxes.size(), 1, solidType.data());
    writeSlab(group, "splitBoxes", H5T_NATIVE_INT, mineSplits, 0, splitBoxes.size(), 2 * SpaceDim, splitCorners.data());
    writeSlab(group, "splitReason", H5T_NATIVE_INT, mineSplits, 0, splitBoxes.size(), 1, splitReason.data());

    H5Gclose(group);
  }

  H5Pclose(transfer);
  H5Fclose(file);

  // A description a viewer can open, pointing at the datasets just written.
  if (procID() == 0) {
    std::ofstream xdmf(stem + ".xmf");

    const std::string base = a_fileName + ".h5";

    xdmf << "<?xml version=\"1.0\" ?>\n";
    xdmf << "<Xdmf Version=\"3.0\">\n  <Domain>\n";
    xdmf << "    <Grid Name=\"" << phaseName << "\" GridType=\"Collection\" CollectionType=\"Spatial\">\n";

    for (int lvl = 0; lvl < numLevels; lvl++) {
      if (totalTriangles[lvl] == 0) {
        continue;
      }

      const std::string name = "level" + std::to_string(lvl);

      xdmf << "      <Grid Name=\"" << name << "\" GridType=\"Uniform\">\n";
      xdmf << "        <Topology TopologyType=\"Triangle\" NumberOfElements=\"" << totalTriangles[lvl] << "\">\n";
      xdmf << "          <DataItem Dimensions=\"" << totalTriangles[lvl] << " 3\" NumberType=\"Int\" Format=\"HDF\">"
           << base << ":/" << name << "/connectivity</DataItem>\n";
      xdmf << "        </Topology>\n";
      xdmf << "        <Geometry GeometryType=\"XYZ\">\n";
      xdmf << "          <DataItem Dimensions=\"" << totalVertices[lvl]
           << " 3\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">" << base << ":/" << name
           << "/vertices</DataItem>\n";
      xdmf << "        </Geometry>\n      </Grid>\n";
    }

    xdmf << "    </Grid>\n  </Domain>\n</Xdmf>\n";
  }
#else
  MayDay::Warning("PolyhedralGeometryShop::writeSurface - built without HDF5, nothing written");
#endif
}

void
PolyhedralGeometryShop::sanityCheck() const
{
  CH_TIME("PolyhedralGeometryShop::sanityCheck");

  if (m_verbose) {
    pout() << "PolyhedralGeometryShop::sanityCheck" << endl;
  }

  this->sanityCheck(m_graphs);
}

void
PolyhedralGeometryShop::sanityCheck(const Vector<RefCountedPtr<PolyhedralEBGraph>>& a_graphs) const
{
  CH_TIME("PolyhedralGeometryShop::sanityCheck(graphs)");

  if (m_verbose) {
    pout() << "PolyhedralGeometryShop::sanityCheck(graphs)" << endl;
  }

  using PolyhedralEB::CutCellBody;
  using PolyhedralEB::CutCellSurface;

  long long numOpen     = 0;
  long long numOverused = 0;
  long long numTouching = 0;
  int       numReported = 0;

  // A regular cell and a covered one cannot share a face: the four nodes of that face belong to both, and they
  // would have to be fluid for one cell and solid for the other. The classification is corner-based, so this holds
  // by construction inside a level -- and the test is here because it is what the edge of the tiled region has to
  // honour as well, where the cells beyond the tiles are classified from their own corners rather than carried:
  // a covered cell at the edge of a tile with a regular cell across from it would need an interface between them
  // that neither cell holds. Cheap, exact and dimension-independent, so it runs on every graph.
  for (int lvl = 0; lvl < a_graphs.size(); lvl++) {
    const PolyhedralEBGraph& graph = *a_graphs[lvl];

    if (!graph.isDefined()) {
      continue;
    }

    const Box& domainBox = graph.getDomain().domainBox();

    const DisjointBoxLayout&               grids  = graph.getGrids();
    const LevelData<BaseFab<signed char>>& states = graph.getCellStates();

    for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box box = grids[dit()];

      const BaseFab<signed char>& state = states[dit()];

      for (BoxIterator bit(box); bit.ok(); ++bit) {
        const IntVect iv = bit();

        for (int dir = 0; dir < SpaceDim; dir++) {
          for (int side = 0; side < 2; side++) {
            const IntVect jv = iv + (2 * side - 1) * BASISV(dir);

            if (!domainBox.contains(jv) || !state.box().contains(jv)) {
              continue;
            }

            const bool touching = (state(iv, 0) == PolyhedralEBGraph::s_regular &&
                                   state(jv, 0) == PolyhedralEBGraph::s_covered) ||
                                  (state(iv, 0) == PolyhedralEBGraph::s_covered &&
                                   state(jv, 0) == PolyhedralEBGraph::s_regular);

            if (touching) {
              numTouching++;

              if (numReported < 10) {
                pout() << "PolyhedralGeometryShop::sanityCheck - level " << lvl << " cell " << iv << " is "
                       << static_cast<int>(state(iv, 0)) << " and its neighbour " << jv << " is "
                       << static_cast<int>(state(jv, 0)) << ", with no cut cell between them" << endl;

                numReported++;
              }
            }
          }
        }
      }
    }
  }

#if CH_SPACEDIM == 3
  // A triangle edge can be shared only by cells that touch, so every edge of a cell's interface must be used
  // exactly twice among the interface triangles of the cell and of the cells in its 3^D neighbourhood -- on this
  // level through the ghost cells the graph holds, and across a level boundary through the finer cells behind a
  // face the finer level describes, reconstructed for the purpose. An edge in the domain boundary is used once;
  // an edge in a face the coarser level describes is checked from that side, where the finer cells are gathered.
  // Vertices are welded onto a lattice a small fraction of the finest spacing wide.
  int finest = -1;

  for (int lvl = 0; lvl < a_graphs.size(); lvl++) {
    if (a_graphs[lvl]->isDefined()) {
      finest = lvl;
    }
  }

  if (finest < 0) {
    return;
  }

  // Two vertices are the same vertex when they are closer than this. Within a level the cells compute a shared
  // vertex from the same nodes and the same edge, so they agree bit for bit; across a level boundary the finer
  // cells behind a face are rebuilt at the finer spacing and agree only to round-off, about 1e-12 of a cell.
  // The distance is what decides, not a lattice: two points on either side of a lattice wall are as close as any
  // other pair, and quantising them apart is how a check of this kind reports a hole that is not there. It is
  // the spacing of the level being checked that sets the scale, since that is the cell the round-off is a
  // fraction of; measuring every level against the finest one asks the coarse levels to agree far closer than
  // the arithmetic that built them can.
  Real tolerance = s_weldSpacing * a_graphs[finest]->getDx();

  using Edge = std::array<int, 2>;

  // The vertices of a neighbourhood's triangles, merged by proximity, and the triangles as vertex numbers.
  std::vector<Real> points;
  std::vector<int>  order;
  std::vector<int>  vertexOf;
  std::vector<Real> vertexPoint;

  auto identify = [&]() -> void {
    const int numPoints = static_cast<int>(points.size() / 3);

    vertexOf.assign(numPoints, -1);

    order.resize(numPoints);

    for (int i = 0; i < numPoints; i++) {
      order[i] = i;
    }

    std::sort(order.begin(), order.end(), [&](const int a, const int b) -> bool {
      return points[3 * a] < points[3 * b];
    });

    int numVertices = 0;

    vertexPoint.clear();

    for (int k = 0; k < numPoints; k++) {
      const int i = order[k];

      // Everything within the tolerance in the first coordinate is adjacent in this order, so the scan back over
      // that window meets every candidate.
      for (int l = k - 1; l >= 0 && points[3 * i] - points[3 * order[l]] <= tolerance; l--) {
        const int j = order[l];

        bool same = true;

        for (int d = 0; d < SpaceDim; d++) {
          same = same && (std::abs(points[3 * i + d] - points[3 * j + d]) <= tolerance);
        }

        if (same) {
          vertexOf[i] = vertexOf[j];

          break;
        }
      }

      if (vertexOf[i] < 0) {
        vertexOf[i] = numVertices++;

        for (int d = 0; d < SpaceDim; d++) {
          vertexPoint.push_back(points[3 * i + d]);
        }
      }
    }
  };

  auto appendPoints = [&](const Vector<Real>& a_facets) -> void {
    for (int i = 0; i < a_facets.size(); i++) {
      points.push_back(a_facets[i]);
    }
  };

  // The triangles of one facet list, as vertex numbers, from the identification above. The facet list must be one
  // of those appendPoints was given, and a_offset its first point.
  auto edgesOf = [&](const int a_offset, const int a_numFacets, std::vector<Edge>& a_edges) -> void {
    for (int t = 0; t < a_numFacets; t++) {
      const int v[3] = {vertexOf[a_offset + 3 * t], vertexOf[a_offset + 3 * t + 1], vertexOf[a_offset + 3 * t + 2]};

      // a triangle two of whose vertices are the same vertex has collapsed to a line and bounds nothing
      if (v[0] == v[1] || v[1] == v[2] || v[2] == v[0]) {
        continue;
      }

      for (int k = 0; k < 3; k++) {
        const int p = std::min(v[k], v[(k + 1) % 3]);
        const int q = std::max(v[k], v[(k + 1) % 3]);

        a_edges.push_back({p, q});
      }
    }
  };

  for (int lvl = 0; lvl < a_graphs.size(); lvl++) {
    const PolyhedralEBGraph& graph = *a_graphs[lvl];

    if (!graph.isDefined()) {
      continue;
    }

    const Real dx     = graph.getDx();
    const Box& domain = graph.getDomain().domainBox();

    tolerance = s_weldSpacing * dx;

    const DisjointBoxLayout&                 grids      = graph.getGrids();
    const LayoutData<IntVectSet>&            cutCells   = graph.getCutCells();
    const LevelData<IVSFAB<CutCellSurface>>& surfaces   = graph.getSurfaces();
    const LevelData<BaseFab<signed char>>&   refined    = graph.getRefinedMask();
    const LevelData<BaseFab<signed char>>&   faceStates = graph.getFaceStates();

    // whether a vertex lies in a face plane of a cell of this level, for reading which plane an edge lies in
    auto inPlane = [&](const int a_vertex, const int a_cellCoordinate, const int a_dir, const int a_side) -> bool {
      const Real plane = m_probLo[a_dir] + dx * static_cast<Real>(a_cellCoordinate + a_side);

      return std::abs(vertexPoint[SpaceDim * a_vertex + a_dir] - plane) <= tolerance;
    };

    for (DataIterator dit(grids); dit.ok(); ++dit) {
      const Box                     box        = grids[dit()];
      const Box                     grown      = grow(box, 1) & domain;
      const IntVectSet&             cut        = cutCells[dit()];
      const IVSFAB<CutCellSurface>& stored     = surfaces[dit()];
      const BaseFab<signed char>&   refinedFab = refined[dit()];
      const BaseFab<signed char>&   faces      = faceStates[dit()];

      // The interface triangles of every cut cell of the box and its one-cell ring that this level writes, built
      // once each and reached through a cell-indexed table.
      BaseFab<int> table(grown, 1);

      table.setVal(-1);

      Vector<Vector<Real>> facets;

      for (BoxIterator cit(grown); cit.ok(); ++cit) {
        const IntVect iv = cit();

        if (!cut.contains(iv) || refinedFab(iv, 0) != 0) {
          continue;
        }

        CutCellBody body;

        this->defineBody(body, graph, stored, refinedFab, iv);

        table(iv, 0) = facets.size();

        facets.push_back(Vector<Real>());

        body.appendInterfaceFacets(facets.back(), iv, m_probLo, dx);
      }

      for (BoxIterator bit(box); bit.ok(); ++bit) {
        const IntVect iv = bit();

        if (table(iv, 0) < 0) {
          continue;
        }

        // The triangles of this cell and of everything in its neighbourhood, gathered before any of them is
        // matched, since a vertex is identified against every other vertex of the neighbourhood at once.
        points.clear();

        std::vector<int> ownOffsets;
        std::vector<int> ownCounts;
        std::vector<int> aroundOffsets;
        std::vector<int> aroundCounts;

        ownOffsets.push_back(0);
        ownCounts.push_back(facets[table(iv, 0)].size() / 9);

        appendPoints(facets[table(iv, 0)]);

        for (BoxIterator nit(grow(Box(iv, iv), 1) & grown); nit.ok(); ++nit) {
          if (table(nit(), 0) >= 0) {
            aroundOffsets.push_back(static_cast<int>(points.size() / 3));
            aroundCounts.push_back(facets[table(nit(), 0)].size() / 9);

            appendPoints(facets[table(nit(), 0)]);
          }
        }

        // The finer cells around this one, reconstructed at their spacing. Every neighbour the finer level
        // carries is taken, not only the six across a face: an interface edge that runs along a cell edge is
        // shared by two facets whose cells meet along that edge alone, so the cell holding the other half of it
        // can be a diagonal neighbour, and leaving those out reports an edge as open that is not.
        {
          Box neighbourhood(iv - IntVect::Unit, iv + IntVect::Unit);

          neighbourhood &= domain;

          for (BoxIterator nit2(neighbourhood); nit2.ok(); ++nit2) {
            const IntVect neighbour = nit2();

            if (neighbour == iv || !refinedFab.box().contains(neighbour) || refinedFab(neighbour, 0) == 0) {
              continue;
            }

            const Box fineBox = refine(Box(neighbour, neighbour), 2);

            BaseFab<Real> fineNodeValues;
            BaseFab<Real> fineIntercept[SpaceDim];

            PolyhedralGeometryShop::fillNodeValues(*m_baseIF, fineNodeValues, fineBox, m_probLo, 0.5 * dx);
            PolyhedralGeometryShop::defineIntercepts(fineIntercept, fineBox);

            for (BoxIterator fit(fineBox); fit.ok(); ++fit) {
              CutCellSurface fineSurface;

              PolyhedralGeometryShop::buildSurface(*m_baseIF,
                                                   fineIntercept,
                                                   fineSurface,
                                                   fineNodeValues,
                                                   fit(),
                                                   m_probLo,
                                                   0.5 * dx);

              if (CutCellBody::classify(fineSurface) != CutCellBody::Kind::Cut) {
                continue;
              }

              CutCellBody fineBody;

              if (!fineBody.define(fineSurface)) {
                continue;
              }

              // The finer level classifies with the volume threshold and the dust rule, and a cell those rules
              // turn into a regular or a covered one writes no surface there. A reconstruction that kept it would
              // hold triangles the level itself does not have, and every edge of them would be reported open.
              if (m_volumeThreshold > 0.0 && fineBody.volumeFraction() < m_volumeThreshold) {
                continue;
              }

              if (PolyhedralGeometryShop::isDust(fineBody)) {
                continue;
              }

              Vector<Real> fineFacets;

              fineBody.appendInterfaceFacets(fineFacets, fit(), m_probLo, 0.5 * dx);

              aroundOffsets.push_back(static_cast<int>(points.size() / 3));
              aroundCounts.push_back(fineFacets.size() / 9);

              appendPoints(fineFacets);
            }
          }
        }

        identify();

        std::vector<Edge> own;
        std::vector<Edge> around;

        for (size_t i = 0; i < ownOffsets.size(); i++) {
          edgesOf(ownOffsets[i], ownCounts[i], own);
        }

        for (size_t i = 0; i < aroundOffsets.size(); i++) {
          edgesOf(aroundOffsets[i], aroundCounts[i], around);
        }

        std::sort(around.begin(), around.end());

        for (const Edge& edge : own) {
          // an edge in a face plane of this cell that is the domain boundary, or that the coarser level
          // describes, is not this cell's to check
          bool exempt = false;

          for (int dir = 0; dir < SpaceDim && !exempt; dir++) {
            for (int side = 0; side < 2 && !exempt; side++) {
              const int state = faces(iv, 2 * dir + side);

              if (state != PolyhedralEBGraph::s_faceBoundary && state != PolyhedralEBGraph::s_faceCoarser) {
                continue;
              }

              exempt = inPlane(edge[0], iv[dir], dir, side) && inPlane(edge[1], iv[dir], dir, side);
            }
          }

          if (exempt) {
            continue;
          }

          const auto range = std::equal_range(around.begin(), around.end(), edge);
          const auto uses  = std::distance(range.first, range.second);

          if (uses == 2) {
            continue;
          }

          if (uses < 2) {
            numOpen++;
          }
          else {
            numOverused++;
          }

          if (numReported < 10) {
            pout() << std::setprecision(17) << "PolyhedralGeometryShop::sanityCheck - level " << lvl << " cell " << iv
                   << ": edge used " << uses << " times:";

            for (int d = 0; d < SpaceDim; d++) {
              pout() << " " << vertexPoint[SpaceDim * edge[0] + d];
            }

            pout() << " ->";

            for (int d = 0; d < SpaceDim; d++) {
              pout() << " " << vertexPoint[SpaceDim * edge[1] + d];
            }

            pout() << endl;

            numReported++;
          }
        }
      }
    }
  }
#endif

  const long long totalOpen     = ParallelOps::sum(numOpen);
  const long long totalOverused = ParallelOps::sum(numOverused);
  const long long totalTouching = ParallelOps::sum(numTouching);

  if (procID() == 0) {
    pout() << "PolyhedralGeometryShop::sanityCheck - " << totalOpen << " interior edges open, " << totalOverused
           << " interior edges used more than twice, " << totalTouching << " regular cells against a covered one"
           << endl;
  }

  if (totalTouching > 0) {
    MayDay::Error("PolyhedralGeometryShop::sanityCheck - a regular cell shares a face with a covered one");
  }

  if (totalOpen > 0 || totalOverused > 0) {
    MayDay::Error("PolyhedralGeometryShop::sanityCheck - the interface is not closed away from the domain boundary");
  }
}

void
PolyhedralGeometryShop::testGraphCopy() const
{
  CH_TIME("PolyhedralGeometryShop::testGraphCopy");

  if (m_verbose) {
    pout() << "PolyhedralGeometryShop::testGraphCopy" << endl;
  }

  // Twice: once onto the tiles split into octants, once onto the tiles as they are; the ranks are shuffled both
  // times, so the second copy moves every box between ranks without changing any box.
  for (int split = 1; split >= 0; split--) {
    const std::string what = (split == 1) ? "octants" : "same tiles";

    Timer timer("PolyhedralGeometryShop::testGraphCopy (" + what + ")");

    Vector<RefCountedPtr<PolyhedralEBGraph>> copies(m_graphs.size());

    timer.startEvent("Copy onto " + what);

    for (int lvl = 0; lvl < m_graphs.size(); lvl++) {
      copies[lvl] = RefCountedPtr<PolyhedralEBGraph>(new PolyhedralEBGraph());

      const PolyhedralEBGraph& graph = *m_graphs[lvl];

      if (!graph.isDefined()) {
        continue;
      }

      // The pieces are given ranks that walk the rank list with a stride coprime to its length, so that a tile's
      // pieces land on different ranks and no rank keeps what it had. The layout covers the same cells as the
      // graph's, which is what the copy needs.
      const DisjointBoxLayout& grids = graph.getGrids();

      Vector<Box> pieces;

      for (LayoutIterator lit = grids.layoutIterator(); lit.ok(); ++lit) {
        const Box tile = grids[lit()];

        if (split == 1) {
          Vector<Box> children;

          domainSplit(tile, children, std::max(1, tile.shortside() / 2), 1);

          pieces.append(children);
        }
        else {
          pieces.push_back(tile);
        }
      }

      int stride = numProc() / 2 + 1;

      while (std::gcd(stride, numProc()) != 1) {
        stride++;
      }

      Vector<int> ranks(pieces.size());

      for (int i = 0; i < pieces.size(); i++) {
        ranks[i] = (numProc() > 1) ? static_cast<int>((static_cast<long>(i) * stride + 1) % numProc()) : 0;
      }

      timer.startEvent("Layout, level " + std::to_string(lvl));
      DisjointBoxLayout shuffled(pieces, ranks, graph.getDomain());

      shuffled.close();
      timer.stopEvent("Layout, level " + std::to_string(lvl));

      timer.startEvent("Copy, level " + std::to_string(lvl));
      copies[lvl]->define(graph, shuffled, *m_baseIF);
      timer.stopEvent("Copy, level " + std::to_string(lvl));
    }

    timer.stopEvent("Copy onto " + what);

    // linked to one another as the originals were, which is what sets the refined masks and the seam faces
    timer.startEvent("Link the copies");

    for (int lvl = 0; lvl + 1 < copies.size(); lvl++) {
      if (copies[lvl]->isDefined() && copies[lvl + 1]->isDefined()) {
        PolyhedralEBGraph::link(*copies[lvl], *copies[lvl + 1]);
      }
    }

    timer.stopEvent("Link the copies");

    // the copies must be a closed surface, and copied back they must be the originals
    timer.startEvent("Check the copies");
    this->sanityCheck(copies);
    timer.stopEvent("Check the copies");

    // back onto the original layouts, linked again, and compared level by level
    Vector<RefCountedPtr<PolyhedralEBGraph>> back(m_graphs.size());

    for (int lvl = 0; lvl < m_graphs.size(); lvl++) {
      back[lvl] = RefCountedPtr<PolyhedralEBGraph>(new PolyhedralEBGraph());

      if (!m_graphs[lvl]->isDefined()) {
        continue;
      }

      timer.startEvent("Copy back, level " + std::to_string(lvl));
      back[lvl]->define(*copies[lvl], m_graphs[lvl]->getGrids(), *m_baseIF);
      timer.stopEvent("Copy back, level " + std::to_string(lvl));
    }

    for (int lvl = 0; lvl + 1 < back.size(); lvl++) {
      if (back[lvl]->isDefined() && back[lvl + 1]->isDefined()) {
        PolyhedralEBGraph::link(*back[lvl], *back[lvl + 1]);
      }
    }

    for (int lvl = 0; lvl < m_graphs.size(); lvl++) {
      const PolyhedralEBGraph& graph = *m_graphs[lvl];

      if (!graph.isDefined()) {
        continue;
      }

      timer.startEvent("Compare, level " + std::to_string(lvl));
      const bool same = back[lvl]->equals(graph);
      timer.stopEvent("Compare, level " + std::to_string(lvl));

      if (!same) {
        pout() << "PolyhedralGeometryShop::testGraphCopy - level " << lvl << " differs after a copy onto " << what
               << " and back" << endl;

        MayDay::Error("PolyhedralGeometryShop::testGraphCopy - the graph did not survive a copy");
      }
    }

    if (procID() == 0) {
      pout() << "PolyhedralGeometryShop::testGraphCopy - every graph is closed on a shuffled layout of " << what
             << " and is unchanged by a copy there and back" << endl;
    }

    if (m_profile) {
      timer.eventReport(pout(), false);
    }
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
PolyhedralGeometryShop::isDust(const PolyhedralEB::CutCellBody& a_body) noexcept
{
  // Exactly, not within a tolerance: a degenerate body has no polygon of positive area at all, so both of these
  // are zero by construction rather than by luck.
  return (1.0 - a_body.volumeFraction() <= 0.0) && (a_body.boundaryArea() <= 0.0);
}

Real
PolyhedralGeometryShop::snappedValue(const BaseIF& a_function, const RealVect& a_point, const Real a_dx) noexcept
{
  const Real value = a_function.value(a_point);

  if (std::abs(value) > s_snapTolerance * a_dx) {
    return value;
  }

  return PolyhedralGeometryShop::resolveTangency(a_function, a_point, a_dx);
}

Real
PolyhedralGeometryShop::resolveTangency(const BaseIF& a_function, const RealVect& a_point, const Real a_dx) noexcept
{
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

  // An edge whose ends are both away from zero and of opposite sign carries one plain root, and Brent's method
  // brackets it in a handful of evaluations where bisection needs tens. An end that reads as exactly zero is a
  // different problem: the function may be zero along a stretch of the edge, and the crossing is then where it
  // leaves that stretch rather than where it reaches it. That is a step in the fluid predicate, which has no
  // continuous root to interpolate, so those edges are walked by bisecting the predicate below.
  if (a_loValue != 0.0 && a_hiValue != 0.0 && ((a_loValue < 0.0) != (a_hiValue < 0.0))) {
    const auto alongEdge = [&](const Real a_t) -> Real {
      RealVect x = lowPoint;

      x[a_dir] += a_dx * a_t;

      return PolyhedralGeometryShop::snappedValue(a_function, x, a_dx);
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
    if (PolyhedralGeometryShop::isDust(body)) {
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
