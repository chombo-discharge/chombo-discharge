/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_PolyhedralEBGraph.cpp
 * @brief  Implementation of CD_PolyhedralEBGraph.H
 * @author Robert Marskar
 */

// Std includes
#include <algorithm>
#include <iomanip>

// Chombo includes
#include <BoxIterator.H>
#include <CH_Timer.H>
#include <Copier.H>
#include <MayDay.H>

// Our includes
#include <CD_PolyhedralEBGraph.H>
#include <CD_PolyhedralGeometryShop.H>
#include <CD_CutCellBody.H>
#include <CD_LoadBalancing.H>
#include <CD_ParallelOps.H>
#include <CD_NamespaceHeader.H>

PolyhedralEBGraph::PolyhedralEBGraph()
  : m_isDefined(false), m_numGhost(0), m_dx(0.0), m_probLo(RealVect::Zero), m_domain(ProblemDomain())
{
  CH_TIME("PolyhedralEBGraph::PolyhedralEBGraph");
}

PolyhedralEBGraph::~PolyhedralEBGraph()
{
  CH_TIME("PolyhedralEBGraph::~PolyhedralEBGraph");
}

void
PolyhedralEBGraph::define(const BaseIF&        a_function,
                          const Vector<Box>&   a_cutTiles,
                          const ProblemDomain& a_domain,
                          const RealVect&      a_probLo,
                          const Real           a_dx,
                          const int            a_numGhost,
                          const Vector<Box>&   a_covered)
{
  CH_TIME("PolyhedralEBGraph::define");

  if (a_dx <= 0.0) {
    MayDay::Error("PolyhedralEBGraph::define - the grid spacing must be positive");
  }
  if (a_numGhost < 1) {
    MayDay::Error("PolyhedralEBGraph::define - at least one ghost cell is needed to see across a box boundary");
  }

  m_domain   = a_domain;
  m_probLo   = a_probLo;
  m_dx       = a_dx;
  m_numGhost = a_numGhost;
  m_covered  = a_covered;

  this->defineGrids(a_cutTiles);
  this->markRefined();

  // one on every cell some tile of this level carries, zero elsewhere, with one ghost cell
  LevelData<BaseFab<signed char>> carried;

  this->markCarried(carried);

  this->defineCells(a_function, carried);
  this->defineOuterFaces(carried);

  m_isDefined = true;
}

void
PolyhedralEBGraph::defineGrids(const Vector<Box>& a_cutTiles)
{
  CH_TIME("PolyhedralEBGraph::defineGrids");

  // Balanced by cell count, as AmrMesh balances its grids before it knows anything better.
  Vector<long> loads(a_cutTiles.size());

  for (int i = 0; i < a_cutTiles.size(); i++) {
    loads[i] = a_cutTiles[i].numPts();
  }

  Vector<int> ranks;

  LoadBalancing::makeBalance(ranks, loads, a_cutTiles);

  m_grids.define(a_cutTiles, ranks, m_domain);
  m_grids.close();

  this->defineData();
}

void
PolyhedralEBGraph::defineData()
{
  CH_TIME("PolyhedralEBGraph::defineData");

  m_cutCells.define(m_grids);
  m_cellStates.define(m_grids, 1, m_numGhost * IntVect::Unit);
  m_faceStates.define(m_grids, 2 * SpaceDim, IntVect::Zero);
  m_refined.define(m_grids, 1, std::max(2, m_numGhost) * IntVect::Unit);

  // Each box's cut-cell set is a bitmap over the box and its ghost ring, which is where its cells come from: a
  // set that starts empty would be a tree, several kilobytes for a few hundred scattered cells, and the surface
  // container keeps a copy of it.
  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    const Box grown = grow(m_grids[dit()], m_numGhost) & m_domain.domainBox();

    m_cutCells[dit()] = IntVectSet(DenseIntVectSet(grown, false));
  }
}

void
PolyhedralEBGraph::findUnresolvedCells(const BaseFab<Real>& a_nodeValues,
                                       const Box&           a_region,
                                       const Vector<Box>&   a_covered,
                                       BaseFab<bool>&       a_unresolved)
{
  CH_TIME("PolyhedralEBGraph::findUnresolvedCells");

  using PolyhedralEB::CutCellBody;
  using PolyhedralEB::CutCellSurface;

  a_unresolved.setVal(false);

  // What the finer level carries here. Such a cell is described up there, so whatever this level makes of it is
  // not what gets used, and filling it would coarsen a feature the mesh has already resolved.
  BaseFab<bool> covered(a_region, 1);

  covered.setVal(false);

  for (int i = 0; i < a_covered.size(); i++) {
    const Box overlap = a_covered[i] & a_region;

    if (!overlap.isEmpty()) {
      covered.setVal(true, overlap, 0);
    }
  }

  // Only the combinatorics matter, so the crossings are placed at the middle of the edges that carry one rather
  // than being solved for: where an edge carries a crossing follows from its ends, and that is all the sheet
  // count reads.
  for (BoxIterator bit(a_region); bit.ok(); ++bit) {
    const IntVect iv = bit();

    if (covered(iv, 0)) {
      continue;
    }

    CutCellSurface surface;

    PolyhedralGeometryShop::fillCorners(surface, a_nodeValues, iv);

    for (int e = 0; e < CutCellSurface::s_numEdges; e++) {
      int low  = 0;
      int high = 0;

      PolyhedralEB::detail::edgeCorners(e, low, high);

      const bool crosses = PolyhedralEB::isFluid(surface.m_corner[low]) !=
                           PolyhedralEB::isFluid(surface.m_corner[high]);

      surface.m_crossing[e] = crosses ? 0.5 : CutCellSurface::s_noCrossing;
    }

    if (CutCellBody::numSheets(surface) > 1) {
      a_unresolved(iv, 0) = true;
    }
  }
}

void
PolyhedralEBGraph::defineCells(const BaseIF& a_function, const LevelData<BaseFab<signed char>>& a_carried)
{
  CH_TIME("PolyhedralEBGraph::defineCells");

  using PolyhedralEB::CutCellBody;
  using PolyhedralEB::CutCellSurface;

  // The surfaces of a box's cut cells are kept in the order a BoxIterator meets the cells, and moved into the
  // level container once every box's cut-cell set is known; the container is defined over those sets.
  LayoutData<Vector<CutCellSurface>> kept(m_grids);

  long long numFilled = 0;

  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    const Box box = m_grids[dit()];

    BaseFab<signed char>& states = m_cellStates[dit()];
    BaseFab<signed char>& faces  = m_faceStates[dit()];
    IntVectSet&           cut    = m_cutCells[dit()];

    Vector<CutCellSurface>& surfaces = kept[dit()];

    // Ghost cells are filled afterwards: by exchange where another tile carries them, from the function where
    // none does. Until then, and outside the domain for good, they are regular so that nothing reads them as a
    // boundary of the fluid.
    states.setVal(s_regular);
    faces.setVal(s_faceClosed);

    const BaseFab<signed char>& refined = m_refined[dit()];

    // Node values once per node and each crossed edge bisected once, shared by the cells of the box, as the
    // generator does when it builds a box.
    BaseFab<Real> nodeValues;
    BaseFab<Real> intercept[SpaceDim];

    // One cell wider than the box, so that every node of the box has all the cells around it to hand: the snap
    // below decides a node from the cells it is a corner of, and a node on the box's boundary has cells on the
    // other side. Both boxes then reach the same verdict for the nodes they share without being told.
    const Box grown = grow(box, 1) & m_domain.domainBox();

    PolyhedralGeometryShop::fillNodeValues(a_function, nodeValues, grown, m_probLo, m_dx);
    PolyhedralGeometryShop::defineIntercepts(intercept, box);

    // Which cells cannot be described at all, decided before any of them is classified.
    BaseFab<bool> unresolved(grown, 1);

    PolyhedralEBGraph::findUnresolvedCells(nodeValues, grown, m_covered, unresolved);

    // Counted over the cells this box owns, not over the ring it also filled: neighbouring boxes reach into
    // one another's ring, and a cell counted there would be reported once per box that reaches it.
    for (BoxIterator bit(box); bit.ok(); ++bit) {
      if (unresolved(bit(), 0)) {
        numFilled++;
      }
    }

    for (BoxIterator bit(box); bit.ok(); ++bit) {
      const IntVect iv = bit();

      CutCellSurface surface;

      // A cell the finer level carries is described there. It is classified from its corners so that its
      // neighbours read a state for it, but no crossing is solved for, no body is built, and it is not a cut cell
      // of this level; its faces stay closed, since what opens onto it is the finer level's to say.
      if (refined(iv, 0) != 0) {
        PolyhedralGeometryShop::fillCorners(surface, nodeValues, iv);

        const CutCellBody::Kind corners = CutCellBody::classify(surface);

        states(iv, 0) = (corners == CutCellBody::Kind::Covered) ? s_covered
                        : (corners == CutCellBody::Kind::Cut)   ? s_cut
                                                                : s_regular;

        continue;
      }

      PolyhedralGeometryShop::buildSurface(a_function, intercept, surface, nodeValues, iv, m_probLo, m_dx);

      // A cell the surface enters as more than one sheet holds a feature thinner than itself, and one body and
      // one interface cannot describe it: read from the nodes, a plate through the middle comes out as two
      // slivers hugging opposite edges, and the fluid runs straight through a barrier that should stop it. The
      // cell is filled, which is the only reading that stays single valued and keeps the barrier a barrier, and
      // it errs toward blocking rather than leaking. Unlike moving a node, this changes no value another level
      // reads, so the children a coarse cell restricts against still agree with it about every edge.
      const CutCellBody::Kind kind = unresolved(iv, 0) ? CutCellBody::Kind::Covered : CutCellBody::classify(surface);

      // A cell filled next door leaves this one with a face onto nothing. A cell whose corners make it regular
      // is then not regular at all: it is full, but that face is closed and the fluid it used to open onto is
      // gone, so what closes the cell there is interface lying in the plane of that face. It is kept as a cut
      // cell holding the whole of itself, which is what the index space does through
      // GeometryShop::fixRegularCellsNextToCovered, and the closure check reads a regular cell against a covered
      // one as a fault for exactly this reason.
      bool nextToFilled = false;

      for (int dir = 0; dir < SpaceDim && !nextToFilled; dir++) {
        for (int side = 0; side < 2 && !nextToFilled; side++) {
          const IntVect other = iv + (2 * side - 1) * BASISV(dir);

          nextToFilled = unresolved.box().contains(other) && unresolved(other, 0);
        }
      }

      int state = s_regular;

      CutCellBody body;

      if (kind == CutCellBody::Kind::Covered) {
        state = s_covered;
      }
      else if (kind == CutCellBody::Kind::Cut || nextToFilled) {
        if (!body.define(surface)) {
          pout() << "PolyhedralEBGraph::defineCells - cell " << iv << " did not close" << endl;

          MayDay::Error("PolyhedralEBGraph::defineCells - a cut cell's body did not close");
        }

        // A body with nothing in it is a regular cell. Nothing else is discarded: a cell's interface is what
        // closes the surface against its neighbours, so a cell dropped for holding little would leave a hole
        // exactly the size of what it held, and no repair on the neighbours' faces can put it back. The volume
        // threshold that keeps such cells out of the index space is applied where the index space is built,
        // which is where its reason -- a solver that would rather not see a cell of no volume -- applies.
        if (PolyhedralGeometryShop::isDust(body) && !nextToFilled) {
          state = s_regular;
        }
        else {
          state = s_cut;
        }
      }

      states(iv, 0) = state;

      // The apertures decide whether a face is open; what lies across an open face is decided afterwards.
      for (int dir = 0; dir < SpaceDim; dir++) {
        for (int side = 0; side < 2; side++) {
          bool open = false;

          const IntVect other = iv + (2 * side - 1) * BASISV(dir);

          const bool ontoFilled = unresolved.box().contains(other) && unresolved(other, 0);

          if (ontoFilled) {
            open = false;
          }
          else if (state == s_regular) {
            open = true;
          }
          else if (state == s_cut) {
            open = body.areaFraction(dir, (side == 0) ? Side::Lo : Side::Hi) > 0.0;
          }

          faces(iv, 2 * dir + side) = open ? s_faceSameLevel : s_faceClosed;
        }
      }

      if (state == s_cut) {
        cut |= iv;

        surfaces.push_back(surface);
      }
    }
  }

  // ghost cells another tile carries take that tile's state, the rest are classified from the function, and the
  // ghost cut cells enter the sets
  m_cellStates.exchange();

  this->defineGhostCells(a_function, a_carried);

  const long long totalFilled = ParallelOps::sum(numFilled);

  if (totalFilled > 0 && procID() == 0) {
    pout() << "PolyhedralEBGraph::defineCells - filled " << totalFilled
           << " cells holding a feature thinner than themselves" << endl;
  }

  m_surfaces.define(m_grids, 1, m_numGhost * IntVect::Unit, IVSFABFactory<CutCellSurface>(m_cutCells));

  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    const Box         box      = m_grids[dit()];
    const IntVectSet& cut      = m_cutCells[dit()];
    const auto&       surfaces = kept[dit()];

    IVSFAB<CutCellSurface>& stored = m_surfaces[dit()];

    int next = 0;

    for (BoxIterator bit(box); bit.ok(); ++bit) {
      const IntVect iv = bit();

      if (cut.contains(iv)) {
        stored(iv, 0) = surfaces[next++];
      }
    }

    CH_assert(next == surfaces.size());
  }

  m_surfaces.exchange();
}

void
PolyhedralEBGraph::defineGhostCells(const BaseIF& a_function, const LevelData<BaseFab<signed char>>& a_carried)
{
  CH_TIME("PolyhedralEBGraph::defineGhostCells");

  using PolyhedralEB::CutCellBody;
  using PolyhedralEB::CutCellSurface;

  const Box& domainBox = m_domain.domainBox();

  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    const Box box   = m_grids[dit()];
    const Box grown = grow(box, m_numGhost) & domainBox;

    BaseFab<signed char>&       states  = m_cellStates[dit()];
    const BaseFab<signed char>& carried = a_carried[dit()];
    const BaseFab<signed char>& refined = m_refined[dit()];

    IntVectSet& cut = m_cutCells[dit()];

    // A ghost cell no tile carries is classified from its corners alone: regular or covered where the tiles end
    // because the surface does, and possibly cut where they end because the coarser level describes the cells
    // beyond -- such a cell has no surface here, and the face onto it says to ask the coarser level. The cells
    // are a shell on the outside of the tiled region, so the function is evaluated at a node only when a cell
    // that needs it comes by, and each node once.
    Box nodeBox = grown;
    nodeBox.surroundingNodes();

    BaseFab<Real>        nodeValues(nodeBox, 1);
    BaseFab<signed char> nodeKnown(nodeBox, 1);

    nodeKnown.setVal(0);

    for (BoxIterator bit(grown); bit.ok(); ++bit) {
      const IntVect iv = bit();

      if (box.contains(iv)) {
        continue;
      }

      if (carried(iv, 0) == 0) {
        CutCellSurface surface;

        for (int c = 0; c < CutCellSurface::s_numCorners; c++) {
          IntVect node = iv;

          for (int d = 0; d < SpaceDim; d++) {
            node[d] += (c >> d) & 1;
          }

          if (nodeKnown(node, 0) == 0) {
            RealVect x = m_probLo;

            for (int d = 0; d < SpaceDim; d++) {
              x[d] += m_dx * static_cast<Real>(node[d]);
            }

            nodeValues(node, 0) = PolyhedralGeometryShop::snappedValue(a_function, x, m_dx);
            nodeKnown(node, 0)  = 1;
          }

          surface.m_corner[c] = nodeValues(node, 0);
        }

        switch (CutCellBody::classify(surface)) {
        case CutCellBody::Kind::Covered: {
          states(iv, 0) = s_covered;

          break;
        }
        case CutCellBody::Kind::Cut: {
          states(iv, 0) = s_cut;

          break;
        }
        default: {
          states(iv, 0) = s_regular;

          break;
        }
        }
      }
      else if (states(iv, 0) == s_cut && refined(iv, 0) == 0) {
        // The surfaces are kept with ghost cells, so that a box holds its neighbours' cut cells' surfaces as
        // well: each box's set takes in the ghost cells that are cut and that some tile carries, which is exactly
        // the set the owning tiles hold in that region, and an exchange fills them.
        cut |= iv;
      }
    }
  }
}

void
PolyhedralEBGraph::define(const PolyhedralEBGraph& a_source, const DisjointBoxLayout& a_grids, const BaseIF& a_function)
{
  CH_TIME("PolyhedralEBGraph::define(copy)");

  using PolyhedralEB::CutCellSurface;

  if (!a_source.isDefined()) {
    MayDay::Error("PolyhedralEBGraph::define - the source graph is not defined");
  }

  // The layouts are global, so the coverage is checked by cell count: a layout that covers other cells would
  // leave cut cells without surfaces, and outer faces pointing the wrong way.
  long long numCells       = 0;
  long long numSourceCells = 0;

  for (LayoutIterator lit = a_grids.layoutIterator(); lit.ok(); ++lit) {
    numCells += a_grids[lit()].numPts();
  }

  for (LayoutIterator lit = a_source.m_grids.layoutIterator(); lit.ok(); ++lit) {
    numSourceCells += a_source.m_grids[lit()].numPts();
  }

  if (numCells != numSourceCells) {
    MayDay::Error("PolyhedralEBGraph::define - the layout does not cover the cells of the source graph");
  }

  m_domain   = a_source.m_domain;
  m_probLo   = a_source.m_probLo;
  m_dx       = a_source.m_dx;
  m_numGhost = a_source.m_numGhost;
  m_covered  = a_source.m_covered;
  m_grids    = a_grids;

  this->defineData();
  this->markRefined();

  LevelData<BaseFab<signed char>> carried;

  this->markCarried(carried);

  // Valid cells by copy, ghost cells by exchange; the states of the ghost cells no tile carries from the
  // function, and the cut-cell sets from the states, before the surfaces have a container to land in.
  const Copier copier(a_source.m_grids, m_grids, m_domain);

  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    m_cellStates[dit()].setVal(s_regular);
    m_faceStates[dit()].setVal(s_faceClosed);
  }

  a_source.m_cellStates.copyTo(Interval(0, 0), m_cellStates, Interval(0, 0), copier);
  a_source.m_faceStates.copyTo(Interval(0, 2 * SpaceDim - 1), m_faceStates, Interval(0, 2 * SpaceDim - 1), copier);

  m_cellStates.exchange();

  // The refined mask is not copied: it reaches into ghost cells no tile of this level carries, which an
  // exchange cannot fill, and it follows from the covered region, which the copy takes from the source and marks
  // from above. A face state copied as finer stays finer, and linking the copy to its finer level checks the mask.

  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    const Box box = m_grids[dit()];

    const BaseFab<signed char>& states  = m_cellStates[dit()];
    const BaseFab<signed char>& refined = m_refined[dit()];

    IntVectSet& cut = m_cutCells[dit()];

    for (BoxIterator bit(box); bit.ok(); ++bit) {
      if (states(bit(), 0) == s_cut && refined(bit(), 0) == 0) {
        cut |= bit();
      }
    }
  }

  this->defineGhostCells(a_function, carried);

  m_surfaces.define(m_grids, 1, m_numGhost * IntVect::Unit, IVSFABFactory<CutCellSurface>(m_cutCells));

  a_source.m_surfaces.copyTo(Interval(0, 0), m_surfaces, Interval(0, 0), copier);

  m_surfaces.exchange();

  m_isDefined = true;
}

bool
PolyhedralEBGraph::equals(const PolyhedralEBGraph& a_other) const
{
  CH_TIME("PolyhedralEBGraph::equals");

  using PolyhedralEB::CutCellSurface;

  if (!m_isDefined || !a_other.m_isDefined || !(m_grids == a_other.m_grids)) {
    MayDay::Error("PolyhedralEBGraph::equals - both graphs must be defined over the same layout");
  }

  int same = 1;

  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    const BaseFab<signed char>& states      = m_cellStates[dit()];
    const BaseFab<signed char>& otherStates = a_other.m_cellStates[dit()];
    const BaseFab<signed char>& faces       = m_faceStates[dit()];
    const BaseFab<signed char>& otherFaces  = a_other.m_faceStates[dit()];
    const BaseFab<signed char>& refined     = m_refined[dit()];
    const BaseFab<signed char>& otherRef    = a_other.m_refined[dit()];

    for (BoxIterator bit(states.box() & m_domain.domainBox()); bit.ok(); ++bit) {
      same = same && (states(bit(), 0) == otherStates(bit(), 0));
    }

    for (BoxIterator bit(faces.box()); bit.ok(); ++bit) {
      for (int comp = 0; comp < 2 * SpaceDim; comp++) {
        same = same && (faces(bit(), comp) == otherFaces(bit(), comp));
      }
    }

    for (BoxIterator bit(refined.box() & m_domain.domainBox()); bit.ok(); ++bit) {
      same = same && (refined(bit(), 0) == otherRef(bit(), 0));
    }

    const IntVectSet& cut      = m_cutCells[dit()];
    const IntVectSet& otherCut = a_other.m_cutCells[dit()];

    same = same && (cut == otherCut);

    if (same) {
      const IVSFAB<CutCellSurface>& surfaces      = m_surfaces[dit()];
      const IVSFAB<CutCellSurface>& otherSurfaces = a_other.m_surfaces[dit()];

      for (IVSIterator ivsit(cut); ivsit.ok(); ++ivsit) {
        const CutCellSurface& a = surfaces(ivsit(), 0);
        const CutCellSurface& b = otherSurfaces(ivsit(), 0);

        for (int e = 0; e < CutCellSurface::s_numEdges; e++) {
          same = same && (a.m_crossing[e] == b.m_crossing[e]);
        }

        for (int c = 0; c < CutCellSurface::s_numCorners; c++) {
          same = same && (a.m_corner[c] == b.m_corner[c]);
        }
      }
    }
  }

  return ParallelOps::min(same) == 1;
}

void
PolyhedralEBGraph::markCarried(LevelData<BaseFab<signed char>>& a_carried) const
{
  CH_TIME("PolyhedralEBGraph::markCarried");

  // One in the valid cells, zero in the ghost cells, and an exchange fills the ghost cells another tile covers. A
  // ghost cell the exchange leaves at zero is carried by no tile of this level.
  a_carried.define(m_grids, 1, m_numGhost * IntVect::Unit);

  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    a_carried[dit()].setVal(0);
    a_carried[dit()].setVal(1, m_grids[dit()], 0);
  }

  a_carried.exchange();
}

void
PolyhedralEBGraph::markRefined()
{
  CH_TIME("PolyhedralEBGraph::markRefined");

  // The covered region is the same list on every rank, so every box marks its own cells and ghost cells from it
  // with no communication, including ghost cells no tile of this level carries.
  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    BaseFab<signed char>& refined = m_refined[dit()];

    refined.setVal(0);

    const Box region = refined.box() & m_domain.domainBox();

    for (int i = 0; i < m_covered.size(); i++) {
      const Box overlap = m_covered[i] & region;

      if (!overlap.isEmpty()) {
        refined.setVal(1, overlap, 0);
      }
    }
  }
}

void
PolyhedralEBGraph::defineOuterFaces(const LevelData<BaseFab<signed char>>& a_carried)
{
  CH_TIME("PolyhedralEBGraph::defineOuterFaces");

  const Box& domainBox = m_domain.domainBox();

  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    const Box                   box    = m_grids[dit()];
    const BaseFab<signed char>& marker = a_carried[dit()];

    BaseFab<signed char>& faces = m_faceStates[dit()];

    for (BoxIterator bit(box); bit.ok(); ++bit) {
      const IntVect iv = bit();

      for (int dir = 0; dir < SpaceDim; dir++) {
        for (int side = 0; side < 2; side++) {
          const int face = 2 * dir + side;

          if (faces(iv, face) == s_faceClosed) {
            continue;
          }

          const IntVect neighbour = iv + (2 * side - 1) * BASISV(dir);

          if (!domainBox.contains(neighbour)) {
            faces(iv, face) = s_faceBoundary;
          }
          else if (marker(neighbour, 0) == 0) {
            faces(iv, face) = s_faceCoarser;
          }
        }
      }
    }
  }
}

void
PolyhedralEBGraph::link(PolyhedralEBGraph& a_coarse, const PolyhedralEBGraph& a_fine)
{
  CH_TIME("PolyhedralEBGraph::link");

  if (!a_coarse.isDefined() || !a_fine.isDefined()) {
    MayDay::Error("PolyhedralEBGraph::link - both graphs must be defined");
  }
  if (refine(a_coarse.m_domain, 2) != a_fine.m_domain) {
    MayDay::Error("PolyhedralEBGraph::link - the fine domain is not the coarse domain refined by two");
  }

  const Box& domainBox = a_coarse.m_domain.domainBox();

  // A marker that is one on the fine tiles, coarsened, copied onto the coarse layout and its ghost cells. Where
  // it lands, the fine level carries the coarse cell.
  DisjointBoxLayout coarsenedFine;

  coarsen(coarsenedFine, a_fine.m_grids, 2);

  LevelData<BaseFab<signed char>> fineMarker(coarsenedFine, 1, IntVect::Zero);
  LevelData<BaseFab<signed char>> coarMarker(a_coarse.m_grids, 1, 2 * IntVect::Unit);

  for (DataIterator dit(coarsenedFine); dit.ok(); ++dit) {
    fineMarker[dit()].setVal(1);
  }

  for (DataIterator dit(a_coarse.m_grids); dit.ok(); ++dit) {
    coarMarker[dit()].setVal(0);
  }

  const Copier copier(coarsenedFine, a_coarse.m_grids, a_coarse.m_domain, 2 * IntVect::Unit);

  fineMarker.copyTo(Interval(0, 0), coarMarker, Interval(0, 0), copier);

  for (DataIterator dit(a_coarse.m_grids); dit.ok(); ++dit) {
    const Box                   box    = a_coarse.m_grids[dit()];
    const BaseFab<signed char>& marker = coarMarker[dit()];

    const BaseFab<signed char>& refined = a_coarse.m_refined[dit()];
    BaseFab<signed char>&       faces   = a_coarse.m_faceStates[dit()];

    // The coarse level was told its covered region when it was defined, from the same tiles the fine level is
    // built over, so the two must agree cell for cell -- over the box and the two ghost cells a neighbour's body on
    // a level boundary reads. A disagreement means the levels were not built from one set of tiles.
    for (BoxIterator bit(grow(box, 2) & domainBox); bit.ok(); ++bit) {
      if ((refined(bit(), 0) != 0) != (marker(bit(), 0) != 0)) {
        pout() << "PolyhedralEBGraph::link - cell " << bit() << " is " << ((refined(bit(), 0) != 0) ? "" : "not ")
               << "in the coarse level's covered region but is " << ((marker(bit(), 0) != 0) ? "" : "not ")
               << "under the fine level's tiles" << endl;

        MayDay::Error("PolyhedralEBGraph::link - the coarse level's covered region disagrees with the fine tiles");
      }
    }

    for (BoxIterator bit(box); bit.ok(); ++bit) {
      const IntVect iv = bit();

      if (marker(iv, 0) != 0) {
        continue;
      }

      // an open face of an unrefined cell onto a refined neighbour is described by the fine level
      for (int dir = 0; dir < SpaceDim; dir++) {
        for (int side = 0; side < 2; side++) {
          const int face = 2 * dir + side;

          if (faces(iv, face) != s_faceSameLevel) {
            continue;
          }

          const IntVect neighbour = iv + (2 * side - 1) * BASISV(dir);

          if (domainBox.contains(neighbour) && marker(neighbour, 0) != 0) {
            faces(iv, face) = s_faceFiner;
          }
        }
      }
    }
  }
}

bool
PolyhedralEBGraph::isDefined() const noexcept
{
  return m_isDefined;
}

const ProblemDomain&
PolyhedralEBGraph::getDomain() const noexcept
{
  return m_domain;
}

Real
PolyhedralEBGraph::getDx() const noexcept
{
  return m_dx;
}

const DisjointBoxLayout&
PolyhedralEBGraph::getGrids() const noexcept
{
  return m_grids;
}

const LayoutData<IntVectSet>&
PolyhedralEBGraph::getCutCells() const noexcept
{
  return m_cutCells;
}

const LevelData<BaseFab<signed char>>&
PolyhedralEBGraph::getCellStates() const noexcept
{
  return m_cellStates;
}

const LevelData<BaseFab<signed char>>&
PolyhedralEBGraph::getFaceStates() const noexcept
{
  return m_faceStates;
}

const LevelData<BaseFab<signed char>>&
PolyhedralEBGraph::getRefinedMask() const noexcept
{
  return m_refined;
}

const LevelData<IVSFAB<PolyhedralEB::CutCellSurface>>&
PolyhedralEBGraph::getSurfaces() const noexcept
{
  return m_surfaces;
}

#include <CD_NamespaceFooter.H>
