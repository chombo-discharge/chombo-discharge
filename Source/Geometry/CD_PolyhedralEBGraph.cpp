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
                          const int            a_numGhost)
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

  this->defineGrids(a_cutTiles);

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
  m_refined.define(m_grids, 1, 2 * IntVect::Unit);

  // Each box's cut-cell set is a bitmap over the box and its ghost ring, which is where its cells come from: a
  // set that starts empty would be a tree, several kilobytes for a few hundred scattered cells, and the surface
  // container keeps a copy of it.
  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    const Box grown = grow(m_grids[dit()], m_numGhost) & m_domain.domainBox();

    m_cutCells[dit()] = IntVectSet(DenseIntVectSet(grown, false));
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

    m_refined[dit()].setVal(0);

    // Node values once per node and each crossed edge bisected once, shared by the cells of the box, as the
    // generator does when it builds a box.
    BaseFab<Real> nodeValues;
    BaseFab<Real> intercept[SpaceDim];

    PolyhedralGeometryShop::fillNodeValues(a_function, nodeValues, box, m_probLo, m_dx);
    PolyhedralGeometryShop::defineIntercepts(intercept, box);

    for (BoxIterator bit(box); bit.ok(); ++bit) {
      const IntVect iv = bit();

      CutCellSurface surface;

      PolyhedralGeometryShop::buildSurface(a_function, intercept, surface, nodeValues, iv, m_probLo, m_dx);

      const CutCellBody::Kind kind = CutCellBody::classify(surface);

      int state = s_regular;

      CutCellBody body;

      if (kind == CutCellBody::Kind::Covered) {
        state = s_covered;
      }
      else if (kind == CutCellBody::Kind::Cut) {
        if (!body.define(surface)) {
          pout() << "PolyhedralEBGraph::defineCells - cell " << iv << " did not close" << endl;

          MayDay::Error("PolyhedralEBGraph::defineCells - a cut cell's body did not close");
        }

        // A body with nothing in it is a regular cell. Nothing else is discarded: a cell's interface is what
        // closes the surface against its neighbours, so a cell dropped for holding little would leave a hole
        // exactly the size of what it held, and no repair on the neighbours' faces can put it back. The volume
        // threshold that keeps such cells out of the index space is applied where the index space is built,
        // which is where its reason -- a solver that would rather not see a cell of no volume -- applies.
        if (PolyhedralGeometryShop::isDust(body)) {
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

          if (state == s_regular) {
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
      else if (states(iv, 0) == s_cut) {
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
  m_grids    = a_grids;

  this->defineData();

  LevelData<BaseFab<signed char>> carried;

  this->markCarried(carried);

  // Valid cells by copy, ghost cells by exchange; the states of the ghost cells no tile carries from the
  // function, and the cut-cell sets from the states, before the surfaces have a container to land in.
  const Copier copier(a_source.m_grids, m_grids, m_domain);

  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    m_cellStates[dit()].setVal(s_regular);
    m_faceStates[dit()].setVal(s_faceClosed);
    m_refined[dit()].setVal(0);
  }

  a_source.m_cellStates.copyTo(Interval(0, 0), m_cellStates, Interval(0, 0), copier);
  a_source.m_faceStates.copyTo(Interval(0, 2 * SpaceDim - 1), m_faceStates, Interval(0, 2 * SpaceDim - 1), copier);

  m_cellStates.exchange();

  // The refined mask is not copied: it reaches into ghost cells no tile of this level carries, which an
  // exchange cannot fill, and it is derived from the finer level in any case. The copy is linked to its finer
  // level afterwards, as the original was, which sets the mask and the faces the finer level describes; a face
  // state copied as finer stays finer, since link only ever marks faces.

  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    const Box box = m_grids[dit()];

    const BaseFab<signed char>& states = m_cellStates[dit()];

    IntVectSet& cut = m_cutCells[dit()];

    for (BoxIterator bit(box); bit.ok(); ++bit) {
      if (states(bit(), 0) == s_cut) {
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

    BaseFab<signed char>& refined = a_coarse.m_refined[dit()];
    BaseFab<signed char>& faces   = a_coarse.m_faceStates[dit()];

    // the mask keeps two ghost cells, so a box knows whether its neighbours' cells are refined, and their
    // neighbours' in turn, which is what rebuilding a neighbour's body on a level boundary asks
    refined.copy(marker, grow(box, 2) & domainBox);

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
