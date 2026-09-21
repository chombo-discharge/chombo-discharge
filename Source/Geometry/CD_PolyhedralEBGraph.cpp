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
                          const Real           a_volumeThreshold)
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
  this->defineCells(a_function, a_volumeThreshold);
  this->defineOuterFaces();

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

  m_cutCells.define(m_grids);
  m_cellStates.define(m_grids, 1, m_numGhost * IntVect::Unit);
  m_faceStates.define(m_grids, 2 * SpaceDim, IntVect::Zero);
  m_refined.define(m_grids, 1, IntVect::Unit);
}

void
PolyhedralEBGraph::defineCells(const BaseIF& a_function, const Real a_volumeThreshold)
{
  CH_TIME("PolyhedralEBGraph::defineCells");

  using PolyhedralEB::CutCellBody;
  using PolyhedralEB::CutCellSurface;

  const Box& domainBox = m_domain.domainBox();

  // The surfaces of a box's cut cells are kept in the order a BoxIterator meets the cells, and moved into the
  // level container once every box's cut-cell set is known; the container is defined over those sets.
  LayoutData<Vector<CutCellSurface>> kept(m_grids);

  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    const Box box   = m_grids[dit()];
    const Box grown = grow(box, m_numGhost) & domainBox;

    BaseFab<int>& states = m_cellStates[dit()];
    BaseFab<int>& faces  = m_faceStates[dit()];
    IntVectSet&   cut    = m_cutCells[dit()];

    Vector<CutCellSurface>& surfaces = kept[dit()];

    // Ghost cells outside the domain have no cell to describe; they are regular so that nothing reads them as a
    // boundary of the fluid.
    states.setVal(s_regular);
    faces.setVal(s_faceClosed);

    m_refined[dit()].setVal(0);

    // Node values once per node and each crossed edge bisected once, shared by the cells of the region, as the
    // generator does when it builds a box.
    BaseFab<Real> nodeValues;
    BaseFab<Real> intercept[SpaceDim];

    PolyhedralGeometryShop::fillNodeValues(a_function, nodeValues, grown, m_probLo, m_dx);
    PolyhedralGeometryShop::defineIntercepts(intercept, grown);

    // A ghost cell is classified from its corners alone, with no crossing bisected and no body built: a cell
    // another tile carries is overwritten with that tile's exact state by the exchange below, and a cell no tile
    // carries is regular or covered by construction, which the corners decide. The valid cells get the full
    // reconstruction, since their bodies decide the thresholds and the apertures.
    for (BoxIterator bit(grown); bit.ok(); ++bit) {
      const IntVect iv = bit();

      CutCellSurface surface;

      if (!box.contains(iv)) {
        PolyhedralGeometryShop::fillCorners(surface, nodeValues, iv);

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

        continue;
      }

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

        // the generator's rules: too little fluid is a covered cell, too little solid a regular one
        if (a_volumeThreshold > 0.0 && body.volumeFraction() < a_volumeThreshold) {
          state = s_covered;
        }
        else if (PolyhedralGeometryShop::isDust(body, a_volumeThreshold)) {
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

  // ghost cells another tile carries take that tile's state
  m_cellStates.exchange();

  m_surfaces.define(m_grids, 1, IntVect::Zero, IVSFABFactory<CutCellSurface>(m_cutCells));

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
}

void
PolyhedralEBGraph::defineOuterFaces()
{
  CH_TIME("PolyhedralEBGraph::defineOuterFaces");

  const Box& domainBox = m_domain.domainBox();

  // A marker that is one on every cell a tile of this level carries and zero elsewhere: one in the valid cells,
  // zero in the ghost cells, and an exchange fills the ghost cells another tile covers. A ghost cell the exchange
  // leaves at zero is carried by no tile of this level.
  LevelData<BaseFab<int>> carried(m_grids, 1, IntVect::Unit);

  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    carried[dit()].setVal(0);
    carried[dit()].setVal(1, m_grids[dit()], 0);
  }

  carried.exchange();

  for (DataIterator dit(m_grids); dit.ok(); ++dit) {
    const Box           box    = m_grids[dit()];
    const BaseFab<int>& marker = carried[dit()];

    BaseFab<int>& faces = m_faceStates[dit()];

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

  LevelData<BaseFab<int>> fineMarker(coarsenedFine, 1, IntVect::Zero);
  LevelData<BaseFab<int>> coarMarker(a_coarse.m_grids, 1, IntVect::Unit);

  for (DataIterator dit(coarsenedFine); dit.ok(); ++dit) {
    fineMarker[dit()].setVal(1);
  }

  for (DataIterator dit(a_coarse.m_grids); dit.ok(); ++dit) {
    coarMarker[dit()].setVal(0);
  }

  const Copier copier(coarsenedFine, a_coarse.m_grids, a_coarse.m_domain, IntVect::Unit);

  fineMarker.copyTo(Interval(0, 0), coarMarker, Interval(0, 0), copier);

  for (DataIterator dit(a_coarse.m_grids); dit.ok(); ++dit) {
    const Box           box    = a_coarse.m_grids[dit()];
    const BaseFab<int>& marker = coarMarker[dit()];

    BaseFab<int>& refined = a_coarse.m_refined[dit()];
    BaseFab<int>& faces   = a_coarse.m_faceStates[dit()];

    // the mask keeps one ghost cell, so a box knows whether its neighbours' cells are refined as well
    refined.copy(marker, grow(box, 1) & domainBox);

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

const LevelData<BaseFab<int>>&
PolyhedralEBGraph::getCellStates() const noexcept
{
  return m_cellStates;
}

const LevelData<BaseFab<int>>&
PolyhedralEBGraph::getFaceStates() const noexcept
{
  return m_faceStates;
}

const LevelData<BaseFab<int>>&
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
