/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_AMRPolyhedralEBGraph.cpp
 * @brief  Implementation of CD_AMRPolyhedralEBGraph.H
 * @author Robert Marskar
 */

// Std includes
#include <algorithm>
#include <numeric>

// Chombo includes
#include <BoxIterator.H>
#include <CH_assert.H>
#include <CH_Timer.H>

// Our includes
#include <CD_AMRPolyhedralEBGraph.H>
#include <CD_PolyhedralEBGraph.H>
#include <CD_NamespaceHeader.H>

namespace {

/**
 * @brief A LayoutData over a layout, for every pair of destination level and geometry level.
 * @param[out] a_data            One per pair.
 * @param[in]  a_grids           Destination layout of every level.
 * @param[in]  a_numSourceLevels Number of geometry levels.
 */
template <class T>
void
allocate(std::vector<std::vector<RefCountedPtr<LayoutData<T>>>>& a_data,
         const std::vector<DisjointBoxLayout>&                   a_grids,
         const int                                               a_numSourceLevels)
{
  a_data.assign(a_grids.size(), {});

  for (std::size_t i = 0; i < a_grids.size(); i++) {
    a_data[i].resize(a_numSourceLevels);

    for (int lvl = 0; lvl < a_numSourceLevels; lvl++) {
      a_data[i][lvl] = RefCountedPtr<LayoutData<T>>(new LayoutData<T>(a_grids[i]));
    }
  }
}

/**
 * @brief Position of a cell in a box, counted as a BoxIterator meets the box's cells.
 * @param[in] a_box  The box.
 * @param[in] a_cell A cell of it.
 * @return The position.
 */
long
linearIndex(const Box& a_box, const IntVect& a_cell)
{
  long position = 0;
  long stride   = 1;

  for (int d = 0; d < SpaceDim; d++) {
    position += stride * (a_cell[d] - a_box.smallEnd(d));
    stride *= a_box.size(d);
  }

  return position;
}

} // namespace

AMRPolyhedralEBGraph::AMRPolyhedralEBGraph() noexcept
{
  m_isDefined       = false;
  m_ghost           = 0;
  m_numSourceLevels = 0;
}

AMRPolyhedralEBGraph::~AMRPolyhedralEBGraph() noexcept
{}

void
AMRPolyhedralEBGraph::define(const std::vector<DisjointBoxLayout>& a_grids,
                             const std::vector<ProblemDomain>&     a_domains,
                             const std::vector<int>&               a_geometryLevels,
                             const int                             a_numSourceLevels,
                             const int                             a_ghost)
{
  CH_TIME("AMRPolyhedralEBGraph::define");

  CH_assert(a_grids.size() == a_domains.size());
  CH_assert(a_grids.size() == a_geometryLevels.size());

  m_grids           = a_grids;
  m_domains         = a_domains;
  m_geometryLevels  = a_geometryLevels;
  m_numSourceLevels = a_numSourceLevels;
  m_ghost           = a_ghost;

  allocate(m_pieces, a_grids, a_numSourceLevels);
  allocate(m_pieceOffsets, a_grids, a_numSourceLevels);
  allocate(m_states, a_grids, a_numSourceLevels);
  allocate(m_surfaceCells, a_grids, a_numSourceLevels);
  allocate(m_surfaces, a_grids, a_numSourceLevels);
  allocate(m_faceOverrides, a_grids, a_numSourceLevels);
  allocate(m_binSize, a_grids, a_numSourceLevels);
  allocate(m_binBox, a_grids, a_numSourceLevels);
  allocate(m_binOffsets, a_grids, a_numSourceLevels);
  allocate(m_binPieces, a_grids, a_numSourceLevels);

  for (std::size_t i = 0; i < a_grids.size(); i++) {
    for (int lvl = 0; lvl < a_numSourceLevels; lvl++) {
      for (DataIterator dit(a_grids[i]); dit.ok(); ++dit) {
        (*m_pieceOffsets[i][lvl])[dit()].assign(1, 0);
      }
    }
  }

  m_isDefined = true;
}

bool
AMRPolyhedralEBGraph::isDefined() const noexcept
{
  return m_isDefined;
}

int
AMRPolyhedralEBGraph::getNumLevels() const noexcept
{
  return static_cast<int>(m_grids.size());
}

int
AMRPolyhedralEBGraph::getNumSourceLevels() const noexcept
{
  return m_numSourceLevels;
}

int
AMRPolyhedralEBGraph::getGhost() const noexcept
{
  return m_ghost;
}

const DisjointBoxLayout&
AMRPolyhedralEBGraph::getGrids(const int a_level) const noexcept
{
  return m_grids[a_level];
}

const ProblemDomain&
AMRPolyhedralEBGraph::getDomain(const int a_level) const noexcept
{
  return m_domains[a_level];
}

int
AMRPolyhedralEBGraph::getGeometryLevel(const int a_level) const noexcept
{
  return m_geometryLevels[a_level];
}

void
AMRPolyhedralEBGraph::getCutCells(std::vector<PolyhedralEB::CutCellDescription>& a_cells,
                                  const int                                      a_level,
                                  const DataIndex&                               a_dit,
                                  const IntVect&                                 a_cell) const
{
  CH_TIME("AMRPolyhedralEBGraph::getCutCells");

  a_cells.clear();

  const int geometryLevel = m_geometryLevels[a_level];

  for (int lvl = 0; lvl < m_numSourceLevels; lvl++) {
    const std::vector<Box>& pieces = this->getPieces(a_level, a_dit, lvl);

    if (pieces.empty()) {
      continue;
    }

    if (lvl == geometryLevel) {
      if (this->getSourceState(a_level, a_dit, lvl, a_cell) == PolyhedralEBGraph::s_cut) {
        this->describe(a_cells, a_level, a_dit, lvl, a_cell);
      }
    }
    else if (lvl < geometryLevel) {
      const IntVect ancestor = coarsen(a_cell, 1 << (geometryLevel - lvl));

      if (this->getSourceState(a_level, a_dit, lvl, ancestor) == PolyhedralEBGraph::s_cut) {
        this->describe(a_cells, a_level, a_dit, lvl, ancestor);
      }
    }
    else {
      // Every cut cell of the finer level under this one. A cut cell holds a surface, and the cells that hold one are
      // kept in order, so only the stretch of them between the corners of the cells under this one is looked at.
      const Box children = refine(Box(a_cell, a_cell), 1 << (lvl - geometryLevel));

      const std::vector<IntVect>& withSurface = this->getSurfaceCells(a_level, a_dit, lvl);

      const auto first = std::lower_bound(withSurface.begin(),
                                          withSurface.end(),
                                          children.smallEnd(),
                                          PolyhedralEB::CutCellFaceOverrides::precedes);

      for (auto it = first; it != withSurface.end(); ++it) {
        if (PolyhedralEB::CutCellFaceOverrides::precedes(children.bigEnd(), *it)) {
          break;
        }

        if (children.contains(*it) && this->getSourceState(a_level, a_dit, lvl, *it) == PolyhedralEBGraph::s_cut) {
          this->describe(a_cells, a_level, a_dit, lvl, *it);
        }
      }
    }
  }
}

int
AMRPolyhedralEBGraph::getState(const int a_level, const DataIndex& a_dit, const IntVect& a_cell) const
{
  const int geometryLevel = m_geometryLevels[a_level];

  bool anyRegular = false;
  bool anyCovered = false;
  bool anyCut     = false;

  const auto take = [&](const int a_state) {
    anyRegular = anyRegular || a_state == PolyhedralEBGraph::s_regular;
    anyCovered = anyCovered || a_state == PolyhedralEBGraph::s_covered;
    anyCut     = anyCut || a_state == PolyhedralEBGraph::s_cut;
  };

  for (int lvl = 0; lvl < m_numSourceLevels; lvl++) {
    const std::vector<Box>& pieces = this->getPieces(a_level, a_dit, lvl);

    if (pieces.empty()) {
      continue;
    }

    if (lvl == geometryLevel) {
      take(this->getSourceState(a_level, a_dit, lvl, a_cell));
    }
    else if (lvl < geometryLevel) {
      take(this->getSourceState(a_level, a_dit, lvl, coarsen(a_cell, 1 << (geometryLevel - lvl))));
    }
    else {
      const Box children = refine(Box(a_cell, a_cell), 1 << (lvl - geometryLevel));

      std::vector<int> near;

      this->piecesNear(near, a_level, a_dit, lvl, children);

      for (const int p : near) {
        const Box overlap = pieces[p] & children;

        if (overlap.isEmpty()) {
          continue;
        }

        for (BoxIterator bit(overlap); bit.ok(); ++bit) {
          take(this->getSourceState(a_level, a_dit, lvl, bit()));
        }
      }
    }
  }

  if (anyCut || (anyRegular && anyCovered)) {
    return PolyhedralEBGraph::s_cut;
  }

  if (anyRegular) {
    return PolyhedralEBGraph::s_regular;
  }

  if (anyCovered) {
    return PolyhedralEBGraph::s_covered;
  }

  return s_unowned;
}

const std::vector<Box>&
AMRPolyhedralEBGraph::getPieces(const int a_level, const DataIndex& a_dit, const int a_sourceLevel) const noexcept
{
  return (*m_pieces[a_level][a_sourceLevel])[a_dit];
}

int
AMRPolyhedralEBGraph::getSourceState(const int        a_level,
                                     const DataIndex& a_dit,
                                     const int        a_sourceLevel,
                                     const IntVect&   a_sourceCell) const
{
  const std::vector<Box>&         pieces  = (*m_pieces[a_level][a_sourceLevel])[a_dit];
  const std::vector<int>&         offsets = (*m_pieceOffsets[a_level][a_sourceLevel])[a_dit];
  const std::vector<signed char>& states  = (*m_states[a_level][a_sourceLevel])[a_dit];

  if (pieces.empty()) {
    return s_unowned;
  }

  const int     binSize    = (*m_binSize[a_level][a_sourceLevel])[a_dit];
  const Box&    binBox     = (*m_binBox[a_level][a_sourceLevel])[a_dit];
  const auto&   binOffsets = (*m_binOffsets[a_level][a_sourceLevel])[a_dit];
  const auto&   binPieces  = (*m_binPieces[a_level][a_sourceLevel])[a_dit];
  const IntVect bin        = coarsen(a_sourceCell, binSize);

  CH_assert(!binOffsets.empty());

  if (!binBox.contains(bin)) {
    return s_unowned;
  }

  const long b = linearIndex(binBox, bin);

  for (int k = binOffsets[b]; k < binOffsets[b + 1]; k++) {
    const int p = binPieces[k];

    if (pieces[p].contains(a_sourceCell)) {
      return states[offsets[p] + linearIndex(pieces[p], a_sourceCell)];
    }
  }

  return s_unowned;
}

const std::vector<IntVect>&
AMRPolyhedralEBGraph::getSurfaceCells(const int a_level, const DataIndex& a_dit, const int a_sourceLevel) const noexcept
{
  return (*m_surfaceCells[a_level][a_sourceLevel])[a_dit];
}

const PolyhedralEB::CutCellSurface*
AMRPolyhedralEBGraph::getSourceSurface(const int        a_level,
                                       const DataIndex& a_dit,
                                       const int        a_sourceLevel,
                                       const IntVect&   a_sourceCell) const
{
  using PolyhedralEB::CutCellFaceOverrides;

  const std::vector<IntVect>& cells = (*m_surfaceCells[a_level][a_sourceLevel])[a_dit];

  const auto it = std::lower_bound(cells.begin(), cells.end(), a_sourceCell, CutCellFaceOverrides::precedes);

  if (it == cells.end() || *it != a_sourceCell) {
    return nullptr;
  }

  return &(*m_surfaces[a_level][a_sourceLevel])[a_dit][it - cells.begin()];
}

const PolyhedralEB::CutCellFaceOverrides&
AMRPolyhedralEBGraph::getFaceOverrides(const int        a_level,
                                       const DataIndex& a_dit,
                                       const int        a_sourceLevel) const noexcept
{
  return (*m_faceOverrides[a_level][a_sourceLevel])[a_dit];
}

PolyhedralEB::CutCellFaceOverrides&
AMRPolyhedralEBGraph::getFaceOverrides(const int a_level, const DataIndex& a_dit, const int a_sourceLevel) noexcept
{
  return (*m_faceOverrides[a_level][a_sourceLevel])[a_dit];
}

void
AMRPolyhedralEBGraph::addPiece(const int          a_level,
                               const DataIndex&   a_dit,
                               const int          a_sourceLevel,
                               const Box&         a_region,
                               const signed char* a_states)
{
  std::vector<Box>&         pieces  = (*m_pieces[a_level][a_sourceLevel])[a_dit];
  std::vector<int>&         offsets = (*m_pieceOffsets[a_level][a_sourceLevel])[a_dit];
  std::vector<signed char>& states  = (*m_states[a_level][a_sourceLevel])[a_dit];

  pieces.push_back(a_region);
  states.insert(states.end(), a_states, a_states + a_region.numPts());
  offsets.push_back(static_cast<int>(states.size()));
}

void
AMRPolyhedralEBGraph::addSurface(const int                           a_level,
                                 const DataIndex&                    a_dit,
                                 const int                           a_sourceLevel,
                                 const IntVect&                      a_sourceCell,
                                 const PolyhedralEB::CutCellSurface& a_surface)
{
  (*m_surfaceCells[a_level][a_sourceLevel])[a_dit].push_back(a_sourceCell);
  (*m_surfaces[a_level][a_sourceLevel])[a_dit].push_back(a_surface);
}

void
AMRPolyhedralEBGraph::piecesNear(std::vector<int>& a_pieces,
                                 const int         a_level,
                                 const DataIndex&  a_dit,
                                 const int         a_sourceLevel,
                                 const Box&        a_region) const
{
  a_pieces.clear();

  const std::vector<Box>& pieces = (*m_pieces[a_level][a_sourceLevel])[a_dit];

  if (pieces.empty()) {
    return;
  }

  const int   binSize    = (*m_binSize[a_level][a_sourceLevel])[a_dit];
  const Box&  binBox     = (*m_binBox[a_level][a_sourceLevel])[a_dit];
  const auto& binOffsets = (*m_binOffsets[a_level][a_sourceLevel])[a_dit];
  const auto& binPieces  = (*m_binPieces[a_level][a_sourceLevel])[a_dit];

  const Box bins = coarsen(a_region, binSize) & binBox;

  for (BoxIterator bit(bins); bit.ok(); ++bit) {
    const long b = linearIndex(binBox, bit());

    for (int k = binOffsets[b]; k < binOffsets[b + 1]; k++) {
      a_pieces.push_back(binPieces[k]);
    }
  }

  std::sort(a_pieces.begin(), a_pieces.end());

  a_pieces.erase(std::unique(a_pieces.begin(), a_pieces.end()), a_pieces.end());
}

void
AMRPolyhedralEBGraph::finish()
{
  CH_TIME("AMRPolyhedralEBGraph::finish");

  using PolyhedralEB::CutCellFaceOverrides;
  using PolyhedralEB::CutCellSurface;

  for (std::size_t i = 0; i < m_grids.size(); i++) {
    for (int lvl = 0; lvl < m_numSourceLevels; lvl++) {
      for (DataIterator dit(m_grids[i]); dit.ok(); ++dit) {
        std::vector<IntVect>&        cells    = (*m_surfaceCells[i][lvl])[dit()];
        std::vector<CutCellSurface>& surfaces = (*m_surfaces[i][lvl])[dit()];

        std::vector<std::size_t> order(cells.size());

        std::iota(order.begin(), order.end(), 0);

        std::sort(order.begin(), order.end(), [&cells](const std::size_t a_first, const std::size_t a_second) -> bool {
          return CutCellFaceOverrides::precedes(cells[a_first], cells[a_second]);
        });

        std::vector<IntVect>        sortedCells(cells.size());
        std::vector<CutCellSurface> sortedSurfaces(surfaces.size());

        for (std::size_t n = 0; n < order.size(); n++) {
          sortedCells[n]    = cells[order[n]];
          sortedSurfaces[n] = surfaces[order[n]];
        }

        cells.swap(sortedCells);
        surfaces.swap(sortedSurfaces);

        // the pieces' index: bins as wide as the widest piece, so a piece reaches at most two bins a direction
        const std::vector<Box>& pieces = (*m_pieces[i][lvl])[dit()];

        int& binSize = (*m_binSize[i][lvl])[dit()];
        Box& binBox  = (*m_binBox[i][lvl])[dit()];

        std::vector<int>& binOffsets = (*m_binOffsets[i][lvl])[dit()];
        std::vector<int>& binPieces  = (*m_binPieces[i][lvl])[dit()];

        binSize = 1;
        binBox  = Box();

        binOffsets.clear();
        binPieces.clear();

        if (pieces.empty()) {
          continue;
        }

        Box hull = pieces[0];

        for (const Box& piece : pieces) {
          hull = minBox(hull, piece);

          for (int d = 0; d < SpaceDim; d++) {
            binSize = std::max(binSize, piece.size(d));
          }
        }

        binBox = coarsen(hull, binSize);

        std::vector<int> counts(binBox.numPts(), 0);

        for (const Box& piece : pieces) {
          for (BoxIterator bit(coarsen(piece, binSize)); bit.ok(); ++bit) {
            counts[linearIndex(binBox, bit())]++;
          }
        }

        binOffsets.assign(counts.size() + 1, 0);

        for (std::size_t b = 0; b < counts.size(); b++) {
          binOffsets[b + 1] = binOffsets[b] + counts[b];
        }

        binPieces.assign(binOffsets.back(), -1);

        std::vector<int> fill(binOffsets.begin(), binOffsets.end() - 1);

        for (std::size_t p = 0; p < pieces.size(); p++) {
          for (BoxIterator bit(coarsen(pieces[p], binSize)); bit.ok(); ++bit) {
            binPieces[fill[linearIndex(binBox, bit())]++] = static_cast<int>(p);
          }
        }
      }
    }
  }
}

void
AMRPolyhedralEBGraph::describe(std::vector<PolyhedralEB::CutCellDescription>& a_cells,
                               const int                                      a_level,
                               const DataIndex&                               a_dit,
                               const int                                      a_sourceLevel,
                               const IntVect&                                 a_sourceCell) const
{
  using PolyhedralEB::CutCellDescription;
  using PolyhedralEB::CutCellFaceOverrides;

  const PolyhedralEB::CutCellSurface* surface = this->getSourceSurface(a_level, a_dit, a_sourceLevel, a_sourceCell);

  CH_assert(surface != nullptr);

  CutCellDescription description;

  description.m_level   = a_sourceLevel;
  description.m_cell    = a_sourceCell;
  description.m_state   = this->getSourceState(a_level, a_dit, a_sourceLevel, a_sourceCell);
  description.m_surface = *surface;

  // the cell's own entry, copied out of the box's overrides through their linear form
  const CutCellFaceOverrides& overrides = this->getFaceOverrides(a_level, a_dit, a_sourceLevel);

  const int entry = overrides.find(a_sourceCell);

  if (entry >= 0) {
    std::vector<char> buffer(overrides.linearSize(entry));

    overrides.linearOut(buffer.data(), entry);

    description.m_overrides.linearIn(buffer.data(), a_sourceCell);
  }

  a_cells.push_back(description);
}

#include <CD_NamespaceFooter.H>
