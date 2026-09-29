/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_CutCellFaceOverrides.cpp
 * @brief  Implementation of CD_CutCellFaceOverrides.H
 * @author Robert Marskar
 */

// Std includes
#include <algorithm>

// Chombo includes
#include <CH_assert.H>

// Our includes
#include <CD_CutCellFaceOverrides.H>
#include <CD_NamespaceHeader.H>

namespace PolyhedralEB {

CutCellFaceOverrides::CutCellFaceOverrides() noexcept
{
  this->clear();
}

void
CutCellFaceOverrides::clear() noexcept
{
  m_cells.clear();
  m_face.clear();
  m_reason.clear();
  m_vertices.clear();
  m_vertexEdges.clear();

  m_cellFaces.assign(1, 0);
  m_facePolygons.assign(1, 0);
  m_polygonVertices.assign(1, 0);
}

void
CutCellFaceOverrides::beginCell(const IntVect& a_cell)
{
  CH_assert(m_cells.empty() || precedes(m_cells.back(), a_cell));

  m_cells.push_back(a_cell);
  m_cellFaces.push_back(m_cellFaces.back());
}

void
CutCellFaceOverrides::beginFace(const int a_face, const int a_reason)
{
  CH_assert(!m_cells.empty());
  CH_assert(a_face >= 0 && a_face < 2 * SpaceDim);
  CH_assert(a_reason == s_finer || a_reason == s_closed);

  m_face.push_back(a_face);
  m_reason.push_back(a_reason);
  m_facePolygons.push_back(m_facePolygons.back());

  m_cellFaces.back()++;
}

void
CutCellFaceOverrides::addPolygon(const RealVect* a_vertices, const int* a_vertexEdges, const int a_numVertices)
{
  CH_assert(!m_face.empty());
  CH_assert(m_reason.back() == s_finer);
  CH_assert(a_numVertices >= 3);

  for (int i = 0; i < a_numVertices; i++) {
    m_vertices.push_back(a_vertices[i]);
    m_vertexEdges.push_back(a_vertexEdges[i]);
  }

  m_polygonVertices.push_back(static_cast<int>(m_vertices.size()));

  m_facePolygons.back()++;
}

void
CutCellFaceOverrides::endCell() noexcept
{
  CH_assert(!m_cells.empty());

  const std::size_t n = m_cells.size();

  if (m_cellFaces[n] == m_cellFaces[n - 1]) {
    m_cells.pop_back();
    m_cellFaces.pop_back();
  }
}

int
CutCellFaceOverrides::numCells() const noexcept
{
  return static_cast<int>(m_cells.size());
}

int
CutCellFaceOverrides::find(const IntVect& a_cell) const noexcept
{
  const auto it = std::lower_bound(m_cells.begin(), m_cells.end(), a_cell, precedes);

  if (it == m_cells.end() || *it != a_cell) {
    return -1;
  }

  return static_cast<int>(it - m_cells.begin());
}

const IntVect&
CutCellFaceOverrides::cell(const int a_index) const noexcept
{
  CH_assert(a_index >= 0 && a_index < this->numCells());

  return m_cells[a_index];
}

void
CutCellFaceOverrides::faces(const int a_index, int& a_begin, int& a_end) const noexcept
{
  CH_assert(a_index >= 0 && a_index < this->numCells());

  a_begin = m_cellFaces[a_index];
  a_end   = m_cellFaces[a_index + 1];
}

int
CutCellFaceOverrides::face(const int a_face) const noexcept
{
  return m_face[a_face];
}

int
CutCellFaceOverrides::reason(const int a_face) const noexcept
{
  return m_reason[a_face];
}

void
CutCellFaceOverrides::polygons(const int a_face, int& a_begin, int& a_end) const noexcept
{
  a_begin = m_facePolygons[a_face];
  a_end   = m_facePolygons[a_face + 1];
}

void
CutCellFaceOverrides::vertices(const int a_polygon, int& a_begin, int& a_end) const noexcept
{
  a_begin = m_polygonVertices[a_polygon];
  a_end   = m_polygonVertices[a_polygon + 1];
}

const RealVect&
CutCellFaceOverrides::vertex(const int a_vertex) const noexcept
{
  return m_vertices[a_vertex];
}

int
CutCellFaceOverrides::vertexEdge(const int a_vertex) const noexcept
{
  return m_vertexEdges[a_vertex];
}

bool
CutCellFaceOverrides::equals(const CutCellFaceOverrides& a_other) const noexcept
{
  if (m_cells != a_other.m_cells || m_cellFaces != a_other.m_cellFaces || m_face != a_other.m_face ||
      m_reason != a_other.m_reason || m_facePolygons != a_other.m_facePolygons ||
      m_polygonVertices != a_other.m_polygonVertices || m_vertexEdges != a_other.m_vertexEdges) {
    return false;
  }

  for (std::size_t i = 0; i < m_vertices.size(); i++) {
    for (int d = 0; d < SpaceDim; d++) {
      if (m_vertices[i][d] != a_other.m_vertices[i][d]) {
        return false;
      }
    }
  }

  return true;
}

bool
CutCellFaceOverrides::precedes(const IntVect& a_first, const IntVect& a_second) noexcept
{
  for (int d = SpaceDim - 1; d >= 0; d--) {
    if (a_first[d] != a_second[d]) {
      return a_first[d] < a_second[d];
    }
  }

  return false;
}

} // namespace PolyhedralEB

#include <CD_NamespaceFooter.H>
