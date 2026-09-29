/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_CutCellFaceOverrides.H
 * @brief  Declaration of the faces of a box's cut cells that their own surfaces do not describe.
 * @author Robert Marskar
 */

#ifndef CD_CUTCELLFACEOVERRIDES_H
#define CD_CUTCELLFACEOVERRIDES_H

// Std includes
#include <vector>

// Chombo includes
#include <IntVect.H>
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

namespace PolyhedralEB {

/**
 * @brief The faces of a box's cut cells that do not come from the cells' own chords.
 * @details A cut cell's body is built from its surface, and most of its faces are exactly what the surface's chords
 * give. Two kinds of face are not. A face onto a cell the finer level carries is described at the finer level's
 * resolution: the chord is replaced by the finer cells' chords, merged into the polygons stored here. A face onto a
 * covered cell the finer level does not carry is closed, and what it gave up becomes interface in the plane of that
 * face; such a face stores no polygon, only that it is closed. With these the body of every cut cell is a function of
 * its own surface and its own entries here, and needs nothing from its neighbours or from the finer level.
 *
 * Stored compressed, in the order the faces are applied: the cells in increasing order, each with its faces, each
 * face with its polygons, each polygon with its vertices. A cell's restricted faces come before its closed ones,
 * and each kind in increasing face order, which is the order a body applies them in. Cells are ordered as a
 * BoxIterator meets them, the highest direction the most significant, and found by bisection in that order.
 */
class CutCellFaceOverrides
{
public:
  /**
   * @brief Why a face is overridden: restricted to the finer level's chords.
   */
  static constexpr int s_finer = 0;

  /**
   * @brief Why a face is overridden: closed onto a covered cell.
   */
  static constexpr int s_closed = 1;

  /**
   * @brief Constructor. Holds no cell.
   */
  CutCellFaceOverrides() noexcept;

  /**
   * @brief Destructor.
   */
  ~CutCellFaceOverrides() noexcept = default;

  /**
   * @brief Drop every entry.
   */
  void
  clear() noexcept;

  /**
   * @brief Open the entries of a cell.
   * @details Cells are added in increasing BoxIterator order, each once, and a cell's faces follow it.
   * @param[in] a_cell The cell.
   */
  void
  beginCell(const IntVect& a_cell);

  /**
   * @brief Open a face of the cell last begun.
   * @param[in] a_face   Face index, 2 * dir + side.
   * @param[in] a_reason s_finer or s_closed.
   */
  void
  beginFace(const int a_face, const int a_reason);

  /**
   * @brief Add a polygon to the face last begun.
   * @param[in] a_vertices    Vertices in circuit order, in the cell's frame.
   * @param[in] a_vertexEdges Cell edge each vertex came from, or -1.
   * @param[in] a_numVertices Number of vertices.
   */
  void
  addPolygon(const RealVect* a_vertices, const int* a_vertexEdges, const int a_numVertices);

  /**
   * @brief Drop the cell last begun if no face was added to it.
   */
  void
  endCell() noexcept;

  /**
   * @brief Number of cells with entries.
   * @return The count.
   */
  int
  numCells() const noexcept;

  /**
   * @brief A cell's position among the cells with entries.
   * @param[in] a_cell The cell.
   * @return Its index, or -1 if it has none.
   */
  int
  find(const IntVect& a_cell) const noexcept;

  /**
   * @brief The cell at an index.
   * @param[in] a_index Index among the cells.
   * @return The cell.
   */
  const IntVect&
  cell(const int a_index) const noexcept;

  /**
   * @brief The faces of a cell, as a range of face indices.
   * @param[in]  a_index Index among the cells.
   * @param[out] a_begin First face.
   * @param[out] a_end   One past the last face.
   */
  void
  faces(const int a_index, int& a_begin, int& a_end) const noexcept;

  /**
   * @brief Which face of its cell a face entry is.
   * @param[in] a_face Face entry.
   * @return 2 * dir + side.
   */
  int
  face(const int a_face) const noexcept;

  /**
   * @brief Why a face entry overrides its face.
   * @param[in] a_face Face entry.
   * @return s_finer or s_closed.
   */
  int
  reason(const int a_face) const noexcept;

  /**
   * @brief The polygons of a face entry, as a range of polygon indices.
   * @param[in]  a_face  Face entry.
   * @param[out] a_begin First polygon.
   * @param[out] a_end   One past the last polygon.
   */
  void
  polygons(const int a_face, int& a_begin, int& a_end) const noexcept;

  /**
   * @brief The vertices of a polygon, as a range of vertex indices.
   * @param[in]  a_polygon Polygon index.
   * @param[out] a_begin   First vertex.
   * @param[out] a_end     One past the last vertex.
   */
  void
  vertices(const int a_polygon, int& a_begin, int& a_end) const noexcept;

  /**
   * @brief A vertex position, in its cell's frame.
   * @param[in] a_vertex Vertex index.
   * @return The position.
   */
  const RealVect&
  vertex(const int a_vertex) const noexcept;

  /**
   * @brief The cell edge a vertex came from.
   * @param[in] a_vertex Vertex index.
   * @return The edge, or -1.
   */
  int
  vertexEdge(const int a_vertex) const noexcept;

  /**
   * @brief Whether two sets of entries are the same, bit for bit.
   * @param[in] a_other The other set.
   * @return True if nothing differs.
   */
  bool
  equals(const CutCellFaceOverrides& a_other) const noexcept;

protected:
  /**
   * @brief Whether one cell comes before another in BoxIterator order.
   * @param[in] a_first  One cell.
   * @param[in] a_second The other.
   * @return True if a_first is met first.
   */
  static bool
  precedes(const IntVect& a_first, const IntVect& a_second) noexcept;

  /**
   * @brief The cells with entries, in increasing order.
   */
  std::vector<IntVect> m_cells;

  /**
   * @brief First face of every cell, and one past the last face of the last.
   */
  std::vector<int> m_cellFaces;

  /**
   * @brief Face index of every face entry.
   */
  std::vector<int> m_face;

  /**
   * @brief Reason of every face entry.
   */
  std::vector<int> m_reason;

  /**
   * @brief First polygon of every face entry, and one past the last polygon of the last.
   */
  std::vector<int> m_facePolygons;

  /**
   * @brief First vertex of every polygon, and one past the last vertex of the last.
   */
  std::vector<int> m_polygonVertices;

  /**
   * @brief Every vertex position.
   */
  std::vector<RealVect> m_vertices;

  /**
   * @brief Cell edge every vertex came from, or -1.
   */
  std::vector<int> m_vertexEdges;
};

} // namespace PolyhedralEB

#include <CD_NamespaceFooter.H>

#endif
