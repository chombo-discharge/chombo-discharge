/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_CutCellSurfaceImplem.H
 * @brief  Implementation of CD_CutCellSurface.H
 * @author Robert Marskar
 */

#ifndef CD_CUTCELLSURFACEIMPLEM_H
#define CD_CUTCELLSURFACEIMPLEM_H

// Std includes
#include <cstring>

// Our includes
#include <CD_CutCellSurface.H>
#include <CD_NamespaceHeader.H>

namespace PolyhedralEB {

inline CutCellSurface::CutCellSurface() noexcept
{
  for (int edge = 0; edge < s_numEdges; edge++) {
    m_crossing[edge] = s_noCrossing;
  }

  for (int corner = 0; corner < s_numCorners; corner++) {
    m_corner[corner] = 0.0;
  }
}

inline bool
CutCellSurface::hasCrossing(const int a_edge) const noexcept
{
  return m_crossing[a_edge] != s_noCrossing;
}

inline int
CutCellSurface::linearSize() const noexcept
{
  return static_cast<int>(sizeof(m_crossing) + sizeof(m_corner));
}

inline void
CutCellSurface::linearOut(void* a_buffer) const noexcept
{
  char* buffer = static_cast<char*>(a_buffer);

  std::memcpy(buffer, m_crossing, sizeof(m_crossing));
  std::memcpy(buffer + sizeof(m_crossing), m_corner, sizeof(m_corner));
}

inline void
CutCellSurface::linearIn(const void* a_buffer) noexcept
{
  const char* buffer = static_cast<const char*>(a_buffer);

  std::memcpy(m_crossing, buffer, sizeof(m_crossing));
  std::memcpy(m_corner, buffer + sizeof(m_crossing), sizeof(m_corner));
}

} // namespace PolyhedralEB

#include <CD_NamespaceFooter.H>

#endif
