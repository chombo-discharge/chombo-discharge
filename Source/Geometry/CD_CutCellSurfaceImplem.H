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

inline void
CutCellSurface::store(Real* a_values, const int a_stride) const noexcept
{
  for (int e = 0; e < s_numEdges; e++) {
    a_values[e * a_stride] = m_crossing[e];
  }

  for (int c = 0; c < s_numCorners; c++) {
    a_values[(s_numEdges + c) * a_stride] = m_corner[c];
  }
}

inline void
CutCellSurface::load(const Real* a_values, const int a_stride) noexcept
{
  for (int e = 0; e < s_numEdges; e++) {
    m_crossing[e] = a_values[e * a_stride];
  }

  for (int c = 0; c < s_numCorners; c++) {
    m_corner[c] = a_values[(s_numEdges + c) * a_stride];
  }
}

} // namespace PolyhedralEB

#include <CD_NamespaceFooter.H>

#endif
