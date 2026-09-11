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

} // namespace PolyhedralEB

#include <CD_NamespaceFooter.H>

#endif
