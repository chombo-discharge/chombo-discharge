/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_CutCellBodyImplem.H
 * @brief  Implementation of CD_CutCellBody.H
 * @author Robert Marskar
 */

#ifndef CD_CUTCELLBODYIMPLEM_H
#define CD_CUTCELLBODYIMPLEM_H

// Our includes
#include <CD_CutCellBody.H>
#include <CD_NamespaceHeader.H>

namespace PolyhedralEB {

/**
 * @brief Visit every cell a body refines into, cutting the subtree once rather than a path per cell.
 * @details Depth-first, so only one set of children is live per level rather than the whole tree:
 * at 4,544 bytes a body that is about 36 kB a level, on the stack of whichever thread is cutting.
 *
 * The visitor is called as a_visit(cell, body, depthRemaining). A cell the interface never enters
 * is handed over once with the depth still to go, because everything under it is the same kind,
 * and enumerating those cells is what would make refinement cost the volume instead of the
 * surface. Only cells the interface does enter are followed to the bottom, where depthRemaining
 * is zero.
 * @param[in] a_body  Body to refine.
 * @param[in] a_cell  Cell the body belongs to, at its own level.
 * @param[in] a_depth How many levels to descend.
 * @param[in] a_visit What to do with each cell reached.
 * @return False if any cut failed, in which case the visit is abandoned part way.
 */
template <typename Visitor>
inline bool
refineSubtree(const CutCellBody& a_body, const IntVect& a_cell, const int a_depth, Visitor&& a_visit) noexcept
{
  // A cell the interface never enters says the same of everything under it, so the visitor is
  // handed it once with the depth still to go rather than the cells it stands for. That is what
  // keeps this proportional to the surface rather than to the volume: enumerating the whole
  // refinement of a whole or empty cell is the cost this exists to avoid.
  if (a_depth <= 0 || a_body.kind() != CutCellBody::Kind::Cut) {
    a_visit(a_cell, a_body, a_depth);

    return true;
  }

  CutCellBody children[1 << SpaceDim];

  if (!a_body.refine(children)) {
    return false;
  }

  for (int c = 0; c < (1 << SpaceDim); c++) {
    IntVect child = 2 * a_cell;

    for (int d = 0; d < SpaceDim; d++) {
      child[d] += (c >> d) & 1;
    }

    if (!refineSubtree(children[c], child, a_depth - 1, a_visit)) {
      return false;
    }
  }

  return true;
}

} // namespace PolyhedralEB

#include <CD_NamespaceFooter.H>

#endif
