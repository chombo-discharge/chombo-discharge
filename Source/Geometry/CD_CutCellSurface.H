/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_CutCellSurface.H
 * @brief  Declaration of the data a cut cell's embedded boundary is reconstructed from.
 * @author Robert Marskar
 */

#ifndef CD_CUTCELLSURFACE_H
#define CD_CUTCELLSURFACE_H

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
 * @brief The values a cut cell's embedded boundary is reconstructed from.
 * @details Holds the implicit function at the cell's corners and the parameter along each cell
 * edge where it changes sign, both in the cell's own [-1/2, 1/2]^SpaceDim frame. Neither value
 * belongs to the cell alone: a corner is shared by the 2^SpaceDim cells meeting there and an
 * edge by the four cells meeting along it, so two neighbouring cells reconstruct the same
 * surface bit for bit rather than to within a tolerance. That is what makes a face aperture
 * agree when it is computed from either side.
 */
class CutCellSurface
{
public:
  /**
   * @brief Number of cell edges a crossing can sit on.
   */
  static constexpr int s_numEdges = 2 * SpaceDim * (SpaceDim - 1);

  /**
   * @brief Number of cell corners.
   */
  static constexpr int s_numCorners = 1 << SpaceDim;

  /**
   * @brief Marks an edge whose endpoints lie on the same side of the interface.
   */
  static constexpr Real s_noCrossing = -1.0;

  /**
   * @brief Constructor. Leaves every edge without a crossing and every corner at zero.
   */
  CutCellSurface() noexcept;

  /**
   * @brief Destructor.
   */
  ~CutCellSurface() noexcept = default;

  /**
   * @brief Whether the edge carries a crossing.
   * @param[in] a_edge Edge index in [0, s_numEdges).
   * @return True if the implicit function changes sign along the edge.
   */
  inline bool
  hasCrossing(const int a_edge) const noexcept;

  /**
   * @brief Position along an edge where the implicit function changes sign.
   * @details Parameterised from the edge's low corner, so that both cells sharing the edge read
   * the same value. Only meaningful where hasCrossing is true.
   */
  Real m_crossing[s_numEdges];

  /**
   * @brief Implicit function at the cell corners, ordered so that corner c has offset
   * ((c >> d) & 1) in direction d.
   */
  Real m_corner[s_numCorners];
};

#include <CD_NamespaceFooter.H>

#include <CD_CutCellSurfaceImplem.H>

#endif
