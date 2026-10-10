/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_GraphCopyHeader.H
 * @brief  Declaration of the record that opens every piece of a polyhedral graph sent to another layout.
 * @author Robert Marskar
 */

#ifndef CD_GRAPHCOPYHEADER_H
#define CD_GRAPHCOPYHEADER_H

// Std includes
#include <type_traits>

// Chombo includes
#include <Box.H>

// Our includes
#include <CD_NamespaceHeader.H>

namespace PolyhedralEB {

/**
 * @brief What a piece of a polyhedral graph is and where it goes, as it travels between ranks.
 * @details A piece is a box of cells of one source level, all held by one tile of that level and all owned by it,
 * that one destination box needs. Plain data, copied into and out of an MPI buffer as bytes: the sender writes one
 * ahead of every piece it sends, and the receiver, which works out the same pieces from the same replicated data,
 * reads them back to learn what arrived and where it belongs. The counts say how much of the piece's payload follows.
 */
struct GraphCopyHeader
{
  /**
   * @brief Level of the destination box in the destination hierarchy.
   */
  int m_destLevel;

  /**
   * @brief Index of the destination box in its level's layout.
   */
  int m_destBox;

  /**
   * @brief Geometry level the cells come from.
   */
  int m_sourceLevel;

  /**
   * @brief Index of the box of the source level's graph layout that holds the cells.
   */
  int m_sourceBox;

  /**
   * @brief Low corner of the piece, on the source level.
   */
  int m_lo[SpaceDim];

  /**
   * @brief High corner of the piece, on the source level.
   */
  int m_hi[SpaceDim];

  /**
   * @brief Number of surfaces the piece carries.
   */
  int m_numSurfaces;

  /**
   * @brief Bytes of face overrides the piece carries.
   */
  int m_overrideBytes;

  /**
   * @brief The piece as a box.
   * @return The box on the source level.
   */
  inline Box
  region() const noexcept
  {
    return Box(IntVect(D_DECL(m_lo[0], m_lo[1], m_lo[2])), IntVect(D_DECL(m_hi[0], m_hi[1], m_hi[2])));
  }
};

// The header travels between ranks as raw bytes.
static_assert(std::is_trivially_copyable<GraphCopyHeader>::value, "GraphCopyHeader must be trivially copyable");

} // namespace PolyhedralEB

#include <CD_NamespaceFooter.H>

#endif
