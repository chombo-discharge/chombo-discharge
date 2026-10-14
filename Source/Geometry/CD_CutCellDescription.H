/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_CutCellDescription.H
 * @brief  Declaration of everything one cut cell of the polyhedral geometry is built from.
 * @author Robert Marskar
 */

#ifndef CD_CUTCELLDESCRIPTION_H
#define CD_CUTCELLDESCRIPTION_H

// Chombo includes
#include <IntVect.H>

// Our includes
#include <CD_CutCellFaceOverrides.H>
#include <CD_CutCellSurface.H>
#include <CD_NamespaceHeader.H>

namespace PolyhedralEB {

/**
 * @brief Everything one cut cell of the polyhedral geometry is built from, and which cell it is.
 * @details The cell's level and index, its surface -- corner values and edge crossings, and the role that says whether
 * its faces are its own -- and the face overrides that replace the faces that are not, holding this cell's entry alone.
 * Self-contained: CutCellBody::define(m_surface, m_overrides, m_cell) builds the cell's polyhedron from it, and
 * nothing else is needed. No interior is stored; the polyhedron, its fan of interface triangles included, is built
 * by whoever wants it.
 */
struct CutCellDescription
{
  /**
   * @brief Geometry level of the cell.
   */
  int m_level;

  /**
   * @brief The cell, on its level.
   */
  IntVect m_cell;

  /**
   * @brief The cell's state, as its level holds it.
   */
  int m_state;

  /**
   * @brief The cell's surface.
   */
  CutCellSurface m_surface;

  /**
   * @brief The cell's face overrides; empty where its role is native.
   */
  CutCellFaceOverrides m_overrides;
};

} // namespace PolyhedralEB

#include <CD_NamespaceFooter.H>

#endif
