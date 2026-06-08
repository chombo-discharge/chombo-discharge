/*
 * SPDX-FileCopyrightText: 2025-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_Rod.H
  @brief  File containing headers for a rod-like geometry.
  @author Robert Marskar
*/

#ifndef CD_ROD_H
#define CD_ROD_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_EBGeometryIF.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Geometry class representing a single cylindrical rod electrode.
*/
class Rod : public ComputationalGeometry
{
public:
  /**
    @brief Full constructor.
    @details This function parses all input parameters and creates the geometry.
  */
  Rod() noexcept;

  /**
    @brief Destructor (does nothing).
  */
  virtual ~Rod() noexcept = default;
};

#include <CD_NamespaceFooter.H>

#endif
