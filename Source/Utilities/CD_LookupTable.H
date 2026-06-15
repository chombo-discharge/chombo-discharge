/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_LookupTable.H
  @brief  Declaration of a table for looking up coefficients etc. 
  @author Robert Marskar
*/

#ifndef CD_LOOKUPTABLE_H
#define CD_LOOKUPTABLE_H

// Std includes
#include <memory>
#include <vector>
#include <type_traits>
#include <iostream>

// Our includes
#include <CD_NamespaceHeader.H>

namespace LookupTable {
  /**
  @brief Strategy for obtaining data when requesting data that is out of range. 
  We can use either constant (i.e. virtually padding the data with the first/last entry) or 
  with linear interpolation (i.e. virtually padding the data with a slope connecting the first
  two or last two entries in the table). 
*/
  enum class OutOfRangeStrategy
  {
    Constant,
    Interpolate
  };

  /**
  @brief Enum for classifying the coordinate system.
*/
  enum class Spacing
  {
    Uniform,
    Exponential
  };

} // namespace LookupTable

#include <CD_NamespaceFooter.H>

#include <CD_LookupTable1D.H>

#endif
