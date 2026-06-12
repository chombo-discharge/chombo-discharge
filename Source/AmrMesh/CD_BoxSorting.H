/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_BoxSorting.H
  @brief  Declaration of an enum for identifying various ways of sorting grid patches. 
  @author Robert Marskar
*/

#ifndef CD_BOXSORTING_H
#define CD_BOXSORTING_H

#include <CD_NamespaceHeader.H>

/**
  @brief Enum for sorting boxes
*/
enum class BoxSorting
{
  None,
  Std,
  Shuffle,
  Morton,
  Hilbert
};

#include <CD_NamespaceFooter.H>

#endif
