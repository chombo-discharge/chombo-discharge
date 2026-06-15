/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_DepositionType.H
  @brief  Declaration of deposition types. 
  @author Robert Marskar
*/

#ifndef CD_DEPOSITIONTYPE_H
#define CD_DEPOSITIONTYPE_H

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Deposition types. 
  @details NGP = Nearest grid point, CIC = Cloud-in-cell, TSC = Triangle-shaped cloud
*/
enum class DepositionType
{
  NGP,
  CIC,
  TSC
};

#include <CD_NamespaceFooter.H>

#endif
