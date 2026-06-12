/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_EBIntersection.H
  @brief  Enum for putting some logic into EB intersection tests. 
  @author Robert Marskar
*/

#ifndef CD_EBINTERSECTION_H
#define CD_EBINTERSECTION_H

#include <CD_NamespaceHeader.H>

/**
  @brief Enum for putting some logic into how we think about intersection between particles and EBs. 
*/
enum class EBIntersection
{
  Bisection,
  Raycast
};

#include <CD_NamespaceFooter.H>

#endif
