/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_EBRepresentation.H
  @brief  Enum for putting some logic into EB intersection tests. 
  @author Robert Marskar
*/

#ifndef CD_EBREPRESENTATION_H
#define CD_EBREPRESENTATION_H

#include <CD_NamespaceHeader.H>

/**
  @brief Enum for putting some logic into how we think about EBs. This is just a simply supporting class for various solvers/algorithms that can distinguish
  between ways of think about the EB: 
*/
enum class EBRepresentation
{
  ImplicitFunction,
  Discrete,
  Voxel,
  Levelset
};

#include <CD_NamespaceFooter.H>

#endif
