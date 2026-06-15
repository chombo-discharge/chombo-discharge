/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_CopyStrategy.H
  @brief  Simple enum for distinguishing copying strategies
  @author Robert Marskar
*/

#ifndef CD_COPYSTRATEGY_H
#define CD_COPYSTRATEGY_H

#include <CD_NamespaceHeader.H>

/**
  @brief Enum for distinguishing how we copy data
  Valid => valid region
  ValidGhost => valid+ghost region
*/
enum class CopyStrategy
{
  Valid,
  ValidGhost
};
#include <CD_NamespaceFooter.H>

#endif
