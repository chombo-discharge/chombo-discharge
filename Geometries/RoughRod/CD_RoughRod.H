/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_RoughRod.H
  @brief  Rod geometry with surface roughness at tip
  @author Robert Marskar
*/

#ifndef CD_ROUGHROD_H
#define CD_ROUGHROD_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Class that describes a rod electrode with surface roughness
*/
class RoughRod : public ComputationalGeometry
{
public:
  RoughRod();
  ~RoughRod();
};

#include <CD_NamespaceFooter.H>

#endif
