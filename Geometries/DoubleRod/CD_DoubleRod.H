/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_DoubleRod.H
  @brief  Double rod geometry.
  @author Robert Marskar
*/

#ifndef CD_DOUBLEROD_H
#define CD_DOUBLEROD_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/**
  @brief This is a class that defines two dielectric rods.
*/
class DoubleRod : public ComputationalGeometry
{
public:
  DoubleRod();
  ~DoubleRod();
};

#include <CD_NamespaceFooter.H>

#endif
