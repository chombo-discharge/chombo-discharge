/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_RegularGeometry.H
  @brief  Regular geometry (no EBs)
  @author Robert Marskar
*/

#ifndef CD_REGULARGEOMETRY_H
#define CD_REGULARGEOMETRY_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/**
  @brief This class implements an empty geometry
*/
class RegularGeometry : public ComputationalGeometry
{
public:
  RegularGeometry();
  ~RegularGeometry();
};

#include <CD_NamespaceFooter.H>

#endif
