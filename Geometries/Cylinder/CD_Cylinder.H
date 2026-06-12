/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_Cylinder.H
  @brief  Cylinder geometry
  @author Robert Marskar
*/

#ifndef CD_CYLINDER_H
#define CD_CYLINDER_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/**
  @brief This is a class that defines a cylinder geometry. 
*/
class Cylinder : public ComputationalGeometry
{
public:
  Cylinder();
  ~Cylinder();
};

#include <CD_NamespaceFooter.H>

#endif
