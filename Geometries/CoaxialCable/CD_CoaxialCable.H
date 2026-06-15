/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_CoaxialCable.H
  @brief  Coaxial cable geometry
  @author Robert Marskar
*/

#ifndef CD_COAXIALCABLE_H
#define CD_COAXIALCABLE_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/**
  @brief This is a class that defines a coaxial cable geometry 
*/
class CoaxialCable : public ComputationalGeometry
{
public:
  CoaxialCable();
  ~CoaxialCable();
};

#include <CD_NamespaceFooter.H>

#endif
