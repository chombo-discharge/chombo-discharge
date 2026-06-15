/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_Aerosol.H
  @brief  Declaration of an aerosol geometry
  @author Robert Marskar
*/

#ifndef CD_AEROSOL_H
#define CD_AEROSOL_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Simple class for making a bunch of (manually constructored) dielectric sphere
*/
class Aerosol : public ComputationalGeometry
{
public:
  /**
    @brief Constructor
  */
  Aerosol();

  /**
    @brief Destructor
  */
  ~Aerosol();
};

#include <CD_NamespaceFooter.H>

#endif
