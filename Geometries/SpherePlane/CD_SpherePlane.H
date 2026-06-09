/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_SpherePlane.H
  @brief  Sphere-plane gap geometry
  @author Robert Marskar
*/

#ifndef CD_SPHEREPLANE_H
#define CD_SPHEREPLANE_H

// our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Sphere-plane gap geometry.
  @details The geometry is always oriented along +y, in both 2D and 3D. 
*/
class SpherePlane : public ComputationalGeometry
{
public:
  /**
    @brief Constructor
  */
  SpherePlane() noexcept;

  /**
    @brief Destructor
  */
  virtual ~SpherePlane() noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
