/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_DiskProfiledPlane.H
  @brief  Declaration of a geometry describing a disk electrode and a profiled dielectric surface.
  @author Robert Marskar
*/

#ifndef CD_DISKPROFILEDPLANE_H
#define CD_DISKPROFILEDPLANE_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Geometry class for disk electrode and profiled surface experiments. 
  @details This geometry consists of a disk electrode oriented along +y (user can shift it), and a profiled
  plane with specified profiles: squares or square channels, spheres, or cylinders. Profiled are organized
  along a user-specified lattice (coinciding with the coordinate directions), and can be shifted. 
*/
class DiskProfiledPlane : public ComputationalGeometry
{
public:
  /**
    @brief Constructor - reads in all parameters
  */
  DiskProfiledPlane() noexcept;

  /**
    @brief Destructor - does nothing
  */
  virtual ~DiskProfiledPlane() = default;

protected:
  /**
    @brief Define the electrode
  */
  virtual void
  defineElectrode() noexcept;

  /**
    @brief Define the dielectric plate.
  */
  virtual void
  defineDielectric() noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
