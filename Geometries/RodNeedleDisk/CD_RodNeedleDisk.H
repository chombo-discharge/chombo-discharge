/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_RodNeedleDisk.H
  @brief  Rod plane electrode geometry, with a needle protruding out of the rod.
  @author Fanny Skirbekk
  @author Robert Marskar
*/

#ifndef CD_RODNEEDLEDISK_H
#define CD_RODNEEDLEDISK_H

//Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/**
  @brief A class that defines a rod electrode with a needle out of the tip, and a plane (electrode) a set distance from the needle tip.
  @details This geometry is an axi-symmetric geometry oriented along one of the axis. All input points are in 1D along the specified axis.
*/
class RodNeedleDisk : public ComputationalGeometry
{
public:
  /**
    @brief Constructor. Puts object in usable state based on the input given in the .options-file
  */
  RodNeedleDisk() noexcept;
  /**
    @brief Destructor. Destroys the object and cleans up when the object goes out of scope.
  */
  virtual ~RodNeedleDisk() noexcept;

protected:
  /**
    @brief Define the rod+needle electrode
  */
  static void
  defineRodNeedle() noexcept;

  /**
    @brief Defines a Plane electrode
  */
  static void
  defineDisk();
};

#include <CD_NamespaceFooter.H>

#endif
