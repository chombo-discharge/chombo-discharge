/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_RodPlaneProfile.H
  @brief  Rod-plane geometry. The dielectric plane can have surface roughness
  @author Robert Marskar
*/

#ifndef CD_RODPLANEPROFILE_H
#define CD_RODPLANEPROFILE_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Geometry class which defines an electrode blade and a dielectric slab with various surface profiles.
*/
class RodPlaneProfile : public ComputationalGeometry
{
public:
  /**
    @brief Enumeration of supported dielectric surface profiles.
  */
  enum class profile
  {
    circle, ///< Circular surface profile.
    square, ///< Square surface profile.
  };

  /**
    @brief Constructor. Defines the geometry.
  */
  RodPlaneProfile();

  /**
    @brief Destructor.
  */
  ~RodPlaneProfile();

protected:
  /**
    @brief Active surface profile type.
  */
  profile m_profile;

  /**
    @brief Get the base implicit function for the selected profile.
    @return Implicit function for the dielectric slab.
  */
  BaseIF*
  getBaseIF();

  /**
    @brief Get the base implicit function for a circular profile.
    @return Implicit function for the circular surface profile.
  */
  static BaseIF*
  getBaseIFCircle();

  /**
    @brief Get the base implicit function for a square profile.
    @return Implicit function for the square surface profile.
  */
  static BaseIF*
  getBaseIFSquare();
};

#include <CD_NamespaceFooter.H>

#endif
