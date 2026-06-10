/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_RodDielectric.H
  @brief  Definition of a electrode rod and dielectric dielectric geometry
  @author Robert Marskar
*/

#ifndef CD_RODDIELECTRIC_H
#define CD_RODDIELECTRIC_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Class that defines an electrode rod and a dielectric (some shapes supported).
*/
class RodDielectric : public ComputationalGeometry
{
public:
  /**
    @brief Constructor. Defines the geometry.
  */
  RodDielectric();

  /**
    @brief Destructor.
  */
  ~RodDielectric();

protected:
  /**
    @brief Define the electrode geometry.
  */
  static void
  defineElectrode();

  /**
    @brief Define the insulator geometry.
  */
  void
  defineInsulator();

  /**
    @brief Get an axis-aligned box implicit function.
    @return Signed distance function for a box.
  */
  static RefCountedPtr<BaseIF>
  getBox();

  /**
    @brief Get a plane implicit function.
    @return Signed distance function for a plane.
  */
  static RefCountedPtr<BaseIF>
  getPlane();

  /**
    @brief Get a sphere implicit function.
    @return Signed distance function for a sphere.
  */
  static RefCountedPtr<BaseIF>
  getSphere();

  /**
    @brief Get a Perlin-noise-perturbed box implicit function.
    @return Signed distance function for a Perlin box.
  */
  static RefCountedPtr<BaseIF>
  getPerlinBox();
};

#include <CD_NamespaceFooter.H>

#endif
