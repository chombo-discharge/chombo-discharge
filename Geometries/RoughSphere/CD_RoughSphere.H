/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_RoughSphere.H
  @brief  Rough sphere geometry
  @author Robert Marskar
*/

#ifndef CD_ROUGHSPHERE_H
#define CD_ROUGHSPHERE_H

#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/**
  @brief This is a class that a rough sphere (either electrode or dielectric)
*/
class RoughSphere : public ComputationalGeometry
{
public:
  RoughSphere();
  ~RoughSphere() override;

  /**
    @brief Copy constructor
    @param[in] a_other Other geometry
  */
  RoughSphere(const RoughSphere& a_other) = default;

  /**
    @brief Copy assignment
    @param[in] a_other Other geometry
    @return Reference to this
  */
  RoughSphere&
  operator=(const RoughSphere& a_other) = default;

  /**
    @brief Move constructor
    @param[in,out] a_other Other geometry
  */
  RoughSphere(RoughSphere&& a_other) = default;

  /**
    @brief Move assignment
    @param[in,out] a_other Other geometry
    @return Reference to this
  */
  RoughSphere&
  operator=(RoughSphere&& a_other) noexcept = default;
};

#include <CD_NamespaceFooter.H>

#endif
