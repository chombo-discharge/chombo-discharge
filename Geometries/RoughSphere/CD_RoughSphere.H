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

  RoughSphere(const RoughSphere&) = default;
  RoughSphere&
  operator=(const RoughSphere&) = default;
  RoughSphere(RoughSphere&&)    = default;
  RoughSphere&
  operator=(RoughSphere&&) noexcept = default;
};

#include <CD_NamespaceFooter.H>

#endif
