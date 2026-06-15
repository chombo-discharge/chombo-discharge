/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_NoisePlane.H
  @brief  Rough sphere geometry
  @author Robert Marskar
*/

#ifndef CD_NOISEPLANE_H
#define CD_NOISEPLANE_H

#include <CD_ComputationalGeometry.H>
#include <CD_NamespaceHeader.H>

/**
  @brief This is a class that a rough sphere (either electrode or dielectric)
*/
class NoisePlane : public ComputationalGeometry
{
public:
  /**
    @brief Constructor.
  */
  NoisePlane();

  /**
    @brief Destructor (does nothing)
  */
  ~NoisePlane() override;

  NoisePlane(const NoisePlane&) = delete;
  NoisePlane&
  operator=(const NoisePlane&) = delete;
  NoisePlane(NoisePlane&&)     = delete;
  NoisePlane&
  operator=(NoisePlane&&) = delete;
};

#include <CD_NamespaceFooter.H>

#endif
