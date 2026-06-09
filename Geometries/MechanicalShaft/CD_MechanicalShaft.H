/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_MechanicalShaft.H
  @brief  Declaration of a mechanical shaft geometry
  @author Robert Marskar
*/

#ifndef CD_MECHANICALSHAFT_H
#define CD_MECHANICALSHAFT_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_EBGeometryIF.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Mechanical shaft class. Consists of a single ring electrode and an infinitely long> shaft with various shapes.
*/
class MechanicalShaft : public ComputationalGeometry
{
public:
  /**
    @brief Constructor. Does nothing
  */
  MechanicalShaft() noexcept;

  /**
    @brief Destructor
  */
  virtual ~MechanicalShaft() noexcept;

protected:
  /**
    @brief Define the electrode
  */
  virtual void
  defineElectrode() noexcept;

  /**
    @brief Define the dielectric
  */
  virtual void
  defineDielectric() noexcept;

  /**
    @brief Get a cylinder oriented along the z-axis.
    @return Implicit function for a cylinder.
  */
  virtual std::shared_ptr<EBGeometry::ImplicitFunction<Real>>
  getSimpleCylinder() const noexcept;

  /**
    @brief Get a polygonal shaft oriented along the z-axis.
    @return Implicit function for a polygonal shaft.
  */
  virtual std::shared_ptr<EBGeometry::ImplicitFunction<Real>>
  getPolygon() const noexcept;

  /**
    @brief Get a profiled shaft oriented along the z-axis.
    @return Implicit function for a circular profiled shaft.
  */
  virtual std::shared_ptr<EBGeometry::ImplicitFunction<Real>>
  getCircularProfiles() const noexcept;
};

#include <CD_NamespaceFooter.H>

#endif
