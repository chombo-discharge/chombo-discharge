/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_SphereSdf.H
  @brief  Declaration of a signed distance function for sphere
  @author Robert Marskar
*/

#ifndef CD_SPHERESDF_H
#define CD_SPHERESDF_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Signed distance function for sphere
*/
class SphereSdf : public BaseIF
{
public:
  /**
    @brief Constructor
    @param[in] a_center      Sphere center
    @param[in] a_radius      Sphere radius
    @param[in] a_fluidInside Fluid on inside or not. 
  */
  SphereSdf(const RealVect& a_center, const Real& a_radius, const bool& a_fluidInside);

  /**
    @brief Copy constructor
    @param[in] a_inputIF Input cylinder
  */
  SphereSdf(const SphereSdf& a_inputIF);

  /**
    @brief Destructor
  */
  virtual ~SphereSdf();

  /**
    @brief Get distance to sphere. 
    @param[in] a_point Physical position. 
    @return Function value at the given point
  */
  virtual Real
  value(const RealVect& a_point) const;

  /**
    @brief IF factory method
    @return New copy of this implicit function
  */
  virtual BaseIF*
  newImplicitFunction() const;

protected:
  /**
    @brief Sphere center
  */
  RealVect m_center;

  /**
    @brief Cylinder radius
  */
  Real m_radius;

  /**
    @brief Inside or outside domain
  */
  bool m_fluidInside;
};

#include <CD_NamespaceFooter.H>

#endif
