/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_TorusSdf.H
  @brief  Declaration of a signed distance function for a torus
  @author Robert Marskar
*/

#ifndef CD_TORUSSDF_H
#define CD_TORUSSDF_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Signed distance function for a torus (oriented along z).
*/
class TorusSdf : public BaseIF
{
public:
  /**
    @brief No weak construction (yet)
  */
  TorusSdf() = delete;

  /**
    @brief Constructor
    @param[in] a_center      Center
    @param[in] a_majorRadius Major radius
    @param[in] a_minorRadius Minor radius
    @param[in] a_fluidInside Fluid on inside of or not
  */
  TorusSdf(const RealVect a_center, const Real a_majorRadius, const Real a_minorRadius, const bool a_fluidInside);

  /**
    @brief Copy constructor
    @param[in] a_inputIF Input cylinder
  */

  TorusSdf(const TorusSdf& a_inputIF);

  /**
    @brief Destructor
  */
  virtual ~TorusSdf();

  /**
    @brief Test whether inside or not
    @param[in] a_point Physical position. 
    @return Function value at the given point
  */
  virtual Real
  value(const RealVect& a_point) const override;

  /**
    @brief IF factory method
    @return New copy of this implicit function
  */
  virtual BaseIF*
  newImplicitFunction() const override;

protected:
  /**
    @brief Sphere center
  */
  RealVect m_center;

  /**
    @brief Cylinder radius
  */
  Real m_majorRadius;

  /**
    @brief Cylinder radius
  */
  Real m_minorRadius;

  /**
    @brief Inside or outside domain
  */
  bool m_fluidInside;
};

#include <CD_NamespaceFooter.H>

#endif
