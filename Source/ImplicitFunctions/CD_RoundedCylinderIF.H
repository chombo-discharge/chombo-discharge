/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_RoundedCylinderIF.H
  @brief  Declaration of a cylindrical implicit function with rounded edges
  @author Robert Marskar
*/

#ifndef CD_ROUNDEDCYLINDERIF_H
#define CD_ROUNDEDCYLINDERIF_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Cylinder object with rounded edges (no endcaps)
*/
class RoundedCylinderIF : public BaseIF
{
public:
  /**
    @brief Constructor
    @param[in] a_center1     Cylinder endpoint
    @param[in] a_center2     Other cylinder endpoint
    @param[in] a_radius      Cylinder radius
    @param[in] a_curv        Radius of curvature
    @param[in] a_fluidInside Domain on inside or outside
  */
  RoundedCylinderIF(const RealVect& a_center1,
                    const RealVect& a_center2,
                    const Real      a_radius,
                    const Real      a_curv,
                    const bool      a_fluidInside);

  /**
    @brief Copy constructor
    @param[in] a_inputIF Input cylinder
  */
  RoundedCylinderIF(const RoundedCylinderIF& a_inputIF);

  /**
    @brief Get distance to object. 
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
    @brief Cylinder center
  */
  RealVect m_center1;

  /**
    @brief Cylinder center
  */
  RealVect m_center2;

  /**
    @brief Cylinder length
  */
  Real m_length;

  /**
    @brief Cylinder radius
  */
  Real m_radius;

  /**
    @brief Radius of curvature
  */
  Real m_curv;

  /**
    @brief Inside or outside domain
  */
  bool m_fluidInside;

  /**
    @brief Implicit function
  */
  RefCountedPtr<BaseIF> m_baseIF;

  /**
    @brief Make the implicit function
  */
  void
  makeBaseIF();

  /**
    @brief Make BaseIF
    @return Newly created implicit function for the rounded cylinder.
  */
#if CH_SPACEDIM == 2
  BaseIF*
  makeBaseIF2D();
#elif CH_SPACEDIM == 3
  BaseIF*
  makeBaseIF3D();
#endif
};

#include <CD_NamespaceFooter.H>

#endif
