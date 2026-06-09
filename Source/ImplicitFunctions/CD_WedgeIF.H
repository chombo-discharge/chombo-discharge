/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_WedgeIF.H
  @brief  Declaration of an implicit function for a wedge
  @author Robert Marskar
*/

#ifndef CD_WEDGEIF_H
#define CD_WEDGEIF_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Implicit function for a wedge. 
*/
class WedgeIF : public BaseIF
{
public:
  /**
    @brief Base constructor
    @param[in] a_dir    Direction, 0=x,1=y,2=z
    @param[in] a_angle  Opening angle
    @param[in] a_curv   Corner curvature
    @param[in] a_point  Center of circle used for curvature
    @param[in] a_inside Fluid on inside or not
  */
  WedgeIF(const int a_dir, const Real a_angle, const Real a_curv, const RealVect a_point, const bool a_inside);

  /**
    @brief Copy constructor
    @param[in] a_inputIF Another polygon rod
  */
  WedgeIF(const WedgeIF& a_inputIF);

  /**
    @brief Destructor
  */
  virtual ~WedgeIF();

  /**
    @brief Get distance to object
    @param[in] a_pos Position
    @return Function value at the given point
  */
  virtual Real
  value(const RealVect& a_pos) const;

  /**
    @brief Factory method.
    @return New copy of this implicit function
  */
  virtual BaseIF*
  newImplicitFunction() const;

protected:
  /**
    @brief Base function
  */
  RefCountedPtr<BaseIF> m_baseIF;
};

#include <CD_NamespaceFooter.H>

#endif
