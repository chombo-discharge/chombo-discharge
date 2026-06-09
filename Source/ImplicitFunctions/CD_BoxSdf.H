/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_BoxSdf.H
  @brief  Declare an interface for defining an axis-aligned box
  @author Robert Marskar
*/

#ifndef CD_BOXSDF_H
#define CD_BOXSDF_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Class for defining a two- or three-dimensional box with arbitrary centroid and orientation.
*/
class BoxSdf : public BaseIF
{
public:
  /**
    @brief Base constructor
    @param[in] a_loCorner Lower left corner
    @param[in] a_hiCorner Upper right corner
    @param[in] a_fluidInside Fluid inside
  */
  BoxSdf(const RealVect& a_loCorner, const RealVect& a_hiCorner, const bool& a_fluidInside);

  /**
    @brief Copy constructor
    @param[in] a_inputIF Another BoxSdf
  */
  BoxSdf(const BoxSdf& a_inputIF);

  /**
    @brief Destructor
  */
  virtual ~BoxSdf();

  /**
    @brief Test if object is on inside or outside
    @return Function value at the given point
    @param[in,out] a_pos Description
  */
  virtual Real
  value(const RealVect& a_pos) const;

  /**
    @brief IF Factory method
    @return New copy of this implicit function
  */
  virtual BaseIF*
  newImplicitFunction() const;

protected:
  /**
    @brief Lower left corner
  */
  RealVect m_loCorner;

  /**
    @brief Upper right corner
  */
  RealVect m_hiCorner;

  /**
    @brief Fluid inside box or not. 
  */
  bool m_fluidInside;
};

#include <CD_NamespaceFooter.H>

#endif
