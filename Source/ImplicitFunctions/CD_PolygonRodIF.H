/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_PolygonRodIF.H
  @brief  Declaration of a n-sided rod
  @author Robert Marskar
*/

#ifndef CD_POLYGONRODIF_H
#define CD_POLYGONRODIF_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Implicit function class describing a rod with a polygonal cross scetion. 
  @note  This is a 3D object only!
*/
class PolygonRodIF : public BaseIF
{
public:
  /**
    @brief Base constructor
    @param[in] a_endPoint1   One endpoint
    @param[in] a_endPoint2   Other endpoint
    @param[in] a_radius      Outer radius
    @param[in] a_cornerCurv  Rounding radius on corners
    @param[in] a_numSides    Number of sides
    @param[in] a_fluidInside Fluid on inside or not
  */
  PolygonRodIF(const RealVect& a_endPoint1,
               const RealVect& a_endPoint2,
               const Real      a_radius,
               const Real      a_cornerCurv,
               const int       a_numSides,
               const bool      a_fluidInside);

  /**
    @brief Copy constructor
    @param[in] a_inputIF Another polygon rod
  */
  PolygonRodIF(const PolygonRodIF& a_inputIF);

  /**
    @brief Destructor
  */
  virtual ~PolygonRodIF();

  /**
    @brief Value function. Returns distance to object (ish).
    @return Function value at the given point
    @param[in,out] a_pos Description
  */
  virtual Real
  value(const RealVect& a_pos) const;

  /**
    @brief Factory method
    @return New copy of this implicit function
  */
  virtual BaseIF*
  newImplicitFunction() const;

protected:
  /**
    @brief Base function
  */
  RefCountedPtr<BaseIF> m_baseif;
};

#include <CD_NamespaceFooter.H>

#endif
