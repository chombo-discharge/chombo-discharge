/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CylinderSdf.H
  @brief  Declaration of a cylindrical signed distance function
  @date   Nov. 2017
*/

#ifndef CD_CylinderSdf_H
#define CD_CylinderSdf_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Declaration of a cylinder IF class
*/
class CylinderSdf : public BaseIF
{
public:
  /*!
    @brief Constructor
    @param[in] a_center1 Cylinder center
    @param[in] a_center2 Cylinder center
    @param[in] a_radius  Cylinder radius
    @param[in] a_inside  Domain on inside or outside
  */
  CylinderSdf(const RealVect& a_center1, const RealVect& a_center2, const Real& a_radius, const bool& a_fluidInside);

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input cylinder
  */
  CylinderSdf(const CylinderSdf& a_inputIF);

  /*!
    @brief Test whether inside or not
    @param[in] a_point Physical position. 
  */
  virtual Real
  value(const RealVect& a_point) const;

  /*!
    @brief IF factory method
  */
  virtual BaseIF*
  newImplicitFunction() const;

protected:
  /*!
    @brief One end poitn
  */
  RealVect m_endPoint1;

  /*!
    @brief The other end point. 
  */
  RealVect m_endPoint2;

  /*!
    @brief Top
  */
  RealVect m_top;

  /*!
    @brief Top
  */
  RealVect m_axis;

  /*!
    @brief Center
  */
  RealVect m_center;

  /*!
    @brief Cylinder length
  */
  Real m_length;

  /*!
    @brief Cylinder radius
  */
  Real m_radius;

  /*!
    @brief Inside or outside domain
  */
  bool m_fluidInside;
};

#include <CD_NamespaceFooter.H>

#endif
