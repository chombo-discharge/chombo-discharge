/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RodIF.H
  @brief  Declaration of a cylinder with spheres at both ends
  @author Robert marskar
*/

#ifndef CD_RodIF_H
#define CD_RodIF_H

// Chombo includes
#include <BaseIF.H>
#include <IntersectionIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Cylinder with rounded caps at its ends. 
*/
class RodIF : public BaseIF
{
public:
  /*!
    @brief Constructor
    @param[in] a_center1 Cylinder center (at center of sphere)
    @param[in] a_center2 Cylinder center (at center of sphere)
    @param[in] a_radius  Cylinder and sphere radius
    @param[in] a_fluidInside Domain on inside or outside
  */
  RodIF(const RealVect& a_center1, const RealVect& a_center2, const Real& a_radius, const bool& a_fluidInside);

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input cylinder
  */

  RodIF(const RodIF& a_inputIF);

  /*!
    @brief Value function. Returns distance to objcet. 
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
  RefCountedPtr<BaseIF> m_baseif;
};

#include <CD_NamespaceFooter.H>

#endif
