/* chombo-discharge
 * Copyright © 2022 NTNU.
 * Copyright © 2022 Fanny Skirbekk. 
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_NeedleIF.H
  @brief  Declaration of needle electrode (cylinder with one pointy end and one flat end)
  @author Fanny Skirbekk
*/

#ifndef CD_NeedleIF_H
#define CD_NeedleIF_H

// Chombo includes
#include <BaseIF.H>
#include <IntersectionIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Cylinder with one pointed end and one flat end
*/
class NeedleIF : public BaseIF
{
public:
  /*!
    @brief Constructor
    @param[in] a_length        Needle length
    @param[in] a_radius        Needle radius
    @param[in] a_fluidInside   Fluid inside cylinder-part of needle, true/false
    @param[in] a_tipRadius     Radius of needle tip curvature           
    @param[in] a_angle         Needle tip (cone) opening angle
    @param[in] a_cornerCurve   Corner curve between rod and needle
  */
  NeedleIF(const Real& a_length,
           const Real& a_radius,
           const bool& a_fluidInside,
           const Real& a_tipRadius,
           const Real& a_angle,
           const Real& a_cornerCurve);

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input Needle
  */
  NeedleIF(const NeedleIF& a_inputIF);

  /*!
    @brief Value function. Returns distance to object.
    @param[in] a_point Physical position.
  */
  virtual Real
  value(const RealVect& a_point) const override;

  /*!
    @brief IF factory method
  */
  virtual BaseIF*
  newImplicitFunction() const override;

protected:
  RefCountedPtr<BaseIF> m_baseif;
  /*!
    @brief curvature-radius of needle tip
  */
  Real m_tipRadius;
};

#include <CD_NamespaceFooter.H>

#endif
