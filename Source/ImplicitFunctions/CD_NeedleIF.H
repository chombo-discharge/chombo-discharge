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
#include <RefCountedPtr.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Cylinder with one pointed end and one circular end. 
*/
class NeedleIF : public BaseIF
{
public:
  /*!
    @brief Constructor
    @param[in] a_length      Needle length
    @param[in] a_radius      Needle radius
    @param[in] a_tipRadius   Radius of needle tip curvature           
    @param[in] a_angle       Needle tip (cone) opening angle
    @param[in] a_cornerCurve Corner curve between rod and needle
    @param[in] flipInside    Flip inside to outside
  */
  NeedleIF(const Real& a_length,
           const Real& a_radius,
           const Real& a_tipRadius,
           const Real& a_angle,
           const Real& a_cornerCurve,
           const bool& a_flipInside);

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
  /*!
    @brief Base implicit function. This is a union between a cylinder and a cone.
  */
  RefCountedPtr<BaseIF> m_implicitFunction;

  /*!
    @brief curvature-radius of needle tip
  */
  Real m_tipRadius;
};

#include <CD_NamespaceFooter.H>

#endif
