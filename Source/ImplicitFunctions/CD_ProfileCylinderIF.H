/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file    CD_ProfileCylinderIF.H
  @brief   Declaration of a cylinder with a surface profile
  @author  Robert Marskar
*/

#ifndef CD_ProfileCylinderIF_H
#define CD_ProfileCylinderIF_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Implicit function which takes the form of a cylinder with a surface profiles (grooves). 
  @details The grooves are torii that are cut out of the cylinder. Use the input parameters to adjust the height, distance, and width of the grooves. 
  @note This is a 3D object only!
*/
class ProfileCylinderIF : public BaseIF
{
public:
  /*!
    @brief No weak construction, use the full constructor.
  */
  ProfileCylinderIF() = delete;

  /*!
    @brief Full constructor
    @param[in] a_endPoint1        Endpoint of cylinder
    @param[in] a_endPoint2        Other endpoint of cylinder
    @param[in] a_cylinderRadius   Cylinder radius
    @param[in] a_torusMajorRadius Major radius of torii
    @param[in] a_torusMajorRadius Minor radius of torii
    @param[in] a_ccDistance       Distance between torii
    @param[in] a_shift            Shift groove pattern along the axis (relative to cylinder center)
    @param[in] a_roundingRadius   Rounding radius on geometry
    @param[in] a_numLeft          Number of torii to the left of the centerpoint
    @param[in] a_numRight         Number of torii to the right of the centerpoint
    @param[in] a_fluidInside      Fluid on inside or not
  */
  ProfileCylinderIF(const RealVect a_endPoint1,
                    const RealVect a_endPoint2,
                    const Real     a_cylinderRadius,
                    const Real     a_torusMajorRadius,
                    const Real     a_torusMinorRadius,
                    const Real     a_ccDistance,
                    const Real     a_shift,
                    const Real     a_roundingRadius,
                    const int      a_numLeft,
                    const int      a_numRight,
                    const bool     a_fluidInside);

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Other cylinder
  */
  ProfileCylinderIF(const ProfileCylinderIF& a_inputIF);

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~ProfileCylinderIF();

  /*!
    @brief Value function. Returns distance to object.
    @param[in] a_pos Position
  */
  virtual Real
  value(const RealVect& a_pos) const override;

  /*!
    @brief Factory method
  */
  virtual BaseIF*
  newImplicitFunction() const override;

protected:
  /*!
    @brief Implicit function. This is an intersection between a cylinder and several torii. 
  */
  RefCountedPtr<BaseIF> m_baseIF;
};

#include <CD_NamespaceFooter.H>

#endif
