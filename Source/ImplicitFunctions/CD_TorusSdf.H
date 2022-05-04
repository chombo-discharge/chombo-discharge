/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TorusSdf.H
  @brief  Declaration of a signed distance function for a torus
  @author Robert Marskar
*/

#ifndef CD_TorusSdf_H
#define CD_TorusSdf_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Signed distance function for a torus (oriented along z).
*/
class TorusSdf : public BaseIF
{
public:
  /*!
    @brief No weak construction (yet)
  */
  TorusSdf() = delete;

  /*!
    @brief Constructor
    @param[in] a_center      Center
    @param[in] a_majorRadius Major radius
    @param[in] a_minorRadius Minor radius
    @param[in] a_fluidInside Fluid on inside of or not
  */
  TorusSdf(const RealVect a_center, const Real a_majorRadius, const Real a_minorRadius, const bool a_fluidInside);

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input cylinder
  */

  TorusSdf(const TorusSdf& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~TorusSdf();

  /*!
    @brief Test whether inside or not
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
    @brief Sphere center
  */
  RealVect m_center;

  /*!
    @brief Cylinder radius
  */
  Real m_majorRadius;

  /*!
    @brief Cylinder radius
  */
  Real m_minorRadius;

  /*!
    @brief Inside or outside domain
  */
  bool m_fluidInside;
};

#include <CD_NamespaceFooter.H>

#endif
