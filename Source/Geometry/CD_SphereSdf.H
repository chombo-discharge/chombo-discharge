/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_SphereSdf.H
  @brief  Signed distance function for sphere
  @author Robert Marskar
*/

#ifndef CD_SphereSdf_H
#define CD_SphereSdf_H

// Chombo includes
#include <BaseIF.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Signed distance function for sphere
*/
class SphereSdf : public BaseIF {
public:

  /*!
    @brief Constructor
    @param[in] a_center      Sphere center
    @param[in] a_radius      Sphere radius
    @param[in] a_fluidInside Fluid on inside of sphere
  */
  SphereSdf(const RealVect& a_center, const Real& a_radius, const bool& a_fluidInside);

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input cylinder
  */
  
  SphereSdf(const SphereSdf& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~SphereSdf();

  /*!
    @brief Test whether inside or not
    @param[in] a_point Physical position. 
  */
  virtual Real value(const RealVect& a_point) const; 

  /*!
    @brief IF factory method
  */
  virtual BaseIF* newImplicitFunction() const;

protected:

  /*!
    @brief Sphere center
  */
  RealVect m_center;

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
