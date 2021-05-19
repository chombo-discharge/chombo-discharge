/*!
  @file new_sphere_if.H
  @brief New SphereIF function, the old one didn't use the correct distance (return distance squared)
  @date Sep. 2017
  @author Robert Marskar
*/

#ifndef _NEW_SPHERE_IF
#define _NEW_SPHERE_IF_

#include <BaseIF.H>
#include <IntersectionIF.H>

#include "CD_NamespaceHeader.H"

/*!
  @brief Declaration of a hemispheric rod IF class
  @details This class is esseentially a cylinder with a hemisphere at one of its ends. 
*/
class new_sphere_if : public BaseIF {
public:

  /*!
    @brief Constructor
    @param[in] a_center      Sphere center
    @param[in] a_radius      Sphere radius
    @param[in] a_fluidInside Fluid on inside of sphere
  */

  new_sphere_if(const RealVect& a_center, const Real& a_radius, const bool& a_fluidInside);

  /*!
    @brief Copy constructor
    @param[in] a_inputIF Input cylinder
  */
  
  new_sphere_if(const new_sphere_if& a_inputIF);

  /*!
    @brief Destructor
  */
  virtual ~new_sphere_if();

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
#include "CD_NamespaceFooter.H"
#endif
