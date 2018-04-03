/*!
  @brief  bv_sphere.H
  @brief  Declaration of a bounding volume sphere
  @author Robert Marskar
  @date   Apr. 2018
  @todo   Segregate implementation
*/

#ifndef _BV_SPHERE_
#define _BV_SPHERE_

#include "RealVect.H"

/*!
  @brief Bv_Sphere class
*/
class bv_sphere {
public:

  /*!
    @brief Weak constructor. Must subsequently call define. 
  */
  bv_sphere();

  /*!
    @brief Full constructor
  */
  bv_sphere(const RealVect a_center, const Real a_radius);

  /*!
    @brief Destructor
  */ 
  virtual ~bv_sphere();

  /*!
    @brief Get sphere center
  */
  virtual RealVect get_center() const;

  /*!
    @brief Return radius
  */
  virtual Real get_radius() const;

protected:

  /*!
    @brief Center of bounding volume sphere
  */
  RealVect m_center;

  /*!
    @brief Bounding volume sphere radius
  */
  Real m_radius;
};
#endif 
