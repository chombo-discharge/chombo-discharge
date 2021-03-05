/*!
  @file   ritter_sphere.H
  @brief  Declaration of class that uses Ritter's method for computing a bounding sphere for a set of points
  @author Robert Marskar
  @date   Ap. 2018
*/

#ifndef _RITTER_SPHERE_
#define _RITTER_SPHERE_

#include "RealVect.H"

/*!
  @brief Class for computing a bounding sphere for a set of points by using Ritter's method
*/
class ritter_sphere {
public:

  /*!
    @brief Weak constructor
  */
  ritter_sphere();

  /*!
    @brief Full constructor
  */
  ritter_sphere(const Vector<RealVect>& a_points);

  /*!
    @brief Destructor
  */
  virtual ~ritter_sphere();

  /*!
    @brief Define sphere
  */
  virtual void define(const Vector<RealVect>& a_points);

  /*!
    @brief Check if point is outside bounding sphere
    @return True if point is outside, false otherwise
  */
  virtual bool inside(const RealVect a_x0) const;

  /*!
    @brief Get distance to ritter sphere
  */
  virtual Real dist(const RealVect a_x0) const;

  /*!
    @brief Get radius
  */
  virtual Real get_radius() const;

  /*!
    @brief Get center
  */
  virtual RealVect get_center() const;

protected:

  /*!
    @brief Sphere radius
  */
  Real m_radius;

  /*!
    @brief Sphere center
  */
  RealVect m_center;
};
#endif
