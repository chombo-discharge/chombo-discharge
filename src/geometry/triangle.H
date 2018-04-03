/*!
  @brief  triangle.H
  @brief  Declaration of a triangle class
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _TRIANGLE_
#define _TRIANGLE_

#include "surface_element.H"

/*!
  @brief Triangle class
*/
class triangle : public surface_element {
public:

  /*!
    @brief Weak constructor. Must subsequently call define. 
  */
  triangle();

  /*!
    @brief Full constructor
  */
  triangle(const RealVect a_x1, const RealVect a_x2, const RealVect a_x3, const RealVect a_normal);

  /*!
    @brief Destructor
  */
  virtual ~triangle();

  /*!
    @brief Define function
  */
  virtual void define(const RealVect a_x1, const RealVect a_x2, const RealVect a_x3, const RealVect a_normal);

  /*!
    @brief Get the signed distance
  */
  virtual Real signed_distance(const RealVect a_x0) const;

  /*!
    @brief Get centroid
  */
  virtual RealVect centroid() const;
  
protected:

  /*!
    @brief Triangle corner
  */
  RealVect m_x1;

  /*!
    @brief Triangle corner
  */
  RealVect m_x2;

  /*!
    @brief Triangle corner
  */
  RealVect m_x3;

  /*!
    @brief Auxiliary vector
  */
  RealVect m_v1;

    /*!
    @brief Auxiliary vector
  */
  RealVect m_v2;

    /*!
    @brief Auxiliary vector
  */
  RealVect m_v3;

  /*!
    @brief Triangle centroid
  */
  RealVect m_centroid;

  /*!
    @brief Triangle normal
  */
  RealVect m_normal;

  /*!
    @brief Use the 3D signed distance method
  */
  virtual Real signed_distance3D(const RealVect a_x0) const;

  /*!
    @brief Use the 2D signed distance method
  */
  virtual Real signed_distance2D(const RealVect a_x0) const;
};
#endif

