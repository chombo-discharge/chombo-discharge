/*!
  @brief  surface_element.H
  @brief  Declaration of a surface_element class
  @author Robert Marskar
  @date   Apr. 2018
  @todo   Segregate implementation
*/

#ifndef _SURFACE_ELEMENT_
#define _SURFACE_ELEMENT_

#include "bv_sphere.H"

#include "RealVect.H"

/*!
  @brief Surface_Element class
*/
class surface_element {
public:

  /*!
    @brief Weak constructor. Must subsequently call define. 
  */
  surface_element();

  /*!
    @brief Destructor
  */ 
  virtual ~surface_element();

  /*!
    @brief Get the signed distance from the surface_element
  */
  virtual Real signed_distance(const RealVect a_x) const = 0;

  /*!
    @brief Get element centroid
  */
  virtual RealVect centroid() const = 0;

  /*!
    @brief Get bounding volume sphere
  */
  virtual bv_sphere get_bv_sphere();
  
protected:

  /*!
    @brief Bounding volume sphere
  */
  bv_sphere m_bv;
};
#endif 
