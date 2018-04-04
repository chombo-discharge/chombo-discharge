/*!
  @file   dcel_triangle.H
  @brief  Declaration of a triangle class for DCEL surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_TRIANGLE_
#define _DCEL_TRIANGLE_

#include "dcel_poly.H"

/*!
  @brief Abstract class for declaring polygons with DCEL surface tesselations
*/
class dcel_triangle : public dcel_poly {
public:
  
  /*!
    @brief Constructor
  */
  dcel_triangle();

  /*!
    @brief Destructor
  */
  ~dcel_triangle();

  /*!
    @brief Signed distance for triangles. 
  */
  virtual Real signed_distance(const RealVect a_x0) const ;
};

#endif
