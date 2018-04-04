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
  @brief Implementation of a triangle class for DCEL surface tesselations. This implements the signed 
  distance function for a triangle. 
  @details For cases in which a point lies closest to a triangle vertex or edge, the signed distance function
  for dcel_edge is called. 
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
    @brief Define function. I'm overriding this because I need something more general
  */
  virtual void define(const RealVect a_normal, const dcel_edge* const a_edge);

  /*!
    @brief Compute triangle area
  */
  virtual void compute_area();

  /*!
    @brief Signed distance for triangles. 
  */
  virtual Real signed_distance(const RealVect a_x0) const;

protected:

  RealVect m_x1;
  RealVect m_x2;
  RealVect m_x3;
  RealVect m_v1v2;
  RealVect m_v2v3;
  RealVect m_v1v3;
  Real     m_len2; // Length of the v2v3 vector
  
};

#endif
