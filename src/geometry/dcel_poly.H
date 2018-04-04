/*!
  @file   dcel_poly.cpp
  @brief  Declaration of a polygon class for DCEL surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_POLY_
#define _DCEL_POLY_

#include "RealVect.H" 

class dcel_vert;
class dcel_edge;

/*!
  @brief Abstract class for declaring polygons with DCEL surface tesselations
*/
class dcel_poly {
public:
  
  /*!
    @brief Constructor
  */
  dcel_poly();

  /*!
    @brief Destructor
  */
  ~dcel_poly();

  /*!
    @brief Define edge
  */
  virtual void define(const RealVect a_normal, const dcel_edge* const a_edge);

  /*!
    @brief Get edge
  */
  virtual const dcel_edge* get_edge() const;

  /*!
    @brief Get face normal vector
  */
  virtual RealVect get_normal() const;

  /*!
    @brief Signed distance. This should be pure in order to do general polygons. 
  */
  virtual Real signed_distance(const RealVect a_x0) const = 0;
  
protected:

  /*!
    @brief One of the edges of the polygon
  */
  const dcel_edge* m_edge;

  /*!
    @brief Face normal
  */
  RealVect m_normal;
};

#endif
