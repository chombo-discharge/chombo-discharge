/*!
  @file   dcel_edge.H
  @brief  Declaration of a dcel_edge class for describing surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_VERT_
#define _DCEL_VERT_

#include "RealVect.H"

class dcel_edge;
class dcel_poly;

/*!
  @brief Vertex class for handling surface tesselations with DCEL meshes
*/
class dcel_vert {
public:

  /*!
    @brief Constructor
  */
  dcel_vert();

  /*!
    @brief Destructor
  */
  ~dcel_vert();

  /*!
    @brief Get vertex
  */
  virtual const dcel_edge* get_edge() const;

  /*!
    @brief Get position
  */
  virtual RealVect get_pos() const;

  /*!
    @brief Get pseudo-normal
  */
  virtual RealVect get_normal() const;
  
  /*!
    @brief Get the signed distance between this edge and a point
  */
  virtual Real signed_distance(const RealVect a_x0) const;

protected:

  /*!
    @brief Vertex position
  */
  RealVect m_pos;

  /*!
    @brief Vertex normal
  */
  RealVect m_normal;

  /*!
    @brief One of the half-edges in the mesh
  */
  dcel_edge* m_edge;
};
#endif
