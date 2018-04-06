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
    @brief Constructor. No information about edge
  */
  dcel_vert(const RealVect a_pos);

  /*!
    @brief Destructor
  */
  ~dcel_vert();

  /*!
    @brief Define function
  */
  virtual void define(const RealVect m_pos, const dcel_edge* const a_edge, const RealVect a_normal = RealVect::Zero);

  /*!
    @brief Set position
  */
  virtual void set_pos(const RealVect a_pos);

  /*!
    @brief Set normal vector
  */
  virtual void set_normal(const RealVect a_normal);

  /*!
    @brief Set edge
  */
  virtual void set_edge(const dcel_edge* const a_edge);

  /*!
    @brief Get reference to polygon
  */
  virtual void add_polygon(const dcel_poly* const a_poly);

  /*!
    @brief Get vertex
  */
  virtual const dcel_edge* get_edge() const;

  /*!
    @brief Get all polygons connected to this vertex
  */
  virtual Vector<const dcel_poly*> get_polygons() const;

  /*!
    @brief Get all cached polygons
  */
  virtual Vector<const dcel_poly*> get_polycache() const;

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
  const dcel_edge* m_edge;

  /*!
    @brief Reference to polygons
  */
  Vector<const dcel_poly*> m_polycache;
};
#endif
