/*!
  @file   dcel_edge.H
  @brief  Declaration of a dcel_edge class for describing surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_VERT_
#define _DCEL_VERT_

#include "RealVect.H"
#include "RefCountedPtr.H"

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
  virtual void define(const RealVect m_pos, const RefCountedPtr<dcel_edge>& a_edge, const RealVect a_normal = RealVect::Zero);

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
  virtual void set_edge(const RefCountedPtr<dcel_edge>& a_edge);

  /*!
    @brief Get reference to polygon
  */
  virtual void add_polygon(const RefCountedPtr<dcel_poly>& a_poly);

  /*!
    @brief Get edge
  */
  virtual const RefCountedPtr<dcel_edge>& get_edge() const;

    /*!
    @brief Get edge
  */
  virtual RefCountedPtr<dcel_edge>& get_edge();

  /*!
    @brief Get all polygons connected to this vertex
  */
  virtual Vector<RefCountedPtr<dcel_poly> > get_polygons();

  /*!
    @brief Get all cached polygons
  */
  virtual Vector<RefCountedPtr<dcel_poly> >& get_polycache();

  /*!
    @brief Get all cached polygons
  */
  virtual const Vector<RefCountedPtr<dcel_poly> >& get_polycache() const;

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
  RefCountedPtr<dcel_edge> m_edge;

  /*!
    @brief Reference to polygons
  */
  Vector<RefCountedPtr<dcel_poly> > m_polycache;
};
#endif
