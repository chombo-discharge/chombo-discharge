/*!
  @file   dcel_edge.H
  @brief  Declaration of a dcel_edge class for describing surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_VERT_
#define _DCEL_VERT_

#include "RealVect.H"

#include <vector>
#include <memory>

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
  dcel_vert(const RealVect& a_pos);

  /*!
    @brief Destructor
  */
  ~dcel_vert();

  /*!
    @brief Define function
  */
  
  void define(const RealVect& m_pos, const std::shared_ptr<dcel_edge>& a_edge, const RealVect a_normal = RealVect::Zero);

  /*!
    @brief Set position
  */
  
  void set_pos(const RealVect& a_pos);

  /*!
    @brief Set normal vector
  */
  
  void set_normal(const RealVect& a_normal);

  /*!
    @brief Set edge
  */
  
  void set_edge(const std::shared_ptr<dcel_edge>& a_edge);

  /*!
    @brief Get reference to polygon
  */
  
  void add_polygon(const std::shared_ptr<dcel_poly>& a_poly);

  /*!
    @brief Get edge
  */
  
  const std::shared_ptr<dcel_edge>& get_edge() const;

    /*!
    @brief Get edge
  */
  
  std::shared_ptr<dcel_edge>& get_edge();

  /*!
    @brief Get all polygons connected to this vertex
  */
  
  std::vector<std::shared_ptr<dcel_poly> > get_polygons();

  /*!
    @brief Get all cached polygons
  */
  
  std::vector<std::shared_ptr<dcel_poly> >& get_polycache();

  /*!
    @brief Get all cached polygons
  */
  
  const std::vector<std::shared_ptr<dcel_poly> >& get_polycache() const;

  /*!
    @brief Get position
  */
  
  const RealVect& get_pos() const;

  /*!
    @brief Get pseudo-normal
  */
  
  const RealVect& get_normal() const;
  
  /*!
    @brief Get the signed distance between this edge and a point
  */
  
  Real signed_distance(const RealVect a_x0) const;

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
  std::shared_ptr<dcel_edge> m_edge;

  /*!
    @brief Reference to polygons
  */
  std::vector<std::shared_ptr<dcel_poly> > m_polycache;
};

#endif
