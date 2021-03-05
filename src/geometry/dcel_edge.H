/*!
  @file   dcel_edge.H
  @brief  Declaration of a dcel_edge class for describing surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _DCEL_EDGE_
#define _DCEL_EDGE_

#include "RealVect.H"

#include <vector>
#include <memory>

namespace dcel {
  
  class vertex;
  class polygon;

  /*!
    @brief Edge class for handling surface tesselations with DCEL meshes
  */
  class edge {
  public:

    /*!
      @brief Constructor
    */
    edge();

    /*!
      @brief Destructor
    */
    ~edge();

    /*!
      @brief Define function
    */
  
    void define(const std::shared_ptr<vertex>& a_vert,
		const std::shared_ptr<edge>& a_pair,
		const std::shared_ptr<edge>& a_next,
		const std::shared_ptr<edge>& a_prev,
		const RealVect                    a_normal = RealVect::Zero);


    /*!
      @brief Set polygon
    */
  
    void set_poly(const std::shared_ptr<polygon>& a_poly);

    /*!
      @brief Set vertex
    */
  
    void set_vert(const std::shared_ptr<vertex>& a_vert);

    /*!
      @brief Set pair
    */
  
    void set_pair(const std::shared_ptr<edge>& a_pair);

    /*!
      @brief Set next edge
    */
  
    void set_next(const std::shared_ptr<edge>& a_next);

    /*!
      @brief Set previous
    */
  
    void set_prev(const std::shared_ptr<edge>& a_prev);
  
    /*!
      @brief Set normal vector
    */
  
    void set_normal(const RealVect a_normal);
  
    /*!
      @brief Get vertex
    */
  
    const std::shared_ptr<vertex>& get_vert() const;

    /*!
      @brief Get vertex
    */
  
    std::shared_ptr<vertex>& get_vert();

    /*!
      @brief Get other vertex
    */
  
    const std::shared_ptr<vertex>& get_other_vert() const;

    /*!
      @brief Get other vertex
    */
  
    std::shared_ptr<vertex>& get_other_vert();

    /*!
      @brief Get the other half edge
    */
  
    const std::shared_ptr<edge>& get_pair() const;

    /*!
      @brief Get the other half edge
    */
  
    std::shared_ptr<edge>& get_pair();

    /*!
      @brief Get the previous half edge
    */
  
    const std::shared_ptr<edge>& get_prev() const;

    /*!
      @brief Get the previous half edge
    */
  
    std::shared_ptr<edge>& get_prev();

    /*!
      @brief Get the next half edge
    */
  
    const std::shared_ptr<edge>& get_next() const;

    /*!
      @brief Get the next half edge
    */
  
    std::shared_ptr<edge>& get_next();

    /*!
      @brief Get polygon
    */
  
    const std::shared_ptr<polygon>& get_poly() const;

    /*!
      @brief Get polygon
    */
  
    std::shared_ptr<polygon>& get_poly();
  
    /*!
      @brief Get pseudo-normal
    */
  
    RealVect get_normal() const;

    /*!
      @brief Get the signed distance between this edge and a point
    */
  
    Real signed_distance(const RealVect a_x0) const;

  protected:

    /*!
      @brief Edge pseudo-normal
    */
    RealVect m_normal;

    /*!
      @brief Vertex at end of edge
    */
    std::shared_ptr<vertex> m_vert; 

    /*!
      @brief Edge pair
    */
    std::shared_ptr<edge> m_pair;

    /*!
      @brief Previous edge
    */
    std::shared_ptr<edge> m_prev; 

    /*!
      @brief Next edge
    */
    std::shared_ptr<edge> m_next;

    /*!
      @brief Reference to polygon
    */
    std::shared_ptr<polygon> m_poly;
  };

}
#endif
