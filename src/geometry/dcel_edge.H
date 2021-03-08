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

    edge();
    edge(const edge& a_otherEdge);
    ~edge();

    void define(const std::shared_ptr<vertex>& a_vertex,
		const std::shared_ptr<edge>&   a_pairEdge,
		const std::shared_ptr<edge>&   a_nextEdge,
		const std::shared_ptr<edge>&   a_previousEdge,
		const RealVect                 a_normal = RealVect::Zero);

    // Set functions
    void setVertex(const std::shared_ptr<vertex>& a_vertex) noexcept;
    void setPairEdge(const std::shared_ptr<edge>& a_pairEdge) noexcept;
    void setNextEdge(const std::shared_ptr<edge>& a_nextEdge) noexcept;
    void setPreviousEdge(const std::shared_ptr<edge>& a_previousEdge) noexcept;
    void setNormal(const RealVect a_normal) noexcept;
    void setPolygon(const std::shared_ptr<polygon>& a_polygon) noexcept;
    void normalizeNormalVector() noexcept;

    // Get functions for vertices
    std::shared_ptr<vertex>& getVertex() noexcept;
    const std::shared_ptr<vertex>& getVertex() const noexcept;

    std::shared_ptr<vertex>& getOtherVertex() noexcept;
    const std::shared_ptr<vertex>& getOtherVertex() const noexcept;

    // Get functions for edges
    std::shared_ptr<edge>& getPairEdge() noexcept;
    const std::shared_ptr<edge>& getPairEdge() const noexcept;

    std::shared_ptr<edge>& getPreviousEdge() noexcept;  
    const std::shared_ptr<edge>& getPreviousEdge() const noexcept;
  
    std::shared_ptr<edge>& getNextEdge() noexcept;
    const std::shared_ptr<edge>& getNextEdge() const noexcept;

    RealVect& getNormal() noexcept;
    const RealVect& getNormal() const noexcept;

    std::shared_ptr<polygon>& getPolygon() noexcept;  
    const std::shared_ptr<polygon>& getPolygon() const noexcept;

    Real signedDistance(const RealVect a_x0) const noexcept;

  protected:

    RealVect m_normal;

    std::shared_ptr<vertex>  m_vertex; 
    std::shared_ptr<edge>    m_pairEdge;
    std::shared_ptr<edge>    m_previousEdge; 
    std::shared_ptr<edge>    m_nextEdge;
    std::shared_ptr<polygon> m_polygon;
  };

}
#endif
