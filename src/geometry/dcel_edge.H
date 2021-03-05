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
    void setPolygon(const std::shared_ptr<polygon>& a_polygon);
    void setVertex(const std::shared_ptr<vertex>& a_vertex);
    void setPairEdge(const std::shared_ptr<edge>& a_pairEdge);
    void setNextEdge(const std::shared_ptr<edge>& a_nextEdge);
    void setPreviousEdge(const std::shared_ptr<edge>& a_previousEdge);
    void setNormal(const RealVect a_normal);

    // Get functions for vertices
    std::shared_ptr<vertex>& getVertex();
    const std::shared_ptr<vertex>& getVertex() const;

    std::shared_ptr<vertex>& getOtherVertex();
    const std::shared_ptr<vertex>& getOtherVertex() const;

    // Get functions for edges
    std::shared_ptr<edge>& getPairEdge();
    const std::shared_ptr<edge>& getPairEdge() const;

    std::shared_ptr<edge>& getPreviousEdge();  
    const std::shared_ptr<edge>& getPreviousEdge() const;
  
    std::shared_ptr<edge>& getNextEdge();
    const std::shared_ptr<edge>& getNextEdge() const;
    

    // Get functions for polygons
    std::shared_ptr<polygon>& getPolygon();  
    const std::shared_ptr<polygon>& getPolygon() const;

    RealVect& getNormal();
    const RealVect& getNormal() const;

    Real signedDistance(const RealVect a_x0) const;

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
