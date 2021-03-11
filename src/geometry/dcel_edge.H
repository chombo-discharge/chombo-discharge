/*!
  @file   dcel_edge.H
  @brief  Declaration of a dcel_edge class for describing surface tesselations
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_EDGE_
#define _DCEL_EDGE_

#include "dcel_vec.H"

#include <vector>
#include <memory>

namespace dcel {
  
  class vertex;
  class face;

  class edge {
  public:

    edge();
    
    edge(const edge& a_otherEdge);
    
    edge(const std::shared_ptr<vertex>& a_vertex);
    
    ~edge();

    inline
    void define(const std::shared_ptr<vertex>& a_vertex,
		const std::shared_ptr<edge>&   a_pairEdge,
		const std::shared_ptr<edge>&   a_nextEdge,
		const std::shared_ptr<edge>&   a_previousEdge,
		const Vec3<double>             a_normal = Vec3<double>::Zero);

    inline
    void setVertex(const std::shared_ptr<vertex>& a_vertex) noexcept;

    inline
    void setPairEdge(const std::shared_ptr<edge>& a_pairEdge) noexcept;

    inline
    void setNextEdge(const std::shared_ptr<edge>& a_nextEdge) noexcept;

    inline
    void setPreviousEdge(const std::shared_ptr<edge>& a_previousEdge) noexcept;

    inline
    void setFace(const std::shared_ptr<face>& a_face) noexcept;

    inline
    void reconcile() noexcept;

    inline
    std::shared_ptr<vertex>& getVertex() noexcept;

    inline
    const std::shared_ptr<vertex>& getVertex() const noexcept;

    inline
    std::shared_ptr<vertex>& getOtherVertex() noexcept;

    inline
    const std::shared_ptr<vertex>& getOtherVertex() const noexcept;

    inline
    std::shared_ptr<edge>& getPairEdge() noexcept;

    inline
    const std::shared_ptr<edge>& getPairEdge() const noexcept;

    inline
    std::shared_ptr<edge>& getPreviousEdge() noexcept;

    inline
    const std::shared_ptr<edge>& getPreviousEdge() const noexcept;

    inline
    std::shared_ptr<edge>& getNextEdge() noexcept;

    inline
    const std::shared_ptr<edge>& getNextEdge() const noexcept;

    inline
    Vec3<double>& getNormal() noexcept;

    inline
    const Vec3<double>& getNormal() const noexcept;

    inline
    std::shared_ptr<face>& getFace() noexcept;

    inline
    const std::shared_ptr<face>& getFace() const noexcept;

    inline
    double signedDistance(const Vec3<double>& a_x0) const noexcept;

    inline
    double unsignedDistance2(const Vec3<double>& a_x0) const noexcept;

  protected:

    Vec3<double> m_normal;

    // Store these for efficient queries. 
    double     m_len2;
    Vec3<double> m_x2x1;

    std::shared_ptr<vertex>  m_vertex; 
    std::shared_ptr<edge>    m_pairEdge;
    std::shared_ptr<edge>    m_previousEdge; 
    std::shared_ptr<edge>    m_nextEdge;
    std::shared_ptr<face>    m_face;

    inline
    double projectPointToEdge(const Vec3<double>& a_x0) const noexcept;

    inline
    void normalizeNormalVector() noexcept;

    inline
    void computeEdgeLength() noexcept;

    inline
    void computeNormal() noexcept;
  };

}

#include "dcel_edgeI.H"

#endif
