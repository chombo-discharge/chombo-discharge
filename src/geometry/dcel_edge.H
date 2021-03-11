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

    using Vec3 = Vec3T<double>;

    edge();
    
    edge(const edge& a_otherEdge);
    
    edge(const std::shared_ptr<vertex>& a_vertex);
    
    ~edge();

    inline
    auto define(const std::shared_ptr<vertex>& a_vertex,
		const std::shared_ptr<edge>&   a_pairEdge,
		const std::shared_ptr<edge>&   a_nextEdge,
		const std::shared_ptr<edge>&   a_previousEdge,
		const Vec3                     a_normal) noexcept -> void;

    inline
    auto setVertex(const std::shared_ptr<vertex>& a_vertex) noexcept -> void;

    inline
    auto setPairEdge(const std::shared_ptr<edge>& a_pairEdge) noexcept -> void;

    inline
    auto setNextEdge(const std::shared_ptr<edge>& a_nextEdge) noexcept -> void;

    inline
    auto setPreviousEdge(const std::shared_ptr<edge>& a_previousEdge) noexcept -> void;

    inline
    auto setFace(const std::shared_ptr<face>& a_face) noexcept -> void;

    inline
    auto reconcile() noexcept -> void;

    inline
    auto getVertex() noexcept -> std::shared_ptr<vertex>&;

    inline
    auto getVertex() const noexcept -> const std::shared_ptr<vertex>&;

    inline
    auto getOtherVertex() noexcept -> std::shared_ptr<vertex>&;

    inline
    auto getOtherVertex() const noexcept -> const std::shared_ptr<vertex>&;

    inline
    auto getPairEdge() noexcept -> std::shared_ptr<edge>&;

    inline
    auto getPairEdge() const noexcept -> const std::shared_ptr<edge>& ;

    inline
    std::shared_ptr<edge>& getPreviousEdge() noexcept;

    inline
    const std::shared_ptr<edge>& getPreviousEdge() const noexcept;

    inline
    std::shared_ptr<edge>& getNextEdge() noexcept;

    inline
    const std::shared_ptr<edge>& getNextEdge() const noexcept;

    inline
    auto getNormal() noexcept -> Vec3&;

    inline
    auto getNormal() const noexcept -> const Vec3&;

    inline
    std::shared_ptr<face>& getFace() noexcept;

    inline
    const std::shared_ptr<face>& getFace() const noexcept;

    inline
    double signedDistance(const Vec3& a_x0) const noexcept;

    inline
    double unsignedDistance2(const Vec3& a_x0) const noexcept;

  protected:

    Vec3   m_normal;
    Vec3   m_x2x1;
    double m_len2;

    std::shared_ptr<vertex>  m_vertex; 
    std::shared_ptr<edge>    m_pairEdge;
    std::shared_ptr<edge>    m_previousEdge; 
    std::shared_ptr<edge>    m_nextEdge;
    std::shared_ptr<face>    m_face;

    inline
    double projectPointToEdge(const Vec3& a_x0) const noexcept;

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
