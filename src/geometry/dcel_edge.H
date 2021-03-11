/*!
  @file   dcel_edge.H
  @brief  Declaration of a dcel_edge class for describing surface tesselations
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_EDGE_
#define _DCEL_EDGE_

#include "dcel_vec.H"
#include "dcel_face.H"

#include <vector>
#include <memory>

namespace dcel {

  template <class T> class vertexT;
  template <class T> class edgeT;
  template <class T> class faceT;
  template <class T> class edgeIteratorT;

  template<class T>
  class edgeT {
  public:

    using Vec3   = Vec3T<T>;
    using vertex = vertexT<T>;
    using edge   = edgeT<T>;
    using face   = faceT<T>;

    using edgeIterator = edgeIteratorT<T>;

    edgeT();
    
    edgeT(const edge& a_otherEdge);
    
    edgeT(const std::shared_ptr<vertex>& a_vertex);
    
    ~edgeT();

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
    auto getPreviousEdge() noexcept -> std::shared_ptr<edge>&;

    inline
    auto getPreviousEdge() const noexcept -> const std::shared_ptr<edge>& ;

    inline
    auto getNextEdge() noexcept -> std::shared_ptr<edge>& ;

    inline
    auto getNextEdge() const noexcept -> const std::shared_ptr<edge>& ;

    inline
    auto getNormal() noexcept -> Vec3&;

    inline
    auto getNormal() const noexcept -> const Vec3&;

    inline
    std::shared_ptr<face>& getFace() noexcept;

    inline
    const std::shared_ptr<face>& getFace() const noexcept;

    inline
    auto signedDistance(const Vec3& a_x0) const noexcept -> T;

    inline
    auto unsignedDistance2(const Vec3& a_x0) const noexcept -> T;

  protected:

    Vec3   m_normal;
    Vec3   m_x2x1;
    T m_len2;

    std::shared_ptr<vertex>  m_vertex; 
    std::shared_ptr<edge>    m_pairEdge;
    std::shared_ptr<edge>    m_previousEdge; 
    std::shared_ptr<edge>    m_nextEdge;
    std::shared_ptr<face>    m_face;

    inline
    auto projectPointToEdge(const Vec3& a_x0) const noexcept -> T;

    inline
    auto normalizeNormalVector() noexcept -> void;

    inline
    auto computeEdgeLength() noexcept -> void;

    inline
    auto computeNormal() noexcept -> void;
  };

}

#include "dcel_edgeI.H"

#endif
