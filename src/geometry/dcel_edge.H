/*!
  @file   dcel_edge.H
  @brief  Declaration of a dcel_edge class for describing surface tesselations
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_EDGE_
#define _DCEL_EDGE_

#include "Vec.H"

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

    using vertexPtr = std::shared_ptr<vertex>;
    using edgePtr   = std::shared_ptr<edge>;
    using facePtr   = std::shared_ptr<face>;

    using edgeIterator = edgeIteratorT<T>;

    edgeT();
    
    edgeT(const edge& a_otherEdge);
    
    edgeT(const vertexPtr& a_vertex);
    
    ~edgeT();

    inline
    decltype(auto) define(const vertexPtr& a_vertex,
			  const edgePtr&   a_pairEdge,
			  const edgePtr&   a_nextEdge,
			  const edgePtr&   a_previousEdge,
			  const Vec3                     a_normal) noexcept;

    inline
    decltype(auto) setVertex(const vertexPtr& a_vertex) noexcept;

    inline
    decltype(auto) setPairEdge(const edgePtr& a_pairEdge) noexcept;

    inline
    decltype(auto) setNextEdge(const edgePtr& a_nextEdge) noexcept;

    inline
    decltype(auto) setPreviousEdge(const edgePtr& a_previousEdge) noexcept;

    inline
    decltype(auto) setFace(const facePtr& a_face) noexcept;

    inline
    decltype(auto) reconcile() noexcept;

    inline
    decltype(auto) getVertex() noexcept;

    inline
    decltype(auto) getVertex() const noexcept;

    inline
    decltype(auto) getOtherVertex() noexcept;

    inline
    decltype(auto) getOtherVertex() const noexcept;

    inline
    decltype(auto) getPairEdge() noexcept;

    inline
    decltype(auto) getPairEdge() const noexcept;

    inline
    decltype(auto) getPreviousEdge() noexcept;

    inline
    decltype(auto) getPreviousEdge() const noexcept;

    inline
    decltype(auto) getNextEdge() noexcept;

    inline
    decltype(auto) getNextEdge() const noexcept;

    inline
    decltype(auto) getNormal() noexcept;

    inline
    decltype(auto) getNormal() const noexcept;

    inline
    decltype(auto) getFace() noexcept;

    inline
    decltype(auto) getFace() const noexcept;

    inline
    decltype(auto) signedDistance(const Vec3& a_x0) const noexcept;

    inline
    decltype(auto) unsignedDistance2(const Vec3& a_x0) const noexcept;

    inline
    decltype(auto) getFlippedEdge() const noexcept;

  protected:

    Vec3   m_normal;

    // These are not necessary and exist because of performance reasons. 
    Vec3 m_x2x1;
    T    m_invLen2;

    vertexPtr  m_vertex; 
    edgePtr    m_pairEdge;
    edgePtr    m_previousEdge; 
    edgePtr    m_nextEdge;
    facePtr    m_face;

    inline
    decltype(auto) projectPointToEdge(const Vec3& a_x0) const noexcept;

    inline
    decltype(auto) normalizeNormalVector() noexcept;

    inline
    decltype(auto) computeEdgeLength() noexcept;

    inline
    decltype(auto) computeNormal() noexcept;
  };

}

#include "dcel_edgeI.H"

#endif
