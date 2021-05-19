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

#include "CD_NamespaceHeader.H"

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
    void define(const vertexPtr& a_vertex,
		const edgePtr&   a_pairEdge,
		const edgePtr&   a_nextEdge,
		const edgePtr&   a_previousEdge,
		const Vec3                     a_normal) noexcept;

    inline
    void setVertex(const vertexPtr& a_vertex) noexcept;

    inline
    void setPairEdge(const edgePtr& a_pairEdge) noexcept;

    inline
    void setNextEdge(const edgePtr& a_nextEdge) noexcept;

    inline
    void setPreviousEdge(const edgePtr& a_previousEdge) noexcept;

    inline
    void setFace(const facePtr& a_face) noexcept;

    inline
    void reconcile() noexcept;

    inline
    vertexPtr& getVertex() noexcept;

    inline
    const vertexPtr& getVertex() const noexcept;

    inline
    vertexPtr& getOtherVertex() noexcept;

    inline
    const vertexPtr& getOtherVertex() const noexcept;

    inline
    edgePtr& getPairEdge() noexcept;

    inline
    const edgePtr& getPairEdge() const noexcept;

    inline
    edgePtr& getPreviousEdge() noexcept;

    inline
    const edgePtr& getPreviousEdge() const noexcept;

    inline
    edgePtr& getNextEdge() noexcept;

    inline
    const edgePtr& getNextEdge() const noexcept;

    inline
    Vec3T<T>& getNormal() noexcept;

    inline
    const Vec3T<T>& getNormal() const noexcept;

    inline
    facePtr& getFace() noexcept;

    inline
    const facePtr& getFace() const noexcept;

    inline
    T signedDistance(const Vec3& a_x0) const noexcept;

    inline
    T unsignedDistance2(const Vec3& a_x0) const noexcept;

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
    T projectPointToEdge(const Vec3& a_x0) const noexcept;

    inline
    void normalizeNormalVector() noexcept;

    inline
    void computeEdgeLength() noexcept;

    inline
    void computeNormal() noexcept;
  };
}
#include "CD_NamespaceFooter.H"

#include "dcel_edgeI.H"

#endif
