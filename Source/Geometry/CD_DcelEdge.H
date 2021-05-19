/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_DcelEdge.H
  @brief  Declaration of an edge class for DCEL tools
  @author Robert Marskar
*/

#ifndef CD_DcelEdge_H
#define CD_DcelEdge_H

// Std includes
#include <vector>
#include <memory>

// Our includes
#include <CD_Vec.H>
#include <CD_DcelFace.H>
#include <CD_NamespaceHeader.H>

namespace Dcel {

  template <class T> class VertexT;
  template <class T> class edgeT;
  template <class T> class FaceT;
  template <class T> class edgeIteratorT;

  template<class T>
  class edgeT {
  public:

    using Vec3   = Vec3T<T>;
    
    using Vertex = VertexT<T>;
    using edge   = edgeT<T>;
    using Face   = FaceT<T>;

    using VertexPtr = std::shared_ptr<Vertex>;
    using edgePtr   = std::shared_ptr<edge>;
    using FacePtr   = std::shared_ptr<Face>;

    using edgeIterator = edgeIteratorT<T>;

    edgeT();
    
    edgeT(const edge& a_otherEdge);
    
    edgeT(const VertexPtr& a_vertex);
    
    ~edgeT();

    inline
    void define(const VertexPtr& a_vertex,
		const edgePtr&   a_pairEdge,
		const edgePtr&   a_nextEdge,
		const edgePtr&   a_previousEdge,
		const Vec3                     a_normal) noexcept;

    inline
    void setVertex(const VertexPtr& a_vertex) noexcept;

    inline
    void setPairEdge(const edgePtr& a_pairEdge) noexcept;

    inline
    void setNextEdge(const edgePtr& a_nextEdge) noexcept;

    inline
    void setPreviousEdge(const edgePtr& a_previousEdge) noexcept;

    inline
    void setFace(const FacePtr& a_face) noexcept;

    inline
    void reconcile() noexcept;

    inline
    VertexPtr& getVertex() noexcept;

    inline
    const VertexPtr& getVertex() const noexcept;

    inline
    VertexPtr& getOtherVertex() noexcept;

    inline
    const VertexPtr& getOtherVertex() const noexcept;

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
    FacePtr& getFace() noexcept;

    inline
    const FacePtr& getFace() const noexcept;

    inline
    T signedDistance(const Vec3& a_x0) const noexcept;

    inline
    T unsignedDistance2(const Vec3& a_x0) const noexcept;

  protected:

    Vec3   m_normal;

    // These are not necessary and exist because of performance reasons. 
    Vec3 m_x2x1;
    T    m_invLen2;

    VertexPtr  m_vertex; 
    edgePtr    m_pairEdge;
    edgePtr    m_previousEdge; 
    edgePtr    m_nextEdge;
    FacePtr    m_face;

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

#include <CD_NamespaceFooter.H>

#include <CD_DcelEdgeImplem.H>

#endif
