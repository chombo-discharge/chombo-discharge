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
  template <class T> class EdgeT;
  template <class T> class FaceT;
  template <class T> class edgeIteratorT;

  template<class T>
  class EdgeT {
  public:

    using Vec3   = Vec3T<T>;
    
    using Vertex = VertexT<T>;
    using Edge   = EdgeT<T>;
    using Face   = FaceT<T>;

    using VertexPtr = std::shared_ptr<Vertex>;
    using EdgePtr   = std::shared_ptr<Edge>;
    using FacePtr   = std::shared_ptr<Face>;

    using edgeIterator = edgeIteratorT<T>;

    EdgeT();
    
    EdgeT(const Edge& a_otherEdge);
    
    EdgeT(const VertexPtr& a_vertex);
    
    ~EdgeT();

    inline
    void define(const VertexPtr& a_vertex,
		const EdgePtr&   a_pairEdge,
		const EdgePtr&   a_nextEdge,
		const EdgePtr&   a_previousEdge,
		const Vec3                     a_normal) noexcept;

    inline
    void setVertex(const VertexPtr& a_vertex) noexcept;

    inline
    void setPairEdge(const EdgePtr& a_pairEdge) noexcept;

    inline
    void setNextEdge(const EdgePtr& a_nextEdge) noexcept;

    inline
    void setPreviousEdge(const EdgePtr& a_previousEdge) noexcept;

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
    EdgePtr& getPairEdge() noexcept;

    inline
    const EdgePtr& getPairEdge() const noexcept;

    inline
    EdgePtr& getPreviousEdge() noexcept;

    inline
    const EdgePtr& getPreviousEdge() const noexcept;

    inline
    EdgePtr& getNextEdge() noexcept;

    inline
    const EdgePtr& getNextEdge() const noexcept;

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
    EdgePtr    m_pairEdge;
    EdgePtr    m_previousEdge; 
    EdgePtr    m_nextEdge;
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
