/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_DcelVertex.H
  @brief  Declaration of a vertex class for describing surface tesselations
  @author Robert Marskar
*/

#ifndef CD_DcelVertex_H
#define CD_DcelVertex_H

// Std includes
#include <vector>
#include <memory>

// Our includes
#include <CD_Vec.H>
#include <CD_DcelEdge.H>
#include <CD_DcelFace.H>
#include <CD_NamespaceHeader.H>

namespace Dcel {

  template <class T> class VertexT;
  template <class T> class EdgeT;
  template <class T> class FaceT;
  template <class T> class edgeIteratorT;

  template <class T>
  class VertexT {
  public:

    using Vec3   = Vec3T<T>;
    
    using Vertex = VertexT<T>;
    using Edge   = EdgeT<T>;
    using Face   = FaceT<T>;

    using VertexPtr = std::shared_ptr<Vertex>;
    using EdgePtr   = std::shared_ptr<Edge>;
    using FacePtr   = std::shared_ptr<Face>;
    
    using edgeIterator = edgeIteratorT<T>;
    
    VertexT();
    VertexT(const Vec3& a_pos);
    VertexT(const Vec3& a_pos, const Vec3& a_normal);
    VertexT(const Vertex& a_otherVertex);
    ~VertexT();

    inline
    void define(const Vec3& m_pos, const EdgePtr& a_edge, const Vec3& a_normal) noexcept;

    inline
    void setPosition(const Vec3& a_pos) noexcept;

    inline
    void setNormal(const Vec3& a_normal) noexcept;

    inline
    void setEdge(const EdgePtr& a_edge) noexcept;

    inline
    void addFace(const FacePtr& a_face) noexcept;

    inline
    void normalizeNormalVector() noexcept;

    inline
    void computeVertexNormalAverage() noexcept;

    inline
    void computeVertexNormalAverage(const std::vector<FacePtr>& a_faces) noexcept;

    inline
    void computeVertexNormalAngleWeighted() noexcept;

    inline
    void computeVertexNormalAngleWeighted(const std::vector<FacePtr>& a_faces) noexcept;

    inline
    Vec3T<T>& getPosition() noexcept;

    inline
    const Vec3T<T>& getPosition() const noexcept;

    inline
    Vec3T<T>& getNormal() noexcept;

    inline
    const Vec3T<T>& getNormal() const noexcept;

    inline
    EdgePtr& getOutgoingEdge() noexcept;

    inline
    const EdgePtr& getOutgoingEdge() const noexcept;

    inline
    std::vector<FacePtr>& getFaces() noexcept;

    inline
    const std::vector<FacePtr>& getFaces() const noexcept;

    inline
    T signedDistance(const Vec3& a_x0) const noexcept;

    inline
    T unsignedDistance2(const Vec3& a_x0) const noexcept;

  protected:
    
    EdgePtr m_outgoingEdge;
    
    Vec3 m_pos;
    Vec3 m_normal;
    
    std::vector<FacePtr > m_faces;
  };
}

#include <CD_NamespaceFooter.H>

#include <CD_DcelVertexImplem.H>

#endif
