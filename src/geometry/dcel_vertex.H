/*!
  @file   vertex.H
  @brief  Declaration of a vertex class for describing surface tesselations
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _VERTEX_
#define _VERTEX_

#include "Vec.H"

#include "dcel_edge.H"
#include "dcel_face.H"

#include <vector>
#include <memory>

namespace dcel {

  template <class T> class vertexT;
  template <class T> class edgeT;
  template <class T> class faceT;
  template <class T> class edgeIteratorT;

  template <class T>
  class vertexT {
  public:

    using Vec3   = Vec3T<T>;
    
    using vertex = vertexT<T>;
    using edge   = edgeT<T>;
    using face   = faceT<T>;

    using vertexPtr = std::shared_ptr<vertex>;
    using edgePtr   = std::shared_ptr<edge>;
    using facePtr   = std::shared_ptr<face>;
    
    using edgeIterator = edgeIteratorT<T>;
    
    vertexT();
    vertexT(const Vec3& a_pos);
    vertexT(const Vec3& a_pos, const Vec3& a_normal);
    vertexT(const vertex& a_otherVertex);
    ~vertexT();

    inline
    void define(const Vec3& m_pos, const edgePtr& a_edge, const Vec3& a_normal) noexcept;

    inline
    void setPosition(const Vec3& a_pos) noexcept;

    inline
    void setNormal(const Vec3& a_normal) noexcept;

    inline
    void setEdge(const edgePtr& a_edge) noexcept;

    inline
    void addFace(const facePtr& a_face) noexcept;

    inline
    void normalizeNormalVector() noexcept;

    inline
    void computeVertexNormalAverage() noexcept;

    inline
    void computeVertexNormalAverage(const std::vector<facePtr>& a_faces) noexcept;

    inline
    void computeVertexNormalAngleWeighted() noexcept;

    inline
    void computeVertexNormalAngleWeighted(const std::vector<facePtr>& a_faces) noexcept;

    inline
    Vec3T<T>& getPosition() noexcept;

    inline
    const Vec3T<T>& getPosition() const noexcept;

    inline
    Vec3T<T>& getNormal() noexcept;

    inline
    const Vec3T<T>& getNormal() const noexcept;

    inline
    edgePtr& getOutgoingEdge() noexcept;

    inline
    const edgePtr& getOutgoingEdge() const noexcept;

    inline
    std::vector<facePtr>& getFaces() noexcept;

    inline
    const std::vector<facePtr>& getFaces() const noexcept;

    inline
    T signedDistance(const Vec3& a_x0) const noexcept;

    inline
    T unsignedDistance2(const Vec3& a_x0) const noexcept;

  protected:
    
    edgePtr m_outgoingEdge;
    
    Vec3 m_pos;
    Vec3 m_normal;
    
    std::vector<facePtr > m_faces;


  };
}

#include "dcel_vertexI.H"

#endif
