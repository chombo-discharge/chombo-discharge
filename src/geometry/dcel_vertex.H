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
    decltype(auto) define(const Vec3& m_pos, const edgePtr& a_edge, const Vec3& a_normal) noexcept;

    inline
    decltype(auto) setPosition(const Vec3& a_pos) noexcept;

    inline
    decltype(auto) setNormal(const Vec3& a_normal) noexcept;

    inline
    decltype(auto) setEdge(const edgePtr& a_edge) noexcept;

    inline
    decltype(auto) addFace(const facePtr& a_face) noexcept;

    inline
    decltype(auto) normalizeNormalVector() noexcept;

    inline
    decltype(auto) computeVertexNormalAverage() noexcept;

    inline
    decltype(auto) computeVertexNormalAngleWeighted() noexcept;

    inline
    decltype(auto) getPosition() noexcept;

    inline
    decltype(auto) getPosition() const noexcept;

    inline
    decltype(auto) getNormal() noexcept;

    inline
    decltype(auto) getNormal() const noexcept;

    inline
    decltype(auto) getOutgoingEdge() noexcept;

    inline
    decltype(auto) getOutgoingEdge() const noexcept;

    inline
    decltype(auto) getFaces() noexcept;

    inline
    decltype(auto) getFaces() const noexcept;

    inline
    decltype(auto) signedDistance(const Vec3& a_x0) const noexcept;

    inline
    decltype(auto) unsignedDistance2(const Vec3& a_x0) const noexcept;

  protected:
    
    edgePtr m_outgoingEdge;
    
    Vec3 m_pos;
    Vec3 m_normal;
    
    std::vector<facePtr > m_faces;


  };
}

#include "dcel_vertexI.H"

#endif
