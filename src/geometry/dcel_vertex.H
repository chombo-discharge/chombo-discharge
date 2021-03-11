/*!
  @file   vertex.H
  @brief  Declaration of a vertex class for describing surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _VERTEX_
#define _VERTEX_

#include "dcel_vec.H"
#include "dcel_edge.H"
#include "dcel_face.H"

#include <vector>
#include <memory>

namespace dcel {

  template <class T> class vertexT;
  template <class T> class edgeT;
  template <class T> class faceT;
  template <class T> class edge_iteratorT;

  template <class T>
  class vertexT {
  public:

    using Vec3   = Vec3T<T>;
    using vertex = vertexT<T>;
    using edge   = edgeT<T>;
    using face   = faceT<T>;
    
    using edge_iterator = edge_iteratorT<T>;

    vertexT();
    vertexT(const Vec3& a_pos);
    vertexT(const Vec3& a_pos, const Vec3& a_normal);
    vertexT(const vertex& a_otherVertex);
    ~vertexT();

    inline
    auto define(const Vec3& m_pos, const std::shared_ptr<edge>& a_edge, const Vec3 a_normal = Vec3::Zero) noexcept -> void;

    inline
    auto setPosition(const Vec3& a_pos) noexcept -> void;

    inline
    auto setNormal(const Vec3& a_normal) noexcept -> void;

    inline
    auto setEdge(const std::shared_ptr<edge>& a_edge) noexcept -> void;

    inline
    auto addFaceToCache(const std::shared_ptr<face>& a_face) noexcept -> void;

    inline
    auto normalizeNormalVector() noexcept -> void;

    inline
    auto computeVertexNormalAverage() noexcept -> void;

    inline
    auto computeVertexNormalAngleWeighted() noexcept -> void;

    inline
    auto clearFaceCache() noexcept -> void;

    inline
    auto getPosition() noexcept -> Vec3&;

    inline
    auto getPosition() const noexcept -> const Vec3&;

    inline
    auto getNormal() noexcept -> Vec3&;

    inline
    auto getNormal() const noexcept -> const Vec3&;

    inline
    auto getEdge() noexcept -> std::shared_ptr<edge>&;

    inline
    auto getEdge() const noexcept -> const std::shared_ptr<edge>& ;

    inline
    auto getFaceCache() noexcept -> std::vector<std::shared_ptr<face> >&;

    inline
    auto getFaceCache() const noexcept -> const std::vector<std::shared_ptr<face> >& ;

    inline
    auto signedDistance(const Vec3& a_x0) const noexcept -> T;

    inline
    auto unsignedDistance2(const Vec3& a_x0) const noexcept -> T;

  protected:

    Vec3 m_pos;
    Vec3 m_normal;
    
    std::shared_ptr<edge> m_edge;
    
    std::vector<std::shared_ptr<face> > m_faceCache;

  };
}

#include "dcel_vertexI.H"

#endif
