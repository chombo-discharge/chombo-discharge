/*!
  @file   vertex.H
  @brief  Declaration of a vertex class for describing surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _VERTEX_
#define _VERTEX_

#include "dcel_vec.H"

#include <vector>
#include <memory>

namespace dcel {

  class edge;
  class face;

  class vertex {
  public:

    using Vec3 = Vec3T<double>;

    vertex();
    vertex(const Vec3& a_pos);
    vertex(const Vec3& a_pos, const Vec3& a_normal);
    vertex(const vertex& a_otherVertex);
    ~vertex();

    inline
    void define(const Vec3& m_pos, const std::shared_ptr<edge>& a_edge, const Vec3 a_normal = Vec3::Zero) noexcept;

    inline
    void setPosition(const Vec3& a_pos) noexcept;

    inline
    void setNormal(const Vec3& a_normal) noexcept;

    inline
    void setEdge(const std::shared_ptr<edge>& a_edge) noexcept;

    inline
    void addFaceToCache(const std::shared_ptr<face>& a_face) noexcept;

    inline
    void normalizeNormalVector() noexcept;

    inline
    void computeVertexNormalAverage() noexcept;

    inline
    void computeVertexNormalAngleWeighted() noexcept;

    inline
    void clearFaceCache() noexcept;

    inline
    Vec3& getPosition() noexcept;

    inline
    const Vec3& getPosition() const noexcept;

    inline
    std::shared_ptr<edge>& getEdge() noexcept;

    inline
    const std::shared_ptr<edge>& getEdge() const noexcept;

    inline
    Vec3& getNormal() noexcept;

    inline
    const Vec3& getNormal() const noexcept;

    inline
    std::vector<std::shared_ptr<face> >& getFaceCache() noexcept;

    inline
    const std::vector<std::shared_ptr<face> >& getFaceCache() const noexcept;

    inline
    double signedDistance(const Vec3& a_x0) const noexcept;

    inline
    double unsignedDistance2(const Vec3& a_x0) const noexcept; 

  protected:

    Vec3 m_pos;
    Vec3 m_normal;
    
    std::shared_ptr<edge> m_edge;
    
    std::vector<std::shared_ptr<face> > m_faceCache;

  };
}

#include "dcel_vertexI.H"
#endif
