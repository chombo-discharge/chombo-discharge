/*!
  @file   vertex.H
  @brief  Declaration of a vertex class for describing surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _VERTEX_
#define _VERTEX_

#include "RealVect.H"

#include <vector>
#include <memory>

namespace dcel {

  class edge;
  class face;

  class vertex {
  public:

    vertex();
    vertex(const RealVect& a_pos);
    vertex(const RealVect& a_pos, const RealVect& a_normal);
    vertex(const vertex& a_otherVertex);
    ~vertex();

    void define(const RealVect& m_pos, const std::shared_ptr<edge>& a_edge, const RealVect a_normal = RealVect::Zero) noexcept;
    void setPosition(const RealVect& a_pos) noexcept;
    void setNormal(const RealVect& a_normal) noexcept;
    void setEdge(const std::shared_ptr<edge>& a_edge) noexcept;
    void addFaceToCache(const std::shared_ptr<face>& a_face) noexcept;
    void normalizeNormalVector() noexcept;
    void computeVertexNormalAverage() noexcept;
    void computeVertexNormalAngleWeighted() noexcept;
    void clearFaceCache() noexcept;

    inline
    RealVect& getPosition() noexcept;

    inline
    const RealVect& getPosition() const noexcept;

    inline
    std::shared_ptr<edge>& getEdge() noexcept;

    inline
    const std::shared_ptr<edge>& getEdge() const noexcept;

    inline
    RealVect& getNormal() noexcept;

    inline
    const RealVect& getNormal() const noexcept;

    std::vector<std::shared_ptr<face> > getFaces() noexcept;

    std::vector<std::shared_ptr<face> >& getFaceCache() noexcept;
    const std::vector<std::shared_ptr<face> >& getFaceCache() const noexcept;
  
    Real signedDistance(const RealVect& a_x0) const noexcept;
    Real unsignedDistance2(const RealVect& a_x0) const noexcept; // Returns squared distance to point

  protected:

    RealVect m_pos;
    RealVect m_normal;
    
    std::shared_ptr<edge> m_edge;
    
    std::vector<std::shared_ptr<face> > m_faceCache;

  };
}

#include "dcel_vertexI.H"
#endif
