/*!
  @file   vertex.H
  @brief  Declaration of a vertex class for describing surface tesselations
  @author Robert Marskar
  @date   Apr. 2018
*/

#ifndef _VERTEX_
#define _VERTEX_

#include "dcel_vec.H"
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

    inline
    void define(const RealVect& m_pos, const std::shared_ptr<edge>& a_edge, const RealVect a_normal = RealVect::Zero) noexcept;

    inline
    void setPosition(const RealVect& a_pos) noexcept;

    inline
    void setNormal(const RealVect& a_normal) noexcept;

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


    inline
    std::vector<std::shared_ptr<face> >& getFaceCache() noexcept;

    inline
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
