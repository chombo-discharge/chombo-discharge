/*!
  @file   dcel_mesh.H
  @brief  Declaration of a mesh class for handling surface tesselations
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_MESH_
#define _DCEL_MESH_

#include "dcel_algorithms.H"

#include <vector>
#include <memory>
#include <functional>
#include <map>

namespace dcel {

  template <class T> class vertexT;
  template <class T> class edgeT;
  template <class T> class faceT;

  template <class T>
  class meshT {
  public:

    using Vec3   = Vec3T<T>;
    
    using vertex = vertexT<T>;
    using edge   = edgeT<T>;
    using face   = faceT<T>;

    using vertexPtr = std::shared_ptr<vertex>;
    using edgePtr   = std::shared_ptr<edge>;
    using facePtr   = std::shared_ptr<face>;
    
    using mesh   = meshT<T>;

    meshT();
    meshT(const mesh& a_otherMesh) = delete;
    meshT(std::vector<facePtr >&   a_faces,
	  std::vector<edgePtr >&   a_edges,
	  std::vector<vertexPtr >& a_vertices);

    ~meshT();

    inline
    void define(std::vector<facePtr >&   a_faces,
		std::vector<edgePtr >&   a_edges,
		std::vector<vertexPtr >& a_vertices) noexcept;
    
    inline
    void sanityCheck() const noexcept;

    inline
    void setSearchAlgorithm(const SearchAlgorithm a_algorithm) noexcept;

    inline
    void setInsideOutsideAlgorithm(const dcel::InsideOutsideAlgorithm a_algorithm) noexcept;

    inline
    void reconcile(VertexNormalWeight a_weight = VertexNormalWeight::Angle) noexcept;

    inline
    std::vector<vertexPtr>& getVertices() noexcept;

    inline
    const std::vector<vertexPtr>& getVertices() const noexcept;

    inline
    std::vector<edgePtr>& getEdges() noexcept;

    inline
    const std::vector<edgePtr>& getEdges() const noexcept;

    inline
    std::vector<facePtr>& getFaces() noexcept;

    inline
    const std::vector<facePtr>& getFaces() const noexcept;

    inline
    T signedDistance(const Vec3& a_x0) const noexcept;

    inline
    T signedDistance(const Vec3& a_x0, SearchAlgorithm a_algorithm) const noexcept;

  protected:

    SearchAlgorithm m_algorithm;

    std::vector<vertexPtr >  m_vertices;
    std::vector<edgePtr >    m_edges;
    std::vector<facePtr >    m_faces;

    inline
    std::vector<Vec3T<T> > getAllVertexCoordinates() const noexcept;

    inline
    void reconcileFaces() noexcept;

    inline
    void reconcileEdges() noexcept;

    inline
    void reconcileVertices(VertexNormalWeight a_weight) noexcept;

    inline
    T DirectSignedDistance(const Vec3& a_point) const noexcept;

    inline
    T DirectSignedDistance2(const Vec3& a_point) const noexcept;

    inline
    void incrementWarning(std::map<std::string, int>& a_warnings, const std::string& a_warn) const noexcept;

    inline
    void printWarnings(const std::map<std::string, int>& a_warnings) const noexcept;
  };
}

#include "dcel_meshI.H"

#endif
