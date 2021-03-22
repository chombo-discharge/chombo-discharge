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
    decltype(auto) define(std::vector<facePtr >&   a_faces,
			  std::vector<edgePtr >&   a_edges,
			  std::vector<vertexPtr >& a_vertices) noexcept;
    
    inline
    decltype(auto) sanityCheck() const noexcept;

    inline
    decltype(auto) setSearchAlgorithm(const SearchAlgorithm a_algorithm) noexcept;

    inline
    decltype(auto) setInsideOutsideAlgorithm(const dcel::InsideOutsideAlgorithm a_algorithm) noexcept;

    inline
    decltype(auto) reconcile(VertexNormalWeight a_weight = VertexNormalWeight::Angle) noexcept;

    inline
    decltype(auto) getVertices() noexcept;

    inline
    decltype(auto) getVertices() const noexcept;

    inline
    decltype(auto) getEdges() noexcept;

    inline
    decltype(auto) getEdges() const noexcept;

    inline
    decltype(auto) getFaces() noexcept;

    inline
    decltype(auto) getFaces() const noexcept;

    inline
    decltype(auto) signedDistance(const Vec3& a_x0) const noexcept;

    inline
    decltype(auto) signedDistance(const Vec3& a_x0, SearchAlgorithm a_algorithm) const noexcept;

  protected:

    SearchAlgorithm m_algorithm;

    std::vector<vertexPtr >  m_vertices;
    std::vector<edgePtr >    m_edges;
    std::vector<facePtr >    m_faces;

    inline
    decltype(auto) getAllVertexCoordinates() const noexcept;

    inline
    decltype(auto) reconcileFaces() noexcept;

    inline
    decltype(auto) reconcileEdges() noexcept;

    inline
    decltype(auto) reconcileVertices(VertexNormalWeight a_weight) noexcept;

    inline
    decltype(auto) DirectSignedDistance(const Vec3& a_point) const noexcept;

    inline
    decltype(auto) DirectSignedDistance2(const Vec3& a_point) const noexcept;

    inline
    decltype(auto) incrementWarning(std::map<std::string, int>& a_warnings, const std::string& a_warn) const noexcept;

    inline
    decltype(auto) printWarnings(const std::map<std::string, int>& a_warnings) const noexcept;
  };
}

#include "dcel_meshI.H"

#endif
