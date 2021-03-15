/*!
  @file   dcel_BVH.H
  @brief  Declaration of boundary volume hierarchies for use with dcel meshes.
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_BVH_H_
#define _DCEL_BVH_H_

#include "dcel_vec.H"

#include <memory>
#include <vector>
#include <functional>

namespace dcel {

  template <class T>
  using SortingFunction = std::function<bool(const T&, const T&, const int)>;

  template <class T>
  using CostFunction = std::function<std::pair<double, int>(const std::vector<std::shared_ptr<T> >&, const int)>;

  template <class T>
  using StopFunction = std::function<bool(const std::vector<std::shared_ptr<T> >&, const int)>;

  template <class T, class BV>
  using BVConstructor = std::function<BV(const std::vector<std::shared_ptr<T> >&)>;
  
  namespace BVH {

    enum class NodeType {
      Regular,
      Leaf,
    };
        
    template <class T, class BV>
    class NodeT {
    public:

      using Vec3 = Vec3T<double>;
      using Node = NodeT<T, BV>;

      NodeT();
      NodeT(std::shared_ptr<Node>& a_parent);
      NodeT(std::vector<std::shared_ptr<T> >& a_data);
      ~NodeT();

      inline
      decltype(auto) setParent(const std::shared_ptr<Node>& a_parent) noexcept;

      inline
      decltype(auto) setLeft(const std::shared_ptr<Node>& a_left) noexcept;

      inline
      decltype(auto) setRight(const std::shared_ptr<Node>& a_right) noexcept;

      inline
      decltype(auto) setNodeType(const NodeType a_nodeType) noexcept;

      inline
      decltype(auto) setDepth(const int a_depth) noexcept;

      inline
      decltype(auto) setObjects(const std::vector<std::shared_ptr<T> >& a_data) noexcept;

      inline
      decltype(auto) computeBoundingVolume(BVConstructor<T, BV> bvFunction) noexcept;

      inline
      decltype(auto) topDownSortAndPartitionObjects(SortingFunction<T>   sortingFunc,
						    CostFunction<T>      costFunc,
						    StopFunction<T>      stopFunc,
						    BVConstructor<T, BV> bvFunc) noexcept -> void;

      inline
      decltype(auto) getParent() noexcept;

      inline
      decltype(auto) getParent() const noexcept;

      inline
      decltype(auto) getLeft() noexcept;

      inline
      decltype(auto) getLeft() const noexcept;

      inline
      decltype(auto) getRight() noexcept;

      inline
      decltype(auto) getRight() const noexcept;

      inline
      decltype(auto) getNodeType() const noexcept;
      
      inline
      decltype(auto) getDepth() const noexcept;

      inline
      decltype(auto) getSignedDistanceToPoint(const Vec3& a_x0) const noexcept;

      inline
      decltype(auto) getBoundingVolume() noexcept;

      inline
      decltype(auto) getBoundingVolume() const noexcept;

      inline
      decltype(auto) getObjects() noexcept;

      inline
      decltype(auto) getObjects() const noexcept;

    protected:

      BV m_bv;

      NodeType m_nodeType;

      int m_depth;
      
      std::vector<std::shared_ptr<T> > m_objects;

      std::shared_ptr<Node> m_parent;
      std::shared_ptr<Node> m_left;
      std::shared_ptr<Node> m_right;
    };
  }
}

#include "dcel_BVHI.H"

#endif
