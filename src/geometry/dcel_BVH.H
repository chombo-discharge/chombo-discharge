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
  namespace BVH {

    template <class Obj, class BV>
    class NodeT;

    template <class Obj>
    using SortingFunction = std::function<bool(const Obj&, const Obj&, const int)>;

    template <class Obj>
    using CostFunction = std::function<std::pair<double, int>(const std::vector<std::shared_ptr<Obj> >&, const int)>;

    template <class Obj, class BV>
    using SortAndPartitionFunction = std::function<int(NodeT<Obj, BV>&)>;

    template <class Obj, class BV>
    using StopFunction = std::function<bool(const NodeT<Obj, BV>&)>;

    template <class Obj, class BV>
    using BVConstructor = std::function<BV(const std::vector<std::shared_ptr<Obj> >&)>;

    enum class NodeType {
      Regular,
      Leaf,
    };

    // Obj is the object type you want to enclose, BV is the bounding volume you use for it. 
    template <class Obj, class BV>
    class NodeT {
    public:

      using Vec3       = Vec3T<double>;
      using Node       = NodeT<Obj, BV>;
      using NodePtr    = std::shared_ptr<Node>;
      using ObjectList = std::vector<std::shared_ptr<Obj> >;

      NodeT();
      NodeT(std::shared_ptr<Node>& a_parent);
      NodeT(const ObjectList& a_data);
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
      decltype(auto) setObjects(const ObjectList& a_data) noexcept;

      inline
      decltype(auto) computeBoundingVolume(const BVConstructor<Obj, BV> bvFunction) noexcept;

      inline
      decltype(auto) topDownSortAndPartitionObjects(SortingFunction<Obj>   sortingFunc,
						    CostFunction<Obj>      costFunc,
						    StopFunction<Obj, BV>  stopFunc,
						    BVConstructor<Obj, BV> bvFunc) noexcept -> void;

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
      
      ObjectList m_objects;

      NodePtr m_parent;
      NodePtr m_left;
      NodePtr m_right;

      inline
      decltype(auto) insertNode(std::shared_ptr<Node>&  a_node, const ObjectList& a_objects) noexcept;

      inline
      decltype(auto) setToRegularNode() noexcept;
      
    };
  }
}

#include "dcel_BVHI.H"

#endif
