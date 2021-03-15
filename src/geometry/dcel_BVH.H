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

    template <class T, class Obj, class BV>
    class NodeT;

    template <class T, class Obj, class BV>
    using PartitionFunctionT = std::function<int(NodeT<T, Obj, BV>&)>;
    
    template <class T, class Obj, class BV>
    using StopFunctionT = std::function<bool(const NodeT<T, Obj, BV>&)>;

    template <class T, class Obj, class BV>
    using BVConstructorT = std::function<BV(const NodeT<T, Obj, BV>&)>;

    enum class NodeType {
      Regular,
      Leaf,
    };

    // T is the precision for Vec3, Obj is the object type you want to enclose, BV is the bounding volume you use for it. 
    template <class T, class Obj, class BV>
    class NodeT {
    public:

      using ObjectList = std::vector<std::shared_ptr<Obj> >;

      using Vec3       = Vec3T<T>;
      
      using Node       = NodeT<T, Obj, BV>;
      using NodePtr    = std::shared_ptr<Node>;
      
      using PartitionFunction = PartitionFunctionT<T, Obj, BV>;
      using StopFunction      = StopFunctionT<T, Obj, BV>;
      using BVConstructor     = BVConstructorT<T, Obj, BV>;

      NodeT();
      NodeT(NodePtr& a_parent);
      NodeT(const ObjectList& a_data);
      ~NodeT();

      inline
      decltype(auto) setNodeType(const NodeType a_nodeType) noexcept;

      inline
      decltype(auto) setDepth(const int a_depth) noexcept;

      inline
      decltype(auto) topDownSortAndPartitionObjects(const StopFunction&      a_stopFunc,
						    const PartitionFunction& a_partFunc,
						    const BVConstructor&     a_bvFunc) noexcept -> void;

      inline
      decltype(auto) getDepth() const noexcept;

      inline
      decltype(auto) getObjects() noexcept;

      inline
      decltype(auto) getObjects() const noexcept;

      inline
      decltype(auto) signedDistance(const Vec3& a_point, T& a_pruneDistance) const noexcept;

    protected:

      BV m_bv;

      NodeType m_nodeType;

      int m_depth;
      
      ObjectList m_objects;

      NodePtr m_parent;
      NodePtr m_left;
      NodePtr m_right;

      inline
      decltype(auto) insertNode(NodePtr&  a_node, const ObjectList& a_objects) noexcept;

      inline
      decltype(auto) setToRegularNode() noexcept;

      inline
      decltype(auto) setObjects(const ObjectList& a_data) noexcept;

      inline
      decltype(auto) computeBoundingVolume(const BVConstructor& bvFunction) noexcept;

      inline
      decltype(auto) getBoundingVolume() noexcept;

      inline
      decltype(auto) getBoundingVolume() const noexcept;

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
      decltype(auto) getDistanceToBV(const Vec3& a_x0) const noexcept;

      inline
      decltype(auto) getDistanceToData(const Vec3& a_x0) const noexcept;
      
      inline
      decltype(auto) getLeafData() const noexcept -> ObjectList;

      inline
      decltype(auto) pruneTree(const Vec3& a_point, T& a_pruneDistance) const noexcept -> ObjectList;

      inline
      decltype(auto) setParent(const NodePtr& a_parent) noexcept;

      inline
      decltype(auto) setLeft(const NodePtr& a_left) noexcept;

      inline
      decltype(auto) setRight(const NodePtr& a_right) noexcept;
      
    };
  }
}

#include "dcel_BVHI.H"

#endif
