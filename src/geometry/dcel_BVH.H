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

    template <class T, class P, class BV>
    class NodeT;

    template <class P>
    using PrimitiveListT = std::vector<std::shared_ptr<P> >;

    template <class T, class P, class BV>
    using PartitionFunctionT = std::function<int(NodeT<T, P, BV>&)>;
    
    template <class T, class P, class BV>
    using StopFunctionT = std::function<bool(const NodeT<T, P, BV>&)>;

    template <class T, class P, class BV>
    using BVConstructorT = std::function<BV(const NodeT<T, P, BV>&)>;

    template <class T, class P>
    //    using DistanceFunctionT = std::function<T(const std::vector<std::shared_ptr<P> >&, const Vec3T<T>& a_point)>;
    using DistanceFunctionT = std::function<T(const PrimitiveListT<P>&, const Vec3T<T>& a_point)>;

    enum class NodeType {
      Regular,
      Leaf,
    };

    // T is the precision for Vec3, P is the primitive type you want to enclose, BV is the bounding volume you use for it. 
    template <class T, class P, class BV>
    class NodeT {
    public:

      using PrimitiveList = PrimitiveListT<P>;
      
      using Vec3       = Vec3T<T>;
      using Node       = NodeT<T, P, BV>;
      using NodePtr    = std::shared_ptr<Node>;
      
      using PartitionFunction = PartitionFunctionT<T, P, BV>;
      using StopFunction      = StopFunctionT<T, P, BV>;
      using BVConstructor     = BVConstructorT<T, P, BV>;
      using DistanceFunction  = DistanceFunctionT<T, P>;

      NodeT();
      NodeT(NodePtr& a_parent);
      NodeT(const PrimitiveList& a_primitives);
      ~NodeT();

      inline
      decltype(auto) setNodeType(const NodeType a_nodeType) noexcept;

      inline
      decltype(auto) setDepth(const int a_depth) noexcept;

      inline
      decltype(auto) topDownSortAndPartitionPrimitives(const StopFunction&      a_stopFunc,
						       const PartitionFunction& a_partFunc,
						       const BVConstructor&     a_bvFunc) noexcept -> void;

      inline
      decltype(auto) getDepth() const noexcept;

      inline
      decltype(auto) getPrimitives() noexcept;

      inline
      decltype(auto) getPrimitives() const noexcept;

      inline
      decltype(auto) pruneBranches(const Vec3& a_point, T& a_closest, const DistanceFunction& a_distFunc) const noexcept -> void;

    protected:

      BV m_bv;

      NodeType m_nodeType;

      int m_depth;
      
      PrimitiveList m_primitives;

      NodePtr m_parent;
      NodePtr m_left;
      NodePtr m_right;

      inline
      decltype(auto) insertNode(NodePtr&  a_node, const PrimitiveList& a_primitives) noexcept;

      inline
      decltype(auto) setToRegularNode() noexcept;

      inline
      decltype(auto) setPrimitives(const PrimitiveList& a_primitives) noexcept;

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
