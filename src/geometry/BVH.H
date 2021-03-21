/*!
  @file   BVH.H
  @brief  Declaration of boundary volume hierarchies. 
  @author Robert Marskar
  @date   March 2021
*/

#ifndef _DCEL_BVH_H_
#define _DCEL_BVH_H_

#include "Vec.H"

#include <memory>
#include <vector>
#include <functional>
#include <queue>

namespace BVH {

  template <class T, class P, class BV>
  class NodeT;

  // Implementation of list of primitives. 
  template <class P>
  using PrimitiveListT = std::vector<std::shared_ptr<const P> >;

  // Partitioning functions. 
  template <class P>
  using PartitionFunctionT = std::function<std::pair<PrimitiveListT<P>, PrimitiveListT<P> >(const PrimitiveListT<P>&)>;
    
  template <class T, class P, class BV>
  using StopFunctionT = std::function<bool(const NodeT<T, P, BV>&)>;

  template <class P, class BV>
  using BVConstructorT = std::function<BV(const PrimitiveListT<P>&)>;

  enum class NodeType {
    Regular,
    Leaf,
  };

  // T is the precision for Vec3, P is the primitive type you want to enclose, BV is the bounding volume you use for it. P must supply a function
  // signedDistance(...) and BV must supply a function getDistance. 
  template <class T, class P, class BV>
  class NodeT {
  public:

    using PrimitiveList = PrimitiveListT<P>;
      
    using Vec3       = Vec3T<T>;
    using Node       = NodeT<T, P, BV>;
    using NodePtr    = std::shared_ptr<Node>;
      
    //      using PartitionFunction = PartitionFunctionT<T, P, BV>;
    using PartitionFunction = PartitionFunctionT<P>;
    using StopFunction      = StopFunctionT<T, P, BV>;
    using BVConstructor     = BVConstructorT<P, BV>;

    NodeT();
    NodeT(NodePtr& a_parent);
    NodeT(const std::vector<std::shared_ptr<P> >& a_primitives);
    NodeT(const std::vector<std::shared_ptr<const P> >& a_primitives);
    ~NodeT();

    inline
    decltype(auto) topDownSortAndPartitionPrimitives(const StopFunction&      a_stopFunc,
						     const PartitionFunction& a_partFunc,
						     const BVConstructor&     a_bvFunc) noexcept -> void;

    inline
    decltype(auto) getDepth() const noexcept;

    inline
    decltype(auto) getPrimitives() const noexcept;

    inline
    decltype(auto) pruneOrdered(const Vec3& a_point) const noexcept;

    inline
    decltype(auto) pruneOrdered2(const Vec3& a_point) const noexcept;

    inline
    decltype(auto) pruneUnordered(const Vec3& a_point) const noexcept;

    inline
    decltype(auto) pruneUnordered2(const Vec3& a_point) const noexcept;

    inline
    decltype(auto) prunePriorityQueue(const Vec3& a_point) const noexcept;

    inline
    decltype(auto) prunePriorityQueue2(const Vec3& a_point) const noexcept;

  protected:

    BV m_bv;

    NodeType m_nodeType;

    int m_depth;
      
    PrimitiveList m_primitives;

    NodePtr m_parent;
    NodePtr m_left;
    NodePtr m_right;

    inline
    decltype(auto) setNodeType(const NodeType a_nodeType) noexcept;

    inline
    decltype(auto) setDepth(const int a_depth) noexcept;

    inline
    decltype(auto) insertNode(NodePtr&  a_node, const PrimitiveList& a_primitives) noexcept;

    inline
    decltype(auto) setToRegularNode() noexcept;

    inline
    decltype(auto) setPrimitives(const PrimitiveList& a_primitives) noexcept;

    inline
    decltype(auto) getBoundingVolume() noexcept;

    inline
    decltype(auto) getBoundingVolume() const noexcept;

    inline
    decltype(auto) getDistanceToBoundingVolume(const Vec3& a_point) const noexcept;

    inline
    decltype(auto) getDistanceToBoundingVolume2(const Vec3& a_point) const noexcept;

    inline
    decltype(auto) getDistanceToPrimitives(const Vec3& a_point) const noexcept;

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
    decltype(auto) getPrimitives() noexcept;

    inline
    decltype(auto) setParent(const NodePtr& a_parent) noexcept;

    inline
    decltype(auto) setLeft(const NodePtr& a_left) noexcept;

    inline
    decltype(auto) setRight(const NodePtr& a_right) noexcept;

    inline
    decltype(auto) pruneOrdered(T& a_closest, const Vec3& a_point) const noexcept -> void;

    inline
    decltype(auto) pruneOrdered2(T& a_minDist2, std::shared_ptr<const P>& a_closest, const Vec3& a_point) const noexcept -> void;

    inline
    decltype(auto) pruneUnordered(T& a_closest, const Vec3& a_point) const noexcept -> void;

    inline
    decltype(auto) pruneUnordered2(T& a_minDist2, std::shared_ptr<const P>& a_closest, const Vec3& a_point) const noexcept -> void;
  };
}

#include "BVHI.H"

#endif
