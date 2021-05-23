/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_BVH.H
  @brief  Declaration of boundary volume hierarchies. 
  @author Robert Marskar
*/

#ifndef CD_BVH_H
#define CD_BVH_H

// Std includes
#include <memory>
#include <vector>
#include <functional>
#include <queue>

// Our includes
#include <CD_Vec.H>
#include <CD_NamespaceHeader.H>

namespace BVH {

  template <class T, class P, class BV>
  class NodeT;

  // Implementation of list of primitives. 
  template <class P>
  using PrimitiveListT = std::vector<std::shared_ptr<const P> >;

  // Stop, partitioning, and bounding volume constructor functions. 
  template <class T, class P, class BV>
  using StopFunctionT = std::function<bool(const NodeT<T, P, BV>&)>;
  
  template <class P>
  using PartitionFunctionT = std::function<std::pair<PrimitiveListT<P>, PrimitiveListT<P> >(const PrimitiveListT<P>&)>;

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
    void topDownSortAndPartitionPrimitives(const StopFunction&      a_stopFunc,
					   const PartitionFunction& a_partFunc,
					   const BVConstructor&     a_bvFunc) noexcept;

    inline
    int getDepth() const noexcept;

    inline
    const PrimitiveList& getPrimitives() const noexcept;

    inline
    T pruneOrdered(const Vec3& a_point) const noexcept;

    inline
    T pruneOrdered2(const Vec3& a_point) const noexcept;

    inline
    T pruneUnordered(const Vec3& a_point) const noexcept;

    inline
    T pruneUnordered2(const Vec3& a_point) const noexcept;

    inline
    T prunePriorityQueue(const Vec3& a_point) const noexcept;

    inline
    T prunePriorityQueue2(const Vec3& a_point) const noexcept;

  protected:

    BV m_bv;

    NodeType m_nodeType;

    int m_depth;
      
    PrimitiveList m_primitives;

    NodePtr m_parent;
    NodePtr m_left;
    NodePtr m_right;

    inline
    void setNodeType(const NodeType a_nodeType) noexcept;

    inline
    void setDepth(const int a_depth) noexcept;

    inline
    void insertNode(NodePtr& a_node, const PrimitiveList& a_primitives) noexcept;

    inline
    void setToRegularNode() noexcept;

    inline
    void setPrimitives(const PrimitiveList& a_primitives) noexcept;

    inline
    BV& getBoundingVolume() noexcept;

    inline
    const BV& getBoundingVolume() const noexcept;

    inline
    T getDistanceToBoundingVolume(const Vec3& a_point) const noexcept;

    inline
    T getDistanceToBoundingVolume2(const Vec3& a_point) const noexcept;

    inline
    T getDistanceToPrimitives(const Vec3& a_point) const noexcept;

    inline
    Node& getParent() noexcept;

    inline
    const Node& getParent() const noexcept;

    inline
    Node& getLeft() noexcept;

    inline
    const Node& getLeft() const noexcept;

    inline
    Node& getRight() noexcept;

    inline
    const Node& getRight() const noexcept;

    inline
    NodeType getNodeType() const noexcept;

    inline
    PrimitiveList& getPrimitives() noexcept;

    inline
    void setParent(const NodePtr& a_parent) noexcept;

    inline
    void setLeft(const NodePtr& a_left) noexcept;

    inline
    void setRight(const NodePtr& a_right) noexcept;

    inline
    void pruneOrdered(T& a_closest, const Vec3& a_point) const noexcept;

    inline
    void pruneOrdered2(T& a_minDist2, std::shared_ptr<const P>& a_closest, const Vec3& a_point) const noexcept;

    inline
    void pruneUnordered(T& a_closest, const Vec3& a_point) const noexcept;

    inline
    void pruneUnordered2(T& a_minDist2, std::shared_ptr<const P>& a_closest, const Vec3& a_point) const noexcept;
  };
}

#include <CD_NamespaceFooter.H>

#include <CD_BVHImplem.H>

#endif
