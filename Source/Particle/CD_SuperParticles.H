/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_SuperParticles.H
  @brief  Namespace containing various particle merging utilities.
  @author Robert Marskar
*/

#ifndef CD_SuperParticles_H
#define CD_SuperParticles_H

// Std includes
#include <memory>
#include <vector>

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Namespace for various particle merging tools. 
*/
namespace SuperParticles {

  // clang-format off
  /*!
    @brief Particle list alias
  */
  template <class P>
  using ParticleList = std::vector<P>;

  /*!
    @brief Partitioning function for dividing a particle list into two subsets.
  */
  template <class P>
  using ParticlePartitioner = std::function<std::pair<ParticleList<P>, ParticleList<P>>(ParticleList<P>& a_particles)>;  

  /*!
    @brief Polymorphic function for ending particle splitting/merging.
    @details Should return true if we should keep partitioning the particles and false otherwise.
  */
  template <class P>
  using StopFunction = std::function<bool(const ParticleList<P>& a_particles)>;
  // clang-format on

  /*!
    @brief Node in a particle-merging KD-tree. 
    @details This node type is used for partitioning particles into spatial subsets and merging/splitting them. 
    @details The template argument P is the particle type - the users decides how to merge/split/partition. 
  */
  template <class P>
  class KDNode
  {
  public:

    KDNode(const KDNode&)  = delete;
    KDNode(const KDNode&&) = delete;

    KDNode&
    operator=(const KDNode&) = delete;
    KDNode&
    operator=(const KDNode&&) = delete;

    /*!
      @brief Default constructor
    */
    KDNode();

    /*!
      @brief Valid constructor. Takes list of primitives.
      @param[in] a_primitive List of primitives to partition later. 
      @note The input particle list is modifiable!
    */
    KDNode(const ParticleList<P>& a_primitives);

    /*!
      @brief Destructor. Does nothing.
    */
    virtual ~KDNode();

    /*!
      @brief Top-down partitioning building a KD-tree.
      @param[in] a_partitionFunction Partitioning function for dividing into subsets.
      @param[in] a_stopFunction      Stop function for ending particle recursion.
    */
    inline void
    topDownPartition(const ParticlePartitioner<P>& a_partitionFunction,
                     const StopFunction<P>&        a_stopFunction) noexcept;

    /*!
      @brief Get primitives in this node
    */
    inline const ParticleList<P>&
    getPrimitives() const noexcept;

    /*!
      @brief Is leaf node or not.
    */
    inline bool
    isLeafNode() const noexcept;

    /*!
      @brief Is leaf node or not.
    */
    inline bool
    isInteriorNode() const noexcept;

    /*!
      @brief Gather primitives further down in the subtree and return all primitives (in the leaf nodes)
    */
    inline ParticleList<P>
    gatherPrimitives() const noexcept;

    /*!
      @brief Move the primitives list further down in the subtree into this vector.
      @note This clears the primitives from the leaf nodes. Use gatherPrimitives if you do not
      want to change the tree contents. 
    */
    inline ParticleList<P>
    movePrimitives() noexcept;

  protected:
    /*!
      @brief Left KD-node
    */
    std::shared_ptr<KDNode> m_left;

    /*!
      @brief Right KD-node
    */
    std::shared_ptr<KDNode> m_right;

    /*!
      @brief List of primitives
    */
    std::vector<P> m_primitives;

    /*!
      @brief Gather primitives further down in the subtree and return all primitives (in the leaf nodes)
      @param[inout] a_primitives List of primitives
    */      
    inline void
    gatherPrimitives(ParticleList<P>& a_primitives) const noexcept;

    /*!
      @brief Move primitives further down in the subtree and return all primitives (in the leaf nodes)
      @param[inout] a_primitives List of primitives
      @note This clears the primitives from the leaf nodes. Use gatherPrimitives if you do not
      want to change the tree contents. 
    */      
    inline void
    movePrimitives(ParticleList<P>& a_primitives) noexcept;        
  };

  /*!
    @brief Class for describing a single node with two children and one parent. This is used together with SuperParticles::Tree for particle handling. 
  */
  template <class T>
  class Node : public std::enable_shared_from_this<Node<T>>
  {
  public:
    /*!
      @brief Default constructor.
    */
    Node();

    /*!
      @brief Constructor. Sets the parent. 
      @param[in] a_parent Parent node
    */
    Node(std::shared_ptr<Node<T>>& a_parent);

    /*!
      @brief Destructor. Destroys data in node
    */
    ~Node();

    /*!
      @brief Can split node or not. 
      @details This will check if there is more than one primitive in the node (in which case we can split) or if the single primitive in the node can be split. 
    */
    inline bool
    canSplit() const;

    /*!
      @brief Check if node is a leaf node
      @details This returns true if this node does not have children nodes. 
    */
    inline bool
    isLeafNode() const;

    /*!
      @brief Check if node is a regular node -- i.e. not a leaf node.
    */
    inline bool
    isRegularNode() const;

    /*!
      @brief Check if node is a root node. 
      @details Returns true if this node does not have a parent node. 
    */
    inline bool
    isRootNode() const;

    /*!
      @brief Set parent node
      @param[in] a_parent Pointer to parent node
    */
    inline void
    setParent(std::shared_ptr<Node<T>>& a_parent);

    /*!
      @brief Set pointer to left child node
      @param[in] a_left Left child node
    */
    inline void
    setLeft(std::shared_ptr<Node<T>>& a_left);

    /*!
      @brief Set pointer to left child node
      @param[in] a_right Right child node
    */
    inline void
    setRight(std::shared_ptr<Node<T>>& a_right);

    /*!
      @brief Set data. This sets the primitives in the node and the total mass in the node.
      @param[in] a_data Data
      @param[in] a_mass Mass of all primitives.
    */
    inline void
    setData(const std::vector<T>& a_data, const Real a_mass);

    /*!
      @brief Set the mass in the node
      @param[in] a_mass Mass of all primitives.
    */
    inline void
    setMass(const Real a_mass);

    /*!
      @brief Split node along the input coordinate.
      @details This will split the node along the input coordinate and make this node into a regular node. The left/right nodes become new leaf nodes.
      @param[in] a_dir Splitting coordinate -- this will compute a coordinate plane and put about half the mass on one side and half the mass on the other. 
    */
    inline void
    split(const int a_dir);

    /*!
      @brief Get the data in this node.
      @return Returns m_data
    */
    inline std::vector<T>&
    getData();

    /*!
      @brief Get the data in this node.
      @return Returns m_data
    */
    inline const std::vector<T>&
    getData() const;

    /*!
      @brief Get the total mass in this node
      @return m_mass
    */
    inline Real
    mass() const;

    /*!
      @brief Get this node's parent node
      @return m_parent
    */
    inline std::shared_ptr<Node<T>>&
    getParent();

    /*!
      @brief Get this node's left child node
      @return m_left
    */
    inline std::shared_ptr<Node<T>>&
    getLeft();

    /*!
      @brief Get this node's right child node
      @return m_right
    */
    inline std::shared_ptr<Node<T>>&
    getRight();

  protected:
    /*!
      @brief Pointer to parent node (can be nullptr if this node is the root node)
    */
    std::shared_ptr<Node<T>> m_parent;

    /*!
      @brief Pointer to left child node (can be nullptr if this node is a leaf node)
    */
    std::shared_ptr<Node<T>> m_left;

    /*!
      @brief Pointer to left child node (can be nullptr if this node is a leaf node)
    */
    std::shared_ptr<Node<T>> m_right;

    /*!
      @brief Total mass in this node
    */
    Real m_mass;

    /*!
      @brief List of primitives contained in this node
    */
    std::vector<T> m_data;
  };

  /*!
    @brief Declaration of a node partitioner function. 
    @param[inout] a_node Node to be split
    @param[in]    a_dir  Coordinate to split along
  */
  template <class T>
  using NodePartitioner = std::function<void(std::shared_ptr<Node<T>>& a_node, const int a_dir)>;

  /*!
    @brief Partitioner function which splits a Node<T> into left/right child nodes with masses differing by at most one physical particle
    @param[inout] a_node Node to be split
    @param[in]    a_dir  Coordinate to split along
    @note Implemented in CD_SuperParticlesImplem.H
  */

  /*!
    @brief Class for making a tree structure with SuperParticles::Node<T> in it. 
    @details This is the main function that is used for partitioning particles ala ItoSolver -- it uses the SuperParticles::Node<T> structure to recursively 
    partition a list of primitives into a new list of primitives. 
    @note This class should work with any Node<T> -- refer to SuperParticles::Node<T> for requirements on the template parameters. To build the tree the user
    will just call buildTree() which selects the first splitting plane and number of leaves in the tree. The user can also input his own partitioning function
    which determines how leaf nodes are split. 
  */
  template <class T>
  class Tree
  {
  public:
    /*!
      @brief Weak constructor. Must subsequently call the define function
    */
    Tree();

    /*!
      @brief Full constructor. Takes a list of primitives to be merged/split.
      @param[in] a_point List of primitives that will be partitioned by the tree
    */
    Tree(const std::vector<T>& a_points);

    /*!
      @brief Destructor (does nothing)
    */
    virtual ~Tree();

    /*!
      @brief Full define function. Takes a list of primitives to be merged/split.
      @param[in] a_point List of primitives that will be partitioned by the tree
    */
    inline void
    define(const std::vector<T>& a_points);

    /*!
      @brief Build tree with desired number of leaves. This will recursively partition the primitives.
      @param[in] a_firstDir        Coordinate of the first splitting direction (i.e. the first plane that we split along). Subsequent planes are cycled x->y->z->x->y->z etc.
      @param[in] a_numTargetLeaves Desired number of leaves when building the tree.
      @param[in] a_partitioner     Partitioning function. 
    */
    inline void
    buildTree(const int a_firstDir, const int a_numTargetLeaves, const NodePartitioner<T>& a_partitioner);

    /*!
      @brief Get the leaf nodes
      @return Returns m_leaves
    */
    inline std::vector<std::shared_ptr<Node<T>>>&
    getLeaves();

    /*!
      @brief Get the leaf nodes
      @return Returns m_leaves
    */
    inline const std::vector<std::shared_ptr<Node<T>>>&
    getLeaves() const;

  protected:
    /*!
      @brief Is defined or not
    */
    bool m_isDefined;

    /*!
      @brief Pointer to root node. 
    */
    std::shared_ptr<Node<T>> m_root;

    /*!
      @brief All leaf nodes. Yes, I know this is O(log N) lookup but why not just keep a list since we ALWAYS need to return it. 
    */
    std::vector<std::shared_ptr<Node<T>>> m_leaves;
  };
} // namespace SuperParticles

#include <CD_NamespaceFooter.H>

#include <CD_SuperParticlesImplem.H>
#include <CD_ParticlePartitioners.H>

#endif
