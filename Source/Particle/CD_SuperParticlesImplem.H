/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_SuperParticlesImplem.H
  @brief  Implementation of CD_SuperParticles.H
  @author Robert Marskar
*/

#ifndef CD_SuperParticlesImplem_H
#define CD_SuperParticlesImplem_H

// Std includes
#include <utility>

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_SuperParticles.H>
#include <CD_NamespaceHeader.H>

#define SuperParticlesVectorCapacity 2048

namespace SuperParticles {

  template <class P>
  inline KDNode<P>::KDNode() : m_left(nullptr), m_right(nullptr), m_weight(0.0)
  {}

  template <class P>
  inline KDNode<P>::KDNode(ParticleList& a_particles)
    : m_left(nullptr), m_right(nullptr), m_weight(0.0), m_particles(std::move(a_particles))
  {}

  template <class P>
  inline KDNode<P>::~KDNode()
  {}

  template <class P>
  inline const typename KDNode<P>::ParticleList&
  KDNode<P>::getParticles() const noexcept
  {
    return m_particles;
  }

  template <class P>
  inline typename KDNode<P>::ParticleList&
  KDNode<P>::getParticles() noexcept
  {
    return m_particles;
  }

  template <class P>
  inline const Real&
  KDNode<P>::weight() const noexcept
  {
    return m_weight;
  }

  template <class P>
  inline Real&
  KDNode<P>::weight() noexcept
  {
    return m_weight;
  }

  template <class P>
  inline bool
  KDNode<P>::isLeafNode() const noexcept
  {
    return (m_left == nullptr) && (m_right == nullptr);
  }

  template <class P>
  inline bool
  KDNode<P>::isInteriorNode() const noexcept
  {
    return !(this->isLeafNode());
  }

  template <class P>
  inline typename KDNode<P>::ParticleList
  KDNode<P>::gatherParticles() const noexcept
  {
    ParticleList primitives;

    this->gatherParticles(primitives);

    return primitives;
  }

  template <class P>
  inline void
  KDNode<P>::gatherParticles(ParticleList& a_particles) const noexcept
  {
    if (this->isLeafNode()) {
      a_particles.reserve(a_particles.size() + m_particles.size());
      a_particles.insert(a_particles.end(), m_particles.begin(), m_particles.end());
    }
    else {
      m_left->gatherParticles(a_particles);
      m_right->gatherParticles(a_particles);
    }
  }

  template <class P>
  inline typename KDNode<P>::ParticleList
  KDNode<P>::moveParticles() noexcept
  {
    ParticleList primitives;

    this->moveParticles(primitives);

    return primitives;
  }

  template <class P>
  inline void
  KDNode<P>::moveParticles(ParticleList& a_particles) noexcept
  {
    if (this->isLeafNode()) {
      a_particles.reserve(a_particles.size() + m_particles.size());

      std::move(m_particles.begin(), m_particles.end(), std::back_inserter(a_particles));
    }
    else {
      m_left->moveParticles(a_particles);
      m_right->moveParticles(a_particles);
    }
  }

  template <class P>
  inline std::shared_ptr<KDNode<P>>&
  KDNode<P>::getLeft() noexcept
  {
    return m_left;
  }

  template <class P>
  inline std::shared_ptr<KDNode<P>>&
  KDNode<P>::getRight() noexcept
  {
    return m_right;
  }

  template <class P, Real& (P::*weight)(), const RealVect& (P::*position)() const>
  typename KDNode<P>::Partitioner PartitionEqualWeight = [](KDNode<P>& a_node) -> void {
    CH_assert(!(a_node.isInteriorNode()));

    constexpr Real splitThresh = 2.0 - std::numeric_limits<Real>::min();

    // Particles and node weight.
    typename KDNode<P>::ParticleList& particles = a_node.getParticles();

    const Real W = a_node.weight();

    // A. Figure out which coordinate direction we should partition and sort
    //    the particles.
    RealVect loCorner = +std::numeric_limits<Real>::max() * RealVect::Unit;
    RealVect hiCorner = -std::numeric_limits<Real>::max() * RealVect::Unit;

    for (auto& p : particles) {
      const RealVect& pos = (p.*position)();
      for (int dir = 0; dir < SpaceDim; dir++) {
        loCorner[dir] = std::min(pos[dir], loCorner[dir]);
        hiCorner[dir] = std::max(pos[dir], hiCorner[dir]);
      }
    }

    const int splitDir = (hiCorner - loCorner).maxDir(true);

    auto sortCrit = [splitDir](const P& p1, const P& p2) -> bool {
      return (p1.*position)()[splitDir] < (p2.*position)()[splitDir];
    };

    std::sort(particles.begin(), particles.end(), sortCrit);

    // B. Determine the "median particle" and start computing the weight in the
    //    two halves.
    size_t id = 0;
    Real   wl = 0.0;
    Real   wr = W - (particles[id].*weight)();

    for (size_t i = 1; i < particles.size(); i++) {
      const Real& w = (particles[id].*weight)();

      if (wl + w < wr) {
        id = i;
        wl += w;
        wr = W - wl - (particles[id].*weight)();
      }
      else {
        break;
      }
    }

    // C. Copy the two particle halves to each subnode.
    P p = particles[id];

    typename KDNode<P>::ParticleList pl;
    typename KDNode<P>::ParticleList pr;

    std::move(particles.begin(), particles.begin() + id, std::back_inserter(pl));
    std::move(particles.begin() + id + 1, particles.end(), std::back_inserter(pr));

    const Real& pw = (p.*weight)();
    const Real  dw = wr - wl;

    CH_assert(wl + wr + pw == W);

    // D. Assign the median particle; split the particle if we can.
    if (pw >= splitThresh && pw >= std::abs(dw)) {
      Real dwl = dw;
      Real dwr = 0.0;
      Real ddw = pw - dw;

      const long long N = llround(ddw);

      if (N > 0LL) {

        const long long Nr = N / 2;
        const long long Nl = N - Nr;

        dwl += (ddw / N) * Nl;
        dwr += (ddw / N) * Nr;
      }

      if (dwl > 0.0) {
        P il(p);

        CH_assert(dwl >= 1.0);

        wl += dwl;
        (il.*weight)() = dwl;
        pl.emplace_back(std::move(il));
      }

      if (dwr > 0.0) {
        P ir(p);

        CH_assert(dwr >= 1.0);

        wr += dwr;
        (ir.*weight)() = dwr;
        pr.emplace_back(std::move(ir));
      }
    }
    else {
      if (wl <= wr) {
        wl += pw;
        pl.emplace_back(std::move(p));
      }
      else {
        wr += pw;
        pr.emplace_back(std::move(p));
      }
    }

    // E. If this breaks, weight is not conserved or we broke the median particle splitting; the weight difference
    //    between the left/right node should be at most one physical particle.
    CH_assert(wl + wr == W);
    CH_assert(std::abs(wl - wr) <= 1.0);

    // F. Instantiate the child nodes.
    particles.resize(0);

    a_node.getLeft()  = std::make_shared<KDNode<P>>(pl);
    a_node.getRight() = std::make_shared<KDNode<P>>(pr);

    a_node.getLeft()->weight()  = wl;
    a_node.getRight()->weight() = wr;

    // G. Debug code; make sure particle weights make sense
#ifndef NDEBUG
    Real WL = 0.0;
    Real WR = 0.0;

    // Note: Not auto& l/r : pl/pr because the particles were into the child node.
    for (auto& l : a_node.getLeft()->getParticles()) {
      CH_assert((l.*weight)() >= 1.0);

      WL += (l.*weight)();
    }

    for (auto& r : a_node.getRight()->getParticles()) {
      CH_assert((r.*weight)() >= 1.0);

      WR += (r.*weight)();
    }

    CH_assert(WL == wl);
    CH_assert(WR == wr);
#endif
  };

  template <class P, const Real& (P::*weight)() const>
  typename KDNode<P>::StopFunction InsufficientWeightForSplitting =
    [](const typename KDNode<P>::ParticleList& a_particles) {
      constexpr Real splitThresh = 2.0 - std::numeric_limits<Real>::min();

      Real W = 0.0;
      for (const auto& p : a_particles) {
        W += (p.*weight)();
      }

      return (W < splitThresh) ? true : false;
    };

  template <class T>
  Node<T>::Node()
  {
    m_parent = nullptr;
    m_left   = nullptr;
    m_right  = nullptr;

    m_data.reserve(SuperParticlesVectorCapacity);
  }

  template <class T>
  Node<T>::~Node()
  {
    m_data.resize(0);
  }

  template <class T>
  Node<T>::Node(std::shared_ptr<Node<T>>& a_parent) : Node<T>()
  {
    this->setParent(a_parent);
  }

  template <class T>
  inline bool
  Node<T>::isLeafNode() const
  {
    return m_left == nullptr && m_right == nullptr;
  }

  template <class T>
  inline bool
  Node<T>::isRegularNode() const
  {
    return !(this->isLeafNode());
  }

  template <class T>
  inline bool
  Node<T>::isRootNode() const
  {
    return m_parent == nullptr;
  }

  template <class T>
  inline void
  Node<T>::setParent(std::shared_ptr<Node<T>>& a_parent)
  {
    m_parent = a_parent;
  }

  template <class T>
  inline void
  Node<T>::setLeft(std::shared_ptr<Node<T>>& a_left)
  {
    m_left = a_left;
  }

  template <class T>
  inline void
  Node<T>::setRight(std::shared_ptr<Node<T>>& a_right)
  {
    m_right = a_right;
  }

  template <class T>
  inline void
  Node<T>::setData(const std::vector<T>& a_data, const Real a_mass)
  {
    m_data = a_data;
    m_mass = a_mass;
  }

  template <class T>
  inline void
  Node<T>::setMass(const Real a_mass)
  {
    m_mass = a_mass;
  }

  template <class T>
  inline std::shared_ptr<Node<T>>&
  Node<T>::getLeft()
  {
    return m_left;
  }

  template <class T>
  inline std::shared_ptr<Node<T>>&
  Node<T>::getRight()
  {
    return m_right;
  }

  template <class T>
  inline std::shared_ptr<Node<T>>&
  Node<T>::getParent()
  {
    return m_parent;
  }

  template <class T>
  inline std::vector<T>&
  Node<T>::getData()
  {
    return m_data;
  }

  template <class T>
  inline const std::vector<T>&
  Node<T>::getData() const
  {
    return m_data;
  }

  template <class T>
  inline bool
  Node<T>::canSplit() const
  {
    bool ret = false;
    if (m_data.size() > 0) {
      if (m_data.size() > 1 || m_data[0].canSplit()) {
        ret = true;
      }
    }

    return m_mass >= 2.0;
    return ret;
  }

  template <class T>
  inline Real
  Node<T>::mass() const
  {
    return m_mass;
  }

  template <class T>
  Tree<T>::Tree()
  {
    m_isDefined = false;
  }

  template <class T>
  Tree<T>::Tree(const std::vector<T>& a_data) : Tree()
  {
    this->define(a_data);
  }

  template <class T>
  Tree<T>::~Tree()
  {
    m_leaves.resize(0);
  }

  template <class T>
  inline void
  Tree<T>::define(const std::vector<T>& a_data)
  {

    // Make the root node if we don't have it.
    if (!m_root) {
      m_root = std::make_shared<Node<T>>();
    }

    // Compute total mass.
    Real mass = 0.0;
    for (const auto& d : a_data) {
      mass += d.mass();
    }

    CH_assert(mass >= 1.0);

    m_root->setData(a_data, mass);

    m_leaves.resize(1);
    m_leaves[0] = m_root;

    m_isDefined = true;
  }

  template <class T>
  inline std::vector<std::shared_ptr<Node<T>>>&
  Tree<T>::getLeaves()
  {
    CH_assert(m_isDefined);

    return m_leaves;
  }

  template <class T>
  inline const std::vector<std::shared_ptr<Node<T>>>&
  Tree<T>::getLeaves() const
  {
    CH_assert(m_isDefined);

    return m_leaves;
  }

  template <class T>
  inline void
  Tree<T>::buildTree(const int a_firstDir, const int a_numTargetLeaves, const NodePartitioner<T>& a_partitioner)
  {
    CH_TIME("SuperParticles::Tree<T>::buildTree");

    CH_assert(m_isDefined);
    CH_assert(a_numTargetLeaves >= 1);
    CH_assert(a_firstDir >= 0);
    CH_assert(a_firstDir <= SpaceDim - 1);

    // Should evaluate to an empty loop and be purged by the compiler if DEBUG=FALSE
    for (const auto& d : m_root->getData()) {
      CH_assert(d.mass() >= 1.0);
    }

    // TLDR: This function splits the primitives recursively until we either have a_numTargetLeaves leaf nodes OR none of the leaf nodes can be
    //       partitioned further.

    // Reset leaves and root. I'm just resizing here so I get the capacity
    // and avoid vector resizing during the splitting process.
    std::vector<std::shared_ptr<Node<T>>> newLeaves;
    m_leaves.resize(a_numTargetLeaves, nullptr);
    newLeaves.resize(a_numTargetLeaves, nullptr);
    m_leaves[0] = m_root;

    // We always begin with one leaf node.
    int numCurrentLeaves = 1;

    int  splitDir      = a_firstDir;
    bool keepSplitting = (a_numTargetLeaves > 1);

    while (keepSplitting) {

      // Flag which determines if we can split any of the nodes further. This will remain false if none of the
      // leaf nodes could be split OR if we didn't need more leaf nodes (i.e. the tree has finished building).
      bool canSplit = false;

#if 0 // I'm not sure this matters so I'm commenting it out. R.M. Nov. 2021. \
      // Sort the leaves by their mass(?). We will split the ones with the largest masses.      
      std::sort(m_leaves.begin(), m_leaves.begin() + numCurrentLeaves,
		[](const std::shared_ptr<Node<T> >& n1,
		   const std::shared_ptr<Node<T> >& n2){
		  return n1->mass() < n2->mass();
		});
#endif

      // These are the number of leaves that we need in order to reach a_numTargetLeaves.
      int leavesNeeded = a_numTargetLeaves - numCurrentLeaves;

      int N = 0;
      for (int i = 0; i < numCurrentLeaves; i++) {

        // If we needed a new leaf AND we could split this node, split it and go further. This will give us two leaf nodes
        // in the "new" list of leaf nodes.
        if (m_leaves[i]->canSplit() && leavesNeeded > 0) {
          a_partitioner(m_leaves[i], splitDir);

          newLeaves[N]     = m_leaves[i]->getLeft();
          newLeaves[N + 1] = m_leaves[i]->getRight();

          N += 2;

          leavesNeeded--;

          canSplit = true;
        }
        else { // We did not split the leaf node so just append the current leaf to the list of "new" leaves.
          newLeaves[N] = m_leaves[i];

          N++;
        }
      }

      // Assign new leaf nodes.
      m_leaves.assign(newLeaves.begin(), newLeaves.begin() + N);

      // Update the current number of leaf nodes.
      numCurrentLeaves = N;

      // Did the tree finish building or are all the nodes unsplittable leaf nodes? If so, exit.
      keepSplitting = (N < a_numTargetLeaves) && canSplit;

      // Update the splitting coordinate -- the next leaves are split along a different coordinate.
      splitDir = (splitDir + 1) % SpaceDim;
    }
  }

  /*!
    @brief Partitioning function which partitions masses equally in the kd-tree. 
  */
  template <class T>
  NodePartitioner<T> NodePartitionEqualMass = [](std::shared_ptr<Node<T>>& a_node, const int a_dir) {
    CH_TIME("NodePartitionEqualMass<T>");

    CH_assert(a_node != nullptr);
    CH_assert(a_dir >= 0);
    CH_assert(a_dir < SpaceDim);

    // Get a handle to the node data.
    std::vector<T>& data = a_node->getData();
    const Real      mass = a_node->mass();

    // Some basic assertions -- the loop should evaluate to an empty loop when DEBUG=TRUE.
    CH_assert(data.size() > 0);
    CH_assert(mass > 1.0);

    for (const auto& n : data) {
      CH_assert(n.mass() >= 1.0);
    }

    // TLDR: This function splits the primitives in the node into two subvolumes. A splitting plane is made along the input coordinate and it is positioned
    //       such that the total mass on the left/right side of the splitting plane differ by at most one computational particle. The guts of the code mainly
    //       consists of tests for how to divide the "median" particle -- there are various cases that we need to cover here (such as the median particle being
    //       a physical particle or a superparticle).

    // 1. Sort data along the input coordinate. We need this because we need to have something sensible to iterate through
    //    when we place the splitting plane.
    std::sort(std::begin(data), std::end(data), [a_dir](const T& p1, const T& p2) { return p1[a_dir] < p2[a_dir]; });

    // 2. Figure out where to place the splitting plane -- this is not the physical coordinate but index in data (the list of primitives). We are
    //    actually looking for the particle on some median where the masses between the two half-planes can be made approximately the same. We begin
    //    by setting the median on i = 0 in the primitives vector (data). We then go through the primitives and check if adding the median particle
    //    ends up placing MORE mass in the left node than the right node. Once that happens we have the inflection point where we need to split. The median
    //    particle will be assigned later because it may be a superparticle and we can equilibriate the masses even further.
    //
    //
    // In all of the below. massLeft and massRight are the masses on each side of the median particle. As remarked above, the median particle could be split
    // later.
    int  medianParticleIndex = 0;
    Real massLeft            = 0.0;                           // Zero mass, as it should be.
    Real massRight = mass - data[medianParticleIndex].mass(); // So, entire mass minus the median particle mass.

    for (int i = 1; i < data.size(); i++) {
      const Real& medianParticleMass = data[medianParticleIndex].mass(); // This is the current median particle mass.

      // If adding the median particle to the mass in the left node still yields a lower mass than the right now, we just shift the median particle index. So, the previous
      // "median particle" goes to the left node but we then have a new median particle.
      if (massLeft + medianParticleMass < massRight) {
        massLeft += medianParticleMass;

        medianParticleIndex = i;

        massRight = mass - massLeft -
                    data[medianParticleIndex].mass(); // Right mass = Total mass - left mass - median particle mass
      }
      else {
        break;
      }
    }

    // 2. Make the two leaves and push data into them. Currently, we know that all particles to the left of the "median" particle go into the left half-plane
    //    and all particles to the right of the "median" particle go into the right half-plane. We do not yet know how to place the median particle -- it could
    //    go into either the left node or the right node, or it could be split into and placed into either node.
    std::shared_ptr<Node<T>>& left  = a_node->getLeft();
    std::shared_ptr<Node<T>>& right = a_node->getRight();

    if (!left)
      left = std::make_shared<Node<T>>();
    if (!right)
      right = std::make_shared<Node<T>>();

    std::vector<T>& leftData  = left->getData();
    std::vector<T>& rightData = right->getData();

    leftData.assign(data.begin(), data.begin() + medianParticleIndex);
    rightData.assign(data.begin() + medianParticleIndex + 1, data.end());

    // 3. We next need to assign the median particle. If we can split up the median point into several ones, we do it. This would be the case if the median point
    //    is a superparticle; in this case we divide the weights as best we can and then assign them to the left/right halves. If we can't then we just assign the
    //    particle to whichever left/right child node has the lowest mass.
    const T& medianParticle = data[medianParticleIndex];

    if (medianParticle.canSplit()) {
      const RealVect& medianPos  = medianParticle.pos();  // Physical position
      const Real&     medianMass = medianParticle.mass(); // Always positive
      const Real&     medianEnergy =
        medianParticle.energy(); // Note: This is average energy so it remain the same for subparticles.
      const Real massDiff = massRight - massLeft; // Not necessarily positive.

      // This hook is for when all mass goes into one part.
      if (medianMass < std::abs(massDiff)) {
        if (massDiff > 0) { // Right mass is bigger than left so put entire particle in left mass.
          massLeft += medianMass;
          leftData.emplace_back(T(medianPos, medianMass, medianEnergy));
        }
        else {                     // Left mass is bigger than right so put everything into the right node.
          massRight += medianMass; // Goes into right
          rightData.emplace_back(T(medianPos, medianMass, medianEnergy));
        }
      }
      else {
        // Distribute mass as best we can. This might seem weird because we can have massLeft > massRight in which
        // case massDiff < 0. This still works since we are simply subtracting the difference from the left leaf first and
        // then divide a larger total mass between the two leaves. Essentially, we are subtracting the part of the median
        // particle that gives massLeft > massRight from the left particle first, and then we reassign a larger total mass.
        Real tmpMassLeft  = massDiff;
        Real tmpMassRight = 0.0;
        Real tmpMassRem   = medianMass - massDiff;

        const long long N = llround(tmpMassRem);
        if (N > 0LL) {
          const long long NR = N / 2;
          const long long NL = N - NR;

          tmpMassLeft += (tmpMassRem / N) * NL;
          tmpMassRight += (tmpMassRem / N) * NR;
        }

        if (tmpMassLeft > 0.0) {
          leftData.emplace_back(T(medianPos, tmpMassLeft, medianEnergy));
          massLeft += tmpMassLeft;
        }

        if (tmpMassRight > 0.0) {
          rightData.emplace_back(T(medianPos, tmpMassRight, medianEnergy));
          massRight += tmpMassRight;
        }
      }
    }
    else {
      // Not a computational particle -- put the particle in the child node with the lowest mass.
      if (massLeft <= massRight) {
        leftData.emplace_back(medianParticle);
        massLeft += medianParticle.mass();
      }
      else {
        rightData.emplace_back(medianParticle);
        massRight += medianParticle.mass();
      }
    }

    // Update the tree structure.
    left->setMass(massLeft);
    right->setMass(massRight);

    left->setParent(a_node);
    right->setParent(a_node);

    data.resize(0);

    // Some assertions in case this code every breaks.
    CH_assert((std::abs(mass - (massLeft + massRight)) <= 1.E-12));
    CH_assert(massLeft >= std::numeric_limits<Real>::min());
    CH_assert(massRight >= std::numeric_limits<Real>::min());

    // These should be evaluated to empty loopps by the compiler if DEBUG=FALSE.
    for (const auto& l : leftData) {
      CH_assert(l.mass() >= std::numeric_limits<Real>::min());
    }

    for (const auto& r : rightData) {
      CH_assert(r.mass() >= std::numeric_limits<Real>::min());
    }
  };

} // namespace SuperParticles

#include <CD_NamespaceFooter.H>

#endif
